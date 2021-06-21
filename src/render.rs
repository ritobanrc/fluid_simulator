use egui_wgpu_backend::RenderPass;
use futures::executor::block_on;
use std::sync::mpsc::{Receiver, RecvError};
use std::time::Instant;
use winit::{
    dpi::PhysicalSize,
    event::*,
    event_loop::{ControlFlow, EventLoop},
    window::WindowBuilder,
};

mod camera;
mod scene;
mod state;

pub use camera::{Camera, CameraController};
pub use scene::{Scene, Vertex};
pub use state::State;

pub fn render_texture(
    image_dir: std::path::PathBuf,
    rx: Receiver<Vec<Vertex>>,
    width: u32,
    height: u32,
    num_frames: usize,
) -> Result<(), RecvError> {
    let mut state = block_on(State::with_texture(width, height));

    let start = Instant::now();

    let verticies = rx.recv()?;
    let mut scene = Scene::new(&verticies, &state.device, (width, height));

    for i in 0..num_frames {
        let verticies = rx.recv().expect("Failed to recieve verts");
        state.update(&mut scene, &verticies);
        state.render(&scene, None).expect("Swapchain error");

        block_on(async {
            let buffer_slice = state.get_buffer().unwrap().slice(..);

            // NOTE: We have to create the mapping THEN device.poll() before await
            // the future. Otherwise the application will freeze.
            let mapping = buffer_slice.map_async(wgpu::MapMode::Read);
            state.device.poll(wgpu::Maintain::Wait);
            mapping.await.unwrap();

            let data = buffer_slice.get_mapped_range();

            use image::{ImageBuffer, Rgba};
            let buffer = ImageBuffer::<Rgba<u8>, _>::from_raw(width, height, data).unwrap();

            let mut path = image_dir.clone();
            path.push(format!("frame{:04}.png", i));
            buffer.save(&path).unwrap();

            if (i + 1) % 20 == 0 {
                let end = Instant::now();
                println!("Completed frame {}. Time elapsed: {:?}", i + 1, end - start);
            }
        });
        state.get_buffer().unwrap().unmap();
    }

    let end = Instant::now();
    println!(
        "Completed simulating {} frames. Time elapsed: {:?}",
        num_frames,
        end - start
    );

    Ok(())
}

/// A custom event type for the winit app.
enum CustomEvent {
    RequestRedraw,
}

/// This is the repaint signal type that egui needs for requesting a repaint from another thread.
/// It sends the custom RequestRedraw event to the winit event loop.
struct RepaintSignal(std::sync::Mutex<winit::event_loop::EventLoopProxy<CustomEvent>>);

impl epi::RepaintSignal for RepaintSignal {
    fn request_repaint(&self) {
        self.0
            .lock()
            .unwrap()
            .send_event(CustomEvent::RequestRedraw)
            .ok();
    }
}

pub fn open_window(rx: Receiver<Vec<Vertex>>) -> Result<(), RecvError> {
    let event_loop = EventLoop::with_user_event();
    let size = PhysicalSize::new(1280, 720);
    let window = WindowBuilder::new()
        .with_inner_size(size)
        .with_title("Fluid Simulation")
        .build(&event_loop)
        .unwrap();

    let mut state = block_on(State::from_window(&window));

    // -------------- EGUI Stuff -----------------
    //
    //
    let mut platform =
        egui_winit_platform::Platform::new(egui_winit_platform::PlatformDescriptor {
            physical_width: size.width,
            physical_height: size.height,
            scale_factor: window.scale_factor(),
            font_definitions: egui::FontDefinitions::default(),
            style: egui::Style::default(),
        });
    let mut egui_render_pass = RenderPass::new(&state.device, wgpu::TextureFormat::Bgra8UnormSrgb);

    let repaint_signal = std::sync::Arc::new(RepaintSignal(std::sync::Mutex::new(
        event_loop.create_proxy(),
    )));

    let mut ui_state = UIState::default();
    // -------------------------------------------

    let size = window.inner_size();
    let verticies = rx.recv()?;
    let mut scene = Scene::new(&verticies, &state.device, (size.width, size.height));

    let start_time = std::time::Instant::now();
    let mut last_render_time = start_time;
    let mut frame_count = 0;

    event_loop.run(move |event, _, control_flow| {
        platform.handle_event(&event);

        match event {
            Event::WindowEvent {
                ref event,
                window_id,
            } if window_id == window.id() => match event {
                WindowEvent::CloseRequested => *control_flow = ControlFlow::Exit,
                WindowEvent::KeyboardInput { input, .. } => match input {
                    KeyboardInput {
                        state: ElementState::Pressed,
                        virtual_keycode: Some(VirtualKeyCode::Escape),
                        ..
                    } => *control_flow = ControlFlow::Exit,
                    _ => {
                        state.input(event, &mut scene);
                    }
                },
                WindowEvent::Resized(physical_size) => {
                    state.resize(*physical_size);
                }
                WindowEvent::ScaleFactorChanged { new_inner_size, .. } => {
                    // new_inner_size is &&mut so we have to dereference it twice
                    state.resize(**new_inner_size);
                }
                window_event => {
                    // TODO: figure out what to do with this bool
                    //       basically, make a better input system in general
                    state.input(window_event, &mut scene);
                }
            },
            Event::RedrawRequested(_) => {
                let now = std::time::Instant::now();
                let dt = now - last_render_time;
                last_render_time = now;
                frame_count += 1;
                if frame_count % 100 == 0 {
                    println!("Running at {:?} FPS", 1_000_000 / dt.as_micros())
                }

                platform.update_time(start_time.elapsed().as_secs_f64());
                // Get swapchain output frame
                //let egui_start = Instant::now();
                platform.begin_frame();
                let mut app_output = epi::backend::AppOutput::default();

                let mut frame = epi::backend::FrameBuilder {
                    info: epi::IntegrationInfo {
                        web_info: None,
                        cpu_usage: None,
                        seconds_since_midnight: None,
                        native_pixels_per_point: Some(window.scale_factor() as _),
                    },
                    tex_allocator: &mut egui_render_pass,
                    output: &mut app_output,
                    repaint_signal: repaint_signal.clone(),
                }
                .build();

                ui_state.egui_update(platform.context(), &mut frame);

                let (_output, paint_commands) = platform.end_frame();
                let paint_jobs = platform.context().tessellate(paint_commands);

                let egui_state = EguiRenderState {
                    platform: &platform,
                    render_pass: &mut egui_render_pass,
                    paint_jobs: &paint_jobs,
                    screen_descriptor: egui_wgpu_backend::ScreenDescriptor {
                        physical_width: size.width,
                        physical_height: size.height,
                        scale_factor: window.scale_factor() as f32,
                    },
                };

                let verticies = rx.recv().expect("Failed to recieve verts");
                state.update(&mut scene, &verticies);
                match state.render(&scene, Some(egui_state)) {
                    Ok(_) => {}
                    // Recreate the swap_chain if lost
                    Err(wgpu::SwapChainError::Lost) => state.resize(state.size),
                    // The system is out of memory, we should probably quit
                    Err(wgpu::SwapChainError::OutOfMemory) => *control_flow = ControlFlow::Exit,
                    // All other errors (Outdated, Timeout) should be resolved by the next frame
                    Err(e) => eprintln!("{:?}", e),
                }

                *control_flow = ControlFlow::Poll;
            }
            Event::MainEventsCleared | Event::UserEvent(CustomEvent::RequestRedraw) => {
                // RedrawRequested will only trigger once, unless we manually
                // request it.
                window.request_redraw();
            }
            _ => {}
        }
    });
}

#[derive(Default)]
pub struct UIState {
    counter: i32,
}

impl UIState {
    fn egui_update(&mut self, ctx: egui::CtxRef, _frame: &mut epi::Frame) {
        egui::Window::new("My Window").show(&ctx, |ui| {
            ui.label("Hello World!");
            //
            // Put the buttons and label on the same row:
            ui.horizontal(|ui| {
                if ui.button("-").clicked() {
                    self.counter -= 1;
                }
                ui.label(self.counter.to_string());
                if ui.button("+").clicked() {
                    self.counter += 1;
                }
            });
        });
    }
}

/// The Egui-associated state that `State::render` needs
pub struct EguiRenderState<'a> {
    platform: &'a egui_winit_platform::Platform,
    render_pass: &'a mut egui_wgpu_backend::RenderPass,
    paint_jobs: &'a [egui::ClippedMesh],
    screen_descriptor: egui_wgpu_backend::ScreenDescriptor,
}
