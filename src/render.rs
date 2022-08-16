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
mod texture;
mod ui;

pub use camera::{Camera, CameraController};
pub use scene::{Scene, Vertex};
pub use state::State;
pub use ui::{EguiRenderState, UIState};

use crate::Simulation;

pub fn render_texture(
    image_dir: std::path::PathBuf,
    rx: Receiver<Vec<Vertex>>,
    width: u32,
    height: u32,
    num_frames: usize,
) -> Result<(), RecvError> {
    let mut state = block_on(State::with_texture(width, height));

    let start = Instant::now();

    let vertices = rx.recv()?;
    let mut scene = Scene::new(&vertices, &state.device, (width, height));

    for i in 0..num_frames {
        let vertices = rx.recv().expect("Failed to recieve verts");
        state.update(&mut scene, 20., &vertices);
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

pub fn open_window(files: Option<Vec<std::path::PathBuf>>) -> eyre::Result<()> {
    let event_loop = EventLoop::with_user_event();
    let mut size = PhysicalSize::new(1280, 720);
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

    let verticies = Vec::new();
    let mut scene = Scene::new(&verticies, &state.device, (size.width, size.height));

    let start_time = std::time::Instant::now();
    let mut last_render_time = start_time;

    let mut rx = None;

    let mut stop_tx = None;

    let mut frame_number = 0;

    event_loop.run(move |event, _, control_flow| {
        platform.handle_event(&event);

        match event {
            Event::DeviceEvent { event, .. } => {
                state.input(&event, &mut scene);
            }
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
                    _ => {}
                },
                WindowEvent::Resized(physical_size) => {
                    size = *physical_size;
                    state.resize(*physical_size);
                }
                WindowEvent::ScaleFactorChanged { new_inner_size, .. } => {
                    size = **new_inner_size;
                    state.resize(**new_inner_size);
                }
                _ => {}
            },
            Event::RedrawRequested(_) => {
                let now = std::time::Instant::now();
                let dt = now - last_render_time;
                last_render_time = now;

                ui_state.render_dt = dt;

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

                let mut should_start_simulation = false;
                ui_state.egui_update(
                    platform.context(),
                    &mut frame,
                    &mut should_start_simulation,
                    &stop_tx,
                );

                if should_start_simulation && rx.is_none() {
                    let (new_stop_tx, new_stop_rx) = std::sync::mpsc::channel();
                    stop_tx = Some(new_stop_tx);
                    rx = Some(start_simulation(&ui_state, new_stop_rx, None));
                }

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

                let vertices = if let Some(files) = &files {
                    let file = &files[frame_number];
                    let expected = format!("{:03}.dat", frame_number);
                    let name = file.file_name().unwrap();
                    if name != &expected[..] {
                        panic!(
                            "File name not expected. Expected {:?}, found {:?}",
                            expected, name
                        );
                    }
                    let verts: Vec<Vertex> = rmp_serde::decode::from_read(
                        std::fs::File::open(file)
                            .expect(&format!("failed to open file: {:?}", file)),
                    )
                    .expect(&format!("failed to parse file: {:?}", file));
                    verts
                } else if let Some(this_rx) = &rx {
                    match this_rx.recv() {
                        Ok(verts) => verts,
                        Err(_) => {
                            stop_tx = None;
                            rx = None;
                            Vec::new()
                        }
                    }
                } else {
                    use crate::initial_condition::InitialCondition;
                    let mut nop = NopSimulation::default();
                    ui_state.initial_condition.add_particles(&mut nop);

                    nop.vertices
                };

                state.update(&mut scene, ui_state.particle_size, &vertices);

                match state.render(&scene, Some(egui_state)) {
                    Ok(_) => {
                        frame_number += 1;
                    }
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
/// An empty simulation only used for displaying the particles' initial condition before the main
/// simulation starts
struct NopSimulation {
    vertices: Vec<Vertex>,
}

impl Simulation for NopSimulation {
    type Parameters = ();

    fn new(_params: Self::Parameters) -> Self {
        NopSimulation {
            vertices: Vec::new(),
        }
    }

    fn simulate_frame(&mut self) -> Vec<Vertex> {
        panic!("NopSimulation::simulate_frame should never be called");
    }

    fn add_particle(&mut self, _mass: crate::Scalar, position: crate::Vec3, velocity: crate::Vec3) {
        let vel = velocity.magnitude_squared() as f32;
        let pos = position.cast();
        self.vertices.push(Vertex {
            position: [pos.x, pos.y, pos.z],
            color: [vel, 0.5 * vel + 0.5, 1.],
        });
    }
}

pub fn start_simulation(
    ui_state: &ui::UIState,
    stop_rx: Receiver<()>,
    last_frame: Option<usize>,
) -> Receiver<Vec<Vertex>> {
    match &ui_state.algorithm {
        ui::Algorithm::Mpm(params) => {
            match &params.constitutive_model {
                ui::ConstituveModels::NeoHookean(nh) => {
                    println!("Simulating MPM w/ NeoHookean Model: {:?}", nh);
                    start_simulation_helper::<crate::MpmSimulation<_>>(
                        // Would be very nice to use type-changing struct update
                        // TODO: Probably better to turn this into a macro
                        // syntax (https://github.com/rust-lang/rfcs/pull/2528) here, but alas, its not stable yet
                        crate::MpmParameters {
                            num_particles: params.num_particles,
                            h: params.h,
                            bounds: params.bounds.clone(),
                            delta_time: params.delta_time,
                            transfer_scheme: params.transfer_scheme,
                            constitutive_model: nh.clone(),
                            cfl: params.cfl,
                        },
                        &ui_state.initial_condition,
                        stop_rx,
                        last_frame,
                    )
                }

                ui::ConstituveModels::NewtonianFluid(nf) => {
                    println!("Simulating MPM w/ Newtonian Fluid Model: {:?}", nf);
                    start_simulation_helper::<crate::MpmSimulation<_>>(
                        crate::MpmParameters {
                            num_particles: params.num_particles,
                            h: params.h,
                            bounds: params.bounds.clone(),
                            delta_time: params.delta_time,
                            transfer_scheme: params.transfer_scheme,
                            constitutive_model: nf.clone(),
                            cfl: params.cfl,
                        },
                        &ui_state.initial_condition,
                        stop_rx,
                        last_frame,
                    )
                }

                ui::ConstituveModels::FixedCorotated(fc) => {
                    println!("Simulating MPM w/ FixedCorotated Model: {:?}", fc);
                    start_simulation_helper::<crate::MpmSimulation<_>>(
                        crate::MpmParameters {
                            num_particles: params.num_particles,
                            h: params.h,
                            bounds: params.bounds.clone(),
                            delta_time: params.delta_time,
                            transfer_scheme: params.transfer_scheme,
                            constitutive_model: fc.clone(),
                            cfl: params.cfl,
                        },
                        &ui_state.initial_condition,
                        stop_rx,
                        last_frame,
                    )
                }

                ui::ConstituveModels::SnowPlasticity(sp) => {
                    println!(
                        "Simulating MPM w/ Drucker-Prager Snow Plastcity Model: {:?}",
                        sp
                    );
                    start_simulation_helper::<crate::MpmSimulation<_>>(
                        crate::MpmParameters {
                            num_particles: params.num_particles,
                            h: params.h,
                            bounds: params.bounds.clone(),
                            delta_time: params.delta_time,
                            transfer_scheme: params.transfer_scheme,
                            constitutive_model: sp.clone(),
                            cfl: params.cfl,
                        },
                        &ui_state.initial_condition,
                        stop_rx,
                        last_frame,
                    )
                }
            }
        }
        ui::Algorithm::Sph(params) => start_simulation_helper::<crate::SphSimulation>(
            params.clone(),
            &ui_state.initial_condition,
            stop_rx,
            last_frame,
        ),
    }
}

fn start_simulation_helper<S: Simulation + 'static>(
    params: S::Parameters,
    initial_condition: &ui::InitialConditions,
    stop_rx: Receiver<()>,
    last_frame: Option<usize>,
) -> Receiver<Vec<Vertex>> {
    let mut s = S::new(params);
    let (tx, rx) = std::sync::mpsc::channel();

    use crate::initial_condition::InitialCondition;
    initial_condition.add_particles(&mut s);

    // TODO: Get this to work
    //println!(
    //"Running simulation with {:?} particles",
    //s.params.num_particles
    //);

    std::thread::spawn(move || {
        let mut num_frames = 0;
        loop {
            let verts = s.simulate_frame();

            tx.send(verts).unwrap();

            match stop_rx.try_recv() {
                Ok(()) => {
                    println!("Stopping simulation.");
                    return;
                }
                Err(std::sync::mpsc::TryRecvError::Disconnected) => {
                    eprintln!("Stop Recieved disconnected.");
                    return;
                }
                Err(std::sync::mpsc::TryRecvError::Empty) => {}
            }

            num_frames += 1;
            if let Some(last_frame) = last_frame {
                if num_frames >= last_frame {
                    println!("Finished simulation w/ {:?} frames.", num_frames);
                    return;
                }
            }
        }
    });

    rx
}
