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
        state.render(&scene).expect("Swapchain error");

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

pub fn open_window(rx: Receiver<Vec<Vertex>>) -> Result<(), RecvError> {
    let event_loop = EventLoop::new();
    let window = WindowBuilder::new()
        .with_inner_size(PhysicalSize::new(1280, 720))
        .with_title("Fluid Simulation")
        .build(&event_loop)
        .unwrap();

    let mut state = block_on(State::from_window(&window));

    let size = window.inner_size();
    let verticies = rx.recv()?;
    let mut scene = Scene::new(&verticies, &state.device, (size.width, size.height));

    let mut last_render_time = std::time::Instant::now();
    let mut frame_count = 0;

    event_loop.run(move |event, _, control_flow| match event {
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

            let verticies = rx.recv().expect("Failed to recieve verts");
            state.update(&mut scene, &verticies);
            match state.render(&scene) {
                Ok(_) => {}
                // Recreate the swap_chain if lost
                Err(wgpu::SwapChainError::Lost) => state.resize(state.size),
                // The system is out of memory, we should probably quit
                Err(wgpu::SwapChainError::OutOfMemory) => *control_flow = ControlFlow::Exit,
                // All other errors (Outdated, Timeout) should be resolved by the next frame
                Err(e) => eprintln!("{:?}", e),
            }
        }
        Event::MainEventsCleared => {
            // RedrawRequested will only trigger once, unless we manually
            // request it.
            window.request_redraw();
        }
        _ => {}
    });
}
