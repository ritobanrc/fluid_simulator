use futures::executor::block_on;
use std::sync::mpsc::{Receiver, RecvError};
use winit::{
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

pub fn open_window(rx: Receiver<Vec<Vertex>>) -> Result<(), RecvError> {
    let event_loop = EventLoop::new();
    let window = WindowBuilder::new().build(&event_loop).unwrap();

    let mut state = block_on(state::State::new(&window));

    let size = window.inner_size();
    let verticies = rx.recv()?;
    let mut scene = Scene::new(&verticies, &state.device, (size.width, size.height));

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
