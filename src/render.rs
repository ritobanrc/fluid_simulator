use cgmath::{vec3, Point3};
use futures::executor::block_on;
use winit::{
    event::*,
    event_loop::{ControlFlow, EventLoop},
    window::{Window, WindowBuilder},
};

mod scene;
mod state;

pub use scene::{Scene, Vertex};
pub use state::State;

pub fn open_window(particles: &[Point3<f32>]) {
    let event_loop = EventLoop::new();
    let window = WindowBuilder::new().build(&event_loop).unwrap();

    let mut state = block_on(state::State::new(&window));

    let verts: Vec<Vertex> = particles
        .into_iter()
        .map(|p| Vertex {
            position: [p.x, p.y, p.z],
            color: [0.1, 0.1, 0.1],
        })
        .collect();
    let mut scene = Scene::new(&state.device, &verts);

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
                _ => {}
            },
            WindowEvent::MouseInput {
                device_id: _device_id,
                state,
                ..
            } => {
                if let winit::event::ElementState::Released = state {
                    scene.clear_color = [
                        rand::random::<f32>(),
                        rand::random::<f32>(),
                        rand::random::<f32>(),
                    ]
                }
            }
            WindowEvent::Resized(physical_size) => {
                state.resize(*physical_size);
            }
            WindowEvent::ScaleFactorChanged { new_inner_size, .. } => {
                // new_inner_size is &&mut so we have to dereference it twice
                state.resize(**new_inner_size);
            }
            _ => {}
        },
        Event::RedrawRequested(_) => {
            state.update();
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
