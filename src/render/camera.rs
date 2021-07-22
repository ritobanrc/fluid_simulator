use na::{Point3, Vector2, Vector3};

#[derive(Debug)]
pub struct Camera {
    eye: Point3<f32>,
    target: Point3<f32>,
    up: Vector3<f32>,
    aspect: f32,
    fovy: f32,
    znear: f32,
    zfar: f32,
}

impl Camera {
    pub fn build_view_projection_matrix(&self) -> na::Matrix4<f32> {
        #[rustfmt::skip]
        pub const OPENGL_TO_WGPU_MATRIX: na::Matrix4<f32> = na::Matrix4::new(
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.5, 0.0,
            0.0, 0.0, 0.5, 1.0,
        );

        let view = na::Matrix4::look_at_rh(&self.eye, &self.target, &self.up);
        let proj = na::Matrix4::new_perspective(self.aspect, self.fovy, self.znear, self.zfar);

        OPENGL_TO_WGPU_MATRIX * proj * view
    }

    /// Creates a new camera with a bunch of default settings.
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            // position the camera one unit up and 2 units back
            // +z is out of the screen
            eye: Point3::new(-2., 1., 1.),
            // have it look at the origin
            target: Point3::new(0.5, 0.5, 1.),
            // which way is "up"
            up: na::Vector3::new(0., 1., 0.),
            aspect: width as f32 / height as f32,
            fovy: 45.0,
            znear: 0.1,
            zfar: 100.0,
        }
    }
}

pub struct CameraController {
    speed: f32,
    zoom: f32,
    middle_mouse: bool,
    shift: bool,
    delta: Vector2<f32>,
}

impl CameraController {
    pub fn new(speed: f32) -> Self {
        Self {
            speed,
            zoom: 0.,
            middle_mouse: false,
            shift: false,
            delta: Vector2::zeros(),
        }
    }

    pub fn process_events(&mut self, event: &winit::event::DeviceEvent) -> bool {
        use winit::event::{DeviceEvent, KeyboardInput, MouseScrollDelta};
        match event {
            DeviceEvent::MouseWheel { delta } => {
                match delta {
                    MouseScrollDelta::LineDelta(_dx, dy) => self.zoom = *dy * 20., // Naively assume each line is 20 pixels
                    MouseScrollDelta::PixelDelta(winit::dpi::PhysicalPosition {
                        x: _dx,
                        y: dy,
                    }) => self.zoom = *dy as f32,
                }
                true
            }
            DeviceEvent::MouseMotion { delta } => {
                self.delta = Vector2::new(delta.0, delta.1).cast();
                true
            }
            DeviceEvent::Button { button: 2, state } => {
                match state {
                    winit::event::ElementState::Pressed => self.middle_mouse = true,
                    winit::event::ElementState::Released => self.middle_mouse = false,
                }
                true
            }
            DeviceEvent::Key(KeyboardInput {
                virtual_keycode: keycode,
                state,
                ..
            }) => match keycode {
                Some(winit::event::VirtualKeyCode::LShift)
                | Some(winit::event::VirtualKeyCode::RShift) => {
                    match state {
                        winit::event::ElementState::Pressed => self.shift = true,
                        winit::event::ElementState::Released => self.shift = false,
                    };
                    true
                }
                _ => false,
            },
            _ => false,
        }
    }

    pub fn update_camera(&mut self, camera: &mut Camera) {
        let forward = camera.target - camera.eye;
        let forward_norm = forward.normalize();
        let forward_mag = forward.magnitude();

        // Zoom
        if self.zoom != 0. && (self.zoom < 0. || forward_mag > self.speed) {
            camera.eye += self.zoom.signum() * forward_norm * self.speed;
            self.zoom = 0.; // We've handled it
        }

        if self.middle_mouse {
            if self.shift {
                let mut shift = Vector3::zeros();
                shift += self.delta.y * camera.up;
                shift -= self.delta.x * camera.up.cross(&-forward);

                camera.eye += 0.004 * shift;
                camera.target += 0.004 * shift;

                self.delta = Vector2::zeros();
            } else {
                // Orbit
                if self.delta.y != 0. {
                    let axis = na::Unit::new_normalize(camera.up.cross(&-forward));
                    camera.eye = camera.target
                        - na::Rotation3::from_axis_angle(&axis, -0.015 * self.delta.y) * forward;

                    self.delta.y = 0.;
                }

                if self.delta.x != 0. {
                    let axis = na::Unit::new_normalize(camera.up);
                    let forward = camera.target - camera.eye; // recalculate in case this was modified by the delta.y path
                    camera.eye = camera.target
                        - na::Rotation3::from_axis_angle(&axis, -0.015 * self.delta.x) * forward;

                    self.delta.x = 0.;
                }
            }
        }
    }
}
