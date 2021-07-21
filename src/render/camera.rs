use na::{Point3, Vector3};

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
}

impl CameraController {
    pub fn new(speed: f32) -> Self {
        Self { speed, zoom: 0. }
    }

    pub fn process_events(&mut self, event: &winit::event::WindowEvent) -> bool {
        use winit::event::{MouseScrollDelta, WindowEvent};
        match event {
            WindowEvent::MouseWheel { delta, .. } => {
                match delta {
                    MouseScrollDelta::LineDelta(_dx, dy) => self.zoom = *dy,
                    MouseScrollDelta::PixelDelta(winit::dpi::PhysicalPosition {
                        x: _dx,
                        y: dy,
                    }) => self.zoom = *dy as f32,
                }
                true
            }
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
        /*

        let right = forward_norm.cross(&camera.up);

        // Redo radius calc in case the up/ down is pressed.
        let forward = camera.target - camera.eye;
        let forward_mag = forward.magnitude();

        if self.is_right_pressed {
            // Rescale the distance between the target and eye so
            // that it doesn't change. The eye therefore still
            // lies on the circle made by the target and eye.
            camera.eye = camera.target - (forward + right * self.speed).normalize() * forward_mag;
        }
        if self.is_left_pressed {
            camera.eye = camera.target - (forward - right * self.speed).normalize() * forward_mag;
        } */
    }
}
