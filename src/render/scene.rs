/// This module contains all of the _data_ that is going to be sent to the GPU. Essentially, all
/// of the data that describes the scene. This includes vertex buffers, uniforms, the clear color,
/// and camera position.
use crate::render::{Camera, CameraController};
use wgpu::util::DeviceExt;

pub struct Scene {
    pub clear_color: [f32; 3],
    pub camera: Camera,
    pub camera_controller: CameraController,
    pub num_particles: u32,
    pub vertex_buffer: wgpu::Buffer,
    pub uniforms: Uniforms,
    pub uniform_state: UniformState,
}

impl Scene {
    pub fn new(verts: &[Vertex], device: &wgpu::Device, screen_size: (u32, u32)) -> Self {
        let camera = Camera::new(screen_size.0, screen_size.1);
        let mut uniforms = Uniforms::default();
        uniforms.update_view_proj(&camera);

        let uniform_state = UniformState::new(device, uniforms);

        Self {
            clear_color: [0., 0., 0.],
            camera,
            camera_controller: CameraController::new(0.2),
            num_particles: verts.len() as u32,
            vertex_buffer: device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Vertex Buffer"),
                contents: bytemuck::cast_slice(verts),
                usage: wgpu::BufferUsage::VERTEX | wgpu::BufferUsage::COPY_DST,
            }),
            uniforms,
            uniform_state,
        }
    }
}

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Vertex {
    pub position: [f32; 3],
    pub color: [f32; 3],
}

impl Vertex {
    pub fn desc() -> wgpu::VertexBufferLayout<'static> {
        const ATTRS: [wgpu::VertexAttribute; 2] =
            wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3];

        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Vertex>() as wgpu::BufferAddress,
            step_mode: wgpu::InputStepMode::Vertex,
            attributes: &ATTRS,
        }
    }
}

// We need this for Rust to store our data correctly for the shaders
#[repr(C)]
// This is so we can store this in a buffer
#[derive(Debug, Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Uniforms {
    // We can't use nalgebra with bytemuck directly so we'll have
    // to convert the Matrix4 into a 4x4 f32 array
    pub(super) view_proj: [[f32; 4]; 4],
    pub(super) u_point_size: f32,
}

impl Default for Uniforms {
    fn default() -> Self {
        Self {
            view_proj: na::Matrix4::identity().into(),
            u_point_size: 20.,
        }
    }
}

impl Uniforms {
    pub fn update_view_proj(&mut self, camera: &Camera) {
        self.view_proj = camera.build_view_projection_matrix().into();
    }
}
/// This contains all of the GPU-associated state relating to uniforms, including the buffer and
/// the bind group
pub struct UniformState {
    pub buffer: wgpu::Buffer,
    pub bind_group: wgpu::BindGroup,
    pub bind_group_layout: wgpu::BindGroupLayout,
}

impl UniformState {
    fn new(device: &wgpu::Device, uniforms: Uniforms) -> Self {
        let buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Uniform Buffer"),
            contents: bytemuck::cast_slice(&[uniforms]),
            usage: wgpu::BufferUsage::UNIFORM | wgpu::BufferUsage::COPY_DST,
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStage::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    min_binding_size: None,
                    has_dynamic_offset: false,
                },
                count: None,
            }],
            label: Some("uniform_bind_group_layout"),
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: buffer.as_entire_binding(),
            }],
            label: Some("uniform_bind_group"),
        });

        UniformState {
            buffer,
            bind_group,
            bind_group_layout,
        }
    }
}
