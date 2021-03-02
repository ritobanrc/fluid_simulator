use wgpu::util::DeviceExt;

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Vertex {
    pub position: [f32; 3],
    pub color: [f32; 3],
}

impl Vertex {
    pub fn desc() -> wgpu::VertexBufferLayout<'static> {
        const ATTRS: [wgpu::VertexAttribute; 2] =
            wgpu::vertex_attr_array![0 => Float3, 1 => Float3];

        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Vertex>() as wgpu::BufferAddress,
            step_mode: wgpu::InputStepMode::Vertex,
            attributes: &ATTRS,
        }
    }
}

pub struct Scene {
    pub clear_color: [f32; 3],
    pub vertex_buffer: wgpu::Buffer,
}

impl Scene {
    pub fn new(device: &wgpu::Device, verts: &[Vertex]) -> Self {
        Self {
            clear_color: [0., 0., 0.],
            vertex_buffer: device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Vertex Buffer"),
                contents: bytemuck::cast_slice(verts),
                usage: wgpu::BufferUsage::VERTEX,
            }),
        }
    }
}
