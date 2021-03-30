use crate::render::{Scene, Vertex};
use wgpu::{
    Adapter, Device, Instance, Queue, RenderPipeline, Surface, SwapChain, SwapChainDescriptor,
};
use winit::event::WindowEvent;
use winit::window::Window;

/// Contains all of the `wgpu` related state
pub struct State {
    pub surface: Surface,
    pub device: Device,
    pub queue: Queue,
    pub sc_desc: SwapChainDescriptor,
    pub swap_chain: SwapChain,
    pub size: winit::dpi::PhysicalSize<u32>,
    pub render_pipeline: Option<RenderPipeline>,
}

fn create_instance() -> Instance {
    Instance::new(wgpu::BackendBit::PRIMARY)
}

async fn create_surface(instance: &Instance, window: &Window) -> (wgpu::Surface, wgpu::Adapter) {
    // SAFETY: I'm not actually sure LMAO, I'm just doing what the docs tell me
    let surface = unsafe { instance.create_surface(window) };
    let adapter = instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: Default::default(),
            compatible_surface: Some(&surface),
        })
        .await
        .unwrap();

    (surface, adapter)
}

async fn create_device_queue(adapter: &Adapter) -> (Device, Queue) {
    adapter
        .request_device(
            &wgpu::DeviceDescriptor {
                label: Some("Default device"),
                features: wgpu::Features::empty(),
                limits: wgpu::Limits::default(),
            },
            None, // Trace path
        )
        .await
        .unwrap()
}

fn create_swapchain(
    surface: &Surface,
    device: &Device,
    size: winit::dpi::PhysicalSize<u32>,
) -> (wgpu::SwapChainDescriptor, wgpu::SwapChain) {
    let sc_desc = SwapChainDescriptor {
        usage: wgpu::TextureUsage::RENDER_ATTACHMENT,
        format: wgpu::TextureFormat::Bgra8UnormSrgb,
        width: size.width,
        height: size.height,
        present_mode: wgpu::PresentMode::Fifo,
    };
    let swap_chain = device.create_swap_chain(&surface, &sc_desc);

    (sc_desc, swap_chain)
}

fn create_render_pipeline(device: &Device, scene: &Scene) -> RenderPipeline {
    let vs_module = device.create_shader_module(&wgpu::include_spirv!("../shader.vert.spv"));
    let fs_module = device.create_shader_module(&wgpu::include_spirv!("../shader.frag.spv"));

    let render_pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("Render Pipeline Layout"),
        bind_group_layouts: &[&scene.uniform_state.bind_group_layout],
        push_constant_ranges: &[],
    });

    let render_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("Render Pipeline"),
        layout: Some(&render_pipeline_layout),
        vertex: wgpu::VertexState {
            module: &vs_module,
            entry_point: "main",
            buffers: &[Vertex::desc()],
        },
        fragment: Some(wgpu::FragmentState {
            // 2.
            module: &fs_module,
            entry_point: "main",
            targets: &[wgpu::ColorTargetState {
                format: wgpu::TextureFormat::Bgra8UnormSrgb, // NOTE: Don't hardcode this
                color_blend: wgpu::BlendState::REPLACE,
                alpha_blend: wgpu::BlendState::REPLACE,
                write_mask: wgpu::ColorWrite::ALL,
            }],
        }),
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::PointList,
            front_face: wgpu::FrontFace::Ccw,
            cull_mode: wgpu::CullMode::Back,
            strip_index_format: None,
            polygon_mode: wgpu::PolygonMode::Fill,
        },
        depth_stencil: None,
        multisample: wgpu::MultisampleState {
            count: 1,
            mask: !0,
            alpha_to_coverage_enabled: false,
        },
    });

    render_pipeline
}

impl State {
    // Creating some of the wgpu types requires async code
    pub async fn new(window: &Window) -> Self {
        // The instance is a handle to our GPU
        // BackendBit::PRIMARY => Vulkan + Metal + DX12 + Browser WebGPU
        let instance = create_instance();
        let (surface, adapter) = create_surface(&instance, window).await;
        let (device, queue) = create_device_queue(&adapter).await;
        let (sc_desc, swap_chain) = create_swapchain(&surface, &device, window.inner_size());

        Self {
            surface,
            device,
            queue,
            sc_desc,
            swap_chain,
            size: window.inner_size(),
            render_pipeline: None,
        }
    }

    pub fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>) {
        self.size = new_size;
        self.sc_desc.width = new_size.width;
        self.sc_desc.height = new_size.height;
        self.swap_chain = self.device.create_swap_chain(&self.surface, &self.sc_desc);
    }

    pub fn input(&mut self, event: &WindowEvent, scene: &mut Scene) -> bool {
        scene.camera_controller.process_events(event)
    }

    pub fn update(&mut self, scene: &mut Scene, verts: &[Vertex]) {
        scene.camera_controller.update_camera(&mut scene.camera);
        scene.uniforms.update_view_proj(&scene.camera);
        self.queue.write_buffer(
            &scene.uniform_state.buffer,
            0,
            bytemuck::cast_slice(&[scene.uniforms]),
        );

        self.queue
            .write_buffer(&scene.vertex_buffer, 0, bytemuck::cast_slice(verts));
    }

    pub fn render(&mut self, scene: &Scene) -> Result<(), wgpu::SwapChainError> {
        let frame = self.swap_chain.get_current_frame()?.output;
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Render Encoder"),
            });

        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("Render Pass Descriptor"),
            color_attachments: &[wgpu::RenderPassColorAttachmentDescriptor {
                attachment: &frame.view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color {
                        r: scene.clear_color[0] as f64,
                        g: scene.clear_color[1] as f64,
                        b: scene.clear_color[2] as f64,
                        a: 1.0,
                    }),
                    store: true,
                },
            }],
            depth_stencil_attachment: None,
        });

        let device = &self.device;
        let render_pipeline = &self
            .render_pipeline
            .get_or_insert_with(|| create_render_pipeline(device, scene));

        render_pass.set_pipeline(&render_pipeline);

        render_pass.set_vertex_buffer(0, scene.vertex_buffer.slice(..));
        render_pass.set_bind_group(0, &scene.uniform_state.bind_group, &[]);
        render_pass.draw(0..scene.num_particles, 0..1);

        drop(render_pass);

        self.queue.submit(std::iter::once(encoder.finish()));

        Ok(())
    }
}
