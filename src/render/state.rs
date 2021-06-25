use crate::render::{Scene, Vertex};
use wgpu::{
    util::DeviceExt, Adapter, Device, Instance, Queue, RenderPipeline, Surface, SwapChain,
    SwapChainDescriptor,
};
use winit::event::WindowEvent;
use winit::window::Window;

use std::convert::TryInto;

pub enum RenderTarget {
    SwapChain {
        surface: Surface,
        sc_desc: SwapChainDescriptor,
        swap_chain: SwapChain,
    },
    Texture {
        texture: wgpu::Texture,
        texture_view: wgpu::TextureView,
        texture_size: winit::dpi::PhysicalSize<u32>,
        buffer: wgpu::Buffer,
    },
}

/// Contains all of the `wgpu` & `winit` related state
pub struct State {
    pub device: Device,
    pub queue: Queue,
    pub render_target: RenderTarget,
    pub size: winit::dpi::PhysicalSize<u32>,
    pub render_pipeline: Option<RenderPipeline>,
}

fn create_instance() -> Instance {
    Instance::new(wgpu::BackendBit::PRIMARY)
}

async fn create_adapter_with_surface(
    instance: &Instance,
    window: &Window,
) -> (wgpu::Surface, wgpu::Adapter) {
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

async fn crate_adapter(instance: &Instance) -> wgpu::Adapter {
    instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: Default::default(),
            compatible_surface: None,
        })
        .await
        .unwrap()
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

fn create_render_texture(
    texture_size: winit::dpi::PhysicalSize<u32>,
    device: &Device,
) -> RenderTarget {
    let texture_desc = wgpu::TextureDescriptor {
        size: wgpu::Extent3d {
            width: texture_size.width,
            height: texture_size.height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Rgba8UnormSrgb,
        usage: wgpu::TextureUsage::COPY_SRC | wgpu::TextureUsage::RENDER_ATTACHMENT,
        label: None,
    };
    let texture = device.create_texture(&texture_desc);
    let texture_view = texture.create_view(&Default::default());

    let u32_size = std::mem::size_of::<u32>() as u32;

    let output_buffer_size =
        (u32_size * texture_size.width * texture_size.height) as wgpu::BufferAddress;
    let output_buffer_desc = wgpu::BufferDescriptor {
        size: output_buffer_size,
        usage: wgpu::BufferUsage::COPY_DST
                // this tells wpgu that we want to read this buffer from the cpu
                | wgpu::BufferUsage::MAP_READ,
        label: None,
        mapped_at_creation: false,
    };
    let output_buffer = device.create_buffer(&output_buffer_desc);

    RenderTarget::Texture {
        texture,
        texture_view,
        texture_size,
        buffer: output_buffer,
    }
}

fn create_render_pipeline(
    device: &Device,
    scene: &Scene,
    texture_format: wgpu::TextureFormat,
) -> RenderPipeline {
    fn dont_validate(
        mut desc: wgpu::ShaderModuleDescriptor<'_>,
    ) -> wgpu::ShaderModuleDescriptor<'_> {
        desc.flags.remove(wgpu::ShaderFlags::VALIDATION);
        desc
    }

    let vs_module = device.create_shader_module(&wgpu::include_spirv!("../shader.vert.spv"));
    let fs_module =
        device.create_shader_module(&&dont_validate(wgpu::include_spirv!("../shader.frag.spv")));

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
                format: texture_format,
                blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                write_mask: wgpu::ColorWrite::ALL,
            }],
        }),
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::PointList,
            front_face: wgpu::FrontFace::Ccw,
            cull_mode: None,
            strip_index_format: None,
            polygon_mode: wgpu::PolygonMode::Fill,
            clamp_depth: false,
            conservative: false,
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
    pub async fn from_window(window: &Window) -> Self {
        // The instance is a handle to our GPU
        // BackendBit::PRIMARY => Vulkan + Metal + DX12 + Browser WebGPU
        let instance = create_instance();
        let (surface, adapter) = create_adapter_with_surface(&instance, window).await;
        let (device, queue) = create_device_queue(&adapter).await;

        let (sc_desc, swap_chain) = create_swapchain(&surface, &device, window.inner_size());
        let render_target = RenderTarget::SwapChain {
            surface,
            sc_desc,
            swap_chain,
        };

        Self {
            device,
            queue,
            render_target,
            size: window.inner_size(),
            render_pipeline: None,
        }
    }

    pub async fn with_texture(width: u32, height: u32) -> Self {
        // The instance is a handle to our GPU
        let instance = create_instance();
        let adapter = crate_adapter(&instance).await;
        let (device, queue) = create_device_queue(&adapter).await;

        // NOTE: This has nothing to do with winit, its literally just using the
        // PhysicalSize struct
        let texture_size = winit::dpi::PhysicalSize::new(width, height);
        let render_target = create_render_texture(texture_size, &device);

        Self {
            device,
            queue,
            render_target,
            size: texture_size,
            render_pipeline: None,
        }
    }

    pub fn get_buffer(&self) -> Option<&wgpu::Buffer> {
        if let RenderTarget::Texture { buffer, .. } = &self.render_target {
            Some(buffer)
        } else {
            None
        }
    }

    pub fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>) {
        self.size = new_size;
        match &mut self.render_target {
            RenderTarget::SwapChain {
                surface,
                sc_desc,
                swap_chain,
            } => {
                sc_desc.width = new_size.width;
                sc_desc.height = new_size.height;
                *swap_chain = self.device.create_swap_chain(&surface, &sc_desc);
            }
            RenderTarget::Texture { .. } => {
                panic!("Cannot resize a texture target.")
            }
        }
    }

    #[allow(clippy::single_match, clippy::collapsible_match)]
    pub fn input(&mut self, event: &WindowEvent, scene: &mut Scene) -> bool {
        use winit::event::{KeyboardInput, VirtualKeyCode};
        match event {
            WindowEvent::KeyboardInput {
                input:
                    KeyboardInput {
                        virtual_keycode: Some(keycode),
                        ..
                    },
                ..
            } => match keycode {
                VirtualKeyCode::C => {
                    println!("{:#?}", scene.camera);
                }
                _ => {}
            },
            _ => {}
        };
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

        if verts.len() != scene.num_particles as usize {
            scene.num_particles = verts.len() as u32;

            scene.vertex_buffer.destroy();
            scene.vertex_buffer =
                self.device
                    .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                        label: Some("Vertex Buffer"),
                        contents: bytemuck::cast_slice(verts),
                        usage: wgpu::BufferUsage::VERTEX | wgpu::BufferUsage::COPY_DST,
                    })
        } else {
            self.queue
                .write_buffer(&scene.vertex_buffer, 0, bytemuck::cast_slice(verts));
        }
    }

    pub fn render(
        &mut self,
        scene: &Scene,
        egui_state: Option<crate::render::EguiRenderState>,
    ) -> Result<(), wgpu::SwapChainError> {
        let swap_chain_texture;
        let frame = match &self.render_target {
            RenderTarget::SwapChain { swap_chain, .. } => {
                swap_chain_texture = swap_chain.get_current_frame()?;
                &swap_chain_texture.output.view
            }
            RenderTarget::Texture {
                ref texture_view, ..
            } => &texture_view,
        };

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Render Encoder"),
            });

        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("Render Pass Descriptor"),
            color_attachments: &[wgpu::RenderPassColorAttachment {
                view: frame,
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
        let format = match self.render_target {
            RenderTarget::SwapChain { .. } => wgpu::TextureFormat::Bgra8UnormSrgb,
            RenderTarget::Texture { .. } => wgpu::TextureFormat::Rgba8UnormSrgb,
        };
        let render_pipeline = &self
            .render_pipeline
            .get_or_insert_with(|| create_render_pipeline(device, scene, format));

        render_pass.set_pipeline(&render_pipeline);

        render_pass.set_vertex_buffer(0, scene.vertex_buffer.slice(..));
        render_pass.set_bind_group(0, &scene.uniform_state.bind_group, &[]);
        render_pass.draw(0..scene.num_particles, 0..1);

        drop(render_pass);

        if let Some(egui_state) = egui_state {
            egui_state.render_pass.update_texture(
                &self.device,
                &self.queue,
                &egui_state.platform.context().texture(),
            );
            egui_state
                .render_pass
                .update_user_textures(&self.device, &self.queue);
            egui_state.render_pass.update_buffers(
                &self.device,
                &self.queue,
                &egui_state.paint_jobs,
                &egui_state.screen_descriptor,
            );

            egui_state.render_pass.execute(
                &mut encoder,
                &frame,
                &egui_state.paint_jobs,
                &egui_state.screen_descriptor,
                None,
            );
        }

        if let RenderTarget::Texture {
            texture,
            buffer,
            texture_size,
            ..
        } = &self.render_target
        {
            let u32_size = std::mem::size_of::<u32>() as u32;
            encoder.copy_texture_to_buffer(
                wgpu::ImageCopyTexture {
                    texture: &texture,
                    mip_level: 0,
                    origin: wgpu::Origin3d::ZERO,
                },
                wgpu::ImageCopyBuffer {
                    buffer: &buffer,
                    layout: wgpu::ImageDataLayout {
                        offset: 0,
                        bytes_per_row: (u32_size * texture_size.width).try_into().ok(),
                        rows_per_image: texture_size.height.try_into().ok(),
                    },
                },
                wgpu::Extent3d {
                    width: texture_size.width,
                    height: texture_size.height,
                    depth_or_array_layers: 1,
                },
            );
        }

        self.queue.submit(std::iter::once(encoder.finish()));

        Ok(())
    }
}
