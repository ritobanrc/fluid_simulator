use eyre::WrapErr;
use std::fmt::Display;

use egui::{DragValue, Ui, Widget};

use crate::{
    initial_condition::{Block, InitialCondition, Sphere},
    mpm::{
        models::{FixedCorotated, IsotropicParameters, NeoHookean, NewtonianFluid, SnowPlasticity},
        parameters::MpmParameters,
    },
    sph::SphParamaters,
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct UIState {
    pub(super) algorithm: Algorithm,
    pub(super) initial_condition: InitialConditions,
    pub(super) json_filename: String,
    pub(super) particle_size: f32,
    #[serde(skip)]
    pub(super) render_dt: std::time::Duration,
}

impl Default for UIState {
    fn default() -> Self {
        UIState {
            algorithm: Default::default(),
            initial_condition: Default::default(),
            json_filename: String::from("./example/simulation_settings.json"),
            particle_size: 20.,
            render_dt: Default::default(),
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(super) enum Algorithm {
    Mpm(MpmParameters<ConstituveModels>),
    Sph(SphParamaters),
}

impl Default for Algorithm {
    fn default() -> Self {
        Algorithm::Mpm(MpmParameters::default())
    }
}

impl PartialEq for Algorithm {
    fn eq(&self, other: &Self) -> bool {
        use Algorithm::*;
        matches!((self, other), (Mpm(_), Mpm(_)) | (Sph(_), Sph(_)))
    }
}

impl EguiInspector for Algorithm {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.heading("Algorithm");

        ui.horizontal(|ui| {
            ui.selectable_value(self, Algorithm::Sph(SphParamaters::default()), "Sph");
            ui.selectable_value(self, Algorithm::Mpm(MpmParameters::default()), "Mpm");
        });

        ui.end_row();

        match self {
            Algorithm::Mpm(params) => params.egui_update(ui),
            Algorithm::Sph(params) => params.egui_update(ui),
        }
    }
}

impl Display for Algorithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Algorithm::Mpm(_) => write!(f, "Mpm"),
            Algorithm::Sph(_) => write!(f, "Sph"),
        }
    }
}

#[derive(serde::Serialize, serde::Deserialize, Debug, Clone)]
pub(super) enum InitialConditions {
    Block(crate::initial_condition::Block),
    Sphere(crate::initial_condition::Sphere),
}

impl InitialCondition for InitialConditions {
    fn add_particles<S: crate::Simulation>(&self, s: &mut S) {
        match self {
            InitialConditions::Block(a) => a.add_particles(s),
            InitialConditions::Sphere(a) => a.add_particles(s),
        }
    }
}

impl Default for InitialConditions {
    fn default() -> Self {
        InitialConditions::Block(Block::default())
    }
}

impl PartialEq for InitialConditions {
    fn eq(&self, other: &Self) -> bool {
        use InitialConditions::*;
        matches!((self, other), (Block(_), Block(_)) | (Sphere(_), Sphere(_)))
    }
}

impl EguiInspector for InitialConditions {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.heading("Initial Condition");

        ui.horizontal(|ui| {
            ui.selectable_value(self, InitialConditions::Block(Block::default()), "Block");
            ui.selectable_value(self, InitialConditions::Sphere(Sphere::default()), "Sphere");
        });
        ui.end_row();

        match self {
            Self::Block(a) => a.egui_update(ui),
            Self::Sphere(a) => a.egui_update(ui),
        }
    }
}

impl EguiInspector for Block {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.label("Size: ");
        ui.end_row();
        self.size.egui_update(ui);

        ui.label("Spacing: ");
        DragValue::new(&mut self.spacing)
            .speed(0.01)
            .clamp_range(0. ..=1.)
            .ui(ui);
        ui.end_row();

        ui.label("Jitter: ");
        self.jitter.egui_update(ui);
        ui.end_row();

        ui.label("Density: ");
        DragValue::new(&mut self.density)
            .speed(0.5)
            .clamp_range(0. ..=1000.)
            .ui(ui);
        ui.end_row();
    }
}

impl EguiInspector for Sphere {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.label("Number of Particles: ");
        DragValue::new(&mut self.num_particles).speed(10).ui(ui);
        ui.end_row();

        ui.label("Center: ");
        self.center.egui_update(ui);
        ui.end_row();

        ui.label("Radius: ");
        DragValue::new(&mut self.radius).speed(0.05).ui(ui);
        ui.end_row();

        ui.label("Density: ");
        DragValue::new(&mut self.density)
            .speed(0.5)
            .clamp_range(0. ..=1000.)
            .ui(ui);
        ui.end_row();
    }
}

impl UIState {
    pub fn egui_update(
        &mut self,
        ctx: egui::CtxRef,
        _frame: &mut epi::Frame,
        start_simulation: &mut bool,
        stop_tx: &Option<std::sync::mpsc::Sender<()>>,
    ) {
        egui::SidePanel::left("My Window", 200.).show(&ctx, |ui| {
            egui::Grid::new("side-panel-grid-1")
                .striped(true)
                .spacing([40., 8.])
                .show(ui, |ui| {
                    ui.end_row();

                    self.initial_condition.egui_update(ui);

                    ui.separator();
                    ui.end_row();

                    self.algorithm.egui_update(ui);

                    let mut fps = 1_000_000 / self.render_dt.as_micros() as u64;
                    ui.add(egui::Label::new("FPS").heading());
                    ui.add(egui::DragValue::new(&mut fps));
                    ui.end_row();

                    ui.add(egui::Label::new("Particle Size").heading());
                    ui.add(
                        egui::DragValue::new(&mut self.particle_size)
                            .clamp_range(0. ..=100.)
                            .speed(0.5),
                    );
                    ui.end_row();
                });

            ui.separator();
            ui.end_row();

            ui.vertical_centered_justified(|ui| {
                if ui.button("Start Simulation").clicked() {
                    *start_simulation = true;
                }

                if let Some(stop_tx) = stop_tx {
                    if ui.button("Stop Simulation").clicked() {
                        stop_tx
                            .send(())
                            .expect("Channel to Simulation failed to send.");
                    }
                }
            });

            ui.separator();
            ui.end_row();

            ui.heading("Import/Export JSON Settings: ");
            egui::Grid::new("side-panel-grid-2")
                .striped(true)
                .spacing([40., 4.])
                .show(ui, |ui| {
                    ui.label("Filename");
                    ui.text_edit_singleline(&mut self.json_filename);
                });

            ui.vertical_centered_justified(|ui| {
                if ui.button("Export JSON Settings").clicked() {
                    let json = serde_json::to_string_pretty(self).expect("Serialization failed");
                    if let Err(e) = std::fs::write(&self.json_filename, json) {
                        eprintln!(
                            "Failed to write json settings to {}: {}",
                            self.json_filename, e
                        );
                    }
                }
                if ui.button("Import JSON Settings").clicked() {
                    let ui_state: Self = std::fs::read(&self.json_filename)
                        .wrap_err_with(|| {
                            format!("Failed to read JSON settings file: {}", &self.json_filename)
                        })
                        .and_then(|json| {
                            serde_json::from_slice(&json)
                                .wrap_err("Serde failed to deserialize JSON.")
                        })
                        .unwrap();

                    *self = ui_state;
                }
            });
        });
    }
}

trait EguiInspector {
    fn egui_update(&mut self, ui: &mut Ui);
}

impl<T: egui::math::Numeric + na::Scalar, D: na::Dim, S: na::storage::StorageMut<T, D>>
    EguiInspector for na::Vector<T, D, S>
{
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            for x in self.iter_mut() {
                // TODO: support different speed values
                ui.add(egui::DragValue::new(x).speed(0.01));
            }
        });
    }
}

impl<T: EguiInspector> EguiInspector for std::ops::Range<T> {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.label("\tMin");
        self.start.egui_update(ui);

        ui.end_row();

        ui.label("\tMax");
        self.end.egui_update(ui);
        ui.end_row();
    }
}

impl<CM: EguiInspector> EguiInspector for MpmParameters<CM> {
    fn egui_update(&mut self, ui: &mut Ui) {
        use crate::mpm::parameters::TransferScheme;
        ui.label("h: ");
        ui.add(egui::Slider::new(&mut self.h, 0. ..=0.1));
        ui.end_row();

        ui.label("dt: ");
        ui.add(egui::Slider::new(&mut self.delta_time, 0. ..=0.05));
        ui.end_row();

        ui.label("Bounds");
        ui.end_row();
        self.bounds.egui_update(ui);

        ui.label("CFL");
        let mut use_cfl = self.cfl.is_some();
        ui.horizontal(|ui| {
            if ui.checkbox(&mut use_cfl, "").changed() {
                match use_cfl {
                    true => self.cfl = Some(0.8),
                    false => self.cfl = None,
                }
            }
            if let Some(ref mut cfl) = self.cfl {
                ui.add(egui::DragValue::new(cfl).clamp_range(0. ..=1.5).speed(0.05));
            }
        });
        ui.end_row();

        ui.label("Boundary Mu");
        ui.add(egui::Slider::new(&mut self.boundary_mu, 0. ..=1.));
        ui.end_row();

        ui.separator();
        ui.end_row();

        ui.heading("Transfer Scheme");
        ui.end_row();
        ui.radio_value(&mut self.transfer_scheme, TransferScheme::PIC, "PIC");
        ui.end_row();
        ui.radio_value(&mut self.transfer_scheme, TransferScheme::FLIP, "FLIP");
        ui.end_row();
        ui.radio_value(&mut self.transfer_scheme, TransferScheme::APIC, "APIC");
        ui.end_row();

        // NOTE: This is a kinda ugly hack to keep the radio button selected even when the blend
        // value is changed. Basically, when the radio button is selected, we're making egui think
        // that clicking the radio button sets it to the blend value it already is.
        let pic_flip_blend_select = if let TransferScheme::PIC_FLIP(blend) = self.transfer_scheme {
            TransferScheme::PIC_FLIP(blend)
        } else {
            TransferScheme::PIC_FLIP(0.95)
        };

        ui.radio_value(
            &mut self.transfer_scheme,
            pic_flip_blend_select,
            "PIC/FLIP Blend",
        );

        if let TransferScheme::PIC_FLIP(ref mut blend) = self.transfer_scheme {
            ui.add(
                egui::DragValue::new(blend)
                    .clamp_range(0. ..=1.)
                    .speed(0.02),
            );
        }
        ui.end_row();

        self.constitutive_model.egui_update(ui);
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(super) enum ConstituveModels {
    NeoHookean(NeoHookean),
    NewtonianFluid(NewtonianFluid),
    FixedCorotated(FixedCorotated),
    SnowPlasticity(SnowPlasticity),
}

impl PartialEq for ConstituveModels {
    fn eq(&self, other: &Self) -> bool {
        use ConstituveModels::*;
        matches!(
            (self, other),
            (NeoHookean(_), NeoHookean(_))
                | (NewtonianFluid(_), NewtonianFluid(_))
                | (FixedCorotated(_), FixedCorotated(_))
                | (SnowPlasticity(_), SnowPlasticity(_))
        )
    }
}

impl Default for ConstituveModels {
    fn default() -> Self {
        ConstituveModels::NeoHookean(NeoHookean::default())
    }
}

impl Display for ConstituveModels {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            ConstituveModels::NeoHookean(_) => write!(f, "NeoHookean"),
            ConstituveModels::NewtonianFluid(_) => write!(f, "Newtonian Fluid"),
            ConstituveModels::FixedCorotated(_) => write!(f, "Fixed Corotated"),
            ConstituveModels::SnowPlasticity(_) => write!(f, "Snow Plasticity"),
        }
    }
}

impl EguiInspector for ConstituveModels {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.separator();
        ui.end_row();

        ui.heading("Constitutive Model");

        egui::ComboBox::from_label("")
            .selected_text(format!("{}", self))
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    self,
                    ConstituveModels::NeoHookean(NeoHookean::default()),
                    "NeoHookean",
                );
                ui.selectable_value(
                    self,
                    ConstituveModels::NewtonianFluid(NewtonianFluid::default()),
                    "Newtonian Fluid",
                );
                ui.selectable_value(
                    self,
                    ConstituveModels::FixedCorotated(FixedCorotated::default()),
                    "Fixed Corotated",
                );
                ui.selectable_value(
                    self,
                    ConstituveModels::SnowPlasticity(SnowPlasticity::default()),
                    "Snow Plasticity",
                );
            });
        ui.end_row();

        match self {
            Self::NeoHookean(a) => a.egui_update(ui),
            Self::NewtonianFluid(a) => a.egui_update(ui),
            Self::FixedCorotated(a) => a.egui_update(ui),
            Self::SnowPlasticity(a) => a.egui_update(ui),
        }
    }
}

impl EguiInspector for IsotropicParameters {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.label("Young's Modulus: ");
        let response_youngs_modulus =
            ui.add(egui::Slider::new(&mut self.youngs_modulus, 0. ..=10_000.));
        ui.end_row();

        ui.label("Poisson's Ratio: ");
        let response_poissons_ratio =
            ui.add(egui::Slider::new(&mut self.poissons_ratio, 0. ..=0.5));
        ui.end_row();

        if response_youngs_modulus.changed() || response_poissons_ratio.changed() {
            self.recalculate_lame_parameters();
        }
    }
}

impl EguiInspector for NeoHookean {
    fn egui_update(&mut self, ui: &mut Ui) {
        self.0.egui_update(ui);
    }
}

impl EguiInspector for FixedCorotated {
    fn egui_update(&mut self, ui: &mut Ui) {
        self.0.egui_update(ui);
    }
}

impl EguiInspector for SnowPlasticity {
    fn egui_update(&mut self, ui: &mut Ui) {
        self.init_isotropic_params.egui_update(ui);

        ui.label("Hardening: ");
        ui.add(egui::Slider::new(&mut self.hardening, 0. ..=20.));
        ui.end_row();

        ui.label("Critical Compression: ");
        ui.add(egui::Slider::new(
            &mut self.critical_compression,
            1.5e-2..=4.0e-2,
        ));
        ui.end_row();

        ui.label("Critical Stretch: ");
        ui.add(egui::Slider::new(
            &mut self.critical_stretch,
            4.0e-3..=1.0e-2,
        ));
        ui.end_row();
    }
}

impl EguiInspector for NewtonianFluid {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.label("k: ");
        ui.add(egui::Slider::new(&mut self.k, 0. ..=10.));
        ui.end_row();

        ui.label("Rest Density: ");
        ui.add(egui::Slider::new(&mut self.rest_density, 0. ..=2000.));
        ui.end_row();
    }
}

impl EguiInspector for SphParamaters {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.label("h: ");
        ui.add(egui::Slider::new(&mut self.h, 0. ..=0.1));
        ui.end_row();

        ui.label("dt: ");
        ui.add(egui::Slider::new(&mut self.delta_time, 0. ..=0.05));
        ui.end_row();

        ui.label("Rest density: ");
        ui.add(egui::Slider::new(&mut self.rest_density, 0. ..=5000.));
        ui.end_row();

        ui.label("k: ");
        ui.add(egui::Slider::new(&mut self.k, 0. ..=10.));
        ui.end_row();

        ui.label("mu (viscosity): ");
        ui.add(egui::Slider::new(&mut self.mu, 0. ..=10.));
        ui.end_row();

        ui.label("gravity: ");
        self.gravity.egui_update(ui);
        ui.end_row();

        ui.label("velocity damping");
        ui.add(egui::Slider::new(&mut self.velocity_damping, 0. ..=1.));
        ui.end_row();

        ui.label("bounds");
        self.bounds.egui_update(ui);
        ui.end_row();
    }
}

/// The Egui-associated state that `State::render` needs
pub struct EguiRenderState<'a> {
    pub(super) platform: &'a egui_winit_platform::Platform,
    pub(super) render_pass: &'a mut egui_wgpu_backend::RenderPass,
    pub(super) paint_jobs: &'a [egui::ClippedMesh],
    pub(super) screen_descriptor: egui_wgpu_backend::ScreenDescriptor,
}
