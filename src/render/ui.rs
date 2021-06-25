use std::fmt::Display;

use egui::{ComboBox, Ui, Widget};

use crate::{
    mpm::{
        parameters::{NeoHookean, NewtonianFluid},
        MpmParameters,
    },
    sph::SphParamaters,
};

#[derive(Default)]
pub struct UIState {
    pub(super) algorithm: Algorithm,
}

#[derive(Clone, Debug)]
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
        match (self, other) {
            (Mpm(_), Mpm(_)) | (Sph(_), Sph(_)) => true,
            _ => false,
        }
    }
}

impl EguiInspector for Algorithm {
    fn egui_update(&mut self, ui: &mut Ui) {
        match self {
            Algorithm::Mpm(params) => params.egui_update(ui),
            Algorithm::Sph(params) => params.egui_update(ui),
        }
    }
}

impl Display for Algorithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            &Algorithm::Mpm(_) => write!(f, "Mpm"),
            &Algorithm::Sph(_) => write!(f, "Sph"),
        }
    }
}

impl UIState {
    pub fn egui_update(
        &mut self,
        ctx: egui::CtxRef,
        _frame: &mut epi::Frame,
        start_simulation: &mut bool,
    ) {
        egui::SidePanel::left("My Window", 200.).show(&ctx, |ui| {
            egui::Grid::new("side-panel-grid")
                .striped(true)
                .spacing([40., 4.])
                .show(ui, |ui| {
                    ComboBox::from_label("")
                        .selected_text(format!("{}", self.algorithm))
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.algorithm,
                                Algorithm::Mpm(MpmParameters::default()),
                                "Mpm",
                            );

                            ui.selectable_value(
                                &mut self.algorithm,
                                Algorithm::Sph(SphParamaters::default()),
                                "Sph",
                            );
                        });

                    ui.end_row();

                    self.algorithm.egui_update(ui);
                });

            if ui.button("Play").clicked() {
                *start_simulation = true;
            }
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
        ui.vertical(|ui| {
            ui.horizontal(|ui| {
                ui.add_space(10.);
                ui.label("Min");
                self.start.egui_update(ui);
            });

            ui.horizontal(|ui| {
                ui.add_space(10.);
                ui.label("Max");
                self.end.egui_update(ui);
            });
        });
        ui.end_row();
    }
}

impl<CM: EguiInspector> EguiInspector for MpmParameters<CM> {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.label("h: ");
        ui.add(egui::Slider::new(&mut self.h, 0. ..=0.1));
        ui.end_row();

        ui.label("dt: ");
        ui.add(egui::Slider::new(&mut self.delta_time, 0. ..=0.05));
        ui.end_row();

        ui.label("Bounds");
        ui.end_row();
        self.bounds.egui_update(ui);

        self.constitutive_model.egui_update(ui);
    }
}

#[derive(Debug, Clone)]
pub(super) enum ConstituveModels {
    NeoHookean(NeoHookean),
    NewtonianFluid(NewtonianFluid),
}

impl PartialEq for ConstituveModels {
    fn eq(&self, other: &Self) -> bool {
        use ConstituveModels::*;
        match (self, other) {
            (NeoHookean(_), NeoHookean(_)) | (NewtonianFluid(_), NewtonianFluid(_)) => true,
            _ => false,
        }
    }
}

impl Default for ConstituveModels {
    fn default() -> Self {
        ConstituveModels::NeoHookean(NeoHookean::default())
    }
}

impl Display for ConstituveModels {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            &ConstituveModels::NeoHookean(_) => write!(f, "NeoHookean"),
            &ConstituveModels::NewtonianFluid(_) => write!(f, "Newtonian Fluid"),
        }
    }
}

impl EguiInspector for ConstituveModels {
    fn egui_update(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
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
        });
        ui.end_row();

        egui::CollapsingHeader::new(format!("{} Constitutive Model", self))
            .default_open(true)
            .show(ui, |ui| match self {
                Self::NeoHookean(a) => a.egui_update(ui),
                Self::NewtonianFluid(a) => a.egui_update(ui),
            });
    }
}

impl EguiInspector for NeoHookean {
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
