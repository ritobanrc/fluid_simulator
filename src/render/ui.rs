use std::fmt::Display;

use egui::{ComboBox, Ui, Widget};

use crate::{mpm::MpmParameters, sph::SphParamaters};

#[derive(Default)]
pub struct UIState {
    pub algorithm: Algorithm,
}

#[derive(Clone, Debug)]
pub enum Algorithm {
    Mpm(MpmParameters),
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
        egui::SidePanel::left("My Window", 300.).show(&ctx, |ui| {
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

            self.algorithm.egui_update(ui);

            if ui.button("Play").clicked() {
                *start_simulation = true;
            }
        });
    }
}

trait EguiInspector {
    fn egui_update(&mut self, ui: &mut Ui);
}

impl EguiInspector for MpmParameters {
    fn egui_update(&mut self, ui: &mut Ui) {
        egui::Slider::new(&mut self.h, 0. ..=0.1).ui(ui);

        ui.label("Bounds");
        ui.horizontal(|ui| {
            ui.label("Min");

            ui.add(egui::DragValue::new(&mut self.bounds.start.x).speed(0.01));
            ui.add(egui::DragValue::new(&mut self.bounds.start.y).speed(0.01));
            ui.add(egui::DragValue::new(&mut self.bounds.start.z).speed(0.01));
        });
        ui.horizontal(|ui| {
            ui.label("Max");

            ui.add(egui::DragValue::new(&mut self.bounds.end.x).speed(0.01));
            ui.add(egui::DragValue::new(&mut self.bounds.end.y).speed(0.01));
            ui.add(egui::DragValue::new(&mut self.bounds.end.z).speed(0.01));
        });
    }
}

impl EguiInspector for SphParamaters {
    fn egui_update(&mut self, ui: &mut Ui) {}
}

/// The Egui-associated state that `State::render` needs
pub struct EguiRenderState<'a> {
    pub(super) platform: &'a egui_winit_platform::Platform,
    pub(super) render_pass: &'a mut egui_wgpu_backend::RenderPass,
    pub(super) paint_jobs: &'a [egui::ClippedMesh],
    pub(super) screen_descriptor: egui_wgpu_backend::ScreenDescriptor,
}
