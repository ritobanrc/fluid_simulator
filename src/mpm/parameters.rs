use super::{Scalar, Vec3};
use nalgebra::Vector3;
use std::ops::Range;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MpmParameters<CM> {
    /// The total number of Lagrangian particles in the simulation
    pub num_particles: usize,
    /// The grid spacing
    pub h: Scalar,
    /// The bounds of the simulation.
    pub bounds: Range<Vec3>,
    /// The size of the time step. Larger time steps will simulate faster, but may be unstable or
    /// innaccurate.
    pub delta_time: Scalar,
    /// Which algorithm to use for transferring between particles and grid nodes (PIC, FLIP, APIC,
    /// or a PIC/FLIP blend)
    pub transfer_scheme: TransferScheme,
    /// Paramters for the Constitutive Model.
    pub constitutive_model: CM,
    /// Whether or not to use the CFL timestep restriction,
    pub cfl: Option<f64>,
}

#[derive(Debug, Copy, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[allow(non_camel_case_types)]
pub enum TransferScheme {
    PIC,
    FLIP,
    APIC,
    PIC_FLIP(f64),
}

impl Default for TransferScheme {
    fn default() -> Self {
        TransferScheme::APIC
    }
}

impl<CM: Default> Default for MpmParameters<CM> {
    fn default() -> Self {
        MpmParameters {
            num_particles: 0,
            h: 0.1,
            bounds: Vector3::zeros()..Vector3::new(2., 2., 2.),
            delta_time: 0.01,
            transfer_scheme: TransferScheme::default(),
            constitutive_model: CM::default(),
            cfl: Some(0.8),
        }
    }
}
