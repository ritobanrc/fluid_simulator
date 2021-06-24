use super::{Scalar, Vec3};
use nalgebra::Vector3;
use std::ops::Range;

#[derive(Debug, Clone)]
pub struct MpmParameters {
    /// The total number of Lagrangian particles in the simulation
    pub num_particles: usize,
    /// The radius of each particle/the grid spacing
    pub h: Scalar,
    /// The bounds of the simulation.
    pub bounds: Range<Vec3>,
    /// The size of the time step. Larger time steps will simulate faster, but may be unstable or
    /// innaccurate.
    pub delta_time: Scalar,
    /// Whether or not to use APIC (as opposed to PIC). TODO: Handle both FLIP and PIC
    pub use_affine: bool,
    /// Paramters for the Neo-Hookean Constitutive Model. TODO: Add support for other models
    pub constitutive_model: NeoHookeanParameters,
}

impl Default for MpmParameters {
    fn default() -> Self {
        MpmParameters {
            num_particles: 0,
            h: 0.05,
            bounds: Vector3::zeros()..Vector3::new(2., 2., 2.),
            delta_time: 0.01,
            use_affine: false,
            constitutive_model: NeoHookeanParameters::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct NeoHookeanParameters {
    pub youngs_modulus: Scalar,
    pub poissons_ratio: Scalar,
}

impl Default for NeoHookeanParameters {
    fn default() -> Self {
        NeoHookeanParameters {
            youngs_modulus: 10000.,
            poissons_ratio: 0.1,
        }
    }
}

impl NeoHookeanParameters {
    /// Calculates and returns the lame parametrs mu and lambda
    pub fn get_lame_parameters(&self) -> (Scalar, Scalar) {
        let mu = self.youngs_modulus / (2. * (1. + self.poissons_ratio));
        let lambda = self.youngs_modulus * self.poissons_ratio
            / ((1. + self.poissons_ratio) * (1. - 2. * self.poissons_ratio));

        (mu, lambda)
    }
}
