use super::{MpmSimulation, Scalar, Vec3};
use nalgebra::{Matrix3, Vector3};
use std::ops::Range;

#[derive(Debug, Clone)]
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
    /// Paramters for the Neo-Hookean Constitutive Model. TODO: Add support for other models
    pub constitutive_model: CM,
}

#[derive(Debug, Copy, Clone, PartialEq)]
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
            h: 0.05,
            bounds: Vector3::zeros()..Vector3::new(2., 2., 2.),
            delta_time: 0.01,
            transfer_scheme: TransferScheme::default(),
            constitutive_model: CM::default(),
        }
    }
}

/// Trait that describes a constitutive model. Note that the default implementations for both the
/// Piola Kirchoff stress and the Cauchy stress use each other, so you must override at least one
/// of them (otherwise they will recurse infinitely).
///
/// The current implementation of MPM uses only the Piola Kichoff stress.
pub trait ConstitutiveModel: Send + Sized {
    #[allow(non_snake_case)]
    fn piola_kirchoff(&self, s: &MpmSimulation<Self>, p: usize) -> Matrix3<Scalar> {
        let F = s.particles.deformation_gradient[p];
        let F_inv_trans = F
            .try_inverse()
            .expect("Deformation gradient is not invertible")
            .transpose();
        let J = F.determinant();
        let sigma = self.cauchy_stress(s, p);

        J * sigma * F_inv_trans
    }

    #[allow(non_snake_case)]
    fn cauchy_stress(&self, s: &MpmSimulation<Self>, p: usize) -> Matrix3<Scalar> {
        let F = s.particles.deformation_gradient[p];
        let F_trans = F.transpose();
        let J = F.determinant();
        let P = self.piola_kirchoff(s, p);

        1. / J * P * F_trans
    }
}

#[derive(Debug, Clone)]
pub struct NeoHookean {
    pub youngs_modulus: Scalar,
    pub poissons_ratio: Scalar,
    pub mu: Scalar,
    pub lambda: Scalar,
}

impl NeoHookean {
    pub fn new(youngs_modulus: Scalar, poissons_ratio: Scalar) -> Self {
        let mu = youngs_modulus / (2. * (1. + poissons_ratio));
        let lambda =
            youngs_modulus * poissons_ratio / ((1. + poissons_ratio) * (1. - 2. * poissons_ratio));

        NeoHookean {
            youngs_modulus,
            poissons_ratio,
            mu,
            lambda,
        }
    }

    pub fn recalculate_lame_parameters(&mut self) {
        self.mu = self.youngs_modulus / (2. * (1. + self.poissons_ratio));
        self.lambda = self.youngs_modulus * self.poissons_ratio
            / ((1. + self.poissons_ratio) * (1. - 2. * self.poissons_ratio));
    }
}

impl Default for NeoHookean {
    fn default() -> Self {
        Self::new(9500., 0.05)
    }
}

impl ConstitutiveModel for NeoHookean {
    fn piola_kirchoff(&self, s: &MpmSimulation<Self>, p: usize) -> Matrix3<Scalar> {
        #![allow(non_snake_case)]
        let F = s.particles.deformation_gradient[p];
        let J = F.determinant();

        if J <= 0. {
            eprintln!("J = {:?} will result in NaNs", J);
        }

        let F_inv_trans = F
            .try_inverse()
            .expect("Deformation gradient is not invertible")
            .transpose();

        let piola_kirchoff = self.mu * (F - F_inv_trans) + self.lambda * J.ln() * F_inv_trans;

        if cfg!(debug_assertions) && piola_kirchoff.iter().any(|x| x.is_nan()) {
            eprintln!("Piola Kirchoff is NaN: {:?}", piola_kirchoff);
            eprintln!("    Deformation Gradient: {:?}", F);
            eprintln!("    J: {:?}", J);
            eprintln!("    Lame Parameters: {:?}, {:?}", self.mu, self.lambda);
        }

        piola_kirchoff
    }
}

#[derive(Debug, Clone)]
pub struct NewtonianFluid {
    pub k: Scalar,
    pub rest_density: Scalar,
}

impl NewtonianFluid {
    pub fn new(k: Scalar, rest_density: Scalar) -> Self {
        NewtonianFluid { k, rest_density }
    }
}

impl Default for NewtonianFluid {
    fn default() -> Self {
        Self::new(4., 500.)
    }
}

impl ConstitutiveModel for NewtonianFluid {
    fn cauchy_stress(&self, s: &MpmSimulation<Self>, p: usize) -> Matrix3<Scalar> {
        #![allow(non_snake_case)]

        let deformation_gradient = s.particles.deformation_gradient[p];
        // Hopefully the optimization is smart enough not to calculate J twice
        let J = deformation_gradient.determinant();
        let initial_volume = s.particles.initial_volume[p];

        let volume = J * initial_volume;
        let density = s.particles.mass[p] / volume;

        let pressure = self.k * (density - self.rest_density);

        Matrix3::from_diagonal_element(-pressure)
    }
}

#[derive(Debug, Clone)]
pub struct FixedCorotated {
    // TODO: Remove this code duplication w/ Neo-Hookean
    pub youngs_modulus: Scalar,
    pub poissons_ratio: Scalar,
    pub mu: Scalar,
    pub lambda: Scalar,
}

impl FixedCorotated {
    pub fn new(youngs_modulus: Scalar, poissons_ratio: Scalar) -> Self {
        let mu = youngs_modulus / (2. * (1. + poissons_ratio));
        let lambda =
            youngs_modulus * poissons_ratio / ((1. + poissons_ratio) * (1. - 2. * poissons_ratio));

        FixedCorotated {
            youngs_modulus,
            poissons_ratio,
            mu,
            lambda,
        }
    }

    pub fn recalculate_lame_parameters(&mut self) {
        self.mu = self.youngs_modulus / (2. * (1. + self.poissons_ratio));
        self.lambda = self.youngs_modulus * self.poissons_ratio
            / ((1. + self.poissons_ratio) * (1. - 2. * self.poissons_ratio));
    }
}

impl Default for FixedCorotated {
    fn default() -> Self {
        Self::new(9500., 0.05)
    }
}

impl ConstitutiveModel for FixedCorotated {
    fn piola_kirchoff(&self, s: &MpmSimulation<Self>, p: usize) -> Matrix3<Scalar> {
        #![allow(non_snake_case)]
        let F = s.particles.deformation_gradient[p];
        let J = F.determinant();

        let svd = F.svd(true, true);
        let R = svd.u.unwrap() * svd.v_t.unwrap();

        let F_inv_trans = F
            .try_inverse()
            .expect("Deformation gradient is not invertible")
            .transpose();

        let mu_term = 2. * self.mu * (F - R);
        let lambda_term = self.lambda * (J - 1.) * J * F_inv_trans;

        mu_term + lambda_term
    }
}
