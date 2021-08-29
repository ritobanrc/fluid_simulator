use super::{MpmSimulation, Scalar};
use nalgebra::Matrix3;

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

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
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

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IsotropicParameters {
    pub youngs_modulus: Scalar,
    pub poissons_ratio: Scalar,
    pub mu: Scalar,
    pub lambda: Scalar,
}

impl IsotropicParameters {
    pub fn new(youngs_modulus: Scalar, poissons_ratio: Scalar) -> Self {
        let mut base = Self {
            youngs_modulus,
            poissons_ratio,
            mu: 0.,
            lambda: 0.,
        };
        base.recalculate_lame_parameters();
        base
    }

    pub fn recalculate_lame_parameters(&mut self) {
        self.mu = self.youngs_modulus / (2. * (1. + self.poissons_ratio));
        self.lambda = self.youngs_modulus * self.poissons_ratio
            / ((1. + self.poissons_ratio) * (1. - 2. * self.poissons_ratio));
    }
}

impl Default for IsotropicParameters {
    fn default() -> Self {
        Self::new(9500., 0.2)
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct NeoHookean(pub IsotropicParameters);

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

        let piola_kirchoff = self.0.mu * (F - F_inv_trans) + self.0.lambda * J.ln() * F_inv_trans;

        if cfg!(debug_assertions) && piola_kirchoff.iter().any(|x| x.is_nan()) {
            eprintln!("Piola Kirchoff is NaN: {:?}", piola_kirchoff);
            eprintln!("    Deformation Gradient: {:?}", F);
            eprintln!("    J: {:?}", J);
            eprintln!("    Lame Parameters: {:?}, {:?}", self.0.mu, self.0.lambda);
        }

        piola_kirchoff
    }
}

#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize)]
pub struct FixedCorotated(pub IsotropicParameters);

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

        let mu_term = 2. * self.0.mu * (F - R);
        let lambda_term = self.0.lambda * (J - 1.) * J * F_inv_trans;

        mu_term + lambda_term
    }
}
