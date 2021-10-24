use super::{MpmParticles, Scalar};
use nalgebra::Matrix3;

/// Trait that describes a constitutive model. Note that the default implementations for both the
/// Piola Kirchoff stress and the Cauchy stress use each other, so you must override at least one
/// of them (otherwise they will recurse infinitely).
///
/// The current implementation of MPM uses only the Piola Kichoff stress.
pub trait ConstitutiveModel: Send + Sized {
    #[allow(non_snake_case)]
    fn piola_kirchoff(&self, particles: &MpmParticles, p: usize) -> Matrix3<Scalar> {
        let F =
            particles.elastic_deformation_gradient[p] * particles.plastic_deformation_gradient[p];
        let F_inv_trans = F
            .try_inverse()
            .expect("Deformation gradient is not invertible")
            .transpose();
        let J = F.determinant();
        let sigma = self.cauchy_stress(particles, p);

        J * sigma * F_inv_trans
    }

    #[allow(non_snake_case)]
    fn cauchy_stress(&self, particles: &MpmParticles, p: usize) -> Matrix3<Scalar> {
        let F = particles.elastic_deformation_gradient[p];
        let F_trans = F.transpose();
        let J = F.determinant();
        let P = self.piola_kirchoff(particles, p);

        1. / J * P * F_trans
    }

    fn handle_plasticity(&self, _particles: &mut MpmParticles, _p: usize) {}
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
    fn cauchy_stress(&self, particles: &MpmParticles, p: usize) -> Matrix3<Scalar> {
        #![allow(non_snake_case)]

        let deformation_gradient = &particles.elastic_deformation_gradient[p];
        // Hopefully the optimization is smart enough not to calculate J twice
        let J = deformation_gradient.determinant();
        let initial_volume = particles.initial_volume[p];

        let volume = J * initial_volume;
        let density = particles.mass[p] / volume;

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
    fn piola_kirchoff(&self, particles: &MpmParticles, p: usize) -> Matrix3<Scalar> {
        #![allow(non_snake_case)]
        let F = &particles.elastic_deformation_gradient[p];
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
    fn piola_kirchoff(&self, particles: &MpmParticles, p: usize) -> Matrix3<Scalar> {
        #![allow(non_snake_case)]
        let F = &particles.elastic_deformation_gradient[p];
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

/// The Snow Plasticity Constitutive Model described in
/// Stomakhin, A., Schroeder, C., Chai, L., Teran, J., and Selle, A., "A material point method for snow simulation," ACM Transactions on Graphics (SIGGRAPH 2013), 32(4), pp. 102:1-102:10 (2013).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SnowPlasticity {
    /// The hardening coefficient $\xi$. Describes how fast the material breaks after plastic
    /// defomraiton starts.
    pub hardening: Scalar,
    /// The critical compression threshold $\theta_c$, for determining when the material starts breaking
    pub critical_compression: Scalar,
    /// The critical stretch threashold $\theta_s$, also determines when the material starts
    /// breaking. Dry and powdery snow has small $\theta_c$ and $\theta_s$, while wet and chunky
    /// snow has the opposite.
    pub critical_stretch: Scalar,
    /// The initial values of $\mu$ and $\lambda$.
    pub init_isotropic_params: IsotropicParameters,
}

impl Default for SnowPlasticity {
    fn default() -> Self {
        Self {
            hardening: 10.,
            critical_compression: 2.5e-2,
            critical_stretch: 7.5e-3,
            init_isotropic_params: IsotropicParameters::new(1.4e5, 0.2),
        }
    }
}

impl ConstitutiveModel for SnowPlasticity {
    fn piola_kirchoff(&self, particles: &MpmParticles, p: usize) -> Matrix3<Scalar> {
        #![allow(non_snake_case)]
        let Fp = particles.plastic_deformation_gradient[p];
        let Jp = Fp.determinant();

        let Fp_inv_trans = Fp
            .try_inverse()
            .expect("Plastic Deformation gradient is not invertible")
            .transpose();

        let param_scale = Scalar::exp(self.hardening * (1. - Jp));

        let current_isotropic_params = IsotropicParameters {
            youngs_modulus: 0.,
            poissons_ratio: 0.,
            mu: param_scale * self.init_isotropic_params.mu,
            lambda: param_scale * self.init_isotropic_params.lambda,
        };
        let corotated = FixedCorotated(current_isotropic_params);

        corotated.piola_kirchoff(particles, p) * Fp_inv_trans
    }

    fn handle_plasticity(&self, particles: &mut MpmParticles, p: usize) {
        // TODO: Rewrite this whole function
        #![allow(non_snake_case)]
        let Fe = &particles.elastic_deformation_gradient[p];
        let Fp = &particles.plastic_deformation_gradient[p];
        let F = Fe * Fp;

        let mut svd = Fe.svd(true, true);

        svd.singular_values = svd
            .singular_values
            .map(|sigma| sigma.clamp(1. - self.critical_compression, 1. + self.critical_stretch));

        particles.elastic_deformation_gradient[p] = svd
            .recompose()
            .expect("Impossible, U and V_T were constructed.");

        let sigma_inv = svd.singular_values.map(|sigma| 1. / sigma);

        particles.plastic_deformation_gradient[p] = svd.v_t.unwrap().transpose()
            * Matrix3::from_diagonal(&sigma_inv)
            * svd.u.unwrap().transpose()
            * F;
    }
}
