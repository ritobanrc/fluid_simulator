use super::{Scalar, Vec3};
use na::Matrix3;

/// Contains all of the particle data: position, velocity, force, etc.
pub struct MpmParticles {
    pub mass: Vec<Scalar>,
    pub position: Vec<Vec3>,
    pub velocity: Vec<Vec3>,
    pub initial_volume: Vec<Scalar>,
    pub elastic_deformation_gradient: Vec<Matrix3<Scalar>>,
    pub plastic_deformation_gradient: Vec<Matrix3<Scalar>>,
    // TODO: Remove the memory overhead of storing the affine matrices when APIC isn't used
    pub affine_matrix: Vec<Matrix3<Scalar>>,
    // Probably other stuff: `mass`, `cell`, `initial_volume`, `density`, `pressure`, etc.
}

impl Default for MpmParticles {
    fn default() -> Self {
        Self {
            mass: Vec::new(),
            position: Vec::new(),
            velocity: Vec::new(),
            initial_volume: Vec::new(),
            elastic_deformation_gradient: Vec::new(),
            plastic_deformation_gradient: Vec::new(),
            affine_matrix: Vec::new(),
        }
    }
}

impl MpmParticles {
    /// Adds a new particle at the given position.
    /// All other data (velocity, force, etc.) set to zero by default
    pub(crate) fn add_particle(&mut self, mass: Scalar, position: Vec3, velocity: Vec3) {
        // TODO: Units -- if mass is in kg, this should probably be a lot smaller, like
        // 0.001
        self.mass.push(mass);
        self.position.push(position);
        self.velocity.push(velocity);
        self.elastic_deformation_gradient.push(Matrix3::identity());

        // TODO: Don't use this memory if not necessary
        self.plastic_deformation_gradient.push(Matrix3::identity());
        self.affine_matrix.push(Matrix3::zeros());

        // NOTE: I am intentionally not populating `initial_volume`, it will be populated
        // on the first time step in `calculate_initial_volumes`
    }

    pub fn total_mass(&self) -> Scalar {
        self.mass.iter().sum()
    }

    pub fn total_momentum(&self) -> Vec3 {
        self.mass
            .iter()
            .zip(&self.velocity)
            .map(|(&m, v)| m * v)
            .sum()
    }
}
