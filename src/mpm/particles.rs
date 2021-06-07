use super::{Scalar, Vec3};
use na::Matrix3;

/// Contains all of the particle data: position, velocity, force, etc.
pub struct MpmParticles {
    pub position: Vec<Vec3>,
    pub velocity: Vec<Vec3>,
    pub force: Vec<Vec3>,
    pub deformation_gradient: Vec<Matrix3<Scalar>>,
    // Probably other stuff: `mass`, `cell`, `initial_volume`, `density`, `pressure`, etc.
}

impl Default for MpmParticles {
    fn default() -> Self {
        Self {
            position: Vec::new(),
            velocity: Vec::new(),
            force: Vec::new(),
            deformation_gradient: Vec::new(),
        }
    }
}

impl MpmParticles {
    /// Adds a new particle at the given position.
    /// All other data (velocity, force, etc.) set to zero by default
    pub(crate) fn add_particle(&mut self, position: Vec3) {
        self.position.push(position);
        self.velocity.push(Vec3::zeros());
        self.force.push(Vec3::zeros());
        self.deformation_gradient.push(Matrix3::zeros());
    }
}
