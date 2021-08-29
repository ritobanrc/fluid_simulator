use crate::MpmSimulation;
use crate::{Scalar, Vec3};
use itertools::izip;

trait SimulationStatistics {
    fn total_time(&self) -> Scalar;
    fn total_mass(&self) -> Scalar;
    fn total_linear_momentum(&self) -> Vec3;
    fn total_angular_momentum(&self) -> Vec3;
    fn total_energy(&self) -> Scalar;
}

impl<CM> SimulationStatistics for MpmSimulation<CM> {
    fn total_time(&self) -> Scalar {
        self.time
    }

    fn total_mass(&self) -> Scalar {
        self.particles.total_mass()
    }

    fn total_linear_momentum(&self) -> Vec3 {
        self.particles
            .mass
            .iter()
            .zip(&self.particles.velocity)
            .map(|(&m, v)| m * v)
            .sum()
    }

    fn total_angular_momentum(&self) -> Vec3 {
        izip!(
            &self.particles.mass,
            &self.particles.velocity,
            &self.particles.position
        )
        .map(|(&m, v, x)| m * v.cross(x))
        .sum()
    }

    fn total_energy(&self) -> Scalar {
        self.particles
            .mass
            .iter()
            .zip(&self.particles.velocity)
            .map(|(&m, v)| m * v.dot(v))
            .sum()
    }
}
