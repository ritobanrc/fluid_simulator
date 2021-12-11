use crate::MpmSimulation;
use crate::SphSimulation;
use crate::{Scalar, Vec3};
use itertools::izip;

pub trait SimulationStatistics {
    fn total_time(&self) -> Scalar;
    fn total_mass(&self) -> Scalar;
    fn total_linear_momentum(&self) -> Vec3;
    fn total_angular_momentum(&self) -> Vec3;
    fn total_energy(&self) -> Scalar;
    fn total_volume(&self) -> Scalar;
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

    fn total_volume(&self) -> Scalar {
        #![allow(non_snake_case)]
        (0..self.params.num_particles)
            .map(|p| {
                let Fe = &self.particles.elastic_deformation_gradient[p];
                let Fp = &self.particles.plastic_deformation_gradient[p];
                let F = Fe * Fp;

                let J = F.determinant();
                let init_volume = self.particles.initial_volume[p];
                J * init_volume
            })
            .sum()
    }
}

impl SimulationStatistics for SphSimulation {
    fn total_time(&self) -> Scalar {
        self.time
    }

    fn total_mass(&self) -> Scalar {
        self.masses.iter().sum()
    }

    fn total_linear_momentum(&self) -> Vec3 {
        self.masses
            .iter()
            .zip(&self.velocities)
            .map(|(&m, v)| m * v)
            .sum()
    }

    fn total_angular_momentum(&self) -> Vec3 {
        izip!(&self.masses, &self.velocities, &self.positions)
            .map(|(&m, v, x)| m * v.cross(x))
            .sum()
    }

    fn total_energy(&self) -> Scalar {
        self.masses
            .iter()
            .zip(&self.velocities)
            .map(|(&m, v)| m * v.dot(v))
            .sum()
    }

    fn total_volume(&self) -> Scalar {
        self.density
            .iter()
            .zip(&self.masses)
            .map(|(rho, m)| m / rho)
            .sum()
    }
}
