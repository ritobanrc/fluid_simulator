mod grid;
mod parameters;
mod particles;

pub use parameters::MpmParameters;

use nalgebra::Vector3;
use std::ops::RangeInclusive as Range;

use crate::render::Vertex;
use crate::Simulation;

use grid::MpmGrid;
use particles::MpmParticles;

type Scalar = f32;
type Vec3 = Vector3<Scalar>;

/// Contains all of the state for the Material Point Method Simulation
pub struct MpmSimulation {
    pub particles: MpmParticles,
    pub grid: MpmGrid,
    pub params: MpmParameters,
}

impl MpmSimulation {
    /// Creates a new simulation with the given parameters.
    pub fn new(params: MpmParameters) -> MpmSimulation {
        MpmSimulation {
            particles: MpmParticles::default(),
            grid: MpmGrid,
            params,
        }
    }

    /// Adds a particle to the simulation.
    pub(crate) fn add_particle(&mut self, position: crate::Vec3) {
        self.params.num_particles += 1;

        // FIXME: Migrate the entire crate over to `nalgebra` so you don't have to deal
        // with crap like this
        let position = Vec3::new(position.x, position.y, position.z);
        self.particles.add_particle(position);
    }
}

impl Simulation for MpmSimulation {
    fn simulate_frame(&mut self) -> Vec<Vertex> {
        self.create_verts()
    }
}

impl MpmSimulation {
    /// Returns an array of `Vertex`es, to be passed to the `render` module.
    /// The color of each vertex is based on the magnitude of the velocity of the particles
    fn create_verts(&self) -> Vec<Vertex> {
        self.particles
            .position
            .iter()
            .zip(self.particles.velocity.iter())
            .map(|(pos, vel)| {
                let vel = vel.magnitude_squared();
                Vertex {
                    position: [pos.x, pos.y, pos.z],
                    color: [vel, 0.5 * vel + 0.5, 1.],
                }
            })
            .collect()
    }
}
