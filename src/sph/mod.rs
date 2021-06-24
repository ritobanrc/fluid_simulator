pub(crate) mod grid;
pub(crate) mod kernels;

use grid::{Coord, Grid};

use crate::render::Vertex;
use crate::{Scalar, Simulation, Vec3};
use kernels::*;

pub struct SphSimulation {
    pub masses: Vec<Scalar>,
    pub positions: Vec<Vec3>,
    pub velocities: Vec<Vec3>,
    pub force: Vec<Vec3>,
    pub grid: Grid,
    pub params: SphParamaters,
}

impl SphSimulation {
    pub(crate) fn position_to_coord(&self, pos: Vec3) -> Coord {
        pos.map(|i| (i / self.params.h) as usize)
    }

    pub(crate) fn coord(&self, index: usize) -> Coord {
        self.position_to_coord(self.positions[index])
    }
}

impl Simulation for SphSimulation {
    type Parameters = SphParamaters;

    fn new(params: SphParamaters) -> Self {
        let grid_bounds = params.bounds + Vec3::new(0.1, 0.1, 0.1);
        SphSimulation {
            masses: Vec::new(),
            positions: Vec::new(),
            velocities: Vec::new(),
            force: Vec::new(),
            grid: Grid::new(
                (grid_bounds.x / params.h) as usize,
                (grid_bounds.y / params.h) as usize,
                (grid_bounds.z / params.h) as usize,
            ),
            params,
        }
    }

    fn simulate_frame(&mut self) -> Vec<Vertex> {
        sph_simulate_frame(self)
    }

    fn add_particle(&mut self, position: Vec3, velocity: Vec3) {
        let index = self.masses.len();
        self.masses.push((10. * self.params.h).powi(3));
        self.positions.push(position);
        self.velocities.push(Vec3::zeros());
        self.force.push(Vec3::zeros());

        self.grid
            .add_particle(self.position_to_coord(position), index);
        self.params.num_particles += 1;
    }
}

fn sph_simulate_frame(s: &mut SphSimulation) -> Vec<Vertex> {
    // This is a separate function because `self` is way too long, I wanna write
    // `s.positions[i]`
    let mut verts = Vec::with_capacity(s.params.num_particles);

    let mut densities: Vec<Scalar> = Vec::with_capacity(s.params.num_particles);

    for i in 0..s.params.num_particles {
        let neighbors = s.grid.get_neighbors(s.coord(i));
        densities.push(
            neighbors
                .map(|j| {
                    s.masses[j] * Poly6Kernel::value(s.positions[i] - s.positions[j], s.params.h)
                })
                .sum(),
        );
    }

    for i in 0..s.params.num_particles {
        let pressure_i: Scalar = s.params.k * (densities[i] - s.params.rest_density);
        let neighbors = s.grid.get_neighbors(s.coord(i));

        let force_pressure = -neighbors
            .clone()
            .map(|j| {
                if i == j {
                    return Vec3::zeros();
                }
                let r_ij = s.positions[i] - s.positions[j];

                if r_ij.magnitude_squared() > s.params.h * s.params.h {
                    return Vec3::zeros();
                }

                let pressure_j = s.params.k * (densities[j] - s.params.rest_density);

                s.masses[j] * (pressure_i + pressure_j) / (2. * densities[j])
                    * SpikyKernel::gradient(r_ij, s.params.h)
            })
            .sum::<Vec3>();

        let force_viscosity = s.params.mu
            * neighbors
                .map(|j| {
                    if i == j {
                        return Vec3::zeros();
                    }
                    let vdiff = s.velocities[j] - s.velocities[i];

                    let r_ij = s.positions[i] - s.positions[j];

                    if r_ij.magnitude_squared() > s.params.h * s.params.h {
                        return Vec3::zeros();
                    }

                    s.masses[j] * vdiff / densities[j]
                        * ViscosityKernel::laplacian(r_ij, s.params.h)
                })
                .sum::<Vec3>();

        let force_gravity = s.params.gravity * densities[i];
        s.force[i] = force_pressure + force_viscosity + force_gravity;
    }

    for i in 0..s.params.num_particles {
        s.velocities[i] = s.velocities[i] + s.params.delta_time / densities[i] * s.force[i];
        let old_coord = s.coord(i);
        s.positions[i] = s.positions[i] + s.params.delta_time * s.velocities[i];

        update_bounds(
            &mut s.positions[i],
            &mut s.velocities[i],
            0.8,
            Vec3::zeros(),
            s.params.bounds,
        );

        let new_coord = s.coord(i);
        if new_coord != old_coord {
            s.grid.update_particle(i, old_coord, new_coord);
        }

        let pos = s.positions[i].cast();
        let vel = s.velocities[i].magnitude_squared() as f32;

        verts.push(Vertex {
            position: [pos.x, pos.y, pos.z],
            color: [vel, 0.5 * vel + 0.5, 1.],
        });
    }

    verts
}

fn update_bounds(
    position: &mut Vec3,
    velocity: &mut Vec3,
    velocity_damping: Scalar,
    bounds_min: Vec3,
    bounds_max: Vec3,
) {
    (0..3).for_each(|i| {
        if position[i] < bounds_min[i] - 0.01 {
            velocity[i] *= -velocity_damping;
            position[i] = bounds_min[i];
        }

        if position[i] > bounds_max[i] + 0.01 {
            velocity[i] *= -velocity_damping;
            position[i] = bounds_max[i];
        }
    })
}

/// A struct containing all of the high-level parameters for the SPH simulation
#[derive(Clone, Debug)]
pub struct SphParamaters {
    pub num_particles: usize,
    /// The time step
    pub delta_time: Scalar,
    /// The radius of the smoothing kernel,
    pub h: Scalar,
    /// The density of the fluid without any forces
    pub rest_density: Scalar,
    /// The ideal gas constant used in the state equation pressure solver
    pub k: Scalar,
    /// The viscosity constant
    pub mu: Scalar,
    /// The force of gravity
    pub gravity: Vec3,
    /// The velocity damping at the boundary
    pub velocity_damping: Scalar,
    /// The bounds of the simulation. Note that this is the maximum bound, the minimum bound
    /// is assumed to be (0, 0, 0)
    pub bounds: Vec3,
}

impl Default for SphParamaters {
    fn default() -> Self {
        Self {
            num_particles: 0,
            delta_time: 0.01,
            h: 0.04,
            rest_density: 1000.,
            k: 4.,
            mu: 8.,
            gravity: Vec3::new(0., -1., 0.),
            velocity_damping: 0.8,
            bounds: Vec3::new(3., 3., 3.),
        }
    }
}
