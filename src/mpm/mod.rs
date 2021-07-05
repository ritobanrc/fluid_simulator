mod grid;
pub mod parameters;
mod particles;

use na::Matrix3;
pub use parameters::{MpmParameters, NeoHookean, NewtonianFluid};

use nalgebra::Vector3;

use crate::mpm::parameters::TransferScheme;
use crate::render::Vertex;
use crate::Simulation;

use grid::MpmGrid;
use particles::MpmParticles;

use self::grid::weights::ParticleGridWeights;
use self::parameters::ConstitutiveModel;

type Scalar = f64;
type Vec3 = Vector3<Scalar>;

/// Contains all of the state for the Material Point Method Simulation
pub struct MpmSimulation<CM> {
    pub particles: MpmParticles,
    pub grid: MpmGrid,
    pub params: MpmParameters<CM>,
    pub weights: ParticleGridWeights,
}

impl<CM> MpmSimulation<CM> {
    fn precompute_weights(&mut self) {
        self.weights
            .precompute(&self.grid.data, &self.particles.position);
    }

    fn particles_to_grid(&mut self) {
        // NOTE: Because Dp is proportional to the identity matrix (See Course Notes Pg. 42)
        //       its stored as just a float
        #![allow(non_snake_case)]
        let Dp = (self.params.h * self.params.h) / 3.;
        let Dp_inv = 1. / Dp;

        // Destructure `self` into its constituent parts because the borrow checker
        // isn't smart enough to figure out that we should be able to access these fields
        // simultaneously
        let MpmSimulation {
            ref particles,
            ref mut grid,
            ref weights,
            ref params,
        } = self;

        for p in 0..self.params.num_particles {
            weights.particle_grid_iterator(p).for_each(|(i, weight)| {
                if cfg!(debug_assertions) && !grid.data.coord_in_grid(i) {
                    // FIXME: This check should not be necessary, ideally,
                    // `particle_grid_iterator` simply should not return grid cells that are
                    // out of bounds. But I'm leaving this in here for now cause debugging
                    // boundary issues is a pain
                    eprintln!(
                        "Particle neighborhood cell out of range. Cell: {:?}, Particle Pos: {:?}",
                        i, particles.position[p]
                    );
                    return;
                }

                // Eqn. 172, JSTSS Sigraphh 2016 Course Notes
                *grid.mass_mut(i).unwrap() += particles.mass[p] * weight;

                // Eqn. 128, Course Notes
                let mut v_adjusted = particles.velocity[p];
                if params.transfer_scheme == TransferScheme::APIC {
                    let xi = grid.data.coord_to_pos(i);
                    let xp = particles.position[p];
                    v_adjusted += particles.affine_matrix[p] * Dp_inv * (xi - xp);
                }

                *grid.momentum_mut(i).unwrap() += weight * particles.mass[p] * v_adjusted;
            })
        }
    }

    fn grid_to_particles(&mut self) {
        let MpmSimulation {
            ref mut particles,
            ref grid,
            ref weights,
            ref params,
        } = self;

        for p in 0..params.num_particles {
            //let initial_particle_velocity = particles.velocity[p];
            //particles.velocity[p] = Vector3::zeros();
            if params.transfer_scheme == TransferScheme::APIC {
                particles.affine_matrix[p] = Matrix3::zeros();
            }

            let mut xpic_next_velocity = Vec3::zeros();
            let mut flip_next_velocity = particles.velocity[p];

            weights.particle_grid_iterator(p).for_each(|(i, weight)| {
                if cfg!(debug_assertions) && !grid.data.coord_in_grid(i) {
                    eprintln!(
                        "Particle neighborhood cell out of range. Cell: {:?}, Particle Pos: {:?}",
                        i, particles.position[p]
                    );
                    return;
                }

                xpic_next_velocity += weight * grid.velocity(i).unwrap();
                flip_next_velocity +=
                    weight * (grid.velocity(i).unwrap() - grid.velocity_prev(i).unwrap());

                if params.transfer_scheme == TransferScheme::APIC {
                    let xi = grid.data.coord_to_pos(i);
                    let xp = particles.position[p];
                    particles.affine_matrix[p] +=
                        weight * grid.velocity(i).unwrap() * (xi - xp).transpose();
                }
            });

            particles.velocity[p] = match self.params.transfer_scheme {
                TransferScheme::PIC | TransferScheme::APIC => xpic_next_velocity,
                TransferScheme::FLIP => flip_next_velocity,
                TransferScheme::PIC_FLIP(blend) => {
                    blend * flip_next_velocity + (1. - blend) * xpic_next_velocity
                }
            };
        }
    }

    fn calculate_initial_volumes(&mut self) {
        let h3 = self.params.h * self.params.h * self.params.h;

        self.particles
            .initial_volume
            .reserve(self.params.num_particles);

        for p in 0..self.params.num_particles {
            let density: Scalar = self
                .weights
                .particle_grid_iterator(p)
                .map(|(i, weight)| self.grid.mass(i).unwrap() * weight / h3)
                .sum();

            self.particles
                .initial_volume
                .push(self.particles.mass[p] / density);
        }
    }
}

impl<CM: ConstitutiveModel> MpmSimulation<CM> {
    fn compute_forces(&mut self) {
        for i in 0..self.grid.data.num_cells {
            // Gravity. Fg = -mg
            self.grid.force[i] = self.grid.mass[i] * Vec3::new(0., -1., 0.);
        }

        for p in 0..self.params.num_particles {
            let piola_kirchoff = self.params.constitutive_model.piola_kirchoff(&self, p);

            let MpmSimulation {
                ref particles,
                ref mut grid,
                ref weights,
                params: _,
            } = self;

            weights
                .particle_grid_iterator_grad(p)
                .for_each(|(i, weight_grad)| {
                    let initial_volume = particles.initial_volume[p];

                    let force_contrib = -1.
                        * initial_volume
                        * piola_kirchoff
                        * particles.deformation_gradient[p].transpose()
                        * weight_grad;

                    if cfg!(debug_assertions) && force_contrib.iter().any(|x| x.is_nan()) {
                        eprintln!("Force Contribution is NaN: {:?}", force_contrib);
                        eprintln!(
                            "    Deformation Gradient: {}",
                            particles.deformation_gradient[p]
                        );
                        eprintln!("    Piola Kirchoff: {}", piola_kirchoff);
                    }

                    *grid.force_mut(i).unwrap() += force_contrib;
                });
        }
    }
}

impl<CM> MpmSimulation<CM> {
    fn update_deformation_gradient(&mut self) {
        for p in 0..self.params.num_particles {
            // Eqn. 181 Course Notes
            let fact = Matrix3::identity()
                + self.params.delta_time * {
                    self.weights
                        .particle_grid_iterator_grad(p)
                        .map(|(i, grad)| self.grid.velocity(i).unwrap() * grad.transpose())
                        .sum::<Matrix3<Scalar>>()
                };

            self.particles.deformation_gradient[p] *= fact;
        }
    }

    fn advect_particles(&mut self) {
        for p in 0..self.params.num_particles {
            self.particles.position[p] += self.params.delta_time * self.particles.velocity[p];

            // TODO: Figure out how to correctly handle these boundary conditions
            update_bounds(
                &mut self.particles.position[p],
                self.params.bounds.start,
                self.params.bounds.end,
            );
        }
    }
}

pub(crate) fn update_bounds(position: &mut Vec3, bounds_min: Vec3, bounds_max: Vec3) {
    (0..3).for_each(|i| {
        if position[i] < bounds_min[i] {
            position[i] = bounds_min[i];
        }

        if position[i] > bounds_max[i] {
            position[i] = bounds_max[i];
        }
    })
}

impl<CM: ConstitutiveModel> Simulation for MpmSimulation<CM> {
    type Parameters = MpmParameters<CM>;

    /// Creates a new simulation with the given parameters.
    fn new(params: Self::Parameters) -> Self {
        MpmSimulation {
            particles: MpmParticles::default(),
            grid: MpmGrid::new(&params),
            params,
            weights: ParticleGridWeights::default(),
        }
    }

    /// Adds a particle to the simulation.
    fn add_particle(&mut self, position: Vec3, velocity: Vec3) {
        self.params.num_particles += 1;

        self.particles.add_particle(position, velocity);
        self.weights.add_particle();
    }

    fn simulate_frame(&mut self) -> Vec<Vertex> {
        self.grid.clear_grid();
        self.precompute_weights();

        self.particles_to_grid();

        if self.particles.initial_volume.is_empty() {
            self.calculate_initial_volumes();
        }

        if cfg!(debug_assertions) {
            println!(
                "Particle Mass: {:?}, Grid Mass: {:?}",
                self.particles.total_mass(),
                self.grid.total_mass()
            );

            println!(
                "Particle Momentum: {:?}, Grid Momentum: {:?}",
                self.particles.total_momentum(),
                self.grid.total_momentum()
            );
        }

        self.grid.compute_velocities();
        self.compute_forces();
        self.grid.velocity_update(self.params.delta_time);

        self.update_deformation_gradient();

        self.grid_to_particles();
        self.advect_particles();

        self.create_verts()
    }
}

impl<CM> MpmSimulation<CM> {
    /// Returns an array of `Vertex`es, to be passed to the `render` module.
    /// The color of each vertex is based on the magnitude of the velocity of the particles
    fn create_verts(&self) -> Vec<Vertex> {
        self.particles
            .position
            .iter()
            .zip(self.particles.velocity.iter())
            .map(|(pos, vel)| {
                let pos = pos.cast();
                let vel = vel.magnitude_squared() as f32;
                Vertex {
                    position: [pos.x, pos.y, pos.z],
                    color: [vel, 0.5 * vel + 0.5, 1.],
                }
            })
            .collect()
    }
}
