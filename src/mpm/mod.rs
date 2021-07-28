mod grid;
pub mod parameters;
mod particles;

use na::Matrix3;
pub use parameters::{FixedCorotated, MpmParameters, NeoHookean, NewtonianFluid};

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
        //
        //  TODO: Support other degree polynomials
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
                    // out of bounds.
                    eprintln!(
                        "Particle neighborhood cell out of range. Cell: {:?}, Particle Pos: {:?}",
                        i, particles.position[p]
                    );
                    return;
                }

                // Eqn. 172, JSTSS Sigraphh 2016 Course Notes
                *grid.mass_mut(i).unwrap() += particles.mass[p] * weight;

                let mut v_affine = particles.velocity[p];
                if params.transfer_scheme == TransferScheme::APIC {
                    let Cp = particles.affine_matrix[p] * Dp_inv;
                    let xi = grid.data.coord_to_pos(i);
                    let xp = particles.position[p];
                    v_affine += Cp * (xi - xp);
                }

                *grid.momentum_mut(i).unwrap() += weight * particles.mass[p] * v_affine;
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
            let mut pic_next_velocity = Vec3::zeros();
            let mut flip_next_velocity = particles.velocity[p];

            if params.transfer_scheme == TransferScheme::APIC {
                particles.affine_matrix[p] = Matrix3::zeros();
            }

            weights.particle_grid_iterator(p).for_each(|(i, weight)| {
                if cfg!(debug_assertions) && !grid.data.coord_in_grid(i) {
                    eprintln!(
                        "Particle neighborhood cell out of range. Cell: {:?}, Particle Pos: {:?}",
                        i, particles.position[p]
                    );
                    return;
                }

                pic_next_velocity += weight * grid.velocity(i).unwrap();
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
                TransferScheme::PIC | TransferScheme::APIC => pic_next_velocity,
                TransferScheme::FLIP => flip_next_velocity,
                TransferScheme::PIC_FLIP(blend) => {
                    blend * flip_next_velocity + (1. - blend) * pic_next_velocity
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
        for &i in &self.grid.valid_grid_indices {
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
    fn update_deformation_gradient(&mut self, delta_time: Scalar) {
        for p in 0..self.params.num_particles {
            // Eqn. 181 Course Notes
            let fact = Matrix3::identity()
                + delta_time * {
                    self.weights
                        .particle_grid_iterator_grad(p)
                        .map(|(i, grad)| self.grid.velocity(i).unwrap() * grad.transpose())
                        .sum::<Matrix3<Scalar>>()
                };

            self.particles.deformation_gradient[p] *= fact;
        }
    }

    fn advect_particles(&mut self, delta_time: Scalar) {
        for p in 0..self.params.num_particles {
            self.particles.position[p] += delta_time * self.particles.velocity[p];

            // TODO: Figure out how to correctly handle these boundary conditions
            update_bounds(
                &mut self.particles.position[p],
                self.params.bounds.start,
                self.params.bounds.end,
            );
        }
    }

    fn calculate_cfl(&self) -> Scalar {
        let mut max_affine_velocity_mag = 0.;
        let d = 3.;
        let gamma_over_k_delta_x = 6. * Scalar::sqrt(d) / self.params.h;

        for p in 0..self.params.num_particles {
            let mut velocity_mag = self.particles.velocity[p].magnitude();
            if self.params.transfer_scheme == TransferScheme::APIC {
                velocity_mag += gamma_over_k_delta_x * self.particles.affine_matrix[p].norm();
            }

            if velocity_mag > max_affine_velocity_mag {
                max_affine_velocity_mag = velocity_mag;
            }
        }

        self.params.h / max_affine_velocity_mag
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
    fn add_particle(&mut self, mass: Scalar, position: Vec3, velocity: Vec3) {
        self.params.num_particles += 1;

        self.particles.add_particle(mass, position, velocity);
        self.weights.add_particle();
    }

    fn simulate_frame(&mut self) -> Vec<Vertex> {
        let mut time_simulated = 0.;
        let mut num_substeps = 0;

        while time_simulated < self.params.delta_time {
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

            let mut cfl_dt = 0.8 * self.calculate_cfl();
            if cfl_dt + time_simulated > self.params.delta_time {
                cfl_dt = self.params.delta_time - time_simulated;
            }
            time_simulated += cfl_dt;
            println!("CFL Time Step: {:?}", cfl_dt);

            self.compute_forces();
            self.grid.velocity_update(cfl_dt);

            self.update_deformation_gradient(cfl_dt);

            self.grid_to_particles();
            self.advect_particles(cfl_dt);

            num_substeps += 1;
        }
        println!("Simulated frame with {:?} substeps. ", num_substeps);

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
