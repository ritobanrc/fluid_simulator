mod grid;
mod parameters;
mod particles;

use na::Matrix3;
pub use parameters::MpmParameters;

use nalgebra::Vector3;

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
            grid: MpmGrid::new(&params),
            params,
        }
    }

    /// Adds a particle to the simulation.
    pub(crate) fn add_particle(&mut self, position: [Scalar; 3], velocity: [Scalar; 3]) {
        self.params.num_particles += 1;

        // FIXME: Migrate the entire crate over to `nalgebra` so you don't have to deal
        // with crap like this
        let position = Vec3::new(position[0], position[1], position[2]);
        let velocity = Vec3::new(velocity[0], velocity[1], velocity[2]);
        self.particles.add_particle(position, velocity);
    }

    fn particles_to_grid(&mut self) {
        // NOTE: Because Dp is proportional to the identity matrix (See Course Notes Pg. 42)
        //       its stored as just a float
        #![allow(non_snake_case)]
        let Dp = (self.params.h * self.params.h) / 3.;
        let Dp_inv = 1. / Dp;

        for p in 0..self.params.num_particles {
            self.grid
                .data
                .clone() // FIXME: This clone is unnecessary, and just a hack to get around the borrow checker
                .particle_grid_iterator(self.particles.position[p])
                .for_each(|(i, weight)| {
                    if cfg!(debug_assertions) && !self.grid.data.coord_in_grid(i) {
                        // FIXME: This check should not be necessary, ideally,
                        // `particle_grid_iterator` simply should not return grid cells that are
                        // out of bounds. But I'm leaving this in here for now cause debugging
                        // boundary issues is a pain
                        eprintln!(
                            "Particle neighborhood cell out of range. Cell: {:?}, Particle Pos: {:?}",
                            i, self.particles.position[p]
                        );
                        return;
                    }

                    // Eqn. 172, JSTSS Sigraphh 2016 Course Notes
                    *self.grid.mass_mut(i).unwrap() += self.particles.mass[p] * weight;
                    // Eqn. 128, Course Notes

                    let mut v_adjusted = self.particles.velocity[p];
                    if self.params.use_affine {
                        let xi = self.grid.data.coord_to_pos(i);
                        let xp = self.particles.position[p];
                        v_adjusted += self.particles.affine_matrix[p] * Dp_inv * (xi - xp);
                    }
                    *self.grid.momentum_mut(i).unwrap() += weight * self.particles.mass[p] * v_adjusted;
                })
        }
    }

    fn grid_to_particles(&mut self) {
        for p in 0..self.params.num_particles {
            self.particles.velocity[p] = Vector3::zeros();
            if self.params.use_affine {
                self.particles.affine_matrix[p] = Matrix3::zeros();
            }

            self
                .grid
                .data
                .clone()
                .particle_grid_iterator(self.particles.position[p])
                .for_each(|(i, weight)| {
                    if cfg!(debug_assertions) && !self.grid.data.coord_in_grid(i) {
                        eprintln!(
                            "Particle neighborhood cell out of range. Cell: {:?}, Particle Pos: {:?}",
                            i, self.particles.position[p]
                        );
                        return;
                    }

                    self.particles.velocity[p] += weight * self.grid.velocity(i).unwrap();

                    if self.params.use_affine {
                        let xi = self.grid.data.coord_to_pos(i);
                        let xp = self.particles.position[p];
                        self.particles.affine_matrix[p] += weight * self.grid.velocity(i).unwrap() * (xi - xp).transpose();
                    }
                });
        }
    }

    fn compute_forces(&mut self) {
        for i in 0..self.grid.data.num_cells {
            // Just gravity, for now. Fg = -mg
            self.grid.force[i] = self.grid.mass[i] * Vec3::new(0., -1., 0.);
        }

        for p in 0..self.params.num_particles {
            let piola_kirchoff =
                self.neo_hookean_piola_kirchoff(self.particles.deformation_gradient[p]);

            self.grid
                .data
                .clone()
                .particle_grid_iterator_grad(self.particles.position[p])
                .for_each(|(i, weight_grad)| {
                    let initial_volume = 1.0; // TODO: Compute this from a density parameter as described in SSCTS 13

                    let force_contrib = -1.
                        * initial_volume
                        * piola_kirchoff
                        * self.particles.deformation_gradient[p].transpose()
                        * weight_grad;

                    if cfg!(debug_assertions) && force_contrib.iter().any(|x| x.is_nan()) {
                        eprintln!("Force Contribution is NaN: {:?}", force_contrib);
                        eprintln!(
                            "    Deformation Gradient: {}",
                            self.particles.deformation_gradient[p]
                        );
                        eprintln!("    Piola Kirchoff: {}", piola_kirchoff);
                    }

                    *self.grid.force_mut(i).unwrap() += force_contrib;
                });
        }
    }

    fn fixed_corotated_piola_kirchoff(
        &self,
        deformation_gradient: Matrix3<Scalar>,
    ) -> Matrix3<Scalar> {
        #![allow(non_snake_case)]

        let F = deformation_gradient;
        let svd = F.svd(true, true);
        let R = svd.u.unwrap() * svd.v_t.unwrap();
        let J = F.determinant();

        let (mu, lambda) = self.params.constitutive_model.get_lame_parameters();
        2. * mu * (F - R) + lambda * (J - 1.) * J * F.try_inverse().unwrap().transpose()
    }

    fn neo_hookean_piola_kirchoff(&self, deformation_gradient: Matrix3<Scalar>) -> Matrix3<Scalar> {
        #![allow(non_snake_case)]
        let F = deformation_gradient;
        let J = F.determinant();

        let (mu, lambda) = self.params.constitutive_model.get_lame_parameters();

        let F_inv_trans = F
            .try_inverse()
            .expect("Deformation gradient is not invertible")
            .transpose();

        // TODO: Figure out what base this logarithm is supposed to be
        let piola_kirchoff = mu * (F - F_inv_trans) + lambda * J.log10() * F_inv_trans;

        if cfg!(debug_assertions) && piola_kirchoff.iter().any(|x| x.is_nan()) {
            eprintln!("Piola Kirchoff is NaN: {:?}", piola_kirchoff);
            eprintln!("    Deformation Gradient: {:?}", F);
            eprintln!("    J: {:?}", J);
            eprintln!("    Lame Parameters: {:?}, {:?}", mu, lambda);
        }

        piola_kirchoff
    }

    fn update_deformation_gradient(&mut self) {
        for p in 0..self.params.num_particles {
            // Eqn. 181 Course Notes
            let fact = Matrix3::identity()
                + self.params.delta_time * {
                    self.grid
                        .data
                        .particle_grid_iterator_grad(self.particles.position[p])
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

impl Simulation for MpmSimulation {
    fn simulate_frame(&mut self) -> Vec<Vertex> {
        self.grid.clear_grid();
        self.particles_to_grid();

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
