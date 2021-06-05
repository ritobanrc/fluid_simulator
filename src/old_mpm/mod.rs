//mod algo;
mod grid;
mod math;

use crate::render::Vertex;
use crate::{Scalar, Simulation, Vec3};
use cgmath::{vec3, InnerSpace, Matrix, Matrix3, SquareMatrix, Zero};
use itertools::iproduct;

use std::ops::RangeInclusive as Range; // NOTE: All ranges should be inclusive

use grid::{Coord, MpmGrid};
use math::*;

/// This struct stores all of the quantities on each particle
pub struct MpmSimulation {
    frame_counter: usize,
    mass: Vec<Scalar>,
    volume: Vec<Scalar>,
    position: Vec<Vec3>,
    velocity: Vec<Vec3>,
    //affine_matrices: Vec<Matrix3<Scalar>>,
    deformation_gradient: Vec<Matrix3<Scalar>>,
    //grid: MpmGrid,
    pub params: MpmParmaters,
}

impl MpmSimulation {
    pub fn new(params: MpmParmaters) -> Self {
        MpmSimulation {
            frame_counter: 0,
            mass: Vec::new(),
            volume: Vec::new(),
            position: Vec::new(),
            velocity: Vec::new(),
            deformation_gradient: Vec::new(),
            //grid: MpmGrid::new(&params),
            //affine_matrices: Vec::new(),
            params,
        }
    }

    pub fn add_particle(&mut self, position: Vec3) {
        self.mass.push(0.1);
        self.position.push(position);
        self.velocity.push(Vec3::zero());
        //self.affine_matrices.push(Matrix3::identity());
        self.deformation_gradient.push(Matrix3::identity());
        self.params.num_particles += 1;
    }
}

// TODO: Get rid of this
macro_rules! unwrap_or_continue {
    ($res:expr) => {
        match $res {
            Some(val) => val,
            None => {
                continue;
            }
        }
    };
}

impl Simulation for MpmSimulation {
    fn simulate_frame(&mut self) -> Vec<Vertex> {
        let mut verts = Vec::with_capacity(self.params.num_particles);

        //let grid_bounds: Coord = (self.params.bounds / self.params.h)
        //.cast()
        //.expect("Failed to cast f32 to usize");

        let mut grid = MpmGrid::new(&self.params);

        // Step 1: APIC Particle to Grid Transfer
        for p in 0..self.params.num_particles {
            let node: Coord = unwrap_or_continue!((self.position[p] / self.params.h).cast()); // FIXME: Get rid of unwrap or continue
                                                                                              //.expect(&format!(
                                                                                              //"Failed to cast f32 to usize: {:?}",
                                                                                              //self.position[p] / self.params.h
                                                                                              //));

            //let mut dp = Matrix3::zero();

            // Instead of doing every possible pair of particle-node, or iterating over all
            // particles for each node, we're going to instead iterate over all _particles_, and
            // add their contribution to nodes in the 3x3x3 neighborhood around them.
            for (ax, ay, az) in iproduct!(0..3, 0..3, 0..3) {
                if any(node + vec3(ax, ay, az), |n| n == 0) {
                    continue;
                }
                let other_node = vec3(node.x + ax - 1, node.y + ay - 1, node.z + az - 1);
                // FIXME: This feels like there's a bug here
                let other_node_pos: Vec3 =
                    other_node.cast().expect("Failed to cast usize to f32") * self.params.h;

                let weight = weight(other_node_pos, self.position[p], self.params.h);
                match grid.mass_mut(other_node) {
                    Some(mass) => {
                        *mass += weight * self.mass[p];
                    }
                    None => continue,
                }

                let momentum = weight
                    * self.mass[p]
                    * (self.velocity[p]/*+ self.affine_matrices[p] * dp_inv * diff*/);

                match grid.momentum_mut(other_node) {
                    Some(p) => *p += momentum,
                    None => continue,
                }
                //let diff = other_node_pos - self.position[p];
                //dp += weight * outer(diff, diff);
            }

            //let dp_inv = dp
            //.invert()
            //.expect(&format!("Dp matrix is non-invertible: {:?}", dp));

            //for (ax, ay, az) in iproduct!(0..3, 0..3, 0..3) {
            //// NOTE: There has to be some way to fix this code duplication.
            //let other_node = vec3(node.x + ax - 1, node.y + ay - 1, node.z + az - 1);
            //let other_node_pos: Vec3 =
            //self.params.h * other_node.cast().expect("Failed to cast usize to f32");

            //let weight = weight(other_node_pos, self.position[p], self.params.h);
            ////let diff = other_node_pos - self.position[p];

            //let momentum = weight
            //* self.mass[p]
            //* (self.velocities[p][>+ self.affine_matrices[p] * dp_inv * diff<]);

            //match grid.momentum_mut(other_node) {
            //Some(p) => *p += momentum,
            //None => continue,
            //}
            //}
        }

        let mut density = vec![0.; self.params.num_particles];

        let mut finding_initial_volume = false;
        if self.volume.is_empty() {
            self.volume.reserve(self.params.num_particles);
            finding_initial_volume = true;
        }

        let h3 = self.params.h * self.params.h * self.params.h;
        let cell_densities: Vec<_> = grid.mass.iter().map(|m| m / h3).collect();

        for p in 0..self.params.num_particles {
            let node: Coord = unwrap_or_continue!((self.position[p] / self.params.h).cast());

            //let mut density = 0.;

            for (ax, ay, az) in iproduct!(0..3, 0..3, 0..3) {
                if any(node + vec3(ax, ay, az), |n| n == 0) {
                    continue;
                }
                let other_node = vec3(node.x + ax - 1, node.y + ay - 1, node.z + az - 1);
                let other_node_pos: Vec3 = self.params.h * unwrap_or_continue!(other_node.cast());
                let idx = grid.coord_to_index(other_node);

                if idx >= cell_densities.len() {
                    continue;
                }

                let weight = weight(other_node_pos, self.position[p], self.params.h);

                density[p] += cell_densities[idx] * weight;
            }

            if finding_initial_volume {
                self.volume.push(self.mass[p] / density[p]);
            }
        }

        // Step 2: Compute Grid Velocities from Momentea
        let mut velocities: Vec<_> = grid
            .momentum
            .iter()
            .zip(grid.mass.iter())
            .map(|(momentum, &mass)| {
                if mass != 0. {
                    momentum / mass
                } else {
                    Vec3::zero()
                }
            })
            .collect();

        // TODO: Step 3, Identify Grid Degrees of Freedom
        // Step 4: Compute explicit grid forces
        //         TODO: Actually calculate the forces using a deformation gradient
        let g = -100. * Vec3::unit_y();
        let mut force: Vec<_> = (0..grid.mass.len())
            .map(|i| {
                if grid.mass[i] != 0. {
                    g * grid.mass[i]
                } else {
                    Vec3::zero()
                }
            })
            .collect();

        // Sum over all the particles
        //      Get the particle _initial_ Volume (constant)
        //      derivative of PE wrt to deformation graident -- i.e. the piola kirchoff stress
        //          P = mu(F - F^-T) + lambda * log(J) F^-T
        //      multiply by transpose of deformation gradient
        //      multiply by gradient of kernel function
        //
        for p in 0..self.params.num_particles {
            let F = self.deformation_gradient[p];
            //let F_trans = F.transpose();
            //let F_inv_trans = F_trans
            //.invert()
            //.expect("Failed to invert deformation gradient");

            //let youngs_modulus = 0.1;
            //let poissons_ratio = 0.475;

            //let mu = youngs_modulus / (2. * (1. + poissons_ratio));
            //let lambda = youngs_modulus * poissons_ratio
            // / ((1. + poissons_ratio) * (1. - 2. * poissons_ratio));
            let J = F.determinant();
            //let piola_kirchoff = mu * (F - F_inv_trans) + lambda * J.ln() * F_inv_trans;
            let initial_volume = self.volume[p];
            //let force_contrib = initial_volume * piola_kirchoff * F_trans;
            let volume = initial_volume * J;
            //let force_contrib = initial_volume * piola_kirchoff * F_trans;
            let pressure = 4. * (density[p] - 500.);
            let cauchy_stress = -pressure * Matrix3::identity();
            let force_contrib = volume * cauchy_stress;

            //let pressure = density[p] - 10.;
            let node: Coord = unwrap_or_continue!((self.position[p] / self.params.h).cast()); // FIXME: Don't use unwrap or continue

            // TODO: Abstract over this
            for (ax, ay, az) in iproduct!(0..3, 0..3, 0..3) {
                if any(node + vec3(ax, ay, az), |n| n == 0) {
                    continue;
                }
                let other_node = vec3(node.x + ax - 1, node.y + ay - 1, node.z + az - 1);
                let other_node_pos: Vec3 =
                    self.params.h * other_node.cast().expect("Failed to cast usize to f32");
                let idx = grid.coord_to_index(other_node);
                // TODO: Don't just ignore particles that get out of bounds
                if idx >= velocities.len() {
                    continue;
                }

                let weight_grad = weight_grad(other_node_pos, self.position[p], self.params.h);
                force[idx] += force_contrib * weight_grad;
            }
        }

        // Step 5: Grid Velocity update
        velocities.iter_mut().enumerate().for_each(|(i, v)| {
            if grid.mass[i] != 0. {
                let delta_v = self.params.delta_time * force[i] / grid.mass[i];
                *v += delta_v;

                let coord = grid.index_to_coord(i);
                // Enforce boundary conditions
                // NOTE: As far as I can tell, these grid level boundary conditions just
                // add unnecessary instability to the system, the particle ones work
                // just as well
                //if coord.y < 3 {
                //v.y = 0.;
                //}
                //if coord.x < 3 {
                //v.x = 0.
                //}
                //if coord.z < 3 {
                //v.z = 0.
                //}

                //if coord.x >= grid.size.x - 2 {
                //v.x = 0.;
                //}
                //if coord.y >= grid.size.y - 2 {
                //v.y = 0.;
                //}
                //if coord.z >= grid.size.z - 2 {
                //v.z = 0.;
                //}
            }
        });

        // Step 6: Update Particle Deformation Gradient
        for p in 0..self.params.num_particles {
            let mut fact = Matrix3::identity();
            let node: Coord = unwrap_or_continue!((self.position[p] / self.params.h).cast());

            for (ax, ay, az) in iproduct!(0..3, 0..3, 0..3) {
                if any(node + vec3(ax, ay, az), |n| n == 0) {
                    continue;
                }
                let other_node = vec3(node.x + ax - 1, node.y + ay - 1, node.z + az - 1);
                let other_node_pos: Vec3 = self.params.h * unwrap_or_continue!(other_node.cast()); // TODO: GET RID OF THIS
                let idx = grid.coord_to_index(other_node);
                // TODO: Don't just ignore particles that get out of bounds
                if idx >= velocities.len() {
                    continue;
                }

                let velocity = velocities[idx];
                let weight = weight_grad(other_node_pos, self.position[p], self.params.h);
                fact += self.params.delta_time * outer(velocity, weight);
            }

            self.deformation_gradient[p] = fact * self.deformation_gradient[p];
        }

        //
        // Step 7: Grid to Particle Transfer
        for p in 0..self.params.num_particles {
            let node: Coord = unwrap_or_continue!((self.position[p] / self.params.h).cast());

            //self.affine_matrices[p] = Matrix3::zero();
            self.velocity[p] = Vec3::zero();

            for (ax, ay, az) in iproduct!(0..3, 0..3, 0..3) {
                if any(node + vec3(ax, ay, az), |n| n == 0) {
                    continue;
                }
                let other_node = vec3(node.x + ax - 1, node.y + ay - 1, node.z + az - 1);
                let other_node_pos: Vec3 = self.params.h * unwrap_or_continue!(other_node.cast());
                let idx = grid.coord_to_index(other_node);

                let weight = weight(other_node_pos, self.position[p], self.params.h);

                //let diff = other_node_pos - self.position[p];

                // TODO: Don't just ignore particles that get out of bounds
                if idx >= velocities.len() {
                    continue;
                }
                self.velocity[p] += weight * velocities[idx];

                //self.affine_matrices[p] += weight * outer(velocities[idx], diff);
            }
            self.position[p] += self.params.delta_time * self.velocity[p];

            let bounds_min = self.params.bounds.start();
            let bounds_max = self.params.bounds.end();

            let position = &mut self.position[p];
            //let velocity = &mut self.velocity[p];

            (0..3).for_each(|i| {
                if position[i] < bounds_min[i] {
                    //velocity[i] *= -0.5;
                    position[i] = bounds_min[i];
                }

                if position[i] > bounds_max[i] {
                    //velocity[i] *= -0.5;
                    position[i] = bounds_max[i];
                }
            });

            // Step 8: Particle Advection

            //let pos = self.position[p];
            let vel = self.velocity[p].magnitude2();
            verts.push(Vertex {
                position: self.position[p].into(),
                color: [vel, 0.5 * vel + 0.5, 1.],
            });
        }

        self.frame_counter += 1;
        //println!("Finished Frame {}", self.frame_counter);

        verts
    }
}

pub struct MpmParmaters {
    /// The total number of Lagrangian particles in the simulation
    pub num_particles: usize,
    /// The radius of each particle
    pub h: Scalar,
    /// The bounds of the simulation.
    pub bounds: Range<Vec3>,
    /// The size of the time step. Larger time steps will simulate faster, but may be unstable or
    /// innaccurate.
    pub delta_time: f32,
}

impl Default for MpmParmaters {
    fn default() -> Self {
        MpmParmaters {
            num_particles: 0,
            h: 0.05,
            bounds: Vec3::zero()..=vec3(0.5, 1.5, 2.),
            delta_time: 0.001,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_outer_product() {
        let v0 = Vec3::new(1., 3., 5.);
        let m = outer(v0, v0);
        assert_eq!(m, Matrix3::new(1., 3., 5., 3., 9., 15., 5., 15., 25.));
    }
}
