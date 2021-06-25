mod data;
pub mod weights;

use na::Vector3;

use super::{MpmParameters, Scalar, Vec3};
use data::GridData;

/// Stores the grid data for the Mpm Simulation
pub struct MpmGrid {
    pub mass: Vec<Scalar>,
    pub velocity: Vec<Vec3>,
    pub momentum: Vec<Vec3>,
    pub force: Vec<Vec3>,
    pub data: GridData,
}

macro_rules! grid_impls {
    ($($field_name:ident | $mut_name:ident : $type:ty),*) => {
        $(impl MpmGrid {
            #[allow(dead_code)]
            pub fn $field_name(&self, coord: Vector3<usize>) -> Option<&$type> {
                self.$field_name.get(self.data.coord_to_index(coord))
            }

            #[allow(dead_code)]
            pub fn $mut_name(&mut self, coord: Vector3<usize>) -> Option<&mut $type> {
                self.$field_name.get_mut(self.data.coord_to_index(coord))
            }
        })*
    };
}

grid_impls!(
    mass | mass_mut: Scalar,
    velocity | velocity_mut: Vec3,
    momentum | momentum_mut: Vec3,
    force | force_mut: Vec3
);

impl MpmGrid {
    pub fn new<CM>(params: &MpmParameters<CM>) -> Self {
        let grid_bounds_start = params.bounds.start - Vec3::from_element(params.h);
        let grid_bounds_end = params.bounds.end + Vec3::from_element(params.h);

        let data = GridData::new(params.h, grid_bounds_start..grid_bounds_end);

        Self {
            mass: vec![0.; data.num_cells],
            velocity: vec![Vec3::zeros(); data.num_cells],
            momentum: vec![Vec3::zeros(); data.num_cells],
            force: vec![Vec3::zeros(); data.num_cells],
            data,
        }
    }

    /// Fills each of the arrays in the grid with zeros.
    pub fn clear_grid(&mut self) {
        self.mass.fill(0.);
        self.velocity.fill(Vector3::zeros());
        self.momentum.fill(Vector3::zeros());
        self.force.fill(Vector3::zeros());
    }

    pub fn total_mass(&self) -> Scalar {
        self.mass.iter().sum()
    }

    pub fn total_momentum(&self) -> Vec3 {
        self.momentum.iter().sum()
    }

    /// Fill the `velocity` array using the `momentum` array
    pub fn compute_velocities(&mut self) {
        for i in 0..self.data.num_cells {
            if self.mass[i] != 0. {
                self.velocity[i] = self.momentum[i] / self.mass[i];
            }
            // Note that we don't need to handle the `else` case because `velocity` has already been zeroed out
        }
    }

    pub fn velocity_update(&mut self, delta_time: Scalar) {
        for i in 0..self.data.num_cells {
            if self.mass[i] == 0. {
                continue;
            }

            let delta_v = delta_time * self.force[i] / self.mass[i];

            self.velocity[i] += delta_v;

            // Enforce boundary conditions
            // TODO: Abstract over different types of BCs
            let coord = self.data.index_to_coord(i);
            if coord.x <= 1 || coord.x >= self.data.size.x - 2 {
                self.velocity[i].x = 0.
            }
            if coord.y <= 1 || coord.y >= self.data.size.y - 2 {
                self.velocity[i].y = 0.
            }
            if coord.z <= 1 || coord.z >= self.data.size.z - 2 {
                self.velocity[i].z = 0.
            }
        }
    }
}
