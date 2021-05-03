pub(crate) mod grid;
pub(crate) mod kernels;

use cgmath::prelude::*;
use cgmath::{point3, vec3, Point3, Vector3};
use grid::{Coord, Grid};

use crate::{Scalar, Vec3};

pub struct Simulation {
    pub masses: Vec<Scalar>,
    pub positions: Vec<Point3<Scalar>>,
    pub velocities: Vec<Vec3>,
    pub force: Vec<Vec3>,
    pub grid: Grid,
    pub h: Scalar,
}

impl Simulation {
    pub fn new(h: Scalar, bounds: Vec3) -> Self {
        Simulation {
            masses: Vec::new(),
            positions: Vec::new(),
            velocities: Vec::new(),
            force: Vec::new(),
            grid: Grid::new(
                (bounds.x / h) as usize,
                (bounds.y / h) as usize,
                (bounds.z / h) as usize,
            ),
            h,
        }
    }

    pub(crate) fn position_to_coord(&self, pos: Vec3) -> Coord {
        pos.map(|i| (i / self.h) as usize)
    }

    pub(crate) fn coord(&self, index: usize) -> Coord {
        self.position_to_coord(self.positions[index].to_vec())
    }

    pub(crate) fn add_particle(&mut self, position: Point3<Scalar>) {
        let index = self.masses.len();
        self.masses.push((10. * self.h).powi(3));
        self.positions.push(position);
        self.velocities.push(Vector3::zero());
        self.force.push(Vector3::zero());

        self.grid
            .add_particle(self.position_to_coord(position.to_vec()), index);
    }
}
