use cgmath::{vec3, Vector3, Zero};

use crate::{Scalar, Vec3};

/// A 3d Coordinate composed of 3 integers.
pub type Coord = Vector3<usize>;

pub struct MpmGrid {
    pub(crate) mass: Vec<Scalar>,
    pub(crate) momentum: Vec<Vec3>,
    pub(crate) bounds: Coord,
}

impl MpmGrid {
    pub fn new(width: usize, height: usize, depth: usize) -> Self {
        let length = width * height * depth;
        MpmGrid {
            mass: vec![0.; length],
            momentum: vec![Vec3::zero(); length],
            bounds: vec3(width, height, depth),
        }
    }

    pub fn coord_to_index(&self, i: Coord) -> usize {
        i.x + self.bounds.x * i.y + self.bounds.x * self.bounds.y * i.z
    }

    pub fn index_to_coord(&self, mut i: usize) -> Coord {
        let z = i / (self.bounds.x * self.bounds.y);
        i -= z * self.bounds.x * self.bounds.y;
        let y = i / self.bounds.x;
        let x = i % self.bounds.x;
        vec3(x, y, z)
    }

    pub(crate) fn mass(&self, c: Coord) -> Option<Scalar> {
        self.mass.get(self.coord_to_index(c)).copied()
    }

    pub(crate) fn mass_mut(&mut self, c: Coord) -> Option<&mut Scalar> {
        let idx = self.coord_to_index(c);
        self.mass.get_mut(idx)
    }

    pub(crate) fn momentum(&self, c: Coord) -> Option<Vec3> {
        self.momentum.get(self.coord_to_index(c)).copied()
    }

    pub(crate) fn momentum_mut(&mut self, c: Coord) -> Option<&mut Vec3> {
        let idx = self.coord_to_index(c);
        self.momentum.get_mut(idx)
    }
}
