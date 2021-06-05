use crate::mpm::MpmParmaters;
use cgmath::{vec3, Array, Vector3, Zero};
use std::ops::RangeInclusive as Range; // NOTE: All ranges should be inclusive

use crate::{Scalar, Vec3};

/// A 3d Coordinate composed of 3 integers.
pub type Coord = Vector3<usize>;

pub struct MpmGrid {
    pub(crate) mass: Vec<Scalar>,
    pub(crate) momentum: Vec<Vec3>,
    pub(crate) size: Coord,
    h: Scalar,
    bounds: Range<Vec3>,
}

impl MpmGrid {
    pub fn new(params: &MpmParmaters) -> Self {
        let size: Coord = ((params.bounds.end() - params.bounds.start()) / params.h)
            .cast()
            .expect(&format!(
                "Failed to cast Vec3 to Coord: {:?}",
                (params.bounds.start() - params.bounds.end()) / params.h
            ));

        let length = size.product();
        MpmGrid {
            mass: vec![0.; length],
            momentum: vec![Vec3::zero(); length],
            size,
            h: params.h, // NOTE: Consider cloning `params` or using `Rc`
            bounds: params.bounds.clone(),
        }
    }

    pub fn coord_to_index(&self, i: Coord) -> usize {
        i.x + self.size.x * i.y + self.size.x * self.size.y * i.z
    }

    pub fn index_to_coord(&self, mut i: usize) -> Coord {
        let z = i / (self.size.x * self.size.y);
        i -= z * self.size.x * self.size.y;
        let y = i / self.size.x;
        let x = i % self.size.x;
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

    /// The size of the `Range` returned by `particle_range`. Corresponds to `degree` in PhysBAM's
    /// `PARTICLE_GRID_WEIGHTS_SPLINE`
    const NEIGHBORHOOD_SIZE: Scalar = 2.;

    /// Returns a box of size NEIGHBORHOOD_SIZE around `pos`. See https://www.desmos.com/calculator/gwhf6t3lqs
    pub fn particle_range(&self, pos: Vec3) -> Range<Vec3> {
        let lower = pos - Vec3::from_value(MpmGrid::NEIGHBORHOOD_SIZE * self.h / 2.);
        let upper = lower + Vec3::from_value(MpmGrid::NEIGHBORHOOD_SIZE + 1.);

        lower..=upper
    }

    /// Given the lower point of a cell, returns the coordinate
    pub fn position_to_coord(&self, pos: Vec3) -> Option<Coord> {
        ((pos - self.bounds.start()) / self.h).cast()
    }
}
