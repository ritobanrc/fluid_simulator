use na::Vector3;

use super::{MpmParameters, Scalar, Vec3};
use std::ops::Range;

/// Stores the metadata associated with the `MpmGrid`
#[derive(Debug, Clone)]
pub struct GridData {
    /// The total number of cells in the grid
    pub num_cells: usize,
    /// The number of nodes in each direction
    pub size: Vector3<usize>,
    /// The grid spacing
    pub h: Scalar,
    //// Reciprocal of the grid spacing, for efficient computation
    pub one_over_h: Scalar,
    /// The bounds of the domain spanned by the grid.
    pub bounds: Range<Vec3>,
}

impl GridData {
    /// Creates a new `GridData`.
    /// Note that the `bounds` in the returned `GridData` may not be the same as the `bounds` passed into the function, depending on rounding.
    pub fn new(h: Scalar, bounds: Range<Vec3>) -> GridData {
        let size = (bounds.end - bounds.start) / h;
        let size = size.map(|x| x.ceil() as usize + 1);
        let num_cells = size.iter().product();

        println!(
            "Using Grid w/ Size: [{}, {}, {}] and {} cells",
            size.x, size.y, size.z, num_cells
        );

        let actual_bounds = bounds.start + h * size.cast::<Scalar>();

        GridData {
            num_cells,
            size,
            h,
            one_over_h: 1. / h,
            bounds: bounds.start..actual_bounds,
        }
    }

    pub fn coord_to_pos(&self, coord: Vector3<usize>) -> Vec3 {
        coord.cast::<Scalar>() * self.h + self.bounds.start
    }

    pub fn coord_to_index(&self, i: Vector3<usize>) -> usize {
        i.x + self.size.x * i.y + self.size.x * self.size.y * i.z
    }

    pub fn index_to_coord(&self, mut i: usize) -> Vector3<usize> {
        let z = i / (self.size.x * self.size.y);
        i -= z * self.size.x * self.size.y;
        let y = i / self.size.x;
        let x = i % self.size.x;
        Vector3::new(x, y, z)
    }

    /// Gets a 4x4 neighborhood of grid nodes around a point
    pub fn particle_neighborhood(&self, p: Vec3) -> Range<Vector3<usize>> {
        // transform the point into "grid space"
        assert!(
            self.bounds.contains_point(p),
            "Particle {:?} not in Grid Bounds {:?}",
            p,
            self.bounds
        );

        let grid_space_pos = (p - self.bounds.start) * self.one_over_h;

        // usize cast shouldn't panic because of assert above
        // TODO: Implement this for other degree kernel functions.
        //       For now, this is just degree 3 (i.e. the diamater of interpolating stencil is 4
        //       grid units).
        //
        //       PhysBAM uses the formula `min = floor(((x - 0.5 * degree * h) - start)/h),
        //                                 max = min + degree + 1`
        //
        //       See `PARTICLE_GRID_WEIGHTS_SPLINE.cpp` line 139, 202, and GRID.h line 177
        let min = grid_space_pos.map(|x| x.floor() as usize) - Vector3::from_element(1);
        let max = grid_space_pos.map(|x| x.ceil() as usize) + Vector3::from_element(1);

        min..max
    }

    pub(crate) fn coord_in_grid(&self, coord: Vector3<usize>) -> bool {
        for i in 0..3 {
            if coord[i] >= self.size[i] {
                return false;
            }
        }
        return true;
    }
}

trait ContainsPoint {
    fn contains_point(&self, a: Vec3) -> bool;
}

impl ContainsPoint for Range<Vec3> {
    /// Checks if the range contains the point. Note that function does _not_ treat it as an
    /// exclusive range, and does return true at `self.end`.
    fn contains_point(&self, a: Vec3) -> bool {
        self.start.x <= a.x
            && self.start.y <= a.y
            && self.start.z <= a.z
            && self.end.x >= a.x
            && self.end.y >= a.y
            && self.end.z >= a.z
    }
}
