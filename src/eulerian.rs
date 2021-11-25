/// Implements a purely Eulerian fluid simulator, similar to the one described in Stam 03 or
/// Bridson 08.
mod conjugate_gradient;

use crate::util::RangeExt;
use serde::{Deserialize, Serialize};
use std::ops::Range;

use crate::fem::math::*;

/// Represents a Grid (Co-located or Staggered) in world space.
#[derive(Serialize, Deserialize)]
pub struct Grid {
    /// Represents the world-space domain of the grid.
    pub domain: Range<TV>,
    /// The number of cells in each direction.
    /// Should always by domain.size() / dx, componentwise.
    pub cells: IV,
    /// Size of each grid cell
    /// Not serialized, since it can be calculated from `domain` and `cells`
    #[serde(skip)]
    pub dx: TV,
    /// Reciprocal of each component of `dx`
    /// Also not serialized, since it can be easily calculated from `dx`.
    #[serde(skip)]
    pub one_over_dx: TV,
}

impl Grid {
    /// Creates a new grid given the number of cells in each direction and the size of the domain.
    pub fn new(cells: IV, domain: Range<TV>) -> Self {
        let dx = domain.size().component_div(&na::convert::<_, TV>(cells));
        let one_over_dx = TV::from_element(1.).component_div(&dx);
        Self {
            domain,
            dx,
            one_over_dx,
            cells,
        }
    }

    /// Recalculates dx and one_over_dx using `domain` and `cells`.
    pub fn recalculate_dx(&mut self) {
        self.dx = self
            .domain
            .size()
            .component_div(&na::convert::<_, TV>(self.cells));
        self.one_over_dx = TV::from_element(1.).component_div(&self.dx);
    }

    /// Returns the number of nodes in each direction
    /// This is one greater than the number of cells.
    pub fn num_nodes(&self) -> IV {
        self.cells + IV::from_element(1)
    }

    /// Returns the number of nodes in each direction
    /// This is one greater than the number of cells.
    pub fn num_cells(&self) -> IV {
        self.cells
    }

    /// Returns the volume (in 3d) or area (in 2d) of a cell
    pub fn cell_size(&self) -> T {
        self.dx.iter().product()
    }

    /// Returns the area (in 3d) or length (in 2d) of each face
    pub fn face_areas(&self) -> TV {
        self.cell_size() * self.one_over_dx
    }

    /// Returns the area (in 3d) or length (in 2d) of faces in a particular direction.
    pub fn face_area(&self, axis: usize) -> T {
        self.cell_size() * self.one_over_dx[axis]
    }

    /// Location of the center of the cell.
    pub fn cell_center(&self, idx: IV) -> TV {
        self.domain.start
            + (na::convert::<_, TV>(idx) + TV::from_element(0.5)).component_mul(&self.dx)
    }

    /// Index of the nearest cell center
    pub fn cell_index(&self, x: TV) -> Option<IV> {
        na::try_convert::<_, IV>(
            (x - self.domain.start)
                .component_mul(&self.one_over_dx)
                .map(T::floor),
        )
    }

    /// Index of the node to the lower-left
    pub fn node_lower(&self, x: TV) -> Option<IV> {
        self.cell_index(x)
    }

    pub fn node_x(&self, node: IV) -> TV {
        self.domain.start + node.cast::<T>().component_mul(&self.dx)
    }

    pub fn nodes<'a>(&self) -> impl Iterator<Item = IV> + 'a {
        RangeIterator::new(IV::zeros()..self.num_nodes())
    }

    pub fn cells<'a>(&self) -> impl Iterator<Item = IV> + 'a {
        RangeIterator::new(IV::zeros()..self.num_cells())
    }

    pub const fn binary_counts<'a>(&self) -> [IV; 1 << DIM] {
        [
            IV::new(0, 0, 0),
            IV::new(1, 0, 0),
            IV::new(0, 1, 0),
            IV::new(1, 1, 0),
            IV::new(0, 0, 1),
            IV::new(1, 0, 1),
            IV::new(0, 1, 1),
            IV::new(1, 1, 1),
        ]
    }
}

pub(crate) fn calculate_strides(cells: IV) -> IV {
    let mut strides = IV::zeros();
    strides[DIM - 1] = 1;
    for i in (1..=DIM - 1).rev() {
        strides[i - 1] = strides[i] * cells[i];
    }
    strides
}

// TODO: handle ghost cells and all kinds of other crap
pub struct RangeIterator {
    index: IV,
    range: Range<IV>,
    done: bool,
}

impl RangeIterator {
    fn new(range: Range<IV>) -> Self {
        Self {
            index: range.start,
            range,
            done: false,
        }
    }
}

impl Iterator for RangeIterator {
    type Item = IV;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let result = self.index;

        let mut i = DIM - 1;
        loop {
            self.index[i] += 1;
            if self.index[i] >= self.range.end[i] {
                self.index[i] = self.range.start[i];
                if i != 0 {
                    i -= 1;
                } else {
                    self.done = true;
                    break;
                }
            } else {
                break;
            }
        }

        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn range_tests() {
        let range = IV::new(0, 0, 0)..IV::new(2, 2, 2);
        let v = RangeIterator::new(range).collect::<Vec<_>>();
        println!("{:?}", v);
        assert_eq!(v.len(), 8);
    }
}
