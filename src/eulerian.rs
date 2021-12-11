/// Implements a purely Eulerian fluid simulator, similar to the one described in Stam 03 or
/// Bridson 08.
mod conjugate_gradient;

use crate::util::{RangeExt, VecExt};
use nalgebra_sparse::{CooMatrix, CsrMatrix};
use serde::{Deserialize, Serialize};
use std::ops::Range;

use crate::math::*;

#[derive(Default)]
struct ArrayNd<T> {
    data: Vec<T>,
    dimension: UV,
    stride: UV,
}

impl<T> ArrayNd<T>
where
    T: num::Zero + Clone,
{
    fn zeros(dim: UV) -> Self {
        let size = dim.iter().product();
        Self {
            data: vec![num::Zero::zero(); size],
            dimension: dim,
            stride: calculate_strides(dim),
        }
    }

    fn zeros_like<U>(other: &ArrayNd<U>) -> Self {
        Self {
            data: vec![num::Zero::zero(); other.data.len()],
            dimension: other.dimension,
            stride: other.stride,
        }
    }
}

pub(crate) fn calculate_strides(cells: UV) -> UV {
    let mut strides = UV::zeros();
    strides[DIM - 1] = 1;
    for i in (1..=DIM - 1).rev() {
        strides[i - 1] = strides[i] * cells[i];
    }
    strides
}

impl<T> std::ops::Index<UV> for ArrayNd<T> {
    type Output = T;
    fn index(&self, idx: UV) -> &Self::Output {
        debug_assert!(idx.all_lt(&self.dimension));
        &self.data[idx.dot(&self.stride)]
    }
}

impl<T> std::ops::IndexMut<UV> for ArrayNd<T> {
    fn index_mut(&mut self, idx: UV) -> &mut Self::Output {
        debug_assert!(idx.all_lt(&self.dimension));
        &mut self.data[idx.dot(&self.stride)]
    }
}

#[derive(Default)]
struct MacArray<T>([ArrayNd<T>; DIM]);

impl<T> std::ops::Index<FaceIndex> for MacArray<T> {
    type Output = T;
    fn index(&self, fi: FaceIndex) -> &Self::Output {
        &self.0[fi.axis][fi.cell]
    }
}

impl<T> std::ops::IndexMut<FaceIndex> for MacArray<T> {
    fn index_mut(&mut self, fi: FaceIndex) -> &mut Self::Output {
        &mut self.0[fi.axis][fi.cell]
    }
}

struct EulerianSimulation {
    grid: Grid,
    density: ArrayNd<T>,
    pressure: ArrayNd<T>,
    velocity: MacArray<T>,
    params: EulerianParams,
}

#[derive(Debug, Clone)]
struct EulerianParams {
    delta_time: T,
    rho: T,
}

impl crate::Simulation for EulerianSimulation {
    type Parameters = EulerianParams;

    fn new(params: Self::Parameters) -> Self {
        Self {
            grid: Grid::new(UV::from_element(10), TV::zeros()..TV::from_element(10.)),
            density: Default::default(),
            pressure: Default::default(),
            velocity: Default::default(),
            params,
        }
    }

    fn add_particle(&mut self, mass: crate::Scalar, position: crate::Vec3, velocity: crate::Vec3) {
        panic!("this makes no sense fix it")
    }

    fn simulate_frame(&mut self) -> Vec<crate::Vertex> {
        // advect density
        // advect other things (temperature, level set, etc.)
        // advect velocity
        // apply gravity
        // apply viscosity
        // build poisson matrix
        // perform pressure projection
        Vec::new()
    }
}

impl EulerianSimulation {
    fn discrete_divergence(&self) -> ArrayNd<T> {
        // TODO: this feels like C code, I should just be able to use a rust iterator and collect
        // into an ArrayNd.
        //
        // TODO: handle boundaries
        let mut divergence = ArrayNd::zeros(self.grid.cells);
        for cell in self.grid.cells() {
            divergence[cell] = (0..DIM)
                .map(|axis| {
                    // TODO: this should be in the grid class
                    let first = FaceIndex::new(cell, axis);
                    let second = FaceIndex::new(cell + UV::ith_axis(axis).into_inner(), axis);
                    self.grid.dx[axis] * (self.velocity[second] - self.velocity[first])
                })
                .sum();
        }
        divergence
    }

    fn build_poison_matrix(&self) -> CsrMatrix<T> {
        let num_cells = self.grid.num_cells().iter().product();
        // TODO: use only actual DOFs here
        let strides = calculate_strides(self.grid.num_cells());
        let idx = |cell: UV| cell.dot(&strides);
        let mut matrix = CooMatrix::new(num_cells, num_cells);

        let scale =
            self.params.delta_time / (self.params.rho * self.grid.dx.max() * self.grid.dx.max());

        let add_diag = |matrix: &mut CooMatrix<T>, cell, value| {
            let idx = idx(cell);
            // we (ab)use the fact that nalgerba sparse adds together duplicate entries when
            // converting to Csr/Csc format
            matrix.push(idx, idx, value)
        };
        let add_plus_axis = |matrix: &mut CooMatrix<T>, axis, cell, value| {
            let row = idx(cell);
            let col = idx(cell + UV::ith_axis(axis).into_inner());

            // TODO: the matrix is symmetric, store only one half of it
            matrix.push(row, col, value);
            matrix.push(col, row, value);
        };

        for cell in self.grid.cells() {
            // if cell is not fluid, just skip it
            // this should be done automatically when we use only actual dofs
            for axis in 0..DIM {
                add_diag(&mut matrix, cell, scale);
                // TODO: do these next two lines only `if cell.add_axis(axis, 1) == fluid`
                add_diag(&mut matrix, cell + UV::ith_axis(axis).into_inner(), scale);
                add_plus_axis(&mut matrix, axis, cell, -scale);
            }
        }

        CsrMatrix::from(&matrix)
    }

    fn apply_pressure(&mut self) {
        let scale = self.params.delta_time / self.params.rho * self.grid.dx.max();

        for cell in self.grid.cells() {
            for axis in 0..DIM {
                self.velocity[FaceIndex::new(cell, axis)] -= scale * self.pressure[cell];
                let neighbor = cell + UV::ith_axis(axis).into_inner();
                self.velocity[FaceIndex::new(neighbor, axis)] += scale * self.pressure[cell];
            }
        }
        // TODO: handle solids
    }
}

fn interpolate_node<Z>(grid: &Grid, z: &ArrayNd<Z>, x: TV) -> Z {
    todo!()
}

fn interpolate_cell<Z>(grid: &Grid, z: &ArrayNd<Z>, x: TV) -> Z {
    todo!()
}

fn interpolate_face<Z>(grid: &Grid, z: &MacArray<Z>, x: TV) -> na::SVector<Z, DIM> {
    todo!()
}

trait InterpolateFace<Z> {
    fn interpolate_face(grid: &Grid, z: &MacArray<Z>, x: TV) -> na::SVector<Z, DIM>;
    fn gradient_face(grid: &Grid, z: &MacArray<Z>, x: TV) -> na::SVector<Z, DIM>;
    fn hessian_face(grid: &Grid, z: &MacArray<Z>, x: TV) -> na::SVector<Z, DIM>;
}

trait InterpolationOutput {
    fn stencil_size() -> usize;
}

impl<T, const N: usize> InterpolationOutput for na::SVector<T, N> {
    fn stencil_size() -> usize {
        N
    }
}

trait InterpolationStencil {
    //const STENCIL_SIZE: usize;
    //type StencilSize: na::Dim;
    type Output: InterpolationOutput;
    fn value(t: T) -> Self::Output;
    fn diff(t: T) -> Self::Output;
    fn hess(t: T) -> Self::Output;
}

struct LinearStencil;

impl InterpolationStencil for LinearStencil {
    //const STENCIL_SIZE: usize = 2;
    type Output = na::Vector2<T>;
    fn value(t: T) -> Self::Output {
        na::vector![1. - t, t]
    }

    fn diff(t: T) -> Self::Output {
        na::vector![-1., 1.]
    }

    fn hess(t: T) -> Self::Output {
        na::vector![-1., 1.]
    }
}

impl<Z, Stencil: InterpolationStencil + Default> InterpolateFace<Z> for Stencil {
    fn interpolate_face(grid: &Grid, z: &MacArray<Z>, x: TV) -> na::SVector<Z, DIM> {
        let stencil = Stencil::default();
        let ss = Stencil::Output::stencil_size();
        let half_stencil_width = 0.5 * Stencil::Output::stencil_size() as T;
        let pos = (x - grid.domain.start).component_mul(&grid.one_over_dx)
            - TV::from_element(half_stencil_width);

        for axis in 0..DIM {
            let pos_axis = pos + TV::ith(axis, 0.5);
            let bottom_left_idx = pos.map(|x| x.floor() as usize); // NOTE: there should be a clamp here see physbam code
        }

        todo!()
    }

    fn gradient_face(grid: &Grid, z: &MacArray<Z>, x: TV) -> na::SVector<Z, DIM> {
        todo!()
    }
    fn hessian_face(grid: &Grid, z: &MacArray<Z>, x: TV) -> na::SVector<Z, DIM> {
        todo!()
    }
}

fn linear_stencil(t: T) -> na::SVector<T, 2> {
    na::vector![1. - t, t]
}

fn semi_lagrangian_advection<Z: num::Zero + Clone>(
    grid: &Grid,
    z: &ArrayNd<Z>,
    vel: &MacArray<T>,
    dt: T,
) -> ArrayNd<Z> {
    let mut z_out = ArrayNd::zeros_like(z);
    for cell in grid.cells() {
        let x = grid.cell_x(cell); // TODO: make a CellIterator struct with a reference to grid
        z_out[cell] = interpolate_cell(
            grid,
            z,
            euler_integrate(x, &|x| interpolate_face(grid, vel, x), -dt),
        );
    }
    z_out
}

fn euler_integrate<V: Fn(TV) -> TV>(x: TV, vel: &V, dt: T) -> TV {
    x + dt * vel(x)
}

fn rk2_integrate<V: Fn(TV) -> TV>(x: TV, vel: &V, dt: T) -> TV {
    let x1 = euler_integrate(x, vel, 0.5 * dt);
    x + dt * vel(x1)
}

/// Represents a Grid (Co-located or Staggered) in world space.
#[derive(Serialize, Deserialize)]
pub struct Grid {
    /// Represents the world-space domain of the grid.
    pub domain: Range<TV>,
    /// The number of cells in each direction.
    /// Should always by domain.size() / dx, componentwise.
    pub cells: UV,
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
    pub fn new(cells: UV, domain: Range<TV>) -> Self {
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
    pub fn num_nodes(&self) -> UV {
        self.cells + UV::from_element(1)
    }

    /// Returns the number of nodes in each direction
    /// This is one greater than the number of cells.
    pub fn num_cells(&self) -> UV {
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
    pub fn cell_center(&self, idx: UV) -> TV {
        self.domain.start
            + (na::convert::<_, TV>(idx) + TV::from_element(0.5)).component_mul(&self.dx)
    }

    /// Index of the nearest cell center
    pub fn cell_index(&self, x: TV) -> Option<UV> {
        na::try_convert::<_, UV>(
            (x - self.domain.start)
                .component_mul(&self.one_over_dx)
                .map(T::floor),
        )
    }

    /// Index of the node to the lower-left
    pub fn node_lower(&self, x: TV) -> Option<UV> {
        self.cell_index(x)
    }

    pub fn node_x(&self, node: UV) -> TV {
        self.domain.start + node.cast::<T>().component_mul(&self.dx)
    }

    pub fn cell_x(&self, cell: UV) -> TV {
        self.domain.start + (cell.cast::<T>() + TV::from_element(0.5)).component_mul(&self.dx)
    }

    pub fn nodes<'a>(&self) -> impl Iterator<Item = UV> + 'a {
        RangeIterator::new(UV::zeros()..self.num_nodes())
    }

    pub fn cells<'a>(&self) -> impl Iterator<Item = UV> + 'a {
        RangeIterator::new(UV::zeros()..self.num_cells())
    }

    pub fn faces<'a>(&'a self) -> impl Iterator<Item = FaceIndex> + 'a {
        (0..DIM).flat_map(move |axis| {
            let start = UV::zeros();
            let end = self.num_cells() + UV::ith_axis(axis).into_inner();
            RangeIterator::new(start..end).map(move |cell| FaceIndex { cell, axis })
        })
    }

    pub const fn binary_counts<'a>(&self) -> [UV; 1 << DIM] {
        [
            UV::new(0, 0, 0),
            UV::new(1, 0, 0),
            UV::new(0, 1, 0),
            UV::new(1, 1, 0),
            UV::new(0, 0, 1),
            UV::new(1, 0, 1),
            UV::new(0, 1, 1),
            UV::new(1, 1, 1),
        ]
    }
}

pub struct FaceIndex {
    pub cell: UV,
    pub axis: usize, // NOTE: use ndarray `Axis` here
}

impl FaceIndex {
    fn new(cell: UV, axis: usize) -> Self {
        Self { cell, axis }
    }
}

// TODO: handle ghost cells and all kinds of other crap
pub struct RangeIterator {
    index: UV,
    range: Range<UV>,
    done: bool,
}

impl RangeIterator {
    fn new(range: Range<UV>) -> Self {
        Self {
            index: range.start,
            range,
            done: false,
        }
    }
}

impl Iterator for RangeIterator {
    type Item = UV;

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
        let range = UV::new(0, 0, 0)..UV::new(2, 2, 2);
        let v = RangeIterator::new(range).collect::<Vec<_>>();
        println!("{:?}", v);
        assert_eq!(v.len(), 8);
    }
}
