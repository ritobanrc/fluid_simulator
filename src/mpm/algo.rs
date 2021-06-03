use super::grid::{Coord, MpmGrid};
use super::MpmSimulation;

pub fn particle_to_grid(s: &MpmSimulation) -> MpmGrid {
    let grid_bounds: Coord = (s.params.bounds / s.params.h)
        .cast()
        .expect("Failed to cast f32 to usize");

    let mut grid = MpmGrid::new(grid_bounds.x, grid_bounds.y, grid_bounds.z);
}
