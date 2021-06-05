use super::grid::{Coord, MpmGrid};
use super::MpmSimulation;

pub fn compute_regions(s: &MpmSimulation) {
    use super::math::RangeExt;
    let particle_regions = s
        .position
        .iter()
        .map(|x| s.grid.particle_range(*x))
        .collect::<Vec<_>>();

    let particle_cells = particle_regions
        .iter()
        .map(|range| s.grid.position_to_coord(*range.start()));

    let _ = particle_regions
        .iter()
        .map(|region| {
            s.position
                .iter()
                .filter(move |&&x| region.contains_point(x))
        })
        .collect::<Vec<_>>();

    todo!()
}

pub fn particle_to_grid(s: &MpmSimulation) -> MpmGrid {
    let mut grid = MpmGrid::new(&s.params);

    grid
}
