use itertools::iproduct;
use na::Vector3;
use smallvec::SmallVec;
use std::ops::{Index, IndexMut};

/// Represents a single grid cell. A grid cell contains a list of the particles within it.
///
/// A `SmallVec` is used to prevent unnecessary allocation.
type GridCell = SmallVec<[usize; 3]>;

/// A 3d Coordinate composed of 3 integers.
pub type Coord = Vector3<usize>;

/// The Grid structure used to speed up SPH neighbor finding.
pub struct Grid {
    pub(crate) grid: Vec<GridCell>,
    bounds: Coord,
}

impl Grid {
    pub fn new(width: usize, height: usize, depth: usize) -> Self {
        Grid {
            // TODO: Figure out if this is wasting a bunch of memory
            grid: vec![GridCell::new(); width * height * depth],
            bounds: Vector3::new(width, height, depth),
        }
    }

    fn coord_to_index(&self, i: Coord) -> usize {
        i.x + self.bounds.x * i.y + self.bounds.x * self.bounds.y * i.z
    }

    pub fn add_particle(&mut self, coord: Coord, index: usize) {
        self[coord].push(index);
    }

    pub fn get_neighbors(&self, coord: Coord) -> impl Iterator<Item = usize> + Clone + '_ {
        iproduct!(
            coord.x.saturating_sub(1)..=coord.x + 1,
            coord.y.saturating_sub(1)..=coord.y + 1,
            coord.z.saturating_sub(1)..=coord.z + 1
        )
        .filter_map(move |(x, y, z)| self.get(Vector3::new(x, y, z)))
        .flat_map(|cell| cell.iter().copied())
    }

    pub fn update_particle(&mut self, index: usize, old_coord: Coord, new_coord: Coord) {
        let old = self
            .get_mut(old_coord)
            .expect("Grid::update_particle old_coord not found");

        let index_in_coord = old
            .iter()
            .position(|&x| x == index)
            .expect("Grid::update_particle called with incorrect old_coord.");

        old.remove(index_in_coord);

        let new = self
            .get_mut(new_coord)
            .expect("Griid::update_particle new_coord out of bounds");

        new.push(index);
    }

    fn get(&self, i: Coord) -> Option<&GridCell> {
        if i.x >= self.bounds.x || i.y >= self.bounds.y || i.z >= self.bounds.z {
            return None;
        }
        let index = self.coord_to_index(i);
        Some(&self.grid[index])
    }

    #[allow(unused)]
    fn get_mut(&mut self, i: Coord) -> Option<&mut GridCell> {
        if i.x >= self.bounds.x || i.y >= self.bounds.y || i.z >= self.bounds.z {
            return None;
        }
        let index = self.coord_to_index(i);
        Some(&mut self.grid[index])
    }

    #[allow(dead_code)]
    pub(crate) fn measure_spilled(&self) -> usize {
        self.grid.iter().filter(|x| x.spilled()).count()
    }
}

impl Index<Coord> for Grid {
    type Output = GridCell;

    fn index(&self, i: Coord) -> &Self::Output {
        if i.x >= self.bounds.x || i.y >= self.bounds.y || i.z >= self.bounds.z {
            panic!("Attempted to get index out of bounds: {:?}", i);
        }
        let index = self.coord_to_index(i);
        &self.grid[index]
    }
}

impl IndexMut<Coord> for Grid {
    fn index_mut(&mut self, i: Coord) -> &mut Self::Output {
        if i.x >= self.bounds.x || i.y >= self.bounds.y || i.z >= self.bounds.z {
            panic!("Attempted to get index out of bounds: {:?}", i);
        }
        let index = self.coord_to_index(i);
        &mut self.grid[index]
    }
}
