use na::Vector3;

use super::{MpmParameters, Scalar, Vec3};
use std::ops::Range;

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
            pub fn $field_name(&self, coord: Vector3<usize>) -> Option<&$type> {
                self.$field_name.get(self.data.coord_to_index(coord))
            }

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
    pub fn new(params: &MpmParameters) -> Self {
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

    pub fn compute_forces(&mut self) {
        for i in 0..self.data.num_cells {
            // Just gravity, for now. Fg = -mg
            self.force[i] = self.mass[i] * Vec3::new(0., -1., 0.);
            // TODO: Add a constitutive model
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
    //pub one_over_h: Scalar,
    /// The bounds of the domain spanned by the grid.
    pub bounds: Range<Vec3>,
}

impl GridData {
    /// Creates a new `GridData`.
    /// Note that the `bounds` in the returned `GridData` may not be the same as the `bounds` passed into the function, depending on rounding.
    fn new(h: Scalar, bounds: Range<Vec3>) -> GridData {
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
            //one_over_h: 1. / h,
            bounds: bounds.start..actual_bounds,
        }
    }

    /// Converts the position in world space to a grid location. Note that this function does
    /// not explicitly round, so avoid passing in coordinates that are not exactly on grid locations.
    fn pos_to_coord(&self, pos: Vec3) -> Result<Vector3<usize>, ParticleGridError> {
        use num::ToPrimitive;
        let coord_f32 = (pos - self.bounds.start) / self.h;

        let mut coord_usize = Vector3::zeros();

        for (i, &v) in coord_f32.into_iter().enumerate() {
            // TODO: Fix this abomination of a function, converting floats to integers
            // shouldn't be this complicated
            coord_usize[i] = match v.to_usize() {
                Some(v) => v,
                None => {
                    return Err(ParticleGridError::FloatToUsize {
                        error: format!("Failed to convert position to coordinate: {:?}", pos),
                        float: v,
                    })
                }
            };
        }

        Ok(coord_usize)
    }

    pub fn coord_to_pos(&self, coord: Vector3<usize>) -> Vec3 {
        coord.cast::<Scalar>() * self.h + self.bounds.start
    }

    pub fn particle_neighborhood(&self, p: Vec3) -> Range<Vector3<usize>> {
        // transform the point into "grid space"
        assert!(
            self.bounds.contains_point(p),
            "Particle {:?} not in Grid Bounds {:?}",
            p,
            self.bounds
        );

        let grid_space_pos = (p - self.bounds.start) / self.h;

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

    pub fn particle_grid_iterator(
        &self,
        p: Vec3,
    ) -> impl Iterator<Item = (Vector3<usize>, Scalar)> + '_ {
        NeighborhoodIter::new(self.particle_neighborhood(p)).map(move |i| {
            let grid_pos = i.cast::<Scalar>() * self.h + self.bounds.start;

            // TODO: Migrate entirely to f64 to avoid the casts
            let weight = self.weight((p - grid_pos).cast());

            (i, weight as f32)
        })
    }

    fn weight(&self, v: Vector3<f64>) -> f64 {
        v.into_iter().map(|x| kernel(x / self.h as f64)).product()
    }

    fn weight_grad(&self, v: Vector3<f64>) -> Vector3<f64> {
        // TODO: Add derivative tests to this
        let mut grad = Vector3::zeros();

        for i in 0..3 {
            grad[i] = 1. / self.h as f64
                * (kernel_derivative(v[i]) * kernel(v[(i + 1) % 3]) * kernel(v[(i + 2) % 3]));
        }
        grad
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

#[derive(Debug, Clone)]
enum ParticleGridError {
    FloatToUsize { error: String, float: f32 },
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

struct NeighborhoodIter {
    neighborhood: Range<Vector3<usize>>,
    current: Vector3<usize>,
}

impl NeighborhoodIter {
    fn new(neighborhood: Range<Vector3<usize>>) -> Self {
        let start = neighborhood.start;
        NeighborhoodIter {
            neighborhood,
            current: start,
        }
    }
}

impl Iterator for NeighborhoodIter {
    type Item = Vector3<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        self.current.x += 1;

        // TODO: Make this algorithm work for arbitrary dimensions
        if self.current.x > self.neighborhood.end.x {
            self.current.x = self.neighborhood.start.x;
            self.current.y += 1;
            if self.current.y > self.neighborhood.end.y {
                self.current.y = self.neighborhood.start.y;
                self.current.z += 1;
                if self.current.z > self.neighborhood.end.z {
                    return None;
                }
            }
        }

        Some(self.current)
    }
}

/// The cubic kernel function N(x) described in Eqn. (122), pg. 33
/// MPM SIGGRAPH Course Notes 2016
pub fn kernel(x: f64) -> f64 {
    let x = x.abs();
    if 0. <= x && x < 1. {
        let x3 = x * x * x;
        let x2 = x * x;

        0.5 * x3 - x2 + 2. / 3.
    } else if 1. <= x && x < 2. {
        let a = 2. - x;
        1. / 6. * a * a * a
    } else {
        0.
    }
}

/// The derivative of the cubic kernel function, N'(x).
///
/// N(x) = 1/2*|x|^3 - |x|^2 + 2/3    0 <= |x| < 1
///        1/6(2 - |x|)^3             1 <= |x| < 2
///        0                          2 <= |x|
///
/// N'(x) =             3/2 |x|^2  - 2|x|       0 <= |x| < 1
///         abs'(x)    -1/2(2 - |x|)^2          1 <= |x| < 2
///                     0                       2 <= |x|
///
/// Where abs' is the derivative of the absolute value function, because chain rule.
pub fn kernel_derivative(x: f64) -> f64 {
    let chain_rule = if x < 0. { -1. } else { 1. };
    let x = x.abs();

    chain_rule
        * if 0. <= x && x < 1. {
            let x2 = x * x;

            1.5 * x2 - 2. * x
        } else if 1. <= x && x < 2. {
            let a = 2. - x;
            -0.5 * a * a
        } else {
            0.
        }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn index_coord_test(i in 0usize..8000) {
            let grid = GridData::new(0.1, Vector3::new(-1., -1., -1.)..Vector3::new(1., 1., 1.));
            let coord = grid.index_to_coord(i);
            let index = grid.coord_to_index(coord);

            assert_eq!(index, i);
        }
    }

    #[derive(Debug)]
    struct DerivativeError {
        x: f64,
        estimated: f64,
        actual: f64,
        error: f64,
    }

    /// Tests if z_prime is the derivative of z, sampling `num_trials` points within the
    /// domain. Returns the Root-Mean-Square-Error if `Ok`, returns a `DerivativeError`
    /// otherwise
    fn test_derivative<F1: Fn(f64) -> f64, F2: Fn(f64) -> f64>(
        z: F1,
        z_prime: F2,
        domain: Range<f64>,
        num_trials: usize,
    ) -> Result<f64, DerivativeError> {
        let delta_x = 1e-5;

        let mut sum_square_error = 0.;

        for x in crate::linspace(domain.start, domain.end, num_trials) {
            let estimated = (z(x + delta_x) - z(x - delta_x)) / (2. * delta_x);

            let actual = (z_prime(x + delta_x) + z_prime(x - delta_x)) / 2.;

            let error = estimated - actual;
            if error.abs() >= 0.01 {
                return Err(DerivativeError {
                    x,
                    estimated,
                    actual,
                    error,
                });
            }
            sum_square_error += error * error;
        }

        Ok((sum_square_error / num_trials as f64).sqrt())
    }

    #[test]
    fn test_derivative_test() {
        assert!(test_derivative(f64::sin, f64::cos, -2. ..2., 100).is_ok());
        assert!(test_derivative(f64::sin, |x| 2. * f64::cos(x), -2. ..2., 100).is_err());

        assert!(test_derivative(|x| x * x, |x| 2. * x, -2. ..2., 100).is_ok());
        assert!(test_derivative(|x| x * x, |x| 3. * x, -2. ..2., 100).is_err());

        assert!(test_derivative(f64::exp, f64::exp, -2. ..2., 100).is_ok());
        assert!(test_derivative(f64::exp, f64::exp2, -2. ..2., 100).is_err());
    }

    #[test]
    fn test_kernel_derivative() {
        assert!(test_derivative(kernel, kernel_derivative, -2. ..2., 100).is_ok());
        assert!(test_derivative(kernel, |x| 2. * kernel_derivative(x), -2. ..2., 100).is_err());
    }
}
