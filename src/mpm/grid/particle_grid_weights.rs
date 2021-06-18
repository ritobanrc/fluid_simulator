use super::data::GridData;
use super::{MpmParameters, Scalar, Vec3};
use na::Vector3;
use std::ops::Range;

impl GridData {
    pub fn particle_grid_iterator(
        &self,
        p: Vec3,
    ) -> impl Iterator<Item = (Vector3<usize>, Scalar)> + '_ {
        NeighborhoodIter::new(self.particle_neighborhood(p)).map(move |i| {
            let grid_pos = self.coord_to_pos(i);

            // TODO: Migrate entirely to f64 to avoid the casts
            let weight = self.weight((p - grid_pos).cast());

            (i, weight as f32)
        })
    }

    pub fn particle_grid_iterator_grad(
        &self,
        p: Vec3,
    ) -> impl Iterator<Item = (Vector3<usize>, Vec3)> + '_ {
        NeighborhoodIter::new(self.particle_neighborhood(p)).map(move |i| {
            let grid_pos = self.coord_to_pos(i);

            // TODO: Migrate entirely to f64 to avoid the casts
            let weight_grad = self.weight_grad((p - grid_pos).cast());

            (i, weight_grad.cast())
        })
    }

    fn weight(&self, v: Vector3<f64>) -> f64 {
        v.into_iter()
            .map(|x| kernel(x * self.one_over_h as f64))
            .product()
    }

    fn weight_grad(&self, mut v: Vector3<f64>) -> Vector3<f64> {
        v *= self.one_over_h as f64;

        let k = v.map(kernel);
        let kd = v.map(kernel_derivative);

        let grad = self.one_over_h as f64
            * Vector3::new(kd.x * k.y * k.z, k.x * kd.y * k.z, k.x * k.y * kd.z);

        grad
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

    #[test]
    fn test_weight_gradient() {
        let h = 1.;
        let domain = Vector3::repeat(-2. * h)..Vector3::repeat(2. * h);
        let num_steps = 10;

        // the bounds shouldn't have any effect on the weight or its gradient
        let grid = GridData::new(h as f32, Vector3::zeros()..Vector3::zeros());

        for y in crate::linspace(domain.start.y, domain.end.y, num_steps) {
            for z in crate::linspace(domain.start.z, domain.end.z, num_steps) {
                let result = test_derivative(
                    |x| grid.weight(Vector3::new(x, y, z)),
                    |x| grid.weight_grad(Vector3::new(x, y, z)).x,
                    domain.start.x..domain.end.x,
                    10,
                );

                assert!(result.is_ok(), "∂/∂x failed. Error: {:#?}", result);
            }
        }

        for x in crate::linspace(domain.start.x, domain.end.x, num_steps) {
            for z in crate::linspace(domain.start.z, domain.end.z, num_steps) {
                let result = test_derivative(
                    |y| grid.weight(Vector3::new(x, y, z)),
                    |y| grid.weight_grad(Vector3::new(x, y, z)).y,
                    domain.start.y..domain.end.y,
                    10,
                );

                assert!(result.is_ok(), "∂/∂y failed. Error: {:#?}", result);
            }
        }

        for x in crate::linspace(domain.start.x, domain.end.x, num_steps) {
            for y in crate::linspace(domain.start.y, domain.end.y, num_steps) {
                let result = test_derivative(
                    |z| grid.weight(Vector3::new(x, y, z)),
                    |z| grid.weight_grad(Vector3::new(x, y, z)).z,
                    domain.start.z..domain.end.z,
                    10,
                );

                assert!(result.is_ok(), "∂/∂z failed. Error: {:#?}", result);
            }
        }
    }
}
