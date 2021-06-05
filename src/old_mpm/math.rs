use crate::{Scalar, Vec3};
use cgmath::{vec3, Array, Matrix3, Vector3};

/// The kernel function N(x) described in Eqn. (122), pg. 33
/// MPM SIGGRAPH Course Notes 2016
pub fn kernel(x: Scalar) -> Scalar {
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
/// Derivation:
/// N(x) = 1/2*|x|^3 - |x|^2 + 2/3    0 <= |x| < 1
///        1/6(2 - |x|)^3             1 <= |x| < 2
///        0                          2 <= |x|
///
/// N'(x) = 3/2 |x|^2  - 2|x|        0 <= |x| < 1
///         -1/2(2 - |x|)^2          1 <= |x| < 2
///          0                       2 <= |x|
pub fn kernel_grad(x: Scalar) -> Scalar {
    let x = x.abs();
    if 0. <= x && x < 1. {
        let x2 = x * x;

        1.5 * x2 - 2. * x
    } else if 1. <= x && x < 2. {
        let a = 2. - x;
        -0.5 * a * a
    } else {
        0.
    }
}

pub fn weight(grid_pos: Vec3, particle_pos: Vec3, h: Scalar) -> Scalar {
    (particle_pos - grid_pos).map(|x| kernel(x / h)).product()
}

pub fn weight_grad(grid_pos: Vec3, particle_pos: Vec3, h: Scalar) -> Vec3 {
    let arg = (particle_pos - grid_pos) / h;
    vec3(
        1. / h * kernel_grad(arg.x) * kernel(arg.y) * kernel(arg.z),
        kernel(arg.x) * 1. / h * kernel_grad(arg.y) * kernel(arg.z),
        kernel(arg.x) * kernel(arg.y) * 1. / h * kernel_grad(arg.z),
    )
}

pub fn outer(a: Vec3, b: Vec3) -> Matrix3<Scalar> {
    Matrix3::from_cols(
        vec3(a.x * b.x, a.x * b.y, a.x * b.z),
        vec3(a.y * b.x, a.y * b.y, a.y * b.z),
        vec3(a.z * b.x, a.z * b.y, a.z * b.z),
    )
}

use std::ops::RangeInclusive;

pub trait RangeExt {
    fn contains_point(&self, x: Vec3) -> bool;
}

impl RangeExt for RangeInclusive<Vec3> {
    fn contains_point(&self, a: Vec3) -> bool {
        self.start().x <= a.x
            && self.start().y <= a.y
            && self.start().z <= a.z
            && self.end().x >= a.x
            && self.end().y >= a.y
            && self.end().z >= a.z
    }
}

pub fn any<T, F: Fn(T) -> bool>(v: Vector3<T>, f: F) -> bool {
    f(v.x) || f(v.y) || f(v.z)
}

pub fn all<T, F: Fn(T) -> bool>(v: Vector3<T>, f: F) -> bool {
    f(v.x) && f(v.y) && f(v.z)
}
