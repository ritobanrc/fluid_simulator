use crate::util::*;
use crate::{Scalar, Vec3};
use std::ops::Range;

/// An object represented implicitly by a signed distance field.
pub trait ImplicitObject {
    fn signed_distance(&self, x: Vec3) -> Scalar;

    fn normal(&self, x: Vec3) -> Vec3;
}

impl ImplicitObject for Range<Vec3> {
    /// Credit to Ingilo Quielez.
    /// https://iquilezles.org/www/articles/distfunctions/distfunctions.htm
    fn signed_distance(&self, x: Vec3) -> Scalar {
        let p = x - self.center();
        let q = p.abs() - 0.5 * self.size();
        q.map(|x| x.max(0.0)).magnitude() + Scalar::min(q.max(), 0.0)
    }

    /// See PhysBAM Public_Library/Core/Math_Tools/RANGE.cpp for source.
    fn normal(&self, x: Vec3) -> Vec3 {
        let min_dists = self.start - x;
        let max_dists = x - self.end;

        if self.contains_point(&x) {
            let phi = min_dists.component_max(&max_dists);
            let (axis, dist) = phi.argmax();
            Vec3::ith(axis, dist.signum())
        } else {
            min_dists
                .zip_map(&max_dists, |min, max| {
                    if max > min && max > 0. {
                        max
                    } else if min > 0. {
                        -min
                    } else {
                        0.
                    }
                })
                .normalize()
        }
    }
}

pub struct Sphere {
    pub center: Vec3,
    pub radius: f64,
}

impl ImplicitObject for Sphere {
    fn signed_distance(&self, x: Vec3) -> Scalar {
        (x - self.center).magnitude() - self.radius
    }

    fn normal(&self, x: Vec3) -> Vec3 {
        (x - self.center).normalize()
    }
}

pub struct Invert<O>(pub O);

impl<O: ImplicitObject> ImplicitObject for Invert<O> {
    fn signed_distance(&self, x: Vec3) -> Scalar {
        -self.0.signed_distance(x)
    }

    fn normal(&self, x: Vec3) -> Vec3 {
        -self.0.normal(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_implicit_object() {
        let block = Vec3::new(-3., 4., -1.)..Vec3::new(6., 5., 7.);
        let sphere = Sphere {
            center: Vec3::from_element(3.),
            radius: 5.,
        };

        let objects: &[&dyn ImplicitObject] = &[&block, &sphere];

        for object in objects {
            diff_test(
                |x| object.signed_distance(x),
                |x| object.normal(x),
                Vec3::from_element(-10.)..Vec3::from_element(10.),
                1e-6,
            );
        }
    }
}
