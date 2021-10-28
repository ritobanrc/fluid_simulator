use crate::{Scalar, Vec3};
use std::ops::Range;

#[cfg(test)]
pub fn diff_test<F, DF>(f: F, df: DF, domain: Range<Vec3>, eps: Scalar)
where
    F: Fn(Vec3) -> Scalar,
    DF: Fn(Vec3) -> Vec3,
{
    for _ in 0..1000 {
        let x0 = Vec3::new_random().component_mul(&domain.size()) - domain.start;
        let dx = (2. * Vec3::new_random() - Vec3::ones()) * eps;

        let x1 = x0 + dx;

        let a0 = f(x0);
        let a1 = f(x1);

        let d0 = df(x0);
        let d1 = df(x1);

        let u = (a1 - a0) / eps;
        let v = (d0 + d1).dot(&dx) / (2. * eps);

        let err = (u - v).abs();
        let rel_err = err / (Scalar::max(u.abs(), v.abs()));

        if rel_err > 100. * eps {
            eprintln!(
                "Diff Test Failed: {:?} -- {:?} vs {:?} at {:?}",
                err, u, v, x0
            );
        }
    }
}

pub trait RangeExt {
    fn size(&self) -> Vec3;

    fn center(&self) -> Vec3;

    fn contains_point(&self, x: &Vec3) -> bool;

    fn thickened(&self, amount: Scalar) -> Self;
}

impl RangeExt for Range<Vec3> {
    fn size(&self) -> Vec3 {
        self.end - self.start
    }

    fn center(&self) -> Vec3 {
        0.5 * (self.start + self.end)
    }

    fn contains_point(&self, x: &Vec3) -> bool {
        self.start.all_lt(x) && self.end.all_gt(x)
    }

    fn thickened(&self, amount: Scalar) -> Self {
        self.start - Vec3::from_element(amount)..self.end + Vec3::from_element(amount)
    }
}

pub trait VecExt {
    fn all_lt(&self, other: &Self) -> bool;

    fn all_gt(&self, other: &Self) -> bool;

    fn component_max(&self, other: &Self) -> Self;

    fn ones() -> Self;
}

impl VecExt for Vec3 {
    fn all_lt(&self, other: &Self) -> bool {
        self.x < other.x && self.y < other.y && self.z < other.z
    }

    fn all_gt(&self, other: &Self) -> bool {
        self.x > other.x && self.y > other.y && self.z > other.z
    }

    fn component_max(&self, other: &Self) -> Self {
        self.zip_map(other, |a, b| a.max(b))
    }

    fn ones() -> Self {
        Self::from_element(1.)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diff_test() {
        diff_test(
            |x| x.x.sin() + x.y * x.y + x.z.exp(),
            |x| Vec3::new(x.x.cos(), 2. * x.y, x.z.exp()),
            Vec3::from_element(-10.)..Vec3::from_element(10.),
            1e-6,
        );
    }
}
