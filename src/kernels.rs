use crate::{Scalar, Vec3};
use cgmath::prelude::*;
use cgmath::Vector3;

pub trait SmoothingKernel {
    fn value(_r: Vec3, _h: Scalar) -> Scalar {
        0.
    }

    fn gradient_mag(_r: Vec3, _h: Scalar) -> Scalar {
        0.
    }

    fn gradient(r: Vec3, h: Scalar) -> Vec3 {
        r.normalize() * Self::gradient_mag(r, h)
    }

    fn laplacian(_r: Vec3, _h: Scalar) -> Scalar {
        0.
    }
}

pub struct SpikyKernel;

impl SmoothingKernel for SpikyKernel {
    fn value(r: Vec3, h: Scalar) -> Scalar {
        let r_mag = r.magnitude();
        if r_mag >= 0. && r_mag <= h {
            let c = 15. / (std::f32::consts::PI * h.powi(6));
            let h_sub_r = h - r_mag;
            c * h_sub_r * h_sub_r * h_sub_r
        } else {
            0.
        }
    }

    fn gradient_mag(r: Vector3<Scalar>, h: Scalar) -> Scalar {
        let r_mag = r.magnitude();
        if r_mag >= 0. && r_mag <= h {
            let c = 15. * -3. / (std::f32::consts::PI * h.powi(6));
            let h_sub_r = h - r_mag;
            c * h_sub_r * h_sub_r
        } else {
            0.
        }
    }
}

pub struct Poly6Kernel;

impl SmoothingKernel for Poly6Kernel {
    fn value(r: Vector3<Scalar>, h: Scalar) -> Scalar {
        let c = 315. / (64. * std::f32::consts::PI * h.powi(9));
        let mag2 = r.magnitude2();
        // NOTE: I smell numerical instability coming from the >= 0, TODO, use an epsilon
        // value
        if mag2 <= h * h && mag2 >= 0. {
            c * (h * h - mag2).powi(3)
        } else {
            0.
        }
    }

    fn gradient_mag(r: Vec3, h: Scalar) -> Scalar {
        let c = 315. / (64. * std::f32::consts::PI * h.powi(9));
        let mag2 = r.magnitude2();
        if mag2 <= h * h && mag2 > 0. {
            c * 3. * -2. * mag2.sqrt() * (h * h - mag2) * (h * h - mag2)
        } else {
            0.
        }
    }
}

pub struct ViscosityKernel;

impl SmoothingKernel for ViscosityKernel {
    fn laplacian(r: Vec3, h: Scalar) -> Scalar {
        let c = 45. / (2. * std::f32::consts::PI * h.powi(6));

        let mag2 = r.magnitude2();
        if mag2 <= h * h && mag2 > 0. {
            c * (h - mag2.sqrt())
        } else {
            0.
        }
    }
}
