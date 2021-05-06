use crate::{Scalar, Vec3};
use cgmath::Matrix3;

struct MpmSimulation {
    masses: Vec<Scalar>,
    volume: Vec<Scalar>,
    position: Vec<Vec3>,
    velocities: Vec<Vec3>,
    affine_matrices: Vec<Matrix3<Scalar>>,
}
