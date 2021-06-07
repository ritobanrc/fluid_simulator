use nalgebra::Vector3;
use std::ops::RangeInclusive as Range;

type Scalar = f32;
type Vec3 = Vector3<Scalar>;

pub struct MpmParameters {
    /// The total number of Lagrangian particles in the simulation
    pub num_particles: usize,
    /// The radius of each particle
    pub h: Scalar,
    /// The bounds of the simulation. Note that this is a `RangeInclusive`. The upper bounds
    /// are included as part of the simulation domain
    pub bounds: Range<Vec3>,
    /// The size of the time step. Larger time steps will simulate faster, but may be unstable or
    /// innaccurate.
    pub delta_time: f32,
}

impl Default for MpmParameters {
    fn default() -> Self {
        MpmParameters {
            num_particles: 0,
            h: 0.05,
            bounds: Vector3::zeros()..=Vector3::new(5., 5., 5.),
            delta_time: 0.01,
        }
    }
}
