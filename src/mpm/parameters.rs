use nalgebra::Vector3;
use std::ops::Range;

type Scalar = f32;
type Vec3 = Vector3<Scalar>;

pub struct MpmParameters {
    /// The total number of Lagrangian particles in the simulation
    pub num_particles: usize,
    /// The radius of each particle/the grid spacing
    pub h: Scalar,
    /// The bounds of the simulation.
    pub bounds: Range<Vec3>,
    /// The size of the time step. Larger time steps will simulate faster, but may be unstable or
    /// innaccurate.
    pub delta_time: f32,
    /// Whether or not to use APIC (as opposed to PIC). TODO: Handle both FLIP and PIC
    pub use_affine: bool,
}

impl Default for MpmParameters {
    fn default() -> Self {
        MpmParameters {
            num_particles: 0,
            h: 0.05,
            bounds: Vector3::zeros()..Vector3::new(2., 2., 2.),
            delta_time: 0.01,
            use_affine: true,
        }
    }
}
