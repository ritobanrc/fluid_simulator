mod render;

use cgmath::prelude::*;
use cgmath::{point3, Point3, Vector3};
use rand::Rng;
use render::Vertex;

trait SPHDiscretization {}

type ParticleIdx = usize;
type Scalar = f32;
type Vec3 = Vector3<Scalar>;

pub struct Simulation {
    pub masses: Vec<Scalar>,
    pub positions: Vec<Point3<Scalar>>,
    pub velocities: Vec<Point3<Scalar>>,
    pub accelerations: Vec<Point3<Scalar>>,
}

trait SmoothingKernel {
    fn value(r: Vec3, h: Scalar) -> Scalar;

    fn gradient_mag(r: Vec3, h: Scalar) -> Scalar;
}

struct SpikyKernel;

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

struct Poly6Kernel;

impl SmoothingKernel for Poly6Kernel {
    fn value(r: Vector3<Scalar>, h: Scalar) -> Scalar {
        let c = 315. / (64. * std::f32::consts::PI * h.powi(9));
        let mag2 = r.magnitude2();
        if mag2 <= h * h && mag2 > 0. {
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

fn main() {
    let num_particles = 1000;
    let mut s = Simulation {
        masses: Vec::new(),
        positions: Vec::new(),
        velocities: Vec::new(),
        accelerations: Vec::new(),
    };

    let mut rng = rand::thread_rng();

    let h = 1.;
    for _i in 0..num_particles {
        s.masses.push(1.0);
        s.positions.push(point3(
            rng.gen::<Scalar>() * 2. - 1.,
            rng.gen::<Scalar>() * 2. - 1.,
            rng.gen::<Scalar>() * 2. - 1.,
        ));
        s.velocities.push(point3(0., 0., 0.));
        s.accelerations.push(point3(0., 0., 0.));
    }

    let mut verts = Vec::with_capacity(num_particles);

    for i in 0..num_particles {
        let density: Scalar = (0..num_particles)
            .map(|j| s.masses[j] * Poly6Kernel::value(s.positions[i] - s.positions[j], h))
            .sum();

        dbg!(density);

        let pos = s.positions[i];
        verts.push(Vertex {
            position: [pos.x, pos.y, pos.z],
            color: [density / 150., 1., density / 150.],
        });
    }

    render::open_window(&verts);
}
