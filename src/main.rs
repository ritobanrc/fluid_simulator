mod kernels;
mod render;

use crate::kernels::{Poly6Kernel, SmoothingKernel, SpikyKernel, ViscosityKernel};

use std::sync::mpsc::channel;

use cgmath::prelude::*;
use cgmath::{point3, vec3, Point3, Vector3};
use num::{Float, FromPrimitive, Zero};
use rand::Rng;
use render::Vertex;
use structopt::StructOpt;

type Scalar = f32;
type Vec3 = Vector3<Scalar>;

#[derive(Default)]
pub struct Simulation {
    pub masses: Vec<Scalar>,
    pub positions: Vec<Point3<Scalar>>,
    pub velocities: Vec<Vec3>,
    pub force: Vec<Vec3>,
}

impl Simulation {
    fn add_particle(&mut self, position: Point3<Scalar>) {
        self.masses.push(1.0);
        self.positions.push(position);
        self.velocities.push(Zero::zero());
        self.force.push(Zero::zero());
    }
}

fn linspace<T>(start: T, stop: T, num_steps: usize) -> impl Iterator<Item = T>
where
    T: Float + FromPrimitive,
{
    let delta: T = (stop - start) / T::from_usize(num_steps - 1).expect("out of range");
    return (0..num_steps).map(move |i| start + T::from_usize(i).expect("out of range") * delta);
}

#[derive(StructOpt, Debug)]
#[structopt(name = "sph_solver")]
struct Opt {
    #[structopt(short, long)]
    window: bool,
    #[structopt(short, long)]
    image_dir: Option<std::path::PathBuf>,
    #[structopt(short, long, default_value = "600")]
    frames: usize,
}

fn main() {
    //let num_particles = 1000;
    let mut s = Simulation::default();
    let mut rng = rand::thread_rng();

    let opt = Opt::from_args();

    let h = 0.1;
    // Water column scenario
    for x in linspace(0., 0.5, 5) {
        for y in linspace(0., 1., 10) {
            for z in linspace(0., 0.5, 5) {
                let jitter_x = rng.gen::<f32>() * h / 8. - h / 16.;
                let jitter_z = rng.gen::<f32>() * h / 8. - h / 16.;
                s.add_particle(point3(x + jitter_x, y, z + jitter_z));
            }
        }
    }

    let num_particles = s.masses.len();
    println!("Running simulation with {:?} particles", num_particles);

    let delta_time = 0.01;
    let rest_density = 1000.;
    let k = 4.;
    let mu = 8.;

    let gravity = -Vec3::unit_y();
    let (tx, rx) = channel::<Vec<Vertex>>();

    std::thread::spawn(move || loop {
        let mut verts = Vec::with_capacity(num_particles);

        let mut densities: Vec<Scalar> = Vec::with_capacity(num_particles);
        for i in 0..num_particles {
            densities.push(
                (0..num_particles)
                    .map(|j| s.masses[j] * Poly6Kernel::value(s.positions[i] - s.positions[j], h))
                    .sum(),
            );
        }

        for i in 0..num_particles {
            let pressure_i: Scalar = k * (densities[i] - rest_density);
            //dbg!(pressure_i);
            let force_pressure = -(0..num_particles)
                .map(|j| {
                    if i == j {
                        return Vec3::zero();
                    }
                    let r_ij = s.positions[i] - s.positions[j];

                    if r_ij.magnitude2() > h * h {
                        return Vec3::zero();
                    }

                    let pressure_j = k * (densities[j] - rest_density);
                    //dbg!(pressure_j);
                    s.masses[j] * (pressure_i + pressure_j) / (2. * densities[j])
                        * SpikyKernel::gradient(r_ij, h)
                })
                .sum::<Vec3>();

            let force_viscosity = mu
                * (0..num_particles)
                    .map(|j| {
                        if i == j {
                            return Vec3::zero();
                        }
                        let vdiff = s.velocities[j] - s.velocities[i];

                        let r_ij = s.positions[i] - s.positions[j];

                        if r_ij.magnitude2() > h * h {
                            return Vec3::zero();
                        }

                        s.masses[j] * vdiff / densities[j] * ViscosityKernel::laplacian(r_ij, h)
                    })
                    .sum::<Vec3>();

            let force_gravity = gravity * densities[i];
            s.force[i] = force_pressure + force_viscosity + force_gravity;
        }

        for i in 0..num_particles {
            s.velocities[i] = s.velocities[i] + delta_time / densities[i] * s.force[i];
            s.positions[i] = s.positions[i] + delta_time * s.velocities[i];

            if s.positions[i].y < -0.01 {
                s.velocities[i].y *= -0.5;
                s.positions[i].y = 0.;
            }

            if s.positions[i].x < -0.03 {
                s.velocities[i].x *= -0.5;
                s.positions[i].x = 0.;
            }

            if s.positions[i].x > 0.53 {
                s.velocities[i].x *= -0.5;
                s.positions[i].x = 0.5;
            }

            if s.positions[i].z < -0.03 {
                s.velocities[i].z *= -0.5;
                s.positions[i].z = 0.;
            }

            if s.positions[i].z > 2. {
                s.velocities[i].z *= -0.5;
                s.positions[i].z = 1.99;
            }

            let pos = s.positions[i];
            let density = densities[i];
            verts.push(Vertex {
                position: [pos.x, pos.y, pos.z],
                color: [density / 150., 1., density / 150.],
            });
        }

        tx.send(verts).unwrap();
    });

    if opt.window {
        println!("Displaying fluid simulation in window.");
        render::open_window(rx).expect("Failed to recieve vertecies");
    } else {
        if let Some(path) = opt.image_dir {
            std::fs::create_dir_all(&path).unwrap();
            render::render_texture(path, rx, 512, 256, opt.frames)
                .expect("Failed to recieve verticies");
        } else {
            eprintln!("Fluid sim is not being displayed or saved anywhere! Did you mean to run with -w or -i?")
        }
    }
}
