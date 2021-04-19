mod grid;
mod kernels;
mod render;

use crate::grid::{Coord, Grid};
use crate::kernels::{Poly6Kernel, SmoothingKernel, SpikyKernel, ViscosityKernel};
use crate::render::Vertex;

use std::sync::mpsc::channel;

use cgmath::prelude::*;
use cgmath::{point3, vec3, Point3, Vector3};
use num::{Float, FromPrimitive, Zero};
use rand::Rng;
use structopt::StructOpt;

type Scalar = f32;
type Vec3 = Vector3<Scalar>;

pub struct Simulation {
    pub masses: Vec<Scalar>,
    pub positions: Vec<Point3<Scalar>>,
    pub velocities: Vec<Vec3>,
    pub force: Vec<Vec3>,
    //pub grid: Grid,
    pub h: Scalar,
}

impl Simulation {
    pub fn new(h: Scalar) -> Self {
        Simulation {
            masses: Vec::new(),
            positions: Vec::new(),
            velocities: Vec::new(),
            force: Vec::new(),
            //grid: Grid::new(
            //(bounds.x / h) as usize,
            //(bounds.y / h) as usize,
            //(bounds.z / h) as usize,
            //),
            h,
        }
    }

    fn position_to_coord(&self, pos: Vec3) -> Coord {
        pos.map(|i| (i / self.h) as usize)
    }

    fn coord(&self, index: usize) -> Coord {
        self.position_to_coord(self.positions[index].to_vec())
    }

    fn add_particle(&mut self, position: Point3<Scalar>) {
        //let index = self.masses.len();
        self.masses.push(0.125);
        self.positions.push(position);
        self.velocities.push(Zero::zero());
        self.force.push(Zero::zero());

        //self.grid
        //.add_particle(self.position_to_coord(position.to_vec()), index);
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
    //let bounds = vec3(0.6, 1.1, 2.1);
    let h = 0.1;

    let mut s = Simulation::new(h);
    let mut rng = rand::thread_rng();

    let opt = Opt::from_args();

    // Water column scenario
    for x in linspace(0., 0.5, (0.5 / h) as usize) {
        for y in linspace(0., 1., (1. / h) as usize) {
            for z in linspace(0., 0.5, (0.5 / h) as usize) {
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
            //let neighbors = s.grid.get_neighbors(s.coord(i));

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
            //let old_coord = s.coord(i);
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

            //let new_coord = s.coord(i);
            //if new_coord != old_coord {
            //s.grid.update_particle(i, old_coord, new_coord);
            //}

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
            render::render_texture(path, rx, 1280, 720, opt.frames)
                .expect("Failed to recieve verticies");
        } else {
            eprintln!("Fluid sim is not being displayed or saved anywhere! Did you mean to run with -w or -i?")
        }
    }
}
