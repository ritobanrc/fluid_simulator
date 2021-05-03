mod render;
mod sph;

use crate::render::Vertex;
use crate::sph::kernels::*;

use std::sync::mpsc::channel;

use cgmath::prelude::*;
use cgmath::{point3, vec3, Point3, Vector3};
use num::{Float, FromPrimitive, Zero};
use rand::Rng;
use structopt::StructOpt;

type Scalar = f32;
type Vec3 = Vector3<Scalar>;

fn linspace<T>(start: T, stop: T, num_steps: usize) -> impl Iterator<Item = T>
where
    T: Float + FromPrimitive,
{
    let delta: T = (stop - start) / T::from_usize(num_steps - 1).expect("out of range");
    return (0..num_steps).map(move |i| start + T::from_usize(i).expect("out of range") * delta);
}

fn update_bounds(
    position: &mut Point3<f32>,
    velocity: &mut Vec3,
    velocity_damping: Scalar,
    bounds_min: Vec3,
    bounds_max: Vec3,
) {
    (0..3).for_each(|i| {
        if position[i] < bounds_min[i] - 0.01 {
            velocity[i] *= -velocity_damping;
            position[i] = bounds_min[i];
        }

        if position[i] > bounds_max[i] + 0.01 {
            velocity[i] *= -velocity_damping;
            position[i] = bounds_max[i];
        }
    })
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
    let bounds = vec3(0.5, 1.5, 2.);
    let grid_bounds = bounds + vec3(0.1, 0.1, 0.1);

    let h = 0.04;

    let mut s = sph::Simulation::new(h, grid_bounds);
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
    println!("Created Grid with {:?} Cells", s.grid.grid.len());

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
            let neighbors = s.grid.get_neighbors(s.coord(i));
            densities.push(
                neighbors
                    .map(|j| s.masses[j] * Poly6Kernel::value(s.positions[i] - s.positions[j], h))
                    .sum(),
            );
        }

        for i in 0..num_particles {
            let pressure_i: Scalar = k * (densities[i] - rest_density);
            let neighbors = s.grid.get_neighbors(s.coord(i));

            let force_pressure = -neighbors
                .clone()
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
                * neighbors
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
            let old_coord = s.coord(i);
            s.positions[i] = s.positions[i] + delta_time * s.velocities[i];

            update_bounds(
                &mut s.positions[i],
                &mut s.velocities[i],
                0.8,
                Vector3::zero(),
                bounds,
            );

            let new_coord = s.coord(i);
            if new_coord != old_coord {
                s.grid.update_particle(i, old_coord, new_coord);
            }

            let pos = s.positions[i];
            let vel = s.velocities[i].magnitude2();
            //let color = s.velocities[i].magnitude2();
            verts.push(Vertex {
                position: [pos.x, pos.y, pos.z],
                color: [vel, 0.5 * vel + 0.5, 1.],
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
