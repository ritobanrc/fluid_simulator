mod mpm;
mod render;
mod sph;

extern crate nalgebra as na;

use crate::mpm::{MpmParameters, MpmSimulation};
use crate::render::Vertex;
use crate::sph::{SphParamaters, SphSimulation};

use std::sync::mpsc::channel;

use cgmath::{point3, vec3, Point3, Vector3};

use num::{Float, FromPrimitive};
use rand::{rngs::StdRng, Rng, SeedableRng};
use structopt::StructOpt;

type Scalar = f32;
type Vec3 = cgmath::Vector3<Scalar>;

fn linspace<T>(start: T, stop: T, num_steps: usize) -> impl Iterator<Item = T>
where
    T: Float + FromPrimitive,
{
    let delta: T = (stop - start) / T::from_usize(num_steps - 1).expect("out of range");
    return (0..num_steps).map(move |i| start + T::from_usize(i).expect("out of range") * delta);
}

pub(crate) fn update_bounds(
    position: &mut Vec3,
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
    output_dir: Option<std::path::PathBuf>,
    #[structopt(short, long, default_value = "600")]
    frames: usize,
}

trait Simulation {
    fn simulate_frame(&mut self) -> Vec<Vertex>;
}

fn main() {
    let params = MpmParameters::default();
    let h = params.h;

    let mut s = MpmSimulation::new(params);
    let mut rng = StdRng::from_seed([0; 32]);

    let opt = Opt::from_args();

    // Block scenario
    for x in linspace(0.5, 1.5, (1. / h) as usize) {
        for y in linspace(0.5, 1.5, (1. / h) as usize) {
            for z in linspace(0.5, 1.5, (1. / h) as usize) {
                let jitter_x = rng.gen::<f64>() * h / 8. - h / 16.;
                let jitter_z = rng.gen::<f64>() * h / 8. - h / 16.;

                let velocity = [0., 0., 0.];
                //rng.gen::<f64>() * 2. - 1.,
                //rng.gen::<f64>() * 2. - 1.,
                //rng.gen::<f64>() * 2. - 1.,
                //];
                s.add_particle([x + jitter_x, y, z + jitter_z], velocity);
            }
        }
    }

    println!(
        "Running simulation with {:?} particles",
        s.params.num_particles
    );

    let (tx, rx) = channel::<Vec<Vertex>>();

    std::thread::spawn(move || loop {
        let verts = s.simulate_frame();

        tx.send(verts).unwrap();
    });

    if opt.window {
        println!("Displaying fluid simulation in window.");
        render::open_window(rx).expect("Failed to recieve vertecies");
    } else if let Some(path) = opt.output_dir {
        std::fs::create_dir_all(&path).unwrap();
        render::render_texture(path, rx, 1920, 1080, opt.frames)
            .expect("Failed to recieve verticies");
    } else {
        eprintln!("Fluid sim is not being displayed or saved anywhere! Did you mean to run with -w or -i?")
    }
}
