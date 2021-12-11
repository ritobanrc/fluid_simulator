mod eulerian;
mod initial_condition;
mod mpm;
mod render;
mod sph;
mod statistics;
mod util;

extern crate nalgebra as na;

use crate::mpm::{MpmParameters, MpmSimulation};
use crate::render::Vertex;
use crate::sph::SphSimulation;

use std::sync::mpsc::channel;

use structopt::StructOpt;

type Scalar = f64;
type Vec3 = na::Vector3<Scalar>;
pub trait Simulation: Send {
    type Parameters;

    fn new(params: Self::Parameters) -> Self;

    fn simulate_frame(&mut self) -> Vec<Vertex>;

    fn add_particle(&mut self, mass: Scalar, position: Vec3, velocity: Vec3);
}

pub mod math {
    pub const DIM: usize = 3;

    pub type Dim = na::Const<DIM>;

    pub type T = f64;
    pub type TV = na::SVector<T, DIM>;
    pub type IV = na::SVector<isize, DIM>;
    pub type UV = na::SVector<usize, DIM>;

    pub type Mat = na::SMatrix<T, DIM, DIM>;
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

fn main() -> eyre::Result<()> {
    let opt = Opt::from_args();

    use crate::initial_condition::InitialCondition;
    if opt.window {
        println!("Displaying fluid simulation in window.");
        render::open_window().expect("Failed to recieve vertecies");
    } else if let Some(path) = opt.output_dir {
        let (tx, rx) = channel::<Vec<Vertex>>();

        let params = MpmParameters::<mpm::NeoHookean>::default();
        let mut s = MpmSimulation::new(params);
        crate::initial_condition::Block::default().add_particles(&mut s);

        println!(
            "Running simulation with {:?} particles",
            s.params.num_particles
        );

        std::thread::spawn(move || loop {
            let verts = s.simulate_frame();

            tx.send(verts).unwrap();
        });

        std::fs::create_dir_all(&path).unwrap();
        render::render_texture(path, rx, 1920, 1080, opt.frames)
            .expect("Failed to recieve verticies");
    } else {
        eprintln!("Fluid sim is not being displayed or saved anywhere! Did you mean to run with -w or -i?")
    }

    Ok(())
}
