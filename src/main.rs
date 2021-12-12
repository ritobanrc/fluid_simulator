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
    input_file: Option<std::path::PathBuf>,
    #[structopt(short, long)]
    output_dir: Option<std::path::PathBuf>,
    #[structopt(short, long, default_value = "600")]
    frames: usize,
}

fn main() -> eyre::Result<()> {
    let opt = Opt::from_args();

    use eyre::WrapErr;

    if let Some(input_file) = opt.input_file {
        if input_file.is_dir() {
            let mut files = std::fs::read_dir(input_file)?
                .filter_map(|entry| Some(entry.ok()?.path()))
                .collect::<Vec<_>>();

            files.sort();

            render::open_window(Some(files)).expect("Failed to recieve vertecies");
        } else if input_file.extension().unwrap() == "json" {
            if let Some(path) = opt.output_dir {
                let ui_state: render::UIState = std::fs::read(&input_file)
                    .wrap_err_with(|| {
                        format!("Failed to read JSON settings file: {:?}", &input_file)
                    })
                    .and_then(|json| {
                        serde_json::from_slice(&json).wrap_err("Serde failed to deserialize JSON.")
                    })?;

                let (stop_tx, stop_rx) = std::sync::mpsc::channel();
                let vert_rx = render::start_simulation(&ui_state, stop_rx);

                for frame in 0..opt.frames {
                    println!("starting Frame: {:?}", frame);
                    let verts = vert_rx.recv()?;
                    let mut path = path.clone();
                    path.push(format!("{:03}.dat", frame));
                    let mut writer = std::fs::File::create(&path)?;
                    rmp_serde::encode::write(&mut writer, &verts)?;
                }

                drop(stop_tx);
            } else {
                return Err(eyre::eyre!("Must specify output directory."));
            }
        }
    } else {
        println!("Displaying fluid simulation in window.");
        render::open_window(None).expect("Failed to recieve vertecies");
    }

    Ok(())
}
