use crate::{Simulation, Vec3};
use num::{Float, FromPrimitive};
use rand::{rngs::StdRng, Rng, SeedableRng};

pub trait InitialCondition {
    fn add_particles<S: Simulation>(&self, s: &mut S);
}

pub struct Block;

fn linspace<T>(start: T, stop: T, num_steps: usize) -> impl Iterator<Item = T>
where
    T: Float + FromPrimitive,
{
    let delta: T = (stop - start) / T::from_usize(num_steps - 1).expect("out of range");
    (0..num_steps).map(move |i| start + T::from_usize(i).expect("out of range") * delta)
}

impl InitialCondition for Block {
    fn add_particles<S: Simulation>(&self, s: &mut S) {
        let mut rng = StdRng::from_seed([0; 32]);
        let h = 0.05; // TODO: Get this form simulation or smth

        // Block scenario
        // TODO: Don't hardcode all these numbers, make them possible to set from the gui
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
                    s.add_particle([x + jitter_x, y, z + jitter_z].into(), velocity.into());
                }
            }
        }
    }
}

pub struct Sphere {
    pub num_particles: usize,
    pub center: Vec3,
    pub radius: f64,
}

impl Default for Sphere {
    fn default() -> Self {
        Sphere {
            num_particles: 5000,
            center: Vec3::new(1., 1., 1.),
            radius: 0.25,
        }
    }
}

impl InitialCondition for Sphere {
    fn add_particles<S: Simulation>(&self, s: &mut S) {
        let mut rng = StdRng::from_seed([0; 32]);

        for _ in 0..self.num_particles {
            let pos = loop {
                let rand: Vec3 = rng.gen::<[f64; 3]>().into();
                let pos = rand * 2. - Vec3::from_element(1.);

                if pos.magnitude_squared() < self.radius * self.radius {
                    break pos + self.center;
                }
            };

            s.add_particle(pos, Vec3::zeros());
        }
    }
}
