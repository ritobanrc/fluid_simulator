use crate::{Scalar, Simulation, Vec3};
use itertools::iproduct;
use na::Vector3;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::ops::Range;

pub trait InitialCondition {
    fn add_particles<S: Simulation>(&self, s: &mut S);
}

pub struct Block {
    pub total_mass: Scalar,
    pub size: Range<Vec3>,
    pub spacing: Scalar,
    pub jitter: Vec3,
}

impl Default for Block {
    fn default() -> Self {
        Block {
            total_mass: 500.,
            size: Vec3::from_element(0.5)..Vec3::from_element(1.5),
            spacing: 0.05,
            jitter: Vec3::from_element(0.05 / 8.),
        }
    }
}

impl InitialCondition for Block {
    fn add_particles<S: Simulation>(&self, s: &mut S) {
        let mut rng = StdRng::from_seed([0; 32]);

        let min = self.size.start;
        let max = self.size.end;

        if self.spacing == 0.0 {
            return;
        }

        let counts = ((max - min) / self.spacing).map(|x| x.ceil() as usize);
        let num_particles: usize = counts.iter().product();
        let mass = self.total_mass / num_particles as Scalar;

        for (i, j, k) in iproduct!(0..counts.x, 0..counts.y, 0..counts.z) {
            let idx = Vector3::new(i, j, k);
            let pos = idx.cast::<Scalar>() * self.spacing + min;

            let rand: Vec3 = rng.gen::<[Scalar; 3]>().into();
            let jitter = rand.component_mul(&self.jitter) - self.jitter / 2.;

            let pos = pos + jitter;

            s.add_particle(mass, pos, Vector3::zeros());
        }
    }
}

pub struct Sphere {
    pub num_particles: usize,
    pub total_mass: Scalar,
    pub center: Vec3,
    pub radius: f64,
}

impl Default for Sphere {
    fn default() -> Self {
        Sphere {
            num_particles: 5000,
            total_mass: 500.,
            center: Vec3::new(1., 1., 1.),
            radius: 0.25,
        }
    }
}

impl InitialCondition for Sphere {
    fn add_particles<S: Simulation>(&self, s: &mut S) {
        let mut rng = StdRng::from_seed([0; 32]);
        let mass = self.total_mass / self.num_particles as Scalar;

        for _ in 0..self.num_particles {
            let pos = loop {
                let rand: Vec3 = rng.gen::<[f64; 3]>().into();
                let pos = rand * 2. - Vec3::from_element(1.);

                if pos.magnitude_squared() < self.radius * self.radius {
                    break pos + self.center;
                }
            };

            s.add_particle(mass, pos, Vec3::zeros());
        }
    }
}
