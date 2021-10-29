mod math {
    pub const DIM: usize = 3;

    pub type Dim = na::Const<DIM>;

    pub type T = f64;
    pub type TV = na::SVector<T, DIM>;
    pub type IV = na::SVector<isize, DIM>;
    pub type UV = na::SVector<usize, DIM>;

    pub type Mat = na::SMatrix<T, DIM, DIM>;
}

use crate::mpm::models::IsotropicParameters;
use math::*;

type SimplexIdx = na::SVector<usize, { DIM + 1 }>;
type Simplex = na::SVector<TV, { DIM + 1 }>;

#[derive(Default)]
struct DeformedMesh {
    deformed_verts: Vec<TV>,
    base_verts: Vec<TV>,
    indices: Vec<SimplexIdx>,
}

impl DeformedMesh {
    fn num_verts(&self) -> usize {
        self.deformed_verts.len()
    }

    fn num_simplexes(&self) -> usize {
        self.indices.len()
    }

    fn get_base_simplex(&self, idx: usize) -> Simplex {
        self.indices[idx].map(|i| self.base_verts[i])
    }

    fn get_deformed_simplex(&self, idx: usize) -> Simplex {
        self.indices[idx].map(|i| self.deformed_verts[i])
    }
}

struct FiniteElements {
    mesh: DeformedMesh,
    mass: Vec<T>,
    velocity: Vec<TV>,
    force: Vec<TV>,
    params: FEMParameters,
}

struct FEMParameters {
    gravity: TV,
    dt: T,
    isotropic: IsotropicParameters,
}

impl Default for FEMParameters {
    fn default() -> Self {
        FEMParameters {
            gravity: TV::new(0., -9.8, 0.),
            dt: 0.01,
            isotropic: IsotropicParameters::default(),
        }
    }
}

trait BasicGeometry {
    type TV;
    fn normal(&self) -> Self::TV;
}

/// Triangle in 3D
impl BasicGeometry for na::Vector3<TV> {
    type TV = na::Vector3<T>;

    fn normal(&self) -> Self::TV {
        (self[1] - self[0]).cross(&(self[2] - self[0])).normalize()
    }
}

/// Line
impl BasicGeometry for na::Vector2<TV> {
    type TV = na::Vector2<T>;

    fn normal(&self) -> Self::TV {
        let a = self[1] - self[0];
        Self::TV::new(-a.y, a.x)
    }
}

impl crate::Simulation for FiniteElements {
    type Parameters = FEMParameters;

    fn new(params: Self::Parameters) -> Self {
        Self {
            mesh: DeformedMesh::default(),
            mass: Vec::new(),
            velocity: Vec::new(),
            force: Vec::new(),
            params,
        }
    }

    fn add_particle(&mut self, mass: crate::Scalar, position: crate::Vec3, velocity: crate::Vec3) {
        panic!("adding particles to FEM does not make sense.")
    }

    #[allow(non_snake_case)]
    fn simulate_frame(&mut self) -> Vec<crate::render::Vertex> {
        let verts = Vec::new();

        for v in 0..self.mesh.num_verts() {
            self.force[v] = self.mass[v] * self.params.gravity;
        }

        for simplex in 0..self.mesh.num_simplexes() {
            let base = self.mesh.get_base_simplex(simplex);
            let deformed = self.mesh.get_deformed_simplex(simplex);

            let mut base_mat = Mat::zeros();
            let mut deformed_mat = Mat::zeros();

            for a in 1..DIM {
                *base_mat.column_mut(a - 0) = *(base[a] - base[0]);
                *deformed_mat.column_mut(a - 1) = *(deformed[a] - deformed[0]);
            }

            let inv_base = base_mat
                .try_inverse()
                .expect(&format!("Element is not invertible: {:?}", base));

            let F = deformed_mat * inv_base;

            let (Q, F_tilde) = (|| {
                let svd = F.svd(true, true);
                let Q = svd.u? * svd.v_t?;
                let F_tilde =
                    svd.v_t?.transpose() * Mat::from_diagonal(&svd.singular_values) * svd.v_t?;

                Some((Q, F_tilde))
            })()
            .expect("failed to compute polar decomposition somehow");

            let strain = 0.5 * (F_tilde * F_tilde.transpose()) - Mat::identity();
            let lambda = self.params.isotropic.lambda;
            let mu = self.params.isotropic.mu;
            let stress = lambda * strain.trace() * Mat::identity() + 2. * mu * strain;

            for (i, &v) in self.mesh.indices[simplex].iter().enumerate() {
                let tri = base.remove_row(i);
                let normal = tri.normal();
                self.force[v] = Q * stress * normal;
            }
        }

        for v in 0..self.mesh.num_verts() {
            self.velocity[v] += self.params.dt * (self.force[v] / self.mass[v]);
            self.mesh.deformed_verts[v] += self.params.dt * self.velocity[v];
        }

        verts
    }
}
