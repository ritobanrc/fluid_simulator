use crate::math::*;
use nalgebra_sparse::csr::CsrMatrix;

#[allow(non_snake_case)]
fn conjugate_gradient(A: &CsrMatrix<T>, b: &na::DVector<T>, tol: T) -> na::DVector<T> {
    let mut x = na::DVector::zeros(b.nrows());
    let m: na::DVector<_> = A * &x;
    let mut r = b - m;
    let mut p = r.clone();

    loop {
        let s = A * &p;
        let alpha = r.dot(&r) / (p.dot(&s)); // step length
        x = x + alpha * &p; // approximate solution
        let r_next = &r - alpha * s; // calculate residual

        if &r.abs().max() < &tol {
            return x;
        }

        let beta = r_next.dot(&r_next) / r.dot(&r); // improvement this step
        p = &r_next + beta * &p; // new search direction

        r = r_next;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_conjugate_gradient() {
        type CGMatrix = na::Matrix4<T>;
        type CGVector = na::Vector4<T>;

        fn solve(A: CGMatrix, b: CGVector) -> CGVector {
            A.try_inverse().unwrap() * b
        }

        fn random_spd() -> CGMatrix {
            let a = CGMatrix::new_random();
            a * a.transpose()
        }

        for _ in 1..1000 {
            let a = random_spd();
            let b = CGVector::new_random();

            let exact = solve(a, b);
            let approx = conjugate_gradient(a, b, 10e-15);

            assert!(
                ((approx - exact).abs().component_div(&exact)).max() < 10e-5,
                "A = {}, b = {}, exact = {}, approx = {}",
                a,
                b,
                exact,
                approx
            );
        }
    }
}
