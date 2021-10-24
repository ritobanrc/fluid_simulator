type Float = f64;

type CGMatrix = na::Matrix4<Float>;
type CGVector = na::Vector4<Float>;

trait KrylovSystem {
    type KrylovVector;

    fn apply(&self, x: Self::KrylovVector) -> Self::KrylovVector;
}

#[allow(non_snake_case)]
fn conjugate_gradient(A: CGMatrix, b: CGVector, tol: Float) -> CGVector {
    let mut x = CGVector::zeros();
    let mut r = b - A * x;
    let mut p = r;

    loop {
        let s = A * p;
        let alpha = r.dot(&r) / (p.dot(&s));
        x = x + alpha * p;
        let r_prev = r;
        r = r - alpha * s;

        if r.abs().max() < tol {
            return x;
        }

        let beta = r.dot(&r) / r_prev.dot(&r_prev);
        p = r + beta * p;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_conjugate_gradient() {
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
