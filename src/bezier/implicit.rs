use crate::polynomial::{Poly1x2d, Poly3x2d};
use cgmath::{BaseFloat, Point2};

/// Returns the determinant of the following matrix:
///
/// ```text
/// | x  y  1 |
/// | xi yi 1 |
/// | xj yj 1 |
/// ```
fn impl_l_det<S>(pi: Point2<S>, pj: Point2<S>) -> Poly1x2d<S>
where
    S: BaseFloat,
{
    let a = Poly1x2d {
        k: S::zero(),
        x: pi.y - pj.y,
        y: S::zero(),
    };
    let b = Poly1x2d {
        k: S::zero(),
        x: S::zero(),
        y: pi.x - pj.x,
    };
    let c = Poly1x2d {
        k: pi.x * pj.y - pj.x * pi.y,
        x: S::zero(),
        y: S::zero(),
    };
    a - b + c
}

fn binom(n: usize, k: usize) -> f64 {
    (1..=k).map(|i| (n + 1 - i) as f64 / i as f64).product()
}

fn impl_l<S, const N: usize>(curve: [Point2<S>; N], i: usize, j: usize) -> Poly1x2d<S>
where
    S: BaseFloat,
{
    let det = impl_l_det(curve[i], curve[j]);
    let n = N - 1; // degree
    det * S::from(binom(n, i) * binom(n, j)).unwrap()
}

/// Expands the determinant of a row major matrix of polynomials.
fn expand_det3<S>(matrix: [[Poly1x2d<S>; 3]; 3]) -> Poly3x2d<S>
where
    S: BaseFloat,
{
    let a = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2]);
    let b = matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[2][0] * matrix[1][2]);
    let c = matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]);
    a - b + c
}

/// Returns an implicit function for a 2D cubic bézier curve.
///
/// The curve is located at f(x, y) = 0.
///
/// # Details
/// Curve implicitization is implemented using the method outlined in chapter 17 of "Computer
/// Aided Geometric Design" by Thomas W. Sederberg.
pub fn implicit_cubic<S>(curve: [Point2<S>; 4]) -> Poly3x2d<S>
where
    S: BaseFloat,
{
    let l32 = impl_l(curve, 3, 2);
    let l31 = impl_l(curve, 3, 1);
    let l30 = impl_l(curve, 3, 0);
    let l21 = impl_l(curve, 2, 1);
    let l20 = impl_l(curve, 2, 0);
    let l10 = impl_l(curve, 1, 0);

    expand_det3([[l32, l31, l30], [l31, l30 + l21, l20], [l30, l20, l10]])
}

#[test]
fn test_implicit_cubic() {
    use super::evaluate;
    use cgmath::{assert_relative_eq, assert_relative_ne};

    let curve1 = [
        Point2::new(1_f64, 0.),
        Point2::new(5., 0.),
        Point2::new(5., 2.),
        Point2::new(4., 3.),
    ];
    let i_curve1 = implicit_cubic(curve1);

    for i in 0..10 {
        let t = (i as f64) / 10.;

        let curve_point = evaluate(&curve1, t);

        // the curve point is on the curve! (implicit equation = 0)
        assert_relative_eq!(
            i_curve1.eval(curve_point.x, curve_point.y),
            0.,
            epsilon = 1e-8
        );

        // try some random points that shouldn't be on the curve
        assert_relative_ne!(
            i_curve1.eval(curve_point.x, curve_point.y + 1.),
            0.,
            epsilon = 1e-8
        );
        assert_relative_ne!(
            i_curve1.eval(curve_point.x + 0.1, curve_point.y),
            0.,
            epsilon = 1e-8
        );
    }
}
