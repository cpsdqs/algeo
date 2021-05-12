use super::{BezierCurve, DerivativeSpace};
use cgmath::num_traits::NumCast;
use std::ops;

/// Returns the derivative of a bézier curve.
///
/// # Panics
/// - if the number of points cannot be represented by the scalar type S
///
/// # Examples
/// ```
/// # use cgmath::{Point2, Vector2};
/// # use cgmath::assert_abs_diff_eq;
/// # use cgmath::EuclideanSpace;
/// # use algeo::bezier::{self, DerivativeSpace};
/// # #[allow(non_upper_case_globals)]
/// # const to_vectors: fn(p: [Point2<f64>; 4]) -> [Vector2<f64>; 4] = DerivativeSpace::from_integral;
/// let control_points = [
///     Point2::new(0., 2.),
///     Point2::new(4., 3.),
///     Point2::new(6., 0.),
///     Point2::new(9., 4.),
/// ];
/// let derivative: [Vector2<f64>; 3] = bezier::derive(&control_points);
///
/// fn cubic_bezier_derivative(p: [Point2<f64>; 4], t: f64) -> Vector2<f64> {
///     let p: [Vector2<f64>; 4] = to_vectors(p);
///     3. * (1. - t) * (1. - t) * (p[1] - p[0])
///         + 6. * (1. - t) * t * (p[2] - p[1])
///         + 3. * t * t * (p[3] - p[2])
/// }
///
/// assert_abs_diff_eq!(
///     bezier::evaluate(&derivative, 0.),
///     cubic_bezier_derivative(control_points, 0.)
/// );
/// ```
pub fn derive<S, P, V, L, R>(points: &L) -> R
where
    L: BezierCurve<P>,
    P: ops::Sub<P, Output = V> + Clone,
    V: ops::Mul<S, Output = V>,
    S: NumCast + Clone,
    R: DerivativeSpace<L::Derivative> + BezierCurve<V>,
{
    let mut derivative = R::from_integral(points.reduced());
    let n: S =
        NumCast::from(points.count() - 1).expect("could not cast point count to scalar type");
    for i in 0..(points.count() - 1) {
        derivative.set(
            i,
            (points.get(i + 1).clone() - points.get(i).clone()) * n.clone(),
        );
    }
    derivative
}

#[test]
fn test_derive() {
    use super::evaluate;
    use cgmath::assert_abs_diff_eq;
    use cgmath::{Point2, Vector2};

    fn quad_bezier_derivative(
        p0: Point2<f64>,
        p1: Point2<f64>,
        p2: Point2<f64>,
        t: f64,
    ) -> Vector2<f64> {
        2. * (1. - t) * (p1 - p0) + 2. * t * (p2 - p1)
    }

    let a = Point2::new(0., 0.);
    let b = Point2::new(0., 1.);
    let c = Point2::new(1., 1.);
    let curve = [a, b, c];
    let derivative: [Vector2<f64>; 2] = derive(&curve);

    for t in [0., 0.3, 0.5, 0.9, 1.].iter() {
        assert_abs_diff_eq!(
            evaluate(&derivative, *t),
            quad_bezier_derivative(a, b, c, *t)
        );
    }
}
