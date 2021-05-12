use super::BezierCurve;
use cgmath::BaseFloat;
use std::ops;

/// Implements De Casteljau’s algorithm to evaluate a bézier curve at the given location.
///
/// Will probably panic if the list of points is empty.
pub fn evaluate<S, P, V, L>(points: &L, t: S) -> P
where
    L: BezierCurve<P>,
    P: ops::Sub<P, Output = V> + ops::Add<V, Output = P> + Clone,
    V: ops::Mul<S, Output = V>,
    S: BaseFloat + Clone,
{
    if points.count() <= 1 {
        points.get(0).clone()
    } else {
        let mut new_points = points.reduced();
        for i in 0..new_points.count() {
            let p = points.get(i).clone() + (points.get(i + 1).clone() - points.get(i).clone()) * t;
            new_points.set(i, p);
        }
        evaluate(&new_points, t)
    }
}

#[test]
fn test_evaluate() {
    use cgmath::assert_abs_diff_eq;
    use cgmath::Point2;

    let a = Point2::new(1., 0.);
    let b = Point2::new(5., 1.);

    let ab_mid = evaluate(&[a, b], 0.5);
    assert_abs_diff_eq!(ab_mid, Point2::new(3., 0.5));

    let a = Point2::new(0., 0.);
    let b = Point2::new(0., 1.);
    let c = Point2::new(1., 1.);

    let abc_start = evaluate(&[a, b, c], 0.);
    let abc_mid = evaluate(&[a, b, c], 0.5);
    let abc_end = evaluate(&[a, b, c], 1.);
    assert_abs_diff_eq!(abc_start, a);
    assert_abs_diff_eq!(abc_mid, Point2::new(0.25, 0.75));
    assert_abs_diff_eq!(abc_end, c);
}
