use super::BezierCurve;
use std::ops;

/// Splits a bézier curve at t to create two new bézier curves.
pub fn subdivide<S, P, V, L>(points: &L, t: S) -> (L, L)
where
    L: BezierCurve<P>,
    P: ops::Sub<P, Output = V> + ops::Add<V, Output = P> + Clone,
    V: ops::Mul<S, Output = V>,
    S: Clone,
{
    let mut points_a = points.clone();
    let mut points_b = points.clone();

    // the control points for the two new curves are simply the intermediate starting points
    // in the de Casteljau algorithm, which is very convenient!
    #[inline]
    fn de_casteljau<S, P, V, L, L2>(points: &L, t: S, points_a: &mut L2, points_b: &mut L2)
    where
        L: BezierCurve<P>,
        L2: BezierCurve<P>,
        P: ops::Sub<P, Output = V> + ops::Add<V, Output = P> + Clone,
        V: ops::Mul<S, Output = V>,
        S: Clone,
    {
        let index = points_a.count() - points.count();
        points_a.set(index, points.get(0).clone());
        points_b.set(
            points_b.count() - index - 1,
            points.get(points.count() - 1).clone(),
        );

        if points.count() > 1 {
            let mut new_points = points.reduced();
            for i in 0..new_points.count() {
                let p = points.get(i).clone()
                    + (points.get(i + 1).clone() - points.get(i).clone()) * t.clone();
                new_points.set(i, p);
            }
            de_casteljau(&new_points, t, points_a, points_b);
        }
    }

    de_casteljau(points, t, &mut points_a, &mut points_b);

    (points_a, points_b)
}

#[test]
fn test_subdiv() {
    use super::evaluate;
    use cgmath::assert_abs_diff_eq;
    use cgmath::Vector2;

    let a = Vector2::new(0., 1.);
    let b = Vector2::new(5., 3.);
    let c = Vector2::new(3., 8.);
    let d = Vector2::new(8., 2.);
    let curve = [a, b, c, d];
    let (split_a, split_b) = subdivide(&curve, 0.5);

    assert_abs_diff_eq!(evaluate(&curve, 0.), evaluate(&split_a, 0.),);
    assert_abs_diff_eq!(evaluate(&curve, 0.25), evaluate(&split_a, 0.5));
    assert_abs_diff_eq!(evaluate(&curve, 0.5), evaluate(&split_a, 1.));
    assert_abs_diff_eq!(evaluate(&curve, 0.5), evaluate(&split_b, 0.));
    assert_abs_diff_eq!(evaluate(&curve, 0.75), evaluate(&split_b, 0.5));
    assert_abs_diff_eq!(evaluate(&curve, 1.), evaluate(&split_b, 1.));
}
