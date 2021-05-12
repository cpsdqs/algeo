use super::BezierCurve;
use cgmath::num_traits::Float;
use cgmath::MetricSpace;
use std::ops;

/// Returns cheap lower and upper bounds for the arc length of a bézier curve.
///
/// Specifically, this will be the length from the start point to the end point as the lower bound,
/// and the polyline through all control points as the upper bound.
pub fn hull_arclen_bounds<S, P, L>(points: &L) -> (S, S)
where
    L: BezierCurve<P>,
    P: MetricSpace<Metric = S> + Clone,
    S: Float + ops::AddAssign<S>,
{
    let lower = points
        .get(0)
        .clone()
        .distance(points.get(points.count() - 1).clone());
    let mut upper = S::zero();
    for i in 0..(points.count() - 1) {
        upper += points.get(i).clone().distance(points.get(i + 1).clone());
    }
    (lower, upper)
}
