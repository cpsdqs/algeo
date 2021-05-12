use crate::polynomial::Poly3;
use cgmath::{BaseFloat, Point2};

/// Returns the parametric form of a 2D cubic bézier curve.
///
/// Returns two polynomials for (x, y).
pub fn parametric_cubic<S>(curve: [Point2<S>; 4]) -> (Poly3<S>, Poly3<S>)
where
    S: BaseFloat,
{
    // (1-t)^3 A + 3t(1-t)^2 B + 3t^2(1-t) C + t^3 D
    // = (-t^3 + 3t^2 - 3t + 1) A
    //   + (3t^3 - 6t^2 + 3t) B
    //   + (-3t^3 + 3t^2) C
    //   + t^3 D

    let n3 = S::from(3).unwrap();
    let n6 = S::from(6).unwrap();

    let [a, b, c, d] = curve;
    let x_ttt = -a.x + n3 * b.x - n3 * c.x + d.x;
    let y_ttt = -a.y + n3 * b.y - n3 * c.y + d.y;
    let x_tt = n3 * a.x - n6 * b.x + n3 * c.x;
    let y_tt = n3 * a.y - n6 * b.y + n3 * c.y;
    let x_t = -n3 * a.x + n3 * b.x;
    let y_t = -n3 * a.y + n3 * b.y;
    let x_k = a.x;
    let y_k = a.y;

    (
        Poly3 {
            k: x_k,
            x: x_t,
            xx: x_tt,
            xxx: x_ttt,
        },
        Poly3 {
            k: y_k,
            x: y_t,
            xx: y_tt,
            xxx: y_ttt,
        },
    )
}

#[test]
fn test_parametric_cubic() {
    use super::evaluate;
    use cgmath::assert_relative_eq;

    let curve1 = [
        Point2::new(1_f64, 0.),
        Point2::new(5., 0.),
        Point2::new(5., 2.),
        Point2::new(4., 3.),
    ];

    let (p_curve1_x, p_curve1_y) = parametric_cubic(curve1);

    for i in 0..10 {
        let t = (i as f64) / 10.;

        let p = evaluate(&curve1, t);
        let p_x = p_curve1_x.eval(t);
        let p_y = p_curve1_y.eval(t);

        assert_relative_eq!(p.x, p_x, epsilon = 1e-8);
        assert_relative_eq!(p.y, p_y, epsilon = 1e-8);
    }
}
