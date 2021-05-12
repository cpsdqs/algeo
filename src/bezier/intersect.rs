use super::{evaluate, implicit_cubic, parametric_cubic};
use cgmath::{BaseFloat, Point2};

/// Finds intersections of two 2D cubic bézier curves.
///
/// Curve intersections are given as a t parameter for the curve given in `a` and a 2D point.
///
/// # Panics
/// - if S is not isomorphic to f64
///
/// # Details
/// Intersection is implemented using the implicitization method outlined in chapter 17 of
/// "Computer Aided Geometric Design" by Thomas. W. Sederberg and in the paper
/// "Implicitization, Inversion, and Intersection of Planar Rational Cubic Curves" by Sederberg,
/// Anderson, and Goldman in *Computer Vision, Graphics, and Image Processing* (1985).
///
/// This method creates an implicit function for the second curve (see [`implicit_cubic`]),
/// substitutes the first curve into the implicit function in parametric form (see
/// [`parametric_cubic`]), and solves for the equation's roots (see [`roots::find_roots_eigen`]).
pub fn intersect_cubic<S>(
    a: [Point2<S>; 4],
    b: [Point2<S>; 4],
) -> impl Iterator<Item = (f64, Point2<S>)>
where
    S: BaseFloat + 'static,
{
    let b_implicit = implicit_cubic(b);
    let (a_x, a_y) = parametric_cubic(a);

    let polynomial = b_implicit
        .subst(a_x.as_slice(), a_y.as_slice())
        .into_iter()
        .map(|s| s.to_f64().unwrap())
        .collect();

    roots::find_roots_eigen(polynomial)
        .into_iter()
        .filter(|t| *t >= 0. && *t <= 1.)
        .map(move |t| (
            t,
            evaluate(&a, S::from(t).unwrap()),
        ))
}

#[test]
fn test_intersect_cubic() {
    use cgmath::assert_relative_eq;

    let curve1 = [
        Point2::new(0., 0.),
        Point2::new(5., 11.),
        Point2::new(7., 2.),
        Point2::new(16., 0.),
    ];
    let curve2 = [
        Point2::new(1., 6.),
        Point2::new(2., 0.),
        Point2::new(14., 10.),
        Point2::new(11., 1.),
    ];
    let mut ips = intersect_cubic(curve1, curve2).collect::<Vec<_>>();
    let mut ips2 = intersect_cubic(curve2, curve1).collect::<Vec<_>>();

    ips.sort_by(|a, b| a.1.x.partial_cmp(&b.1.x).unwrap());
    ips2.sort_by(|a, b| a.1.x.partial_cmp(&b.1.x).unwrap());

    let p_ref = vec![
        Point2::new(2.43, 4.11),
        Point2::new(7.12, 4.54),
        Point2::new(11.26, 1.88),
    ];

    assert_eq!(p_ref.len(), ips.len());
    assert_eq!(ips.len(), ips2.len());

    for (i, (j, k)) in p_ref.iter().zip(ips.iter().zip(ips2.iter())) {
        assert_relative_eq!(*i, j.1, epsilon = 0.02);
        assert_relative_eq!(j.1, k.1, epsilon = 1e-5);
    }
}
