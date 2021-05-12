//! Utilities for dealing with polynomials

use cgmath::Zero;
use std::ops;

/// Third degree polynomial
#[derive(Debug, Clone, Copy, Default)]
#[repr(C)]
pub struct Poly3<S> {
    // SAFETY: do not modify layout!
    pub k: S,
    pub x: S,
    pub xx: S,
    pub xxx: S,
}

/// First degree polynomial in 2 dimensions
#[derive(Debug, Clone, Copy, Default)]
pub struct Poly1x2d<S> {
    pub k: S,
    pub x: S,
    pub y: S,
}

/// Second degree polynomial in 2 dimensions
#[derive(Debug, Clone, Copy, Default)]
pub struct Poly2x2d<S> {
    pub k: S,
    pub x: S,
    pub y: S,
    pub xy: S,
    pub xx: S,
    pub yy: S,
}

/// Third degree polynomial in 2 dimensions
#[derive(Debug, Clone, Copy, Default)]
pub struct Poly3x2d<S> {
    pub k: S,
    pub x: S,
    pub y: S,
    pub xy: S,
    pub xx: S,
    pub yy: S,
    pub xxy: S,
    pub xyy: S,
    pub xxx: S,
    pub yyy: S,
}

impl<S> ops::Add for Poly1x2d<S>
where
    S: ops::Add<S, Output = S>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Poly1x2d {
            k: self.k + rhs.k,
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl<S> ops::Sub for Poly1x2d<S>
where
    S: ops::Sub<S, Output = S>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Poly1x2d {
            k: self.k - rhs.k,
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<S> ops::Mul<S> for Poly1x2d<S>
where
    S: ops::Mul<S, Output = S> + Copy,
{
    type Output = Self;
    fn mul(self, rhs: S) -> Self {
        Poly1x2d {
            k: self.k * rhs,
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl<S> ops::Mul<Self> for Poly1x2d<S>
where
    S: ops::Mul<S, Output = S> + Copy + Zero,
{
    type Output = Poly2x2d<S>;
    fn mul(self, rhs: Self) -> Self::Output {
        Poly2x2d {
            k: self.k * rhs.k,
            x: self.k * rhs.x + self.x * rhs.k,
            y: self.k * rhs.y + self.y * rhs.k,
            xx: self.x * rhs.x,
            xy: self.x * rhs.y + self.y * rhs.x,
            yy: self.y * rhs.y,
        }
    }
}

impl<S> ops::Mul<Poly2x2d<S>> for Poly1x2d<S>
where
    S: ops::Mul<S, Output = S> + Copy + Zero,
{
    type Output = Poly3x2d<S>;
    fn mul(self, rhs: Poly2x2d<S>) -> Self::Output {
        let Poly1x2d { k, x, y } = self;
        Poly3x2d {
            k: rhs.k * k,
            x: rhs.k * x + rhs.x * k,
            y: rhs.k * y + rhs.y * k,
            xx: rhs.x * x + rhs.xx * k,
            xy: rhs.x * y + rhs.y * x + rhs.xy * k,
            yy: rhs.y * y + rhs.yy * k,
            xxy: rhs.xx * y + rhs.xy * x,
            xyy: rhs.xy * y + rhs.yy * x,
            xxx: rhs.xx * x,
            yyy: rhs.yy * y,
        }
    }
}

impl<S> ops::Add for Poly2x2d<S>
where
    S: ops::Add<S, Output = S>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Poly2x2d {
            k: self.k + rhs.k,
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            xx: self.xx + rhs.xx,
            xy: self.xy + rhs.xy,
            yy: self.yy + rhs.yy,
        }
    }
}

impl<S> ops::Sub for Poly2x2d<S>
where
    S: ops::Sub<S, Output = S>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Poly2x2d {
            k: self.k - rhs.k,
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            xx: self.xx - rhs.xx,
            xy: self.xy - rhs.xy,
            yy: self.yy - rhs.yy,
        }
    }
}

impl<S> ops::Add for Poly3x2d<S>
where
    S: ops::Add<S, Output = S>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Poly3x2d {
            k: self.k + rhs.k,
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            xx: self.xx + rhs.xx,
            xy: self.xy + rhs.xy,
            yy: self.yy + rhs.yy,
            xxy: self.xxy + rhs.xxy,
            xyy: self.xyy + rhs.xyy,
            xxx: self.xxx + rhs.xxx,
            yyy: self.yyy + rhs.yyy,
        }
    }
}

impl<S> ops::Sub for Poly3x2d<S>
where
    S: ops::Sub<S, Output = S>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Poly3x2d {
            k: self.k - rhs.k,
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            xx: self.xx - rhs.xx,
            xy: self.xy - rhs.xy,
            yy: self.yy - rhs.yy,
            xxy: self.xxy - rhs.xxy,
            xyy: self.xyy - rhs.xyy,
            xxx: self.xxx - rhs.xxx,
            yyy: self.yyy - rhs.yyy,
        }
    }
}

impl<S> Poly3<S>
where
    S: ops::Add<S, Output = S> + ops::Mul<S, Output = S> + Copy,
{
    pub fn eval(&self, x: S) -> S {
        self.k + self.x * x + self.xx * x * x + self.xxx * x * x * x
    }

    pub fn as_slice(&self) -> &[S; 4] {
        // SAFETY: struct is repr(C) so the layout is the same
        unsafe { &*(self as *const Self as *const [S; 4]) }
    }
}

impl<S> Poly3x2d<S>
where
    S: ops::Add<S, Output = S> + ops::Mul<S, Output = S> + Copy,
{
    pub fn eval(&self, x: S, y: S) -> S {
        self.k
            + self.x * x
            + self.y * y
            + self.xx * x * x
            + self.xy * x * y
            + self.yy * y * y
            + self.xxy * x * x * y
            + self.xyy * x * y * y
            + self.xxx * x * x * x
            + self.yyy * y * y * y
    }
}

impl<S> Poly3x2d<S>
where
    S: ops::Add<S, Output = S> + ops::Mul<S, Output = S> + Zero + Copy,
{
    pub fn subst(&self, x: &[S], y: &[S]) -> Vec<S> {
        let mut out = vec![self.k];
        let mut add_fac = |deg: usize, a: S| {
            if out.len() < deg + 1 {
                out.resize(deg + 1, S::zero());
            }
            out[deg] = out[deg] + a;
        };
        let d = |x: &S| *x;
        for (deg, xa) in x.iter().map(d).enumerate() {
            add_fac(deg, self.x * xa);
            for (deg2, xb) in x.iter().map(d).enumerate() {
                add_fac(deg + deg2, self.xx * xa * xb);
                for (deg3, xc) in x.iter().map(d).enumerate() {
                    add_fac(deg + deg2 + deg3, self.xxx * xa * xb * xc);
                }
            }
        }
        for (deg, ya) in y.iter().map(d).enumerate() {
            add_fac(deg, self.y * ya);
            for (deg2, yb) in y.iter().map(d).enumerate() {
                add_fac(deg + deg2, self.yy * ya * yb);
                for (deg3, yc) in y.iter().map(d).enumerate() {
                    add_fac(deg + deg2 + deg3, self.yyy * ya * yb * yc);
                }
            }
        }
        for (deg, xa) in x.iter().map(d).enumerate() {
            for (deg2, ya) in y.iter().map(d).enumerate() {
                add_fac(deg + deg2, self.xy * xa * ya);
                for (deg3, xb) in x.iter().map(d).enumerate() {
                    add_fac(deg + deg2 + deg3, self.xxy * xa * xb * ya);
                }
                for (deg3, yb) in y.iter().map(d).enumerate() {
                    add_fac(deg + deg2 + deg3, self.xyy * xa * ya * yb);
                }
            }
        }
        out
    }
}

#[test]
fn test_poly3x2d_subst() {
    use cgmath::assert_abs_diff_eq;

    let p = Poly3x2d {
        k: 3.,
        x: 2.,
        y: 5.,
        xy: -5.,
        xx: -24.,
        yy: 3.,
        xxy: 9.,
        xyy: -16.,
        xxx: -44.,
        yyy: 1.,
    };
    let px = Poly3 {
        k: -4.,
        x: 1.,
        xx: 2.,
        xxx: 3.,
    };
    let py = Poly3 {
        k: 4.,
        x: -2.,
        xx: 7.,
        xxx: 1.,
    };

    let f_ref = |t: f64| p.eval(px.eval(t), py.eval(t));
    let f_subst = p.subst(px.as_slice(), py.as_slice());

    fn eval_poly(f: &[f64], t: f64) -> f64 {
        let mut out = 0.;
        for (i, a) in f.iter().enumerate() {
            out += a * t.powf(i as f64);
        }
        out
    }

    for i in 0..10 {
        let t = i as f64;
        assert_abs_diff_eq!(f_ref(t), eval_poly(&f_subst, t));
    }
}
