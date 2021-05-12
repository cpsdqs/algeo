use cgmath::{BaseNum, EuclideanSpace};

/// A finite list of items (points or vectors) representing a bézier curve. Generalizes over [T; n]
pub trait BezierCurve<T> {
    /// Derivative type.
    type Derivative: BezierCurve<T>;
    /// Returns the number of items.
    fn count(&self) -> usize;
    /// Gets an item at the specified index.
    fn get(&self, i: usize) -> &T;
    /// Sets an item at the specified index.
    fn set(&mut self, i: usize, v: T);
    /// Returns the same curve.
    fn clone(&self) -> Self;
    /// Returns the same list, but shorter by one element.
    fn reduced(&self) -> Self::Derivative;
}

/// Type conversion trait: converts a point type to a vector type
///
/// TODO: use GAT in FiniteList instead when they’re stable
pub trait DerivativeSpace<I> {
    fn from_integral(this: I) -> Self;
}

macro_rules! id_derivative_space {
    ($($ty:ty),+) => {
        $(impl DerivativeSpace<$ty> for $ty {
            fn from_integral(this: $ty) -> Self {
                this
            }
        })+
    }
}
id_derivative_space!(i8, u8, i16, u16, i32, u32, i64, u64, isize, usize, f32, f64);
impl<S> DerivativeSpace<cgmath::Vector1<S>> for cgmath::Vector1<S> {
    fn from_integral(this: Self) -> Self {
        this
    }
}

impl<S> DerivativeSpace<cgmath::Vector2<S>> for cgmath::Vector2<S> {
    fn from_integral(this: Self) -> Self {
        this
    }
}

impl<S> DerivativeSpace<cgmath::Vector3<S>> for cgmath::Vector3<S> {
    fn from_integral(this: Self) -> Self {
        this
    }
}

impl<S> DerivativeSpace<cgmath::Vector4<S>> for cgmath::Vector4<S> {
    fn from_integral(this: Self) -> Self {
        this
    }
}

impl<S> DerivativeSpace<cgmath::Point1<S>> for cgmath::Vector1<S>
where
    S: BaseNum,
{
    fn from_integral(this: cgmath::Point1<S>) -> Self {
        this.to_vec()
    }
}

impl<S> DerivativeSpace<cgmath::Point2<S>> for cgmath::Vector2<S>
where
    S: BaseNum,
{
    fn from_integral(this: cgmath::Point2<S>) -> Self {
        this.to_vec()
    }
}

impl<S> DerivativeSpace<cgmath::Point3<S>> for cgmath::Vector3<S>
where
    S: BaseNum,
{
    fn from_integral(this: cgmath::Point3<S>) -> Self {
        this.to_vec()
    }
}

impl<T> BezierCurve<T> for [T; 0] {
    type Derivative = Self;
    fn count(&self) -> usize {
        0
    }
    fn get(&self, i: usize) -> &T {
        &self[i]
    }
    fn set(&mut self, _i: usize, _v: T) {}
    fn clone(&self) -> Self {
        []
    }
    fn reduced(&self) -> Self::Derivative {
        []
    }
}

impl<T, U> DerivativeSpace<[U; 0]> for [T; 0]
where
    T: DerivativeSpace<U>,
{
    fn from_integral(_: [U; 0]) -> Self {
        []
    }
}

macro_rules! impl_bezier_curve {
    (m $($n:expr),+) => {
        $(impl_bezier_curve!($n);)+
    };
    ($N:expr) => {
        impl<T> BezierCurve<T> for [T; $N] where T: Clone + Copy {
            type Derivative = [T; $N - 1];

            fn count(&self) -> usize {
                $N
            }

            fn get(&self, i: usize) -> &T {
                &self[i]
            }

            fn set(&mut self, i: usize, v: T) {
                self[i] = v;
            }

            fn clone(&self) -> Self {
                *self
            }

            fn reduced(&self) -> Self::Derivative {
                let mut res = [self[0]; $N - 1];
                for i in 1..($N - 1) {
                    res[i] = self[i];
                }
                res
            }
        }
        impl<T, U> DerivativeSpace<[U; $N]> for [T; $N]
        where
            T: DerivativeSpace<U> + Copy,
            U: Copy,
        {
            fn from_integral(this: [U; $N]) -> Self {
                let mut res = [T::from_integral(this[0]); $N];
                for i in 1..$N {
                    res[i] = T::from_integral(this[i]);
                }
                res
            }
        }
    };
}
impl_bezier_curve!(m 1, 2, 3, 4, 5, 6);
