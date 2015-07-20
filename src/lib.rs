extern crate num;

extern "C" {
    /// Regularized incomplete beta function.
    fn incbet(a: f64, b: f64, x: f64) -> f64;
    /// Inverse of incomplete beta integral.
    fn incbi(a: f64, b: f64, y: f64) -> f64;
    /// Regularized incomplete gamma integral.
    fn igam(a: f64, x: f64) -> f64;
    /// Complemented incomplete gamma integral.
    fn igamc(a: f64, x: f64) -> f64;
    /// Inverse of complemented incomplete gamma integral.
    fn igami(a: f64, p: f64) -> f64;
    /// Normal distribution function.
    fn ndtr(x: f64) -> f64;
    /// Inverse of Normal distribution function.
    fn ndtri(x: f64) -> f64;
    /// Bessel function of non-integer order.
    fn jv(v: f64, x: f64) -> f64;

    /// Regularized incomplete beta function.
    fn incbetf(a: f32, b: f32, x: f32) -> f32;
    /// Inverse of incomplete beta integral.
    fn incbif(a: f32, b: f32, y: f32) -> f32;
    /// Regularized incomplete gamma integral.
    fn igamf(a: f32, x: f32) -> f32;
    /// Complemented incomplete gamma integral.
    fn igamcf(a: f32, x: f32) -> f32;
    /// Inverse of complemented incomplete gamma integral.
    fn igamif(a: f32, p: f32) -> f32;
    /// Normal distribution function.
    fn ndtrf(x: f32) -> f32;
    /// Inverse of Normal distribution function.
    fn ndtrif(x: f32) -> f32;
    /// Bessel function of non-integer order.
    fn jvf(v: f32, x: f32) -> f32;
}

/// Special functions on primitive floating point numbers.
/// These are essential for most statistical applications.
pub trait FloatSpecial {
    /// Regularized incomplete beta function.
    fn betainc(self, a: Self, b: Self) -> Self;
    /// Inverse of incomplete beta integral.
    fn betainc_inv(self, a: Self, b: Self) -> Self;
    /// Regularized incomplete gamma integral.
    fn gammainc(self, a: Self) -> Self;
    /// Complemented incomplete gamma integral.
    fn gammac(self, a: Self) -> Self;
    /// Inverse of complemented incomplete gamma integral.
    fn gammac_inv(self, a: Self) -> Self;
    /// Normal distribution function.
    fn norm(self) -> Self;
    /// Inverse of Normal distribution function.
    fn norm_inv(self) -> Self;
    /// Bessel function of non-integer order of the first kind.
    fn besselj(self, v: Self) -> Self;
}

impl FloatSpecial for f64 {
    fn betainc(self, a: f64, b: f64) -> f64 {
        unsafe { incbet(a, b, self) }
    }
    fn betainc_inv(self, a: f64, b: f64) -> f64 {
        unsafe { incbi(a, b, self) }
    }
    fn gammainc(self, a: f64) -> f64 {
        unsafe { igam(a, self) }
    }
    fn gammac(self, a: f64) -> f64 {
        unsafe { igamc(a, self) }
    }
    fn gammac_inv(self, a: f64) -> f64 {
        unsafe { igami(a, self) }
    }
    fn norm(self) -> f64 {
        unsafe { ndtr(self) }
    }
    fn norm_inv(self) -> f64 {
        unsafe { ndtri(self) }
    }
    fn besselj(self, v: f64) -> f64 {
        unsafe { jv(v, self) }
    }
}

impl FloatSpecial for f32 {
    fn betainc(self, a: f32, b: f32) -> f32 {
        unsafe { incbetf(a, b, self) }
    }
    fn betainc_inv(self, a: f32, b: f32) -> f32 {
        unsafe { incbif(a, b, self) }
    }
    fn gammainc(self, a: f32) -> f32 {
        unsafe { igamf(a, self) }
    }
    fn gammac(self, a: f32) -> f32 {
        unsafe { igamcf(a, self) }
    }
    fn gammac_inv(self, a: f32) -> f32 {
        unsafe { igamif(a, self) }
    }
    fn norm(self) -> f32 {
        unsafe { ndtrf(self) }
    }
    fn norm_inv(self) -> f32 {
        unsafe { ndtrif(self) }
    }
    fn besselj(self, v: f32) -> f32 {
        unsafe { jvf(v, self) }
    }
}

#[cfg(test)]
mod test {
    use ::std::fmt::Debug;
    use ::num::traits::{Float, FromPrimitive};

    fn assert_almost_eq<T: Float + FromPrimitive + Debug>(a: T, b: T) {
        let tol: T = FromPrimitive::from_f32(1e-6).unwrap();
        if (a - b).abs() > tol {
            panic!("{:?} vs. {:?}", a, b);
        }
    }

    mod double {
        use super::super::FloatSpecial;
        use super::assert_almost_eq;

        #[test]
        fn beta() {
            assert_almost_eq(0.5f64.betainc(2.0, 3.0), 0.6875);
        }

        #[test]
        fn gamma() {
            assert_almost_eq(4.0f64.gammainc(2.0), 0.90842180555632912);
            assert_almost_eq(1.0 - 4.0f64.gammainc(2.0), 4.0f64.gammac(2.0));
            assert_almost_eq(4.0f64.gammac(2.0).gammac_inv(2.0), 4.0);
        }

        #[test]
        fn norm() {
            assert_almost_eq(2.0f64.norm(), 0.9772499);
            assert_almost_eq(0.9f64.norm_inv(), 1.281552);
        }

        #[test]
        fn bessel() {
            assert_almost_eq(10.0f64.besselj(2.0), 0.25463031368512062);
        }
    }

    mod single {
        use super::super::FloatSpecial;
        use super::assert_almost_eq;

        #[test]
        fn beta() {
            assert_almost_eq(0.5f32.betainc(2.0, 3.0), 0.6875);
        }

        #[test]
        fn gamma() {
            assert_almost_eq(4.0f32.gammainc(2.0), 0.90842180555632912);
            assert_almost_eq(1.0 - 4.0f32.gammainc(2.0), 4.0f32.gammac(2.0));
            assert_almost_eq(4.0f32.gammac(2.0).gammac_inv(2.0), 4.0);
        }

        #[test]
        fn norm() {
            assert_almost_eq(2.0f32.norm(), 0.9772499);
            assert_almost_eq(0.9f32.norm_inv(), 1.281552);
        }

        #[test]
        fn bessel() {
            assert_almost_eq(10.0f32.besselj(2.0), 0.25463031368512062);
        }
    }
}
