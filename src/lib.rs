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

#[cfg(test)]
mod test {
    use ::std;
    use ::num::traits::{Float, FromPrimitive};

    use super::FloatSpecial;

    fn assert_almost_eq<T: Float + FromPrimitive + std::fmt::Debug>(a: T, b: T) {
        let tol: T = FromPrimitive::from_f32(1e-6).unwrap();
        if (a - b).abs() > tol {
            panic!("{:?} vs. {:?}", a, b);
        }
    }

    #[test]
    fn test_beta() {
        assert_almost_eq(0.5f64.betainc(2.0, 3.0), 0.6875);
    }

    #[test]
    fn test_gamma() {
        assert_almost_eq(4.0f64.gammainc(2.0), 0.90842180555632912);
        assert_almost_eq(1.0 - 4.0f64.gammainc(2.0), 4.0f64.gammac(2.0));
        assert_almost_eq(4.0f64.gammac(2.0).gammac_inv(2.0), 4.0);
    }

    #[test]
    fn test_norm() {
        assert_almost_eq(2.0f64.norm(), 0.9772499);
        assert_almost_eq(0.9f64.norm_inv(), 1.281552);
    }

    #[test]
    fn test_bessel() {
        assert_almost_eq(10.0f64.besselj(2.0), 0.25463031368512062);
    }
}
