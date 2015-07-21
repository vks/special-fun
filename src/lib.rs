extern crate num;

extern "C" {
    // Floating point numeric utilities
    /// Round to nearest or event integer valued f64.
    fn round(x: f64) -> f64;
    /// Largest integer less than or equal to x.
    fn floor(x: f64) -> f64;
    /// Smallest integer greater than or equal to x.
    fn ceil(x: f64) -> f64;
    /// Return the significand between 0.5 and 1. Write exponent to expnt.
    /// x = y * 2**expn
    fn frexp(x: f64, expnt: &mut i32) -> f64;
    /// Multiply x by 2**n.
    fn ldexp(x: f64, n: i32) -> f64;
    /// Absolute value.
    fn fabs(x: f64) -> f64;
    /// Return 1 if the sign bit of x is 1, else 0.
    fn signbit(x: f64) -> i32;
    /// Return 1 if x is NaN, else 0.
    fn isnan(x: f64) -> i32;
    /// Return 1 if x is finite, else 0.
    fn isfinite(x: f64) -> i32;

    // Roots
    /// Cube root.
    fn cbrt(x: f64) -> f64;
    /// Square root.
    fn sqrt(x: f64) -> f64;
    /// Integer square root.
    fn lsqrt(x: i64) -> i64;

    // Exponential functions
    /// Exponential function.
    fn exp(x: f64) -> f64;
    /// Base 10 exponential function.
    fn exp10(x: f64) -> f64;
    /// Base 2 exponential function.
    fn exp2(x: f64) -> f64;
    /// Exponential of squared argument.
    fn expx2(x: f64, sign: i32) -> f64;
    /// Exponential integral.
    fn ei(x: f64) -> f64;
    /// Error function.
    fn erf(x: f64) -> f64;
    /// Complementary error function.
    fn erfc(x: f64) -> f64;
    /// Power function.
    fn pow(x: f64, y: f64) -> f64;
    /// Integer power function.
    fn powi(x: f64, n: i32) -> f64;

    // Logarithmic functions
    /// Natural logarithm.
    fn log(x: f64) -> f64;
    /// Common logarithm.
    fn log10(x: f64) -> f64;
    /// Base 2 logarithm.
    fn log2(x: f64) -> f64;
    /// Dilogarithm (Spence's function).
    fn spence(x: f64) -> f64;

    // Trigonometric functions
    /// Circular sine.
    fn sin(x: f64) -> f64;
    /// Circular cosine.
    fn cos(x: f64) -> f64;
    /// Circular tangent.
    fn tan(x: f64) -> f64;
    /// Inverse circular sine.
    fn asin(x: f64) -> f64;
    /// Inverse circular cosine.
    fn acos(x: f64) -> f64;
    /// Inverse circular tangent.
    fn atan(x: f64) -> f64;
    /// Quadrant-correct inverse circular tangent.
    fn atan2(y: f64, x: f64) -> f64;
    /// Sine and cosine integrals.
    fn sici(x: f64, si: &mut f64, ci: &mut f64) -> f64;

    // Hyperbolic functions
    /// Hyperbolic sine.
    fn sinh(x: f64) -> f64;
    /// Hyperbolic cosine.
    fn cosh(x: f64) -> f64;
    /// Hyperbolic tangent.
    fn tanh(x: f64) -> f64;
    /// Inverse hyperbolic sine.
    fn asinh(x: f64) -> f64;
    /// Inverse hyperbolic cosine.
    fn acosh(x: f64) -> f64;
    /// Inverse hyperbolic tangent.
    fn atanh(x: f64) -> f64;
    /// Hyperbolic sine and cosine integrals.
    fn shichi(x: f64, chi: &mut f64, shi: &mut f64);

    // Beta functions
    /// Beta function.
    fn beta(a: f64, b: f64) -> f64;
    /// Regularized incomplete beta function.
    fn incbet(a: f64, b: f64, x: f64) -> f64;
    /// Inverse of incomplete beta integral.
    fn incbi(a: f64, b: f64, y: f64) -> f64;

    // Gamma functions
    /// Gamma function.
    fn gamma(x: f64) -> f64;
    /// Reciprocal gamma function.
    fn rgamma(x: f64) -> f64;
    /// Natural logarithm of gamma function.
    fn lgam(x: f64) -> f64;
    /// Regularized incomplete gamma integral.
    fn igam(a: f64, x: f64) -> f64;
    /// Complemented incomplete gamma integral.
    fn igamc(a: f64, x: f64) -> f64;
    /// Inverse of complemented incomplete gamma integral.
    fn igami(a: f64, p: f64) -> f64;
    /// Psi (digamma) function.
    fn psi(x: f64) -> f64;
    /// Factorial function.
    fn fac(i: i32) -> f64;

    // Bessel functions
    /// Bessel function of order zero.
    fn j0(x: f64) -> f64;
    /// Bessel function of order one.
    fn j1(x: f64) -> f64;
    /// Bessel function of integer order.
    fn jn(n: i32, x: f64) -> f64;
    /// Bessel function of real order.
    fn jv(n: f64, x: f64) -> f64;

    /// Bessel function of the second kind, order zero.
    fn y0(x: f64) -> f64;
    /// Bessel function of the second kind, order one.
    fn y1(x: f64) -> f64;
    /// Bessel function of the second kind, integer order.
    fn yn(n: i32, x: f64) -> f64;
    /// Bessel function of the second kind, real order.
    fn yv(v: f64, x: f64) -> f64;

    /// Modified Bessel function of order zero.
    fn i0(x: f64) -> f64;
    /// Modified Bessel function of order zero, exponentially scaled.
    fn i0e(x: f64) -> f64;
    /// Modified Bessel function of order one.
    fn i1(x: f64) -> f64;
    /// Modified Bessel function of order one, exponentially scaled.
    fn i1e(x: f64) -> f64;
    /// Modified Bessel function of real order.
    fn iv(v: f64, x: f64) -> f64;

    /// Modified Bessel function of the third kind, order zero.
    fn k0(x: f64) -> f64;
    /// Modified Bessel function of the third kind, order zero,
    /// exponentially scaled.
    fn k0e(x: f64) -> f64;
    /// Modified Bessel function of the third kind, order one.
    fn k1(x: f64) -> f64;
    /// Modified Bessel function of the third kind, order one,
    /// exponentially scaled.
    fn k1e(x: f64) -> f64;
    /// Modified Bessel function of the third kind, integer order.
    fn kn(n: i32, x: f64) -> f64;

    // Elliptic functions
    /// Incomplete elliptic integral of the first kind.
    fn ellik(phi: f64, m: f64) -> f64;
    /// Incomplete elliptic integral of the second kind.
    fn ellie(phi: f64, m: f64) -> f64;
    /// Complete elliptic integral of the first kind.
    fn ellpk(m1: f64) -> f64;
    /// Complete elliptic integral of the second kind.
    fn ellpe(m1: f64) -> f64;
    /// Jacobian elliptic function.
    fn ellpj(u: f64, m: f64, sn: &mut f64, cn: &mut f64, dn: &mut f64, phi: &mut f64) -> i32;

    // Hypergeometric functions
    /// Gauss hypergeometric function 2F1.
    fn hyp2f1(a: f64, b: f64, c: f64, x: f64) -> f64;
    /// Confluent hypergeometric function.
    fn hyperg(a: f64, b: f64, x: f64) -> f64;

    // Distributions
    /// Binomial distribution.
    fn bdtr(k: i32, n: i32, p: f64) -> f64;
    /// Complemented binomial distribution.
    fn bdtrc(k: i32, n: i32, p: f64) -> f64;
    /// Inverse of binomial distribution.
    fn bdtri(k: i32, n: i32, y: f64) -> f64;

    /// Negative binomial distribution.
    fn nbdtr(k: i32, n: i32, p: f64) -> f64;
    /// Complemented negative binomial distribution.
    fn nbdtrc(k: i32, n: i32, p: f64) -> f64;
    /// Inverse of negative binomial distribution.
    fn nbdtri(k: i32, n: i32, p: f64) -> f64;

    /// Beta distribution.
    fn btdtr(a: f64, b: f64, x: f64) -> f64;

    /// Chi-square distribution.
    fn chdtr(df: f64, x: f64) -> f64;
    /// Complemented chi-square distribution.
    fn chdtrc(v: f64, x: f64) -> f64;
    /// Inverse of complemented chi-square distribution.
    fn chdtri(df: f64, y: f64) -> f64;

    /// F distribution.
    fn fdtr(df1: i32, df2: i32, x: f64) -> f64;
    /// Complemented F distribution.
    fn fdtrc(df1: i32, df2: i32, x: f64) -> f64;
    /// Inverse of complemented F distribution.
    fn fdtri(df1: i32, df2: i32, p: f64) -> f64;

    /// Gamma distribution.
    fn gdtr(a: f64, b: f64, x: f64) -> f64;
    /// Complemented gamma distribution.
    fn gdtrc(a: f64, b: f64, x: f64) -> f64;

    /// Normal distribution.
    fn ndtr(x: f64) -> f64;
    /// Inverse of normal distribution.
    fn ndtri(y: f64) -> f64;

    /// Poisson distribution.
    fn pdtr(k: i32, m: f64) -> f64;
    /// Complemented Poisson distribution.
    fn pdtrk(k: i32, m: f64) -> f64;
    /// Inverse of Poisson distribution.
    fn pdtri(k: i32, y: f64) -> f64;

    /// Student's t distribution.
    fn stdtr(k: i16, t: f64) -> f64;
    /// Inverse of Student's t distribution.
    fn stdtri(k: i32, p: f64) -> f64;

    // Misc special functions
    /// Airy function.
    fn airy(x: f64, ai: &mut f64, aip: &mut f64, bi: &mut f64, bip: &mut f64) -> i32;
    /// Dawson's integral.
    fn dawsn(x: f64) -> f64;
    /// Fresnel integral.
    fn fresnl(x: f64, s: &mut f64, c: &mut f64);
    /// Integral of Planck's black body radiation formula.
    fn plancki(lambda: f64, temperature: f64) -> f64;
    /// Struve function.
    fn struve(v: f64, x: f64) -> f64;
    /// Riemann zeta function.
    fn zetac(x: f64) -> f64;
    /// Riemann zeta function of two arguments.
    fn zeta(x: f64, q: f64) -> f64;
}

/// Special functions on primitive floating point numbers.
pub trait FloatSpecial {
    /// Beta function.
    fn beta(self, b: Self) -> Self;
    /// Regularized incomplete beta function.
    fn betainc(self, a: Self, b: Self) -> Self;
    /// Inverse of incomplete beta integral.
    fn betainc_inv(self, a: Self, b: Self) -> Self;

    /// Factorial.
    fn factorial(self) -> Self;
    /// Gamma function.
    fn gamma(self) -> Self;
    /// Reciprocal gamma function.
    fn rgamma(self) -> Self;
    /// Logarithm of gamma function.
    fn loggamma(self) -> Self;

    /// Regularized incomplete gamma integral.
    fn gammainc(self, a: Self) -> Self;
    /// Complemented incomplete gamma integral.
    fn gammac(self, a: Self) -> Self;
    /// Inverse of complemented incomplete gamma integral.
    fn gammac_inv(self, a: Self) -> Self;

    /// Digamma function.
    fn digamma(self) -> Self;

    /// Normal distribution function.
    fn norm(self) -> Self;
    /// Inverse of Normal distribution function.
    fn norm_inv(self) -> Self;

    /// Bessel function of real order of the first kind.
    fn besselj(self, v: Self) -> Self;
    /// Bessel function of real order of the second kind.
    fn bessely(self, v: Self) -> Self;
    /// Modified bessel function of real order of the first kind.
    fn besseli(self, v: Self) -> Self;
    /// Modified bessel function of integer order of the second kind.
    fn besselk(self, v: i32) -> Self;

    /// Riemann zeta function.
    fn riemann_zeta(self) -> Self;
    /// Hurwitz zeta function.
    fn hurwitz_zeta(self, q: Self) -> Self;
}

impl FloatSpecial for f64 {
    fn beta(self, b: f64) -> f64 {
        unsafe { beta(self, b) }
    }
    fn betainc(self, a: f64, b: f64) -> f64 {
        unsafe { incbet(a, b, self) }
    }
    fn betainc_inv(self, a: f64, b: f64) -> f64 {
        unsafe { incbi(a, b, self) }
    }

    fn factorial(self) -> f64 {
        unsafe { gamma(self + 1.0) }
    }
    fn gamma(self) -> f64 {
        unsafe { gamma(self) }
    }
    fn rgamma(self) -> f64 {
        unsafe { rgamma(self) }
    }
    fn loggamma(self) -> f64 {
        unsafe { lgam(self) }
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

    fn digamma(self) -> f64 {
        unsafe { psi(self) }
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
    fn bessely(self, v: f64) -> f64 {
        unsafe { yv(v, self) }
    }
    fn besseli(self, v: f64) -> f64 {
        unsafe { iv(v, self) }
    }
    fn besselk(self, n: i32) -> f64 {
        unsafe { kn(n, self) }
    }

    fn riemann_zeta(self) -> f64 {
        unsafe { 1.0 + zetac(self) }
    }
    fn hurwitz_zeta(self, q: f64) -> f64 {
        unsafe { zeta(self, q) }
    }
}

/*
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
*/

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
            assert_almost_eq(5.0f64.beta(2.0), 0.03333333333333);
            assert_almost_eq(1.5f64.beta(2.0), 0.26666666666666);
            assert_almost_eq(0.5f64.betainc(2.0, 3.0), 0.6875);
            assert_almost_eq(0.6875f64.betainc_inv(2.0, 3.0), 0.5);
        }

        #[test]
        fn factorial() {
            assert_almost_eq(0.5f64.factorial(), 0.886226925452758);
        }

        #[test]
        fn gamma() {
            assert_almost_eq(3.5f64.gamma(), 3.32335097044784);
            assert_eq!(1.0f64.rgamma(), 1.0);
            assert_almost_eq(4.0f64.rgamma(), 0.1666666666666666666);
            assert_almost_eq(13.2f64.loggamma(), 20.49400419456603678498394);
            assert_almost_eq(4.0f64.gammainc(2.0), 0.90842180555632912);
            assert_almost_eq(1.0 - 4.0f64.gammainc(2.0), 4.0f64.gammac(2.0));
            assert_almost_eq(4.0f64.gammac(2.0).gammac_inv(2.0), 4.0);
        }

        #[test]
        fn digamma() {
            assert_almost_eq(1.0f64.digamma(), -0.5772156649015328606065121);
        }

        #[test]
        fn norm() {
            assert_almost_eq(2.0f64.norm(), 0.9772499);
            assert_almost_eq(0.9f64.norm_inv(), 1.281552);
        }

        #[test]
        fn besselj() {
            assert_almost_eq(10.0f64.besselj(2.0), 0.25463031368512062);
            assert_almost_eq(1000.0f64.besselj(2.0), -0.024777229528606);
            assert_almost_eq(0.75f64.besselj(4.0), 0.000801070086542314);
        }

        #[test]
        fn bessely() {
            assert_almost_eq(
                3.141592653589793f64.bessely(1.0), 0.3588729167767189594679827);
        }

        #[test]
        fn besseli() {
            assert_eq!(0.0f64.besseli(0.0), 1.0);
            assert_eq!(0.0f64.besseli(1.0), 0.0);
            assert_almost_eq(1.0f64.besseli(0.0), 1.266065877752008335598245);
        }

        #[test]
        fn besselk() {
            assert_almost_eq(1.0f64.besselk(0), 0.4210244382407083333356274);
            assert_almost_eq(100.0f64.besselk(0), 4.656628229175902018939005e-45);
        }

        #[test]
        fn riemann_zeta() {
            assert_almost_eq(2.0f64.riemann_zeta(), 1.64493406684822);
            assert_eq!(0.0f64.riemann_zeta(), -0.5);
            //assert_almost_eq(-1.0f64.riemann_zeta(), -0.0833333333333);
            assert_almost_eq(50.0f64.riemann_zeta(), 1.0);
        }

        #[test]
        fn hurwitz_zeta() {
            assert_almost_eq(2.0f64.hurwitz_zeta(3.0), 0.3949340668482);
            //assert_almost_eq(0.0f64.hurwitz_zeta(10.0), -9.5);
        }
    }

    /*
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
    */
}
