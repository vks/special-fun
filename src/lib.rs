extern crate num;

use num::traits::Float;

// double precision
#[allow(dead_code)]
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
    /// Compute accurately exponential of squared argument.
    fn expm1(x: f64) -> f64;
    /// Compute accurately exp(x) - 1 for x close to 0.
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
    /// Compute accurately log(1 + x) for x close to 0.
    fn log1p(x: f64) -> f64;
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
    /// Compute accurately cos(x) - 1 for x close to 0.
    fn cosm1(x: f64) -> f64;
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
    /// Confluent hypergeometric function 1F1.
    fn hyperg(a: f64, b: f64, x: f64) -> f64;
    /// Hypergeometric function 1F2.
    fn onef2(a: f64, b: f64, c: f64, x: f64, err: &mut f64) -> f64;
    /// Gauss hypergeometric function 2F1.
    fn hyp2f1(a: f64, b: f64, c: f64, x: f64) -> f64;
    /// Hypergeometric function 3F0.
    fn threef0(a: f64, b: f64, c: f64, x: f64, err: &mut f64) -> f64;

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

// single precision
#[allow(dead_code)]
extern "C" {
    // Floating point numeric utilities
    /// Round to nearest or event integer valued f32.
    fn roundf(x: f32) -> f32;
    /// Largest integer less than or equal to x.
    fn floorf(x: f32) -> f32;
    /// Smallest integer greater than or equal to x.
    fn ceilf(x: f32) -> f32;
    /// Return the significand between 0.5 and 1. Write exponent to expnt.
    /// x = y * 2**expn
    fn frexpf(x: f32, expnt: &mut i32) -> f32;
    /// Multiply x by 2**n.
    fn ldexpf(x: f32, n: i32) -> f32;
    /// Absolute value.
    fn fabsf(x: f32) -> f32;
    /// Return 1 if the sign bit of x is 1, else 0.
    fn signbitf(x: f32) -> i32;
    /// Return 1 if x is NaN, else 0.
    fn isnanf(x: f32) -> i32;
    /// Return 1 if x is finite, else 0.
    fn isfinitef(x: f32) -> i32;

    // Roots
    /// Cube root.
    fn cbrtf(x: f32) -> f32;
    /// Square root.
    fn sqrtf(x: f32) -> f32;
    /// Integer square root.
    fn lsqrtf(x: i64) -> i64;

    // Exponential functions
    /// Exponential function.
    fn expf(x: f32) -> f32;
    /// Base 10 exponential function.
    fn exp10f(x: f32) -> f32;
    /// Base 2 exponential function.
    fn exp2f(x: f32) -> f32;
    /// Compute accurately exponential of squared argument.
    fn expm1f(x: f32) -> f32;
    /// Compute accurately exp(x) - 1 for x close to 0.
    fn expx2f(x: f32, sign: i32) -> f32;
    /// Exponential integral.
    fn eif(x: f32) -> f32;
    /// Error function.
    fn erff(x: f32) -> f32;
    /// Complementary error function.
    fn erfcf(x: f32) -> f32;
    /// Power function.
    fn powf(x: f32, y: f32) -> f32;
    /// Integer power function.
    fn powif(x: f32, n: i32) -> f32;

    // Logarithmic functions
    /// Natural logarithm.
    fn logf(x: f32) -> f32;
    /// Common logarithm.
    fn log10f(x: f32) -> f32;
    /// Base 2 logarithm.
    fn log2f(x: f32) -> f32;
    /// Compute accurately log(1 + x) for x close to 0.
    fn log1pf(x: f32) -> f32;
    /// Dilogarithm (Spence's function).
    fn spencef(x: f32) -> f32;

    // Trigonometric functions
    /// Circular sine.
    fn sinf(x: f32) -> f32;
    /// Circular cosine.
    fn cosf(x: f32) -> f32;
    /// Circular tangent.
    fn tanf(x: f32) -> f32;
    /// Inverse circular sine.
    fn asinf(x: f32) -> f32;
    /// Inverse circular cosine.
    fn acosf(x: f32) -> f32;
    /// Inverse circular tangent.
    fn atanf(x: f32) -> f32;
    /// Quadrant-correct inverse circular tangent.
    fn atan2f(y: f32, x: f32) -> f32;
    /// Compute accurately cos(x) - 1 for x close to 0.
    fn cosm1f(x: f32) -> f32;
    /// Sine and cosine integrals.
    fn sicif(x: f32, si: &mut f32, ci: &mut f32) -> f32;

    // Hyperbolic functions
    /// Hyperbolic sine.
    fn sinhf(x: f32) -> f32;
    /// Hyperbolic cosine.
    fn coshf(x: f32) -> f32;
    /// Hyperbolic tangent.
    fn tanhf(x: f32) -> f32;
    /// Inverse hyperbolic sine.
    fn asinhf(x: f32) -> f32;
    /// Inverse hyperbolic cosine.
    fn acoshf(x: f32) -> f32;
    /// Inverse hyperbolic tangent.
    fn atanhf(x: f32) -> f32;
    /// Hyperbolic sine and cosine integrals.
    fn shichif(x: f32, chi: &mut f32, shi: &mut f32);

    // Beta functions
    /// Beta function.
    fn betaf(a: f32, b: f32) -> f32;
    /// Regularized incomplete beta function.
    fn incbetf(a: f32, b: f32, x: f32) -> f32;
    /// Inverse of incomplete beta integral.
    fn incbif(a: f32, b: f32, y: f32) -> f32;

    // Gamma functions
    /// Gamma function.
    fn gammaf(x: f32) -> f32;
    /// Reciprocal gamma function.
    fn rgammaf(x: f32) -> f32;
    /// Natural logarithm of gamma function.
    fn lgamf(x: f32) -> f32;
    /// Regularized incomplete gamma integral.
    fn igamf(a: f32, x: f32) -> f32;
    /// Complemented incomplete gamma integral.
    fn igamcf(a: f32, x: f32) -> f32;
    /// Inverse of complemented incomplete gamma integral.
    fn igamif(a: f32, p: f32) -> f32;
    /// Psi (digamma) function.
    fn psif(x: f32) -> f32;
    /// Factorial function.
    fn facf(i: i32) -> f32;

    // Bessel functions
    /// Bessel function of order zero.
    fn j0f(x: f32) -> f32;
    /// Bessel function of order one.
    fn j1f(x: f32) -> f32;
    /// Bessel function of integer order.
    fn jnf(n: i32, x: f32) -> f32;
    /// Bessel function of real order.
    fn jvf(n: f32, x: f32) -> f32;

    /// Bessel function of the second kind, order zero.
    fn y0f(x: f32) -> f32;
    /// Bessel function of the second kind, order one.
    fn y1f(x: f32) -> f32;
    /// Bessel function of the second kind, integer order.
    fn ynf(n: i32, x: f32) -> f32;
    /// Bessel function of the second kind, real order.
    fn yvf(v: f32, x: f32) -> f32;

    /// Modified Bessel function of order zero.
    fn i0f(x: f32) -> f32;
    /// Modified Bessel function of order zero, exponentially scaled.
    fn i0ef(x: f32) -> f32;
    /// Modified Bessel function of order one.
    fn i1f(x: f32) -> f32;
    /// Modified Bessel function of order one, exponentially scaled.
    fn i1ef(x: f32) -> f32;
    /// Modified Bessel function of real order.
    fn ivf(v: f32, x: f32) -> f32;

    /// Modified Bessel function of the third kind, order zero.
    fn k0f(x: f32) -> f32;
    /// Modified Bessel function of the third kind, order zero,
    /// exponentially scaled.
    fn k0ef(x: f32) -> f32;
    /// Modified Bessel function of the third kind, order one.
    fn k1f(x: f32) -> f32;
    /// Modified Bessel function of the third kind, order one,
    /// exponentially scaled.
    fn k1ef(x: f32) -> f32;
    /// Modified Bessel function of the third kind, integer order.
    fn knf(n: i32, x: f32) -> f32;

    // Elliptic functions
    /// Incomplete elliptic integral of the first kind.
    fn ellikf(phi: f32, m: f32) -> f32;
    /// Incomplete elliptic integral of the second kind.
    fn ellief(phi: f32, m: f32) -> f32;
    /// Complete elliptic integral of the first kind.
    fn ellpkf(m1: f32) -> f32;
    /// Complete elliptic integral of the second kind.
    fn ellpef(m1: f32) -> f32;
    /// Jacobian elliptic function.
    fn ellpjf(u: f32, m: f32, sn: &mut f32, cn: &mut f32, dn: &mut f32, phi: &mut f32) -> i32;

    // Hypergeometric functions
    /// Confluent hypergeometric function 1F1.
    fn hypergf(a: f32, b: f32, x: f32) -> f32;
    /// Hypergeometric function 1F2.
    fn onef2f(a: f32, b: f32, c: f32, x: f32, err: &mut f32) -> f32;
    /// Gauss hypergeometric function 2F1.
    fn hyp2f1f(a: f32, b: f32, c: f32, x: f32) -> f32;
    /// Hypergeometric function 3F0.
    fn threef0f(a: f32, b: f32, c: f32, x: f32, err: &mut f32) -> f32;

    // Distributions
    /// Binomial distribution.
    fn bdtrf(k: i32, n: i32, p: f32) -> f32;
    /// Complemented binomial distribution.
    fn bdtrcf(k: i32, n: i32, p: f32) -> f32;
    /// Inverse of binomial distribution.
    fn bdtrif(k: i32, n: i32, y: f32) -> f32;

    /// Negative binomial distribution.
    fn nbdtrf(k: i32, n: i32, p: f32) -> f32;
    /// Complemented negative binomial distribution.
    fn nbdtrcf(k: i32, n: i32, p: f32) -> f32;
    /// Inverse of negative binomial distribution.
    fn nbdtrif(k: i32, n: i32, p: f32) -> f32;

    /// Beta distribution.
    fn btdtrf(a: f32, b: f32, x: f32) -> f32;

    /// Chi-square distribution.
    fn chdtrf(df: f32, x: f32) -> f32;
    /// Complemented chi-square distribution.
    fn chdtrcf(v: f32, x: f32) -> f32;
    /// Inverse of complemented chi-square distribution.
    fn chdtrif(df: f32, y: f32) -> f32;

    /// F distribution.
    fn fdtrf(df1: i32, df2: i32, x: f32) -> f32;
    /// Complemented F distribution.
    fn fdtrcf(df1: i32, df2: i32, x: f32) -> f32;
    /// Inverse of complemented F distribution.
    fn fdtrif(df1: i32, df2: i32, p: f32) -> f32;

    /// Gamma distribution.
    fn gdtrf(a: f32, b: f32, x: f32) -> f32;
    /// Complemented gamma distribution.
    fn gdtrcf(a: f32, b: f32, x: f32) -> f32;

    /// Normal distribution.
    fn ndtrf(x: f32) -> f32;
    /// Inverse of normal distribution.
    fn ndtrif(y: f32) -> f32;

    /// Poisson distribution.
    fn pdtrf(k: i32, m: f32) -> f32;
    /// Complemented Poisson distribution.
    fn pdtrkf(k: i32, m: f32) -> f32;
    /// Inverse of Poisson distribution.
    fn pdtrif(k: i32, y: f32) -> f32;

    /// Student's t distribution.
    fn stdtrf(k: i16, t: f32) -> f32;
    /// Inverse of Student's t distribution.
    fn stdtrif(k: i32, p: f32) -> f32;

    // Misc special functions
    /// Airy function.
    fn airyf(x: f32, ai: &mut f32, aip: &mut f32, bi: &mut f32, bip: &mut f32) -> i32;
    /// Dawson's integral.
    fn dawsnf(x: f32) -> f32;
    /// Fresnel integral.
    fn fresnlf(x: f32, s: &mut f32, c: &mut f32);
    /// Integral of Planck's black body radiation formula.
    fn planckif(lambda: f32, temperature: f32) -> f32;
    /// Struve function.
    fn struvef(v: f32, x: f32) -> f32;
    /// Riemann zeta function.
    fn zetacf(x: f32) -> f32;
    /// Riemann zeta function of two arguments.
    fn zetaf(x: f32, q: f32) -> f32;
}

/// Special functions on primitive floating point numbers.
pub trait FloatSpecial: Float {
    /// Beta function.
    fn beta(self, b: Self) -> Self;
    /// Logarithm of beta function.
    fn logbeta(self, b: Self) -> Self {
        self.loggamma() + b.loggamma() - (self + b).loggamma()
    }
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

    /// Confluent hypergeometric function 1F1.
    fn hyp1f1(self, a: Self, b: Self) -> Self;
    /// Hypergeometric function 1F2.
    fn hyp1f2(self, a: Self, b: Self, c: Self) -> Self;
    /// Gauss hypergeometric function 2F1.
    fn hyp2f1(self, a: Self, b: Self, c: Self) -> Self;
    /// Hypergeometric function 3F0.
    fn hyp3f0(self, a: Self, b: Self, c: Self) -> Self;

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

    fn hyp1f1(self, a: f64, b: f64) -> f64 {
        unsafe { hyperg(a, b, self) }
    }
    fn hyp1f2(self, a: f64, b: f64, c: f64) -> f64 {
        let mut err = 0.0;
        unsafe { onef2(a, b, c, self, &mut err) }
    }
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> f64 {
        unsafe { hyp2f1(a, b, c, self) }
    }
    fn hyp3f0(self, a: f64, b: f64, c: f64) -> f64 {
        let mut err = 0.0;
        unsafe { threef0(a, b, c, self, &mut err) }
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

impl FloatSpecial for f32 {
    fn beta(self, b: f32) -> f32 {
        unsafe { betaf(self, b) }
    }
    fn betainc(self, a: f32, b: f32) -> f32 {
        unsafe { incbetf(a, b, self) }
    }
    fn betainc_inv(self, a: f32, b: f32) -> f32 {
        unsafe { incbif(a, b, self) }
    }

    fn factorial(self) -> f32 {
        unsafe { gammaf(self + 1.0) }
    }
    fn gamma(self) -> f32 {
        unsafe { gammaf(self) }
    }
    fn rgamma(self) -> f32 {
        unsafe { rgammaf(self) }
    }
    fn loggamma(self) -> f32 {
        unsafe { lgamf(self) }
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

    fn digamma(self) -> f32 {
        unsafe { psif(self) }
    }

    fn hyp1f1(self, a: f32, b: f32) -> f32 {
        unsafe { hypergf(a, b, self) }
    }
    fn hyp1f2(self, a: f32, b: f32, c: f32) -> f32 {
        let mut err = 0.0;
        unsafe { onef2f(a, b, c, self, &mut err) }
    }
    fn hyp2f1(self, a: f32, b: f32, c: f32) -> f32 {
        unsafe { hyp2f1f(a, b, c, self) }
    }
    fn hyp3f0(self, a: f32, b: f32, c: f32) -> f32 {
        let mut err = 0.0;
        unsafe { threef0f(a, b, c, self, &mut err) }
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
    fn bessely(self, v: f32) -> f32 {
        unsafe { yvf(v, self) }
    }
    fn besseli(self, v: f32) -> f32 {
        unsafe { ivf(v, self) }
    }
    fn besselk(self, n: i32) -> f32 {
        unsafe { knf(n, self) }
    }

    fn riemann_zeta(self) -> f32 {
        unsafe { 1.0 + zetacf(self) }
    }
    fn hurwitz_zeta(self, q: f32) -> f32 {
        unsafe { zetaf(self, q) }
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
            assert_almost_eq(5.0f64.beta(2.0), 0.03333333333333);
            assert_almost_eq(1.5f64.beta(2.0), 0.26666666666666);
            assert_almost_eq(1.5f64.logbeta(3.7), -2.1763500732197696);
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
        fn hypergeometric() {
            assert_almost_eq(1.5f64.hyp1f1(1.5, 3.0), 2.269381460919952778587441);
            assert_almost_eq(10.0f64.hyp1f2(1.5, 3.0, 2.25), 7.0792797035649206);
            assert_eq!(1.0f64.hyp2f1(-2.5, 3.5, 1.5), 0.0);
            assert_almost_eq(1.0f64.hyp2f1(-2.5, 3.0, 4.0), 0.06926406926406926406926407);
            assert_almost_eq(0.2f64.hyp3f0(-0.5, 0.6, -0.7), 1.04370133264);
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

    mod single {
        use super::super::FloatSpecial;
        use super::assert_almost_eq;

        #[test]
        fn beta() {
            assert_almost_eq(5.0f32.beta(2.0), 0.03333333333333);
            assert_almost_eq(1.5f32.beta(2.0), 0.26666666666666);
            assert_almost_eq(1.5f32.logbeta(3.7), -2.1763500732197696);
            assert_almost_eq(0.5f32.betainc(2.0, 3.0), 0.6875);
            assert_almost_eq(0.6875f32.betainc_inv(2.0, 3.0), 0.5);
        }

        #[test]
        fn factorial() {
            assert_almost_eq(0.5f32.factorial(), 0.886226925452758);
        }

        #[test]
        fn gamma() {
            assert_almost_eq(3.5f32.gamma(), 3.32335097044784);
            assert_eq!(1.0f32.rgamma(), 1.0);
            assert_almost_eq(4.0f32.rgamma(), 0.1666666666666666666);
            assert_almost_eq(13.2f32.loggamma(), 20.49400419456603678498394);
            assert_almost_eq(4.0f32.gammainc(2.0), 0.90842180555632912);
            assert_almost_eq(1.0 - 4.0f32.gammainc(2.0), 4.0f32.gammac(2.0));
            assert_almost_eq(4.0f32.gammac(2.0).gammac_inv(2.0), 4.0);
        }

        #[test]
        fn digamma() {
            assert_almost_eq(1.0f32.digamma(), -0.5772156649015328606065121);
        }

        #[test]
        fn hypergeometric() {
            assert_almost_eq(1.5f32.hyp1f1(1.5, 3.0), 2.269381460919952778587441);
            assert_almost_eq(10.0f32.hyp1f2(1.5, 3.0, 2.25), 7.0792797035649206);
            assert_eq!(1.0f32.hyp2f1(-2.5, 3.5, 1.5), 0.0);
            assert_almost_eq(1.0f32.hyp2f1(-2.5, 3.0, 4.0), 0.06926406926406926406926407);
            assert_almost_eq(0.2f32.hyp3f0(-0.5, 0.6, -0.7), 1.04370133264);
        }

        #[test]
        fn norm() {
            assert_almost_eq(2.0f32.norm(), 0.9772499);
            assert_almost_eq(0.9f32.norm_inv(), 1.281552);
        }

        #[test]
        fn besselj() {
            assert_almost_eq(10.0f32.besselj(2.0), 0.25463031368512062);
            assert_almost_eq(1000.0f32.besselj(2.0), -0.024777229528606);
            assert_almost_eq(0.75f32.besselj(4.0), 0.000801070086542314);
        }

        #[test]
        fn bessely() {
            assert_almost_eq(
                3.141592653589793f32.bessely(1.0), 0.3588729167767189594679827);
        }

        #[test]
        fn besseli() {
            assert_eq!(0.0f32.besseli(0.0), 1.0);
            assert_eq!(0.0f32.besseli(1.0), 0.0);
            assert_almost_eq(1.0f32.besseli(0.0), 1.266065877752008335598245);
        }

        #[test]
        fn besselk() {
            assert_almost_eq(1.0f32.besselk(0), 0.4210244382407083333356274);
            assert_almost_eq(100.0f32.besselk(0), 4.656628229175902018939005e-45);
        }

        #[test]
        fn riemann_zeta() {
            assert_almost_eq(2.0f32.riemann_zeta(), 1.64493406684822);
            assert_eq!(0.0f32.riemann_zeta(), -0.5);
            //assert_almost_eq(-1.0f32.riemann_zeta(), -0.0833333333333);
            assert_almost_eq(50.0f32.riemann_zeta(), 1.0);
        }

        #[test]
        fn hurwitz_zeta() {
            assert_almost_eq(2.0f32.hurwitz_zeta(3.0), 0.3949340668482);
            //assert_almost_eq(0.0f32.hurwitz_zeta(10.0), -9.5);
        }
    }
}
