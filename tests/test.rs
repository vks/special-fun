extern crate special_fun;
extern crate num_traits;

use special_fun::FloatSpecial;

use ::std::fmt::Debug;
use ::num_traits::{Float, FromPrimitive};

fn assert_almost_eq<T: Float + FromPrimitive + Debug>(a: T, b: T) {
    let tol: T = FromPrimitive::from_f32(1e-6).unwrap();
    if (a - b).abs() > tol {
        panic!("{:?} vs. {:?}", a, b);
    }
}

mod double {
    use super::FloatSpecial;
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
    fn error() {
        assert_eq!(0.0f64.erf(), 0.0);
        assert_almost_eq(1.0f64.erf(), 0.842700792949715);
        assert_almost_eq(-1.0f64.erf(), -0.842700792949715);
        assert_almost_eq(1.0f64.erfc(), 0.157299207050285);
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
    use super::FloatSpecial;
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
    fn error() {
        assert_eq!(0.0f64.erf(), 0.0);
        assert_almost_eq(1.0f64.erf(), 0.842700792949715);
        assert_almost_eq(-1.0f64.erf(), -0.842700792949715);
        assert_almost_eq(1.0f64.erfc(), 0.157299207050285);
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
