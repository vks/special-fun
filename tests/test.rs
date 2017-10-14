extern crate special_fun;
extern crate num_traits;

use ::std::fmt::Debug;
use ::num_traits::{Float, FromPrimitive};

fn assert_almost_eq<T: Float + FromPrimitive + Debug>(a: T, b: T) {
    let tol: T = FromPrimitive::from_f32(1e-6).unwrap();
    if (a - b).abs() > tol {
        panic!("{:?} vs. {:?}", a, b);
    }
}

mod cephes_double {
    use special_fun::cephes_double;
    use super::assert_almost_eq;
    use std::f64::consts::PI;
    #[test]
    fn exp10() {
        assert_almost_eq(cephes_double::exp10(0.0f64), 1.0);
        assert_almost_eq(cephes_double::exp10(1.0f64), 10.0);
        assert_almost_eq(cephes_double::exp10(1.2f64), 15.848931924611133);
        assert_almost_eq(cephes_double::exp10(1.5f64), 31.622776601683793);
        assert_almost_eq(cephes_double::exp10(308.2547155599167f64), 1.7976931348620926e308);
    }

    #[test]
    fn expm1() {
        assert_almost_eq(cephes_double::expm1(0.0f64), 0.0);
        assert_almost_eq(cephes_double::expm1(0.001f64), 0.0010005001667083846);
        assert_almost_eq(cephes_double::expm1(0.1f64), 0.10517091807564771);
        assert_almost_eq(cephes_double::expm1(1.0f64), 1.718281828459045);
    }

    #[test]
    fn expx2(){
        assert_almost_eq(cephes_double::expx2(2.0f64, 1), 54.598150033144236);
        assert_almost_eq(cephes_double::expx2(2.0f64, 2), 54.598150033144236);
        assert_almost_eq(cephes_double::expx2(2.0f64, -1), 1.0/54.598150033144236);
        assert_almost_eq(cephes_double::expx2(2.0f64, -2), 1.0/54.598150033144236);
    }

    #[test]
    fn ei(){
        assert_almost_eq(cephes_double::ei(1.0f64), 1.89511781);
        assert_almost_eq(cephes_double::ei(0.0f64), 0.0);
        assert_almost_eq(cephes_double::ei(0.6f64), 0.7698812899373594);
    }

    #[test]
    fn erf() {
        assert_almost_eq(cephes_double::erf(0.0f64), 0.0);
        assert_almost_eq(cephes_double::erf(0.5f64), 0.52049987781304654);
        assert_almost_eq(cephes_double::erf(-0.5f64), -0.52049987781304654);
        assert_almost_eq(cephes_double::erf(1.0f64), 0.84270079294971487);
        assert_almost_eq(cephes_double::erf(-1.0f64), -0.84270079294971487);
        assert_almost_eq(cephes_double::erf(10.0f64), 0.99999999999999999);
        assert_almost_eq(cephes_double::erf(-10.0f64), -0.99999999999999999);
    }

    #[test]
    fn erfc() {
        assert_almost_eq(cephes_double::erfc(0.0f64), 1.0 - 0.0);
        assert_almost_eq(cephes_double::erfc(0.5f64), 1.0 - 0.52049987781304654);
        assert_almost_eq(cephes_double::erfc(-0.5f64), 1.0 + 0.52049987781304654);
        assert_almost_eq(cephes_double::erfc(1.0f64), 1.0 - 0.84270079294971487);
        assert_almost_eq(cephes_double::erfc(-1.0f64), 1.0 + 0.84270079294971487);
        assert_almost_eq(cephes_double::erfc(10.0f64), 1.0 - 0.99999999999999999);
        assert_almost_eq(cephes_double::erfc(-10.0f64), 1.0 + 0.99999999999999999);
    }

    #[test]
    fn log1p() {
        assert_almost_eq(cephes_double::log1p(0.0f64), 0.0);
        assert_almost_eq(cephes_double::log1p(1.0e-5_f64), 0.000009999950000333332);
        assert_almost_eq(cephes_double::log1p(1.0f64), 0.6931471805599453);
    }

    #[test]
    fn spence() {
        assert_almost_eq(cephes_double::spence(0.0f64), PI*PI/6.0f64);
        assert_almost_eq(cephes_double::spence(1.0f64), 0.0);
        assert_almost_eq(cephes_double::spence(2.0f64), -PI*PI/12.0f64);
    }

    #[test]
    fn cosm1() {
        assert_almost_eq(cephes_double::cosm1(0.0f64), 0.0);
        assert_almost_eq(cephes_double::cosm1(1.0e-5_f64), -0.00000000004999999999958334);
        assert_almost_eq(cephes_double::cosm1(0.1f64), -0.004995834721974235);
    }

    #[test]
    fn sici() {
        let mut si = 0.0_f64;
        let mut ci = 0.0_f64;

        cephes_double::sici(0.5_f64, &mut si, &mut ci);
        assert_almost_eq(si, 0.49310741804306674);
        assert_almost_eq(ci, -0.17778407880661287);

        cephes_double::sici(0.1_f64, &mut si, &mut ci);
        assert_almost_eq(si, 0.09994446110827694);
        assert_almost_eq(ci, -1.7278683866572966);

        cephes_double::sici(0.0_f64, &mut si, &mut ci);
        assert_almost_eq(si, 0.0);
    }

    #[test]
    fn shichi() {
        let mut shi = 0.0_f64;
        let mut chi = 0.0_f64;

        cephes_double::shichi(0.5_f64, &mut shi, &mut chi);
        assert_almost_eq(shi, 0.5069967498196671);
        assert_almost_eq(chi, -0.05277684495649357);

        cephes_double::shichi(0.1_f64, &mut shi, &mut chi);
        assert_almost_eq(shi, 0.10005557222505701);
        assert_almost_eq(chi, -1.7228683861943335);

        cephes_double::shichi(0.0_f64, &mut shi, &mut chi);
        assert_almost_eq(shi, 0.0);
    }

    #[test]
    fn beta() {
        assert_almost_eq(cephes_double::beta(5.0f64, 2.0f64), 0.03333333333333);
        assert_almost_eq(cephes_double::beta(1.5f64, 2.0f64), 0.26666666666666);
        assert_almost_eq(cephes_double::beta(1.0f64, 1.0f64), 1.0);
    }

    #[test]
    fn incbet() {
        assert_almost_eq(cephes_double::incbet(5.0f64, 2.0f64, 1.0f64), 1.0);
        assert_almost_eq(cephes_double::incbet(5.0f64, 2.0f64, 0.5f64), 0.109375f64);
        assert_almost_eq(cephes_double::incbet(1.0f64, 2.0f64, 0.2f64), 0.36f64);
    }

    #[test]
    fn incbi() {
        assert_almost_eq(cephes_double::incbi(5.0f64, 2.0f64, 1.0f64), 1.0);
        assert_almost_eq(cephes_double::incbi(5.0f64, 2.0f64, 0.5f64), 0.73555001670434);
        assert_almost_eq(cephes_double::incbi(1.0f64, 2.0f64, 0.2f64), 0.10557280900008412);
        assert_almost_eq(cephes_double::incbi(2.0, 3.0, 0.6875f64), 0.5);
    }

    #[test]
    fn gamma() {
        assert_almost_eq(cephes_double::gamma(0.0f64), 1.0);
        assert_almost_eq(cephes_double::gamma(0.5f64), PI.sqrt());
        assert_almost_eq(cephes_double::gamma(1.0f64), 1.0);
        assert_almost_eq(cephes_double::gamma(1.5f64), 0.5*PI.sqrt());
        assert_almost_eq(cephes_double::gamma(2.0f64), 1.0);
        assert_almost_eq(cephes_double::gamma(2.5f64), 0.75*PI.sqrt());
        assert_almost_eq(cephes_double::gamma(3.0f64), 2.0);
        assert_almost_eq(cephes_double::gamma(3.5f64), 15.0f64/8.0f64*PI.sqrt());
        assert_almost_eq(cephes_double::gamma(4.0f64), 6.0);
    }

    #[test]
    fn rgamma() {
        assert_almost_eq(cephes_double::rgamma(0.0f64), 0.0);
        assert_almost_eq(cephes_double::rgamma(0.5f64), 1.0/PI.sqrt());
        assert_almost_eq(cephes_double::rgamma(1.0f64), 1.0);
        assert_almost_eq(cephes_double::rgamma(1.5f64), 1.0/(0.5*PI.sqrt()));
        assert_almost_eq(cephes_double::rgamma(2.0f64), 1.0);
        assert_almost_eq(cephes_double::rgamma(2.5f64), 1.0/(0.75*PI.sqrt()));
        assert_almost_eq(cephes_double::rgamma(3.0f64), 1.0/2.0);
        assert_almost_eq(cephes_double::rgamma(3.5f64), 1.0/(15.0f64/8.0f64*PI.sqrt()));
        assert_almost_eq(cephes_double::rgamma(4.0f64), 1.0/6.0);
    }

    #[test]
    fn lgam() {
        assert_almost_eq(cephes_double::lgam(0.5f64), PI.sqrt().ln());
        assert_almost_eq(cephes_double::lgam(1.0f64), 0.0);
        assert_almost_eq(cephes_double::lgam(1.5f64), (0.5*PI.sqrt()).ln());
        assert_almost_eq(cephes_double::lgam(2.0f64), 0.0);
        assert_almost_eq(cephes_double::lgam(2.5f64), (0.75f64*PI.sqrt()).ln());
        assert_almost_eq(cephes_double::lgam(3.0f64), 2.0f64.ln());
        assert_almost_eq(cephes_double::lgam(3.5f64), (15.0f64/8.0f64*PI.sqrt()).ln());
        assert_almost_eq(cephes_double::lgam(4.0f64), 6.0f64.ln());
    }

    #[test]
    fn igam() {
        assert_almost_eq(cephes_double::igam(1.0f64, 1.0f64), 1.0f64 - (-1.0f64).exp());
        assert_almost_eq(cephes_double::igam(1.0f64, 1.5f64), 1.0f64 - (-1.5f64).exp());
        assert_almost_eq(cephes_double::igam(1.0f64, 2.0f64), 1.0f64 - (-2.0f64).exp());
        assert_almost_eq(cephes_double::igam(0.5f64, 1.0f64),
                         PI.sqrt() * cephes_double::erf(1.0f64.sqrt()) / cephes_double::gamma(0.5f64));
        assert_almost_eq(cephes_double::igam(0.5f64, 1.5f64),
                         PI.sqrt() * cephes_double::erf(1.5f64.sqrt()) / cephes_double::gamma(0.5f64));
        assert_almost_eq(cephes_double::igam(0.5f64, 2.0f64),
                         PI.sqrt() * cephes_double::erf(2.0f64.sqrt()) / cephes_double::gamma(0.5f64));
    }

    #[test]
    fn igamc() {
        assert_almost_eq(cephes_double::igamc(1.0f64, 1.0f64),  1.0 - (1.0f64 - (-1.0f64).exp()));
        assert_almost_eq(cephes_double::igamc(1.0f64, 1.5f64), 1.0 - (1.0f64 - (-1.5f64).exp()));
        assert_almost_eq(cephes_double::igamc(1.0f64, 2.0f64), 1.0 - (1.0f64 - (-2.0f64).exp()));
        assert_almost_eq(cephes_double::igamc(0.5f64, 1.0f64),
                         1.0 - (PI.sqrt() * cephes_double::erf(1.0f64.sqrt()) / cephes_double::gamma(0.5f64)));
        assert_almost_eq(cephes_double::igamc(0.5f64, 1.5f64),
                         1.0 - (PI.sqrt() * cephes_double::erf(1.5f64.sqrt()) / cephes_double::gamma(0.5f64)));
        assert_almost_eq(cephes_double::igamc(0.5f64, 2.0f64),
                         1.0 - (PI.sqrt() * cephes_double::erf(2.0f64.sqrt()) / cephes_double::gamma(0.5f64)));
    }

    #[test]
    fn igami() {
        assert_almost_eq(cephes_double::igami(0.5f64, 0.1f64), 1.3527717270477);
        assert_almost_eq(cephes_double::igami(0.5f64, 0.2f64), 0.8211872075749);
        assert_almost_eq(cephes_double::igami(0.2f64, 0.1f64), 0.6049023209866);
    }

    #[test]
    fn psi() {
        let euler = 0.5772156649015329;
        assert_almost_eq(cephes_double::psi(1.0f64), -euler);
        assert_almost_eq(cephes_double::psi(0.5f64), -2.0*2.0f64.ln() - euler);
        assert_almost_eq(cephes_double::psi(20.0f64), 2.970523992242149);
    }

    #[test]
    fn fac() {
        assert_almost_eq(cephes_double::fac(0i32), 1.0);
        assert_almost_eq(cephes_double::fac(1i32), 1.0);
        assert_almost_eq(cephes_double::fac(2i32), 2.0);
        assert_almost_eq(cephes_double::fac(3i32), 6.0);
        assert_almost_eq(cephes_double::fac(4i32), 24.0);
        assert_almost_eq(cephes_double::fac(5i32), 120.0);
        assert_almost_eq(cephes_double::fac(6i32), 720.0);
        assert_almost_eq(cephes_double::fac(20i32), 2432902008176640000.0);
    }
}

mod double {
    use special_fun::FloatSpecial;
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
    use special_fun::FloatSpecial;
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
