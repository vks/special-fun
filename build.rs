extern crate cc;

use std::path::PathBuf;

fn build(dir: &str, files: &str) {
    let mut build = cc::Build::new();
    build.include(dir);
    for file in files.split(' ') {
        let mut path = PathBuf::from(dir);
        path.push(file);
        build.file(path);
    }
    build.flag("-Wall");
    build.flag_if_supported("-fno-builtin");
    build.flag_if_supported("/Oi-");
    build.compile(dir);
}

fn main() {
    let dir_double = "cephes-double";
    let dir_single = "cephes-single";
    let files_double =
        "acosh.c airy.c asin.c asinh.c atan.c atanh.c bdtr.c beta.c \
         btdtr.c cbrt.c chbevl.c chdtr.c clog.c cmplx.c const.c \
         cosh.c dawsn.c drand.c ei.c ellie.c ellik.c ellpe.c ellpj.c ellpk.c \
         exp.c exp10.c exp2.c expn.c expx2.c fabs.c fac.c fdtr.c \
         fresnl.c gamma.c gdtr.c hyp2f1.c hyperg.c i0.c i1.c igami.c incbet.c \
         incbi.c igam.c isnan.c iv.c j0.c j1.c jn.c jv.c k0.c k1.c kn.c kolmogorov.c \
         log.c log2.c log10.c lrand.c nbdtr.c ndtr.c ndtri.c pdtr.c planck.c \
         polevl.c polmisc.c polylog.c polyn.c pow.c powi.c psi.c rgamma.c round.c \
         shichi.c sici.c sin.c sindg.c sinh.c spence.c stdtr.c struve.c \
         tan.c tandg.c tanh.c unity.c yn.c zeta.c zetac.c \
         sqrt.c floor.c setprec.c mtherr.c";
    let files_single =
        "acoshf.c airyf.c asinf.c asinhf.c atanf.c \
         atanhf.c bdtrf.c betaf.c cbrtf.c chbevlf.c chdtrf.c \
         clogf.c cmplxf.c constf.c coshf.c dawsnf.c ellief.c \
         ellikf.c ellpef.c ellpkf.c ellpjf.c expf.c exp2f.c \
         exp10f.c expnf.c expx2f.c facf.c fdtrf.c floorf.c fresnlf.c \
         gammaf.c gdtrf.c hypergf.c hyp2f1f.c igamf.c igamif.c \
         incbetf.c incbif.c i0f.c i1f.c ivf.c j0f.c j1f.c \
         jnf.c jvf.c k0f.c k1f.c knf.c logf.c log2f.c \
         log10f.c nbdtrf.c ndtrf.c ndtrif.c pdtrf.c polynf.c \
         powif.c powf.c psif.c rgammaf.c shichif.c sicif.c \
         sindgf.c sinf.c sinhf.c spencef.c sqrtf.c stdtrf.c \
         struvef.c tandgf.c tanf.c tanhf.c ynf.c zetaf.c \
         zetacf.c polevlf.c setprecf.c mtherrf.c";

    build(dir_double, files_double);
    build(dir_single, files_single);
}
