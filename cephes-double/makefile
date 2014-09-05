# Double precision Cephes library
# Makefile for unix or GCC

CC = gcc
#CC = /a/gnu/linux/gcc/xgcc -B/a/gnu/linux/gcc/ -I/a/gnu/linux/gcc/include
#CFLAGS = -g -O2 -Wall -fno-builtin # -ffloat-store
CFLAGS = -g -O2 -Wall -fno-builtin -fPIC
AR = ar
RANLIB = ranlib
INCS = mconf.h
AS = as

#ifndef $(OUT_DIR)
#	OUT_DIR = .
#endif

OBJS = $(OUT_DIR)/acosh.o $(OUT_DIR)/airy.o $(OUT_DIR)/asin.o $(OUT_DIR)/asinh.o $(OUT_DIR)/atan.o $(OUT_DIR)/atanh.o $(OUT_DIR)/bdtr.o $(OUT_DIR)/beta.o \
$(OUT_DIR)/btdtr.o $(OUT_DIR)/cbrt.o $(OUT_DIR)/chbevl.o $(OUT_DIR)/chdtr.o $(OUT_DIR)/clog.o $(OUT_DIR)/cmplx.o $(OUT_DIR)/const.o \
$(OUT_DIR)/cosh.o $(OUT_DIR)/dawsn.o $(OUT_DIR)/drand.o $(OUT_DIR)/ei.o $(OUT_DIR)/ellie.o $(OUT_DIR)/ellik.o $(OUT_DIR)/ellpe.o $(OUT_DIR)/ellpj.o $(OUT_DIR)/ellpk.o \
$(OUT_DIR)/exp.o $(OUT_DIR)/exp10.o $(OUT_DIR)/exp2.o $(OUT_DIR)/expn.o $(OUT_DIR)/expx2.o $(OUT_DIR)/fabs.o $(OUT_DIR)/fac.o $(OUT_DIR)/fdtr.o \
$(OUT_DIR)/fresnl.o $(OUT_DIR)/gamma.o $(OUT_DIR)/gdtr.o $(OUT_DIR)/hyp2f1.o $(OUT_DIR)/hyperg.o $(OUT_DIR)/i0.o $(OUT_DIR)/i1.o $(OUT_DIR)/igami.o $(OUT_DIR)/incbet.o \
$(OUT_DIR)/incbi.o $(OUT_DIR)/igam.o $(OUT_DIR)/isnan.o $(OUT_DIR)/iv.o $(OUT_DIR)/j0.o $(OUT_DIR)/j1.o $(OUT_DIR)/jn.o $(OUT_DIR)/jv.o $(OUT_DIR)/k0.o $(OUT_DIR)/k1.o $(OUT_DIR)/kn.o $(OUT_DIR)/kolmogorov.o \
$(OUT_DIR)/log.o $(OUT_DIR)/log2.o $(OUT_DIR)/log10.o $(OUT_DIR)/lrand.o $(OUT_DIR)/nbdtr.o $(OUT_DIR)/ndtr.o $(OUT_DIR)/ndtri.o $(OUT_DIR)/pdtr.o $(OUT_DIR)/planck.o \
$(OUT_DIR)/polevl.o $(OUT_DIR)/polmisc.o $(OUT_DIR)/polylog.o $(OUT_DIR)/polyn.o $(OUT_DIR)/pow.o $(OUT_DIR)/powi.o $(OUT_DIR)/psi.o $(OUT_DIR)/rgamma.o $(OUT_DIR)/round.o \
$(OUT_DIR)/shichi.o $(OUT_DIR)/sici.o $(OUT_DIR)/sin.o $(OUT_DIR)/sindg.o $(OUT_DIR)/sinh.o $(OUT_DIR)/spence.o $(OUT_DIR)/stdtr.o $(OUT_DIR)/struve.o \
$(OUT_DIR)/tan.o $(OUT_DIR)/tandg.o $(OUT_DIR)/tanh.o $(OUT_DIR)/unity.o $(OUT_DIR)/yn.o $(OUT_DIR)/zeta.o $(OUT_DIR)/zetac.o \
$(OUT_DIR)/sqrt.o $(OUT_DIR)/floor.o $(OUT_DIR)/setprec.o $(OUT_DIR)/mtherr.o

all: $(OUT_DIR)/libmd.a $(OUT_DIR)/mtst $(OUT_DIR)/dtestvec $(OUT_DIR)/monot $(OUT_DIR)/dcalc $(OUT_DIR)/paranoia # stamp-timing

lib: $(OUT_DIR)/libmd.a

stamp-timing: $(OUT_DIR)/libmd.a $(OUT_DIR)/mtst $(OUT_DIR)/time-it
	$(OUT_DIR)/time-it "$(OUT_DIR)/mtst > /dev/null"
	touch stamp-timing

time-it: $(OUT_DIR)/time-it.o
	$(CC) -o $(OUT_DIR)/time-it $(OUT_DIR)/time-it.o

$(OUT_DIR)/time-it.o: time-it.c
	$(CC) -O2 -c time-it.c -o $(OUT_DIR)/time-it.o

$(OUT_DIR)/dcalc: $(OUT_DIR)/dcalc.o $(OUT_DIR)/libmd.a
	$(CC) -o $(OUT_DIR)/dcalc $(OUT_DIR)/dcalc.o $(OUT_DIR)/libmd.a
#	aout2exe mtst

$(OUT_DIR)/mtst: $(OUT_DIR)/mtst.o $(OUT_DIR)/libmd.a
	$(CC) -v -o $(OUT_DIR)/mtst $(OUT_DIR)/mtst.o $(OUT_DIR)/libmd.a
#	gcc -Wl,-verbose -b i486-linuxaout -v -o mtst mtst.o libmd.a
#	coff2exe mtst

$(OUT_DIR)/mtst.o: mtst.c
	$(CC) -O2 -Wall -c mtst.c -o $(OUT_DIR)/mtst.o

$(OUT_DIR)/dtestvec: $(OUT_DIR)/dtestvec.o $(OUT_DIR)/libmd.a
	$(CC) -o $(OUT_DIR)/dtestvec $(OUT_DIR)/dtestvec.o $(OUT_DIR)/libmd.a

$(OUT_DIR)/dtestvec.o: dtestvec.c
	$(CC) -g -c dtestvec.c -o $(OUT_DIR)/dtestvec.o

$(OUT_DIR)/monot: $(OUT_DIR)/monot.o $(OUT_DIR)/libmd.a
	$(CC) -o $(OUT_DIR)/monot $(OUT_DIR)/monot.o $(OUT_DIR)/libmd.a

$(OUT_DIR)/monot.o: monot.c
	$(CC) -g -c monot.c -o $(OUT_DIR)/monot.o

$(OUT_DIR)/paranoia: $(OUT_DIR)/paranoia.o $(OUT_DIR)/setprec.o $(OUT_DIR)/libmd.a
	$(CC) -o $(OUT_DIR)/paranoia $(OUT_DIR)/paranoia.o $(OUT_DIR)/setprec.o $(OUT_DIR)/libmd.a

$(OUT_DIR)/paranoia.o: paranoia.c
	$(CC) -Wno-implicit -c paranoia.c -o $(OUT_DIR)/paranoia.o

$(OUT_DIR)/libmd.a: $(OBJS) $(INCS)
# for real Unix:
	$(AR) crus $(OUT_DIR)/libmd.a $(OBJS)
# for djgcc MSDOS:
#	>libmd.rf -rv libmd.a $(OBJS)
#	$(AR) @libmd.rf
	$(RANLIB) $(OUT_DIR)/libmd.a

# If the following are all commented out, the C versions
# will be used by default.

# IBM PC:
#sqrt.o: sqrt.387
#	$(AS) -o sqrt.o sqrt.387
#
#floor.o: floor.387
#	$(AS) -o floor.o floor.387
#
#setprec.o: setprec.387
#	$(AS) -o setprec.o setprec.387

# ELF versions for linux (no underscores)
#sqrt.o: sqrtelf.387
#	$(AS) -o sqrt.o sqrtelf.387
#
#floor.o: floorelf.387
#	$(AS) -o floor.o floorelf.387
#
#setprec.o: setprelf.387
#	$(AS) -o setprec.o setprelf.387

# Motorola 68881. Caution, subroutine return conventions vary.
#sqrt.o: sqrt.688
#	$(AS) -o sqrt.o sqrt.688
#
#setprec.o: setprec.688
#	$(AS) -o setprec.o setprec.688

# SPARC:
#sqrt.o: sqrt.spa
#	$(AS) -o sqrt.o sqrt.spa

clean:
	rm -f $(OUT_DIR)/*.o
	rm -f $(OUT_DIR)/mtst
	rm -f $(OUT_DIR)/paranoia
	rm -f $(OUT_DIR)/dcalc
	rm -f $(OUT_DIR)/dtestvec
	rm -f $(OUT_DIR)/monot
	rm -f $(OUT_DIR)/libmd.a
	rm -f $(OUT_DIR)/time-it
	rm -f $(OUT_DIR)/dtestvec

