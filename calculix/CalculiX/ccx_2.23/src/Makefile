# ============================================================
# CalculiX CCX 2.23 Makefile - Windows (MSYS2 gfortran)
# ============================================================

CFLAGS = -O2 -I ../SPOOLES.2.2 -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -Wno-implicit-int -Wno-return-mismatch -Wno-maybe-uninitialized -Wno-unused
FFLAGS = -O2 -cpp -fdefault-real-8 -fdefault-double-8 -std=legacy -Wno-unused -Wno-conversion

CC=gcc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $< -o $@
.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

include Makefile.inc

SCCXMAIN = ccx_2.23.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

SPOOLES_LIB = ../SPOOLES.2.2/spooles.a

LDFLAGS = -L/mingw64/lib -larpack -lopenblas -lpthread -lgfortran -lquadmath -lm

ccx_2.23.exe: $(OCCXMAIN) ccx_2.23.a $(SPOOLES_LIB)
	$(CC) $(CFLAGS) -c ccx_2.23.c
	$(FC) $(FFLAGS) -o $@ ccx_2.23.o ccx_2.23.a $(SPOOLES_LIB) $(LDFLAGS) -fopenmp

ccx_2.23.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

clean:
	rm -f *.o *.a ccx_2.23.exe
