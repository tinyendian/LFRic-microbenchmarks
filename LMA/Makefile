#!-----------------------------------------------------------------------------
#! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
#! For further details please refer to the file LICENCE.original which you
#! should have received as part of this distribution.
#!-----------------------------------------------------------------------------

# FC=gfortran
# DEBUGFLAGS = # -g
# FFLAGS=-O3 -fopenmp -march=skylake-avx512 -mtune=skylake-avx512 -mprefer-vector-width=512
# LDFLAGS=-O3 -fopenmp -march=skylake-avx512 -mtune=skylake-avx512 -mprefer-vector-width=512

FC=ifort
DEBUGFLAGS = -g -qopt-report -traceback -check bounds -warn all
FFLAGS=-O3 -xHost -qopt-zmm-usage=high -qopenmp -qopenmp-simd
LDFLAGS=-O3 -xHost -qopt-zmm-usage=high -qopenmp -qopenmp-simd

# FC=ftn
# DEBUGFLAGS = # -g
# FFLAGS=-O3 -hscalar3 -hvector3 -hcache3 -hnoacc -hpreferred_vector_width=512
# LDFLAGS=-O3 -hscalar3 -hvector3 -hcache3 -hnoacc -hpreferred_vector_width=512

# FC=nvfortran
# DEBUGFLAGS= # -g
# FFLAGS=-O3 -Mvect=simd:512 -tp=skylake -mp -target=multicore
# LDFLAGS=-O3 -Mvect=simd:512 -tp=skylake -mp -target=multicore

# FC=nvfortran
# DEBUGFLAGS=
# FFLAGS=-g -O2 -Minfo=accel -acc -gpu=cc60
# LDFLAGS=-g -O2 -Minfo=accel -acc -gpu=cc60

VPATH=../common:.

%: %.o
	$(FC) $(FFLAGS) $(DEBUGFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.F90 
	$(FC) $(FFLAGS) $(DEBUGFLAGS) -c $<

all : lma_driver

clean:
	rm -f *.s *.o *.mod *.optrpt *.lst lma_driver dino.saur.gz

tar:
	tar -cvf dino.saur Makefile *.F90
	gzip dino.saur

# dependencies
constants_mod.o:
dino_mod.o: constants_mod.o
compare_mod.o: constants_mod.o
argument_mod.o: 
kernel_mod.o:
colourist_mod.o: constants_mod.o
matrix_vector_kernel_mod.o: constants_mod.o argument_mod.o kernel_mod.o
compute_loop_mod.o : constants_mod.o matrix_vector_kernel_mod.o
lma_driver.o: constants_mod.o dino_mod.o compute_loop_mod.o compare_mod.o colourist_mod.o
lma_driver: lma_driver.o constants_mod.o argument_mod.o dino_mod.o kernel_mod.o compute_loop_mod.o matrix_vector_kernel_mod.o compare_mod.o colourist_mod.o

