ifeq ("$(COMPILER)","gfortran")
  FC=gfortran
  DEBUGFLAGS = -g
  FFLAGS=-O3 -fopenmp -march=skylake-avx512 -mtune=skylake-avx512 -mprefer-vector-width=512
  LDFLAGS=-O3 -fopenmp -march=skylake-avx512 -mtune=skylake-avx512 -mprefer-vector-width=512
else ifeq ("$(COMPILER)","ifort")
  FC=ifort
  DEBUGFLAGS = -g -qopt-report -traceback -check bounds -warn all
  FFLAGS=-O3 -xHost -qopt-zmm-usage=high -qopenmp -qopenmp-simd
  LDFLAGS=-O3 -xHost -qopt-zmm-usage=high -qopenmp -qopenmp-simd
else ifeq ("$(COMPILER)","ftn")
  FC=ftn
  DEBUGFLAGS = -g
  FFLAGS=-O3 -hscalar3 -hvector3 -hcache3 -hnoacc -hpreferred_vector_width=512
  LDFLAGS=-O3 -hscalar3 -hvector3 -hcache3 -hnoacc -hpreferred_vector_width=512
else ifeq ("$(COMPILER)","nvfortran-skl")
  FC=nvfortran
  DEBUGFLAGS= -g
  FFLAGS=-O3 -Mvect=simd:512 -tp=skylake -mp -target=multicore
  LDFLAGS=-O3 -Mvect=simd:512 -tp=skylake -mp -target=multicore
else ifeq ("$(COMPILER)","nvfortran-p100-acc")
  FC=nvfortran
  DEBUGFLAGS= -g
  FFLAGS=-O3 -Minfo=accel -acc -gpu=cc60
  LDFLAGS=-O3 -Minfo=accel -acc -gpu=cc60
else ifeq ("$(COMPILER)","nvfortran-a100-acc")
  FC=nvfortran
  DEBUGFLAGS=-g
  FFLAGS=-O3 -march=x86-64-v3 -Minfo=accel -acc -gpu=cc80
  LDFLAGS=-O3 -march=x86-64-v3 -Minfo=accel -acc -gpu=cc80
else ifeq ("$(COMPILER)","nvfortran-a100-omp")
  FC=nvfortran
  DEBUGFLAGS= -g
  FFLAGS=-O3 -march=x86-64-v3 -Minfo=all -mp=gpu -gpu=cc80
  LDFLAGS=-O3 -march=x86-64-v3 -Minfo=all -mp=gpu -gpu=cc80
else ifneq ("$(MAKECMDGOALS)","clean")
  $(error Variable COMPILER not set, or unknown compiler)
endif

%: %.o
	$(FC) $(FFLAGS) $(DEBUGFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.F90
	$(FC) $(FFLAGS) $(DEBUGFLAGS) -c $<

all: mem_access_pattern

clean:
	rm -f *.s *.o *.mod *.optrpt *.lst *.bin *.fat *.ptx *.gpu *.cubin memtest

mem_access_pattern: mem_access_pattern.o
