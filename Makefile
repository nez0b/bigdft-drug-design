
MKL=/share/apps/intel2021-base/mkl/2021.4.0

BigDFT=/mnt/c/Users/NTUQ-RA/Documents/pojen/Qchem/bigdft-suite/build/install

INCLUDE=-I${BigDFT}/include


LIBRARY=-L${BigDFT}/lib

FC=mpifort
CC=gcc 
FCFLAGS=-O2 -Wno-error -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow -fopenmp -m64 -g -Wl,--no-as-needed
#FC=gfortran
#FCFLAGS=-O2 -fopenmp

all: toy_model.x
	${MAKE} toy_model.x

toy_model.x: toy_model.f90
	$(FC) -o $@ $(FCFLAGS) $(INCLUDE) toy_model.f90 $(LIBRARY)

.PHONY: clean
clean :
	rm -f toy_model.x
