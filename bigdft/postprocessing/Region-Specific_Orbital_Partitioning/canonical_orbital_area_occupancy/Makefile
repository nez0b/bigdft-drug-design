MKL=/aracbox/intel/compilers_and_libraries_2019.4.227/linux/mkl

BigDFT=/home2/r09222059/BigDFT_workspace/bigdft-suite/build_2025/install

# filename=toy_model_truncate_test
filename=toy_model_truncate_test
filename=toy_model_truncate_test_github
filename=my_localization
# filename=mat_dia_test
filename=plot_based_indices_input
filename=canonical_orbital_area_occupancy


INCLUDE=-I${BigDFT}/include 

LIBRARY=-L${BigDFT}/lib -L/home2/r09222059/.conda/envs/bigdft_new/lib -ldl -lrt -lstdc++ -lbigdft-1 -labinit -lxcf90 -lxc -lGaIn -lPSolver-1 -lCheSS-1 -lfutile-1 \
		-latlab-1 -lyaml -ldicts -lNTPoly -lorbs -lfmalloc-1 -llapack -lblas

FC=mpiifort
CC=gcc
FCFLAGS=-O2 -fopenmp -fbounds-check -fbacktrace -Wl,--no-as-needed

all: ${filename}.x
	${MAKE} ${filename}.x

${filename}.x: ${filename}.f90
	$(FC) -o $@ $(FCFLAGS) $(INCLUDE) ${filename}.f90 $(LIBRARY) 

.PHONY: clean
clean :
	rm -f ${filename}.x ${filename}.out ${filename}.log
