# How to compile BigDFT-suite

The following procedure has been tested on ubuntu machines, and ubuntu-WSL2

## Prerequisites

Download BigDFT source code:
```bash
git clone https://gitlab.com/l_sim/bigdft-suite
git checkout 1.9.4 # switch to stable version
```

Dependencies:
```bash
sudo apt-get install gfortran cmake autoconf pkg-config #fortran compiler, cmake
sudo apt-get install libblas-dev liblapack-dev #blas and lapack
sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev # open-mpi
```

## Install

```bash
cd $(BigDFT-src-location)
mkdir build
cd build
#python ../Installer.py autogen
python ../Installer.py build -f ../rcfiles/ubuntu_MPI.rc
```


## GPU Version
Enable Intel MKL and MPI Libraries
```bash
source /aracbox/intel/bin/compilervars.sh -arch intel64 -platform linux
source /aracbox/intel/compilers_and_libraries_2019.3.199/linux/mpi/intel64/bin/mpivars.sh
# check whether mpiifort and $MKLROOT exist
which mpiifort
echo $MKLROOT
```
Compile CUDA version BigDFT
```bash
cd $(BigDFT-src-location)
mkdir build-cuda
cd build-cuda
python ../Installer.py build -f ../rcfiles/spring-cuda.rc
```
Modify `input.yaml` to enable CUDA in Poisson Solver:
```yaml
+ psolver:
+   setup:
+     accel: CUDA
```

And in `toy_model.f90`, 
```fortran
! line 180
- use dictionaries  
+ use dictionaries, dict_set => set
```
And modify the setting before initialize the pkernel:
```fortran
! line 550
+ call dict_set(dict//'setup'//'accel', 'CUDA')

pkernel=pkernel_init(iproc,nproc,dict,&
      dom,(/Lzd%Glr%nboxi(2,1), Lzd%Glr%nboxi(2,2), Lzd%Glr%nboxi(2,3)/),&
       (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/))
```


