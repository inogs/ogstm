module purge
module load profile/advanced
module load autoload

module load netcdf/4.7.3--spectrum_mpi--10.3.1--binary
module load netcdff/4.5.2--spectrum_mpi--10.3.1--binary
module load cmake
#module load gnu
#module load spectrum_mpi
module load numpy/1.19.4--python--3.8.2
export FC=gfortran
export CC=gcc 
#module load petsc/3.13.1--spectrum_mpi--10.3.1--binary
