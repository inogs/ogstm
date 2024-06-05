module load cmake
module load gcc
module load netcdf-fortran

export CC=gcc
export NETCDF_CFLAGS=$(nc-config --cflags)
export NETCDF_CLIBS=$(nc-config --libs)
export NETCDF_LIB=$(nc-config --libdir)
export NETCDF_INC=$(nc-config --includedir)

export FC=gfortran
export NETCDF_FFLAGS=$(nf-config --fflags)
export NETCDF_FLIBS=$(nf-config --flibs)
export NETCDFF_LIB=$(nf-config --prefix)/lib
export NETCDFF_INC=$(nf-config --includedir)
