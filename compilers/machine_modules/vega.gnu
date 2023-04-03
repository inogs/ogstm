module load netCDF-Fortran/4.6.0-gompi-2022a
module load CMake/3.23.1-GCCcore-11.3.0
module load MKL

# modules define these $.._HOME env variable
export NETCDF_INC=$(nc-config --includedir)
export NETCDF_LIB=$(nc-config --libdir)
export NETCDFF_INC=$EBROOTNETCDFMINFORTRAN/include
export NETCDFF_LIB=$EBROOTNETCDFMINFORTRAN/lib

# debug
#echo NETCDF_INC=$NETCDF_INC
#echo NETCDF_LIB=$NETCDF_LIB
#echo NETCDFF_INC=$NETCDFF_INC
#echo NETCDFF_LIB=$NETCDFF_LIB

