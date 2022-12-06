# Include file for compilation with gfortran
# you may also need the following netcdf libraries 
# -L/usr/local/lib -lhdf5_hl -lhdf5 -lz

#  NCLIB= -L${NETCDFF_LIB} -lnetcdff -L${NETCDF_LIB} -lnetcdf
#
#  FC=mpif90
#  LD=mpif90
#
#  AR=ar
#  ARFLAGS=-r
#  
#  LDFLAGS= ${NCLIB}

FFLAGS=-cpp -Og -g3 -ggdb -Wall -Wextra -fbounds-check -fimplicit-none -pedantic -fsanitize=address -ffpe-trap=invalid,zero,overflow -I${NETCDF_INC} -I${NETCDFF_INC} --free-line-length-160
LDFLAGS=-L${NETCDFF_LIB} -lnetcdff -L${NETCDF_LIB} -lnetcdf
