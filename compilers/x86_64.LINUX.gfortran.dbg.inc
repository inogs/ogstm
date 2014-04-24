    CC=mpicc
    FC=mpif90  $(OPENMP_FLAG)
    LD=mpif90  $(OPENMP_FLAG)
	AR=ar
	ARFLAGS=-r
	CPPFLAGS=
# F Flags
        FFLAGS_DBG= -O0 -g -Wall -Wextra -fbounds-check -fimplicit-none -pedantic -ffixed-line-length-132
	FFLAGS= $(FFLAGS_DBG)  -I$(NETCDF_INC)  -I$(OPENMPI_HOME)/include   -I$(BFM_INC) -fno-underscoring
	LDFLAGS= -L $(NETCDF_LIB) -lnetcdff -L $(MPI_LIB)/lib -lmpi_f90 -L$(BFM_LIB) -lbfm -fno-underscoring
# C Flags
        CFLAGS_DBG=
        CFLAGS=  ${CFLAGS_DBG} -I$(NETCDF_INC) -I$(OPENMPI_HOME)/include
