    CC=mpicc
    FC=mpif90  $(OPENMP_FLAG)
    LD=mpif90  $(OPENMP_FLAG)
	AR=ar
	ARFLAGS=-r
	CPPFLAGS=
    DA_INC=/home/epascolo/TryCompile/3d_var
    DA_LIB=$(DA_INC)
# F Flags
        FFLAGS_DBG= -O0 -g -Wall -Wextra -fbounds-check -fimplicit-none -ffpe-trap=invalid,overflow -pedantic -ffixed-line-length-132
	FFLAGS= $(FFLAGS_DBG)  -I$(NETCDF_INC)  -I$(OPENMPI_HOME)/include   -I$(BFM_INC)  -I$(DA_INC) -fno-underscoring
	LDFLAGS= -L $(NETCDF_LIB) -lnetcdff -L $(MPI_LIB)/lib -lmpi_f90 -L$(BFM_LIB) -lbfm -L$(DA_LIB) -lvar_3d -L$(FE_LIB) -lf_exit -fno-underscoring
# C Flags
        CFLAGS_DBG=
        CFLAGS=  ${CFLAGS_DBG} -I$(NETCDF_INC) -I$(OPENMPI_HOME)/include
