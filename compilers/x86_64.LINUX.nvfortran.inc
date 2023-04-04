    CC=mpicc
    FC=mpif90  $(OPENMP_FLAG)
    LD=mpif90  $(OPENMP_FLAG)
	AR=ar
	ARFLAGS=-r
	CPPFLAGS=
# F Flags
#        FFLAGS_OPT= -O2  -fimplicit-none  -ffixed-line-length-132
        FFLAGS_OPT= -O2 -g -pg -acc=gpu -gpu=cc70
	FFLAGS= $(FFLAGS_OPT)  -I$(NETCDF_INC)  -I$(OPENMPI_HOME)/include  -I$(BFM_INC) -fno-underscoring
	LDFLAGS= -L $(NETCDF_LIB) -lnetcdff -L $(MPI_LIB)/lib -lmpi_f90 -L$(BFM_LIB) -lbfm -fno-underscoring
# C Flags
        CFLAGS_OPT= -g -pg
        CFLAGS=  ${CFLAGS_OPT} -I$(NETCDF_INC) -I$(OPENMPI_HOME)/include
