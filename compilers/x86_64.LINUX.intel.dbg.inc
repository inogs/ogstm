        CC=mpicc
        FC=mpif90 $(OPENMP_FLAG)
        LD=mpif90 $(OPENMP_FLAG) 
	AR=ar
	ARFLAGS=-r
	CPPFLAGS=
        MPILIB= -L $(MPI_LIB)/lib -lmpi_f90
# F Flags
        FFLAGS_DBG=  -O0 -g -traceback -fp-stack-check -check bounds -fpe0
	FFLAGS= $(FFLAGS_DBG)  -I$(NETCDF_INC)  -I$(OPENMPI_HOME)/include -I$(BFM_INC) -extend_source
	LDFLAGS= -L $(NETCDF_LIB) -lnetcdff -lnetcdf -L$(BFM_LIB) -lbfm 
# C Flags
        CFLAGS_DBG=
        CFLAGS=  ${CFLAGS_DBG} -I$(NETCDF_INC) -I$(OPENMPI_HOME)/include
