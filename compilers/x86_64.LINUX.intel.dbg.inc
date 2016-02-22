        CC=mpiicc
        FC=mpiifort $(OPENMP_FLAG)
        LD=mpiifort $(OPENMP_FLAG)
	DA_INC=/pico/scratch/userexternal/gbolzon0/TRANSITION/carbonatic_part/Code/3d_var
  DA_LIB=$(DA_INC)
  FE_LIB=$(DA_INC)/libfexit
  AR=ar
	ARFLAGS=-r
	CPPFLAGS=
        MPILIB= -L $(MPI_LIB)/lib -lmpi_f90
# F Flags
        FFLAGS_DBG=  -O0 -g -traceback -fp-stack-check -check bounds -fpe0
	FFLAGS= $(FFLAGS_DBG)  -I$(NETCDF_INC)  -I$(OPENMPI_HOME)/include -I$(BFM_INC) -I$(DA_INC) -extend_source
	LDFLAGS= -L $(NETCDF_LIB) -lnetcdff -lnetcdf -L$(BFM_LIB) -lbfm -L$(DA_LIB) -lvar_3d -L$(FE_LIB) -lf_exit
# C Flags
        CFLAGS_DBG=
        CFLAGS=  ${CFLAGS_DBG} -I$(NETCDF_INC) -I$(OPENMPI_HOME)/include
