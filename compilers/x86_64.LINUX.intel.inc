    CC=mpicc 
    FC=mpif90 $(OPENMP_FLAG) 
    LD=mpif90 $(OPENMP_FLAG)
    AR=ar
    ARFLAGS=-r
    CPPFLAGS=
# F Flags
    MPILIB= -L $(MPI_LIB)/lib -lmpi_f90    
    FFLAGS_OPT= -xSSE4.2 -fno-math-errno -unroll=3 -opt-subscript-in-range -align all -axSSE4.2 -cpp -heap-arrays
    FFLAGS= $(FFLAGS_OPT)  -I$(NETCDF_INC)  -I$(OPENMPI_HOME)/include -I$(BFM_INC) -extend_source
    LDFLAGS= -L $(NETCDF_LIB) -lnetcdff -lnetcdf -L $(BFM_LIB) -lbfm
# C Flags
    CFLAGS_DBG=
    CFLAGS=  ${CFLAGS_DBG} -I$(NETCDF_INC) -I$(OPENMPI_HOME)/include
