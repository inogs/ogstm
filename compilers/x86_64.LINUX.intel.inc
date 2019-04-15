    CC=mpicc 
    FC=mpiifort $(OPENMP_FLAG) 
    LD=mpiifort $(OPENMP_FLAG)
    AR=ar
    ARFLAGS=-r
    CPPFLAGS=
    DA_LIB=$(DA_INC)
    FE_LIB=$(DA_INC)/libfexit 
# F Flags
    MPILIB= -L $(MPI_LIB)/lib -lmpi_f90    
    FFLAGS_OPT= -fno-math-errno -unroll=3 -qopt-subscript-in-range -align all  -cpp -heap-arrays
    FFLAGS= $(FFLAGS_OPT)  -I$(NETCDF_INC) -I$(NETCDFF_INC)  -I$(OPENMPI_HOME)/include -I$(BFM_INC) -I$(DA_INC) -extend_source
    LDFLAGS= -L $(NETCDF_LIB) -L $(NETCDFF_LIB)  -lnetcdff -lnetcdf -L $(BFM_LIB) -lbfm  -L$(DA_LIB) -lvar_3d  -L$(FE_LIB) -lf_exit
# C Flags
    CFLAGS_DBG=
    CFLAGS=  ${CFLAGS_DBG} -I$(NETCDF_INC) -I$(OPENMPI_HOME)/include
