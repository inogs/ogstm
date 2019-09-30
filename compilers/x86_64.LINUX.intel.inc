    CC=mpiicc 
    FC=mpiifort $(OPENMP_FLAG) 
    LD=mpiifort $(OPENMP_FLAG)
    AR=ar
    ARFLAGS=-r
    CPPFLAGS=
# F Flags
    MPILIB= -L $(MPI_LIB)/lib -lmpi_f90    
    FFLAGS_KEY= -Dkey_INCLUDE_BFM_PELCO2 -Dkey_kef -Dkey_mpp -Dkey_mpp_mpi -Dkey_trahdfbilap -Dkey_trahdfcoef1d -Dkey_trc_bfm -Dkey_trc_dmp -Dkey_trc_hdfbilap -Dkey_trc_sed -Dkey_trc_smolar 
    FFLAGS_OPT= -std=f2003 -fno-math-errno -O2 -xAVX -qopt-report5 -g -cpp -align array64byte
    FFLAGS= $(FFLAGS_KEY) $(FFLAGS_OPT)  -I$(NETCDF_INC) -I$(NETCDFF_INC)  -I$(OPENMPI_HOME)/include -I$(BFM_INC) -extend_source
    LDFLAGS= -L $(NETCDF_LIB) -lnetcdff -L $(NETCDFF_LIB) -lnetcdf -L $(BFM_LIB) -lbfm
# C Flags
    CFLAGS_DBG=
    CFLAGS=  ${CFLAGS_DBG} -I$(NETCDF_INC) -I$(OPENMPI_HOME)/include
