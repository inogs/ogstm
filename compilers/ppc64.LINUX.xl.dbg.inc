#   LOC=/bgsys/drivers/ppcfloor/comm/xl.ndebug/bin    
    CC=mpicc
    FC=/bgsys/drivers/ppcfloor/comm/xl.ndebug/bin/mpixlf90_r
    LD=/bgsys/drivers/ppcfloor/comm/xl.ndebug/bin/mpixlf90_r 
    AR=ar
    ARFLAGS=-r
    CPPFLAGS=
    DA_LIB=$(DA_INC)
    FE_LIB=$(DA_INC)/libfexit
# F Flags
    FFLAGS_OPT=-O3 -qstrict -qarch=qp -qtune=qp -g -qdpc=e -qalias=noaryovrlp:nointptr -qalign=linuxppc -qsmp=noauto:omp
    FFLAGS_DBG=-O0 -g -qsmp=noauto:omp -qflttrap=nanq
    FFLAGS=$(FFLAGS_DBG) -I$(NETCDF_INC) -I$(BFM_INC) -qfixed=130
    LDFLAGS= -L$(NETCDF_LIB) -lnetcdff -lnetcdf -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -L$(SZIP_LIB) -lsz -L$(ZLIB_LIB) -lz -L$(BFM_LIB) -lbfm -L$(DA_LIB) -lvar_3d -L$(FE_LIB) -lf_exit -qsmp=noauto:omp
