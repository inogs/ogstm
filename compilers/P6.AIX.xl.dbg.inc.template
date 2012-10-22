	CC=/usr/bin/mpcc_r -c  
	FC=/usr/bin/mpxlf90_r -g  -qsmp=omp:noauto -cpp
	LD=/usr/bin/mpxlf90_r -g  -qsmp=omp:noauto
	AR=ar
	ARFLAGS=-r
        BFM_LIB= -L @@BFM_LIB@@ -lbfm  
	MPI_TRACE_LIB= -L/usr/lpp/ppe.hpct/lib64/ -lmpitrace 
# C Flags
	CPPFLAGS=
        CFLAGS_DBG=
	CFLAGS=  ${CFLAGS_DBG} -I$(NETCDF_INC)
#  Fortran Flags
        FFLAGS_DBG= -C -O0 -qflttrap=ov:zero:inv:en -qsigtrap=xl__trcedump -qinitauto=7FF7FFFF
	FFLAGS= ${FFLAGS_DBG} -I$(NETCDF_INC) -qsuppress=1501-245 -qcache=auto -qsource -g -c -qnolm  -q64 -qfixed=120 -qmaxmem=-1  -qdpc 
	LDFLAGS= ${MPI_TRACE_LIB} -L ${NETCDF_LIB} -lnetcdff ${BFM_LIB} ${FFLAGS_DBG} -qnolm -q64 -bmaxdata:0x400000000 -qfixed=120 -qsuppress=1501-245 -qcache=auto -qsource -g -qdpc 
