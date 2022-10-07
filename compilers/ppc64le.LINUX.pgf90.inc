FFLAGS=-fast -O3 -acc=strict -ta=tesla:cc70,host -I$(NETCDFF_INC)
CFLAGS=-fast -O3 -I$(NETCDF_INC)
LDFLAGS=-lnetcdf -lnetcdff -lmpi_f90 -L$(BFM_LIB) -lbfm -fno-underscoring

