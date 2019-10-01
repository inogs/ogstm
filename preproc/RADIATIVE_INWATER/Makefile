# Generic Makefile to compile c-postproc utils
# for Code_Saturne
#
# Arnau Miro - UPC (ESEIAAT)

# Configuration
#
OGSTM_DIR = ../../../ogstm
BFM_DIR   = ../../../bfm
OGSTM_BLD = ../../../OGSTM_BUILD
BUILDDIR  = build
OPTRTDIR  = optrpt
COMPILERS = x86_64.LINUX.intel

# Environment
#
OGSTM_SRC = $(OGSTM_DIR)/src
OGSTM_CMP = $(OGSTM_DIR)/compilers
BFM_INC   = ${BFM_DIR}/include
BFM_LIB   = ${BFM_DIR}/lib

# Sources and Objects
#
SOURCES := $(shell find $(OGSTM_SRC)  -name '*.f90' -or  -name '*.F90' -or -name '*.F' )
SRC_OUT := $(shell find "$(OGSTM_SRC)/BC/src/"  -name '*.f90' -or  -name '*.F90' -or -name '*.F')
SRC_OUT += $(shell find $(OGSTM_SRC)  -name '*.f03' -or -name '*.pf' -or -name '*.backup' )
SOURCES := $(filter-out $(SRC_OUT),$(SOURCES))
# Modules
SRC_MOD := $(shell grep -Ril "end module" $(OGSTM_SRC)/*)
SRC_MOD := $(filter-out $(SRC_OUT),$(SRC_MOD))
SOURCES := $(filter-out $(SRC_MOD),$(SOURCES))
# Convert .f90, .F90 and .F to .o
OBJECTS := $(SOURCES:.f90=.o)
OBJECTS := $(OBJECTS:.F90=.o)
OBJECTS := $(OBJECTS:.F=.o)
# Change to build dir
OBJECTS := $(foreach obj,$(OBJECTS),$(patsubst $(shell dirname $(obj))/%,$(BUILDDIR)/%,$(obj)))
# By hand just because it won't pick the dependencies -.-
MODULES := build/stringop.o build/calendar.o build/bc_data.o build/parameter.o build/bc.o build/Timers.o build/memory.o build/bc_aux.o build/hard_open.o build/rivers.o build/sponge.o build/nudging.o build/TimeManager.o build/bc_handle.o build/bc_set.o build/closed.o build/OPT_mem.o build/SLAE.o build/BIO_mem.o build/FN_mem.o build/DIA_mem.o build/SED_mem.o build/get_mem_mod.o build/ADV_mem.o build/BC_mem.o build/IO_mem.o build/bc_atm.o build/bc_co2.o build/ogstm_mpi.o build/step.o build/ZDF_mem.o build/ogstm.o

# Compilers
#
include $(OGSTM_CMP)/$(COMPILERS).inc

# Default
#
# One rule to compile them all, one rule to find them,
# One rule to bring them all and in the compiler link them.
default: dirs libs apps

# Dirs
#
.PHONY: dirs
dirs:
	@mkdir -p $(BUILDDIR) $(OPTRTDIR)

# Applications
#
# The first application links with the library created during the compilation
# of OGSTM-BFM. The second links with the library created by this Makefile
.PHONY: apps
apps: compute.xx compute2.xx

compute.xx: libogstm_lib.a compute.F90
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS) -L./
compute2.xx: libogstm.a compute.F90
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS) -L./

# Libraries
#
.PHONY:libs
libs: libogstm_lib.a libogstm.a

libogstm_lib.a: $(OGSTM_BLD)/libogstm_lib.a
	@ln -s $< $@
libogstm.a: dirs modules objects
	ar qc $@ $(MODULES) $(OBJECTS)
	ranlib $@

# Fortran modules
#
.PHONY: modules
modules: $(MODULES)

build/stringop.o: $(OGSTM_SRC)/General/stringop.F
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/calendar.o: $(OGSTM_SRC)/General/calendar.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_data.o: $(OGSTM_SRC)/BC/bc_data.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/parameter.o: $(OGSTM_SRC)/General/parameter.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)	
build/bc.o: $(OGSTM_SRC)/BC/bc.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/Timers.o: $(OGSTM_SRC)/General/Timers.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/memory.o: $(OGSTM_SRC)/General/memory.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_aux.o: $(OGSTM_SRC)/BC/bc_aux.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/hard_open.o: $(OGSTM_SRC)/BC/hard_open.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/rivers.o: $(OGSTM_SRC)/BC/rivers.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/sponge.o: $(OGSTM_SRC)/BC/sponge.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/nudging.o: $(OGSTM_SRC)/BC/nudging.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/TimeManager.o: $(OGSTM_SRC)/General/TimeManager.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_handle.o: $(OGSTM_SRC)/BC/bc_handle.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_set.o: $(OGSTM_SRC)/BC/bc_set.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/closed.o: $(OGSTM_SRC)/BC/closed.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/OPT_mem.o: $(OGSTM_SRC)/BIO-OPTICS/OPT_mem.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/SLAE.o: $(OGSTM_SRC)/BIO-OPTICS/SLAE.F90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/BIO_mem.o: $(OGSTM_SRC)/BIO/BIO_mem.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/FN_mem.o: $(OGSTM_SRC)/BIO/FN_mem.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/DIA_mem.o: $(OGSTM_SRC)/IO/DIA_mem.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/SED_mem.o: $(OGSTM_SRC)/BIO/SED_mem.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/get_mem_mod.o: $(OGSTM_SRC)/General/get_mem_mod.F90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/ADV_mem.o: $(OGSTM_SRC)/PHYS/ADV_mem.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/BC_mem.o: $(OGSTM_SRC)/IO/BC_mem.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/IO_mem.o: $(OGSTM_SRC)/IO/IO_mem.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_atm.o: $(OGSTM_SRC)/IO/bc_atm.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_co2.o: $(OGSTM_SRC)/IO/bc_co2.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/ogstm_mpi.o: $(OGSTM_SRC)/MPI/ogstm_mpi.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/step.o: $(OGSTM_SRC)/General/step.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/ZDF_mem.o: $(OGSTM_SRC)/PHYS/ZDF_mem.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)
build/ogstm.o: $(OGSTM_SRC)/General/ogstm.f90
	$(FC) $(FFLAGS) -I build/ -c $<
	@mv *.o *.mod $(BUILDDIR)
	@mv *.optrpt $(OPTRTDIR)

# Objects
#
.PHONY: objects
objects: $(OBJECTS)

$(OBJECTS): $(SOURCES)
	$(FC) $(FFLAGS) -I build/ -c -o $@ $<

# Clean
#
.PHONY: clean
clean: 
	@rm -rf *.o *.optrpt $(BUILDDIR) $(OPTRTDIR)

.PHONY: clean_libs
clean_libs:
	@rm -f *.a

.PHONY: clean_apps
clean_apps: 
	@rm -rf *.xx

.PHONY: cleanall
cleanall: clean clean_libs clean_apps

