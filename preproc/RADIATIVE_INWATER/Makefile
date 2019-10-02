# Generic Makefile to compile c-postproc utils
# for Code_Saturne
#
# Arnau Miro - UPC (ESEIAAT)

# Configuration
#
OGSTM_DIR = ../../../ogstm
BFM_DIR   = ../../../bfm
OGSTM_BLD = ../../../OGSTM_BUILD
FOLDERS   = BIO BIO-OPTICS General IO MPI namelists PHYS BC
BUILDDIR  = build
INCDIR    = include
OPTRTDIR  = optrpt
COMPILERS = x86_64.LINUX.intel

# Environment
#
OGSTM_SRC = $(OGSTM_DIR)/src
OGSTM_CMP = $(OGSTM_DIR)/compilers
OGSTM_INC = $(OGSTM_BLD)
BFM_INC   = ${BFM_DIR}/include
BFM_LIB   = ${BFM_DIR}/lib

# Sources and Objects
#
# Provide the Makefile a search path for all the sources
VPATH = $(shell echo $(addprefix $(OGSTM_SRC)/,$(FOLDERS)) | sed "s/ \+/:/g")
# Actively look for the sources to compile
SOURCES := $(foreach folder,$(FOLDERS),$(shell find $(OGSTM_SRC)/$(folder)  -name '*.f90' -or  -name '*.F90' -or -name '*.F' -or -name '*.c'))
SRC_OUT := $(shell find "$(OGSTM_SRC)/BC/src/"  -name '*.f90' -or  -name '*.F90' -or -name '*.F')
SRC_OUT += $(shell find $(OGSTM_SRC)  -name '*.f03' -or -name '*.pf' -or -name '*.backup' )
SOURCES := $(filter-out $(SRC_OUT),$(SOURCES))
# Actively look for the modules and take them out from the sources
SRC_MOD := $(shell grep -Ril "end module" $(OGSTM_SRC)/*)
SRC_MOD := $(filter-out $(SRC_OUT),$(SRC_MOD))
SOURCES := $(filter-out $(SRC_MOD),$(SOURCES))
# Generate objects
OBJECTS := $(SOURCES:.c=.o)
OBJECTS := $(OBJECTS:.f90=.o)
OBJECTS := $(OBJECTS:.F90=.o)
OBJECTS := $(OBJECTS:.F=.o)
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
default: dirs modules objects libs apps

# Dirs
#
.PHONY: dirs
dirs:
	@mkdir -p $(BUILDDIR) $(INCDIR) $(OPTRTDIR)

# Applications
#
# The first application links with the library created during the compilation
# of OGSTM-BFM. The second links with the library created by this Makefile
.PHONY: apps
apps: compute.xx compute2.xx

compute.xx: compute.F90 libogstm_lib.a
	$(FC) $(FFLAGS) -I $(OGSTM_INC) -o $@ $^ $(LDFLAGS)
compute2.xx: compute.F90 libogstm.a
	$(FC) $(FFLAGS) -I $(INCDIR) -o $@ $^ $(LDFLAGS)

# Libraries
#
.PHONY:libs
libs: dirs libogstm_lib.a libogstm.a

libogstm_lib.a: $(OGSTM_BLD)/libogstm_lib.a
	@ln -s $< $@
libogstm.a: $(MODULES) $(OBJECTS)
	ar qc $@ $^
	ranlib $@

# Objects
#
.PHONY: objects
objects: dirs modules $(OBJECTS)

$(BUILDDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/%.o: %.F
	$(FC) $(FFLAGS) -I $(INCDIR) -c -o $@ $<
	@mv $(@:.o=.optrpt) $(OPTRTDIR)

$(BUILDDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c -o $@ $<
	@mv $(@:.o=.optrpt) $(OPTRTDIR)

$(BUILDDIR)/%.o: %.F90
	$(FC) $(FFLAGS) -I $(INCDIR) -c -o $@ $<
	@mv $(@:.o=.optrpt) $(OPTRTDIR)

# Fortran modules
#
.PHONY: dirs modules
modules: $(MODULES)

build/stringop.o: stringop.F
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/calendar.o: calendar.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_data.o: bc_data.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/parameter.o: parameter.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc.o: bc.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/Timers.o: Timers.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/memory.o: memory.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_aux.o: bc_aux.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/hard_open.o: hard_open.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/rivers.o: rivers.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/sponge.o: sponge.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/nudging.o: nudging.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/TimeManager.o: TimeManager.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_handle.o: bc_handle.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_set.o: bc_set.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/closed.o: closed.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/OPT_mem.o: OPT_mem.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/SLAE.o: SLAE.F90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/BIO_mem.o: BIO_mem.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/FN_mem.o: FN_mem.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/DIA_mem.o: DIA_mem.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/SED_mem.o: SED_mem.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/get_mem_mod.o: get_mem_mod.F90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/ADV_mem.o: ADV_mem.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/BC_mem.o: BC_mem.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/IO_mem.o: IO_mem.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_atm.o: bc_atm.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/bc_co2.o: bc_co2.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/ogstm_mpi.o: ogstm_mpi.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/step.o: step.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/ZDF_mem.o: ZDF_mem.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)
build/ogstm.o: ogstm.f90
	$(FC) $(FFLAGS) -I $(INCDIR) -c $<
	@mv *.o $(BUILDDIR)
	@mv *.mod $(INCDIR)
	@mv *.optrpt $(OPTRTDIR)

# Clean
#
.PHONY: clean
clean: 
	@rm -rf *.o *.optrpt $(BUILDDIR) $(INCDIR) $(OPTRTDIR)

.PHONY: clean_libs
clean_libs:
	@rm -f *.a

.PHONY: clean_apps
clean_apps: 
	@rm -rf *.xx

.PHONY: cleanall
cleanall: clean clean_libs clean_apps

