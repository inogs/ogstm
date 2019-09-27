# Generic Makefile to compile c-postproc utils
# for Code_Saturne
#
# Arnau Miro - UPC (ESEIAAT)

# Config
#
OGSTM_SRC  = /galileo/home/userexternal/plazzari/CODE/ogstm/src
OGSTM_COMP = /galileo/home/userexternal/plazzari/CODE/ogstm/compilers
BUILDDIR   = build
COMPILERS  = x86_64.LINUX.intel

# Sources and Objects
#
SOURCES := $(shell find $(OGSTM_SRC)  -name '*.f90' -or  -name '*.F90' -or -name '*.F' )
SRC_OUT := $(shell find "$(OGSTM_SRC)/BC/src/"  -name '*.f90' -or  -name '*.F90' -or -name '*.F')
SRC_OUT += $(OGSTM_SRC)/General/ogstm.f90
SOURCES := $(filter-out $(SRC_OUT),$(SOURCES))
# Convert .f90, .F90 and .F to .o
OBJECTS := $(SOURCES:.f90=.o)
OBJECTS := $(OBJECTS:.F90=.o)
OBJECTS := $(OBJECTS:.F=.o)
# Change to build dir
OBJECTS := $(foreach obj,$(OBJECTS),$(patsubst $(shell dirname $(obj))/%,$(BUILDDIR)/%,$(obj)))

# Compilers
#
include $(OGSTM_COMP)/$(COMPILERS).inc

LIBS = -logstm

# Default
#
default: ogstm.a compute.xx

# Applications
#
compute.xx: compute.o
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

# Libraries
#
ogstm.a: $(OBJECTS)
	ar rc $@ $<
	ranlib $@

# Generic object makers
#
$(OBJECTS): $(SOURCES)
	@mkdir -p $(BUILDDIR)
	$(FC) $(FFLAGS) -c -o $@ $<
%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
%.o: %.F
	$(FC) $(FFLAGS) -c -o $@ $<

# Clean
#
remove_libs:
	@rm -f *.a

clean: remove_libs
	@rm -f build
