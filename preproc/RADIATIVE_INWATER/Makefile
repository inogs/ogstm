# Generic Makefile to compile c-postproc utils
# for Code_Saturne
#
# Arnau Miro - UPC (ESEIAAT)

# Compilers
#
CC       = icc
CFLAGS   = -Ofast -xHost -Wall -fPIC
#-std=gnu11
CXX      = icpc
CXXFLAGS = -Ofast -xHost -Wall -fPIC
#-std=c++11
BUILDDIR = BUILD
OGSTM_PATH = /galileo/home/userexternal/plazzari/CODE/ogstm
GENERAL_PATH = $(OGSTM_PATH)/General
BIO_OPTICS_PATH = $(OGSTM_PATH)/BIO-OPTICS
MPI_PATH = $(OGSTM_PATH)/MPI

SOURCES := $(shell find $(OGSTM_PATH)  -name '*.f90' -or  -name '*.F90' -or -name '*.F' )
# Get list of object files, with paths
OBJECTS := $(patsubst $(SOURCES)/%.f90,$(BUILDDIR)/%.o,$(SOURCES))
#OBJECTS += $(patsubst $(BUILDDIR)/,$(SOURCES:%.F90=%.o))
#OBJECTS += $(patsubst $(BUILDDIR)/,$(SOURCES:%.F=%.o))



default:
	@echo "Compile with make <option> to generate a program."
	@echo "Options:"

# One rule to compile them all, one rule to find them,
# One rule to bring them all and in the compiler link them.
all: libs search_neg join_case join_probe \
	time_average cs-toolslib
	@echo "All programs compiled successfully"

# Deprecated binaries in favor of the cs-toolslib
bin: export_field
	@echo "All programs compiled successfully"

# Applications
#
# Libraries
#
libs: ogstm.a
	@echo "Libraries compiled successfully"
general.o: $(wildcard $(GENERAL_PATH)/*.F90)
	@echo $(OBJECTS)

ogstm.a: general.o
	ar rc $@ $<
	ranlib $@
	@mv $@ $(LPATH)

# Generic object makers
#
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $< $(DFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(DFLAGS)

%.o: %.F90
	$(FC) $(FCLAGS) -c -o $@ $< $(DFLAGS)

%.o: %.f90
	$(FC) $(FCLAGS) -c -o $@ $< $(DFLAGS)

%.o: %.F
	$(FC) $(FCLAGS) -c -o $@ $< $(DFLAGS)

# Clean
#
clean:
	@rm -f *.o */*.o

remove_libs:
	@rm -f *.a */*.a
