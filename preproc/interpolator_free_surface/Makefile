include $(HOME)/DEV/compilers/compiler.inc

INSTDIR	= .
LDFLAGS += -L$(LIBDIR) -lnc-medlevel -lf_exit
OBJS	= modulo16.o quartic_solver.o div16.o wzv16.o comp_ssh16.o stringop.o calendar.o ForcingGenerator.o

EXE	= ForcingGenerator

target: $(EXE)

ForcingGenerator:$(OBJS) 
		$(LD) $(OBJS)  $(LDFLAGS) -o $@


clean:	
	@rm -f $(OBJS) *.mod

erase:	clean
	@rm -f $(EXE)

install:	$(EXE)
	cp -p $(EXE) $(INSTDIR)

uninstall:	$(EXE)
	rm -f $(INSTDIR)/$(EXE)

.SUFFIXES: .f90 .F90 .F .o
.f90.o:
	$(FC) $(FFLAGS) $< -c -g

.F90.o:
	$(FC) $(FFLAGS) $< -c -g

.F.o:
	$(FC) $(FFLAGS) $< -c -extend_source

