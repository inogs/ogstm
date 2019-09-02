CPPOPTS=  -C -P -I. -Dsun -Ucray 

FPP= 	fpp $(CPPOPTS)
#FTRAN= gfortran -g -DRAD2HR
FTRAN= ifort -g -DRAD2HR -I $(NETCDFF_INC)

LOUT=	clrtrans.o lidata.o light.o light_verbose.o ocalbedo.o rdatopt.o IOnc.o setdate.o \
	setlte.o sfcirr.o sfcsolz.o slingo.o slingo_verbose.o sunmod.o findbinxy2.o \
	monavgpar.o getarea.o monavgswr.o monavgday0m.o 

monrad: monrad.F $(LOUT) 
	$(FTRAN) -L $(NETCDFF_LIB)  -lnetcdff -o monrad monrad.F $(LOUT) 

clrtrans.o: clrtrans.F
	$(FTRAN) -c clrtrans.F

lidata.o: lidata.F
	$(FTRAN) -c lidata.F

light.o: light.F
	$(FTRAN) -c light.F

light_verbose.o: light_verbose.F
	$(FTRAN) -c light_verbose.F

ocalbedo.o: ocalbedo.F
	$(FTRAN) -c ocalbedo.F

rdatopt.o: rdatopt.F
	$(FTRAN) -c rdatopt.F

IOnc.o: IOnc.f90
	$(FTRAN) -c IOnc.f90

rstfle.o: rstfle.F
	$(FTRAN) -c rstfle.F

setdate.o: setdate.F
	$(FTRAN) -c setdate.F

setlte.o: setlte.F
	$(FTRAN) -c setlte.F

sfcirr.o: sfcirr.F
	$(FTRAN) -c sfcirr.F

sfcsolz.o: sfcsolz.F
	$(FTRAN) -c sfcsolz.F

slingo.o: slingo.F
	$(FTRAN) -c slingo.F

slingo_verbose.o: slingo_verbose.F
	$(FTRAN) -c slingo_verbose.F

sunmod.o: sunmod.F
	$(FTRAN) -c sunmod.F

findbinxy2.o: findbinxy2.F
	$(FTRAN) -c findbinxy2.F

monavgswr.o: monavgswr.F
	$(FTRAN) -c monavgswr.F

monavgpar.o: monavgpar.F
	$(FTRAN) -c monavgpar.F

monavgday0m.o: monavgday0m.F
	$(FTRAN) -c monavgday0m.F

getarea.o: getarea.F
	$(FTRAN) -c getarea.F


