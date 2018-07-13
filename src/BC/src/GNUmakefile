.PHONY: clean
ifneq ($(BASEMK_INCLUDED),YES)
include $(PFUNIT)/include/base.mk
endif

SRC_DIR_GENERAL=/home/marco/MHPC/OGS_code/devel/ModelBuild/ogstm/src/General

SUT: bc_data$(OBJ_EXT)

bc_data$(OBJ_EXT): bc_data.f03 calendar$(OBJ_EXT)
	$(FC) -c $(FFLAGS) -std=f2003 $<

calendar$(OBJ_EXT): $(SRC_DIR_GENERAL)/calendar.f90 stringop$(OBJ_EXT)
	$(FC) -c $(FFLAGS) $<

stringop$(OBJ_EXT): $(SRC_DIR_GENERAL)/stringop.F
	$(FC) -c $(FFLAGS) $<

clean: local-E1-clean

local-E1-clean:
	$(RM) *$(OBJ_EXT) *.mod

