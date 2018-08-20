.PHONY: clean
ifneq ($(BASEMK_INCLUDED),YES)
include $(PFUNIT)/include/base.mk
endif

SRC_DIR_GENERAL=/home/marco/MHPC/OGS_code/devel/ModelBuild/ogstm/src/General

SUT: sponge$(OBJ_EXT) rivers$(OBJ_EXT) nudging$(OBJ_EXT) bc$(OBJ_EXT) bc_data$(OBJ_EXT)

nudging$(OBJ_EXT): nudging.f03 bc$(OBJ_EXT)
	$(FC) -c -I$(NETCDFF_INCLUDE) $(FFLAGS) -std=f2003 $<

sponge$(OBJ_EXT): sponge.f03 bc$(OBJ_EXT) bc_aux_testing$(OBJ_EXT)
	$(FC) -c $(FFLAGS) -std=f2003 $<

rivers$(OBJ_EXT): rivers.f03 bc$(OBJ_EXT) bc_aux_testing$(OBJ_EXT)
	$(FC) -c $(FFLAGS) -std=f2003 $<

bc$(OBJ_EXT): bc.f03 bc_data$(OBJ_EXT)
	$(FC) -c $(FFLAGS) -std=f2003 $<

bc_data$(OBJ_EXT): bc_data.f03 calendar$(OBJ_EXT)
	$(FC) -c $(FFLAGS) -std=f2003 $<

calendar$(OBJ_EXT): $(SRC_DIR_GENERAL)/calendar.f90 stringop$(OBJ_EXT)
	$(FC) -c $(FFLAGS) $<

stringop$(OBJ_EXT): $(SRC_DIR_GENERAL)/stringop.F
	$(FC) -c $(FFLAGS) $<

bc_aux_testing$(OBJ_EXT): bc_aux_testing.f90
	$(FC) -c -I$(NETCDFF_INCLUDE) $(FFLAGS) $<

clean: local-E1-clean

local-E1-clean:
	$(RM) *$(OBJ_EXT) *.mod

