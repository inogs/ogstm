.PHONY: clean tests
ifneq ($(BASEMK_INCLUDED),YES)
include $(PFUNIT)/include/base.mk
endif

SRCS = $(wildcard *.pf)
OBJS = $(SRCS:.pf=$(OBJ_EXT))

# Uncomment and replace only for debugging
# tests: $(OBJS) test_sponge_default.F90 test_rivers_default.F90
tests: $(OBJS)

testSuites.inc: $(SRCS)

FFLAGS += -I$(SRC_DIR) -I$(PFUNIT)/mod

ifndef SKIP_INTENTIONALLY_BROKEN
  FFLAGS += -DINCLUDE_INTENTIONALLY_BROKEN
endif

%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $< $@

%$(OBJ_EXT): %.F90
	$(FC) -c $(FFLAGS) $(FPPFLAGS) $<

clean: local-E1-clean

local-E1-clean:
	$(RM) *$(OBJ_EXT) *.mod
