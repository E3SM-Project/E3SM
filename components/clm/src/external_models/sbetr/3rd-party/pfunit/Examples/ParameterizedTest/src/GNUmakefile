.PHONY: clean
ifneq ($(BASEMK_INCLUDED),YES)
include $(PFUNIT)/include/base.mk
endif

SUT: SphericalCoordinates$(OBJ_EXT)

%$(OBJ_EXT): %.F90
	$(FC) -c $(FFLAGS) $<

clean: local-E1-clean

local-E1-clean: 
	$(RM) *$(OBJ_EXT) *.mod

