# ===================================
# Check if building with LifeV, Albany, and/or PHG external libraries

BUILD_INTERFACE=false  # This will become true if any of the external libraries are being used.

## LifeV can solve L1L2 or FO
#ifeq "$(LIFEV)" "true"
#	EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_L1L2
#	EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_FIRSTORDER
#	BUILD_INTERFACE = true
#endif # LIFEV IF
#
## Albany can only solve FO at present
#ifeq "$(ALBANY)" "true"
#	EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_FIRSTORDER
#	BUILD_INTERFACE = true
#endif # ALBANY IF
#
## Currently LifeV AND Albany is not allowed
#ifeq "$(LIFEV)" "true"
#ifeq "$(ALBANY)" "true"
#	$(error Compiling with both LifeV and Albany is not allowed at this time.)
#endif
#endif
#
## PHG currently requires LifeV
#ifeq "$(PHG)" "true"
#ifneq "$(LIFEV)" "true"
#	$(error Compiling with PHG requires LifeV at this time.)
#endif
#endif
## PHG can only Stokes at present
#ifeq "$(PHG)" "true"
#	EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_STOKES
#	BUILD_INTERFACE = true
#endif # PHG IF

#override CPPFLAGS += $(EXTERNAL_DYCORE_FLAG)
# ===================================

.SUFFIXES: .F .o .cpp
.PHONY: mode_forward shared

all: core_landice shared mode_forward

core_landice: 
#	ar -ru libdycore.a mode_forward/*.o
#	ar -ru libdycore.a shared/*.o

#core_reg:
#	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.landice )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.landice stream_list.landice. listed )

gen_includes:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml
	(if [ ! -d inc ]; then mkdir -p inc; fi) # To generate *.inc files
	(cd inc; $(REG_PARSE) < ../Registry_processed.xml )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )

shared: 
	(cd shared; $(MAKE))

mode_forward: 
	(cd mode_forward; $(MAKE))

clean:
	$(RM) *.o *.mod *.f90 libdycore.a
	$(cd shared; $(MAKE) clean)
	$(cd mode_forward; $(MAKE) clean)
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	$(RM) -r default_inputs

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $*.cpp $(CXINCLUDES) $(CPPINCLUDES) -lmpi_cxx -lstdc++ $(CPPFLAGS)
