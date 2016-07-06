
.SUFFIXES: .F .o .cpp
.PHONY: mode_forward shared analysis_members subglacial_hydro

SHARED_INCLUDES  = -I$(PWD)/../framework -I$(PWD)/../external/esmf_time_f90 -I$(PWD)/../operators
SHARED_INCLUDES += -I$(PWD)/shared -I$(PWD)/analysis_members -I$(PWD)/mode_forward -I$(PWD)/subglacial_hydro

all: core_landice

shared: 
	(cd shared; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(SHARED_INCLUDES)")

analysis_members: shared
	(cd analysis_members; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(SHARED_INCLUDES)")

subglacial_hydro: shared
	(cd subglacial_hydro; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(SHARED_INCLUDES)")

mode_forward: shared analysis_members subglacial_hydro
	(cd mode_forward; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(SHARED_INCLUDES)")

core_landice: mode_forward shared analysis_members subglacial_hydro
	ar -ru libdycore.a `find . -type f -name "*.o"`

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.landice )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.landice stream_list.landice. listed )

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

gen_includes:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml
	(if [ ! -d inc ]; then mkdir -p inc; fi) # To generate *.inc files
	(cd inc; $(REG_PARSE) < ../Registry_processed.xml )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )

clean:
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	$(RM) -r default_inputs
	(cd shared; $(MAKE) clean)
	(cd mode_forward; $(MAKE) clean)
	(cd analysis_members; $(MAKE) clean)
	(cd subglacial_hydro; $(MAKE) clean)
