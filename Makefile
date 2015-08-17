
.SUFFIXES: .F .o .cpp
.PHONY: mode_forward shared

all: shared mode_forward lib_landice core_landice

shared: 
	(cd shared; $(MAKE))

mode_forward: shared 
	(cd mode_forward; $(MAKE))

lib_landice: shared mode_forward

core_landice: lib_landice
	ar -ru libdycore.a mode_forward/*.o

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
