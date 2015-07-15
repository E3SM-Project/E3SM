.SUFFIXES: .F .o
.PHONY: forward_model analysis_members

all: column_package forward_model analysis_members core_cice

core_cice:
	ar -ru libdycore.a column/*.o
	ar -ru libdycore.a forward_model/*.o
	ar -ru libdycore.a analysis_members/*.o

gen_includes:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml
	(if [ ! -d inc ]; then mkdir -p inc; fi) # To generate *.inc files
	(cd inc; $(REG_PARSE) < ../Registry_processed.xml )

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.cice )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.cice stream_list.cice. listed )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )

column_package:
	(cd column; $(MAKE))

forward_model: column_package
	(cd forward_model; $(MAKE))

analysis_members: column_package forward_model
	(cd analysis_members; $(MAKE))

clean:
	$(RM) *.o *.mod *.f90 libdycore.a
	(cd column; $(MAKE) clean)
	(cd forward_model; $(MAKE) clean)
	(cd analysis_members; $(MAKE) clean)
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	$(RM) -r default_inputs
