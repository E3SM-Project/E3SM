.SUFFIXES: .F .c .o


OCEAN_SHARED_INCLUDES = -I$(PWD)/../framework -I$(PWD)/../external/esmf_time_f90 -I$(PWD)/../operators
OCEAN_SHARED_INCLUDES += -I$(PWD)/BGC -I$(PWD)/shared -I$(PWD)/analysis_members -I$(PWD)/cvmix -I$(PWD)/mode_forward -I$(PWD)/mode_analysis -I$(PWD)/mode_init

all: shared libcvmix analysis_members libBGC
	(cd mode_forward; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(OCEAN_SHARED_INCLUDES)" all )
	(cd mode_analysis; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(OCEAN_SHARED_INCLUDES)" all )
	(cd mode_init; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(OCEAN_SHARED_INCLUDES)" all )
	(cd driver; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(OCEAN_SHARED_INCLUDES)" all )
	if [ -e libdycore.a ]; then \
		($(RM) libdycore.a) \
	fi
	ar -ru libdycore.a `find . -type f -name "*.o"`

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.ocean )
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.ocean.forward mode=forward )
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.ocean.analysis mode=analysis )
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.ocean.init mode=init )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.ocean stream_list.ocean. mutable )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.ocean.forward stream_list.ocean.forward. mutable mode=forward )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.ocean.analysis stream_list.ocean.analysis. mutable mode=analysis )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.ocean.init stream_list.ocean.init. mutable mode=init )

gen_includes:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml
	(if [ ! -d inc ]; then mkdir -p inc; fi) # To generate *.inc files
	(cd inc; $(REG_PARSE) < ../Registry_processed.xml )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cp $(ROOT_DIR)/default_inputs/namelist.ocean $(ROOT_DIR)/namelist.ocean )
	( cp $(ROOT_DIR)/default_inputs/namelist.ocean.forward $(ROOT_DIR)/namelist.ocean.forward )
	( cp $(ROOT_DIR)/default_inputs/namelist.ocean.analysis $(ROOT_DIR)/namelist.ocean.analysis )
	( cp $(ROOT_DIR)/default_inputs/namelist.ocean.init $(ROOT_DIR)/namelist.ocean.init )
	( cp $(ROOT_DIR)/default_inputs/streams.ocean $(ROOT_DIR)/streams.ocean )
	( cp $(ROOT_DIR)/default_inputs/streams.ocean.forward $(ROOT_DIR)/streams.ocean.forward )
	( cp $(ROOT_DIR)/default_inputs/streams.ocean.analysis $(ROOT_DIR)/streams.ocean.analysis )
	( cp $(ROOT_DIR)/default_inputs/streams.ocean.init $(ROOT_DIR)/streams.ocean.init )

cvmix_source: get_cvmix.sh
	(/bin/bash ./get_cvmix.sh)
	(cd cvmix)

BGC_source: get_BGC.sh
	(/bin/bash ./get_BGC.sh)
	(cd BGC)

libcvmix: cvmix_source
	if [ -d cvmix ]; then \
		(cd cvmix; make all FC="$(FC)" FCFLAGS="$(FFLAGS)" FINCLUDES="$(FINCLUDES)") \
	else \
		(exit 1) \
	fi

libBGC: BGC_source
	if [ -d BGC ]; then \
		(cd BGC; make all FC="$(FC)" FCFLAGS="$(FFLAGS)" FINCLUDES="$(FINCLUDES)") \
	else \
		(exit 1) \
	fi

shared: libcvmix libBGC
	(cd shared; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(OCEAN_SHARED_INCLUDES)")

analysis_members: libcvmix shared
	( cd analysis_members; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(OCEAN_SHARED_INCLUDES)" CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 

clean:
	if [ -d cvmix ]; then \
		(cd cvmix; make clean) \
	fi
	if [ -d inc ]; then \
		($(RM) -r inc) \
	fi
	if [ -d BGC ]; then \
		(cd BGC; make clean) \
	fi
	(cd mode_forward; $(MAKE) clean)
	(cd mode_analysis; $(MAKE) clean)
	(cd mode_init; $(MAKE) clean)
	(cd driver; $(MAKE) clean)
	(cd analysis_members; $(MAKE) clean)
	(cd shared; $(MAKE) clean)
	($(RM) *.mod libdycore.a Registry_processed.xml)
	$(RM) -r default_inputs
