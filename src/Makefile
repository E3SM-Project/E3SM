.SUFFIXES: .F .c .o

ifneq "$(ESM)" ""

include Makefile.in.$(ESM)

else

ifeq "$(AUTOCLEAN)" "true"
AUTOCLEAN_DEPS=clean_shared
else
AUTOCLEAN_DEPS=
endif

all: mpas

mpas: $(AUTOCLEAN_DEPS) externals frame ops dycore drver
	$(LINKER) $(LDFLAGS) -o $(EXE_NAME) driver/*.o -L. -ldycore -lops -lframework $(LIBS) -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time

externals: $(AUTOCLEAN_DEPS)
	( cd external; $(MAKE) FC="$(FC)" SFC="$(SFC)" CC="$(CC)" SCC="$(SCC)" FFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)" CPP="$(CPP)" NETCDF="$(NETCDF)" CORE="$(CORE)" all )

drver:  $(AUTOCLEAN_DEPS) externals frame ops dycore
	( cd driver; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
endif

build_tools: externals
	(cd tools; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(SCC)" CFLAGS="$(CFLAGS)")

frame: $(AUTOCLEAN_DEPS) externals
	( cd framework; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf framework/libframework.a libframework.a

ops: $(AUTOCLEAN_DEPS) externals frame
	( cd operators; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf operators/libops.a libops.a

dycore: $(AUTOCLEAN_DEPS) build_tools externals frame ops
	( cd core_$(CORE); $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" REG_PARSE="$(PWD)/tools/registry/parse" gen_includes )
	( cd core_$(CORE); $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" NL_GEN="$(PWD)/tools/input_gen/namelist_gen" ST_GEN="$(PWD)/tools/input_gen/streams_gen" core_input_gen )
	( cd core_$(CORE); $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf core_$(CORE)/libdycore.a libdycore.a


clean: clean_shared clean_core

clean_core:
	if [ -d core_$(CORE) ] ; then \
	   ( cd core_$(CORE); $(MAKE) clean ) \
	fi;

clean_shared:
ifeq "$(AUTOCLEAN)" "true"
	@echo ""
	@echo "*********************************************************************************************"
	@echo "The MPAS infrastructure is currently built for a core different from $(CORE)."
	@echo "The infrastructure will be cleaned and re-built for the $(CORE) core."
	@echo "*********************************************************************************************"
	@echo ""
endif
	$(RM) libframework.a libops.a libdycore.a lib$(CORE).a *.o
	( cd tools; $(MAKE) clean )
	( cd external; $(MAKE) clean )
	( cd framework; $(MAKE) clean )
	( cd operators; $(MAKE) clean )
	( cd driver; $(MAKE) clean )
