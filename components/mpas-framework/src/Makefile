.SUFFIXES: .F .c .o

ifneq "$(ESM)" ""

include Makefile.in.$(ESM)

else

ifeq "$(AUTOCLEAN)" "true"
AUTOCLEAN_DEPS=clean_shared
else
AUTOCLEAN_DEPS=
endif

all: build_tools externals frame ops

externals: $(AUTOCLEAN_DEPS)
	( cd external; $(MAKE) FC="$(FC)" SFC="$(SFC)" CC="$(CC)" SCC="$(SCC)" FFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)" CPP="$(CPP)" NETCDF="$(NETCDF)" CORE="$(CORE)" all )

endif

build_tools: externals
	(cd tools; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(SCC)" CFLAGS="$(CFLAGS)")

frame: $(AUTOCLEAN_DEPS) externals
	( cd framework; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" FCINCLUDES="$(FCINCLUDES)" MOAB_PATH=$(MOAB_PATH) all ) 
	ln -sf framework/libframework.a libframework.a

ops: $(AUTOCLEAN_DEPS) externals frame
	( cd operators; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf operators/libops.a libops.a

clean: clean_shared

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
