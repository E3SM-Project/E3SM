.SUFFIXES: .F .c .o

CVMIX_REPO_ADDRESS=http://cvmix.googlecode.com/svn/trunk/src/shared

OCEAN_SHARED_INCLUDES=-I../shared -I../analysis_members -I../cvmix -I../../framework -I../../external/esmf_time_f90 -I../../operators

OCEAN_LIBRARIES=cvmix/*.o analysis_members/*.o shared/*.o

ifdef MODE

ifeq ($(wildcard ./mode_$(MODE)), ) # CHECK FOR EXISTENCE OF MODE DIRECTORY
all: exit

core_reg: exit

error_msg: error_head
	@echo "$(MODE) is not a valid build mode for the ocean core"

else # IFEQ ($(wildcard....

all: shared libcvmix analysis_members
	(cd mode_$(MODE); $(MAKE) FCINCLUDES="$(FCINCLUDES) $(OCEAN_SHARED_INCLUDES)" )
	if [ -e libdycore.a ]; then \
		($(RM) libdycore.a) \
	fi
	ar -ru libdycore.a $(OCEAN_LIBRARIES) mode_$(MODE)/*.o

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

endif # IFEQ ($(wildcard....

else # IFDEF MODE

all: exit

core_reg: exit

error_msg: error_head
	@echo "The ocean core requires a build mode."

endif # IFDEF MODE

libcvmix:
	if [ ! -d cvmix ]; then \
		(svn checkout $(CVMIX_REPO_ADDRESS) cvmix) \
	fi
	if [ -d cvmix ]; then \
		(cd cvmix; svn update; make all FC="$(FC)" FFLAGS="$(FFLAGS)" FINCLUDES="$(FINCLUDES)") \
	fi

shared: libcvmix
	(cd shared; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(OCEAN_SHARED_INCLUDES)")

analysis_members: libcvmix shared
	( cd analysis_members; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(OCEAN_SHARED_INCLUDES)" CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 

error_head:
	@echo ""
	@echo ""
	@echo "*************************************"
	@echo "ERROR"

error_tail: error_head error_msg
	@echo "Available build modes are:"
	@ls -d mode_* | grep ".*" | sed "s/mode_/    /g"
	@echo ""
	@echo "Please specify at build time as follows:"
	@echo "    make target CORE=ocean MODE=build_mode"
	@echo "*************************************"
	@echo ""
	@echo ""

exit: error_head error_msg error_tail
	@exit 1



clean:
	if [ -d cvmix ]; then \
		(cd cvmix; make clean) \
	fi
	(cd mode_forward; $(MAKE) clean)
	(cd mode_analysis; $(MAKE) clean)
	(cd analysis_members; $(MAKE) clean)
	(cd shared; $(MAKE) clean)
	($(RM) *.mod libdycore.a Registry_processed.xml)
