# Top level Makefile for MiniCAM
# William.Gustafson@pnl.gov; 9-Sep-2008

export BUILDPATH = $(PWD)
include ./configure.gcam

OBJDIR = ./objs

all: gcam

gcam-prof: export PROF = -pg
gcam-prof: clean libgcam.a main_dir

gcam: libgcam.a main_dir

gcamccsm: libgcam.a ccsmcpl_dir

libgcam.a: dirs
	$(AR) libgcam.a $(OBJDIR)/*.o

dirs : containers_dir demographics_dir emissions_dir marketplace_dir resources_dir sectors_dir solution_solvers_dir solution_util_dir technologies_dir util_base_dir util_logger_dir util_database_dir util_curves_dir consumers_dir investment_dir reporting_dir climate_dir functions_dir target_finder_dir land_allocator_dir ccarbon_model_dir policy_dir #main_dir agLU_dir magicc_dir 

## special case patterns first
# util has subdirs
util_%_dir:
	$(MAKE) -C ../../util/$*/source $@
	ln -sf $(PWD)/../../util/$*/source/*.o $(OBJDIR)

#solution also has subdirs
solution_%_dir:
	$(MAKE) -C ../../solution/$*/source $@
	ln -sf $(PWD)/../../solution/$*/source/*.o $(OBJDIR)

#general pattern
%_dir:
	$(MAKE) -C ../../$*/source $@
	ln -sf $(PWD)/../../$*/source/*.o $(OBJDIR)

# main has additional instructions and doesn't do the softlink
main_dir :
	@ echo '----------------------------------------------------------------'
	rm -f ../../main/source/gcam.exe
	$(MAKE) -C ../../main/source  BUILDPATH=$(BUILDPATH) main_dir 
	cp ../../main/source/gcam.exe ../../../../exe/


# special because it has a funny directory name
agLU_dir :
	$(MAKE) -C ../../agLU/fortran_source agLU_dir
	ln -sf $(PWD)/../../agLU/fortran_source/*.o $(OBJDIR)

clean :
	-$(RM) libgcamlib.a
	-$(RM) objs/*.o
	-$(MAKE) -C ../../agLU/fortran_source  clean 
	-$(MAKE) -C ../../magicc/source  clean 
	-$(MAKE) -C ../../containers/source  clean 
	-$(MAKE) -C ../../demographics/source  clean 
	-$(MAKE) -C ../../emissions/source  clean 
	-$(MAKE) -C ../../marketplace/source  clean 
	-$(MAKE) -C ../../resources/source  clean 
	-$(MAKE) -C ../../sectors/source  clean 
	-$(MAKE) -C ../../solution/solvers/source  clean 
	-$(MAKE) -C ../../solution/util/source  clean 
	-$(MAKE) -C ../../technologies/source  clean 
	-$(MAKE) -C ../../util/base/source  clean 
	-$(MAKE) -C ../../util/logger/source  clean 
	-$(MAKE) -C ../../util/database/source  clean 
	-$(MAKE) -C ../../util/curves/source  clean 
	-$(MAKE) -C ../../consumers/source  clean 
	-$(MAKE) -C ../../investment/source  clean 
	-$(MAKE) -C ../../reporting/source  clean 
	-$(MAKE) -C ../../climate/source  clean 
	-$(MAKE) -C ../../functions/source  clean 
	-$(MAKE) -C ../../target_finder/source  clean 
	-$(MAKE) -C ../../land_allocator/source  clean 
	-$(MAKE) -C ../../ccarbon_model/source  clean 
	-$(MAKE) -C ../../policy/source  clean 
	-$(MAKE) -C ../../main/source  clean 
	-$(MAKE) -C ../../ccsmcpl/source  clean 
