PWD=$(shell pwd)
ifeq "$(MODE)" "analysis"
	EXE_NAME=ocean_analysis_model
	NAMELIST_SUFFIX=ocean_analysis
	FCINCLUDES += -I$(PWD)/src/core_ocean/mode_analysis -I$(PWD)/src/core_ocean/shared -I$(PWD)/src/core_ocean/analysis_members -I$(PWD)/src/core_ocean/cvmix
else ifeq "$(MODE)" "forward"
	EXE_NAME=ocean_forward_model
	NAMELIST_SUFFIX=ocean_forward
	FCINCLUDES += -I$(PWD)/src/core_ocean/mode_forward -I$(PWD)/src/core_ocean/shared -I$(PWD)/src/core_ocean/analysis_members -I$(PWD)/src/core_ocean/cvmix
else
	EXE_NAME=ocean*_model
	NAMELIST_SUFFIX=ocean*
endif

report_builds:
	@echo "CORE=ocean MODE=analysis"
	@echo "CORE=ocean MODE=forward"
