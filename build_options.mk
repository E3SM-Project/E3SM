PWD=$(shell pwd)
EXE_NAME=ocean_model
NAMELIST_SUFFIX=ocean
FCINCLUDES += -I$(PWD)/src/core_ocean/driver -I$(PWD)/src/core_ocean/mode_forward -I$(PWD)/src/core_ocean/mode_analysis -I$(PWD)/src/core_ocean/shared -I$(PWD)/src/core_ocean/analysis_members -I$(PWD)/src/core_ocean/cvmix

report_builds:
	@echo "CORE=ocean"
