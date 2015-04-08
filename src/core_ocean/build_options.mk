ifeq "$(ROOT_DIR)" ""
	ROOT_DIR=$(shell pwd)/src
endif
ifeq "$(MODE)" "analysis"
	EXE_NAME=ocean_analysis_model
	NAMELIST_SUFFIX=ocean_analysis
	FCINCLUDES += -I$(ROOT_DIR)/core_ocean/mode_analysis -I$(ROOT_DIR)/core_ocean/shared -I$(ROOT_DIR)/core_ocean/analysis_members -I$(ROOT_DIR)/core_ocean/cvmix
else ifeq "$(MODE)" "forward"
	EXE_NAME=ocean_forward_model
	NAMELIST_SUFFIX=ocean_forward
	FCINCLUDES += -I$(ROOT_DIR)/core_ocean/mode_forward -I$(ROOT_DIR)/core_ocean/shared -I$(ROOT_DIR)/core_ocean/analysis_members -I$(ROOT_DIR)/core_ocean/cvmix
else
	EXE_NAME=ocean*_model
	NAMELIST_SUFFIX=ocean*
endif

report_builds:
	@echo "CORE=ocean MODE=analysis"
	@echo "CORE=ocean MODE=forward"
