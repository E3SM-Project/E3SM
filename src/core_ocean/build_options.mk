ifeq "$(ROOT_DIR)" ""
	ROOT_DIR=$(shell pwd)/src
endif
EXE_NAME=ocean_model
NAMELIST_SUFFIX=ocean
FCINCLUDES += -I$(ROOT_DIR)/core_ocean/driver -I$(ROOT_DIR)/core_ocean/mode_forward -I$(ROOT_DIR)/core_ocean/mode_analysis -I$(ROOT_DIR)/core_ocean/shared -I$(ROOT_DIR)/core_ocean/analysis_members -I$(ROOT_DIR)/core_ocean/cvmix
override CPPFLAGS += -DCORE_OCEAN

report_builds:
	@echo "CORE=ocean"
