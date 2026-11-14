ifeq "$(ROOT_DIR)" ""
	ROOT_DIR=$(shell pwd)/src
endif
EXE_NAME=ocean_model
NAMELIST_SUFFIX=ocean
# The registry parser writes the generated .inc files here
CPPINCLUDES += -I$(ROOT_DIR)/inc
FCINCLUDES += -I$(ROOT_DIR)/inc
FCINCLUDES += -I$(ROOT_DIR)/driver
FCINCLUDES += -I$(ROOT_DIR)/mode_forward -I$(ROOT_DIR)/mode_analysis -I$(ROOT_DIR)/mode_init
FCINCLUDES += -I$(ROOT_DIR)/shared -I$(ROOT_DIR)/analysis_members
FCINCLUDES += -I$(ROOT_DIR)/cvmix/src/shared
FCINCLUDES += -I$(ROOT_DIR)/BGC
FCINCLUDES += -I$(ROOT_DIR)/MARBL/include
FCINCLUDES += -I$(ROOT_DIR)/gotm/build/modules
FCINCLUDES += -I$(ROOT_DIR)/ppr/src
override CPPFLAGS += -DCORE_OCEAN

report_builds:
	@echo "CORE=ocean"
