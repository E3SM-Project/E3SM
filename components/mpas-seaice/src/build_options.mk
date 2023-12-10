ifeq "$(ROOT_DIR)" ""
   ROOT_DIR=$(shell pwd)/src
endif
EXE_NAME=seaice_model
NAMELIST_SUFFIX=seaice
FCINCLUDES += -I$(ROOT_DIR)/column -I$(ROOT_DIR)/shared -I$(ROOT_DIR)/analysis_members -I$(ROOT_DIR)/model_forward
override CPPFLAGS += -DCORE_SEAICE -DUSE_SNICARHC
ifneq "$(ESM)" ""
override CPPFLAGS += -Dcoupled -DCCSMCOUPLED
endif

report_builds:
	@echo "CORE=seaice"
