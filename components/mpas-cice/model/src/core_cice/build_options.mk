ifeq "$(ROOT_DIR)" ""
   ROOT_DIR=$(shell pwd)/src
endif
EXE_NAME=cice_model
NAMELIST_SUFFIX=cice
FCINCLUDES += -I$(ROOT_DIR)/core_cice/column -I$(ROOT_DIR)/core_cice/shared -I$(ROOT_DIR)/core_cice/analysis_members -I$(ROOT_DIR)/core_cice/model_forward
override CPPFLAGS += -DCORE_CICE
ifneq "$(ESM)" ""
override CPPFLAGS += -Dcoupled -DCCSMCOUPLED
endif

report_builds:
	@echo "CORE=cice"
