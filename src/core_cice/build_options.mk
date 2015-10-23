PWD=$(shell pwd)
EXE_NAME=cice_model
NAMELIST_SUFFIX=cice
FCINCLUDES += -I$(PWD)/src/core_cice/column -I$(PWD)/src/core_cice/forward_model -I$(PWD)/src/core_cice/analysis_members
override CPPFLAGS += -DCORE_CICE
ifneq "$(ESM)" ""
override CPPFLAGS += -Dcoupled -DCCSMCOUPLED
endif

report_builds:
	@echo "CORE=cice"
