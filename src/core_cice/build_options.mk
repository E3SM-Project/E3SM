PWD=$(shell pwd)
EXE_NAME=cice_model
NAMELIST_SUFFIX=cice
FCINCLUDES += -I$(PWD)/src/core_cice/column -I$(PWD)/src/core_cice/shared -I$(PWD)/src/core_cice/analysis_members -I$(PWD)/src/core_cice/model_forward
override CPPFLAGS += -DCORE_CICE

report_builds:
	@echo "CORE=cice"
