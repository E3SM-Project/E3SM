PWD=$(shell pwd)
EXE_NAME=cice_model
NAMELIST_SUFFIX=cice
override CPPFLAGS += -DCORE_CICE

report_builds:
	@echo "CORE=cice"
