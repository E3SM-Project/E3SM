PWD=$(shell pwd)
EXE_NAME=landice_model
NAMELIST_SUFFIX=landice
override CPPFLAGS += -DCORE_LANDICE

report_builds:
	@echo "CORE=landice"
