PWD=$(shell pwd)
EXE_NAME=init_atmosphere_model
NAMELIST_SUFFIX=init_atmosphere
override CPPFLAGS += -DCORE_INIT_ATMOSPHERE

report_builds:
	@echo "CORE=init_atmosphere"
