PWD=$(shell pwd)
EXE_NAME=atmosphere_model
NAMELIST_SUFFIX=atmosphere
override CPPFLAGS += -DCORE_ATMOSPHERE

report_builds:
	@echo "CORE=atmosphere"
