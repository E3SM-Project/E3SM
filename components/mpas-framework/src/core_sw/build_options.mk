PWD=$(shell pwd)
EXE_NAME=sw_model
NAMELIST_SUFFIX=sw
override CPPFLAGS += -DCORE_SW

report_builds:
	@echo "CORE=sw"
