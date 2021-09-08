PWD=$(shell pwd)
EXE_NAME=test_model
NAMELIST_SUFFIX=test
override CPPFLAGS += -DCORE_TEST

report_builds:
	@echo "CORE=test"
