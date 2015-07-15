ifeq "$(ROOT_DIR)" ""
        ROOT_DIR=$(shell pwd)/src
endif
EXE_NAME=landice_model
NAMELIST_SUFFIX=landice
FCINCLUDES += -I$(ROOT_DIR)/core_landice
override CPPFLAGS += -DCORE_LANDICE

report_builds:
	@echo "CORE=landice"
