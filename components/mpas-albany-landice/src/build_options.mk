ifeq "$(ROOT_DIR)" ""
        ROOT_DIR=$(shell pwd)/src
endif
EXE_NAME=landice_model
NAMELIST_SUFFIX=landice
FCINCLUDES += -I$(ROOT_DIR)/mode_forward -I$(ROOT_DIR)/shared -I$(ROOT_DIR)/analysis_members
override CPPFLAGS += -DCORE_LANDICE

# ===================================
# Check if building with LifeV, Albany, and/or PHG external libraries

# LifeV can solve L1L2 or FO
ifeq "$(LIFEV)" "true"
    EXTERNAL_DYCORE_FLAG += -DLIFEV
    EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_L1L2
    EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_FIRSTORDER
    EXTERNAL_DYCORE_FLAG += -DMPAS_LI_BUILD_INTERFACE
endif # LIFEV IF

# Albany can only solve FO at present
ifeq "$(ALBANY)" "true"
    EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_FIRSTORDER
    EXTERNAL_DYCORE_FLAG += -DMPAS_LI_BUILD_INTERFACE
endif # ALBANY IF

# Currently LifeV AND Albany is not allowed
ifeq "$(LIFEV)" "true"
ifeq "$(ALBANY)" "true"
    $(error Compiling with both LifeV and Albany is not allowed at this time.)
endif
endif

# PHG currently requires LifeV
ifeq "$(PHG)" "true"
ifneq "$(LIFEV)" "true"
    $(error Compiling with PHG requires LifeV at this time.)
endif
endif
# PHG can only Stokes at present
ifeq "$(PHG)" "true"
    EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_STOKES
    EXTERNAL_DYCORE_FLAG += -DMPAS_LI_BUILD_INTERFACE
endif # PHG IF

override CPPFLAGS += $(EXTERNAL_DYCORE_FLAG)
# ===================================

# Optional Sea Level model
ifeq "$(SLM)" "true"
    override CPPFLAGS += -DUSE_SEALEVELMODEL
endif

# ===================================
report_builds:
	@echo "CORE=landice"
