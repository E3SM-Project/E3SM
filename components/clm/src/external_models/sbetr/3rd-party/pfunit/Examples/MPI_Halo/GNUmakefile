.PHONY: clean distclean tests all
.DEFAULT_GOAL = tests

TOP_DIR := $(shell pwd)
SRC_DIR=$(TOP_DIR)/src
TEST_DIR=$(TOP_DIR)/tests

VPATH = . $(SRC_DIR) $(TEST_DIR)

include $(PFUNIT)/include/base.mk

ifeq ($(USEMPI),YES)
   MPIF90=mpif90
endif

EXE = tests$(EXE_EXT)

all: $(EXE)
	mpirun -np 4 ./$(EXE)

SUT:
	make -C $(SRC_DIR) SUT
	make -C $(TEST_DIR) tests

tests: all

$(EXE): testSuites.inc myTests.pf Halo.F90 SUT
	$(MPIF90) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/include -Itests $(PFUNIT)/include/driver.F90 $(TEST_DIR)/*$(OBJ_EXT) $(SRC_DIR)/*$(OBJ_EXT) -L$(PFUNIT)/lib -lpfunit $(LIBS) -DUSE_MPI $(FFLAGS) $(FPPFLAGS)

echo: local-E0-echo

local-E0-echo:
	@echo PFUNIT: $(PFUNIT)
	@echo BUILDROBUST: \"$(BUILDROBUST)\"
	@echo COMPILER: $(COMPILER)
	@echo FFLAGS: $(FFLAGS)
	@echo FPPFLAGS: $(FPPFLAGS)
	@echo FC: $(FC)
	@echo MPI: $(MPI)
	@echo CPPFLAGS: $(CPPFLAGS)
	@echo LIBS: $(LIBS)

clean: local-E0-clean

local-E0-clean:
	make -C $(SRC_DIR) clean
	make -C $(TEST_DIR) clean
	rm -f $(EXE) *$(OBJ_EXT)

export FC
export FPPFLAGS
export FFLAGS
export SRC_DIR
export TEST_DIR
export OBJ_EXT
export LIB_EXT
export EXE_EXT
export MPIF90
