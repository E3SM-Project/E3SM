.PHONY: clean distclean tests all
.DEFAULT_GOAL = tests

TOP_DIR := $(shell pwd)
TEST_DIR=$(TOP_DIR)/tests

VPATH = . $(TEST_DIR)

include $(PFUNIT)/include/base.mk

ifeq ($(USEMPI),YES)
   MPIF90=mpif90
endif

EXE = tests$(EXE_EXT)

all: $(EXE)
	mpirun -np 4 ./$(EXE)

SUT:
	make -C $(TEST_DIR) tests

tests: all

$(EXE): testSuites.inc parameterizedTests.pf SUT
	$(MPIF90) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/include -Itests $(PFUNIT)/include/driver.F90 $(TEST_DIR)/*$(OBJ_EXT) -L$(PFUNIT)/lib -lpfunit $(LIBS) -DUSE_MPI $(FFLAGS) $(FPPFLAGS)

clean: local-E0-clean

local-E0-clean:
	make -C $(TEST_DIR) clean
	rm -f $(EXE) *$(OBJ_EXT)

export FC
export FPPFLAGS
export FFLAGS
export TEST_DIR
export OBJ_EXT
export LIB_EXT
export EXE_EXT
export MPIF90
