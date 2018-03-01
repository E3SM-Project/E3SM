.PHONY: clean tests all
.DEFAULT_GOAL = tests

TOP_DIR := $(shell pwd)
SRC_DIR=$(TOP_DIR)/src
TEST_DIR=$(TOP_DIR)/tests

VPATH = . $(SRC_DIR) $(TEST_DIR)

include $(PFUNIT)/include/base.mk

ifeq ($(USEMPI),YES)
   FC=mpif90
endif

EXE = tests$(EXE_EXT)

# NOTE: use the -robust flag
all: $(EXE)
ifeq ($(USEMPI),YES)
	echo  Warning:  Testing -robust flag with MPI, robust will override.
	echo            driver.F90 will call mpirun on itself with robust set.
	echo            "Top" level mpirun not needed, may be problemmatic.
	./$(EXE) -robust
	# mpirun -np 1 ./$(EXE) -robust
else
	./$(EXE) -robust -max-timeout-duration 0.16 -max-launch-duration 4.72
endif

SUT:
	make -C $(SRC_DIR) SUT
	make -C $(TEST_DIR) tests

tests: all

$(EXE): testSuites.inc myTests.pf add.F90 SUT
	$(FC) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/include -Itests $(PFUNIT)/include/driver.F90 $(TEST_DIR)/*$(OBJ_EXT) $(SRC_DIR)/*$(OBJ_EXT) -L$(PFUNIT)/lib -lpfunit $(LIBS) $(FFLAGS) $(FPPFLAGS)

clean: local-E0-clean

local-E0-clean:
	make -C $(SRC_DIR) clean
	make -C $(TEST_DIR) clean
	rm -f tests.x

export FC
export FPPFLAGS
export FFLAGS
export SRC_DIR
export TEST_DIR
export OBJ_EXT
export LIB_EXT
export EXE_EXT
