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

all: $(EXE)
ifeq ($(USEMPI),YES)
	mpirun -np 1 ./$(EXE)
else
	./$(EXE)
endif


SUT:
	make -C $(SRC_DIR) SUT
	make -C $(TEST_DIR) tests

tests: all

$(EXE): testSuites.inc SUT
	$(FC) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/include -Itests $(PFUNIT)/include/driver.F90 $(TEST_DIR)/*$(OBJ_EXT) $(SRC_DIR)/*$(OBJ_EXT) -L$(PFUNIT)/lib -lpfunit $(FFLAGS) $(FPPFLAGS)

clean: local-E0-clean

local-E0-clean:
	$(MAKE) -C $(SRC_DIR) clean
	$(MAKE) -C $(TEST_DIR) clean
	$(RM) -f $(EXE) *$(OBJ_EXT)

export FC
export FPPFLAGS
export FFLAGS
export SRC_DIR
export TEST_DIR
export OBJ_EXT
export LIB_EXT
export EXE_EXT
