# -*- mode: makefile -*-
#
# makefile controling generation and execution of unit tests with pFUnit
#

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

SRC_DIR = ..
PFUNIT_DIR = ../../../tpls/pfunit

include $(PFUNIT_DIR)/include/extensions.mk

TEST_SRCS = $(wildcard *.pf)
TEST_OBJS = $(TEST_SRCS:.pf=$(OBJ_EXT))

# this list should only include the modules in the test suite (and
# their dependancies)! Probably need a more complex .o dependancy
# setup like the main makefile...
PFLOTRAN_OBJS = \
	$(SRC_DIR)/eos_gas.o \
	$(SRC_DIR)/eos_water.o \
	$(SRC_DIR)/input_aux.o \
	$(SRC_DIR)/geometry.o \
	$(SRC_DIR)/logging.o \
	$(SRC_DIR)/option.o \
	$(SRC_DIR)/option_flow.o \
	$(SRC_DIR)/option_transport.o \
	$(SRC_DIR)/pflotran_constants.o \
    $(SRC_DIR)/characteristic_curves.o \
	$(SRC_DIR)/string.o \
	$(SRC_DIR)/units.o \
	$(SRC_DIR)/utility.o

FFLAGS = -g -I. -I$(SRC_DIR) -I$(PFUNIT_DIR)/source -I$(PFUNIT_DIR)/include

all : test

test : pfunit unittests
	@echo "----------------------------------------------------------------------"
	@echo "Running pflotran unit tests :"
	@./unittests
	@echo "----------------------------------------------------------------------"

unittests : testSuites.inc $(TEST_OBJS) $(PFLOTRAN_OBJS)
	$(FC) -o $@ $(PFUNIT_DIR)/include/driver.F90 $(FFLAGS) \
	$(TEST_OBJS) $(PFLOTRAN_OBJS) -L$(PFUNIT_DIR)/source -lpfunit  ${PETSC_LIB}

#
# generic rules for converting from .pf --> .F90 --> .o
#
%.F90 : %.pf
	$(PYTHON) $(PFUNIT_DIR)/bin/pFUnitParser.py $< $@

%$(OBJ_EXT) : %.F90
	$(FC) -c $(FC_FLAGS) $(FFLAGS) ${FCPPFLAGS} $<

#
# build pfunit if necessary
#
pfunit : $(PFUNIT_DIR)/source/libpfunit.a

$(PFUNIT_DIR)/source/libpfunit.a :
	./build-pfunit.sh -c $(FC) -d $(PFUNIT_DIR)

#
# cleanup: note, this removes F90 files!
#
clean-tests : clean-pfunit
	-rm -f unittests *~ *.F90 *.o *.mod
	-rm -rf unittests.dSYM

clean-pfunit :
	./build-pfunit.sh -c $(FC) -d $(PFUNIT_DIR) -n

# uncomment the PRECIOUS rule to keep the intermediate F90 files
#.PRECIOUS : %.F90

