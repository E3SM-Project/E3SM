TOP = ..
include ../include/base.mk
.PHONY: all generate clean_generated_code testCodeMacs

all: generated.inc dependencies.inc libpfunit$(LIB_EXT)

# "GENERATED_CODE" is defined in "generated.inc".
# GENERATED_CODE = AssertReal.F90 AssertComplex.F90 AssertInteger1.F90
# 
generated.inc:
	$(PYTHON) ./GenerateAssertsOnArrays.py --maxRank $(MAX_RANK)

dependencies.inc: generated.inc dependenciesSource.tmp
	./MakeDependenciesInc dependenciesSource.tmp

# Would prefer to have this in base.mk, but haven't figured out the rules to do so.
# gmake looks at ../source/dependencies.tmp, instead of making a new one at ".".
#

#		$(FC) -E $(FFLAGS) $(CPPFLAGS) $$i > $${i}_cpp.f90;
dependenciesSource.tmp: generated.inc
	for i in *.F90; do \
		$(F90) $(F90_PP_ONLY) $(FFLAGS) $(CPPFLAGS) $$i $(F90_PP_OUTPUT) $${i}_cpp.f90; \
		$(TOP)/tools/MakeDependenciesTmp $$i $${i}_cpp.f90 \
			$@ "PrivateException|MockSUT|SUT"; \
		$(RM) $${i}_cpp.f90; done

# 	@echo "# This file automatically generated. Do not modify." > ./generated.inc

AssertArrays.fh: generated.inc

Assert.F90: AssertArrays.fh

include ./generated.inc

# $(GENERATED_CODE):
#	./GenerateAssertsOnArrays.py

testCodeMacs:
	$(PYTHON) $(TOP)/bin/codeMacs.py --inFile testCM.cm --outFile testCM.F90

# %.F90: %.cm
#	$(PYTHON) $(TOP)/bin/codeMacs.py $< $@

clean_generated_code:
	rm -f generated.inc AssertArrays.fh
	rm -f $(GENERATED_CODE)

distclean: local-build1-distclean

local-build1-distclean: clean_generated_code


ifneq ($(BUILDROBUST),YES)
EXCLUDE += UnixProcess.F90 RobustRunner.F90
endif
SRCS = $(filter-out $(EXCLUDE),$(wildcard *.F90))
OBJS = $(SRCS:%.F90=%$(OBJ_EXT))

EXCLUDE_MPI_SRCS = MpiContext.F90 MpiStubs.F90 MpiTestCase.F90 MpiTestMethod.F90

ifneq ($(USEMPI),YES)
   SRCS_SERIAL=$(filter-out $(EXCLUDE_MPI_SRCS),$(SRCS))
   OBJS_SERIAL=$(SRCS_SERIAL:%.F90=%$(OBJ_EXT))
endif

ifeq ($(USEMPI),YES)
   libpfunit$(LIB_EXT): $(OBJS)
else
   libpfunit$(LIB_EXT): $(OBJS_SERIAL)
endif
	$(AR)$@ $?
	-@$(RANLIB) $@

# Only include the dependencies file if we are not doing a clean.
ifeq (,$(findstring distclean,$(MAKECMDGOALS)))
include ./dependencies.inc
endif

