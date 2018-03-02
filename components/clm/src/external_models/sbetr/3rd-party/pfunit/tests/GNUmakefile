
# Gnu makefile for $(TOP)/tests/

TOP = ..
include ../include/base.mk
.PHONY: all tests
EXE = tests$(EXE_EXT)
REMOTE = remote$(EXE_EXT)

all: dependencies.inc $(EXE) $(REMOTE)

VPATH = . ../source
tests: dependencies.inc $(EXE)

#generated.inc:
#	touch generated.inc

# include ./generated.inc

dependencies.inc: dependenciesTests.tmp
	./MakeDependenciesInc dependenciesTests.tmp

# Would prefer to have this in base.mk, but haven't figured out the rules to do so.
# gmake looks at ../source/dependencies.tmp, instead of making a new one at ".".
#

#		$(FC) -E $(FFLAGS) $(CPPFLAGS) $$i > $${i}_cpp.f90;
dependenciesTests.tmp: 
	for i in *.F90; do \
		$(F90) $(F90_PP_ONLY) $(FFLAGS) $(CPPFLAGS) $$i $(F90_PP_OUTPUT) $${i}_cpp.f90; \
		$(TOP)/tools/MakeDependenciesTmp $$i $${i}_cpp.f90 \
			$@ "PrivateException|MockSUT|SUT"; \
		$(RM) $${i}_cpp.f90; done

SRCS_BASE = $(wildcard *.F90)
PROGRAMS = selfTests.F90 RemoteProgram.F90
EXCLUDE = $(PROGRAMS)
ifneq ($(BUILDROBUST),YES)
EXCLUDE += RobustRunner.F90 Test_RobustRunner.F90 Test_UnixProcess.F90
endif


SRCS = $(filter-out $(EXCLUDE),$(SRCS_BASE))
OBJS = $(SRCS:%.F90=%$(OBJ_EXT))

ifneq ($(USEMPI),YES)
   SRCS_SERIAL=$(filter-out Test_Mpi%,$(SRCS))
   OBJS_SERIAL=$(SRCS_SERIAL:%.F90=%$(OBJ_EXT))
endif

FFLAGS +=$M$(MOD_DIR)

ifneq ($(UNAME),Windows)
LIBS = $L$(LIB_DIR) -lpfunit

ifeq ($(USEMPI),YES)
   LIBS += $(LIBMPI)
endif
else
LIBS = $(LIB_DIR)/libpfunit$(LIB_EXT)
ifeq ($(USEMPI),YES)
# TODO: This must still be addressed for MPI
   LIBS += $(LIBMPI)
endif
endif

$(EXE) $(REMOTE): libpfunit$(LIB_EXT)

# Some platforms need a space between the OUTPUT_FLAG and $@.  Others don't.
# So we put an optional space into OUTPUT_FLAG.  See GNUmakefile in the top 
# directory.
ifeq ($(USEMPI),YES)

$(EXE): selfTests$(OBJ_EXT) $(OBJS)
	$(FC) $(OUTPUT_FLAG)$@ $^ $(FFLAGS) $(LIBS) -g
$(REMOTE): RemoteProgram$(OBJ_EXT) $(OBJS)
	$(FC) $(OUTPUT_FLAG)$@ $^ $(FFLAGS) $(LIBS) -g

else

$(EXE): selfTests$(OBJ_EXT) $(OBJS_SERIAL)
	$(FC) $(OUTPUT_FLAG)$@ $^ $(FFLAGS) $(LIBS) $(DEBUG_FLAGS)
$(REMOTE): RemoteProgram$(OBJ_EXT) $(OBJS_SERIAL)
	$(FC) $(OUTPUT_FLAG)$@ $^ $(FFLAGS) $(LIBS) $(DEBUG_FLAGS)

endif

# Only include the dependencies file if we are not doing a clean, otherwise the
# Files are generated, only to be deleted again.
ifeq (,$(findstring distclean,$(MAKECMDGOALS)))
include ./dependencies.inc
endif

