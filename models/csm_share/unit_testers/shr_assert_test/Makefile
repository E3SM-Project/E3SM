# Makefile for tests of shr_assert_mod

CSM_SHARE_DIR ?= ../../

VPATH := .:$(CSM_SHARE_DIR)/shr

# Default to intel on yellowstone for now.
COMPILER ?= INTEL
OS ?= LINUX

# Shared flags
CPPFLAGS += -DCPR$(COMPILER) -D$(OS)
INCFLAGS += -I. -I$(CSM_SHARE_DIR)/shr
CFLAGS   += $(INCFLAGS)
FCFLAGS  += $(INCFLAGS)

ifeq ($(COMPILER),INTEL)
	CC = icc
	FC = ifort
	LD = $(FC)
	CPPFLAGS:= $(CPPFLAGS) -DFORTRANUNDERSCORE
else ifeq ($(COMPILER),GNU)
	CC = gcc
	FC = gfortran
	LD = $(FC)
	CPPFLAGS:= $(CPPFLAGS) -DFORTRANUNDERSCORE
	FCFLAGS += -fno-range-check
else ifeq ($(COMPILER),PGI)
	CC = pgcc
	FC = pgfortran
	LD = $(FC)
	CPPFLAGS:= $(CPPFLAGS) -DFORTRANUNDERSCORE
else ifeq ($(COMPILER),NAG)
	CC = gcc
	FC = nagfor
	LD = $(FC)
	CPPFLAGS:= $(CPPFLAGS) -DFORTRANUNDERSCORE
	FCFLAGS += -gline -C=all
else ifeq ($(COMPILER),LAHEY)
	CC = gcc
	FC = lf95
	LD = $(FC)
	CPPFLAGS:= $(CPPFLAGS) -DFORTRANUNDERSCORE
endif

# Let some paths/names be changed by user.
GENF90 ?=  $(CSM_SHARE_DIR)/../../tools/cprnc/genf90/genf90.pl
MKDEPENDS ?= $(CSM_SHARE_DIR)/../../scripts/ccsm_utils/Machines/mkDepends
EXENAME ?= test_shr_assert

# Make commands
.PHONY: all test clean distclean

all:  Depends $(EXENAME)

test: all
	$(EXENAME)

clean:
	$(RM) $(EXENAME) *.mod *.o Filepath Srcfiles Depends

distclean: clean
	$(RM) *.log

SOURCES := shr_kind_mod.F90 shr_log_mod.F90 shr_sys_mod.F90 \
	shr_isnan.c shr_test_infnan_mod.F90 shr_infnan_mod.F90 \
	shr_assert_mod.F90 test_utils.F90 logical_tests.F90 \
	domain_tests.F90 test_shr_assert.F90

OBJS := $(addsuffix .o, $(basename $(SOURCES)))

# Regenerate using genf90.pl
$(CSM_SHARE_DIR)/shr/shr_assert_mod.F90: shr_assert_mod.F90.in
	$(GENF90) $< >$@

# Don't use implicit rules.
.SUFFIXES:

%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

%.o: %.F90
	$(FC) -c $(CPPFLAGS) $(FCFLAGS) $<

# Linking
$(EXENAME): $(OBJS)
	$(LD) -o $@ $(CPPFLAGS) $(FCFLAGS) $(LDFLAGS) $^ $(LIBS)

# Dependency autogeneration
Depends: Filepath Srcfiles $(SOURCES)
	$(MKDEPENDS) Filepath Srcfiles >$@

space := $(null) $(null)
Srcfiles:
	@echo -e "$(subst $(space),\n,$(SOURCES))" >$@

Filepath:
	@echo -e "$(subst :,\n,$(VPATH))" >$@

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include Depends
endif
endif
