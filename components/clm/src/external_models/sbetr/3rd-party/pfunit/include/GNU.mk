F90 ?= gfortran

I=-I
M=-I
L=-L

FFLAGS += -g -O0 -fbacktrace
FFLAGS += -fbounds-check -fcheck=mem
FPPFLAGS += -DGNU

# The ramifications across all GNUish configurations of eliding CPPFLAGS here are not known. MLR 2013-1104
CPPFLAGS += -DGNU

F90_PP_ONLY = -E
F90_PP_OUTPUT = >

ifeq ($(USEOPENMP),YES)
FFLAGS += -fopenmp
LIBS += -lgomp
endif

