F90 ?= pgfortran

I=-I
M=-I
L=-L

FFLAGS += -O0 -g -traceback -Mallocatable=03 -Mbounds -Mchkfpstk -Mchkstk -DPGI

ifeq ($(USEOPENMP),YES)
FFLAGS += -mp
endif

FPPFLAGS += -DPGI
CPPFLAGS += -DPGI -Mpreprocess

F90_PP_ONLY = -E
F90_PP_OUTPUT = >

ifeq ($(DSO),YES)
  FFLAGS +=-PIC
endif

LDFLAGS+= -ldl




