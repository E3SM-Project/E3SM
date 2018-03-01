
# GNU makefile include for NAG compilers. Ver. 2014-0428-1 MLR

F90 ?= nagfor

I=-I
M=-I
L=-L

FFLAGS += -g -O0 -f2008 -w=uda -mismatch_all -fpp -C=present

ifeq ($(USEOPENMP),YES)
FFLAGS += -openmp
else
FFLAGS += -gline
endif

CPPFLAGS += -DNAG

F90_PP_ONLY = -F
F90_PP_OUTPUT = -o

# For OS X Mavericks (i.e. Apple LLVM version 6.0...), bring your own CPP.
# CPP = /opt/local/bin/cpp-mp-4.9 -traditional -C
# CPP =cpp -traditional -C

ifeq ($(DSO),YES)
  FFLAGS +=-PIC
endif

LDFLAGS+= -ldl
