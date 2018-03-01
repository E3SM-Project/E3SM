# Makefile created by mkmf $Id: mkmf,v 18.0.18.4 2012/12/04 15:24:15 Seth.Underwood Exp $

ifeq ($(COMPILER),intel)
 FC := ifort
 LD := ifort
 LDFLAGS += -mkl
 CC := icc

# Sandy Bridge/Ivy Bridge

 FFLAGS := -O3 -xHost -fp-model fast -mkl -no-prec-div -no-prec-sqrt -override-limits

 CFLAGS := -O3 -xHost -fp-model fast -std=c99

 CPPDEFS = -DINTEL_MKL -DHAVE_SSE2


# Knights Corner

# FFLAGS := -mmic -O3 -qopt-report=5 -fp-model fast -no-prec-div -no-prec-sqrt -I./

# Haswell

# FFLAGS := -O3 -xCORE-AVX2 -no-prec-div -no-prec-sqrt  -I./ -DCPRINTEL

endif

ifeq ($(COMPILER),pgi)
 FC := pgfortran
 LD := pgfortran
 CC := pgcc

 FFLAGS := -fastsse
 CFLAGS := -fastsse
 CPPDEFS =
endif

ifeq ($(COMPILER),gnu)
 FC := gfortran
 LD := gfortran
 CC := gcc

 FFLAGS := -Ofast -march=native
 CFLAGS := -Ofast -march=native -std=gnu99
 CPPDEFS = -DHAVE_SSE2
endif

ifeq ($(COMPILER),nag)
 FC := nagfor
 LD := nagfor
 CC := gcc

 FFLAGS := -O4
 CFLAGS := -Ofast -march=native -std=gnu99
 CPPDEFS = -DHAVE_SSE2
endif

ifeq ($(COMPILER),ibm)
 FC := xlf2003
 LD := xlf2003
 CC := xlc

 FFLAGS := -O4
 CFLAGS := -O4
 CPPDEFS =
endif

CPPDEFS := $(CPPDEFS) -DDSFMT_MEXP=19937

ifeq ($(COMPILER),ibm)
  cpre = $(null)-WF,-D$(null)
  FPPDEFS := $(patsubst -D%,$(cpre)%,$(CPPDEF))
else
  FPPDEFS := $(CPPDEFS)
endif

FFLAGS := $(FFLAGS) $(FPPDEFS) -I../../include -I./
CFLAGS := $(CFLAGS) $(CPPDEFS) -I../../include

.DEFAULT:
	-echo $@ does not exist.
all: ./shr_RandNum.exe
dSFMT.o: ../../src/dsfmt_f03/dSFMT.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c ../../src/dsfmt_f03/dSFMT.c
dSFMT_interface.o: ../../src/dsfmt_f03/dSFMT_interface.F90
	$(FC) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c ../../src/dsfmt_f03/dSFMT_interface.F90
dSFMT_utils.o: ../../src/dsfmt_f03/dSFMT_utils.c ../../include/dSFMT.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c ../../src/dsfmt_f03/dSFMT_utils.c
kissvec_mod.o: ../../src/kissvec/kissvec_mod.F90
	$(FC) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c ../../src/kissvec/kissvec_mod.F90
kissvec.o: ../../src/kissvec/kissvec.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c ../../src/kissvec/kissvec.c
test_shr_RandNum.o: ./test_shr_RandNum.F90 shr_RandNum_mod.o
	$(FC) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c ./test_shr_RandNum.F90
mersennetwister_mod.o: ../../src/mt19937/mersennetwister_mod.F90
	$(FC) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c ../../src/mt19937/mersennetwister_mod.F90
shr_RandNum_mod.o: ../../src/shr_RandNum_mod.F90 kissvec_mod.o mersennetwister_mod.o dSFMT_interface.o
	$(FC) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c ../../src/shr_RandNum_mod.F90
SRC = ./test_shr_RandNum.F90 ../kissvec/kissvec_mod.F90 ../mt19937/mersennetwister_mod.F90 ../dsfmt_f03/dSFMT.c ../dsfmt_f03/dSFMT_interface.F90 ../dsfmt_f03/dSFMT_utils.c ../shr_RandNum_mod.F90 ../../include/dSFMT.h
OBJ = test_shr_RandNum.o kissvec_mod.o mersennetwister_mod.o dSFMT.o dSFMT_interface.o dSFMT_utils.o shr_RandNum_mod.o kissvec.o
clean:
	-rm -f .shr_RandNum.exe.cppdefs $(OBJ) *.mod ./shr_RandNum.exe *.s
shr_RandNum.exe: $(OBJ)
	$(LD) $(OBJ) -o shr_RandNum.exe  $(LDFLAGS)
