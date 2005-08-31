
#
# for debugging
# -DCHECKS -DHANDLE_INFO
#


include ../Makefile.conf

#
# The values used from Makefile.conf
#

# ALLCFLAGS= -DFORTRAN_UNDERSCORE_
# ALLCFLAGS= -DFORTRAN_SAME
# ALLCFLAGS= -DFORTRAN_CAPS

# FC=pgf90
# AR=ar rv
# CC=cc

###############################


.SUFFIXES: .F90 .c


.F90.o:
	$(FC) $(FFLAGS) -c -o $@ $<

.c.o:
	$(CC) $(ALLCFLAGS) -c -o $@ $<


###############################

LIBOBJS= mpi.o send.o recv.o collective.o req.o list.o fort.o handles.o

LIB=libmpi.a

$(LIB)(%.o): %.o
	$(AR) $(LIB) $%


LIB_MEMBERS= $(foreach file, $(LIBOBJS), $(LIB)($(file)))


#
# Default target is library

$(LIB): $(LIB_MEMBERS)
	@echo "Done building $@"


###############################

#
# test programs
#


all: ctest ftest


ctest: $(LIB) ctest.c
	$(CC) $(ALLCFLAGS) -o $@ ctest.c -L. -lmpi

ftest: $(LIB) ftest.F90
	$(FC) -o $@ ftest.F90 -L. -lmpi


clean:
	/bin/rm -f *.o ctest ftest $(LIB)

