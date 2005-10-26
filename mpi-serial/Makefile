SHELL		= /bin/sh

# SOURCE FILES

MODULE		= mpi

SRCS_F90	= fort.F90

SRCS_C		= mpi.c \
		  send.c \
		  recv.c \
		  collective.c \
		  req.c \
		  list.c \
		  handles.c \
                  comm.c \
                  group.c


OBJS_ALL	= $(SRCS_C:.c=.o) \
		  $(SRCS_F90:.F90=.o)


###############################

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

# TARGETS

default: lib

examples: ctest ftest


MPIFH= mpif$(BITS).h

fort.o: mpif.h

lib:
	@if [ ! $(BITS) ] ; \
           then echo "Please setenv BITS (32 or 64)"; \
                exit 1; fi
	@if [ ! -r $(MPIFH) ] ; \
           then echo "Error: there is no $(MPIFH) -" \
                      "check the environment value of BITS" ; \
                exit 1; fi
	cp -f mpif$(BITS).h mpif.h
	chmod -w mpif.h
	$(MAKE) $(LIB)



lib$(MODULE).a: $(OBJS_ALL)
	echo $(OBJS_ALL)
	$(RM) $@
	$(AR) $@ $(OBJS_ALL)


LIB	= lib$(MODULE).a



###############################
#RULES

.SUFFIXES:
.SUFFIXES: .F90 .c .o

$(CRULE):
	$(CC) -c $(ALLCFLAGS) $*.c

$(F90RULE):
	$(FC) -c $(INCPATH) $(DEFS) $(FCFLAGS) $(F90FLAGS) $(MPEUFLAGS) $*.F90

$(F90RULECPP):
	$(FPP) $(DEFS) $(FPPFLAGS) $*.F90 $*.f90 
	$(FC) -c $(INCPATH) $(FCFLAGS) $(F90FLAGS) $(MPEUFLAGS) $*.f90
	$(RM) $*.f90


clean:
	/bin/rm -f *.o ctest ftest $(LIB) mpif.h


install: lib
	$(MKINSTALLDIRS) $(libdir) $(includedir)
	$(INSTALL) lib$(MODULE).a -m 644 $(libdir)
	$(INSTALL) mpi.h -m 644 $(includedir)
	$(INSTALL) mpif.h -m 644 $(includedir)

###############################

#
# test programs
#


ctest: lib ctest.c
	$(CC) $(ALLCFLAGS) -o $@ ctest.c -L. -lmpi

ftest: lib ftest.F90
	$(FC) -o $@ ftest.F90 -L. -lmpi

