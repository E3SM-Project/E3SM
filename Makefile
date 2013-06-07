SHELL		= /bin/sh
###############################
include ../Makefile.conf

VPATH=$(SRCDIR)/mpi-serial
# SOURCE FILES

MODULE		= mpi-serial

SRCS_F90	= fort.F90

SRCS_C		= mpi.c \
		  send.c \
		  recv.c \
		  collective.c \
		  req.c \
		  list.c \
		  handles.c \
                  comm.c \
                  group.c \
                  time.c \
                  pack.c


OBJS_ALL	= $(SRCS_C:.c=.o) \
		  $(SRCS_F90:.F90=.o)


INCPATH:= $(INCFLAG)$(SRCDIR)/mpi-serial $(INCFLAG). $(INCFLAG)../ $(INCPATH)

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

default: lib$(MODULE).a

examples: ctest ftest


MPIFH= mpif.$(FORT_SIZE).h


mpif.h: $(MPIFH)
	cp -f $< $@

fort.o: mpif.h

lib:
	@if [ ! "$(FORT_SIZE)" ] ; \
           then echo "Please set FORT_SIZE (e.g. real4double8 or real8double16) when you do the main MCT configure"; \
                exit 1; fi
	@if [ ! -r $(MPIFH) ] ; \
           then echo "Error: there is no $(MPIFH) -" \
                      "check the value of FORT_SIZE in the main MCT configure" ; \
                exit 1; fi
	cp -f $(MPIFH) mpif.h
	chmod -w mpif.h
	$(MAKE) $(LIB)



lib$(MODULE).a: $(OBJS_ALL)
	echo $(OBJS_ALL)
	$(RM) $@
	$(AR) $@ $(OBJS_ALL)
	$(RANLIB) $@


LIB	= lib$(MODULE).a


###############################
#RULES

.SUFFIXES:
.SUFFIXES: .F90 .c .o 

.c.o:
	$(CC) -c $(INCPATH) $(CPPDEFS) $(CFLAGS) $<

.F90.o:
	$(FC) -c $(INCFLAG). $(INCPATH) $(FPPDEFS) $(FCFLAGS) $(MPEUFLAGS) $<

MYF90FLAGS=$(INCPATH) $(DEFS) $(FCFLAGS)  $(MPEUFLAGS)

clean:
	/bin/rm -f *.o ctest ftest $(LIB) mpif.h


install: lib
	$(MKINSTALLDIRS) $(libdir) $(includedir)
	$(INSTALL) lib$(MODULE).a -m 644 $(libdir)
	$(INSTALL) mpi.h -m 644 $(includedir)
	$(INSTALL) mpif.h -m 644 $(includedir)



###############################
#
# Create mpif.realXdoubleY.h filesfrom mpif.master.h template
#

mpif:
	make-mpif 4 8
	make-mpif 8 8
	make-mpif 8 16


###############################

#
# test programs
#


ctest: lib ctest.c
	$(CC) $(ALLCFLAGS) -o $@ ctest.c -L. -lmpi-serial

ftest: lib ftest.F90
	$(FC) $(MYF90FLAGS) -o $@ ftest.F90 -L. -lmpi-serial

