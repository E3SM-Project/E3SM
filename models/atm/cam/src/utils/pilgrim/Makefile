#-------------------------------------------------------------------------
#         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
#-------------------------------------------------------------------------
#
# !ROUTINE: Makefile
#
# !DESCRIPTION:
#
#     Makefile for GFIO Library and documentation.
#     You must enter
#
#       ./configure
#
#     before attempting to make anything in this Makefile.
#
# !SEE ALSO:
#
#     configure
#     Makefile.conf.*
#
#
# !REVISION HISTORY:
#
#  14Mar02  Sawyer    Adapted previous Makefile for FVDAS
#
#
#-----------------------------------------------------------------------
#
#                System Dependent Parameters
#                ---------------------------
 
COREROOT = ../../..
COREBIN  = $(COREROOT)/bin
CORELIB  = $(COREROOT)/lib
COREINC  = $(COREROOT)/include
COREETC  = $(COREROOT)/etc

CP		= /bin/cp -p
RANLIB		= touch
MKDIR		= mkdir -p

PILGRIMLIB      =  libpilgrim.a
SEQSRCS         =  debugutilitiesmodule.F90 \
                   decompmodule.F90 ghostmodule.F90
PUMSRCS         =  mod_comm.F90 parutilitiesmodule.F90 \
                   redistributemodule.F90 puminterfaces.F90
CSRCS           =  memstuff.c mp_assign_to_cpu.c

SEQOBJS         = ${SEQSRCS:.F90=.o}
PUMOBJS         = ${PUMSRCS:.F90=.o} 
COBJS           = ${CSRCS:.c=.o} 

OPTIONS         = $(_IMPI)

.SUFFIXES:	.F90 .c .o
# -------------------------------------------------------------------

include Makefile.conf

All:	mpi

all:	mpi

mpi:	spmdclean library parutil

mlp:	spmdclean library mlpparutil

tests:	mpi
	(cd unit_testers; make tests)

library:	$(SEQOBJS)
	-@${RM} -f $(PILGRIMLIB)
	$(AR) $(PILGRIMLIB) $(SEQOBJS)

parutil:	$(PUMOBJS)
	$(AR) $(PILGRIMLIB) $(PUMOBJS)

mlpparutil:	$(COBJS)
	$(FC) $(FFLAGS) $(DFLAGS) $(OPTIONS) -DUSE_MLP -c $(PUMSRCS)
	$(AR) $(PILGRIMLIB) $(PUMOBJS) $(COBJS)

# Export library (fvDAS specific)
# -------------------------------
export: mpi
	$(MKDIR) $(CORELIB) $(COREINC)/pilgrim
	$(CP) libpilgrim.a $(CORELIB)
	$(CP) *.[Mm][Oo][Dd]   $(COREINC)/pilgrim
	$(RANLIB) $(CORELIB)/libpilgrim.a

.F90.o:
	$(FC) $(FFLAGS) $(DFLAGS) $(OPTIONS) -c $<

.c.o:
	${CC} ${INCLUDE} ${CFLAGS} -DUSE_MLP $(OPTIONS) -c $<

clean:
	-@${RM} -f $(PILGRIMLIB) $(SEQOBJS) $(PUMOBJS) $(COBJS) *~ *.mod *.MOD
	-@${RM} -rf rii_files

libclean:
	-@${RM} -f $(PILGRIMLIB)

spmdclean:
	-@${RM} -f $(PUMOBJS)
