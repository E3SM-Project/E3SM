#!/bin/make

MODULE		= mpi0

SRCS_mks	= Makefile			\
		  Makefile.dependency		\
		  Makefile.rules		\
		  configure			\
		  unconfigure			\
		  README

# Single PE replacements of a minimum subset of libmpi.a, if there is
# no any "message-passing" call.

SRCS_F90	= mpi_abort.F90			\
		  mpi_comm_dup.F90		\
		  mpi_comm_free.F90		\
		  mpi_comm_rank.F90		\
		  mpi_comm_size.F90		\
		  mpi_error_string.F90		\
		  mpi_finalize.F90		\
		  mpi_init.F90			\
		  mpi_initialized.F90		\
		  mpi_wtick.F90			\
		  mpi_wtime.F90			\
		  mpi_barrier.F90		\
		  mpi_bcast.F90			\
		  mpi_gather.F90		\
		  mpi0_copy.F90

SRCS_F		= m_mpi0.F

SRCS_H		= Const_MPItype.H	\
		  mpif.H

SRCS_ALL	= $(SRCS_f)		\
		  $(SRCS_F)		\
		  $(SRCS_f90)		\
		  $(SRCS_F90)		\
		  $(SRCS_C)

OBJS_ALL	= $(SRCS_f:.f=.o)	\
		  $(SRCS_F:.F=.o)	\
		  $(SRCS_f90:.f90=.o)	\
		  $(SRCS_F90:.F90=.o)	\
		  $(SRCS_C:.c=.o)

ARCH_TARGETS	= $(SRCS_mks) $(SRCS_ALL) $(SRCS_H)
ALL_TARGETS	= lib$(MODULE).a
CLEAN_FILES	= mpif.h

include Makefile.rules
include Makefile.conf

# targets

lib$(MODULE).a:	$(OBJS_ALL)
	$(RM) $@
	@ echo $(AR) $@ $(OBJS_ALL)
	@ $(AR) $@ $(OBJS_ALL)

#   A special dependency, if there is not MPI on your system, defined
#
#	LOCAL_mpif_H=mpif.h
#
#   in your Makefile.conf.

m_mpi0.o: $(LOCAL_mpif_H)

# dependencies:

include Makefile.dependency

#.
