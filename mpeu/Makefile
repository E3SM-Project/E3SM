#!/bin/make

MODULE		= mpeu

SRCS_mks	= Makefile			\
		  Makefile.conf.in              \
		  Makefile.dependency		\
		  Makefile.rules		\
		  configure			\
		  unconfigure			\
		  UCupdate.sh			\
		  LCupdate.sh			\
		  externals.mk			\
		  README

SRCS_C		= get_zeits.c

SRCS_F		= m_mpif.F

SRCS_F90	= m_IndexBin_char.F90		\
		  m_IndexBin_integer.F90	\
		  m_IndexBin_logical.F90	\
		  m_List.F90			\
		  m_MergeSorts.F90		\
		  m_Filename.F90		\
		  m_Permuter.F90		\
		  m_SortingTools.F90		\
		  m_String.F90			\
		  m_StrTemplate.F90		\
		  m_chars.F90			\
		  m_die.F90			\
		  m_dropdead.F90		\
                  m_FileResolv.F90		\
		  m_flow.F90			\
		  m_inpak90.F90			\
		  m_ioutil.F90			\
		  m_mall.F90			\
		  m_mpif90.F90			\
		  m_mpout.F90			\
		  m_rankMerge.F90		\
		  m_realkinds.F90		\
		  m_stdio.F90			\
		  m_TraceBack.F90		\
		  m_zeit.F90			

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

ARCH_TARGETS	= $(SRCS_mks) $(SRCS_ALL)
ALL_TARGETS	= libeu.a lib$(MODULE).a $(MODULE)
SUBD_TARGETS	= mpi0
CLEAN_FILES	= mpif.h

M		= mod
MODS_OBJ        = m_chars.o		\
		  m_die.o		\
		  m_flow.o		\
		  m_inpak90.o 		\
		  m_ioutil.o 		\
		  m_List.o 		\
		  m_mall.o		\
		  m_mpif90.o		\
		  m_mpout.o		\
		  m_Filename.o		\
		  m_Permuter.o		\
		  m_realkinds.o		\
		  m_SortingTools.o	\
		  m_stdio.o		\
		  m_String.o		\
		  m_StrTemplate.o	\
		  m_zeit.o              \
                  m_FileResolv.o        \
                  m_MergeSorts.o	\
                  m_TraceBack.o		\
		  assert.H

include Makefile.rules
include Makefile.conf

# targets

libeu.a: $(OBJS_ALL) mpi0/libmpi0.a
	$(RM) $@
	cp -p mpi0/libmpi0.a $@
	@ echo $(AR) $@ $(OBJS_ALL)
	@ $(AR) $@ $(OBJS_ALL)

mpi0/libmpi0.a:
	@(cd mpi0; make libmpi0.a)


lib$(MODULE).a:	$(OBJS_ALL)
	$(RM) $@
	@ echo $(AR) $@ $(OBJS_ALL)
	@ $(AR) $@ $(OBJS_ALL)

$(MODULE): $(MODS_OBJ)
	env M=$M $(MUPDATE) $@ $?

#   A special dependency, if there is not MPI on your system, defined
#
#	LOCAL_mpif_H=mpif.h
#
#   in your Makefile.conf.

m_mpif.o: $(LOCAL_mpif_H)

# dependencies:

include Makefile.dependency

#.
