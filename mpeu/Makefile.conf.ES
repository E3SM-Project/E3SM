#!/bin/make
#-----------------------------------------------------------------------
#   Math and Computer Science Division, Argonne National Laboratory    |
#-----------------------------------------------------------------------
# CVS $Id$
# CVS $Name$  
#-----------------------------------------------------------------------
#BOP
#
# !ROUTINE: Makefile.conf.ES-UX
#
# !DESCRIPTION:
#       An include file of configurations for the NEC SX platform.
#   See Makefile for additional `make' procedures.  This include file is
#   used to overide the default environments (including commands) defined
#   in file Makefile.
#
# !INTERFACE:
#
#       % ln -s Makefile.conf.ES Makefile.conf
#       % make ...
#
# !SEE ALSO: Makefile
#
# !REVISION HISTORY: (`cvs log <>' for more)
#       25MAY02 - R. L. Jacob - initial version
#
#EOP
#-----------------------------------------------------------------------

_Dsys		= -DsysSUPER_UX
_Dnod		= 

# Environments
# ============

# If use the vendor implementation

_LMPI		= 
_lMPI		= 
_IMPI		= 

#_R		= 

_Iconf		= $(_IMPI)
_I		= $(_IMPI)

SHELL		= /usr/bin/sh

CPP		= /usr/lib/cpp

F77		= esmpif90
F90		= esmpif90
FC		= $(F90)
CC		= escc


FFLAGS          = $(_R) -C vopt -Wf"-pvctl fullmsg -L fmtlist transform map"

CFLAGS		= -O $(_Dsys)

ZIP		= /usr/bin/compress
ZCAT		= /usr/bin/zcat
TAR		= /usr/bin/tar

LD		= esmpif90
LDFLAGS		= $(FFLAGS)

AR		= esar cqs
#RM		= rm -f
MV		= /usr/bin/mv

LIBprof		= 
LIBsys	 	= 


CPPFLAGS	=$(_Dsys)
FPPFLAGS	= -EP $(CPPFLAGS)

.F.o:
	$(FC) $(FPPFLAGS) $(FFLAGS) $(_I) -c $*.F

.F90.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $(_I) $*.F90

#.----------------------------------------------------------------------

MUPDATE         = sh LCupdate.sh
