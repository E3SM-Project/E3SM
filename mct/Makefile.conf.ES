#!/bin/make
#-----------------------------------------------------------------------
#   Math and Computer Science Division, Argonne National Laboratory    |
#-----------------------------------------------------------------------
# CVS $Id$
# CVS $Name$ 
#-----------------------------------------------------------------------
#BOP
#
# !ROUTINE: Makefile.conf.SUPER-UX
#
# !DESCRIPTION:
#       An include file of configurations for the NEC SX platform.
#   See Makefile for additional `make' procedures.  This include file is
#   used to overide the default environments (including commands) defined
#   in file Makefile.
#
# !INTERFACE:
#
#       % ln -s Makefile.conf.SUPER-UX Makefile.conf
#       % make ...
#
# !SEE ALSO: Makefile
#
# !REVISION HISTORY: (`cvs log <>' for more)
#       25MAY02 - R. L. Jacob - initial version
#       23Oct03 - R. Jacob - incorporate mods from Yoshi.
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

_Iconf		= $(_IMPI)

# This module requires ../mpeu module

_Lmpeu		= -L../mpeu
_lmpeu		= -lmpeu
_Impeu		= -I../mpeu

_I		= -I. $(_Impeu) $(_IMPI)

SHELL		= /usr/bin/sh

CPP		= /usr/lib/cpp

F77		= esmpif90
F90		= esmpif90
FC		= $(F90)

#_Ftrace                = -ftrace
#_Ftrace                =  -vampir 
#_Ftrace                = -Wf,"-vampir" -vampir

#FFLAGS         = -C vsafe -Wf"-pvctl fullmsg -L fmtlist transform map" $(_Ftrace)
FFLAGS          = -C vopt -Wf"-pvctl fullmsg -L fmtlist transform map"

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


CPPFLAGS	= $(_Dsys)
FPPFLAGS	= -EP $(CPPFLAGS)

.F.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $(_I) $*.F

.F90.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $(_I) $*.F90

#.----------------------------------------------------------------------

