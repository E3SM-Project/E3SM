#!/bin/make
#-----------------------------------------------------------------------
#   Math and Computer Science Division, Argonne National Laboratory    |
#-----------------------------------------------------------------------
# CVS $Id$
# CVS $Name$  
#-----------------------------------------------------------------------
#BOP
#
# !ROUTINE: Makefile.conf.ES
# 
# !DESCRIPTION: 
#	An include file of configurations for the NEC SX platform.
#   See Makefile for additional `make' procedures.  This include file is 
#   used to overide the default environments (including commands) defined 
#   in file Makefile.
#
# !INTERFACE:
#
#	% ln -s Makefile.conf.ES Makefile.conf
#	% make ...
#
# !SEE ALSO: Makefile
#
# !REVISION HISTORY: (`cvs log <>' for more)
#	25MAY02 - R. L. Jacob - initial version
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
_IMPI		= -I/usr/include

#_R		= -ew

_Iconf		= $(_IMPI)
_I		= $(_IMPI)

SHELL		= /usr/bin/sh

CPP		= /usr/lib/cpp

F77		= mpif90
F90		= mpif90
FC		= $(F90)

_Ftrace                = -ftrace
#_Ftrace                = -Wf,"-vampir" -vampir
#_Ftrace                = -vampir

# orignal
#FFLAGS		= $(_R) -C hopt -Wf"-pvctl noassoc"
FFLAGS		= $(_R) -C vopt -Wf"-pvctl fullmsg -L fmtlist" $(_Ftrace)
CFLAGS		= -O $(_Dsys)

ZIP		= /usr/bin/compress
ZCAT		= /usr/bin/zcat
TAR		= /usr/bin/tar

LD		= mpif90
LDFLAGS		= $(FFLAGS)

#AR		= ar cqs
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

