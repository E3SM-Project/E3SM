#!/bin/make

MODULE		= mct

SRCS_mks	=
include SRCS_mks.mk

SRCS_C		=

SRCS_f		=
include SRCS_f.mk

SRCS_F		=
include SRCS_F.mk

SRCS_f90	=
include SRCS_f90.mk

SRCS_F90	=
include SRCS_F90.mk

TEXFILES	=
include SRCS_tex.mk

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
ALL_TARGETS	= lib$(MODULE).a

#	Try to define _Ddft, _Dsys and _Dnod in proper Makefiles.
#
#		1) _D for the make() time macros, as in
#		make _D="-Dxxx ..."
#
#		2) _Ddft for the default global macros defined
#		in the main Makefile.
#
#		3) _Dsys for the OS dependent macros defined
#		in the main Makefile and can be redefined in
#		the Makefile.conf file.
#
#		4) _Dnod for the OS or the node dependent
#		macros defined in the Makefile.conf file.

_Ddft	=
_Dsys	= -Dsys`uname -s`

_DS	= $(_Ddft) $(_Dsys) $(_Dnod) $(_D)

include Makefile.rules
include Makefile.conf

# targets

lib$(MODULE).a:	$(OBJS_ALL)
	$(RM) $@
	@ echo $(AR) $@ $(OBJS_ALL)
	@ $(AR) $@ $(OBJS_ALL)

doc:	$(TEXFILES)

do-it-anyway:

#________________________________________
# Regular targets

#________________________________________
# Just dependencies

include Makefile.dependency

#.
