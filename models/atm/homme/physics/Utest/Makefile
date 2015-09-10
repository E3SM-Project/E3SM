ARCH =$(shell uname -s)
include ../../Params.inc
include ../../Makefile.$(ARCH)

#NETCDFINC=/usr/local/include
#NETCDFLIB=/usr/local/lib32/r4i4
#NCARG=yes


VPATH=../emanuel/ ../../src ../crcp
SRC = kinds.F90 physical_constants.F90 \
      control_mod.F90 convect43c.F90 physicsUtest.F90

OBJ = $(SRC:.F90=.o)

all: physicsUtest 

#FFLAGS = -g -qfullpath -O2 -qarch=auto -WF,-DCONVECT43C
#,-DUSE43B

ifeq ('$(NCARG)','yes')
  F90 = ncargf90 -qsuffix=cpp=F90
  F77 = ncargf77
  FFLAGS += -WF,-DNCARG
else
  F90 = xlf90_r -qsuffix=cpp=F90
  F77 = xlf_r
endif
FC=$(F90)

ifneq ('$(NETCDFINC)','')
  FFLAGS += -I$(NETCDFINC) -I$(NETCDFLIB) -WF,-DNETCDF
  LDFLAGS += -L$(NETCDFLIB) -lnetcdf
endif



physicsUtest: $(OBJ)
	$(F90) $(FFLAGS) $(OBJ) -o $@


%.o: %.F90
	$(F90) $(FFLAGS) -c $<


$(ORIGOBJ): %.o: %.F
	$(F77) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJ) $(ORIGOBJ)

origtest: $(ORIGOBJ)
	$(F77) $(FFLAGS) $(ORIGOBJ) -o $@