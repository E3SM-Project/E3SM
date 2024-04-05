EXEDIR = .
EXENAME = cube_to_target
RM = rm

.SUFFIXES:
.SUFFIXES: .F90 .o

FC = ifort
DEBUG = FALSE


# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /gpfs/fs1/soft/chrysalis/spack/opt/spack/linux-centos8-x86_64/intel-20.0.4/netcdf-fortran-4.4.4-rdxohvp/lib
#/global/common/software/nersc/pm-2023q1/spack-stacks-1/views/climate-utils/lib
#/public/software/mathlib/netcdf/4.3.2/intel/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /gpfs/fs1/soft/chrysalis/spack/opt/spack/linux-centos8-x86_64/intel-20.0.4/netcdf-fortran-4.4.4-rdxohvp/include
#/global/common/software/nersc/pm-2023q1/spack-stacks-1/views/climate-utils/include
#/public/software/mathlib/netcdf/4.3.2/intel/include
endif

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

#------------------------------------------------------------------------
# LF95
#------------------------------------------------------------------------
#
# setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/lib
#
ifeq ($(FC),lf95)
#
# Tramhill
#
  INC_NETCDF :=/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/include
  LIB_NETCDF :=/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/lib

  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff -lcurl -lhdf5 -lhdf5_hl -mcmodel=medium 
  FFLAGS   := -c --trace --trap --wide -CcdRR8 -I$(INC_NETCDF)
  ifeq ($(DEBUG),TRUE)
#   FFLAGS += --chk aesu  -Cpp --trace
    FFLAGS += -g --chk a,e,s,u --pca
  else
    FFLAGS += -O
  endif

endif



.F90.o:
	$(FC) $(FFLAGS) $<


#------------------------------------------------------------------------
# AIX
# #------------------------------------------------------------------------
#
 #ifeq ($(UNAMES),AIX)
 FC = ifort #xlf90
 #FFLAGS = -c -I$(INC_NETCDF) -I/BIGDATA1/iapcas_mhzhang_xiejinbo/topo_tool/cube_to_target/functional/  -convert big_endian

 FFLAGS = -c -I$(INC_NETCDF)   -convert big_endian  -traceback
 #FFLAGS        := -c -I$(INC_NETCDF)  -no-prec-div -traceback -convert big_endian -fp-model source  -assume byterecl -ftz   -m64  -mcmodel=large -safe-cray-ptr
 LDFLAGS = -L$(LIB_NETCDF)  -lnetcdff
 #LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff -m64 -static-intel
 .F90.o:
	$(FC) $(FFLAGS) -qsuffix=f=F90 $<
#         #endif


.F90.o:
	$(FC) $(FFLAGS) $<










#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

#OBJS := reconstruct.o remap.o cube_to_target.o shr_kind_mod.o
OBJS := reconstruct.o remap.o shr_kind_mod.o transform.o  sub_xjb.o  cube_to_target.o 
#OBJS := reconstruct.o remap.o cube_to_target.o sub.o shr_kind_mod.o
#sub.o  shr_kind_mod.o

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) -I$(INC_NETCDF) $(LDFLAGS) 

clean:
	$(RM) -f $(OBJS)  *.mod $(EXEDIR)/$(EXENAME)

cube_to_target.o: shr_kind_mod.o remap.o reconstruct.o   sub_xjb.o  transform.o
remap.o: 
reconstruct.o: remap.o
#reconstruct.o : shr_kind_mod.o
