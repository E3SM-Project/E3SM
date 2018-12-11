#------------------------------------------------------------------------
# Makefile for objects/ccsmcpl
# tcraig 19-Apr-2010
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       =  CCSM_GCAM_interface.o CCSM_GCAM_interface_wrapper.o

MODULES    = gcam_comp_mod.o 

ccsmcpl_dir:   ${OBJS} ${MODULES}


clean:
	rm *.o *.mod
