#===============================================================================
# $Id: Makefile 46158 2013-04-19 18:41:34Z mlevy@ucar.edu $
# $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_141106/check_maps/src/Makefile $
#
# This is the makefile for the RegridWeightGen application in ESMF
# Altered for use with CESM's gen_atm_ocn_maps.sh tool
#===============================================================================

################################################################################
### Finding and including esmf.mk                                         ######
### NOTE: On Cheyenne, ESMF is only available with the intel compilers    ######
###       And require running the following:                              ######
###           $ module load esmf_libs/7.0.0                               ######
###           $ module load esmf-7.0.0-ncdfio-uni-O                       ######
###       The executable should run on the cheyenne login node            ######
################################################################################

# Don't require ESMF to be loaded to run "make clean":
USE_ESMF=TRUE
ifeq ($(MAKECMDGOALS),clean)
  USE_ESMF=FALSE
endif

ifeq ($(USE_ESMF),TRUE)
  ifneq ($(origin ESMFMKFILE), environment)
    $(error Environment variable ESMFMKFILE was not set.)
  endif
  include $(ESMFMKFILE)
endif
TARGET = ESMF_RegridWeightGenCheck

# To compile with verbose output, run
# $ gmake VERBOSE=TRUE
#
# To get back to concise output, run
# $ gmake VERBOSE=FALSE
# or touch the source file prior to running gmake again
#
ifneq ($(strip $(VERBOSE)),)
  ifneq ($(VERBOSE),FALSE)
    # Technically any string other than VERBOSE=FALSE will compile verbose
    # ESMF_F90COMPILECPPFLAGS += -WF,-DVERBOSE=1
    ESMF_F90COMPILECPPFLAGS += -DVERBOSE
  endif
  $(shell rm -f $(TARGET).o)
endif


################################################################################
### Target / Linking ###########################################################

../ESMF_RegridWeightGenCheck: $(TARGET).o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) \
  $(ESMF_F90LINKRPATHS) -o $@ $(TARGET).o $(ESMF_F90ESMFLINKLIBS)

.PHONY: verbose clean
verbose:
	$(MAKE) VERBOSE=TRUE

clean:
	rm -f ../ESMF_RegridWeightGenCheck ../PET*.Log ../hostfile *.o

###############################################################################
### Compiler rules using ESMF_ variables supplied by esmf.mk ##################

.SUFFIXES: .F90 .o

.F90.o:
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) \
  $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $<

