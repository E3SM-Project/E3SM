#
# Defined externally in Makefile.conf
#
#   INCLUDES
#   LIBS
#   MPICC
#   MPIF90
#   COPTS
#   FOPTS
#   CFLAGS
#   FFLAGS
#   AWK
#   AR

RM=/bin/rm -f
FDEPENDS=fdepends.awk
DEPSUF = .d
MODSUF = .mod
CPPSUF = .f90


ifeq (,$(PIOARCH))
  PIOARCH=conf
endif
include Makefile.$(PIOARCH)
export PIOARCH

SRCS_C = 

SRCS_F90 =  pio.F90 \
            pio_kinds.F90  \
            nf_mod.F90     \
            ionf_mod.F90 \
            pio_types.F90  \
            piolib_mod.F90 \
	    pio_mpi_utils.F90 \
            pio_nf_utils.F90 \
	    pio_utils.F90 \
            pio_quicksort.F90 
           

TEMPLATES_F90 = pionfatt_mod.F90.in \
	        pionfread_mod.F90.in \
                pionfwrite_mod.F90.in \
                pionfput_mod.F90.in \
                pionfget_mod.F90.in \
	        alloc_mod.F90.in \
                box_rearrange.F90.in \
                rearrange.F90.in \
                pio_support.F90.in \
	        mct_rearrange.F90.in \
	        iompi_mod.F90.in  \
	        piodarray.F90.in \
	        pio_spmd_utils.F90.in

TEMPSRCF90 = $(TEMPLATES_F90:.in=)
SRCS_F90 += $(TEMPSRCF90)

OBJS=  $(SRCS_C:.c=.o) \
       $(SRCS_F90:.F90=.o)


MODFILES := $(SRCS_F90:.F90=$(MODSUF))

# File foo.DEPSUF will contain make dependency for foo.F90
DEPENDS := $(SRCS_F90:.F90=$(DEPSUF))
PERL = /usr/bin/perl

LIB= libpio.a

LEVEL=0


all: depends
	@$(MAKE) LEVEL=1 $(LIB)



depends: $(DEPENDS)
	@echo "Done updating dependencies"


#
# Suffix rules
#

.SUFFIXES:
.SUFFIXES: .o .c .F90 $(DEPSUF)


# For each file foo.F90, make the depency file foo.d

%$(DEPSUF): %.F90
	@echo 'Making dependencies for' $< '-->' $@
	@$(AWK) -f $(FDEPENDS) -v NAME=$(basename $<) -v SUF=$(suffix $<) $< > $@


ifeq ($(EXPLICIT_CPP),yes)
SRCS_CPP= $(SRCS_F90:.F90=$(CPPSUF))
.F90.o:
	@if [ -w $*.f90 ] ; then echo "ERROR: file $*.f90 is writable - the .f90 suffix is reserved for temporary cpp output" ; exit 1; fi
	$(RM) $*.f90
	$(CPP) $(CPPFLAGS) $(CFLAGS) $(COPTS) $(INCLUDES) -o $*.f90 $*.F90 
	chmod a-w $*.f90
	$(MPIF90) -c $(FFLAGS) $(FOPTS) $(INCLUDES) $*.f90
else
SRCS_CPP=
.F90.o:
	$(MPIF90) -c $(FFLAGS) $(FOPTS) $(INCLUDES) $*.F90
endif

.c.o:
	$(MPICC) -c $(CFLAGS) $(COPTS) $(INCLUDES) $*.c


$(TEMPSRCF90): $(TEMPLATE_F90)
	$(PERL) genf90.pl $@.in > $*.F90



$(LIB): $(OBJS)
	$(RM) $@
	$(AR) $(ARFLAGS) $@ $(OBJS)


predist: predistclean $(TEMPSRCF90)
     
clean:
	$(RM) $(LIB) $(OBJS) $(MODFILES) $(DEPENDS) $(SRCS_CPP)

predistclean: clean
	$(RM) $(TEMPSRCF90)

#
# Automatically generated module dependencies
#

ifneq (0,$(LEVEL))
  include $(DEPENDS)
endif

alloc_mod.o : alloc_mod.F90
box_rearrange.o : box_rearrange.F90
iompi_mod.o : iompi_mod.F90
mct_rearrange.o : mct_rearrange.F90
piodarray.o : piodarray.F90
pionfatt_mod.o : pionfatt_mod.F90
pionfget_mod.o : pionfget_mod.F90
pionfput_mod.o : pionfput_mod.F90
pionfread_mod.o : pionfread_mod.F90
pionfwrite_mod.o : pionfwrite_mod.F90
pio_support.o : pio_support.F90
rearrange.o : rearrange.F90

alloc_mod.F90 : alloc_mod.F90.in
box_rearrange.F90 : box_rearrange.F90.in
iompi_mod.F90 : iompi_mod.F90.in
mct_rearrange.F90 : mct_rearrange.F90.in
piodarray.F90 : piodarray.F90.in
pionfatt_mod.F90 : pionfatt_mod.F90.in
pionfget_mod.F90 : pionfget_mod.F90.in
pionfput_mod.F90 : pionfput_mod.F90.in
pionfread_mod.F90 : pionfread_mod.F90.in
pionfwrite_mod.F90 : pionfwrite_mod.F90.in
pio_support.F90 : pio_support.F90.in
rearrange.F90 : rearrange.F90.in
