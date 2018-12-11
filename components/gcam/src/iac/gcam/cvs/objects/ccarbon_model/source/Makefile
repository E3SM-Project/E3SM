#------------------------------------------------------------------------
# Makefile for objects/ccarbon_model
# William.Gustafson@pnl.gov; 8-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = carbon_model_utils.o \
             asimple_carbon_calc.o \
             land_carbon_densities.o

ccarbon_model_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
