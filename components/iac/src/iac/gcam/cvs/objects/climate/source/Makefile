#------------------------------------------------------------------------
# Makefile for objects/climate
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = magicc_model.o \
             ObjECTS_MAGICC_others.o \
             MAGICC_array.o \
             MAGICC_IO_helpers.o \
             ObjECTS_MAGICC.o

climate_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
