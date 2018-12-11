#------------------------------------------------------------------------
# Makefile for objects/agLU
# William.Gustafson@pnl.gov; 8-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = Ag2Calibrate.o \
             Ag2CH4N2O.o \
             Ag2control.o \
             Ag2Demand.o \
             Ag2Emiss.o \
             Ag2GetData.o \
             Ag2init.o \
             Ag2Link.o \
             Ag2model.o \
             Ag2Output.o \
             Ag2Supply.o

MODULES    = Ag2Global8.o

$(OBJS): $(MODULES)

agLU_dir: ${MODULES} ${OBJS}

clean:
	rm *.o *.mod
