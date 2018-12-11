#------------------------------------------------------------------------
# Makefile for objects/util/curves
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = curve.o \
             data_point.o \
             explicit_point_set.o \
             point_set.o \
             point_set_curve.o \
             xy_data_point.o

util_curves_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
