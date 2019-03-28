#------------------------------------------------------------------------
# Makefile for objects/util/database
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = output_helper.o

util_database_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
