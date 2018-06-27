#------------------------------------------------------------------------
# Makefile for objects/util/logger
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = logger.o \
             logger_factory.o \
             plain_text_logger.o \
             xml_logger.o

util_logger_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
