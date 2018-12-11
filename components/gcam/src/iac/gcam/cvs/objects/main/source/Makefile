#------------------------------------------------------------------------
# Makefile for objects/main
# William.Gustafson@pnl.gov; 8-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = main.o

main_dir: ${OBJS} gcam

gcam : 
	$(RANLIB) ${PATHOFFSET}/build/linux/libgcamlib.a
	$(LD) -o gcam $(LDFLAGS) main.o $(LIB) ${PATHOFFSET}/build/linux/libgcamlib.a

clean:
	rm *.o
