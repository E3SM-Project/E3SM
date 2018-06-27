#------------------------------------------------------------------------
# Makefile for objects/resources
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = accumulated_grade.o \
             accumulated_post_grade.o \
             grade.o \
             renewable_subresource.o \
             resource.o \
             smooth_renewable_subresource.o \
             subresource.o \
             unlimited_resource.o \
             depleting_fixed_resource.o

resources_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
