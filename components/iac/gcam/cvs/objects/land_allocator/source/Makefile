#------------------------------------------------------------------------
# Makefile for objects/land_allocator_finder
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = aland_allocator_item.o \
             land_leaf.o \
             land_node.o \
             land_use_history.o \
             land_allocator.o \
             carbon_land_leaf.o \
             set_carbon_density.o \
             unmanaged_land_leaf.o

land_allocator_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
