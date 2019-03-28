#------------------------------------------------------------------------
# Makefile for objects/target_finder
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = bisecter.o \
             concentration_target.o \
             emissions_stabalization_target.o \
             forcing_target.o \
             policy_target_runner.o \
             simple_policy_target_runner.o \
             target_factory.o \
             secanter.o \
             kyoto_forcing_target.o \
             temperature_target.o

target_finder_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
