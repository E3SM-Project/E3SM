#------------------------------------------------------------------------
# Makefile for objects/util/base
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = atom.o \
             atom_registry.o \
             configuration.o \
             input_finder.o \
             model_time.o \
             summary.o \
             supply_demand_curve.o \
             timer.o \
             calibrate_share_weight_visitor.o \
             calibrate_resource_visitor.o \
             interpolation_rule.o \
             interpolation_function_factory.o \
             fixed_interpolation_function.o \
             linear_interpolation_function.o \
             s_curve_interpolation_function.o \
             util.o

util_base_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
