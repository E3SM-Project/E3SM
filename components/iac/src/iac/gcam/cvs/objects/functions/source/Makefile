#------------------------------------------------------------------------
# Makefile for objects/functions
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = ademand_function.o \
             aproduction_function.o \
             building_demand_input.o \
             ces_production_function.o \
             demand_input.o \
             efficiency.o \
             energy_input.o \
             function_manager.o \
             function_utils.o \
             government_demand_function.o \
             household_demand_function.o \
             input_subsidy.o \
             input_tax.o \
             intensity.o \
             investment_demand_function.o \
             leontief_production_function.o \
             minicam_input.o \
             minicam_leontief_production_function.o \
             minicam_price_elasticity_function.o \
             non_energy_input.o \
             production_input.o \
             renewable_input.o \
             sgm_input.o \
             trade_demand_function.o \
             utility_demand_function.o \
             trade_input.o \
             node_input.o \
             nested_ces_production_function.o \
             logit_production_function.o \
             input_OM_fixed.o \
             input_OM_var.o \
             input_capital.o \
             leaf_input_finder.o

functions_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
