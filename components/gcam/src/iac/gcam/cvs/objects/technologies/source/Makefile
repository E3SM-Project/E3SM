#------------------------------------------------------------------------
# Makefile for objects/technologies
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = base_technology.o \
             building_generic_dmd_technology.o \
             cal_data_output.o \
             cal_data_output_percap.o \
             capture_component_factory.o \
             default_technology.o \
             expenditure.o \
             fixed_production_state.o \
             ag_production_technology.o \
             global_technology_database.o \
             intermittent_technology.o \
             internal_gains.o \
             marginal_profit_calculator.o \
             non_energy_use_capture_component.o \
             nuke_fuel_technology.o \
             output_factory.o \
             phased_shutdown_decider.o \
             power_plant_capture_component.o \
             primary_output.o \
             production_state_factory.o \
             production_technology.o \
             profit_shutdown_decider.o \
             residue_biomass_output.o \
             retired_production_state.o \
             secondary_output.o \
             sgm_output.o \
             shutdown_decider_factory.o \
             solar_technology.o \
             standard_capture_component.o \
             standard_technical_change_calc.o \
             technology.o \
             technology_type.o \
             tran_technology.o \
             unmanaged_land_technology.o \
             variable_production_state.o \
             vintage_production_state.o \
             wind_technology.o \
             technology_container.o \
             res_secondary_output.o \
             stub_technology_container.o \
             s_curve_shutdown_decider.o \
             empty_technology.o \
             generic_output.o \
             input_factory.o

technologies_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
