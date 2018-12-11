#------------------------------------------------------------------------
# Makefile for objects/sectors
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = backup_calculator_factory.o \
             building_final_demand.o \
             cal_quantity_tabulator.o \
             capacity_limit_backup_calculator.o \
             CSP_backup_calculator.o \
             energy_final_demand.o \
             export_sector.o \
             factor_supply.o \
             final_demand_sector.o \
             ag_supply_sector.o \
             ag_supply_subsector.o \
             more_sector_info.o \
             production_sector.o \
             sector.o \
             sector_utils.o \
             subsector.o \
             subsector_add_techcosts.o \
             supply_sector.o \
             tran_sector.o \
             tran_subsector.o \
             wind_backup_calculator.o

sectors_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
