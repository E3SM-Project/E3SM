#------------------------------------------------------------------------
# Makefile for objects/emissions
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = acomplex_emissions.o \
             aemissions_coef.o \
             aggr_emissions_coef.o \
             aghg.o \
             by_product.o \
             co2_emissions.o \
             emissions_driver_factory.o \
             emissions_summer.o \
             generic_emissions.o \
             ghg_factory.o \
             ghg_mac.o \
             input_driver.o \
             input_emissions_coef.o \
             input_output_driver.o \
             luc_emissions_summer.o \
             output_driver.o \
             read_emissions_coef.o \
             so2_emissions.o \
             total_sector_emissions.o

emissions_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
