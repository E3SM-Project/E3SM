#------------------------------------------------------------------------
# Makefile for objects/reporting
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = batch_csv_outputter.o \
             demand_components_table.o \
             govt_results.o \
             graph_printer.o \
             indirect_emissions_calculator.o \
             input_output_table.o \
             land_allocator_printer.o \
             rcp_emissions_visitor.o \
             sector_report.o \
             sector_results.o \
             sgm_gen_table.o \
             social_accounting_matrix.o \
             storage_table.o \
             energy_balance_table.o \
             get_glm_data.o \
             xml_db_outputter.o

reporting_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
