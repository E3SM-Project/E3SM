#------------------------------------------------------------------------
# Makefile for objects/investment
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = accelerator.o \
             investment_growth_calculator.o \
             investment_utils.o \
             levelized_cost_calculator.o \
             market_based_investment.o \
             output_growth_calculator.o \
             rate_logit_distributor.o \
             simple_expected_profit_calculator.o \
             get_distributed_investment_visitor.o \
             investable_counter_visitor.o \
             output_share_levelized_cost_calculator.o \
             set_share_weight_visitor.o

investment_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
