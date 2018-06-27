#------------------------------------------------------------------------
# Makefile for objects/consumers
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = consumer.o \
             govt_consumer.o \
             household_consumer.o \
             invest_consumer.o \
             trade_consumer.o \
             calc_capital_good_price_visitor.o

consumers_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
