#------------------------------------------------------------------------
# Makefile for objects/marketplace
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = calibration_market.o \
             demand_market.o \
             inverse_calibration_market.o \
             market.o \
             market_locator.o \
             market_subsidy.o \
             market_tax.o \
             marketplace.o \
             normal_market.o \
             price_market.o \
             cached_market.o \
             market_RES.o \
             trial_value_market.o

marketplace_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
