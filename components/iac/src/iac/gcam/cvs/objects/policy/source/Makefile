#------------------------------------------------------------------------
# Makefile for objects/policy
# William.Gustafson@pnl.gov; 8-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = policy_ghg.o \
             policy_portfolio_standard.o

policy_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
