#------------------------------------------------------------------------
# Makefile for objects/demographics
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = age_cohort.o \
             demographic.o \
             female.o \
             gender.o \
             male.o \
             population.o \
             population_mini_cam.o \
             population_sgm_fixed.o \
             population_sgm_rate.o

demographics_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
