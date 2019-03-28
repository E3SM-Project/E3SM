#------------------------------------------------------------------------
# Makefile for objects/solution/solvers
# William.Gustafson@pnl.gov; 5-Sep-2008
#------------------------------------------------------------------------

#PATHOFFSET = path to objects directory
PATHOFFSET = ../../..
include ${PATHOFFSET}/build/linux/configure.gcam

OBJS       = bisect_all.o \
             bisect_one.o \
             bisect_policy.o \
             bisect_policy_nr_solver.o \
             bisection_nr_solver.o \
             log_newton_raphson.o \
             log_newton_raphson_sd.o \
             solver_factory.o \
             user_configurable_solver.o \
             solver_component_factory.o \
             solver_component.o

solution_solvers_dir: ${OBJS}

-include $(DEPS)

clean:
	rm *.o
