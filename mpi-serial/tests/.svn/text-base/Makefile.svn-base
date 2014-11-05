###############################

#
# test programs Makefile
#

# Parent dir Makefile.conf has all necessary vars
include ../Makefile.conf

TINC    = -I..
LDFLAGS = -L..
LIBS    = -l$(MODULE)

runtests: all
	./ctest
	@echo
	./ftest
	@echo

all: ctest ftest ctest2 ftest2

ctest: ctest.c
	$(CC) $(DEFS) $(TINC) $(ALLCFLAGS) -o $@ ctest.c $(LDFLAGS) $(LIBS)

ftest: ftest.F90
	$(FC) $(DEFS) $(TINC) $(MYF90FLAGS) -o $@ ftest.F90 $(LDFLAGS) $(LIBS)

ctest2: ctest_old.c
	$(CC) $(DEFS) $(TINC) $(ALLCFLAGS) -o $@ ctest_old.c $(LDFLAGS) $(LIBS)

ftest2: ftest_old.F90
	$(FC) $(DEFS) $(TINC) $(MYF90FLAGS) -o $@ ftest_old.F90 $(LDFLAGS) $(LIBS)

stest: stest.F90 stest2.o
	$(FC) $(DEFS) $(TINC) $(MYF90FLAGS) -o $@ stest.F90 stest2.o $(LDFLAGS) $(LIBS)


clean:
	rm -f ctest ftest ctest2 ftest2
	rm -f *.o
