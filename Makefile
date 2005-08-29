
FC= pgf90

# -DCHECKS -DHANDLE_INFO
CFLAGS= -O3 


.SUFFIXES: .F90

.F90.o:
	$(FC) $(FFLAGS) -c -o $@ $<



LIBOBJS= mpi.o send.o recv.o collective.o req.o list.o fort.o handles.o

LIB=libmpi.a

$(LIB)(%.o): %.o
	$(AR) $(ARFLAGS) $(LIB) $%

LIB_MEMBERS= $(foreach file, $(LIBOBJS), $(LIB)($(file)))

$(LIB): $(LIB_MEMBERS)
	@echo "Done builing $@"



ctest: $(LIB) ctest.c
	cc -o $@ ctest.c -L. -lmpi

ftest: $(LIB) ftest.F90
	$(FC) -o $@ ftest.F90 -L. -lmpi


clean:
	/bin/rm -f *.o ctest ftest $(LIB)

