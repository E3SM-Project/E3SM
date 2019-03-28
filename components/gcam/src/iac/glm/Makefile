include ./configure.glm

BUILDPATH = $(PWD)

SOURCES=glm.future.iesm.c ccsm_glm_interface.c glm_comp_mod.F90
LIBOBJECTS=$(SOURCES:.c=.o) glm_comp_mod.o
EXECUTABLE=glm
LDFLAGS += -lm -fpic -lnetcdf


all: $(LIBOBJECTS)

glmlib: $(LIBOBJECTS)

glm:
	${CXX} -c ${CXXFLAGS} ${ARCHFLAGS} ${INCLUDE} -DMAIN glm.future.iesm.c
	$(LD) $(LDFLAGS) glm.future.iesm.o -L$(LIB_NETCDF) -lnetcdf -o $@
clean:
	rm -f *.o $(EXECUTABLE) *.a *.mod *.test *.txt *latest

