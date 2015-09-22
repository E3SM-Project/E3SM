if TRILINOS

##lib_LIBRARIES = libglimmertrilinos.a
lib_LTLIBRARIES = libglimmertrilinos.la

variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

BUILT_SOURCES = 
CLEANFILES = $(BUILT_SOURCES) 

libglimmertrilinos_la_SOURCES = \
                                trilinosLinearSolver.cpp \
                                matrixInterface.cpp \
                                trilinosNoxSolver.cpp \
                                trilinosModelEvaluator.cpp

libglimmertrilinos_la_DEPENDENCIES = matrixInterface.hpp \
                                 trilinosModelEvaluator.hpp



#ALL_SRC= $(libglimmertrilinos_la_SOURCES)

AM_FCFLAGS = $(CXXFLAGS) $(NETCDF_CPPFLAGS) $(HDF5_CPPFLAGS)
AM_CXXFLAGS = $(CXXFLAGS)

LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) 

AM_CFLAGS=$(CFORTRAN_CFLAGS)

endif TRILINOS
