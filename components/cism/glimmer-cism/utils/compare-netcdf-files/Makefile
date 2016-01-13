# Set the environment variable NETCDF to the path of your NETCDF installation.

CXX=g++
CXXFLAGS= -O3 -m64
CXXLIBS = -I$(NETCDF)/include -L$(NETCDF)/lib -lnetcdf_c++

all:
	$(CXX) compare.cpp $(CXXLIBS) $(CXXFLAGS) -o compare

clean:
	rm compare

