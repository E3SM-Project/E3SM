#ifndef SIMPLE_NETCDF_HPP
#define SIMPLE_NETCDF_HPP
// Provide a simple netCDF interface for RRTMGP that will read data in as YAKL
// arrays. This replaces the simple netCDF interface that is included with YAKL
// that requires netcdf-cxx4. We use it here instead of SCORPIO because we use
// the I/O routines from RRTMGP++ mostly unchanged, and they used the YAKL I/O
// for convenience. This small library mimics the functionality of the YAKL I/O
// and provides a nearly 100% compatible API, but uses netCDF-C directly rather
// than netcdf-cxx4. This library is provided specifically to work with
// existing RRTMGP I/O routines; all other codes developed for SCREAM should
// use SCORPIO for I/O.
#include <netcdf.h>
#include "YAKL.h"

using namespace yakl;
namespace simple_netcdf {

    class SimpleNetCDF {

        protected:

            int ncid;

        public:

            // Constructor
            SimpleNetCDF() {};

            // Destructor
            ~SimpleNetCDF() {
                //close();
            };

            void close() {
                handle_error(nc_close(ncid));
            }

            void create(std::string filename, int mode=NC_CLOBBER) {
                handle_error(nc_create(filename.c_str(), mode, &ncid));
            };

            void open(std::string filename, int mode=NC_NOWRITE) {
                handle_error(nc_open(filename.c_str(), mode, &ncid));
            };

            void open(char *filename) {
                handle_error(nc_open(filename, NC_NOWRITE, &ncid));
            }

            // NetCDF routines return an integer error code. Define a function
            // here to abort program execution and throw an error code if we
            // encounter a non-zero NetCDF return code. We will wrap our
            // NetCDF calls with this function to handle these errors in a
            // consistent way
            void handle_error(int err) {
                if (err) {
                    std::cout << "ERROR: " << nc_strerror(err) << std::endl;
                    abort();
                }
            }

            void handle_error(int err, const char *file, int line) {
                if (err) {
                    std::cout << "ERROR: " << nc_strerror(err) << " at line " << line << " in " << file << std::endl;
                    abort();
                }
            }

            // Read a netCDF array to a YAKL array
            template <class T, int rank, int myMem, int myStyle> void read(Array<T,rank,myMem,myStyle> &arr, std::string varName) {

                // Get variable ID
                int varid;
                handle_error(nc_inq_varid(ncid, varName.c_str(), &varid), __FILE__, __LINE__);

                // Get variable dimension sizes
                int ndims;
                int dimids[NC_MAX_VAR_DIMS];
                nc_type vtype;
                handle_error(nc_inq_var(ncid, varid, NULL, &vtype, &ndims, dimids, NULL), __FILE__, __LINE__);
                std::vector<int> dimSizes(ndims);
                size_t dimsize;
                for (int i = 0; i < ndims; i++) {
                    handle_error(nc_inq_dimlen(ncid, dimids[i], &dimsize), __FILE__, __LINE__); 
                    dimSizes[i] = dimsize;
                }

                // If style is fortran, we need to reverse array dims
                if (myStyle == styleFortran) {
                    std::reverse(dimSizes.begin(), dimSizes.end());
                }

                // Allocate (or reshape) the yakl array
                arr = Array<T,rank,myMem,myStyle>(varName.c_str(),dimSizes);

                // Read variable data
                if (myMem == memDevice) {
                    auto arrHost = arr.createHostCopy();
                    if (std::is_same<T,bool>::value) {
                        // Create boolean array from integer arrays
                        Array<int,rank,memHost,myStyle> tmp("tmp",dimSizes);
                        handle_error(nc_get_var(ncid, varid, tmp.data()), __FILE__, __LINE__);
                        for (int i=0; i < arr.totElems(); i++) { arrHost.myData[i] = tmp.myData[i] == 1; }
                    } else {
                        // Need to be careful with floats; nc_get_var is overloaded on type, but we need
                        // to make sure we read floats from file with the float procedure, and doubles
                        // with that for doubles. The danger is if the user passes a yakl array here
                        // with type double, but tries to read type float from file.
                        // TODO: why does the YAKL implementation for this work fine, but this version
                        // calling nc_get_var directly does not?
                        if (vtype == NC_FLOAT) {
                            Array<float,rank,memHost,myStyle> tmp("tmp",dimSizes);
                            handle_error(nc_get_var(ncid, varid, tmp.data()), __FILE__, __LINE__);
                            for (size_t i=0; i < arr.totElems(); i++) { arrHost.myData[i] = tmp.myData[i]; }
                        } else {
                            handle_error(nc_get_var(ncid, varid, arrHost.data()), __FILE__, __LINE__);
                        }
                    }
                    arrHost.deep_copy_to(arr);
                } else {
                    if (std::is_same<T,bool>::value) {
                        // Create boolean array from integer arrays
                        Array<int,rank,memHost,myStyle> tmp("tmp",dimSizes);
                        handle_error(nc_get_var(ncid, varid, tmp.data()), __FILE__, __LINE__);
                        for (size_t i=0; i < arr.totElems(); i++) { arr.myData[i] = tmp.myData[i] == 1; }
                    } else {
                        if (vtype == NC_FLOAT) {
                            Array<float,rank,memHost,myStyle> tmp("tmp",dimSizes);
                            handle_error(nc_get_var(ncid, varid, tmp.data()), __FILE__, __LINE__);
                            for (size_t i=0; i < arr.totElems(); i++) { arr.myData[i] = tmp.myData[i]; }
                        } else {
                            handle_error(nc_get_var(ncid, varid, arr.data()), __FILE__, __LINE__);
                        }
                    }
                }

            }

            // Read a scalar type
            template <class T> void read(T &arr , std::string varName) {
                // Get variable ID
                int varid;
                handle_error(nc_inq_varid(ncid, varName.c_str(), &varid), __FILE__, __LINE__);

                // Read data
                handle_error(nc_get_var(ncid, varid, &arr), __FILE__, __LINE__);
            }

            // Check if variable exists in file
            bool varExists (std::string varName) {
                int varid;
                int ncerr = nc_inq_varid(ncid, varName.c_str(), &varid);
                if (ncerr == 0) {
                    return true;
                } else {
                    return false;
                }
            }

            bool dimExists (std::string dimName) {
                int dimid;
                int ncerr = nc_inq_dimid(ncid, dimName.c_str(), &dimid);
                if (ncerr == 0) {
                    return true;
                } else {
                    return false;
                }
            }

            size_t getDimSize(std::string dimName) {
                // Get dimension ID
                int dimid;
                handle_error(nc_inq_dimid(ncid, dimName.c_str(), &dimid));

                // Get dimension size
                size_t dimSize;
                handle_error(nc_inq_dimlen(ncid, dimid, &dimSize));

                return dimSize;
            }

            void addDim(std::string dimName, int dimSize, int *dimid) {
                // Put file into define mode
                int ncerr = nc_redef(ncid);
                if ((ncerr != NC_NOERR) and (ncerr != NC_EINDEFINE)) {
                    handle_error(ncerr, __FILE__, __LINE__);
                }

                // Define dimension
                handle_error(nc_def_dim(ncid, dimName.c_str(), dimSize, dimid), __FILE__, __LINE__);

                // End define mode
                handle_error(nc_enddef(ncid), __FILE__, __LINE__);
            }

            void addVar(std::string varName, nc_type varType, int ndims, int dimids[], int *varid) {
                // Put file into define mode
                int ncerr = nc_redef(ncid);
                if ((ncerr != NC_NOERR) and (ncerr != NC_EINDEFINE)) {
                    handle_error(ncerr, __FILE__, __LINE__);
                }

                // Define variable
                handle_error(nc_def_var(ncid, varName.c_str(), varType, ndims, dimids, varid), __FILE__, __LINE__);

                // End define mode
                handle_error(nc_enddef(ncid), __FILE__, __LINE__);
            }

            template <class T> void putVar(T const &arr, std::string varName) {
                // Make sure file is not in define mode
                int ncerr = nc_enddef(ncid);
                if ((ncerr != NC_NOERR) and (ncerr != NC_ENOTINDEFINE)) {
                    handle_error(ncerr, __FILE__, __LINE__);
                }

                // Get variable Id
                int varid;
                handle_error(nc_inq_varid(ncid, varName.c_str(), &varid), __FILE__, __LINE__);

                // Write variable data
                handle_error(nc_put_var(ncid, varid, arr), __FILE__, __LINE__);
            }

            template <class T, int rank, int myMem, int myStyle> 
            void write(Array<T,rank,myMem,myStyle> const &arr, std::string varName, std::vector<std::string> dimNames) {

                // Make sure length of dimension names is equal to rank of array
                if (rank != dimNames.size()) { yakl_throw("dimNames.size() != Array rank"); }

                // Get dimension sizes
                // Define dimensions if they do not exist and get dimension IDs
                int dimids[rank];
                size_t dimSize;
                size_t idim;
                for (size_t i = 0; i < dimNames.size(); i++) {
                    // If style is Fortran, dimension ordering is reversed
                    if (myStyle == styleC) {
                        idim = i;
                    } else {
                        idim = rank - 1 - i;
                    }
                    int ncerr = nc_inq_dimid(ncid, dimNames[i].c_str(), &dimids[idim]);
                    if (ncerr == NC_NOERR) {
                        // check that size is correct
                        handle_error(nc_inq_dimlen(ncid, dimids[idim], &dimSize), __FILE__, __LINE__);
                        if (dimSize != arr.dimension[i]) {
                            yakl_throw("dimSize != arr.dimension[i]");
                        }
                    } else {
                        addDim(dimNames[i], arr.dimension[i], &dimids[idim]);
                    }
                }

                // Add variable if it does not exist
                if (!varExists(varName)) {
                    int varid;
                    addVar(varName, getType<T>(), rank, dimids, &varid);
                }

                // Write data to file
                if (myMem == memDevice) {
                    putVar(arr.createHostCopy().data(), varName);
                } else {
                    putVar(arr.data(), varName);
                }
            }

            template <class T> void write(T arr, std::string varName) {
                // If variable does not exist, try to add it
                if (!varExists(varName)) {
                    int dimids[1] = {0};
                    int varid;
                    addVar(varName, getType<T>(), 0, dimids, &varid);
                }
                // Write to file
                putVar(&arr, varName);
            }

            // Determine nc_type corresponding to intrinsic type
            template <class T> nc_type getType() const {
                     if ( std::is_same<T,          char>::value ) { return NC_CHAR;   }
                else if ( std::is_same<T,unsigned  char>::value ) { return NC_UBYTE;  }
                else if ( std::is_same<T,         short>::value ) { return NC_SHORT;  }
                else if ( std::is_same<T,unsigned short>::value ) { return NC_USHORT; }
                else if ( std::is_same<T,           int>::value ) { return NC_INT;    }
                else if ( std::is_same<T,unsigned   int>::value ) { return NC_UINT;   }
                else if ( std::is_same<T,          long>::value ) { return NC_INT64;  }
                else if ( std::is_same<T,unsigned  long>::value ) { return NC_UINT64; }
                else if ( std::is_same<T,         float>::value ) { return NC_FLOAT;  }
                else if ( std::is_same<T,        double>::value ) { return NC_DOUBLE; }
                else if ( std::is_same<T,std::string   >::value ) { return NC_STRING; }
                else { yakl_throw("Invalid type"); }
                return -1;
            }

    };  // class SimpleNetCDF

} // namespace simple_netcdf
#endif
