#ifndef OMEGA_IO_H
#define OMEGA_IO_H
//===-- base/IO.h - basic IO utilities --------------------------*- C++ -*-===//
//
/// \file
/// \brief Defines basic IO functions for use by IO streams
///
/// The routines here provide the interfaces to the parallel IO environment,
/// currently the SCORPIO parallel IO library.
/// These functions should generally be accessed through the IOStreams
/// interfaces and not used directly. A few parameters are defined via
/// the Omega input configuration:
/// \ConfigInput
/// # Basic parallel IO configuration
/// IO:
///    # Number of MPI tasks to use for IO
///    # Default is 1 for safety but should be set to a subset of the
///    #  total MPI tasks for efficiency, eg to map efficiently to the
///    #  underlying hardware (network interfaces, NVM, etc)
///    IOTasks:  1
///    # The stride in MPI tasks when the number of IO tasks is less
///    #  than the total number of MPI tasks
///    IOStride: 1
///    # When using parallel IO, the data must be rearranged to match
///    #  the IO task decomposition. Choices are box and subset. Box is
///    #  the default and generally preferred (ensures each IO task has
///    #  a contiguous chunk of data). Subset is available for exploring
///    #  the most efficient approach. See SCORPIO documentation for details.
///    IORearranger: box
///    # The IO supports a number of file formats. We specify a default
///    #  here, but this value can be overridden on a file-by-file basis
///    #  through the streams interface. Choices include all the various
///    #  netCDF file formats as well as the ADIOS format.
///    IODefaultFormat: NetCDF4
/// \EndConfigInput
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "mpi.h"
#include "pio.h"

#include <map>
#include <string>

namespace OMEGA {

// Define a number of enum classes to enumerate various options
/// Choice of parallel IO rearranger algorithm
enum IORearranger {
   IORearrBox     = PIO_REARR_BOX,    ///< box rearranger (default)
   IORearrSubset  = PIO_REARR_SUBSET, ///< subset rearranger
   IORearrDefault = PIO_REARR_BOX,    ///< default value (Box)
   IORearrUnknown = -1                ///< unknown or undefined
};

/// Supported file formats
enum IOFileFmt {
   IOFmtNetCDF3  = PIO_IOTYPE_NETCDF,   ///< NetCDF3 classic format
   IOFmtPnetCDF  = PIO_IOTYPE_PNETCDF,  ///< Parallel NetCDF3
   IOFmtNetCDF4c = PIO_IOTYPE_NETCDF4C, ///< NetCDF4 (HDF5-compatible) cmpressed
   IOFmtNetCDF4p = PIO_IOTYPE_NETCDF4P, ///< NetCDF4 (HDF5-compatible) parallel
   IOFmtNetCDF4  = PIO_IOTYPE_NETCDF4P, ///< NetCDF4 (HDF5-compatible) parallel
   IOFmtHDF5     = PIO_IOTYPE_HDF5,     ///< native HDF5 format
   IOFmtADIOS    = PIO_IOTYPE_ADIOS,    ///< ADIOS format
   IOFmtUnknown  = -1,                  ///< Unknown or undefined
   IOFmtDefault  = PIO_IOTYPE_NETCDF4C, ///< NetCDF4 is default
};

/// File operations
enum IOMode {
   IOModeUnknown, /// Unknown or undefined
   IOModeRead,    /// Read  (input)
   IOModeWrite,   /// Write (output)
};

/// Behavior (for output files) when a file already exists
enum class IOIfExists {
   Fail,    /// Fail with an error
   Replace, /// Replace the file
   Append,  /// Append to the existing file
};

/// Floating point precision to allow reduced precision to save space
enum class IOPrecision {
   Double, /// Maintain full double precision (64-bit) for reals
   Single, /// Reduce all floating point variables to 32-bit reals
};

/// Data types for PIO corresponding to Omega types
enum IODataType {
   IOTypeI4      = PIO_INT,    /// 32-bit integer
   IOTypeI8      = PIO_INT64,  /// 64-bit integer
   IOTypeR4      = PIO_REAL,   /// 32-bit real
   IOTypeR8      = PIO_DOUBLE, /// 64-bit real
   IOTypeChar    = PIO_CHAR,   /// Character/string
   IOTypeLogical = PIO_INT     /// Logicals are converted to ints for IO
};

/// The IO system id, defined on IO initialization and used by all
/// IO functions
extern int IOSysID;

/// The default file format set on initialization. This value will
/// be assumed if not overridden during file open.
extern IOFileFmt IODefaultFileFmt;

// Utilities

/// Converts string choice for PIO rearranger to an enum
IORearranger
IORearrFromString(const std::string &Rearr ///< [in] choice of IO rearranger
);
/// Converts string choice for File Format to an enum
IOFileFmt
IOFileFmtFromString(const std::string &Format ///< [in] choice of IO file format
);
/// Converts string choice for IO read/write mode to an enum
IOMode IOModeFromString(
    const std::string &Mode ///< [in] choice of IO mode (read/write)
);
/// Converts string choice for existence behavior to an enum
IOIfExists IOIfExistsFromString(
    const std::string &IfExists ///< [in] choice of behavior on file existence
);
/// Converts string choice for floating point precision to an enum
IOPrecision IOPrecisionFromString(
    const std::string &Precision ///< [in] choice of floating point precision
);

// Methods

/// Initializes the IO system based on configuration inputs and
/// default MPI communicator
int IOInit(const MPI_Comm &InComm ///< [in] MPI communicator to use
);

/// This routine opens a file for reading or writing, depending on the
/// Mode argument. The filename with full path must be supplied and
/// a FileID is returned to be used by other IO functions.
/// The format of the file is assumed to be the default defined on init
/// but can be optionally changed through this open function.
/// For files to be written, optional arguments govern the behavior to be
/// used if the file already exists, and the precision of any floating point
/// variables. Returns an error code.
int IOFileOpen(
    int &FileID,                 ///< [out] returned fileID for this file
    const std::string &Filename, ///< [in] name (incl path) of file to open
    IOMode Mode,                 ///< [in] mode (read or write)
    IOFileFmt Format      = IOFmtDefault,     ///< [in] (optional) file format
    IOIfExists IfExists   = IOIfExists::Fail, ///< [in] behavior if file exists
    IOPrecision Precision = IOPrecision::Double ///< [in] precision of floats
);

/// Closes an open file using the fileID, returns an error code
int IOFileClose(int &FileID /// [in] ID of the file to be closed
);

/// Retrieves a dimension length from an input file, given the name
/// of the dimension. Returns the length if exists, but returns a negative
/// value if a dimension of that length is not found in the file.
int IOGetDimLength(int FileID, ///< [in] ID of the file containing dim
                   const std::string &DimName ///< [in] name of dimension
);

/// Creates a PIO decomposition description to describe the layout of
/// a distributed array of given type.
int IOCreateDecomp(
    int &DecompID,                      ///< [out] ID for the new decomposition
    IODataType VarType,                 ///< [in] data type of array
    int NDims,                          ///< [in] number of array dimensions
    const std::vector<int> &DimLengths, ///< [in] global dimension lengths
    int Size,                           ///< [in] local size of array
    const std::vector<int> &GlobalIndx, ///< [in] global indx for each loc indx
    IORearranger Rearr                  ///< [in] rearranger method to use
);

/// Removes a PIO decomposition to free memory.
int IODestroyDecomp(int &DecompID ///< [inout] ID for decomp to remove
);

/// Reads a distributed array. This overloaded interface is for integer
/// arrays. Also, all arrays are assumed to be in contiguous storage so
/// the arrays of any dimension are treated as a 1-d array with the full
/// local size.
int IOReadArray(int *Array,                 ///< [out] array to be read
                int Size,                   ///< [in] local size of array
                const std::string &VarName, ///< [in] name of variable to read
                int FileID,  ///< [in] ID of open file to read from
                int DecompID ///< [in] decomposition ID for this var
);

/// Reads a distributed array. This overloaded interface is for 32-bit real
/// arrays. Also, all arrays are assumed to be in contiguous storage so
/// the arrays of any dimension are treated as a 1-d array with the full
/// local size.
int IOReadArray(float *Array,               ///< [out] array to be read
                int Size,                   ///< [in] local size of array
                const std::string &VarName, ///< [in] name of variable to read
                int FileID,  ///< [in] ID of open file to read from
                int DecompID ///< [in] decomposition ID for this var
);

/// Reads a distributed array. This overloaded interface is for 64-bit real
/// arrays. Also, all arrays are assumed to be in contiguous storage so
/// the arrays of any dimension are treated as a 1-d array with the full
/// local size.
int IOReadArray(double *Array,              ///< [out] array to be read
                int Size,                   ///< [in] local size of array
                const std::string &VarName, ///< [in] name of variable to read
                int FileID,  ///< [in] ID of open file to read from
                int DecompID ///< [in] decomposition ID for this var
);

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_IO_H
