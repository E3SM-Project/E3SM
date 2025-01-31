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
namespace IO {

// Define a number of parameters and enum classes to enumerate various options
/// ID for global metadata (ie metadata not associated with a variable)
constexpr int GlobalID = PIO_GLOBAL;

/// Length for unlimited dimensions
constexpr int Unlimited = PIO_UNLIMITED;

/// Choice of parallel IO rearranger algorithm
enum Rearranger {
   RearrBox     = PIO_REARR_BOX,    ///< box rearranger (default)
   RearrSubset  = PIO_REARR_SUBSET, ///< subset rearranger
   RearrDefault = PIO_REARR_BOX,    ///< default value (Box)
   RearrUnknown = -1                ///< unknown or undefined
};

/// Supported file formats
enum FileFmt {
   FmtNetCDF3  = PIO_IOTYPE_NETCDF,   ///< NetCDF3 classic format
   FmtPnetCDF  = PIO_IOTYPE_PNETCDF,  ///< Parallel NetCDF3
   FmtNetCDF4c = PIO_IOTYPE_NETCDF4C, ///< NetCDF4 (HDF5-compatible) cmpressed
   FmtNetCDF4p = PIO_IOTYPE_NETCDF4P, ///< NetCDF4 (HDF5-compatible) parallel
   FmtNetCDF4  = PIO_IOTYPE_NETCDF4P, ///< NetCDF4 (HDF5-compatible) parallel
   FmtHDF5     = PIO_IOTYPE_HDF5,     ///< native HDF5 format
   FmtADIOS    = PIO_IOTYPE_ADIOS,    ///< ADIOS format
   FmtUnknown  = -1,                  ///< Unknown or undefined
   FmtDefault  = PIO_IOTYPE_NETCDF4C, ///< NetCDF4 is default
};

/// File operations
enum Mode {
   ModeUnknown, /// Unknown or undefined
   ModeRead,    /// Read  (input)
   ModeWrite,   /// Write (output)
};

/// Behavior (for output files) when a file already exists
enum class IfExists {
   Fail,    /// Fail with an error
   Replace, /// Replace the file
   Append,  /// Append to the existing file
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
extern int SysID;

/// The default file format set on initialization. This value will
/// be assumed if not overridden during file open.
extern FileFmt DefaultFileFmt;

/// The default rearranger set on initialization. This value will
/// be assumed if not overridden during file open.
extern Rearranger DefaultRearr;

// Utilities

/// Converts string choice for PIO rearranger to an enum
Rearranger
RearrFromString(const std::string &Rearr ///< [in] choice of IO rearranger
);
/// Converts string choice for File Format to an enum
FileFmt
FileFmtFromString(const std::string &Format ///< [in] choice of IO file format
);
/// Converts string choice for IO read/write mode to an enum
Mode ModeFromString(
    const std::string &Mode ///< [in] choice of IO mode (read/write)
);
/// Converts string choice for existence behavior to an enum
IfExists IfExistsFromString(
    const std::string &IfExists ///< [in] choice of behavior on file existence
);

// Methods

/// Initializes the IO system based on configuration inputs and
/// default MPI communicator
int init(const MPI_Comm &InComm ///< [in] MPI communicator to use
);

/// This routine opens a file for reading or writing, depending on the
/// Mode argument. The filename with full path must be supplied and
/// a FileID is returned to be used by other IO functions.
/// The format of the file is assumed to be the default defined on init
/// but can be optionally changed through this open function.
/// For files to be written, optional arguments govern the behavior to be
/// used if the file already exists, and the precision of any floating point
/// variables. Returns an error code.
int openFile(
    int &FileID,                    ///< [out] returned fileID for this file
    const std::string &Filename,    ///< [in] name (incl path) of file to open
    Mode Mode,                      ///< [in] mode (read or write)
    FileFmt Format    = FmtDefault, ///< [in] (optional) file format
    IfExists IfExists = IfExists::Fail ///< [in] behavior if file exists
);

/// Closes an open file using the fileID, returns an error code
int closeFile(int &FileID /// [in] ID of the file to be closed
);

/// Retrieves a dimension length from an input file, given the name
/// of the dimension. Returns the length if exists, but returns a negative
/// value if a dimension of that length is not found in the file.
int getDimFromFile(int FileID, ///< [in] ID of the file containing dim
                   const std::string &DimName, ///< [in] name of dimension
                   int &DimID,    ///< [out] ID assigned to this dimension
                   int &DimLength ///< [out] global length of the dimension
);

/// Defines a dimension for an output file. Returns a dimension id to
/// be used in future field definitions as well as an error flag.
int defineDim(int FileID, ///< [in] ID of the file containing dim
              const std::string &DimName, ///< [in] name of dimension
              int Length,                 ///< [in] length of dimension
              int &DimID                  ///< [out] dimension id assigned
);

/// Writes metadata (name, value) associated with a variable and/or file.
/// The variable ID can be GlobalID for global file and/or simulation
/// attributes/metadata. Specific interfaces for each data type are aliased
/// to a generic interface form.
int writeMeta(const std::string &MetaName, ///< [in] name of metadata
              I4 MetaValue,                ///< [in] value of metadata
              int FileID,                  ///< [in] ID of the file for writing
              int VarID ///< [in] ID for variable associated with metadata
);
int writeMeta(const std::string &MetaName, ///< [in] name of metadata
              I8 MetaValue,                ///< [in] value of metadata
              int FileID,                  ///< [in] ID of the file for writing
              int VarID ///< [in] ID for variable associated with metadata
);
int writeMeta(const std::string &MetaName, ///< [in] name of metadata
              R4 MetaValue,                ///< [in] value of metadata
              int FileID,                  ///< [in] ID of the file for writing
              int VarID ///< [in] ID for variable associated with metadata
);
int writeMeta(const std::string &MetaName, ///< [in] name of metadata
              R8 MetaValue,                ///< [in] value of metadata
              int FileID,                  ///< [in] ID of the file for writing
              int VarID ///< [in] ID for variable associated with metadata
);
int writeMeta(const std::string &MetaName,  ///< [in] name of metadata
              const std::string &MetaValue, ///< [in] value of metadata
              int FileID,                   ///< [in] ID of the file for writing
              int VarID ///< [in] ID for variable associated with metadata
);

/// Reads metadata (name, value) associated with a variable and/or file.
/// The variable ID can be GlobalID for global file and/or simulation
/// metadata. Specific interfaces for each supported type are aliased
/// to a generic interface form.
int readMeta(const std::string &MetaName, ///< [in] name of metadata
             I4 &MetaValue,               ///< [out] value of metadata
             int FileID,                  ///< [in] ID of the file for writing
             int VarID ///< [in] ID for variable associated with metadata
);
int readMeta(const std::string &MetaName, ///< [in] name of metadata
             I8 &MetaValue,               ///< [out] value of metadata
             int FileID,                  ///< [in] ID of the file for writing
             int VarID ///< [in] ID for variable associated with metadata
);
int readMeta(const std::string &MetaName, ///< [in] name of metadata
             R4 &MetaValue,               ///< [out] value of metadata
             int FileID,                  ///< [in] ID of the file for writing
             int VarID ///< [in] ID for variable associated with metadata
);
int readMeta(const std::string &MetaName, ///< [in] name of metadata
             R8 &MetaValue,               ///< [out] value of metadata
             int FileID,                  ///< [in] ID of the file for writing
             int VarID ///< [in] ID for variable associated with metadata
);
int readMeta(const std::string &MetaName, ///< [in] name of metadata
             std::string &MetaValue,      ///< [out] value of metadata
             int FileID,                  ///< [in] ID of the file for writing
             int VarID ///< [in] ID for variable associated with metadata
);

/// Defines a variable for an output file. The name and dimensions of
/// the variable must be supplied. An ID is assigned to the variable
/// for later use in the writing of the variable.
int defineVar(int FileID, ///< [in] ID of the file containing dim
              const std::string &VarName, ///< [in] name of variable
              IODataType VarType,         ///< [in] data type for the variable
              int NDims,                  ///< [in] number of dimensions
              int *DimIDs, ///< [in] vector of NDims dimension IDs
              int &VarID   ///< [out] id assigned to this variable
);

/// Ends define mode signifying all field definitions and metadata
/// have been written and the larger data sets can now be written
int endDefinePhase(int FileID ///< [in] ID of the file being written
);

/// Creates a PIO decomposition description to describe the layout of
/// a distributed array of given type.
int createDecomp(
    int &DecompID,                      ///< [out] ID for the new decomposition
    IODataType VarType,                 ///< [in] data type of array
    int NDims,                          ///< [in] number of array dimensions
    const std::vector<int> &DimLengths, ///< [in] global dimension lengths
    int Size,                           ///< [in] local size of array
    const std::vector<int> &GlobalIndx, ///< [in] global indx for each loc indx
    Rearranger Rearr                    ///< [in] rearranger method to use
);

/// Removes a PIO decomposition to free memory.
int destroyDecomp(int &DecompID ///< [inout] ID for decomp to remove
);

/// Reads a distributed array. We use a void pointer here to create
/// a generic interface for all types. Arrays are assumed to be in contiguous
/// storage so the arrays of any dimension are treated as a 1-d array with
/// the full local size. The routine returns the array as well as the id
/// assigned to the array should that be needed to retrieve variable metadata.
int readArray(void *Array,                ///< [out] array to be read
              int Size,                   ///< [in] local size of array
              const std::string &VarName, ///< [in] name of variable to read
              int FileID,                 ///< [in] ID of open file to read from
              int DecompID, ///< [in] decomposition ID for this var
              int &VarID    ///< [out] variable ID in case metadata needed
);

/// Reads a non-distributed variable. We use a void pointer here to create
/// a generic interface for all types. Arrays are assumed to be in contiguous
/// storage so the arrays of any dimension are treated as a 1-d array with
/// the full local size. The routine returns the variable as well as the id
/// assigned to the variable should that be needed later.
int readNDVar(void *Variable,             ///< [out] variable to be read
              const std::string &VarName, ///< [in] name of variable to read
              int FileID,                 ///< [in] ID of open file to read from
              int &VarID ///< [out] variable ID in case metadata needed
);

/// Writes a distributed array. A void pointer is used to create a generic
/// interface. Arrays are assumed to be in contiguous storage and the variable
/// must have a valid ID assigned by the defineVar function. A void pointer
/// to a scalar FillValue is also required to fill missing values.
int writeArray(void *Array,     ///< [in] array to be written
               int Size,        ///< [in] size of array to be written
               void *FillValue, ///< [in] value to use for missing entries
               int FileID,      ///< [in] ID of open file to write to
               int DecompID,    ///< [in] decomposition ID for this var
               int VarID        ///< [in] variable ID assigned by defineVar
);

/// Writes a non-distributed variable. A void pointer is used for a generic
/// interface. Arrays are assumed to be in contiguous storage and the variable
/// must have a valid ID assigned by the defineVar function. A void pointer
/// to a scalar FillValue is also required to fill missing values.
int writeNDVar(void *Variable, ///< [in] variable to be written
               int FileID,     ///< [in] ID of open file to write to
               int VarID       ///< [in] variable ID assigned by defineVar
);

} // end namespace IO
} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_IO_H
