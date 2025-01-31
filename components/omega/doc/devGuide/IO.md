(omega-dev-IO)=

## Parallel IO (IO)

For input and output of data needed for Omega, we use the Software for
Caching Output and Reads for Parallel I/O
([SCORPIO](https://github.com/E3SM-Project/scorpio)) library. This
library supports parallel reading and writing of distributed arrays in various
self-describing
formats like [NetCDF](https://docs.unidata.ucar.edu/netcdf/),
[HDF5](https://www.hdfgroup.org/solutions/hdf5/),
and [ADIOS](https://csmd.ornl.gov/adios).
SCORPIO is responsible for rearranging data from the data decomposition used
for computation to a different rearrangement for parallel IO. In particular,
for optimal IO, the user would specify a different number of IO tasks than
for computation. For example, the user could specify the number of IO tasks to
match underlying hardware like the network interfaces (NICs) on a node.
Users and developers will generally access IO via
[IOStreams](#omega-user-iostreams)
and other interfaces. The base IO layer described here only provides
an Omega-aware wrapper around SCORPIO calls.

The base interfaces provide functions for file operations (open/close),
reading and writing of metadata, and reading and writing of data arrays.
Interfaces at this level utilize raw pointers to data and assume
contiguous storage for arrays. SCORPIO utilizes integer handles
to various files and data types so these must often be defined or
retrieved for many operations. Although there is no IO class, we encapsulate
the IO routines within the Omega and IO namespaces.

Before using any IO functions, the parallel IO system must be initialized
using:
```c++
   int Err = IO::init(Comm);
```
where Comm is an MPI communicator and should in most cases be the communicator
from the Omega default MachEnv (see [MachEnv](#omega-dev-mach-env)). This
function also extracts the user-defined variables from the model configuration,
include the number of IO tasks, the IO task stride, the default data
rearranger method, and the default file format
(see [User Guide](#omega-user-IO)).

As mentioned above, most I/O operations will take place within the IOStreams
module, but the base IO functions can be accessed directly. To open and close
files for reading/writing, use:
```c++
   int Err = IO::openFile(FileID, Filename, Mode);
   Err = IO::closeFile(FileID);
```
or
```c++
   int Err = IO::openFile(FileID, Filename, Mode, Format, IfExists);
   Err = IO::closeFile(FileID);
```
In both case, an integer FileID is assigned to the open file which is then
used by all subsequent operations. The Filename is a ``std::string`` that
should include the full path and name of the file. The mode can be either
``OMEGA::IO::ModeRead`` or ``OMEGA::IO::ModeWrite``. The shorter interface
uses default values for the file format, overwrite behavior and precision.
In the longer interface, the user must supply these values using the enums:
```c++
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

/// Behavior (for output files) when a file already exists
enum class IfExists {
   Fail,    /// Fail with an error
   Replace, /// Replace the file
   Append,  /// Append to the existing file
};

```
with the defaults for each being ``IO::FmtDefault``, and
``IO::IfExists::Fail``. The E3SM standard format is currently
NetCDF4. Earlier NetCDF formats should be avoided, but are provided in
case an input file is in an earlier format.

Once the file is open, data is read/written using one of two interfaces,
depending on whether the array is decomposed across MPI tasks or not. For
large decomposed arrays, the interface is:
```c++
int Err = IO::readArray (&Array, Size, VariableName, FileID, DecompID, VarID);
int Err = IO::writeArray(&Array, Size, &FillValue,   FileID, DecompID, VarID);
```
where the pointer to the data array is passed and data is assumed to be
contiguous with the full local size of the array passed as Size. The FileID is
the integer ID for the open file. The DecompID is a defined data decomposition
as described below. For reading, the variable name (as a ``std::string``) is
supplied and the variable ID (VarID) is returned in case it is needed for
reading of variable metadata. For writing, a FillValue is supplied to fill
undefined locations in an array and the variable ID must have been assigned
in a prior defineVar call prior to the write as described below.

For arrays or scalars that are not distributed, the non-distributed variable
interface must be used:
```c++
int Err = IO::readNDVar(&Array, VariableName, FileID, VarID);
int Err = IO::writeNDVar(&Array, FileID, VarID);
```
with arguments similar to the distributed array calls above. Note that
when defining dimensions for these fields, the dimensions must be
non-distributed. For scalars, the number of dimensions should be zero.

The IO subsystem must know how the data is laid out in the parallel
decomposition. Both the dimensions of the array and the decomposition
across tasks must be defined. For each dimension, a dimension must be
defined using:
```c++
   int Err = IO::defineDim(FileID, DimName, Length, DimID);
```
where FileID is the ID of an open file, the DimName is a ``std::string``
with the dimension name (eg NCells, NEdges, NVertices, NVertLevels or
NTracers), length is the length of the full global array and DimID is
the ID assigned to this dimension. Note that for reading a file, we
supply the function:
```c++
   int Err = IO::getDimFromFile(FileID, DimName, DimID, DimLength);
```
that can be used to inquire about the dimension length and retrieve the
dimension ID from the file.

Once the dimensions are defined, the decomposition of an array must
be defined using:
```c++
   int Err = createDecomp(DecompID, IODataType, NDims, DimLengths,
                          Size, GlobalIndx, Rearr);
```
where DecompID is the ID assigned to the newly created decomposition,
NDims is the number of dimensions for the array, DimLengths is an
integer ``std::vector`` of the NDims local dimension lengths, Size is the
full size of the local array, and Rearr is the data rearranger method
to use to map to the number of IO tasks. The Rearr can be set to
``OMEGA::IO::DefaultRearr`` to make use of the overall default defined
when the IO system was initialized, but can also be set explicitly to
``OMEGA::IO::RearrBox`` or ``OMEGA::IO::RearrSubset``. The box rearranger
is generally preferred (see [UserGuide](#omega-user-IO)). The GlobalIndx
array describes the global location (as a zero-based offset) of each
local array entry. This can be computed from the Omega Default Decomp
arrays. For example, an array dimensioned (NCellsAll,NVertLevels) would
have an offset computed using:
```c++
   std::vector<int> OffsetCell(NCellsAll*NVertLevels,-1);
   int Add = 0;
   for (int Cell; Cell = 0; Cell < NCellsOwned){
      for (int k; k = 0; k < NVertLevels){
         OffsetCell[Add] = (CellIDH(Cell)-1)*NVertLevels + k;
         Add++;
      }
   }
```
Note that we exclude Halo layers by assigning an offset of -1. Finally,
the data type of the array must be supplied. To map the standard Omega
data types to the data types used in the IO subsystem, we define:
```c++
enum IODataType {
   IOTypeI4      = PIO_INT,    /// 32-bit integer
   IOTypeI8      = PIO_INT64,  /// 64-bit integer
   IOTypeR4      = PIO_REAL,   /// 32-bit real
   IOTypeR8      = PIO_DOUBLE, /// 64-bit real
   IOTypeChar    = PIO_CHAR,   /// Character/string
   IOTypeLogical = PIO_INT     /// Logicals are converted to ints for IO
};
```
so that in the above interface, we would supply for example ``IO::IOTypeI4``
for an Omega I4 data type.

Now that dimensions and decompositions have been defined, a variable can
be defined (this is required for writing only) using:
```c++
   int Err = IO::defineVar(FileID, VarName, IODataType, NDims, DimIDs, VarID);
```
where VarID is the ID assigned to the variable, FileID is the usual ID of
the data file, VarName is a ``std::string`` holding the variable name,
IODataType is the data type described above, NDims are the number of dimensions,
and DimIDs are an integer ``std::vector`` holding the dimension IDs for each
dimension. For scalar variables, NDims should be set to zero and a null pointer
should be used in place of the DimID argument. Once defined, the variable ID
is used in all IO calls related to this variable.

In addition to data in a file, we can also read and write metadata. As with
the data itself, metadata is typically managed by the IOStreams and Metadata
interfaces, but the base IO module contains interfaces for reading and
writing metadata associated with either an array or the file or simulation
itself. To read/write metadata, use:
```c++
   int Err = IO::writeMeta(MetaName, MetaValue, FileID, VarID);
   int Err = IO::readMeta (MetaName, MetaValue, FileID, VarID);
```
where MetaName is a ``std::string`` holding the name of the metadata and
the MetaValue is the value of the MetaData. All supported Omega data types are
allowed except boolean which must be converted to an integer type. The FileID
is once again the ID of the open data file and VarID is the variable to which
this metadata is attached. For global file and simulation metadata not attached
to a variable, the ID ``IO::GlobalID`` is used to denote global metadata.

For an example of the full read/write process, the IO unit test contains
the full reading and writing of a data file and associated metadata. We
will repeat again that users and developers are not expected to learn and
use the interfaces above, but should access the IO system from the IOStreams
interfaces.
