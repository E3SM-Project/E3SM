(omega-design-IO)=

# Input/Output (IO)

## 1 Overview

The OMEGA model must be able to read and write data to a filesystem.
For performance at high resolution, much of this IO must occur in parallel
and interact with a high-performance filesystem. We describe here an IO
layer that provides interfaces to the underlying SCORPIO parallel I/O
library used in E3SM. It primarily provides a translation layer to
read/write OMEGA metadata and Kokkos arrays. It works together with the
OMEGA Metadata and IOStreams capabilities. Users will interact with
IO primarily through the IOStreams and should not need to access this
layer directly.

## 2 Requirements

### 2.1 Requirement: Read fields and metadata

The model must be able to read desired data and metadata from
an input file into the internal model storage/decomposition.
Metadata requirements are described in the Metadata design document.

### 2.2 Requirement: Write fields and metadata

The model must be able to write desired data and metadata
to an output file for later use. Metadata requirements are
described in the Metadata design document.

### 2.3 Requirement: Self-describing formats

All files must (eventually) be in a self-describing
format like HDF5 or netCDF. Note that alternative formats
can be supported as long as scripts are provided to convert
to netCDF/HDF5 efficiently as part of model run scripts.
SCORPIO supports a number of underlying formats and the
user must be able to select the appropriate format (eg
multiple netcdf formats and ADIOS format).

### 2.4 Requirement: multiple files/streams

The model must support reading from and writing to an
arbitrary number of files with different properties
(eg precision, time frequencies, contents). There will
be a Streams capability described in a separate
design document that will manage many of these properties,
but the underlying I/O layer here must be able to
support multiple open files, each with unique
properties.

### 2.5 Requirement: parallel I/O

Performance at high resolution requires a parallel I/O
solution in which multiple processors can be writing/reading
data to maximize bandwidth to the filesystem and minimize
time spent in I/O.

### 2.6 Requirement: parallel I/O tuning

Some properties of the parallel I/O must be configurable
to tune the parallel I/O for performance on a particular
architecture and filesystem. At a minimum, the user must
be able to specify the number of I/O tasks and the stride
of those tasks.

### 2.7 Requirement: Data types and type conversion

The I/O system must be able to read/write all supported
data types for metadata and all supported Kokkos array types.
In some cases, output files at reduced precision are required
so an option to convert data to reduced precision is needed.

### 2.8 Requirement: Kokkos arrays and host/device support

Distributed data in OMEGA is stored as Kokkos array types.
We must be able to read/write Kokkos arrays and be able
to move data between host and device as needed.

### 2.9 Requirement: Modes on file existence

If an output file exists, the user must be able to specify
whether the model should overwrite the existing file, exit
with an error message or append data to the existing file.
On input, if a file does not exist, the model must exit with
an appropriate error code.

### 2.10 Desired: Asynchronous I/O

For performance, it would be desireable to enable the model
to launch I/O tasks and resume computation while the actual
writing to file takes place. This will mask I/O costs by
overlapping I/O functions with model computations.

### 2.11 Desired: File compression

It may be desireable to apply data compression while
performing I/O to minimize data file sizes. This compression
may be lossless or not depending on the use for each
file.


## 3 Algorithmic Formulation

Parallel I/O systems use various algorithms for rearranging
data in parallel decompositions to the decomposition required
for the I/O processor layout. These algorithms are described
in the I/O library (scorpio) documentation and related
publications.

## 4 Design

The OMEGA model I/O will be built on top of the SCORPIO
parallel I/O library used across E3SM components. The
I/O interfaces here generally provide wrappers for
translating internal OMEGA metadata representations and
Kokkos array types to the form required by SCORPIO. OMEGA
users and developers will generally interact with I/O
through the IO Streams layer that manages all files and
associated file contents.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

There will be a section in the input OMEGA configuration file
for managing overall parallel IO options. Currently, this
will include three variables:

```yaml
    IO:
       ioTasks:  1
       ioStride: 1
       ioRearranger: box
```

The default for IO tasks and stride will be 1 for safety, but
should always be changed to appropriate values for a particular
simulation and machine architecture. For example, the user
should set the stride such that the number of cores performing
I/O is limited to the optimal number for a given multi-processor
node. Scorpio supports a few rearranger methods for moving the
data to IO tasks (see Scorpio documentation) but the box
rearranger is often preferred and is the default. The rearranger
is ignored for the ADIOS back-end since it writes separate files
from each task and rearranges the data in a post-processing step.
These input variables correspond to variables of the
same name.

We define a number of enums to support various options, including
the rearranger method, file format, operation (read/write) and
behavior when file exists:

```c++
   enum class IORearranger {
      box,    /// box rearranger (default)
      subset, /// subset rearranger
      unknown /// unknown or undefined
      };

   enum class IOFileFormat {
      undefined,
      default,
      netcdf,
      pnetcdf,
      netcdf4,
      pnetcdf5,
      hdf5,
      adios,
      };

   enum class IOMode {
      read,
      write,
      };

   enum class IOIfExists {
      overwrite,
      fail,
      };

   enum class IOPrecision {
      single,
      double,
      };
```


#### 4.1.2 Class/structs/data types

The PIO routines require some information about the decomposition,
so there will be one class containing environment information including
array decompositions for each data type, mesh location and array size.

```c++
   class IOEnv {

      private:

         /// track and store all defined environments internally
         static std::vector<std::shared_ptr<IOEnv>> definedEnvs;

         /// name for this environment
         std::string name;

         /// default file format and rearranger method
         IOFormat format;
         IORearranger rearranger;

         /// ptr to defined Scorpio IO system with MPI communicator info
         int* iosystem;

         /// pointers to defined decompositions for every multi-dim
         /// array, data type and cell location
         int* decompCell1DI4;
         int* decompCell1DR4;
         int* decompCell1DR8;
         int* decompCell2DI4;
         int* decompCell2DR4;
         int* decompCell2DR8;
         // continue for dimensions up to 5 and edge, vertex
         //   mesh locations
         // note that Scorpio does not support logical arrays
         //   or I8 types so these will be converted to an
         //   an appropriate type during read/write

         // Methods below will be private and only accessed
         // via the OMEGA IOstreams interfaces with IOstreams
         // as a friend class

      public:
         // Users will not access methods or vars directly
         // but only through the IOstreams interfaces

      friend class IOStreams;
   }
```

There will also be an IOField class to combine the Metadata
and a pointer to the array holding data. This will be used
to specify the fields in an IOStream.

```c++
   class IOField {

      private:

         /// track and store all defined fields internally
         static std::vector<std::shared_ptr<IOField>> definedFields;

         std::shared_ptr<Metadata> metaPtr;

         /// only one of the pointers will be defined based on array type
         std::shared_ptr<Array1DI4> data1DI4;
         [replicated up to 5D arrays of I4,I8,R4,R8]
         std::shared_ptr<HostArray1DI4> dataHost1DI4;
         [replicated up to 5D arrays of I4,I8,R4,R8]


      public:
         // See methods below

      friend class IOStreams;
   }
```


### 4.2 Methods

As noted above, these methods are actually private and accessed
only through the IOStreams interfaces.

#### 4.2.1 Environment constructor

There will be a constructor that initializes the IO system
and decompositions. Note that the defaults for file format
and rearranger can be overridden on a per-file basis.

```c++
   IOEnv(const std::string name, /// name for this env
         const MachEnv omegaEnv, /// machine env with MPI info
         const Decomp  omegaDecomp, /// mesh decomposition
         const int ioTasks,      /// number of IO tasks
         const int ioStride,     /// stride in MPI ranks for io tasks
         const IOFormat format,  /// default IO format
         const IORearranger rearranger, /// default rearranger method
         );
```

#### 4.2.2 File open/close

There will be interfaces for opening and closing files for reading and writing.
For the file open, a file id will be returned.

```c++
   /// Opens a file for reading or writing, returns an error code.
   /// If the IOFileFormat argument is `default`, the default
   /// format from IOEnv will be used. `ifexists` will be ignored for
   /// reads.
   int IOFileOpen(int& fileID, /// returned fileID for this file
                  const std::string filename,  /// name of file to open
                  const IOEnv& myEnv,          /// overall IO environment
                  const IOMode mode,           /// mode (read or write)
                  const IOPrecision precision, /// precision of floats
                  const IOIfExists ifexists,   /// behavior if file exists
                  const IOFileFormat format,   /// file format
                 );

   /// closes an open file using the fileID, returns an error code
   int IOFileClose(int& fileID  /// ID of the file to be closed
                  );

```
#### 4.2.3 Write operations

Because the IOStreams and metadata manager have aggregated all
metadata, dims and pointers to data arrays, a file can be written
with a single call using the aggregated data from the IOStream:

```c++
   /// writes all data associated with the file
   int IOWrite(const int fileID, /// id for the open file
               const std::vector<std::shared_ptr<IOField>> contents);
```

#### 4.2.4 Read operations

Unlike the write function, not all data within a file may need to
be read, so we read each field individually:

```c++
   /// reads a field the file
   int IORead(const int fileID, /// id for the open file
              std::shared_ptr<IOField> field);
```

#### 4.2.5 Defining IO fields

Each field meant to be read or written should be defined, typically
at the same time a module defines the associated metadata. The IOField
is then constructed first with a pointer to the Metadata:

```c++
   /// IOField constructor
   IOField(std::shared_ptr<Metadata> fieldMeta);
```

Then a pointer to the array is attached using an attach function
(aliased to the various array types):

```c++
   /// IOField attach data
   IOField::attachData(std::shared_ptr<Array1DI4>);
   [replicate generic interface for all supported array types]
```


## 5 Verification and Testing

### 5.1 Test via IOStreams

Because the functions here are only accessed via the IOStreams interfaces,
the testing of these routines will be part of the IOStreams unit test.
  - tests requirement 2.1-2.9
