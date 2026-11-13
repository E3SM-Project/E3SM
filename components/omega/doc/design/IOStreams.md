
(omega-design-IOStreams)=

# IOStreams

## 1 Overview

OMEGA must be able to read and write data throughout a simulation.
Data can be read or written either once or periodically. Files
will contain different fields to read/write. We consider each file
or file sequence an IO stream with its own unique set of properties.
We describe here the requirements and design of these IO streams.
The design relies on companion designs for managing metadata and a
lower-level IO functions for writing metadata and data stored in
Kokkos arrays.

## 2 Requirements

### 2.1 Requirement: Multiple streams

An arbitrary number of input and output streams, each with its own
properties and precision, must be supported.

### 2.2 Requirement: Reduced precision

For some output, it is not necessary to retain the full precision
that the model supports. Reduced precision for each stream must
be supported to reduce file size when full precision is not
required.

### 2.3 Requirement: Contents for each stream

Users must be able to supply a list of fields to be included in
each stream. This contents list should be checked against
available fields (from OMEGA defined metadata and metadata groups)
and exit if a field is not available.

### 2.4 Requirement: Exact restart

There must be one stream for checkpointing and restarting and
must include all fields needed to exactly restart the model,
maintaining bit-for-bit when compared to a non-restarted simulation
over the same time interval. This also generally requires full
precision for all fields.

### 2.5 Requirement: Restart pointer

To avoid the need to modify configuration files with a new
input filename on every restart, a mechanism for determining
the last successful restart is needed. This is often implemented
by writing the name of the last successful restart file to a
standard location (pointer file) so that the model can read the
name of the last restart and continue the simulation.

### 2.6 Requirement: Time intervals

The user must be able to specify a time interval for any
repeating input and output streams (eg monthly output).
The time interval can be any interval supported by the OMEGA
time manager.

### 2.7 Requirement: Optional start/stop times

For some streams, we may wish to supply a start and stop
time to sample the simulation only over a particular time
period. For example, requesting high-frequency output over a
short interval. The user must be able to optionally request
a start and stop time for the stream.

### 2.8 Requirement: Filenames

The user must be able to supply a filename associated with the
stream. For repeating I/O, the filename will be a template
that will specify how the filename will be modified to reflect
the time associated with the output. The template should
be able to support any conventions (eg E3SM) required. Filenames
should include the full path.

### 2.9 Requirement: Multiple time slices per file

Streams must be able to support multiple time instances for
the requested fields in a single file. For example, monthly
forcing inputs may be provided in a single input file. Similarly,
for higher-frequency output, it may be desirable to include
several intervals in a single file. An appropriate filename
convention/template must also be provided for this situation.

### 2.10 Requirement: Behavior on file existence

To support multiple time slices per file and to support cases
where the model is re-run over the same interval (eg during
testing or repeating a failed run), the user must specify the
desired behavior when an output file of the same name already
exists. A minimum set of options should include replacing
the file, overwriting portions of a file (eg specific time
slices), appending to the existing file, or exiting with an
error message.

### 2.11 Desired: Time averaging

In some cases, it is desirable to accumulate the time average of
selected fields during the simulation in order to sample every
timestep (or other frequency higher than the output interval). While
the time-averaging capability itself will be implemented elsewhere
in the code, time-averaging over intervals longer than the typical
checkpoint/restart interval would require additional restart
capabilities to maintain these time averages and a filename
template should describe the averaging interval in some format.

### 2.12 Desired: Data compression

In the future, it may be possible to apply lossy or loss-less
compression on data to save space. An option to compress
data would be desirable if it can be supported. Like
reduced-precision, this option would probably be on a per-stream
basis.

## 3 Algorithmic Formulation

No algorithms are needed. Most functions carried out by other
modules and libraries.

## 4 Design

IOStreams will be the primary mechanism used for all input and
output. The streams will be defined in a streams section of the input
Omega configuration file. Because YAML uses some common characters
to denote lists, etc., filenames should be in quotes. For example:

```yaml
   IOStreams:

      mesh:
         mode: read
         name: '/mypath/meshFileName'
         freqUnits: initial
         freq: 1
         precision: double
         contents:
         - GRPmesh
         - [other fields?]

      restart-read:
         mode: read
         name: '/mypath/omega-restart.$Y-$M-$D_$h.$m.$s.nc'
         pointerFilename: '/mypath/pointerFileName'
         freqUnits: initial
         freq: 1
         precision: double
         contents:
         - GRPrestart
         - [other fields or groups as needed]

      restart-write:
         mode: write
         name: '/mypath/omega-restart.$Y-$M-$D_$h.$m.$s.nc'
         pointerFilename: '/mypath/pointerFileName'
         freqUnits: years
         freq: 1
         precision: double
         ifExists: overwrite
         contents:
         - GRPrestart
         - [other fields or groups as needed]

      history:
         mode: write
         name: '/mypath/omega-history.$Y-$M-$D_$h.$m.$s.nc'
         freqUnits: months
         freq: 1
         precision: single
         ifExists: overwrite
         contents:
         - temperature
         - salinity
         - normalVelocity
         - [other fields or groups as needed]

      # 10-day high-freq data stored in an annual file
      highFreq:
         mode: write
         name: '/mypath/omega-highfreq.$Y-$M.nc'
         freqUnits: days
         freq: 10
         precision: single
         ifExists: append
         startDate: yyyy-mm-dd
         endDate: yyyy-mm-dd
         contents:
         - temperature
         - salinity
         - normalVelocity
         - [other fields or groups as needed]
```

Note that there is a stdlib iostreams so our use of IOStreams runs
the risk of inadvertent conflict. Use of both the Omega namespace and
the capitalized IOStreams will be used to distinguish the two.

### 4.1 Data types and parameters


#### 4.1.1 Parameters

There are no parameters specific to the IOStreams, though they will
utilize some parameters (like precision and other options) from a
lower-level IO module.

#### 4.1.2 Class/structs/data types

The main class is an IOStream that carries all the information
related to each input or output stream. Like other classes, all
instantiations of streams are managed within the class.

```c++
   class IOStream {

      private:
         std::string name;      /// name of stream
         std::string filename;  /// filename or template

         IOMode mode;           /// mode (read or write)
         IOPrecision precision; /// precision for floating point vars
         Alarm sAlarm;          /// time mgr alarm for read/write

         bool usePointer;       /// flag for using a pointer file
         std::string ptrFilename; /// name of pointer file

         bool useStartEnd;      /// flag for using start, end times
         TimeInstant startTime; /// start time for stream
         TimeInstant endTime;   /// end   time for stream

         /// Contents of stream in the form of a vector of
         /// pointers to defined IOFields
         std::vector<std::shared_ptr<IOField>> contents;

         /// Store and maintain all defined streams in this vector
         static int numStreams;
         static std::vector<std::shared_ptr<IOStream>> allStreams;

      public:
         [see methods described below]
   };
```

### 4.2 Methods

### 4.2.1 Create/Destroy

There will be a constructor that creates a stream and fills the
scalar variables. And a destructor that deletes the stream and removes
it from the list of defined streams.

```c++
   IOStream(int&  streamID, /// id of created stream
            const std::string name,     /// name of stream
            const std::string filename, /// file name or template
            const IOmode      mode,     /// read/write mode
            const IOPrecision precision, /// precision for floats
            const IOIfExists  ifExists,  /// action if file exists
            const std::string freqUnits, /// time frequency for I/O
            const int         freq,      /// freq in above units for I/O
            const std::string pointerFile, /// pointer filename if used
            const std::string startDate, /// optional start date string
            const std::string endDate.   /// optional end date string
            );

   ~IOStream();
```

### 4.2.2 Add fields

Once an IOStream is created, the contents must be added using one of
two interfaces. If the streamId is still available, the first interface
can be used and is faster (the streamID is the index into the vector
of defined streams). Otherwise, an interface that takes the
stream name can be used. An error code is returned.

```c++
   int IOStreams::AddField(
                const int streamID, /// id of stream to be modified
                const std::string fieldName /// name of field
                );

   int IOStreams::AddField(
                const std::string streamName, /// name of stream
                const std::string fieldName /// name of field
                );
```

### 4.2.3 Write streams

Generally, we want to write all streams at the end of the timestep
if it's time to write, so there will be one interface to simply
check the alarms for all streams and write the data if it's time.
After writing, the alarm will be reset.

```c++
   int IOStreams::WriteAll();
```

However, if the user needs to write a stream from elsewhere in
the code, there will be an interface to write a single stream
by either streamID or by stream name after checking the alarm to
see if it's time to write. The alarm will be reset once the writing
is complete.

```c++
   int IOStreams::Write(const int streamID, /// id of stream to be modified
                       );

   int IOStreams::Write(const std::string streamName, /// name of stream
                       );
```

### 4.2.3 Read streams

Unlike writing, each input stream is typically read from the
appropriate module, so only single-stream reads are provided and
with only the name interface since the id may no longer be
accessible by the calling routine.

```c++
   int IOStreams::Read(const std::string streamName, /// name of stream
                      );

```


## 5 Verification and Testing

### 5.1 Test All

A test driver and configuration file will create a number of streams
with fields, frequencies and other options that attempt to span all
possible configurations (though may not be able to test all file
formats). This driver will march in time to write streams at the
requested frequencies. At the end of the driver, input streams that
mirror the output streams will be read and the fields compared to
determine if they match.
  - tests requirements 2.1-2.9
