(omega-dev-iostreams)=

## IO Streams (IOStream)

Most input and output for Omega occurs through IOStreams. Each stream
defines a file, the contents to be read/written and the time frequency
for reading and writing. Defining streams via the input configuration
file is described in the [User Guide](#omega-user-iostreams). IOStreams
are built on top of the parallel IO infrastructure described in the
[IO Section](#omega-dev-IO) and the Field and Metadata described in the
[Field Section](#omega-dev-Field). Here we describe the classes and functions
used to implement IOStreams. Any module accessing an IOStream instance
or related functions must include the ``IOStream.h`` header file.

All IOStreams are initialized in a two-step process. A call to the
init routine should take place early in the Omega initialization after
the ModelClock has been initialized using:
```c++
   int Err = IOStream::init(ModelClock);
```
This routine extracts all the stream definitions from the input configuration
file and creates all the Streams. This initialization also defines the
contents of each Stream but does not yet validate those contents against all
the defined Fields. The contents of all streams should be validated at the
end of initialization (when all Fields have been defined) using the call:
```c++
   bool AllValidate = IOStream::validateAll();
```
However, if a stream is needed (eg a read stream) during initialization
before the validateAll call, a single stream can be validated using
```c++
   bool Validated = MyStream.validate();
```
and the validation status can be checked with
```c++
   bool Validate = MyStream.isValidated();
```
All streams must be validated before use to make sure the Fields have
been defined and the relevant data arrays have been attached to Fields and
are available to access.  At the end of a simulation, IOStreams must be
finalized using
```c++
   int Err = IOStream::finalize(ModelClock);
```
so that any final writes can take place for the OnShutdown streams and to
deallocate all defined streams and arrays. If a stream needs to be removed
before that time, an erase function is provided:
```c++
   IOStream::erase(StreamName);
```

For most output streams, we provide a writeAll interface that should be placed
at an appropriate time during the time step loop:
```c++
   int Err = IOStream::writeAll(ModelClock);
```
This function checks each write stream and writes the file if it is time, based
on a time manager alarm that is defined during initialization for each stream
based on the time frequency in the streams configuration. After writing the
file, the alarm is reset for the next write time. If a file must be written
outside of this routine, a single-stream write can take place using:
```c++
   int Err = IOStream::write(StreamName, ModelClock);
```

Reading files (eg for initialization, restart or forcing) does not often
take place all at once, so no readAll interface is provided. Instead, each
input stream is read using:
```c++
   int Err = IOStream::read(StreamName, ModelClock, ReqMetadata);
```
where ReqMetadata is a variable of type Metadata (defined in Field but
essentially a ``std::map<std::string, std::any>`` for the name/value pair).
This variable should incude the names of global metadata that are desired
from the input file. For example, if a time string is needed to verify the
input file corresponds to a desired time, the required metadata can be
initialized with
```c++
   Metadata ReqMetadata;
   ReqMetadata["ForcingTime"] = "";
```
The Metadata corresponding to ForcingTime will then be read from the file
and inserted as the Metadata value. If no metadata is to be read from the
file, then an empty ReqMetadata variable can be passed.

As described in the [User Guide](#omega-user-iostreams), all streams are
defined in the input configuration file and most other IOStream functions
are associated either with that initialization or to support the read/write
functions above.
