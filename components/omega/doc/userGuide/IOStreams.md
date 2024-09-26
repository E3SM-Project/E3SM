(omega-user-iostreams)=

## IO Streams (IOStream)

IO Streams are the primary mechanism for users to specify input and output
for Omega. An IOStream can be defined for any number of fields and at desired
time frequencies (including one-time or at startup/shutdown). IOStreams are
defined in the Omega input configuration file in an IOStreams section:

```yaml
Omega:
  # other config options removed for brevity
  IOStreams:
    InitialState:
      UsePointerFile: false
      Filename: OmegaMesh.nc
      Mode: read
      Precision: double
      Freq: 1
      FreqUnits: OnStartup
      UseStartEnd: false
      Contents:
        - Restart
    RestartWrite:
      UsePointerFile: true
      PointerFilename: ocn.pointer
      Filename: ocn.restart.$Y-$M-$D_$h.$m.$s
      Mode: write
      IfExists: replace
      Precision: double
      Freq: 6
      FreqUnits: months
      UseStartEnd: false
      Contents:
        - Restart
    History:
      UsePointerFile: false
      Filename: ocn.hist.$SimTime
      Mode: write
      IfExists: replace
      Precision: double
      Freq: 1
      FreqUnits: months
      UseStartEnd: false
      Contents:
        - Tracers
    Highfreq:
      UsePointerFile: false
      Filename: ocn.hifreq.$Y-$M-$D_$h.$m.$s
      Mode: write
      IfExists: replace
      Precision: single
      Freq: 10
      FreqUnits: days
      UseStartEnd: true
      StartTime: 0001-06-01_00.00.00
      EndTime: 0001-06-30_00.00.00
      Contents:
        - Tracers
```

Each stream has a number of required and optional parameters for customizing
input and output. These options are indented below the stream name as shown
in the sample YAML entries above. They include:
- **UsePointerFile:** A required flag that is either true or false. A pointer
file is used for cases like restart files where the last file written can
be stored for the next job submission so that the configuration file does
not need to be edited between job submissions.
- **PointerFilename:** Only required if UsePointerFile is true and should
be set to the full filename (with path) for the pointer file. Each stream
using a pointer file must define a unique pointer file name.
- **Filename:** Required in all cases except input streams using a pointer
file. This is the complete name (with path) of the file to be read or written.
A filename template is also supported in which simulation (or real) time
can be used in the file name. As the examples above show, accepted keys for
a template can be:
  - $SimTime for the current simulation time in a standard time string (note
     that this time string may include colon separators that can be a problem
     for filenames so using the individual keys below is preferred).
  - $Y for the current simulation year
  - $M for the current simulation month
  - $D for the current simulation day
  - $h for the current simulation hour
  - $m for the current simulation minute
  - $s for the current simulation second
  - $WallTime for the time IRL for use when you might need the actual time for
    a debug time stamp
- **Mode:** A required field that is either read or write. There is no
   readwrite option (eg for restarts) so a separate stream should be
   defined for such cases as in the examples above.
- **IfExists:** A required field for write streams that determines behavior
   if the file already exists. Acceptable options are:
   - Fail if you want the code to exit with an error
   - Replace if you want to replace the existing file with the new file
   - Append if you want to append (eg multiple time slices) to the existing
     file (this option is not currently supported).
- **Precision:** A field that determines whether floating point numbers are
   written in full (double) precision or reduced (single). Acceptable values
   are double or single. If not present, double is assumed, but a warning
   message will be generated so it is best to explicitly include it.
- **Freq:** A required integer field that determines the frequency of
   input/output in units determined by the next FreqUnits entry.
- **FreqUnits:** A required field that, combined with the integer frequency,
   determines the frequency of input/output. Acceptable values include:
   - OnStartup for files read/written once on startup
   - OnShutdown for files read/written only once on model exit
   - AtTime or OnTime or Time or TimeInstant for a one-time read or write
     at the time specified in the StartTime entry
   - Years for a frequency every Freq years (*not* Freq times per year)
   - Months for a frequency every Freq months (*not* Freq times per month)
   - Days for a frequency every Freq days (*not* Freq times per day)
   - Hours for a frequency every Freq hours (*not* Freq times per hour)
   - Minutes for a frequency every Freq minutes (*not* Freq times per minute)
   - Seconds for a frequency every Freq seconds (*not* Freq times per seconds)
- **UseStartEnd:** A required entry that is true or false and is used if the
   I/O is desired only within a certain time interval. An example might be
   for specifying high-frequency output within a certain period of a simulation.
- **StartTime:** A field only required when UseStartEnd is true or if
   the FreqUnits request a one-time read/write. The StartTime must be a time
   string of the format YYYY-MM-DD_hh.mm.ss (though the delimiters can be
   any non-numeric character). The year entry is the integer year and can be
   four or more digits. The StartTime is inclusive - the I/O will occur at or
   after that date/time.
- **EndTime:** A field that is only required when UseStartEnd is true. It
   requires the same format as StartTime but unlike StartTime, the EndTime
   is not inclusive and I/O only occurs for times before the EndTime. If a
   file is desired at the EndTime, the user should specify an EndTime slightly
   later (less than a time step) than the desired end time.
- **Contents:** This is a required field that contains an itemized list of
   each Field or FieldGroup that is desired in the output. The name must
   match a name of a defined Field or Group within Omega. Group names are
   preferred to keep the list of fields short so Omega will define convenient
   FieldGroups like Restart, State, Tracers that will include all members
   of the group. If only a subset of Fields from a Group is desired, the
   individual Field names should be specified and not the Group name.

This streams configuration should be sufficient to define all input and output
from the model and provides a relatively simple interface for a typical user.
However, if necessary (eg before streams have been defined), the specific
interfaces in the lower level [IO](#omega-user-IO) module can be used.
