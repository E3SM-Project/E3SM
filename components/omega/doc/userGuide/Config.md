(omega-user-config)=

## Model Configuration (Config)

Model configuration refers to all of the necessary variables to configure
and run the model. It contains all the input parameters, choices of methods,
physical coefficients for parameterizations and various options for I/O.
Typically, this is read as a single input file that can serve as
a description of the configuration for provenance as well.

In Omega, the input configuration file is in [YAML](https://yaml.org) format
and must be a file named ``omega.yml`` in the same directory as the executable.
However, Unix soft links can be used to point a link with that name to a
specific configuration file stored elsewhere.
In YAML format, most of the configuration variables are typically represented
as key-value pairs called maps, where the name and the value are separated
by a colon. However, lists of variables (eg contents of an IO stream) can
be used and are called sequences. Sequences can be in the form of a list
with an entry on each line starting with a dash and space. They can also
be formatted as a vector with a comma-delimited list within square brackets.
Each of these can be nested to build up the full configuration with indentation
used to separate the nests.
Each map or sequence is called a node in YAML.  A typical omega configration
file will start with the main omega map node with sub-maps associated with
various modules in omega. For example, a file might look like this:

```yaml
omega:
   TimeManagement:
      RestartOn: false
      RestartTimestampName: restartTimestamp
      StartTime: 0001-01-01_00:00:00
      StopTime: none
      RunDuration: 0010_00:00:00
      CalendarType: noleap

   [Other config options in a similar way]

   Hmix:
      HmixScaleWithMesh: false
      MaxMeshDensity: -1.0
      HmixUseRefWidth: false
      HmixRefWidth: 30.0e3

   [more config options]

   MyVector: [1, 2, 3, 4, 5]

   Streams:

      Mesh:
         Type: input
         FilenameTemplate: mesh.nc
         InputInterval: initial_only

      Output:
         Type: output
         FilenameTemplate: output/output.$Y-$M-$D_$h.$m.$s.nc
         FilenameInterval: 01-00-00_00:00:00
         ReferenceTime: 0001-01-01_00:00:00
         Precision: single
         OutputInterval: 0001_00:00:00
         Contents:
         - tracers
         - layerThickness
         - ssh
         - kineticEnergyCell
         - relativeVorticityCell
         - [other fields]

      [other streams in similar form]
```

The variable names in the example above may not be actual omega variables and
are shown only as an example.
We intend to create a script that will generate a default input file that
users can modify for standalone ocean experiments but this capability is not
currently available. Similarly, standard configurations for E3SM compsets will
be provided as part of E3SM releases. The configuration supports all Omega
data types (bool, string, I4, I8, R4, R8). Because YAML internally represents
all data as strings, configuration variables are converted by omega into the
correct type and this conversion may slightly change floating point results
in the last digits, but in a consistent manner based on C++ coercion rules.
This can be minimized by supplying all digits for the precision desired.

Users will configure the model primarily by editing or modifying this
input YAML-formatted file.  Details of the implementation within omega
can be found in the [Developer's Guide](#omega-dev-config) and the actual
interfaces for extracting configuration variables into the modules that
own them are described there as well.
