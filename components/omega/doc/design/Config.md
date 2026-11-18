(omega-design-config)=
# Config

## 1 Overview

Model configuration refers to all of the necessary variables to configure
and run the model. It contains all the input parameters, choices of methods,
physical coefficients for parameterizations and various options for I/O.
Typically, this is read as a single input file that can also serve as
a description of the configuration for provenance as well.


## 2 Requirements

### 2.1 Requirement: Human readability

File or files containing the configuration must be easily understood and
modified by non-expert users. It is also desirable to have minimum markup
that would interfere with readability.

### 2.2 Requirement: Standard format

Configuration files must conform to a standard format for ease in parsing
and for potential interoperability. Examples of standards typically used
include YAML, JSON, XML.

### 2.3 Requirement: Archiving and provenance

The entire model configuration must be able to be saved and archived to
act as model provenance.

### 2.4 Requirement: Internal accessibility

Although most input parameters will only be relevant to a particular
module or class, it is inevitable that some parameterizations or process
models will have dependencies on other configuration choices. All
configuration variables must therefore be accessible to all other
modules/classes.

### 2.5 Requirement: Support of data types

Model configuration must support parameters in the standard data types,
including logical (boolean), integer, float/double, and string data types.
Limited support for vectors or arrays of the above types may be desirable.

### 2.6 Requirement: Efficiency and parallelism

Configuration should never consume significant resources, in either run
time or storage. For parallel execution, the configuration
will need to be replicated across MPI ranks, though in many cases, the
variables will be manipulated into a different form (eg string variables
describing choices will be converted to enums). This may happen before
or after replication.

### 2.8 Desired: Single configuration input

For simplicity, we desire a single file for model configuration where users
can define all model configuration details and provenance can be maintained
more easily. If multiple files are required, there should be obvious links
or includes in the main configuration file to make the dependency clear.

### 2.9 Required: Hierarchy or grouping

For both readability and for easier encapsulation by modules/classes, the
configuration should be organized in a logical way to make it clear where
specific configuration variables are to be used. Parameters primarily or
specifically associated with a module or class should be in a grouping
of those parameters with the configuration.

### 2.10 Requirement: Language support

Because the configuration will also be used for provenance, it is likely
that the configuration input file will need to be read by other languages
outside of the OMEGA C++ model. This requires an ability to parse a
configuration input file from other common languages (eg python).

### 2.11 Requirement: Optional or missing values

If a configuration variable is missing from the configuration input file,
there must be an option to either supply a default value or throw an error
and exit, depending on the user choice. If a default value is supplied,
then the default value must be added to the configuration so that it is
included in the output and become part of the provenance.

### 2.12 Requirement: Extra values

If extra or unexpected values are encountered, they will be ignored.

### 2.13 Desired: Automated generation of default input and error checking

While the source code defines the configuration variables, it would
be desirable to have a means to extract from the source code what the
code is expecting into a default input config file. This would also
enable some external error checking for missing or extra entries.

### 2.14 Desired: Acceptable values

For users modifying an input configuration, it would be desirable
to document the acceptable values or range of values that each
variable can be assigned.

## 3 Algorithmic Formulation

There are no specific algorithms needed for this other than those
used in the anticipated parsing/storage packages associated with
the standard format chosen.


## 4 Design

We select the YAML format that meets the above requirements with
improved readability over other standard forms. Within the OMEGA
model, we use the yaml-cpp library, a third-party implementation for
parsing YAML input and efficiently storing/retrieving information
as YAML nodes. Many of the configuration variables will be stored as
maps in YAML, the package features a map syntax similar to the C++
standard template library map type.

Because extracting variables from a large and complex config
structure is less efficient, we plan for a single reading of the
full configuration. Each initialization routine in OMEGA will then
extract needed variables and manipulate them as needed for the later
forward integration. In this implementation, we will read/write the
configuration from a master rank and each initialization routine
within OMEGA will manipulate and broadcast necessary variables
across ranks for parallel execution. Because not all variables will
need to be broadcast and many others will be converted to more
efficient types (eg string options to logical or enums), we believe
this model will be more efficient than broadcasting the full config
structure and manipulating afterward.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

There are no global parameters or shared constants.

#### 4.1.2 Class/structs/data types

We define a Config type, which is actually an alias of YAML::node:

```c++
using Config = YAML::node;
```

from the yaml-cpp library. A YAML node is more fully and accurately
defined in the YAML specification, but for the purposes of this
document, our configuration is represented in YAML as a set of nested
map nodes, where a map is simply a keyword-value pair. At the lowest
level, these nodes are the simple `variable-name: value` maps. The next
level up is a map of the module name to the collection of maps associated
with the module. The root node corresponds to the full model configuration
and is simply a collection of all those module maps. I/O stream/file
configuration will be part of this configuration in a design TBD later
but will be a similar hierarchy under the full omega config node.

An example YAML input file might then look like:

```yaml
omega:
   timeManagement:
      doRestart: false
      restartTimestampName: restartTimestamp
      startTime: 0001-01-01_00:00:00
      stopTime: none
      runDuration: 0010_00:00:00
      calendarType: noleap

   [Other config options in a similar way]

   hmix:
      hmixScaleWithMesh: false
      maxMeshDensity: -1.0
      hmixUseRefWidth: false
      hmixRefWidth: 30.0e3

   [more config options]

   streams:

      mesh:
         type: input
         filenameTemplate: mesh.nc
         inputInterval: initial_only

      output:
         type: output
         filenameTemplate: output/output.$Y-$M-$D_$h.$m.$s.nc
         filenameInterval: 01-00-00_00:00:00
         referenceTime: 0001-01-01_00:00:00
         clobberMode: truncate
         precision: single
         outputInterval: 0001_00:00:00
         contents:
         - tracers
         - layerThickness
         - ssh
         - kineticEnergyCell
         - relativeVorticityCell
         - [other fields]

      [other streams in similar form]
```

### 4.2 Methods

All of the methods in the YAML::Node class are obviously supported,
but we will alias or wrap some of the most common in the OMEGA context
to be associated with Config.

#### 4.2.1 File read and master config

The most common use case should be creating a Config by reading a
YAML configuration file using:

```c++
Config omegaConfig = ConfigRead("omega.yml");
```

where the argument is the name for the YAML input file. In OMEGA,
we will retain this master configuration throughout the initialization
as omegaConfig.

#### 4.2.2 Get/Retrieval

Once the configuration has been read, we will need to retrieve variables
from the Config. Because our config is a hierarchy, there are really no
variables at the top level and we need to first retrieve the sub-config
associated with the local module/group. In the sample above, if we need
to retrieve a variable from the hmix group, we first retrieve the hmix
config and then the variable using:

```c++
Config hmixConfig = ConfigGet(omegaConfig,"hmix",iErr);
Real refWidth{0.0};
bool useRefWidth{false};
refWidth    = ConfigGet(hmixConfig, "hmixRefWidth",    iErr);
useRefWidth = ConfigGet(hmixConfig, "hmixUseRefWidth", iErr);
```

where there is a retrieval function for all supported Omega data types:
bool, I4, I8, R4, R8, Real, std::string. These retrievals are just
overloaded wrappers around the YAML form: `configName["varName"].as<type>`
with some error checking and reporting. Rather than a templated form,
we use simple overloading to keep a cleaner interface. If the variable
or config is missing, these functions will print an error message and
return with a non-zero error argument.

Another interface will allow the setting of a default value if the
variable is missing from the input config.  This interface simply
adds the default value as an additional argument, for example:

```c++
refWidth = ConfigGet(hmixConfig, "hmixRefWidth", defaultVal, iErr);
```

In this case, if the variable does not exist, it will not only
use the default value but print a warning that the default is
being used because the entry is missing.


#### 4.2.3 Change an existing value

While the intent is for all config variables to be set using the config
file read interface above, the capability modify a value is also
required. The syntax is essentially the inverse of the get/retrieval
above. Similar to that case, the sub-group will need to be retrieved first.

```c++
Config hmixConfig = ConfigGet(omegaConfig, "hmix", iErr);
ConfigSet(hmixConfig, "hmixRefWidth", 10.0e3, iErr);
ConfigSet(hmixConfig, "hmixUseRefWidth", true, iErr);
```

There will be overloaded interfaces for each supported type. For
literals (as in the example above), they will be cast to an appropriate
type according to C++ default type conversion and will be converted
to the desired type on retrieval (YAML internal storage is ignorant of
the type and only performs the type cast on retrieval).

#### 4.2.3 Adding new entries

It may be necessary to build up a configuration that does not yet
exist or add new entries to an existing group. We provide an Add
interface to distinguish this case from the Set case above.

```c++
// For an existing subgroup:
Config hmixConfig = ConfigGet(omegaConfig, "hmix", iErr);
ConfigAdd(hmixConfig, "hmixRefWidth", 10.0e3, iErr);
ConfigAdd(hmixConfig, "hmixUseRefWidth", true, iErr);

// To add a new subgroup:
Config hmixConfig; // empty Config constructor
ConfigAdd(hmixConfig, "hmixRefWidth", 10.0e3, iErr); // build subgroup
ConfigAdd(hmixConfig, "hmixUseRefWidth", true, iErr);
ConfigAdd(omegaConfig, hmixConfig); // add new subgroup to parent
```

There will be overloaded interfaces for each supported type. For
literals (as in the example above), they will be cast to an appropriate
type according to C++ default type conversion and will be converted
to the desired type on retrieval (YAML internal storage is ignorant of
the type and only performs the type cast on retrieval).

#### 4.2.4 Existence

It is not expected that a user would test the existence since the
Get/Set functions will perform the test internally. However, to
satisfy requirement 2.11, we will add a function to test the
existence of an entry, given a config or sub-config.
Using the hmix example again:

```c++
if (ConfigExists(hmixConfig,"hmixRefWidth") {
   // variable exists, do stuff
}
```

Note that this can also be used to test the existence of a complete
sub-group as well:

```c++
bool hmixExists = ConfigExists(omegaConfig, "hmix");
```

#### 4.2.5 File write

While we may decide to save provenance a different way, a write interface
is supplied to write a configuration to an output YAML file:

```c++
err = ConfigWrite(myConfig, "outputFileName");
```

#### 4.2.6 Constructor/destructor

A default constructor for an empty Config and destructor will be provided:

```c++
Config myConfig;
delete myConfig;
```

The destructor may be important to free up space since the Config is
likely to only be used during the init phase in the current plan.

### 4.3 Documentation and Auto-generation of Default Inputs

In order to keep the source code, input files and documentation consistent
and avoid missing or extra entries, we propose inserting a block within
each source code header (where other interfaces will be documented).
This block would follow Doxygen-like format and look something like:

```c++
/// \ConfigInput
/// #
/// # Group description (eg. hmix: the horizontal mix configuration)
/// #
/// groupName:
///    #
///    # Parameter description (eg horizontal mixing coeff)
///    #    more description (units, acceptable values or range)
///    #
///    varName1: defaultValue1
///    #
///    # Parameter description
///    #
///    varName2: defaultValue2
///    [ continue for remaining vars in this block]
///
/// \EndConfigInput (might not be necessary?)
```

The block between the ConfigInput lines could be extracted verbatim
and written (with proper indenting) into a fully documented yaml
default input file or could be extracted and the `#`-delimited
comments stripped for a more concise yaml input file.
In addition, we may be able to similarly extract the same info for
the User and Developer Guides.

## 5 Verification and Testing

The selection of YAML automatically satisfies Requirements 2.1, 2.2 and
2.10. Requirement 2.8 will be enforced by Omega development. The other
requirements will be tested with a parallel unit test driver that performs
the following tests in order and will output the test name and a
PASS/FAIL for use with CTest or other testing frameworks.


### 5.1 Test default constructor

Create an empty config on master rank using default constructor.
  - tests constructor needed to satisfy other requirements

### 5.2 Test Set function

Add configuration variables to the empty config with at least one of
each supported type with a few extras to test behavior for other
requirements. Also add multiple levels of hierarchy and a large
enough number of subgroups and parameters to simulate a full omega
config.
  - tests function needed for requirements 2.3, 2.5, 2.9, 2.11, 2.12

### 5.3 Test Get function

Retrieve all the variables set in the test above and verify they
return identical values.
  - completes test of 2.4, 2.5, 2.9

Broadcast variables from master after retrieval to provide test
of performance for this choice of parallel implementation.
  - tests 2.6

### 5.4 Test for missing variables

Inquire for a missing variable and output a PASS if the variable
was not found. Also add the missing variable and assign a default
value, the retrieve to test this behavior.
  - tests 2.11

### 5.5 Test write and re-read

Write a YAML output file using the above constructed Config. Then
read in the new file. Add an extra variable to the new file.
Verify that the newly read Config matches the original, ignoring
the extra variable.
   - tests 2.3, 2.12
