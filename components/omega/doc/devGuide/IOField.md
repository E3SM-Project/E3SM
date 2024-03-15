(omega-dev-iofield)=

## IO Fields (IOField)

To determine which fields are available for input or output (I/O), OMEGA
includes an IOField class. An IOField is simply a class that includes a
pointer to field metadata and a pointer to a data array containing the data.
A module registers an IO field by first defining the MetaData
(see [MetaData](#omega-dev-metadata) then calling:
```c++
int Err = OMEGA::IOField::define("FieldName");
```
where the FieldName must match the field name for the previously defined
MetaData. Once the IOField is registered, users can select these fields
and include them in IOStreams (see [IOStreams](#omega-dev-iostreams)).
The actual field data must be in one of the supported array data types
(see [DataTypes](#omega-dev-data-types)) and can be attached to the IOField
using:
```c++
int Err = OMEGA::IOField::attachData<ArrayType>("FieldName", DataArray);
```
where ArrayType is the data type of the DataArray (eg. Array2DI4). The
array can be either a host or device array and the IO system will determine
whether any data transfer is required based on the type. The
DataArray must not be a local temporary because IOStreams will need to
retrieve this array when it comes time to write/read the stream. The
"attach" function is separate from field definition because some data arrays
may not have been initialized at the time of field definition and also
because some array locations may change in time (e.g. different time index or
array for different time levels for prognostic variables).

When it is time to write an IOStream, the IOStream module can extract
the metadata and data for all selected IOFields in the contents using:
```c++
std::shared_ptr<OMEGA::MetaData> MyFieldMeta =
    OMEGA::IOField::getMetaData("MyField");

Array2DI4 MyFieldData = OMEGA::IOField::getData<Array2DI4>("MyField");
```
where Array2DI4 is just an example of a supported array type. Note that
all the array type copies are just shallow copies that copy pointers
and adjust the reference counts for the field.
As is probably clear from the calls above, the data attach/getData functions
are implemented as templates for each of the supported Array data types.
Similarly, the MetaData are represented by shared pointers that also
enable reference counting.

Additional utility functions can erase (remove) an IOField or clear all
IOFields (this must be done before the program exits):
```c++
OMEGA::IOField::erase("MyFieldName");
OMEGA::IOField::clear();
```
Finally, the class includes a query function to determine whether a
field has been defined:
```c++
OMEGA::IOField::isDefined("MyField");
```
