(omega-user-field)=

## Fields and Metadata (Field)

Within OMEGA, each module registers fields that it owns together with the
metadata associated with that field. Metadata must be compliant with the
Climate and Forecast [(CF) metadata conventions.](http://cfconventions.org/)
This information can then be used by IOStreams or other modules that may
require metadata information. Two special fields are defined internally as
CodeMeta and FieldMeta to hold any global metadata not associated with a
data field, but may be required for provenance. Required metadata are enforced
by the Field create method but additional metadata can be added after field
creation. The data array holding field values can be attached/retrieved and
updated whenever needed.

Groups of fields can also be defined to provide a shortcut to a list of fields
commonly used together, like the model state or tracer groups. This reduces the
size of the input YAML configuration. The user can select which fields or field
groups they wish to include via the streams section of the input YAML
configuration file. This is described further in {ref} `omega-user-iostreams`.
There are no user-configurable options for the Field class itself since all
fields are defined and updated internally, but further details on the methods
and implementation can be found in the {ref}`omega-dev-field` section of the
Developer's Guide.
