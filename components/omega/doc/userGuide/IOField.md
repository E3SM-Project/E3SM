(omega-user-iofield)=

## IO Fields (IOField)

Within Omega, each module registers which fields are available for
input and output.  The user then selects which fields they wish to
include via the streams section of the input configuration. This is
described further in the [IOStreams](#omega-user-iostreams) section.
An IOField contains pointers to both the metadata and the data arrays
associated with the field. There are no user-configurable options
but further details on the methods and implementation can be found
in the [Developer's Guide](#omega-dev-iofield).
