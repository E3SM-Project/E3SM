(omega-user-IO)=

## Parallel IO (IO)

For input and output of data needed for Omega, we use the Software for
Caching Output and Reads for Parallel I/O
([SCORPIO](https://github.com/E3SM-Project/scorpio)) library. This
library supports parallel reading and writing of distributed arrays in various
self-describing
formats like [NetCDF](https://docs.unidata.ucar.edu/netcdf/),
[HDF5](https://www.hdfgroup.org/solutions/hdf5/),
and [ADIOS](https://csmd.ornl.gov/adios).
SCORPIO is responsible for rearranging data from data decomposition used
for computation to a different rearrangement for parallel IO. In particular,
for optimal IO, the user would specify a different number of IO tasks than
for computation. For example, the user could specify the number of IO tasks to
match underlying hardware like the network interfaces (NICs) on a node.
Users will generally access IO via the [IOStreams](#omega-user-iostreams)
and other interfaces and the base IO layer described here only provides
an Omega-aware wrapper around SCORPIO calls.

There are some general parameters that must be set for IO performance and
formatting via the input configuration file. These are:
```yaml
IO:
    IOTasks:  1
    IOStride: 1
    IOBaseTask: 0
    IORearranger: box
    IODefaultFormat: pnetcdf
```
where ``IOTasks`` is the total number of IOTasks to assign to reading
and writing. The default is 1 (serial IO) for safety but this number
should be set appropriately for the underlying hardware. A simple
starting point might be one IOTask per node or socket. The ``IOStride``
is set to spread the IOTasks across the total number of MPI tasks so
that every IOStride task (starting with the root task IOBaseTask) is an
IOTask. The product of IOTasks and IOStride should equal the total number
of MPI Tasks.

These options apply to standalone Omega runs, where Omega creates its own
SCORPIO IO system. In coupled E3SM runs the IO system is created by the
coupled driver, so Omega reuses the driver's instance rather than
creating its own. The above IO options cannot be overridden in
``user_nl_omega`` and are set to ``-9999`` in the generated ``omega.yml``
to show they are not used.

Instead, the IO task layout, rearranger and file format are set by case
variables in ``env_run.xml`` in the case directory:

| Omega option          | Case variable             |
|-----------------------|---------------------------|
| ``IOTasks``           | ``PIO_NUMTASKS_OCN``      |
| ``IOStride``          | ``PIO_STRIDE_OCN``        |
| ``IOBaseTask``        | ``PIO_ROOT_OCN``          |
| ``IORearranger``      | ``PIO_REARRANGER_OCN``    |
| ``IODefaultFormat``   | ``PIO_TYPENAME_OCN``      |
| (netCDF file version) | ``PIO_NETCDF_FORMAT_OCN`` |

Each variable has a separate value for each model component, and Omega uses
the ocean (``OCN``) value. These can be changed with ``./xmlchange``, for
example ``./xmlchange PIO_STRIDE_OCN=56``. Searching ``env_run.xml`` for a
variable name requires omitting the component suffix (e.g. ``PIO_STRIDE``),
since the component is given as an attribute of each value.

In coupled runs, ``PIO_NETCDF_FORMAT_OCN`` sets the netCDF file version
(for example classic, ``64bit_offset`` or ``64bit_data``) for files written
with the ``netcdf`` and ``pnetcdf`` formats. It has no equivalent in the Omega
configuration, where ``pnetcdf`` files are always written with the
``64bit_data`` version.

When using parallel IO, the data must be rearranged to match the IO task
decomposition. There are two algorithms for rearranging data available
in SCORPIO: box and subset. Box is the default and preferred for most
configurations. It ensures each IO task has a contiguous chunk of data to
read/write but may draw from the full range of MPI tasks. In the subset
rearranger, each IO task is assigned a subset of the MPI tasks and only
communicates with those tasks. It no longer is guaranteed to have a
contiguous chunk of data to read/write. This option can sometimes result
in more optimal communication patterns but at the cost of less efficient
I/O as the underlying library must handle the non-contiguous data. The
subset option is available for exploring the most efficient approach.
See SCORPIO documentation for details.

Finally, the user can specify the file format. The SCORPIO library
supports all the various NetCDF formats, though the use of older
NetCDF formats is strongly discouraged. PNetCDF is currently the
default for Omega. In addition, SCORPIO supports
native HDF5 files (note that the latest NetCDF formats are all implemented
with HDF5 and are mostly compatible with HDF5) The ADIOS format is different
since ADIOS writes each chunk of data to a separate file and these files must
be combined later. E3SM provides the necessary tools to combine/convert
ADIOS files to netcdf as part of pre/post-processing. For standalone Omega,
users would need to become familiar with these tools if they wish to use
ADIOS.
