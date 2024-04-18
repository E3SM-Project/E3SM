User's Guide
============

Guidance for using E3SM is available from [E3SM's public web site](<https://e3sm.org/model/running-e3sm/e3sm-quick-start>).

**Configuring MPAS-seaice**
---------------------------

MPAS-seaice is controlled using namelist options.

- Default namelist values are found in
``E3SM/components/mpas-seaice/bld/namelist_files/namelist_defaults_mpassi.xml``.
- Namelist options are defined in
``E3SM/components/mpas-seaice/bld/namelist_files/namelist_definitions_mpassi.xml``,
including type, category (``seaice_model``), group, valid values and a brief description. Each namelist variable is defined in an ``entry`` element.  The content of the element is the documentation of how the variable is used.  Other aspects of the variable's definition are expressed as attributes of the ``entry`` element.
- Some namelist values or combinations are not allowed and will generate warnings and often abort the code.  The consistency checks for using MPAS-seaice within E3SM are in ``mpas_seaice_initialize`` (subroutines ``seaice_check_configs_coupled``, ``seaice_check_constants_coupled``), and those specific to Icepack can be found in subroutine ``check_column_package_configs`` in ``mpas_seaice_icepack.F``.

Related namelist variables are grouped according to their application.

| Namelist Groups     | Relevant application |
| ------------------- | -------------------- |
| ``seaice_model``    |  general options     |
| ``io``              |  input/output        |
| ``decomposition``   | mesh parallelization |
| ``restart``         | restarting the code  |
| ``dimensions``      | column physics dimensions (layers, categories)  |
| ``initialize``      | initialization       |
| ``use_sections``    | turn entire parameterizations on and off |
| ``forcing``         | forcing for standalone configurations    |
| ``velocity_solver`` | algorithms for solving the dynamics (velocity and stress) equations |
| ``advection``       | advection                                |
| ``column_package``  | general column package software configurations  |
| ``biogeochemistry`` | biogeochemistry                          |
| ``shortwave``       | radiation                                |
| ``snow``            | advanced snow physics                    |
| ``meltponds``       | melt pond parameterization flags and parameters |
| ``thermodynamics``  | basic thermodynamics                     |
| ``itd``             | ice thickness distribution               |
| ``floesize``        | floe size distribution                   |
| ``ridging``         | mechanical redistribution                |
| ``atmosphere``      | atmospheric boundary layer and coupling  |
| ``ocean``           | oceanic boundary layer and coupling      |
| ``diagnostics``     | diagnostic output                        |
| ``prescribed_ice``  | for testing atmosphere simulations       |

**Icepack**
-----------

The Icepack software has replaced the original ``colpkg`` column physics code in MPAS-seaice. The ``column_package`` option is still available but is no longer being supported in MPAS-seaice.

Full documentation for E3SM's version of Icepack can be found in [E3SM's Icepack readthedocs](<https://e3sm-icepack.readthedocs.io/en/latest>).  The most up-to-date documentation from the CICE Consortium's main Icepack repository is [here](<https://cice-consortium-icepack.readthedocs.io/en/main>).

The MPAS-seaice driver for Icepack is

``E3SM/components/mpas-seaice/src/shared/mpas_seaice_icepack.f``

and the mapping between the names of Icepack's namelist options and those in MPAS-seaice can be found in subroutine ``init_icepack_package_configs`` (see the argument list for ``call subroutine icepack_init_parameters`` and comments at the end of ``init_icepack_package_configs``.

**Configuring Model Input and Output**
--------------------------------------

The reading and writing of model fields in MPAS is handled by user-configurable streams. A stream represents a fixed set of model fields, together with dimensions and attributes, that are all written or read together to or from the same file or set of files. Each MPAS model core may define its own set of default streams that it typically uses for reading initial conditions, for writing and reading restart fields, and for writing additional model history fields. Besides these default streams, users may define new streams to, e.g., write certain diagnostic fields at a higher temporal frequency than the usual model history fields.

Streams are defined in XML configuration files that are created at build time for each model core. The name of this XML file is simply ‘streams.’ suffixed with the name of the core. For example, the streams for the sw (shallow-water) core are defined in a file named ‘streams.sw’. An XML stream file may further reference other text files that contain lists of the model fields that are read or written in each of the streams defined in the XML stream file.

Changes to the XML stream configuration file will take effect the next time an MPAS core is run; there is no need to re-compile after making modifications to the XML files. As described in the next section, it is therefore possible, e.g., to change the interval at which a stream is written, the template for the filenames associated with a stream, or the set of fields that are written to a stream, without the need to re-compile any code.

Two classes of streams exist in MPAS: immutable streams and mutable streams. Immutable streams are those for which the set of fields that belong to the stream may not be modified at model run-time; however, it is possible to modify the interval at which the stream is read or written, the filename template describing the files containing the stream on disk, and several other parameters of the stream. In contrast, all aspects of mutable streams, including the set of fields that belong to the stream, may be modified at run-time. The motivation for the creation of two stream classes is the idea that an MPAS core may not function correctly if certain fields are not read in upon model start-up or written to restart files, and it is therefore not reasonable for users to modify this set of required fields at run-time. An MPAS core developer may choose to implement such streams as immutable streams. Since fields may not be added to an immutable stream at run-time, new immutable streams may not be defined at run-time, and the only type of new stream that may be defined at run-time is the mutable stream type.
