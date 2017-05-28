.. _data-model-introduction:

Introduction
============

--------
Overview
--------
The CIME data models perform the basic function of reading external data files, modifying those data, and then sending the data to the driver via the CIME coupling interfaces.
The fields sent to the driver are the same as those that would be sent by an active component.
This takes advantage of the fact that the driver and other models have no fundamental knowledge of whether another component is fully active or just a data model.
So, for example, the data atmosphere model (datm) sends the same fields as the prognostic Community Atmosphere Model (CAM).
However, rather than determining these fields prognostically, most data models simply read prescribed data.

The data models typically read gridded data from observations or reanalysis products.
Out of the box, they often provide a few possible data sources and/or time periods that you can choose from when setting up a case.
However, data models can also be configured to read output from a previous coupled run.
For example, you can perform a fully-coupled run in which you ask for particular extra output streams; you can then use these saved "driver history" files as inputs to datm to run a later land-only spinup.

In some cases, data models have prognostic functionality, that is, they also receive and use data sent by the driver.
However, in most cases, the data models are not running prognostically and have no need to receive any data from the driver.

The CIME data models have parallel capability and share significant amounts of source code. 
Methods for reading and interpolating data have been established and can easily be reused:
The data model calls strdata ("stream data") methods which then call stream methods. 
The stream methods are responsible for managing lists of input data files and their time axis. 
The information is then passed up to the strdata methods where the data is read and interpolated in space and time. 
The interpolated data is passed up to the data model where final fields are derived, packed, and returned to the driver.

------
Design
------
Data models function by reading in different streams of input data and interpolating those data both spatially and temporally to the appropriate final model grid and model time. 

- **Each data model**

  - communicates with the driver with fields on only the data model model grid.

  - can be associated with multiple streams

  - is associated with only one datamode value (specified in the ``shr_strdata_nml`` namelist group)

  - has an xml variable in ``env_run.xml`` that specifies its mode. These are:
    ``DATM_MODE``, ``DICE_MODE``, ``DLND_MODE``, ``DOCN_MODE``, ``DROF_MODE``, ``DWAV_MODE``.

- **Each ``DXXX_MODE`` xml variable variable specfies 2 things:**

  - the list of streams that are associated with the data model.

  - a ``datamode`` namelist variable that is associated with each data model and that determines if additional operations need to be performed on on the input streams before returning to the driver.

    at a minimum, all data models support ``datamode`` values of ``NULL`` and ``COPYALL``.

    - ``NULL`` - turns off the data model as a provider of data to the coupler.

    - ``COPYALL`` - copies all fields directly from the input data streams. Any required fields not found on an input stream will be set to zero.
 
- **Each data model stream**

  - can be associated with multiple stream input files (specified in the ``shr_strdata_nml`` namelist group).

- **Each stream input file** 

  - can contain data on a unique grid and unique temporal time stamps.

  - is interpolated to a single model grid and the present model time.

More details of the data model design are covered in :ref:`design details<design-details>`. 

--------------
Namelist Input
--------------

Each data model has two namelist groups in its input namelist file: a **stream-dependent** and a **stream-independent** namelist group. 

The stream-dependent namelist group (``shr_strdata_nml``) specifies the data model mode, stream description text files, and interpolation options. 
The stream description files will be provided as separate input files and contain the files and fields that need to be read.
The stream-independent namelist group (one of ``[datm_nml, dice_nml, dlnd_nml, docn_nml, drof_nml, dwav_nml]``) contains namelist input such as the data model decomposition, etc.

From a user perspective, for any data model, it is important to know what modes are supported and the internal field names in the data model.
That information will be used in the strdata namelist and stream input files.

Users will primarily setup different data model configurations through namelist settings.
**The strdata and stream input options and format are identical for all data models**. 
The data model specific namelist has significant overlap between data models, but each data model has a slightly different set of input namelist variables and each model reads that namelist from a unique filename.
The detailed namelist options for each data model will be described later, but each model will specify a filename or filenames for strdata namelist input and each strdata namelist will specify a set of stream input files.

The following example illustrates the basic set of namelist inputs::

   &dlnd_nml
      decomp = '1d'
   /
   &shr_strdata_nml
     dataMode   = 'CPLHIST'
     domainFile = 'grid.nc'
     streams    = 'streama', 'streamb', 'streamc'
     mapalgo    = 'interpa', 'interpb', 'interpc'
   /

As mentioned above, the ``dataMode`` namelist variable that is associated with each data model specifies if there is any additional operations that need to be performed on that data model's input streams before return to the driver.
At a minimum, all data models support ``datamode`` values of ``NULL`` and ``COPYALL``.

- ``NULL`` - turns off the data model as a provider of data to the coupler.

- ``COPYALL`` - copies all fields directly from the input data streams. Any required fields not found on an input stream will be set to zero.

Three stream description files are then expected to be available, ``streama``, ``streamb`` and ``streamc``.
Those files specify the input data filenames, input data grids, and input fields that are expected, among other things. 
The stream files are **not** Fortran namelist format.
Their format and options will be described later.
As an example, one of the stream description files might look like
::

   <stream>
      <dataSource>
         GENERIC
      </dataSource>
      <fieldInfo>
         <variableNames>
            dn10  dens
            slp_  pslv
            q_10  shum
            t_10  tbot
            u_10  u
            v_10  v
         </variableNames>
         <filePath>
            /glade/proj3/cseg/inputdata/atm/datm7/NYF
         </filePath>
         <offset>
            0
         </offset>
         <fileNames>
            nyf.ncep.T62.050923.nc
         </fileNames>
      </fieldInfo>
      <domainInfo>
         <variableNames>
            time   time
            lon    lon
            lat    lat
            area   area
            mask   mask
         </variableNames>
         <filePath>
            /glade/proj3/cseg/inputdata/atm/datm7/NYF
         </filePath>
         <fileNames>
            nyf.ncep.T62.050923.nc
         </fileNames>
      </domainInfo>
   </stream>


In general, these examples of input files are not complete, but they do show the general hierarchy and feel of the data model input.

-------------
Next Sections
-------------
In the next sections, more details will be presented including a full description of the science modes and namelist settings for the data atmosphere, data land, data runoff, data ocean, and data ice models; namelist settings for the strdata namelist input; a description of the format and options for the stream description input files; and a list of internal field names for each of the data components.
The internal data model field names are important because they are used to setup the stream description files and to map the input data fields to the internal data model field names.


