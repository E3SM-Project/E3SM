.. _input-namelists:

Input Namelists
===============

Each data model has two namelist groups in its input namelist file: a **stream-dependent** and a **stream-independent** namelist group.

The stream-dependent namelist group (``shr_strdata_nml``) specifies the data model mode, stream description text files, and interpolation options.
The stream description files will be provided as separate input files and contain the files and fields that need to be read.
The stream-independent namelist group (one of ``[datm_nml, dice_nml, dlnd_nml, docn_nml, drof_nml, dwav_nml]``) contains namelist input such as the data model decomposition, etc.

For users wanting to introduce new data sources for any data model, it is important to know what modes are supported and the internal field names in the data model.
That information will be used in the ``shr_strdata_nml`` namelist and stream input files.

Users will primarily set up different data model configurations through namelist settings.
**The stream input options and format are identical for all data models**.
The data model-specific namelist has significant overlap between data models, but each data model has a slightly different set of input namelist variables and each model reads that namelist from a unique filename.
The detailed namelist options for each data model will be described later, but each model will specify a filename or filenames for stream namelist input and each ``shr_strdata_nml`` namelist will specify a set of stream input files.

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

As mentioned above, the ``dataMode`` namelist variable that is associated with each data model specifies if there are any additional operations that need to be performed on that data model's input streams before return to the driver.
At a minimum, all data models support ``datamode`` values of ``NULL`` and ``COPYALL``.

- ``NULL`` - turns off the data model as a provider of data to the coupler.

- ``COPYALL`` - copies all fields directly from the input data streams. Any required fields not found on an input stream will be set to zero.

Three stream description files (see :ref:`input streams<input-streams>`) are then expected to be available, ``streama``, ``streamb`` and ``streamc``.
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
