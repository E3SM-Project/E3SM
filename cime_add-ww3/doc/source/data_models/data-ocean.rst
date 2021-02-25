.. _data-ocean:

===================
Data Ocean (DOCN)
===================

Data ocean can be run both as a prescribed component, simply reading
in SST data from a stream, or as a prognostic slab ocean model
component.

The data ocean component (DOCN) always returns SSTs to the driver.
The atmosphere/ocean fluxes are computed in the coupler.  Therefore,
the data ocean model does not compute fluxes like the data ice (DICE)
model.  DOCN has two distinct modes of operation.  DOCN can run as a
pure data model, reading in ocean SSTs (normally climatological) from
input datasets, performing time/spatial interpolations, and passing
these to the coupler.  Alternatively, DOCN can compute updated SSTs by
running as a slab ocean model where bottom ocean heat flux convergence
and boundary layer depths are read in and used with the
atmosphere/ocean and ice/ocean fluxes obtained from the driver.

DOCN running in prescribed mode assumes that the only field in the
input stream is SST and also that SST is in Celsius and must be
converted to Kelvin.  All other fields are set to zero except for
ocean salinity, which is set to a constant reference salinity value.
Normally the ice fraction data (used for prescribed CICE) is found in
the same data files that provide SST data to the data ocean model
since SST and ice fraction data are derived from the same
observational data sets and are consistent with each other.  For DOCN
prescribed mode, default yearly climatological datasets are provided
for various model resolutions.

DOCN running as a slab ocean model is used in conjunction with active
ice mode running in full prognostic mode (e.g. CICE for CESM).  This
mode computes a prognostic sea surface temperature and a freeze/melt
potential (surface Q-flux) used by the sea ice model.  This
calculation requires an external SOM forcing data file that includes
ocean mixed layer depths and bottom-of-the-slab Q-fluxes.
Scientifically appropriate bottom-of-the-slab Q-fluxes are normally
ocean resolution dependent and are derived from the ocean model output
of a fully coupled CCSM run.  Note that this mode no longer runs out
of the box, the default testing SOM forcing file is not scientifically
appropriate and is provided for testing and development purposes only.
Users must create scientifically appropriate data for their particular
application or use one of the standard SOM forcing files from the full
prognostic control runs.  For CESM, some of these are available in the
`inputdata repository
<https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/ocn/docn7/SOM/>`_.
The user then modifies the ``$DOCN_SOM_FILENAME`` variable in
env_run.xml to point to the appropriate SOM forcing dataset.

.. note:: A tool is available to derive valid `SOM forcing <http://www.cesm.ucar.edu/models/cesm1.2/data8/doc/SOM.pdf>`_ and more information on creating the SOM forcing is also available.

.. _docn-xml-vars:

-------------
xml variables
-------------

The following are xml variables that CIME supports for DOCN.  These
variables are defined in
``$CIMEROOT/src/components/data_comps/docn/cime_config/config_component.xml``.
These variables will appear in ``env_run.xml`` and are used by the
DOCN ``cime_config/buildnml`` script to generate the DOCN namelist
file ``docn_in`` and the required associated stream files for the
case.

.. note:: These xml variables are used by the the docn's **cime_config/buildnml** script in conjunction with docn's **cime_config/namelist_definition_docn.xml** file to generate the namelist file ``docn_in``.

.. csv-table:: "DOCN xml variables"
   :header: "xml variable", "description"
   :widths: 15, 85

   "DOCN_MODE", "Data mode"
   "", "Valid values are: null, prescribed, som, interannual, ww3"
   "DOCN_SOM_FILENAME", "Sets SOM forcing data filename for pres runs, only used in D and E compset"
   "SSTICE_STREAM", "Prescribed SST and ice coverage stream name."
   "", "Sets SST and ice coverage stream name for prescribed runs."
   "SSTICE_DATA_FILENAME", "Prescribed SST and ice coverage data file name."
   "", "Sets SST and ice coverage data file name for DOCN prescribed runs."
   "SSTICE_YEAR_ALIGN", "The model year that corresponds to SSTICE_YEAR_START on the data file."
   "", "Prescribed SST and ice coverage data will be aligned so that the first year of"
   "", "data corresponds to SSTICE_YEAR_ALIGN in the model. For instance, if the first"
   "", "year of prescribed data is the same as the first year of the model run, this" 
   "", "should be set to the year given in RUN_STARTDATE."
   "", "If SSTICE_YEAR_ALIGN is later than the model's starting year, or if the model is"
   "", "run after the prescribed data ends (as determined by SSTICE_YEAR_END), the"
   "", "default behavior is to assume that the data from SSTICE_YEAR_START to SSTICE_YEAR_END"
   "", "cyclically repeats. This behavior is controlled by the *taxmode* stream option"
   "SSTICE_YEAR_START", "The first year of data to use from SSTICE_DATA_FILENAME."
   "", "This is the first year of prescribed SST and ice coverage data to use. For"
   "", "example, if a data file has data for years 0-99, and SSTICE_YEAR_START is 10,"
   "", "years 0-9 in the file will not be used."
   "SSTICE_YEAR_END", "The last year of data to use from SSTICE_DATA_FILENAME."
   "", "This is the last year of prescribed SST and ice coverage data to use. For"
   "", "example, if a data file has data for years 0-99, and value is 49,"
   "", "years 50-99 in the file will not be used."

.. note:: For multi-year runs requiring AMIP datasets of sst/ice_cov fields, you need to set the xml variables for ``DOCN_SSTDATA_FILENAME``, ``DOCN_SSTDATA_YEAR_START``, and ``DOCN_SSTDATA_YEAR_END``. CICE in prescribed mode also uses these values.

.. _docn-datamodes:

---------------
datamode values
---------------

The xml variable ``DOCN_MODE`` (described in :ref:`docn_mode`) sets the streams that are associated with DOCN and also sets the namelist variable ``datamode``.
``datamode`` (which appears in ``shr_strdata_nml``) specifies what additional operations need to be done by DOCN on the streams before returning to the driver.

Each data model has its own set of supported ``datamode`` values. The following are the supported DOCN ``datamode`` values, as defined in the file ``namelist_definition_docn.xml``.

.. csv-table:: "Valid values for datamode namelist variable"
   :header: "datamode variable", "description"
   :widths: 20, 80

   "NULL", "Turns off the data model as a provider of data to the coupler.  The ocn_present flag will be set to false and the coupler will assume no exchange of data to or from the data model."
   "COPYALL", "The default science mode of the data model is the COPYALL mode. This mode will examine the fields found in all input data streams; if any input field names match the field names used internally, they are copied into the export array and passed directly to the coupler without any special user code.  Any required fields not found on an input stream will be set to zero."
   "SSTDATA", "assumes the only field in the input stream is SST. It also assumes the SST is in Celsius and must be converted to Kelvin.  All other fields are set to zero except for ocean salinity, which is set to a constant reference salinity value. Normally the ice fraction data is found in the same data files that provide SST data to the data ocean model. They are normally found in the same file because the SST and ice fraction data are derived from the same observational data sets and are consistent with each other. They are normally found in the same file because the SST and ice fraction data are derived from the same observational data sets and are consistent with each other."
   "IAF", "is the interannually varying version of SSTDATA"
   "SOM", "(slab ocean model) mode is a prognostic mode.  This mode computes a prognostic sea surface temperature and a freeze/melt potential (surface Q-flux) used by the sea ice model.  This calculation requires an external SOM forcing data file that includes ocean mixed layer depths and bottom-of-the-slab Q-fluxes. Scientifically appropriate bottom-of-the-slab Q-fluxes are normally ocean resolution dependent and are derived from the ocean model output of a fully coupled CCSM run.  Note that while this mode runs out of the box, the default SOM forcing file is not scientifically appropriate and is provided for testing and development purposes only. Users must create scientifically appropriate data for their particular application.  A tool is available to derive valid SOM forcing."

.. _docn_mode:

-------------------------------
DOCN_MODE, datamode and streams
-------------------------------

The following table describes the valid values of ``DOCN_MODE`` (defined in the ``config_component.xml`` file for DOCN), and how they relate to the associated input streams and the ``datamode`` namelist variable.
CIME will generate a value of ``DOCN_MODE`` based on the compset.

.. csv-table:: "Relationship between DOCN_MODE, datamode and streams"
   :header: "DOCN_MODE, "description-streams-datamode"
   :widths: 20, 80

   "null", "null mode"
   "", "streams: none"
   "", "datamode: null"
   "prescribed","run with prescribed climatological SST and ice-coverage"
   "","streams:  prescribed"
   "","datamode: SSTDATA"
   "interannual", "run with interannual SST and ice-coverage"
   "","streams:  prescribed"
   "","datamode: SSTDATA"
   "som", "run in slab ocean mode"
   "","streams:  som"
   "","datamode: SOM"
   "ww3", "ww3 mode"
   "", "streams: ww3"
   "", "datamode: COPYALL"

.. _docn-namelists:

---------
Namelists
---------

As is the case for all data models, DOCN namelists can be separated into two groups, stream-independent and stream-dependent. 

The namelist file for DOCN is ``docn_in`` (or ``docn_in_NNN`` for multiple instances).

The stream dependent group is :ref:`shr_strdata_nml<input-streams>` .

As part of the stream dependent namelist input, DOCN supports two science modes, ``SSTDATA`` (prescribed mode) and ``SOM`` (slab ocean mode). 

.. _docn-stream-independent-namelists:

The stream-independent group is ``docn_nml`` and the DOCN stream-independent namelist variables are:

=====================  ======================================================
decomp                 decomposition strategy (1d, root)
    
                       1d => vector decomposition, root => run on master task
restfilm               master restart filename 
restfils               stream restart filename 
force_prognostic_true  TRUE => force prognostic behavior
=====================  ======================================================

To change the namelist settings in docn_in, edit the file user_nl_docn. 

.. _docn-mode-independent-streams:

---------------------------------
Datamode independent streams
---------------------------------

There are no datamode independent streams for DOCN.

.. _docn-fields:

-----------
Field names
-----------

DOCN defines a set of pre-defined internal field names as well as mappings for how those field names map to the fields sent to the coupler.

.. note:: In general, the stream input file should translate the stream input variable names into the ``docn_fld`` names below for use within the data ocn model.

.. csv-table:: "DOCN internal field names"
   :header: "docn_fld (avifld)", "driver_fld (avofld)"
   :widths: 30, 30

   "t", "So_t" 
   "u", "So_u" 
   "v", "So_v" 
   "dhdx", "So_dhdx" 
   "dhdy", "So_dhdy" 
   "s", "So_s" 
   "h", "strm_h (internal to docn_comp_mod only)"
   "qbot", "strm_qbot (internal to docn_comp_mod only)"

.. _creating-sstdata-input-from-prognostic-run:

---------------------------------------------------------------------
Creating SSTDATA mode input from a fully prognostic run (CESM only)
---------------------------------------------------------------------

The following outlines the steps you would take to create monthly averages of SST and ice coverage from a previous fully prognostic run that can then be read as as stream data by DOCN.

As an example, the following uses an f09_g16 CESM B-configuration simulation using CAM5 physics and with cosp enabled. The procedure to create the SST/ICE file is as follows:

1. Save monthly averaged 'aice' information from cice code (this is the default).

2. Save monthly averaged SST information from pop2. To do this, copy $SRCROOT/pop2/input_templates/gx1v6_tavg_contents to $CASEROOT/SourceMods/src.pop2 and change the 2 in front of SST to 1 for monthly frequency.

3. Extract (using ncrcat) SST from monthly pop2 history files and form a single netcdf file containing just SST; change SST to SST_cpl.
   ::

      > ncrcat -v SST case.pop.h.*.nc temp.nc
      > ncrename -v SST,SST_cpl temp.nc sst_cpl.nc

4. Extract aice from monthly cice history files and form a single netcdf file containing aice; change aice to ice_cov; divide values by 100 (to convert from percent to fraction).
   ::

      > ncrcat -v aice case.cice.h.*.nc temp.nc
      > ncrename -v aice,ice_cov temp.nc temp2.nc
      > ncap2 -s 'ice_cov=ice_cov/100.' temp2.nc ice_cov.nc

5. Modify fill values in the sst_cpl file (which are over land points) to have value -1.8 and remove fill and missing value designators; change coordinate lengths and names: to accomplish this, first run ncdump, then replace _ with -1.8 in SST_cpl, then remove lines with _FillValue and missing_value. 
   (Note: although it might be possible to merely change the fill value to -1.8, this is conforming to other SST/ICE files, which have SST_cpl explicitly set to -1.8 over land.) 
   To change coordinate lengths and names, replace nlon by lon, nlat by lat, TLONG by lon, TLAT by lat. 
   The last step is to run ncgen. Note: when using ncdump followed by ncgen, precision will be lost; however, one can specify -d 9,17 to maximize precision - as in the following example:
   ::

      > ncdump -d 9,17 old.nc > old
      > ncgen -o new.nc new

6. Modify fill values in the ice_cov file (which are over land points) to have value 1 and remove fill and missing value designators; change coordinate lengths and names; patch longitude and latitude to replace missing values.
   To accomplish this, first run ncdump, then replace _ with 1 in ice_cov, then remove lines with _FillValue and missing_value. 
   To change coordinate lengths and names, replace ni by lon, nj by lat, TLON by lon, TLAT by lat. 
   To patch longitude and latitude arrays, replace values of those arrays with those in sst_cpl file. 
   The last step is to run ncgen. 
   (Note: the replacement of longitude and latitude missing values by actual values should not be necessary but is safer.)

7. Combine (using ncks) the two netcdf files.
   ::

      > ncks -v ice_cov ice_cov.nc sst_cpl.nc

   Rename the file to ssticetemp.nc. 
   The time variable will refer to the number of days at the end of each month, counting from year 0, whereas the actual simulation began at year 1. 
   However, we want time values to be in the middle of each month, referenced to the first year of the simulation (first time value equals 15.5).
   Extract (using ncks) time variable from existing amip sst file (for correct number of months - 132 in this example) into working netcdf file.
   ::

      > ncks -d time,0,131 -v time amipsst.nc ssticetemp.nc

   Add date variable: ncdump date variable from existing amip sst file; modify first year to be year 0 instead of 1949 (do not including leading zeroes or it will interpret as octal) and use correct number of months; ncgen to new netcdf file; extract date (using ncks) and place in working netcdf file.
   ::

      > ncks -v date datefile.nc ssticetemp.nc

   Add datesec variable: extract (using ncks) datesec (correct number of months) from existing amip sst file and place in working netcdf file.
   ::

      > ncks -d time,0,131 -v datesec amipsst.nc ssticetemp.nc

8. At this point, you have an SST/ICE file in the correct format. 

9. Due to CAM's linear interpolation between mid-month values, you need to apply a procedure to assure that the computed monthly means are consistent with the input data. 
   To do this, invoke ``$SRCROOT/components/cam/tools/icesst/bcgen`` and following the following steps:

   a. Rename SST_cpl to SST, and ice_cov to ICEFRAC in the current SST/ICE file:
      ::

	 > ncrename -v SST_cpl,SST -v ice_cov,ICEFRAC ssticetemp.nc

   b. In driver.f90, sufficiently expand the lengths of variables prev_history and history (16384 should be sufficient); also comment out the test that the climate year be between 1982 and 2001 (lines 152-158).

   c. In bcgen.f90 and setup_outfile.f90, change the dimensions of xlon and ???TODO xlat to (nlon,nlat); this is to accommodate use of non-cartesian ocean grid.

   d. In setup_outfile.f90, modify the 4th and 5th ???TODO arguments in the calls to wrap_nf_def_var for *lon* and *lat* to be *2* and *dimids*; this is to accommodate use of non-cartesian ocean grid.

   e. Adjust Makefile to have proper path for LIB_NETCDF and INC_NETCDF.

   f. Modify namelist accordingly.

   g. Make bcgen and execute per instructions. The resulting sstice_ts.nc file is the desired ICE/SST file.

9. Place the new SST/ICE file in desired location and modify ``env_run.xml`` to have :

   a. ``SSTICE_DATA_FILENAME`` point to the complete path of your SST/ICE file.

   b. ``SSTICE_GRID_FILENAME`` correspond to full path of (in this case) gx1v6 grid file.

   c. ``SSTICE_YEAR_START`` set to 0

   d. ``SSTICE_YEAR_END`` to one less than the total number of years
      
   e. ``SSTICE_YEAR_ALIGN`` to 1 (for CESM, since CESM starts counting at year 1).
