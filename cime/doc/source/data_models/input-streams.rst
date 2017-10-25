.. _input-streams:

Input Streams
=============

--------
Overview
--------
An *input data stream* is a time-series of input data files where all the fields in the stream are located in the same data file and all share the same spatial and temporal coordinates (ie. are all on the same grid and share the same time axis). Normally a time axis has a uniform dt, but this is not a requirement.

The data models can have multiple input streams.

The data for one stream may be all in one file or may be spread over several files. For example, 50 years of monthly average data might be contained all in one data file or it might be spread over 50 files, each containing one year of data.

The data models can *loop* over stream data -- i.e., repeatedly cycle over some subset of an input stream's time axis. When looping, the models can only loop over whole years. For example, an input stream might have SST data for years 1950 through 2000, but a model could loop over the data for years 1960 through 1980. A model *cannot* loop over partial years, for example, from 1950-Feb-10 through 1980-Mar-15.

The input data must be in a netcdf file and the time axis in that file must be CF-1.0 compliant.

There are two main categories of information that the data models need to know about a stream:

- data that describes what a user wants -- what streams to use and how to use them -- things that can be changed by a user.
- data that describes the stream data -- meta-data about the inherent properties of the data itself -- things that cannot be changed by a user.

Generally, information about what streams a user wants to use and how to use them is input via the strdata ("stream data") Fortran namelist, while meta-data that describes the stream data itself is found in an xml-like text file called a "stream description file."

--------------------------------------------------
Stream Data and shr_strdata_nml namelists
--------------------------------------------------
The stream data (referred to as *strdata*) input is set via a Fortran namelist called ``shr_strdata_nml``.
That namelist, the associated strdata datatype, and the methods are contained in the share source code file, ``shr_strdata_mod.F90``.
In general, strdata input defines an array of input streams and operations to perform on those streams.
Therefore, many namelist inputs are arrays of character strings.
Different variables of the same index are associated. For instance, mapalgo(1) spatial interpolation will be performed between streams(1) and the target domain.

Each data model has an associated input namelist file, ``xxx_in``, where ``xxx=[datm,dlnd,dice,docn,drof,dwav]``.

The input namelist file for each data model has a stream dependent namelist group, ``shr_strdata_nml``, and a stream independent namelist group.
The  ``shr_strdata_nml`` namelist variables **are the same for all data models**.

=========== ==========================================================================================================================
File        Namelist Groups
=========== ==========================================================================================================================
datm_in     datm_nml, shr_strdata_nml
dice_in     dice_nml, shr_strdata_nml
dlnd_in     dlnd_nml, shr_strdata_nml
docn_in     docn_nml, shr_strdata_nml
drof_in     drof_nml, shr_strdata_nml
dwav_in     dwav_nml, shr_strdata_nml
=========== ==========================================================================================================================

.. _shr-strdata-nml:

The following table summaries the ``shr_strdata_nml`` entries.

=========== ==========================================================================================================================
Namelist    Description
=========== ==========================================================================================================================
dataMode    component specific mode.

            Each CIME data model has its own datamode values as described below:

            :ref:`datm dataMode<datm-datamodes>`

	    :ref:`dice dataMode<dice-datamodes>`

	    :ref:`dlnd dataMode<dlnd-datamodes>`

	    :ref:`docn dataMode<docn-datamodes>`

	    :ref:`drof dataMode<drof-datamodes>`

	    :ref:`dwav dataMode<dwav-datamodes>`

domainFile  component domain (all streams will be mapped to this domain).

            Spatial gridfile associated with the strdata.  grid information will
	    be read from this file and that grid will serve as the target grid
	    for all input data for this strdata input.
	    If the value is **null** then the domain of the first stream
	    will be used as the component domain

	    default="null"

streams     character array (up to 30 elements) of input stream filenames and associated years of data.

            Each array entry consists of a stream_input_filename year_align year_first year_last.
            The  stream_input_filename is a stream text input file and the format and options are described elsewhere.
	    The year_align, year_first, and year_last provide information about  the time axis of the file and how to relate
	    the input time axis to the model time axis.

	    default="null".

fillalgo    array (up to 30 elements) of fill algorithms associated with the array of streams.

            Valid options are just copy (ie. no fill), special value, nearest neighbor, nearest neighbor in "i" direction,
	    or nearest neighbor in "j" direction.

	    valid values:  'copy','spval','nn','nnoni','nnonj'

	    default value='nn'

fillmask    array (up to 30 elements) of fill masks.

            valid values: "nomask,srcmask,dstmask,bothmask"

            default="nomask"

fillread    array (up to 30 elements) fill mapping files to read. Secifies the weights file to read in instead of
            computing the weights on the fly for the fill operation.  If this is set, fillalgo and fillmask are ignored.

	    default='NOT_SET'

fillwrite   array of fill mapping file to write

	    default='NOT_SET'

mapalgo     array of spatial interpolation algorithms

            default="bilinear"

mapmask     array of spatial interpolation mask

	    default='NOT_SET'

mapread     array of spatial interpolation mapping files to read (optional)

	    default='NOT_SET'

mapwrite    array (up to 30 elements) of spatial interpolation mapping files to write (optional). Specifies the weights file
            to generate after weights are computed on the fly for the mapping (interpolation) operation, thereby allowing
	    users to save and reuse a set of weights later.
	    default='NOT_SET'

tintalgo    array (up to 30 elements) of time interpolation algorithm options associated with the array of streams.

            valid values: lower,upper,nearest,linear,coszen
            lower   = Use lower time-value

            upper   = Use upper time-value

            nearest = Use the nearest time-value

            linear  = Linearly interpolate between the two time-values

            coszen  = Scale according to the cosine of the solar zenith angle (for solar)

            default="linear"

taxMode     array (up to 30 elements) of time interpolation modes.

            Time axis interpolation modes are associated with the array of streams for
	    handling data outside the specified stream time axis.

	    Valid options are to cycle the data based on the first, last, and align
	    settings associated with the stream dataset, to extend the first and last
	    valid value indefinitely, or to limit the interpolated data to fall only between
	    the least and greatest valid value of the time array.

	    valid values: cycle,extend,limit

	    extend = extrapolate before and after the period by using the first or last value.

	    cycle  = cycle between the range of data

	    limit  = restrict to the period for which the data is valid

	    default="cycle"

dtlimit     array (up to 30 elements) of setting delta time axis limit.

            Specifies delta time ratio limits placed on the time interpolation
	    associated with the array of streams.  Causes the model to stop if
	    the ratio of the running maximum delta time divided by the minimum delta time
	    is greater than the dtlimit for that stream.  For instance, with daily data,
	    the delta time should be exactly one day throughout the dataset and
	    the computed maximum divided by minimum delta time should always be 1.0.
	    For monthly data, the delta time should be between 28 and 31 days and the
            maximum ratio should be about 1.1.  The running value of the delta
            time is computed as data is read and any wraparound or cycling is also
            included.  this input helps trap missing data or errors in cycling.
            to turn off trapping, set the value to 1.0e30 or something similar.

            default=1.5

vectors     paired vector field names
=========== ==========================================================================================================================


``shr_strdata_nml`` contains a namelist variable, ``streams``, that specifies a list of input stream description files and for each file what years of data to use, and how to align the input stream time axis with the model run time axis.

The general input format for the ``streams`` namelist variable is:
::

   &shr_strdata_nml
     streams = 'stream1.txt year_align year_first year_last ',
               'stream2.txt year_align year_first year_last ',
                ...
               'streamN.txt year_align year_first year_last '
   /

where:

.. code-block:: none

   streamN.txt
      the stream description file, a plain text file containing details about the input stream (see below)
   year_first
      the first year of data that will be used
   year_last
      the last year of data that will be used
   year_align
      a model year that will be aligned with data for year_first

---------------------
Details on year_align
---------------------

The ``year_align`` value gives the simulation year corresponding to
``year_first``. A common usage is to set this to the year of
``RUN_STARTDATE``. With this setting, the forcing in the first year of
the run will be the forcing of year ``year_first``. Another use case is
to align the calendar of transient forcing with the model calendar. For
example, setting ``year_align`` = ``year_first`` will lead to the
forcing calendar being the same as the model calendar. The forcing for a
given model year would be the forcing of the same year. This would be
appropriate in transient runs where the model calendar is setup to span
the same year range as the forcing data.

For some data model modes, ``year_align`` can be set via an xml variable
whose name ends with ``YR_ALIGN`` (there are a few such xml variables,
each pertaining to a particular data model mode).

An example of this is land-only historical simulations in which we run
the model for 1850 to 2010 using atmospheric forcing data that is only
available for 1901 to 2010. In this case, we want to run the model for
years 1850 (so ``RUN_STARTDATE`` has year 1850) through 1900 by looping
over the forcing data for 1901-1920, and then run the model for years
1901-2010 using the forcing data from 1901-2010. To do this, we
initially set::

  ./xmlchange DATM_CLMNCEP_YR_ALIGN=1901
  ./xmlchange DATM_CLMNCEP_YR_START=1901
  ./xmlchange DATM_CLMNCEP_YR_END=1920

When the model has completed year 1900, then we set::

  ./xmlchange DATM_CLMNCEP_YR_ALIGN=1901
  ./xmlchange DATM_CLMNCEP_YR_START=1901
  ./xmlchange DATM_CLMNCEP_YR_END=2010

With this setup, the correlation between model run year and forcing year
looks like this::

  RUN   Year : 1850 ... 1860 1861 ... 1870 ... 1880 1881 ... 1890 ... 1900 1901 ... 2010
  FORCE Year : 1910 ... 1920 1901 ... 1910 ... 1920 1901 ... 1910 ... 1920 1901 ... 2010

Setting ``DATM_CLMNCEP_YR_ALIGN`` to 1901 tells the code that you want
to align model year 1901 with forcing data year 1901, and then it
calculates what the forcing year should be if the model starts in
year 1850.

--------------------------------------------------
Customizing shr_strdata_nml values
--------------------------------------------------

The contents of ``shr_strdata_nml`` are automatically generated by that data model's **cime_config/buildnml** script.
These contents are easily customizable for your target experiment.
As an example we refer to the following ``datm_in`` contents (that would appear in both ``$CASEROOT/CaseDocs`` and ``$RUNDIR``):
::

   \&shr_strdata_nml
      datamode   = 'CLMNCEP'
      domainfile = '/glade/proj3/cseg/inputdata/share/domains/domain.lnd.fv1.9x2.5_gx1v6.090206.nc'
      dtlimit    = 1.5,1.5,1.5,1.5
      fillalgo   = 'nn','nn','nn','nn'
      fillmask   = 'nomask','nomask','nomask','nomask'
      mapalgo    = 'bilinear','bilinear','bilinear','bilinear'
      mapmask    = 'nomask','nomask','nomask','nomask'
      streams    = "datm.streams.txt.CLM_QIAN.Solar  1895 1948 1972  ",
                   "datm.streams.txt.CLM_QIAN.Precip 1895 1948 1972  ",
                   "datm.streams.txt.CLM_QIAN.TPQW   1895 1948 1972  ",
                   "datm.streams.txt.presaero.trans_1850-2000 1849 1849 2006"
      taxmode    = 'cycle','cycle','cycle','cycle'
      tintalgo   = 'coszen','nearest','linear','linear'
      vectors    = 'null'
   /


As is discussed in the :ref:`CIME User's Guide<running-a-case>`, to change the contents of ``datm_in``, you must edit ``$CASEROOT/user_nl_datm``.
In the above example, you can to this to change any of the above settings **except for the names**

.. code-block:: none

   datm.streams.txt.CLM_QIAN.Solar
   datm.streams.txt.CLM_QIAN.Precip
   datm.streams.txt.CLM_QIAN.TPQW
   datm.streams.txt.presaero.trans_1850-2000

Other than these names, any namelist variable from ``shr_strdata_nml`` can be modified by adding the appropriate keyword/value pairs to ``user_nl_datm``.

As an example, the following could be the contents of ``$CASEROOT/user_nl_datm``:
::

   !------------------------------------------------------------------------
   ! Users should ONLY USE user_nl_datm to change namelists variables
   ! Users should add all user specific namelist changes below in the form of
   ! namelist_var = new_namelist_value
   ! Note that any namelist variable from shr_strdata_nml and datm_nml can
   ! be modified below using the above syntax
   ! User preview_namelists to view (not modify) the output namelist in the
   ! directory $CASEROOT/CaseDocs
   ! To modify the contents of a stream txt file, first use preview_namelists
   ! to obtain the contents of the stream txt files in CaseDocs, and then
   ! place a copy of the  modified stream txt file in $CASEROOT with the string
   ! user_ prepended.
   !------------------------------------------------------------------------
   streams    = "datm.streams.txt.CLM_QIAN.Solar  1895 1948 1900  ",
                "datm.streams.txt.CLM_QIAN.Precip 1895 1948 1900  ",
                "datm.streams.txt.CLM_QIAN.TPQW   1895 1948 1900  ",
                "datm.streams.txt.presaero.trans_1850-2000 1849 1849 2006"

and the contents of ``shr_strdata_nml`` (in both ``$CASEROOT/CaseDocs`` and ``$RUNDIR``) would be
::

   datamode   = 'CLMNCEP'
   domainfile = '/glade/proj3/cseg/inputdata/share/domains/domain.lnd.fv1.9x2.5_gx1v6.090206.nc'
   dtlimit    = 1.5,1.5,1.5,1.5
   fillalgo   = 'nn','nn','nn','nn'
   fillmask   = 'nomask','nomask','nomask','nomask'
   mapalgo    = 'bilinear','bilinear','bilinear','bilinear'
   mapmask    = 'nomask','nomask','nomask','nomask'
   streams    = "datm.streams.txt.CLM_QIAN.Solar  1895 1948 1900  ",
                "datm.streams.txt.CLM_QIAN.Precip 1895 1948 1900  ",
                "datm.streams.txt.CLM_QIAN.TPQW   1895 1948 1900  ",
                "datm.streams.txt.presaero.trans_1850-2000 1849 1849 2006"
   taxmode    = 'cycle','cycle','cycle','cycle'
   tintalgo   = 'coszen','nearest','linear','linear'
   vectors    = 'null'

As is discussed in the :ref:`CIME User's Guide<running-a-case>`, you should use **preview_namelists** to view (not modify) the output namelist in ``CaseDocs``.


.. _stream_description_file:

-----------------------
Stream Description File
-----------------------
The *stream description file* is not a Fortran namelist, but a locally built xml-like parsing implementation.
Sometimes it is called a "stream dot-text file" because it has a ".txt." in the filename.
Stream description files contain data that specifies the names of the fields in the stream, the names of the input data files, and the file system directory where the data files are located.

The data elements found in the stream description file are:

``dataSource``
  A comment about the source of the data -- always set to GENERIC and is there only for backwards compatibility.

``domainInfo``
  Information about the domain data for this stream specified by the following 3 sub elements.

  ``variableNames``
    A list of the domain variable names. This is a paired list with the name of the variable in the netCDF file on the left and the name of the corresponding model variable on the right. This data models require five variables in this list. The names of model's variables (names on the right) must be: "time," "lon," "lat," "area," and "mask."

  ``filePath``
    The file system directory where the domain data file is located.

  ``fileNames``
    The name of the domain data file. Often the domain data is located in the same file as the field data (above), in which case the name of the domain file could simply be the name of the first field data file. Sometimes the field data files don't contain the domain data required by the data models, in this case, one new file can be created that contains the required data.


``fieldInfo``
  Information about the stream data for this stream specified by the following 3 required sub elements and optional offset element.

  ``variableNames``
    A list of the field variable names. This is a paired list with the name of the variable in the netCDF file on the left and the name of the corresponding model variable on the right. This is the list of fields to read in from the data file; there may be other fields in the file which are not read in (i.e., they won't be used).

  ``filePath``
    The file system directory where the data files are located.

  ``fileNames``
    The list of data files to use. If there is more than one file, the files must be in chronological order, that is, the dates in time axis of the first file are before the dates in the time axis of the second file.

  ``offset``
    The offset allows a user to shift the time axis of a data stream by a fixed and constant number of seconds. For instance, if a data set contains daily average data with timestamps for the data at the end of the day, it might be appropriate to shift the time axis by 12 hours so the data is taken to be at the middle of the day instead of the end of the day. This feature supports only simple shifts in seconds as a way of correcting input data time axes without having to modify the input data time axis manually. This feature does not support more complex shifts such as end of month to mid-month. But in conjunction with the time interpolation methods in the strdata input, hopefully most user needs can be accommodated with the two settings. Note that a positive offset advances the input data time axis forward by that number of seconds.

The data models advance in time discretely.
At a given time, they read/derive fields from input files.
Those input files have data on a discrete time axis as well.
Each data point in the input files is associated with a discrete time (as opposed to a time interval).
Depending on whether you pick lower, upper, nearest, linear, or coszen, the data in the input file will be "interpolated" to the time in the model.

The offset shifts the time axis of the input data the given number of seconds.
So if the input data is at 0, 3600, 7200, 10800 seconds (hourly) and you set an offset of 1800, then the input data will be set at times 1800, 5400, 9000, and 12600.
So a model at time 3600 using linear interpolation would have data at "n=2" with offset of 0 will have data at "n=(2+3)/2" with an offset of 1800.
n=2 is the 2nd data in the time list 0, 3600, 7200, 10800 in this example.
n=(2+3)/2 is the average of the 2nd and 3rd data in the time list 0, 3600, 7200, 10800.
offset can be positive or negative.

Actual example:
::

   <stream>
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
    <fieldInfo>
      <variableNames>
         dn10  dens
         slp_  pslv
         q10   shnum
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
   </stream>

--------------------------------------------------
Customizing stream description files
--------------------------------------------------

Each data model's **cime-config/buildnml** utility automatically generates the required stream description files for the case.
The directory contents of each data model will look like the following (using DATM as an example)
::

   $CIMEROOT/components/data_comps/datm/cime_config/buildnml
   $CIMEROOT/components/data_comps/datm/cime_config/namelist_definition_datm.xml

The ``namelist_definition_datm.xml`` file defines and sets default values for all the namelist variables and associated groups and also provides out-of-the box settings for the target data model and target stream.
**buildnml** utilizes these two files to construct the stream files for the given compset settings. You can modify the generated stream files for your particular needs by doing the following:


1. Copy the relevant description file from ``$CASEROOT/CaseDocs`` to ``$CASEROOT`` and pre-pend a "\user_"\ string to the filename. Change the permission of the file to write. For example, assuming you are in **$CASEROOT**
   ::

      cp $CASEROOT/CaseDocs/datm.streams.txt.CLM_QIAN.Solar  $CASEROOT/user_datm.streams.txt.CLM_QIAN.Solar
      chmod u+w $CASEROOT/user_datm.streams.txt.CLM_QIAN.Solar

2. Edit ``$CASEROOT/user_datm.streams.txt.CLM_QIAN.Solar`` with your desired changes.

   - *Be sure not to put any tab characters in the file: use spaces instead*.

   - In contrast to other user_nl_xxx files, be sure to set all relevant data model settings in the xml files, issue the **preview_namelist** command and THEN edit the ``user_datm.streams.txt.CLM_QIAN.Solar`` file.

   - **Once you have created a user_xxx.streams.txt.* file, further modifications to the relevant data model settings in the xml files will be ignored.**

   - If you later realize that you need to change some settings in an xml file, you should remove the user_xxx.streams.txt.* file(s), make the modifications in the xml file, rerun **preview_namelists**, and then reintroduce your modifications into a new user_xxx.streams.txt.* stream file(s).

3. Call **preview_namelists** and verify that your changes do indeed appear in the resultant stream description file appear in ``CaseDocs/datm.streams.txt.CLM_QIAN.Solar``. These changes will also appear in ``$RUNDIR/datm.streams.txt.CLM_QIAN.Solar``.
