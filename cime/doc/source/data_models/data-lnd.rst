.. _data-lnd:

Data Land (DLND)
================

The land model is unique because it supports land data and snow data (*lnd and sno*) almost as if they were two separate components, but they are in fact running in one component model through one interface.
The lnd (land) data consist of fields sent to the atmosphere.
This set of data is used when running DLND with an active atmosphere.
In general this is not a mode that is used or supported.
The sno (snow) data consist of fields sent to the glacier model. This set of data is used when running dlnd with an active glacier model (TG compsets). Both sets of data are assumed to be on the same grid.

.. _dlnd-xml-vars:

---------------
xml variables
---------------

The following are xml variables that CIME supports for DLND.
These variables are defined in ``$CIMEROOT/src/components/data_comps/dlnd/cime_config/config_component.xml``.
These variables will appear in ``env_run.xml`` and are used by the DLND ``cime_config/buildnml`` script to generate the DLND namelist file ``dlnd_in`` and the required associated stream files for the case.

.. note:: These xml variables are used by the the dlnd's **cime_config/buildnml** script in conjunction with dlnd's **cime_config/namelist_definition_dlnd.xml** file to generate the namelist file ``dlnd_in``.

.. csv-table:: "DLND xml variables"
   :header: "xml variable", "description"
   :widths: 15, 85

   "DLND_MODE", "Mode for data land component"
   "", "Valid values are: NULL, CPLHIST, GLC_CPLHIST"

   "DLND_CPLHIST_DOMAIN_FILE", "Coupler history forcing data mode - full pathname of model domain file"
   "DLND_CPLHIST_CASE", "Coupler history forcing data mode - case name"
   "DLND_CPLHIST_DIR", "Coupler history forcing data mode - directory containing coupler history data"
   "DLND_CPLHIST_YR_ALIGN",  "Coupler history forcing data mode - simulation year corresponding to DLND_CPLHIST_YR_START"
   "DLND_CPLHIST_YR_START", "Coupler history forcing data mode - starting year to loop data over"
   "DLND_CPLHIST_YR_END", "Coupler history forcing data mode - ending year to loop data over"

.. note:: If ``DLND_MODE`` is set to ``CPLHIST``, it is normally assumed that the model domain will be identical to **all** of the stream domains. To ensure this, the xml variables ``LND_DOMAIN_PATH`` and ``LND_DOMAIN_FILE`` are ignored and a valid setting **must be given** for ``DLND_CPLHIST_DOMAIN_FILE``. If ``DLND_CPLHIST_DOMAIN_FILE`` is set to ``null``, then the dlnd component domain information is read in from the first coupler history file in the target stream and  it is assumed that the first coupler stream file that is pointed to contains the domain  information for that stream. Alternatively, ``DLND_CPLHIST_DOMAIN_FILE`` can be set to ``$LND_DOMAIN_PATH/$LND_DOMAIN_FILE``.

.. _dlnd-datamodes:

--------------------
datamode values
--------------------

The xml variable ``DLND_MODE`` (described in :ref:`dlnd_mode`) sets the streams that are associated with DLND and also sets the namelist variable ``datamode``.
``datamode`` (which appears in ``shr_strdata_nml``) specifies what additional operations need to be done by DLND on the streams before returning to the driver.

Each data model has its own set of supported ``datamode`` values. The following are the supported DLND ``datamode`` values, as defined in the file ``namelist_definition_dlnd.xml``.

.. csv-table:: "Valid values for datamode namelist variable"
   :header: "datamode variable", "description"
   :widths: 20, 80

   "NULL", "Turns off the data model as a provider of data to the coupler.  The ``lnd_present`` flag will be set to false and the coupler will assume no exchange of data to or from the data model."
   "COPYALL", "The default science mode of the data model is the COPYALL mode. This mode will examine the fields found in all input data streams; if any input field names match the field names used internally, they are copied into the export array and passed directly to the coupler without any special user code.  Any required fields not found on an input stream will be set to zero."

.. _dlnd_mode:

-------------------------------
DLND_MODE, datamode and streams
-------------------------------

The following table describes the valid values of ``DLND_MODE`` (defined in the ``config_component.xml`` file for DLND), and how they relate to the associated input streams and the ``datamode`` namelist variable.
CIME will generate a value of ``DLND_MODE`` based on the compset.

.. csv-table:: "Relationship between DLND_MODE, datamode and streams"
   :header: "DLND_MODE", "description-streams-datamode"
   :widths: 20, 80

   "NULL", "null mode"
   "", "streams: none"
   "", "datamode: null"
   "CPLHIST", "land forcing data (e.g. produced by CESM/CLM) from a previous model run are read in from coupler history files"
   "", "streams: lnd.cplhist"
   "", "datamode: COPYALL"
   "GLC_CPLHIST", "glc coupling fields (e.g. produced by CESM/CLM) from a previous model run are read in from coupler history files"
   "", "streams: sno.cplhist"
   "", "datamode: COPYALL"

---------
Namelists
---------

The namelist file for DLND is ``dlnd_in`` (or ``dlnd_in_NNN`` for multiple instances).

As is the case for all data models, DLND namelists can be separated into two groups, stream-independent and stream-dependent.

The stream dependent group is :ref:`shr_strdata_nml<input-streams>`.

.. _dlnd-stream-independent-namelists:

The stream-independent group is ``dlnd_nml`` and the DLND stream-independent namelist variables are:

=====================  ======================================================
decomp                 decomposition strategy (1d, root)

                       1d => vector decomposition, root => run on master task
restfilm               master restart filename
restfils               stream restart filename
force_prognostic_true  TRUE => force prognostic behavior
=====================  ======================================================

To change the namelist settings in ``dlnd_in``, edit the file ``user_nl_dlnd``.

.. _dlnd-mode-independent-streams:

--------------------------------------
Streams independent of DLND_MODE value
--------------------------------------

There are no datamode independent streams for DLND.

.. _dlnd-fields:

-----------
Field names
-----------

DLND defines a set of pre-defined internal field names as well as mappings for how those field names map to the fields sent to the coupler.
In general, the stream input file should translate the stream input variable names into the ``dlnd_fld`` names below for use within the data land model.

.. csv-table:: "DLND internal field names"
   :header: "dlnd_fld (avifld)", "driver_fld (avofld)"
   :widths: 30, 30

   "t", "Sl_t"
   "tref", "Sl_tref"
   "qref", "Sl_qref"
   "avsdr", "Sl_avsdr"
   "anidr", "Sl_anidr"
   "avsdf", "Sl_avsdf"
   "anidf", "Sl_anidf"
   "snowh", "Sl_snowh"
   "taux", "Fall_taux"
   "tauy", "Fall_tauy"
   "lat", "Fall_lat"
   "sen", "Fall_sen"
   "lwup", "Fall_lwup"
   "evap", "Fall_evap"
   "swnet", "Fall_swnet"
   "lfrac", "Sl_landfrac"
   "fv", "Sl_fv"
   "ram1", "Sl_ram1"
   "flddst1", "Fall_flxdst1"
   "flxdst2", "Fall_flxdst2"
   "flxdst3", "Fall_flxdst3"
   "flxdst4", "Fall_flxdst4"
   "tsrfNN", "Sl_tsrfNN"
   "topoNN", "Sl_topoNN"
   "qiceNN",  "Flgl_qiceNN"

where NN = (00,01,02,..., ``glc_nec``), and ``glc_nec`` is the number of glacier elevation classes.
Note that the number of elevation classes on the input files must be the same as in the run.
