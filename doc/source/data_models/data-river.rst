.. _data-river:

=================
Data River (DROF)
=================

The data river model (DROF) provides river runoff data primarily to be used by the prognostic ocean component.
This data can either be observational (climatological or interannual river data) or data from a previous model run that is output to coupler history files and then read back in by DROF.

.. _drof-xml-vars:

-------------
xml variables
-------------

The following are xml variables that CIME supports for DROF.
These variables are defined in ``$CIMEROOT/src/components/data_comps/drof/cime_config/config_component.xml``.
These variables will appear in ``env_run.xml`` and are used by the DROF ``cime_config/buildnml`` script to generate the DROF namelist file ``drof_in`` and the required associated stream files for the case.

.. note:: These xml variables are used by the the drof's **cime_config/buildnml** script in conjunction with drof's **cime_config/namelist_definition_drof.xml** file to generate the namelist file ``drof_in``.

.. csv-table:: "DROF xml variables"
   :header: "xml variable", "description"
   :widths: 15, 85

   "DROF_MODE",                "Data mode"
   "",                         "Valid values are: NULL,CPLHIST,DIATREN_ANN_RX1,DIATREN_IAF_RX1,IAF_JRA"

   "DROF_CPLHIST_DOMAIN_FILE", "Coupler history forcing data mode - full pathname of model domain file "
   "DROF_CPLHIST_CASE",        "Coupler history forcing data mode - case name"
   "DROF_CPLHIST_DIR",         "Coupler history forcing data mode - directory containing coupler history forcing data"
   "DROF_CPLHIST_YR_ALIGN",    "Coupler history forcing data mode - simulation year corresponding to DROF_CPLHIST_YR_START"
   "DROF_CPLHIST_YR_START",    "Coupler history forcing data mode - starting year to loop forcing data over"
   "DROF_CPLHIST_YR_END",      "Coupler history forcing data mode - ending year to loop forcing data over"

.. note:: If ``DROF_MODE`` is set to ``CPLHIST``, it is normally assumed that the model domain will be identical to **all** of the stream domains. To ensure this, the xml variables ``ROF_DOMAIN_PATH`` and ``ROF_DOMAIN_FILE`` are ignored and a valid setting **must be given** for ``DROF_CPLHIST_DOMAIN_FILE``. If ``DROF_CPLHIST_DOMAIN_FILE`` is set to ``null``, then the drof component domain information is read in from the first coupler history file in the target stream and  it is assumed that the first coupler stream file that is pointed to contains the domain  information for that stream. This is the default mode that should be used for this mode. Alternatively, ``DROF_CPLHIST_DOMAIN_FILE`` can be set to ``$ROF_DOMAIN_PATH/$ROF_DOMAIN_FILE`` in a non-default configuration.

.. _drof-datamodes:

--------------------
datamode values
--------------------

The xml variable ``DROF_MODE`` (described in :ref:`drof_mode`) sets the streams that are associated with DROF and also sets the namelist variable ``datamode``.
``datamode`` (which appears in ``shr_strdata_nml``) specifies what additional operations need to be done by DROF on the streams before returning to the driver.

Each data model has its own set of supported ``datamode`` values. The following are the supported DROF ``datamode`` values, as defined in the file ``namelist_definition_drof.xml``.

.. csv-table:: "Valid values for datamode namelist variable"
   :header: "datamode variable", "description"
   :widths: 20, 80

   "NULL", "Turns off the data model as a provider of data to the coupler.  The rof_present flag will be set to false and the coupler will assume no exchange of data to or from the data model."
   "COPYALL", "Copies all fields directly from the input data streams Any required fields not found on an input stream will be set to zero."

---------
Namelists
---------

The data river runoff model (DROF) provides data river input to prognostic components such as the ocean.

The namelist file for DROF is ``drof_in``.

As is the case for all data models, DROF namelists can be separated into two groups, stream-independent and stream-dependent.
The stream dependent group is :ref:`shr_strdata_nml<input-streams>`.
The stream-independent group is ``drof_nml`` and the DROF stream-independent namelist variables are:

.. _drof-stream-independent-namelists:

=====================  ======================================================
decomp                 decomposition strategy (1d, root)

                       1d => vector decomposition, root => run on master task
restfilm               master restart filename
restfils               stream restart filename
force_prognostic_true  TRUE => force prognostic behavior
=====================  ======================================================

To change the namelist settings in ``drof_in``, edit the file ``user_nl_drof`` in your case directory.

.. _drof_mode:

-------------------------------
DROF_MODE, datamode and streams
-------------------------------

The following table describes the valid values of ``DROF_MODE`` (defined in the ``config_component.xml`` file for DROF), and how they relate to the associated input streams and the ``datamode`` namelist variable.
CIME will generate a value of ``DROF_MODE`` based on the compset.

.. csv-table:: "Relationship between DROF_MODE, datamode and streams"
   :header: "DROF_MODE", "description-streams-datamode"
   :widths: 15, 85

   "NULL", "null mode"
   "", "streams: none"
   "", "datamode: NULL"
   "DIATREN_ANN_RX1", "Reads in annual forcing river data used for CORE2 forcing runs."
   "", "streams: rof.diatren_ann_rx1"
   "", "datamode: COPYALL"
   "DIATREN_IAF_RX1", "Reads in intra-annual forcing river data used for CORE2 forcing runs."
   "", "streams: rof.diatren_iaf_rx1"
   "", "datamode: COPYALL"
   "CPLHIST", "Reads in data from coupler history files generated by a previous run."
   "", "streams: rof.cplhist"
   "", "datamode: COPYALL"
   "IAF_JRA", "Reads in intra-annual forcing river data used for JRA-55 forcing runs."
   "", "streams: rof.iaf_jra"
   "", "datamode: COPYALL"

.. _drof-mode-independent-streams:

------------------------------------------
Streams independent of DROF_MODE value
------------------------------------------

There are no datamode independent streams for DROF.

.. _drof-fields:

----------------
DROF Field names
----------------

DROF defines a set of pre-defined internal field names as well as mappings for how those field names map to the fields sent to the coupler.
In general, the stream input file should translate the stream input variable names into the ``drof_fld`` names for use within the data rofosphere model.

.. csv-table:: "DROF internal field names"
   :header: "drof_fld (avifld)", "driver_fld (avofld)"
   :widths: 30, 30

   "roff", "Forr_rofl"
   "ioff", "Forr_rofi"
