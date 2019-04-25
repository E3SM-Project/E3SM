.. _data-atm:

Data Atmosphere (DATM)
======================

DATM is normally used to provide observational forcing data (or forcing data produced by a previous run using active components) to drive prognostic components.
In the case of CESM, these would be: CLM (I compset), POP2 (C compset), and POP2/CICE (G compset).
As a result, DATM variable settings are specific to the compset that will be targeted.
As examples, CORE2_NYF (CORE2 normal year forcing) is the DATM mode used in C and G compsets.
CLM_QIAN, CLMCRUNCEP, CLMGSWP3 and CLM1PT are DATM modes using observational data for forcing CLM in I compsets.

.. _datm-xml-vars:

------------------
xml variables
------------------
The following are ``$CASEROOT`` xml variables that CIME supports for DATM.
These variables are defined in ``$CIMEROOT/src/components/data_comps/datm/cime_config/config_component.xml``.
These variables will appear in ``env_run.xml`` and the resulting values are compset dependent.

.. note:: These xml variables are used by the the datm's **cime_config/buildnml** script in conjunction with datm's **cime_config/namelist_definition_datm.xml** file to generate the namelist file ``datm_in``.

.. csv-table:: "DATM xml variables"
   :header: "xml variable", "description"
   :widths: 20, 80

   "DATM_MODE",                "Mode for atmospheric component"
   "",                         "Valid values are: CORE2_NYF,CORE2_IAF,CLM_QIAN,CLM_QIAN_WISO,CLM1PT,CLMCRUNCEP,"
   "",                         "CLMCRUNCEP_V5,CLMGSWP3,WW3,CPLHIST,CORE_IAF_JRA"

   "DATM_PRESAERO",            "Optional prescribed aerosol forcing"
   "DATM_TOPO",                "Optional Surface topography"
   "DATM_CO2_TSERIES",         "Optional CO2 time series type"

   "DATM_CPLHIST_DOMAIN_FILE", "Coupler history forcing data mode - full pathname of model domain file "
   "DATM_CPLHIST_CASE",        "Coupler history forcing data mode - case name"
   "DATM_CPLHIST_DIR",         "Coupler history forcing data mode - directory containing coupler history data"
   "DATM_CPLHIST_YR_ALIGN",    "Coupler history forcing data mode - simulation year corresponding to DATM_CPLHIST_YR_START"
   "DATM_CPLHIST_YR_START",    "Coupler history forcing data mode - starting year to loop data over"
   "DATM_CPLHIST_YR_END",      "Coupler history forcing data mode - ending year to loop data over"

   "DATM_CLMNCEP_YR_ALIGN",    "I compsets only - simulation year corresponding to data starting year"
   "DATM_CLMNCEP_YR_START",    "I compsets only - data model starting year to loop data over"
   "DATM_CLMNCEP_YR_END",      "I compsets only - data model ending year to loop data over"

.. note:: If ``DATM_MODE`` is set to ``CPLHIST``, it is normally assumed that the model domain will be identical to **all** of the stream domains. To ensure this, the xml variables ``ATM_DOMAIN_PATH`` and ``ATM_DOMAIN_FILE`` are ignored and a valid setting **must be given** for ``DATM_CPLHIST_DOMAIN_FILE``. If ``DATM_CPLHIST_DOMAIN_FILE`` is set to ``null``, then the datm component domain information is read in from the first coupler history file in the target stream and it is assumed that the first coupler stream file that is pointed to contains the domain information for that stream. This is the default that should be used for this mode. Alternatively, ``DATM_CPLHIST_DOMAIN_FILE`` can be set to ``$ATM_DOMAIN_PATH/$ATM_DOMAIN_FILE`` in a non-default configuration.

.. _datm-datamodes:

--------------------
datamode values
--------------------

The xml variable ``DATM_MODE`` (described in :ref:`datm_mode`) sets the streams that are associated with DATM and also sets the namelist variable ``datamode``.
``datamode`` (which appears in ``shr_strdata_nml``) specifies what additional operations need to be done by DATM on the streams before returning to the driver.

Each data model has its own set of supported ``datamode`` values. The following are the supported DATM ``datamode`` values, as defined in the file ``namelist_definition_datm.xml``.

.. csv-table:: "Valid values for datamode namelist variable"
   :header: "datamode variable", "description"
   :widths: 20, 80

   "NULL", "This mode turns off the data model as a provider of data to the coupler. The ``atm_present`` flag will be set to ``false`` and the coupler assumes no exchange of data to or from the data model."
   "COPYALL", "The default science mode of the data model is the COPYALL mode. This mode will examine the fields found in all input data streams; if any input field names match the field names used internally, they are copied into the export array and passed directly to the coupler without any special user code.  Any required fields not found on an input stream will be set to zero except for aerosol deposition fields which will be set to a special value. "
   "CLMNCEP", "In conjunction with NCEP climatological atmosphere data, provides the atmosphere forcing favored by the Land Model Working Group when coupling an active land model with observed atmospheric forcing. This  mode replicates code previously found in CLM (circa 2005), before the LMWG started using the CIME coupling infrastructure and data models to do active-land-only simulations."
   "CORE2_NYF", "Coordinated Ocean-ice Reference Experiments (CORE) Version 2 Normal Year Forcing."
   "CORE2_IAF", "In conjunction with CORE Version 2 atmospheric forcing data, provides the atmosphere forcing favored by the Ocean Model Working Group when coupling an active ocean model with observed atmospheric forcing. This mode and associated data sets implement the CORE-IAF Version 2 forcing data, as developed by Large and Yeager (2008) at NCAR.  Note that CORE2_NYF and CORE2_IAF work exactly the same way."
   "CORE_IAF_JRA", "In conjunction with JRA-55 Project, provides the atmosphere forcing when coupling an active ocean model with observed atmospheric forcing. This mode and associated data sets implement the JRA-55 v1.3 forcing data."

.. _datm_mode:

-------------------------------
DATM_MODE, datamode and streams
-------------------------------

The following table describes the valid values of ``DATM_MODE`` (defined in the ``config_component.xml`` file for DATM), and how they relate to the associated input streams and the ``datamode`` namelist variable.
CIME will generate a value of ``DATM_MODE`` based on the compset.

.. csv-table:: "Relationship between DATM_MODE, datamode and streams"
   :header: "DATM_MODE", "description-streams-datamode"
   :widths: 15, 85

   "NULL", "null mode"
   "", "streams: none"
   "", "datamode: NULL"
   "CORE2_NYF","CORE2 normal year forcing (C ang G compsets)"
   "", "streams: CORE2_NYF.GISS,CORE2_NYF.GXGXS,CORE2_NYF.NCEP"
   "", "datamode: CORE2_NYF"
   "CORE2_IAF","CORE2 interannual year forcing (C ang G compsets)"
   "", "streams: CORE2_IAF.GCGCS.PREC,CORE2_IAF.GISS.LWDN,CORE2_IAF.GISS.SWDN,CORE2_IAF.GISS.SWUP,"
   "", "CORE2_IAF.NCEP.DN10,CORE2_IAF.NCEP.Q_10,CORE2_IAF.NCEP.SLP_,CORE2_IAF.NCEP.T_10,CORE2_IAF.NCEP.U_10,"
   "", "CORE2_IAF.NCEP.V_10,CORE2_IAF.CORE2.ArcFactor"
   "", "datamode: CORE2_IAF"
   "CORE_IAF_JRA",JRA-55 intra-annual year forcing(C ang G compsets)"
   "", "streams: CORE_IAF_JRA.PREC,CORE_IAF_JRA.LWDN,CORE_IAF_JRA.SWDN,"
   "", "CORE_IAF_JRA.Q_10,CORE_IAF_JRA.SLP_,CORE_IAF_JRA.T_10,CORE_IAF_JRA.U_10,"
   "", "CORE_IAF_JRA.V_10,CORE_IAF_JRA.CORE2.ArcFactor"
   "", "datamode: CORE_IAF_JRA"
   "CLM_QIAN_WISO","QIAN atm input data with water isotopes (I compsets)"
   "", "streams: CLM_QIAN_WISO.Solar,CLM_QIAN_WISO.Precip,CLM_QIAN_WISO.TPQW"
   "", "datamode: CLMNCEP"
   "CLM_QIAN", "QIAN atm input data (I compsets)"
   "", "streams: CLM_QIAN.Solar,CLM_QIAN.Precip,CLM_QIAN.TPQW"
   "", "datamode: CLMNCEP"
   "CLMCRUNCEP","CRUNCEP atm input data (I compsets)"
   "", "streams: CLMCRUNCEP.Solar,CLMCRUNCEP.Precip,CLMCRUNCEP.TPQW"
   "", "datamode: CLMNCEP"
   "CLMCRUNCEP_V5","CRUNCEP atm input data (I compsets)"
   "","streams: CLMCRUNCEP_V5.Solar,CLMCRUNCEP_V5.Precip,CLMCRUNCEP_V5.TPQW"
   "","datamode: CLMNCEP"
   "CLMGSWP3","GSWP3 atm input data (I compsets)"
   "","streams: CLMGSWP3.Solar,CLMGSWP3.Precip,CLMGSWP3.TPQW"
   "","datamode: CLMNCEP"
   "CLM1PT", "single point tower site atm input data"
   "","streams: CLM1PT.$ATM_GRID"
   "","datamode: CLMNCEP"
   "CPLHIST","user generated forcing data from using coupler history files used to spinup relevant prognostic components (for CESM this is CLM, POP and CISM)"
   "","streams: CPLHISTForcing.Solar,CPLHISTForcing.nonSolarFlux,"
   "","CPLHISTForcing.State3hr,CPLHISTForcing.State1hr"
   "","datamode: CPLHIST"
   "WW3","WW3 wave watch data from a short period of hi WW3 wave watch data from a short period of hi temporal frequency COREv2 data"
   "","streams: WW3"
   "","datamode: COPYALL"

--------------
Namelists
--------------

The DATM namelist file is ``datm_in`` (or ``datm_in_NNN`` for multiple instances). DATM namelists can be separated into two groups: *stream-independent* namelist variables that are specific to the DATM model and *stream-specific* namelist variables whose names are common to all the data models.

Stream dependent input is in the namelist group ``"shr_strdata_nml`` which is discussed in :ref:`input streams <input-streams>` and is the same for all data models.

.. _datm-stream-independent-namelists:

The stream-independent group is ``datm_nml`` and the DATM stream-independent namelist variables are:

=====================  =============================================================================================
datm_nml vars          description
=====================  =============================================================================================
decomp                 decomposition strategy (1d, root)

                       1d => vector decomposition, root => run on master task
restfilm               master restart filename
restfils               stream restart filename
force_prognostic_true  TRUE => force prognostic behavior
bias_correct           if set, include bias correction streams in namelist
anomaly_forcing        if set, includ anomaly forcing streams in namelist
factorfn               filename containing correction factors for use in CORE2 modes (CORE2_IAF and CORE2_NYF)
presaero               if true, prescribed aerosols are sent from datm
iradsw                 frequency to update radiation in number of time steps (of hours if negative)
wiso_datm              if true, turn on water isotopes
=====================  =============================================================================================

.. _datm-mode-independent-streams:

------------------------------------------
Streams independent of DATM_MODE value
------------------------------------------

In general, each ``DATM_MODE`` xml variable is identified with a unique set of streams.
However, there are several streams in DATM that can accompany any ``DATM_MODE`` setting.
Currently, these are streams associated with prescribed aerosols, co2 time series, topography, anomoly forcing and bias correction.
These mode-independent streams are activated different, depending on the stream.

- ``prescribed aerosol stream:``
  To add this stream, set ``$DATM_PRESAERO`` to a supported value other than ``none``.

- ``co2 time series stream``:
  To add this stream, set ``$DATM_CO2_TSERIES`` to a supported value other than ``none``.

- ``topo stream``:
  To add this stream, set ``$DATM_TOPO`` to a supported value other than ``none``.

- ``anomaly forcing stream:``
  To add this stream, you need to add any of the following keywword/value pair to the end of ``user_nl_datm``:
  ::

    Anomaly.Forcing.Precip = <filename>
    Anomaly.Forcing.Temperature = <filename>
    Anomaly.Forcing.Pressure = <filename>
    Anomaly.Forcing.Humidity = <filename>
    Anomaly.Forcing.Uwind = <filename>
    Anomaly.Forcing.Vwind = <filename>
    Anomaly.Forcing.Shortwave = <filename>
    Anomaly.Forcing.Longwave = <filename>

- ``bias_correct stream:``
  To add this stream, you need to add any of the following keywword/value pair to the end of ``user_nl_datm``:
  ::

   BC.QIAN.CMAP.Precip = <filename>
   BC.QIAN.GPCP.Precip = <filename>
   BC.CRUNCEP.CMAP.Precip = <filename>
   BC.CRUNCEP.GPCP.Precip = <filename>

.. _datm-fields:

----------------
DATM Field names
----------------

DATM defines a set of pre-defined internal field names as well as mappings for how those field names map to the fields sent to the coupler.
In general, the stream input file should translate the stream input variable names into the ``datm_fld`` names for use within the data atmosphere model.

.. csv-table:: "DATM internal field names"
   :header: "datm_fld (avifld)", "driver_fld (avofld)"
   :widths: 30, 30

    "z",         "Sa_z"
    "topo",      "Sa_topo"
    "u",         "Sa_u"
    "v",         "Sa_v"
    "tbot",      "Sa_tbot"
    "ptem",      "Sa_ptem"
    "shum",      "Sa_shum"
    "dens",      "Sa_dens"
    "pbot",      "Sa_pbot"
    "pslv",      "Sa_pslv"
    "lwdn",      "Faxa_lwdn"
    "rainc",     "Faxa_rainc"
    "rainl",     "Faxa_rainl"
    "snowc",     "Faxa_snowc"
    "snowl",     "Faxa_snowl"
    "swndr",     "Faxa_swndr"
    "swvdr",     "Faxa_swvdr"
    "swndf",     "Faxa_swndf"
    "swvdf",     "Faxa_swvdf"
    "swnet",     "Faxa_swnet"
    "co2prog",   "Sa_co2prog"
    "co2diag",   "Sa_co2diag"
    "bcphidry",  "Faxa_bcphidry"
    "bcphodry",  "Faxa_bcphodry"
    "bcphiwet",  "Faxa_bcphiwet"
    "ocphidry",  "Faxa_ocphidry"
    "ocphodry",  "Faxa_ocphodry"
    "ocphiwet",  "Faxa_ocphiwet"
    "dstwet1",   "Faxa_dstwet1"
    "dstwet2",   "Faxa_dstwet2"
    "dstwet3",   "Faxa_dstwet3"
    "dstwet4",   "Faxa_dstwet4"
    "dstdry1",   "Faxa_dstdry1"
    "dstdry2",   "Faxa_dstdry2"
    "dstdry3",   "Faxa_dstdry3"
    "dstdry4",   "Faxa_dstdry4"
    "tref",      "Sx_tref"
    "qref",      "Sx_qref"
    "avsdr",     "Sx_avsdr"
    "anidr",     "Sx_anidr"
    "avsdf",     "Sx_avsdf"
    "anidf",     "Sx_anidf"
    "ts",        "Sx_t"
    "to",        "So_t"
    "snowhl",    "Sl_snowh"
    "lfrac",     "Sf_lfrac"
    "ifrac",     "Sf_ifrac"
    "ofrac",     "Sf_ofrac"
    "taux",      "Faxx_taux"
    "tauy",      "Faxx_tauy"
    "lat",       "Faxx_lat"
    "sen",       "Faxx_sen"
    "lwup",      "Faxx_lwup"
    "evap",      "Faxx_evap"
    "co2lnd",    "Fall_fco2_lnd"
    "co2ocn",    "Faoo_fco2_ocn"
    "dms",       "Faoo_fdms_ocn"
    "precsf",    "Sa_precsf"
    "prec_af",   "Sa_prec_af"
    "u_af",      "Sa_u_af"
    "v_af",      "Sa_v_af"
    "tbot_af",   "Sa_tbot_af"
    "pbot_af",   "Sa_pbot_af"
    "shum_af",   "Sa_shum_af"
    "swdn_af",   "Sa_swdn_af"
    "lwdn_af",   "Sa_lwdn_af"
    "rainc_18O", "Faxa_rainc_18O"
    "rainc_HDO", "Faxa_rainc_HDO"
    "rainl_18O", "Faxa_rainl_18O"
    "rainl_HDO", "Faxa_rainl_HDO"
    "snowc_18O", "Faxa_snowc_18O"
    "snowc_HDO", "Faxa_snowc_HDO"
    "snowl_18O", "Faxa_snowl_18O"
    "snowl_HDO", "Faxa_snowl_HDO"
    "shum_16O",  "Sa_shum_16O"
    "shum_18O",  "Sa_shum_18O"
