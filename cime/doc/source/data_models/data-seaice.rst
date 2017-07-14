.. _data-seaice:

Data Ice (DICE)
================

DICE is a combination of a data model and a prognostic model. 
The data functionality reads in ice coverage. 
The prognostic functionality calculates the ice/atmosphere and ice/ocean fluxes. 
DICE receives the same atmospheric input from the coupler as the active CICE model (i.e., atmospheric  states, shortwave fluxes, and ocean ice melt flux) and acts very similarly to CICE running in prescribed mode. 
Currently, this component is only used to drive POP in "C" compsets.

.. _dice-xml-vars:

---------------
xml variables
---------------
The following are xml variables that CIME supports for DICE. 
These variables are defined in ``$CIMEROOT/src/components/data_comps/dice/cime_config/config_component.xml``.
These variables will appear in ``env_run.xml`` and are used by the DICE ``cime_config/buildnml`` script to generate the DICE namelist file ``dice_in`` and the required associated stream files for the case.

.. note:: These xml variables are used by the the dice's **cime_config/buildnml** script in conjunction with dice's **cime_config/namelist_definition_dice.xml** file to generate the namelist file ``dice_in``.

.. csv-table:: "DICE xml variables"
   :header: "xml variable", "description"
   :widths: 15, 85

   "DICE_MODE", "Mode for sea-ice component"
   "","Valid values are: null, prescribed, ssmi, ssmi_iaf, ww3 "


.. _dice-datamodes:

--------------------
datamode values
--------------------

The xml variable ``DICE_MODE`` sets the streams that are associated with DICE and also sets the namelist variable ``datamode`` that specifies what additional operations need to be done by DICE on the streams before returning to the driver.
One of the variables in ``shr_strdata_nml`` is ``datamode``, whose value is a character string.  Each data model has a unique set of ``datamode`` values that it supports. 
The valid values for ``datamode`` are set in the file ``namelist_definition_dice.xml`` using the xml variable ``DICE_MODE`` in the ``config_component.xml`` file for DICE. 
CIME will generate a value ``datamode`` that is compset dependent. 

The following are the supported DICE datamode values and their relationship to the ``DICE_MODE`` xml variable value.

.. csv-table:: "Valid values for datamode namelist variable"
   :header: "datamode variable", "description"
   :widths: 20, 80

   "NULL", "Turns off the data model as a provider of data to the coupler.  The ice_present flag will be set to false and the coupler will assume no exchange of data to or from the data model."
   "COPYALL", "The default science mode of the data model is the COPYALL mode. This mode will examine the fields found in all input data streams, if any input field names match the field names used internally, they are copied into the export array and passed directly to the coupler without any special user code.  Any required fields not found on an input stream will be set to zero."
   "SSTDATA","Is a prognostic mode. It requires data be sent to the ice model. Ice fraction (extent) data is read from an input stream, atmosphere state variables are received from the coupler, and then an atmosphere-ice surface flux is computed and sent to the coupler. It is called ``SSTDATA`` mode because normally the ice fraction data is found in the same data files that provide SST data to the data ocean model. They are normally found in the same file because the SST and ice fraction data are derived from the same observational data sets and are consistent with each other. "

-------------------------------
DICE_MODE, datamode and streams
-------------------------------

The following tabe describes the valid values of ``DICE_MODE``, and how it relates to the associated input streams and the ``datamode`` namelist variable.

.. csv-table:: "Relationship between DICE_MODE, datamode and streams"
   :header: "DICE_MODE, "description-streams-datamode"
   :widths: 20, 80

   "null", "null mode"
   "", "streams: none"
   "", "datamode: null"
   "prescribed","prognostic mode - requires data to be sent to DICE"
   "","streams:  prescribed"
   "","datamode: SSTDATA"
   "ssmi", "Special Sensor Microwave Imager climatological data"
   "","streams: SSMI"
   "","datamode: SSTDATA"
   "ssmi", "Special Sensor Microwave Imager inter-annual forcing data"
   "","streams: SSMI_IAF"
   "","datamode: SSTDATA"
   "ww3", "ww3 mode"
   "", "streams: ww3"
   "", "datamode: COPYALL"

NIf DICE_MODE is set to ``ssmi``, ``ssmi_iaf`` or ``prescribed``, it is a prognostic mode and requires data be sent to the ice model.
Ice fraction (extent) data is read from an input stream, atmosphere state variables are received from the coupler, and then an atmosphere-ice surface flux is computed and sent to the coupler. 
Normally the ice fraction data is found in the same data files that provide SST data to the data ocean model. 
They are normally found in the same file because the SST and ice fraction data are derived from the same observational data sets and are consistent with each other.

.. _dice-namelists:

---------
Namelists
---------

The namelist file for DICE is ``dice_in`` (or ``dice_in_NNN`` for multiple instances).

As is the case for all data models, DICE namelists can be separated into two groups, stream-independent and stream-dependent. 

The stream dependent group is :ref:`shr_strdata_nml<input-streams>`. 

.. _dice-stream-independent-namelists:

The stream-independent group is ``dice_nml`` and the DICE stream-independent namelist variables are:

=====================  ======================================================
decomp                 decomposition strategy (1d, root)
    
                       1d => vector decomposition, root => run on master task
flux_qacc              activates water accumulation/melt wrt Q
flux_qacc0             initial water accumulation value
flux_qmin              bound on melt rate
flux_swpf              short-wave penetration factor
restfilm               master restart filename 
restfils               stream restart filename 
force_prognostic_true  TRUE => force prognostic behavior
=====================  ======================================================

To change the namelist settings in ``dice_in``, edit the file ``user_nl_dice``. 

.. _dice-mode-independent-streams:

--------------------------------------
Streams independent of DICE_MODE value
--------------------------------------

There are no datamode independent streams for DICE.

.. _dice-fields:

-----------
Field names
-----------

DICE defines a set of pre-defined internal field names as well as mappings for how those field names map to the fields sent to the coupler.
In general, the stream input file should translate the stream input variable names into the ``dice_fld`` names below for use within the data ice model.

.. csv-table:: "DICE internal field names"
   :header: "dice_fld (avifld)", "driver_fld (avofld)"
   :widths: 30, 30

   "to",    "So_t"	       
   "s",	    "So_s"	       
   "uo",    "So_u"	       
   "vo",    "So_v"	       
   "dhdx",  "So_dhdx"      
   "dhdy",  "So_dhdy"      
   "q",	    "Fioo_q"       
   "z",	    "Sa_z"	       
   "ua",    "Sa_u"	       
   "va",    "Sa_v"	       
   "ptem",  "Sa_ptem"      
   "tbot",  "Sa_tbot"      
   "shum",  "Sa_shum"      
   "dens",  "Sa_dens"      
   "swndr", "Faxa_swndr"   
   "swvdr", "Faxa_swvdr"   
   "swndf", "Faxa_swndf"   
   "swvdf", "Faxa_swvdf"   
   "lwdn",  "Faxa_lwdn"    
   "rain",  "Faxa_rain"    
   "snow",  "Faxa_snow"    
   "t",	    "Si_t"	       
   "tref",  "Si_tref"      
   "qref",  "Si_qref"      
   "ifrac", "Si_ifrac"     
   "avsdr", "Si_avsdr"     
   "anidr", "Si_anidr"     
   "avsdf", "Si_avsdf"     
   "anidf", "Si_anidf"     
   "tauxa", "Faii_taux"    
   "tauya", "Faii_tauy"    
   "lat",   "Faii_lat"     
   "sen",   "Faii_sen"     
   "lwup",  "Faii_lwup"    
   "evap",  "Faii_evap"    
   "swnet", "Faii_swnet"   
   "swpen", "Fioi_swpen"   
   "melth", "Fioi_melth"   
   "meltw", "Fioi_meltw"   
   "salt",  "Fioi_salt"    
   "tauxo", "Fioi_taux"    
   "tauyo", "Fioi_tauy"    











































