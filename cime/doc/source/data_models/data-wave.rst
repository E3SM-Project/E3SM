.. _data-wave:

=================
Data Wave (DWAV)
=================

The data wave model (DWAV) provides data wave forcing primarily to be used by the prognostic ocean component.
Currently, this data is climatological.

.. _dwav-xml-vars:

-------------
xml variables
-------------

The following are XML variables that CIME supports for DWAV.
These variables will appear in ``env_run.xml`` and are used by the DWAV ``cime_config/buildnml`` script to generate the DWAV namelist file ``dwav_in`` and the required associated stream files for the case.

.. note:: These XML variables are used by the the DWAV **cime_config/buildnml** script in conjunction with the DWAV **cime_config/namelist_definition_dwav.xml** file to generate the namelist file ``dwav_in``.

.. csv-table:: "DROF xml variables"
   :header: "xml variable", "description"
   :widths: 15, 85

   "DWAV_MODE", "Data mode"
   "", "Valid values are: NULL, CLIMO"

.. _dwav-datamodes:

--------------------
DWAV datamode values
--------------------

One of the variables in ``shr_strdata_nml`` is ``datamode``, whose value is a character string.
Each data model has a unique set of ``datamode`` values that it supports.

The valid values for ``datamode`` are set by the XML variable ``DWAV_MODE`` in the ``config_component.xml`` file for DWAV.
CIME will generate a value ``datamode`` that is compset dependent.

The following are the supported DWAV datamode values and their relationship to the ``DWAV_MODE`` xml variable value.

.. csv-table:: Relationship between ``DWAV_MODE`` xml variables and ``datamode`` namelist variables
   :header: "DWAV_MODE (xml)", "datamode (namelist)"
   :widths: 15, 90

   "NULL", "NULL"	    
   "", "This mode turns off the data model as a provider of data to the coupler. "
   "", "The ``wav_present`` flag will be set to ``false`` and the coupler assumes no exchange of data to or from the data model."
   "CLIMO", "COPYALL"
   "", "Examines the fields found in all input data streams and if any input field names match the field names used internally, "
   "", "they  are copied into the export array and passed directly to the coupler  without any special user code."


.. _dwav-namelists:

---------
Namelists
---------

The data wave model (DWAV) provides data wave input to prognostic components such as the ocean.

The namelist file for DWAV is ``dwav_in``.

As is the case for all data models, DWAV namelists can be separated into two groups, stream-independent and stream-dependent.
The stream dependent group is :ref:`shr_strdata_nml<input-streams>`.
The stream-independent group is ``dwav_nml`` and the DWAV stream-independent namelist variables are:

.. _dwav-stream-independent-namelists:

=====================  ======================================================
decomp                 decomposition strategy (1d, root)

                       1d => vector decomposition, root => run on master task
restfilm               master restart filename
restfils               stream restart filename
force_prognostic_true  TRUE => force prognostic behavior
=====================  ======================================================

To change the namelist settings in ``dwav_in``, edit the file ``user_nl_dwav`` in your case directory.

.. _dwav_mode:

-------------------------------
DWAV_MODE, datamode and streams
-------------------------------

The following table describes the valid values of ``DWAV_MODE`` (defined in the ``config_component.xml`` file for DWAV), and how they relate to the associated input streams and the ``datamode`` namelist variable.
CIME will generate a value of ``DWAV_MODE`` based on the compset.

.. csv-table:: "Relationship between DWAV_MODE, datamode and streams"
   :header: "DWAV_MODE", "description-streams-datamode"
   :widths: 15, 85

   "NULL", "null mode"
   "", "streams: none"
   "", "datamode: NULL"


.. _dwav-mode-independent-streams:

--------------------------------------
Streams independent of DWAV_MODE value
--------------------------------------

There are no datamode independent streams for DWAV.

.. _dwav-fields:

----------------
Field names
----------------

DWAV defines a set of pre-defined internal field names as well as mappings for how those field names map to the fields sent to the coupler.
In general, the stream input file should translate the stream input variable names into the ``dwav_fld`` names below for use within the data wave model.

.. csv-table:: "DWAV internal field names"
   :header: "dwav_fld (avifld)", "driver_fld (avofld)"
   :widths: 30, 30

   "lamult", "Sw_lamult" 
   "ustokes","Sw_ustokes"
   "vstokes", "Sw_vstokes"





