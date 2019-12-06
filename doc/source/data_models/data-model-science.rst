.. _data-model-science:

Data Model Science
==================

When a given data model is run, the user must specify which *science mode* it will run in.
Each data model has a fixed set of fields that it must send to the coupler, but it is the choice of mode that specifies how that set of fields is to be computed. 
Each mode activates various assumptions about what input fields are available from the input data streams, what input fields are available from the the coupler, and how to use this input data to compute the output fields sent to the coupler.

In general, a mode might specify...

- that fields be set to a time invariant constant (so that no input data is needed)
- that fields be taken directly from input data files (the input streams)
- that fields be computed using data read in from input files
- that fields be computed using data received from the coupler
- some combination of the above.

If a science mode is chosen that is not consistent with the input data provided, the model may abort (perhaps with a "missing data" error message), or the model may send erroneous data to the coupler (for example, if a mode assumes an input stream has temperature in Kelvin, but it really has temperature in Celsius).
Such an error is unlikely unless a user has edited the run scripts to specify either non-standard input data or a non-standard science mode. 
When editing the run scripts to use non-standard stream data or modes, users must be careful that the input data is consistent with the science mode and should verify that the data model is providing data to the coupler as expected.

The data model mode is a character string that is set in the namelist variable ``datamode`` in the namelist group ``shr_strdata_nml``. Although each data model, 
``datm``, ``dlnd``, ``drof``, ``docn``, ``dice`` and ``dwav`` has its own set of valid datamode values, two modes are common to all data models: ``COPYALL`` and ``NULL``.

``dataMode = "COPYALL"``
  The default mode is ``COPYALL`` -- the model will assume *all* the data that must be sent to the coupler will be found in the input data streams, and that this data can be sent to the coupler, unaltered, except for spatial and temporal interpolation.

``dataMode = "NULL"``
  ``NULL`` mode turns off the data model as a provider of data to the coupler. The ``model_present`` flag (eg. ``atm_present``) will be set to false and the coupler will assume no exchange of data to or from the data model.
