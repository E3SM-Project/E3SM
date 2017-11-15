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
For example, you can perform a fully-coupled run in which you ask for particular extra output streams; you can then use these saved "coupler history" files as inputs to datm to run a later land-only spinup.

In some cases, data models have prognostic functionality, that is, they also receive and use data sent by the driver.
However, in most cases, the data models are not running prognostically and have no need to receive any data from the driver.

The CIME data models have parallel capability and share significant amounts of source code.
Methods for reading and interpolating data have been established and can easily be reused:
The data model calls strdata ("stream data") methods which then call stream methods.
The stream methods are responsible for managing lists of input data files and their time axes.
The information is then passed up to the strdata methods where the data is read and interpolated in space and time.
The interpolated data is passed up to the data model where final fields are derived, packed, and returned to the driver.

------
Design
------
Data models function by reading in different streams of input data and interpolating those data both spatially and temporally to the appropriate final model grid and model time.
The strdata implementation does the following:

1. determines nearest lower and upper bound data from the input dataset 
2. if that is new data then read lower and upper bound data
3. fill lower and upper bound data
4. spatially map lower and upper bound data to model grid
5. time interpolate lower and upper bound data to model time
6. return fields to data model

The two timestamps of input data that bracket the present model time are read first.
These are called the lower and upper bounds of data and will change as the model advances. 
Those two sets of inputdata are first filled based on the user setting of the namelist variables ``str_fillalgo`` and ``str_fillmask``. 
That operation occurs on the input data grid.
The lower and upper bound data are then spatially mapped to the model grid based upon the user setting of the namelist variables ``str_mapalgo`` and ``str_mapmask``. 
Spatial interpolation only occurs if the input data grid and model grid are not identical, and this is determined in the strdata module automatically.
Time interpolation is the final step and is done using a time interpolation method specified by the user in namelist (via the ``shr_strdata_nml`` namelist variable ``tintalgo``). 
A final set of fields is then available to the data model on the model grid and for the current model time.
(See the :ref:`stream data namelist section <shr-strdata-nml>` for details on these and other namelist variables.)

**Each data model**

- communicates with the driver with fields on only the data model grid

- can be associated with multiple :ref:`streams<input-streams>`

- has an xml variable in ``env_run.xml`` that specifies its mode.
  These are: ``DATM_MODE``, ``DICE_MODE``, ``DLND_MODE``, ``DOCN_MODE``, ``DROF_MODE``, ``DWAV_MODE``.
  Each data model mode specifies the streams that are associated with that data model.

- has two :ref:`namelist<input-namelists>` groups in its input namelist file: a **stream-dependent** and a **stream-independent** namelist group.

- is associated with only one stream-independent namelist variable ``datamode`` (specified in the ``shr_strdata_nml`` namelist group) that determines if additional operations need to be performed on the input streams before returning to the driver.


**Each** ``DXXX_MODE`` **xml variable variable specfies 2 things:**

- the list of streams that are associated with the data model.

- a ``datamode`` namelist variable that is associated with each data model and that determines if additional operations need to be performed on the input streams before returning to the driver.

  At a minimum, all data models support ``datamode`` values of ``NULL`` and ``COPYALL``.

  - ``NULL`` - turns off the data model as a provider of data to the coupler.

  - ``COPYALL`` - copies all fields directly from the input data streams. Any required fields not found on an input stream will be set to zero.


**Each data model stream**

- can be associated with multiple stream input files (specified in the ``shr_strdata_nml`` namelist group).


**Each stream input file**

- can contain data on a unique grid and unique temporal time stamps.

- is interpolated to a single model grid and the present model time.

More details of the data model design are covered in :ref:`design details<design-details>`.

-------------
Next Sections
-------------
In the next sections, more details will be presented, including a full description of the science modes and namelist settings for the data atmosphere, data land, data runoff, data ocean, and data ice models; namelist settings for the strdata namelist input; a description of the format and options for the stream description input files; and a list of internal field names for each of the data components.
The internal data model field names are important because they are used to setup the stream description files and to map the input data fields to the internal data model field names.
