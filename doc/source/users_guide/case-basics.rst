.. _case-basics:

*********************************
The basics of CIME cases 
*********************************

Two concepts to understand before working with CIME are component sets and model grids.

- *Component sets*, which are usually referred to as "compsets," define both individual model components and any component-specific namelist or configuration settings that are used in a case.

- *Model grids* specify the grid or resolution for each component of the model.

Creating a CIME experiment or *case* requires, at a minimum, specifying a compset and a model grid.

Out-of-the-box compsets and model grids each have two names: a *longname* and an *alias* name. Examples of both follow.

Aliases are used for convenience. *Compset aliases* are unique; each is associated with one and only one compset. *Grid aliases*, on the other hand, are overloaded; the same grid alias may result in a different grid depending on the associated compset. Always confirm that the *compset longname* and the *grid longname* are correct when using aliases to create a case.

================
 Component sets
================

A compset longname has this form::

  TIME_ATM[%phys]_LND[%phys]_ICE[%phys]_OCN[%phys]_ROF[%phys]_GLC[%phys]_WAV[%phys]_ESP[_BGC%phys]

Supported values for each element of the longname::

  TIME = model time period (e.g. 1850, 2000, 20TR, RCP8...)

  CIME supports the following values for ATM,LND,ICE,OCN,ROF,GLC,WAV and ESP.
  ATM  = [DATM, SATM, XATM]
  LND  = [DLND, SLND, XLND]
  ICE  = [DICE, SICE, SICE]
  OCN  = [DOCN, SOCN, XOCN]
  ROF  = [DROF, SROF, XROF]
  GLC  = [SGLC, XGLC]
  WAV  = [SWAV, XWAV]
  ESP  = [SESP]

A CIME-driven model may have other options available.  Use **query_config** to determine the available options.

The OPTIONAL %phys attributes specify sub-modes of the given system.
For example, DOCN%DOM is the DOCN data ocean (rather than slab-ocean) mode.
ALL the possible %phys choices for each component are listed by
calling **query_case** with the --compsets all argument.  ALL data models have
a %phys option that corresponds to the data model mode.

As an example, this actual CESM compset longname refers to running a pre-industrial control with active CESM components CAM, CLM, CICE, POP2, MOSART, CISM2 and WW3 in a BDRD BGC coupling scenario::

   1850_CAM60_CLM50%BGC_CICE_POP2%ECO_MOSART_CISM2%NOEVOLVE_WW3_BGC%BDRD

The alias for this compset is B1850.

Either a compset longname or a compset alias can be input to **create_newcase**. You can also create your own custom compset. See *How do I create my own compset?* in the FAQ.

