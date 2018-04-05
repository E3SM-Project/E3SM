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

===============================
 Model Grids
===============================

A model grid longname has the form::

  a%name_l%name_oi%name_r%name_m%mask_g%name_w%name

For reference::

  a%  = atmosphere grid
  l%  = land grid
  oi% = ocean/sea-ice grid (must be the same)
  r%  = river grid
  m%  = ocean mask grid
  g%  = internal land-ice grid
  w%  = wave component grid

The ocean mask grid determines land/ocean boundaries in the model.
On the ocean grid, a grid cell is assumed to be either all ocean or all land.
The land mask on the land grid is obtained by mapping the ocean mask
(using first-order conservative mapping) from the ocean grid to the land grid.

From the point of view of model coupling, the glc grid is assumed to
be identical to the land grid. The internal land-ice grid can be different,
however, and is specified by the g% value.

As an example, examine this actual grid longname::

   a%ne30np4_l%ne30np4_oi%gx1v7_r%r05_m%gx1v7_g%null_w%null

It refers to a model grid with a ne30np4 spectral element (approximately 1-degree) atmosphere and land grids, gx1v7 Greenland pole, 1-degree ocean and sea-ice grids, a 1/2 degree river routing grid, null wave and internal cism grids, and an gx1v7 ocean mask.
The alias for this grid is ne30_g16.

CIME also permits users to introduce their own :ref:`user-defined grids <adding-a-grid>`.

Component grids are denoted by the following naming convention:

- "[dlat]x[dlon]" are regular lon/lat finite volume grids where dlat and dlon are the approximate grid spacing. The shorthand convention is "fnn" where nn generally is a pair of numbers indicating the resolution. An example is 1.9x2.5 or f19 for the approximately "2-degree" finite-volume grid. Note that CAM uses an [nlat]x[nlon] naming convention internally for this grid.

- "Tnn" are spectral lon/lat grids where nn is the spectral truncation value for the resolution. The shorthand name is identical. Example: T85.

- "ne[X]np[Y]" are cubed sphere resolutions where X and Y are integers. The short name generally is ne[X]. Examples: ne30np4 or ne30.

- "pt1" is a single grid point.

- "gx[D]v[n]" is a POP displaced pole grid where D is the approximate resolution in degrees and n is the grid version. The short name generally is g[D][n]. An example is gx1v7 or g17 for a grid of approximately 1-degree resolution.

- "tx[D]v[n]" is a POP tripole grid where D is the approximate resolution in degrees and n is the grid version.

- "oRSS[x]to[y]" is an MPAS grid with grid spacing from x to y kilometers.

- "oEC[x]to[y]" is an MPAS grid with grid spacing from x to y kilometers.
