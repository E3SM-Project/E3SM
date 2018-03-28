.. _grids:

========================
Model grids
========================

CIME looks at the xml node ``GRIDS_SPEC_FILE`` in the **config_files.xml** file to identify supported out-of-the-box model grids for the target model. The node has the following contents:
::

   <entry id="GRIDS_SPEC_FILE">
     <type>char</type>
     <default_value>$CIMEROOT/cime_config/$MODEL/config_grids.xml</default_value>
     <group>case_last</group>
     <file>env_case.xml</file>
     <desc>file containing specification of all supported model grids, domains and mapping files (for documentation only - DO NOT EDIT)</desc>
     <schema>$CIMEROOT/cime_config/xml_schemas/config_grids_v2.xsd</schema>
   </entry>

Grid longname
-------------

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
