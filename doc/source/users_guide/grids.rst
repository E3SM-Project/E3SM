.. _grids:

========================
Model grids
========================

CIME looks at the xml node ``GRIDS_SPEC_FILE`` in  **$CIMEROOT/config/$models/config_files.xml** file to identify supported out-of-the-box model grids for the target model.

The node has the following contents:
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

CIME model grids generally are associated with a specific combination of atmosphere, land, land-ice, river-runoff and ocean/ice grids. The naming convention for these grids uses only atmosphere, land, and ocean/ice grid specifications.

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

.. _adding-cases:

Adding grids
-------------

.. _adding-a-grid:

CIME supports numerous out-of-the box model resolutions. To see the grids that are supported, call `query_config <../Tools_user/query_config.html>`_ as shown below.
   ::

      > query_config --grids

The most common resolutions have the atmosphere and land components on one grid and the ocean and ice on a second grid. The following overview assumes that this is the case.
The naming convention looks like *f19_g17*, where the f19 indicates that the atmosphere and land are on the 1.9x2.5 (finite volume dycore) grid while the g17 means the ocean and ice are on the gx1v6 one-degree displaced pole grid.

CIME enables users to add their own component grid combinations.
The steps for adding a new component grid to the model system follow. This process can be simplified if the atmosphere and land are running on the same grid.

1. The first step is to generate SCRIP grid files for the atmosphere, land, ocean, land-ice, river and wave component grids that will comprise your model grid.
   If you are introducing just one new grid, you can leverage SCRIP grid files that are already in place for the other components.
   There is no supported functionality for creating the SCRIP format file.

2. Build the **check_map** utility by following the instructions in **$CCSMROOT/mapping/check_maps/INSTALL**. Also confirm that the `ESMF <http://www.cesm.ucar.edu/models2.0/external-link-here>`_ toolkit is installed on your machine.

   When you add new user-defined grid files, you also need to generate a set of mapping files so the coupler can send data from a component on one grid to a component on another grid.
   There is an ESMF tool that tests the mapping file by comparing a mapping of a smooth function to its true value on the destination grid.
   We have tweaked this utility to test a suite of smooth functions, as well as ensure conservation (when the map is conservative).
   Before generating mapping functions it is *highly recommended* that you build this utility.

3. Generate these mapping files:
   ::

     atm <-> ocn
     atm <-> wav
     lnd <-> rof
     lnd <-> glc
     ocn <-> wav
     rof -> ocn

  Using the SCRIP grid files from Step 1, generate a set of conservative (area-averaged) and non-conservative (patch and bilinear) mapping files.

  You can do this by calling **gen_cesm_maps.sh** in ``$CCSMROOT/tools/mapping/gen_mapping_files/``.
  This script generates all the mapping files needed except ``rof -> ocn``, which is discussed below.
  This script uses the ESMF offline weight generation utility, which you must build *prior* to running **gen_cesm_maps.sh**.

  The **README** file in the **gen_mapping_files/** directory describes how to run **gen_cesm_maps.sh**. The basic usage is shown here:
   ::

    > cd $CCSMROOT/mapping/gen_mapping_files
    > ./gen_cesm_maps.sh \
       --fileocn  <input SCRIP ocn_grid full pathname>  \
       --fileatm  <input SCRIP atm grid full pathname>  \
       --filelnd  <input SCRIP lnd grid full pathname>  \
       --filertm  <input SCRIP rtm grid full pathname>  \
       --nameocn  <ocnname in output mapping file> \
       --nameatm  <atmname in output mapping file> \
       --namelnd  <lndname in output mapping file> \
       --namertm  <rtmname in output mapping file>

  This command generates the following mapping files:
   ::

     map_atmname_TO_ocnname_aave.yymmdd.nc
     map_atmname_TO_ocnname_blin.yymmdd.nc
     map_atmname_TO_ocnname_patc.yymmdd.nc
     map_ocnname_TO_atmname_aave.yymmdd.nc
     map_ocnname_TO_atmname_blin.yymmdd.nc
     map_atmname_TO_lndname_aave.yymmdd.nc
     map_atmname_TO_lndname_blin.yymmdd.nc
     map_lndname_TO_atmname_aave.yymmdd.nc
     map_ocnname_TO_lndname_aave.yymmdd.nc
     map_lndname_TO_rtmname_aave.yymmdd.nc
     map_rtmname_TO_lndname_aave.yymmdd.nc

   .. note:: You do not need to specify all four grids. For example, if you are running with the atmosphere and land on the same grid, then you do not need to specify the land grid (and atm<->rtm maps will be generated).
                   If you also omit the runoff grid, then only the 5 atm<->ocn maps will be generated.

   .. note:: ESMF_RegridWeightGen runs in parallel, and the ``gen_cesm_maps.sh`` script has been written to run on yellowstone.
                   To run on any other machine, you may need to add some environment variables to ``$CCSMROOT/mapping/gen_mapping_files/gen_ESMF_mapping_file/create_ESMF_map.sh`` -- search for hostname to see where to edit the file.

4. Generate atmosphere, land and ocean / ice domain files.

   Using the conservative ocean to land and ocean to atmosphere mapping files created in the previous step, you can create domain files for the atmosphere, land, and ocean; these are basically grid files with consistent masks and fractions.
   You make these files by calling **gen_domain** in **$CCSMROOT/mapping/gen_domain_files**.
   The **INSTALL** file in the **gen_domain_files/** directory describes how to build the **gen_domain** executable. The **README** file in the same directory explains how to use the tool. The basic usage is:
   ::

      > ./gen_domain -m ../gen_mapping_files/map_ocnname_TO_lndname_aave.yymmdd.nc -o ocnname -l lndname
      > ./gen_domain -m ../gen_mapping_files/map_ocnname_TO_atmname_aave.yymmdd.nc -o ocnname -l atmname

   These commands generate the following domain files:
   ::

      domain.lnd.lndname_ocnname.yymmdd.nc
      domain.ocn.lndname_ocnname.yymmdd.nc
      domain.lnd.atmname_ocnname.yymmdd.nc
      domain.ocn.atmname_ocnname.yymmdd.nc
      domain.ocn.ocnname.yymmdd.nc

   .. note:: The input atmosphere grid is assumed to be unmasked (global). Land cells whose fraction is zero will have land mask = 0.

   .. note:: If the ocean and land grids *are identical* then the mapping file will simply be unity and the land fraction will be one minus the ocean fraction.

5. If you are adding a new ocn or rtm grid, create a new rtm->ocn mapping file. (Otherwise you can skip this step.)
   The process for mapping from the runoff grid to the ocean grid is currently undergoing many changes.
   At this time, if you are running with a new ocean or runoff grid, please contact Michael Levy (mlevy_AT_ucar_DOT_edu) for assistance. If you are running with standard ocean and runoff grids, the mapping file should already exist and you do not need to generate it.


6. CESM specific: If you are adding a new atmosphere grid, this means you are also generating a new land grid, and you will need to create a new CLM surface dataset. (Otherwise you can skip this step).
   You need to first generate mapping files for CLM surface dataset (since this is a non-standard grid).
   ::

      > cd $CCSMROOT/models/lnd/clm/tools/mkmapdata
      > ./mkmapdata.sh --gridfile <lnd SCRIP grid file> --res <atm resolution name> --gridtype global

    These mapping files are then used to generate CLM surface dataset. Below is an example for a current day surface dataset (model year 2000).

    ::

       > cd  $CCSMROOT/models/lnd/clm/tools/mksurfdata_map
       > ./mksurfdata.pl -res usrspec -usr_gname <atm resolution name> -usr_gdate yymmdd -y 2000

7. Create grid file needed for create_newcase.
   The next step is to add the necessary new entries in the appropriate ``config_grids.xml`` file.
   You will need to modify ``$CIMEROOT/config/cesm/config_grids.xml`` or ``$CIMEROOT/config/e3sm/config_grids.xml`` depending on the value of ``$CIME_MODEL``.
   You will need to:

   - add a single  ``<model_grid>`` entry
   - add possibly multiple ``<domain>`` entries for  every new component grid that you have added
   - add possibly multiple ``<gridmap>`` entries for all the new component combinations that require new mapping files

8. Test new grid.

   Below assume that the new grid is an atmosphere grid.
   ::

      Test the new grid with all data components.
      (write an example)
      Test the new grid with CAM(newgrid), CLM(newgrid), DOCN(gx1v6), DICE(gx1v6)
      (write an example)
