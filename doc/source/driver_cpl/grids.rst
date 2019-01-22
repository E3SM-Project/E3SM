=====================
Grids
=====================

----------------------------
Standard Grid Configurations
----------------------------
The standard implementation for grids in CIME has been that the atmosphere and land models are run on identical grids and the ocean and sea ice model are run on identical grids.
The ocean model mask is used to derive a complementary mask for the land grid such that for any given combination of atmosphere/land and ocean/ice grids, there is a unique land mask. 
This approach for dealing with grids is still used a majority of the time.
But there is a new capability, called ``trigrid`` that allows the atmosphere and land grids to be unique. 
A typical grid is the finite volume "1 degree" atmosphere/land grid matched with the  "1 degree" ocean/ice grid. 
The runoff grid is generally unique to runoff and the land ice grid is coupled on the land grid with interpolation carried out to a unique land ice grid inside that component.

Historically, the ocean grid has been the higher resolution grid in CIME model configurations.
While that is no longer always the case, the current driver implementation largely reflects that presumption. 
The atmosphere/ocean fluxes in the coupler are computed on the ocean grid.
The driver namelist variable ``aoflux_grid`` allows users to specify the atmosphere/ocean flux computation grid in the coupler in the future. 
In addition, the default mapping approach used also reflects the presumption that the ocean is generally higher resolution.
Fluxes are always mapped using a locally conservative area average methods to preserve conservation. 
However, states are mapped using bilinear interpolation from the atmosphere grid to the ocean grid to better preserve gradients, while they are mapped using a locally conservative area average approach from the ocean grid to the atmosphere grid.
These choices are based on the presumption that the ocean grid is higher resolution.

There has always been an option that all grids (atmosphere, land, ocean, and ice) could be identical, and this is still supported. 
There are a couple of namelist variables, ``samegrid_ao``, ``samegrid_al``, and ``samegrid_ro`` that tell the coupler whether to expect that the following grids; atmosphere/ocean, atmosphere/land, and runoff/ocean respectively are identical. 
These are set automaticaly in the driver namelist depending on the grid chosen and impact mapping as well as domain checking.

----------------------
Trigrid Configurations
----------------------
Grid configurations are allowed where the atmosphere and land grids are unique. 

The trigrid implementation introduces an ambiguity in the definition of the mask.
This ambiguity is associated with an inability to define an absolutely consistent ocean/land mask across all grids in the system. 
A summary of trigrid support follows:
- The land mask is defined on the atmosphere grid as the complement of the ocean mask mapped conservatively to the atmosphere grid. 
- Then the land and ocean masks are exactly complementary on the atmosphere grid where conservative merging are critical.
- No precise land fraction needs to be defined in the land grid. 
- The only requirement is that the land model compute data on a masked grid such that when mapped to the atmosphere grid, all atmosphere grid points that contain some fraction of land have valid values computed in the land model.
- There are an infinite number of land fraction masks that can accomplish this including a fraction field that is exactly one at every grid cell. 
- In the land model, all land fraction masks produce internally conservative results.
- The only place where the land fraction becomes important is mapping the land model output to the runoff model. 
- In that case, the land fraction on the land grid is applied to the land to runoff mapping.

---------
Fractions
---------
The component grid fractions in the coupler are defined and computed in ``$CIMEROOT/driver_cpl/driver/seq_frac_mct``.
A slightly modified version of the notes from this file is pasted below.
Just to clarify some of the terms. 
- fractions_a, fractions_l, fractions_i, and fractions_o are the fractions on the atmosphere, land, ice, and ocean grids. 
- afrac, lfrac, ifrac, and ofrac are the atmosphere, land, ice, and ocean fractions on those grids. 
So fractions_a(lfrac) is the land fraction on the atmosphere grid.
lfrin in the land fraction defined in the land model. 
This can be different from lfrac because of the trigrid implementation.
lfrac is the land fraction consistent with the ocean mask and lfrin is the land fraction in the land model. 
ifrad and ofrad are fractions at the last radiation timestep.
These fractions preserve conservation of heat in the net shortwave calculation because the net shortwave calculation is one timestep behind the ice fraction evolution in the system. 
When the variable "dom" is mentioned below, that refers to a field sent from a component at initialization.
::

   !  the fractions fields are now afrac, ifrac, ofrac, lfrac, and lfrin.
   !    afrac = fraction of atm on a grid
   !    lfrac = fraction of lnd on a grid
   !    ifrac = fraction of ice on a grid
   !    ofrac = fraction of ocn on a grid
   !    lfrin = land fraction defined by the land model
   !    ifrad = fraction of ocn on a grid at last radiation time
   !    ofrad = fraction of ice on a grid at last radiation time
   !      afrac, lfrac, ifrac, and ofrac are the self-consistent values in the
   !      system.  lfrin is the fraction on the land grid and is allowed to
   !      vary from the self-consistent value as descibed below.  ifrad
   !      and ofrad are needed for the swnet calculation.
   !  the fractions fields are defined for each grid in the fraction bundles as
   !    needed as follows.
   !    character(*),parameter :: fraclist_a = 'afrac:ifrac:ofrac:lfrac:lfrin'
   !    character(*),parameter :: fraclist_o = 'afrac:ifrac:ofrac:ifrad:ofrad'
   !    character(*),parameter :: fraclist_i = 'afrac:ifrac:ofrac'
   !    character(*),parameter :: fraclist_l = 'afrac:lfrac:lfrin'
   !    character(*),parameter :: fraclist_g = 'gfrac'
   !
   !  we assume ocean and ice are on the same grids, same masks
   !  we assume ocn2atm and ice2atm are masked maps
   !  we assume lnd2atm is a global map
   !  we assume that the ice fraction evolves in time but that
   !    the land model fraction does not.  the ocean fraction then
   !    is just the complement of the ice fraction over the region
   !    of the ocean/ice mask.
   !  we assume that component domains are filled with the total
   !    potential mask/fraction on that grid, but that the fractions
   !    sent at run time are always the relative fraction covered.
   !    for example, if an atm cell can be up to 50% covered in
   !    ice and 50% land, then the ice domain should have a fraction
   !    value of 0.5 at that grid cell.  at run time though, the ice
   !    fraction will be between 0.0 and 1.0 meaning that grid cells
   !    is covered with between 0.0 and 0.5 by ice.  the "relative" fractions
   !    sent at run-time are corrected by the model to be total fractions
   !    such that
   !  in general, on every grid,
   !              fractions_*(afrac) = 1.0
   !              fractions_*(ifrac) + fractions_*(ofrac) + fractions_*(lfrac) = 1.0
   !  where fractions_* are a bundle of fractions on a particular grid and
   !    *frac (ie afrac) is the fraction of a particular component in the bundle.
   !
   !  the fractions are computed fundamentally as follows (although the
   !    detailed implementation might be slightly different)
   !  initialization (frac_init):
   !    afrac is set on all grids
   !      fractions_a(afrac) = 1.0
   !      fractions_o(afrac) = mapa2o(fractions_a(afrac))
   !      fractions_i(afrac) = mapa2i(fractions_a(afrac))
   !      fractions_l(afrac) = mapa2l(fractions_a(afrac))
   !    initially assume ifrac on all grids is zero
   !      fractions_*(ifrac) = 0.0
   !    fractions/masks provided by surface components
   !      fractions_o(ofrac) = dom_o(frac)  ! ocean "mask"
   !      fractions_l(lfrin) = dom_l(frac)  ! land model fraction
   !    then mapped to the atm model
   !      fractions_a(ofrac) = mapo2a(fractions_o(ofrac))
   !      fractions_a(lfrin) = mapl2a(fractions_l(lfrin))
   !    and a few things are then derived
   !      fractions_a(lfrac) = 1.0 - fractions_a(ofrac)
   !        this is truncated to zero for very small values (< 0.001)
   !        to attempt to preserve non-land gridcells.
   !      fractions_l(lfrac) = mapa2l(fractions_a(lfrac))
   !    one final term is computed
   !      dom_a(ascale) = fractions_a(lfrac)/fractions_a(lfrin)
   !      dom_l(ascale) = mapa2l(dom_a(ascale))
   !        these are used to correct land fluxes in budgets and lnd2rtm coupling
   !        and are particularly important when the land model is running on
   !        a different grid than the atm model.  in the old system, this term
   !        was treated as effectively 1.0 since there was always a check that
   !        fractions_a(lfrac) ~ fractions_a(lfrin), namely that the land model
   !        provided a land frac that complemented the ocean grid.  this is
   !        no longer a requirement in this new system and as a result, the
   !        ascale term can be thought of as a rescaling of the land fractions
   !        in the land model to be exactly complementary to the ocean model
   !        on whatever grid it may be running.
   !  run-time (frac_set):
   !    update fractions on ice grid
   !      fractions_i(ifrac) = i2x_i(Si_ifrac)  ! ice frac from ice model
   !      fractions_i(ofrac) = 1.0 - fractions_i(ifrac)
   !        note: the relative fractions are corrected to total fractions
   !      fractions_o(ifrac) = mapi2o(fractions_i(ifrac))
   !      fractions_o(ofrac) = mapi2o(fractions_i(ofrac))
   !      fractions_a(ifrac) = mapi2a(fractions_i(ifrac))
   !      fractions_a(ofrac) = mapi2a(fractions_i(ofrac))
   !
   !  fractions used in merging are as follows
   !    mrg_x2a uses fractions_a(lfrac,ofrac,ifrac)
   !    mrg_x2o needs to use fractions_o(ofrac,ifrac) normalized to one
   !      normalization happens in mrg routine
   !
   !  fraction corrections in mapping are as follows
   !    mapo2a uses *fractions_o(ofrac) and /fractions_a(ofrac)
   !    mapi2a uses *fractions_i(ifrac) and /fractions_a(ifrac)
   !    mapl2a uses *fractions_l(lfrin) and /fractions_a(lfrin)
   !    mapa2* should use *fractions_a(afrac) and /fractions_*(afrac) but this
   !      has been defered since the ratio always close to 1.0
   !
   !  budgets use the standard afrac, ofrac, ifrac, and lfrac to compute
   !    quantities except in the land budget which uses lfrin multiplied
   !    by the scale factor, dom_l(ascale) to compute budgets.
   !
   !  fraction and domain checks
   !    initialization:
   !      dom_i = mapo2i(dom_o)  ! lat, lon, mask, area
   !      where fractions_a(lfrac) > 0.0, fractions_a(lfrin) is also > 0.0
   !         this ensures the land will provide data everywhere the atm needs it
   !         and allows the land frac to be subtlely different from the
   !         land fraction specified in the atm.
   !      dom_a = mapl2a(dom_l)  ! if atm/lnd same grids
   !      dom_a = mapo2a(dom_o)  ! if atm/ocn same grids
   !      dom_a = mapi2a(dom_i)  ! if atm/ocn same grids
   !      0.0-eps < fractions_*(*) < 1.0+eps
   !      fractions_l(lfrin) = fractions_l(lfrac)
   !        only if atm/lnd same grids (but this is not formally required)
   !        this is needed until dom_l(ascale) is sent to the land model
   !        as an additional field for use in l2r mapping.
   !    run time:
   !      fractions_a(lfrac) + fractions_a(ofrac) + fractions_a(ifrac) ~ 1.0
   !      0.0-eps < fractions_*(*) < 1.0+eps
   
---------------
Domain Checking
---------------
Domain checking is a very important initialization step in the system.
Domain checking verifies that the longitudes, latitudes, areas, masks, and fractions of different grids are consistent with each other.
The subroutine that carries out domain checking is in ``$CIMEROOT/driver_cpl/driver/seq_domain_mct``.
Tolerances for checking the domains can be set in the drv_in driver namelist via the namelist variables, ``eps_frac``, ``eps_amask``, ``eps_agrid``, ``eps_aarea``, ``eps_omask``, ``eps_ogrid``, and ``eps_oarea``. 
These values are derived in the coupler namelist from the script env variables, EPS_FRAC, EPS_AMASK, EPS_AGRID, EPS_AAREA, EPS_OMASK, EPS_OGRID, and EPS_OAREA in the env_run.xml.
If an error is detected in the domain checking, the model will write an error message and abort.

The domain checking is dependent on the grids and in particular, the samegrid input namelist settings. But it basically does the following,
::

   ocean/ice grid comparison:
      verifies the grids are the same size
      verifies the difference in longitudes and latitudes is less than eps_ogrid.
      verifies the difference in masks is less than eps_omask
      verifies the difference in areas is less than eps_oarea

   atmosphere/land grid comparison (if samegrid_al):
      verifies the grids are the same size
      verifies the difference in longitudes and latitudes is less than eps_agrid.
      verifies the difference in masks is less than eps_amask
      verifies the difference in areas is less than eps_aarea

   atmosphere/ocean grid comparison (if samegrid_ao):
      verifies the grids are the same size
      verifies the difference in longitudes and latitudes is less than eps_agrid.
      verifies the difference in masks is less than eps_amask
      verifies the difference in areas is less than eps_aarea

   fractions
      verifies that the land fraction on the atmosphere grid and the ocean fraction 
      on the atmosphere grid add to one within a tolerance of eps_frac.

There are a number of subtle aspects in the domain checking like whether to check over masked grid cells, but these issues are less important than recognizing that errors in the domain checking should be treated seriously.
It is easy to make the errors go away by changing the tolerances, but by doing so, critical grid errors that can impact conservation and consistency in a simulation might be overlooked.


-----------------------
Mapping (Interpolation)
-----------------------
Mapping files to support interpolation of fields between grids are computed offline.
This is done using the ESMF offline regridding utility.
First, note that historically, the ocean grid has been the higher resolution grid.
While that is no longer always the case, the current implementation largely reflects that presumption.
In general, mapping of fluxes is done using a locally conservative area average approach to preserve conservation. 
State fields are generally mapped using bilinear interpolation from the atmosphere grid to the ocean grid to better preserve gradients, but state fields are generally mapped using the conservative area average approach from the ocean grid to the atmosphere grid.
But this is not a requirement of the system. 
The individual state and flux mapping files are specified at runtime using the ``seq_maps.rc`` input file, and any valid mapping file using any mapping approach can be specified in that input file.

The ``seq_maps.c`` file contains information about the mapping files as well as the mapping type. 
There are currently two types of mapping implementations, "X" and "Y".

- "X" mapping rearranges the source data to the destination grid decomposition and then a local mapping is done from the source to the destination grid on the destination decomposition. In the "X" type, the source grid is rearranged. 
- "Y" mapping does a local mapping from the source grid to the destination grid on the source grid decomposition. That generates a partial sum of the destination values which are then rearranged to the destination decomposition and summed. Both options produce reasonable results, although they may differ in value by "roundoff" due to differences in order or operations. The type chosen impacts performance. In both implementations, the number of flops is basically identical. The difference is the communication. In the "Y" type, the destination grid is rearranged. 

Since historically, the ocean grid is higher resolution than the atmosphere grid, "X" mapping is used for atmosphere to ocean/ice mapping and "Y" mapping is used from ocean/ice to atmosphere mapping to optimize mapping performance.

Mapping corrections are made in some cases in the polar region.
In particular, the current bilinear and area conservative mapping approaches introduce relatively large errors in mapping vector fields around the pole. 
The current coupler can correct the interpolated surface wind velocity near the pole when mapping from the atmosphere to the ocean and ice grids.
There are several options that correct the vector mapping and these are set in the env variable VECT_MAP. 
The npfix option only affects ocean and ice grid cells that are northward of the last latitude line of the atmospheric grid.
The algorithm is contained in the file models/drv/driver/map_atmocn_mct.F90 and is only valid when the atmosphere grid is a longitude/latitude grid. 
This feature is generally on by default.
The other alternative is the cart3d option which converts the surface u and v velocity to 3d x,y,z vectors then maps those three vectors before coverting back to u and v east and north directions on the surface. 
Both vector mapping methods introduce errors of different degrees but are generally much better than just mapping vector fields as if they were individual scalars. 
The ``vect_map`` namelist input is set in the ``drv_in`` file.

The input mapping files are assumed to be valid for grids with masks of value zero or one where grid points with a mask of zero are never considered in the mapping.
Well defined, locally conservative area mapping files as well as bilinear mapping files can be generated using this masked approach. 
However, there is another issue which is that a grid fraction in an active cell might actually change over time.
This is not the case for land fraction, but it is the case for relative ice and ocean fractions.
The ice fraction is constantly evolving in the system in general.
To improve the accuracy of the ice and ocean mapping, the ocean/ice fields are scaled by the local fraction before mapping and unscaled by the mapped fraction after mapping. 
The easiest way to demonstate this is via an example.
Consider a case where two ice cells of equal area underlie a single atmosphere cell completely. 
The mapping weight of each ice cell generated offline would be 0.5 in this case and if ice temperatures of -1.0 and -2.0 in the two cells respectively were mapped to the atmosphere grid, a resulting ice temperature on the atmosphere grid of -1.5 would result. 
Consider the case where one cell has an ice fraction of 0.3 and the other has a fraction of 0.5. 
Mapping the ice fraction to the atmospheric cell results in a value of 0.4. 
If the same temperatures are mapped in the same way, a temperature of -1.5 results which is reasonable, but not entirely accurate. 
Because of the relative ice fractions, the weight of the second cell should be greater than the weight of the first cell. 
Taking this into account properly results in a fraction weighted ice temperature of -1.625 in this example. 
This is the fraction correction that is carried out whenever ocean and ice fields are mapped to the atmosphere grid. 
Time varying fraction corrections are not required in other mappings to improve accuracy because their relative fractions remain static.

-------------------------
Area Correction of Fluxes
-------------------------
To improve conservation in the system, all fluxes sent to and received from components are corrected for the area differences between the components.
There are many reasonable ways to compute an area of a grid cell, but they are not generally consistent. 
One assumption with respect to conservation of fluxes is that the area acting upon the flux is well defined.
Differences in area calculations can result in differences of areas up to a few percent and if these are not corrected, will impact overall mass and heat conservation. 
Areas are extracted for each grid from the mapping files.
In this implementation, it is assumed that the areas in all mapping files are computed reasonably and consistently for each grid and on different grids. 
Those mapping areas are used to correct the fluxes for each component by scaling the fluxes sent to and received by the component by the ratio of the mapping area and the component area.
The areas from the components are provided to the coupler by the component at initialization. 
The minimum and maximum value of each area corrections is written to the coupler log file at initialization.
One critical point is that if mapping files are generated by different tools offline and used in the driver, an error could be introduced that is related to inconsistent areas provided by different mapping files.


