module Hydrology2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculation of soil/snow hydrology.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl,   only : iulog, use_vichydro, use_cn
  use abortutils,   only : endrun
  use decompMod   , only: bounds_type

  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: Hydrology2NoDrainage        ! Calculates soil/snow hydrology without drainage
  public :: Hydrology2Drainage        ! Calculates soil/snow hydrolog: drainage
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Hydrology2NoDrainage(bounds, &
       num_nolakec, filter_nolakec, &
       num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, &
       num_snowc, filter_snowc, &
       num_nosnowc, filter_nosnowc)
    !
    ! !DESCRIPTION:
    ! This is the main subroutine to execute the calculation of soil/snow
    ! hydrology
    ! Calling sequence is:
    !  Hydrology2:                 surface hydrology driver
    !    -> SnowWater:             change of snow mass and snow water onto soil
    !    -> SurfaceRunoff:         surface runoff
    !    -> Infiltration:          infiltration into surface soil layer
    !    -> SoilWater:             soil water movement between layers
    !          -> Tridiagonal      tridiagonal matrix solution
    !    -> Drainage:              subsurface runoff
    !    -> SnowCompaction:        compaction of snow layers
    !    -> CombineSnowLayers:     combine snow layers that are thinner than minimum
    !    -> DivideSnowLayers:      subdivide snow layers that are thicker than maximum
    !
    ! !USES:
    use clmtype
    use clm_atmlnd      , only : clm_a2l
    use clm_varcon      , only : denh2o, denice, istice, istwet, istsoil, istice_mec, spval, &
         icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall, &
         icol_shadewall, istdlak, tfrz, hfus, grav
    use clm_varcon      , only : istcrop
    use clm_varpar      , only : nlevgrnd, nlevsno, nlevsoi, nlevurb
    use SnowHydrologyMod, only : SnowCompaction, CombineSnowLayers, DivideSnowLayers, &
         SnowWater, BuildSnowFilter
    use SoilHydrologyMod, only : Infiltration, SoilWater, Drainage, SurfaceRunoff, WaterTable
    use clm_time_manager, only : get_step_size, get_nstep
    use CLMVICMapMod    , only : CLMVICMap
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(:)   ! column filter for non-lake points
    integer, intent(in) :: num_hydrologyc              ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:)! column filter for soil points
    integer, intent(in) :: num_urbanc                  ! number of column urban points in column filter
    integer, intent(in) :: filter_urbanc(:)    ! column filter for urban points
    integer  :: num_snowc                  ! number of column snow points
    integer  :: filter_snowc(:)    ! column filter for snow points
    integer  :: num_nosnowc                ! number of column non-snow points
    integer  :: filter_nosnowc(:)  ! column filter for non-snow points
    !
    ! !LOCAL VARIABLES:
    integer  :: g,l,c,j,fc                 ! indices
    integer  :: nstep                      ! time step number
    real(r8) :: dtime                      ! land model time step (sec)
    real(r8) :: psi,vwc,fsattmp,psifrz     ! temporary variables for soilpsi calculation
    real(r8) :: watdry                     ! temporary
    real(r8) :: rwat(bounds%begc:bounds%endc)              ! soil water wgted by depth to maximum depth of 0.5 m
    real(r8) :: swat(bounds%begc:bounds%endc)              ! same as rwat but at saturation
    real(r8) :: rz(bounds%begc:bounds%endc)                ! thickness of soil layers contributing to rwat (m)
    real(r8) :: tsw                        ! volumetric soil water to 0.5 m
    real(r8) :: stsw                       ! volumetric soil water to 0.5 m at saturation
    real(r8) :: snowmass                   ! liquid+ice snow mass in a layer [kg/m2]
    real(r8) :: snowcap_scl_fct            ! temporary factor used to correct for snow capping
    real(r8) :: fracl                      ! fraction of soil layer contributing to 10cm total soil water
    real(r8) :: s_node                     ! soil wetness (-)
    real(r8) :: icefrac(bounds%begc:bounds%endc,1:nlevsoi)
    !-----------------------------------------------------------------------

   associate(& 
   ityplun              => lun%itype               , & ! Input:  [integer (:)]  landunit type                            
   urbpoi               => lun%urbpoi              , & ! Input:  [logical (:)]  true => landunit is an urban point       
   snow_depth           => cps%snow_depth          , & ! Input:  [real(r8) (:)] snow height of snow covered area (m)     
   snowdp               => cps%snowdp              , & ! Input:  [real(r8) (:)]  gridcell averaged snow height (m)       
   frac_sno_eff         => cps%frac_sno_eff        , & ! Input:  [real(r8) (:)] eff.  snow cover fraction (col) [frc]    
   qflx_evap_soi        => pwf_a%qflx_evap_soi     , & ! Input:  [real(r8) (:)]  soil evaporation                        
   h2osfc               => cws%h2osfc              , & ! Input:  [real(r8) (:)]  surface water (mm)                      
   frac_h2osfc          => cps%frac_h2osfc         , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   t_h2osfc             => ces%t_h2osfc            , & ! Input:  [real(r8) (:)]  surface water temperature               
   qflx_drain_perched   => cwf%qflx_drain_perched  , & ! Input:  [real(r8) (:)]  sub-surface runoff from perched zwt (mm H2O /s)
   qflx_floodg          => clm_a2l%forc_flood      , & ! Input:  [real(r8) (:)]  gridcell flux of flood water from RTM   
   qflx_h2osfc_surf     => cwf%qflx_h2osfc_surf    , & ! Input:  [real(r8) (:)] surface water runoff (mm/s)              
   cactive              => col%active              , & ! Input:  [logical (:)]  true=>do computations on this column 
   cgridcell            => col%gridcell            , & ! Input:  [integer (:)]  column's gridcell                        
   clandunit            => col%landunit            , & ! Input:  [integer (:)]  column's landunit                        
   ctype                => col%itype               , & ! Input:  [integer (:)]  column type                              
   snl                  => cps%snl                 , & ! Input:  [integer (:)]  number of snow layers                    
   t_grnd               => ces%t_grnd              , & ! Output: [real(r8) (:)]  ground temperature (Kelvin)             
   h2ocan               => pws_a%h2ocan            , & ! Input:  [real(r8) (:)]  canopy water (mm H2O)                   
   h2osno               => cws%h2osno              , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                     
   wf                   => cps%wf                  , & ! Output: [real(r8) (:)]  soil water as frac. of whc for top 0.05 m added by F. Li and S. Levis 
   wf2                  => cps%wf2                 , & ! Output: [real(r8) (:)]  soil water as frac. of whc for top 0.17 m added by F. Li and S. Levis 
   snowice              => cws%snowice             , & ! Output: [real(r8) (:)]  average snow ice lens                   
   snowliq              => cws%snowliq             , & ! Output: [real(r8) (:)]  average snow liquid water               
   zwt                  => cws%zwt                 , & ! Input:  [real(r8) (:)]  water table depth (m)                   
   fcov                 => cws%fcov                , & ! Input:  [real(r8) (:)]  fractional impermeable area             
   fsat                 => cws%fsat                , & ! Input:  [real(r8) (:)]  fractional area with water table at surface
   wa                   => cws%wa                  , & ! Input:  [real(r8) (:)]  water in the unconfined aquifer (mm)    
   qcharge              => cws%qcharge             , & ! Input:  [real(r8) (:)]  aquifer recharge rate (mm/s)            
   watsat               => cps%watsat              , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)
   sucsat               => cps%sucsat              , & ! Input:  [real(r8) (:,:)]  minimum soil suction (mm)             
   bsw                  => cps%bsw                 , & ! Input:  [real(r8) (:,:)]  Clapp and Hornberger "b"              
   z                    => cps%z                   , & ! Input:  [real(r8) (:,:)]  layer depth  (m)                      
   dz                   => cps%dz                  , & ! Input:  [real(r8) (:,:)]  layer thickness depth (m)             
   zi                   => cps%zi                  , & ! Input:  [real(r8) (:,:)]  interface depth (m)                   
   t_soisno             => ces%t_soisno            , & ! Output: [real(r8) (:,:)]  soil temperature (Kelvin)             
   h2osoi_ice           => cws%h2osoi_ice          , & ! Output: [real(r8) (:,:)]  ice lens (kg/m2)                      
   h2osoi_liq           => cws%h2osoi_liq          , & ! Output: [real(r8) (:,:)]  liquid water (kg/m2)                  
   h2osoi_vol           => cws%h2osoi_vol          , & ! Output: [real(r8) (:,:)]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   t_soi_10cm           => ces%t_soi_10cm          , & ! Output: [real(r8) (:)]  soil temperature in top 10cm of soil (Kelvin)
   tsoi17               => ces%tsoi17              , & ! Output: [real(r8) (:)]  soil temperature in top 17cm of soil (Kelvin) added by F. Li and S. Levis
   h2osoi_liqice_10cm   => cws%h2osoi_liqice_10cm  , & ! Output: [real(r8) (:)]  liquid water + ice lens in top 10cm of soil (kg/m2)
   qflx_evap_tot        => pwf_a%qflx_evap_tot     , & ! Input:  [real(r8) (:)]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
   qflx_drain           => cwf%qflx_drain          , & ! Output: [real(r8) (:)]  sub-surface runoff (mm H2O /s)          
   qflx_surf            => cwf%qflx_surf           , & ! Output: [real(r8) (:)]  surface runoff (mm H2O /s)              
   qflx_infl            => cwf%qflx_infl           , & ! Output: [real(r8) (:)]  infiltration (mm H2O /s)                
   qflx_qrgwl           => cwf%qflx_qrgwl          , & ! Output: [real(r8) (:)]  qflx_surf at glaciers, wetlands, lakes  
   endwb                => cwbal%endwb             , & ! Output: [real(r8) (:)]  water mass end of the time step         
   begwb                => cwbal%begwb             , & ! Input:  [real(r8) (:)]  water mass begining of the time step    
   soilpsi              => cps%soilpsi             , & ! Output: [real(r8) (:,:)]  soil water potential in each soil layer (MPa)
   smp_l                => cws%smp_l               , & ! Input:  [real(r8) (:,:)]  soil matrix potential [mm]            
   hk_l                 => cws%hk_l                , & ! Input:  [real(r8) (:,:)]  hydraulic conductivity (mm/s)         
   qflx_rsub_sat        => cwf%qflx_rsub_sat       , & ! Input:  [real(r8) (:)]  soil saturation excess [mm h2o/s]       
   qflx_runoff          => cwf%qflx_runoff         , & ! Output: [real(r8) (:)]  total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
   qflx_runoff_u        => cwf%qflx_runoff_u       , & ! Output: [real(r8) (:)]  Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
   qflx_runoff_r        => cwf%qflx_runoff_r       , & ! Output: [real(r8) (:)]  Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
   t_grnd_u             => ces%t_grnd_u            , & ! Output: [real(r8) (:)]  Urban ground temperature (Kelvin)       
   t_grnd_r             => ces%t_grnd_r            , & ! Output: [real(r8) (:)]  Rural ground temperature (Kelvin)       
   snot_top             => cps%snot_top            , & ! Output: [real(r8) (:)]  snow temperature in top layer (col) [K] 
   dTdz_top             => cps%dTdz_top            , & ! Output: [real(r8) (:)]  temperature gradient in top layer (col) [K m-1]
   snw_rds              => cps%snw_rds             , & ! Output: [real(r8) (:,:)]  effective snow grain radius (col,lyr) [microns, m^-6]
   snw_rds_top          => cps%snw_rds_top         , & ! Output: [real(r8) (:)]  effective snow grain size, top layer(col) [microns]
   sno_liq_top          => cps%sno_liq_top         , & ! Output: [real(r8) (:)]  liquid water fraction in top snow layer (col) [frc]
   frac_sno             => cps%frac_sno            , & ! Output: [real(r8) (:)]  snow cover fraction (col) [frc]         
   h2osno_top           => cps%h2osno_top          , & ! Output: [real(r8) (:)]  mass of snow in top layer (col) [kg]    
   mss_bcpho            => cps%mss_bcpho           , & ! Output: [real(r8) (:,:)]  mass of hydrophobic BC in snow (col,lyr) [kg]
   mss_bcphi            => cps%mss_bcphi           , & ! Output: [real(r8) (:,:)]  mass of hydrophillic BC in snow (col,lyr) [kg]
   mss_bctot            => cps%mss_bctot           , & ! Output: [real(r8) (:,:)]  total mass of BC (pho+phi) (col,lyr) [kg]
   mss_bc_col           => cps%mss_bc_col          , & ! Output: [real(r8) (:)]  total mass of BC in snow column (col) [kg]
   mss_bc_top           => cps%mss_bc_top          , & ! Output: [real(r8) (:)]  total mass of BC in top snow layer (col) [kg]
   mss_cnc_bcphi        => cps%mss_cnc_bcphi       , & ! Output: [real(r8) (:,:)]  mass concentration of BC species 1 (col,lyr) [kg/kg]
   mss_cnc_bcpho        => cps%mss_cnc_bcpho       , & ! Output: [real(r8) (:,:)]  mass concentration of BC species 2 (col,lyr) [kg/kg]
   mss_ocpho            => cps%mss_ocpho           , & ! Output: [real(r8) (:,:)]  mass of hydrophobic OC in snow (col,lyr) [kg]
   mss_ocphi            => cps%mss_ocphi           , & ! Output: [real(r8) (:,:)]  mass of hydrophillic OC in snow (col,lyr) [kg]
   mss_octot            => cps%mss_octot           , & ! Output: [real(r8) (:,:)]  total mass of OC (pho+phi) (col,lyr) [kg]
   mss_oc_col           => cps%mss_oc_col          , & ! Output: [real(r8) (:)]  total mass of OC in snow column (col) [kg]
   mss_oc_top           => cps%mss_oc_top          , & ! Output: [real(r8) (:)]  total mass of OC in top snow layer (col) [kg]
   mss_cnc_ocphi        => cps%mss_cnc_ocphi       , & ! Output: [real(r8) (:,:)]  mass concentration of OC species 1 (col,lyr) [kg/kg]
   mss_cnc_ocpho        => cps%mss_cnc_ocpho       , & ! Output: [real(r8) (:,:)]  mass concentration of OC species 2 (col,lyr) [kg/kg]
   mss_dst1             => cps%mss_dst1            , & ! Output: [real(r8) (:,:)]  mass of dust species 1 in snow (col,lyr) [kg]
   mss_dst2             => cps%mss_dst2            , & ! Output: [real(r8) (:,:)]  mass of dust species 2 in snow (col,lyr) [kg]
   mss_dst3             => cps%mss_dst3            , & ! Output: [real(r8) (:,:)]  mass of dust species 3 in snow (col,lyr) [kg]
   mss_dst4             => cps%mss_dst4            , & ! Output: [real(r8) (:,:)]  mass of dust species 4 in snow (col,lyr) [kg]
   mss_dsttot           => cps%mss_dsttot          , & ! Output: [real(r8) (:,:)]  total mass of dust in snow (col,lyr) [kg]
   mss_dst_col          => cps%mss_dst_col         , & ! Output: [real(r8) (:)]  total mass of dust in snow column (col) [kg]
   mss_dst_top          => cps%mss_dst_top         , & ! Output: [real(r8) (:)]  total mass of dust in top snow layer (col) [kg]
   mss_cnc_dst1         => cps%mss_cnc_dst1        , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 1 (col,lyr) [kg/kg]
   mss_cnc_dst2         => cps%mss_cnc_dst2        , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 2 (col,lyr) [kg/kg]
   mss_cnc_dst3         => cps%mss_cnc_dst3        , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 3 (col,lyr) [kg/kg]
   mss_cnc_dst4         => cps%mss_cnc_dst4        , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 4 (col,lyr) [kg/kg]
   do_capsnow           => cps%do_capsnow          , & ! Output: [logical (:)]  true => do snow capping                  
   qflx_snwcp_ice       => pwf_a%qflx_snwcp_ice    , & ! Output: [real(r8) (:)]  excess snowfall due to snow capping (mm H2O /s) [+]`
   smpmin               => cps%smpmin              , & ! Input:  [real(r8) (:)]  restriction for min of soil potential (mm)
   snow_persistence     => cps%snow_persistence      & ! Output: [real(r8) (:)]  counter for length of time snow-covered
   )

    ! Determine time step and step size

    nstep = get_nstep()
    dtime = get_step_size()

    ! Determine initial snow/no-snow filters (will be modified possibly by
    ! routines CombineSnowLayers and DivideSnowLayers below

    call BuildSnowFilter(bounds, num_nolakec, filter_nolakec, &
         num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! Determine the change of snow mass and the snow water onto soil

    call SnowWater(bounds, num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! Determine soil hydrology
    if (use_vichydro) then
       ! mapping soilmoist from CLM to VIC layers for runoff calculations
       call CLMVICMap(bounds, num_hydrologyc, filter_hydrologyc)
    end if

    call SurfaceRunoff(bounds, num_hydrologyc, filter_hydrologyc, &
                       num_urbanc, filter_urbanc)

    call Infiltration(bounds,  num_hydrologyc, filter_hydrologyc, &
                      num_urbanc, filter_urbanc)

    call SoilWater(bounds, num_hydrologyc, filter_hydrologyc, &
                   num_urbanc, filter_urbanc)

    if (use_vichydro) then
       ! mapping soilmoist from CLM to VIC layers for runoff calculations
       call CLMVICMap(bounds, num_hydrologyc, filter_hydrologyc)
    end if

    call WaterTable(bounds, num_hydrologyc, filter_hydrologyc, &
                  num_urbanc, filter_urbanc)
                  
    ! Natural compaction and metamorphosis.

    call SnowCompaction(bounds, num_snowc, filter_snowc)
    
    ! Combine thin snow elements
    
    call CombineSnowLayers(bounds, num_snowc, filter_snowc)
    
    ! Divide thick snow elements
    
    call DivideSnowLayers(bounds, num_snowc, filter_snowc)

    ! Set empty snow layers to zero

    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j <= snl(c) .and. snl(c) > -nlevsno) then
             h2osoi_ice(c,j) = 0._r8
             h2osoi_liq(c,j) = 0._r8
             t_soisno(c,j) = 0._r8
             dz(c,j) = 0._r8
             z(c,j) = 0._r8
             zi(c,j-1) = 0._r8
          end if
       end do
    end do

    ! Build new snow filter

    call BuildSnowFilter(bounds, num_nolakec, filter_nolakec, &
         num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! For columns where snow exists, accumulate 'time-covered-by-snow' counters.
    ! Otherwise, re-zero counter, since it is bareland
    do fc = 1, num_snowc
      c = filter_snowc(fc)
      snow_persistence(c) = snow_persistence(c) + dtime
    end do
    do fc = 1, num_nosnowc
      c = filter_nosnowc(fc)
      snow_persistence(c) = 0._r8
    enddo

    ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
    ! over all snow layers for history output

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       snowice(c) = 0._r8
       snowliq(c) = 0._r8
    end do

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             snowice(c) = snowice(c) + h2osoi_ice(c,j)
             snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
          end if
       end do
    end do

    ! Calculate column average snow depth
    do c = bounds%begc,bounds%endc
       snowdp(c) = snow_depth(c) * frac_sno_eff(c)
    end do

    ! Determine ground temperature, ending water balance and volumetric soil water
    ! Calculate soil temperature and total water (liq+ice) in top 10cm of soil
    ! Calculate soil temperature and total water (liq+ice) in top 17cm of soil
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       if (.not. urbpoi(l)) then
          t_soi_10cm(c) = 0._r8
          tsoi17(c) = 0._r8
          h2osoi_liqice_10cm(c) = 0._r8
       end if
    end do
    do j = 1, nlevsoi
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (.not. urbpoi(l)) then
             ! soil T at top 17 cm added by F. Li and S. Levis
            if (zi(c,j) <= 0.17_r8) then
              fracl = 1._r8
              tsoi17(c) = tsoi17(c) + t_soisno(c,j)*dz(c,j)*fracl
            else
              if (zi(c,j) > 0.17_r8 .and. zi(c,j-1) .lt. 0.17_r8) then 
                fracl = (0.17_r8 - zi(c,j-1))/dz(c,j)
                tsoi17(c) = tsoi17(c) + t_soisno(c,j)*dz(c,j)*fracl
               end if
            end if

            if (zi(c,j) <= 0.1_r8) then
              fracl = 1._r8
              t_soi_10cm(c) = t_soi_10cm(c) + t_soisno(c,j)*dz(c,j)*fracl
              h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
              (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                                       fracl
            else
              if (zi(c,j) > 0.1_r8 .and. zi(c,j-1) .lt. 0.1_r8) then
                 fracl = (0.1_r8 - zi(c,j-1))/dz(c,j)
                 t_soi_10cm(c) = t_soi_10cm(c) + t_soisno(c,j)*dz(c,j)*fracl
                 h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
                 (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                                          fracl
              end if
            end if
          end if
       end do
    end do

    do fc = 1, num_nolakec
       
       c = filter_nolakec(fc)
       l = clandunit(c)

       ! t_grnd is weighted average of exposed soil and snow
       if (snl(c) < 0) then
          t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
               + (1 - frac_sno_eff(c)- frac_h2osfc(c)) * t_soisno(c,1) &
               + frac_h2osfc(c) * t_h2osfc(c)
       else
          t_grnd(c) = (1 - frac_h2osfc(c)) * t_soisno(c,1) + frac_h2osfc(c) * t_h2osfc(c)
       endif

       if (urbpoi(l)) then
          t_grnd_u(c) = t_soisno(c,snl(c)+1)
       else
          t_soi_10cm(c) = t_soi_10cm(c)/0.1_r8
           tsoi17(c) =  tsoi17(c)/0.17_r8         ! F. Li and S. Levis
       end if
       if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
         t_grnd_r(c) = t_soisno(c,snl(c)+1)
       end if

    end do

    do j = 1, nlevgrnd
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
               .or. ctype(c) == icol_roof) .and. j > nlevurb) then
          else
            h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
          end if
       end do
    end do

    if (use_cn) then
       ! Update soilpsi.
       ! ZMS: Note this could be merged with the following loop updating smp_l in the future.
       do j = 1, nlevgrnd
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)

             if (h2osoi_liq(c,j) > 0._r8) then

                vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)

                ! the following limit set to catch very small values of 
                ! fractional saturation that can crash the calculation of psi

                ! use the same contants used in the supercool so that psi for frozen soils is consistent
                fsattmp = max(vwc/watsat(c,j), 0.001_r8)
                psi = sucsat(c,j) * (-9.8e-6_r8) * (fsattmp)**(-bsw(c,j))  ! Mpa
                soilpsi(c,j) = min(max(psi,-15.0_r8),0._r8)

             else 
                soilpsi(c,j) = -15.0_r8
             end if
          end do
       end do
    end if

    ! Update smp_l for history and for ch4Mod.
    ! ZMS: Note, this form, which seems to be the same as used in SoilWater, DOES NOT distinguish between
    ! ice and water volume, in contrast to the soilpsi calculation above. It won't be used in ch4Mod if
    ! t_soisno <= tfrz, though.
    do j = 1, nlevgrnd
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
          s_node = min(1.0_r8, s_node)

          smp_l(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
          smp_l(c,j) = max(smpmin(c), smp_l(c,j))
       end do
    end do

    if (use_cn) then
       ! Available soil water up to a depth of 0.05 m.
       ! Potentially available soil water (=whc) up to a depth of 0.05 m.
       ! Water content as fraction of whc up to a depth of 0.05 m.
       
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          rwat(c) = 0._r8
          swat(c) = 0._r8
          rz(c)   = 0._r8
       end do

       do j = 1, nlevgrnd
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             !if (z(c,j)+0.5_r8*dz(c,j) <= 0.5_r8) then
             if (z(c,j)+0.5_r8*dz(c,j) <= 0.05_r8) then
                watdry = watsat(c,j) * (316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))
                rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
                swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
                rz(c) = rz(c) + dz(c,j)
             end if
          end do
       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (rz(c) /= 0._r8) then
             tsw  = rwat(c)/rz(c)
             stsw = swat(c)/rz(c)
          else
             watdry = watsat(c,1) * (316230._r8/sucsat(c,1)) ** (-1._r8/bsw(c,1))
             tsw = h2osoi_vol(c,1) - watdry
             stsw = watsat(c,1) - watdry
          end if
          wf(c) = tsw/stsw
       end do

       do j = 1, nlevgrnd
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             if (z(c,j)+0.5_r8*dz(c,j) <= 0.17_r8) then
                watdry = watsat(c,j) * (316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))
                rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
                swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
                rz(c) = rz(c) + dz(c,j)
             end if
          end do
       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (rz(c) /= 0._r8) then
             tsw  = rwat(c)/rz(c)
             stsw = swat(c)/rz(c)
          else
             watdry = watsat(c,1) * (316230._r8/sucsat(c,1)) ** (-1._r8/bsw(c,1))
             tsw = h2osoi_vol(c,1) - watdry
             stsw = watsat(c,1) - watdry
          end if
          wf2(c) = tsw/stsw
       end do
    end if

    !  Calculate column-integrated aerosol masses, and
    !  mass concentrations for radiative calculations and output
    !  (based on new snow level state, after SnowFilter is rebuilt.
    !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there 
    !  can be zero snow layers but an active column in filter)

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       ! Zero column-integrated aerosol mass before summation
       mss_bc_col(c)  = 0._r8
       mss_oc_col(c)  = 0._r8
       mss_dst_col(c) = 0._r8

       do j = -nlevsno+1, 0

          ! layer mass of snow:
          snowmass = h2osoi_ice(c,j)+h2osoi_liq(c,j)

          ! Correct the top layer aerosol mass to account for snow capping. 
          ! This approach conserves the aerosol mass concentration
          ! (but not the aerosol amss) when snow-capping is invoked

          if (j == snl(c)+1) then
             if (do_capsnow(c)) then
                snowcap_scl_fct = snowmass / (snowmass+(qflx_snwcp_ice(c)*dtime))

                mss_bcpho(c,j) = mss_bcpho(c,j)*snowcap_scl_fct
                mss_bcphi(c,j) = mss_bcphi(c,j)*snowcap_scl_fct
                mss_ocpho(c,j) = mss_ocpho(c,j)*snowcap_scl_fct
                mss_ocphi(c,j) = mss_ocphi(c,j)*snowcap_scl_fct
                
                mss_dst1(c,j)  = mss_dst1(c,j)*snowcap_scl_fct
                mss_dst2(c,j)  = mss_dst2(c,j)*snowcap_scl_fct
                mss_dst3(c,j)  = mss_dst3(c,j)*snowcap_scl_fct
                mss_dst4(c,j)  = mss_dst4(c,j)*snowcap_scl_fct 
             endif
          endif

          if (j >= snl(c)+1) then
             mss_bctot(c,j)     = mss_bcpho(c,j) + mss_bcphi(c,j)
             mss_bc_col(c)      = mss_bc_col(c)  + mss_bctot(c,j)
             mss_cnc_bcphi(c,j) = mss_bcphi(c,j) / snowmass
             mss_cnc_bcpho(c,j) = mss_bcpho(c,j) / snowmass

             mss_octot(c,j)     = mss_ocpho(c,j) + mss_ocphi(c,j)
             mss_oc_col(c)      = mss_oc_col(c)  + mss_octot(c,j)
             mss_cnc_ocphi(c,j) = mss_ocphi(c,j) / snowmass
             mss_cnc_ocpho(c,j) = mss_ocpho(c,j) / snowmass
             
             mss_dsttot(c,j)    = mss_dst1(c,j)  + mss_dst2(c,j) + mss_dst3(c,j) + mss_dst4(c,j)
             mss_dst_col(c)     = mss_dst_col(c) + mss_dsttot(c,j)
             mss_cnc_dst1(c,j)  = mss_dst1(c,j)  / snowmass
             mss_cnc_dst2(c,j)  = mss_dst2(c,j)  / snowmass
             mss_cnc_dst3(c,j)  = mss_dst3(c,j)  / snowmass
             mss_cnc_dst4(c,j)  = mss_dst4(c,j)  / snowmass
         
          else
             !set variables of empty snow layers to zero
             snw_rds(c,j)       = 0._r8

             mss_bcpho(c,j)     = 0._r8
             mss_bcphi(c,j)     = 0._r8
             mss_bctot(c,j)     = 0._r8
             mss_cnc_bcphi(c,j) = 0._r8
             mss_cnc_bcpho(c,j) = 0._r8

             mss_ocpho(c,j)     = 0._r8
             mss_ocphi(c,j)     = 0._r8
             mss_octot(c,j)     = 0._r8
             mss_cnc_ocphi(c,j) = 0._r8
             mss_cnc_ocpho(c,j) = 0._r8

             mss_dst1(c,j)      = 0._r8
             mss_dst2(c,j)      = 0._r8
             mss_dst3(c,j)      = 0._r8
             mss_dst4(c,j)      = 0._r8
             mss_dsttot(c,j)    = 0._r8
             mss_cnc_dst1(c,j)  = 0._r8
             mss_cnc_dst2(c,j)  = 0._r8
             mss_cnc_dst3(c,j)  = 0._r8
             mss_cnc_dst4(c,j)  = 0._r8
          endif
       enddo
       
       ! top-layer diagnostics
       h2osno_top(c)  = h2osoi_ice(c,snl(c)+1) + h2osoi_liq(c,snl(c)+1)
       mss_bc_top(c)  = mss_bctot(c,snl(c)+1)
       mss_oc_top(c)  = mss_octot(c,snl(c)+1)
       mss_dst_top(c) = mss_dsttot(c,snl(c)+1)
    enddo
    
    ! Zero mass variables in columns without snow
    do fc = 1, num_nosnowc
       c = filter_nosnowc(fc)
            
       h2osno_top(c)      = 0._r8
       snw_rds(c,:)       = 0._r8

       mss_bc_top(c)      = 0._r8
       mss_bc_col(c)      = 0._r8    
       mss_bcpho(c,:)     = 0._r8
       mss_bcphi(c,:)     = 0._r8
       mss_bctot(c,:)     = 0._r8
       mss_cnc_bcphi(c,:) = 0._r8
       mss_cnc_bcpho(c,:) = 0._r8

       mss_oc_top(c)      = 0._r8
       mss_oc_col(c)      = 0._r8    
       mss_ocpho(c,:)     = 0._r8
       mss_ocphi(c,:)     = 0._r8
       mss_octot(c,:)     = 0._r8
       mss_cnc_ocphi(c,:) = 0._r8
       mss_cnc_ocpho(c,:) = 0._r8

       mss_dst_top(c)     = 0._r8
       mss_dst_col(c)     = 0._r8
       mss_dst1(c,:)      = 0._r8
       mss_dst2(c,:)      = 0._r8
       mss_dst3(c,:)      = 0._r8
       mss_dst4(c,:)      = 0._r8
       mss_dsttot(c,:)    = 0._r8
       mss_cnc_dst1(c,:)  = 0._r8
       mss_cnc_dst2(c,:)  = 0._r8
       mss_cnc_dst3(c,:)  = 0._r8
       mss_cnc_dst4(c,:)  = 0._r8

       ! top-layer diagnostics (spval is not averaged when computing history fields)
       snot_top(c)        = spval
       dTdz_top(c)        = spval
       snw_rds_top(c)     = spval
       sno_liq_top(c)     = spval
    enddo

    end associate 

  end subroutine Hydrology2NoDrainage

  !-----------------------------------------------------------------------
  subroutine Hydrology2Drainage(bounds, &
       num_nolakec, filter_nolakec, &
       num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, &
       num_do_smb_c, filter_do_smb_c)
    !
    ! !DESCRIPTION:
    ! This is the main subroutine to execute the calculation of soil/snow
    ! hydrology
    ! muszala March 7 2013 - description needs updating after Jinyun mods.
    ! Calling sequence is:
    !  Hydrology2:                 surface hydrology driver
    !    -> Drainage:              subsurface runoff
    !
    ! !USES:
    use clmtype
    use clm_atmlnd      , only : clm_a2l, a2l_downscaled_col
    use clm_varcon      , only : istice, istwet, istsoil, istice_mec, spval, &
         icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall, &
         icol_shadewall, denh2o, denice
    use clm_varcon      , only : istcrop, secspday
    use clm_varctl      , only : glc_dyn_runoff_routing, glc_snow_persistence_max_days
    use clm_varctl      , only : use_vichydro
    use clm_varpar      , only : nlevgrnd, nlevurb
    use SoilHydrologyMod, only : Drainage
    use clm_time_manager, only : get_step_size, get_nstep
    use CLMVICMapMod    , only : CLMVICMap

    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds    ! bounds
    integer, intent(in) :: num_nolakec         ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(:)   ! column filter for non-lake points
    integer, intent(in) :: num_hydrologyc      ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:)! column filter for soil points
    integer, intent(in) :: num_urbanc          ! number of column urban points in column filter
    integer, intent(in) :: filter_urbanc(:)    ! column filter for urban points
    integer, intent(in) :: num_do_smb_c        ! number of bareland columns in which SMB is calculated, in column filter    
    integer, intent(in) :: filter_do_smb_c(:)  ! column filter for bare land SMB columns      
    !
    ! !LOCAL VARIABLES:
    integer  :: g,l,c,j,fc                 ! indices
    integer  :: nstep                      ! time step number
    real(r8) :: dtime                      ! land model time step (sec)
    real(r8) :: vol_liq(bounds%begc:bounds%endc,1:nlevgrnd)! partial volume of liquid water in layer
    real(r8) :: hk(bounds%begc:bounds%endc,1:nlevgrnd)     ! hydraulic conductivity (mm h2o/s)
    !-----------------------------------------------------------------------

   associate(&    
   forc_rain                           =>    a2l_downscaled_col%forc_rain , & ! Input:  [real(r8) (:)]  rain rate [mm/s]                                  
   forc_snow                           =>    a2l_downscaled_col%forc_snow , & ! Input:  [real(r8) (:)]  snow rate [mm/s]                                  
   ityplun                             =>    lun%itype                    , & ! Input:  [integer (:)]  landunit type                                      
   urbpoi                              =>   lun%urbpoi                    , & ! Output: [logical (:)]  true => landunit is an urban point                 
   dz                                  =>    cps%dz                       , & ! Input:  [real(r8) (:,:)]  layer thickness depth (m)                       
   h2osfc                              =>    cws%h2osfc                   , & ! Input:  [real(r8) (:)]  surface water (mm)                                
   qflx_drain_perched                  =>    cwf%qflx_drain_perched       , & ! Input:  [real(r8) (:)]  sub-surface runoff from perched zwt (mm H2O /s)   
   qflx_floodg                         =>    clm_a2l%forc_flood           , & ! Input:  [real(r8) (:)]  gridcell flux of flood water from RTM             
   cactive                             =>    col%active                   , & ! Input:  [logical (:)]  true=>do computations on this column 
   qflx_h2osfc_surf                    =>    cwf%qflx_h2osfc_surf         , & ! Input:  [real(r8) (:)] surface water runoff (mm/s)                        
   cgridcell                           =>   col%gridcell                  , & ! Input:  [integer (:)]  column's gridcell                                  
   clandunit                           =>   col%landunit                  , & ! Input:  [integer (:)]  column's landunit                                  
   ctype                               =>    col%itype                    , & ! Input:  [integer (:)]  column type                                        
   h2ocan                              =>    pws_a%h2ocan                 , & ! Input:  [real(r8) (:)]  canopy water (mm H2O)                             
   h2osno                              =>    cws%h2osno                   , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                               
   fcov                                =>    cws%fcov                     , & ! Input:  [real(r8) (:)]  fractional impermeable area                       
   fsat                                =>    cws%fsat                     , & ! Input:  [real(r8) (:)]  fractional area with water table at surface       
   wa                                  =>    cws%wa                       , & ! Input:  [real(r8) (:)]  water in the unconfined aquifer (mm)              
   qcharge                             =>    cws%qcharge                  , & ! Input:  [real(r8) (:)]  aquifer recharge rate (mm/s)                      
   h2osoi_ice                          =>    cws%h2osoi_ice               , & ! Output: [real(r8) (:,:)]  ice lens (kg/m2)                                
   h2osoi_liq                          =>    cws%h2osoi_liq               , & ! Output: [real(r8) (:,:)]  liquid water (kg/m2)                            
   h2osoi_vol                          =>    cws%h2osoi_vol               , & ! Output: [real(r8) (:,:)]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   qflx_evap_tot                       =>    pwf_a%qflx_evap_tot          , & ! Input:  [real(r8) (:)]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg     
   qflx_drain                          =>    cwf%qflx_drain               , & ! Output: [real(r8) (:)]  sub-surface runoff (mm H2O /s)                    
   qflx_surf                           =>    cwf%qflx_surf                , & ! Output: [real(r8) (:)]  surface runoff (mm H2O /s)                        
   qflx_infl                           =>    cwf%qflx_infl                , & ! Output: [real(r8) (:)]  infiltration (mm H2O /s)                          
   qflx_qrgwl                          =>    cwf%qflx_qrgwl               , & ! Output: [real(r8) (:)]  qflx_surf at glaciers, wetlands, lakes            
   qflx_irrig                          =>    pwf_a%qflx_irrig             , & ! Input:  [real(r8) (:)]  irrigation flux (mm H2O /s)                       
   endwb                               =>    cwbal%endwb                  , & ! Output: [real(r8) (:)]  water mass end of the time step                   
   begwb                               =>    cwbal%begwb                  , & ! Input:  [real(r8) (:)]  water mass begining of the time step              
   qflx_rsub_sat                       =>    cwf%qflx_rsub_sat            , & ! Input:  [real(r8) (:)]  soil saturation excess [mm h2o/s]                 
   qflx_runoff                         =>    cwf%qflx_runoff              , & ! Output: [real(r8) (:)]  total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
   qflx_runoff_u                       =>    cwf%qflx_runoff_u            , & ! Output: [real(r8) (:)]  Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
   qflx_runoff_r                       =>    cwf%qflx_runoff_r            , & ! Output: [real(r8) (:)]  Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
   qflx_snwcp_ice                      =>    pwf_a%qflx_snwcp_ice         , & ! Output: [real(r8) (:)]  excess snowfall due to snow capping (mm H2O /s) [+]`
   qflx_glcice                         =>    cwf%qflx_glcice              , & ! Output: [real(r8) (:)]  flux of new glacier ice (mm H2O /s)               
   qflx_glcice_frz                     =>    cwf%qflx_glcice_frz          , & ! Output: [real(r8) (:)]  ice growth (positive definite) (mm H2O/s)
   qflx_glcice_melt                    =>    cwf%qflx_glcice_melt         , & ! Input:  [real(r8) (:)]  ice melt (positive definite) (mm H2O/s)      
   snow_persistence                    =>    cps%snow_persistence           & ! Input:  [real(r8) (:)]  counter for length of time snow-covered
   )
    
    ! Determine time step and step size

    nstep = get_nstep()
    dtime = get_step_size()
    
    if (use_vichydro) then
       call CLMVICMap(bounds, num_hydrologyc, filter_hydrologyc)
    endif
    call Drainage(bounds, num_hydrologyc, filter_hydrologyc, &
                  num_urbanc, filter_urbanc)


    do j = 1, nlevgrnd
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
               .or. ctype(c) == icol_roof) .and. j > nlevurb) then
          else
            h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
          end if
       end do
    end do


    do fc = 1, num_nolakec
       
       c = filter_nolakec(fc)
       l = clandunit(c)

       if (ctype(c) == icol_roof .or. ctype(c) == icol_sunwall &
          .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_road_imperv) then
          endwb(c) = h2ocan(c) + h2osno(c)
       else
          ! add h2osfc to water balance
          endwb(c) = h2ocan(c) + h2osno(c) + h2osfc(c) + wa(c)
       end if

    end do

    do j = 1, nlevgrnd
      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
            .or. ctype(c) == icol_roof) .and. j > nlevurb) then
               
         else
            endwb(c) = endwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
         end if
      end do
    end do


    ! Prior to summing up wetland/ice hydrology, calculate land ice contributions/sinks
    ! to this hydrology.
    ! 1) Generate SMB from capped-snow amount.  This is done over istice_mec
    !    columns, and also any other columns included in do_smb_c filter, where
    !    perennial snow has remained for at least snow_persistence_max.
    ! 2) If using glc_dyn_runoff_routing, zero qflx_snwcp_ice: qflx_snwcp_ice is the flux
    !    sent to ice runoff, but for glc_dyn_runoff_routing, we do NOT want this to be
    !    sent to ice runoff (instead it is sent to CISM).
    do c = bounds%begc,bounds%endc
       qflx_glcice_frz(c) = 0._r8
    end do 
    do fc = 1,num_do_smb_c
       c = filter_do_smb_c(fc)
       l = clandunit(c)
       ! In the following, we convert glc_snow_persistence_max_days to r8 to avoid overflow
       if ( (snow_persistence(c) >= (real(glc_snow_persistence_max_days, r8) * secspday)) &
             .or. lun%itype(l) == istice_mec) then
         qflx_glcice_frz(c) = qflx_snwcp_ice(c)  
         qflx_glcice(c) = qflx_glcice(c) + qflx_glcice_frz(c)
         if (glc_dyn_runoff_routing) qflx_snwcp_ice(c) = 0._r8
       end if     
    end do
    
    ! Determine wetland and land ice hydrology (must be placed here
    ! since need snow updated from CombineSnowLayers)
    
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       g = cgridcell(c)
       if (ityplun(l)==istwet .or. ityplun(l)==istice      &
                              .or. ityplun(l)==istice_mec) then
          qflx_drain(c)         = 0._r8
          qflx_drain_perched(c) = 0._r8
          qflx_h2osfc_surf(c)   = 0._r8
          qflx_surf(c)          = 0._r8
          qflx_infl(c)          = 0._r8
          qflx_qrgwl(c) = forc_rain(c) + forc_snow(c) + qflx_floodg(g) - qflx_evap_tot(c) - qflx_snwcp_ice(c) - &
               (endwb(c)-begwb(c))/dtime
               
          ! With glc_dyn_runoff_routing = false (the less realistic way, typically used
          ! when NOT coupling to CISM), excess snow immediately runs off, whereas melting
          ! ice stays in place and does not run off. The reverse is true with
          ! glc_dyn_runoff_routing = true: in this case, melting ice runs off, and excess
          ! snow is sent to CISM, where it is converted to ice. These corrections are
          ! done here: 
          if (glc_dyn_runoff_routing .and. ityplun(l)==istice_mec) then
             ! If using glc_dyn_runoff_routing, add meltwater from istice_mec ice columns to the runoff.
             !    Note: The meltwater contribution is computed in PhaseChanges (part of Biogeophysics2)
             qflx_qrgwl(c) = qflx_qrgwl(c) + qflx_glcice_melt(c)
             ! Also subtract the freezing component of qflx_glcice: this ice is added to
             ! CISM's ice column rather than running off. (This is analogous to the
             ! subtraction of qflx_snwcp_ice from qflx_qrgwl above, which accounts for
             ! snow that should be put into ice runoff rather than liquid runoff. But for
             ! glc_dyn_runoff_routing=true, qflx_snwcp_ice has been zeroed out, and has
             ! been put into qflx_glcice_frz.)
             qflx_qrgwl(c) = qflx_qrgwl(c) - qflx_glcice_frz(c)
          endif
          
          fcov(c)       = spval
          fsat(c)       = spval
          qcharge(c)    = spval
          qflx_rsub_sat(c) = spval
       else if (urbpoi(l) .and. ctype(c) /= icol_road_perv) then
          fcov(c)               = spval
          fsat(c)               = spval
          qflx_drain_perched(c) = 0._r8
          qflx_h2osfc_surf(c)   = 0._r8
          qcharge(c)            = spval
          qflx_rsub_sat(c)      = spval
       end if

       qflx_runoff(c) = qflx_drain(c) + qflx_surf(c)  + qflx_h2osfc_surf(c) + qflx_qrgwl(c) + qflx_drain_perched(c)

       if ((ityplun(l)==istsoil .or. ityplun(l)==istcrop) &
           .and. cactive(c)) then
          qflx_runoff(c) = qflx_runoff(c) - qflx_irrig(c)
       end if
       if (urbpoi(l)) then
         qflx_runoff_u(c) = qflx_runoff(c)
       else if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
         qflx_runoff_r(c) = qflx_runoff(c)
       end if

    end do
      
  end associate
    
  end subroutine Hydrology2Drainage

end module Hydrology2Mod
