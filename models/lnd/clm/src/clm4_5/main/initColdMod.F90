module initColdMod

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  ! Create cold start initial conditions  
  ! This is always called now - and will be overwritten by spun-up initial
  ! coniditions if either finidat or finidat_interp_source is set to
  ! non-blank fields
  !---------------------------------------------------------------------------
  !
  implicit none
  SAVE
  private                              ! By default make data private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: initCold   ! Make arbitrary initial conditions
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: initColdDefault
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initCold(bounds)
    !
    ! DESCRIPTION:
    ! Create cold start initial conditions 
    ! This is always called - and will be overwritten by spun-up initial
    ! coniditions if either finidat or finidat_interp_source is set to
    ! non-blank fields
    !
    ! USES:
    use clm_varctl   , only : use_cn, use_cndv, use_lch4, iulog
    use UrbanInitMod , only : initColdUrban
    use SLakeInitMod , only : initColdSLake
    use CNInitMod    , only : initColdCN
    use CNDVInitMod  , only : initColdCNDV
    use ch4InitMod   , only : initColdCH4
    use spmdMod      , only : masterproc
    use decompMod    , only : bounds_type
    !
    ! ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds   ! bounds
    !-----------------------------------------------------------------------
  
    if ( masterproc )then
        write(iulog,*) 'Setting initial data to non-spun up values for all clm model components'
    end if

    ! Create cold start initial data fields

    call initColdDefault(bounds)

    call initColdSLake(bounds)

    call initColdUrban(bounds)

    if (use_cn)   call initColdCN(bounds)

    if (use_cndv) call initColdCNDV(bounds)

    if (use_lch4) call initColdCH4(bounds)

  end subroutine initCold

  !-----------------------------------------------------------------------
  subroutine initColdDefault(bounds)
    !
    ! !DESCRIPTION:
    ! Initializes default cold start conditions - some of these might
    ! be overwritten in later calls to initColdXXX subroutines
    !
    ! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype      , only : pft, pps, pes, pef, pws  
    use clmtype      , only : col, cps, ces, cws, pws_a, cwf
    use clmtype      , only : lun, lps
    use clmtype      , only : grc
    use shr_const_mod, only : SHR_CONST_TKFRZ
    use clm_varpar   , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use clm_varcon   , only : bdsno, istice, istwet, istsoil, zlnd
    use clm_varcon   , only : denice, denh2o, spval, sb, icol_road_perv
    use clm_varcon   , only : icol_road_imperv, icol_roof, icol_sunwall
    use clm_varcon   , only : icol_shadewall, istcrop, istice_mec, h2osno_max
    use clm_varcon   , only : col_itype_to_icemec_class
    use clm_varctl   , only : iulog, use_vancouver, use_mexicocity
    use clm_varsur   , only : topo_glc_mec
    use spmdMod      , only : masterproc
    use decompMod    , only : bounds_type
    use domainMod    , only : ldomain
    use SNICARMod    , only : snw_rds_min
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds   ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,p,g    ! indices
    integer  :: icemec_class ! current icemec class (1..maxpatch_glcmec)
    integer  :: nlevs        ! number of levels
    real(r8) :: snowbd       ! temporary calculation of snow bulk density (kg/m3)
    real(r8) :: fmelt        ! snowbd/100
    !-----------------------------------------------------------------------

    ! cps%watsat              Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
    ! cps%sucsat              Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       

    ! lps%sabs_roof_dir       Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
    ! lps%sabs_roof_dif       Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
    ! lps%sabs_sunwall_dir    Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
    ! lps%sabs_sunwall_dif    Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
    ! lps%sabs_shadewall_dir  Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
    ! lps%sabs_shadewall_dif  Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
    ! lps%sabs_improad_dir    Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
    ! lps%sabs_improad_dif    Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
    ! lps%sabs_perroad_dir    Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
    ! lps%sabs_perroad_dif    Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux

    ! cps%z                   Output: [real(r8) (:,:) ]  layer depth  (m) OVER SNOW ONLY                 
    ! cps%zi                  Output: [real(r8) (:,:) ]  interface level below a "z" level (m) OVER SNOW ONLY           
    ! cps%dz                  Output: [real(r8) (:,:) ]  layer thickness depth (m) OVER SNOW ONLY           
    ! cps%frac_h2osfc         Output: [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
    ! cps%frac_sno            Output: [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
    ! cps%snl                 Output: [integer (:)    ]  number of snow layers                              
    ! cps%soilpsi             Output: [real(r8) (:,:) ]  soil water potential in each soil layer (MPa)   
    ! cps%snow_depth          Output: [real(r8) (:)   ]  snow height (m)                                   
    ! cps%snw_rds             Output: [real(r8) (:,:) ]  effective snow grain radius (col,lyr) [microns, m^-6]
    ! cps%snw_rds_top         Output: [real(r8) (:)   ]  snow grain size, top (col) [microns]              
    ! cps%sno_liq_top         Output: [real(r8) (:)   ]  liquid water fraction (mass) in top snow layer (col) [frc]
    ! cps%snow_persistence    Output: [real(r8) (:)   ]  perennial snow counter

    ! cps%mss_bcpho           Output: [real(r8) (:,:) ]  mass of hydrophobic BC in snow (col,lyr) [kg]   
    ! cps%mss_bcphi           Output: [real(r8) (:,:) ]  mass of hydrophillic BC in snow (col,lyr) [kg]  
    ! cps%mss_bctot           Output: [real(r8) (:,:) ]  total mass of BC (pho+phi) (col,lyr) [kg]       
    ! cps%mss_cnc_bcphi       Output: [real(r8) (:,:) ]  mass concentration of BC species 1 (col,lyr) [kg/kg]
    ! cps%mss_cnc_bcpho       Output: [real(r8) (:,:) ]  mass concentration of BC species 2 (col,lyr) [kg/kg]
    ! cps%mss_ocpho           Output: [real(r8) (:,:) ]  mass of hydrophobic OC in snow (col,lyr) [kg]   
    ! cps%mss_ocphi           Output: [real(r8) (:,:) ]  mass of hydrophillic OC in snow (col,lyr) [kg]  
    ! cps%mss_octot           Output: [real(r8) (:,:) ]  total mass of OC (pho+phi) (col,lyr) [kg]       
    ! cps%mss_oc_col          Output: [real(r8) (:)   ]  total mass of OC in snow column (col) [kg]        
    ! cps%mss_oc_top          Output: [real(r8) (:)   ]  total mass of OC in top snow layer (col) [kg]     
    ! cps%mss_cnc_ocphi       Output: [real(r8) (:,:) ]  mass concentration of OC species 1 (col,lyr) [kg/kg]
    ! cps%mss_cnc_ocpho       Output: [real(r8) (:,:) ]  mass concentration of OC species 2 (col,lyr) [kg/kg]
    ! cps%mss_dst1            Output: [real(r8) (:,:) ]  mass of dust species 1 in snow (col,lyr) [kg]   
    ! cps%mss_dst2            Output: [real(r8) (:,:) ]  mass of dust species 2 in snow (col,lyr) [kg]   
    ! cps%mss_dst3            Output: [real(r8) (:,:) ]  mass of dust species 3 in snow (col,lyr) [kg]   
    ! cps%mss_dst4            Output: [real(r8) (:,:) ]  mass of dust species 4 in snow (col,lyr) [kg]   
    ! cps%mss_dsttot          Output: [real(r8) (:,:) ]  total mass of dust in snow (col,lyr) [kg]       
    ! cps%mss_dst_col         Output: [real(r8) (:)   ]  total mass of dust in snow column (col) [kg]      
    ! cps%mss_dst_top         Output: [real(r8) (:)   ]  total mass of dust in top snow layer (col) [kg]   
    ! cps%mss_cnc_dst1        Output: [real(r8) (:,:) ]  mass concentration of dust species 1 (col,lyr) [kg/kg]
    ! cps%mss_cnc_dst2        Output: [real(r8) (:,:) ]  mass concentration of dust species 2 (col,lyr) [kg/kg]
    ! cps%mss_cnc_dst3        Output: [real(r8) (:,:) ]  mass concentration of dust species 3 (col,lyr) [kg/kg]
    ! cps%mss_cnc_dst4        Output: [real(r8) (:,:) ]  mass concentration of dust species 4 (col,lyr) [kg/kg]

    ! cps%albgrd              Output: [real(r8) (:,:) ]  ground albedo (direct)                
    ! cps%albgri              Output: [real(r8) (:,:) ]  ground albedo (diffuse)               
    ! cps%albsod              Output: [real(r8) (:,:) ]  direct-beam soil albedo (col,bnd) [frc]
    ! cps%albsoi              Output: [real(r8) (:,:) ]  diffuse soil albedo (col,bnd) [frc]   
    ! cps%albsnd_hst          Output: [real(r8) (:,:) ]  snow albedo, direct, for history files (col,bnd) [frc]
    ! cps%albsni_hst          Output: [real(r8) (:,:) ]  snow ground albedo, diffuse, for history files (col,bnd) [frc]
    ! cps%albgrd_pur(:,:)     Output: [real(r8) (:,:) ]  pure snow ground direct albedo (numrad)
    ! cps%albgri_pur(:,:)     Output: [real(r8) (:,:) ]  pure snow ground diffuse albedo (numrad)
    ! cps%albgrd_bc(:,:)      Output: [real(r8) (:,:) ]  ground direct albedo without BC  (numrad)
    ! cps%albgri_bc(:,:)      Output: [real(r8) (:,:) ]  ground diffuse albedo without BC (numrad)
    ! cps%albgrd_oc(:,:)      Output: [real(r8) (:,:) ]  ground direct albedo without OC  (numrad)
    ! cps%albgri_oc(:,:)      Output: [real(r8) (:,:) ]  ground diffuse albedo without OC (numrad)
    ! cps%albgrd_dst(:,:)     Output: [real(r8) (:,:) ]  ground direct albedo without dust  (numrad)
    ! cps%albgri_dst(:,:)     Output: [real(r8) (:,:) ]  ground diffuse albedo without dust (numrad)

    ! cws%h2osfc              Output: [real(r8) (:)   ]  surface water (mm)                                
    ! cws%zwt_perched         Output: [real(r8) (:)   ]  perched water table depth (m)                     
    ! cws%int_snow            Output: [real(r8) (:)   ]  integrated snowfall                               
    ! cws%h2osoi_ice          Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
    ! cws%h2osoi_liq          Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
    ! cws%h2osoi_vol          Output: [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    ! cws%h2osno              Output: [real(r8) (:)   ]  snow water (mm H2O)                               
    ! cws%wa                  Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
    ! cws%zwt                 Output: [real(r8) (:)   ]  water table depth (m)                             
    ! cws%fsat                Output: [real(r8) (:)   ] fractional area with water table at surface        
    ! cws%frost_table         Output: [real(r8) (:)   ]  frost table depth (m)                             
    ! cwf%qflx_h2osfc_surf    Output: [real(r8) (:)   ] surface water runoff (mm/s)                        
    ! cwf%qflx_snow_melt      Output: [real(r8) (:)   ]  snow melt (net)                                   

    ! ces%t_h2osfc            Output: [real(r8) (:)   ]  surface water temperature                         
    ! ces%t_soisno            Output: [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    ! ces%t_lake              Output: [real(r8) (:,:) ]  lake temperature (Kelvin)  (1:nlevlak)          
    ! ces%t_grnd              Output: [real(r8) (:)   ]  ground temperature (Kelvin)                       
    ! ces%tsoi17              Output: [real(r8) (:)   ]  soil T for top 0.17 m                             

    ! pps%n_irrig_steps_left  Output: [integer (:)    ]  number of time steps for which we still need to irrigate today (if 0, ignore irrig_rate)
    ! pps%irrig_rate          Output: [real(r8) (:)   ]  current irrigation rate [mm/s]                    
    ! pps%albd                Output: [real(r8) (:,:) ]  surface albedo (direct)               
    ! pps%albi                Output: [real(r8) (:,:) ]  surface albedo (diffuse)              
    ! pps%fabd                Output: [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
    ! pps%fabi                Output: [real(r8) (:,:) ]  flux absorbed by canopy per unit diffuse flux
    ! pps%fabd_sun            Output: [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit direct flux
    ! pps%fabd_sha            Output: [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit direct flux
    ! pps%fabi_sun            Output: [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit diffuse flux
    ! pps%fabi_sha            Output: [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit diffuse flux
    ! pps%ftdd                Output: [real(r8) (:,:) ]  down direct flux below canopy per unit direct flux
    ! pps%ftid                Output: [real(r8) (:,:) ]  down diffuse flux below canopy per unit direct flux
    ! pps%ftii                Output: [real(r8) (:,:) ]  down diffuse flux below canopy per unit diffuse flu

    ! pws%h2ocan              Output: [real(r8) (:)   ]  canopy water (mm H2O) (pft-level)                 
    ! pws_a%h2ocan            Output: [real(r8) (:)   ]  canopy water (mm H2O) (column-level)              

    ! pes%t_veg               Output: [real(r8) (:)   ]  vegetation temperature (Kelvin)                   
    ! pes%t_ref2m             Output: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)       
    ! pes%t_ref2m_u           Output: [real(r8) (:)   ]  Urban 2 m height surface air temperature (Kelvin) 
    ! pes%t_ref2m_r           Output: [real(r8) (:)   ]  Rural 2 m height surface air temperature (Kelvin) 
    ! pef%eflx_lwrad_out      Output: [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)    

    if ( masterproc )then
        write(iulog,*) 'Setting initial data to non-spun up values for non-lake points '
    end if

    associate(snl => cps%snl) ! Output: [integer (:)    ]  number of snow layers   

    ! -----------------------------------------------------------------
    ! initialize grid cell-level quantities
    ! -----------------------------------------------------------------

    do g = bounds%begg, bounds%endg
       ! glcmask (from a file) provides a rough guess of the icemask (from CISM); thus, in
       ! initialization, set icemask equal to glcmask; icemask will later get updated at
       ! the start of the run loop, as soon as we have data from CISM
       grc%icemask(g) = ldomain%glcmask(g)
    end do
    
    ! -----------------------------------------------------------------
    ! initialize glc_topo
    ! -----------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       g = col%gridcell(c)
       if (lun%itype(l) == istice_mec) then
          ! For ice_mec landunits, initialize glc_topo based on surface dataset; this
          ! will get overwritten in the run loop by values sent from CISM
          icemec_class = col_itype_to_icemec_class(col%itype(c))
          cps%glc_topo(c) = topo_glc_mec(g, icemec_class)
       else
          ! For other landunits, arbitrarily initialize glc_topo to 0 m; for landunits
          ! where this matters, this will get overwritten in the run loop by values sent
          ! from CISM
          cps%glc_topo(c) = 0._r8
       end if
    end do

    ! -----------------------------------------------------------------
    ! initialize h2osfc, frac_h2osfc, t_h2osfc, qflx_snow_melt
    ! -----------------------------------------------------------------

    do c = bounds%begc,bounds%endc
       cws%h2osfc(c)           = 0._r8
       cps%frac_h2osfc(c)      = 0._r8
       cwf%qflx_h2osfc_surf(c) = 0._r8
       cwf%qflx_snow_melt(c)   = 0._r8
       ces%t_h2osfc(c)         = 274._r8
    enddo

    ! -----------------------------------------------------------------
    ! Initialize canopy water, snow water and snow depth
    ! -----------------------------------------------------------------

    ! NOTE: h2ocan, h2osno, and snow_depth has valid values everywhere

    ! canopy water (pft level and column)
    do p = bounds%begp, bounds%endp
       pws%h2ocan(p) = 0._r8
    end do
    do c = bounds%begc,bounds%endc
       pws_a%h2ocan(c) = 0._r8
    end do

    do c = bounds%begc,bounds%endc
       ! snow water
       ! Note: Glacier_mec columns are initialized with half the maximum snow cover.
       ! This gives more realistic values of qflx_glcice sooner in the simulation
       ! for columns with net ablation, at the cost of delaying ice formation
       ! in columns with net accumulation.
       l = col%landunit(c)
       g = col%gridcell(c)

       if (lun%itype(l)==istice) then
          cws%h2osno(c) = h2osno_max
       elseif (lun%itype(l)==istice_mec .or. &
              (lun%itype(l)==istsoil .and. ldomain%glcmask(g) > 0._r8)) then 
          ! Initialize a non-zero snow thickness where the ice sheet can/potentially operate.
          ! Using glcmask to capture all potential vegetated points around GrIS (ideally
          ! we would use icemask from CISM, but that isn't available until after initialization.)
          cws%h2osno(c) = 0.5_r8 * h2osno_max  
       else
          cws%h2osno(c) = 0._r8
       endif

       ! integrated snowfall
       cws%int_snow(c) = cws%h2osno(c)

       cps%snow_depth(c) = cws%h2osno(c) / bdsno
       cps%snow_persistence(c) = 0._r8
    end do

    ! -----------------------------------------------------------------
    ! initialize snow levels and interfaces (lake and non-lake points)
    ! Note that cps%zi(0) is set in routine iniTimeConst.
    ! -----------------------------------------------------------------

     do c = bounds%begc,bounds%endc
        cps%dz(c,-nlevsno+1: 0) = 1.e36_r8
        cps%z (c,-nlevsno+1: 0) = 1.e36_r8
        cps%zi(c,-nlevsno  :-1) = 1.e36_r8

        l = col%landunit(c)
        if (.not. lun%lakpoi(l)) then
           if (cps%snow_depth(c) < 0.01_r8) then
              snl(c)             = 0
              cps%dz(c,-nlevsno+1:0) = 0._r8
              cps%z (c,-nlevsno+1:0) = 0._r8
              cps%zi(c,-nlevsno+0:0) = 0._r8
           else
              if ((cps%snow_depth(c) >= 0.01_r8) .and. (cps%snow_depth(c) <= 0.03_r8)) then
                 snl(c)  = -1
                 cps%dz(c,0) = cps%snow_depth(c)
              else if ((cps%snow_depth(c) > 0.03_r8) .and. (cps%snow_depth(c) <= 0.04_r8)) then
                 snl(c)   = -2
                 cps%dz(c,-1) = cps%snow_depth(c)/2._r8
                 cps%dz(c, 0) = cps%dz(c,-1)
              else if ((cps%snow_depth(c) > 0.04_r8) .and. (cps%snow_depth(c) <= 0.07_r8)) then
                 snl(c)   = -2
                 cps%dz(c,-1) = 0.02_r8
                 cps%dz(c, 0) = cps%snow_depth(c) - cps%dz(c,-1)
              else if ((cps%snow_depth(c) > 0.07_r8) .and. (cps%snow_depth(c) <= 0.12_r8)) then
                 snl(c)   = -3
                 cps%dz(c,-2) = 0.02_r8
                 cps%dz(c,-1) = (cps%snow_depth(c) - 0.02_r8)/2._r8
                 cps%dz(c, 0) = cps%dz(c,-1)
              else if ((cps%snow_depth(c) > 0.12_r8) .and. (cps%snow_depth(c) <= 0.18_r8)) then
                 snl(c)   = -3
                 cps%dz(c,-2) = 0.02_r8
                 cps%dz(c,-1) = 0.05_r8
                 cps%dz(c, 0) = cps%snow_depth(c) - cps%dz(c,-2) - cps%dz(c,-1)
              else if ((cps%snow_depth(c) > 0.18_r8) .and. (cps%snow_depth(c) <= 0.29_r8)) then
                 snl(c)   = -4
                 cps%dz(c,-3) = 0.02_r8
                 cps%dz(c,-2) = 0.05_r8
                 cps%dz(c,-1) = (cps%snow_depth(c) - cps%dz(c,-3) - cps%dz(c,-2))/2._r8
                 cps%dz(c, 0) = cps%dz(c,-1)
              else if ((cps%snow_depth(c) > 0.29_r8) .and. (cps%snow_depth(c) <= 0.41_r8)) then
                 snl(c)   = -4
                 cps%dz(c,-3) = 0.02_r8
                 cps%dz(c,-2) = 0.05_r8
                 cps%dz(c,-1) = 0.11_r8
                 cps%dz(c, 0) = cps%snow_depth(c) - cps%dz(c,-3) - cps%dz(c,-2) - cps%dz(c,-1)
              else if ((cps%snow_depth(c) > 0.41_r8) .and. (cps%snow_depth(c) <= 0.64_r8)) then
                 snl(c)   = -5
                 cps%dz(c,-4) = 0.02_r8
                 cps%dz(c,-3) = 0.05_r8
                 cps%dz(c,-2) = 0.11_r8
                 cps%dz(c,-1) = (cps%snow_depth(c) - cps%dz(c,-4) - cps%dz(c,-3) - cps%dz(c,-2))/2._r8
                 cps%dz(c, 0) = cps%dz(c,-1)
              else if (cps%snow_depth(c) > 0.64_r8) then
                 snl(c)   = -5
                 cps%dz(c,-4) = 0.02_r8
                 cps%dz(c,-3) = 0.05_r8
                 cps%dz(c,-2) = 0.11_r8
                 cps%dz(c,-1) = 0.23_r8
                 cps%dz(c, 0) = cps%snow_depth(c)-cps%dz(c,-4)-cps%dz(c,-3)-cps%dz(c,-2)-cps%dz(c,-1)
              endif
           end if
           do j = 0, snl(c)+1, -1
              cps%z(c,j)    = cps%zi(c,j) - 0.5_r8*cps%dz(c,j)
              cps%zi(c,j-1) = cps%zi(c,j) - cps%dz(c,j)
           end do
        else !lake
           snl(c)             = 0
           cps%dz(c,-nlevsno+1:0) = 0._r8
           cps%z (c,-nlevsno+1:0) = 0._r8
           cps%zi(c,-nlevsno+0:0) = 0._r8
        end if
     end do

    ! -----------------------------------------------------------------
    ! Set snow/soil temperature, note:
    ! t_soisno only has valid values over non-lake
    ! t_lake   only has valid values over lake
    ! t_grnd has valid values over all land
    ! t_veg  has valid values over all land
    ! -----------------------------------------------------------------

    do c = bounds%begc,bounds%endc

       ces%t_soisno(c,-nlevsno+1:nlevgrnd) = spval
       ces%t_lake(c,1:nlevlak) = spval

       l = col%landunit(c)
       if (.not. lun%lakpoi(l)) then  !not lake
          ces%t_soisno(c,-nlevsno+1:0) = spval
          if (snl(c) < 0) then    !snow layer temperatures
             do j = snl(c)+1, 0
                ces%t_soisno(c,j) = 250._r8
             enddo
          endif
          if (lun%itype(l)==istice .or. lun%itype(l)==istice_mec) then
             do j = 1, nlevgrnd
                ces%t_soisno(c,j) = 250._r8
             end do
          else if (lun%itype(l) == istwet) then
             do j = 1, nlevgrnd
                ces%t_soisno(c,j) = 277._r8
             end do
          else if (lun%urbpoi(l)) then
             if (use_vancouver) then
                if (col%itype(c) == icol_road_perv .or. col%itype(c) == icol_road_imperv) then 
                   ! Set road top layer to initial air temperature and interpolate other
                   ! layers down to 20C in bottom layer
                   do j = 1, nlevgrnd
                      ces%t_soisno(c,j) = 297.56 - (j-1) * ((297.56-293.16)/(nlevgrnd-1)) 
                   end do
                   ! Set wall and roof layers to initial air temperature
                else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall .or. col%itype(c) == icol_roof) then
                   ces%t_soisno(c,1:nlevurb) = 297.56
                else
                   ces%t_soisno(c,1:nlevgrnd) = 283._r8
                end if
             else if (use_mexicocity) then
                if (col%itype(c) == icol_road_perv .or. col%itype(c) == icol_road_imperv) then 
                   ! Set road top layer to initial air temperature and interpolate other
                   ! layers down to 22C in bottom layer
                   do j = 1, nlevgrnd
                      ces%t_soisno(c,j) = 289.46 - (j-1) * ((289.46-295.16)/(nlevgrnd-1)) 
                   end do
                   ! Set wall and roof layers to initial air temperature
                else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall .or. col%itype(c) == icol_roof) then
                   ces%t_soisno(c,1:nlevurb) = 289.46
                else
                   ces%t_soisno(c,1:nlevgrnd) = 283._r8
                end if
             else
                if (col%itype(c) == icol_road_perv .or. col%itype(c) == icol_road_imperv) then 
                   ces%t_soisno(c,1:nlevgrnd) = 274._r8
                   ! Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
                   ! shock from large heating/air conditioning flux
                else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
                     .or. col%itype(c) == icol_roof) then
                   ces%t_soisno(c,1:nlevurb) = 292._r8
                end if
             end if
          else
             ces%t_soisno(c,1:nlevgrnd) = 274._r8
          endif
          ces%t_grnd(c) = ces%t_soisno(c,snl(c)+1)
       else                     !lake
          ces%t_lake(c,1:nlevlak) = 277._r8
          ces%t_grnd(c) = ces%t_lake(c,1)
       endif
       ces%tsoi17(c) = ces%t_grnd(c)

    end do

    ! -----------------------------------------------------------------
    ! Initialize pft level irrigation variables
    ! -----------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       c = pft%column(p)
       l = pft%landunit(p)

       ! Initialize Irrigation to zero
       if (lun%itype(l)==istsoil) then
          pps%n_irrig_steps_left(p) = 0
          pps%irrig_rate(p)         = 0.0_r8
       end if

       if (use_vancouver) then
          pes%t_veg(p) = 297.56
          pes%t_ref2m(p) = 297.56
          if (lun%urbpoi(l)) then
             pes%t_ref2m_u(p) = 297.56
          else
             pes%t_ref2m_u(p) = spval
          end if
          if (lun%ifspecial(l)) then
             pes%t_ref2m_r(p) = spval
          else
             pes%t_ref2m_r(p) = 297.56
          end if
       else if (use_mexicocity) then
          pes%t_veg(p) = 289.46
          pes%t_ref2m(p) = 289.46
          if (lun%urbpoi(l)) then
             pes%t_ref2m_u(p) = 289.46
          else
             pes%t_ref2m_u(p) = spval
          end if
          if (lun%ifspecial(l)) then
             pes%t_ref2m_r(p) = spval
          else
             pes%t_ref2m_r(p) = 289.46
          end if
       else
          pes%t_veg(p) = 283._r8
          pes%t_ref2m(p) = 283._r8
          if (lun%urbpoi(l)) then
             pes%t_ref2m_u(p) = 283._r8
          else
             pes%t_ref2m_u(p) = spval
          end if
          if (lun%ifspecial(l)) then
             pes%t_ref2m_r(p) = spval
          else
             pes%t_ref2m_r(p) = 283._r8
          end if
       end if
       pef%eflx_lwrad_out(p) = sb * (ces%t_grnd(c))**4
    end do

    ! -----------------------------------------------------------------
    ! Initialize surface albedos and absorbed fluxes to reasonable values
    ! -----------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%urbpoi(l)) then
          ! From Bonan 1996 (LSM technical note)
          cps%frac_sno(c) = min( cps%snow_depth(c)/0.05_r8, 1._r8)
       else
          cps%frac_sno(c) = 0._r8
          ! snow cover fraction as in Niu and Yang 2007
          if(cps%snow_depth(c) .gt. 0.0)  then
             snowbd   = min(400._r8, cws%h2osno(c)/cps%snow_depth(c)) !bulk density of snow (kg/m3)
             fmelt    = (snowbd/100.)**1.
             ! 100 is the assumed fresh snow density; 1 is a melting factor that could be
             ! reconsidered, optimal value of 1.5 in Niu et al., 2007
             cps%frac_sno(c) = tanh( cps%snow_depth(c) /(2.5 * zlnd * fmelt) )
          endif
       end if
    end do

    cps%albgrd(bounds%begc:bounds%endc, :)     = 0.2_r8
    cps%albgri(bounds%begc:bounds%endc, :)     = 0.2_r8
    cps%albsod(bounds%begc:bounds%endc, :)     = 0.2_r8
    cps%albsoi(bounds%begc:bounds%endc, :)     = 0.2_r8
    cps%albsnd_hst(bounds%begc:bounds%endc, :) = 0.6_r8
    cps%albsni_hst(bounds%begc:bounds%endc, :) = 0.6_r8
    pps%albd(bounds%begp:bounds%endp, :)       = 0.2_r8
    pps%albi(bounds%begp:bounds%endp, :)       = 0.2_r8

    cps%albgrd_pur(bounds%begc:bounds%endc, :) = 0.2_r8
    cps%albgri_pur(bounds%begc:bounds%endc, :) = 0.2_r8
    cps%albgrd_bc(bounds%begc:bounds%endc, :)  = 0.2_r8
    cps%albgri_bc(bounds%begc:bounds%endc, :)  = 0.2_r8
    cps%albgrd_oc(bounds%begc:bounds%endc, :)  = 0.2_r8
    cps%albgri_oc(bounds%begc:bounds%endc, :)  = 0.2_r8
    cps%albgrd_dst(bounds%begc:bounds%endc, :) = 0.2_r8
    cps%albgri_dst(bounds%begc:bounds%endc, :) = 0.2_r8

    pps%fabi(bounds%begp:bounds%endp, :)       = 0.0_r8
    pps%fabd(bounds%begp:bounds%endp, :)       = 0.0_r8
    pps%fabi_sun(bounds%begp:bounds%endp, :)   = 0.0_r8
    pps%fabd_sun(bounds%begp:bounds%endp, :)   = 0.0_r8
    pps%fabd_sha(bounds%begp:bounds%endp, :)   = 0.0_r8
    pps%fabi_sha(bounds%begp:bounds%endp, :)   = 0.0_r8
    pps%ftdd(bounds%begp:bounds%endp, :)       = 1.0_r8
    pps%ftid(bounds%begp:bounds%endp, :)       = 0.0_r8
    pps%ftii(bounds%begp:bounds%endp, :)       = 1.0_r8

    ! -----------------------------------------------------------------
    ! Initialize frost table
    ! -----------------------------------------------------------------

    cws%wa(bounds%begc:bounds%endc)  = 5000._r8
    cws%zwt(bounds%begc:bounds%endc) = 0._r8

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       if (.not. lun%lakpoi(l)) then  !not lake
          if (lun%urbpoi(l)) then
             if (col%itype(c) == icol_road_perv) then
                cws%wa(c)  = 4800._r8
                cws%zwt(c) = (25._r8 + cps%zi(c,nlevsoi)) - cws%wa(c)/0.2_r8 /1000._r8  ! One meter below soil column
             else
                cws%wa(c)  = spval
                cws%zwt(c) = spval
             end if
             ! initialize frost_table, zwt_perched
             cws%zwt_perched(c) = spval
             cws%frost_table(c) = spval
          else
             cws%wa(c)  = 4000._r8
             cws%zwt(c) = (25._r8 + cps%zi(c,nlevsoi)) - cws%wa(c)/0.2_r8 /1000._r8  ! One meter below soil column
             ! initialize frost_table, zwt_perched to bottom of soil column
             cws%zwt_perched(c) = cps%zi(c,nlevsoi)
             cws%frost_table(c) = cps%zi(c,nlevsoi)
          end if
       end if
    end do

    ! -----------------------------------------------------------------
    ! Set soil water
    ! -----------------------------------------------------------------

    ! volumetric water is set first and liquid content and ice lens are obtained
    ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
    ! and urban pervious road (other urban columns have zero soil water)

    cws%h2osoi_vol(bounds%begc:bounds%endc,         1:) = spval
    cws%h2osoi_liq(bounds%begc:bounds%endc,-nlevsno+1:) = spval
    cws%h2osoi_ice(bounds%begc:bounds%endc,-nlevsno+1:) = spval

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       if (.not. lun%lakpoi(l)) then  !not lake

          ! volumetric water
          if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (j > nlevsoi) then
                   cws%h2osoi_vol(c,j) = 0.0_r8
                else
                   cws%h2osoi_vol(c,j) = 0.15_r8
                endif
             end do
          else if (lun%urbpoi(l)) then
             if (col%itype(c) == icol_road_perv) then
               nlevs = nlevgrnd
               do j = 1, nlevs
                  if (j <= nlevsoi) then
                     cws%h2osoi_vol(c,j) = 0.3_r8
                  else
                     cws%h2osoi_vol(c,j) = 0.0_r8
                  end if
               end do
             else if (col%itype(c) == icol_road_imperv) then
               nlevs = nlevgrnd
               do j = 1, nlevs
                  cws%h2osoi_vol(c,j) = 0.0_r8
               end do
             else
               nlevs = nlevurb
               do j = 1, nlevs
                  cws%h2osoi_vol(c,j) = 0.0_r8
               end do
             end if
          else if (lun%itype(l) == istwet) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (j > nlevsoi) then
                   cws%h2osoi_vol(c,j) = 0.0_r8
                else
                   cws%h2osoi_vol(c,j) = 1.0_r8
                endif
             end do
          else if (lun%itype(l) == istice .or. lun%itype(l) == istice_mec) then
             nlevs = nlevgrnd 
             do j = 1, nlevs
                cws%h2osoi_vol(c,j) = 1.0_r8
             end do
          endif
          do j = 1, nlevs
             cws%h2osoi_vol(c,j) = min(cws%h2osoi_vol(c,j),cps%watsat(c,j))
        
             ! soil layers
             if (ces%t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                cws%h2osoi_ice(c,j) = cps%dz(c,j)*denice*cws%h2osoi_vol(c,j)
                cws%h2osoi_liq(c,j) = 0._r8
             else
                cws%h2osoi_ice(c,j) = 0._r8
                cws%h2osoi_liq(c,j) = cps%dz(c,j)*denh2o*cws%h2osoi_vol(c,j)
             endif
          end do

       end if
    end do

    do j = -nlevsno+1, 0
       do c = bounds%begc,bounds%endc
          l = col%landunit(c)
          if (.not. lun%lakpoi(l)) then  !not lake
             if (j > snl(c)) then
                cws%h2osoi_ice(c,j) = cps%dz(c,j)*250._r8
                cws%h2osoi_liq(c,j) = 0._r8
             end if
          end if
       end do
    end do

    ! -----------------------------------------------------------------
    ! initialize SNICAR fields:
    ! -----------------------------------------------------------------

    do c = bounds%begc,bounds%endc
       cps%mss_bctot(c,:)     = 0._r8
       cps%mss_bcpho(c,:)     = 0._r8
       cps%mss_bcphi(c,:)     = 0._r8
       cps%mss_cnc_bcphi(c,:) = 0._r8
       cps%mss_cnc_bcpho(c,:) = 0._r8

       cps%mss_octot(c,:)     = 0._r8
       cps%mss_ocpho(c,:)     = 0._r8
       cps%mss_ocphi(c,:)     = 0._r8
       cps%mss_cnc_ocphi(c,:) = 0._r8
       cps%mss_cnc_ocpho(c,:) = 0._r8
       
       cps%mss_dst1(c,:)      = 0._r8
       cps%mss_dst2(c,:)      = 0._r8
       cps%mss_dst3(c,:)      = 0._r8
       cps%mss_dst4(c,:)      = 0._r8
       cps%mss_dsttot(c,:)    = 0._r8
       cps%mss_cnc_dst1(c,:)  = 0._r8
       cps%mss_cnc_dst2(c,:)  = 0._r8
       cps%mss_cnc_dst3(c,:)  = 0._r8
       cps%mss_cnc_dst4(c,:)  = 0._r8
       
       if (snl(c) < 0) then
          cps%snw_rds(c,snl(c)+1:0) = snw_rds_min
          cps%snw_rds(c,-nlevsno+1:snl(c)) = 0._r8
          cps%snw_rds_top(c) = snw_rds_min
          cps%sno_liq_top(c) = cws%h2osoi_liq(c,snl(c)+1) / &
               (cws%h2osoi_liq(c,snl(c)+1)+cws%h2osoi_ice(c,snl(c)+1))
       elseif (cws%h2osno(c) > 0._r8) then
          cps%snw_rds(c,0)             = snw_rds_min
          cps%snw_rds(c,-nlevsno+1:-1) = 0._r8
          cps%snw_rds_top(c)           = spval
          cps%sno_liq_top(c)           = spval
       else
          cps%snw_rds(c,:)   = 0._r8
          cps%snw_rds_top(c) = spval
          cps%sno_liq_top(c) = spval
       endif
    enddo

  end associate
  end subroutine initColdDefault

end module initColdMod
