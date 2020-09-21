module FireMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for fire dynamics 
  ! created in Nov, 2012  and revised in Apr, 2013 by F. Li and S. Levis
  ! based on Li et al. (2012a,b; 2013)
  ! revised in Apr, 2014 according Li et al.(2014)
  ! Fire-related parameters were calibrated or tuned in Apr, 2013 based on the 
  ! 20th Century transient simulations at f19_g16 with (newfire05_clm45sci15_clm4_0_58) 
  ! a CLM4.5 version, Qian et al. (2006) atmospheric forcing, and
  ! climatological lightning data.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_const_mod          , only : SHR_CONST_PI,SHR_CONST_TKFRZ
  use shr_infnan_mod         , only : shr_infnan_isnan
  use shr_strdata_mod        , only : shr_strdata_type, shr_strdata_create, shr_strdata_print
  use shr_strdata_mod        , only : shr_strdata_advance
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use elm_varctl             , only : iulog
  use elm_varpar             , only : nlevdecomp, ndecomp_pools
  use elm_varcon             , only : dzsoi_decomp
  use pftvarcon              , only : fsr_pft, fd_pft, noveg
  use spmdMod                , only : masterproc, mpicom, comp_id
  use fileutils              , only : getavu, relavu
  use controlMod             , only : NLFilename
  use decompMod              , only : gsmap_lnd_gdc2glo
  use domainMod              , only : ldomain
  use abortutils             , only : endrun
  use decompMod              , only : bounds_type
  use subgridAveMod          , only : p2c
  use CNDecompCascadeConType , only : decomp_cascade_con
  use VegetationPropertiesType         , only : veg_vp
  use atm2lndType            , only : atm2lnd_type
  use CNStateType            , only : cnstate_type
  use CNCarbonFluxType       , only : carbonflux_type
  use CNCarbonStateType      , only : carbonstate_type
  use CNNitrogenFluxType     , only : nitrogenflux_type
  use CNNitrogenStateType    , only : nitrogenstate_type
  use EnergyFluxType         , only : energyflux_type
  use SoilHydrologyType      , only : soilhydrology_type  
  use TemperatureType        , only : temperature_type
  use WaterstateType         , only : waterstate_type
  use GridcellType           , only : grc_pp
  use TopounitDataType       , only : top_as, top_af ! atmospheric state and flux variables  
  use ColumnType             , only : col_pp
  use ColumnDataType         , only : col_es, col_ws, col_cs, col_cf, col_ns, col_nf 
  use ColumnDataType         , only : col_ps, col_pf 
  use VegetationType         , only : veg_pp                
  use VegetationDataType     , only : veg_cs, veg_cf, veg_ns, veg_nf  
  use VegetationDataType     , only : veg_ps, veg_pf  
  use mct_mod
  use PhosphorusFluxType     , only : phosphorusflux_type
  use PhosphorusStateType    , only : phosphorusstate_type
  use SharedParamsMod      , only : ParamsShareInst
  use elm_varctl             , only : nu_com
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: FireInit    ! Initialization of CNFire
  public :: FireInterp  ! Interpolate fire data
  public :: FireArea    ! Calculate fire area
  public :: FireFluxes  ! Calculate fire fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: hdm_init    ! position datasets for dynamic human population density
  private :: hdm_interp  ! interpolates between two years of human pop. density file data
  private :: lnfm_init   ! position datasets for Lightning
  private :: lnfm_interp ! interpolates between two years of Lightning file data

  ! !PRIVATE MEMBER DATA:
  real(r8), pointer     :: forc_lnfm(:)        ! Lightning frequency
  real(r8), pointer     :: forc_hdm(:)         ! Human population density
  real(r8), parameter   :: secsphr = 3600._r8  ! Seconds in an hour
  real(r8), parameter   :: borealat = 40._r8   ! Latitude for boreal peat fires

  type(shr_strdata_type) :: sdat_hdm    ! Human population density input data stream
  type(shr_strdata_type) :: sdat_lnfm   ! Lightning input data stream
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine FireInit( bounds )
    !
    ! !DESCRIPTION:
    ! Initialize CN Fire module
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !-----------------------------------------------------------------------
#ifndef CPL_BYPASS
    call hdm_init(bounds)
    call hdm_interp(bounds)
    call lnfm_init(bounds)
    call lnfm_interp(bounds)
#endif
  end subroutine FireInit

  !-----------------------------------------------------------------------
  subroutine FireInterp(bounds)
    !
    ! !DESCRIPTION:
    ! Interpolate CN Fire datasets
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !-----------------------------------------------------------------------
#ifndef CPL_BYPASS
    call hdm_interp(bounds)
    call lnfm_interp(bounds)
#endif
  end subroutine FireInterp

  !-----------------------------------------------------------------------
  subroutine FireArea (bounds, &
       num_soilc, filter_soilc, &
       num_soilp, filter_soilp, &
       atm2lnd_vars, temperature_vars, energyflux_vars, soilhydrology_vars, waterstate_vars, &
       cnstate_vars, carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! Computes column-level burned area 
    !
    ! !USES:
    use clm_time_manager     , only: get_step_size, get_days_per_year, get_curr_date, get_nstep
    use elm_varpar           , only: max_patch_per_col
    use elm_varcon           , only: secspday
    use elm_varctl           , only: use_nofire, spinup_state, spinup_mortality_factor
    use dynSubgridControlMod , only: run_has_transient_landcover
    use pftvarcon            , only: nc4_grass, nc3crop, ndllf_evr_tmp_tree
    use pftvarcon            , only: nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree, nbrdlf_evr_shrub
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds 
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(energyflux_type)    , intent(in)    :: energyflux_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter  :: lfuel=75._r8    ! lower threshold of fuel mass (gC/m2) for ignition, Li et al.(2014)
    real(r8), parameter  :: ufuel=1050._r8  ! upper threshold of fuel mass(gC/m2) for ignition 
    real(r8), parameter  :: g0=0.05_r8      ! g(W) when W=0 m/s
    !
    ! a1 parameter for cropland fire in (Li et. al., 2014), but changed from
    ! /timestep to /hr
    real(r8), parameter :: cropfire_a1 = 0.3_r8
    !
    ! c parameter for peatland fire in Li et. al. (2013)
    ! boreal peat fires (was different in paper),changed from /timestep to /hr
    real(r8), parameter :: boreal_peatfire_c = 4.2e-5_r8
    !
    ! non-boreal peat fires (was different in paper)
    real(r8), parameter :: non_boreal_peatfire_c = 0.001_r8
    !
    integer  :: g,t,l,c,p,pi,j,fc,fp,kyr, kmo, kda, mcsec   ! index variables
    real(r8) :: dt       ! time step variable (s)
    real(r8) :: m        ! top-layer soil moisture (proportion)
    real(r8) :: dayspyr  ! days per year
    real(r8) :: cli      ! effect of climate on deforestation fires (0-1)
    real(r8), parameter ::cli_scale = 0.035_r8   !global constant for deforestation fires (/d)
    real(r8) :: cri      ! thresholds used for cli, (mm/d), see Eq.(7) in Li et al.(2013)
    real(r8) :: fb       ! availability of fuel for regs A and C
    real(r8) :: fhd      ! impact of hd on agricultural fire
    real(r8) :: fgdp     ! impact of gdp on agricultural fire
    real(r8) :: fire_m   ! combustability of fuel for fire occurrence
    real(r8) :: spread_m ! combustability of fuel for fire spread
    real(r8) :: Lb_lf    ! length-to-breadth ratio added by Lifang 
    integer  :: i_cwd    ! cwd pool
    real(r8) :: lh       ! anthro. ignitions (count/km2/hr)
    real(r8) :: fs       ! hd-dependent fires suppression (0-1)
    real(r8) :: ig       ! total ignitions (count/km2/hr)
    real(r8) :: hdmlf    ! human density
    real(r8) :: btran_col(bounds%begc:bounds%endc)
    logical  :: transient_landcover  ! whether this run has any prescribed transient landcover
    !-----------------------------------------------------------------------

    associate(                                                                & 
         is_cwd             =>    decomp_cascade_con%is_cwd                 , & ! Input:  [logical  (:)     ]  TRUE => pool is a cwd pool                         
         
         forc_rh            =>    top_as%rhbot                              , & ! Input:  [real(r8) (:)     ]  relative humidity                                 
         forc_wind          =>    top_as%windbot                            , & ! Input:  [real(r8) (:)     ]  atmospheric wind speed (m/s)                       
         forc_t             =>    top_as%tbot                               , & ! Input:  [real(r8) (:)     ]  atmospheric temperature (Kelvin)                  
         forc_rain          =>    top_af%rain                               , & ! Input:  [real(r8) (:)     ]  rain rate (kg H2O/m**2/s, or mm liquid H2O/s)                                              
         forc_snow          =>    top_af%snow                               , & ! Input:  [real(r8) (:)     ]  snow rate (kg H2O/m**2/s, or mm liquid H2O/s)                                            
#ifdef CPL_BYPASS
         forc_hdm           =>    atm2lnd_vars%forc_hdm                     , & ! Input:  [real(r8) (:)     ]  population density
         forc_lnfm          =>    atm2lnd_vars%forc_lnfm                    , & ! Input:  [real(r8) (:)     ]  ligntning data 
#endif 
         prec60             =>    top_af%prec60d                            , & ! Input:  [real(r8) (:)     ]  60-day running mean of tot. precipitation, mm liquid H2O/s 
         prec10             =>    top_af%prec10d                            , & ! Input:  [real(r8) (:)     ]  10-day running mean of tot. precipitation, mm liquid H2O/s 
         
         tsoi17             =>    col_es%t_soi17cm            , & ! Input:  [real(r8) (:)     ]  soil T for top 0.17 m                             
         
         btran2             =>    energyflux_vars%btran2_patch              , & ! Input:  [real(r8) (:)     ]  root zone soil wetness                            
         
         fsat               =>    soilhydrology_vars%fsat_col               , & ! Input:  [real(r8) (:)     ]  fractional area with water table at surface       
         
         wf                 =>    col_ws%wf                    , & ! Input:  [real(r8) (:)     ]  soil water as frac. of whc for top 0.05 m         
         wf2                =>    col_ws%wf2                   , & ! Input:  [real(r8) (:)     ]  soil water as frac. of whc for top 0.17 m         
        
         lfpftd             =>    cnstate_vars%lfpftd_patch                 , & ! Input:  [real(r8) (:)     ]  decrease of pft weight (0-1) on the col. for dt
         cropf_col          =>    cnstate_vars%cropf_col                    , & ! Input:  [real(r8) (:)     ]  cropland fraction in veg column                   
         gdp_lf             =>    cnstate_vars%gdp_lf_col                   , & ! Input:  [real(r8) (:)     ]  gdp data                                          
         peatf_lf           =>    cnstate_vars%peatf_lf_col                 , & ! Input:  [real(r8) (:)     ]  peatland fraction data                            
         abm_lf             =>    cnstate_vars%abm_lf_col                   , & ! Input:  [integer  (:)     ]  prescribed crop fire time                          
         baf_crop           =>    cnstate_vars%baf_crop_col                 , & ! Output: [real(r8) (:)     ]  burned area fraction for cropland (/sec)  
         baf_peatf          =>    cnstate_vars%baf_peatf_col                , & ! Output: [real(r8) (:)     ]  burned area fraction for peatland (/sec)  
         burndate           =>    cnstate_vars%burndate_patch               , & ! Output: [integer  (:)     ]  burn date for crop                                 
         fbac               =>    cnstate_vars%fbac_col                     , & ! Output: [real(r8) (:)     ]  total burned area out of conversion (/sec)
         fbac1              =>    cnstate_vars%fbac1_col                    , & ! Output: [real(r8) (:)     ]  burned area out of conversion region due to land use fire
         farea_burned       =>    cnstate_vars%farea_burned_col             , & ! Output: [real(r8) (:)     ]  total fractional area burned (/sec)
         nfire              =>    cnstate_vars%nfire_col                    , & ! Output: [real(r8) (:)     ]  fire counts (count/km2/sec), valid only in Reg. C
         fsr_col            =>    cnstate_vars%fsr_col                      , & ! Output: [real(r8) (:)     ]  fire spread rate at column level
         fd_col             =>    cnstate_vars%fd_col                       , & ! Output: [real(r8) (:)     ]  fire duration rate at column level
         lgdp_col           =>    cnstate_vars%lgdp_col                     , & ! Output: [real(r8) (:)     ]  gdp limitation factor for nfire                   
         lgdp1_col          =>    cnstate_vars%lgdp1_col                    , & ! Output: [real(r8) (:)     ]  gdp limitation factor for baf per fire            
         lpop_col           =>    cnstate_vars%lpop_col                     , & ! Output: [real(r8) (:)     ]  pop limitation factor for baf per fire            
         lfwt               =>    cnstate_vars%lfwt_col                     , & ! Output: [real(r8) (:)     ]  fractional coverage of non-crop and non-bare-soil Patches
         trotr1_col         =>    cnstate_vars%trotr1_col                   , & ! Output: [real(r8) (:)     ]  pft weight of BET on the gridcell (0-1)           
         trotr2_col         =>    cnstate_vars%trotr2_col                   , & ! Output: [real(r8) (:)     ]  pft weight of BDT on the gridcell (0-1)           
         dtrotr_col         =>    cnstate_vars%dtrotr_col                   , & ! Output: [real(r8) (:)     ]  decreased frac. coverage of BET+BDT on grid for dt
         lfc                =>    cnstate_vars%lfc_col                      , & ! Output: [real(r8) (:)     ]  conversion area frac. of BET+BDT that haven't burned before
         wtlf               =>    cnstate_vars%wtlf_col                     , & ! Output: [real(r8) (:)     ]  fractional coverage of non-crop Patches              
         
         deadcrootc         =>    veg_cs%deadcrootc         , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
         deadcrootc_storage =>    veg_cs%deadcrootc_storage , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
         deadcrootc_xfer    =>    veg_cs%deadcrootc_xfer    , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
         frootc             =>    veg_cs%frootc             , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C                               
         frootc_storage     =>    veg_cs%frootc_storage     , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
         frootc_xfer        =>    veg_cs%frootc_xfer        , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
         livecrootc         =>    veg_cs%livecrootc         , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
         livecrootc_storage =>    veg_cs%livecrootc_storage , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
         livecrootc_xfer    =>    veg_cs%livecrootc_xfer    , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
         totvegc            =>    veg_cs%totvegc            , & ! Input:  [real(r8) (:)     ]  (gC/m2) total vegetation carbon, excluding cpool  
         leafc              =>    veg_cs%leafc              , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
         leafc_storage      =>    veg_cs%leafc_storage      , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
         leafc_xfer         =>    veg_cs%leafc_xfer         , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
         deadstemc          =>    veg_cs%deadstemc          , & ! Input:[real(r8) (:)       ]  (gC/m2) dead stem C             
         totlitc            =>    col_cs%totlitc              , & ! Input:  [real(r8) (:)     ]  (gC/m2) total lit C (column-level mean)           
         decomp_cpools_vr   =>    col_cs%decomp_cpools_vr     , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  VR decomp. (litter, cwd, soil)
         rootc_col          =>    col_cs%rootc                , & ! Output: [real(r8) (:)     ]  root carbon                                       
         totvegc_col        =>    col_cs%totvegc              , & ! Output: [real(r8) (:)     ]  totvegc at column level                           
         leafc_col          =>    col_cs%leafc                , & ! Output: [real(r8) (:)     ]  leaf carbon at column level                       
         deadstemc_col      =>    col_cs%deadstemc            , & ! Output: [real(r8) (:)     ] dead stem carbon at column level
         fuelc              =>    col_cs%fuelc                , & ! Output: [real(r8) (:)     ]  fuel avalability factor for Reg.C                 
         fuelc_crop         =>    col_cs%fuelc_crop             & ! Output: [real(r8) (:)     ]  fuel avalability factor for Reg.A                 
         )
 
      transient_landcover = run_has_transient_landcover()

      !pft to column average 
      call p2c(bounds, num_soilc, filter_soilc, &
           totvegc(bounds%begp:bounds%endp), &
           totvegc_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           leafc(bounds%begp:bounds%endp), &
           leafc_col(bounds%begc:bounds%endc))
     
      call p2c(bounds, num_soilc, filter_soilc, &
           deadstemc(bounds%begp:bounds%endp), &
           deadstemc_col(bounds%begc:bounds%endc))

     call get_curr_date (kyr, kmo, kda, mcsec)
     dayspyr = get_days_per_year()
     ! Get model step size
     dt      = real( get_step_size(), r8 )
     !
     ! On first time-step, just set area burned to zero and exit
     !
     if ( get_nstep() == 0 )then
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           farea_burned(c) = 0._r8
           baf_crop(c)     = 0._r8
           baf_peatf(c)    = 0._r8
           fbac(c)         = 0._r8
           fbac1(c)        = 0._r8
           cropf_col(c)    = 0._r8 
        end do
        return
     end if
     !
     ! Calculate fraction of crop (cropf_col) and non-crop and non-bare-soil 
     ! vegetation (lfwt) in vegetated column
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        cropf_col(c) = 0._r8 
        lfwt(c)      = 0._r8   
     end do
     do pi = 1,max_patch_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           if (pi <=  col_pp%npfts(c)) then
              p = col_pp%pfti(c) + pi - 1
              ! For crop veg types
              if( veg_pp%itype(p) > nc4_grass )then
                 cropf_col(c) = cropf_col(c) + veg_pp%wtcol(p)
              end if
              ! For natural vegetation (non-crop and non-bare-soil)
              if( veg_pp%itype(p) >= ndllf_evr_tmp_tree .and. veg_pp%itype(p) <= nc4_grass )then
                 lfwt(c) = lfwt(c) + veg_pp%wtcol(p)
              end if
           end if
        end do
     end do
     ! 
     ! Calculate crop fuel   
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        fuelc_crop(c)=0._r8
     end do
     do pi = 1,max_patch_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           if (pi <=  col_pp%npfts(c)) then
              p = col_pp%pfti(c) + pi - 1
              ! For crop PFTs, fuel load includes leaf and litter; only
              ! column-level litter carbon
              ! is available, so we use leaf carbon to estimate the
              ! litter carbon for crop PFTs
              if( veg_pp%itype(p) > nc4_grass .and. veg_pp%wtcol(p) > 0._r8 .and. leafc_col(c) > 0._r8 )then
                 fuelc_crop(c)=fuelc_crop(c) + (leafc(p) + leafc_storage(p) + &
                      leafc_xfer(p))*veg_pp%wtcol(p)/cropf_col(c)     + &
                      totlitc(c)*leafc(p)/leafc_col(c)*veg_pp%wtcol(p)/cropf_col(c)
              end if
           end if
        end do
     end do
     !   
     ! Calculate noncrop column variables
     ! 5/22/2018, PET: switched the use of column-weight-on-gridcell, to column-weight-on-topounit
     ! in the calculation of summed weight for tropical trees
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        fsr_col(c)   = 0._r8
        fd_col(c)    = 0._r8
        rootc_col(c) = 0._r8
        lgdp_col(c)  = 0._r8
        lgdp1_col(c) = 0._r8
        lpop_col(c)  = 0._r8
        btran_col(c) = 0._r8
        wtlf(c)      = 0._r8
        trotr1_col(c)= 0._r8
        trotr2_col(c)= 0._r8
        if (transient_landcover) then    !true when landuse data is used
           dtrotr_col(c)=0._r8
        end if
     end do
     do pi = 1,max_patch_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           g = col_pp%gridcell(c)
           if (pi <=  col_pp%npfts(c)) then
              p = col_pp%pfti(c) + pi - 1

              ! For non-crop -- natural vegetation and bare-soil
              if( veg_pp%itype(p)  <  nc3crop .and. cropf_col(c)  <  1.0_r8 )then
                 if( .not. shr_infnan_isnan(btran2(p))) then
                    if (btran2(p)  <=  1._r8 ) then
                       btran_col(c) = btran_col(c)+btran2(p)*veg_pp%wtcol(p)
                       wtlf(c)      = wtlf(c)+veg_pp%wtcol(p)
                    end if
                 end if
                 if( veg_pp%itype(p) == nbrdlf_evr_trp_tree .and. veg_pp%wtcol(p)  >  0._r8 )then
                    trotr1_col(c)=trotr1_col(c)+veg_pp%wtcol(p)*col_pp%wttopounit(c)
                 end if
                 if( veg_pp%itype(p) == nbrdlf_dcd_trp_tree .and. veg_pp%wtcol(p)  >  0._r8 )then
                    trotr2_col(c)=trotr2_col(c)+veg_pp%wtcol(p)*col_pp%wttopounit(c)
                 end if
                 if (transient_landcover) then    !true when landuse data is used
                    if( veg_pp%itype(p) == nbrdlf_evr_trp_tree .or. veg_pp%itype(p) == nbrdlf_dcd_trp_tree )then
                       if(lfpftd(p) > 0._r8)then
                          dtrotr_col(c)=dtrotr_col(c)+lfpftd(p)*col_pp%wttopounit(c)
                       end if
                    end if
                 end if
                 rootc_col(c) = rootc_col(c) + (frootc(p) + frootc_storage(p) + &
                      frootc_xfer(p) + deadcrootc(p) +                &
                      deadcrootc_storage(p) + deadcrootc_xfer(p) +    &
                      livecrootc(p)+livecrootc_storage(p) +           &
                      livecrootc_xfer(p))*veg_pp%wtcol(p)

                 fsr_col(c) = fsr_col(c) + fsr_pft(veg_pp%itype(p))*veg_pp%wtcol(p)/(1.0_r8-cropf_col(c))

                 if( lfwt(c)  /=  0.0_r8 )then    
                    hdmlf=forc_hdm(g)

                    ! all these constants are in Li et al. BG (2012a,b;2013)

                    if( hdmlf  >  0.1_r8 )then            
                       ! For NOT bare-soil
                       if( veg_pp%itype(p)  /=  noveg )then
                          ! For shrub and grass (crop already excluded above)
                          if( veg_pp%itype(p)  >=  nbrdlf_evr_shrub )then      !for shurb and grass
                             lgdp_col(c)  = lgdp_col(c) + (0.1_r8 + 0.9_r8*    &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (gdp_lf(c)/8._r8)**0.5_r8))*veg_pp%wtcol(p) &
                                  /(1.0_r8 - cropf_col(c))
                             lgdp1_col(c) = lgdp1_col(c) + (0.2_r8 + 0.8_r8*   &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (gdp_lf(c)/7._r8)))*veg_pp%wtcol(p)/lfwt(c)
                             lpop_col(c)  = lpop_col(c) + (0.2_r8 + 0.8_r8*    &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (hdmlf/450._r8)**0.5_r8))*veg_pp%wtcol(p)/lfwt(c)
                          else   ! for trees
                             if( gdp_lf(c)  >  20._r8 )then
                                lgdp_col(c)  =lgdp_col(c)+0.39_r8*veg_pp%wtcol(p)/(1.0_r8 - cropf_col(c))
                             else    
                                lgdp_col(c) = lgdp_col(c)+veg_pp%wtcol(p)/(1.0_r8 - cropf_col(c))
                             end if
                             if( gdp_lf(c)  >  20._r8 )then   
                                lgdp1_col(c) = lgdp1_col(c)+0.62_r8*veg_pp%wtcol(p)/lfwt(c)
                             else
                                if( gdp_lf(c)  >  8._r8 ) then
                                   lgdp1_col(c)=lgdp1_col(c)+0.83_r8*veg_pp%wtcol(p)/lfwt(c)
                                else
                                   lgdp1_col(c)=lgdp1_col(c)+veg_pp%wtcol(p)/lfwt(c)
                                end if
                             end if
                             lpop_col(c) = lpop_col(c) + (0.4_r8 + 0.6_r8*    &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (hdmlf/125._r8)))*veg_pp%wtcol(p)/lfwt(c) 
                          end if
                       end if
                    else
                       lgdp_col(c)  = lgdp_col(c)+veg_pp%wtcol(p)/(1.0_r8 - cropf_col(c))
                       lgdp1_col(c) = lgdp1_col(c)+veg_pp%wtcol(p)/lfwt(c)
                       lpop_col(c)  = lpop_col(c)+veg_pp%wtcol(p)/lfwt(c)
                    end if
                 end if

                 fd_col(c) = fd_col(c) + fd_pft(veg_pp%itype(p)) * veg_pp%wtcol(p) * secsphr / (1.0_r8-cropf_col(c))         
              end if
           end if
        end do
     end do

     ! estimate annual decreased fractional coverage of BET+BDT
     ! land cover conversion in CLM4.5 is the same for each timestep except for the beginning  

     if (transient_landcover) then    !true when landuse data is used
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           if( dtrotr_col(c)  >  0._r8 )then
              if( kmo == 1 .and. kda == 1 .and. mcsec == 0)then
                 lfc(c) = 0._r8
              end if
              if( kmo == 1 .and. kda == 1 .and. mcsec == dt)then
                 lfc(c) = dtrotr_col(c)*dayspyr*secspday/dt
              end if
           else
              lfc(c)=0._r8
           end if
        end do
     end if
     !
     ! calculate burned area fraction in cropland
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        baf_crop(c)=0._r8
     end do

     do fp = 1,num_soilp
        p = filter_soilp(fp)  
        if( kmo == 1 .and. kda == 1 .and. mcsec == 0 )then
           burndate(p) = 10000 ! init. value; actual range [0 365]
        end if
     end do

     do pi = 1,max_patch_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           g = col_pp%gridcell(c)
           t = col_pp%topounit(c)
           hdmlf=forc_hdm(g)
           if (pi <=  col_pp%npfts(c)) then
              p = col_pp%pfti(c) + pi - 1
              ! For crop
              if( forc_t(t)  >=  SHR_CONST_TKFRZ .and. veg_pp%itype(p)  >  nc4_grass .and.  &
                   kmo == abm_lf(c) .and. forc_rain(t)+forc_snow(t) == 0._r8  .and. &
                   burndate(p) >= 999 .and. veg_pp%wtcol(p)  >  0._r8 )then ! catch  crop burn time

                 ! calculate human density impact on ag. fire
                 fhd = 0.04_r8+0.96_r8*exp(-1._r8*SHR_CONST_PI*(hdmlf/350._r8)**0.5_r8)

                 ! calculate impact of GDP on ag. fire
                 fgdp = 0.01_r8+0.99_r8*exp(-1._r8*SHR_CONST_PI*(gdp_lf(c)/10._r8))

                 ! calculate burned area
                 fb   = max(0.0_r8,min(1.0_r8,(fuelc_crop(c)-lfuel)/(ufuel-lfuel)))

                 ! crop fire only for generic crop types at this time
                 ! managed crops are treated as grasses if crop model is turned on
                 baf_crop(c) = baf_crop(c) + cropfire_a1/secsphr*fb*fhd*fgdp*veg_pp%wtcol(p)
                 if( fb*fhd*fgdp*veg_pp%wtcol(p)  >  0._r8)then
                    burndate(p)=kda
                 end if
              end if
           end if
        end do
     end do
     !
     ! calculate peatland fire
     !
     do fc = 1, num_soilc
        c = filter_soilc(fc)
        t = col_pp%topounit(c) 
        g = col_pp%gridcell(c)
        if(grc_pp%latdeg(g) < borealat )then
           baf_peatf(c) = non_boreal_peatfire_c/secsphr*max(0._r8, &
                min(1._r8,(4.0_r8-prec60(t)*secspday)/ &
                4.0_r8))**2*peatf_lf(c)*(1._r8-fsat(c))
        else
           baf_peatf(c) = boreal_peatfire_c/secsphr*exp(-SHR_CONST_PI*(max(wf2(c),0._r8)/0.3_r8))* &
                max(0._r8,min(1._r8,(tsoi17(c)-SHR_CONST_TKFRZ)/10._r8))*peatf_lf(c)* &
                (1._r8-fsat(c))
        end if
     end do
     !
     ! calculate other fires
     !

     ! Set the number of timesteps for e-folding.
     ! When the simulation has run fewer than this number of steps,
     ! re-scale the e-folding time to get a stable early estimate.

     ! find which pool is the cwd pool
     i_cwd = 0
     do l = 1, ndecomp_pools
        if ( is_cwd(l) ) then
           i_cwd = l
        endif
     end do

     !
     ! begin column loop to calculate fractional area affected by fire
     !
     do fc = 1, num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)
        t = col_pp%topounit(c)
        hdmlf=forc_hdm(g)
        if( cropf_col(c)  <  1.0 )then
           if (trotr1_col(c)+trotr2_col(c)>0.6_r8) then
              farea_burned(c)=min(1.0_r8,baf_crop(c)+baf_peatf(c))
           else
              fuelc(c) = totlitc(c)+totvegc_col(c)-rootc_col(c)-fuelc_crop(c)*cropf_col(c)
              if (spinup_state == 1) fuelc(c) = fuelc(c) + ((spinup_mortality_factor - 1._r8)*deadstemc_col(c))
              do j = 1, nlevdecomp 
                if (spinup_state == 1 .and. kyr < 40) then 
                  fuelc(c) = fuelc(c)+decomp_cpools_vr(c,j,i_cwd) *dzsoi_decomp(j) * &
                    decomp_cascade_con%spinup_factor(i_cwd) 
                else if (spinup_state == 1 .and. kyr >= 40) then 
                  fuelc(c) = fuelc(c)+decomp_cpools_vr(c,j,i_cwd) *dzsoi_decomp(j) * &
                    decomp_cascade_con%spinup_factor(i_cwd) / cnstate_vars%scalaravg_col(c,j)
                else  
                  fuelc(c) = fuelc(c)+decomp_cpools_vr(c,j,i_cwd) * dzsoi_decomp(j)
                end if 
              end do
              fuelc(c) = fuelc(c)/(1._r8-cropf_col(c))
              fb       = max(0.0_r8,min(1.0_r8,(fuelc(c)-lfuel)/(ufuel-lfuel)))
              m        = max(0._r8,wf(c))
              fire_m   = exp(-SHR_CONST_PI *(m/0.69_r8)**2)*(1.0_r8 - max(0._r8, &
                   min(1._r8,(forc_rh(t)-30._r8)/(80._r8-30._r8))))*  &
                   min(1._r8,exp(SHR_CONST_PI*(forc_t(t)-SHR_CONST_TKFRZ)/10._r8))
              lh       = 0.0035_r8*6.8_r8*hdmlf**(0.43_r8)/30._r8/24._r8
              fs       = 1._r8-(0.01_r8+0.98_r8*exp(-0.025_r8*hdmlf))
              ig       = (lh+forc_lnfm(g)/(5.16_r8+2.16_r8*cos(3._r8*grc_pp%lat(g)))*0.25_r8)*(1._r8-fs)*(1._r8-cropf_col(c)) 
              nfire(c) = ig/secsphr*fb*fire_m*lgdp_col(c) !fire counts/km2/sec
              Lb_lf    = 1._r8+10.0_r8*(1._r8-EXP(-0.06_r8*forc_wind(t)))
              if ( wtlf(c) > 0.0_r8 )then
                 spread_m = (1.0_r8 - max(0._r8,min(1._r8,(btran_col(c)/wtlf(c)-0.3_r8)/ &
                      (0.7_r8-0.3_r8))))*(1.0_r8-max(0._r8, &
                      min(1._r8,(forc_rh(t)-30._r8)/(80._r8-30._r8))))
              else
                 spread_m = 0.0_r8
              end if
              farea_burned(c) = min(1._r8,(g0*spread_m*fsr_col(c)* &
                   fd_col(c)/1000._r8)**2*lgdp1_col(c)* &
                   lpop_col(c)*nfire(c)*SHR_CONST_PI*Lb_lf+ &
                   baf_crop(c)+baf_peatf(c))  ! fraction (0-1) per sec
           end if
           !
           ! if landuse change data is used, calculate deforestation fires and 
           ! add it in the total of burned area fraction
           !
           if (transient_landcover) then    !true when landuse change data is used
              if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 )then
                 if(( kmo == 1 .and. kda == 1 .and. mcsec == 0) .or. &
                      dtrotr_col(c) <=0._r8 )then
                    fbac1(c)        = 0._r8
                    farea_burned(c) = baf_crop(c)+baf_peatf(c)
                 else
                    cri = (4.0_r8*trotr1_col(c)+1.8_r8*trotr2_col(c))/(trotr1_col(c)+trotr2_col(c))
                    cli = (max(0._r8,min(1._r8,(cri-prec60(t)*secspday)/cri))**0.5)* &
                         (max(0._r8,min(1._r8,(cri-prec10(t)*secspday)/cri))**0.5)* &
                         max(0.0005_r8,min(1._r8,19._r8*dtrotr_col(c)*dayspyr*secspday/dt-0.001_r8))* &
                         max(0._r8,min(1._r8,(0.25_r8-(forc_rain(t)+forc_snow(t))*secsphr)/0.25_r8))
                    farea_burned(c) = cli*(cli_scale/secspday)+baf_crop(c)+baf_peatf(c)
                    ! burned area out of conversion region due to land use fire
                    fbac1(c) = max(0._r8,cli*(cli_scale/secspday) - 2.0_r8*lfc(c)/dt)   
                 end if
                 ! total burned area out of conversion 
                 fbac(c) = fbac1(c)+baf_crop(c)+baf_peatf(c) 
              else
                 fbac(c) = farea_burned(c)
              end if
           end if

        else
           farea_burned(c) = min(1._r8,baf_crop(c)+baf_peatf(c))
        end if

        if (use_nofire) then
           ! zero out the fire area if NOFIRE flag is on

           farea_burned(c) = 0._r8
           baf_crop(c)     = 0._r8
           baf_peatf(c)    = 0._r8
           fbac(c)         = 0._r8
           fbac1(c)        = 0._r8
           ! with NOFIRE, tree carbon is still removed in landuse change regions by the
           ! landuse code
        end if

     end do  ! end of column loop

   end associate

 end subroutine FireArea

 !-----------------------------------------------------------------------
 subroutine FireFluxes (num_soilc, filter_soilc, num_soilp, filter_soilp, &
      cnstate_vars, carbonstate_vars, nitrogenstate_vars, &
      carbonflux_vars,nitrogenflux_vars,phosphorusstate_vars,phosphorusflux_vars)
   !
   ! !DESCRIPTION:
   ! Fire effects routine for coupled carbon-nitrogen code (CN).
   ! Relies primarily on estimate of fractional area burned, from FireArea().
   !
   ! Total fire carbon emissions (g C/m2 land area/yr) 
   !  =avg(COL_FIRE_CLOSS)*seconds_per_year + avg(SOMC_FIRE)*seconds_per_year + 
   !   avg(LF_CONV_CFLUX)*seconds_per_year*min(1.0,avg(LFC2)*seconds_per_year)*0.8
   ! where avg means the temporal average in a year
   ! seconds_per_year is the number of seconds in a year.
   !
   ! !USES:
   use pftvarcon            , only: cc_leaf,cc_lstem,cc_dstem,cc_other,fm_leaf,fm_lstem,fm_other,fm_root,fm_lroot,fm_droot
   use pftvarcon            , only: nc3crop,lf_flab,lf_fcel,lf_flig,fr_flab,fr_fcel,fr_flig
   use clm_time_manager     , only: get_step_size,get_days_per_year,get_curr_date
   use elm_varpar           , only: max_patch_per_col
   use elm_varctl           , only: spinup_state, spinup_mortality_factor
   use dynSubgridControlMod , only: get_flanduse_timeseries
   use elm_varcon           , only: secspday
   use dynSubgridControlMod , only: run_has_transient_landcover
   !
   ! !ARGUMENTS:
   integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
   integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
   integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
   type(cnstate_type)       , intent(inout) :: cnstate_vars
   type(carbonstate_type)   , intent(inout) :: carbonstate_vars
   type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
   type(carbonflux_type)    , intent(inout) :: carbonflux_vars
   type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
  
   ! ! adding phosphorus state and flux variables
   type(phosphorusstate_type) , intent(in)    :: phosphorusstate_vars
   type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
   !
   ! !LOCAL VARIABLES:
   integer :: g,c,p,j,l,pi,kyr, kmo, kda, mcsec   ! indices
   real(r8):: m_veg                ! speedup factor for accelerated decomp
   integer :: fp,fc                ! filter indices
   real(r8):: f                    ! rate for fire effects (1/s)
   real(r8):: dt                   ! time step variable (s)
   real(r8):: dayspyr              ! days per year

   real(r8), pointer :: m_leafn_to_litter_fire                  (:)
   real(r8), pointer :: m_leafn_storage_to_litter_fire          (:)
   real(r8), pointer :: m_leafn_xfer_to_litter_fire             (:)
   real(r8), pointer :: m_livestemn_to_litter_fire              (:)
   real(r8), pointer :: m_livestemn_storage_to_litter_fire      (:)
   real(r8), pointer :: m_livestemn_xfer_to_litter_fire         (:)
   real(r8), pointer :: m_livestemn_to_deadstemn_fire           (:)
   real(r8), pointer :: m_deadstemn_to_litter_fire              (:)
   real(r8), pointer :: m_deadstemn_storage_to_litter_fire      (:)
   real(r8), pointer :: m_deadstemn_xfer_to_litter_fire         (:)
   real(r8), pointer :: m_frootn_to_litter_fire                 (:)
   real(r8), pointer :: m_frootn_storage_to_litter_fire         (:)
   real(r8), pointer :: m_frootn_xfer_to_litter_fire            (:)
   real(r8), pointer :: m_livecrootn_to_litter_fire             (:)
   real(r8), pointer :: m_livecrootn_storage_to_litter_fire     (:)
   real(r8), pointer :: m_livecrootn_xfer_to_litter_fire        (:)
   real(r8), pointer :: m_livecrootn_to_deadcrootn_fire         (:)
   real(r8), pointer :: m_deadcrootn_to_litter_fire             (:)
   real(r8), pointer :: m_deadcrootn_storage_to_litter_fire     (:)
   real(r8), pointer :: m_deadcrootn_xfer_to_litter_fire        (:)
   real(r8), pointer :: m_retransn_to_litter_fire               (:)
   real(r8), pointer :: m_npool_to_litter_fire                  (:)
   real(r8), pointer :: m_decomp_npools_to_fire_vr              (:,:,:)
   real(r8), pointer :: m_n_to_litr_met_fire                    (:,:)
   real(r8), pointer :: m_n_to_litr_cel_fire                    (:,:)
   real(r8), pointer :: m_n_to_litr_lig_fire                    (:,:)

   real(r8), pointer :: m_leafp_to_litter_fire                  (:)
   real(r8), pointer :: m_leafp_storage_to_litter_fire          (:)
   real(r8), pointer :: m_leafp_xfer_to_litter_fire             (:)
   real(r8), pointer :: m_livestemp_to_litter_fire              (:)
   real(r8), pointer :: m_livestemp_storage_to_litter_fire      (:)
   real(r8), pointer :: m_livestemp_xfer_to_litter_fire         (:)
   real(r8), pointer :: m_livestemp_to_deadstemp_fire           (:)
   real(r8), pointer :: m_deadstemp_to_litter_fire              (:)
   real(r8), pointer :: m_deadstemp_storage_to_litter_fire      (:)
   real(r8), pointer :: m_deadstemp_xfer_to_litter_fire         (:)
   real(r8), pointer :: m_frootp_to_litter_fire                 (:)
   real(r8), pointer :: m_frootp_storage_to_litter_fire         (:)
   real(r8), pointer :: m_frootp_xfer_to_litter_fire            (:)
   real(r8), pointer :: m_livecrootp_to_litter_fire             (:)
   real(r8), pointer :: m_livecrootp_storage_to_litter_fire     (:)
   real(r8), pointer :: m_livecrootp_xfer_to_litter_fire        (:)
   real(r8), pointer :: m_livecrootp_to_deadcrootp_fire         (:)
   real(r8), pointer :: m_deadcrootp_to_litter_fire             (:)
   real(r8), pointer :: m_deadcrootp_storage_to_litter_fire     (:)
   real(r8), pointer :: m_deadcrootp_xfer_to_litter_fire        (:)
   real(r8), pointer :: m_retransp_to_litter_fire               (:)
   real(r8), pointer :: m_ppool_to_litter_fire                  (:)
   real(r8), pointer :: m_decomp_ppools_to_fire_vr              (:,:,:)
   real(r8), pointer :: m_p_to_litr_met_fire                    (:,:)
   real(r8), pointer :: m_p_to_litr_cel_fire                    (:,:)
   real(r8), pointer :: m_p_to_litr_lig_fire                    (:,:)

   logical           :: transient_landcover  ! whether this run has any prescribed transient landcover

   !-----------------------------------------------------------------------

   ! NOTE: VR      = Vertically Resolved
   !       conv.   = conversion
   !       frac.   = fraction
   !       BAF     = Burned Area Fraction
   !       ann.    = annual
   !       GC      = gridcell
   !       dt      = timestep
   !       C       = Carbon
   !       N       = Nitrogen
   !       emis.   = emissions
   !       decomp. = decomposing

   associate(                                                                                                   & 
        is_cwd                              =>    decomp_cascade_con%is_cwd                                   , & ! Input:  [logical  (:)     ]  TRUE => pool is a cwd pool                         
        is_litter                           =>    decomp_cascade_con%is_litter                                , & ! Input:  [logical  (:)     ]  TRUE => pool is a litter pool                      

        woody                               =>    veg_vp%woody                                            , & ! Input:  [real(r8) (:)     ]  woody lifeform (1=woody, 0=not woody)             

        cropf_col                           =>    cnstate_vars%cropf_col                                      , & ! Input:  [real(r8) (:)     ]  cropland fraction in veg column                   
        croot_prof                          =>    cnstate_vars%croot_prof_patch                               , & ! Input:  [real(r8) (:,:)   ]  (1/m) profile of coarse roots                   
        stem_prof                           =>    cnstate_vars%stem_prof_patch                                , & ! Input:  [real(r8) (:,:)   ]  (1/m) profile of stems                          
        froot_prof                          =>    cnstate_vars%froot_prof_patch                               , & ! Input:  [real(r8) (:,:)   ]  (1/m) profile of fine roots                     
        leaf_prof                           =>    cnstate_vars%leaf_prof_patch                                , & ! Input:  [real(r8) (:,:)   ]  (1/m) profile of leaves                         
        farea_burned                        =>    cnstate_vars%farea_burned_col                               , & ! Input:  [real(r8) (:)     ]  fractional area burned (/sec)
        lfc                                 =>    cnstate_vars%lfc_col                                        , & ! Input:  [real(r8) (:)     ]  conv. area frac. of BET+BDT that haven't burned before
        lfc2                                =>    cnstate_vars%lfc2_col                                       , & ! Output: [real(r8) (:)     ]  conv. area frac. of BET+BDT burned this dt (/sec)
        fbac1                               =>    cnstate_vars%fbac1_col                                      , & ! Input:  [real(r8) (:)     ]  burned area out of conv. region due to LU fire 
        fbac                                =>    cnstate_vars%fbac_col                                       , & ! Input:  [real(r8) (:)     ]  total burned area out of conversion (/sec)
        baf_crop                            =>    cnstate_vars%baf_crop_col                                   , & ! Input:  [real(r8) (:)     ]  BAF for cropland                                  
        baf_peatf                           =>    cnstate_vars%baf_peatf_col                                  , & ! Input:  [real(r8) (:)     ]  BAF for peatlabd                                  
        trotr1_col                          =>    cnstate_vars%trotr1_col                                     , & ! Input:  [real(r8) (:)     ]  pft weight of BET on the gridcell (0-1)           
        trotr2_col                          =>    cnstate_vars%trotr2_col                                     , & ! Input:  [real(r8) (:)     ]  pft weight of BDT on the gridcell (0-1)           
        dtrotr_col                          =>    cnstate_vars%dtrotr_col                                     , & ! Input:  [real(r8) (:)     ]  ann. decreased frac. coverage of BET+BDT (0-1) on GC
        
        decomp_cpools_vr                    =>    col_cs%decomp_cpools_vr                       , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  VR decomp. (litter, cwd, soil)
        totsomc                             =>    col_cs%totsomc                                , & ! Input:  [real(r8) (:)     ]  (gC/m2) total soil organic matter C
        leafcmax                            =>    veg_cs%leafcmax                             , & ! Output: [real(r8) (:)     ]  (gC/m2) ann max leaf C                            
        leafc                               =>    veg_cs%leafc                                , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
        leafc_storage                       =>    veg_cs%leafc_storage                        , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
        leafc_xfer                          =>    veg_cs%leafc_xfer                           , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
        livestemc                           =>    veg_cs%livestemc                            , & ! Input:  [real(r8) (:)     ]  (gC/m2) live stem C                               
        livestemc_storage                   =>    veg_cs%livestemc_storage                    , & ! Input:  [real(r8) (:)     ]  (gC/m2) live stem C storage                       
        livestemc_xfer                      =>    veg_cs%livestemc_xfer                       , & ! Input:  [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
        deadstemc                           =>    veg_cs%deadstemc                            , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead stem C                               
        deadstemc_storage                   =>    veg_cs%deadstemc_storage                    , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
        deadstemc_xfer                      =>    veg_cs%deadstemc_xfer                       , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
        frootc                              =>    veg_cs%frootc                               , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C                               
        frootc_storage                      =>    veg_cs%frootc_storage                       , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
        frootc_xfer                         =>    veg_cs%frootc_xfer                          , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
        livecrootc                          =>    veg_cs%livecrootc                           , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
        livecrootc_storage                  =>    veg_cs%livecrootc_storage                   , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
        livecrootc_xfer                     =>    veg_cs%livecrootc_xfer                      , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
        deadcrootc                          =>    veg_cs%deadcrootc                           , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
        deadcrootc_storage                  =>    veg_cs%deadcrootc_storage                   , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
        deadcrootc_xfer                     =>    veg_cs%deadcrootc_xfer                      , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
        gresp_storage                       =>    veg_cs%gresp_storage                        , & ! Input:  [real(r8) (:)     ]  (gC/m2) growth respiration storage                
        gresp_xfer                          =>    veg_cs%gresp_xfer                           , & ! Input:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
        cpool                               =>    veg_cs%cpool                                , & ! Input:  [real(r8) (:)     ]  (gC/m2) C pool        

        decomp_npools_vr                    =>    col_ns%decomp_npools_vr                     , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  VR decomp. (litter, cwd, soil)
        leafn                               =>    veg_ns%leafn                              , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N                                    
        leafn_storage                       =>    veg_ns%leafn_storage                      , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
        leafn_xfer                          =>    veg_ns%leafn_xfer                         , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N transfer                           
        livestemn                           =>    veg_ns%livestemn                          , & ! Input:  [real(r8) (:)     ]  (gN/m2) live stem N                               
        livestemn_storage                   =>    veg_ns%livestemn_storage                  , & ! Input:  [real(r8) (:)     ]  (gN/m2) live stem N storage                       
        livestemn_xfer                      =>    veg_ns%livestemn_xfer                     , & ! Input:  [real(r8) (:)     ]  (gN/m2) live stem N transfer                      
        deadstemn                           =>    veg_ns%deadstemn                          , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead stem N                               
        deadstemn_storage                   =>    veg_ns%deadstemn_storage                  , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead stem N storage                       
        deadstemn_xfer                      =>    veg_ns%deadstemn_xfer                     , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead stem N transfer                      
        frootn                              =>    veg_ns%frootn                             , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N                               
        frootn_storage                      =>    veg_ns%frootn_storage                     , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
        frootn_xfer                         =>    veg_ns%frootn_xfer                        , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N transfer                      
        livecrootn                          =>    veg_ns%livecrootn                         , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N                        
        livecrootn_storage                  =>    veg_ns%livecrootn_storage                 , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
        livecrootn_xfer                     =>    veg_ns%livecrootn_xfer                    , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N transfer               
        deadcrootn                          =>    veg_ns%deadcrootn                         , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead coarse root N                        
        deadcrootn_storage                  =>    veg_ns%deadcrootn_storage                 , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead coarse root N storage                
        deadcrootn_xfer                     =>    veg_ns%deadcrootn_xfer                    , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead coarse root N transfer               
        retransn                            =>    veg_ns%retransn                           , & ! Input:  [real(r8) (:)     ]  (gN/m2) plant pool of retranslocated N            
        npool                               =>    veg_ns%npool                              , & ! Input:  [real(r8) (:)     ]  (gN/m2) plant pool of stored N  
        ! add phosphorus state variables - X.YANG
        decomp_ppools_vr                    =>    col_ps%decomp_ppools_vr                     , & ! Input:  [real(r8) (:,:,:) ]  (gP/m3)  VR decomp. (litter, cwd, soil)
        leafp                               =>    veg_ps%leafp                              , & ! Input:  [real(r8) (:)     ]  (gP/m2) leaf P                                    
        leafp_storage                       =>    veg_ps%leafp_storage                      , & ! Input:  [real(r8) (:)     ]  (gP/m2) leaf P storage                            
        leafp_xfer                          =>    veg_ps%leafp_xfer                         , & ! Input:  [real(r8) (:)     ]  (gP/m2) leaf P transfer                           
        livestemp                           =>    veg_ps%livestemp                          , & ! Input:  [real(r8) (:)     ]  (gP/m2) live stem P                               
        livestemp_storage                   =>    veg_ps%livestemp_storage                  , & ! Input:  [real(r8) (:)     ]  (gP/m2) live stem P storage                       
        livestemp_xfer                      =>    veg_ps%livestemp_xfer                     , & ! Input:  [real(r8) (:)     ]  (gP/m2) live stem P transfer                      
        deadstemp                           =>    veg_ps%deadstemp                          , & ! Input:  [real(r8) (:)     ]  (gP/m2) dead stem P                               
        deadstemp_storage                   =>    veg_ps%deadstemp_storage                  , & ! Input:  [real(r8) (:)     ]  (gP/m2) dead stem P storage                       
        deadstemp_xfer                      =>    veg_ps%deadstemp_xfer                     , & ! Input:  [real(r8) (:)     ]  (gP/m2) dead stem P transfer                      
        frootp                              =>    veg_ps%frootp                             , & ! Input:  [real(r8) (:)     ]  (gP/m2) fine root P                               
        frootp_storage                      =>    veg_ps%frootp_storage                     , & ! Input:  [real(r8) (:)     ]  (gP/m2) fine root P storage                       
        frootp_xfer                         =>    veg_ps%frootp_xfer                        , & ! Input:  [real(r8) (:)     ]  (gP/m2) fine root P transfer                      
        livecrootp                          =>    veg_ps%livecrootp                         , & ! Input:  [real(r8) (:)     ]  (gP/m2) live coarse root P                        
        livecrootp_storage                  =>    veg_ps%livecrootp_storage                 , & ! Input:  [real(r8) (:)     ]  (gP/m2) live coarse root P storage                
        livecrootp_xfer                     =>    veg_ps%livecrootp_xfer                    , & ! Input:  [real(r8) (:)     ]  (gP/m2) live coarse root P transfer               
        deadcrootp                          =>    veg_ps%deadcrootp                         , & ! Input:  [real(r8) (:)     ]  (gP/m2) dead coarse root P                        
        deadcrootp_storage                  =>    veg_ps%deadcrootp_storage                 , & ! Input:  [real(r8) (:)     ]  (gP/m2) dead coarse root P storage                
        deadcrootp_xfer                     =>    veg_ps%deadcrootp_xfer                    , & ! Input:  [real(r8) (:)     ]  (gP/m2) dead coarse root P transfer               
        retransp                            =>    veg_ps%retransp                           , & ! Input:  [real(r8) (:)     ]  (gP/m2) plant pool of retranslocated P            
        ppool                               =>    veg_ps%ppool                              , & ! Input:  [real(r8) (:)     ]  (gP/m2) plant pool of storage P            

        fire_mortality_c_to_cwdc            =>    col_cf%fire_mortality_c_to_cwdc                , & ! Input:  [real(r8) (:,:)   ]  C flux fire mortality to CWD (gC/m3/s)
        somc_fire                           =>    col_cf%somc_fire                               , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) fire C emissions due to peat burning
        m_leafc_to_fire                     =>    veg_cf%m_leafc_to_fire                       , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) fire C emis. from leafc	    
        m_leafc_storage_to_fire             =>    veg_cf%m_leafc_storage_to_fire               , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) fire C emis. from leafc_storage   
        m_leafc_xfer_to_fire                =>    veg_cf%m_leafc_xfer_to_fire                  , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) fire C emis. from leafc_xfer	    
        m_livestemc_to_fire                 =>    veg_cf%m_livestemc_to_fire                   , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) fire C emis. from livestemc	         
        m_livestemc_storage_to_fire         =>    veg_cf%m_livestemc_storage_to_fire           , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. livestemc_storage	       
        m_livestemc_xfer_to_fire            =>    veg_cf%m_livestemc_xfer_to_fire              , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. livestemc_xfer	       
        m_deadstemc_to_fire                 =>    veg_cf%m_deadstemc_to_fire                   , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadstemc_xfer	       
        m_deadstemc_storage_to_fire         =>    veg_cf%m_deadstemc_storage_to_fire           , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadstemc_storage	       
        m_deadstemc_xfer_to_fire            =>    veg_cf%m_deadstemc_xfer_to_fire              , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadstemc_xfer	       
        m_frootc_to_fire                    =>    veg_cf%m_frootc_to_fire                      , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. frootc		       
        m_frootc_storage_to_fire            =>    veg_cf%m_frootc_storage_to_fire              , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. frootc_storage	       
        m_frootc_xfer_to_fire               =>    veg_cf%m_frootc_xfer_to_fire                 , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. frootc_xfer		       
        m_livecrootc_to_fire                =>    veg_cf%m_livecrootc_to_fire                  , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. livecrootc		    	
        m_livecrootc_storage_to_fire        =>    veg_cf%m_livecrootc_storage_to_fire          , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. livecrootc_storage	       
        m_livecrootc_xfer_to_fire           =>    veg_cf%m_livecrootc_xfer_to_fire             , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. livecrootc_xfer	       
        m_deadcrootc_to_fire                =>    veg_cf%m_deadcrootc_to_fire                  , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadcrootc		    	
        m_deadcrootc_storage_to_fire        =>    veg_cf%m_deadcrootc_storage_to_fire          , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadcrootc_storage	       
        m_deadcrootc_xfer_to_fire           =>    veg_cf%m_deadcrootc_xfer_to_fire             , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadcrootc_xfer	       
        m_gresp_storage_to_fire             =>    veg_cf%m_gresp_storage_to_fire               , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. gresp_storage	
        m_gresp_xfer_to_fire                =>    veg_cf%m_gresp_xfer_to_fire                  , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. gresp_xfer           
        m_cpool_to_fire                     =>    veg_cf%m_cpool_to_fire                       , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. cpool     
  
        fire_mortality_n_to_cwdn            =>    col_nf%fire_mortality_n_to_cwdn              , & ! Input:  [real(r8) (:,:)   ]  N flux fire mortality to CWD (gN/m3/s)
        m_leafn_to_fire                     =>    veg_nf%m_leafn_to_fire                     , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. leafn		  
        m_leafn_storage_to_fire             =>    veg_nf%m_leafn_storage_to_fire             , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. leafn_storage	  
        m_leafn_xfer_to_fire                =>    veg_nf%m_leafn_xfer_to_fire                , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. leafn_xfer	       
        m_livestemn_to_fire                 =>    veg_nf%m_livestemn_to_fire                 , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. livestemn	       
        m_livestemn_storage_to_fire         =>    veg_nf%m_livestemn_storage_to_fire         , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. livestemn_s	       
        m_livestemn_xfer_to_fire            =>    veg_nf%m_livestemn_xfer_to_fire            , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. livestemn_xfer       
        m_deadstemn_to_fire                 =>    veg_nf%m_deadstemn_to_fire                 , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadstemn	       
        m_deadstemn_storage_to_fire         =>    veg_nf%m_deadstemn_storage_to_fire         , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadstemn_storage    
        m_deadstemn_xfer_to_fire            =>    veg_nf%m_deadstemn_xfer_to_fire            , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadstemn_xfer       
        m_frootn_to_fire                    =>    veg_nf%m_frootn_to_fire                    , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. frootn	       
        m_frootn_storage_to_fire            =>    veg_nf%m_frootn_storage_to_fire            , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. frootn_storage       
        m_frootn_xfer_to_fire               =>    veg_nf%m_frootn_xfer_to_fire               , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. frootn_xfer	       
        m_livecrootn_to_fire                =>    veg_nf%m_livecrootn_to_fire                , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. m_livecrootn_to_fire 
        m_livecrootn_storage_to_fire        =>    veg_nf%m_livecrootn_storage_to_fire        , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. livecrootn_storage   
        m_livecrootn_xfer_to_fire           =>    veg_nf%m_livecrootn_xfer_to_fire           , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. livecrootn_xfer      
        m_deadcrootn_to_fire                =>    veg_nf%m_deadcrootn_to_fire                , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadcrootn	       
        m_deadcrootn_storage_to_fire        =>    veg_nf%m_deadcrootn_storage_to_fire        , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadcrootn_storage   
        m_deadcrootn_xfer_to_fire           =>    veg_nf%m_deadcrootn_xfer_to_fire           , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadcrootn_xfer      
        m_retransn_to_fire                  =>    veg_nf%m_retransn_to_fire                  , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. retransn             
        m_npool_to_fire                     =>    veg_nf%m_npool_to_fire                     , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. npooln             
       

        ! add phosphorus fluxes with fire - these will be adding to 
        ! -.XYANG 
        fire_mortality_p_to_cwdp            =>    col_pf%fire_mortality_p_to_cwdp              , & ! Input:  [real(r8) (:,:)   ]  P flux fire mortality to CWD (gP/m3/s)
        m_leafp_to_fire                     =>    veg_pf%m_leafp_to_fire                     , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. leafp		  
        m_leafp_storage_to_fire             =>    veg_pf%m_leafp_storage_to_fire             , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. leafp_storage	  
        m_leafp_xfer_to_fire                =>    veg_pf%m_leafp_xfer_to_fire                , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. leafp_xfer	       
        m_livestemp_to_fire                 =>    veg_pf%m_livestemp_to_fire                 , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. livestemp	       
        m_livestemp_storage_to_fire         =>    veg_pf%m_livestemp_storage_to_fire         , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. livestemp_s	       
        m_livestemp_xfer_to_fire            =>    veg_pf%m_livestemp_xfer_to_fire            , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. livestemp_xfer       
        m_deadstemp_to_fire                 =>    veg_pf%m_deadstemp_to_fire                 , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. deadstemp	       
        m_deadstemp_storage_to_fire         =>    veg_pf%m_deadstemp_storage_to_fire         , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. deadstemp_storage    
        m_deadstemp_xfer_to_fire            =>    veg_pf%m_deadstemp_xfer_to_fire            , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. deadstemp_xfer       
        m_frootp_to_fire                    =>    veg_pf%m_frootp_to_fire                    , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. frootp	       
        m_frootp_storage_to_fire            =>    veg_pf%m_frootp_storage_to_fire            , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. frootp_storage       
        m_frootp_xfer_to_fire               =>    veg_pf%m_frootp_xfer_to_fire               , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. frootp_xfer	       
        m_livecrootp_to_fire                =>    veg_pf%m_livecrootp_to_fire                , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. m_livecrootp_to_fire 
        m_livecrootp_storage_to_fire        =>    veg_pf%m_livecrootp_storage_to_fire        , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. livecrootp_storage   
        m_livecrootp_xfer_to_fire           =>    veg_pf%m_livecrootp_xfer_to_fire           , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. livecrootp_xfer      
        m_deadcrootp_to_fire                =>    veg_pf%m_deadcrootp_to_fire                , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. deadcrootp	       
        m_deadcrootp_storage_to_fire        =>    veg_pf%m_deadcrootp_storage_to_fire        , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. deadcrootp_storage   
        m_deadcrootp_xfer_to_fire           =>    veg_pf%m_deadcrootp_xfer_to_fire           , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. deadcrootp_xfer      
        m_retransp_to_fire                  =>    veg_pf%m_retransp_to_fire                  , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. retransp             
        m_ppool_to_fire                     =>    veg_pf%m_ppool_to_fire                     , & ! Input:  [real(r8) (:)     ]  (gP/m2/s) P emis. ppool     
        
        m_leafc_to_litter_fire              =>    veg_cf%m_leafc_to_litter_fire                , & ! Output: [real(r8) (:)     ]                                                    
        m_leafc_storage_to_litter_fire      =>    veg_cf%m_leafc_storage_to_litter_fire        , & ! Output: [real(r8) (:)     ]                                                    
        m_leafc_xfer_to_litter_fire         =>    veg_cf%m_leafc_xfer_to_litter_fire           , & ! Output: [real(r8) (:)     ]                                                    
        m_livestemc_to_litter_fire          =>    veg_cf%m_livestemc_to_litter_fire            , & ! Output: [real(r8) (:)     ]                                                    
        m_livestemc_storage_to_litter_fire  =>    veg_cf%m_livestemc_storage_to_litter_fire    , & ! Output: [real(r8) (:)     ]                                                    
        m_livestemc_xfer_to_litter_fire     =>    veg_cf%m_livestemc_xfer_to_litter_fire       , & ! Output: [real(r8) (:)     ]                                                    
        m_livestemc_to_deadstemc_fire       =>    veg_cf%m_livestemc_to_deadstemc_fire         , & ! Output: [real(r8) (:)     ]                                                    
        m_deadstemc_to_litter_fire          =>    veg_cf%m_deadstemc_to_litter_fire            , & ! Output: [real(r8) (:)     ]                                                    
        m_deadstemc_storage_to_litter_fire  =>    veg_cf%m_deadstemc_storage_to_litter_fire    , & ! Output: [real(r8) (:)     ]                                                    
        m_deadstemc_xfer_to_litter_fire     =>    veg_cf%m_deadstemc_xfer_to_litter_fire       , & ! Output: [real(r8) (:)     ]                                                    
        m_frootc_to_litter_fire             =>    veg_cf%m_frootc_to_litter_fire               , & ! Output: [real(r8) (:)     ]                                                    
        m_frootc_storage_to_litter_fire     =>    veg_cf%m_frootc_storage_to_litter_fire       , & ! Output: [real(r8) (:)     ]                                                    
        m_frootc_xfer_to_litter_fire        =>    veg_cf%m_frootc_xfer_to_litter_fire          , & ! Output: [real(r8) (:)     ]                                                    
        m_livecrootc_to_litter_fire         =>    veg_cf%m_livecrootc_to_litter_fire           , & ! Output: [real(r8) (:)     ]                                                    
        m_livecrootc_storage_to_litter_fire =>    veg_cf%m_livecrootc_storage_to_litter_fire   , & ! Output: [real(r8) (:)     ]                                                    
        m_livecrootc_xfer_to_litter_fire    =>    veg_cf%m_livecrootc_xfer_to_litter_fire      , & ! Output: [real(r8) (:)     ]                                                    
        m_livecrootc_to_deadcrootc_fire     =>    veg_cf%m_livecrootc_to_deadcrootc_fire       , & ! Output: [real(r8) (:)     ]                                                    
        m_deadcrootc_to_litter_fire         =>    veg_cf%m_deadcrootc_to_litter_fire           , & ! Output: [real(r8) (:)     ]                                                    
        m_deadcrootc_storage_to_litter_fire =>    veg_cf%m_deadcrootc_storage_to_litter_fire   , & ! Output: [real(r8) (:)     ]                                                    
        m_deadcrootc_xfer_to_litter_fire    =>    veg_cf%m_deadcrootc_xfer_to_litter_fire      , & ! Output: [real(r8) (:)     ]                                                   
        m_gresp_storage_to_litter_fire      =>    veg_cf%m_gresp_storage_to_litter_fire        , & ! Output: [real(r8) (:)     ]                                                    
        m_gresp_xfer_to_litter_fire         =>    veg_cf%m_gresp_xfer_to_litter_fire           , & ! Output: [real(r8) (:)     ]                             
        m_cpool_to_litter_fire              =>    veg_cf%m_cpool_to_litter_fire                , & ! Output: [real(r8) (:)     ]                       
        m_decomp_cpools_to_fire_vr          =>    col_cf%m_decomp_cpools_to_fire_vr              , & ! Output: [real(r8) (:,:,:) ]  (gC/m3/s) VR decomp. C fire loss
        m_c_to_litr_met_fire                =>    col_cf%m_c_to_litr_met_fire                    , & ! Output: [real(r8) (:,:)   ]                                                  
        m_c_to_litr_cel_fire                =>    col_cf%m_c_to_litr_cel_fire                    , & ! Output: [real(r8) (:,:)   ]                                                  
        m_c_to_litr_lig_fire                =>    col_cf%m_c_to_litr_lig_fire                      & ! Output: [real(r8) (:,:)   ]

        )

     m_leafn_to_litter_fire              =>    veg_nf%m_leafn_to_litter_fire              
     m_leafn_storage_to_litter_fire      =>    veg_nf%m_leafn_storage_to_litter_fire      
     m_leafn_xfer_to_litter_fire         =>    veg_nf%m_leafn_xfer_to_litter_fire         
     m_livestemn_to_litter_fire          =>    veg_nf%m_livestemn_to_litter_fire          
     m_livestemn_storage_to_litter_fire  =>    veg_nf%m_livestemn_storage_to_litter_fire  
     m_livestemn_xfer_to_litter_fire     =>    veg_nf%m_livestemn_xfer_to_litter_fire     
     m_livestemn_to_deadstemn_fire       =>    veg_nf%m_livestemn_to_deadstemn_fire       
     m_deadstemn_to_litter_fire          =>    veg_nf%m_deadstemn_to_litter_fire          
     m_deadstemn_storage_to_litter_fire  =>    veg_nf%m_deadstemn_storage_to_litter_fire  
     m_deadstemn_xfer_to_litter_fire     =>    veg_nf%m_deadstemn_xfer_to_litter_fire     
     m_frootn_to_litter_fire             =>    veg_nf%m_frootn_to_litter_fire             
     m_frootn_storage_to_litter_fire     =>    veg_nf%m_frootn_storage_to_litter_fire     
     m_frootn_xfer_to_litter_fire        =>    veg_nf%m_frootn_xfer_to_litter_fire        
     m_livecrootn_to_litter_fire         =>    veg_nf%m_livecrootn_to_litter_fire         
     m_livecrootn_storage_to_litter_fire =>    veg_nf%m_livecrootn_storage_to_litter_fire 
     m_livecrootn_xfer_to_litter_fire    =>    veg_nf%m_livecrootn_xfer_to_litter_fire    
     m_livecrootn_to_deadcrootn_fire     =>    veg_nf%m_livecrootn_to_deadcrootn_fire     
     m_deadcrootn_to_litter_fire         =>    veg_nf%m_deadcrootn_to_litter_fire         
     m_deadcrootn_storage_to_litter_fire =>    veg_nf%m_deadcrootn_storage_to_litter_fire 
     m_deadcrootn_xfer_to_litter_fire    =>    veg_nf%m_deadcrootn_xfer_to_litter_fire    
     m_retransn_to_litter_fire           =>    veg_nf%m_retransn_to_litter_fire           
     m_npool_to_litter_fire              =>    veg_nf%m_npool_to_litter_fire
     m_decomp_npools_to_fire_vr          =>    col_nf%m_decomp_npools_to_fire_vr            
     m_n_to_litr_met_fire                =>    col_nf%m_n_to_litr_met_fire                  
     m_n_to_litr_cel_fire                =>    col_nf%m_n_to_litr_cel_fire                  
     m_n_to_litr_lig_fire                =>    col_nf%m_n_to_litr_lig_fire                  


     m_leafp_to_litter_fire              =>    veg_pf%m_leafp_to_litter_fire              
     m_leafp_storage_to_litter_fire      =>    veg_pf%m_leafp_storage_to_litter_fire      
     m_leafp_xfer_to_litter_fire         =>    veg_pf%m_leafp_xfer_to_litter_fire         
     m_livestemp_to_litter_fire          =>    veg_pf%m_livestemp_to_litter_fire          
     m_livestemp_storage_to_litter_fire  =>    veg_pf%m_livestemp_storage_to_litter_fire  
     m_livestemp_xfer_to_litter_fire     =>    veg_pf%m_livestemp_xfer_to_litter_fire     
     m_livestemp_to_deadstemp_fire       =>    veg_pf%m_livestemp_to_deadstemp_fire       
     m_deadstemp_to_litter_fire          =>    veg_pf%m_deadstemp_to_litter_fire          
     m_deadstemp_storage_to_litter_fire  =>    veg_pf%m_deadstemp_storage_to_litter_fire  
     m_deadstemp_xfer_to_litter_fire     =>    veg_pf%m_deadstemp_xfer_to_litter_fire     
     m_frootp_to_litter_fire             =>    veg_pf%m_frootp_to_litter_fire             
     m_frootp_storage_to_litter_fire     =>    veg_pf%m_frootp_storage_to_litter_fire     
     m_frootp_xfer_to_litter_fire        =>    veg_pf%m_frootp_xfer_to_litter_fire        
     m_livecrootp_to_litter_fire         =>    veg_pf%m_livecrootp_to_litter_fire         
     m_livecrootp_storage_to_litter_fire =>    veg_pf%m_livecrootp_storage_to_litter_fire 
     m_livecrootp_xfer_to_litter_fire    =>    veg_pf%m_livecrootp_xfer_to_litter_fire    
     m_livecrootp_to_deadcrootp_fire     =>    veg_pf%m_livecrootp_to_deadcrootp_fire     
     m_deadcrootp_to_litter_fire         =>    veg_pf%m_deadcrootp_to_litter_fire         
     m_deadcrootp_storage_to_litter_fire =>    veg_pf%m_deadcrootp_storage_to_litter_fire 
     m_deadcrootp_xfer_to_litter_fire    =>    veg_pf%m_deadcrootp_xfer_to_litter_fire    
     m_retransp_to_litter_fire           =>    veg_pf%m_retransp_to_litter_fire           
     m_ppool_to_litter_fire              =>    veg_pf%m_ppool_to_litter_fire

     m_decomp_ppools_to_fire_vr          =>    col_pf%m_decomp_ppools_to_fire_vr            
     m_p_to_litr_met_fire                =>    col_pf%m_p_to_litr_met_fire                  
     m_p_to_litr_cel_fire                =>    col_pf%m_p_to_litr_cel_fire                  
     m_p_to_litr_lig_fire                =>    col_pf%m_p_to_litr_lig_fire                  

     transient_landcover = run_has_transient_landcover()

     ! Get model step size
     ! calculate burned area fraction per sec
     dt = real( get_step_size(), r8 )

     dayspyr = get_days_per_year()
     !
     ! patch loop
     !
     do fp = 1,num_soilp
        p = filter_soilp(fp)
        c = veg_pp%column(p)

        if( veg_pp%itype(p) < nc3crop .and. cropf_col(c) < 1.0_r8)then
           ! For non-crop (bare-soil and natural vegetation)
           if (transient_landcover) then    !true when landuse data is used
              f = (fbac(c)-baf_crop(c))/(1.0_r8-cropf_col(c))
           else
              f = (farea_burned(c)-baf_crop(c))/(1.0_r8-cropf_col(c))
           end if
        else
           ! For crops
           if(cropf_col(c) > 0._r8)then
             f = baf_crop(c) /cropf_col(c)
           else
             f = 0._r8
           end if
        end if

        ! apply this rate to the pft state variables to get flux rates
        ! biomass burning
        ! carbon fluxes

        m_veg = 1.0_r8
        if (spinup_state == 1) m_veg = spinup_mortality_factor
        m_leafc_to_fire(p)               =  leafc(p)              * f * cc_leaf(veg_pp%itype(p))
        m_leafc_storage_to_fire(p)       =  leafc_storage(p)      * f * cc_other(veg_pp%itype(p))
        m_leafc_xfer_to_fire(p)          =  leafc_xfer(p)         * f * cc_other(veg_pp%itype(p))
        m_livestemc_to_fire(p)           =  livestemc(p)          * f * cc_lstem(veg_pp%itype(p))
        m_livestemc_storage_to_fire(p)   =  livestemc_storage(p)  * f * cc_other(veg_pp%itype(p))
        m_livestemc_xfer_to_fire(p)      =  livestemc_xfer(p)     * f * cc_other(veg_pp%itype(p))
        m_deadstemc_to_fire(p)           =  deadstemc(p)          * m_veg * f * cc_dstem(veg_pp%itype(p))
        m_deadstemc_storage_to_fire(p)   =  deadstemc_storage(p)  * f * cc_other(veg_pp%itype(p))
        m_deadstemc_xfer_to_fire(p)      =  deadstemc_xfer(p)     * f * cc_other(veg_pp%itype(p))
        m_frootc_to_fire(p)              =  frootc(p)             * f * 0._r8
        m_frootc_storage_to_fire(p)      =  frootc_storage(p)     * f * cc_other(veg_pp%itype(p)) 
        m_frootc_xfer_to_fire(p)         =  frootc_xfer(p)        * f * cc_other(veg_pp%itype(p))
        m_livecrootc_to_fire(p)          =  livecrootc(p)         * f * 0._r8
        m_livecrootc_storage_to_fire(p)  =  livecrootc_storage(p) * f * cc_other(veg_pp%itype(p)) 
        m_livecrootc_xfer_to_fire(p)     =  livecrootc_xfer(p)    * f * cc_other(veg_pp%itype(p)) 
        m_deadcrootc_to_fire(p)          =  deadcrootc(p)         * m_veg * f * 0._r8
        m_deadcrootc_storage_to_fire(p)  =  deadcrootc_storage(p) * f*  cc_other(veg_pp%itype(p)) 
        m_deadcrootc_xfer_to_fire(p)     =  deadcrootc_xfer(p)    * f * cc_other(veg_pp%itype(p)) 
        m_gresp_storage_to_fire(p)       =  gresp_storage(p)      * f * cc_other(veg_pp%itype(p))
        m_gresp_xfer_to_fire(p)          =  gresp_xfer(p)         * f * cc_other(veg_pp%itype(p))
        m_cpool_to_fire(p)               =  cpool(p)              * f * cc_other(veg_pp%itype(p))

        ! nitrogen fluxes
        m_leafn_to_fire(p)               =  leafn(p)              * f * cc_leaf(veg_pp%itype(p))
        m_leafn_storage_to_fire(p)       =  leafn_storage(p)      * f * cc_other(veg_pp%itype(p))
        m_leafn_xfer_to_fire(p)          =  leafn_xfer(p)         * f * cc_other(veg_pp%itype(p))
        m_livestemn_to_fire(p)           =  livestemn(p)          * f * cc_lstem(veg_pp%itype(p))
        m_livestemn_storage_to_fire(p)   =  livestemn_storage(p)  * f * cc_other(veg_pp%itype(p))
        m_livestemn_xfer_to_fire(p)      =  livestemn_xfer(p)     * f * cc_other(veg_pp%itype(p))
        m_deadstemn_to_fire(p)           =  deadstemn(p)          * m_veg * f * cc_dstem(veg_pp%itype(p))
        m_deadstemn_storage_to_fire(p)   =  deadstemn_storage(p)  * f * cc_other(veg_pp%itype(p))
        m_deadstemn_xfer_to_fire(p)      =  deadstemn_xfer(p)     * f * cc_other(veg_pp%itype(p))
        m_frootn_to_fire(p)              =  frootn(p)             * f * 0._r8
        m_frootn_storage_to_fire(p)      =  frootn_storage(p)     * f * cc_other(veg_pp%itype(p))
        m_frootn_xfer_to_fire(p)         =  frootn_xfer(p)        * f * cc_other(veg_pp%itype(p))
        m_livecrootn_to_fire(p)          =  livecrootn(p)         * f * 0._r8 
        m_livecrootn_storage_to_fire(p)  =  livecrootn_storage(p) * f * cc_other(veg_pp%itype(p)) 
        m_livecrootn_xfer_to_fire(p)     =  livecrootn_xfer(p)    * f * cc_other(veg_pp%itype(p))
        m_deadcrootn_to_fire(p)          =  deadcrootn(p)         * m_veg * f * 0._r8
        m_deadcrootn_xfer_to_fire(p)     =  deadcrootn_xfer(p)    * f * cc_other(veg_pp%itype(p)) 
        m_deadcrootn_storage_to_fire(p)  =  deadcrootn_storage(p) * f * cc_other(veg_pp%itype(p))
        m_retransn_to_fire(p)            =  retransn(p)           * f * cc_other(veg_pp%itype(p))
        m_npool_to_fire(p)               =  npool(p)              * f * cc_other(veg_pp%itype(p))

        ! phosphorus fluxes
        m_leafp_to_fire(p)               =  leafp(p)              * f * cc_leaf(veg_pp%itype(p))
        m_leafp_storage_to_fire(p)       =  leafp_storage(p)      * f * cc_other(veg_pp%itype(p))
        m_leafp_xfer_to_fire(p)          =  leafp_xfer(p)         * f * cc_other(veg_pp%itype(p))
        m_livestemp_to_fire(p)           =  livestemp(p)          * f * cc_lstem(veg_pp%itype(p))
        m_livestemp_storage_to_fire(p)   =  livestemp_storage(p)  * f * cc_other(veg_pp%itype(p))
        m_livestemp_xfer_to_fire(p)      =  livestemp_xfer(p)     * f * cc_other(veg_pp%itype(p))
        !m_deadstemp_to_fire(p)           =  deadstemp(p)          * m_veg * f * cc_dstem(veg_pp%itype(p))
        m_deadstemp_to_fire(p)           =  deadstemp(p)          * f * cc_dstem(veg_pp%itype(p))
        m_deadstemp_storage_to_fire(p)   =  deadstemp_storage(p)  * f * cc_other(veg_pp%itype(p))
        m_deadstemp_xfer_to_fire(p)      =  deadstemp_xfer(p)     * f * cc_other(veg_pp%itype(p))
        m_frootp_to_fire(p)              =  frootp(p)             * f * 0._r8
        m_frootp_storage_to_fire(p)      =  frootp_storage(p)     * f * cc_other(veg_pp%itype(p))
        m_frootp_xfer_to_fire(p)         =  frootp_xfer(p)        * f * cc_other(veg_pp%itype(p))
        m_livecrootp_to_fire(p)          =  livecrootp(p)         * f * 0._r8 
        m_livecrootp_storage_to_fire(p)  =  livecrootp_storage(p) * f * cc_other(veg_pp%itype(p)) 
        m_livecrootp_xfer_to_fire(p)     =  livecrootp_xfer(p)    * f * cc_other(veg_pp%itype(p))
        m_deadcrootp_to_fire(p)          =  deadcrootp(p)         * m_veg * f * 0._r8
        m_deadcrootp_xfer_to_fire(p)     =  deadcrootp_xfer(p)    * f * cc_other(veg_pp%itype(p)) 
        m_deadcrootp_storage_to_fire(p)  =  deadcrootp_storage(p) * f * cc_other(veg_pp%itype(p))
        m_retransp_to_fire(p)            =  retransp(p)           * f * cc_other(veg_pp%itype(p))
        m_ppool_to_fire(p)               =  ppool(p)              * f * cc_other(veg_pp%itype(p))


        ! mortality due to fire
        ! carbon pools
        m_leafc_to_litter_fire(p)                   =  leafc(p) * f * &
             (1._r8 - cc_leaf(veg_pp%itype(p))) * &
             fm_leaf(veg_pp%itype(p))
        m_leafc_storage_to_litter_fire(p)           =  leafc_storage(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_leafc_xfer_to_litter_fire(p)              =  leafc_xfer(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_livestemc_to_litter_fire(p)               =  livestemc(p) * f * &
             (1._r8 - cc_lstem(veg_pp%itype(p))) * &
             fm_droot(veg_pp%itype(p))    
        m_livestemc_storage_to_litter_fire(p)       =  livestemc_storage(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_livestemc_xfer_to_litter_fire(p)          =  livestemc_xfer(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p)) 
        m_livestemc_to_deadstemc_fire(p)            =  livestemc(p) * f * &
             (1._r8 - cc_lstem(veg_pp%itype(p))) * &
             (fm_lstem(veg_pp%itype(p))-fm_droot(veg_pp%itype(p)))
        m_deadstemc_to_litter_fire(p)               =  deadstemc(p) * m_veg * f * &
             (1._r8 - cc_dstem(veg_pp%itype(p))) * &
             fm_droot(veg_pp%itype(p))    
        m_deadstemc_storage_to_litter_fire(p)       =  deadstemc_storage(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_deadstemc_xfer_to_litter_fire(p)          =  deadstemc_xfer(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_frootc_to_litter_fire(p)                  =  frootc(p)             * f * &
             fm_root(veg_pp%itype(p))
        m_frootc_storage_to_litter_fire(p)          =  frootc_storage(p)     * f * &
             fm_other(veg_pp%itype(p))
        m_frootc_xfer_to_litter_fire(p)             =  frootc_xfer(p)        * f * &
             fm_other(veg_pp%itype(p))
        m_livecrootc_to_litter_fire(p)              =  livecrootc(p)         * f * &
             fm_droot(veg_pp%itype(p))
        m_livecrootc_storage_to_litter_fire(p)      =  livecrootc_storage(p) * f * &
             fm_other(veg_pp%itype(p)) 
        m_livecrootc_xfer_to_litter_fire(p)         =  livecrootc_xfer(p)    * f * &
             fm_other(veg_pp%itype(p)) 
        m_livecrootc_to_deadcrootc_fire(p)          =  livecrootc(p)         * f * &
             (fm_lroot(veg_pp%itype(p))-fm_droot(veg_pp%itype(p)))
        m_deadcrootc_to_litter_fire(p)              =  deadcrootc(p)   * m_veg * f * &
             fm_droot(veg_pp%itype(p))
        m_deadcrootc_storage_to_litter_fire(p)      =  deadcrootc_storage(p) * f * &
             fm_other(veg_pp%itype(p))
        m_deadcrootc_xfer_to_litter_fire(p)         =  deadcrootc_xfer(p)    * f * &
             fm_other(veg_pp%itype(p))      
        m_gresp_storage_to_litter_fire(p)           =  gresp_storage(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))  
        m_gresp_xfer_to_litter_fire(p)              =  gresp_xfer(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p)) 
        m_cpool_to_litter_fire(p)                   =  cpool(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))

        ! nitrogen pools    
        m_leafn_to_litter_fire(p)                  =  leafn(p) * f * &
             (1._r8 - cc_leaf(veg_pp%itype(p))) * &
             fm_leaf(veg_pp%itype(p))
        m_leafn_storage_to_litter_fire(p)          =  leafn_storage(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))  
        m_leafn_xfer_to_litter_fire(p)             =  leafn_xfer(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_livestemn_to_litter_fire(p)              =  livestemn(p) * f * &
             (1._r8 - cc_lstem(veg_pp%itype(p))) * &
             fm_droot(veg_pp%itype(p))
        m_livestemn_storage_to_litter_fire(p)      =  livestemn_storage(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))   
        m_livestemn_xfer_to_litter_fire(p)         =  livestemn_xfer(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_livestemn_to_deadstemn_fire(p)           =  livestemn(p) * f * &
             (1._r8 - cc_lstem(veg_pp%itype(p))) * &
             (fm_lstem(veg_pp%itype(p))-fm_droot(veg_pp%itype(p)))
        m_frootn_to_litter_fire(p)                 =  frootn(p)             * f * &
             fm_root(veg_pp%itype(p))
        m_deadstemn_to_litter_fire(p)              =  deadstemn(p) * m_veg * f *  &
             (1._r8 - cc_dstem(veg_pp%itype(p))) * &
             fm_droot(veg_pp%itype(p))
        m_deadstemn_storage_to_litter_fire(p)       =  deadstemn_storage(p) * f *  &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_deadstemn_xfer_to_litter_fire(p)          =  deadstemn_xfer(p) * f *  &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_frootn_storage_to_litter_fire(p)         =  frootn_storage(p)     * f * &
             fm_other(veg_pp%itype(p))
        m_frootn_xfer_to_litter_fire(p)            =  frootn_xfer(p)        * f * &
             fm_other(veg_pp%itype(p))
        m_livecrootn_to_litter_fire(p)             =  livecrootn(p)         * f * &
             fm_droot(veg_pp%itype(p))
        m_livecrootn_storage_to_litter_fire(p)     =  livecrootn_storage(p) * f * &
             fm_other(veg_pp%itype(p))
        m_livecrootn_xfer_to_litter_fire(p)        =  livecrootn_xfer(p)    * f * &
             fm_other(veg_pp%itype(p)) 
        m_livecrootn_to_deadcrootn_fire(p)         =  livecrootn(p)         * f * &
             (fm_lroot(veg_pp%itype(p))-fm_droot(veg_pp%itype(p)))
        m_deadcrootn_to_litter_fire(p)             =  deadcrootn(p)    * m_veg * f * &
             fm_droot(veg_pp%itype(p))
        m_deadcrootn_storage_to_litter_fire(p)     =  deadcrootn_storage(p) * f * &
             fm_other(veg_pp%itype(p))
        m_deadcrootn_xfer_to_litter_fire(p)        =  deadcrootn_xfer(p)    * f * &
             fm_other(veg_pp%itype(p))
        m_retransn_to_litter_fire(p)               =  retransn(p)           * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p)) 
        m_npool_to_litter_fire(p)                  =  npool(p)              * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))

        ! phosphorus fluxes   
        m_leafp_to_litter_fire(p)                  =  leafp(p) * f * &
             (1._r8 - cc_leaf(veg_pp%itype(p))) * &
             fm_leaf(veg_pp%itype(p))
        m_leafp_storage_to_litter_fire(p)          =  leafp_storage(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))  
        m_leafp_xfer_to_litter_fire(p)             =  leafp_xfer(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_livestemp_to_litter_fire(p)              =  livestemp(p) * f * &
             (1._r8 - cc_lstem(veg_pp%itype(p))) * &
             fm_droot(veg_pp%itype(p))
        m_livestemp_storage_to_litter_fire(p)      =  livestemp_storage(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))   
        m_livestemp_xfer_to_litter_fire(p)         =  livestemp_xfer(p) * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_livestemp_to_deadstemp_fire(p)           =  livestemp(p) * f * &
             (1._r8 - cc_lstem(veg_pp%itype(p))) * &
             (fm_lstem(veg_pp%itype(p))-fm_droot(veg_pp%itype(p)))
        m_frootp_to_litter_fire(p)                 =  frootp(p)             * f * &
             fm_root(veg_pp%itype(p))
        m_deadstemp_to_litter_fire(p)               =  deadstemp(p) * m_veg * f *  &
             (1._r8 - cc_dstem(veg_pp%itype(p))) * &
             fm_droot(veg_pp%itype(p))
        m_deadstemp_storage_to_litter_fire(p)       =  deadstemp_storage(p) * f *  &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_deadstemp_xfer_to_litter_fire(p)          =  deadstemp_xfer(p) * f *  &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))
        m_frootp_storage_to_litter_fire(p)         =  frootp_storage(p)     * f * &
             fm_other(veg_pp%itype(p))
        m_frootp_xfer_to_litter_fire(p)            =  frootp_xfer(p)        * f * &
             fm_other(veg_pp%itype(p))
        m_livecrootp_to_litter_fire(p)             =  livecrootp(p)         * f * &
             fm_droot(veg_pp%itype(p))
        m_livecrootp_storage_to_litter_fire(p)     =  livecrootp_storage(p) * f * &
             fm_other(veg_pp%itype(p))
        m_livecrootp_xfer_to_litter_fire(p)        =  livecrootp_xfer(p)    * f * &
             fm_other(veg_pp%itype(p)) 
        m_livecrootp_to_deadcrootp_fire(p)         =  livecrootp(p)         * f * &
             (fm_lroot(veg_pp%itype(p))-fm_droot(veg_pp%itype(p)))
        m_deadcrootp_to_litter_fire(p)             =  deadcrootp(p)         * f * &
             fm_droot(veg_pp%itype(p)) ! * m_veg
        m_deadcrootp_storage_to_litter_fire(p)     =  deadcrootp_storage(p) * f * &
             fm_other(veg_pp%itype(p))
        m_deadcrootp_xfer_to_litter_fire(p)        =  deadcrootp_xfer(p)    * f * &
             fm_other(veg_pp%itype(p))
        m_retransp_to_litter_fire(p)               =  retransp(p)           * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p)) 
        m_ppool_to_litter_fire(p)                  =  ppool(p)              * f * &
             (1._r8 - cc_other(veg_pp%itype(p))) * &
             fm_other(veg_pp%itype(p))

     end do  ! end of patches loop  

     ! fire-induced transfer of carbon and nitrogen pools to litter and cwd
     ! add phosphorus transfer fluxes -X.YANG 

     do j = 1,nlevdecomp
        do pi = 1,max_patch_per_col
           do fc = 1,num_soilc
              c = filter_soilc(fc)
              if (pi <=  col_pp%npfts(c)) then
                 p = col_pp%pfti(c) + pi - 1
                 if ( veg_pp%active(p) ) then

                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_deadstemc_to_litter_fire(p) * veg_pp%wtcol(p) * stem_prof(p,j)
                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_deadcrootc_to_litter_fire(p) * veg_pp%wtcol(p) * croot_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_deadstemn_to_litter_fire(p) * veg_pp%wtcol(p) * stem_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_deadcrootn_to_litter_fire(p) * veg_pp%wtcol(p) * croot_prof(p,j)
                    ! add phosphorus
                    fire_mortality_p_to_cwdp(c,j) = fire_mortality_p_to_cwdp(c,j) + &
                         m_deadstemp_to_litter_fire(p) * veg_pp%wtcol(p) * stem_prof(p,j)
                    fire_mortality_p_to_cwdp(c,j) = fire_mortality_p_to_cwdp(c,j) + &
                         m_deadcrootp_to_litter_fire(p) * veg_pp%wtcol(p) * croot_prof(p,j)


                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_livestemc_to_litter_fire(p) * veg_pp%wtcol(p) * stem_prof(p,j)
                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_livecrootc_to_litter_fire(p) * veg_pp%wtcol(p) * croot_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_livestemn_to_litter_fire(p) * veg_pp%wtcol(p) * stem_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_livecrootn_to_litter_fire(p) * veg_pp%wtcol(p) * croot_prof(p,j)
                    ! add phosphorus
                    fire_mortality_p_to_cwdp(c,j) = fire_mortality_p_to_cwdp(c,j) + &
                         m_livestemp_to_litter_fire(p) * veg_pp%wtcol(p) * stem_prof(p,j)
                    fire_mortality_p_to_cwdp(c,j) = fire_mortality_p_to_cwdp(c,j) + &
                         m_livecrootp_to_litter_fire(p) * veg_pp%wtcol(p) * croot_prof(p,j)


                    m_c_to_litr_met_fire(c,j)=m_c_to_litr_met_fire(c,j) + &
                         ((m_leafc_to_litter_fire(p)*lf_flab(veg_pp%itype(p)) &
                         +m_leafc_storage_to_litter_fire(p) + &
                         m_leafc_xfer_to_litter_fire(p) + m_cpool_to_litter_fire(p) + &
                         m_gresp_storage_to_litter_fire(p) &
                         +m_gresp_xfer_to_litter_fire(p))*leaf_prof(p,j) + &
                         (m_frootc_to_litter_fire(p)*fr_flab(veg_pp%itype(p)) &
                         +m_frootc_storage_to_litter_fire(p) + &
                         m_frootc_xfer_to_litter_fire(p))*froot_prof(p,j) &
                         +(m_livestemc_storage_to_litter_fire(p) + &
                         m_livestemc_xfer_to_litter_fire(p) &
                         +m_deadstemc_storage_to_litter_fire(p) + &
                         m_deadstemc_xfer_to_litter_fire(p))* stem_prof(p,j)&
                         +(m_livecrootc_storage_to_litter_fire(p) + &
                         m_livecrootc_xfer_to_litter_fire(p) &
                         +m_deadcrootc_storage_to_litter_fire(p) + &
                         m_deadcrootc_xfer_to_litter_fire(p))* croot_prof(p,j))* veg_pp%wtcol(p)    
                    m_c_to_litr_cel_fire(c,j)=m_c_to_litr_cel_fire(c,j) + &
                         (m_leafc_to_litter_fire(p)*lf_fcel(veg_pp%itype(p))*leaf_prof(p,j) + &
                         m_frootc_to_litter_fire(p)*fr_fcel(veg_pp%itype(p))*froot_prof(p,j))* veg_pp%wtcol(p) 
                    m_c_to_litr_lig_fire(c,j)=m_c_to_litr_lig_fire(c,j) + &
                         (m_leafc_to_litter_fire(p)*lf_flig(veg_pp%itype(p))*leaf_prof(p,j) + &
                         m_frootc_to_litter_fire(p)*fr_flig(veg_pp%itype(p))*froot_prof(p,j))* veg_pp%wtcol(p)  

                    m_n_to_litr_met_fire(c,j)=m_n_to_litr_met_fire(c,j) + &
                         ((m_leafn_to_litter_fire(p)*lf_flab(veg_pp%itype(p)) &
                         +m_leafn_storage_to_litter_fire(p) + m_npool_to_litter_fire(p) + &
                         m_leafn_xfer_to_litter_fire(p)+m_retransn_to_litter_fire(p)) &
                         *leaf_prof(p,j) +(m_frootn_to_litter_fire(p)*fr_flab(veg_pp%itype(p)) &
                         +m_frootn_storage_to_litter_fire(p) + &
                         m_frootn_xfer_to_litter_fire(p))*froot_prof(p,j) &
                         +(m_livestemn_storage_to_litter_fire(p) + &
                         m_livestemn_xfer_to_litter_fire(p) &
                         +m_deadstemn_storage_to_litter_fire(p) + &
                         m_deadstemn_xfer_to_litter_fire(p))* stem_prof(p,j)&
                         +(m_livecrootn_storage_to_litter_fire(p) + &
                         m_livecrootn_xfer_to_litter_fire(p) &
                         +m_deadcrootn_storage_to_litter_fire(p) + &
                         m_deadcrootn_xfer_to_litter_fire(p))* croot_prof(p,j))* veg_pp%wtcol(p)    
                    m_n_to_litr_cel_fire(c,j)=m_n_to_litr_cel_fire(c,j) + &
                         (m_leafn_to_litter_fire(p)*lf_fcel(veg_pp%itype(p))*leaf_prof(p,j) + &
                         m_frootn_to_litter_fire(p)*fr_fcel(veg_pp%itype(p))*froot_prof(p,j))* veg_pp%wtcol(p) 
                    m_n_to_litr_lig_fire(c,j)=m_n_to_litr_lig_fire(c,j) + &
                         (m_leafn_to_litter_fire(p)*lf_flig(veg_pp%itype(p))*leaf_prof(p,j) + &
                         m_frootn_to_litter_fire(p)*fr_flig(veg_pp%itype(p))*froot_prof(p,j))* veg_pp%wtcol(p) 

                    ! add phosphorus
                    m_p_to_litr_met_fire(c,j)=m_p_to_litr_met_fire(c,j) + &
                         ((m_leafp_to_litter_fire(p)*lf_flab(veg_pp%itype(p)) &
                         +m_leafp_storage_to_litter_fire(p) + m_ppool_to_litter_fire(p) + &
                         m_leafp_xfer_to_litter_fire(p)+m_retransp_to_litter_fire(p)) &
                         *leaf_prof(p,j) +(m_frootp_to_litter_fire(p)*fr_flab(veg_pp%itype(p)) &
                         +m_frootp_storage_to_litter_fire(p) + &
                         m_frootp_xfer_to_litter_fire(p))*froot_prof(p,j) &
                         +(m_livestemp_storage_to_litter_fire(p) + &
                         m_livestemp_xfer_to_litter_fire(p) &
                         +m_deadstemp_storage_to_litter_fire(p) + &
                         m_deadstemp_xfer_to_litter_fire(p))* stem_prof(p,j)&
                         +(m_livecrootp_storage_to_litter_fire(p) + &
                         m_livecrootp_xfer_to_litter_fire(p) &
                         +m_deadcrootp_storage_to_litter_fire(p) + &
                         m_deadcrootp_xfer_to_litter_fire(p))* croot_prof(p,j))* veg_pp%wtcol(p)    
                    m_p_to_litr_cel_fire(c,j)=m_p_to_litr_cel_fire(c,j) + &
                         (m_leafp_to_litter_fire(p)*lf_fcel(veg_pp%itype(p))*leaf_prof(p,j) + &
                         m_frootp_to_litter_fire(p)*fr_fcel(veg_pp%itype(p))*froot_prof(p,j))* veg_pp%wtcol(p) 
                    m_p_to_litr_lig_fire(c,j)=m_p_to_litr_lig_fire(c,j) + &
                         (m_leafp_to_litter_fire(p)*lf_flig(veg_pp%itype(p))*leaf_prof(p,j) + &
                         m_frootp_to_litter_fire(p)*fr_flig(veg_pp%itype(p))*froot_prof(p,j))* veg_pp%wtcol(p) 

                 end if
              end if
           end do
        end do
     end do
     !
     ! vertically-resolved decomposing C/N fire loss   
     ! column loop
     ! add phosphorus 
     do fc = 1,num_soilc
        c = filter_soilc(fc)

        f = farea_burned(c) 

        ! change CC for litter from 0.4_r8 to 0.5_r8 and CC for CWD from 0.2_r8
        ! to 0.25_r8 according to Li et al.(2014) 
        do j = 1, nlevdecomp
           ! carbon fluxes
           do l = 1, ndecomp_pools
              if ( is_litter(l) ) then
                 m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * f * 0.5_r8 
              end if
              if ( is_cwd(l) ) then
                 m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * &
                      (f-baf_crop(c)) * 0.25_r8
                 call get_curr_date(kyr, kmo, kda, mcsec)
                 if (spinup_state == 1) then 
                   m_decomp_cpools_to_fire_vr(c,j,l) = m_decomp_cpools_to_fire_vr(c,j,l) * &
                     decomp_cascade_con%spinup_factor(l) 
                   if (kyr >= 40) m_decomp_cpools_to_fire_vr(c,j,l) = &
                     m_decomp_cpools_to_fire_vr(c,j,l) / cnstate_vars%scalaravg_col(c,j)
                 end if
              end if
           end do

           ! nitrogen fluxes
           do l = 1, ndecomp_pools
              if ( is_litter(l) ) then
                 m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * f * 0.5_r8
              end if
              if ( is_cwd(l) ) then
                 m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * &
                      (f-baf_crop(c)) * 0.25_r8
                 if (spinup_state == 1) then 
                   m_decomp_npools_to_fire_vr(c,j,l) = m_decomp_npools_to_fire_vr(c,j,l) * &
                     decomp_cascade_con%spinup_factor(l) 
                   if (kyr >= 40) m_decomp_npools_to_fire_vr(c,j,l) = &
                     m_decomp_npools_to_fire_vr(c,j,l) / cnstate_vars%scalaravg_col(c,j)
                 end if             
             end if
           end do

           ! phosphorus fluxes - loss due to fire
           do l = 1, ndecomp_pools
              if ( is_litter(l) ) then
                 m_decomp_ppools_to_fire_vr(c,j,l) = decomp_ppools_vr(c,j,l) * f * 0.5_r8
              end if
              if ( is_cwd(l) ) then
                 m_decomp_ppools_to_fire_vr(c,j,l) = decomp_ppools_vr(c,j,l) * &
                      (f-baf_crop(c)) * 0.25_r8
              end if
           end do

        end do
     end do  ! end of column loop

     ! carbon loss due to deforestation fires

     if (transient_landcover) then    !true when landuse data is used
        call get_curr_date (kyr, kmo, kda, mcsec)
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           lfc2(c)=0._r8
           if( .not. (kmo == 1 .and. kda == 1 .and. mcsec == 0) )then
              if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 .and. dtrotr_col(c) > 0._r8 .and. &
                   lfc(c) > 0._r8 .and. fbac1(c) == 0._r8) then
                 lfc2(c) = max(0._r8, min(lfc(c), (farea_burned(c)-baf_crop(c) - &
                      baf_peatf(c))/2.0*dt))/(dtrotr_col(c)*dayspyr*secspday/dt)/dt
                 lfc(c)  = lfc(c) - max(0._r8, min(lfc(c), (farea_burned(c)-baf_crop(c) - &
                      baf_peatf(c))*dt/2.0_r8))
              end if
           end if
        end do
     end if
     !
     ! Carbon loss due to peat fires
     !
     ! somc_fire is not connected to clm45 soil carbon pool, ie does not decrease
     ! soil carbon b/c clm45 soil carbon was very low in several peatland grids
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)
        if( grc_pp%latdeg(g)  <  borealat)then
           somc_fire(c)= totsomc(c)*baf_peatf(c)*6.0_r8/33.9_r8
        else
           somc_fire(c)= baf_peatf(c)*2.2e3_r8
        end if
     end do

     ! Fang Li has not added aerosol and trace gas emissions due to fire, yet
     ! They will be added here in proportion to the carbon emission
     ! Emission factors differ for various fire types

   end associate 

 end subroutine FireFluxes

 !-----------------------------------------------------------------------
 subroutine hdm_init( bounds )
   !
   ! !DESCRIPTION:
   ! Initialize data stream information for population density.
   !
   ! !USES:
   use elm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use ndepStreamMod    , only : clm_domain_mct
   use histFileMod      , only : hist_addfld1d
   !
   ! !ARGUMENTS:
   implicit none
   type(bounds_type), intent(in) :: bounds  
   !
   ! !LOCAL VARIABLES:
   integer            :: stream_year_first_popdens   ! first year in pop. dens. stream to use
   integer            :: stream_year_last_popdens    ! last year in pop. dens. stream to use
   integer            :: model_year_align_popdens    ! align stream_year_first_hdm with 
   integer            :: nu_nml                      ! unit for namelist file
   integer            :: nml_error                   ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm                     ! domain information 
   character(len=CL)  :: stream_fldFileName_popdens  ! population density streams filename
   character(len=CL)  :: popdensmapalgo = 'bilinear' ! mapping alogrithm for population density
   character(*), parameter :: subName = "('hdmdyn_init')"
   character(*), parameter :: F00 = "('(hdmdyn_init) ',4a)"
   !-----------------------------------------------------------------------

   namelist /popd_streams/          &
        stream_year_first_popdens,  &
        stream_year_last_popdens,   &
        model_year_align_popdens,   &
        popdensmapalgo,             &
        stream_fldFileName_popdens

   ! Allocate pop dens forcing data
   allocate( forc_hdm(bounds%begg:bounds%endg) )

   ! Default values for namelist
   stream_year_first_popdens  = 1       ! first year in stream to use
   stream_year_last_popdens   = 1       ! last  year in stream to use
   model_year_align_popdens   = 1       ! align stream_year_first_popdens with this model year
   stream_fldFileName_popdens = ' '

   ! Read popd_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'popd_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=popd_streams,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading popd_streams namelist'//errMsg(__FILE__, __LINE__))
         end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_popdens, mpicom)
   call shr_mpi_bcast(stream_year_last_popdens, mpicom)
   call shr_mpi_bcast(model_year_align_popdens, mpicom)
   call shr_mpi_bcast(stream_fldFileName_popdens, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'popdens_streams settings:'
      write(iulog,*) '  stream_year_first_popdens  = ',stream_year_first_popdens  
      write(iulog,*) '  stream_year_last_popdens   = ',stream_year_last_popdens   
      write(iulog,*) '  model_year_align_popdens   = ',model_year_align_popdens   
      write(iulog,*) '  stream_fldFileName_popdens = ',stream_fldFileName_popdens
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(sdat_hdm,name="clmhdm",     &
        pio_subsystem=pio_subsystem,                   & 
        pio_iotype=shr_pio_getiotype(inst_name),       &
        mpicom=mpicom, compid=comp_id,                 &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,        &
        nxg=ldomain%ni, nyg=ldomain%nj,                &
        yearFirst=stream_year_first_popdens,           &
        yearLast=stream_year_last_popdens,             &
        yearAlign=model_year_align_popdens,            &
        offset=0,                                      &
        domFilePath='',                                &
        domFileName=trim(stream_fldFileName_popdens),  &
        domTvarName='time',                            &
        domXvarName='lon' ,                            &
        domYvarName='lat' ,                            &  
        domAreaName='area',                            &
        domMaskName='mask',                            &
        filePath='',                                   &
        filename=(/trim(stream_fldFileName_popdens)/) , &
        fldListFile='hdm',                             &
        fldListModel='hdm',                            &
        fillalgo='none',                               &
        mapalgo=popdensmapalgo,                        &
        calendar=get_calendar(),                       &
        tintalgo='nearest',                            &
        taxmode='extend'                           )

   if (masterproc) then
      call shr_strdata_print(sdat_hdm,'population density data')
   endif

   ! Add history fields
   call hist_addfld1d (fname='HDM', units='counts/km^2',      &
         avgflag='A', long_name='human population density',   &
         ptr_lnd=forc_hdm, default='inactive')

end subroutine hdm_init
  
!-----------------------------------------------------------------------
subroutine hdm_interp(bounds)
  !
  ! !DESCRIPTION:
  ! Interpolate data stream information for population density.
  !
  ! !USES:
  use clm_time_manager, only : get_curr_date
  !
  ! !ARGUMENTS:
  type(bounds_type), intent(in) :: bounds  
  !
  ! !LOCAL VARIABLES:
  integer :: g, ig
  integer :: year    ! year (0, ...) for nstep+1
  integer :: mon     ! month (1, ..., 12) for nstep+1
  integer :: day     ! day of month (1, ..., 31) for nstep+1
  integer :: sec     ! seconds into current date for nstep+1
  integer :: mcdate  ! Current model date (yyyymmdd)
  !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(sdat_hdm, mcdate, sec, mpicom, 'hdmdyn')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      forc_hdm(g) = sdat_hdm%avs(1)%rAttr(1,ig)
   end do
   
end subroutine hdm_interp

!-----------------------------------------------------------------------
subroutine lnfm_init( bounds )
  !
  ! !DESCRIPTION:
  !
  ! Initialize data stream information for Lightning.
  !
  ! !USES:
  use elm_varctl       , only : inst_name
  use clm_time_manager , only : get_calendar
  use ncdio_pio        , only : pio_subsystem
  use shr_pio_mod      , only : shr_pio_getiotype
  use clm_nlUtilsMod   , only : find_nlgroup_name
  use ndepStreamMod    , only : clm_domain_mct
  use histFileMod      , only : hist_addfld1d
  !
  ! !ARGUMENTS:
  implicit none
  type(bounds_type), intent(in) :: bounds  
  !
  ! !LOCAL VARIABLES:
  integer            :: stream_year_first_lightng  ! first year in Lightning stream to use
  integer            :: stream_year_last_lightng   ! last year in Lightning stream to use
  integer            :: model_year_align_lightng   ! align stream_year_first_lnfm with 
  integer            :: nu_nml                     ! unit for namelist file
  integer            :: nml_error                  ! namelist i/o error flag
  type(mct_ggrid)    :: dom_clm                    ! domain information 
  character(len=CL)  :: stream_fldFileName_lightng ! lightning stream filename to read
  character(len=CL)  :: lightngmapalgo = 'bilinear'! Mapping alogrithm
  character(*), parameter :: subName = "('lnfmdyn_init')"
  character(*), parameter :: F00 = "('(lnfmdyn_init) ',4a)"
  !-----------------------------------------------------------------------

   namelist /light_streams/         &
        stream_year_first_lightng,  &
        stream_year_last_lightng,   &
        model_year_align_lightng,   &
        lightngmapalgo,             &
        stream_fldFileName_lightng

   ! Allocate lightning forcing data
   allocate( forc_lnfm(bounds%begg:bounds%endg) )

   ! Default values for namelist
    stream_year_first_lightng  = 1      ! first year in stream to use
    stream_year_last_lightng   = 1      ! last  year in stream to use
    model_year_align_lightng   = 1      ! align stream_year_first_lnfm with this model year
    stream_fldFileName_lightng = ' '

   ! Read light_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'light_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=light_streams,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading light_streams namelist'//errMsg(__FILE__, __LINE__))
         end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_lightng, mpicom)
   call shr_mpi_bcast(stream_year_last_lightng, mpicom)
   call shr_mpi_bcast(model_year_align_lightng, mpicom)
   call shr_mpi_bcast(stream_fldFileName_lightng, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'light_stream settings:'
      write(iulog,*) '  stream_year_first_lightng  = ',stream_year_first_lightng  
      write(iulog,*) '  stream_year_last_lightng   = ',stream_year_last_lightng   
      write(iulog,*) '  model_year_align_lightng   = ',model_year_align_lightng   
      write(iulog,*) '  stream_fldFileName_lightng = ',stream_fldFileName_lightng
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(sdat_lnfm,name="clmlnfm",  &
        pio_subsystem=pio_subsystem,                  & 
        pio_iotype=shr_pio_getiotype(inst_name),      &
        mpicom=mpicom, compid=comp_id,                &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,       &
        nxg=ldomain%ni, nyg=ldomain%nj,               &
        yearFirst=stream_year_first_lightng,          &
        yearLast=stream_year_last_lightng,            &
        yearAlign=model_year_align_lightng,           &
        offset=0,                                     &
        domFilePath='',                               &
        domFileName=trim(stream_fldFileName_lightng), &
        domTvarName='time',                           &
        domXvarName='lon' ,                           &
        domYvarName='lat' ,                           &  
        domAreaName='area',                           &
        domMaskName='mask',                           &
        filePath='',                                  &
        filename=(/trim(stream_fldFileName_lightng)/),&
        fldListFile='lnfm',                           &
        fldListModel='lnfm',                          &
        fillalgo='none',                              &
        mapalgo=lightngmapalgo,                       &
        calendar=get_calendar(),                      &
        taxmode='cycle'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_lnfm,'Lightning data')
   endif

   ! Add history fields
   call hist_addfld1d (fname='LNFM', units='counts/km^2/hr',  &
         avgflag='A', long_name='Lightning frequency',        &
         ptr_lnd=forc_lnfm, default='inactive')

end subroutine lnfm_init
  
!-----------------------------------------------------------------------
subroutine lnfm_interp(bounds )
  !
  ! !DESCRIPTION:
  ! Interpolate data stream information for Lightning.
  !
  ! !USES:
  use clm_time_manager, only : get_curr_date
  !
  ! !ARGUMENTS:
  type(bounds_type), intent(in) :: bounds  
  !
  ! !LOCAL VARIABLES:
  integer :: g, ig
  integer :: year    ! year (0, ...) for nstep+1
  integer :: mon     ! month (1, ..., 12) for nstep+1
  integer :: day     ! day of month (1, ..., 31) for nstep+1
  integer :: sec     ! seconds into current date for nstep+1
  integer :: mcdate  ! Current model date (yyyymmdd)
  !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(sdat_lnfm, mcdate, sec, mpicom, 'lnfmdyn')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      forc_lnfm(g) = sdat_lnfm%avs(1)%rAttr(1,ig)
   end do
   
end subroutine lnfm_interp

end module FireMod
