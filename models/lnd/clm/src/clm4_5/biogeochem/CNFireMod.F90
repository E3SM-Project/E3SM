module CNFireMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for fire dynamics 
  ! created in Nov, 2012  and revised in Apr, 2013 by F. Li and S. Levis
  ! based on Li et al. (2012a,b; 2013)"
  ! Fire-related parameters were calibrated or tuned in Apr, 2013 based on the 
  ! 20th Century transient simulations at f19_g16 with (newfire05_clm45sci15_clm4_0_58) 
  ! a CLM4.5 version, Qian et al. (2006) atmospheric forcing, and
  ! climatological lightning data.
  !
  ! !USES:
  use shr_kind_mod   , only: r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_const_mod  , only: SHR_CONST_PI,SHR_CONST_TKFRZ
  use shr_strdata_mod, only: shr_strdata_type, shr_strdata_create, shr_strdata_print
  use shr_strdata_mod, only: shr_strdata_advance
  use shr_log_mod    , only: errMsg => shr_log_errMsg
  use pft2colMod     , only: p2c
  use clm_varctl     , only: iulog
  use clm_varpar     , only: nlevdecomp, ndecomp_pools
  use clm_varcon     , only: dzsoi_decomp
  use pftvarcon      , only: fsr_pft, fd_pft, noveg
  use spmdMod        , only: masterproc, mpicom, comp_id
  use fileutils      , only: getavu, relavu
  use controlMod     , only: NLFilename
  use decompMod      , only: gsmap_lnd_gdc2glo
  use domainMod      , only: ldomain
  use abortutils     , only: endrun
  use decompMod      , only: bounds_type
  use mct_mod
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNFireInit    ! Initialization of CNFire
  public :: CNFireInterp  ! Interpolate fire data
  public :: CNFireArea    ! Calculate fire area
  public :: CNFireFluxes  ! Calculate fire fluxes
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
  subroutine CNFireInit( bounds )
    !
    ! !DESCRIPTION:
    ! Initialize CN Fire module
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !-----------------------------------------------------------------------

    call hdm_init(bounds)
    call hdm_interp(bounds)
    call lnfm_init(bounds)
    call lnfm_interp(bounds)

  end subroutine CNFireInit


  !-----------------------------------------------------------------------
  subroutine CNFireInterp(bounds)
    !
    ! !DESCRIPTION:
    ! Interpolate CN Fire datasets
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !-----------------------------------------------------------------------

    call hdm_interp(bounds)
    call lnfm_interp(bounds)

  end subroutine CNFireInterp

  !-----------------------------------------------------------------------
  subroutine CNFireArea (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! Computes column-level burned area in each timestep
    !
    ! !USES:
    use clmtype
    use clm_time_manager, only: get_step_size, get_days_per_year, get_curr_date, get_nstep
    use clm_varpar      , only: max_pft_per_col
    use clm_varcon      , only: secspday
    use clm_atmlnd      , only: clm_a2l, a2l_downscaled_col
    use shr_infnan_mod  , only: shr_infnan_isnan
    use clm_varctl      , only: fpftdyn, use_nofire
    use pftvarcon       , only: nc4_grass, nc3crop, ndllf_evr_tmp_tree, &
                                nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree, nbrdlf_evr_shrub
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:) ! filter for soil columns
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter  :: lfuel=110._r8   ! lower threshold of fuel mass (gC/m2) for ignition
    real(r8), parameter  :: ufuel=1050._r8  ! upper threshold of fuel mass(gC/m2) for ignition 
    real(r8), parameter  :: g0=0.05_r8      ! g(W) when W=0 m/s
    ! a1 parameter for cropland fire in Li et. al. 2013 (was different in paper)
    real(r8), parameter :: cropfire_a1 = 0.153_r8
    ! c parameter for peatland fire in Li et. al. 2013
    ! boreal peat fires (was different in paper)
    real(r8), parameter :: boreal_peatfire_c = 2.1d-5
    ! non-boreal peat fires (was different in paper)
    real(r8), parameter :: non_boreal_peatfire_c = 0.0005d00

    integer :: g,l,c,p,pi,j,fc,fp,kyr, kmo, kda, mcsec   ! index variables
    real(r8):: dt        ! time step variable (s)
    real(r8):: m         ! top-layer soil moisture (proportion)
    real(r8):: dayspyr   ! days per year
    real(r8) ::cli       !
    real(r8), parameter ::cli_scale = 1370.0_r8
    real(r8) ::cri       !
    real(r8):: fb        ! availability of fuel 
    real(r8):: fhd       ! impact of hd on agricultural fire
    real(r8):: fgdp      ! impact of gdp on agricultural fire
    real(r8):: fire_m    ! combustability of fuel on fire occurrence
    real(r8):: spread_m  ! combustability of fuel on fire spread
    real(r8):: Lb_lf     ! length-to-breadth ratio added by Lifang 
    integer :: i_cwd     ! cwd pool
    real(r8) :: lh       !
    real(r8) :: fs       !
    real(r8) :: ig       !
    real(r8) :: hdmlf    ! human density
    !-----------------------------------------------------------------------

   associate(& 
   prec60                              =>    pps%prec60                   , & ! Input:  [real(r8) (:)]  60-day running mean of tot. precipitation         
   prec10                              =>    pps%prec10                   , & ! Input:  [real(r8) (:)]  10-day running mean of tot. precipitation         
   deadcrootc                          =>    pcs%deadcrootc               , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   deadcrootc_storage                  =>    pcs%deadcrootc_storage       , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   deadcrootc_xfer                     =>    pcs%deadcrootc_xfer          , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   frootc                              =>    pcs%frootc                   , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
   frootc_storage                      =>    pcs%frootc_storage           , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   frootc_xfer                         =>    pcs%frootc_xfer              , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   livecrootc                          =>    pcs%livecrootc               , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   livecrootc_storage                  =>    pcs%livecrootc_storage       , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   livecrootc_xfer                     =>    pcs%livecrootc_xfer          , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   totvegc                             =>    pcs%totvegc                  , & ! Input:  [real(r8) (:)]  (gC/m2) total vegetation carbon, excluding cpool  
   btran2                              =>    pps%btran2                   , & ! Input:  [real(r8) (:)]  root zone soil wetness                            
   leafc                               =>    pcs%leafc                    , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
   leafc_storage                       =>    pcs%leafc_storage            , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   leafc_xfer                          =>    pcs%leafc_xfer               , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   lfpftd                              =>    pps%lfpftd                   , & ! Input:  [real(r8) (:)]  decrease of pft weight (0-1) on the col. for timestep
   burndate                            =>    pps%burndate                 , & ! Input:  [integer (:)]  burn date for crop                                 
   wf                                  =>    cps%wf                       , & ! Input:  [real(r8) (:)]  soil water as frac. of whc for top 0.05 m         
   wf2                                 =>    cps%wf2                      , & ! Input:  [real(r8) (:)]  soil water as frac. of whc for top 0.17 m         
   tsoi17                              =>    ces%tsoi17                   , & ! Input:  [real(r8) (:)]  soil T for top 0.17 m                             
   farea_burned                        =>    cps%farea_burned             , & ! InOut:  [real(r8) (:)]  fractional area burned in this timestep           
   baf_crop                            =>    cps%baf_crop                 , & ! Input:  [real(r8) (:)]  burned area fraction for cropland                 
   baf_peatf                           =>    cps%baf_peatf                , & ! Input:  [real(r8) (:)]  burned area fraction for peatland                 
   fbac                                =>    cps%fbac                     , & ! Input:  [real(r8) (:)]  total burned area out of conversion               
   fbac1                               =>    cps%fbac1                    , & ! Input:  [real(r8) (:)]  burned area out of conversion region due to land use fire
   cropf_col                           =>    cps%cropf_col                , & ! Input:  [real(r8) (:)]  cropland fraction in veg column                   
   gdp_lf                              =>    cps%gdp_lf                   , & ! Input:  [real(r8) (:)]  gdp data                                          
   peatf_lf                            =>    cps%peatf_lf                 , & ! Input:  [real(r8) (:)]  peatland fraction data                            
   abm_lf                              =>    cps%abm_lf                   , & ! Input:  [integer (:)]  proscribed crop fire time                          
   nfire                               =>    cps%nfire                    , & ! InOut:  [real(r8) (:)]  fire counts (count/km2/timestep), valid only in Reg. C
   totlitc                             =>    ccs%totlitc                  , & ! Input:  [real(r8) (:)]  (gC/m2) total lit C (column-level mean)           
   fsr_col                             =>    cps%fsr_col                  , & ! Input:  [real(r8) (:)]  fire spread rate at column level                  
   fd_col                              =>    cps%fd_col                   , & ! Input:  [real(r8) (:)]  fire duration at column level                     
   rootc_col                           =>    ccs%rootc_col                , & ! Input:  [real(r8) (:)]  root carbon                                       
   totvegc_col                         =>    ccs%totvegc_col              , & ! Input:  [real(r8) (:)]  totvegc at column level                           
   leafc_col                           =>    ccs%leafc_col                , & ! Input:  [real(r8) (:)]  leaf carbon at column level                       
   lgdp_col                            =>    cps%lgdp_col                 , & ! Input:  [real(r8) (:)]  gdp limitation factor for nfire                   
   lgdp1_col                           =>    cps%lgdp1_col                , & ! Input:  [real(r8) (:)]  gdp limitation factor for baf per fire            
   lpop_col                            =>    cps%lpop_col                 , & ! Input:  [real(r8) (:)]  pop limitation factor for baf per fire            
   fuelc                               =>    ccs%fuelc                    , & ! Input:  [real(r8) (:)]  fuel avalability factor for Reg.C                 
   fuelc_crop                          =>    ccs%fuelc_crop               , & ! Input:  [real(r8) (:)]  fuel avalability factor for Reg.A                 
   btran_col                           =>    cps%btran_col                , & ! Input:  [real(r8) (:)]  transpiration wetness factor (0 to 1)             
   wtlf                                =>    cps%wtlf                     , & ! Input:  [real(r8) (:)]  fractional coverage of non-crop PFTs              
   lfwt                                =>    cps%lfwt                     , & ! Input:  [real(r8) (:)]  fractional coverage of non-crop and non-bare-soil PFTs
   trotr1_col                          =>    cps%trotr1_col               , & ! Input:  [real(r8) (:)]  pft weight of BET on the gridcell (0-1)           
   trotr2_col                          =>    cps%trotr2_col               , & ! Input:  [real(r8) (:)]  pft weight of BDT on the gridcell (0-1)           
   dtrotr_col                          =>    cps%dtrotr_col               , & ! Input:  [real(r8) (:)]  annual decreased fraction coverage of BET+BDT on gridcell
   prec60_col                          =>    cps%prec60_col               , & ! Input:  [real(r8) (:)]  60-day running mean of tot. precipitation         
   prec10_col                          =>    cps%prec10_col               , & ! Input:  [real(r8) (:)]  10-day running mean of tot. precipitation         
   lfc                                 =>    cps%lfc                      , & ! Input:  [real(r8) (:)]  conversion area frac. of BET+BDT that haven't burned before
   fsat                                =>    cws%fsat                     , & ! Input:  [real(r8) (:)]  fractional area with water table at surface       
   is_cwd                              =>    decomp_cascade_con%is_cwd    , & ! InOut:  [logical (:)]  TRUE => pool is a cwd pool                         
   decomp_cpools_vr                    =>    ccs%decomp_cpools_vr         , & ! Input:  [real(r8) (:,:,:)]  (gC/m3)  vert.-resolved decomposing c pools   
   forc_rh                             =>    clm_a2l%forc_rh              , & ! Input:  [real(r8) (:)]  relative humidity                                 
   forc_wind                           =>    clm_a2l%forc_wind            , & ! Input:  [real(r8) (:)] atmospheric wind speed (m/s)                       
   forc_t                              =>    a2l_downscaled_col%forc_t    , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)                  
   forc_rain                           =>    a2l_downscaled_col%forc_rain , & ! Input:  [real(r8) (:)]  rain                                              
   forc_snow                           =>    a2l_downscaled_col%forc_snow   & ! Input:  [real(r8) (:)]  snow                                              
   )
 
     !pft to column average 
     call p2c(bounds, num_soilc, filter_soilc, &
          prec10(bounds%begp:bounds%endp), &
          prec10_col(bounds%begc:bounds%endc))
     call p2c(bounds, num_soilc, filter_soilc, &
          prec60(bounds%begp:bounds%endp), &
          prec60_col(bounds%begc:bounds%endc))
     call p2c(bounds, num_soilc, filter_soilc, &
          totvegc(bounds%begp:bounds%endp), &
          totvegc_col(bounds%begc:bounds%endc))
     call p2c(bounds, num_soilc, filter_soilc, &
          leafc(bounds%begp:bounds%endp), &
          leafc_col(bounds%begc:bounds%endc))

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
     do pi = 1,max_pft_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           if (pi <=  col%npfts(c)) then
              p = col%pfti(c) + pi - 1
              ! For crop veg types
              if( pft%itype(p) > nc4_grass )then
                 cropf_col(c) = cropf_col(c) + pft%wtcol(p)
              end if
              ! For natural vegetation (non-crop)
              if( pft%itype(p) >= ndllf_evr_tmp_tree .and. pft%itype(p) <= nc4_grass )then
                 lfwt(c) = lfwt(c) + pft%wtcol(p)
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
     do pi = 1,max_pft_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           if (pi <=  col%npfts(c)) then
              p = col%pfti(c) + pi - 1
              ! For crop PFT's
              if( pft%itype(p) > nc4_grass .and. pft%wtcol(p) > 0._r8 .and. leafc_col(c) > 0._r8 )then
                 fuelc_crop(c)=fuelc_crop(c) + (leafc(p) + leafc_storage(p) + &
                      leafc_xfer(p))*pft%wtcol(p)/cropf_col(c)     + &
                      totlitc(c)*leafc(p)/leafc_col(c)*pft%wtcol(p)/cropf_col(c)
              end if
           end if
        end do
     end do
     !   
     ! Calculate noncrop column variables
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
        if (fpftdyn /= ' ') then    !true when landuse data is used
           trotr1_col(c)=0._r8
           trotr2_col(c)=0._r8
           dtrotr_col(c)=0._r8
        end if
     end do
     do pi = 1,max_pft_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           g = col%gridcell(c)
           if (pi <=  col%npfts(c)) then
              p = col%pfti(c) + pi - 1
              ! For non-crop -- natural vegetation and bare-soil
              if( pft%itype(p) .lt. nc3crop .and. cropf_col(c) .lt. 1.0_r8 )then
                 if( .not. shr_infnan_isnan(btran2(p)) .and. btran2(p) .le. 1._r8 )then
                    btran_col(c) = btran_col(c)+btran2(p)*pft%wtcol(p)
                    wtlf(c)      = wtlf(c)+pft%wtcol(p)
                 end if
                 if ( fpftdyn /= ' ' ) then    !true when landuse data is used
                    if( pft%itype(p) == nbrdlf_evr_trp_tree .and. pft%wtcol(p) .gt. 0._r8 )then
                       trotr1_col(c)=trotr1_col(c)+pft%wtcol(p)*col%wtgcell(c)
                    end if
                    if( pft%itype(p) == nbrdlf_dcd_trp_tree .and. pft%wtcol(p) .gt. 0._r8 )then
                       trotr2_col(c)=trotr2_col(c)+pft%wtcol(p)*col%wtgcell(c)
                    end if
                    if( pft%itype(p) == nbrdlf_evr_trp_tree .or. pft%itype(p) == nbrdlf_dcd_trp_tree )then
                       if(lfpftd(p).gt.0._r8)then
                          dtrotr_col(c)=dtrotr_col(c)+lfpftd(p)*col%wtgcell(c)
                       end if
                    end if
                 end if
                 rootc_col(c) = rootc_col(c) + (frootc(p) + frootc_storage(p) + &
                      frootc_xfer(p) + deadcrootc(p) +                &
                      deadcrootc_storage(p) + deadcrootc_xfer(p) +    &
                      livecrootc(p)+livecrootc_storage(p) +           &
                      livecrootc_xfer(p))*pft%wtcol(p)

                 fsr_col(c) = fsr_col(c) + fsr_pft(pft%itype(p))*pft%wtcol(p)/(1.0_r8-cropf_col(c))

                 if( lfwt(c) .ne. 0.0_r8 )then    
                    hdmlf=forc_hdm(g)

                    ! all these constants are in Li et al. BG (2012a,b;2013)

                    if( hdmlf .gt. 0.1_r8 )then            
                       ! For NOT bare-soil
                       if( pft%itype(p) .ne. noveg )then
                          ! For shrub and grass (crop already excluded above)
                          if( pft%itype(p) .ge. nbrdlf_evr_shrub )then      !for shurb and grass
                             lgdp_col(c)  = lgdp_col(c) + (0.1_r8 + 0.9_r8*    &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (gdp_lf(c)/8._r8)**0.5_r8))*pft%wtcol(p) &
                                  /(1.0_r8 - cropf_col(c))
                             lgdp1_col(c) = lgdp1_col(c) + (0.2_r8 + 0.8_r8*   &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (gdp_lf(c)/7._r8)))*pft%wtcol(p)/lfwt(c)
                             lpop_col(c)  = lpop_col(c) + (0.2_r8 + 0.8_r8*    &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (hdmlf/450._r8)**0.5_r8))*pft%wtcol(p)/lfwt(c)
                          else   ! for trees
                             if( gdp_lf(c) .gt. 20._r8 )then
                                lgdp_col(c)  =lgdp_col(c)+0.39_r8*pft%wtcol(p)/(1.0_r8 - cropf_col(c))
                             else    
                                lgdp_col(c) = lgdp_col(c)+pft%wtcol(p)/(1.0_r8 - cropf_col(c))
                             end if
                             if( gdp_lf(c) .gt. 20._r8 )then   
                                lgdp1_col(c) = lgdp1_col(c)+0.62_r8*pft%wtcol(p)/lfwt(c)
                             else
                                if( gdp_lf(c) .gt. 8._r8 ) then
                                   lgdp1_col(c)=lgdp1_col(c)+0.83_r8*pft%wtcol(p)/lfwt(c)
                                else
                                   lgdp1_col(c)=lgdp1_col(c)+pft%wtcol(p)/lfwt(c)
                                end if
                             end if
                             lpop_col(c) = lpop_col(c) + (0.4_r8 + 0.6_r8*    &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (hdmlf/125._r8)))*pft%wtcol(p)/lfwt(c) 
                          end if
                       end if
                    else
                       lgdp_col(c)  = lgdp_col(c)+pft%wtcol(p)/(1.0_r8 - cropf_col(c))
                       lgdp1_col(c) = lgdp1_col(c)+pft%wtcol(p)/lfwt(c)
                       lpop_col(c)  = lpop_col(c)+pft%wtcol(p)/lfwt(c)
                    end if
                 end if

                 fd_col(c) = fd_col(c) + fd_pft(pft%itype(p))*pft%wtcol(p)*secsphr/(1.0_r8-cropf_col(c))         
              end if
           end if
        end do
     end do

     if (fpftdyn /= ' ') then    !true when landuse data is used
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           if( dtrotr_col(c) .gt. 0._r8 )then
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

     do pi = 1,max_pft_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           g= col%gridcell(c)
           hdmlf=forc_hdm(g)
           if (pi <=  col%npfts(c)) then
              p = col%pfti(c) + pi - 1
              ! For crop
              if( forc_t(c) .ge. SHR_CONST_TKFRZ .and. pft%itype(p) .gt. nc4_grass .and.  &
                   kmo == abm_lf(c) .and. forc_rain(c)+forc_snow(c) .eq. 0._r8  .and. &
                   burndate(p) >= 999 .and. pft%wtcol(p) .gt. 0._r8 )then ! catch  crop burn time
                 ! calculate human density impact on ag. fire
                 fhd = 0.04_r8+0.96_r8*exp(-1._r8*SHR_CONST_PI*(hdmlf/350._r8)**0.5_r8)
                 ! calculate impact of GDP on ag. fire
                 fgdp = 0.01_r8+0.99_r8*exp(-1._r8*SHR_CONST_PI*(gdp_lf(c)/10._r8))
                 ! calculate burned area
                 fb   = max(0.0_r8,min(1.0_r8,(fuelc_crop(c)-lfuel)/(ufuel-lfuel)))
                 ! crop fire only for generic crop types at this time
                 ! managed crops are treated as grasses if crop model is turned on
                 ! NOTE: THIS SHOULD TAKE INTO ACCOUNT THE TIME-STEP AND CURRENTLY DOES NOT!
                 !       As such results are only valid for a time-step of a half-hour.
                 baf_crop(c) = baf_crop(c) + cropfire_a1*fb*fhd*fgdp*pft%wtcol(p)
                 if( fb*fhd*fgdp*pft%wtcol(p) .gt. 0._r8)then
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
        g= col%gridcell(c)
        ! NOTE: THIS SHOULD TAKE INTO ACCOUNT THE TIME-STEP AND CURRENTLY DOES NOT!
        !       As such results are only valid for a time-step of a half-hour.
        if(grc%latdeg(g).lt.borealat )then
           baf_peatf(c) = non_boreal_peatfire_c*max(0._r8, &
                min(1._r8,(4.0_r8-prec60_col(c)*secspday)/ &
                4.0_r8))**2*peatf_lf(c)*(1._r8-fsat(c))
        else
           baf_peatf(c) = boreal_peatfire_c*exp(-SHR_CONST_PI*(max(wf2(c),0._r8)/0.3_r8))* &
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
        g = col%gridcell(c)
        hdmlf=forc_hdm(g)
        if( cropf_col(c) .lt. 1.0 )then
           fuelc(c) = totlitc(c)+totvegc_col(c)-rootc_col(c)-fuelc_crop(c)*cropf_col(c)
           do j = 1, nlevdecomp  
              fuelc(c) = fuelc(c)+decomp_cpools_vr(c,j,i_cwd) * dzsoi_decomp(j)
           end do
           fuelc(c) = fuelc(c)/(1._r8-cropf_col(c))
           fb       = max(0.0_r8,min(1.0_r8,(fuelc(c)-lfuel)/(ufuel-lfuel)))
           m        = max(0._r8,wf(c))
           fire_m   = exp(-SHR_CONST_PI *(m/0.69_r8)**2)*(1.0_r8 - max(0._r8, &
                min(1._r8,(forc_rh(g)-30._r8)/(70._r8-30._r8))))*  &
                min(1._r8,exp(SHR_CONST_PI*(forc_t(c)-SHR_CONST_TKFRZ)/10._r8))
           lh       = 0.0035_r8*6.8_r8*hdmlf**(0.43_r8)/30._r8/24._r8
           fs       = 1._r8-(0.01_r8+0.98_r8*exp(-0.025_r8*hdmlf))
           ig       = (lh+forc_lnfm(g)*0.25_r8)*(1._r8-fs)*(1._r8-cropf_col(c)) 
           nfire(c) = ig/secsphr*dt*fb*fire_m*lgdp_col(c) !fire counts/km2/timestep
           Lb_lf    = 1._r8+10.0_r8*(1._r8-EXP(-0.06_r8*forc_wind(g)))
           if ( wtlf(c) > 0.0_r8 )then
              spread_m = (1.0_r8 - max(0._r8,min(1._r8,(btran_col(c)/wtlf(c)-0.3_r8)/ &
                   (0.7_r8-0.3_r8))))*(1.0-max(0._r8, &
                   min(1._r8,(forc_rh(g)-30._r8)/(70._r8-30._r8))))
           else
              spread_m = 0.0_r8
           end if
           farea_burned(c) = min(1._r8,(g0*spread_m*fsr_col(c)* &
                fd_col(c)/1000._r8)**2*lgdp1_col(c)* &
                lpop_col(c)*nfire(c)*SHR_CONST_PI*Lb_lf+ &
                baf_crop(c)+baf_peatf(c))  ! fraction (0-1) per timestep

           !
           ! if landuse change data is used, calculate deforestation fires and 
           ! add it in the total of burned area fraction
           !
           if (fpftdyn /= ' ') then    !true when landuse change data is used
              if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 )then
                 if(( kmo == 1 .and. kda == 1 .and. mcsec == 0) .or. &
                      dtrotr_col(c) <=0._r8 )then
                    fbac1(c)        = 0._r8
                    farea_burned(c) = baf_crop(c)+baf_peatf(c)
                 else
                    cri = (4.0_r8*trotr1_col(c)+1.8_r8*trotr2_col(c))/(trotr1_col(c)+trotr2_col(c))
                    cli = (max(0._r8,min(1._r8,(cri-prec60_col(c)*secspday)/cri))**0.5)* &
                         (max(0._r8,min(1._r8,(cri-prec10_col(c)*secspday)/cri))**0.5)* &
                         max(0.0005_r8,min(1._r8,19._r8*dtrotr_col(c)*dayspyr*secspday/dt-0.001_r8))* &
                         max(0._r8,min(1._r8,(0.25_r8-(forc_rain(c)+forc_snow(c))*secsphr)/0.25_r8))
                    ! NOTE: THIS SHOULD TAKE INTO ACCOUNT THE TIME-STEP AND CURRENTLY DOES NOT!
                    !       As such results are only valid for a time-step of a half-hour.
                    farea_burned(c) = cli/cli_scale +baf_crop(c)+baf_peatf(c)
                    ! burned area out of conversion region due to land use fire
                    fbac1(c) = max(0._r8,cli/cli_scale - 2.0_r8*lfc(c))   
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
 end subroutine CNFireArea

 !-----------------------------------------------------------------------
 subroutine CNFireFluxes (num_soilc, filter_soilc, num_soilp, filter_soilp)
   !
   ! !DESCRIPTION:
   ! Fire effects routine for coupled carbon-nitrogen code (CN).
   ! Relies primarily on estimate of fractional area burned in this
   ! timestep, from CNFireArea().
   !
   ! Total fire carbon emissions (g C/m2 land area/yr) 
   !  =avg(COL_FIRE_CLOSS)*seconds_per_year + avg(SOMC_FIRE)*seconds_per_year + 
   !   avg(LF_CONV_CFLUX)*seconds_per_year*min(1.0,avg(LFC2)/dt*seconds_per_year)*0.8
   ! where dt is the time step size (sec),avg means the temporal average in a year
   ! seconds_per_year is the number of seconds in a year.
   !
   ! !USES:
   use clmtype
   use pftvarcon, only: cc_leaf,cc_lstem,cc_dstem,cc_other,fm_leaf,fm_lstem,fm_other,fm_root,fm_lroot,fm_droot
   use pftvarcon, only: nc3crop,lf_flab,lf_fcel,lf_flig,fr_flab,fr_fcel,fr_flig
   use clm_time_manager, only: get_step_size,get_days_per_year,get_curr_date
   use clm_varpar  , only: max_pft_per_col
   use clm_varctl  , only: fpftdyn, use_cndv
   use clm_varcon  , only: secspday
   !
   ! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
   !
   ! !LOCAL VARIABLES:
   integer :: g,c,p,j,l,pi,kyr, kmo, kda, mcsec   ! indices
   integer :: fp,fc                ! filter indices
   real(r8):: f                    ! rate for fire effects (1/s)
   real(r8):: dt                   ! time step variable (s)
   real(r8):: dayspyr              ! days per year
   !-----------------------------------------------------------------------

   associate(& 
   nind                                =>    pdgvs%nind                                  , & ! Input:  [real(r8) (:)]  number of individuals (#/m2)                      
   woody                               =>    pftcon%woody                                , & ! Input:  [real(r8) (:)]  woody lifeform (1=woody, 0=not woody)             
   farea_burned                        =>    cps%farea_burned                            , & ! Input:  [real(r8) (:)]  timestep fractional area burned (proportion)      
   fire_mortality_c_to_cwdc            =>    ccf%fire_mortality_c_to_cwdc                , & ! Input:  [real(r8) (:,:)]  C fluxes associated with fire mortality to CWD pool (gC/m3/s)
   fire_mortality_n_to_cwdn            =>    cnf%fire_mortality_n_to_cwdn                , & ! Input:  [real(r8) (:,:)]  N fluxes associated with fire mortality to CWD pool (gN/m3/s)
   lfc                                 =>    cps%lfc                                     , & ! Input:  [real(r8) (:)]  conversion area frac. of BET+BDT that haven't burned before
   lfc2                                =>    cps%lfc2                                    , & ! Input:  [real(r8) (:)]  conversion area frac. of BET+BDT that burned this timestep
   fbac1                               =>    cps%fbac1                                   , & ! Input:  [real(r8) (:)]  burned area out of conversion region due to land use fire
   fbac                                =>    cps%fbac                                    , & ! Input:  [real(r8) (:)]  total burned area out of conversion               
   baf_crop                            =>    cps%baf_crop                                , & ! Input:  [real(r8) (:)]  baf for cropland                                  
   baf_peatf                           =>    cps%baf_peatf                               , & ! Input:  [real(r8) (:)]  baf for peatlabd                                  
   leafcmax                            =>    pcs%leafcmax                                , & ! Input:  [real(r8) (:)]  (gC/m2) ann max leaf C                            
   trotr1_col                          =>    cps%trotr1_col                              , & ! Input:  [real(r8) (:)]  pft weight of BET on the gridcell (0-1)           
   trotr2_col                          =>    cps%trotr2_col                              , & ! Input:  [real(r8) (:)]  pft weight of BDT on the gridcell (0-1)           
   dtrotr_col                          =>    cps%dtrotr_col                              , & ! Input:  [real(r8) (:)]  annual decreased fraction coverage of BET+BDT (0-1) on the gridcell
   somc_fire                           =>    ccf%somc_fire                               , & ! Input:  [real(r8) (:)]  (gC/m2/s)fire carbon emissions due to peat burning
   totsomc                             =>    ccs%totsomc                                 , & ! Input:  [real(r8) (:)]  (gC/m2) total soil organic matter carbon          
   decomp_cpools_vr                    =>    ccs%decomp_cpools_vr                        , & ! Input:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   decomp_npools_vr                    =>    cns%decomp_npools_vr                        , & ! Input:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   leafc                               =>    pcs%leafc                                   , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
   leafc_storage                       =>    pcs%leafc_storage                           , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   leafc_xfer                          =>    pcs%leafc_xfer                              , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   livestemc                           =>    pcs%livestemc                               , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
   livestemc_storage                   =>    pcs%livestemc_storage                       , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   livestemc_xfer                      =>    pcs%livestemc_xfer                          , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   deadstemc                           =>    pcs%deadstemc                               , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
   deadstemc_storage                   =>    pcs%deadstemc_storage                       , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   deadstemc_xfer                      =>    pcs%deadstemc_xfer                          , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   frootc                              =>    pcs%frootc                                  , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
   frootc_storage                      =>    pcs%frootc_storage                          , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   frootc_xfer                         =>    pcs%frootc_xfer                             , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   livecrootc                          =>    pcs%livecrootc                              , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   livecrootc_storage                  =>    pcs%livecrootc_storage                      , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   livecrootc_xfer                     =>    pcs%livecrootc_xfer                         , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   deadcrootc                          =>    pcs%deadcrootc                              , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   deadcrootc_storage                  =>    pcs%deadcrootc_storage                      , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   deadcrootc_xfer                     =>    pcs%deadcrootc_xfer                         , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   gresp_storage                       =>    pcs%gresp_storage                           , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   gresp_xfer                          =>    pcs%gresp_xfer                              , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   leafn                               =>    pns%leafn                                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
   leafn_storage                       =>    pns%leafn_storage                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
   leafn_xfer                          =>    pns%leafn_xfer                              , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
   livestemn                           =>    pns%livestemn                               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
   livestemn_storage                   =>    pns%livestemn_storage                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
   livestemn_xfer                      =>    pns%livestemn_xfer                          , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
   deadstemn                           =>    pns%deadstemn                               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
   deadstemn_storage                   =>    pns%deadstemn_storage                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
   deadstemn_xfer                      =>    pns%deadstemn_xfer                          , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
   frootn                              =>    pns%frootn                                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
   frootn_storage                      =>    pns%frootn_storage                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
   frootn_xfer                         =>    pns%frootn_xfer                             , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
   livecrootn                          =>    pns%livecrootn                              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
   livecrootn_storage                  =>    pns%livecrootn_storage                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
   livecrootn_xfer                     =>    pns%livecrootn_xfer                         , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
   deadcrootn                          =>    pns%deadcrootn                              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
   deadcrootn_storage                  =>    pns%deadcrootn_storage                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
   deadcrootn_xfer                     =>    pns%deadcrootn_xfer                         , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
   retransn                            =>    pns%retransn                                , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   pactive                             =>    pft%active                                  , & ! Input:  [logical (:)]  true=>do computations on this pft 

   m_leafc_to_fire                     =>    pcf%m_leafc_to_fire                         , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from leafc             
   m_leafc_storage_to_fire             =>    pcf%m_leafc_storage_to_fire                 , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from leafc_storage     
   m_leafc_xfer_to_fire                =>    pcf%m_leafc_xfer_to_fire                    , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from leafc_xfer        
   m_livestemc_to_fire                 =>    pcf%m_livestemc_to_fire                     , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from livestemc         
   m_livestemc_storage_to_fire         =>    pcf%m_livestemc_storage_to_fire             , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from livestemc_storage 
   m_livestemc_xfer_to_fire            =>    pcf%m_livestemc_xfer_to_fire                , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from livestemc_xfer    
   m_deadstemc_to_fire                 =>    pcf%m_deadstemc_to_fire                     , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from deadstemc_xfer    
   m_deadstemc_storage_to_fire         =>    pcf%m_deadstemc_storage_to_fire             , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from deadstemc_storage 
   m_deadstemc_xfer_to_fire            =>    pcf%m_deadstemc_xfer_to_fire                , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from deadstemc_xfer    
   m_frootc_to_fire                    =>    pcf%m_frootc_to_fire                        , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from frootc            
   m_frootc_storage_to_fire            =>    pcf%m_frootc_storage_to_fire                , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from frootc_storage    
   m_frootc_xfer_to_fire               =>    pcf%m_frootc_xfer_to_fire                   , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from frootc_xfer       
   m_livecrootc_to_fire                =>    pcf%m_livecrootc_to_fire                    , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from livecrootc        
   m_livecrootc_storage_to_fire        =>    pcf%m_livecrootc_storage_to_fire            , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from livecrootc_storage
   m_livecrootc_xfer_to_fire           =>    pcf%m_livecrootc_xfer_to_fire               , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from livecrootc_xfer   
   m_deadcrootc_to_fire                =>    pcf%m_deadcrootc_to_fire                    , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from deadcrootc        
   m_deadcrootc_storage_to_fire        =>    pcf%m_deadcrootc_storage_to_fire            , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from deadcrootc_storage
   m_deadcrootc_xfer_to_fire           =>    pcf%m_deadcrootc_xfer_to_fire               , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from deadcrootc_xfer   
   m_gresp_storage_to_fire             =>    pcf%m_gresp_storage_to_fire                 , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from gresp_storage     
   m_gresp_xfer_to_fire                =>    pcf%m_gresp_xfer_to_fire                    , & ! Input:  [real(r8) (:)]  (gC/m2/s) fire C emissions from gresp_xfer        

   m_leafn_to_fire                     =>    pnf%m_leafn_to_fire                         , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from leafn             
   m_leafn_storage_to_fire             =>    pnf%m_leafn_storage_to_fire                 , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from leafn_storage     
   m_leafn_xfer_to_fire                =>    pnf%m_leafn_xfer_to_fire                    , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from leafn_xfer        
   m_livestemn_to_fire                 =>    pnf%m_livestemn_to_fire                     , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from livestemn         
   m_livestemn_storage_to_fire         =>    pnf%m_livestemn_storage_to_fire             , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from livestemn_storage 
   m_livestemn_xfer_to_fire            =>    pnf%m_livestemn_xfer_to_fire                , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from livestemn_xfer    
   m_deadstemn_to_fire                 =>    pnf%m_deadstemn_to_fire                     , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from deadstemn         
   m_deadstemn_storage_to_fire         =>    pnf%m_deadstemn_storage_to_fire             , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from deadstemn_storage 
   m_deadstemn_xfer_to_fire            =>    pnf%m_deadstemn_xfer_to_fire                , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from deadstemn_xfer    
   m_frootn_to_fire                    =>    pnf%m_frootn_to_fire                        , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from frootn            
   m_frootn_storage_to_fire            =>    pnf%m_frootn_storage_to_fire                , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from frootn_storage    
   m_frootn_xfer_to_fire               =>    pnf%m_frootn_xfer_to_fire                   , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from frootn_xfer       
   m_livecrootn_to_fire                =>    pnf%m_livecrootn_to_fire                    , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from m_livecrootn_to_fire
   m_livecrootn_storage_to_fire        =>    pnf%m_livecrootn_storage_to_fire            , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from livecrootn_storage
   m_livecrootn_xfer_to_fire           =>    pnf%m_livecrootn_xfer_to_fire               , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from livecrootn_xfer   
   m_deadcrootn_to_fire                =>    pnf%m_deadcrootn_to_fire                    , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from deadcrootn        
   m_deadcrootn_storage_to_fire        =>    pnf%m_deadcrootn_storage_to_fire            , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from deadcrootn_storage
   m_deadcrootn_xfer_to_fire           =>    pnf%m_deadcrootn_xfer_to_fire               , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from deadcrootn_xfer   
   m_retransn_to_fire                  =>    pnf%m_retransn_to_fire                      , & ! Input:  [real(r8) (:)]  (gN/m2/s) fire N emissions from retransn          

   m_leafc_to_litter_fire              =>    pcf%m_leafc_to_litter_fire                  , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_storage_to_litter_fire      =>    pcf%m_leafc_storage_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_xfer_to_litter_fire         =>    pcf%m_leafc_xfer_to_litter_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_litter_fire          =>    pcf%m_livestemc_to_litter_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_storage_to_litter_fire  =>    pcf%m_livestemc_storage_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_xfer_to_litter_fire     =>    pcf%m_livestemc_xfer_to_litter_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_deadstemc_fire       =>    pcf%m_livestemc_to_deadstemc_fire           , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_to_litter_fire          =>    pcf%m_deadstemc_to_litter_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_storage_to_litter_fire  =>    pcf%m_deadstemc_storage_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_xfer_to_litter_fire     =>    pcf%m_deadstemc_xfer_to_litter_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_to_litter_fire             =>    pcf%m_frootc_to_litter_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_storage_to_litter_fire     =>    pcf%m_frootc_storage_to_litter_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_xfer_to_litter_fire        =>    pcf%m_frootc_xfer_to_litter_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_litter_fire         =>    pcf%m_livecrootc_to_litter_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_storage_to_litter_fire =>    pcf%m_livecrootc_storage_to_litter_fire     , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_xfer_to_litter_fire    =>    pcf%m_livecrootc_xfer_to_litter_fire        , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_deadcrootc_fire     =>    pcf%m_livecrootc_to_deadcrootc_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_to_litter_fire         =>    pcf%m_deadcrootc_to_litter_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_storage_to_litter_fire =>    pcf%m_deadcrootc_storage_to_litter_fire     , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_xfer_to_litter_fire    =>    pcf%m_deadcrootc_xfer_to_litter_fire        , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_storage_to_litter_fire      =>    pcf%m_gresp_storage_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_xfer_to_litter_fire         =>    pcf%m_gresp_xfer_to_litter_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_decomp_cpools_to_fire_vr          =>    ccf%m_decomp_cpools_to_fire_vr              , & ! Input:  [real(r8) (:,:,:)]  (gC/m3/s) vertically-resolved decomposing C fire loss
   m_c_to_litr_met_fire                =>    ccf%m_c_to_litr_met_fire                    , & ! Input:  [real(r8) (:,:)]                                                  
   m_c_to_litr_cel_fire                =>    ccf%m_c_to_litr_cel_fire                    , & ! Input:  [real(r8) (:,:)]                                                  
   m_c_to_litr_lig_fire                =>    ccf%m_c_to_litr_lig_fire                    , & ! Input:  [real(r8) (:,:)]                                                  

   m_leafn_to_litter_fire              =>    pnf%m_leafn_to_litter_fire                  , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_storage_to_litter_fire      =>    pnf%m_leafn_storage_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_xfer_to_litter_fire         =>    pnf%m_leafn_xfer_to_litter_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_to_litter_fire          =>    pnf%m_livestemn_to_litter_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_storage_to_litter_fire  =>    pnf%m_livestemn_storage_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_xfer_to_litter_fire     =>    pnf%m_livestemn_xfer_to_litter_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_to_deadstemn_fire       =>    pnf%m_livestemn_to_deadstemn_fire           , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_to_litter_fire          =>    pnf%m_deadstemn_to_litter_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_storage_to_litter_fire  =>    pnf%m_deadstemn_storage_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_xfer_to_litter_fire     =>    pnf%m_deadstemn_xfer_to_litter_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_to_litter_fire             =>    pnf%m_frootn_to_litter_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_storage_to_litter_fire     =>    pnf%m_frootn_storage_to_litter_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_xfer_to_litter_fire        =>    pnf%m_frootn_xfer_to_litter_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_to_litter_fire         =>    pnf%m_livecrootn_to_litter_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_storage_to_litter_fire =>    pnf%m_livecrootn_storage_to_litter_fire     , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_xfer_to_litter_fire    =>    pnf%m_livecrootn_xfer_to_litter_fire        , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_to_deadcrootn_fire     =>    pnf%m_livecrootn_to_deadcrootn_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_to_litter_fire         =>    pnf%m_deadcrootn_to_litter_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_storage_to_litter_fire =>    pnf%m_deadcrootn_storage_to_litter_fire     , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_xfer_to_litter_fire    =>    pnf%m_deadcrootn_xfer_to_litter_fire        , & ! Input:  [real(r8) (:)]                                                    
   m_retransn_to_litter_fire           =>    pnf%m_retransn_to_litter_fire               , & ! Input:  [real(r8) (:)]                                                    
   m_decomp_npools_to_fire_vr          =>    cnf%m_decomp_npools_to_fire_vr              , & ! Input:  [real(r8) (:,:,:)]  vertically-resolved decomposing N fire loss (gN/m3/s)
   m_n_to_litr_met_fire                =>    cnf%m_n_to_litr_met_fire                    , & ! Input:  [real(r8) (:,:)]                                                  
   m_n_to_litr_cel_fire                =>    cnf%m_n_to_litr_cel_fire                    , & ! Input:  [real(r8) (:,:)]                                                  
   m_n_to_litr_lig_fire                =>    cnf%m_n_to_litr_lig_fire                    , & ! Input:  [real(r8) (:,:)]                                                  
   
   is_cwd                              =>    decomp_cascade_con%is_cwd                   , & ! Input:  [logical (:)]  TRUE => pool is a cwd pool                         
   is_litter                           =>    decomp_cascade_con%is_litter                , & ! Input:  [logical (:)]  TRUE => pool is a litter pool                      
   croot_prof                          =>    pps%croot_prof                              , & ! Input:  [real(r8) (:,:)]  (1/m) profile of coarse roots                   
   stem_prof                           =>    pps%stem_prof                               , & ! Input:  [real(r8) (:,:)]  (1/m) profile of stems                          
   froot_prof                          =>    pps%froot_prof                              , & ! Input:  [real(r8) (:,:)]  (1/m) profile of fine roots                     
   leaf_prof                           =>    pps%leaf_prof                                 & ! Input:  [real(r8) (:,:)]  (1/m) profile of leaves                         
   )

     ! Get model step size
     ! calculate burned area fraction per sec
     dt = real( get_step_size(), r8 )

     dayspyr = get_days_per_year()
     !
     ! pft loop
     !
     do fp = 1,num_soilp
        p = filter_soilp(fp)
        c = pft%column(p)

        ! get the column-level fractional area burned for this timestep
        ! and convert to a rate per second
        ! For non-crop (bare-soil and natural vegetation)
        if( pft%itype(p) .lt. nc3crop )then
           if (fpftdyn /= ' ') then    !true when landuse data is used
              f = (fbac(c)-baf_crop(c))/ dt
           else
              f = (farea_burned(c))/ dt
           end if
        else
           f = baf_crop(c) / dt
        end if

        ! apply this rate to the pft state variables to get flux rates
        ! biomass burning
        ! carbon fluxes
        m_leafc_to_fire(p)               =  leafc(p)              * f * cc_leaf(pft%itype(p))
        m_leafc_storage_to_fire(p)       =  leafc_storage(p)      * f * cc_other(pft%itype(p))
        m_leafc_xfer_to_fire(p)          =  leafc_xfer(p)         * f * cc_other(pft%itype(p))
        m_livestemc_to_fire(p)           =  livestemc(p)          * f * cc_lstem(pft%itype(p))
        m_livestemc_storage_to_fire(p)   =  livestemc_storage(p)  * f * cc_other(pft%itype(p))
        m_livestemc_xfer_to_fire(p)      =  livestemc_xfer(p)     * f * cc_other(pft%itype(p))
        m_deadstemc_to_fire(p)           =  deadstemc(p)          * f * cc_dstem(pft%itype(p))
        m_deadstemc_storage_to_fire(p)   =  deadstemc_storage(p)  * f * cc_other(pft%itype(p))
        m_deadstemc_xfer_to_fire(p)      =  deadstemc_xfer(p)     * f * cc_other(pft%itype(p))
        m_frootc_to_fire(p)              =  frootc(p)             * f * 0._r8
        m_frootc_storage_to_fire(p)      =  frootc_storage(p)     * f * cc_other(pft%itype(p)) 
        m_frootc_xfer_to_fire(p)         =  frootc_xfer(p)        * f * cc_other(pft%itype(p))
        m_livecrootc_to_fire(p)          =  livecrootc(p)         * f * 0._r8
        m_livecrootc_storage_to_fire(p)  =  livecrootc_storage(p) * f * cc_other(pft%itype(p)) 
        m_livecrootc_xfer_to_fire(p)     =  livecrootc_xfer(p)    * f * cc_other(pft%itype(p)) 
        m_deadcrootc_to_fire(p)          =  deadcrootc(p)         * f * 0._r8
        m_deadcrootc_storage_to_fire(p)  =  deadcrootc_storage(p) * f*  cc_other(pft%itype(p)) 
        m_deadcrootc_xfer_to_fire(p)     =  deadcrootc_xfer(p)    * f * cc_other(pft%itype(p)) 
        m_gresp_storage_to_fire(p)       =  gresp_storage(p)      * f * cc_other(pft%itype(p))
        m_gresp_xfer_to_fire(p)          =  gresp_xfer(p)         * f * cc_other(pft%itype(p))


        ! nitrogen fluxes
        m_leafn_to_fire(p)               =  leafn(p)              * f * cc_leaf(pft%itype(p))
        m_leafn_storage_to_fire(p)       =  leafn_storage(p)      * f * cc_other(pft%itype(p))
        m_leafn_xfer_to_fire(p)          =  leafn_xfer(p)         * f * cc_other(pft%itype(p))
        m_livestemn_to_fire(p)           =  livestemn(p)          * f * cc_lstem(pft%itype(p))
        m_livestemn_storage_to_fire(p)   =  livestemn_storage(p)  * f * cc_other(pft%itype(p))
        m_livestemn_xfer_to_fire(p)      =  livestemn_xfer(p)     * f * cc_other(pft%itype(p))
        m_deadstemn_to_fire(p)           =  deadstemn(p)          * f * cc_dstem(pft%itype(p))
        m_deadstemn_storage_to_fire(p)   =  deadstemn_storage(p)  * f * cc_other(pft%itype(p))
        m_deadstemn_xfer_to_fire(p)      =  deadstemn_xfer(p)     * f * cc_other(pft%itype(p))
        m_frootn_to_fire(p)              =  frootn(p)             * f * 0._r8
        m_frootn_storage_to_fire(p)      =  frootn_storage(p)     * f * cc_other(pft%itype(p))
        m_frootn_xfer_to_fire(p)         =  frootn_xfer(p)        * f * cc_other(pft%itype(p))
        m_livecrootn_to_fire(p)          =  livecrootn(p)         * f * 0._r8 
        m_livecrootn_storage_to_fire(p)  =  livecrootn_storage(p) * f * cc_other(pft%itype(p)) 
        m_livecrootn_xfer_to_fire(p)     =  livecrootn_xfer(p)    * f * cc_other(pft%itype(p))
        m_deadcrootn_to_fire(p)          =  deadcrootn(p)         * f * 0._r8
        m_deadcrootn_xfer_to_fire(p)     =  deadcrootn_xfer(p)    * f * cc_other(pft%itype(p)) 
        m_deadcrootn_storage_to_fire(p)  =  deadcrootn_storage(p) * f * cc_other(pft%itype(p))
        m_retransn_to_fire(p)            =  retransn(p)           * f * cc_other(pft%itype(p))

        ! mortality due to fire
        ! carbon bool 
        m_leafc_to_litter_fire(p)                   =  leafc(p) * f * &
             (1._r8 - cc_leaf(pft%itype(p))) * &
             fm_leaf(pft%itype(p))
        m_leafc_storage_to_litter_fire(p)           =  leafc_storage(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p))
        m_leafc_xfer_to_litter_fire(p)              =  leafc_xfer(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p))
        m_livestemc_to_litter_fire(p)               =  livestemc(p) * f * &
             (1._r8 - cc_lstem(pft%itype(p))) * &
             fm_droot(pft%itype(p))    
        m_livestemc_storage_to_litter_fire(p)       =  livestemc_storage(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p))
        m_livestemc_xfer_to_litter_fire(p)          =  livestemc_xfer(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p)) 
        m_livestemc_to_deadstemc_fire(p)            =  livestemc(p) * f * &
             (1._r8 - cc_lstem(pft%itype(p))) * &
             (fm_lstem(pft%itype(p))-fm_droot(pft%itype(p)))
        m_deadstemc_to_litter_fire(p)               =  deadstemc(p) * f * &
             (1._r8 - cc_dstem(pft%itype(p))) * &
             fm_droot(pft%itype(p))    
        m_deadstemc_storage_to_litter_fire(p)       =  deadstemc_storage(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p))
        m_deadstemc_xfer_to_litter_fire(p)          =  deadstemc_xfer(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p))
        m_frootc_to_litter_fire(p)                  =  frootc(p)             * f * &
             fm_root(pft%itype(p))
        m_frootc_storage_to_litter_fire(p)          =  frootc_storage(p)     * f * &
             fm_other(pft%itype(p))
        m_frootc_xfer_to_litter_fire(p)             =  frootc_xfer(p)        * f * &
             fm_other(pft%itype(p))
        m_livecrootc_to_litter_fire(p)              =  livecrootc(p)         * f * &
             fm_droot(pft%itype(p))
        m_livecrootc_storage_to_litter_fire(p)      =  livecrootc_storage(p) * f * &
             fm_other(pft%itype(p)) 
        m_livecrootc_xfer_to_litter_fire(p)         =  livecrootc_xfer(p)    * f * &
             fm_other(pft%itype(p)) 
        m_livecrootc_to_deadcrootc_fire(p)          =  livecrootc(p)         * f * &
             (fm_lroot(pft%itype(p))-fm_droot(pft%itype(p)))
        m_deadcrootc_to_litter_fire(p)              =  deadcrootc(p)         * f * &
             fm_droot(pft%itype(p))
        m_deadcrootc_storage_to_litter_fire(p)      =  deadcrootc_storage(p) * f * &
             fm_other(pft%itype(p))
        m_deadcrootc_xfer_to_litter_fire(p)         =  deadcrootc_xfer(p)    * f * &
             fm_other(pft%itype(p))      
        m_gresp_storage_to_litter_fire(p)           =  gresp_storage(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p))  
        m_gresp_xfer_to_litter_fire(p)              =  gresp_xfer(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p)) 



        ! nitrogen bool    
        m_leafn_to_litter_fire(p)                  =  leafn(p) * f * &
             (1._r8 - cc_leaf(pft%itype(p))) * &
             fm_leaf(pft%itype(p))
        m_leafn_storage_to_litter_fire(p)          =  leafn_storage(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p))  
        m_leafn_xfer_to_litter_fire(p)             =  leafn_xfer(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p))
        m_livestemn_to_litter_fire(p)              =  livestemn(p) * f * &
             (1._r8 - cc_lstem(pft%itype(p))) * &
             fm_droot(pft%itype(p))
        m_livestemn_storage_to_litter_fire(p)      =  livestemn_storage(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p))   
        m_livestemn_xfer_to_litter_fire(p)         =  livestemn_xfer(p) * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p))
        m_livestemn_to_deadstemn_fire(p)           =  livestemn(p) * f * &
             (1._r8 - cc_lstem(pft%itype(p))) * &
             (fm_lstem(pft%itype(p))-fm_droot(pft%itype(p)))
        m_frootn_to_litter_fire(p)                 =  frootn(p)             * f * &
             fm_root(pft%itype(p))
        m_frootn_storage_to_litter_fire(p)         =  frootn_storage(p)     * f * &
             fm_other(pft%itype(p))
        m_frootn_xfer_to_litter_fire(p)            =  frootn_xfer(p)        * f * &
             fm_other(pft%itype(p))
        m_livecrootn_to_litter_fire(p)             =  livecrootn(p)         * f * &
             fm_droot(pft%itype(p))
        m_livecrootn_storage_to_litter_fire(p)     =  livecrootn_storage(p) * f * &
             fm_other(pft%itype(p))
        m_livecrootn_xfer_to_litter_fire(p)        =  livecrootn_xfer(p)    * f * &
             fm_other(pft%itype(p)) 
        m_livecrootn_to_deadcrootn_fire(p)         =  livecrootn(p)         * f * &
             (fm_lroot(pft%itype(p))-fm_droot(pft%itype(p)))
        m_deadcrootn_to_litter_fire(p)             =  deadcrootn(p)         * f * &
             fm_droot(pft%itype(p))
        m_deadcrootn_storage_to_litter_fire(p)     =  deadcrootn_storage(p) * f * &
             fm_other(pft%itype(p))
        m_deadcrootn_xfer_to_litter_fire(p)        =  deadcrootn_xfer(p)    * f * &
             fm_other(pft%itype(p))
        m_retransn_to_litter_fire(p)               =  retransn(p)           * f * &
             (1._r8 - cc_other(pft%itype(p))) * &
             fm_other(pft%itype(p)) 

        if (use_cndv) then
           if ( woody(pft%itype(p)) == 1._r8 )then
              if ( livestemc(p)+deadstemc(p) > 0._r8 )then
                 nind(p) = nind(p)*(1._r8-1._r8*fm_droot(pft%itype(p))*f) 
              else
                 nind(p) = 0._r8
              end if
           end if
           leafcmax(p) = max(leafc(p)-m_leafc_to_fire(p)*dt, leafcmax(p))
           if (pft%itype(p) == noveg) leafcmax(p) = 0._r8
        end if

     end do  ! end of pfts loop  
     !
     ! fire-affected carbon to litter and cwd
     !
     do j = 1,nlevdecomp
        do pi = 1,max_pft_per_col
           do fc = 1,num_soilc
              c = filter_soilc(fc)
              if (pi <=  col%npfts(c)) then
                 p = col%pfti(c) + pi - 1
                 if ( pactive(p) ) then

                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_deadstemc_to_litter_fire(p) * pft%wtcol(p) * stem_prof(p,j)
                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_deadcrootc_to_litter_fire(p) * pft%wtcol(p) * croot_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_deadstemn_to_litter_fire(p) * pft%wtcol(p) * stem_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_deadcrootn_to_litter_fire(p) * pft%wtcol(p) * croot_prof(p,j)


                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_livestemc_to_litter_fire(p) * pft%wtcol(p) * stem_prof(p,j)
                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_livecrootc_to_litter_fire(p) * pft%wtcol(p) * croot_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_livestemn_to_litter_fire(p) * pft%wtcol(p) * stem_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_livecrootn_to_litter_fire(p) * pft%wtcol(p) * croot_prof(p,j)


                    m_c_to_litr_met_fire(c,j)=m_c_to_litr_met_fire(c,j) + &
                         ((m_leafc_to_litter_fire(p)*lf_flab(pft%itype(p)) &
                         +m_leafc_storage_to_litter_fire(p) + &
                         m_leafc_xfer_to_litter_fire(p) + &
                         m_gresp_storage_to_litter_fire(p) &
                         +m_gresp_xfer_to_litter_fire(p))*leaf_prof(p,j) + &
                         (m_frootc_to_litter_fire(p)*fr_flab(pft%itype(p)) &
                         +m_frootc_storage_to_litter_fire(p) + &
                         m_frootc_xfer_to_litter_fire(p))*froot_prof(p,j) &
                         +(m_livestemc_storage_to_litter_fire(p) + &
                         m_livestemc_xfer_to_litter_fire(p) &
                         +m_deadstemc_storage_to_litter_fire(p) + &
                         m_deadstemc_xfer_to_litter_fire(p))* stem_prof(p,j)&
                         +(m_livecrootc_storage_to_litter_fire(p) + &
                         m_livecrootc_xfer_to_litter_fire(p) &
                         +m_deadcrootc_storage_to_litter_fire(p) + &
                         m_deadcrootc_xfer_to_litter_fire(p))* croot_prof(p,j))* pft%wtcol(p)    
                    m_c_to_litr_cel_fire(c,j)=m_c_to_litr_cel_fire(c,j) + &
                         (m_leafc_to_litter_fire(p)*lf_fcel(pft%itype(p))*leaf_prof(p,j) + &
                         m_frootc_to_litter_fire(p)*fr_fcel(pft%itype(p))*froot_prof(p,j))* pft%wtcol(p) 
                    m_c_to_litr_lig_fire(c,j)=m_c_to_litr_lig_fire(c,j) + &
                         (m_leafc_to_litter_fire(p)*lf_flig(pft%itype(p))*leaf_prof(p,j) + &
                         m_frootc_to_litter_fire(p)*fr_flig(pft%itype(p))*froot_prof(p,j))* pft%wtcol(p)  

                    m_n_to_litr_met_fire(c,j)=m_n_to_litr_met_fire(c,j) + &
                         ((m_leafn_to_litter_fire(p)*lf_flab(pft%itype(p)) &
                         +m_leafn_storage_to_litter_fire(p) + &
                         m_leafn_xfer_to_litter_fire(p)+m_retransn_to_litter_fire(p)) &
                         *leaf_prof(p,j) +(m_frootn_to_litter_fire(p)*fr_flab(pft%itype(p)) &
                         +m_frootn_storage_to_litter_fire(p) + &
                         m_frootn_xfer_to_litter_fire(p))*froot_prof(p,j) &
                         +(m_livestemn_storage_to_litter_fire(p) + &
                         m_livestemn_xfer_to_litter_fire(p) &
                         +m_deadstemn_storage_to_litter_fire(p) + &
                         m_deadstemn_xfer_to_litter_fire(p))* stem_prof(p,j)&
                         +(m_livecrootn_storage_to_litter_fire(p) + &
                         m_livecrootn_xfer_to_litter_fire(p) &
                         +m_deadcrootn_storage_to_litter_fire(p) + &
                         m_deadcrootn_xfer_to_litter_fire(p))* croot_prof(p,j))* pft%wtcol(p)    
                    m_n_to_litr_cel_fire(c,j)=m_n_to_litr_cel_fire(c,j) + &
                         (m_leafn_to_litter_fire(p)*lf_fcel(pft%itype(p))*leaf_prof(p,j) + &
                         m_frootn_to_litter_fire(p)*fr_fcel(pft%itype(p))*froot_prof(p,j))* pft%wtcol(p) 
                    m_n_to_litr_lig_fire(c,j)=m_n_to_litr_lig_fire(c,j) + &
                         (m_leafn_to_litter_fire(p)*lf_flig(pft%itype(p))*leaf_prof(p,j) + &
                         m_frootn_to_litter_fire(p)*fr_flig(pft%itype(p))*froot_prof(p,j))* pft%wtcol(p) 
                 end if
              end if
           end do
        end do
     end do
     !
     ! vertically-resolved decomposing C/N fire loss   
     ! column loop
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)

        f = farea_burned(c) / dt

        ! apply this rate to the column state variables to get flux rates

        do j = 1, nlevdecomp
           ! carbon fluxes
           do l = 1, ndecomp_pools
              if ( is_litter(l) ) then
                 m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * f * 0.4_r8
              end if
              if ( is_cwd(l) ) then
                 m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * &
                      (f-baf_crop(c)/dt) * 0.2_r8
              end if
           end do

           ! nitrogen fluxes
           do l = 1, ndecomp_pools
              if ( is_litter(l) ) then
                 m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * f * 0.4_r8
              end if
              if ( is_cwd(l) ) then
                 m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * &
                      (f-baf_crop(c)/ dt) * 0.2_r8
              end if
           end do

        end do
     end do  ! end of column loop

     !
     ! carbon loss due to deforestation fires
     !
     if (fpftdyn /= ' ') then    !true when landuse data is used
        call get_curr_date (kyr, kmo, kda, mcsec)
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           lfc2(c)=0._r8
           if( .not. (kmo == 1 .and. kda == 1 .and. mcsec == 0) )then
              if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 .and. dtrotr_col(c) > 0._r8 .and. &
                   lfc(c) > 0._r8 .and. fbac1(c) == 0._r8) then
                 lfc2(c) = max(0._r8,min(lfc(c),(farea_burned(c)-baf_crop(c) - &
                      baf_peatf(c))/2.0))/(dtrotr_col(c)*dayspyr*secspday/dt)
                 lfc(c)  = lfc(c)-max(0._r8,min(lfc(c),(farea_burned(c)-baf_crop(c) - &
                      baf_peatf(c))/2.0_r8))
              end if
           end if
        end do
     end if
     !
     ! Carbon loss due to peat fires
     !
     ! somc_fire is not connected to clm45 soil carbon pool, ie does not decrease
     ! soil carbon b/c clm4 soil carbon was very low in peatland areas
     ! Fang Li has not checked clm45 soil carbon in peatland areas
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col%gridcell(c)
        if( grc%latdeg(g) .lt. borealat)then
           somc_fire(c)= totsomc(c)*baf_peatf(c)/dt*6.0_r8/33.9_r8
        else
           somc_fire(c)= baf_peatf(c)/dt*2.2e3_r8
        end if
     end do

     ! Fang Li has not added aerosol and trace gas emissions due to fire, yet
     ! They will be added here in proportion to the carbon emission
     ! Emission factors differ for various fire types

   end associate 
 end subroutine CNFireFluxes

 !-----------------------------------------------------------------------
 subroutine hdm_init( bounds )
   !
   ! !DESCRIPTION:
   !
   ! Initialize data stream information for population density.
   !
   ! !USES:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use ndepStreamMod    , only : clm_domain_mct
   use histFileMod      , only : hist_addfld1d
   !
   ! !ARGUMENTS:
   implicit none
   type(bounds_type), intent(in) :: bounds  ! bounds
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
        filename=(/trim(stream_fldFileName_popdens)/), &
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
  implicit none
  type(bounds_type), intent(in) :: bounds  ! bounds
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
  use clm_varctl       , only : inst_name
  use clm_time_manager , only : get_calendar
  use ncdio_pio        , only : pio_subsystem
  use shr_pio_mod      , only : shr_pio_getiotype
  use clm_nlUtilsMod   , only : find_nlgroup_name
  use ndepStreamMod    , only : clm_domain_mct
  use histFileMod      , only : hist_addfld1d
  !
  ! !ARGUMENTS:
  implicit none
  type(bounds_type), intent(in) :: bounds  ! bounds
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
  implicit none
  type(bounds_type), intent(in) :: bounds  ! bounds
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

end module CNFireMod
