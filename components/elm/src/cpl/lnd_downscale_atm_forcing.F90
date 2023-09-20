module lnd_downscale_atm_forcing

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Disaggregates gridcell quantities to topounits
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_megan_mod  , only : shr_megan_mechcomps_n
  use elm_varpar     , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use elm_varcon     , only : rair, grav, cpair, hfus, tfrz, spval
  use elm_varctl     , only : iulog, use_c13, use_cn, use_lch4, iulog, precip_downscaling_method
  use elm_cpl_indices
  use seq_drydep_mod , only : n_drydep, drydep_method, DD_XLND
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use atm2lndType    , only : atm2lnd_type
  use lnd2atmType    , only : lnd2atm_type
  use glc2lndMod     , only : glc2lnd_type
  use GridcellType   , only : grc_pp
  use TopounitType   , only : top_pp
  use TopounitDataType, only: top_as, top_es, top_af        
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp  
  use spmdmod          , only: masterproc
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: downscale_atm_forcing_to_topounit    ! Calls downscaling subroutines of forcing fields from gridcell to topounit
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: downscale_atm_state_to_topounit    ! Downscale atmosperic state fields from gridcell to topounit
  private :: downscale_longwave_to_topounit      ! Downscale longwave radiation field from gridcell to topounit
  private :: downscale_precip_to_topounit_FNM    ! Downscale precipitation field from gridcell to topounit using Froude number method (FNM)
  private :: downscale_precip_to_topounit_ERMM   ! Downscale precipitation field from gridcell to topounit using elevation ration with maximum elevation method (ERMM)
  private :: build_normalization                 ! Compute normalization factors so that downscaled fields are conservative
!  private :: check_downscale_consistency ! Check consistency of downscaling
  !-----------------------------------------------------------------------

contains

    !-----------------------------------------------------------------------
  subroutine downscale_atm_forcing_to_topounit(g, i, x2l, lnd2atm_vars)
    !
    ! !DESCRIPTION:
    ! Downscale fields from gridcell to topounit
    !
    ! Downscaling is done over topounits if the number of topounits > 1.
    !
    ! !USES:
    use elm_time_manager, only : get_nstep
    use elm_varcon      , only : rair, cpair, grav, lapse_glcmec
    use elm_varcon      , only : glcmec_rain_snow_threshold, o2_molar_const
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use landunit_varcon , only : istice_mec 
    use elm_varctl      , only : glcmec_downscale_rain_snow_convert
    use domainMod       , only : ldomain
    use QsatMod         , only : Qsat
    !
    ! !ARGUMENTS:
    integer                    , intent(in)    :: g  
    integer                    , intent(in)    :: i
    real(r8)                         , intent(in)    :: x2l(:,:)
    type(lnd2atm_type)               , intent(in)    :: lnd2atm_vars

    !
    ! !LOCAL VARIABLES:
    integer :: t, l, c, fc, t2         ! indices
    integer :: clo, cc
    integer :: numt_pg                ! Number of topounits per grid	
	
    ! temporaries for topo downscaling
    real(r8) :: rain_g, snow_g
    real(r8) :: mxElv              ! Maximum elevation value per grid    
    real(r8) :: uovern_t           ! Froude Number
    real(r8) :: grdElv             ! Grid elevation
    real(r8) :: topoElv            ! Topounit elevation
    real(r8) :: max_tpuElv         ! Maximum topounit elevation for calculating elevation range
    real(r8) :: min_tpuElv         ! Minimum topounit elevation for calculating elevation range
    real(r8) :: crnt_temp_t             ! Current downscaled topounit temperature
    real(r8) :: temp_r             ! Temporary topounit rainfall
    real(r8) :: temp_s             ! Temporary topounit snowfall
    real(r8) :: t_th               ! Temperature threshold for snowfall
    real(r8) :: Ta_th1              ! Temperature at 0.5 Celsius in K 273.65
    real(r8) :: Ta_th2              ! Temperature at 2.0 Celsius in K 275.15
    real(r8) :: Ta_th3              ! Temperature at 2.5 Celsius in K 275.65
    real(r8) :: tmp_Snow_frc        ! Current snow fraction

    real(r8) :: e
    real(r8) :: qvsat

    real(r8) :: sum_qbot_g    ! weighted sum of column-level lwrad
    real(r8) :: sum_wtsq_g      ! sum of weights that contribute to sum_lwrad_g
    real(r8) :: qbot_norm_g   ! normalization factors
    real(r8) :: sum_lwrad_g    ! weighted sum of column-level lwrad
    real(r8) :: sum_wtslw_g      ! sum of weights that contribute to sum_lwrad_g
    real(r8) :: lwrad_norm_g   ! normalization factors
    real(r8) :: esatw                ! saturation vapor pressure over water (Pa)
    real(r8) :: esati                ! saturation vapor pressure over ice (Pa)
    real(r8) :: a0,a1,a2,a3,a4,a5,a6 ! coefficients for esat over water
    real(r8) :: b0,b1,b2,b3,b4,b5,b6 ! coefficients for esat over ice
    real(r8) :: tdc, temp            ! Kelvins to Celcius function and its input
    real(r8) :: vp                   ! water vapor pressure (Pa)
    
    real(r8), allocatable :: deltaRain(:)      ! Deviation of subgrid rain from grid rain
    real(r8), allocatable :: deltaSnow(:)      ! Deviation of subgrid snow from grid snow
    real(r8) :: deltaR              ! Temporary deltaRain
    real(r8) :: deltaS              ! Temporary deltaSnow
    real(r8) :: sum_of_hrise        ! Sum of height rise of air parcel of all subgrids of a grid
    real(r8) :: hrise               ! Temporary height rise
    real(r8) :: elvrnge             ! Elevation range between lowest and highest tpu
    real(r8) :: ave_elv
    integer :: elv_flag             ! Elevation flag to trac grids with +ve grid elevation and -ve tpu elevation

    integer :: uaflag = 0
    integer :: precip_dwn = 0         ! Used to turn on/off the downscaling of precipitation  0 = on; 1 = off
    integer :: other_forcing_dwn = 0  ! Used to turn on/off the downscaling of other forcing  0 = on; 1 = off

    character(len=*), parameter :: subname = 'downscale_atm_forcing_to_topounit'
    !----------------------------------------------------------------------------------------
    
    ! Constants to compute vapor pressure
    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
         a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
         a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
         a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
         b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
         b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
         b6=1.838826904e-10_r8)

    !
    ! function declarations
    !
    tdc(temp) = min( 50._r8, max(-50._r8,(temp-SHR_CONST_TKFRZ)) )                       ! Taken from lnd_import_export.F90
    esatw(temp) = 100._r8*(a0+temp*(a1+temp*(a2+temp*(a3+temp*(a4+temp*(a5+temp*a6)))))) ! Taken from lnd_import_export.F90
    esati(temp) = 100._r8*(b0+temp*(b1+temp*(b2+temp*(b3+temp*(b4+temp*(b5+temp*b6)))))) ! Taken from lnd_import_export.F90
    !-----------------------------------------------------------------------
    ! Get required inputs
    numt_pg     = grc_pp%ntopounits(g)	         ! Number of topounits per grid
    grdElv      = grc_pp%elevation(g)            ! Grid level sfc elevation
    mxElv       = grc_pp%MaxElevation(g)         ! Maximum src elevation per grid obtained from the highest elevation topounit
    uovern_t    = x2l(index_x2l_Sa_uovern,i)     ! Froude Number
    
    snow_g      = x2l(index_x2l_Faxa_snowc,i) + x2l(index_x2l_Faxa_snowl,i)
    rain_g      = x2l(index_x2l_Faxa_rainc,i) + x2l(index_x2l_Faxa_rainl,i)
	
    sum_qbot_g = 0._r8
    sum_wtsq_g = 0._r8
    sum_lwrad_g = 0._r8
    sum_wtslw_g = 0._r8
    
    sum_of_hrise = 0._r8
    t_th = 273.15_r8   ! Freezing temperature in K
    Ta_th1 = 273.65_r8 ! Lowest threshold for snow calculation
    Ta_th2 = 275.15_r8 ! Middle threshold for rain/snow partitioning
    Ta_th3 = 275.65_r8 ! Highest threshold for rain/snow partitioning
    tmp_Snow_frc = 0._r8   ! Snow fraction
    elv_flag = 0
    elvrnge = 0._r8
    ave_elv = 0._r8
    max_tpuElv = 0._r8
    min_tpuElv = 0._r8
    if (numt_pg > 1) then          !downscaling is done only if a grid has more than 1 topounits 
       if (precip_downscaling_method == 'FNM') then
          allocate(deltaRain(numt_pg))
          deltaRain(:) = 0._r8
          allocate(deltaSnow(numt_pg))
          deltaSnow(:) = 0._r8
          hrise = 0._r8
       end if
       
       ! calculate elevation range and track grids with +ve elevation but have -ve tpu elevation
       max_tpuElv = -100000._r8
       min_tpuElv = 100000._r8
       do t = grc_pp%topi(g), grc_pp%topf(g)  ! Check occurence of grid elevation is +ve while tpu elevation is -ve
          topoElv  = top_pp%elevation(t)
           if (topoElv > max_tpuElv) then 
              max_tpuElv = topoElv
           end if
           if (topoElv < min_tpuElv) then 
              min_tpuElv = topoElv
           end if
       end do
       elvrnge = max_tpuElv - min_tpuElv
       
       do t = grc_pp%topi(g), grc_pp%topf(g)       
          t2 = t - grc_pp%topi(g) + 1
          topoElv  = top_pp%elevation(t)             ! Topounit sfc elevation  
          
          ! Downscale precipitation
          if (mxElv == 0._r8 .or. precip_dwn == 1) then  ! avoid dividing by 0
             top_af%rain(t) = rain_g
             top_af%snow(t) = snow_g
          else
             if (precip_downscaling_method == 'FNM') then
                call downscale_precip_to_topounit_FNM(mxElv,uovern_t,grdElv,topoElv,rain_g,snow_g,deltaR,deltaS,hrise) !Use FNM method
                deltaRain(t2) = deltaR 
                deltaSnow(t2) = deltaS             
                sum_of_hrise = sum_of_hrise + hrise              
             else
                call downscale_precip_to_topounit_ERMM(t,mxElv,grdElv,topoElv,rain_g, snow_g,elv_flag,elvrnge) ! Use ERMM method                
             end if
          end if 

          ! Downscale other fluxes
          if (other_forcing_dwn == 1) then  ! flag to turn on or off downscaling of other forcing
             top_af%rain(t) = rain_g 
             top_af%snow(t) = snow_g 
             top_af%lwrad(t) = x2l(index_x2l_faxa_lwdn,i)       	  
             ! Update top_as
             top_as%tbot(t)    = x2l(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
             top_as%thbot(t)   = x2l(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
             top_as%pbot(t)    = x2l(index_x2l_Sa_pbot,i)      ! ptcmxy    Atm state Pa
             top_as%qbot(t)    = x2l(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
             top_as%ubot(t)    = x2l(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
             top_as%vbot(t)    = x2l(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
             top_as%zbot(t)    = x2l(index_x2l_Sa_z,i)         ! zgcmxy    Atm state m
		 
             ! assign the state forcing fields derived from other inputs
             ! Horizontal windspeed (m/s)
             top_as%windbot(t) = sqrt(top_as%ubot(t)**2 + top_as%vbot(t)**2)
             ! Relative humidity (percent)
             if (top_as%tbot(t) > SHR_CONST_TKFRZ) then
                e = esatw(tdc(top_as%tbot(t)))
             else
                e = esati(tdc(top_as%tbot(t)))
             end if
             qvsat = 0.622_r8*e / (top_as%pbot(t) - 0.378_r8*e)
             top_as%rhbot(t) = 100.0_r8*(top_as%qbot(t) / qvsat)
             ! partial pressure of oxygen (Pa)
             top_as%po2bot(t) = o2_molar_const * top_as%pbot(t)
             ! air density (kg/m**3) - uses a temporary calculation
             ! of water vapor pressure (Pa)
             vp = top_as%qbot(t) * top_as%pbot(t)  / (0.622_r8 + 0.378_r8 * top_as%qbot(t))
             top_as%rhobot(t) = (top_as%pbot(t) - 0.378_r8 * vp) / (rair * top_as%tbot(t))

             top_af%solad(t,2) = x2l(index_x2l_Faxa_swndr,i)
             top_af%solad(t,1) = x2l(index_x2l_Faxa_swvdr,i)
             top_af%solai(t,2) = x2l(index_x2l_Faxa_swndf,i)
             top_af%solai(t,1) = x2l(index_x2l_Faxa_swvdf,i)
             ! derived flux forcings
             top_af%solar(t) = top_af%solad(t,2) + top_af%solad(t,1) + &
       		       top_af%solai(t,2) + top_af%solai(t,1)
          else
             call downscale_atm_state_to_topounit(g, i, t, x2l, lnd2atm_vars, uaflag)
             call downscale_longwave_to_topounit(g, i, t, x2l, lnd2atm_vars, uaflag)
             top_as%ubot(t)    = x2l(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
             top_as%vbot(t)    = x2l(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
             top_as%zbot(t)    = x2l(index_x2l_Sa_z,i)         ! zgcmxy    Atm state m
          
             sum_qbot_g = sum_qbot_g + top_pp%wtgcell(t)*top_as%qbot(t)
             sum_wtsq_g = sum_wtsq_g + top_pp%wtgcell(t)
                
             ! assign the state forcing fields derived from other inputs
             ! Horizontal windspeed (m/s)
             top_as%windbot(t) = sqrt(top_as%ubot(t)**2 + top_as%vbot(t)**2)
             ! partial pressure of oxygen (Pa)
             top_as%po2bot(t) = o2_molar_const * top_as%pbot(t)
             ! air density (kg/m**3) - uses a temporary calculation
             ! of water vapor pressure (Pa)
             vp = top_as%qbot(t) * top_as%pbot(t)  / (0.622_r8 + 0.378_r8 * top_as%qbot(t))
             top_as%rhobot(t) = (top_as%pbot(t) - 0.378_r8 * vp) / (rair * top_as%tbot(t))
             top_af%solad(t,2) = x2l(index_x2l_Faxa_swndr,i)
             top_af%solad(t,1) = x2l(index_x2l_Faxa_swvdr,i)
             top_af%solai(t,2) = x2l(index_x2l_Faxa_swndf,i)
             top_af%solai(t,1) = x2l(index_x2l_Faxa_swvdf,i)
            ! derived flux forcings
             top_af%solar(t) = top_af%solad(t,2) + top_af%solad(t,1) + &
             top_af%solai(t,2) + top_af%solai(t,1)

            ! Keep track of the gridcell-level weighted sum for later normalization.
            !
            ! This gridcell-level weighted sum just includes points for which we do the
            ! downscaling (e.g., glc_mec points). Thus the contributing weights
            ! generally do not add to 1. So to do the normalization properly, we also
            ! need to keep track of the weights that have contributed to this sum.
             sum_lwrad_g = sum_lwrad_g + top_pp%wtgcell(t)*top_af%lwrad(t)
             sum_wtslw_g = sum_wtslw_g + top_pp%wtgcell(t)
          end if 
       end do
       if (precip_downscaling_method == 'FNM') then
          do t = grc_pp%topi(g), grc_pp%topf(g)
             t2 = t - grc_pp%topi(g) + 1
             if (mxElv == 0.) then  ! avoid dividing by 0
                top_af%rain(t) = rain_g
                top_af%snow(t) = snow_g
             else
                top_af%rain(t) = rain_g + (deltaRain(t2) - (rain_g/mxElv)*(sum_of_hrise/numt_pg))
                top_af%snow(t) = snow_g + (deltaSnow(t2) - (snow_g/mxElv)*(sum_of_hrise/numt_pg))
             end if

          end do
          deallocate(deltaRain)
          deallocate(deltaSnow)
       end if
       
	   ! Precipitation partitioning using simple method following Jordan (1991)
       do t = grc_pp%topi(g), grc_pp%topf(g)
            crnt_temp_t = top_as%tbot(t)
            if (crnt_temp_t > Ta_th3) then            ! No snow or all rain
                temp_r = top_af%rain(t) + top_af%snow(t)
                temp_s = 0._r8
            else if (crnt_temp_t >= Ta_th2 .and. crnt_temp_t <= Ta_th3) then
                temp_s = (top_af%snow(t) + top_af%rain(t))*0.6_r8              ! 0.6 fraction of the total precip is snow
                temp_r = (top_af%snow(t) + top_af%rain(t)) - temp_s
            else if (crnt_temp_t < Ta_th2 .and. crnt_temp_t > Ta_th1) then
                tmp_Snow_frc = (crnt_temp_t-Ta_th2)*(0.4_r8/(Ta_th1-Ta_th2))+0.6_r8  ! Snow fraction value 1-0.6 = 0.4 snowfrc is 1 at t = Ta_th1
                temp_s = (top_af%snow(t) + top_af%rain(t))*tmp_Snow_frc        ! 0.6 is snow fraction at t = Ta_th2  
                temp_r = (top_af%snow(t) + top_af%rain(t)) - temp_s
            else   ! crnt_temp_t <= Ta_th1 ==> all snow
                temp_s = top_af%snow(t) + top_af%rain(t)
                temp_r = 0._r8                
            
            end if
            
            top_af%rain(t) = temp_r
            top_af%snow(t) = temp_s
            
       end do

    else !grid has a single topounit
       ! update top_af using grid level values
       t = grc_pp%topi(g)
       top_af%rain(t) = rain_g 
       top_af%snow(t) = snow_g 
       top_af%lwrad(t) = x2l(index_x2l_faxa_lwdn,i)         
       ! Update top_as
       top_as%tbot(t)    = x2l(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
       top_as%thbot(t)   = x2l(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
       top_as%pbot(t)    = x2l(index_x2l_Sa_pbot,i)      ! ptcmxy    Atm state Pa
       top_as%qbot(t)    = x2l(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
       top_as%ubot(t)    = x2l(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
       top_as%vbot(t)    = x2l(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
       top_as%zbot(t)    = x2l(index_x2l_Sa_z,i)         ! zgcmxy    Atm state m
		 
       ! assign the state forcing fields derived from other inputs
       ! Horizontal windspeed (m/s)
       top_as%windbot(t) = sqrt(top_as%ubot(t)**2 + top_as%vbot(t)**2)
       ! Relative humidity (percent)
       if (top_as%tbot(t) > SHR_CONST_TKFRZ) then
          e = esatw(tdc(top_as%tbot(t)))
       else
          e = esati(tdc(top_as%tbot(t)))
       end if
       qvsat = 0.622_r8*e / (top_as%pbot(t) - 0.378_r8*e)
       top_as%rhbot(t) = 100.0_r8*(top_as%qbot(t) / qvsat)
       ! partial pressure of oxygen (Pa)
       top_as%po2bot(t) = o2_molar_const * top_as%pbot(t)
       ! air density (kg/m**3) - uses a temporary calculation
       ! of water vapor pressure (Pa)
       vp = top_as%qbot(t) * top_as%pbot(t)  / (0.622_r8 + 0.378_r8 * top_as%qbot(t))
       top_as%rhobot(t) = (top_as%pbot(t) - 0.378_r8 * vp) / (rair * top_as%tbot(t))

       top_af%solad(t,2) = x2l(index_x2l_Faxa_swndr,i)
       top_af%solad(t,1) = x2l(index_x2l_Faxa_swvdr,i)
       top_af%solai(t,2) = x2l(index_x2l_Faxa_swndf,i)
       top_af%solai(t,1) = x2l(index_x2l_Faxa_swvdf,i)
       ! derived flux forcings
       top_af%solar(t) = top_af%solad(t,2) + top_af%solad(t,1) + &
       		       top_af%solai(t,2) + top_af%solai(t,1)
 
    end if   
		 
    if (numt_pg > 1) then

          ! Normalize forc_lwrad_c(c) to conserve energy

          call build_normalization(orig_field=x2l(index_x2l_Sa_shum,i), &
              sum_field=sum_qbot_g, sum_wts=sum_wtsq_g, norms=qbot_norm_g)

          do t = grc_pp%topi(g), grc_pp%topf(g)
             top_as%qbot(t) = top_as%qbot(t) * qbot_norm_g

             ! Relative humidity (percent)
             if (top_as%tbot(t) > SHR_CONST_TKFRZ) then
                e = esatw(tdc(top_as%tbot(t)))
             else
                e = esati(tdc(top_as%tbot(t)))
             end if
             qvsat = 0.622_r8*e / (top_as%pbot(t) - 0.378_r8*e)
             top_as%rhbot(t) = 100.0_r8*(top_as%qbot(t) / qvsat)
             ! partial pressure of oxygen (Pa)
             top_as%po2bot(t) = o2_molar_const * top_as%pbot(t)
             ! air density (kg/m**3) - uses a temporary calculation
             ! of water vapor pressure (Pa)
             vp = top_as%qbot(t) * top_as%pbot(t)  / (0.622_r8 + 0.378_r8 * top_as%qbot(t))
             top_as%rhobot(t) = (top_as%pbot(t) - 0.378_r8 * vp) / (rair * top_as%tbot(t))

          end do


       ! Normalize forc_lwrad_c(c) to conserve energy

       call build_normalization(orig_field=x2l(index_x2l_faxa_lwdn,i), &
              sum_field=sum_lwrad_g, sum_wts=sum_wtslw_g, norms=lwrad_norm_g)

       do t = grc_pp%topi(g), grc_pp%topf(g)
          top_af%lwrad(t) = top_af%lwrad(t) * lwrad_norm_g
       end do

    end if

  end subroutine downscale_atm_forcing_to_topounit
  
  !-----------------------------------------------------------------------
  ! Downscale other atmospheric state variables
  !-----------------------------------------------------------------------
  subroutine downscale_atm_state_to_topounit(g, i, t, x2l, lnd2atm_vars, method)
    !
    ! !DESCRIPTION:
    ! Downscale atmospheric forcing fields from gridcell to topounit
    !
    ! Downscaling is done over topounits.
    !
    ! !USES:
    use elm_time_manager, only : get_nstep
    use elm_varcon      , only : rair, cpair, grav, lapse_glcmec
    use elm_varcon      , only : glcmec_rain_snow_threshold
    use landunit_varcon , only : istice_mec 
    use elm_varctl      , only : glcmec_downscale_rain_snow_convert
    use domainMod       , only : ldomain
    use QsatMod         , only : Qsat
    !
    ! !ARGUMENTS:
    integer                    , intent(in)    :: g  
    integer                    , intent(in)    :: i
    integer                    , intent(in)    :: t
    real(r8)                         , intent(in)    :: x2l(:,:)
    type(lnd2atm_type)               , intent(in)    :: lnd2atm_vars
    integer                    , intent(in)    :: method
    !
    ! !LOCAL VARIABLES:
    integer :: l, c, fc         ! indices
    integer :: clo, cc
    integer :: nstep

    ! temporaries for topo downscaling
    real(r8) :: hsurf_g,hsurf_t,Hbot
    real(r8) :: zbot_g, tsfc_g, tbot_g, pbot_g, thbot_g, qbot_g, qs_g, es_g
    real(r8) :: zbot_t, tsfc_t, tbot_t, pbot_t, thbot_t, qbot_t, qs_t, es_t
    real(r8) :: egcm_t, rhos_t
    real(r8) :: dum1,   dum2

    character(len=*), parameter :: subname = 'downscale_atm_state_to_topounit'
    !-----------------------------------------------------------------------

    nstep = get_nstep()
    
    ! Downscale forc_t, forc_th, forc_q, forc_pbot, and forc_rho to columns.
    ! For glacier_mec columns the downscaling is based on surface elevation.
    ! For other columns the downscaling is a simple copy (above).

         ! This is a simple downscaling procedure 
         ! Note that forc_hgt, forc_u, and forc_v are not downscaled.

    hsurf_g = grc_pp%elevation(g)  ! gridcell sfc elevation
    hsurf_t = top_pp%elevation(t)                  ! topounit sfc elevation
    tbot_g  = x2l(index_x2l_Sa_tbot,i)              ! atm sfc temp
    thbot_g = x2l(index_x2l_Sa_ptem,i)              ! atm sfc pot temp
    tsfc_g  = lnd2atm_vars%t_rad_grc(g)             ! sfc rad temp
    qbot_g  = x2l(index_x2l_Sa_shum,i)              ! atm sfc spec humid
    pbot_g  = x2l(index_x2l_Sa_pbot,i)              ! atm sfc pressure
    zbot_g  = x2l(index_x2l_Sa_z,i)                 ! atm ref height

    zbot_t  = zbot_g
    tsfc_t = top_es%t_rad(t)
    if (method == 0 .or. (method == 1 .and. nstep == 0)) then
       tbot_t  = tbot_g-lapse_glcmec*(hsurf_t-hsurf_g) ! sfc temp for column
    else            
       tbot_t = tsfc_t + tbot_g - tsfc_g ! tsfc is from previous time step              
    end if
 
    Hbot    = rair*0.5_r8*(tbot_g+tbot_t)/grav      ! scale ht at avg temp
    pbot_t  = pbot_g*exp(-(hsurf_t-hsurf_g)/Hbot)   ! column sfc press

    ! Derivation of potential temperature calculation:
    ! 
    ! The textbook definition would be:
    ! thbot_c = tbot_c * (p0/pbot_c)^(rair/cpair)
    ! 
    ! Note that pressure is related to scale height as:
    ! pbot_c = p0 * exp(-zbot_c/H)
    !
    ! Using Hbot in place of H, we get:
    ! pbot_c = p0 * exp(-zbot_c/Hbot)
    !
    ! Plugging this in to the textbook definition, then manipulating, we get:
    ! thbot_c = tbot_c * (p0/(p0*exp(-zbot_c/Hbot)))^(rair/cpair)
    !         = tbot_c * (1/exp(-zbot_c/Hbot))^(rair/cpair)
    !         = tbot_c * (exp(zbot_c/Hbot))^(rair/cpair)
    !         = tbot_c * exp((zbot_c/Hbot) * (rair/cpair))

    thbot_t= tbot_t*exp((zbot_t/Hbot)*(rair/cpair))  ! pot temp calc

    call Qsat(tbot_g,pbot_g,es_g,dum1,qs_g,dum2)
    call Qsat(tbot_t,pbot_t,es_t,dum1,qs_t,dum2)

    qbot_t = qbot_g*(qs_t/qs_g)
    egcm_t = qbot_t*pbot_t/(0.622_r8+0.378_r8*qbot_t)
    rhos_t = (pbot_t-0.378_r8*egcm_t) / (rair*tbot_t)

    top_as%tbot(t) = tbot_t
    top_as%thbot(t) = thbot_t
    top_as%qbot(t) = qbot_t
    top_as%pbot(t) = pbot_t

!    call check_downscale_consistency(bounds, atm2lnd_vars)

  end subroutine downscale_atm_state_to_topounit
    
  !-------------------------------------------------------
  ! Downscale longwave radiation  place holder
  subroutine downscale_longwave_to_topounit(g, i, t, x2l, lnd2atm_vars, method)
      
    ! !DESCRIPTION:
    ! Downscale longwave radiation from gridcell to column
    ! Must be done AFTER temperature downscaling
    
    ! !USES:
    use elm_time_manager, only : get_nstep
    use domainMod       , only : ldomain
    use landunit_varcon , only : istice_mec
    use elm_varcon      , only : lapse_glcmec
    use elm_varctl      , only : glcmec_downscale_longwave
    
    ! !ARGUMENTS:
    integer                    , intent(in)    :: g
    integer                    , intent(in)    :: i
    integer                    , intent(in)    :: t
    real(r8)                         , intent(in)    :: x2l(:,:)
    type(lnd2atm_type)               , intent(in)    :: lnd2atm_vars
    integer                    , intent(in)    :: method
    
    ! !LOCAL VARIABLES:
    integer  :: c,l,fc     ! indices
    integer  :: nstep
    real(r8) :: hsurf_t      ! column-level elevation (m)
    real(r8) :: hsurf_g      ! gridcell-level elevation (m)

    real(r8) :: tair_g         ! original gridcell mean air temperature
    real(r8) :: tair_t         ! downscaled topounit air temperature
    real(r8) :: tsfc_g         ! original gridcell surface temperature
    real(r8) :: tsfc_t         ! downscaled topounit surface temperature
    real(r8) :: lwrad_g        ! original gridcell mean LW radiation
    real(r8) :: lwrad_t        ! downscaled topounit LW radiation
    real(r8) :: newsum_lwrad_g ! weighted sum of column-level lwrad after normalization

    character(len=*), parameter :: subname = 'downscale_longwave_to_topounit'
    !-----------------------------------------------------------------------

    nstep = get_nstep()
    
    ! Do the downscaling
    hsurf_g = grc_pp%elevation(g) 
    hsurf_t = top_pp%elevation(t)

    ! Here we assume that deltaLW = (dLW/dT)*(dT/dz)*deltaz
    ! We get dLW/dT = 4*eps*sigma*T^3 = 4*LW/T from the Stefan-Boltzmann law,
    ! evaluated at the mean temp.
    ! We assume the same temperature lapse rate as above.

    tair_g = x2l(index_x2l_Sa_tbot,i)
    tair_t = top_as%tbot(t)
    tsfc_g = lnd2atm_vars%t_rad_grc(g)
    lwrad_g = x2l(index_x2l_Faxa_lwdn,i)
    if (method == 0 .or. (method == 1 .and. nstep == 0)) then
       lwrad_t = lwrad_g - &
          4.0_r8 * lwrad_g/(0.5_r8*(tair_t+tair_g)) * &
          lapse_glcmec * (hsurf_t - hsurf_g)
    else
       
       lwrad_t = lwrad_g + lnd2atm_vars%eflx_lwrad_out_grc(g) * &
          4._r8 * (tair_t - tair_g) / tsfc_g
       
    end if
       top_af%lwrad(t) = lwrad_t

  end subroutine downscale_longwave_to_topounit
  
  !-------------------------------------------------------
  ! The Froude Number Method (FNM) utilizes information on height rise 
  ! (equivalent to Froude Number), which is calculated based on wind speed and static stability, over 
  ! topographically heterogenous grids whenever a stable orographic regime is identified within the target grid 
  ! and uses elevation range when the orographic regime is not stable.
  ! It is can be turned on for coupled E3SM model configurations only.
  subroutine downscale_precip_to_topounit_FNM(mxElv,uovern_t,grdElv,topoElv,rain_g,snow_g,deltaR,deltaS,hrise) 
   
    ! !ARGUMENTS:
    real(r8) :: rain_g, snow_g
    real(r8) :: mxElv                    ! Maximum elevation value per grid
    real(r8) :: uovern_t                 ! Froude Number
    real(r8) :: grdElv                   ! Grid elevation
    real(r8) :: topoElv                  ! Topounit elevation
    real(r8), intent(inout) :: deltaR
    real(r8), intent(inout) :: deltaS
    real(r8), intent(inout) :: hrise    
    
    ! !LOCAL VARIABLES:
    real(r8) :: r
    real(r8) :: topoEl
    real(r8) :: mxEl 
    real(r8) :: grdEl
    character(len=*), parameter :: subname = 'downscale_precip_to_topounit_FNM'
    !--------------------------------------------------------------------------------------
    
    if (mxElv < 0) then
       mxEl = abs(mxElv) ! avoid situations where -ve r can force lower areas more rain than higher elevation areas.
       topoEl = abs(topoElv)
       grdEl = abs(grdElv)
    else 
       mxEl = mxElv ! avoid situations where -ve r can force lower areas more rain than higher elevation areas.
       topoEl = topoElv
       grdEl = grdElv
    end if    
    
    r = 0._r8  
    r = topoEl - grdEl
    if (uovern_t < 0) then
       hrise = r
    else
       hrise = min(uovern_t,r)
    end if    
    
    if (rain_g < 0) then
       deltaR = 0._r8
    else          
       deltaR = (rain_g*(hrise/mxEl))
    end if

    if (snow_g < 0) then
       deltaS = 0._r8 
    else
       deltaS = (snow_g*(hrise/mxEl))
    end if
  end subroutine downscale_precip_to_topounit_FNM  
  
  !-------------------------------------------------------
  ! Elevation Range with Maximum elevation Method (ERMM) utilizes only topographic characteristics 
  ! of the grid and topounits to disaggregate grid-level precipitation into the topounits of the grid. 
  ! It is can be turned on in both offline and coupled E3SM model configurations.
  subroutine downscale_precip_to_topounit_ERMM(t,mxElv,grdElv,topoElv,rain_g, snow_g,elv_flag,elvrnge) 
   
    ! !ARGUMENTS:
    real(r8) :: rain_g, snow_g
    real(r8) :: mxElv              ! Maximum elevation value per grid
    real(r8) :: grdElv             ! Grid elevation
    real(r8) :: topoElv            ! Topounit elevation    
    real(r8) :: elvrnge            ! Elevation range 
    integer :: elv_flag
    integer, intent(in)    :: t    
    
    ! !LOCAL VARIABLES:
    real(r8) :: r
    real(r8) :: topoEl
    real(r8) :: mxEl 
    real(r8) :: grdEl
    
    character(len=*), parameter :: subname = 'downscale_precip_to_topounit_ERMM'
    !--------------------------------------------------------------------------------------

    topoEl = topoElv
    grdEl = grdElv
    r = 0._r8  
    r = topoEl - grdEl
    if (mxElv <= 0._r8 .or. mxElv < elvrnge) then  ! Use elvrnge to downscale grid precipitation across the topounits
       if (rain_g < 0._r8) then
           write(iulog,*) 'WARNING negative grid level rainfall = ', rain_g
           top_af%rain(t) = 0._r8
        else if (elvrnge == 0._r8) then  !flat topography
           top_af%rain(t) = rain_g
        else
           top_af%rain(t) = rain_g + (rain_g*(r/elvrnge))
        end if
        if (snow_g < 0._r8) then
           write(iulog,*) 'WARNING negative grid level snowfall = ', snow_g
           top_af%snow(t) = 0._r8
        else if (elvrnge == 0._r8) then  !flat topography
           top_af%snow(t) = snow_g
        else
           top_af%snow(t) = snow_g + (snow_g*(r/elvrnge))
        end if    
    else  ! mxElv > 0.0
        top_af%rain(t) = rain_g + (rain_g*(r/mxElv))
        top_af%snow(t) = snow_g + (snow_g*(r/mxElv))
    end if 

    if (top_af%rain(t) < 0._r8 .or. top_af%snow(t) < 0._r8) then 
        write(iulog,*) 'WARNING negative tpu precipitation: top_af%rain(t) + top_af%snow(t) = ', (top_af%rain(t) + top_af%snow(t)) 
        write(iulog,*) ' top_af%rain(t) = ', top_af%rain(t)
        write(iulog,*) ' top_af%snow(t) = ', top_af%snow(t)
        write(iulog,*) ' mxElv =  ', mxElv
        write(iulog,*) ' grdEl =  ', grdEl
        write(iulog,*) ' topoEl =  ', topoEl
        write(iulog,*) ' elvrnge =  ', elvrnge
        write(iulog,*) ' r =  ', r
        write(iulog,*) ' rain_g =  ', rain_g
        write(iulog,*) ' snow_g =  ', snow_g
        write(iulog,*) ' (rain_g*(r/mxElv)) =  ', (rain_g*(r/mxElv))
        write(iulog,*) ' (snow_g*(r/mxElv)) =  ', (snow_g*(r/mxElv))         
     end if

  end subroutine downscale_precip_to_topounit_ERMM  
  
 !-----------------------------------------------------------------------
  subroutine build_normalization(orig_field, sum_field, sum_wts, norms)
    !
    ! !DESCRIPTION:
    ! Build an array of normalization factors that can be applied to a downscaled forcing
    ! field, in order to force the mean of the new field to be the same as the mean of
    ! the old field (for conservation).
    !
    ! This allows for the possibility that only a subset of columns are downscaled. Only
    ! the columns that are adjusted should be included in the weighted sum, sum_field;
    ! sum_wts gives the sum of contributing weights on the grid cell level. 

    ! For example, if a grid cell has an original forcing value of 1.0, and contains 4
    ! columns with the following weights on the gridcell, and the following values after
    ! normalization:
    !
    !       col #:    1     2     3     4
    !      weight:  0.1   0.2   0.3   0.4
    ! downscaled?:  yes   yes    no    no
    !       value:  0.9   1.1   1.0   1.0
    !
    ! Then we would have:
    ! orig_field(g) = 1.0
    ! sum_field(g) = 0.1*0.9 + 0.2*1.1 = 0.31
    ! sum_wts(g) = 0.1 + 0.2 = 0.3
    ! norms(g) = 1.0 / (0.31 / 0.3) = 0.9677
    !
    ! The field can then be normalized as:
    !              forc_lwrad_c(c) = forc_lwrad_c(c) * lwrad_norm_g(g)
    !   where lwrad_norm_g is the array of norms computed by this routine

    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: orig_field  ! the original field, at the grid cell level
    real(r8), intent(in)  :: sum_field   ! the new weighted sum across columns (dimensioned by grid cell)
    real(r8), intent(in)  :: sum_wts     ! sum of the weights used to create sum_field (dimensioned by grid cell)
    real(r8), intent(out) :: norms       ! computed normalization factors
    !-----------------------------------------------------------------------

    if (sum_wts == 0._r8) then
       ! Avoid divide by zero; if sum_wts is 0, then the normalization doesn't matter,
       ! because the adjusted values won't affect the grid cell mean.
       norms = 1.0_r8

    else if (sum_field == 0._r8) then
       ! Avoid divide by zero; this should only happen if the gridcell-level value is 0,
       ! in which case the normalization doesn't matter
       norms = 1.0_r8

    else
       ! The standard case
       norms = orig_field / (sum_field / sum_wts)

    end if

  end subroutine build_normalization

end module lnd_downscale_atm_forcing
