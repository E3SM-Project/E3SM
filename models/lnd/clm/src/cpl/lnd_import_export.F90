module lnd_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use abortutils   , only: endrun
  use decompmod    , only: bounds_type
  use lnd2atmType  , only: lnd2atm_type
  use lnd2glcMod   , only: lnd2glc_type
  use atm2lndType  , only: atm2lnd_type
  use glc2lndMod   , only: glc2lnd_type 
  use clm_cpl_indices
  !
  implicit none
  !===============================================================================

contains

  !===============================================================================
  subroutine lnd_import( bounds, x2l, atm2lnd_vars, glc2lnd_vars)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the coupler to the land model 
    !
    ! !USES:
    use clm_varctl      , only: co2_type, co2_ppmv, iulog, use_c13, create_glacier_mec_landunit
    use clm_varcon      , only: rair, o2_molar_const, c13ratio
    use shr_const_mod   , only: SHR_CONST_TKFRZ
    use domainMod       , only: ldomain
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    real(r8)           , intent(in)    :: x2l(:,:) ! driver import state to land model
    type(atm2lnd_type) , intent(inout) :: atm2lnd_vars      ! clm internal input data type
    type(glc2lnd_type) , intent(inout) :: glc2lnd_vars      ! clm internal input data type
    !
    ! !LOCAL VARIABLES:
    integer  :: g,i,nstep,ier        ! indices, number of steps, and error code
    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: e                    ! vapor pressure (Pa)
    real(r8) :: qsat                 ! saturation specific humidity (kg/kg)
    real(r8) :: forc_t               ! atmospheric temperature (Kelvin)
    real(r8) :: forc_q               ! atmospheric specific humidity (kg/kg)
    real(r8) :: forc_pbot            ! atmospheric pressure (Pa)
    real(r8) :: forc_rainl           ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc           ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl           ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag        ! temporary
    real(r8) :: co2_ppmv_prog        ! temporary
    real(r8) :: co2_ppmv_val         ! temporary
    integer  :: co2_type_idx         ! integer flag for co2_type options
    real(r8) :: esatw                ! saturation vapor pressure over water (Pa)
    real(r8) :: esati                ! saturation vapor pressure over ice (Pa)
    real(r8) :: a0,a1,a2,a3,a4,a5,a6 ! coefficients for esat over water
    real(r8) :: b0,b1,b2,b3,b4,b5,b6 ! coefficients for esat over ice
    real(r8) :: tdc, t               ! Kelvins to Celcius function and its input
    integer  :: num                  ! counter
    character(len=32), parameter :: sub = 'lnd_import_mct'

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
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
    !---------------------------------------------------------------------------

    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic' )
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic' )
    end if

    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.


    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)

       ! Determine flooding input, sign convention is positive downward and
       ! hierarchy is atm/glc/lnd/rof/ice/ocn.  so water sent from rof to land is negative,
       ! change the sign to indicate addition of water to system.

       atm2lnd_vars%forc_flood_grc(g)   = -x2l(index_x2l_Flrr_flood,i)  

       atm2lnd_vars%volr_grc(g)   = x2l(index_x2l_Flrr_volr,i) * (ldomain%area(g) * 1.e6_r8)

       ! Determine required receive fields

       atm2lnd_vars%forc_hgt_grc(g)                  = x2l(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
       atm2lnd_vars%forc_u_grc(g)                    = x2l(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
       atm2lnd_vars%forc_v_grc(g)                    = x2l(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
       atm2lnd_vars%forc_solad_grc(g,2)              = x2l(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
       atm2lnd_vars%forc_solad_grc(g,1)              = x2l(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
       atm2lnd_vars%forc_solai_grc(g,2)              = x2l(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
       atm2lnd_vars%forc_solai_grc(g,1)              = x2l(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

       atm2lnd_vars%forc_th_not_downscaled_grc(g)    = x2l(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
       atm2lnd_vars%forc_q_not_downscaled_grc(g)     = x2l(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
       atm2lnd_vars%forc_pbot_not_downscaled_grc(g)  = x2l(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
       atm2lnd_vars%forc_t_not_downscaled_grc(g)     = x2l(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
       atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = x2l(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2

       forc_rainc                                    = x2l(index_x2l_Faxa_rainc,i)   ! mm/s
       forc_rainl                                    = x2l(index_x2l_Faxa_rainl,i)   ! mm/s
       forc_snowc                                    = x2l(index_x2l_Faxa_snowc,i)   ! mm/s
       forc_snowl                                    = x2l(index_x2l_Faxa_snowl,i)   ! mm/s

       ! atmosphere coupling, for prognostic/prescribed aerosols
       atm2lnd_vars%forc_aer_grc(g,1)                = x2l(index_x2l_Faxa_bcphidry,i)
       atm2lnd_vars%forc_aer_grc(g,2)                = x2l(index_x2l_Faxa_bcphodry,i)
       atm2lnd_vars%forc_aer_grc(g,3)                = x2l(index_x2l_Faxa_bcphiwet,i)
       atm2lnd_vars%forc_aer_grc(g,4)                = x2l(index_x2l_Faxa_ocphidry,i)
       atm2lnd_vars%forc_aer_grc(g,5)                = x2l(index_x2l_Faxa_ocphodry,i)
       atm2lnd_vars%forc_aer_grc(g,6)                = x2l(index_x2l_Faxa_ocphiwet,i)
       atm2lnd_vars%forc_aer_grc(g,7)                = x2l(index_x2l_Faxa_dstwet1,i)
       atm2lnd_vars%forc_aer_grc(g,8)                = x2l(index_x2l_Faxa_dstdry1,i)
       atm2lnd_vars%forc_aer_grc(g,9)                = x2l(index_x2l_Faxa_dstwet2,i)
       atm2lnd_vars%forc_aer_grc(g,10)               = x2l(index_x2l_Faxa_dstdry2,i)
       atm2lnd_vars%forc_aer_grc(g,11)               = x2l(index_x2l_Faxa_dstwet3,i)
       atm2lnd_vars%forc_aer_grc(g,12)               = x2l(index_x2l_Faxa_dstdry3,i)
       atm2lnd_vars%forc_aer_grc(g,13)               = x2l(index_x2l_Faxa_dstwet4,i)
       atm2lnd_vars%forc_aer_grc(g,14)               = x2l(index_x2l_Faxa_dstdry4,i)

       ! Determine optional receive fields

       if (index_x2l_Sa_co2prog /= 0) then
          co2_ppmv_prog = x2l(index_x2l_Sa_co2prog,i)   ! co2 atm state prognostic
       else
          co2_ppmv_prog = co2_ppmv
       end if

       if (index_x2l_Sa_co2diag /= 0) then
          co2_ppmv_diag = x2l(index_x2l_Sa_co2diag,i)   ! co2 atm state diagnostic
       else
          co2_ppmv_diag = co2_ppmv
       end if

       if (index_x2l_Sa_methane /= 0) then
          atm2lnd_vars%forc_pch4_grc(g) = x2l(index_x2l_Sa_methane,i)
       endif

       ! Determine derived quantities for required fields

       forc_t = atm2lnd_vars%forc_t_not_downscaled_grc(g)
       forc_q = atm2lnd_vars%forc_q_not_downscaled_grc(g)
       forc_pbot = atm2lnd_vars%forc_pbot_not_downscaled_grc(g)
       
       atm2lnd_vars%forc_hgt_u_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of wind [m]
       atm2lnd_vars%forc_hgt_t_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of temperature [m]
       atm2lnd_vars%forc_hgt_q_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of humidity [m]
       atm2lnd_vars%forc_vp_grc(g)    = forc_q * forc_pbot  / (0.622_r8 + 0.378_r8 * forc_q)
       atm2lnd_vars%forc_rho_not_downscaled_grc(g) = &
            (forc_pbot - 0.378_r8 * atm2lnd_vars%forc_vp_grc(g)) / (rair * forc_t)
       atm2lnd_vars%forc_po2_grc(g)   = o2_molar_const * forc_pbot
       atm2lnd_vars%forc_wind_grc(g)  = sqrt(atm2lnd_vars%forc_u_grc(g)**2 + atm2lnd_vars%forc_v_grc(g)**2)
       atm2lnd_vars%forc_solar_grc(g) = atm2lnd_vars%forc_solad_grc(g,1) + atm2lnd_vars%forc_solai_grc(g,1) + &
                                        atm2lnd_vars%forc_solad_grc(g,2) + atm2lnd_vars%forc_solai_grc(g,2)

       atm2lnd_vars%forc_rain_not_downscaled_grc(g)  = forc_rainc + forc_rainl
       atm2lnd_vars%forc_snow_not_downscaled_grc(g)  = forc_snowc + forc_snowl

       if (forc_t > SHR_CONST_TKFRZ) then
          e = esatw(tdc(forc_t))
       else
          e = esati(tdc(forc_t))
       end if
       qsat           = 0.622_r8*e / (forc_pbot - 0.378_r8*e)
       atm2lnd_vars%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)
       ! Make sure relative humidity is properly bounded
       ! atm2lnd_vars%forc_rh_grc(g) = min( 100.0_r8, atm2lnd_vars%forc_rh_grc(g) )
       ! atm2lnd_vars%forc_rh_grc(g) = max(   0.0_r8, atm2lnd_vars%forc_rh_grc(g) )

       ! Determine derived quantities for optional fields
       ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
       ! Note that forc_pbot is in Pa

       if (co2_type_idx == 1) then
          co2_ppmv_val = co2_ppmv_prog
       else if (co2_type_idx == 2) then
          co2_ppmv_val = co2_ppmv_diag 
       else
          co2_ppmv_val = co2_ppmv
       end if
       atm2lnd_vars%forc_pco2_grc(g)   = co2_ppmv_val * 1.e-6_r8 * forc_pbot 
       if (use_c13) then
          atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
       end if

       ! glc coupling 

       if (create_glacier_mec_landunit) then
          do num = 0,glc_nec
             glc2lnd_vars%frac_grc(g,num)  = x2l(index_x2l_Sg_frac(num),i)
             glc2lnd_vars%topo_grc(g,num)  = x2l(index_x2l_Sg_topo(num),i)
             glc2lnd_vars%hflx_grc(g,num)  = x2l(index_x2l_Flgg_hflx(num),i)
          end do
          glc2lnd_vars%icemask_grc(g)  = x2l(index_x2l_Sg_icemask,i)
          glc2lnd_vars%icemask_coupled_fluxes_grc(g)  = x2l(index_x2l_Sg_icemask_coupled_fluxes,i)
       end if

    end do

  end subroutine lnd_import

  !===============================================================================

  subroutine lnd_export( bounds, lnd2atm_vars, lnd2glc_vars, l2x)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the clm model to the coupler 
    ! 
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use clm_varctl         , only : iulog, create_glacier_mec_landunit
    use clm_time_manager   , only : get_nstep, get_step_size  
    use seq_drydep_mod     , only : n_drydep
    use shr_megan_mod      , only : shr_megan_mechcomps_n
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds  ! bounds
    type(lnd2atm_type), intent(inout) :: lnd2atm_vars ! clm land to atmosphere exchange data type
    type(lnd2glc_type), intent(inout) :: lnd2glc_vars ! clm land to atmosphere exchange data type
    real(r8)          , intent(out)   :: l2x(:,:)! land to coupler export state on land grid
    !
    ! !LOCAL VARIABLES:
    integer  :: g,i   ! indices
    integer  :: ier   ! error status
    integer  :: nstep ! time step index
    integer  :: dtime ! time step   
    integer  :: num   ! counter
    !---------------------------------------------------------------------------

    ! cesm sign convention is that fluxes are positive downward

    l2x(:,:) = 0.0_r8

    do g = bounds%begg,bounds%endg
       i = 1 + (g-bounds%begg)
       l2x(index_l2x_Sl_t,i)        =  lnd2atm_vars%t_rad_grc(g)
       l2x(index_l2x_Sl_snowh,i)    =  lnd2atm_vars%h2osno_grc(g)
       l2x(index_l2x_Sl_avsdr,i)    =  lnd2atm_vars%albd_grc(g,1)
       l2x(index_l2x_Sl_anidr,i)    =  lnd2atm_vars%albd_grc(g,2)
       l2x(index_l2x_Sl_avsdf,i)    =  lnd2atm_vars%albi_grc(g,1)
       l2x(index_l2x_Sl_anidf,i)    =  lnd2atm_vars%albi_grc(g,2)
       l2x(index_l2x_Sl_tref,i)     =  lnd2atm_vars%t_ref2m_grc(g)
       l2x(index_l2x_Sl_qref,i)     =  lnd2atm_vars%q_ref2m_grc(g)
       l2x(index_l2x_Sl_u10,i)      =  lnd2atm_vars%u_ref10m_grc(g)
       l2x(index_l2x_Fall_taux,i)   = -lnd2atm_vars%taux_grc(g)
       l2x(index_l2x_Fall_tauy,i)   = -lnd2atm_vars%tauy_grc(g)
       l2x(index_l2x_Fall_lat,i)    = -lnd2atm_vars%eflx_lh_tot_grc(g)
       l2x(index_l2x_Fall_sen,i)    = -lnd2atm_vars%eflx_sh_tot_grc(g)
       l2x(index_l2x_Fall_lwup,i)   = -lnd2atm_vars%eflx_lwrad_out_grc(g)
       l2x(index_l2x_Fall_evap,i)   = -lnd2atm_vars%qflx_evap_tot_grc(g)
       l2x(index_l2x_Fall_swnet,i)  =  lnd2atm_vars%fsa_grc(g)
       if (index_l2x_Fall_fco2_lnd /= 0) then
          l2x(index_l2x_Fall_fco2_lnd,i) = -lnd2atm_vars%nee_grc(g)  
       end if

       ! Additional fields for DUST, PROGSSLT, dry-deposition and VOC
       ! These are now standard fields, but the check on the index makes sure the driver handles them
       if (index_l2x_Sl_ram1      /= 0 )  l2x(index_l2x_Sl_ram1,i)     =  lnd2atm_vars%ram1_grc(g)
       if (index_l2x_Sl_fv        /= 0 )  l2x(index_l2x_Sl_fv,i)       =  lnd2atm_vars%fv_grc(g)
       if (index_l2x_Sl_soilw     /= 0 )  l2x(index_l2x_Sl_soilw,i)    =  lnd2atm_vars%h2osoi_vol_grc(g,1)
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x(index_l2x_Fall_flxdst1,i)= -lnd2atm_vars%flxdst_grc(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x(index_l2x_Fall_flxdst2,i)= -lnd2atm_vars%flxdst_grc(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x(index_l2x_Fall_flxdst3,i)= -lnd2atm_vars%flxdst_grc(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x(index_l2x_Fall_flxdst4,i)= -lnd2atm_vars%flxdst_grc(g,4)


       ! for dry dep velocities
       if (index_l2x_Sl_ddvel     /= 0 )  then
          l2x(index_l2x_Sl_ddvel:index_l2x_Sl_ddvel+n_drydep-1,i) = &
               lnd2atm_vars%ddvel_grc(g,:n_drydep)
       end if

       ! for MEGAN VOC emis fluxes
       if (index_l2x_Fall_flxvoc  /= 0 ) then
          l2x(index_l2x_Fall_flxvoc:index_l2x_Fall_flxvoc+shr_megan_mechcomps_n-1,i) = &
               -lnd2atm_vars%flxvoc_grc(g,:shr_megan_mechcomps_n)
       end if

       if (index_l2x_Fall_methane /= 0) then
          l2x(index_l2x_Fall_methane,i) = -lnd2atm_vars%flux_ch4_grc(g) 
       endif

       ! sign convention is positive downward with 
       ! hierarchy of atm/glc/lnd/rof/ice/ocn.  so water sent from land to rof is positive

       l2x(index_l2x_Flrl_rofl,i) = lnd2atm_vars%qflx_rofliq_grc(g)
       l2x(index_l2x_Flrl_rofi,i) = lnd2atm_vars%qflx_rofice_grc(g)

       ! glc coupling

       if (create_glacier_mec_landunit) then
          do num = 0,glc_nec
             l2x(index_l2x_Sl_tsrf(num),i)   = lnd2glc_vars%tsrf_grc(g,num)
             l2x(index_l2x_Sl_topo(num),i)   = lnd2glc_vars%topo_grc(g,num)
             l2x(index_l2x_Flgl_qice(num),i) = lnd2glc_vars%qice_grc(g,num)
          end do
       end if

    end do

  end subroutine lnd_export

end module lnd_import_export
