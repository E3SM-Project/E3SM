module lnd_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use abortutils   , only: endrun
  use clm_atmlnd   , only: lnd2atm_type
  use clm_glclnd   , only: lnd2glc_type
  use decompmod    , only: bounds_type
  use clm_cpl_indices
  use clmtype
  implicit none

contains

  !===============================================================================

  subroutine lnd_import( bounds, x2l, a2l, a2l_not_downscaled_gcell, x2s)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the coupler to the land model 
    !
    ! !USES:
    use clm_atmlnd      , only: atm2lnd_type, atm2lnd_downscaled_fields_type
    use clm_glclnd      , only: glc2lnd_type
    use clm_varctl      , only: co2_type, co2_ppmv, iulog, use_c13, create_glacier_mec_landunit
    use clm_varcon      , only: rair, o2_molar_const, c13ratio
    use shr_const_mod   , only: SHR_CONST_TKFRZ
    use domainMod       , only: ldomain
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds                 ! bounds
    real(r8)                             , intent(in)    :: x2l(:,:)               ! driver import state to land model
    type(atm2lnd_type)                   , intent(inout) :: a2l                    ! clm internal input data type
    type(atm2lnd_downscaled_fields_type) , intent(inout) :: a2l_not_downscaled_gcell ! clm internal input data type
    type(glc2lnd_type)                   , intent(inout) :: x2s                    ! clm internal input data type
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

       a2l%forc_flood(g)   = -x2l(index_x2l_Flrr_flood,i)  

       a2l%volr(g)   = x2l(index_x2l_Flrr_volr,i) * (ldomain%area(g) * 1.e6_r8)

       ! Determine required receive fields

       a2l%forc_hgt(g)     = x2l(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
       a2l%forc_u(g)       = x2l(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
       a2l%forc_v(g)       = x2l(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
       a2l%forc_solad(g,2) = x2l(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
       a2l%forc_solad(g,1) = x2l(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
       a2l%forc_solai(g,2) = x2l(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
       a2l%forc_solai(g,1) = x2l(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

       a2l_not_downscaled_gcell%forc_th(g)      = x2l(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
       a2l_not_downscaled_gcell%forc_q(g)       = x2l(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
       a2l_not_downscaled_gcell%forc_pbot(g)    = x2l(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
       a2l_not_downscaled_gcell%forc_t(g)       = x2l(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
       a2l_not_downscaled_gcell%forc_lwrad(g)   = x2l(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2

       forc_rainc          = x2l(index_x2l_Faxa_rainc,i)   ! mm/s
       forc_rainl          = x2l(index_x2l_Faxa_rainl,i)   ! mm/s
       forc_snowc          = x2l(index_x2l_Faxa_snowc,i)   ! mm/s
       forc_snowl          = x2l(index_x2l_Faxa_snowl,i)   ! mm/s

       !
       ! anomaly forcing code 
       !
       ! bias correct atmospheric input fields if streams exist
       if (index_x2l_Sa_precsf /= 0) then
          a2l%bc_precip(g)    = x2l(index_x2l_Sa_precsf,i)
          a2l%bc_precip(g)    = min(1.e2_r8,a2l%bc_precip(g))
          forc_rainc          = forc_rainc * a2l%bc_precip(g)
          forc_rainl          = forc_rainl * a2l%bc_precip(g)
          forc_snowc          = forc_snowc * a2l%bc_precip(g)
          forc_snowl          = forc_snowl * a2l%bc_precip(g)
       endif

       ! adjust atmospheric input fields if anomaly forcing streams exist
       if (index_x2l_Sa_u_af /= 0) then
          a2l%af_uwind(g)    = x2l(index_x2l_Sa_u_af,i)
          a2l%forc_u(g) = a2l%forc_u(g) + a2l%af_uwind(g)
       endif

       if (index_x2l_Sa_v_af /= 0) then
          a2l%af_vwind(g)    = x2l(index_x2l_Sa_v_af,i)
          a2l%forc_v(g) = a2l%forc_v(g) + a2l%af_vwind(g)
       endif

       if (index_x2l_Sa_shum_af /= 0) then
          a2l%af_shum(g)    = x2l(index_x2l_Sa_shum_af,i)
          a2l_not_downscaled_gcell%forc_q(g) = a2l_not_downscaled_gcell%forc_q(g) + a2l%af_shum(g)
          ! avoid possible negative q values
          if(a2l_not_downscaled_gcell%forc_q(g) < 0._r8) then 
             a2l_not_downscaled_gcell%forc_q(g) = 1.e-6_r8
          endif
       endif

       if (index_x2l_Sa_pbot_af /= 0) then
          a2l%af_pbot(g)    = x2l(index_x2l_Sa_pbot_af,i)
          a2l_not_downscaled_gcell%forc_pbot(g) =  a2l_not_downscaled_gcell%forc_pbot(g) + a2l%af_pbot(g)
       endif

       if (index_x2l_Sa_tbot_af /= 0) then
          a2l%af_tbot(g)    = x2l(index_x2l_Sa_tbot_af,i)
          a2l_not_downscaled_gcell%forc_t(g) =  a2l_not_downscaled_gcell%forc_t(g) + a2l%af_tbot(g)
       endif

       if (index_x2l_Sa_lwdn_af /= 0) then
          a2l%af_lwdn(g)    = x2l(index_x2l_Sa_lwdn_af,i)
          a2l_not_downscaled_gcell%forc_lwrad(g) = a2l_not_downscaled_gcell%forc_lwrad(g) * a2l%af_lwdn(g)
       endif

       if (index_x2l_Sa_prec_af /= 0) then
          a2l%af_precip(g)    = x2l(index_x2l_Sa_prec_af,i)
          forc_rainc = forc_rainc * a2l%af_precip(g)
          forc_rainl = forc_rainl * a2l%af_precip(g)
          forc_snowc = forc_snowc * a2l%af_precip(g)
          forc_snowl = forc_snowl * a2l%af_precip(g)
       endif

       if (index_x2l_Sa_swdn_af /= 0) then
          a2l%af_swdn(g)    = x2l(index_x2l_Sa_swdn_af,i)
          a2l%forc_solad(g,2) = a2l%forc_solad(g,2) * a2l%af_swdn(g)
          a2l%forc_solad(g,1) = a2l%forc_solad(g,1) * a2l%af_swdn(g)
          a2l%forc_solai(g,2) = a2l%forc_solai(g,2) * a2l%af_swdn(g)
          a2l%forc_solai(g,1) = a2l%forc_solai(g,1) * a2l%af_swdn(g)
       endif
       !
       ! anomaly forcing code ends
       !

       ! atmosphere coupling, for prognostic/prescribed aerosols
       a2l%forc_aer(g,1)  =  x2l(index_x2l_Faxa_bcphidry,i)
       a2l%forc_aer(g,2)  =  x2l(index_x2l_Faxa_bcphodry,i)
       a2l%forc_aer(g,3)  =  x2l(index_x2l_Faxa_bcphiwet,i)
       a2l%forc_aer(g,4)  =  x2l(index_x2l_Faxa_ocphidry,i)
       a2l%forc_aer(g,5)  =  x2l(index_x2l_Faxa_ocphodry,i)
       a2l%forc_aer(g,6)  =  x2l(index_x2l_Faxa_ocphiwet,i)
       a2l%forc_aer(g,7)  =  x2l(index_x2l_Faxa_dstwet1,i)
       a2l%forc_aer(g,8)  =  x2l(index_x2l_Faxa_dstdry1,i)
       a2l%forc_aer(g,9)  =  x2l(index_x2l_Faxa_dstwet2,i)
       a2l%forc_aer(g,10) =  x2l(index_x2l_Faxa_dstdry2,i)
       a2l%forc_aer(g,11) =  x2l(index_x2l_Faxa_dstwet3,i)
       a2l%forc_aer(g,12) =  x2l(index_x2l_Faxa_dstdry3,i)
       a2l%forc_aer(g,13) =  x2l(index_x2l_Faxa_dstwet4,i)
       a2l%forc_aer(g,14) =  x2l(index_x2l_Faxa_dstdry4,i)

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
          a2l%forc_pch4(g) = x2l(index_x2l_Sa_methane,i)
       endif

       ! Determine derived quantities for required fields

       forc_t = a2l_not_downscaled_gcell%forc_t(g)
       forc_q = a2l_not_downscaled_gcell%forc_q(g)
       forc_pbot = a2l_not_downscaled_gcell%forc_pbot(g)
       
       a2l%forc_hgt_u(g) = a2l%forc_hgt(g)    !observational height of wind [m]
       a2l%forc_hgt_t(g) = a2l%forc_hgt(g)    !observational height of temperature [m]
       a2l%forc_hgt_q(g) = a2l%forc_hgt(g)    !observational height of humidity [m]
       a2l%forc_vp(g)    = forc_q * forc_pbot &
            / (0.622_r8 + 0.378_r8 * forc_q)
       a2l_not_downscaled_gcell%forc_rho(g)   = (forc_pbot - 0.378_r8 * a2l%forc_vp(g)) &
            / (rair * forc_t)
       a2l%forc_po2(g)   = o2_molar_const * forc_pbot
       a2l%forc_wind(g)  = sqrt(a2l%forc_u(g)**2 + a2l%forc_v(g)**2)
       a2l%forc_solar(g) = a2l%forc_solad(g,1) + a2l%forc_solai(g,1) + &
            a2l%forc_solad(g,2) + a2l%forc_solai(g,2)

       a2l_not_downscaled_gcell%forc_rain(g)  = forc_rainc + forc_rainl
       a2l_not_downscaled_gcell%forc_snow(g)  = forc_snowc + forc_snowl

       if (forc_t > SHR_CONST_TKFRZ) then
          e = esatw(tdc(forc_t))
       else
          e = esati(tdc(forc_t))
       end if
       qsat           = 0.622_r8*e / (forc_pbot - 0.378_r8*e)
       a2l%forc_rh(g) = 100.0_r8*(forc_q / qsat)
       ! Make sure relative humidity is properly bounded
       ! a2l%forc_rh(g) = min( 100.0_r8, a2l%forc_rh(g) )
       ! a2l%forc_rh(g) = max(   0.0_r8, a2l%forc_rh(g) )

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
       a2l%forc_pco2(g)   = co2_ppmv_val * 1.e-6_r8 * forc_pbot 
       if (use_c13) then
          a2l%forc_pc13o2(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
       end if

       ! glc coupling 

       if (create_glacier_mec_landunit) then
          do num = 0,glc_nec
             x2s%frac(g,num)  = x2l(index_x2l_Sg_frac(num),i)
             x2s%topo(g,num)  = x2l(index_x2l_Sg_topo(num),i)
             x2s%hflx(g,num)  = x2l(index_x2l_Flgg_hflx(num),i)
          end do
          x2s%icemask(g)  = x2l(index_x2l_Sg_icemask,i)
       end if

    end do

  end subroutine lnd_import

  !===============================================================================

  subroutine lnd_export( bounds, clm_l2a, clm_s2x, l2x)

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
    type(lnd2atm_type), intent(inout) :: clm_l2a ! clm land to atmosphere exchange data type
    type(lnd2glc_type), intent(inout) :: clm_s2x ! clm land to atmosphere exchange data type
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
       l2x(index_l2x_Sl_t,i)        =  clm_l2a%t_rad(g)
       l2x(index_l2x_Sl_snowh,i)    =  clm_l2a%h2osno(g)
       l2x(index_l2x_Sl_avsdr,i)    =  clm_l2a%albd(g,1)
       l2x(index_l2x_Sl_anidr,i)    =  clm_l2a%albd(g,2)
       l2x(index_l2x_Sl_avsdf,i)    =  clm_l2a%albi(g,1)
       l2x(index_l2x_Sl_anidf,i)    =  clm_l2a%albi(g,2)
       l2x(index_l2x_Sl_tref,i)     =  clm_l2a%t_ref2m(g)
       l2x(index_l2x_Sl_qref,i)     =  clm_l2a%q_ref2m(g)
       l2x(index_l2x_Sl_u10,i)      =  clm_l2a%u_ref10m(g)
       l2x(index_l2x_Fall_taux,i)   = -clm_l2a%taux(g)
       l2x(index_l2x_Fall_tauy,i)   = -clm_l2a%tauy(g)
       l2x(index_l2x_Fall_lat,i)    = -clm_l2a%eflx_lh_tot(g)
       l2x(index_l2x_Fall_sen,i)    = -clm_l2a%eflx_sh_tot(g)
       l2x(index_l2x_Fall_lwup,i)   = -clm_l2a%eflx_lwrad_out(g)
       l2x(index_l2x_Fall_evap,i)   = -clm_l2a%qflx_evap_tot(g)
       l2x(index_l2x_Fall_swnet,i)  =  clm_l2a%fsa(g)
       if (index_l2x_Fall_fco2_lnd /= 0) then
          l2x(index_l2x_Fall_fco2_lnd,i) = -clm_l2a%nee(g)  
       end if

       ! Additional fields for DUST, PROGSSLT, dry-deposition and VOC
       ! These are now standard fields, but the check on the index makes sure the driver handles them
       if (index_l2x_Sl_ram1      /= 0 )  l2x(index_l2x_Sl_ram1,i) = clm_l2a%ram1(g)
       if (index_l2x_Sl_fv        /= 0 )  l2x(index_l2x_Sl_fv,i)   = clm_l2a%fv(g)
       if (index_l2x_Sl_soilw     /= 0 )  l2x(index_l2x_Sl_soilw,i)   = clm_l2a%h2osoi_vol(g,1)
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x(index_l2x_Fall_flxdst1,i)= -clm_l2a%flxdst(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x(index_l2x_Fall_flxdst2,i)= -clm_l2a%flxdst(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x(index_l2x_Fall_flxdst3,i)= -clm_l2a%flxdst(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x(index_l2x_Fall_flxdst4,i)= -clm_l2a%flxdst(g,4)


       ! for dry dep velocities
       if (index_l2x_Sl_ddvel     /= 0 )  then
          l2x(index_l2x_Sl_ddvel:index_l2x_Sl_ddvel+n_drydep-1,i) = &
               clm_l2a%ddvel(g,:n_drydep)
       end if

       ! for MEGAN VOC emis fluxes
       if (index_l2x_Fall_flxvoc  /= 0 ) then
          l2x(index_l2x_Fall_flxvoc:index_l2x_Fall_flxvoc+shr_megan_mechcomps_n-1,i) = &
               -clm_l2a%flxvoc(g,:shr_megan_mechcomps_n)
       end if

       if (index_l2x_Fall_methane /= 0) then
          l2x(index_l2x_Fall_methane,i) = -clm_l2a%flux_ch4(g) 
       endif

       ! sign convention is positive downward with 
       ! hierarchy of atm/glc/lnd/rof/ice/ocn.  so water sent from land to rof is positive

       l2x(index_l2x_Flrl_rofl,i) = clm_l2a%rofliq(g)
       l2x(index_l2x_Flrl_rofi,i) = clm_l2a%rofice(g)

       ! glc coupling

       if (create_glacier_mec_landunit) then
          do num = 0,glc_nec
             l2x(index_l2x_Sl_tsrf(num),i)   = clm_s2x%tsrf(g,num)
             l2x(index_l2x_Sl_topo(num),i)   = clm_s2x%topo(g,num)
             l2x(index_l2x_Flgl_qice(num),i) = clm_s2x%qice(g,num)
          end do
       end if

    end do

  end subroutine lnd_export

end module lnd_import_export
