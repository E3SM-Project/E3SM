module clm_atmlnd

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_atmlnd
!
! !DESCRIPTION:
! Handle atm2lnd, lnd2atm mapping
!
! !USES:
  use clm_varpar  , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon  , only : rair, grav, cpair, hfus, tfrz
  use clm_varctl  , only : iulog, use_cn
  use decompMod   , only : get_proc_bounds
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun
  use seq_drydep_mod, only : n_drydep, drydep_method, DD_XLND
  use shr_megan_mod,  only : shr_megan_mechcomps_n
!
! !PUBLIC TYPES:
  implicit none
!----------------------------------------------------
! atmosphere -> land variables structure
!----------------------------------------------------
  type atm2lnd_type
     real(r8), pointer :: forc_t(:)       !atmospheric temperature (Kelvin)
     real(r8), pointer :: forc_u(:)       !atm wind speed, east direction (m/s)
     real(r8), pointer :: forc_v(:)       !atm wind speed, north direction (m/s)
     real(r8), pointer :: forc_wind(:)    !atmospheric wind speed   
     real(r8), pointer :: forc_q(:)       !atmospheric specific humidity (kg/kg)
     real(r8), pointer :: forc_hgt(:)     !atmospheric reference height (m)
     real(r8), pointer :: forc_hgt_u(:)   !obs height of wind [m] (new)
     real(r8), pointer :: forc_hgt_t(:)   !obs height of temperature [m] (new)
     real(r8), pointer :: forc_hgt_q(:)   !obs height of humidity [m] (new)
     real(r8), pointer :: forc_pbot(:)    !atmospheric pressure (Pa)
     real(r8), pointer :: forc_th(:)      !atm potential temperature (Kelvin)
     real(r8), pointer :: forc_vp(:)      !atmospheric vapor pressure (Pa) 
     real(r8), pointer :: forc_rho(:)     !density (kg/m**3)
     real(r8), pointer :: forc_rh(:)      !atmospheric relative humidity (%)
     real(r8), pointer :: forc_psrf(:)    !surface pressure (Pa)
     real(r8), pointer :: forc_pco2(:)    !CO2 partial pressure (Pa)
     real(r8), pointer :: forc_lwrad(:)   !downwrd IR longwave radiation (W/m**2)
     real(r8), pointer :: forc_solad(:,:) !direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll )
     real(r8), pointer :: forc_solai(:,:) !diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)
     real(r8), pointer :: forc_solar(:)   !incident solar radiation
     real(r8), pointer :: forc_rain(:)    !rain rate [mm/s]
     real(r8), pointer :: forc_snow(:)    !snow rate [mm/s]
     real(r8), pointer :: forc_ndep(:)    !nitrogen deposition rate (gN/m2/s)
     real(r8), pointer :: rainf(:)        !ALMA rain+snow [mm/s]
     real(r8), pointer :: forc_pc13o2(:)  !C13O2 partial pressure (Pa)
     real(r8), pointer :: forc_po2(:)     !O2 partial pressure (Pa)
     real(r8), pointer :: forc_flood(:)   ! rof flood (mm/s)
     real(r8), pointer :: volr(:)   ! rof volr (m3)
     real(r8), pointer :: forc_aer(:,:)   ! aerosol deposition array
  end type atm2lnd_type

!----------------------------------------------------
! land -> atmosphere variables structure
!----------------------------------------------------
  type lnd2atm_type
     real(r8), pointer :: t_rad(:)        !radiative temperature (Kelvin)
     real(r8), pointer :: t_ref2m(:)      !2m surface air temperature (Kelvin)
     real(r8), pointer :: q_ref2m(:)      !2m surface specific humidity (kg/kg)
     real(r8), pointer :: u_ref10m(:)     !10m surface wind speed (m/sec)
     real(r8), pointer :: h2osno(:)       !snow water (mm H2O)
     real(r8), pointer :: albd(:,:)       !(numrad) surface albedo (direct)
     real(r8), pointer :: albi(:,:)       !(numrad) surface albedo (diffuse)
     real(r8), pointer :: taux(:)         !wind stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy(:)         !wind stress: n-s (kg/m/s**2)
     real(r8), pointer :: eflx_lh_tot(:)  !total latent HF (W/m**2)  [+ to atm]
     real(r8), pointer :: eflx_sh_tot(:)  !total sensible HF (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lwrad_out(:) !IR (longwave) radiation (W/m**2)
     real(r8), pointer :: qflx_evap_tot(:)!qflx_evap_soi + qflx_evap_can + qflx_tran_veg
     real(r8), pointer :: fsa(:)          !solar rad absorbed (total) (W/m**2)
     real(r8), pointer :: nee(:)          !net CO2 flux (kg CO2/m**2/s) [+ to atm]
     real(r8), pointer :: ram1(:)         !aerodynamical resistance (s/m)
     real(r8), pointer :: fv(:)           !friction velocity (m/s) (for dust model)
     real(r8), pointer :: h2osoi_vol(:,:) !volumetric soil water (0~watsat, m3/m3, nlevgrnd) (for dust model)
     real(r8), pointer :: rofliq(:)       ! rof liq forcing
     real(r8), pointer :: rofice(:)       ! rof ice forcing
     real(r8), pointer :: flxdst(:,:)     !dust flux (size bins)
     real(r8), pointer :: ddvel(:,:)      !dry deposition velocities
     real(r8), pointer :: flxvoc(:,:)     ! VOC flux (size bins)
  end type lnd2atm_type
  
  type(atm2lnd_type),public,target :: clm_a2l      ! a2l fields on clm grid
  type(lnd2atm_type),public,target :: clm_l2a      ! l2a fields on clm grid

! !PUBLIC MEMBER FUNCTIONS:
  public :: init_atm2lnd_type
  public :: init_lnd2atm_type
  public :: clm_map2gcell
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein and Tony Craig, 2006-01-10
!
! !PRIVATE MEMBER FUNCTIONS:

!EOP
!----------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_atm2lnd_type
!
! !INTERFACE:
  subroutine init_atm2lnd_type(beg, end, a2l)
!
! !DESCRIPTION:
! Initialize atmospheric variables required by the land
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (atm2lnd_type), intent(inout):: a2l
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Modified by T Craig, 11/01/05 for finemesh project
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  allocate(a2l%forc_t(beg:end))
  allocate(a2l%forc_u(beg:end))
  allocate(a2l%forc_v(beg:end))
  allocate(a2l%forc_wind(beg:end))
  allocate(a2l%forc_q(beg:end))
  allocate(a2l%forc_rh(beg:end))
  allocate(a2l%forc_hgt(beg:end))
  allocate(a2l%forc_hgt_u(beg:end))
  allocate(a2l%forc_hgt_t(beg:end))
  allocate(a2l%forc_hgt_q(beg:end))
  allocate(a2l%forc_pbot(beg:end))
  allocate(a2l%forc_th(beg:end))
  allocate(a2l%forc_vp(beg:end))
  allocate(a2l%forc_rho(beg:end))
  allocate(a2l%forc_psrf(beg:end))
  allocate(a2l%forc_pco2(beg:end))
  allocate(a2l%forc_lwrad(beg:end))
  allocate(a2l%forc_solad(beg:end,numrad))
  allocate(a2l%forc_solai(beg:end,numrad))
  allocate(a2l%forc_solar(beg:end))
  allocate(a2l%forc_rain(beg:end))
  allocate(a2l%forc_snow(beg:end))
  allocate(a2l%forc_ndep(beg:end))
  allocate(a2l%rainf(beg:end))
  allocate(a2l%forc_pc13o2(beg:end))
  allocate(a2l%forc_po2(beg:end))
  allocate(a2l%forc_flood(beg:end))
  allocate(a2l%volr(beg:end))
  allocate(a2l%forc_aer(beg:end,14))

  ! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  a2l%forc_t(beg:end) = ival
  a2l%forc_u(beg:end) = ival
  a2l%forc_v(beg:end) = ival
  a2l%forc_wind(beg:end) = ival
  a2l%forc_q(beg:end) = ival
  a2l%forc_rh(beg:end) = ival
  a2l%forc_hgt(beg:end) = ival
  a2l%forc_hgt_u(beg:end) = ival
  a2l%forc_hgt_t(beg:end) = ival
  a2l%forc_hgt_q(beg:end) = ival
  a2l%forc_pbot(beg:end) = ival
  a2l%forc_th(beg:end) = ival
  a2l%forc_vp(beg:end) = ival
  a2l%forc_rho(beg:end) = ival
  a2l%forc_psrf(beg:end) = ival
  a2l%forc_pco2(beg:end) = ival
  a2l%forc_lwrad(beg:end) = ival
  a2l%forc_solad(beg:end,1:numrad) = ival
  a2l%forc_solai(beg:end,1:numrad) = ival
  a2l%forc_solar(beg:end) = ival
  a2l%forc_rain(beg:end) = ival
  a2l%forc_snow(beg:end) = ival
  a2l%forc_ndep(beg:end) = ival
  a2l%rainf(beg:end) = nan
  a2l%forc_pc13o2(beg:end) = ival
  a2l%forc_po2(beg:end) = ival
  a2l%forc_flood(beg:end) = ival
  a2l%volr(beg:end) = ival
  a2l%forc_aer(beg:end,:) = ival

end subroutine init_atm2lnd_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lnd2atm_type
!
! !INTERFACE:
  subroutine init_lnd2atm_type(beg, end, l2a)
!
! !DESCRIPTION:
! Initialize land variables required by the atmosphere
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (lnd2atm_type), intent(inout):: l2a
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Modified by T Craig, 11/01/05 for finemesh project
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  allocate(l2a%t_rad(beg:end))
  allocate(l2a%t_ref2m(beg:end))
  allocate(l2a%q_ref2m(beg:end))
  allocate(l2a%u_ref10m(beg:end))
  allocate(l2a%h2osno(beg:end))
  allocate(l2a%albd(beg:end,1:numrad))
  allocate(l2a%albi(beg:end,1:numrad))
  allocate(l2a%taux(beg:end))
  allocate(l2a%tauy(beg:end))
  allocate(l2a%eflx_lwrad_out(beg:end))
  allocate(l2a%eflx_sh_tot(beg:end))
  allocate(l2a%eflx_lh_tot(beg:end))
  allocate(l2a%qflx_evap_tot(beg:end))
  allocate(l2a%fsa(beg:end))
  allocate(l2a%nee(beg:end))
  allocate(l2a%ram1(beg:end))
  allocate(l2a%fv(beg:end))
  allocate(l2a%h2osoi_vol(beg:end,1:nlevgrnd))
  allocate(l2a%rofliq(beg:end))
  allocate(l2a%rofice(beg:end))
  allocate(l2a%flxdst(beg:end,1:ndst))
  if (shr_megan_mechcomps_n>0) then
     allocate(l2a%flxvoc(beg:end,1:shr_megan_mechcomps_n))
  endif
  if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
     allocate(l2a%ddvel(beg:end,1:n_drydep))
  end if

  ! ival = nan   ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  l2a%t_rad(beg:end) = ival
  l2a%t_ref2m(beg:end) = ival
  l2a%q_ref2m(beg:end) = ival
  l2a%u_ref10m(beg:end) = ival
  l2a%h2osno(beg:end) = ival
  l2a%albd(beg:end,1:numrad) = ival
  l2a%albi(beg:end,1:numrad) = ival
  l2a%taux(beg:end) = ival
  l2a%tauy(beg:end) = ival
  l2a%eflx_lwrad_out(beg:end) = ival
  l2a%eflx_sh_tot(beg:end) = ival
  l2a%eflx_lh_tot(beg:end) = ival
  l2a%qflx_evap_tot(beg:end) = ival
  l2a%fsa(beg:end) = ival
  l2a%nee(beg:end) = ival
  l2a%ram1(beg:end) = ival
  l2a%fv(beg:end) = ival
  l2a%h2osoi_vol(beg:end,1:nlevgrnd) = ival
  l2a%rofliq(beg:end) = ival
  l2a%rofice(beg:end) = ival
  l2a%flxdst(beg:end,1:ndst) = ival
  if (shr_megan_mechcomps_n>0) then
     l2a%flxvoc(beg:end,1:shr_megan_mechcomps_n) = ival
  endif
  if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
     l2a%ddvel(beg:end, : ) = ival
  end if

end subroutine init_lnd2atm_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_map2gcell
!
! !INTERFACE: subroutine clm_map2gcell(init)
subroutine clm_map2gcell(init)
!
! !DESCRIPTION:
! Compute l2a component of gridcell derived type
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
  use subgridAveMod
  use clm_varcon  , only : sb
  use clm_varpar  , only : numrad
!
! !ARGUMENTS:
  implicit none
  save
  logical, optional, intent(in) :: init  ! if true=>only set a subset of arguments
!
! !REVISION HISTORY:
! Mariana Vertenstein: created 03/10-25
!
!
! !LOCAL VARIABLES:
  integer :: begp, endp      ! per-proc beginning and ending pft indices
  integer :: begc, endc      ! per-proc beginning and ending column indices
  integer :: begl, endl      ! per-proc beginning and ending landunit indices
  integer :: begg, endg      ! per-proc gridcell ending gridcell indices
!
! !USES:
!
! !REVISION HISTORY:
! 03-04-27 : Created by Mariana Vertenstein
! 03-08-25 : Updated to vector data structure (Mariana Vertenstein)
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: g                           ! indices
  type(gridcell_type), pointer :: gptr   ! pointer to gridcell derived subtype
  type(landunit_type), pointer :: lptr   ! pointer to landunit derived subtype
  type(column_type)  , pointer :: cptr   ! pointer to column derived subtype
  type(pft_type)     , pointer :: pptr   ! pointer to pft derived subtype
  integer             :: n               ! Loop index over nmap
  real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
  real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
  real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
  ! The following converts g of C to kg of CO2
  real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)

!------------------------------------------------------------------------

  ! Set pointers into derived type

  gptr => grc
  lptr => lun
  cptr => col
  pptr => pft

  ! Determine processor bounds

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

  ! Compute gridcell averages. 

  if (present(init)) then

     call c2g(begc, endc, begl, endl, begg, endg, &
          cws%h2osno, clm_l2a%h2osno, &
          c2l_scale_type= 'urbanf', l2g_scale_type='unity')
     do g = begg,endg
        clm_l2a%h2osno(g) = clm_l2a%h2osno(g)/1000._r8
     end do
     
      call c2g(begc, endc, begl, endl, begg, endg, nlevgrnd, &
          cws%h2osoi_vol, clm_l2a%h2osoi_vol, &
          c2l_scale_type= 'urbanf', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pps%albd, clm_l2a%albd,&
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pps%albi, clm_l2a%albi,&
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pef%eflx_lwrad_out, clm_l2a%eflx_lwrad_out,&
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')
     do g = begg,endg
        clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
     end do

  else

     call c2g(begc, endc, begl, endl, begg, endg, cws%h2osno, clm_l2a%h2osno,&
          c2l_scale_type= 'urbanf', l2g_scale_type='unity')
     do g = begg,endg
        clm_l2a%h2osno(g) = clm_l2a%h2osno(g)/1000._r8
     end do

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pps%albd, clm_l2a%albd, &
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pps%albi, clm_l2a%albi, &
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pes%t_ref2m, clm_l2a%t_ref2m, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pes%q_ref2m, clm_l2a%q_ref2m, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pps%u10_clm, clm_l2a%u_ref10m, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pmf%taux, clm_l2a%taux, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pmf%tauy, clm_l2a%tauy, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pef%eflx_lh_tot, clm_l2a%eflx_lh_tot, &
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

!DML note: use new array: gef%eflx_sh_totg
!     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
!          pef%eflx_sh_tot, clm_l2a%eflx_sh_tot, &
!          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

     do g = begg,endg
        clm_l2a%eflx_sh_tot(g) = gef%eflx_sh_totg(g)
     end do

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pwf%qflx_evap_tot, clm_l2a%qflx_evap_tot, & 
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pef%fsa, clm_l2a%fsa, &
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')
                  
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pef%eflx_lwrad_out, clm_l2a%eflx_lwrad_out, &
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')
                  
     if (use_cn) then
        call c2g(begc, endc, begl, endl, begg, endg, &
             ccf%nee, clm_l2a%nee, &
             c2l_scale_type= 'unity', l2g_scale_type='unity')
     else
        call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
             pcf%fco2, clm_l2a%nee, &
             p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
        ! Note that fco2 in is umolC/m2/sec so units need to be changed to gC/m2/sec
        do g = begg,endg
           clm_l2a%nee(g) = clm_l2a%nee(g)*12.011e-6_r8
        end do
     end if

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pps%fv, clm_l2a%fv, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pps%ram1, clm_l2a%ram1, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     do g = begg,endg
        clm_l2a%rofliq(g) = gwf%qflx_runoffg(g)
        clm_l2a%rofice(g) = gwf%qflx_snwcp_iceg(g)
     end do

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, ndst, &
          pdf%flx_mss_vrt_dst, clm_l2a%flxdst, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     if (shr_megan_mechcomps_n>0) then
        call p2g(begp, endp, begc, endc, begl, endl, begg, endg, shr_megan_mechcomps_n, &
             pvf%vocflx, clm_l2a%flxvoc, &
             p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
     endif

     if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
        call p2g(begp, endp, begc, endc, begl, endl, begg, endg, n_drydep, &
             pdd%drydepvel, clm_l2a%ddvel, &
             p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
     endif

     ! Convert from gC/m2/s to kgCO2/m2/s
     do g = begg,endg
        clm_l2a%nee(g) = clm_l2a%nee(g)*convertgC2kgCO2
     end do

     do g = begg,endg
        clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
     end do

  end if

end subroutine clm_map2gcell

!------------------------------------------------------------------------
!------------------------------------------------------------------------
end module clm_atmlnd

