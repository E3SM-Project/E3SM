#ifndef CAM
#include "config.h"

module p3phys

use control_mod,          only: theta_hydrostatic_mode, dt_remap_factor
use dimensions_mod,       only: np, nlev, nlevp , qsize, qsize_d, nelemd
use element_mod,          only: element_t
use eos
use element_ops,          only: get_pottemp, get_temperature, get_hydro_pressure
use element_state,        only: nt=>timelevels
use hybrid_mod,           only: hybrid_t
use kinds,                only: rl=>real_kind, iulog
use parallel_mod,         only: abortmp,iam
use reduction_mod,        only: parallelmax, parallelmin

use physical_constants,   only: latvap, latice, &
                                gravit=>g, p0, &
                                rdry=>rgas, cpdry=>cp, cpv=>cpwater_vapor, cl, cvdry, cvv, &
                                rvapor=>rwater_vapor, rho_liquidH20, pi => dd_pi
use hybvcoord_mod,        only: hvcoord_t
use time_mod,             only: time_at, TimeLevel_t
use micro_p3,             only: p3_main, p3_init
use wv_sat_scream,        only: qv_sat

use time_mod,             only: tstep

use micro_p3_utils,       only: micro_p3_utils_init

implicit none

public :: interface_to_p3

contains


subroutine interface_to_p3(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl)

  use micro_p3_utils, only: avg_diameter, &
                            rho_h2o, &
                            rho_h2os, &
                            qsmall, &
                            mincld, &
                            inv_cp

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  type(TimeLevel_t),  intent(in)            :: tl                       ! time level structure

  real(rl), dimension(np,np,nlev) :: u,v,w,T,p,dp,hommerho,zm,q1,q2,q3,q4,q5,q6,q7,q8,q9
  real(rl), dimension(np,np,nlev) :: pottemp, pnh, hommeexner, hommetemp, hommepres
  real(rl), dimension(np,np)      :: ps
  real(rl), dimension(np,np,nlevp):: zi, hommemu

#if 1

    !INTERNAL VARIABLES
    real(rl) :: dz(nlev)        !geometric layer thickness              m
    real(rl) :: cldliq(nlev)     !cloud liquid water mixing ratio        kg/kg
    real(rl) :: numliq(nlev)     !cloud liquid water drop concentraiton  #/kg
    real(rl) :: rain(nlev)       !rain water mixing ratio                kg/kg
    real(rl) :: numrain(nlev)    !rain water number concentration        #/kg
    real(rl) :: qv(nlev)         !water vapor mixing ratio               kg/kg
    real(rl) :: ice(nlev)        !total ice water mixing ratio           kg/kg
    real(rl) :: qm(nlev)      !rime ice mixing ratio                  kg/kg
    real(rl) :: numice(nlev)     !total ice crystal number concentration #/kg
    real(rl) :: rimvol(nlev)     !rime volume mixing ratio               m3/kg
    real(rl) :: temp(nlev)       !temperature copy needed for tendency   K
    real(rl) :: th(nlev)         !potential temperature                  K

    real(rl) :: precip_liq_surf         !precipitation rate, liquid             m s-1
    real(rl) :: precip_ice_surf         !precipitation rate, solid              m s-1

!out
    real(rl) :: rho_qi(nlev)  !bulk density of ice                    kg m-1

    real(rl) :: pres(nlev)       !pressure at midlevel                   hPa
    real(rl) :: pdel(nlev)

!out
    real(rl) :: qv2qi_depos_tend(nlev)

!out
    real(rl) :: precip_liq_flux(nlev+1)     !grid-box average rain flux (kg m^-2s^-1) nlevp
    real(rl) :: precip_ice_flux(nlev+1)     !grid-box average ice/snow flux (kg m^-2s^-1) nlevp

!out
    real(rl) :: rflx(nlev+1)     !grid-box average rain flux (kg m^-2s^-1) nlevp
    real(rl) :: sflx(nlev+1)     !grid-box average ice/snow flux (kg m^-2s^-1) nlevp
    real(rl) :: cflx(nlev+1)     !grid-box average cloud flux (kg m^-2s^-1) nlevp

    real(rl) :: exner(nlev)      !exner formula for converting between potential and normal temp

!in, set to 0
    real(rl) :: cld_frac_r(nlev)      !rain cloud fraction
    real(rl) :: cld_frac_l(nlev)      !liquid cloud fraction
    real(rl) :: cld_frac_i(nlev)      !ice cloud fraction

!out
    real(rl) :: tend_out(nlev,49) !microphysical tendencies

!out
    real(rl), dimension(nlev) :: liq_ice_exchange ! sum of liq-ice phase change tendenices
    real(rl), dimension(nlev) :: vap_liq_exchange ! sum of vap-liq phase change tendenices
    real(rl), dimension(nlev) :: vap_ice_exchange ! sum of vap-ice phase change tendenices

!out
    real(rl), dimension(nlev) :: diag_equiv_reflectivity,diag_ze_rain,diag_ze_ice ! equivalent reflectivity [dBz]

    !Prescribed CCN concentration
!in, set to 0
    real(rl), dimension(nlev) :: nccn_prescribed

    ! PBUF Variables
!in set to 0
    real(rl) :: ni_activated(nlev)     ! ice nucleation number
    real(rl) :: npccn(nlev)    ! liquid activation number tendency
    
!in, set to 0
!actually it is used as gamma(relvar), so it cannot be zero
    real(rl) :: relvar(nlev)    ! cloud liquid relative variance [-]

!out
    real(rl) :: qr_evap_tend(nlev) ! precipitation evaporation rate

!in, put into state, set to qv for now
    real(rl) :: qv_prev(nlev)   ! qv from previous p3_main call
    real(rl) :: t_prev(nlev)    ! t from previous p3_main call

    !! wetdep 
!out
    real(rl) :: precip_total_tend(nlev)        ! Total precipitation (rain + snow)
    real(rl) :: nevapr(nlev)       ! Evaporation of total precipitation (rain + snow)

    !! COSP simulator
!out
    real(rl) :: rel(nlev)          ! Liquid effective drop radius (microns)
    real(rl) :: rei(nlev)          ! Ice effective drop size (microns)

    !! radiation 
!not used?
    real(rl) :: dei(nlev)          ! Ice effective diameter (um)
!out
    real(rl) :: mu(nlev)           ! Size distribution shape parameter for radiation
    real(rl) :: lambdac(nlev)      ! Size distribution slope parameter for radiation

    ! Derived Variables
!in
    real(rl) :: col_location(3)  ! Array of column lon (index 1) and lat (index 2)

    ! variables for the CNT primary / heterogeneous freezing
    !real(rl), pointer :: frzimm(nlev)
    !real(rl), pointer :: frzcnt(nlev)
    !real(rl), pointer :: frzdep(nlev)
!in, set to 0
    real(rl) :: frzimm_in(nlev)
    real(rl) :: frzcnt_in(nlev)
    real(rl) :: frzdep_in(nlev)

    integer :: it                      !timestep counter    
    integer :: kts                     !closest level to TOM                   -
    integer :: kte                     !near surface level                     -
    integer :: ii, jj, ie, qind
    real(rl) :: dtime

!OPTIONS, remove them out later
!most of them are here https://docs.e3sm.org/E3SM/EAM/user-guide/namelist_parameters/#predicted-particle-properties 
    logical, parameter :: do_predict_nc     = .true.        !prognostic droplet concentration or not?
    logical, parameter :: do_subgrid_clouds = .false.       !use subgrid cloudiness in tendency calculations?
    logical, parameter :: do_prescribed_CCN = .false.       !see eam and micro_p3_init when this variable is true, there is a lot of init
    logical, parameter :: precip_off = .false.       
    real(rl), parameter ::         micro_nccons = 1.0 ! did not find this one anywhere
    real(rl), parameter ::         p3_autocon_coeff    = 30500.0  ! IN  autoconversion coefficient
    real(rl), parameter ::         p3_accret_coeff     = 117.25   ! IN  accretion coefficient
    real(rl), parameter ::         p3_qc_autocon_expon = 3.19     ! IN  autoconversion qc exponent
    real(rl), parameter ::         p3_nc_autocon_expon = -1.1     ! IN  autoconversion nc exponent
    real(rl), parameter ::         p3_qc_accret_expon  = 1.15     ! IN  autoconversion coefficient
    real(rl), parameter ::         p3_wbf_coeff        = 1.0      ! IN  WBF process coefficient
    real(rl), parameter ::         p3_mincdnc          = 20000000.0     ! IN  imposing minimal Nc 
    real(rl), parameter ::         p3_max_mean_rain_size  = 0.005 ! IN  max mean rain size
    real(rl), parameter ::         p3_embryonic_rain_size = 0.000025 ! IN  embryonic rain size for autoconversion


!should depend on a column, but we set to 0 for all of colns
    frzimm_in = 0.0
    frzcnt_in = 0.0
    frzdep_in = 0.0
    nccn_prescribed = 0.0
    ni_activated = 0.0
    npccn = 0.0 
    relvar = 0.1
    cld_frac_i = 1.0
    cld_frac_l = 1.0
    cld_frac_r = 1.0

    !precip_liq_surf = 0.0
    !precip_ice_surf = 0.0
    vap_liq_exchange = 0.0

    frzimm_in = 0.0
    frzcnt_in = 0.0
    frzdep_in = 0.0

    kte=nlev; kts=1;

!this may be wrong counter, do we call forcing only at remap?
!yes for "not flexible" routine, and the flexible one is disabled for now
    it = tl%nstep;
    dtime = tstep*dt_remap_factor ! figure it later

    !col_location uninited

    do ie = nets, nete

    !call get_state(u,v,w,T,p,dp,ps,hommerho,zm,zi,gravit,elem(ie),hvcoord,nt,ntQ)
    dp = elem(ie)%state%dp3d(:,:,:,nt)
    zi = elem(ie)%state%phinh_i(:,:,:,nt)/gravit

    ! get mixing ratios
    ! use qind to avoid compiler warnings when qsize_d<5
    qind=1;  q1  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp
    qind=2;  q2  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp
    qind=3;  q3  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp
    qind=4;  q4  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp
    qind=5;  q5  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp
    qind=6;  q6  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp
    qind=7;  q7  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp
    qind=8;  q8  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp
    qind=9;  q9  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp

    !get_temp uses Q for Rstar, not Qdp, so does not need ntQ
    call get_temperature(elem(ie),hommetemp,hvcoord,nt)

    call get_pottemp(elem(ie),pottemp,hvcoord,nt,ntQ)
    call get_hydro_pressure(hommepres,dp,hvcoord)
    call pnh_and_exner_from_eos(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,nt),&
          dp,elem(ie)%state%phinh_i(:,:,:,nt),pnh,hommeexner,hommemu)

    ! ensure positivity
    ! do we need this?
    where(q1<0); q1=0; endwhere
    where(q3<0); q3=0; endwhere

    do ii=1,np; do jj=1,np;

    cldliq  = q2(ii,jj,:)
    numliq  = q3(ii,jj,:)
    rain    = q4(ii,jj,:)
    numrain = q5(ii,jj,:)
    qv      = q1(ii,jj,:)  ! <!!!!! ---
    ice     = q6(ii,jj,:)
    qm      = q7(ii,jj,:)
    numice  = q8(ii,jj,:)
    rimvol  = q9(ii,jj,:)

    dz = zi(ii,jj,1:nlev) - zi(ii,jj,2:nlevp)
    !this is only virtual pottemp
    !th = elem(ie)%state%vtheta_dp(ii,jj,:,nt)/dp(ii,jj,:)
    th = pottemp(ii,jj,:)
!beyond me
    exner = 1./hommeexner(ii,jj,:)
    pdel = dp(ii,jj,:)
    pres = hommepres(ii,jj,:)

!temporary, set 'previous' vars to the current ones
!INVESTIGATE
!WRONG
    qv_prev = qv
    T_prev = hommetemp(ii,jj,:)

    call p3_main( &
         cldliq(kts:kte),     & ! INOUT  cloud, mass mixing ratio         kg kg-1
         numliq(kts:kte),     & ! INOUT  cloud, number mixing ratio       #  kg-1
         rain(kts:kte),       & ! INOUT  rain, mass mixing ratio          kg kg-1
         numrain(kts:kte),    & ! INOUT  rain, number mixing ratio        #  kg-1
         th(kts:kte),         & ! INOUT  potential temperature            K
         qv(kts:kte),         & ! INOUT  water vapor mixing ratio         kg kg-1
         dtime,                       & ! IN     model time step                  s
         ice(kts:kte),        & ! INOUT  ice, total mass mixing ratio     kg kg-1
         qm(kts:kte),      & ! INOUT  ice, rime mass mixing ratio      kg kg-1
         numice(kts:kte),     & ! INOUT  ice, total number mixing ratio   #  kg-1
         rimvol(kts:kte),     & ! INOUT  ice, rime volume mixing ratio    m3 kg-1
         pres(kts:kte),       & ! IN     pressure at cell midpoints       Pa
         dz(kts:kte),        & ! IN     vertical grid spacing            m
         npccn(kts:kte),      & ! IN ccn activation number tendency kg-1 s-1
         nccn_prescribed(kts:kte), & ! IN ccn prescribed concentration
         ni_activated(kts:kte),    & ! IN activated ice nuclei concentration kg-1
         frzimm_in(kts:kte), &
         frzcnt_in(kts:kte), & ! IN     CNT coupling
         frzdep_in(kts:kte), &
         relvar(kts:kte),     & ! IN cloud liquid relative variance
         it,                          & ! IN     time step counter NOTE: starts at 1 for first time step
         precip_liq_surf,            & ! OUT    surface liquid precip rate       m s-1
         precip_ice_surf,            & ! OUT    surface frozen precip rate       m s-1
         kts,                         & ! IN     vertical index lower bound       -
         kte,                         & ! IN     vertical index upper bound       -
         rel(kts:kte),        & ! OUT    effective radius, cloud          m
         rei(kts:kte),        & ! OUT    effective radius, ice            m
         rho_qi(kts:kte),  & ! OUT    bulk density of ice              kg m-3
         do_predict_nc,               & ! IN     .true.=prognostic Nc, .false.=specified Nc
         do_prescribed_CCN,           & ! IN
         p3_autocon_coeff,            & ! IN  autoconversion coefficient
         p3_accret_coeff,             & ! IN  accretion coefficient
         p3_qc_autocon_expon,         & ! IN  autoconversion qc exponent
         p3_nc_autocon_expon,         & ! IN  autoconversion nc exponent
         p3_qc_accret_expon,          & ! IN  autoconversion coefficient
         p3_wbf_coeff,                & ! IN  WBF process coefficient
         p3_mincdnc,                  & ! IN  imposing minimal Nc 
         p3_max_mean_rain_size,       & ! IN  max mean rain size
         p3_embryonic_rain_size,      & ! IN  embryonic rain size for autoconversion
         ! AaronDonahue new stuff
         pdel(kts:kte), & ! IN pressure level thickness for computing total mass
         exner(kts:kte),      & ! IN exner values
         qv2qi_depos_tend(kts:kte),    & ! OUT Deposition/sublimation rate of cloud ice 
         precip_total_tend(kts:kte),      & ! OUT Total precipitation (rain + snow)
         nevapr(kts:kte),     & ! OUT evaporation of total precipitation (rain + snow)
         qr_evap_tend(kts:kte),  & ! OUT rain evaporation
         precip_liq_flux(kts:kte+1),     & ! OUT grid-box average rain flux (kg m^-2s^-1) nlevp 
         precip_ice_flux(kts:kte+1),     & ! OUT grid-box average ice/snow flux (kgm^-2 s^-1) nlevp
         rflx(kts:kte+1),     & ! OUT grid-box average rain flux (kg m^-2 s^-1) nlevp 
         sflx(kts:kte+1),     & ! OUT grid-box average ice/snow flux (kgm^-2 s^-1) nlevp
         cflx(kts:kte+1),     & ! OUT grid-box average cld droplet flux (kgm^-2 s^-1) nlevp
         cld_frac_r(kts:kte),      & ! IN rain cloud fraction
         cld_frac_l(kts:kte),      & ! IN liquid cloud fraction
         cld_frac_i(kts:kte),      & ! IN ice cloud fraction
         tend_out(kts:kte,:), & ! OUT p3 microphysics tendencies
         mu(kts:kte),         & ! OUT Size distribution shape parameter for radiation
         lambdac(kts:kte),    & ! OUT Size distribution slope parameter for radiation
         liq_ice_exchange(kts:kte),& ! OUT sum of liq-ice phase change tendenices   
         vap_liq_exchange(kts:kte),& ! OUT sun of vap-liq phase change tendencies
         vap_ice_exchange(kts:kte),& ! OUT sum of vap-ice phase change tendencies
         qv_prev(kts:kte),         & ! IN  qv at end of prev p3_main call   kg kg-1
         t_prev(kts:kte),          & ! IN  t at end of prev p3_main call    K
         col_location(:3),         & ! IN column locations
         precip_off,                       & ! IN Option to turn precip (liquid) off
         micro_nccons,                     & ! IN Option for constant droplet concentration
         diag_equiv_reflectivity(kts:kte), & !OUT equivalent reflectivity (rain + ice) [dBz]
         diag_ze_rain(kts:kte),diag_ze_ice(kts:kte)) !OUT equivalent reflectivity for rain and ice [dBz]

#endif
#if 1
    temp(:nlev) = th(:nlev)/exner(:nlev)
    elem(ie)%derived%FT(ii,jj,:) = (temp(:) - hommetemp(ii,jj,:))/dtime
    ! from above
    !cldliq  = q2(ii,jj,:)
    !numliq  = q3(ii,jj,:)
    !rain    = q4(ii,jj,:)
    !numrain = q5(ii,jj,:)
    !qv      = q1(ii,jj,:)  ! <!!!!! ---
    !ice     = q6(ii,jj,:)
    !qm      = q7(ii,jj,:)
    !numice  = q8(ii,jj,:)
    !rimvol  = q9(ii,jj,:)

    qind=1; elem(ie)%derived%FQ(ii,jj,:,qind)  = pdel(:)* ( qv(:)      - q1(ii,jj,:) )/dtime
    qind=2; elem(ie)%derived%FQ(ii,jj,:,qind)  = pdel(:)*( cldliq(:)  - q2(ii,jj,:) )/dtime
    qind=3; elem(ie)%derived%FQ(ii,jj,:,qind)  = pdel(:)*( numliq(:)  - q3(ii,jj,:) )/dtime
    qind=4; elem(ie)%derived%FQ(ii,jj,:,qind)  = pdel(:)*( rain(:)    - q4(ii,jj,:) )/dtime
    qind=5; elem(ie)%derived%FQ(ii,jj,:,qind)  = pdel(:)*( numrain(:) - q5(ii,jj,:) )/dtime
    qind=6; elem(ie)%derived%FQ(ii,jj,:,qind)  = pdel(:)*( ice(:)     - q6(ii,jj,:) )/dtime
    qind=7; elem(ie)%derived%FQ(ii,jj,:,qind)  = pdel(:)*( qm(:)      - q7(ii,jj,:) )/dtime
    qind=8; elem(ie)%derived%FQ(ii,jj,:,qind)  = pdel(:)*( numice(:)  - q8(ii,jj,:) )/dtime
    qind=9; elem(ie)%derived%FQ(ii,jj,:,qind)  = pdel(:)*( rimvol(:)  - q9(ii,jj,:) )/dtime

    !convert flux m/s into [dp] units, which is [kg / m / s^2]
    elem(ie)%accum%precip_mass(ii,jj) = (precip_liq_surf + precip_ice_surf)*dtime*rho_liquidH20*gravit


!if( (precip_liq_surf > 0.0) .or. (precip_ice_surf > 0.0) ) then
!print *, 'precip', precip_liq_surf, precip_ice_surf
!end if
#if 0

!if (ie == 1) then
!print *, 'precip', precip_liq_surf, precip_ice_surf
!print *, 'FT', elem(ie)%derived%FT(ii,jj,:)
!do qind = 1, 9
!!if(maxval(abs( elem(ie)%derived%FQ(ii,jj,:,qind) ),1)> 0.0)then
!print *, qind, maxval(abs(elem(ie)%derived%FQ(ii,jj,:,qind)),1)
!!endif
!enddo
!endif
!stop
#endif




#endif

    enddo; enddo; ! ii, jj
    elem(ie)%derived%FM(:,:,:,:) = 0.0
    elem(ie)%derived%FQps(:,:) = 0.0

    enddo ! ie

end subroutine interface_to_p3






end module p3phys
#endif
                                                                                                                                                                        





