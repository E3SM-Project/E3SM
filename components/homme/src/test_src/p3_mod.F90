#ifndef CAM
#include "config.h"

module p3phys

! Implementation of the dcmip2012 dycore tests for the preqx dynamics target

use control_mod,          only: theta_hydrostatic_mode,&
                   case_planar_bubble, bubble_prec_type, bubble_rj_cpstar_hy, bubble_rj_cpstar_nh, &
                   bubble_rj_eamcpdry, bubble_rj_eamcpstar
use dimensions_mod,       only: np, nlev, nlevp , qsize, qsize_d, nelemd
use element_mod,          only: element_t
use element_state,        only: nt=>timelevels
use hybrid_mod,           only: hybrid_t
use kinds,                only: rl=>real_kind, iulog
use parallel_mod,         only: abortmp,iam
use reduction_mod,        only: parallelmax, parallelmin

use physical_constants,   only: bubble_const1, bubble_const2, bubble_const3, bubble_const4, &
                                bubble_t0_const, bubble_epsilo, bubble_e0, &
                                latvap, latice, &
                                gravit=>g, p0, &
                                rdry=>rgas, cpdry=>cp, cpv=>cpwater_vapor, cl, cvdry, cvv, &
                                rvapor=>rwater_vapor, rhow
use hybvcoord_mod,        only: hvcoord_t
use time_mod,             only: time_at, TimeLevel_t

implicit none

public :: interface_to_p3

contains

subroutine interface_to_p3(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl)

    use micro_p3,       only: p3_main
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

#if 0

    !INTERNAL VARIABLES
    real(rtype) :: dz(pcols,pver)        !geometric layer thickness              m
    real(rtype) :: cldliq(pcols,pver)     !cloud liquid water mixing ratio        kg/kg
    real(rtype) :: numliq(pcols,pver)     !cloud liquid water drop concentraiton  #/kg
    real(rtype) :: rain(pcols,pver)       !rain water mixing ratio                kg/kg
    real(rtype) :: numrain(pcols,pver)    !rain water number concentration        #/kg
    real(rtype) :: qv(pcols,pver)         !water vapor mixing ratio               kg/kg
    real(rtype) :: ice(pcols,pver)        !total ice water mixing ratio           kg/kg
    real(rtype) :: qm(pcols,pver)      !rime ice mixing ratio                  kg/kg
    real(rtype) :: numice(pcols,pver)     !total ice crystal number concentration #/kg
    real(rtype) :: rimvol(pcols,pver)     !rime volume mixing ratio               m3/kg
    real(rtype) :: temp(pcols,pver)       !temperature copy needed for tendency   K
    real(rtype) :: th(pcols,pver)         !potential temperature                  K
    real(rtype) :: precip_liq_surf(pcols)         !precipitation rate, liquid             m s-1
    real(rtype) :: precip_ice_surf(pcols)         !precipitation rate, solid              m s-1

    real(rtype) :: rho_qi(pcols,pver)  !bulk density of ice                    kg m-1
    real(rtype) :: pres(pcols,pver)       !pressure at midlevel                   hPa
    real(rtype) :: qv2qi_depos_tend(pcols,pver)
    real(rtype) :: precip_liq_flux(pcols,pver+1)     !grid-box average rain flux (kg m^-2s^-1) pverp
    real(rtype) :: precip_ice_flux(pcols,pver+1)     !grid-box average ice/snow flux (kg m^-2s^-1) pverp
    real(rtype) :: rflx(pcols,pver+1)     !grid-box average rain flux (kg m^-2s^-1) pverp
    real(rtype) :: sflx(pcols,pver+1)     !grid-box average ice/snow flux (kg m^-2s^-1) pverp
    real(rtype) :: cflx(pcols,pver+1)     !grid-box average cloud flux (kg m^-2s^-1) pverp
    real(rtype) :: exner(pcols,pver)      !exner formula for converting between potential and normal temp
    real(rtype) :: cld_frac_r(pcols,pver)      !rain cloud fraction
    real(rtype) :: cld_frac_l(pcols,pver)      !liquid cloud fraction
    real(rtype) :: cld_frac_i(pcols,pver)      !ice cloud fraction
    real(rtype) :: tend_out(pcols,pver,49) !microphysical tendencies
    real(rtype), dimension(pcols,pver) :: liq_ice_exchange ! sum of liq-ice phase change tendenices
    real(rtype), dimension(pcols,pver) :: vap_liq_exchange ! sum of vap-liq phase change tendenices
    real(rtype), dimension(pcols,pver) :: vap_ice_exchange ! sum of vap-ice phase change tendenices
    real(rtype), dimension(pcols,pver) :: diag_equiv_reflectivity,diag_ze_rain,diag_ze_ice ! equivalent reflectivity [dBz]

    !Prescribed CCN concentration
    real(rtype), dimension(pcols,pver) :: nccn_prescribed

    ! PBUF Variables
    real(rtype), pointer :: ast(:,:)      ! Relative humidity cloud fraction
    real(rtype), pointer :: ni_activated(:,:)     ! ice nucleation number
    real(rtype), pointer :: npccn(:,:)    ! liquid activation number tendency
    real(rtype), pointer :: cmeliq(:,:)
    !!
    real(rtype), pointer :: prec_str(:)    ! [Total] Sfc flux of precip from stratiform [ m/s ]
    real(rtype), pointer :: prec_sed(:)    ! Surface flux of total cloud water from sedimentation
    real(rtype), pointer :: prec_pcw(:)    ! Sfc flux of precip from microphysics [ m/s ]
    real(rtype), pointer :: snow_str(:)    ! [Total] Sfc flux of snow from stratiform   [ m/s ]
    real(rtype), pointer :: snow_pcw(:)    ! Sfc flux of snow from microphysics [ m/s ]
    real(rtype), pointer :: snow_sed(:)    ! Surface flux of cloud ice from sedimentation
    real(rtype), pointer :: relvar(:,:)    ! cloud liquid relative variance [-]
    real(rtype), pointer :: cldo(:,:)      ! Old cloud fraction
    real(rtype), pointer :: qr_evap_tend(:,:) ! precipitation evaporation rate
    real(rtype), pointer :: qv_prev(:,:)   ! qv from previous p3_main call
    real(rtype), pointer :: t_prev(:,:)    ! t from previous p3_main call
    !! wetdep 
    real(rtype), pointer :: qme(:,:)
    real(rtype), pointer :: precip_total_tend(:,:)        ! Total precipitation (rain + snow)
    real(rtype), pointer :: nevapr(:,:)       ! Evaporation of total precipitation (rain + snow)
    !! COSP simulator
    real(rtype), pointer :: rel(:,:)          ! Liquid effective drop radius (microns)
    real(rtype), pointer :: rei(:,:)          ! Ice effective drop size (microns)
    real(rtype), pointer :: flxprc(:,:)     ! P3 grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
    real(rtype), pointer :: flxsnw(:,:)     ! P3 grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
    real(rtype), pointer :: reffrain(:,:)   ! P3 diagnostic rain effective radius (um)
    real(rtype), pointer :: reffsnow(:,:)   ! P3 diagnostic snow effective radius (um)
    real(rtype), pointer :: cvreffliq(:,:)    ! convective cloud liquid effective radius (um)
    real(rtype), pointer :: cvreffice(:,:)    ! convective cloud ice effective radius (um)
    !! radiation 
    real(rtype), pointer :: dei(:,:)          ! Ice effective diameter (um)
    real(rtype), pointer :: mu(:,:)           ! Size distribution shape parameter for radiation
    real(rtype), pointer :: lambdac(:,:)      ! Size distribution slope parameter for radiation

    ! Derived Variables
    real(rtype) :: icimrst(pcols,pver) ! stratus ice mixing ratio - on grid
    real(rtype) :: icwmrst(pcols,pver) ! stratus water mixing ratio - on grid
    real(rtype) :: rho(pcols,pver)
    real(rtype) :: drout2(pcols,pver)
    real(rtype) :: reff_rain(pcols,pver)
    real(rtype) :: col_location(pcols,3),tmp_loc(pcols)  ! Array of column lon (index 1) and lat (index 2)
    integer     :: tmpi_loc(pcols) ! Global column index temp array

    ! Variables used for microphysics output
    real(rtype) :: aqrain(pcols,pver)
    real(rtype) :: anrain(pcols,pver)
    real(rtype) :: nfice(pcols,pver)
    real(rtype) :: efcout(pcols,pver)
    real(rtype) :: efiout(pcols,pver)
    real(rtype) :: ncout(pcols,pver)
    real(rtype) :: niout(pcols,pver)
    real(rtype) :: freqr(pcols,pver)
    real(rtype) :: freql(pcols,pver)
    real(rtype) :: freqi(pcols,pver)
    real(rtype) :: cdnumc(pcols)
    real(rtype) :: icinc(pcols,pver)
    real(rtype) :: icwnc(pcols,pver)

    ! variables for the CNT primary / heterogeneous freezing
    !real(rtype), pointer :: frzimm(:,:)
    !real(rtype), pointer :: frzcnt(:,:)
    !real(rtype), pointer :: frzdep(:,:)
    real(rtype) :: frzimm_in(pcols,pver)
    real(rtype) :: frzcnt_in(pcols,pver)
    real(rtype) :: frzdep_in(pcols,pver)

    integer :: it                      !timestep counter    
    integer :: kts                     !closest level to TOM                   -
    integer :: kte                     !near surface level                     -

    logical :: do_predict_nc           !prognostic droplet concentration or not?
    logical :: do_subgrid_clouds       !use subgrid cloudiness in tendency calculations?
    integer :: icol, ncol, k
    integer :: psetcols, lchnk
    integer :: itim_old


      frzimm_in = 0.0_rtype
      frzcnt_in = 0.0_rtype
      frzdep_in = 0.0_rtype

    !compute exner here
    !compute dz, theta

    kte=nlev; kts=1;
    it = tl%nstep;
    !col_location uninited
!lets do wv=1
!2=cldliq, 3=numliq; 4=rain; 5=numrain;6=cldice; 7=numice; 8=cldrim; 9=rimvol

    cldliq  = state%q(:,:,2)
    numliq  = state%q(:,:,3)
    rain    = state%q(:,:,4)
    numrain = state%q(:,:,5)
    qv      = state%q(:,:,1)
    ice     = state%q(:,:,6)
    qm      = state%q(:,:,8) !Aaron, changed ixqm to ixcldrim to match Kai's code
    numice  = state%q(:,:,7)
    rimvol  = state%q(:,:,9)

    ! Determine the cloud fraction and precip cover
    cld_frac_i(:) = 1.0_rtype
    cld_frac_l(:) = 1.0_rtype
    cld_frac_r(:) = 1.0_rtype

    precip_liq_surf = 0.0_rtype
    precip_ice_surf = 0.0_rtype
    prec_pcw = 0.0_rtype
    snow_pcw = 0.0_rtype
    vap_liq_exchange = 0.0_rtype

    call p3_main( &
         cldliq(its:ite,kts:kte),     & ! INOUT  cloud, mass mixing ratio         kg kg-1
         numliq(its:ite,kts:kte),     & ! INOUT  cloud, number mixing ratio       #  kg-1
         rain(its:ite,kts:kte),       & ! INOUT  rain, mass mixing ratio          kg kg-1
         numrain(its:ite,kts:kte),    & ! INOUT  rain, number mixing ratio        #  kg-1
         th(its:ite,kts:kte),         & ! INOUT  potential temperature            K
         qv(its:ite,kts:kte),         & ! INOUT  water vapor mixing ratio         kg kg-1
         dtime,                       & ! IN     model time step                  s
         ice(its:ite,kts:kte),        & ! INOUT  ice, total mass mixing ratio     kg kg-1
         qm(its:ite,kts:kte),      & ! INOUT  ice, rime mass mixing ratio      kg kg-1
         numice(its:ite,kts:kte),     & ! INOUT  ice, total number mixing ratio   #  kg-1
         rimvol(its:ite,kts:kte),     & ! INOUT  ice, rime volume mixing ratio    m3 kg-1
         pres(its:ite,kts:kte),       & ! IN     pressure at cell midpoints       Pa
         dz(its:ite,kts:kte),        & ! IN     vertical grid spacing            m
         npccn(its:ite,kts:kte),      & ! IN ccn activation number tendency kg-1 s-1
         nccn_prescribed(its:ite,kts:kte), & ! IN ccn prescribed concentration
         ni_activated(its:ite,kts:kte),    & ! IN activated ice nuclei concentration kg-1
         frzimm_in(its:ite,kts:kte), &
         frzcnt_in(its:ite,kts:kte), & ! IN     CNT coupling
         frzdep_in(its:ite,kts:kte), &
         relvar(its:ite,kts:kte),     & ! IN cloud liquid relative variance
         it,                          & ! IN     time step counter NOTE: starts at 1 for first time step
         precip_liq_surf(its:ite),            & ! OUT    surface liquid precip rate       m s-1
         precip_ice_surf(its:ite),            & ! OUT    surface frozen precip rate       m s-1
         its,                         & ! IN     horizontal index lower bound     -
         ite,                         & ! IN     horizontal index upper bound     -
         kts,                         & ! IN     vertical index lower bound       -
         kte,                         & ! IN     vertical index upper bound       -
         rel(its:ite,kts:kte),        & ! OUT    effective radius, cloud          m
         rei(its:ite,kts:kte),        & ! OUT    effective radius, ice            m
         rho_qi(its:ite,kts:kte),  & ! OUT    bulk density of ice              kg m-3
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
         state%pdel(its:ite,kts:kte), & ! IN pressure level thickness for computing total mass
         exner(its:ite,kts:kte),      & ! IN exner values
         qv2qi_depos_tend(its:ite,kts:kte),    & ! OUT Deposition/sublimation rate of cloud ice 
         precip_total_tend(its:ite,kts:kte),      & ! OUT Total precipitation (rain + snow)
         nevapr(its:ite,kts:kte),     & ! OUT evaporation of total precipitation (rain + snow)
         qr_evap_tend(its:ite,kts:kte),  & ! OUT rain evaporation
         precip_liq_flux(its:ite,kts:kte+1),     & ! OUT grid-box average rain flux (kg m^-2s^-1) pverp 
         precip_ice_flux(its:ite,kts:kte+1),     & ! OUT grid-box average ice/snow flux (kgm^-2 s^-1) pverp
         rflx(its:ite,kts:kte+1),     & ! OUT grid-box average rain flux (kg m^-2 s^-1) pverp 
         sflx(its:ite,kts:kte+1),     & ! OUT grid-box average ice/snow flux (kgm^-2 s^-1) pverp
         cflx(its:ite,kts:kte+1),     & ! OUT grid-box average cld droplet flux (kgm^-2 s^-1) pverp
         cld_frac_r(its:ite,kts:kte),      & ! IN rain cloud fraction
         cld_frac_l(its:ite,kts:kte),      & ! IN liquid cloud fraction
         cld_frac_i(its:ite,kts:kte),      & ! IN ice cloud fraction
         tend_out(its:ite,kts:kte,:), & ! OUT p3 microphysics tendencies
         mu(its:ite,kts:kte),         & ! OUT Size distribution shape parameter for radiation
         lambdac(its:ite,kts:kte),    & ! OUT Size distribution slope parameter for radiation
         liq_ice_exchange(its:ite,kts:kte),& ! OUT sum of liq-ice phase change tendenices   
         vap_liq_exchange(its:ite,kts:kte),& ! OUT sun of vap-liq phase change tendencies
         vap_ice_exchange(its:ite,kts:kte),& ! OUT sum of vap-ice phase change tendencies
         qv_prev(its:ite,kts:kte),         & ! IN  qv at end of prev p3_main call   kg kg-1
         t_prev(its:ite,kts:kte),          & ! IN  t at end of prev p3_main call    K
         col_location(its:ite,:3),         & ! IN column locations
         precip_off,                       & ! IN Option to turn precip (liquid) off
         micro_nccons,                     & ! IN Option for constant droplet concentration
         diag_equiv_reflectivity(its:ite,kts:kte), & !OUT equivalent reflectivity (rain + ice) [dBz]
         diag_ze_rain(its:ite,kts:kte),diag_ze_ice(its:ite,kts:kte)) !OUT equivalent reflectivity for rain and ice [dBz]




    temp(:ncol,:pver) = th(:ncol,:pver)/exner(:ncol,:pver)
    ptend%s(:ncol,:pver)           = ( temp(:ncol,:pver) - state%t(:ncol,:pver) )/dtime
    ptend%q(:ncol,:pver,1)         = ( max(0._rtype,qv(:ncol,:pver)     ) - state%q(:ncol,:pver,1)         )/dtime
    ptend%q(:ncol,:pver,ixcldliq)  = ( max(0._rtype,cldliq(:ncol,:pver) ) - state%q(:ncol,:pver,ixcldliq)  )/dtime
    ptend%q(:ncol,:pver,ixnumliq)  = ( max(0._rtype,numliq(:ncol,:pver) ) - state%q(:ncol,:pver,ixnumliq)  )/dtime
    ptend%q(:ncol,:pver,ixrain)    = ( max(0._rtype,rain(:ncol,:pver)   ) - state%q(:ncol,:pver,ixrain)    )/dtime
    ptend%q(:ncol,:pver,ixnumrain) = ( max(0._rtype,numrain(:ncol,:pver)) - state%q(:ncol,:pver,ixnumrain) )/dtime
    ptend%q(:ncol,:pver,ixcldice)  = ( max(0._rtype,ice(:ncol,:pver)    ) - state%q(:ncol,:pver,ixcldice)  )/dtime
    ptend%q(:ncol,:pver,ixnumice)  = ( max(0._rtype,numice(:ncol,:pver) ) - state%q(:ncol,:pver,ixnumice)  )/dtime
    ptend%q(:ncol,:pver,ixcldrim)  = ( max(0._rtype,qm(:ncol,:pver)  ) - state%q(:ncol,:pver,ixcldrim)  )/dtime
    ptend%q(:ncol,:pver,ixrimvol)  = ( max(0._rtype,rimvol(:ncol,:pver) ) - state%q(:ncol,:pver,ixrimvol)  )/dtime


!disable everything
#endif

end subroutine interface_to_p3






end module p3phys
#endif
                                                                                                                                                                        





