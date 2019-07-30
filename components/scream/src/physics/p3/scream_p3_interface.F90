#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module scream_p3_interface_mod

  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool,C_NULL_CHAR
  use micro_p3_utils, only: rtype

  public :: p3_init_f90
  public :: p3_main_f90
  public :: p3_finalize_f90

  integer :: test

  integer :: pcols
  integer :: pver

contains

  !====================================================================!
  subroutine p3_init_f90 () bind(c)
    use micro_p3,       only: p3_init
    use micro_p3_utils, only: micro_p3_utils_init

    character(len=16)  :: micro_p3_tableversion = "4"
    character(len=100) :: micro_p3_lookup_dir = "/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops"
    real               :: cpair  =    1004.64000000000
    real               :: rair   =    287.042311365049
    real               :: rh2o   =    461.504639820160
    real               :: rhoh2o =    1000.00000000000
    real               :: mwh2o  =    18.0160000000000
    real               :: mwdry  =    28.9660000000000
    real               :: gravit =    9.80616000000000
    real               :: latvap =    2501000.00000000
    real               :: latice =    333700.000000000
    real               :: cpliq  =    4188.00000000000
    real               :: tmelt  =    273.150000000000
    real               :: pi     =    3.14159265358979

    call micro_p3_utils_init(cpair,rair,rh2o,rhoh2o,mwh2o,mwdry,gravit,latvap,latice, &
             cpliq,tmelt,pi,0,.false.)
    call p3_init(micro_p3_lookup_dir,micro_p3_tableversion)

    pcols = 1
    pver  = 72
    test = 0
    print *, 'P3 init = ', test, pcols, pver

  end subroutine p3_init_f90
  !====================================================================!
  subroutine p3_main_f90 () bind(c)
    use micro_p3,       only: p3_main

    real :: dtime ! Timestep, should eventually be an input
    !INTERNAL VARIABLES
    real :: th_old(pcols,pver)     !potential temperature from last step   K
    real :: th(pcols,pver)         !potential temperature  K
    real :: qv_old(pcols,pver)     !water vapor from last step             kg/kg
    real :: ssat(pcols,pver)       !supersaturated mixing ratio            kg/kg
    real :: dzq(pcols,pver)        !geometric layer thickness              m
    real :: cldliq(pcols,pver)     !cloud liquid water mixing ratio        kg/kg
    real :: numliq(pcols,pver)     !cloud liquid water drop concentraiton  #/kg
    real :: rain(pcols,pver)       !rain water mixing ratio                kg/kg
    real :: numrain(pcols,pver)    !rain water number concentration        #/kg
    real :: qv(pcols,pver)         !water vapor mixing ratio               kg/kg
    real :: ice(pcols,pver)        !total ice water mixing ratio           kg/kg
    real :: qirim(pcols,pver)      !rime ice mixing ratio                  kg/kg
    real :: numice(pcols,pver)     !total ice crystal number concentration #/kg
    real :: rimvol(pcols,pver)     !rime volume mixing ratio               m3/kg
    real :: temp(pcols,pver)       !potential temperature                  K
    real :: rim(pcols,pver)        !rime mixing ratio                      kg/kg
    real :: prt_liq(pcols)         !precipitation rate, liquid             m s-1
    real :: prt_sol(pcols)         !precipitation rate, solid              m s-1
    real :: diag_ze(pcols,pver)    !equivalent reflectivity                dBZ
    real :: diag_effc(pcols,pver)  !effective radius, cloud                m
    real :: diag_effi(pcols,pver)  !effective radius, ice                  m
    real :: diag_vmi(pcols,pver)   !mass-weighted fall speed of ice        m s-1
    real :: diag_di(pcols,pver)    !mean diameter of ice                   m
    real :: diag_rhoi(pcols,pver)  !bulk density of ice                    kg m-1
    real :: pres(pcols,pver)       !pressure at midlevel                   hPa
    real :: rflx(pcols,pver+1)     !grid-box average rain flux (kg m^-2s^-1) pverp
    real :: sflx(pcols,pver+1)     !grid-box average ice/snow flux (kg m^-2s^-1) pverp
    real :: exner(pcols,pver)      !exner formula for converting between potential and normal temp
    real :: rcldm(pcols,pver)      !rain cloud fraction
    real :: lcldm(pcols,pver)      !liquid cloud fraction
    real :: icldm(pcols,pver)      !ice cloud fraction
    real :: tend_out(pcols,pver,49) !microphysical tendencies
    real :: npccn(pcols,pver)      !liq. activation number tendency
    real :: naai(pcols,pver)       !ice nucleation number
    real :: rel(pcols,pver)        !liq. effective drop radius (microns)
    real :: rei(pcols,pver)        !ice effective drop radius (microns)
    real :: pdel(pcols,pver)       !pressure thickness
    real :: cmeiout(pcols,pver)    !deposition/sublimation rate of cloud ice
    real :: prain(pcols,pver)      !total precip
    real :: nevapr(pcols,pver)     !evap. of total precip
    real :: prer_evap(pcols,pver)  !rain evaporation
    real :: pratot(pcols,pver)     !accretion of cloud by rain
    real :: prctot(pcols,pver)     !autoconversion of cloud by rain
    real :: mu(pcols,pver)         !Size distribution shape parameter for radiation
    real :: lambdac(pcols,pver)    !Size distribution slope parameter for radiation

    ! For rrtmg optics. specified distribution.
    real, parameter :: dcon   = 25.e-6_rtype      ! Convective size distribution effective radius (um)
    real, parameter :: mucon  = 5.3_rtype         ! Convective size distribution shape parameter
    real, parameter :: deicon = 50._rtype         ! Convective ice effective diameter (um)

    integer :: icol, k, ncol

    integer :: it, its, ite, kts, kte
    logical :: log_predictNc = .true.

    ncol = pcols
    dtime = 600.0_rtype

    ! DEAL WITH SSAT
    !==============
    !ssat (supersaturated mixing ratio measured in kg/kg) can be prognosed
    !or diagnosed in p3 depending on p3's hardcoded log_predictSsat parameter.
    !ssat is an intent inout variable, but when log_predictSsat is false 
    !(as we will set it initially), ssat is overwritten rather than used by p3.
    !Thus it shouldn't matter what ssat is, so we give it -999.

    ssat(:,:) = -999.0_rtype

    ! COMPUTE GEOMETRIC THICKNESS OF GRID
    !==============
    exner(:ncol,:pver) = 1.0_rtype!1._rtype/((state%pmid(:ncol,:pver)*1.e-5_rtype)**(rair*inv_cp))
    do icol = 1,ncol
       do k = 1,pver
! Note: dzq is calculated in the opposite direction that pdel is calculated,
! thus when considering any dp/dz calculation we must also change the sign.
          dzq(icol,k) = 1.0_rtype   !state%zi(icol,k) - state%zi(icol,k+1)
          th(icol,k)  = 270.0_rtype !state%t(icol,k)*exner(icol,k) !/(state%pmid(icol,k)*1.e-5)**(rd*inv_cp) 
          pdel(icol,k)  = 100.0_rtype ! should be changed to come from model state.
       end do
    end do

    ! MAKE LOCAL COPIES OF VARS MODIFIED BY P3
    !==============
    !local copies are needed because state is passed into this routine as intent=in
    !while P3 seeks to modify state variables in-place. Also, we need a copy of 
    !old values in order to back out ptend values later. Traditionally, a local copy 
    !is created by copying the whole state. It is much cheaper to just copy the 
    !variables we need. 
    
    cldliq  = 1.0e-6_rtype!state%q(:,:,ixcldliq)
    numliq  = 1.0e6_rtype!state%q(:,:,ixnumliq)
    rain    = 1.0e-5_rtype!state%q(:,:,ixrain)
    numrain = 1.0e5_rtype!state%q(:,:,ixnumrain)
    qv      = 1.0e-4_rtype!state%q(:,:,1)
    ice     = 1.0e-7_rtype!state%q(:,:,ixcldice)
    qirim   = 1.0e-8_rtype!state%q(:,:,ixcldrim) !Aaron, changed ixqirim to ixcldrim to match Kai's code
    numice  = 1.0e5_rtype!state%q(:,:,ixnumice)
    rimvol  = 1.0e4_rtype!state%q(:,:,ixrimvol)
    its     = 1
    ite     = ncol
    kts     = 1
    kte     = pver
    pres    = 1.0_rtype!state%pmid(:,:)
    ! Initialize the raidation dependent variables.
    mu      = mucon
    lambdac = (mucon + 1._rtype)/dcon
    dei     = deicon
    ! Determine the cloud fraction and precip cover
    icldm(:,:) = 1.0_rtype
    lcldm(:,:) = 1.0_rtype
    rcldm(:,:) = 1.0_rtype
    ! CALL P3
    !==============
    call p3_main( &
         cldliq(its:ite,kts:kte),     & ! INOUT  cloud, mass mixing ratio         kg kg-1
         numliq(its:ite,kts:kte),     & ! INOUT  cloud, number mixing ratio       #  kg-1
         rain(its:ite,kts:kte),       & ! INOUT  rain, mass mixing ratio          kg kg-1
         numrain(its:ite,kts:kte),    & ! INOUT  rain, number mixing ratio        #  kg-1
         th(its:ite,kts:kte),         & ! INOUT  potential temperature            K
         qv(its:ite,kts:kte),         & ! INOUT  water vapor mixing ratio         kg kg-1
         dtime,                       & ! IN     model time step                  s
         ice(its:ite,kts:kte),        & ! INOUT  ice, total mass mixing ratio     kg kg-1
         qirim(its:ite,kts:kte),      & ! INOUT  ice, rime mass mixing ratio      kg kg-1
         numice(its:ite,kts:kte),     & ! INOUT  ice, total number mixing ratio   #  kg-1
         rimvol(its:ite,kts:kte),     & ! INOUT  ice, rime volume mixing ratio    m3 kg-1
         pres(its:ite,kts:kte),       & ! IN     pressure at cell midpoints       Pa
         dzq(its:ite,kts:kte),        & ! IN     vertical grid spacing            m
         npccn(its:ite,kts:kte),      & ! IN ccn activation number tendency kg-1 s-1
         naai(its:ite,kts:kte),       & ! IN activated ice nuclei concentration kg-1
         it,                          & ! IN     time step counter NOTE: starts at 1 for first time step
         prt_liq(its:ite),            & ! OUT    surface liquid precip rate       m s-1
         prt_sol(its:ite),            & ! OUT    surface frozen precip rate       m s-1
         its,                         & ! IN     horizontal index lower bound     -
         ite,                         & ! IN     horizontal index upper bound     -
         kts,                         & ! IN     vertical index lower bound       -
         kte,                         & ! IN     vertical index upper bound       -
         diag_ze(its:ite,kts:kte),    & ! OUT    equivalent reflectivity          dBZ  UNUSED?
         rel(its:ite,kts:kte),        & ! OUT    effective radius, cloud          m
         rei(its:ite,kts:kte),        & ! OUT    effective radius, ice            m
         diag_vmi(its:ite,kts:kte),   & ! OUT    mass-weighted fall speed of ice  m s-1
         diag_di(its:ite,kts:kte),    & ! OUT    mean diameter of ice             m
         diag_rhoi(its:ite,kts:kte),  & ! OUT    bulk density of ice              kg m-3
         log_predictNc,               & ! IN     .true.=prognostic Nc, .false.=specified Nc
         ! AaronDonahue new stuff
         pdel(its:ite,kts:kte), & ! IN pressure level thickness for computing total mass
         exner(its:ite,kts:kte),      & ! IN exner values
         cmeiout(its:ite,kts:kte),    & ! OUT Deposition/sublimation rate of cloud ice 
         prain(its:ite,kts:kte),      & ! OUT Total precipitation (rain + snow)
         nevapr(its:ite,kts:kte),     & ! OUT evaporation of total precipitation (rain + snow)
         prer_evap(its:ite,kts:kte),  & ! OUT rain evaporation
         rflx(its:ite,kts:kte+1),     & ! OUT grid-box average rain flux (kg m^-2s^-1) pverp 
         sflx(its:ite,kts:kte+1),     & ! OUT grid-box average ice/snow flux (kgm^-2 s^-1) pverp
         rcldm(its:ite,kts:kte),      & ! IN rain cloud fraction
         lcldm(its:ite,kts:kte),      & ! IN liquid cloud fraction
         icldm(its:ite,kts:kte),      & ! IN ice cloud fraction
         pratot(its:ite,kts:kte),     & ! OUT accretion of cloud by rain
         prctot(its:ite,kts:kte),     & ! OUT autoconversion of cloud by rain
         tend_out(its:ite,kts:kte,:), & ! OUT p3 microphysics tendencies
         mu(its:ite,kts:kte),         & ! OUT Size distribution shape parameter for radiation
         lambdac(its:ite,kts:kte)     & ! OUT Size distribution slope parameter for radiation
         )

    test = test + 1
    print *, 'P3 run = ', test

  end subroutine p3_main_f90
  !====================================================================!
  subroutine p3_finalize_f90 () bind(c)

    test = -999
    print *, 'P3 final = ', test

  end subroutine p3_finalize_f90
  !====================================================================!

end module scream_p3_interface_mod
