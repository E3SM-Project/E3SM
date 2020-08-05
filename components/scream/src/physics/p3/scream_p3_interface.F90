#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module scream_p3_interface_mod

  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool,C_NULL_CHAR, c_float
  use physics_utils, only: rtype

  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

  public :: p3_init_f90
  public :: p3_main_f90
  public :: p3_finalize_f90

  real   :: test

  integer(kind=c_int) :: pcols = 32
  integer(kind=c_int) :: pver  = 72
  integer(kind=c_int) :: ncol
  integer(kind=c_int) :: nlev
  integer(kind=c_int) :: qsize = 9

  character(len=16)   :: micro_p3_tableversion = "4"
  character(len=100)  :: micro_p3_lookup_dir = "./data"
  real(kind=c_real) :: cpair  !=    1004.64000000000
  real(kind=c_real) :: rair   !=    287.042311365049
  real(kind=c_real) :: rh2o   !=    461.504639820160
  real(kind=c_real) :: rhoh2o !=    1000.00000000000
  real(kind=c_real) :: mwh2o  !=    18.0160000000000
  real(kind=c_real) :: mwdry  !=    28.9660000000000
  real(kind=c_real) :: gravit !=    9.80616000000000
  real(kind=c_real) :: latvap !=    2501000.00000000
  real(kind=c_real) :: latice !=    333700.000000000
  real(kind=c_real) :: cpliq  !=    4188.00000000000
  real(kind=c_real) :: tmelt  !=    273.150000000000
  real(kind=c_real) :: pi     !=    3.14159265358979

contains

  !====================================================================!
  subroutine p3_init_f90 () bind(c)
    use micro_p3,       only: p3_init
    use micro_p3_utils, only: micro_p3_utils_init

    call p3_init(micro_p3_lookup_dir,micro_p3_tableversion)

    test = 0.0
    print '(a15,f16.8,4i8)', 'P3 init = ', test, pcols, pver, ncol, nlev

  end subroutine p3_init_f90
  !====================================================================!
  subroutine p3_standalone_init_f90 (q,T,zi,pmid,pdel,ast,naai,ncnuc) bind(c)
    use micro_p3,       only: p3_init
    use micro_p3_utils, only: micro_p3_utils_init

    real(kind=c_real), intent(inout) :: q(pcols,pver,qsize)  ! State array  kg/kg Pa
    real(kind=c_real), intent(inout) :: T(pcols,pver)        !
    real(kind=c_real), intent(inout) :: zi(pcols,pver+1)     !
    real(kind=c_real), intent(inout) :: pmid(pcols,pver)     !
    real(kind=c_real), intent(inout) :: pdel(pcols,pver)     !
    real(kind=c_real), intent(inout) :: ast(pcols,pver)      !
    real(kind=c_real), intent(inout) :: naai(pcols,pver)     ! ice nucleation number
    real(kind=c_real), intent(inout) :: ncnuc(pcols,pver)    ! liquid activation number tendency

    character(len=100) :: case_title

    integer(kind=c_int) :: i, k
    logical(kind=c_bool) :: masterproc

    ! READ inputs from SCM for p3-stand-alone:
    open(unit=981,file='./data/p3_universal_constants.inp',status='old',action='read')
    read(981,'(A)') case_title
    read(981,'(2I8)') ncol, nlev
    if (ncol.gt.pcols.or.nlev.gt.pver) then
       print *, 'ERROR (P3-Init): inconsistency between array dimensions'
       close(981)
       return
    end if
    read(981,'(12E16.8)') cpair,rair,rh2o,rhoh2o,mwh2o,mwdry,gravit,latvap,latice,cpliq,tmelt,pi
    do i = 1,ncol
      do k = 1,nlev
        read(981,'(16E16.8)') ast(i,k), naai(i,k), ncnuc(i,k), pmid(i,k), zi(i,k), T(i,k), &
                         q(i,k,1), q(i,k,2), q(i,k,3), q(i,k,4), q(i,k,5), q(i,k,6), &
                         q(i,k,7), q(i,k,8), q(i,k,9), pdel(i,k)
      end do
      read(981,'(E16.8)') zi(i,nlev+1)
    end do
    close(981)

    masterproc = .false.
    call micro_p3_utils_init(cpair,rair,rh2o,rhoh2o,mwh2o,mwdry,gravit,latvap,latice, &
             cpliq,tmelt,pi,0,masterproc)
    print *, 'P3-Standalone-Init Finished'
  end subroutine p3_standalone_init_f90
  !====================================================================!
  subroutine p3_main_f90 (dtime,zi,pmid,pdel,ast,naai,ncnuc,q,FQ,T) bind(c)
    use micro_p3,       only: p3_main

!    real, intent(in) :: q(pcols,pver,9) ! Tracer mass concentrations from SCREAM      kg/kg
    real(kind=c_real), intent(in)    :: dtime ! Timestep
    real(kind=c_real), intent(inout) :: q(pcols,pver,qsize) ! Tracer mass concentrations from SCREAM kg/kg
    real(kind=c_real), intent(inout) :: FQ(pcols,4,pver)    ! Tracer mass tendency for physics
    real(kind=c_real), intent(inout) :: T(pcols,pver)       ! temperature
    real(kind=c_real), intent(in)    :: zi(pcols,pver+1)    ! vertical level interfaces
    real(kind=c_real), intent(in)    :: pmid(pcols,pver)    ! pressure mid-levels
    real(kind=c_real), intent(in)    :: pdel(pcols,pver)    ! pressure thickness
    real(kind=c_real), intent(in)    :: ast(pcols,pver)     ! cloud fraction
    real(kind=c_real), intent(in)    :: naai(pcols,pver)    ! ice nucleation number
    real(kind=c_real), intent(in)    :: ncnuc(pcols,pver)   ! liquid activation number tendency
    !INTERNAL VARIABLES
    real(kind=c_real) :: th(pcols,pver)         !potential temperature  K
    real(kind=c_real) :: dzq(pcols,pver)        !geometric layer thickness              m
    real(kind=c_real) :: cldliq(pcols,pver)     !cloud liquid water mixing ratio        kg/kg
    real(kind=c_real) :: numliq(pcols,pver)     !cloud liquid water drop concentraiton  #/kg
    real(kind=c_real) :: rain(pcols,pver)       !rain water mixing ratio                kg/kg
    real(kind=c_real) :: numrain(pcols,pver)    !rain water number concentration        #/kg
    real(kind=c_real) :: qv(pcols,pver)         !water vapor mixing ratio               kg/kg
    real(kind=c_real) :: ice(pcols,pver)        !total ice water mixing ratio           kg/kg
    real(kind=c_real) :: qirim(pcols,pver)      !rime ice mixing ratio                  kg/kg
    real(kind=c_real) :: numice(pcols,pver)     !total ice crystal number concentration #/kg
    real(kind=c_real) :: rimvol(pcols,pver)     !rime volume mixing ratio               m3/kg
    ! real(kind=c_real) :: temp(pcols,pver)       !potential temperature                  K
    ! real(kind=c_real) :: rim(pcols,pver)        !rime mixing ratio                      kg/kg
    real(kind=c_real) :: prt_liq(pcols)         !precipitation rate, liquid             m s-1
    real(kind=c_real) :: prt_sol(pcols)         !precipitation rate, solid              m s-1
    real(kind=c_real) :: diag_ze(pcols,pver)    !equivalent reflectivity                dBZ
    ! real(kind=c_real) :: diag_effc(pcols,pver)  !effective radius, cloud                m
    ! real(kind=c_real) :: diag_effi(pcols,pver)  !effective radius, ice                  m
    real(kind=c_real) :: diag_vmi(pcols,pver)   !mass-weighted fall speed of ice        m s-1
    real(kind=c_real) :: diag_di(pcols,pver)    !mean diameter of ice                   m
    real(kind=c_real) :: diag_rhoi(pcols,pver)  !bulk density of ice                    kg m-1
    real(kind=c_real) :: rflx(pcols,pver+1)     !grid-box average rain flux (kg m^-2s^-1) pverp
    real(kind=c_real) :: sflx(pcols,pver+1)     !grid-box average ice/snow flux (kg m^-2s^-1) pverp
    real(kind=c_real) :: exner(pcols,pver)      !exner formula for converting between potential and normal temp
    real(kind=c_real) :: rcldm(pcols,pver)      !rain cloud fraction
    real(kind=c_real) :: lcldm(pcols,pver)      !liquid cloud fraction
    real(kind=c_real) :: icldm(pcols,pver)      !ice cloud fraction
    real(kind=c_real) :: tend_out(pcols,pver,49) !microphysical tendencies
    real(kind=c_real) :: rel(pcols,pver)        !liq. effective drop radius (microns)
    real(kind=c_real) :: rei(pcols,pver)        !ice effective drop radius (microns)
    real(kind=c_real) :: cmeiout(pcols,pver)    !deposition/sublimation rate of cloud ice
    real(kind=c_real) :: prain(pcols,pver)      !total precip
    real(kind=c_real) :: nevapr(pcols,pver)     !evap. of total precip
    real(kind=c_real) :: prer_evap(pcols,pver)  !rain evaporation
    real(kind=c_real) :: pratot(pcols,pver)     !accretion of cloud by rain
    real(kind=c_real) :: prctot(pcols,pver)     !autoconversion of cloud by rain
    real(kind=c_real) :: mu(pcols,pver)         !Size distribution shape parameter for radiation
    real(kind=c_real) :: lambdac(pcols,pver)    !Size distribution slope parameter for radiation
    real(kind=c_real) :: dei(pcols,pver)        !Diameter for ice
    real(kind=c_real) :: liq_ice_exchange(pcols,pver) ! sum of liq-ice phase change tendenices
    real(kind=c_real) :: vap_liq_exchange(pcols,pver) ! sum of vap-liq phase change tendenices
    real(kind=c_real) :: vap_ice_exchange(pcols,pver) ! sum of vap-ice phase change tendenices
    real(kind=c_real) :: qc_relvar(pcols,pver)        ! 1/(var(qc)/mean(qc)**2) for P3 subgrid qc.
    real(kind=c_real) :: inv_cp

    real(kind=c_real) :: col_location(pcols,3)

    ! For rrtmg optics. specified distribution.
    real(kind=c_real), parameter :: dcon   = 25.e-6_rtype      ! Convective size distribution effective radius (um)
    real(kind=c_real), parameter :: mucon  = 5.3_rtype         ! Convective size distribution shape parameter
    real(kind=c_real), parameter :: deicon = 50._rtype         ! Convective ice effective diameter (um)

    integer(kind=c_int) :: icol, k
    integer :: i

    integer(kind=c_int) :: it, its, ite, kts, kte
    logical(kind=c_bool) :: log_predictNc = .true.
    character(len=16) :: precip_frac_method = 'max_overlap'  ! AaronDonahue, Hard-coded for now, should be fixed in the future

    real(kind=c_real) :: qtest

    inv_cp = 1.0_rtype/cpair

    qtest = sum(q)

    ! WHAT DOES P3 NEED FROM THE OUTSIDE WORLD?
    ! Q                      tracer concentrations
    ! pres                   vertical pressure profile
    ! T                      temperature profile
    ! zi                     vertical height of layer interfaces.  Note this could be backed out from pres and rho using the hydrostatic approximation
    ! pdel                   pressure layer thickness, again can be gotten from pres
    ! lcdlm, icldm, rcldm    cloud fractions
    ! exner                  exner expression.  Can be backed out from pressure.
    ! ncnuc and naai         activation arrays.


    ! MAKE LOCAL COPIES OF VARS MODIFIED BY P3
    !==============
    !local copies are needed because state is passed into this routine as intent=in
    !while P3 seeks to modify state variables in-place. Also, we need a copy of
    !old values in order to back out ptend values later. Traditionally, a local copy
    !is created by copying the whole state. It is much cheaper to just copy the
    !variables we need.

    its     = 1
    ite     = ncol
    kts     = 1
    kte     = nlev

    do i = its,ite
      do k = kts,kte
        qv(:,k)      = q(i,k,1) !1.0e-4_rtype!state%q(:,:,1)
        cldliq(:,k)  = q(i,k,2) !1.0e-6_rtype!state%q(:,:,ixcldliq)
        ice(:,k)     = q(i,k,3) !1.0e-7_rtype!state%q(:,:,ixcldice)
        numliq(:,k)  = q(i,k,4) !1.0e6_rtype!state%q(:,:,ixnumliq)
        numice(:,k)  = q(i,k,5) !1.0e5_rtype!state%q(:,:,ixnumice)
        rain(:,k)    = q(i,k,6) !1.0e-5_rtype!state%q(:,:,ixrain)
        numrain(:,k) = q(i,k,7) !1.0e5_rtype!state%q(:,:,ixnumrain)
        qirim(:,k)   = q(i,k,8) !1.0e-8_rtype!state%q(:,:,ixcldrim) !Aaron, changed ixqirim to ixcldrim to match Kai's code
        rimvol(:,k)  = q(i,k,9) !1.0e4_rtype!state%q(:,:,ixrimvol)
      end do
      col_location(i,:) = real(i)
    end do

!    do k = kte,kts,-1
!      pres(:,k)    = 1e3_rtype - (1e3_rtype-0.1)/real(pver)!state%pmid(:,:)
!    end do
    ! COMPUTE GEOMETRIC THICKNESS OF GRID
    !==============
    exner(:ncol,:nlev) = 1._rtype/((pmid(:ncol,:nlev)*1.e-5_rtype)**(rair*inv_cp))
    do icol = 1,ncol
       do k = 1,nlev
! Note: dzq is calculated in the opposite direction that pdel is calculated,
! thus when considering any dp/dz calculation we must also change the sign.
          dzq(icol,k) = zi(icol,k) - zi(icol,k+1) !100.0_rtype   !state%zi(icol,k) - state%zi(icol,k+1)
          th(icol,k)  = t(icol,k)*exner(icol,k) !/(state%pmid(icol,k)*1.e-5)**(rd*inv_cp)
!          pdel(icol,k)  = (1e3_rtype-0.1)/real(pver) ! should be changed to come from model state.
       end do
    end do
    ! Initialize the radiation dependent variables.
    mu      = mucon
    lambdac = (mucon + 1._rtype)/dcon
    dei     = deicon
    ! Determine the cloud fraction and precip cover
    call get_cloud_fraction(its,ite,kts,kte,ast(its:ite,kts:kte),cldliq(its:ite,kts:kte), &
            rain(its:ite,kts:kte),ice(its:ite,kts:kte),precip_frac_method, &
            icldm(its:ite,kts:kte),lcldm(its:ite,kts:kte),rcldm(its:ite,kts:kte))
    !    icldm(:,:) = 1.0_rtype
    !    lcldm(:,:) = 1.0_rtype
    !    rcldm(:,:) = 1.0_rtype

    ! Hack qc_relvar (should get more thoughtful value later):
    qc_relvar(:,:) = 1.0_rtype

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
         pmid(its:ite,kts:kte),       & ! IN     pressure at cell midpoints       Pa
         dzq(its:ite,kts:kte),        & ! IN     vertical grid spacing            m
         ncnuc(its:ite,kts:kte),      & ! IN ccn activation number tendency kg-1 s-1
         naai(its:ite,kts:kte),       & ! IN activated ice nuclei concentration kg-1
         qc_relvar(its:ite,kts:kte),  & ! IN 1/(var(qc)/mean(qc)**2) used in P3.
         it,                          & ! IN     time step counter NOTE: starts at 1 for first time step
         prt_liq(its:ite),            & ! OUT    surface liquid precip rate       m s-1
         prt_sol(its:ite),            & ! OUT    surface frozen precip rate       m s-1
         its,                         & ! IN     horizontal index lower bound     -
         ite,                         & ! IN     horizontal index upper bound     -
         kts,                         & ! IN     vertical index lower bound       -
         kte,                         & ! IN     vertical index upper bound       -
         rel(its:ite,kts:kte),        & ! OUT    effective radius, cloud          m
         rei(its:ite,kts:kte),        & ! OUT    effective radius, ice            m
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
         tend_out(its:ite,kts:kte,:), & ! OUT p3 microphysics tendencies
         mu(its:ite,kts:kte),         & ! OUT Size distribution shape parameter for radiation
         lambdac(its:ite,kts:kte),     & ! OUT Size distribution slope parameter for radiation
         liq_ice_exchange(its:ite,kts:kte),& ! OUT sum of liq-ice phase change tendenices
         vap_liq_exchange(its:ite,kts:kte),& ! OUT sun of vap-liq phase change tendencies
         vap_ice_exchange(its:ite,kts:kte),& ! OUT sum of vap-ice phase change tendencies
         col_location(its:ite,3)           & ! IN location of columns
         )
    do i = its,ite
      do k = kts,kte
        FQ(i,k,1) = FQ(i,k,1) + (qv(i,k)     - q(i,k,1))/dtime
        FQ(i,k,2) = FQ(i,k,2) + (cldliq(i,k) - q(i,k,2))/dtime
        FQ(i,k,3) = FQ(i,k,3) + (ice(i,k)    - q(i,k,3))/dtime
        FQ(i,k,4) = FQ(i,k,4) + (rain(i,k)   - q(i,k,6))/dtime
        q(i,k,1) = qv(i,k)
        q(i,k,2) = cldliq(i,k)
        q(i,k,3) = ice(i,k)
        q(i,k,4) = numliq(i,k)
        q(i,k,5) = numice(i,k)
        q(i,k,6) = rain(i,k)
        q(i,k,7) = numrain(i,k)
        q(i,k,8) = qirim(i,k)
        q(i,k,9) = rimvol(i,k)
      end do
    end do

    print '(a15,f16.8,4e16.8)', 'P3 run = ', test, qtest, sum(q), sum(qv), sum(FQ(:,:,:))
  end subroutine p3_main_f90
  !====================================================================!
  subroutine p3_finalize_f90 () bind(c)

    print *, "P3 finalization..."

  end subroutine p3_finalize_f90
  !================================================================================================
    subroutine get_cloud_fraction(its,ite,kts,kte,ast,qc,qr,qitot,method, &
                  icldm,lcldm,rcldm)

       use micro_p3_utils, only: mincld, qsmall

       integer,intent(in)                                 :: its,ite,kts,kte
       real(rtype),dimension(its:ite,kts:kte),intent(in)  :: ast, qc, qr, qitot
       character(len=16),intent(in)                       :: method
       real(rtype),dimension(its:ite,kts:kte),intent(out) :: icldm, lcldm, rcldm
       real(rtype),dimension(its:ite,kts:kte)             :: cldm

       integer  :: i,k
       integer  :: ktop, kbot, kdir

       ktop = kts        !k of top level
       kbot = kte        !k of bottom level
       kdir = -1         !(k: 1=top, nk=bottom)

       cldm(:,:)  = mincld
       icldm(:,:) = mincld
       lcldm(:,:) = mincld
       do k = kbot,ktop,kdir
          do i=its,ite
             cldm(i,k)  = max(ast(i,k), mincld)
             icldm(i,k) = max(ast(i,k), mincld)
             lcldm(i,k) = max(ast(i,k), mincld)
          end do
       end do

       DO k = ktop,kbot,-kdir  !AaronDonahue TODO: Check to make sure this is correct.  Are we going the correct direction?
          DO i=its,ite
       !!
       !! precipitation fraction
       !!
          rcldm(i,k) = cldm(i,k)
          IF (trim(method) == 'in_cloud') THEN
             IF (k /= ktop) THEN
                IF (qc(i,k) .lt. qsmall .and. qitot(i,k) .lt. qsmall) THEN
                   rcldm(i,k) = rcldm(i,k+kdir)
                END IF
             END IF
          ELSE IF (trim(method) == 'max_overlap') THEN
          ! calculate precip fraction based on maximum overlap assumption

          ! IF rain or snow mix ratios are smaller than threshold,
          ! then leave rcldm as cloud fraction at current level
             IF (k /= ktop) THEN
                IF (qr(i,k+kdir) .ge. qsmall .or. qitot(i,k+kdir) .ge. qsmall) THEN
                   rcldm(i,k) = max(cldm(i,k+kdir),rcldm(i,k))
                END IF
             END IF
          END IF
          END DO ! i
       END DO    ! k

       return
    end subroutine get_cloud_fraction

  !================================================================================================
  !====================================================================!

end module scream_p3_interface_mod
