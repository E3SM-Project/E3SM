!__________________________________________________________________________________________
! This module contains the Predicted Particle Property (P3) bulk microphysics scheme.      !
!                                                                                          !
! This code was originally written by H. Morrison,  MMM Division, NCAR (Dec 2012).         !
! Modifications were made by J. Milbrandt, RPN, Environment Canada (July 2014).            !
! Peter Caldwell (caldwell19@llnl.gov) further modified this code to remove multiple       !
! ice categories and to clean up/simplify the code for conversion to C++ (9/11/18)         !
!                                                                                          !
! Three configurations of the P3 scheme are currently available:                           !
!  1) specified droplet number (i.e. 1-moment cloud water), 1 ice category                 !
!  2) predicted droplet number (i.e. 2-moment cloud water), 1 ice category                 !
!                                                                                          !
!  The  2-moment cloud version is based on a specified aerosol distribution and            !
!  does not include a subgrid-scale vertical velocity for droplet activation. Hence,       !
!  this version should only be used for high-resolution simulations that resolve           !
!  vertical motion driving droplet activation.                                             !
!                                                                                          !
! For details see: Morrison and Milbrandt (2015) [J. Atmos. Sci., 72, 287-311]             !
!                  Milbrandt and Morrison (2016) [J. Atmos. Sci., 73, 975-995]             !
!                                                                                          !
! For questions or bug reports, please contact:                                            !
!    Hugh Morrison   (morrison@ucar.edu), or                                               !
!    Jason Milbrandt (jason.milbrandt@canada.ca)                                           !
!__________________________________________________________________________________________!
!                                                                                          !
! Version:       2.8.2.4 + Peter/Aaron's fixes                                             !
! Last updated:  2018-02-04                                                                !
!__________________________________________________________________________________________!
! Comments from Aaron Donahue:                                                             !
! 1) Need to change the dz coordinate system in sedimentation to be consistent             !
! with E3SM's pressure based coordinate system, i.e. dp.                                   !
! 2) Move all physical constants into a micro_p3_util module and match them to             !
! universal constants in E3SM for consistency.                                            !
! 3) Need to include extra in/out values which correspond with microphysics PBUF           !
! variables and outputs expected in E3SM.                                                  !
!__________________________________________________________________________________________!

! Include bit-for-bit math macros.
#include "bfb_math.inc"

module micro_p3

   ! get real kind from utils
   use physics_utils, only: rtype,rtype8,btype

   ! physical and mathematical constants
   use micro_p3_utils, only: rhosur,rhosui,ar,br,f1r,f2r,rhow,kr,kc,aimm,mi0,nccnst,  &
       eci,eri,bcn,cpw,cons1,cons3,cons4,cons5,cons6,cons7,         &
       inv_rhow,inv_dropmass,qsmall,nsmall,cp,g,rd,rv,ep_2,inv_cp,   &
       thrd,sxth,piov6,rho_rimeMin,     &
       rho_rimeMax,inv_rho_rimeMax,max_total_Ni,dbrk,nmltratio,clbfact_sub,  &
       clbfact_dep,iparam, isize, densize, rimsize, rcollsize, tabsize, colltabsize, &
       get_latent_heat, zerodegc, pi=>pi_e3sm, dnu, &
       rainfrze, icenuct, homogfrze, iulog=>iulog_e3sm, &
       masterproc=>masterproc_e3sm, calculate_incloud_mixingratios, mu_r_constant, &
       lookup_table_1a_dum1_c, &
       p3_QcAutoCon_Expon, p3_QcAccret_Expon

  ! Bit-for-bit math functions.
#ifdef SCREAM_CONFIG_IS_CMAKE
  use physics_share_f2c, only: cxx_pow, cxx_sqrt, cxx_cbrt, cxx_gamma, cxx_log, &
                                 cxx_log10, cxx_exp
#endif

  implicit none
  save

  public  :: p3_init,p3_main

  ! protected items should be treated as private for everyone except tests

  real(rtype), protected, dimension(densize,rimsize,isize,tabsize) :: itab   !ice lookup table values

  !ice lookup table values for ice-rain collision/collection
  real(rtype), protected, dimension(densize,rimsize,isize,rcollsize,colltabsize) :: itabcoll

  ! lookup table values for rain shape parameter mu_r
  real(rtype), protected, dimension(150) :: mu_r_table

  ! lookup table values for rain number- and mass-weighted fallspeeds and ventilation parameters
  real(rtype), protected, dimension(300,10) :: vn_table,vm_table,revap_table

  type realptr
     real(rtype), dimension(:), pointer :: p
  end type realptr

contains

  SUBROUTINE p3_init(lookup_file_dir,version_p3)
    !------------------------------------------------------------------------------------------!
    ! This subroutine initializes all physical constants and parameters needed by the P3       !
    ! scheme, including reading in two lookup table files and creating a third.                !
    ! 'P3_INIT' be called at the first model time step, prior to first call to 'P3_MAIN'.      !
    !------------------------------------------------------------------------------------------!

    implicit none

    ! Passed arguments:
    character*(*), intent(in)    :: lookup_file_dir                !directory of the lookup tables
    character(len=16), intent(in) :: version_p3  !version number of P3 package

    if (masterproc) write(iulog,*) ''
    if (masterproc) write(iulog,*) ' P3 microphysics: v',version_p3

    call p3_init_a(lookup_file_dir, version_p3)
    call p3_init_b()

    if (masterproc) write(iulog,*) '   P3_INIT DONE.'
    if (masterproc) write(iulog,*) ''

  END SUBROUTINE p3_init

  SUBROUTINE p3_init_a(lookup_file_dir,version_p3)

    use scream_abortutils, only : endscreamrun

    ! Passed arguments:
    character*(*), intent(in)     :: lookup_file_dir       !directory of the lookup tables

    character(len=16), intent(in) :: version_p3            !version number of P3 package
    character(len=1024)           :: lookup_file_1         !lookup table, maini
    character(len=1024)           :: version_header_table_1             !version number read from header, table 1
    integer                       :: i,j,ii,jj
#ifdef SCREAM_CONFIG_IS_CMAKE
    real(rtype8)                  :: dum,dumk1,dumk2
    real(rtype8), dimension(12)   :: dumk
#else
    real(rtype)                   :: dum,dumk1,dumk2
    real(rtype), dimension(12)    :: dumk
#endif
    integer                       :: dumi
    character(len=1024)           :: dumstr

    !------------------------------------------------------------------------------------------!

    lookup_file_1 = trim(lookup_file_dir)//'/'//'p3_lookup_table_1.dat-v'//trim(version_p3)

    !------------------------------------------------------------------------------------------!
    ! read in ice microphysics table

    if (masterproc) write(iulog,*) '   P3_INIT (reading/creating look-up tables) ...'

    open(unit=10,file=lookup_file_1, status='old', action='read')

    read(10,*) dumstr, version_header_table_1
    if (trim(version_p3) /= trim(version_header_table_1)) then
       print*
       print*, '***********   WARNING in P3_INIT   *************'
       print*, ' Loading lookupTable_1: v',trim(version_header_table_1)
       print*, ' P3 is intended to use lookupTable_1: v', trim(version_p3)
       print*, '               -- ABORTING -- '
       print*, '************************************************'
       print*
       call endscreamrun()
    end if

    itab(:,:,:,:) = 0.
    itabcoll(:,:,:,:,:) = 0.
    do jj = 1,densize
       do ii = 1,rimsize
          do i = 1,isize
             read(10,*) dumi,dumi,dum,dum,dumk(1),dumk(2),           &
                  dumk(3),dumk(4),dumk(5),dumk(6),dumk(7),dumk(8),dum,                 &
                  dumk(9),dumk(10),dumk(11),dumk(12)
             itab(jj,ii,i,:) = dumk(:)
          enddo
          ! read in table for ice-rain collection
          do i = 1,isize
             do j = 1,rcollsize
                read(10,*) dumi,dumi,dum,dum,dum,dumk1,dumk2,dum
                itabcoll(jj,ii,i,j,1) = dlog10(dumk1)
                itabcoll(jj,ii,i,j,2) = dlog10(dumk2)
             enddo
          enddo
       enddo
    enddo

    ! hm add fix to prevent end-of-file error in nested runs, 3/28/14
    close(10)

    !PMC: deleted ice-ice collision lookup table here b/c only used for nCat>1.
    ! So there is no need to fill lookup values for lookup table 2.

   return

  END SUBROUTINE p3_init_a

  subroutine p3_get_tables(mu_r_user, revap_user, vn_user, vm_user)
    ! This can be called after p3_init_b.
    implicit none
    real(rtype), dimension(150), intent(out) :: mu_r_user
    real(rtype), dimension(300,10), intent(out) :: vn_user, vm_user, revap_user
    mu_r_user(:) = mu_r_table(:)
    revap_user(:,:) = revap_table(:,:)
    vn_user(:,:) = vn_table(:,:)
    vm_user(:,:) = vm_table(:,:)

   return

  end subroutine p3_get_tables

  subroutine p3_set_tables(mu_r_user, revap_user, vn_user, vm_user)
    ! This can be called instead of p3_init_b.
    implicit none
    real(rtype), dimension(150), intent(in) :: mu_r_user
    real(rtype), dimension(300,10), intent(in) :: vn_user, vm_user, revap_user
    mu_r_table(:) = mu_r_user(:)
    revap_table(:,:) = revap_user(:,:)
    vn_table(:,:) = vn_user(:,:)
    vm_table(:,:) = vm_user(:,:)

   return

  end subroutine p3_set_tables

  SUBROUTINE p3_init_b()
    implicit none
    integer                      :: i,ii,jj,kk
    real(rtype)                         :: lamr,mu_r,dm,dum1,dum2,dum3,dum4,dum5,  &
         dd,amg,vt,dia

    !------------------------------------------------------------------------------------------!

    ! Generate lookup table for rain shape parameter mu_r
    ! this is very fast so it can be generated at the start of each run
    ! make a 150x1 1D lookup table, this is done in parameter
    ! space of a scaled mean size proportional qr/Nr -- initlamr

    !write(iulog,*) '   Generating rain lookup-table ...'

    ! AaronDonahue: Switching to table ver 4 means switching to a constand mu_r,
    ! so this section is commented out.
    do i = 1,150              ! loop over lookup table values
!       initlamr = 1./((real(i)*2.)*1.e-6 + 250.e-6)
!
!       ! iterate to get mu_r
!       ! mu_r-lambda relationship is from Cao et al. (2008), eq. (7)
!
!       ! start with first guess, mu_r = 0
!
!       mu_r = 0.
!
!       do ii=1,50
!          lamr = initlamr*((mu_r+3.)*(mu_r+2.)*(mu_r+1.)/6.)**thrd
!
!          ! new estimate for mu_r based on lambda
!          ! set max lambda in formula for mu_r to 20 mm-1, so Cao et al.
!          ! formula is not extrapolated beyond Cao et al. data range
!          dum  = min(20.,lamr*1.e-3)
!          mu_r = max(0.,-0.0201*dum**2+0.902*dum-1.718)
!
!          ! if lambda is converged within 0.1%, then exit loop
!          if (ii.ge.2) then
!             if (abs((lamold-lamr)/lamr).lt.0.001) goto 111
!          end if
!
!          lamold = lamr
!
!       enddo
!
!111    continue
!
!       ! assign lookup table values
       mu_r_table(i) = mu_r_constant

    enddo

    !.......................................................................
    ! Generate lookup table for rain fallspeed and ventilation parameters
    ! the lookup table is two dimensional as a function of number-weighted mean size
    ! proportional to qr/Nr and shape parameter mu_r

    mu_r_loop: do ii = 1,10   !** change 10 to 9, since range of mu_r is 0-8  CONFIRM
       !mu_r_loop: do ii = 1,9   !** change 10 to 9, since range of mu_r is 0-8

!       mu_r = real(ii-1)  ! values of mu
       mu_r = mu_r_constant

       ! loop over number-weighted mean size
       meansize_loop: do jj = 1,300

          if (jj.le.20) then
             dm = (real(jj)*10._rtype-5._rtype)*1.e-6_rtype      ! mean size [m]
          elseif (jj.gt.20) then
             dm = (real(jj-20)*30._rtype+195._rtype)*1.e-6_rtype ! mean size [m]
          endif

          lamr = (mu_r+1._rtype)/dm

          ! do numerical integration over PSD

          dum1 = 0._rtype ! numerator,   number-weighted fallspeed
          dum2 = 0._rtype ! denominator, number-weighted fallspeed
          dum3 = 0._rtype ! numerator,   mass-weighted fallspeed
          dum4 = 0._rtype ! denominator, mass-weighted fallspeed
          dum5 = 0._rtype ! term for ventilation factor in evap
          dd   = 2._rtype

          ! loop over PSD to numerically integrate number and mass-weighted mean fallspeeds
          do kk = 1,10000

             dia = (real(kk)*dd-dd*0.5_rtype)*1.e-6_rtype  ! size bin [m]
             amg = piov6*997._rtype*dia**3           ! mass [kg]
             amg = amg*1000._rtype                   ! convert [kg] to [g]

             !get fallspeed as a function of size [m s-1]
             if (dia*1.e+6_rtype.le.134.43_rtype)      then
                vt = 4.5795e+3_rtype*amg**(2._rtype*thrd)
             elseif (dia*1.e+6_rtype.lt.1511.64_rtype) then
                vt = 4.962e+1_rtype*amg**thrd
             elseif (dia*1.e+6_rtype.lt.3477.84_rtype) then
                vt = 1.732e+1_rtype*amg**sxth
             else
                vt = 9.17_rtype
             endif

             !note: factor of 4.*mu_r is non-answer changing and only needed to
             !      prevent underflow/overflow errors, same with 3.*mu_r for dum5
             dum1 = dum1 + vt*10._rtype**(mu_r*log10(dia)+4._rtype*mu_r)*exp(-lamr*dia)*dd*1.e-6_rtype
             dum2 = dum2 + 10._rtype**(mu_r*log10(dia)+4._rtype*mu_r)*exp(-lamr*dia)*dd*1.e-6_rtype
             dum3 = dum3 + vt*10._rtype**((mu_r+3._rtype)*log10(dia)+4._rtype*mu_r)*exp(-lamr*dia)*dd*1.e-6_rtype
             dum4 = dum4 + 10._rtype**((mu_r+3._rtype)*log10(dia)+4._rtype*mu_r)*exp(-lamr*dia)*dd*1.e-6_rtype
             dum5 = dum5 + (vt*dia)**0.5*10.**((mu_r+1.)*log10(dia)+3.*mu_r)*exp(-lamr*dia)*dd*1.e-6

          enddo ! kk-loop (over PSD)

          dum2 = max(dum2, 1.e-30_rtype)  !to prevent divide-by-zero below
          dum4 = max(dum4, 1.e-30_rtype)  !to prevent divide-by-zero below
          dum5 = max(dum5, 1.e-30_rtype)  !to prevent log10-of-zero below

          vn_table(jj,ii)    = dum1/dum2
          vm_table(jj,ii)    = dum3/dum4
          revap_table(jj,ii) = 10._rtype**(log10(dum5)+(mu_r+1._rtype)*log10(lamr)-(3._rtype*mu_r))

       enddo meansize_loop

    enddo mu_r_loop

  END SUBROUTINE p3_init_b

  !==========================================================================================!

  SUBROUTINE p3_main_part1(kts, kte, kbot, ktop, kdir, log_predictNc, dt, &
       pres, pdel, dzq, ncnuc, exner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, xxlv, xxls, xlf, &
       t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, acn, qv, th, qc, nc, qr, nr, &
       qitot, nitot, qirim, birim, qc_incld, qr_incld, qitot_incld, qirim_incld, &
       nc_incld, nr_incld, nitot_incld, birim_incld, log_nucleationPossible, log_hydrometeorsPresent)

    implicit none

    ! args

    integer, intent(in) :: kts, kte, kbot, ktop, kdir
    logical(btype), intent(in) :: log_predictNc
    real(rtype), intent(in) :: dt

    real(rtype), intent(in), dimension(kts:kte) :: pres, pdel, dzq, ncnuc, exner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, xxlv, xxls, xlf

    real(rtype), intent(inout), dimension(kts:kte) :: t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, &
         acn, qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, qc_incld, qr_incld, qitot_incld, &
         qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld

    logical(btype), intent(out) :: log_nucleationPossible, log_hydrometeorsPresent

    ! locals
    integer :: k
    real(rtype) :: dum

    log_nucleationPossible = .false.
    log_hydrometeorsPresent = .false.

    k_loop_1: do k = kbot,ktop,kdir
       !calculate some time-varying atmospheric variables
       !AaronDonahue - changed "rho" to be defined on nonhydrostatic
       !assumption, consistent with pressure based coordinate system
       !             - moved latent heat calculation to above.  Latent
       !heat is determined by calling a p3_util function so that it
       !can be made consistent with E3SM definition of latent heat
       rho(k)     = pdel(k)/dzq(k)/g  ! pres(k)/(rd*t(k))
       inv_rho(k) = 1._rtype/rho(k)
       qvs(k)     = qv_sat(t(k),pres(k),0)
       qvi(k)     = qv_sat(t(k),pres(k),1)

       supi(k)    = qv(k)/qvi(k)-1._rtype

       rhofacr(k) = bfb_pow(rhosur*inv_rho(k), 0.54_rtype)
       rhofaci(k) = bfb_pow(rhosui*inv_rho(k), 0.54_rtype)
       dum        = 1.496e-6_rtype * bfb_pow(t(k), 1.5_rtype) / (t(k)+120._rtype)  ! this is mu
       acn(k)     = g*rhow/(18._rtype*dum)  ! 'a' parameter for droplet fallspeed (Stokes' law)

       if ((t(k).lt.zerodegc .and. supi(k).ge.-0.05_rtype)) log_nucleationPossible = .true.

       if (qc(k).lt.qsmall) then
      !--- apply mass clipping if mass is sufficiently small
      !    (implying all mass is expected to evaporate/sublimate in one time step)
          qv(k) = qv(k) + qc(k)
          th(k) = th(k) - exner(k)*qc(k)*xxlv(k)*inv_cp
          qc(k) = 0._rtype
          nc(k) = 0._rtype
       else
          log_hydrometeorsPresent = .true.    ! updated further down
      !--- Apply droplet activation here (before other microphysical processes) for consistency with qc increase by saturation
      !    adjustment already applied in macrophysics. If prescribed drop number is used, this is also a good place to
      !    prescribe that value
          if (.not.(log_predictNc)) then
            nc(k) = nccnst*inv_rho(k)
          else
            nc(k) = max(nc(k) + ncnuc(k) * dt,0.0_rtype)
          endif
       endif

       if (qr(k).lt.qsmall) then
          qv(k) = qv(k) + qr(k)
          th(k) = th(k) - exner(k)*qr(k)*xxlv(k)*inv_cp
          qr(k) = 0._rtype
          nr(k) = 0._rtype
       else
          log_hydrometeorsPresent = .true.    ! updated further down
       endif

       if (qitot(k).lt.qsmall .or. (qitot(k).lt.1.e-8_rtype .and.             &
            supi(k).lt.-0.1_rtype)) then
          qv(k) = qv(k) + qitot(k)
          th(k) = th(k) - exner(k)*qitot(k)*xxls(k)*inv_cp
          qitot(k) = 0._rtype
          nitot(k) = 0._rtype
          qirim(k) = 0._rtype
          birim(k) = 0._rtype
       else
          log_hydrometeorsPresent = .true.    ! final update
       endif

       if (qitot(k).ge.qsmall .and. qitot(k).lt.1.e-8_rtype .and.             &
            t(k).ge.zerodegc) then
          qr(k) = qr(k) + qitot(k)
          th(k) = th(k) - exner(k)*qitot(k)*xlf(k)*inv_cp
          qitot(k) = 0._rtype
          nitot(k) = 0._rtype
          qirim(k) = 0._rtype
          birim(k) = 0._rtype
       endif

       t(k) = th(k) * inv_exner(k)

       call calculate_incloud_mixingratios(qc(k),qr(k),qitot(k),qirim(k),nc(k),nr(k),nitot(k),birim(k), &
            inv_lcldm(k),inv_icldm(k),inv_rcldm(k), &
            qc_incld(k),qr_incld(k),qitot_incld(k),qirim_incld(k),nc_incld(k),nr_incld(k),nitot_incld(k),birim_incld(k))

    enddo k_loop_1

  END SUBROUTINE p3_main_part1

  SUBROUTINE p3_main_part2(kts, kte, kbot, ktop, kdir, log_predictNc, dt, odt, &
       pres, pdel, dzq, ncnuc, exner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, naai, qc_relvar, icldm, lcldm, rcldm,&
       t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, acn, qv, th, qc, nc, qr, nr, qitot, nitot, &
       qirim, birim, xxlv, xxls, xlf, qc_incld, qr_incld, qitot_incld, qirim_incld, nc_incld, nr_incld, &
       nitot_incld, birim_incld, mu_c, nu, lamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, cmeiout, prain, &
       nevapr, prer_evap, vap_liq_exchange, vap_ice_exchange, liq_ice_exchange, pratot, &
       prctot, p3_tend_out, log_hydrometeorsPresent)

    implicit none

    ! args

    integer, intent(in) :: kts, kte, kbot, ktop, kdir
    logical(btype), intent(in) :: log_predictNc
    real(rtype), intent(in) :: dt, odt

    real(rtype), intent(in), dimension(kts:kte) :: pres, pdel, dzq, ncnuc, exner, inv_exner, inv_lcldm, inv_icldm,   &
         inv_rcldm, naai, qc_relvar, icldm, lcldm, rcldm

    real(rtype), intent(inout), dimension(kts:kte) :: t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, acn,        &
         qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls, xlf, qc_incld, qr_incld,                    &
         qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld, mu_c, nu, lamc, cdist, cdist1,      &
         cdistr, mu_r, lamr, logn0r, cmeiout, prain, nevapr, prer_evap, vap_liq_exchange,                            &
         vap_ice_exchange, liq_ice_exchange, pratot, prctot

    real(rtype), intent(inout), dimension(kts:kte,49) :: p3_tend_out ! micro physics tendencies

    logical(btype), intent(out) :: log_hydrometeorsPresent

    ! -------- locals ------- !

    ! liquid-phase microphysical process rates:
    !  (all Q process rates in kg kg-1 s-1)
    !  (all N process rates in # kg-1)

    real(rtype) :: qcacc   ! cloud droplet accretion by rain
    real(rtype) :: qcaut   ! cloud droplet autoconversion to rain
    real(rtype) :: ncacc   ! change in cloud droplet number from accretion by rain
    real(rtype) :: ncautc  ! change in cloud droplet number from autoconversion
    real(rtype) :: ncslf   ! change in cloud droplet number from self-collection  (Not in paper?)
    real(rtype) :: nrslf   ! change in rain number from self-collection  (Not in paper?)
    real(rtype) :: qrevp   ! rain evaporation
    real(rtype) :: nrevp   ! change in rain number from evaporation
    real(rtype) :: ncautr  ! change in rain number from autoconversion of cloud water
! Is is assumed that macrophysics handles condensation/evaporation of qc and
! that there is no condensation of rain. Thus qccon, qrcon and qcevp have
! been removed from the original P3-WRF for P3-SCREAM.
!    real(rtype) :: qrcon   ! rain condensation   (Not in paper?)
!    real(rtype) :: qccon   ! cloud droplet condensation
!    real(rtype) :: qcevp   ! cloud droplet evaporation

    ! ice-phase microphysical process rates:
    !  (all Q process rates in kg kg-1 s-1)
    !  (all N process rates in # kg-1)

    real(rtype) :: qccol     ! collection of cloud water by ice
    real(rtype) :: qwgrth    ! wet growth rate
    real(rtype) :: qidep     ! vapor deposition
    real(rtype) :: qrcol     ! collection rain mass by ice
    real(rtype) :: qinuc     ! deposition/condensation freezing nuc
    real(rtype) :: nccol     ! change in cloud droplet number from collection by ice
    real(rtype) :: nrcol     ! change in rain number from collection by ice
    real(rtype) :: ninuc     ! change in ice number from deposition/cond-freezing nucleation
    real(rtype) :: qisub     ! sublimation of ice
    real(rtype) :: qimlt     ! melting of ice
    real(rtype) :: nimlt     ! melting of ice
    real(rtype) :: nisub     ! change in ice number from sublimation
    real(rtype) :: nislf     ! change in ice number from collection within a category (Not in paper?)
    real(rtype) :: qcheti    ! immersion freezing droplets
    real(rtype) :: qrheti    ! immersion freezing rain
    real(rtype) :: ncheti    ! immersion freezing droplets
    real(rtype) :: nrheti    ! immersion freezing rain
    real(rtype) :: nrshdr    ! source for rain number from collision of rain/ice above freezing and shedding
    real(rtype) :: qcshd     ! source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding
    real(rtype) :: rhorime_c ! density of rime (from cloud)
    real(rtype) :: ncshdc    ! source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)
    real(rtype) :: qiberg    ! Bergeron process

    real(rtype)    :: f1pr02   ! mass-weighted fallspeed              See lines  731 -  808  ums
    real(rtype)    :: f1pr03   ! ice collection within a category     See lines  809 -  928  nagg
    real(rtype)    :: f1pr04   ! collection of cloud water by ice     See lines  929 - 1009  nrwat
    real(rtype)    :: f1pr05   ! melting                              See lines 1212 - 1279  vdep
    real(rtype)    :: f1pr07   ! collection of rain number by ice     See lines 1010 - 1209  nrrain
    real(rtype)    :: f1pr08   ! collection of rain mass by ice       See lines 1010 - 1209  qrrain
    real(rtype)    :: f1pr09   ! minimum ice number (lambda limiter)  See lines  704 -  705  nlarge
    real(rtype)    :: f1pr10   ! maximum ice number (lambda limiter)  See lines  704 -  705  nsmall
    real(rtype)    :: f1pr14   ! melting (ventilation term)           See lines 1212 - 1279  vdep1

    real(rtype)    :: mu,dv,sc,dqsdt,ab,kap,epsr,epsc,epsi,epsi_tot, &
         dum1,dum3,dum4,dum5,dum6,dqsidt,abi,rhop,vtrmi1,eii

    integer :: dumi,k,dumj,dumii,dumjj,dumzz

    logical(btype) :: log_exitlevel, log_wetgrowth

   rhorime_c = 400._rtype
   log_hydrometeorsPresent = .false.

   !------------------------------------------------------------------------------------------!
   !   main k-loop (for processes):
   k_loop_main: do k = kbot,ktop,kdir

      ! if relatively dry and no hydrometeors at this level, skip to end of k-loop (i.e. skip this level)
      log_exitlevel = .true.
      if (qc(k).ge.qsmall .or. qr(k).ge.qsmall) log_exitlevel = .false.

      if (qitot(k).ge.qsmall) log_exitlevel = .false.
      !enddo
      if (log_exitlevel .and.                                                           &
         (t(k).lt.zerodegc .and. supi(k).lt.-0.05_rtype)) goto 555   !i.e. skip all process rates

      ! All microphysics tendencies will be computed as IN-CLOUD, they will be mapped back to cell-average later.

      ! initialize warm-phase process rates
      qcacc   = 0._rtype;     qrevp   = 0._rtype;     qcaut   = 0._rtype;
      ncacc   = 0._rtype;     ncslf   = 0._rtype;
      ncautc  = 0._rtype;     nrslf   = 0._rtype;
      nrevp   = 0._rtype;     ncautr  = 0._rtype

      ! initialize ice-phase  process rates
      qisub   = 0._rtype;     nrshdr  = 0._rtype
      qcheti  = 0._rtype;     qrcol   = 0._rtype;     qcshd   = 0._rtype
      qimlt   = 0._rtype;     qccol   = 0._rtype
      qrheti  = 0._rtype;     qinuc   = 0._rtype;     nimlt   = 0._rtype
      nccol   = 0._rtype;     ncshdc  = 0._rtype
      ncheti  = 0._rtype;     nrcol   = 0._rtype;     nislf   = 0._rtype
      ninuc   = 0._rtype;     qidep   = 0._rtype;     qiberg  = 0._rtype
      nrheti  = 0._rtype;     nisub   = 0._rtype;     qwgrth  = 0._rtype

      log_wetgrowth = .false.

      ! skip micro process calculations except nucleation/acvtivation if there no hydrometeors are present
      log_exitlevel = .true.
      if (qc_incld(k).ge.qsmall .or. qr_incld(k).ge.qsmall) log_exitlevel = .false.
      if (qitot_incld(k).ge.qsmall) log_exitlevel=.false.
      if (log_exitlevel) goto 444   !i.e. skip to nucleation

      !time/space varying physical variables
      call get_time_space_phys_variables( &
           t(k),pres(k),rho(k),xxlv(k),xxls(k),qvs(k),qvi(k), &
           mu,dv,sc,dqsdt,dqsidt,ab,abi,kap,eii)

      call get_cloud_dsd2(qc_incld(k),nc_incld(k),mu_c(k),rho(k),nu(k),dnu,lamc(k),     &
           cdist(k),cdist1(k),lcldm(k))
      nc(k) = nc_incld(k)*lcldm(k)

      call get_rain_dsd2(qr_incld(k),nr_incld(k),mu_r(k),lamr(k),   &
           cdistr(k),logn0r(k),rcldm(k))
      nr(k) = nr_incld(k)*rcldm(k)

      ! initialize inverse supersaturation relaxation timescale for combined ice categories
      epsi_tot = 0._rtype

      call impose_max_total_Ni(nitot_incld(k),max_total_Ni,inv_rho(k))

      if (qitot_incld(k).ge.qsmall) then

         !impose lower limits to prevent taking log of # < 0
         nitot_incld(k) = max(nitot_incld(k),nsmall)
         nr_incld(k)    = max(nr_incld(k),nsmall)

         call calc_bulkRhoRime(qitot_incld(k),qirim_incld(k),birim_incld(k),rhop)

         ! if (.not. tripleMoment_on) zitot(k) = diag_mom6(qitot_incld(k),nitot_incld(k),rho(k))
         call find_lookupTable_indices_1a(dumi,dumjj,dumii,dumzz,dum1,dum4,          &
              dum5,dum6,isize,rimsize,densize,                &
              qitot_incld(k),nitot_incld(k),qirim_incld(k),      &
              rhop)
         !qirim_incld(k),zitot(k),rhop)
         call find_lookupTable_indices_1b(dumj,dum3,rcollsize,qr_incld(k),nr_incld(k))

         ! call to lookup table interpolation subroutines to get process rates
         call access_lookup_table(dumjj,dumii,dumi, 2,dum1,dum4,dum5,f1pr02)
         call access_lookup_table(dumjj,dumii,dumi, 3,dum1,dum4,dum5,f1pr03)
         call access_lookup_table(dumjj,dumii,dumi, 4,dum1,dum4,dum5,f1pr04)
         call access_lookup_table(dumjj,dumii,dumi, 5,dum1,dum4,dum5,f1pr05)
         call access_lookup_table(dumjj,dumii,dumi, 7,dum1,dum4,dum5,f1pr09)
         call access_lookup_table(dumjj,dumii,dumi, 8,dum1,dum4,dum5,f1pr10)
         call access_lookup_table(dumjj,dumii,dumi,10,dum1,dum4,dum5,f1pr14)

         ! ice-rain collection processes
         if (qr_incld(k).ge.qsmall) then
            call access_lookup_table_coll(dumjj,dumii,dumj,dumi,1,dum1,dum3,dum4,dum5,f1pr07)
            call access_lookup_table_coll(dumjj,dumii,dumj,dumi,2,dum1,dum3,dum4,dum5,f1pr08)
         else
            f1pr07 = 0._rtype
            f1pr08 = 0._rtype
         endif

         ! adjust Ni if needed to make sure mean size is in bounds (i.e. apply lambda limiters)
         ! note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
         nitot_incld(k) = min(nitot_incld(k),f1pr09*nitot_incld(k))
         nitot_incld(k) = max(nitot_incld(k),f1pr10*nitot_incld(k))

      endif   ! qitot > qsmall

      !----------------------------------------------------------------------
      ! Begin calculations of microphysical processes

      !......................................................................
      ! ice processes
      !......................................................................

      !.......................
      ! collection of droplets
      call ice_cldliq_collection(rho(k),t(k),rhofaci(k),&
           f1pr04,qitot_incld(k),qc_incld(k),nitot_incld(k),nc_incld(k),&
           qccol,nccol,qcshd,ncshdc)

      !....................
      ! collection of rain
      call ice_rain_collection(rho(k),t(k),rhofaci(k),&
           logn0r(k),f1pr07,f1pr08,qitot_incld(k),nitot_incld(k),qr_incld(k),&
           qrcol,nrcol)
      !...................................
      ! collection between ice categories

      !PMC nCat deleted lots of stuff here.

      !.............................................
      ! self-collection of ice
      call ice_self_collection(rho(k),rhofaci(k),&
           f1pr03,eii,qirim_incld(k),qitot_incld(k),nitot_incld(k),&
           nislf)

      !............................................................
      ! melting
      call ice_melting(rho(k),t(k),pres(k),rhofaci(k),&
           f1pr05,f1pr14,xxlv(k),xlf(k),dv,sc,mu,kap,&
           qv(k),qitot_incld(k),nitot_incld(k),&
           qimlt,nimlt)

      !............................................................
      ! calculate wet growth
      call ice_cldliq_wet_growth(rho(k),t(k),pres(k),rhofaci(k),&
           f1pr05,f1pr14,xxlv(k),xlf(k),dv,kap,mu,sc,&
           qv(k),qc_incld(k),qitot_incld(k),nitot_incld(k),qr_incld(k),log_wetgrowth,&
           qrcol,qccol,qwgrth,nrshdr,qcshd)

      !-----------------------------
      ! calcualte total inverse ice relaxation timescale combined for all ice categories
      ! note 'f1pr' values are normalized, so we need to multiply by N
      call calc_ice_relaxation_timescale(rho(k),t(k),rhofaci(k),&
           f1pr05,f1pr14,dv,mu,sc,qitot_incld(k),nitot_incld(k),&
           epsi,epsi_tot)

      !.........................
      ! calculate rime density
      call calc_rime_density(t(k),rhofaci(k),&
           f1pr02,acn(k),lamc(k),mu_c(k),qc_incld(k),qccol,&
           vtrmi1,rhorime_c)
      !............................................................
      ! contact and immersion freezing droplets
      call cldliq_immersion_freezing(t(k),&
           lamc(k),mu_c(k),cdist1(k),qc_incld(k),qc_relvar(k),&
           qcheti,ncheti)

      !............................................................
      ! immersion freezing of rain
      ! for future: get rid of log statements below for rain freezing
      call rain_immersion_freezing(t(k),&
           lamr(k),mu_r(k),cdistr(k),qr_incld(k),&
           qrheti,nrheti)

      !......................................
      ! rime splintering (Hallet-Mossop 1974)
      !PMC comment: Morrison and Milbrandt 2015 part 1 and 2016 part 3 both say
      !that Hallet-Mossop should be neglected if 1 category to compensate for
      !artificial smearing out of ice DSD

      !................................................
      ! condensation/evaporation/deposition/sublimation
      !   (use semi-analytic formulation)

      ! calculate rain evaporation including ventilation
      call calc_liq_relaxation_timescale(rho(k),f1r,f2r,     &
           dv,mu,sc,mu_r(k),lamr(k),cdistr(k),cdist(k),qr_incld(k),qc_incld(k), &
           epsr,epsc)

      call evaporate_sublimate_precip(qr_incld(k),qc_incld(k),nr_incld(k),qitot_incld(k), &
           lcldm(k),rcldm(k),qvs(k),ab,epsr,qv(k), &
           qrevp,nrevp)

      call ice_deposition_sublimation(qitot_incld(k), nitot_incld(k), t(k), &
           qvs(k),qvi(k),epsi,abi,qv(k), &
           qidep,qisub,nisub,qiberg)

444   continue

      !................................................................
      ! deposition/condensation-freezing nucleation
      call ice_nucleation(t(k),inv_rho(k),&
           nitot(k),naai(k),supi(k),odt,log_predictNc,&
           qinuc, ninuc)

      !................
      ! cloud water autoconversion
      ! NOTE: cloud_water_autoconversion must be called before droplet_self_collection
      call cloud_water_autoconversion(rho(k),qc_incld(k),nc_incld(k),qc_relvar(k),&
           qcaut,ncautc,ncautr)

      !............................
      ! self-collection of droplets
      call droplet_self_collection(rho(k),inv_rho(k),qc_incld(k),&
           mu_c(k),nu(k),ncautc,ncslf)

      !............................
      ! accretion of cloud by rain
      call cloud_rain_accretion(rho(k),inv_rho(k),&
           qc_incld(k),nc_incld(k), qr_incld(k),qc_relvar(k),&
           qcacc, ncacc)

      !.....................................
      ! self-collection and breakup of rain
      ! (breakup following modified Verlinde and Cotton scheme)
      call rain_self_collection(rho(k),qr_incld(k),nr_incld(k),&
           nrslf)

      ! Here we map the microphysics tendency rates back to CELL-AVERAGE quantities for updating
      ! cell-average quantities.
      call back_to_cell_average(lcldm(k), rcldm(k), icldm(k), qcacc, qrevp, qcaut,&
           ncacc, ncslf, ncautc, nrslf, nrevp, ncautr, qisub, nrshdr, qcheti,&
           qrcol, qcshd, qimlt, qccol, qrheti, nimlt, nccol, ncshdc, ncheti, nrcol, nislf,&
           qidep, nrheti, nisub, qinuc, ninuc, qiberg)

      !.................................................................
      ! conservation of water
      !.................................................................

      ! The microphysical process rates are computed above, based on the environmental conditions.
      ! The rates are adjusted here (where necessary) such that the sum of the sinks of mass cannot
      ! be greater than the sum of the sources, thereby resulting in overdepletion.
      !-- Limit ice process rates to prevent overdepletion of sources such that
      !   the subsequent adjustments are done with maximum possible rates for the
      !   time step.  (note: most ice rates are adjusted here since they must be done
      !   simultaneously (outside of iice-loops) to distribute reduction proportionally
      !   amongst categories.
      !PMC - might need to rethink above statement since only one category now.
      ! AaronDonahue: Do we need the below checks for the new definition of
      ! how qidep and qisub are derived?
      ! AaronDonahue: UPDATE, if we are following the implementation of MG
      ! then the answer appears to be YES.  There is a similar check in MG
      ! microphysics which limits qidep and qinuc, but does not limit qisub.
      ! So similar but slightly different.  The overall answer would be that
      ! qidep does need some limit.  The next questions are,
      !   1) Should we be taking qinuc into consideration too?
      !   2) Is MG correct in NOT limiting qisub?

      call prevent_ice_overdepletion(pres(k), t(k), qv(k), xxls(k), odt, qidep, qisub)

      ! vapor -- not needed, since all sinks already have limits imposed and the sum, therefore,
      !          cannot possibly overdeplete qv

      ! cloud
      call cloud_water_conservation(qc(k), dt, qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep)

      ! rain
      call rain_water_conservation(qr(k), qcaut, qcacc, qimlt, qcshd, dt, qrevp, qrcol, qrheti)

      ! ice
      call ice_water_conservation(qitot(k), qidep, qinuc, qiberg, qrcol, qccol, qrheti, qcheti, dt, qisub, qimlt)

      !---------------------------------------------------------------------------------
      ! update prognostic microphysics and thermodynamics variables
      !---------------------------------------------------------------------------------

      !-- ice-phase dependent processes:
      call update_prognostic_ice(qcheti, qccol, qcshd, &
           nccol, ncheti, ncshdc, &
           qrcol, nrcol,  qrheti, nrheti, nrshdr, &
           qimlt, nimlt, qisub, qidep, qinuc, ninuc, nislf, nisub, qiberg, &
           exner(k), xxls(k), xlf(k), &
           log_predictNc, log_wetgrowth, dt, nmltratio, rhorime_c, &
           th(k), qv(k), qitot(k), nitot(k), qirim(k), birim(k), qc(k), nc(k), qr(k), nr(k) )

      !-- warm-phase only processes:
      call update_prognostic_liquid(qcacc, ncacc, qcaut, ncautc, ncautr, ncslf,  &
           qrevp, nrevp, nrslf,                                                  &
           log_predictNc, inv_rho(k), exner(k), xxlv(k), dt,                     &
           th(k), qv(k), qc(k), nc(k), qr(k), nr(k))

      !==
      ! AaronDonahue - Add extra variables needed from microphysics by E3SM:
      cmeiout(k) = qidep - qisub + qinuc
      prain(k)   = ( qcacc + qcaut + qcshd + qccol )
      nevapr(k)  = qisub + qrevp
      prer_evap(k) = qrevp
      vap_ice_exchange(k) = qidep - qisub + qinuc
      vap_liq_exchange(k) = - qrevp
      liq_ice_exchange(k) = qcheti + qrheti - qimlt + qiberg + qccol + qrcol

      ! clipping for small hydrometeor values
      if (qc(k).lt.qsmall) then
         qv(k) = qv(k) + qc(k)
         th(k) = th(k) - exner(k)*qc(k)*xxlv(k)*inv_cp
         qc(k) = 0._rtype
         nc(k) = 0._rtype
      else
         log_hydrometeorsPresent = .true.
      endif

      if (qr(k).lt.qsmall) then
         qv(k) = qv(k) + qr(k)
         th(k) = th(k) - exner(k)*qr(k)*xxlv(k)*inv_cp
         qr(k) = 0._rtype
         nr(k) = 0._rtype
      else
         log_hydrometeorsPresent = .true.
      endif

      if (qitot(k).lt.qsmall) then
         qv(k) = qv(k) + qitot(k)
         th(k) = th(k) - exner(k)*qitot(k)*xxls(k)*inv_cp
         qitot(k) = 0._rtype
         nitot(k) = 0._rtype
         qirim(k) = 0._rtype
         birim(k) = 0._rtype
      else
         log_hydrometeorsPresent = .true.
      endif

      call impose_max_total_Ni(nitot(k),max_total_Ni,inv_rho(k))

      ! Record microphysics tendencies for output:
      ! warm-phase process rates
      p3_tend_out(k, 2) = qcacc     ! cloud droplet accretion by rain
      p3_tend_out(k, 3) = qcaut     ! cloud droplet autoconversion to rain
      p3_tend_out(k, 4) = ncacc     ! change in cloud droplet number from accretion by rain
      p3_tend_out(k, 5) = ncautc    ! change in cloud droplet number from autoconversion
      p3_tend_out(k, 6) = ncslf     ! change in cloud droplet number from self-collection  (Not in paper?)
      p3_tend_out(k, 7) = nrslf     ! change in rain number from self-collection  (Not in paper?)
      p3_tend_out(k,11) = qrevp     ! rain evaporation
      p3_tend_out(k,13) = nrevp     ! change in rain number from evaporation
      p3_tend_out(k,14) = ncautr    ! change in rain number from autoconversion of cloud water
      ! ice-phase  process rates
      p3_tend_out(k,15) = qccol     ! collection of cloud water by ice
      p3_tend_out(k,16) = qwgrth    ! wet growth rate
      p3_tend_out(k,17) = qidep     ! vapor deposition
      p3_tend_out(k,18) = qrcol     ! collection rain mass by ice
      p3_tend_out(k,19) = qinuc     ! deposition/condensation freezing nuc
      p3_tend_out(k,20) = nccol     ! change in cloud droplet number from collection by ice
      p3_tend_out(k,21) = nrcol     ! change in rain number from collection by ice
      p3_tend_out(k,22) = ninuc     ! change in ice number from deposition/cond-freezing nucleation
      p3_tend_out(k,23) = qisub     ! sublimation of ice
      p3_tend_out(k,24) = qimlt     ! melting of ice
      p3_tend_out(k,25) = nimlt     ! melting of ice
      p3_tend_out(k,26) = nisub     ! change in ice number from sublimation
      p3_tend_out(k,27) = nislf     ! change in ice number from collection within a category (Not in paper?)
      p3_tend_out(k,28) = qcheti    ! immersion freezing droplets
      p3_tend_out(k,29) = qrheti    ! immersion freezing rain
      p3_tend_out(k,30) = ncheti    ! immersion freezing droplets
      p3_tend_out(k,31) = nrheti    ! immersion freezing rain
      p3_tend_out(k,32) = nrshdr    ! source for rain number from collision of rain/ice above freezing and shedding
      p3_tend_out(k,33) = qcshd     ! source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding
      p3_tend_out(k,34) = 0._rtype  ! used to be qcmul, but that has been removed.  Kept at 0.0 as placeholder.
      p3_tend_out(k,35) = ncshdc    ! source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)
      ! Outputs associated with aerocom comparison:
      pratot(k) = qcacc ! cloud drop accretion by rain
      prctot(k) = qcaut ! cloud drop autoconversion to rain
      !---------------------------------------------------------------------------------

      ! Recalculate in-cloud values for sedimentation
      call calculate_incloud_mixingratios(qc(k),qr(k),qitot(k),qirim(k),nc(k),nr(k),nitot(k),birim(k), &
           inv_lcldm(k),inv_icldm(k),inv_rcldm(k), &
           qc_incld(k),qr_incld(k),qitot_incld(k),qirim_incld(k),nc_incld(k),nr_incld(k),nitot_incld(k),birim_incld(k))

555   continue

   enddo k_loop_main

 END SUBROUTINE p3_main_part2

 subroutine p3_main_part3(kts, kte, kbot, ktop, kdir, &
      exner, lcldm, rcldm, &
      rho, inv_rho, rhofaci, qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls, &
      mu_c, nu, lamc, mu_r, lamr, vap_liq_exchange, &
      ze_rain, ze_ice, diag_vmi, diag_effi, diag_di, diag_rhoi, diag_ze, diag_effc)

   implicit none

   ! args

   integer, intent(in) :: kts, kte, kbot, ktop, kdir

   real(rtype), intent(in), dimension(kts:kte) :: exner, lcldm, rcldm

   real(rtype), intent(inout), dimension(kts:kte) :: rho, inv_rho, rhofaci, &
        qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls, &
        mu_c, nu, lamc, mu_r, &
        lamr, vap_liq_exchange, &
        ze_rain, ze_ice, diag_vmi, diag_effi, diag_di, diag_rhoi, diag_ze, diag_effc

   ! locals
   integer :: k, dumi, dumii, dumjj, dumzz
   real(rtype) :: tmp1, tmp2, dum1, dum4, dum5, dum6, rhop
   real(rtype)    :: f1pr02   ! mass-weighted fallspeed              See lines  731 -  808  ums
   real(rtype)    :: f1pr06   ! effective radius                     See lines 1281 - 1356  eff
   real(rtype)    :: f1pr09   ! minimum ice number (lambda limiter)  See lines  704 -  705  nlarge
   real(rtype)    :: f1pr10   ! maximum ice number (lambda limiter)  See lines  704 -  705  nsmall
   real(rtype)    :: f1pr13   ! reflectivity                         See lines  731 -  808  refl
   real(rtype)    :: f1pr15   ! mass-weighted mean diameter          See lines 1212 - 1279  dmm
   real(rtype)    :: f1pr16   ! mass-weighted mean particle density  See lines 1212 - 1279  rhomm

   k_loop_final_diagnostics:  do k = kbot,ktop,kdir

      ! cloud:
      if (qc(k).ge.qsmall) then
         call get_cloud_dsd2(qc(k),nc(k),mu_c(k),rho(k),nu(k),dnu,lamc(k),  &
              tmp1,tmp2,lcldm(k))
         diag_effc(k) = 0.5_rtype*(mu_c(k)+3._rtype)/lamc(k)
      else
         qv(k) = qv(k)+qc(k)
         th(k) = th(k)-exner(k)*qc(k)*xxlv(k)*inv_cp
         vap_liq_exchange(k) = vap_liq_exchange(k) - qc(k)
         qc(k) = 0._rtype
         nc(k) = 0._rtype
      endif

      ! rain:
      if (qr(k).ge.qsmall) then

         call get_rain_dsd2(qr(k),nr(k),mu_r(k),lamr(k),tmp1,tmp2,rcldm(k))

         ze_rain(k) = nr(k)*(mu_r(k)+6._rtype)*(mu_r(k)+5._rtype)*(mu_r(k)+4._rtype)*           &
              (mu_r(k)+3._rtype)*(mu_r(k)+2._rtype)*(mu_r(k)+1._rtype)/bfb_pow(lamr(k), 6._rtype)
         ze_rain(k) = max(ze_rain(k),1.e-22_rtype)
      else
         qv(k) = qv(k)+qr(k)
         th(k) = th(k)-exner(k)*qr(k)*xxlv(k)*inv_cp
         vap_liq_exchange(k) = vap_liq_exchange(k) - qr(k)
         qr(k) = 0._rtype
         nr(k) = 0._rtype
      endif

      ! ice:

      call impose_max_total_Ni(nitot(k),max_total_Ni,inv_rho(k))

      qi_not_small:  if (qitot(k).ge.qsmall) then

         !impose lower limits to prevent taking log of # < 0
         nitot(k) = max(nitot(k),nsmall)
         nr(k)    = max(nr(k),nsmall)

         call calc_bulkRhoRime(qitot(k),qirim(k),birim(k),rhop)

         ! if (.not. tripleMoment_on) zitot(k) = diag_mom6(qitot(k),nitot(k),rho(k))
         call find_lookupTable_indices_1a(dumi,dumjj,dumii,dumzz,dum1,dum4,          &
              dum5,dum6,isize,rimsize,densize,     &
              qitot(k),nitot(k),           &
              qirim(k),rhop)
         !qirim(k),zitot(k),rhop)

         call access_lookup_table(dumjj,dumii,dumi, 2,dum1,dum4,dum5,f1pr02)
         call access_lookup_table(dumjj,dumii,dumi, 6,dum1,dum4,dum5,f1pr06)
         call access_lookup_table(dumjj,dumii,dumi, 7,dum1,dum4,dum5,f1pr09)
         call access_lookup_table(dumjj,dumii,dumi, 8,dum1,dum4,dum5,f1pr10)
         call access_lookup_table(dumjj,dumii,dumi, 9,dum1,dum4,dum5,f1pr13)
         call access_lookup_table(dumjj,dumii,dumi,11,dum1,dum4,dum5,f1pr15)
         call access_lookup_table(dumjj,dumii,dumi,12,dum1,dum4,dum5,f1pr16)

         ! impose mean ice size bounds (i.e. apply lambda limiters)
         ! note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
         nitot(k) = min(nitot(k),f1pr09*nitot(k))
         nitot(k) = max(nitot(k),f1pr10*nitot(k))

         !--this should already be done in s/r 'calc_bulkRhoRime'
         if (qirim(k).lt.qsmall) then
            qirim(k) = 0._rtype
            birim(k) = 0._rtype
         endif
         !==

         ! note that reflectivity from lookup table is normalized, so we need to multiply by N
         diag_vmi(k)   = f1pr02*rhofaci(k)
         diag_effi(k)  = f1pr06 ! units are in m
         diag_di(k)    = f1pr15
         diag_rhoi(k)  = f1pr16
         ! note factor of air density below is to convert from m^6/kg to m^6/m^3
         ze_ice(k) = ze_ice(k) + 0.1892_rtype*f1pr13*nitot(k)*rho(k)   ! sum contribution from each ice category (note: 0.1892 = 0.176/0.93)
         ze_ice(k) = max(ze_ice(k),1.e-22_rtype)

      else

         qv(k) = qv(k) + qitot(k)
         th(k) = th(k) - exner(k)*qitot(k)*xxls(k)*inv_cp
         qitot(k) = 0._rtype
         nitot(k) = 0._rtype
         qirim(k) = 0._rtype
         birim(k) = 0._rtype
         diag_di(k) = 0._rtype

      endif qi_not_small

      ! sum ze components and convert to dBZ
      diag_ze(k) = 10._rtype*bfb_log10((ze_rain(k) + ze_ice(k))*1.e18_rtype)

      ! if qr is very small then set Nr to 0 (needs to be done here after call
      ! to ice lookup table because a minimum Nr of nsmall will be set otherwise even if qr=0)
      if (qr(k).lt.qsmall) then
         nr(k) = 0._rtype
      endif

   enddo k_loop_final_diagnostics

 end subroutine p3_main_part3

  !==========================================================================================!

  SUBROUTINE p3_main(qc,nc,qr,nr,th,qv,dt,qitot,qirim,nitot,birim,   &
       pres,dzq,ncnuc,naai,qc_relvar,it,prt_liq,prt_sol,its,ite,kts,kte,diag_ze,diag_effc,     &
       diag_effi,diag_vmi,diag_di,diag_rhoi,log_predictNc, &
       pdel,exner,cmeiout,prain,nevapr,prer_evap,rflx,sflx,rcldm,lcldm,icldm,  &
       pratot,prctot,p3_tend_out,mu_c,lamc,liq_ice_exchange,vap_liq_exchange, &
       vap_ice_exchange,col_location)

    !----------------------------------------------------------------------------------------!
    !                                                                                        !
    ! This is the main subroutine for the P3 microphysics scheme.  It is called from the     !
    ! wrapper subroutine ('MP_P3_WRAPPER') and is passed i,k slabs of all prognostic         !
    ! variables -- hydrometeor fields, potential temperature, and water vapor mixing ratio.  !
    ! Microphysical process rates are computed first.  These tendencies are then used to     !
    ! computed updated values of the prognostic variables.  The hydrometeor variables are    !
    ! then updated further due to sedimentation.                                             !
    !                                                                                        !
    ! Several diagnostic values are also computed and returned to the wrapper subroutine,    !
    ! including precipitation rates.                                                         !
    !                                                                                        !
    !----------------------------------------------------------------------------------------!

    use debug_info, only: get_debug_column_id
    implicit none

    !----- Input/ouput arguments:  ----------------------------------------------------------!

    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qc         ! cloud, mass mixing ratio         kg kg-1
    ! note: Nc may be specified or predicted (set by log_predictNc)
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: nc         ! cloud, number mixing ratio       #  kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qr         ! rain, mass mixing ratio          kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: nr         ! rain, number mixing ratio        #  kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qitot      ! ice, total mass mixing ratio     kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qirim      ! ice, rime mass mixing ratio      kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: nitot      ! ice, total number mixing ratio   #  kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: birim      ! ice, rime volume mixing ratio    m3 kg-1

    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qv         ! water vapor mixing ratio         kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: th         ! potential temperature            K
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: pres       ! pressure                         Pa
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: dzq        ! vertical grid spacing            m
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: ncnuc      ! IN ccn activated number tendency kg-1 s-1
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: naai       ! IN actived ice nuclei concentration  1/kg
    real(rtype), intent(in)                                     :: dt         ! model time step                  s

    real(rtype), intent(out),   dimension(its:ite)              :: prt_liq    ! precipitation rate, liquid       m s-1
    real(rtype), intent(out),   dimension(its:ite)              :: prt_sol    ! precipitation rate, solid        m s-1
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_ze    ! equivalent reflectivity          dBZ
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_effc  ! effective radius, cloud          m
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_effi  ! effective radius, ice            m
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_vmi   ! mass-weighted fall speed of ice  m s-1
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_di    ! mean diameter of ice             m
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_rhoi  ! bulk density of ice              kg m-3
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: mu_c       ! Size distribution shape parameter for radiation
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: lamc       ! Size distribution slope parameter for radiation

    integer, intent(in)                                  :: its,ite    ! array bounds (horizontal)
    integer, intent(in)                                  :: kts,kte    ! array bounds (vertical)
    integer, intent(in)                                  :: it         ! time step counter NOTE: starts at 1 for first time step

    logical(btype), intent(in)                           :: log_predictNc ! .T. (.F.) for prediction (specification) of Nc

    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: pdel       ! pressure thickness               Pa
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: exner      ! Exner expression

    ! OUTPUT for PBUF variables used by other parameterizations
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: cmeiout    ! qitend due to deposition/sublimation
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: prain      ! Total precipitation (rain + snow)
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: nevapr     ! evaporation of total precipitation (rain + snow)
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: prer_evap  ! evaporation of rain
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: rflx       ! grid-box average rain flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: sflx       ! grid-box average ice/snow flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: pratot     ! accretion of cloud by rain
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: prctot     ! autoconversion of cloud to rain
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: liq_ice_exchange ! sum of liq-ice phase change tendenices
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: vap_liq_exchange ! sum of vap-liq phase change tendenices
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: vap_ice_exchange ! sum of vap-ice phase change tendenices
    ! INPUT needed for PBUF variables used by other parameterizations

    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: icldm, lcldm, rcldm ! Ice, Liquid and Rain cloud fraction
    ! AaronDonahue, the following variable (p3_tend_out) is a catch-all for passing P3-specific variables outside of p3_main
    ! so that they can be written as ouput.  NOTE TO C++ PORT: This variable is entirely optional and doesn't need to be
    ! included in the port to C++, or can be changed if desired.
    real(rtype), intent(out),   dimension(its:ite,kts:kte,49)   :: p3_tend_out ! micro physics tendencies
    real(rtype), intent(in),    dimension(its:ite,3)            :: col_location
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: qc_relvar

    !
    !----- Local variables and parameters:  -------------------------------------------------!
    !

    real(rtype), dimension(its:ite,kts:kte) :: mu_r  ! shape parameter of rain
    real(rtype), dimension(its:ite,kts:kte) :: t     ! temperature at the beginning of the microhpysics step [K]

    ! 2D size distribution and fallspeed parameters:

    real(rtype), dimension(its:ite,kts:kte) :: lamr
    real(rtype), dimension(its:ite,kts:kte) :: logn0r

    real(rtype), dimension(its:ite,kts:kte) :: nu
    real(rtype), dimension(its:ite,kts:kte) :: cdist
    real(rtype), dimension(its:ite,kts:kte) :: cdist1
    real(rtype), dimension(its:ite,kts:kte) :: cdistr

    ! Variables needed for in-cloud calculations
    real(rtype), dimension(its:ite,kts:kte) :: inv_icldm, inv_lcldm, inv_rcldm ! Inverse cloud fractions (1/cld)
    real(rtype), dimension(its:ite,kts:kte) :: qc_incld, qr_incld, qitot_incld, qirim_incld ! In cloud mass-mixing ratios
    real(rtype), dimension(its:ite,kts:kte) :: nc_incld, nr_incld, nitot_incld, birim_incld ! In cloud number concentrations

    real(rtype), dimension(its:ite,kts:kte)      :: inv_dzq,inv_rho,ze_ice,ze_rain,prec,rho,       &
         rhofacr,rhofaci,acn,xxls,xxlv,xlf,qvs,qvi,supi,       &
         tmparr1,inv_exner

    ! -- scalar locals -- !

    real(rtype) :: odt, timeScaleFactor

    integer :: ktop,kbot,kdir,i

    logical(btype) :: log_nucleationPossible, log_hydrometeorsPresent

    !--These will be added as namelist parameters in the future
    logical(btype), parameter :: debug_ON     = .true.  !.true. to switch on debugging checks/traps throughout code  TODO: Turn this back off as default once the tlay error is found.
    logical(btype), parameter :: debug_ABORT  = .false.  !.true. will result in forced abort in s/r 'check_values'

    real(rtype),dimension(its:ite,kts:kte) :: qc_old, nc_old, qr_old, nr_old, qitot_old, nitot_old, qv_old, th_old

    !-----------------------------------------------------------------------------------!
    !  End of variables/parameters declarations
    !-----------------------------------------------------------------------------------!

    ! direction of vertical leveling:
    !PMC got rid of 'model' option so we could just replace ktop with kts everywhere...
    ktop = kts        !k of top level
    kbot = kte        !k of bottom level
    kdir = -1         !(k: 1=top, nk=bottom)

    !PMC deleted 'threshold size difference' calculation for multicategory here

    inv_dzq    = 1._rtype/dzq  ! inverse of thickness of layers
    odt        = 1._rtype/dt   ! inverse model time step

    ! Compute time scale factor over which to apply soft rain lambda limiter
    ! note: '1./max(30.,dt)' = '1.*min(1./30., 1./dt)'
    timeScaleFactor = min(1._rtype/120._rtype, odt)

    prt_liq   = 0._rtype
    prt_sol   = 0._rtype
    pratot    = 0._rtype
    prctot    = 0._rtype
    prec      = 0._rtype
    mu_r      = 0._rtype
    diag_ze   = -99._rtype

    ze_ice    = 1.e-22_rtype
    ze_rain   = 1.e-22_rtype
    diag_effc = 10.e-6_rtype ! default value
    diag_effi = 25.e-6_rtype ! default value
    diag_vmi  = 0._rtype
    diag_di   = 0._rtype
    diag_rhoi = 0._rtype

    cmeiout = 0._rtype
    prain   = 0._rtype
    nevapr  = 0._rtype
    rflx    = 0._rtype
    sflx    = 0._rtype
    p3_tend_out = 0._rtype

    inv_icldm = 1.0_rtype/icldm
    inv_lcldm = 1.0_rtype/lcldm
    inv_rcldm = 1.0_rtype/rcldm

    prer_evap = 0._rtype
    liq_ice_exchange = 0._rtype
    vap_liq_exchange = 0._rtype
    vap_ice_exchange = 0._rtype

    mu_c = 0.0_rtype
    lamc = 0.0_rtype
    ! AaronDonahue added exner term to replace all instances of th(i,k)/t(i,k), since th(i,k) is updated but t(i,k) is not, and this was
    ! causing energy conservation errors.
    inv_exner = 1._rtype/exner        !inverse of Exner expression, used when converting potential temp to temp
    t       = th    *inv_exner    !compute temperature from theta (value at beginning of microphysics step)
    qv      = max(qv,0._rtype)        !clip water vapor to prevent negative values passed in (beginning of microphysics)
    ! AaronDonahue added this load of latent heat to be consistent with E3SM, since the inconsistentcy was causing water conservation errors.
    call get_latent_heat(its,ite,kts,kte,xxlv,xxls,xlf)

   ! initialize microphysics processes tendency output
    qc_old = qc         ! Liq. microphysics tendency, initialize
    nc_old = nc         ! Liq. # microphysics tendency, initialize
    qr_old = qr         ! Rain microphysics tendency, initialize
    nr_old = nr         ! Rain # microphysics tendency, initialize
    qitot_old = qitot   ! Ice  microphysics tendency, initialize
    nitot_old = nitot   ! Ice  # microphysics tendency, initialize
    qv_old = qv         ! Vapor  microphysics tendency, initialize
    th_old = th         ! Pot. Temp. microphysics tendency, initialize

    !==
    !-----------------------------------------------------------------------------------!
    i_loop_main: do i = its,ite  ! main i-loop (around the entire scheme)

       !update column in the debug_info module
       call get_debug_column_id(i)

!      if (debug_ON) call check_values(qv,T,i,it,debug_ABORT,100,col_location)

       call p3_main_part1(kts, kte, kbot, ktop, kdir, log_predictNc, dt, &
            pres(i,:), pdel(i,:), dzq(i,:), ncnuc(i,:), exner(i,:), inv_exner(i,:), inv_lcldm(i,:), inv_icldm(i,:), inv_rcldm(i,:), xxlv(i,:), xxls(i,:), xlf(i,:), &
            t(i,:), rho(i,:), inv_rho(i,:), qvs(i,:), qvi(i,:), supi(i,:), rhofacr(i,:), rhofaci(i,:), acn(i,:), qv(i,:), th(i,:), qc(i,:), nc(i,:), qr(i,:), nr(i,:), &
            qitot(i,:), nitot(i,:), qirim(i,:), birim(i,:), qc_incld(i,:), qr_incld(i,:), qitot_incld(i,:), qirim_incld(i,:), &
            nc_incld(i,:), nr_incld(i,:), nitot_incld(i,:), birim_incld(i,:), log_nucleationPossible, log_hydrometeorsPresent)

!      if (debug_ON) then
!         tmparr1(i,:) = th(i,:)*inv_exner(i,:)!(pres(i,:)*1.e-5)**(rd*inv_cp)
!         call check_values(qv,tmparr1,i,it,debug_ABORT,200,col_location)
!      endif

       !jump to end of i-loop if log_nucleationPossible=.false.  (i.e. skip everything)
       if (.not. (log_nucleationPossible .or. log_hydrometeorsPresent)) goto 333

       call p3_main_part2(kts, kte, kbot, ktop, kdir, log_predictNc, dt, odt, &
            pres(i,:), pdel(i,:), dzq(i,:), ncnuc(i,:), exner(i,:), inv_exner(i,:), inv_lcldm(i,:), inv_icldm(i,:), inv_rcldm(i,:), naai(i,:), qc_relvar(i,:), icldm(i,:), lcldm(i,:), rcldm(i,:),&
            t(i,:), rho(i,:), inv_rho(i,:), qvs(i,:), qvi(i,:), supi(i,:), rhofacr(i,:), rhofaci(i,:), acn(i,:), qv(i,:), th(i,:), qc(i,:), nc(i,:), qr(i,:), nr(i,:), qitot(i,:), nitot(i,:), &
            qirim(i,:), birim(i,:), xxlv(i,:), xxls(i,:), xlf(i,:), qc_incld(i,:), qr_incld(i,:), qitot_incld(i,:), qirim_incld(i,:), nc_incld(i,:), nr_incld(i,:), &
            nitot_incld(i,:), birim_incld(i,:), mu_c(i,:), nu(i,:), lamc(i,:), cdist(i,:), cdist1(i,:), cdistr(i,:), mu_r(i,:), lamr(i,:), logn0r(i,:), cmeiout(i,:), prain(i,:), &
            nevapr(i,:), prer_evap(i,:), vap_liq_exchange(i,:), vap_ice_exchange(i,:), liq_ice_exchange(i,:), pratot(i,:), &
            prctot(i,:), p3_tend_out(i,:,:), log_hydrometeorsPresent)

       ! measure microphysics processes tendency output
       p3_tend_out(i,:,42) = ( qc(i,:)    - qc_old(i,:) ) * odt    ! Liq. microphysics tendency, measure
       p3_tend_out(i,:,43) = ( nc(i,:)    - nc_old(i,:) ) * odt    ! Liq. # microphysics tendency, measure
       p3_tend_out(i,:,44) = ( qr(i,:)    - qr_old(i,:) ) * odt    ! Rain microphysics tendency, measure
       p3_tend_out(i,:,45) = ( nr(i,:)    - nr_old(i,:) ) * odt    ! Rain # microphysics tendency, measure
       p3_tend_out(i,:,46) = ( qitot(i,:) - qitot_old(i,:) ) * odt ! Ice  microphysics tendency, measure
       p3_tend_out(i,:,47) = ( nitot(i,:) - nitot_old(i,:) ) * odt ! Ice  # microphysics tendency, measure
       p3_tend_out(i,:,48) = ( qv(i,:)    - qv_old(i,:) ) * odt    ! Vapor  microphysics tendency, measure
       p3_tend_out(i,:,49) = ( th(i,:)    - th_old(i,:) ) * odt    ! Pot. Temp. microphysics tendency, measure
       !NOTE: At this point, it is possible to have negative (but small) nc, nr, nitot.  This is not
       !      a problem; those values get clipped to zero in the sedimentation section (if necessary).
       !      (This is not done above simply for efficiency purposes.)

       !      if (debug_ON) then
       !         tmparr1(i,:) = th(i,:)*inv_exner(i,:)!(pres(i,:)*1.e-5)**(rd*inv_cp)
       !         call check_values(qv,tmparr1,i,it,debug_ABORT,300,col_location)
       !      endif

       if (.not. log_hydrometeorsPresent) goto 333

       !------------------------------------------------------------------------------------------!
       ! End of main microphysical processes section
       !==========================================================================================!

       !==========================================================================================!
       ! Sedimentation:

       !------------------------------------------------------------------------------------------!
       ! Cloud sedimentation:  (adaptive substepping)
       p3_tend_out(i,:,36) = qc(i,:) ! Liq. sedimentation tendency, initialize
       p3_tend_out(i,:,37) = nc(i,:) ! Liq. # sedimentation tendency, initialize

       call cloud_sedimentation(kts,kte,ktop,kbot,kdir, &
         qc_incld(i,:),rho(i,:),inv_rho(i,:),lcldm(i,:),acn(i,:),inv_dzq(i,:), &
         dt,odt,dnu,log_predictNc, &
         qc(i,:),nc(i,:),nc_incld(i,:),mu_c(i,:),lamc(i,:),prt_liq(i),p3_tend_out(i,:,36),p3_tend_out(i,:,37))

       !------------------------------------------------------------------------------------------!
       ! Rain sedimentation:  (adaptive substepping)
       p3_tend_out(i,:,38) = qr(i,:) ! Rain sedimentation tendency, initialize
       p3_tend_out(i,:,39) = nr(i,:) ! Rain # sedimentation tendency, initialize

       call rain_sedimentation(kts,kte,ktop,kbot,kdir, &
         qr_incld(i,:),rho(i,:),inv_rho(i,:),rhofacr(i,:),rcldm(i,:),inv_dzq(i,:),dt,odt, &
         qr(i,:),nr(i,:),nr_incld(i,:),mu_r(i,:),lamr(i,:),prt_liq(i),rflx(i,:),p3_tend_out(i,:,38),p3_tend_out(i,:,39))

       !------------------------------------------------------------------------------------------!
       ! Ice sedimentation:  (adaptive substepping)
       p3_tend_out(i,:,40) = qitot(i,:) ! Ice sedimentation tendency, initialize
       p3_tend_out(i,:,41) = nitot(i,:) ! Ice # sedimentation tendency, initialize

       call ice_sedimentation(kts,kte,ktop,kbot,kdir,    &
         rho(i,:),inv_rho(i,:),rhofaci(i,:),icldm(i,:),inv_dzq(i,:),dt,odt, &
         qitot(i,:),qitot_incld(i,:),nitot(i,:),qirim(i,:),qirim_incld(i,:),birim(i,:),birim_incld(i,:),nitot_incld(i,:), &
         prt_sol(i),p3_tend_out(i,:,40),p3_tend_out(i,:,41))

       !.......................................
       ! homogeneous freezing of cloud and rain

       call homogeneous_freezing(kts,kte,ktop,kbot,kdir,t(i,:),exner(i,:),xlf(i,:),  &
         qc(i,:),nc(i,:),qr(i,:),nr(i,:),qitot(i,:),nitot(i,:),qirim(i,:),birim(i,:),th(i,:))

       !...................................................
       ! final checks to ensure consistency of mass/number
       ! and compute diagnostic fields for output
       call p3_main_part3(kts, kte, kbot, ktop, kdir, &
            exner(i,:), lcldm(i,:), rcldm(i,:), &
            rho(i,:), inv_rho(i,:), rhofaci(i,:), qv(i,:), th(i,:), qc(i,:), nc(i,:), qr(i,:), nr(i,:), qitot(i,:), nitot(i,:), &
            qirim(i,:), birim(i,:), xxlv(i,:), xxls(i,:), &
            mu_c(i,:), nu(i,:), lamc(i,:), mu_r(i,:), lamr(i,:), vap_liq_exchange(i,:), &
            ze_rain(i,:), ze_ice(i,:), diag_vmi(i,:), diag_effi(i,:), diag_di(i,:), diag_rhoi(i,:), diag_ze(i,:), diag_effc(i,:))
       !   if (debug_ON) call check_values(qv,Ti,it,debug_ABORT,800,col_location)

       !..............................................
       ! merge ice categories with similar properties

       !   note:  this should be relocated to above, such that the diagnostic
       !          ice properties are computed after merging

       !PMC nCat deleted nCat>1 stuff

       !.....................................................

333    continue

       if (debug_ON) then
          tmparr1(i,:) = th(i,:)*inv_exner(i,:)!(pres(i,:)*1.e-5)**(rd*inv_cp)
          call check_values(qv(i,:),tmparr1(i,:),kts,kte,it,debug_ABORT,900,col_location(i,:))
       endif

       !.....................................................

    enddo i_loop_main

    !PMC deleted "if WRF" stuff
    !PMC deleted typeDiags optional output stuff

    !=== (end of section for diagnostic hydrometeor/precip types)

    ! end of main microphysics routine

    return

  END SUBROUTINE p3_main

  !==========================================================================================!

  SUBROUTINE access_lookup_table(dumjj,dumii,dumi,index,dum1,dum4,dum5,proc)

    implicit none

    real(rtype)    :: dum1,dum4,dum5,proc,iproc1,gproc1,tmp1,tmp2
    integer :: dumjj,dumii,dumi,index

    ! get value at current density index

    ! first interpolate for current rimed fraction index
    iproc1 = itab(dumjj,dumii,dumi,index)+(dum1-real(dumi))*(itab(dumjj,dumii,       &
         dumi+1,index)-itab(dumjj,dumii,dumi,index))

    ! linearly interpolate to get process rates for rimed fraction index + 1

    gproc1 = itab(dumjj,dumii+1,dumi,index)+(dum1-real(dumi))*(itab(dumjj,dumii+1,   &
         dumi+1,index)-itab(dumjj,dumii+1,dumi,index))

    tmp1   = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

    ! get value at density index + 1

    ! first interpolate for current rimed fraction index

    iproc1 = itab(dumjj+1,dumii,dumi,index)+(dum1-real(dumi))*(itab(dumjj+1,dumii,   &
         dumi+1,index)-itab(dumjj+1,dumii,dumi,index))

    ! linearly interpolate to get process rates for rimed fraction index + 1

    gproc1 = itab(dumjj+1,dumii+1,dumi,index)+(dum1-real(dumi))*(itab(dumjj+1,       &
         dumii+1,dumi+1,index)-itab(dumjj+1,dumii+1,dumi,index))

    tmp2   = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

    ! get final process rate
    proc   = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)
    return
  END SUBROUTINE access_lookup_table

  !------------------------------------------------------------------------------------------!
  SUBROUTINE access_lookup_table_coll(dumjj,dumii,dumj,dumi,index,dum1,dum3,          &
       dum4,dum5,proc)

    implicit none

    real(rtype)    :: dum1,dum3,dum4,dum5,proc,dproc1,dproc2,iproc1,gproc1,tmp1,tmp2
    integer :: dumjj,dumii,dumj,dumi,index

    ! This subroutine interpolates lookup table values for rain/ice collection processes

    ! current density index

    ! current rime fraction index
    dproc1  = itabcoll(dumjj,dumii,dumi,dumj,index)+(dum1-real(dumi))*                &
         (itabcoll(dumjj,dumii,dumi+1,dumj,index)-itabcoll(dumjj,dumii,dumi,    &
         dumj,index))

    dproc2  = itabcoll(dumjj,dumii,dumi,dumj+1,index)+(dum1-real(dumi))*             &
         (itabcoll(dumjj,dumii,dumi+1,dumj+1,index)-itabcoll(dumjj,dumii,dumi,  &
         dumj+1,index))

    iproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)

    ! rime fraction index + 1

    dproc1  = itabcoll(dumjj,dumii+1,dumi,dumj,index)+(dum1-real(dumi))*             &
         (itabcoll(dumjj,dumii+1,dumi+1,dumj,index)-itabcoll(dumjj,dumii+1,     &
         dumi,dumj,index))

    dproc2  = itabcoll(dumjj,dumii+1,dumi,dumj+1,index)+(dum1-real(dumi))*           &
         (itabcoll(dumjj,dumii+1,dumi+1,dumj+1,index)-itabcoll(dumjj,dumii+1,   &
         dumi,dumj+1,index))

    gproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)
    tmp1    = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

    ! density index + 1

    ! current rime fraction index

    dproc1  = itabcoll(dumjj+1,dumii,dumi,dumj,index)+(dum1-real(dumi))*             &
         (itabcoll(dumjj+1,dumii,dumi+1,dumj,index)-itabcoll(dumjj+1,dumii,     &
         dumi,dumj,index))

    dproc2  = itabcoll(dumjj+1,dumii,dumi,dumj+1,index)+(dum1-real(dumi))*           &
         (itabcoll(dumjj+1,dumii,dumi+1,dumj+1,index)-itabcoll(dumjj+1,dumii,   &
         dumi,dumj+1,index))

    iproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)

    ! rime fraction index + 1

    dproc1  = itabcoll(dumjj+1,dumii+1,dumi,dumj,index)+(dum1-real(dumi))*           &
         (itabcoll(dumjj+1,dumii+1,dumi+1,dumj,index)-itabcoll(dumjj+1,dumii+1, &
         dumi,dumj,index))

    dproc2  = itabcoll(dumjj+1,dumii+1,dumi,dumj+1,index)+(dum1-real(dumi))*         &
         (itabcoll(dumjj+1,dumii+1,dumi+1,dumj+1,index)-itabcoll(dumjj+1,       &
         dumii+1,dumi,dumj+1,index))

    gproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)
    tmp2    = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

    ! interpolate over density to get final values
    proc    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

    return
  END SUBROUTINE access_lookup_table_coll

  !==========================================================================================!

  real(rtype) function MurphyKoop_svp(t, i_type)

    use scream_abortutils, only : endscreamrun

    implicit none

    !-------------------------------------------------------------------
    ! Compute saturation vapor pressure (returned in units of pa)
    ! Inputs:
    ! "t", units [K]
    !  i_type refers to saturation with respect to liquid (0) or ice (1)
    !--------------------------------------------------------------------

    !Murphy & Koop (2005)
    real(rtype), intent(in) :: t
    integer, intent(in)     :: i_type

    !local vars
    character(len=1000) :: err_msg
    real(rtype)         :: logt, tmp

    !parameters for ice saturation eqn
    real(rtype), parameter :: ic(4)  =(/9.550426_rtype, 5723.265_rtype, 3.53068_rtype, &
         0.00728332_rtype/)

    !parameters for liq saturation eqn
    real(rtype), parameter :: lq(10) = (/54.842763_rtype, 6763.22_rtype, 4.210_rtype, &
         0.000367_rtype, 0.0415_rtype, 218.8_rtype, 53.878_rtype, 1331.22_rtype,       &
         9.44523_rtype, 0.014025_rtype /)

    logt = bfb_log(t)

    if (i_type .eq. 1 .and. t .lt. zerodegc) then

       !(good down to 110 K)
       MurphyKoop_svp = bfb_exp(ic(1) - (ic(2) / t) + (ic(3) * logt) - (ic(4) * t))

    elseif (i_type .eq. 0 .or. t .ge. zerodegc) then

       ! (good for 123 < T < 332 K)
       !For some reason, we cannot add line breaks if we use "bfb_exp", storing experssion in "tmp"
       tmp = lq(1) - (lq(2) / t) - (lq(3) * logt) + (lq(4) * t) + &
            (tanh(lq(5) * (t - lq(6))) * (lq(7) - (lq(8) / t) - &
            (lq(9) * logt) + lq(10) * t))
       MurphyKoop_svp = bfb_exp(tmp)
    else

       write(err_msg,*)'Error: Either MurphyKoop_svp i_type is not 0 or 1 or t=NaN. itype= ', &
            i_type,' and temperature t=',t,' in file: ',__FILE__,' at line:',__LINE__
       call endscreamrun(err_msg)
    endif

    return
  end function MurphyKoop_svp

  !_rtype
  real(rtype) function polysvp1(t,i_type)

    !-------------------------------------------
    !  COMPUTE SATURATION VAPOR PRESSURE
    !  POLYSVP1 RETURNED IN UNITS OF PA.
    !  T IS INPUT IN UNITS OF K.
    !  i_type REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
    !-------------------------------------------

    use scream_abortutils, only : endscreamrun

    use debug_info, only: report_error_info
    implicit none

    real(rtype), intent(in) :: t
    integer, intent(in)     :: i_type

    ! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

    !local variables
    character(len=1000) :: err_msg

    ! ice
    real(rtype) a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
    data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
         6.11147274_rtype,     0.503160820_rtype,     0.188439774e-1_rtype, &
         0.420895665e-3_rtype, 0.615021634e-5_rtype,  0.602588177e-7_rtype, &
         0.385852041e-9_rtype, 0.146898966e-11_rtype, 0.252751365e-14_rtype/

    ! liquid
    real(rtype) a0,a1,a2,a3,a4,a5,a6,a7,a8

    ! V1.7
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
         6.11239921_rtype,      0.443987641_rtype,     0.142986287e-1_rtype, &
         0.264847430e-3_rtype,  0.302950461e-5_rtype,  0.206739458e-7_rtype, &
         0.640689451e-10_rtype,-0.952447341e-13_rtype,-0.976195544e-15_rtype/
    real(rtype) dt

    !-------------------------------------------

    if (i_type.eq.1 .and. t.lt.zerodegc) then
       ! ICE

       !       Flatau formulation:
       dt       = max(-80._rtype,t-273.15_rtype)
       polysvp1 = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+       &
            a8i*dt)))))))
       polysvp1 = polysvp1*100._rtype

       !       Goff-Gratch formulation:
       !        POLYSVP1 = 10.**(-9.09718*(273.16/T-1.)-3.56654*                 &
       !          log10(273.16/T)+0.876793*(1.-T/273.16)+                        &
       !          log10(6.1071))*100.


    elseif (i_type.eq.0 .or. t.ge.zerodegc) then
       ! LIQUID

       !       Flatau formulation:
       dt       = max(-80._rtype,t-273.15_rtype)
       polysvp1 = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
       polysvp1 = polysvp1*100._rtype

       !       Goff-Gratch formulation:
       !        POLYSVP1 = 10.**(-7.90298*(373.16/T-1.)+                         &
       !             5.02808*log10(373.16/T)-                                    &
       !             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+                  &
       !             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+                &
       !             log10(1013.246))*100.

    !PMC added error checking
    else

       call report_error_info('Something went wrong', 'polysvp1')
       write(err_msg,*)'** polysvp1 i_type must be 0 or 1 but is: ', &
            i_type,' temperature is:',t,' in file: ',__FILE__,' at line:',__LINE__

       call endscreamrun(err_msg)
    endif

   return

  end function polysvp1

  subroutine check_temp(t, subname)
    !Check if temprature values are in legit range
    use scream_abortutils, only : endscreamrun
    use ieee_arithmetic,   only : ieee_is_finite, ieee_is_nan

    implicit none

    real(rtype),      intent(in) :: t
    character(len=*), intent(in) :: subname

    character(len=1000) :: err_msg

    if(t <= 0.0_rtype) then
       write(err_msg,*)'Error: Called from:',trim(adjustl(subname)),'; Temperature is:',t,' which is <= 0._r8 in file:',__FILE__, &
            ' at line:',__LINE__
       call endscreamrun(err_msg)
    elseif(.not. ieee_is_finite(t)) then
       write(err_msg,*)'Error: Called from:',trim(adjustl(subname)),'; Temperature is:',t,' which is not finite in file:', &
            __FILE__,' at line:',__LINE__
       call endscreamrun(err_msg)
    elseif(ieee_is_nan(t)) then
       write(err_msg,*)'Error: Called from:',trim(adjustl(subname)),'; Temperature is:',t,' which is NaN in file:',__FILE__, &
             'at line:',__LINE__
       call endscreamrun(err_msg)
    endif

    return
  end subroutine check_temp

  !------------------------------------------------------------------------------------------!

  !======================================================================================!

  subroutine find_lookupTable_indices_1a(dumi,dumjj,dumii,dumzz,dum1,dum4,dum5,dum6,      &
       isize,rimsize,densize,qitot,nitot,qirim,   &
       rhop)

    !------------------------------------------------------------------------------------------!
    ! Finds indices in 3D ice (only) lookup table
    !------------------------------------------------------------------------------------------!

    implicit none

    ! arguments:
    integer, intent(out) :: dumi,dumjj,dumii,dumzz
    real(rtype),    intent(out) :: dum1,dum4,dum5,dum6
    integer, intent(in)  :: isize,rimsize,densize
    real(rtype),    intent(in)  :: qitot,nitot,qirim,rhop

    !------------------------------------------------------------------------------------------!
    ! find index for qi (normalized ice mass mixing ratio = qitot/nitot)
    !             dum1 = (log10(qitot)+16.)/0.70757  !orig
    !             dum1 = (log10(qitot)+16.)*1.41328
    ! we are inverting this equation from the lookup table to solve for i:
    ! qitot/nitot=261.7**((i+10)*0.1)*1.e-18
    dum1 = (bfb_log10(qitot/nitot)+18._rtype)*lookup_table_1a_dum1_c-10._rtype ! For computational efficiency
    dumi = int(dum1)
    ! set limits (to make sure the calculated index doesn't exceed range of lookup table)
    dum1 = min(dum1,real(isize,rtype))
    dum1 = max(dum1,1._rtype)
    dumi = max(1,dumi)
    dumi = min(isize-1,dumi)

    ! find index for rime mass fraction
    dum4  = (qirim/qitot)*3._rtype + 1._rtype
    dumii = int(dum4)
    ! set limits
    dum4  = min(dum4,real(rimsize,rtype))
    dum4  = max(dum4,1._rtype)
    dumii = max(1,dumii)
    dumii = min(rimsize-1,dumii)

    ! find index for bulk rime density
    ! (account for uneven spacing in lookup table for density)
    if (rhop.le.650._rtype) then
       dum5 = (rhop-50._rtype)*0.005_rtype + 1._rtype
    else
       dum5 =(rhop-650._rtype)*0.004_rtype + 4._rtype
    endif
    dumjj = int(dum5)
    ! set limits
    dum5  = min(dum5,real(densize,rtype))
    dum5  = max(dum5,1._rtype)
    dumjj = max(1,dumjj)
    dumjj = min(densize-1,dumjj)

    dum6  = -99
    dumzz = -99

    return
  end subroutine find_lookupTable_indices_1a

  !======================================================================================!

  subroutine find_lookupTable_indices_1b(dumj,dum3,rcollsize,qr,nr)

    !------------------------------------------------------------------------------------------!
    ! Finds indices in 3D rain lookup table
    !------------------------------------------------------------------------------------------!

    implicit none

    ! arguments:
    integer, intent(out) :: dumj
    real(rtype),    intent(out) :: dum3
    integer, intent(in)  :: rcollsize
    real(rtype),    intent(in)  :: qr,nr

    ! local variables:
    real(rtype)                 :: dumlr
    real(rtype)                 :: real_rcollsize

    !------------------------------------------------------------------------------------------!
    real_rcollsize = real(rcollsize)
    ! find index for scaled mean rain size
    ! if no rain, then just choose dumj = 1 and do not calculate rain-ice collection processes
    if (qr.ge.qsmall .and. nr.gt.0._rtype) then
       ! calculate scaled mean size for consistency with ice lookup table
       dumlr = bfb_cbrt(qr/(pi*rhow*nr))
       dum3  = (bfb_log10(1._rtype*dumlr)+5._rtype)*10.70415_rtype
       dumj  = int(dum3)
       ! set limits
       dum3  = min(dum3,real_rcollsize)
       dum3  = max(dum3,1._rtype)
       dumj  = max(1,dumj)
       dumj  = min(rcollsize-1,dumj)
    else
       dumj  = 1
       dum3  = 1._rtype
    endif

    return

  end subroutine find_lookupTable_indices_1b

  !PMC removed find_lookupTable_indices_2 because it was used for multi-category

  !======================================================================================!
  subroutine find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3,mu_r,lamr)

    !------------------------------------------------------------------------------------------!
    ! Finds indices in rain lookup table (3)
    !------------------------------------------------------------------------------------------!

    implicit none

    ! arguments:
    integer, intent(out) :: dumii,dumjj
    real(rtype),    intent(out) :: dum1,rdumii,rdumjj,inv_dum3
    real(rtype),    intent(in)  :: mu_r,lamr

    !------------------------------------------------------------------------------------------!

    ! find location in scaled mean size space
    dum1 = (mu_r+1._rtype)/lamr
    if (dum1.le.195.e-6_rtype) then
       inv_dum3  = 0.1_rtype
       rdumii = (dum1*1.e6_rtype+5._rtype)*inv_dum3
       rdumii = max(rdumii, 1._rtype)
       rdumii = min(rdumii,20._rtype)
       dumii  = int(rdumii)
       dumii  = max(dumii, 1)
       dumii  = min(dumii,20)
    elseif (dum1.gt.195.e-6_rtype) then
       inv_dum3  = thrd*0.1_rtype            !i.e. 1/30
       rdumii = (dum1*1.e+6_rtype-195._rtype)*inv_dum3 + 20._rtype
       rdumii = max(rdumii, 20._rtype)
       rdumii = min(rdumii,300._rtype)
       dumii  = int(rdumii)
       dumii  = max(dumii, 20)
       dumii  = min(dumii,299)
    endif

    ! find location in mu_r space
    rdumjj = mu_r+1._rtype
    rdumjj = max(rdumjj,1._rtype)
    rdumjj = min(rdumjj,10._rtype)
    dumjj  = int(rdumjj)
    dumjj  = max(dumjj,1)
    dumjj  = min(dumjj,9)

   return

  end subroutine find_lookupTable_indices_3


  !===========================================================================================
  subroutine get_cloud_dsd2(qc,nc,mu_c,rho,nu,dnu,lamc,cdist,cdist1,lcldm)

    implicit none

    !arguments:
    real(rtype), dimension(:), intent(in)  :: dnu
    real(rtype),     intent(in)            :: qc,rho,lcldm
    real(rtype),     intent(inout)         :: nc
    real(rtype),     intent(out)           :: mu_c,nu,lamc,cdist,cdist1

    !local variables
    real(rtype)                            :: lammin,lammax
    integer                         :: dumi

    !--------------------------------------------------------------------------

    nu = 0.0_rtype

    if (qc.ge.qsmall) then

       ! set minimum nc to prevent floating point error
       nc   = max(nc,nsmall)
       mu_c = 0.0005714_rtype*(nc*1.e-6_rtype*rho)+0.2714_rtype
       mu_c = 1._rtype/(mu_c*mu_c)-1._rtype
       mu_c = max(mu_c,2._rtype)
       mu_c = min(mu_c,15._rtype)

       ! interpolate for mass distribution spectral shape parameter (for SB warm processes)
       if (iparam.eq.1) then
          dumi = int(mu_c)
          nu   = dnu(dumi)+(dnu(dumi+1)-dnu(dumi))*(mu_c-dumi)
       endif

       ! calculate lamc
       lamc = bfb_cbrt(cons1*nc*(mu_c+3._rtype)*(mu_c+2._rtype)*(mu_c+1._rtype)/qc)

       ! apply lambda limiters
       lammin = (mu_c+1._rtype)*2.5e+4_rtype   ! min: 40 micron mean diameter
       lammax = (mu_c+1._rtype)*1.e+6_rtype    ! max:  1 micron mean diameter

       if (lamc.lt.lammin) then
          lamc = lammin
          nc   = 6._rtype*(lamc*lamc*lamc)*qc/(pi*rhow*(mu_c+3._rtype)*(mu_c+2._rtype)*(mu_c+1._rtype))
       elseif (lamc.gt.lammax) then
          lamc = lammax
          nc   = 6._rtype*(lamc*lamc*lamc)*qc/(pi*rhow*(mu_c+3._rtype)*(mu_c+2._rtype)*(mu_c+1._rtype))
       endif

       cdist  = nc*(mu_c+1._rtype)/lamc
       cdist1 = nc*lcldm/bfb_gamma(mu_c+1._rtype)

    else

       lamc   = 0._rtype
       cdist  = 0._rtype
       cdist1 = 0._rtype

    endif

   return

  end subroutine get_cloud_dsd2


  !===========================================================================================
  subroutine get_rain_dsd2(qr,nr,mu_r,lamr,cdistr,logn0r,rcldm)

    ! Computes and returns rain size distribution parameters

    implicit none

    !arguments:
    real(rtype),     intent(in)            :: qr,rcldm
    real(rtype),     intent(inout)         :: nr
    real(rtype),     intent(out)           :: lamr,mu_r,cdistr,logn0r

    !local variables:
    real(rtype)                            :: inv_dum,lammax,lammin

    !--------------------------------------------------------------------------

    if (qr.ge.qsmall) then

       ! use lookup table to get mu
       ! mu-lambda relationship is from Cao et al. (2008), eq. (7)

       ! find spot in lookup table
       ! (scaled N/q for lookup table parameter space_
       nr      = max(nr,nsmall)
       inv_dum = bfb_cbrt(qr/(cons1*nr*6._rtype))

       ! Apply constant mu_r:  Recall the switch to v4 tables means constant mu_r
       mu_r = mu_r_constant
       lamr   = bfb_cbrt(cons1*nr*(mu_r+3._rtype)*(mu_r+2._rtype)*(mu_r+1._rtype)/(qr))  ! recalculate slope based on mu_r
       lammax = (mu_r+1._rtype)*1.e+5_rtype   ! check for slope
       lammin = (mu_r+1._rtype)*1250._rtype   ! set to small value since breakup is explicitly included (mean size 0.8 mm)

       ! apply lambda limiters for rain
       if (lamr.lt.lammin) then
          lamr = lammin
          nr   = bfb_exp(3._rtype*bfb_log(lamr)+bfb_log(qr)+bfb_log(bfb_gamma(mu_r+1._rtype))-bfb_log(bfb_gamma(mu_r+4._rtype)))/(cons1)
       elseif (lamr.gt.lammax) then
          lamr = lammax
          nr   = bfb_exp(3._rtype*bfb_log(lamr)+bfb_log(qr)+bfb_log(bfb_gamma(mu_r+1._rtype))-bfb_log(bfb_gamma(mu_r+4._rtype)))/(cons1)
       endif

       cdistr  = nr*rcldm/bfb_gamma(mu_r+1._rtype)
       logn0r  = bfb_log10(nr)+(mu_r+1._rtype)*bfb_log10(lamr)-bfb_log10(bfb_gamma(mu_r+1._rtype)) !note: logn0r is calculated as log10(n0r)

    else

       lamr   = 0._rtype
       cdistr = 0._rtype
       logn0r = 0._rtype

    endif

   return

  end subroutine get_rain_dsd2


  !===========================================================================================
  subroutine calc_bulkRhoRime(qi_tot,qi_rim,bi_rim,rho_rime)

    !--------------------------------------------------------------------------------
    !  Calculates and returns the bulk rime density from the prognostic ice variables
    !  and adjusts qirim and birim appropriately.
    !--------------------------------------------------------------------------------

    implicit none

    !arguments:
    real(rtype), intent(in)    :: qi_tot
    real(rtype), intent(inout) :: qi_rim,bi_rim
    real(rtype), intent(out)   :: rho_rime

    !--------------------------------------------------------------------------

    if (bi_rim.ge.1.e-15_rtype) then
       rho_rime = qi_rim/bi_rim
       !impose limits on rho_rime;  adjust bi_rim if needed
       if (rho_rime.lt.rho_rimeMin) then
          rho_rime = rho_rimeMin
          bi_rim   = qi_rim/rho_rime
       elseif (rho_rime.gt.rho_rimeMax) then
          rho_rime = rho_rimeMax
          bi_rim   = qi_rim/rho_rime
       endif
    else
       qi_rim   = 0._rtype
       bi_rim   = 0._rtype
       rho_rime = 0._rtype
    endif

    !set upper constraint qi_rim <= qi_tot
    if (qi_rim.gt.qi_tot .and. rho_rime.gt.0._rtype) then
       qi_rim = qi_tot
       bi_rim = qi_rim/rho_rime
    endif

    !impose consistency
    if (qi_rim.lt.qsmall) then
       qi_rim = 0._rtype
       bi_rim = 0._rtype
    endif

   return

  end subroutine calc_bulkRhoRime


  !===========================================================================================
  subroutine impose_max_total_Ni(nitot_local,max_total_Ni,inv_rho_local)

    !--------------------------------------------------------------------------------
    ! Impose maximum total ice number concentration (total of all ice categories).
    ! If the sum of all nitot(:) exceeds maximum allowable, each category to preserve
    ! ratio of number between categories.
    !--------------------------------------------------------------------------------

    implicit none

    !arguments:
    real(rtype), intent(inout)               :: nitot_local      !PMC - scalar now that nCat deleted.
    real(rtype), intent(in)                  :: max_total_Ni,inv_rho_local

    !local variables:
    real(rtype)                              :: dum

    if (nitot_local.ge.1.e-20_rtype) then
       dum = max_total_Ni*inv_rho_local/nitot_local
       nitot_local = nitot_local*min(dum,1._rtype)
    endif

  end subroutine impose_max_total_Ni


  !===========================================================================================

  real(rtype) function qv_sat(t_atm,p_atm,i_wrt)

    !------------------------------------------------------------------------------------
    ! Calls polysvp1 to obtain the saturation vapor pressure, and then computes
    ! and returns the saturation mixing ratio, with respect to either liquid or ice,
    ! depending on value of 'i_wrt'
    !------------------------------------------------------------------------------------

    implicit none

    !Calling parameters:
    real(rtype)    :: t_atm  !temperature [K]
    real(rtype)    :: p_atm  !pressure    [Pa]
    integer :: i_wrt  !index, 0 = w.r.t. liquid, 1 = w.r.t. ice

    !Local variables:
    real(rtype)            :: e_pres         !saturation vapor pressure [Pa]

    !------------------
    !Check if temprature is within legitimate range
    call check_temp(t_atm, "qv_sat")

    !e_pres = polysvp1(t_atm,i_wrt)
    e_pres = MurphyKoop_svp(t_atm,i_wrt)
    qv_sat = ep_2*e_pres/max(1.e-3_rtype,(p_atm-e_pres))

    return

  end function qv_sat

  !===========================================================================================

  subroutine check_values(Qv,T,kts,kte,timestepcount,force_abort,source_ind,col_loc)

    !------------------------------------------------------------------------------------
    ! Checks current values of prognotic variables for reasonable values and
    ! stops and prints values if they are out of specified allowable ranges.
    !
    ! 'check_consistency' means include trap for inconsistency in moments;
    ! otherwise, only trap for Q, T, and negative Qx, etc.  This option is here
    ! to allow for Q<qsmall.and.N>nsmall or Q>qsmall.and.N<small which can be produced
    ! at the leading edges due to sedimentation and whose values are accpetable
    ! since lambda limiters are later imposed after SEDI (so one does not necessarily
    ! want to trap for inconsistency after sedimentation has been called).
    !
    ! The value 'source_ind' indicates the approximate location in 'p3_main'
    ! from where 'check_values' was called before it resulted in a trap.
    !
    !------------------------------------------------------------------------------------

    use scream_abortutils, only : endscreamrun

    implicit none

    !Calling parameters:
    real(rtype), dimension(kts:kte), intent(in) :: Qv, T
    real(rtype), dimension(3), intent(in) :: col_loc
    integer,                intent(in) :: source_ind,timestepcount,kts,kte
    logical(btype),                intent(in) :: force_abort         !.TRUE. = forces abort if value violation is detected

    !Local variables:
    real(rtype), parameter :: T_low  = 160._rtype !173._rtype
    real(rtype), parameter :: T_high = 355._rtype !323._rtype
    real(rtype), parameter :: Q_high = 40.e-3_rtype
    real(rtype), parameter :: N_high = 1.e+20_rtype
    real(rtype), parameter :: B_high = Q_high*1.e-3_rtype
    real(rtype), parameter :: x_high = 1.e+30_rtype
    real(rtype), parameter :: x_low  = 0._rtype
    integer         :: k
    logical(btype)         :: trap,badvalue_found
    character(len=1000)    :: err_msg

    trap = .false.

    k_loop: do k = kts, kte

       ! check unrealistic values or NANs for T and Qv
       if (.not.(T(k)>T_low .and. T(k)<T_high)) then
          write(iulog,'(a60,i5,a2,i8,a2,f8.4,a2,f8.4,a2,i4,a2,i8,a2,e16.8)') &
             '** WARNING IN P3_MAIN -- src, gcol, lon, lat, lvl, tstep, T:',source_ind,', ',int(col_loc(1)),', ',col_loc(2),', ',col_loc(3),', ',k,', ',timestepcount,', ',T(k)
          trap = .true.
       endif
       if (.not.(Qv(k)>=0. .and. Qv(k)<Q_high)) then
          write(iulog,'(a60,i5,a2,i8,a2,f8.4,a2,f8.4,a2,i4,a2,i8,a2,e16.8)') &
             '** WARNING IN P3_MAIN -- src, gcol, lon, lat, lvl, tstep, Qv:',source_ind,', ',int(col_loc(1)),', ',col_loc(2),', ',col_loc(3),', ',k,', ',timestepcount,', ',Qv(k)
          !trap = .true.  !note, tentatively no trap, since Qv could be negative passed in to mp
       endif

       ! check NANs for mp variables:
       badvalue_found = .false.

    enddo k_loop

    if (trap .and. force_abort) then
       print*
       print*,'** DEBUG TRAP IN P3_MAIN, s/r CHECK_VALUES -- source: ',source_ind
       print*
       if (source_ind/=100) then
          write(err_msg,*)'Source_ind should be 100, source_ind is:', &
               source_ind,' in file:',__FILE__,' at line:',__LINE__
          call endscreamrun(err_msg)
       endif
    endif

   return

  end subroutine check_values

  subroutine ice_cldliq_collection(rho,t,rhofaci,    &
  f1pr04,qitot_incld,qc_incld,nitot_incld,nc_incld,    &
             qccol,nccol,qcshd,ncshdc)

   !.......................
   ! collection of droplets

   ! here we multiply rates by air density, air density fallspeed correction
   ! factor, and collection efficiency since these parameters are not
   ! included in lookup table calculations
   ! for T < 273.15, assume collected cloud water is instantly frozen
   ! note 'f1pr' values are normalized, so we need to multiply by N

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: t
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: f1pr04  ! collection of cloud water by ice
   real(rtype), intent(in) :: qitot_incld
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: nitot_incld
   real(rtype), intent(in) :: nc_incld


   real(rtype), intent(out) :: qccol
   real(rtype), intent(out) :: nccol
   real(rtype), intent(out) :: qcshd
   real(rtype), intent(out) :: ncshdc

   if (qitot_incld .ge.qsmall .and. qc_incld .ge.qsmall) then
      if  (t .le.zerodegc) then
         qccol = rhofaci*f1pr04*qc_incld*eci*rho*nitot_incld
         nccol = rhofaci*f1pr04*nc_incld*eci*rho*nitot_incld
      else if (t .gt. zerodegc) then
         ! for T > 273.15, assume cloud water is collected and shed as rain drops
         ! sink for cloud water mass and number, note qcshed is source for rain mass
         qcshd = rhofaci*f1pr04*qc_incld*eci*rho*nitot_incld
         nccol = rhofaci*f1pr04*nc_incld*eci*rho*nitot_incld
         ! source for rain number, assume 1 mm drops are shed
         ncshdc = qcshd*inv_dropmass
      end if
   end if

   return

  end subroutine ice_cldliq_collection


  subroutine ice_rain_collection(rho,t,rhofaci,    &
  logn0r,f1pr07,f1pr08,qitot_incld,nitot_incld,qr_incld,    &
  qrcol, nrcol)

   !....................
   ! collection of rain

   ! here we multiply rates by air density, air density fallspeed correction
   ! factor, collection efficiency, and n0r since these parameters are not
   ! included in lookup table calculations

   ! for T < 273.15, assume all collected rain mass freezes
   ! note this is a sink for rain mass and number and a source
   ! for ice mass

   ! note 'f1pr' values are normalized, so we need to multiply by N

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: t
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: logn0r
   real(rtype), intent(in) :: f1pr07 !collection of rain number by ice
   real(rtype), intent(in) :: f1pr08 !collection of rain mass by ice
   real(rtype), intent(in) :: qitot_incld
   real(rtype), intent(in) :: nitot_incld
   real(rtype), intent(in) :: qr_incld

   real(rtype), intent(out) :: qrcol
   real(rtype), intent(out) :: nrcol

   if (qitot_incld.ge.qsmall .and. qr_incld.ge.qsmall) then
      if (t.le.zerodegc) then
         ! note: f1pr08 and logn0r are already calculated as log_10
         qrcol = bfb_pow(10._rtype,(f1pr08+logn0r))*rho*rhofaci*eri*nitot_incld
         nrcol = bfb_pow(10._rtype,(f1pr07+logn0r))*rho*rhofaci*eri*nitot_incld
      else if (t .gt. zerodegc) then
         ! rain number sink due to collection
         ! for T > 273.15, assume collected rain number is shed as
         ! 1 mm drops
         ! note that melting of ice number is scaled to the loss
         ! rate of ice mass due to melting
         ! collection of rain above freezing does not impact total rain mass
         nrcol  = bfb_pow(10._rtype,(f1pr07 + logn0r))*rho*rhofaci*eri*nitot_incld
         ! for now neglect shedding of ice collecting rain above freezing, since snow is
         ! not expected to shed in these conditions (though more hevaily rimed ice would be
         ! expected to lead to shedding)
      endif
   endif

   return

  end subroutine ice_rain_collection

  subroutine ice_self_collection(rho,rhofaci,    &
  f1pr03,eii,qirim_incld,qitot_incld,nitot_incld,    &
             nislf)

   ! self-collection of ice

   ! here we multiply rates by collection efficiency, air density,
   ! and air density correction factor since these are not included
   ! in the lookup table calculations
   ! note 'f1pr' values are normalized, so we need to multiply by N

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: f1pr03 ! ice collection within a category
   real(rtype), intent(in) :: eii
   real(rtype), intent(in) :: qirim_incld
   real(rtype), intent(in) :: qitot_incld
   real(rtype), intent(in) :: nitot_incld

   real(rtype), intent(out) :: nislf

   real(rtype) :: tmp1, Eii_fact

   if (qitot_incld.ge.qsmall) then
      ! Determine additional collection efficiency factor to be applied to ice-ice collection.
      ! The computed values of qicol and nicol are multipiled by Eii_fact to gradually shut off collection
      ! if ice is highly rimed.
      if (qirim_incld>0._rtype) then
         tmp1 = qirim_incld/qitot_incld   !rime mass fraction
         if (tmp1.lt.0.6_rtype) then
            Eii_fact=1._rtype
         else if (tmp1.ge.0.6_rtype.and.tmp1.lt.0.9_rtype) then
            ! linear ramp from 1 to 0 for Fr between 0.6 and 0.9
            Eii_fact = 1._rtype-(tmp1-0.6_rtype)/0.3_rtype
         else
            Eii_fact = 0._rtype
         endif
      else
         Eii_fact = 1._rtype
      endif

      nislf = f1pr03*rho*eii*Eii_fact*rhofaci*nitot_incld
   endif

   return

end subroutine ice_self_collection

!PMC note - indentation pattern changes here.

subroutine ice_melting(rho,t,pres,rhofaci,    &
f1pr05,f1pr14,xxlv,xlf,dv,sc,mu,kap,qv,qitot_incld,nitot_incld,    &
           qimlt,nimlt)
   ! melting
   ! need to add back accelerated melting due to collection of ice mass by rain (pracsw1)
   ! note 'f1pr' values are normalized, so we need to multiply by N
   ! currently enhanced melting from collision is neglected
   ! include RH dependence

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: t
   real(rtype), intent(in) :: pres
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: f1pr05 ! melting
   real(rtype), intent(in) :: f1pr14 ! melting (ventilation term)
   real(rtype), intent(in) :: xxlv
   real(rtype), intent(in) :: xlf
   real(rtype), intent(in) :: dv
   real(rtype), intent(in) :: sc
   real(rtype), intent(in) :: mu
   real(rtype), intent(in) :: kap
   real(rtype), intent(in) :: qv
   real(rtype), intent(in) :: qitot_incld
   real(rtype), intent(in) :: nitot_incld

   real(rtype), intent(out) :: qimlt
   real(rtype), intent(out) :: nimlt

   real(rtype) :: qsat0

   if (qitot_incld .ge.qsmall .and. t.gt.zerodegc) then
      qsat0 = qv_sat( zerodegc,pres,0 )

      qimlt = ((f1pr05+f1pr14*bfb_cbrt(sc)*bfb_sqrt(rhofaci*rho/mu))*((t-   &
      zerodegc)*kap-rho*xxlv*dv*(qsat0-qv))*2._rtype*pi/xlf)*nitot_incld

      qimlt = max(qimlt,0._rtype)
      nimlt = qimlt*(nitot_incld/qitot_incld)

   endif

   return

end subroutine ice_melting


subroutine ice_cldliq_wet_growth(rho,t,pres,rhofaci,    &
f1pr05,f1pr14,xxlv,xlf,dv,kap,mu,sc,    &
qv,qc_incld,qitot_incld,nitot_incld,qr_incld,    &
           log_wetgrowth,qrcol,qccol,qwgrth,nrshdr,qcshd)

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: t
   real(rtype), intent(in) :: pres
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: f1pr05 ! melting
   real(rtype), intent(in) :: f1pr14 ! melting (ventilation term)
   real(rtype), intent(in) :: xxlv
   real(rtype), intent(in) :: xlf
   real(rtype), intent(in) :: dv
   real(rtype), intent(in) :: kap
   real(rtype), intent(in) :: mu
   real(rtype), intent(in) :: sc
   real(rtype), intent(in) :: qv
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: qitot_incld
   real(rtype), intent(in) :: nitot_incld
   real(rtype), intent(in) :: qr_incld

   logical(btype), intent(inout) :: log_wetgrowth
   real(rtype), intent(inout) :: qrcol
   real(rtype), intent(inout) :: qccol
   real(rtype), intent(inout) :: qwgrth
   real(rtype), intent(inout) :: nrshdr
   real(rtype), intent(inout) :: qcshd

   real(rtype) :: qsat0, dum, dum1

   if (qitot_incld.ge.qsmall .and. qc_incld+qr_incld.ge.1.e-6_rtype .and. t.lt.zerodegc) then
      qsat0=qv_sat( zerodegc,pres,0 )

      qwgrth = ((f1pr05 + f1pr14*bfb_cbrt(sc)*bfb_sqrt(rhofaci*rho/mu))*       &
      2._rtype*pi*(rho*xxlv*dv*(qsat0-qv)-(t-zerodegc)*           &
      kap)/(xlf+cpw*(t-zerodegc)))*nitot_incld

      qwgrth = max(qwgrth,0._rtype)
      dum    = max(0._rtype,(qccol+qrcol)-qwgrth)
      if (dum.ge.1.e-10_rtype) then
         nrshdr = nrshdr + dum*1.923e+6_rtype   ! 1/5.2e-7, 5.2e-7 is the mass of a 1 mm raindrop
         if ((qccol+qrcol).ge.1.e-10_rtype) then
            dum1  = 1._rtype/(qccol+qrcol)
            qcshd = qcshd + dum*qccol*dum1
            qccol = qccol - dum*qccol*dum1
            qrcol = qrcol - dum*qrcol*dum1
         endif
         ! densify due to wet growth
         log_wetgrowth = .true.
      endif


   end if

   return

end subroutine ice_cldliq_wet_growth


subroutine calc_ice_relaxation_timescale(rho,t,rhofaci,     &
f1pr05,f1pr14,dv,mu,sc,qitot_incld,nitot_incld,    &
epsi,epsi_tot)

   !-----------------------------
   ! calcualte total inverse ice relaxation timescale combined for all ice categories
   ! note 'f1pr' values are normalized, so we need to multiply by N
   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: t
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: f1pr05 ! melting
   real(rtype), intent(in) :: f1pr14 ! melting (ventilation term)
   real(rtype), intent(in) :: dv
   real(rtype), intent(in) :: mu
   real(rtype), intent(in) :: sc
   real(rtype), intent(in) :: qitot_incld
   real(rtype), intent(in) :: nitot_incld

   real(rtype), intent(out) :: epsi
   real(rtype), intent(inout) :: epsi_tot

   if (qitot_incld.ge.qsmall .and. t.lt.zerodegc) then
      epsi = ((f1pr05+f1pr14*bfb_cbrt(sc)*bfb_sqrt(rhofaci*rho/mu))*2._rtype*pi* &
      rho*dv)*nitot_incld
      epsi_tot   = epsi_tot + epsi
   else
      epsi = 0._rtype
   endif

   return

end subroutine calc_ice_relaxation_timescale


subroutine calc_liq_relaxation_timescale(rho,f1r,f2r,     &
dv,mu,sc,mu_r,lamr,cdistr,cdist,qr_incld,qc_incld, &
epsr,epsc)

   implicit none

   real(rtype), intent(in)  :: rho
   real(rtype), intent(in)  :: f1r
   real(rtype), intent(in)  :: f2r
   real(rtype), intent(in)  :: dv
   real(rtype), intent(in)  :: mu
   real(rtype), intent(in)  :: sc
   real(rtype), intent(in)  :: mu_r
   real(rtype), intent(in)  :: lamr
   real(rtype), intent(in)  :: cdistr
   real(rtype), intent(in)  :: cdist
   real(rtype), intent(in)  :: qr_incld
   real(rtype), intent(in)  :: qc_incld
   real(rtype), intent(out) :: epsr
   real(rtype), intent(out) :: epsc

   integer     :: dumii, dumjj
   real(rtype) :: rdumii, rdumjj
   real(rtype) :: dum, dum1, dum2, inv_dum3

   if (qr_incld.ge.qsmall) then
      call find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3,mu_r,lamr)
      !interpolate value at mu_r
      dum1 = revap_table(dumii,dumjj)+(rdumii-real(dumii))*                            &
             (revap_table(dumii+1,dumjj)-revap_table(dumii,dumjj))

      !interoplate value at mu_r+1
      dum2 = revap_table(dumii,dumjj+1)+(rdumii-real(dumii))*                          &
             (revap_table(dumii+1,dumjj+1)-revap_table(dumii,dumjj+1))
      !final interpolation
      dum  = dum1+(rdumjj-real(dumjj))*(dum2-dum1)

      epsr = 2._rtype*pi*cdistr*rho*dv*                                                &
             (f1r*bfb_gamma(mu_r+2._rtype)/lamr +                                    &
              f2r*bfb_sqrt(rho/mu)*bfb_cbrt(sc)*dum)
   else
      epsr = 0._rtype
   endif

   if (qc_incld.ge.qsmall) then
      epsc = 2._rtype*pi*rho*dv*cdist
   else
      epsc = 0._rtype
   endif

   return

end subroutine calc_liq_relaxation_timescale


subroutine calc_rime_density(t,rhofaci,    &
f1pr02,acn,lamc, mu_c,qc_incld,qccol,    &
           vtrmi1,rhorime_c)

   !.........................
   ! calculate rime density

   !     FUTURE:  Add source term for birim (=qccol/rhorime_c) so that all process rates calculations
   !              are done together, before conservation.

   ! NOTE: Tc (ambient) is assumed for the surface temperature.  Technically,
   ! we should diagose graupel surface temperature from heat balance equation.
   ! (but the ambient temperature is a reasonable approximation; tests show
   ! very little sensitivity to different assumed values, Milbrandt and Morrison 2012).

   ! Compute rime density: (based on parameterization of Cober and List, 1993 [JAS])
   ! for simplicty use mass-weighted ice and droplet/rain fallspeeds

   implicit none

   real(rtype), intent(in) :: t
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: f1pr02 !mass-weighted fallspeed
   real(rtype), intent(in) :: acn
   real(rtype), intent(in) :: lamc
   real(rtype), intent(in) :: mu_c
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: qccol

   real(rtype), intent(out) :: vtrmi1
   real(rtype), intent(out) :: rhorime_c

   real(rtype) :: iTc = 0.0_rtype
   real(rtype) :: Vt_qc = 0.0_rtype
   real(rtype) :: D_c  = 0.0_rtype
   real(rtype) :: V_impact = 0.0_rtype
   real(rtype) :: Ri = 0.0_rtype

   ! if (qitot_incld(i,k).ge.qsmall .and. t(i,k).lt.zerodegc) then
   !  NOTE:  condition applicable for cloud only; modify when rain is added back
   if (qccol.ge.qsmall .and. t.lt.zerodegc) then
      ! get mass-weighted mean ice fallspeed
      vtrmi1 = f1pr02*rhofaci
      iTc   = 1._rtype/min(-0.001_rtype,t-zerodegc)

             ! cloud:
      if (qc_incld.ge.qsmall) then
         ! droplet fall speed
         ! (use Stokes' formulation (thus use analytic solution)
         Vt_qc = acn*bfb_gamma(4._rtype+bcn+mu_c)/(bfb_pow(lamc,bcn)*bfb_gamma(mu_c+4._rtype))
         ! use mass-weighted mean size
         D_c = (mu_c+4._rtype)/lamc
         V_impact  = abs(vtrmi1-Vt_qc)
         Ri        = -0.5e+6_rtype * D_c * V_impact * iTc
         !               Ri        = max(1.,min(Ri,8.))
         Ri        = max(1._rtype,min(Ri,12._rtype))
         if (Ri.le.8.) then
            rhorime_c  = (0.051_rtype + 0.114_rtype*Ri - 0.0055_rtype*bfb_square(Ri))*1000._rtype
         else
            ! for Ri > 8 assume a linear fit between 8 and 12,
            ! rhorime = 900 kg m-3 at Ri = 12
            ! this is somewhat ad-hoc but allows a smoother transition
            ! in rime density up to wet growth
            rhorime_c  = 611._rtype+72.25_rtype*(Ri-8._rtype)
         endif
      else
         rhorime_c = 400._rtype
      endif    !if qc>qsmall
   else
      vtrmi1 = 0._rtype ! no velocity if no ice
      rhorime_c = 400._rtype
   endif ! qi > qsmall and T < 273.15

   return

end subroutine calc_rime_density

function subgrid_variance_scaling(relvar, expon) result(res)
  ! Finds a coefficient for process rates based on the inverse relative variance
  ! of cloud water.

  real(rtype), intent(in) :: relvar
  real(rtype), intent(in) :: expon
  real(rtype) :: res

  res = bfb_gamma(relvar+expon)/(bfb_gamma(relvar)*bfb_pow(relvar,expon))

end function subgrid_variance_scaling

subroutine cldliq_immersion_freezing(t,lamc,mu_c,cdist1,qc_incld,qc_relvar,    &
           qcheti,ncheti)

   !............................................................
   ! contact and immersion freezing droplets

   implicit none

   real(rtype), intent(in) :: t
   real(rtype), intent(in) :: lamc
   real(rtype), intent(in) :: mu_c
   real(rtype), intent(in) :: cdist1
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: qc_relvar

   real(rtype), intent(out) :: qcheti
   real(rtype), intent(out) :: ncheti

   real(rtype) :: dum1, dum2, Q_nuc, N_nuc, sbgrd_var_coef

   if (qc_incld.ge.qsmall .and. t.le.rainfrze) then
      ! for future: calculate gamma(mu_c+4) in one place since its used multiple times  !AaronDonahue, TODO
      dum1 = bfb_exp(aimm*(zerodegc-t))
      dum2 = bfb_cube(1._rtype/lamc)
      sbgrd_var_coef = subgrid_variance_scaling(qc_relvar, 2._rtype)
      Q_nuc = sbgrd_var_coef*cons6*cdist1*bfb_gamma(7._rtype+mu_c)*dum1*bfb_square(dum2)
      N_nuc = cons5*cdist1*bfb_gamma(mu_c+4._rtype)*dum1*dum2
      qcheti = Q_nuc
      ncheti = N_nuc
   endif

   return

end subroutine cldliq_immersion_freezing

subroutine rain_immersion_freezing(t,    &
lamr, mu_r, cdistr, qr_incld,    &
qrheti, nrheti)

   !............................................................
   ! immersion freezing of rain
   ! for future: get rid of log statements below for rain freezing

   implicit none

   real(rtype), intent(in) :: t
   real(rtype), intent(in) :: mu_r
   real(rtype), intent(in) :: lamr
   real(rtype), intent(in) :: cdistr
   real(rtype), intent(in) :: qr_incld

   real(rtype), intent(out) :: qrheti
   real(rtype), intent(out) :: nrheti

   real(rtype) :: Q_nuc, N_nuc

   if (qr_incld.ge.qsmall .and. t.le.rainfrze) then

      Q_nuc = cons6*bfb_exp(bfb_log(cdistr) + bfb_log(bfb_gamma(7._rtype+mu_r)) - 6._rtype*bfb_log(lamr))*bfb_exp(aimm*(zerodegc-t))
      N_nuc = cons5*bfb_exp(bfb_log(cdistr) + bfb_log(bfb_gamma(mu_r+4._rtype)) - 3._rtype*bfb_log(lamr))*bfb_exp(aimm*(zerodegc-t))

      qrheti = Q_nuc
      nrheti = N_nuc

   endif

   return

end subroutine rain_immersion_freezing


subroutine ice_nucleation(t,inv_rho,nitot,naai,supi,odt,log_predictNc,    &
   qinuc,ninuc)

   !................................................................
   ! deposition/condensation-freezing nucleation
   ! allow ice nucleation if < -15 C and > 5% ice supersaturation
   ! use CELL-AVERAGE values, freezing of vapor

   implicit none

   real(rtype), intent(in) :: t
   real(rtype), intent(in) :: inv_rho
   real(rtype), intent(in) :: nitot
   real(rtype), intent(in) :: naai
   real(rtype), intent(in) :: supi
   real(rtype), intent(in) :: odt
   logical(btype), intent(in) :: log_predictNc

   real(rtype), intent(inout) :: qinuc
   real(rtype), intent(inout) :: ninuc


   real(rtype) :: dum, N_nuc, Q_nuc

   if ( t .lt.icenuct .and. supi.ge.0.05_rtype) then
      if(.not. log_predictNc) then
!         ! dum = exp(-0.639+0.1296*100.*supi(i,k))*1000.*inv_rho(i,k)  !Meyers et al. (1992)
         dum = 0.005_rtype*bfb_exp(0.304_rtype*(zerodegc-t))*1000._rtype*inv_rho   !Cooper (1986)
         dum = min(dum,100.e3_rtype*inv_rho)
         N_nuc = max(0._rtype,(dum-nitot)*odt)
         if (N_nuc.ge.1.e-20_rtype) then
            Q_nuc = max(0._rtype,(dum-nitot)*mi0*odt)
            qinuc = Q_nuc
            ninuc = N_nuc
         endif
      else
      ! Ice nucleation predicted by aerosol scheme
         ninuc = max(0._rtype, (naai - nitot)*odt)
         qinuc = ninuc * mi0
      endif
   endif

end subroutine

subroutine droplet_self_collection(rho,inv_rho,qc_incld,mu_c,nu,ncautc,    &
   ncslf)

   !............................
   ! self-collection of droplets

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: inv_rho
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: mu_c
   real(rtype), intent(in) :: nu
   real(rtype), intent(in) :: ncautc

   real(rtype), intent(out) :: ncslf

   if (qc_incld.ge.qsmall) then

      if (iparam.eq.1) then
         !Seifert and Beheng (2001)
         ncslf = -kc*(1.e-3_rtype*rho*qc_incld)**2*(nu+2._rtype)/(nu+1._rtype)*         &
              1.e+6_rtype*inv_rho+ncautc
      elseif (iparam.eq.2) then
         !Beheng (994)
         ncslf = -5.5e+16_rtype*inv_rho*mu_c**(-0.63_rtype)*(1.e-3_rtype*rho*qc_incld)**2
      elseif (iparam.eq.3) then
         !Khroutdinov and Kogan (2000)
         ncslf = 0._rtype
      endif

   endif

end subroutine droplet_self_collection

subroutine cloud_rain_accretion(rho,inv_rho,qc_incld,nc_incld,qr_incld,qc_relvar,    &
   qcacc,ncacc)

  !............................
  ! accretion of cloud by rain

  implicit none

  real(rtype), intent(in) :: rho
  real(rtype), intent(in) :: inv_rho
  real(rtype), intent(in) :: qc_incld
  real(rtype), intent(in) :: nc_incld
  real(rtype), intent(in) :: qr_incld
  real(rtype), intent(in) :: qc_relvar

  real(rtype), intent(out) :: qcacc
  real(rtype), intent(out) :: ncacc

  real(rtype) :: dum, dum1, sbgrd_var_coef

  if (qr_incld.ge.qsmall .and. qc_incld.ge.qsmall) then

     if (iparam.eq.1) then
        !Seifert and Beheng (2001)
        dum   = 1._rtype-qc_incld/(qc_incld+qr_incld)
        dum1  = (dum/(dum+5.e-4_rtype))**4
        qcacc = kr*rho*0.001_rtype*qc_incld*qr_incld*dum1
        ncacc = qcacc*rho*0.001_rtype*(nc_incld*rho*1.e-6_rtype)/(qc_incld*rho*   &
             0.001_rtype)*1.e+6_rtype*inv_rho
     elseif (iparam.eq.2) then
        !Beheng (994)
        qcacc = 6._rtype*rho*(qc_incld*qr_incld)
        ncacc = qcacc*rho*1.e-3_rtype*(nc_incld*rho*1.e-6_rtype)/(qc_incld*rho*1.e-3_rtype)* &
             1.e+6_rtype*inv_rho
     elseif (iparam.eq.3) then
        !Khroutdinov and Kogan (2000)
        !print*,'p3_QcAccret_Expon = ',p3_QcAccret_Expon
        sbgrd_var_coef = subgrid_variance_scaling(qc_relvar, 1.15_rtype ) !p3_QcAccret_Expon
        qcacc = sbgrd_var_coef*67._rtype*bfb_pow(qc_incld*qr_incld, 1.15_rtype) !p3_QcAccret_Expon
        ncacc = qcacc*nc_incld/qc_incld
     endif

     if (qcacc.eq.0._rtype) ncacc = 0._rtype
     if (ncacc.eq.0._rtype) qcacc = 0._rtype

  endif

end subroutine cloud_rain_accretion

subroutine rain_self_collection(rho,qr_incld,nr_incld,    &
   nrslf)

   !.....................................
   ! self-collection and breakup of rain
   ! (breakup following modified Verlinde and Cotton scheme)

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: qr_incld
   real(rtype), intent(in) :: nr_incld
   real(rtype), intent(out) :: nrslf

   real(rtype) :: dum, dum1, dum2

   if (qr_incld.ge.qsmall) then

      ! include breakup
      dum1 = 280.e-6_rtype

      ! use mass-mean diameter (do this by using
      ! the old version of lambda w/o mu dependence)
      ! note there should be a factor of 6^(1/3), but we
      ! want to keep breakup threshold consistent so 'dum'
      ! is expressed in terms of lambda rather than mass-mean D

      dum2 = bfb_cbrt(qr_incld/(pi*rhow*nr_incld))
      if (dum2.lt.dum1) then
         dum = 1._rtype
      else if (dum2.ge.dum1) then
         dum = 2._rtype-bfb_exp(2300._rtype*(dum2-dum1))
      endif

      if (iparam.eq.1) then
         nrslf = dum*kr*1.e-3_rtype*qr_incld*nr_incld*rho
      elseif (iparam.eq.2 .or. iparam.eq.3) then
         nrslf = dum*5.78_rtype*nr_incld*qr_incld*rho
      endif

   endif

end subroutine rain_self_collection


subroutine cloud_water_autoconversion(rho,qc_incld,nc_incld,qc_relvar,    &
   qcaut,ncautc,ncautr)

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: nc_incld
   real(rtype), intent(in) :: qc_relvar

   real(rtype), intent(out) :: qcaut
   real(rtype), intent(out) :: ncautc
   real(rtype), intent(out) :: ncautr

   real(rtype) :: sbgrd_var_coef

   qc_not_small: if (qc_incld.ge.1.e-8_rtype) then

      !Khroutdinov and Kogan (2000)
      !print*,'p3_QcAutoCon_Expon = ',p3_QcAutoCon_Expon
      sbgrd_var_coef = subgrid_variance_scaling(qc_relvar, 2.47_rtype)
      qcaut = sbgrd_var_coef*1350._rtype*bfb_pow(qc_incld,2.47_rtype)*bfb_pow(nc_incld*1.e-6_rtype*rho,-1.79_rtype)
      ! note: ncautr is change in Nr; ncautc is change in Nc
      ncautr = qcaut*cons3
      ncautc = qcaut*nc_incld/qc_incld

      if (qcaut .eq.0._rtype) ncautc = 0._rtype
      if (ncautc.eq.0._rtype) qcaut  = 0._rtype

   endif qc_not_small

end subroutine cloud_water_autoconversion

subroutine back_to_cell_average(lcldm,rcldm,icldm,                         &
   qcacc,qrevp,qcaut,                                                      &
   ncacc,ncslf,ncautc,nrslf,nrevp,ncautr,qisub,nrshdr,qcheti,              &
   qrcol,qcshd,qimlt,qccol,qrheti,nimlt,nccol,ncshdc,ncheti,nrcol,nislf,   &
   qidep,nrheti,nisub,qinuc,ninuc,qiberg)

   ! Here we map the microphysics tendency rates back to CELL-AVERAGE quantities for updating
   ! cell-average quantities.

   implicit none

   ! Intersection of cloud fractions for combination of ice (i), rain (r) and liquid (l)
   real(rtype), intent(in) :: lcldm
   real(rtype), intent(in) :: rcldm
   real(rtype), intent(in) :: icldm

   real(rtype), intent(inout) :: qcacc, qrevp, qcaut, ncacc, ncslf, ncautc, nrslf, nrevp, ncautr
   real(rtype), intent(inout) :: qisub, nrshdr, qcheti, qrcol, qcshd, qimlt, qccol, qrheti, nimlt, nccol, ncshdc, ncheti, nrcol, nislf, qidep
   real(rtype), intent(inout) :: nrheti, nisub, qinuc, ninuc, qiberg

   real(rtype) :: ir_cldm, il_cldm, lr_cldm

   ir_cldm = min(icldm,rcldm)  ! Intersection of ICE and RAIN cloud
   il_cldm = min(icldm,lcldm)  ! Intersection of ICE and LIQUID cloud
   lr_cldm = min(lcldm,rcldm)  ! Intersection of LIQUID and RAIN cloud

   ! Some process rates take place within the intersection of liquid, rain and ice cloud fractions.
   ! We calculate the intersection as the minimum between combinations of cloud fractions and use
   ! these values to map back to cell-average quantities where applicable.

          ! map warm-phase process rates to cell-avg
   qcacc   = qcacc*lr_cldm     ! Accretion of liquid to rain
   qrevp   = qrevp*rcldm       ! Evaporation of rain
   qcaut   = qcaut*lcldm       ! Autoconversion of liquid
   ncacc   = ncacc*lr_cldm     ! Number change due to accretion
   ncslf   = ncslf*lcldm       ! Self collection occurs locally in liq. cloud
   ncautc  = ncautc*lcldm      ! Impact of autoconversion on number
   nrslf   = nrslf*rcldm       ! Self collection occurs locally in rain cloud
   nrevp   = nrevp*rcldm       ! Change in rain number due to evaporation
   ncautr  = ncautr*lr_cldm    ! Autoconversion of rain drops within rain/liq cloud

   ! map ice-phase  process rates to cell-avg
   qisub   = qisub*icldm       ! Sublimation of ice in ice cloud
   nrshdr  = nrshdr*il_cldm    ! Rain # increase due to shedding from rain-ice collisions, occurs when ice and liquid interact
   qcheti  = qcheti*il_cldm    ! Immersion freezing of cloud drops
   qrcol   = qrcol*ir_cldm     ! Collection of rain mass by ice
   qcshd   = qcshd*il_cldm     ! Rain mass growth due to shedding of fain drops after collisions with ice, occurs when ice and liquid interact
   qimlt   = qimlt*icldm       ! Melting of ice
   qccol   = qccol*il_cldm     ! Collection of water by ice
   qrheti  = qrheti*rcldm      ! Immersion freezing of rain
   nimlt   = nimlt*icldm       ! Change in number due to melting
   nccol   = nccol*il_cldm     ! Cloud # change due to collection of cld water by ice
   ncshdc  = ncshdc*il_cldm    ! Number change due to shedding, occurs when ice and liquid interact
   ncheti  = ncheti*lcldm      ! Number change associated with freexzing of cld drops
   nrcol   = nrcol*ir_cldm     ! Rain number change due to collection from ice
   nislf   = nislf*icldm       ! Ice self collection
   qidep   = qidep*icldm       ! Vapor deposition to ice phase
   nrheti  = nrheti*rcldm      ! Change in number due to immersion freezing of rain
   nisub   = nisub*icldm       ! Number change due to sublimation of ice
   qiberg  = qiberg*il_cldm    ! Bergeron process
     ! AaronDonahue: These variables are related to aerosol activation and their usage will be changed in a later PR.
   qinuc   = qinuc             ! Deposition and condensation-freezing nucleation, already cell-averaged
   ninuc   = ninuc             ! Number change due to deposition and condensation-freezing, already cell-averaged

end subroutine back_to_cell_average

subroutine prevent_ice_overdepletion(pres,t,qv,xxls,odt,    &
   qidep,qisub)

   !-- Limit ice process rates to prevent overdepletion of sources such that
   !   the subsequent adjustments are done with maximum possible rates for the
   !   time step.  (note: most ice rates are adjusted here since they must be done
   !   simultaneously (outside of iice-loops) to distribute reduction proportionally
   !   amongst categories.
   !PMC - might need to rethink above statement since only one category now.

   implicit none

   real(rtype), intent(in) :: pres
   real(rtype), intent(in) :: t
   real(rtype), intent(in) :: qv
   real(rtype), intent(in) :: xxls
   real(rtype), intent(in) :: odt

   real(rtype), intent(inout) :: qidep
   real(rtype), intent(inout) :: qisub

   real(rtype) :: dumqvi, qdep_satadj

   dumqvi = qv_sat(t,pres,1)
   qdep_satadj = (qv-dumqvi)/(1._rtype+bfb_square(xxls)*dumqvi/(cp*rv*bfb_square(t)))*odt
   qidep  = qidep*min(1._rtype,max(0._rtype, qdep_satadj)/max(qidep, 1.e-20_rtype))
   qisub  = qisub*min(1._rtype,max(0._rtype,-qdep_satadj)/max(qisub, 1.e-20_rtype))

end subroutine prevent_ice_overdepletion

subroutine cloud_water_conservation(qc,dt,    &
   qcaut,qcacc,qccol,qcheti,qcshd,qiberg,qisub,qidep)

   implicit none

   real(rtype), intent(in) :: qc, dt
   real(rtype), intent(inout) :: qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep
   real(rtype) :: sinks, ratio

   sinks   = (qcaut+qcacc+qccol+qcheti+qcshd+qiberg)*dt
   if (sinks .gt. qc .and. sinks.ge.1.e-20_rtype) then
      ratio  = qc/sinks
      qcaut  = qcaut*ratio
      qcacc  = qcacc*ratio
      qccol  = qccol*ratio
      qcheti = qcheti*ratio
      qcshd  = qcshd*ratio
      qiberg = qiberg*ratio
   else
      ratio = 1.0 ! If not limiting sinks on qc then most likely did not run out of qc
   endif

   !PMC: ratio is also frac of step w/ liq. thus we apply qiberg for
   !"ratio" of timestep and vapor deposition and sublimation  for the
   !remaining frac of the timestep.  Only limit if there will be cloud
   !water to begin with.
   if (qc .gt. 1.e-20_rtype) then
      qidep  = qidep*(1._rtype-ratio)
      qisub  = qisub*(1._rtype-ratio)
   end if


end subroutine cloud_water_conservation

subroutine rain_water_conservation(qr,qcaut,qcacc,qimlt,qcshd,dt,    &
   qrevp,qrcol,qrheti)

   implicit none

   real(rtype), intent(in) :: qr, qcaut, qcacc, qimlt, qcshd, dt
   real(rtype), intent(inout) :: qrevp, qrcol, qrheti

   real(rtype) :: sinks, sources, ratio

   sinks   = (qrevp+qrcol+qrheti)*dt
   sources = qr + (qcaut+qcacc+qimlt+qcshd)*dt
   if (sinks.gt.sources .and. sinks.ge.1.e-20_rtype) then
      ratio  = sources/sinks
      qrevp  = qrevp*ratio
      qrcol  = qrcol*ratio
      qrheti = qrheti*ratio
   endif

end subroutine rain_water_conservation

subroutine ice_water_conservation(qitot,qidep,qinuc,qiberg,qrcol,qccol,qrheti,qcheti,dt,    &
   qisub,qimlt)

   implicit none

   real(rtype), intent(in) :: qitot, qidep, qinuc, qrcol, qccol, qrheti, qcheti, qiberg, dt
   real(rtype), intent(inout) :: qisub, qimlt
   real(rtype) :: sinks, sources, ratio

   sinks   = (qisub+qimlt)*dt
   sources = qitot + (qidep+qinuc+qrcol+qccol+  &
        qrheti+qcheti+qiberg)*dt
   if (sinks.gt.sources .and. sinks.ge.1.e-20_rtype) then
      ratio = sources/sinks
      qisub = qisub*ratio
      qimlt = qimlt*ratio
   endif

end subroutine ice_water_conservation


subroutine update_prognostic_ice(qcheti,qccol,qcshd,    &
   nccol,ncheti,ncshdc,    &
   qrcol,nrcol,qrheti,nrheti,nrshdr,    &
   qimlt,nimlt,qisub,qidep,qinuc,ninuc,nislf,nisub,qiberg,    &
   exner,xxls,xlf,    &
   log_predictNc,log_wetgrowth,dt,nmltratio,rhorime_c,    &
   th,qv,qitot,nitot,qirim,birim,qc,nc,qr,nr)

   !-- ice-phase dependent processes:
   implicit none

   real(rtype), intent(in) :: qcheti
   real(rtype), intent(in) :: qccol
   real(rtype), intent(in) :: qcshd
   real(rtype), intent(in) :: nccol
   real(rtype), intent(in) :: ncheti
   real(rtype), intent(in) :: ncshdc

   real(rtype), intent(in) :: qrcol
   real(rtype), intent(in) :: nrcol
   real(rtype), intent(in) :: qrheti
   real(rtype), intent(in) :: nrheti
   real(rtype), intent(in) :: nrshdr

   real(rtype), intent(in) :: qimlt
   real(rtype), intent(in) :: nimlt
   real(rtype), intent(in) :: qisub
   real(rtype), intent(in) :: qidep
   real(rtype), intent(in) :: qinuc
   real(rtype), intent(in) :: ninuc
   real(rtype), intent(in) :: nislf
   real(rtype), intent(in) :: nisub
   real(rtype), intent(in) :: qiberg
   real(rtype), intent(in) :: exner
   real(rtype), intent(in) :: xlf
   real(rtype), intent(in) :: xxls

   logical(btype), intent(in) :: log_predictNc
   logical(btype), intent(in) :: log_wetgrowth
   real(rtype), intent(in) :: dt
   real(rtype), intent(in) :: nmltratio
   real(rtype), intent(in) :: rhorime_c

   real(rtype), intent(inout) :: th
   real(rtype), intent(inout) :: qv
   real(rtype), intent(inout) :: qc
   real(rtype), intent(inout) :: nc
   real(rtype), intent(inout) :: qr
   real(rtype), intent(inout) :: nr
   real(rtype), intent(inout) :: qitot
   real(rtype), intent(inout) :: nitot
   real(rtype), intent(inout) :: qirim
   real(rtype), intent(inout) :: birim

   real(rtype) :: dum

   qc = qc + (-qcheti-qccol-qcshd-qiberg)*dt
   if (log_predictNc) then
      nc = nc + (-nccol-ncheti)*dt
   endif
   qr = qr + (-qrcol+qimlt-qrheti+qcshd)*dt

   ! apply factor to source for rain number from melting of ice, (ad-hoc
   ! but accounts for rapid evaporation of small melting ice particles)
   nr = nr + (-nrcol-nrheti+nmltratio*nimlt+nrshdr+ncshdc)*dt

   if (qitot.ge.qsmall) then
      ! add sink terms, assume density stays constant for sink terms
      birim = birim - ((qisub+qimlt)/qitot)*dt*birim
      qirim = qirim - ((qisub+qimlt)*qirim/qitot)*dt
      qitot = qitot - (qisub+qimlt)*dt
   endif

   dum = (qrcol+qccol+qrheti+qcheti)*dt
   qitot = qitot + (qidep+qinuc+qiberg)*dt + dum
   qirim = qirim + dum

   birim = birim + (qrcol*inv_rho_rimeMax+qccol/rhorime_c+(qrheti+ &
        qcheti)*inv_rho_rimeMax)*dt

   nitot = nitot + (ninuc-nimlt-nisub-nislf+nrheti+ncheti)*dt

   !PMC nCat deleted interactions_loop

   if (qirim.lt.0._rtype) then
      qirim = 0._rtype
      birim = 0._rtype
   endif

   ! densify under wet growth
   ! -- to be removed post-v2.1.  Densification automatically happens
   !    during wet growth due to parameterized rime density --
   if (log_wetgrowth) then
      qirim = qitot
      birim = qirim*inv_rho_rimeMax
   endif

   ! densify in above freezing conditions and melting
   ! -- future work --
   !   Ideally, this will be treated with the predicted liquid fraction in ice.
   !   Alternatively, it can be simplified by tending qirim -- qitot
   !   and birim such that rho_rim (qirim/birim) --> rho_liq during melting.
   ! ==

   qv = qv + (-qidep+qisub-qinuc)*dt

   th = th + exner*((qidep-qisub+qinuc)*xxls*inv_cp +(qrcol+qccol+   &
       qcheti+qrheti-qimlt+qiberg)* xlf*inv_cp)*dt
end subroutine update_prognostic_ice

subroutine update_prognostic_liquid(qcacc,ncacc,qcaut,ncautc,ncautr,ncslf,    &
    qrevp,nrevp,nrslf,                                                        &
    log_predictNc,inv_rho,exner,xxlv,dt,                                      &
    th,qv,qc,nc,qr,nr)

   !-- warm-phase only processes:
   implicit none

   real(rtype), intent(in) :: qcacc
   real(rtype), intent(in) :: ncacc
   real(rtype), intent(in) :: qcaut
   real(rtype), intent(in) :: ncautc
   real(rtype), intent(in) :: ncautr
   real(rtype), intent(in) :: ncslf
   real(rtype), intent(in) :: qrevp
   real(rtype), intent(in) :: nrevp
   real(rtype), intent(in) :: nrslf


   logical(btype), intent(in) :: log_predictNc
   real(rtype), intent(in) :: inv_rho
   real(rtype), intent(in) :: exner
   real(rtype), intent(in) :: xxlv
   real(rtype), intent(in) :: dt

   real(rtype), intent(inout) :: th
   real(rtype), intent(inout) :: qv
   real(rtype), intent(inout) :: qc
   real(rtype), intent(inout) :: nc
   real(rtype), intent(inout) :: qr
   real(rtype), intent(inout) :: nr

   qc = qc + (-qcacc-qcaut)*dt
   qr = qr + (qcacc+qcaut-qrevp)*dt

   if (log_predictNc) then
      nc = nc + (-ncacc-ncautc+ncslf)*dt
   else
      nc = nccnst*inv_rho
   endif
   if (iparam.eq.1 .or. iparam.eq.2) then
      nr = nr + (0.5_rtype*ncautc-nrslf-nrevp)*dt
   else
      nr = nr + (ncautr-nrslf-nrevp)*dt
   endif

   qv = qv + qrevp*dt
   th = th + exner*(-qrevp*xxlv*    &
        inv_cp)*dt

end subroutine update_prognostic_liquid



subroutine ice_deposition_sublimation(qitot_incld,nitot_incld,t,    &
qvs,qvi,epsi,abi,qv,    &
qidep,qisub,nisub,qiberg)

   implicit none

   real(rtype), intent(in)  :: qitot_incld
   real(rtype), intent(in)  :: nitot_incld
   real(rtype), intent(in)  :: t
   real(rtype), intent(in)  :: qvs
   real(rtype), intent(in)  :: qvi
   real(rtype), intent(in)  :: epsi
   real(rtype), intent(in)  :: abi
   real(rtype), intent(in)  :: qv
   real(rtype), intent(out) :: qidep
   real(rtype), intent(out) :: qisub
   real(rtype), intent(out) :: nisub
   real(rtype), intent(out) :: qiberg

   real(rtype) :: oabi

   oabi = 1._rtype/abi
   if (qitot_incld>=qsmall) then
      !Compute deposition/sublimation
      qidep = epsi * oabi * (qv - qvi)
      !Split into deposition or sublimation.
      if (t < zerodegc .and. qidep>0._rtype) then
         qisub=0._rtype
      else
      ! make qisub positive for consistency with other evap/sub processes
         qisub=-min(qidep,0._rtype)
         qidep=0._rtype
      end if
      !sublimation occurs @ any T. Not so for berg.
      if (t < zerodegc) then
         !Compute bergeron rate assuming cloud for whole step.
         qiberg = max(epsi*oabi*(qvs - qvi), 0._rtype)
      else !T>frz
         qiberg=0._rtype
      end if !T<frz
      nisub = qisub*(nitot_incld/qitot_incld)
   else
      qiberg = 0._rtype
      qidep  = 0._rtype
      qisub  = 0._rtype
      nisub  = 0._rtype
   end if

   return

end subroutine ice_deposition_sublimation


subroutine evaporate_sublimate_precip(qr_incld,qc_incld,nr_incld,qitot_incld,    &
lcldm,rcldm,qvs,ab,epsr,qv,    &
qrevp,nrevp)

   implicit none

   real(rtype), intent(in)  :: qr_incld
   real(rtype), intent(in)  :: qc_incld
   real(rtype), intent(in)  :: nr_incld
   real(rtype), intent(in)  :: qitot_incld
   real(rtype), intent(in)  :: lcldm
   real(rtype), intent(in)  :: rcldm
   real(rtype), intent(in)  :: qvs
   real(rtype), intent(in)  :: ab
   real(rtype), intent(in)  :: epsr
   real(rtype), intent(in)  :: qv
   real(rtype), intent(out) :: qrevp
   real(rtype), intent(out) :: nrevp

   real(rtype) :: qclr, cld

   ! It is assumed that macrophysics handles condensation/evaporation of qc and
   ! that there is no condensation of rain. Thus qccon, qrcon and qcevp have
   ! been removed from the original P3-WRF.

   ! Determine temporary cloud fraction, set to zero if cloud water + ice is
   ! very small.  This will ensure that evap/subl of precip occurs over entire
   ! grid cell, since min cloud fraction is specified otherwise.
   if (qc_incld + qitot_incld < 1.e-6_rtype) then
      cld = 0._rtype
   else
      cld = lcldm
   end if

   ! Only calculate if there is some rain fraction > cloud fraction
   qrevp = 0.0_rtype
   if (rcldm > cld) then
      ! calculate q for out-of-cloud region
      qclr = (qv-cld*qvs)/(1._rtype-cld)

      ! rain evaporation
      if (qr_incld.ge.qsmall) then
         qrevp = epsr * (qclr-qvs)/ab
      end if

      ! only evap in out-of-cloud region
      qrevp = -min(qrevp*(rcldm-cld),0._rtype)
      qrevp = qrevp/rcldm
   end if ! rcld>cld
   if (qr_incld.gt.qsmall)  nrevp = qrevp*(nr_incld/qr_incld)

   return

end subroutine evaporate_sublimate_precip

subroutine get_time_space_phys_variables( &
t,pres,rho,xxlv,xxls,qvs,qvi, &
mu,dv,sc,dqsdt,dqsidt,ab,abi,kap,eii)

   implicit none

   real(rtype), intent(in)  :: t
   real(rtype), intent(in)  :: pres
   real(rtype), intent(in)  :: rho
   real(rtype), intent(in)  :: xxlv
   real(rtype), intent(in)  :: xxls
   real(rtype), intent(in)  :: qvs
   real(rtype), intent(in)  :: qvi
   real(rtype), intent(out) :: mu
   real(rtype), intent(out) :: dv
   real(rtype), intent(out) :: sc
   real(rtype), intent(out) :: dqsdt
   real(rtype), intent(out) :: dqsidt
   real(rtype), intent(out) :: ab
   real(rtype), intent(out) :: abi
   real(rtype), intent(out) :: kap
   real(rtype), intent(out) :: eii

   real(rtype) :: dum

   !time/space varying physical variables
   mu     = 1.496e-6_rtype*bfb_pow(t,1.5_rtype)/(t+120._rtype)
   dv     = 8.794e-5_rtype*bfb_pow(t,1.81_rtype)/pres
   sc     = mu/(rho*dv)
   dum    = 1._rtype/(rv*bfb_square(t))
   dqsdt  = xxlv*qvs*dum
   dqsidt = xxls*qvi*dum
   ab     = 1._rtype+dqsdt*xxlv*inv_cp
   abi    = 1._rtype+dqsidt*xxls*inv_cp
   kap    = 1.414e+3_rtype*mu

   ! very simple temperature dependent aggregation efficiency
   if (t.lt.253.15_rtype) then
      eii=0.1_rtype
   else if (t.ge.253.15_rtype.and.t.lt.268.15_rtype) then
      eii=0.1_rtype+(t-253.15_rtype)/15._rtype*0.9_rtype  ! linear ramp from 0.1 to 1 between 253.15 and 268.15 K
   else
      eii=1._rtype
   end if

   return

end subroutine get_time_space_phys_variables

subroutine cloud_sedimentation(kts,kte,ktop,kbot,kdir,   &
   qc_incld,rho,inv_rho,lcldm,acn,inv_dzq,&
   dt,odt,dnu,log_predictNc, &
   qc, nc, nc_incld,mu_c,lamc,prt_liq,qc_tend,nc_tend)

   implicit none
   integer, intent(in) :: kts, kte
   integer, intent(in) :: ktop, kbot, kdir

   real(rtype), intent(in), dimension(kts:kte) :: rho
   real(rtype), intent(in), dimension(kts:kte) :: inv_rho
   real(rtype), intent(in), dimension(kts:kte) :: lcldm
   real(rtype), intent(in), dimension(kts:kte) :: acn
   real(rtype), intent(in), dimension(kts:kte) :: inv_dzq
   real(rtype), intent(in) :: dt
   real(rtype), intent(in) :: odt
   real(rtype), dimension(:), intent(in) :: dnu
   logical(btype), intent(in) :: log_predictNc

   real(rtype), intent(inout), dimension(kts:kte), target :: qc
   real(rtype), intent(inout), dimension(kts:kte), target :: nc
   real(rtype), intent(inout), dimension(kts:kte) :: qc_incld
   real(rtype), intent(inout), dimension(kts:kte) :: nc_incld
   real(rtype), intent(inout), dimension(kts:kte) :: mu_c
   real(rtype), intent(inout), dimension(kts:kte) :: lamc
   real(rtype), intent(inout) :: prt_liq
   real(rtype), intent(inout), dimension(kts:kte) :: qc_tend
   real(rtype), intent(inout), dimension(kts:kte) :: nc_tend

   logical(btype) :: log_qxpresent
   integer :: k
   integer :: k_qxtop, k_qxbot
   integer, parameter :: num_arrays = 2
   type(realptr), dimension(num_arrays) :: vs, fluxes, qnr

   real(rtype) :: dt_left
   real(rtype) :: prt_accum
   real(rtype) :: Co_max
   real(rtype) :: nu
   real(rtype), dimension(kts:kte), target :: V_qc
   real(rtype), dimension(kts:kte), target :: V_nc
   real(rtype), dimension(kts:kte), target :: flux_qx
   real(rtype), dimension(kts:kte), target :: flux_nx

   real(rtype) :: tmp1, tmp2, dum

   k_qxtop = kbot
   log_qxpresent = .false.

   vs(1)%p => V_qc
   vs(2)%p => V_nc
   fluxes(1)%p => flux_qx
   fluxes(2)%p => flux_nx
   qnr(1)%p => qc
   qnr(2)%p => nc

   !find top, determine qxpresent
   do k = ktop,kbot,-kdir
      if (qc(k).ge.qsmall) then
         log_qxpresent = .true.
         k_qxtop = k
         exit
      endif
   enddo

   qc_present: if (log_qxpresent) then

      dt_left   = dt  !time remaining for sedi over full model (mp) time step
      prt_accum = 0._rtype  !precip rate for individual category

      !find bottom
      do k = kbot,k_qxtop,kdir
         if (qc(k).ge.qsmall) then
            k_qxbot = k
            exit
         endif
      enddo

      two_moment: if (log_predictNc) then  !2-moment cloud:
         substep_sedi_c2: do while (dt_left.gt.1.e-4_rtype)

            Co_max = 0._rtype
            V_qc = 0._rtype
            V_nc = 0._rtype

            kloop_sedi_c2: do k = k_qxtop,k_qxbot,-kdir

               qc_notsmall_c2: if (qc_incld(k)>qsmall) then
                  !-- compute Vq, Vn
                  call get_cloud_dsd2(qc_incld(k),nc_incld(k),mu_c(k),rho(k),nu,dnu,   &
                       lamc(k),tmp1,tmp2,lcldm(k))

                  !get_cloud_dsd2 keeps the drop-size distribution within reasonable
                  !bounds by modifying nc_incld. The next line maintains consistency
                  !between nc_incld and nc
                  nc(k) = nc_incld(k)*lcldm(k)
                  
                  dum = 1._rtype / bfb_pow(lamc(k), bcn)
                  V_qc(k) = acn(k)*bfb_gamma(4._rtype+bcn+mu_c(k))*dum/(bfb_gamma(mu_c(k)+4._rtype))
                  V_nc(k) = acn(k)*bfb_gamma(1._rtype+bcn+mu_c(k))*dum/(bfb_gamma(mu_c(k)+1._rtype))

               endif qc_notsmall_c2
               Co_max = max(Co_max, V_qc(k)*dt_left*inv_dzq(k))

            enddo kloop_sedi_c2

            call generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, inv_dzq, inv_rho, rho, num_arrays, vs, fluxes, qnr)

            !Update _incld values with end-of-step cell-ave values
            !Note that lcldm is set in interface to have min of mincld=1e-4
            !so dividing by it is fine.
            qc_incld(:) = qc(:)/lcldm(:)
            nc_incld(:) = nc(:)/lcldm(:)

         enddo substep_sedi_c2
      else
         substep_sedi_c1: do while (dt_left.gt.1.e-4_rtype)

            Co_max  = 0._rtype
            V_qc = 0._rtype

            kloop_sedi_c1: do k = k_qxtop,k_qxbot,-kdir
               qc_notsmall_c1: if (qc_incld(k)>qsmall) then
                  call get_cloud_dsd2(qc_incld(k),nc_incld(k),mu_c(k),rho(k),nu,dnu,   &
                       lamc(k),tmp1,tmp2,lcldm(k))
                  
                  !get_cloud_dsd2 keeps the drop-size distribution within reasonable
                  !bounds by modifying nc_incld. The next line maintains consistency
                  !between nc_incld and nc
                  nc(k) = nc_incld(k)*lcldm(k)
                  
                  dum = 1._rtype / bfb_pow(lamc(k), bcn)
                  V_qc(k) = acn(k)*bfb_gamma(4._rtype+bcn+mu_c(k))*dum/(bfb_gamma(mu_c(k)+4._rtype))
               endif qc_notsmall_c1

               Co_max = max(Co_max, V_qc(k)*dt_left*inv_dzq(k))
            enddo kloop_sedi_c1

            call generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, inv_dzq, inv_rho, rho, 1, vs, fluxes, qnr)

            !Update _incld values with end-of-step cell-ave values
            !Note that lcldm is set in interface to have min of mincld=1e-4
            !so dividing by it is fine.
            qc_incld(:) = qc(:)/lcldm(:)
            nc_incld(:) = nc(:)/lcldm(:)
            
         enddo substep_sedi_c1

      endif two_moment

      ! JGF: Is prt_liq intended to be inout or just out? Inconsistent with rain and ice sed.
      prt_liq = prt_accum*inv_rhow*odt  !note, contribution from rain is added below

   endif qc_present

   qc_tend(:) = ( qc(:) - qc_tend(:) ) * odt ! Liq. sedimentation tendency, measure
   nc_tend(:) = ( nc(:) - nc_tend(:) ) * odt ! Liq. # sedimentation tendency, measure

end subroutine cloud_sedimentation

subroutine rain_sedimentation(kts,kte,ktop,kbot,kdir,   &
   qr_incld,rho,inv_rho,rhofacr,rcldm,inv_dzq,dt,odt,  &
   qr,nr,nr_incld,mu_r,lamr,prt_liq,rflx,qr_tend,nr_tend)

   implicit none
   integer, intent(in) :: kts, kte
   integer, intent(in) :: ktop, kbot, kdir


   real(rtype), intent(in), dimension(kts:kte) :: rho
   real(rtype), intent(in), dimension(kts:kte) :: inv_rho
   real(rtype), intent(in), dimension(kts:kte) :: rhofacr
   real(rtype), intent(in), dimension(kts:kte) :: rcldm
   real(rtype), intent(in), dimension(kts:kte) :: inv_dzq
   real(rtype), intent(in) :: dt
   real(rtype), intent(in) :: odt

   real(rtype), intent(inout), target, dimension(kts:kte) :: qr
   real(rtype), intent(inout), target, dimension(kts:kte) :: nr
   real(rtype), intent(inout), dimension(kts:kte) :: qr_incld
   real(rtype), intent(inout), dimension(kts:kte) :: nr_incld
   real(rtype), intent(inout), dimension(kts:kte) :: mu_r
   real(rtype), intent(inout), dimension(kts:kte) :: lamr
   real(rtype), intent(inout) :: prt_liq
   real(rtype), intent(inout), dimension(kts:kte+1) :: rflx
   real(rtype), intent(inout), dimension(kts:kte) :: qr_tend
   real(rtype), intent(inout), dimension(kts:kte) :: nr_tend

   logical(btype) :: log_qxpresent
   integer :: k
   integer :: k_qxtop, k_qxbot
   integer, parameter :: num_arrays = 2
   type(realptr), dimension(num_arrays) :: vs, fluxes, qnr

   real(rtype) :: dt_left
   real(rtype) :: prt_accum
   real(rtype) :: Co_max
   real(rtype), dimension(kts:kte), target :: V_qr
   real(rtype), dimension(kts:kte), target :: V_nr
   real(rtype), dimension(kts:kte), target :: flux_qx
   real(rtype), dimension(kts:kte), target :: flux_nx

   vs(1)%p => V_qr
   vs(2)%p => V_nr
   fluxes(1)%p => flux_qx
   fluxes(2)%p => flux_nx
   qnr(1)%p => qr
   qnr(2)%p => nr

   k_qxtop = kbot
   log_qxpresent = .false.

   !find top, determine qxpresent
   do k = ktop,kbot,-kdir
      if (qr(k).ge.qsmall) then
         log_qxpresent = .true.
         k_qxtop = k
         exit
      endif !
   enddo

   qr_present: if (log_qxpresent) then

      dt_left   = dt  !time remaining for sedi over full model (mp) time step
      prt_accum = 0._rtype  !precip rate for individual category

      !find bottom
      do k = kbot,k_qxtop,kdir
         if (qr(k).ge.qsmall) then
            k_qxbot = k
            exit
         endif
      enddo

      substep_sedi_r: do while (dt_left.gt.1.e-4_rtype)

         Co_max = 0._rtype
         V_qr = 0._rtype
         V_nr = 0._rtype

         kloop_sedi_r1: do k = k_qxtop,k_qxbot,-kdir

            qr_notsmall_r1: if (qr_incld(k)>qsmall) then

               call compute_rain_fall_velocity(qr_incld(k), rcldm(k), rhofacr(k), nr_incld(k), &
                    mu_r(k), lamr(k), V_qr(k), V_nr(k))

               !in compute_rain_fall_velocity, get_rain_dsd2 keeps the drop-size
               !distribution within reasonable bounds by modifying nr_incld. 
               !The next line maintains consistency between nr_incld and nr.
               nr = nr_incld*rcldm 

            endif qr_notsmall_r1

            Co_max = max(Co_max, V_qr(k)*dt_left*inv_dzq(k))
            !            Co_max = max(Co_max, max(V_nr(k),V_qr(k))*dt_left*inv_dzq(i,k))

         enddo kloop_sedi_r1

         call generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, inv_dzq, inv_rho, rho, num_arrays, vs, fluxes, qnr)

         !-- AaronDonahue, rflx output
         do k = k_qxbot,k_qxtop,kdir
            rflx(k+1) = rflx(k+1) + flux_qx(k) ! AaronDonahue
         enddo

         !Update _incld values with end-of-step cell-ave values
         !Note that rcldm is set in interface to have min of mincld=1e-4
         !so dividing by it is fine.
         qr_incld(:) = qr(:)/rcldm(:)
         nr_incld(:) = nr(:)/rcldm(:)
         
      enddo substep_sedi_r

      prt_liq = prt_liq + prt_accum*inv_rhow*odt

   endif qr_present

   qr_tend(:) = ( qr(:) - qr_tend(:) ) * odt ! Rain sedimentation tendency, measure
   nr_tend(:) = ( nr(:) - nr_tend(:) ) * odt ! Rain # sedimentation tendency, measure

end subroutine rain_sedimentation

subroutine compute_rain_fall_velocity(qr_incld, rcldm, rhofacr, nr_incld, mu_r, lamr, V_qr, V_nr)

   real(rtype), intent(in) :: qr_incld
   real(rtype), intent(in) :: rcldm
   real(rtype), intent(in) :: rhofacr
   real(rtype), intent(inout) :: nr_incld
   real(rtype), intent(out) :: mu_r
   real(rtype), intent(out) :: lamr
   real(rtype), intent(out) :: V_qr
   real(rtype), intent(out) :: V_nr

   real(rtype) :: tmp1, tmp2, dum1, dum2, inv_dum3, rdumii, rdumjj
   integer :: dumii, dumjj

   !Compute Vq, Vn:

   call get_rain_dsd2(qr_incld,nr_incld,mu_r,lamr,tmp1,tmp2,rcldm)

   call find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3,mu_r,lamr)

   !mass-weighted fall speed:

   dum1 = vm_table(dumii,dumjj)+(rdumii-real(dumii))*                       &
      (vm_table(dumii+1,dumjj)-vm_table(dumii,dumjj))       !at mu_r
   dum2 = vm_table(dumii,dumjj+1)+(rdumii-real(dumii))*                     &
      (vm_table(dumii+1,dumjj+1)-vm_table(dumii,dumjj+1))   !at mu_r+1

   V_qr = dum1 + (rdumjj-real(dumjj))*(dum2-dum1)         !interpolated
   V_qr = V_qr*rhofacr                                    !corrected for air density

   ! number-weighted fall speed:
   dum1 = vn_table(dumii,dumjj)+(rdumii-real(dumii))*                       &
      (vn_table(dumii+1,dumjj)-vn_table(dumii,dumjj))       !at mu_r
   dum2 = vn_table(dumii,dumjj+1)+(rdumii-real(dumii))*                     &
      (vn_table(dumii+1,dumjj+1)-vn_table(dumii,dumjj+1))   !at mu_r+1

   V_nr = dum1+(rdumjj-real(dumjj))*(dum2-dum1)            !interpolated
   V_nr = V_nr*rhofacr               !corrected for air density
end subroutine compute_rain_fall_velocity

subroutine ice_sedimentation(kts,kte,ktop,kbot,kdir,    &
   rho,inv_rho,rhofaci,icldm,inv_dzq,dt,odt,  &
   qitot,qitot_incld,nitot,qirim,qirim_incld,birim,birim_incld,nitot_incld,prt_sol,qi_tend,ni_tend)

   implicit none
   integer, intent(in) :: kts, kte
   integer, intent(in) :: ktop, kbot, kdir

   real(rtype), intent(in), dimension(kts:kte) :: rho
   real(rtype), intent(in), dimension(kts:kte) :: inv_rho
   real(rtype), intent(in), dimension(kts:kte) :: rhofaci
   real(rtype), intent(in), dimension(kts:kte) :: icldm
   real(rtype), intent(in), dimension(kts:kte) :: inv_dzq
   real(rtype), intent(in) :: dt
   real(rtype), intent(in) :: odt

   real(rtype), intent(inout), dimension(kts:kte), target :: qitot
   real(rtype), intent(inout), dimension(kts:kte) :: qitot_incld
   real(rtype), intent(inout), dimension(kts:kte), target :: nitot
   real(rtype), intent(inout), dimension(kts:kte) :: nitot_incld
   real(rtype), intent(inout), dimension(kts:kte), target :: qirim
   real(rtype), intent(inout), dimension(kts:kte) :: qirim_incld
   real(rtype), intent(inout), dimension(kts:kte), target :: birim
   real(rtype), intent(inout), dimension(kts:kte) :: birim_incld

   real(rtype), intent(inout) :: prt_sol
   real(rtype), intent(inout), dimension(kts:kte) :: qi_tend
   real(rtype), intent(inout), dimension(kts:kte) :: ni_tend

   logical(btype) :: log_qxpresent
   integer :: k
   integer :: k_qxtop, k_qxbot
   integer, parameter :: num_arrays = 4
   type(realptr), dimension(num_arrays) :: vs, fluxes, qnr

   real(rtype) :: dt_left
   real(rtype) :: prt_accum
   real(rtype) :: Co_max
   real(rtype) :: rhop
   real(rtype), dimension(kts:kte), target :: V_qit
   real(rtype), dimension(kts:kte), target :: V_nit
   real(rtype), dimension(kts:kte), target :: flux_nit
   real(rtype), dimension(kts:kte), target :: flux_bir
   real(rtype), dimension(kts:kte), target :: flux_qir
   real(rtype), dimension(kts:kte), target :: flux_qit
   real(rtype) :: f1pr01 ! number-weighted fallspeed            See lines  731 -  808  uns
   real(rtype) :: f1pr02 ! mass-weighted fallspeed              See lines  731 -  808  ums
   real(rtype) :: f1pr09 ! minimum ice number (lambda limiter)  See lines  704 -  705  nlarge
   real(rtype) :: f1pr10 ! maximum ice number (lambda limiter)  See lines  704 -  705  nsmall

   real(rtype) :: dum1, dum4, dum5, dum6
   integer dumi, dumii, dumjj, dumzz

   log_qxpresent = .false.  !note: this applies to ice category 'iice' only
   k_qxtop       = kbot

   vs(1)%p => V_qit
   vs(2)%p => V_nit
   vs(3)%p => V_qit
   vs(4)%p => V_qit
   fluxes(1)%p => flux_qit
   fluxes(2)%p => flux_nit
   fluxes(3)%p => flux_qir
   fluxes(4)%p => flux_bir
   qnr(1)%p => qitot
   qnr(2)%p => nitot
   qnr(3)%p => qirim
   qnr(4)%p => birim

   !find top, determine qxpresent
   do k = ktop,kbot,-kdir
      if (qitot(k).ge.qsmall) then
         log_qxpresent = .true.
         k_qxtop = k
         exit
      endif !
   enddo  !k-loop

   qi_present: if (log_qxpresent) then

      dt_left   = dt  !time remaining for sedi over full model (mp) time step
      prt_accum = 0._rtype  !precip rate for individual category

      !find bottom
      do k = kbot,k_qxtop,kdir
         if (qitot(k).ge.qsmall) then
            k_qxbot = k
            exit
         endif
      enddo

      substep_sedi_i: do while (dt_left.gt.1.e-4_rtype)

         Co_max = 0._rtype
         V_qit = 0._rtype
         V_nit = 0._rtype

         kloop_sedi_i1: do k = k_qxtop,k_qxbot,-kdir

            !-- compute Vq, Vn (get values from lookup table)
            qi_notsmall_i1: if (qitot_incld(k)>qsmall) then

               !--Compute Vq, Vn:
               nitot_incld(k) = max(nitot_incld(k),nsmall) !impose lower limits to prevent log(<0)
               call calc_bulkRhoRime(qitot_incld(k),qirim_incld(k),birim_incld(k),rhop)
               !if (.not. tripleMoment_on) zitot(i,k) = diag_mom6(qitot(i,k),nitot(i,k),rho(i,k))
               call find_lookupTable_indices_1a(dumi,dumjj,dumii,dumzz,dum1,dum4,    &
                    dum5,dum6,isize,rimsize,densize,          &
                    qitot_incld(k),nitot_incld(k),qirim_incld(k),&
                    rhop)
               call access_lookup_table(dumjj,dumii,dumi, 1,dum1,dum4,dum5,f1pr01)
               call access_lookup_table(dumjj,dumii,dumi, 2,dum1,dum4,dum5,f1pr02)
               call access_lookup_table(dumjj,dumii,dumi, 7,dum1,dum4,dum5,f1pr09)
               call access_lookup_table(dumjj,dumii,dumi, 8,dum1,dum4,dum5,f1pr10)
               !-impose mean ice size bounds (i.e. apply lambda limiters)
               ! note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
               nitot_incld(k) = min(nitot_incld(k),f1pr09*nitot_incld(k))
               nitot_incld(k) = max(nitot_incld(k),f1pr10*nitot_incld(k))
               nitot(k) = nitot_incld(k)*icldm(k)
               !zitot(i,k) = min(zitot(i,k),f1pr020)  !adjust Zi if needed to make sure mu_i is in bounds
               !zitot(i,k) = max(zitot(i,k),f1pr021)
               V_qit(k) = f1pr02*rhofaci(k)     !mass-weighted  fall speed (with density factor)
               V_nit(k) = f1pr01*rhofaci(k)     !number-weighted    fall speed (with density factor)
               !==

            endif qi_notsmall_i1

            Co_max = max(Co_max, V_qit(k)*dt_left*inv_dzq(k))

         enddo kloop_sedi_i1

         call generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, inv_dzq, inv_rho, rho, num_arrays, vs, fluxes, qnr)

         !update _incld variables
         !Note that icldm is set in interface to have min of mincld=1e-4
         !so dividing by it is fine.
         qitot_incld(:) = qitot(:)/icldm(:)
         nitot_incld(:) = nitot(:)/icldm(:)
         qirim_incld(:) = qirim(:)/icldm(:)
         birim_incld(:) = birim(:)/icldm(:)
         
      enddo substep_sedi_i

      prt_sol = prt_sol + prt_accum*inv_rhow*odt

   endif qi_present

   qi_tend(:) = ( qitot(:) - qi_tend(:) ) * odt ! Ice sedimentation tendency, measure
   ni_tend(:) = ( nitot(:) - ni_tend(:) ) * odt ! Ice # sedimentation tendency, measure

end subroutine ice_sedimentation

subroutine generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, inv_dzq, inv_rho, rho, &
     num_arrays, vs, fluxes, qnx)

   implicit none

   integer, intent(in) :: kts, kte, kdir, k_qxtop, kbot, num_arrays
   integer, intent(inout) :: k_qxbot
   real(rtype), intent(in) :: Co_max
   real(rtype), intent(inout) :: dt_left, prt_accum
   real(rtype), dimension(kts:kte), intent(in) :: inv_dzq
   real(rtype), dimension(kts:kte), intent(in) :: inv_rho
   real(rtype), dimension(kts:kte), intent(in) :: rho

   type(realptr), intent(in), dimension(num_arrays), target :: vs, fluxes, qnx

   integer :: tmpint1, k_temp, i
   real(rtype) :: dt_sub

   !-- compute dt_sub
   tmpint1 = int(Co_max+1._rtype)    !number of substeps remaining if dt_sub were constant
   dt_sub  = min(dt_left, dt_left/float(tmpint1))

   ! -- Move bottom cell down by 1 if not at ground already
   if (k_qxbot.eq.kbot) then
      k_temp = k_qxbot
   else
      k_temp = k_qxbot-kdir
   endif

   call calc_first_order_upwind_step(kts, kte, kdir, k_temp, k_qxtop, dt_sub, rho, inv_rho, inv_dzq, num_arrays, fluxes, vs, qnx)

   !accumulated precip during time step
   if (k_qxbot.eq.kbot) prt_accum = prt_accum + fluxes(1)%p(kbot)*dt_sub

   dt_left = dt_left - dt_sub  !update time remaining for sedimentation
   if (k_qxbot.ne.kbot) k_qxbot = k_qxbot - kdir

end subroutine generalized_sedimentation

subroutine calc_first_order_upwind_step(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dzq, num_arrays, fluxes, vs, qnx)

  implicit none

  integer, intent(in) :: kts, kte, kdir, kbot, k_qxtop, num_arrays
  real(rtype), intent(in) :: dt_sub
  real(rtype), dimension(kts:kte), intent(in) :: rho, inv_rho, inv_dzq
  type(realptr), intent(in), dimension(num_arrays), target :: fluxes, vs, qnx

  integer :: i, k
  real(rtype) :: fluxdiv

  !-- calculate fluxes
  do k = kbot,k_qxtop,kdir
     do i = 1, num_arrays
        fluxes(i)%p(k) = vs(i)%p(k) * qnx(i)%p(k) * rho(k)
     end do
  enddo

  do i = 1, num_arrays
     k = k_qxtop

     !--- for top level only (since flux is 0 above)

     !- compute flux divergence
     fluxdiv = -fluxes(i)%p(k) * inv_dzq(k)
     !- update prognostic variables
     qnx(i)%p(k) = qnx(i)%p(k) + fluxdiv*dt_sub*inv_rho(k)

     do k = k_qxtop-kdir,kbot,-kdir
        !-- compute flux divergence
        fluxdiv = (fluxes(i)%p(k+kdir) - fluxes(i)%p(k))*inv_dzq(k)
        !-- update prognostic variables
        qnx(i)%p(k) = qnx(i)%p(k) + fluxdiv*dt_sub*inv_rho(k)
     end do
  end do

end subroutine calc_first_order_upwind_step

subroutine homogeneous_freezing(kts,kte,ktop,kbot,kdir,t,exner,xlf,    &
   qc,nc,qr,nr,qitot,nitot,qirim,birim,th)

   !.......................................
   ! homogeneous freezing of cloud and rain

   implicit none
   integer, intent(in) :: kts, kte
   integer, intent(in) :: ktop, kbot, kdir
   real(rtype), intent(in), dimension(kts:kte) :: t
   real(rtype), intent(in), dimension(kts:kte) :: exner
   real(rtype), intent(in), dimension(kts:kte) :: xlf

   real(rtype), intent(inout), dimension(kts:kte) :: qc
   real(rtype), intent(inout), dimension(kts:kte) :: nc
   real(rtype), intent(inout), dimension(kts:kte) :: qr
   real(rtype), intent(inout), dimension(kts:kte) :: nr

   real(rtype), intent(inout), dimension(kts:kte) :: qitot
   real(rtype), intent(inout), dimension(kts:kte) :: nitot
   real(rtype), intent(inout), dimension(kts:kte) :: qirim
   real(rtype), intent(inout), dimension(kts:kte) :: birim
   real(rtype), intent(inout), dimension(kts:kte) :: th

   real(rtype) :: Q_nuc
   real(rtype) :: N_nuc
   integer :: k

   k_loop_fz:  do k = kbot,ktop,kdir
      if (qc(k).ge.qsmall .and. t(k).lt.homogfrze) then
         Q_nuc = qc(k)
         N_nuc = max(nc(k),nsmall)

         qirim(k) = qirim(k) + Q_nuc
         qitot(k) = qitot(k) + Q_nuc
         birim(k) = birim(k) + Q_nuc*inv_rho_rimeMax
         nitot(k) = nitot(k) + N_nuc
         th(k) = th(k) + exner(k)*Q_nuc*xlf(k)*inv_cp
         qc(k) = 0._rtype
         nc(k) = 0._rtype

      endif

      if (qr(k).ge.qsmall .and. t(k).lt.homogfrze) then
         Q_nuc = qr(k)
         N_nuc = max(nr(k),nsmall)

         qirim(k) = qirim(k) + Q_nuc
         qitot(k) = qitot(k) + Q_nuc
         birim(k) = birim(k) + Q_nuc*inv_rho_rimeMax
         nitot(k) = nitot(k) + N_nuc
         th(k) = th(k) + exner(k)*Q_nuc*xlf(k)*inv_cp
         qr(k) = 0._rtype
         nr(k) = 0._rtype
      endif

   enddo k_loop_fz

end subroutine homogeneous_freezing

end module micro_p3
