!__________________________________________________________________________________________
! This module contains the Predicted Particle Property (P3) bulk microphysics scheme.      !
!                                                                                          !
! This code was originally written by H. Morrison,  MMM Division, NCAR (Dec 2012).         !
! Modifications were made by J. Milbrandt, RPN, Environment Canada (July 2014).            !
! Peter Caldwell (caldwell19@llnl.gov) further modified this code to remove multiple       !
! ice categories and to clean up/simplify the code for conversion to C++ (9/11/18)         !
! Jacob Shpund (jacob.shpund@pnnl.gov) implemented and further modified/cleaned            !
! the scheme in E3SMv2                                                                     ! 
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
! Improvements from Yunpeng Shan, Jiwen Fan, Kai Zhang, and Jacob Shpund                   !
! 1) Impose low Nc limit.                                                                  !
! 2) Fix the Non-BFB issue in the calc_rime_density subroutine.                            !
! 3) Fix the bug in nc_conservation subroutine.                                            !
! 4) Output the in-cloud properties (qi, ni, qc, nc) after microphysics.                   !
! 5) Make exponential parameters in accretion/autoconversion equation consistent with      !
!    heterogeneity factor.                                                                 !
!__________________________________________________________________________________________!

! Include bit-for-bit math macros.
! #include "bfb_math.inc"

module micro_p3

   ! get real kind from utils
   use physics_utils, only: rtype,rtype8,btype

   use phys_control,  only: use_hetfrz_classnuc

   ! physical and mathematical constants
   use micro_p3_utils, only: rho_1000mb,rho_600mb,ar,br,f1r,f2r,rho_h2o,kr,kc,aimm,mi0, &
       eci,eri,bcn,cpw,cons1,cons3,cons4,cons5,cons6,cons7,         &
       inv_rho_h2o,inv_dropmass,qsmall,nsmall,cp,g,rd,rv,ep_2,inv_cp,   &
       thrd,sxth,piov6,rho_rimeMin,     &
       rho_rimeMax,inv_rho_rimeMax,max_total_ni,dbrk,nmltratio,clbfact_sub,  &
       clbfact_dep,iparam, isize, densize, rimsize, rcollsize, ice_table_size, collect_table_size, &
       get_latent_heat, T_zerodegc, pi=>pi_e3sm, dnu, &
       T_rainfrz, T_icenuc, T_homogfrz, iulog=>iulog_e3sm, &
       masterproc=>masterproc_e3sm, calculate_incloud_mixingratios, mu_r_constant, &
       lookup_table_1a_dum1_c, rho_h2o, &
       do_Cooper_inP3

   use wv_sat_scream, only:qv_sat

  ! Bit-for-bit math functions.
#ifdef SCREAM_CONFIG_IS_CMAKE
  use physics_share_f2c, only: scream_pow, scream_sqrt, scream_cbrt, scream_gamma, scream_log, &
                               scream_log10, scream_exp, scream_expm1, scream_tanh
#endif

  implicit none
  save

  public  :: p3_init,p3_main

  ! protected items should be treated as private for everyone except tests

  real(rtype), protected, dimension(densize,rimsize,isize,ice_table_size) :: ice_table_vals   !ice lookup table values

  !ice lookup table values for ice-rain collision/collection
  real(rtype), protected, dimension(densize,rimsize,isize,rcollsize,collect_table_size) :: collect_table_vals

  ! lookup table values for rain shape parameter mu_r
  real(rtype), protected, dimension(150) :: mu_r_table_vals

  ! lookup table values for rain number- and mass-weighted fallspeeds and ventilation parameters
  real(rtype), protected, dimension(300,10) :: vn_table_vals,vm_table_vals,revap_table_vals

  type realptr
     real(rtype), dimension(:), pointer :: p
  end type realptr

contains

pure function bfb_pow(base, exp) result(res)
   implicit none
   real(rtype), intent(in) :: base, exp
   real(rtype) :: res

   res = base**exp
end function bfb_pow

pure function bfb_cbrt(base) result(res)
   implicit none
   real(rtype), intent(in) :: base
   real(rtype) :: res, loc_thrd
   
   loc_thrd = 1._rtype/3._rtype
   res = base**loc_thrd
end function bfb_cbrt

pure function bfb_gamma(val) result(res)
   implicit none
   real(rtype), intent(in) :: val
   real(rtype) :: res
   
   res = gamma(val)
end function bfb_gamma

pure function bfb_log(val) result(res)
   implicit none
   real(rtype), intent(in) :: val
   real(rtype) :: res
   
   res = log(val)
end function bfb_log

pure function bfb_log10(val) result(res)
   implicit none
   real(rtype), intent(in) :: val
   real(rtype) :: res
   
   res = log10(val)
end function bfb_log10

pure function bfb_exp(val) result(res)
   implicit none
   real(rtype), intent(in) :: val
   real(rtype) :: res
   
   res = exp(val)
end function bfb_exp

pure function bfb_sqrt(val) result(res)
   implicit none
   real(rtype), intent(in) :: val
   real(rtype) :: res
   
   res = sqrt(val)
end function bfb_sqrt

pure function bfb_square(val) result(res)
   implicit none
   real(rtype), intent(in) :: val
   real(rtype) :: res
   
   res = val*val
end function bfb_square

pure function bfb_cube(val) result(res)
   implicit none
   real(rtype), intent(in) :: val
   real(rtype) :: res
   
   res = val*val*val
end function bfb_cube

pure function bfb_expm1(val) result(res)
   implicit none
   real(rtype), intent(in) :: val
   real(rtype) :: res
   
   res = exp(val) - 1.0_rtype
end function bfb_expm1

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

    ice_table_vals(:,:,:,:) = 0.
    collect_table_vals(:,:,:,:,:) = 0.
    do jj = 1,densize
       do ii = 1,rimsize
          do i = 1,isize
             read(10,*) dumi,dumi,dum,dum,dumk(1),dumk(2),              &
                  dumk(3),dumk(4),dumk(5),dumk(6),dumk(7),dumk(8),dum,  &
                  dumk(9),dumk(10),dumk(11),dumk(12)
             ice_table_vals(jj,ii,i,:) = dumk(:)
          enddo
          ! read in table for ice-rain collection
          do i = 1,isize
             do j = 1,rcollsize
                read(10,*) dumi,dumi,dum,dum,dum,dumk1,dumk2,dum
                collect_table_vals(jj,ii,i,j,1) = dlog10(dumk1)
                collect_table_vals(jj,ii,i,j,2) = dlog10(dumk2)
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
    mu_r_user(:) = mu_r_table_vals(:)
    revap_user(:,:) = revap_table_vals(:,:)
    vn_user(:,:) = vn_table_vals(:,:)
    vm_user(:,:) = vm_table_vals(:,:)

   return

  end subroutine p3_get_tables

  subroutine p3_set_tables(mu_r_user, revap_user, vn_user, vm_user)
    ! This can be called instead of p3_init_b.
    implicit none
    real(rtype), dimension(150), intent(in) :: mu_r_user
    real(rtype), dimension(300,10), intent(in) :: vn_user, vm_user, revap_user
    mu_r_table_vals(:) = mu_r_user(:)
    revap_table_vals(:,:) = revap_user(:,:)
    vn_table_vals(:,:) = vn_user(:,:)
    vm_table_vals(:,:) = vm_user(:,:)

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
       mu_r_table_vals(i) = mu_r_constant

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

          vn_table_vals(jj,ii)    = dum1/dum2
          vm_table_vals(jj,ii)    = dum3/dum4
          revap_table_vals(jj,ii) = 10._rtype**(log10(dum5)+(mu_r+1._rtype)*log10(lamr)-(3._rtype*mu_r))

       enddo meansize_loop

    enddo mu_r_loop

  END SUBROUTINE p3_init_b

  !==========================================================================================!

  SUBROUTINE p3_main_part1(kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt,                &
       pres, dpres, dz, nc_nuceat_tend, nccn_prescribed, exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i,  &
       inv_cld_frac_r, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion,                           &
       t_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn, qv, th_atm,           &
       qc, nc, qr, nr,                                                                                      &
       qi, ni, qm, bm, qc_incld, qr_incld, qi_incld, qm_incld,                                              &
       nc_incld, nr_incld, ni_incld, bm_incld, is_nucleat_possible, is_hydromet_present, nccnst)

    implicit none

    ! args

    integer, intent(in) :: kts, kte, kbot, ktop, kdir
    logical(btype), intent(in) :: do_predict_nc
    real(rtype), intent(in) :: dt, nccnst

    real(rtype), intent(in), dimension(kts:kte) :: pres, dpres, dz, nc_nuceat_tend, exner, inv_exner, &
         inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, nccn_prescribed

    real(rtype), intent(inout), dimension(kts:kte) :: t_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, &
         acn, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, qc_incld, qr_incld, qi_incld, &
         qm_incld, nc_incld, nr_incld, ni_incld, bm_incld

    logical(btype), intent(out) :: is_nucleat_possible, is_hydromet_present
    logical(btype), intent(in) :: do_prescribed_CCN

    ! locals
    integer :: k
    real(rtype) :: dum

    is_nucleat_possible = .false.
    is_hydromet_present = .false.

    k_loop_1: do k = kbot,ktop,kdir
       !calculate some time-varying atmospheric variables
       !AaronDonahue - changed "rho" to be defined on nonhydrostatic
       !assumption, consistent with pressure based coordinate system
       !             - moved latent heat calculation to above.  Latent
       !heat is determined by calling a p3_util function so that it
       !can be made consistent with E3SM definition of latent heat
       rho(k)     = dpres(k)/dz(k)/g  ! pres(k)/(rd*t(k))
       inv_rho(k) = 1._rtype/rho(k)
       qv_sat_l(k)     = qv_sat(t_atm(k),pres(k),0)
       qv_sat_i(k)     = qv_sat(t_atm(k),pres(k),1)

       qv_supersat_i(k)    = qv(k)/qv_sat_i(k)-1._rtype

       rhofacr(k) = bfb_pow(rho_1000mb*inv_rho(k), 0.54_rtype)
       rhofaci(k) = bfb_pow(rho_600mb*inv_rho(k), 0.54_rtype)
       dum        = 1.496e-6_rtype * bfb_pow(t_atm(k), 1.5_rtype) / (t_atm(k)+120._rtype)  ! this is mu
       acn(k)     = g*rho_h2o/(18._rtype*dum)  ! 'a' parameter for droplet fallspeed (Stokes' law)

       if ((t_atm(k).lt.T_zerodegc .and. qv_supersat_i(k).ge.-0.05_rtype)) is_nucleat_possible = .true.

       if (qc(k).lt.qsmall) then
      !--- apply mass clipping if mass is sufficiently small
      !    (implying all mass is expected to evaporate/sublimate in one time step)
          qv(k) = qv(k) + qc(k)
          th_atm(k) = th_atm(k) - exner(k)*qc(k)*latent_heat_vapor(k)*inv_cp
          qc(k) = 0._rtype
          nc(k) = 0._rtype
       else
          is_hydromet_present = .true.    ! updated further down
      !--- Apply droplet activation here (before other microphysical processes) for consistency with qc increase by saturation
      !    adjustment already applied in macrophysics. If prescribed drop number is used, this is also a good place to
      !    prescribe that value

          if (do_prescribed_CCN) then
             nc(k) = max(nc(k),nccn_prescribed(k))
          else if (do_predict_nc) then
             nc(k) = max(nc(k) + nc_nuceat_tend(k) * dt,0.0_rtype)
          else
             nc(k) = nccnst*inv_rho(k)
          endif
       endif

       if (qr(k).lt.qsmall) then
          qv(k) = qv(k) + qr(k)
          th_atm(k) = th_atm(k) - exner(k)*qr(k)*latent_heat_vapor(k)*inv_cp
          qr(k) = 0._rtype
          nr(k) = 0._rtype
       else
          is_hydromet_present = .true.    ! updated further down
       endif

       if (qi(k).lt.qsmall .or. (qi(k).lt.1.e-8_rtype .and.             &
            qv_supersat_i(k).lt.-0.1_rtype)) then
          qv(k) = qv(k) + qi(k)
          th_atm(k) = th_atm(k) - exner(k)*qi(k)*latent_heat_sublim(k)*inv_cp
          qi(k) = 0._rtype
          ni(k) = 0._rtype
          qm(k) = 0._rtype
          bm(k) = 0._rtype
       else
          is_hydromet_present = .true.    ! final update
       endif

       if (qi(k).ge.qsmall .and. qi(k).lt.1.e-8_rtype .and.             &
            t_atm(k).ge.T_zerodegc) then
          qr(k) = qr(k) + qi(k)
          th_atm(k) = th_atm(k) - exner(k)*qi(k)*latent_heat_fusion(k)*inv_cp
          qi(k) = 0._rtype
          ni(k) = 0._rtype
          qm(k) = 0._rtype
          bm(k) = 0._rtype
       endif

       t_atm(k) = th_atm(k) * inv_exner(k)

       call calculate_incloud_mixingratios(qc(k),qr(k),qi(k),qm(k),nc(k),nr(k),ni(k),bm(k), &
            inv_cld_frac_l(k),inv_cld_frac_i(k),inv_cld_frac_r(k), &
            qc_incld(k),qr_incld(k),qi_incld(k),qm_incld(k),nc_incld(k),nr_incld(k),ni_incld(k),bm_incld(k))

    enddo k_loop_1

  END SUBROUTINE p3_main_part1

  SUBROUTINE p3_main_part2(kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt, inv_dt, &
       p3_autocon_coeff,p3_accret_coeff,p3_qc_autocon_expon,p3_nc_autocon_expon,p3_qc_accret_expon, &
       p3_wbf_coeff,p3_embryonic_rain_size, p3_max_mean_rain_size, &
       pres, dpres, dz, nc_nuceat_tend, exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, ni_activated, &
       inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r, qv_prev, t_prev, &
       t_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn, qv, th_atm, qc, nc, qr, nr, qi, ni, &
       qm, bm, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, &
       ni_incld, bm_incld, mu_c, nu, lamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, qv2qi_depos_tend, precip_total_tend, &
       nevapr, qr_evap_tend, vap_liq_exchange, vap_ice_exchange, liq_ice_exchange, pratot, &
       prctot, frzimm, frzcnt, frzdep, p3_tend_out, is_hydromet_present, do_precip_off, nccnst)

    implicit none

    ! args

    integer, intent(in) :: kts, kte, kbot, ktop, kdir
    logical(btype), intent(in) :: do_predict_nc, do_prescribed_CCN, do_precip_off
    real(rtype), intent(in) :: dt, inv_dt, nccnst
    real(rtype), intent(in) :: p3_autocon_coeff, p3_accret_coeff, p3_qc_autocon_expon, p3_nc_autocon_expon, p3_qc_accret_expon, &
         p3_wbf_coeff, p3_embryonic_rain_size, p3_max_mean_rain_size

    real(rtype), intent(in), dimension(kts:kte) :: pres, dpres, dz, nc_nuceat_tend, exner, inv_exner, inv_cld_frac_l,      &
         inv_cld_frac_i, inv_cld_frac_r, ni_activated, inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r, qv_prev, t_prev, &
         frzimm, frzcnt, frzdep

    real(rtype), intent(inout), dimension(kts:kte) :: t_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn,   &
         qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, qc_incld, qr_incld,    &
         qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld, mu_c, nu, lamc, cdist, cdist1,                                    &
         cdistr, mu_r, lamr, logn0r, qv2qi_depos_tend, precip_total_tend, nevapr, qr_evap_tend, vap_liq_exchange,                      &
         vap_ice_exchange, liq_ice_exchange, pratot, prctot

    real(rtype), intent(inout), dimension(kts:kte,49) :: p3_tend_out ! micro physics tendencies

    logical(btype), intent(out) :: is_hydromet_present

    ! -------- locals ------- !

    ! liquid-phase microphysical process rates:
    !  (all Q process rates in kg kg-1 s-1)
    !  (all N process rates in # kg-1)

    real(rtype) :: qc2qr_accret_tend   ! cloud droplet accretion by rain
    real(rtype) :: qc2qr_autoconv_tend   ! cloud droplet autoconversion to rain
    real(rtype) :: nc_accret_tend   ! change in cloud droplet number from accretion by rain
    real(rtype) :: nc2nr_autoconv_tend  ! change in cloud droplet number from autoconversion
    real(rtype) :: nc_selfcollect_tend   ! change in cloud droplet number from self-collection  (Not in paper?)
    real(rtype) :: nr_selfcollect_tend   ! change in rain number from self-collection  (Not in paper?)
    real(rtype) :: qr2qv_evap_tend   ! rain evaporation
    real(rtype) :: nr_evap_tend   ! change in rain number from evaporation
    real(rtype) :: ncautr  ! change in rain number from autoconversion of cloud water
! Is is assumed that macrophysics handles condensation/evaporation of qc and
! that there is no condensation of rain. Thus qccon, qrcon and qcevp have
! been removed from the original P3-WRF for P3-SCREAM.
!    real(rtype) :: qrcon   ! rain condensation   (Not in paper?)
!    real(rtype) :: qccon   ! cloud droplet condensation
!    real(rtype) :: qcevp   ! cloud droplet evaporation

    ! ice-phase microphysical process rates:
    !  (all Q process rates in kg kg-1 s-1)
    !  (all N process rates in # kg-1 s-1)

    real(rtype) :: qccol     ! collection of cloud water by ice
    real(rtype) :: qwgrth    ! wet growth rate
    real(rtype) :: qidep     ! vapor deposition
    real(rtype) :: qrcol     ! collection rain mass by ice
    real(rtype) :: qinuc     ! deposition/condensation freezing nuc
    real(rtype) :: nc_collect_tend     ! change in cloud droplet number from collection by ice
    real(rtype) :: nr_collect_tend     ! change in rain number from collection by ice
    real(rtype) :: ni_nucleat_tend     ! change in ice number from deposition/cond-freezing nucleation
    real(rtype) :: qi2qv_sublim_tend     ! sublimation of ice
    real(rtype) :: qi2qr_melt_tend     ! melting of ice
    real(rtype) :: ni2nr_melt_tend     ! melting of ice
    real(rtype) :: ni_sublim_tend     ! change in ice number from sublimation
    real(rtype) :: ni_selfcollect_tend     ! change in ice number from collection within a category (Not in paper?)
    real(rtype) :: qc2qi_hetero_freeze_tend    ! immersion freezing droplets
    real(rtype) :: qr2qi_immers_freeze_tend    ! immersion freezing rain
    real(rtype) :: nc2ni_immers_freeze_tend    ! immersion freezing droplets
    real(rtype) :: nr2ni_immers_freeze_tend    ! immersion freezing rain
    real(rtype) :: nr_ice_shed_tend    ! source for rain number from collision of rain/ice above freezing and shedding
    real(rtype) :: qc2qr_ice_shed_tend     ! source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding
    real(rtype) :: rho_qm_cloud ! density of rime (from cloud)
    real(rtype) :: ncshdc    ! source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)
    real(rtype) :: qiberg    ! Bergeron process

    real(rtype)    :: table_val_qi_fallspd   ! mass-weighted fallspeed              See lines  731 -  808  ums
    real(rtype)    :: table_val_ni_self_collect   ! ice collection within a category     See lines  809 -  928  nagg
    real(rtype)    :: table_val_qc2qi_collect   ! collection of cloud water by ice     See lines  929 - 1009  nrwat
    real(rtype)    :: table_val_qi2qr_melting   ! melting                              See lines 1212 - 1279  vdep
    real(rtype)    :: table_val_nr_collect   ! collection of rain number by ice     See lines 1010 - 1209  nrrain
    real(rtype)    :: table_val_qr2qi_collect   ! collection of rain mass by ice       See lines 1010 - 1209  qrrain
    real(rtype)    :: table_val_ni_lammax   ! minimum ice number (lambda limiter)  See lines  704 -  705  nlarge
    real(rtype)    :: table_val_ni_lammin   ! maximum ice number (lambda limiter)  See lines  704 -  705  nsmall
    real(rtype)    :: table_val_qi2qr_vent_melt   ! melting (ventilation term)           See lines 1212 - 1279  vdep1

    real(rtype)    :: mu,dv,sc,dqsdt,ab,kap,epsr,epsc,epsi,epsi_tot, &
         dum1,dum3,dum4,dum5,dum6,dqsidt,abi,rhop,vtrmi1,eii
      
    real(rtype) :: ncheti_cnt,qcheti_cnt,nicnt,qicnt,ninuc_cnt,qinuc_cnt,qi_wetDepos     

    integer :: dumi,k,dumj,dumii,dumjj,dumzz

    logical(btype) :: log_exitlevel, log_wetgrowth

   rho_qm_cloud = 400._rtype
   is_hydromet_present = .false.

   !------------------------------------------------------------------------------------------!
   !   main k-loop (for processes):
   k_loop_main: do k = kbot,ktop,kdir

      ! if relatively dry and no hydrometeors at this level, skip to end of k-loop (i.e. skip this level)
      log_exitlevel = .true.
      if (qc(k).ge.qsmall .or. qr(k).ge.qsmall) log_exitlevel = .false.

      if (qi(k).ge.qsmall) log_exitlevel = .false.
      !enddo
      if (log_exitlevel .and.                                                           &
         (t_atm(k).lt.T_zerodegc .and. qv_supersat_i(k).lt.-0.05_rtype)) goto 555   !i.e. skip all process rates

      ! All microphysics tendencies will be computed as IN-CLOUD, they will be mapped back to cell-average later.

      ! initialize warm-phase process rates
      qc2qr_accret_tend   = 0._rtype;     qr2qv_evap_tend   = 0._rtype;     qc2qr_autoconv_tend   = 0._rtype;
      nc_accret_tend   = 0._rtype;     nc_selfcollect_tend   = 0._rtype;
      nc2nr_autoconv_tend  = 0._rtype;     nr_selfcollect_tend   = 0._rtype;
      nr_evap_tend   = 0._rtype;     ncautr  = 0._rtype

      ! initialize ice-phase  process rates
      qi2qv_sublim_tend   = 0._rtype;     nr_ice_shed_tend  = 0._rtype
      qc2qi_hetero_freeze_tend  = 0._rtype;     qrcol   = 0._rtype;     qc2qr_ice_shed_tend   = 0._rtype
      qi2qr_melt_tend   = 0._rtype;     qccol   = 0._rtype
      qr2qi_immers_freeze_tend  = 0._rtype;     qinuc   = 0._rtype;     ni2nr_melt_tend   = 0._rtype
      nc_collect_tend   = 0._rtype;     ncshdc  = 0._rtype
      nc2ni_immers_freeze_tend  = 0._rtype;     nr_collect_tend   = 0._rtype;     ni_selfcollect_tend   = 0._rtype
      ni_nucleat_tend   = 0._rtype;     qidep   = 0._rtype;     qiberg  = 0._rtype
      nr2ni_immers_freeze_tend  = 0._rtype;     ni_sublim_tend   = 0._rtype;     qwgrth  = 0._rtype

      ncheti_cnt = 0.0_rtype; qcheti_cnt = 0.0_rtype; nicnt = 0.0_rtype 
      qicnt = 0.0_rtype;      ninuc_cnt = 0.0_rtype;  qinuc_cnt = 0.0_rtype
      qi_wetDepos = 0.0_rtype


      log_wetgrowth = .false.

      ! skip micro process calculations except nucleation/acvtivation if there no hydrometeors are present
      log_exitlevel = .true.
      if (qc_incld(k).ge.qsmall .or. qr_incld(k).ge.qsmall) log_exitlevel = .false.
      if (qi_incld(k).ge.qsmall) log_exitlevel=.false.
      if (log_exitlevel) goto 444   !i.e. skip to nucleation

      !time/space varying physical variables
      call get_time_space_phys_variables( &
           t_atm(k),pres(k),rho(k),latent_heat_vapor(k),latent_heat_sublim(k),qv_sat_l(k),qv_sat_i(k), &
           mu,dv,sc,dqsdt,dqsidt,ab,abi,kap,eii)

      call get_cloud_dsd2(qc_incld(k),nc_incld(k),mu_c(k),rho(k),nu(k),dnu,lamc(k),     &
           cdist(k),cdist1(k))
      nc(k) = nc_incld(k)*cld_frac_l(k)

      call get_rain_dsd2(qr_incld(k),p3_max_mean_rain_size,nr_incld(k),mu_r(k),lamr(k),   &
           cdistr(k),logn0r(k))
      nr(k) = nr_incld(k)*cld_frac_r(k)

      ! initialize inverse supersaturation relaxation timescale for combined ice categories
      epsi_tot = 0._rtype

      call impose_max_total_ni(ni_incld(k),max_total_ni,inv_rho(k))

      if (qi_incld(k).ge.qsmall) then

         !impose lower limits to prevent taking log of # < 0
         ni_incld(k) = max(ni_incld(k),nsmall)
         nr_incld(k)    = max(nr_incld(k),nsmall)

         call calc_bulkRhoRime(qi_incld(k),qm_incld(k),bm_incld(k),rhop)
         qm(k)=qm_incld(k)*cld_frac_i(k)
         bm(k)=bm_incld(k)*cld_frac_i(k)

         ! if (.not. tripleMoment_on) zitot(k) = diag_mom6(qi_incld(k),ni_incld(k),rho(k))
         call find_lookupTable_indices_1a(dumi,dumjj,dumii,dumzz,dum1,dum4,   &
              dum5,dum6,isize,rimsize,densize,                                &
              qi_incld(k),ni_incld(k),qm_incld(k),                            &
              rhop)
         !qm_incld(k),zitot(k),rhop)
         call find_lookupTable_indices_1b(dumj,dum3,rcollsize,qr_incld(k),nr_incld(k))

         ! call to lookup table interpolation subroutines to get process rates
         call access_lookup_table(dumjj,dumii,dumi, 2,dum1,dum4,dum5,table_val_qi_fallspd)
         call access_lookup_table(dumjj,dumii,dumi, 3,dum1,dum4,dum5,table_val_ni_self_collect)
         call access_lookup_table(dumjj,dumii,dumi, 4,dum1,dum4,dum5,table_val_qc2qi_collect)
         call access_lookup_table(dumjj,dumii,dumi, 5,dum1,dum4,dum5,table_val_qi2qr_melting)
         call access_lookup_table(dumjj,dumii,dumi, 7,dum1,dum4,dum5,table_val_ni_lammax)
         call access_lookup_table(dumjj,dumii,dumi, 8,dum1,dum4,dum5,table_val_ni_lammin)
         call access_lookup_table(dumjj,dumii,dumi,10,dum1,dum4,dum5,table_val_qi2qr_vent_melt)

         ! ice-rain collection processes
         if (qr_incld(k).ge.qsmall) then
            call access_lookup_table_coll(dumjj,dumii,dumj,dumi,1,dum1,dum3,dum4,dum5,table_val_nr_collect)
            call access_lookup_table_coll(dumjj,dumii,dumj,dumi,2,dum1,dum3,dum4,dum5,table_val_qr2qi_collect)
         else
            table_val_nr_collect = 0._rtype
            table_val_qr2qi_collect = 0._rtype
         endif

         ! adjust Ni if needed to make sure mean size is in bounds (i.e. apply lambda limiters)
         ! note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
         ni_incld(k) = min(ni_incld(k),table_val_ni_lammax*ni_incld(k))
         ni_incld(k) = max(ni_incld(k),table_val_ni_lammin*ni_incld(k))

      endif   ! qi > qsmall

      !----------------------------------------------------------------------
      ! Begin calculations of microphysical processes

      !......................................................................
      ! ice processes
      !......................................................................

      !.......................
      ! collection of droplets
      call ice_cldliq_collection(rho(k),t_atm(k),rhofaci(k),&
           table_val_qc2qi_collect,qi_incld(k),qc_incld(k),ni_incld(k),nc_incld(k),&
           qccol,nc_collect_tend,qc2qr_ice_shed_tend,ncshdc)

      !....................
      ! collection of rain
      call ice_rain_collection(rho(k),t_atm(k),rhofaci(k),&
           logn0r(k),table_val_nr_collect,table_val_qr2qi_collect,qi_incld(k),ni_incld(k),qr_incld(k),&
           qrcol,nr_collect_tend)
      !...................................
      ! collection between ice categories

      !PMC nCat deleted lots of stuff here.

      !.............................................
      ! self-collection of ice
      call ice_self_collection(rho(k),rhofaci(k),&
           table_val_ni_self_collect,eii,qm_incld(k),qi_incld(k),ni_incld(k),&
           ni_selfcollect_tend)

      !............................................................
      ! melting
      call ice_melting(rho(k),t_atm(k),pres(k),rhofaci(k),&
           table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor(k),latent_heat_fusion(k),dv,sc,mu,kap,&
           qv(k),qi_incld(k),ni_incld(k),&
           qi2qr_melt_tend,ni2nr_melt_tend)

      !............................................................
      ! calculate wet growth
      call ice_cldliq_wet_growth(rho(k),t_atm(k),pres(k),rhofaci(k),&
           table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor(k),latent_heat_fusion(k),dv,kap,mu,sc,&
           qv(k),qc_incld(k),qi_incld(k),ni_incld(k),qr_incld(k),log_wetgrowth,&
           qrcol,qccol,qwgrth,nr_ice_shed_tend,qc2qr_ice_shed_tend)

      !-----------------------------
      ! calcualte total inverse ice relaxation timescale combined for all ice categories
      ! note 'f1pr' values are normalized, so we need to multiply by N
      call calc_ice_relaxation_timescale(rho(k),t_atm(k),rhofaci(k),&
           table_val_qi2qr_melting,table_val_qi2qr_vent_melt,dv,mu,sc,qi_incld(k),ni_incld(k),&
           epsi,epsi_tot)

      !.........................
      ! calculate rime density    
      call calc_rime_density(t_atm(k),rhofaci(k),&
           table_val_qi_fallspd,acn(k),lamc(k),mu_c(k),qc_incld(k),qccol,&
           vtrmi1,rho_qm_cloud)                      
      !............................................................
      ! contact and immersion freezing droplets
      if(.not. use_hetfrz_classnuc)then     
         call cldliq_immersion_freezing(t_atm(k),&
                                       lamc(k),mu_c(k),cdist1(k),qc_incld(k),inv_qc_relvar(k),&
                                       qc2qi_hetero_freeze_tend,nc2ni_immers_freeze_tend)
      else
         call CNT_couple (frzimm(k),frzcnt(k),frzdep(k),rho(k),qc_incld(k),nc_incld(k),1, & 
                           ncheti_cnt,qcheti_cnt,nicnt,qicnt,ninuc_cnt,qinuc_cnt,inv_dt)
         call CNT_couple (frzimm(k),frzcnt(k),frzdep(k),rho(k),qc_incld(k),nc_incld(k),2, & 
                           ncheti_cnt,qcheti_cnt,nicnt,qicnt,ninuc_cnt,qinuc_cnt,inv_dt)
      endif

      !............................................................
      ! immersion freezing of rain
      ! for future: get rid of log statements below for rain freezing
      call rain_immersion_freezing(t_atm(k),&
           lamr(k),mu_r(k),cdistr(k),qr_incld(k),&
           qr2qi_immers_freeze_tend,nr2ni_immers_freeze_tend)

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

      call evaporate_rain(qr_incld(k),qc_incld(k),nr_incld(k),qi_incld(k), &
           cld_frac_l(k),cld_frac_r(k),qv(k),qv_prev(k),qv_sat_l(k),qv_sat_i(k), &
           ab,abi,epsr,epsi_tot,T_atm(k),t_prev(k),latent_heat_sublim(k),dqsdt,dt,&
           qr2qv_evap_tend,nr_evap_tend)

      call ice_deposition_sublimation(qi_incld(k), ni_incld(k), t_atm(k), &
           qv_sat_l(k),qv_sat_i(k),epsi,abi,qv(k), p3_wbf_coeff, &
           qidep,qi2qv_sublim_tend,ni_sublim_tend,qiberg)

444   continue

      !................................................................
      ! deposition/condensation-freezing nucleation
      call ice_nucleation(t_atm(k),inv_rho(k),&
           ni(k),ni_activated(k),qv_supersat_i(k),inv_dt,do_predict_nc, do_prescribed_CCN, &
           qinuc, ni_nucleat_tend)

      !................
      ! cloud water autoconversion
      ! NOTE: cloud_water_autoconversion must be called before droplet_self_collection
      call cloud_water_autoconversion(rho(k),qc_incld(k),nc_incld(k),inv_qc_relvar(k),&
           p3_autocon_coeff,p3_qc_autocon_expon,p3_nc_autocon_expon,p3_embryonic_rain_size,&
           do_precip_off,qc2qr_autoconv_tend,nc2nr_autoconv_tend,ncautr)

      !............................
      ! self-collection of droplets
      call droplet_self_collection(rho(k),inv_rho(k),qc_incld(k),&
           mu_c(k),nu(k),nc2nr_autoconv_tend,nc_selfcollect_tend)

      !............................
      ! accretion of cloud by rain
      call cloud_rain_accretion(rho(k),inv_rho(k),&
           qc_incld(k),nc_incld(k), qr_incld(k),inv_qc_relvar(k),&
           p3_accret_coeff,p3_qc_accret_expon,&
           qc2qr_accret_tend, nc_accret_tend)

      !.....................................
      ! self-collection and breakup of rain
      ! (breakup following modified Verlinde and Cotton scheme)
      call rain_self_collection(rho(k),qr_incld(k),nr_incld(k),&
           nr_selfcollect_tend)

      ! Here we map the microphysics tendency rates back to CELL-AVERAGE quantities for updating
      ! cell-average quantities.
      call back_to_cell_average(cld_frac_l(k),cld_frac_r(k),cld_frac_i(k), qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend, &
           nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend, nr_selfcollect_tend, nr_evap_tend, ncautr,                     &
           qi2qv_sublim_tend, nr_ice_shed_tend, qc2qi_hetero_freeze_tend,                                                           &
           qrcol, qc2qr_ice_shed_tend, qi2qr_melt_tend, qccol, qr2qi_immers_freeze_tend, ni2nr_melt_tend, nc_collect_tend,          &
           ncshdc, nc2ni_immers_freeze_tend, nr_collect_tend, ni_selfcollect_tend,                                                  &
           qidep, nr2ni_immers_freeze_tend, ni_sublim_tend, qinuc, ni_nucleat_tend, qiberg,                                         &
           ncheti_cnt, qcheti_cnt, nicnt, qicnt, ninuc_cnt, qinuc_cnt)

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
      ! how qidep and qi2qv_sublim_tend are derived?
      ! AaronDonahue: UPDATE, if we are following the implementation of MG
      ! then the answer appears to be YES.  There is a similar check in MG
      ! microphysics which limits qidep and qinuc, but does not limit qi2qv_sublim_tend.
      ! So similar but slightly different.  The overall answer would be that
      ! qidep does need some limit.  The next questions are,
      !   1) Should we be taking qinuc into consideration too?
      !   2) Is MG correct in NOT limiting qi2qv_sublim_tend?



      
      ! cloud mass
      call cloud_water_conservation(qc(k), dt, qc2qr_autoconv_tend, qc2qr_accret_tend, qccol, qc2qi_hetero_freeze_tend, &
           qc2qr_ice_shed_tend, qiberg, qi2qv_sublim_tend, qidep, qcheti_cnt, qicnt)

      ! rain mass
      call rain_water_conservation(qr(k), qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, dt, &
           qr2qv_evap_tend, qrcol, qr2qi_immers_freeze_tend)

      ! ice mass
      call ice_water_conservation(qi(k), qidep, qinuc, qiberg, qrcol, qccol, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, &
           dt, qinuc_cnt, qcheti_cnt, qicnt, qi2qv_sublim_tend, qi2qr_melt_tend)

      ! cloud number     
      call nc_conservation(nc(k), nc_selfcollect_tend, dt, nc_collect_tend, nc2ni_immers_freeze_tend, &
           nc_accret_tend, nc2nr_autoconv_tend, ncheti_cnt, nicnt)

      ! rain number     
      call nr_conservation(nr(k),ni2nr_melt_tend,nr_ice_shed_tend,ncshdc,nc2nr_autoconv_tend,dt,nr_collect_tend,nmltratio, &
           nr2ni_immers_freeze_tend,nr_selfcollect_tend,nr_evap_tend)
      
      ! ice number     
      call ni_conservation(ni(k),ni_nucleat_tend,nr2ni_immers_freeze_tend,nc2ni_immers_freeze_tend,ncheti_cnt,nicnt,ninuc_cnt,dt,ni2nr_melt_tend,&
           ni_sublim_tend,ni_selfcollect_tend)

      call prevent_ice_overdepletion(pres(k), t_atm(k), qv(k), latent_heat_vapor(k), latent_heat_sublim(k), inv_dt, dt, qidep, qinuc, qinuc_cnt, qi2qv_sublim_tend, qr2qv_evap_tend)

      !call water_vapor_conservation(qv(k), qidep, qinuc, qi2qv_sublim_tend, qr2qv_evap_tend, qinuc_cnt, dt)
           
      call ice_supersat_conservation(qidep, qinuc, qi2qv_sublim_tend, qr2qv_evap_tend, qinuc_cnt, cld_frac_i(k), qv(k), qv_sat_i(k), latent_heat_sublim(k), th_atm(k)/exner(k), dt)

      !---------------------------------------------------------------------------------
      ! update prognostic microphysics and thermodynamics variables
      !---------------------------------------------------------------------------------

      !-- ice-phase dependent processes:
      call update_prognostic_ice(qc2qi_hetero_freeze_tend, qccol, qc2qr_ice_shed_tend,       &
           nc_collect_tend, nc2ni_immers_freeze_tend, ncshdc,                                &
           qrcol, nr_collect_tend,  qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend, nr_ice_shed_tend,            &
           qi2qr_melt_tend, ni2nr_melt_tend, qi2qv_sublim_tend, qidep, qinuc, ni_nucleat_tend, ni_selfcollect_tend,  &
           ni_sublim_tend, qiberg, exner(k), latent_heat_sublim(k), latent_heat_fusion(k),   &
           do_predict_nc, log_wetgrowth, dt, nmltratio, rho_qm_cloud,                        &
           ncheti_cnt, nicnt, ninuc_cnt, qcheti_cnt, qicnt, qinuc_cnt,                       &
           th_atm(k), qv(k), qi(k), ni(k), qm(k), bm(k), qc(k), nc(k), qr(k), nr(k), qi_wetDepos)

      !-- warm-phase only processes:
      call update_prognostic_liquid(qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, &
           nc_selfcollect_tend, qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend,           &
           do_predict_nc, nccnst, do_prescribed_CCN, inv_rho(k), exner(k), latent_heat_vapor(k), dt,                     &
           th_atm(k), qv(k), qc(k), nc(k), qr(k), nr(k))

      !==
      ! AaronDonahue - Add extra variables needed from microphysics by E3SM:
      if(.not. use_hetfrz_classnuc)then     
         qv2qi_depos_tend(k) = qidep - qi2qv_sublim_tend + qinuc
      else
         qv2qi_depos_tend(k) = qidep - qi2qv_sublim_tend + qinuc + qinuc_cnt
      endif
      precip_total_tend(k)   = ( qc2qr_accret_tend + qc2qr_autoconv_tend + qc2qr_ice_shed_tend + qi_wetDepos )
      nevapr(k)  = qi2qv_sublim_tend + qr2qv_evap_tend
      qr_evap_tend(k) = qr2qv_evap_tend
      vap_ice_exchange(k) = qidep - qi2qv_sublim_tend + qinuc
      vap_liq_exchange(k) = - qr2qv_evap_tend
      liq_ice_exchange(k) = qc2qi_hetero_freeze_tend + qr2qi_immers_freeze_tend - qi2qr_melt_tend + qiberg + qccol + qrcol

      ! clipping for small hydrometeor values
      if (qc(k).lt.qsmall) then
         qv(k) = qv(k) + qc(k)
         th_atm(k) = th_atm(k) - exner(k)*qc(k)*latent_heat_vapor(k)*inv_cp
         qc(k) = 0._rtype
         nc(k) = 0._rtype
      else
         is_hydromet_present = .true.
      endif

      if (qr(k).lt.qsmall) then
         qv(k) = qv(k) + qr(k)
         th_atm(k) = th_atm(k) - exner(k)*qr(k)*latent_heat_vapor(k)*inv_cp
         qr(k) = 0._rtype
         nr(k) = 0._rtype
      else
         is_hydromet_present = .true.
      endif

      if (qi(k).lt.qsmall) then
         qv(k) = qv(k) + qi(k)
         th_atm(k) = th_atm(k) - exner(k)*qi(k)*latent_heat_sublim(k)*inv_cp
         qi(k) = 0._rtype
         ni(k) = 0._rtype
         qm(k) = 0._rtype
         bm(k) = 0._rtype
      else
         is_hydromet_present = .true.
      endif

      !impose_max_total_ni is meant to operate on in-cloud vals. ni_incld is an output of
      !calculate_incloud_mixingratios below but we need to generate it earlier for impose_max_total_ni
      ni_incld(k)=ni(k)/cld_frac_i(k)
      call impose_max_total_ni(ni_incld(k),max_total_ni,inv_rho(k))
      ni(k)=ni_incld(k)*cld_frac_i(k)

      ! Record microphysics tendencies for output:
      ! warm-phase process rates
      p3_tend_out(k, 2) = qc2qr_accret_tend           ! cloud droplet accretion by rain
      p3_tend_out(k, 3) = qc2qr_autoconv_tend         ! cloud droplet autoconversion to rain
      p3_tend_out(k, 4) = nc_accret_tend              ! change in cloud droplet number from accretion by rain
      p3_tend_out(k, 5) = nc2nr_autoconv_tend         ! change in cloud droplet number from autoconversion
      p3_tend_out(k, 6) = nc_selfcollect_tend         ! change in cloud droplet number from self-collection  (Not in paper?)
      p3_tend_out(k, 7) = nr_selfcollect_tend         ! change in rain number from self-collection  (Not in paper?)
      !p3_tend_out(k, 8) = nc_nuceat_tend             ! [in 'part1'] change in cld droplet number concentration from macrophysics
      p3_tend_out(k,11) = qr2qv_evap_tend             ! rain evaporation
      p3_tend_out(k,13) = nr_evap_tend                ! change in rain number from evaporation
      p3_tend_out(k,14) = ncautr                      ! change in rain number from autoconversion of cloud water
      ! ice-phase  process rates
      p3_tend_out(k,15) = qccol                       ! collection of cloud water by ice
      p3_tend_out(k,16) = qwgrth                      ! wet growth rate
      p3_tend_out(k,17) = qidep                       ! vapor deposition
      p3_tend_out(k,18) = qrcol                       ! collection rain mass by ice
      p3_tend_out(k,19) = qinuc                       ! deposition/condensation freezing nuc
      p3_tend_out(k,20) = nc_collect_tend             ! change in cloud droplet number from collection by ice
      p3_tend_out(k,21) = nr_collect_tend             ! change in rain number from collection by ice
      p3_tend_out(k,22) = ni_nucleat_tend             ! change in ice number from deposition/cond-freezing nucleation
      p3_tend_out(k,23) = qi2qv_sublim_tend           ! sublimation of ice
      p3_tend_out(k,24) = qi2qr_melt_tend             ! melting of ice
      p3_tend_out(k,25) = ni2nr_melt_tend             ! melting of ice
      p3_tend_out(k,26) = ni_sublim_tend              ! change in ice number from sublimation
      p3_tend_out(k,27) = ni_selfcollect_tend         ! change in ice number from collection within a category (Not in paper?)
      p3_tend_out(k,28) = qc2qi_hetero_freeze_tend    ! immersion freezing droplets
      p3_tend_out(k,29) = qr2qi_immers_freeze_tend    ! immersion freezing rain
      p3_tend_out(k,30) = nc2ni_immers_freeze_tend    ! immersion freezing droplets
      p3_tend_out(k,31) = nr2ni_immers_freeze_tend    ! immersion freezing rain
      p3_tend_out(k,32) = nr_ice_shed_tend            ! source for rain number from collision of rain/ice above freezing and shedding
      p3_tend_out(k,33) = qc2qr_ice_shed_tend         ! source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding
      p3_tend_out(k,34) = 0._rtype                    ! used to be qcmul, but that has been removed.  Kept at 0.0 as placeholder.
      p3_tend_out(k,35) = ncshdc                      ! source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)
      ! Outputs associated with aerocom comparison:
      pratot(k) = qc2qr_accret_tend                   ! cloud drop accretion by rain
      prctot(k) = qc2qr_autoconv_tend                 ! cloud drop autoconversion to rain
      !---------------------------------------------------------------------------------

      ! Recalculate in-cloud values for sedimentation
      call calculate_incloud_mixingratios(qc(k),qr(k),qi(k),qm(k),nc(k),nr(k),ni(k),bm(k), &
           inv_cld_frac_l(k),inv_cld_frac_i(k),inv_cld_frac_r(k), &
           qc_incld(k),qr_incld(k),qi_incld(k),qm_incld(k),nc_incld(k),nr_incld(k),ni_incld(k),bm_incld(k))

555   continue

   enddo k_loop_main

 END SUBROUTINE p3_main_part2

 subroutine p3_main_part3(kts, kte, kbot, ktop, kdir, p3_max_mean_rain_size,                                                &
      exner, cld_frac_l, cld_frac_r, cld_frac_i,                                                                           &
      rho, inv_rho, rhofaci, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim,            &
      mu_c, nu, lamc, mu_r, lamr, vap_liq_exchange,                                                                        &
      ze_rain, ze_ice, diag_vm_qi, diag_eff_radius_qi, diag_diam_qi, rho_qi, diag_equiv_reflectivity, diag_eff_radius_qc,  & 
      diag_ze_rain,diag_ze_ice,mincdnc)

   implicit none

   ! args

   integer, intent(in) :: kts, kte, kbot, ktop, kdir

   real(rtype), intent(in) :: p3_max_mean_rain_size
   real(rtype), intent(in) :: mincdnc 
   real(rtype), intent(in), dimension(kts:kte) :: exner, cld_frac_l, cld_frac_r, cld_frac_i

   real(rtype), intent(inout), dimension(kts:kte) :: rho, inv_rho, rhofaci, &
        qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, &
        mu_c, nu, lamc, mu_r, &
        lamr, vap_liq_exchange, &
        ze_rain, ze_ice, diag_vm_qi, diag_eff_radius_qi, diag_diam_qi, rho_qi, diag_equiv_reflectivity, diag_eff_radius_qc, &
        diag_ze_rain, diag_ze_ice

   ! locals
   integer :: k, dumi, dumii, dumjj, dumzz
   real(rtype) :: tmp1, tmp2, dum1, dum4, dum5, dum6, rhop
   real(rtype)    :: table_val_qi_fallspd   ! mass-weighted fallspeed              See lines  731 -  808  ums
   real(rtype)    :: table_val_ice_eff_radius   ! effective radius                     See lines 1281 - 1356  eff
   real(rtype)    :: table_val_ni_lammax   ! minimum ice number (lambda limiter)  See lines  704 -  705  nlarge
   real(rtype)    :: table_val_ni_lammin   ! maximum ice number (lambda limiter)  See lines  704 -  705  nsmall
   real(rtype)    :: table_val_ice_reflectivity   ! reflectivity                         See lines  731 -  808  refl
   real(rtype)    :: table_val_ice_mean_diam   ! mass-weighted mean diameter          See lines 1212 - 1279  dmm
   real(rtype)    :: table_val_ice_bulk_dens   ! mass-weighted mean particle density  See lines 1212 - 1279  rhomm

   real(rtype)    :: qc_incld     !in-cloud qi
   real(rtype)    :: nc_incld     !in-cloud ni
   real(rtype)    :: qr_incld     !in-cloud qi
   real(rtype)    :: nr_incld     !in-cloud ni
   real(rtype)    :: qi_incld     !in-cloud qi
   real(rtype)    :: ni_incld     !in-cloud ni
   real(rtype)    :: qm_incld     !in-cloud qm
   real(rtype)    :: bm_incld     !in-cloud bm

   k_loop_final_diagnostics:  do k = kbot,ktop,kdir

      ! cloud:
      if (qc(k).ge.qsmall) then
         qc_incld = qc(k)/cld_frac_l(k)
         nc_incld = nc(k)/cld_frac_l(k)
         call get_cloud_dsd2(qc_incld,nc_incld,mu_c(k),rho(k),nu(k),dnu,lamc(k),  &
              tmp1,tmp2)
       if (mincdnc.gt.0._rtype) nc_incld = max(nc_incld,mincdnc/rho(k)) 
         diag_eff_radius_qc(k) = 0.5_rtype*(mu_c(k)+3._rtype)/lamc(k)
         nc(k) = nc_incld*cld_frac_l(k) !limiters in dsd2 may change nc_incld. Enforcing consistency here.
      else
         diag_eff_radius_qc(k) = 0.0_rtype
         qv(k) = qv(k)+qc(k)
         th_atm(k) = th_atm(k)-exner(k)*qc(k)*latent_heat_vapor(k)*inv_cp
         vap_liq_exchange(k) = vap_liq_exchange(k) - qc(k)
         qc(k) = 0._rtype
         nc(k) = 0._rtype
      endif

      ! rain:
      if (qr(k).ge.qsmall) then
         qr_incld = qr(k)/cld_frac_r(k)
         nr_incld = nr(k)/cld_frac_r(k)
         call get_rain_dsd2(qr_incld,p3_max_mean_rain_size,nr_incld,mu_r(k),lamr(k),tmp1,tmp2)
         nr(k) = nr_incld*cld_frac_r(k) !limiters might change nc_incld... enforcing consistency

         !Note that integrating over the drop-size PDF as done here should only be done to in-cloud
         !quantities but radar reflectivity is likely meant to be a cell ave. Thus nr in the next line
         !really should be cld_frac_r * nr/cld_frac_r. Not doing that since cld_frac_r cancels out.
         ze_rain(k) = nr(k)*(mu_r(k)+6._rtype)*(mu_r(k)+5._rtype)*(mu_r(k)+4._rtype)*           &
              (mu_r(k)+3._rtype)*(mu_r(k)+2._rtype)*(mu_r(k)+1._rtype)/bfb_pow(lamr(k), 6._rtype)
         ze_rain(k) = max(ze_rain(k),1.e-22_rtype)
         diag_ze_rain(k) = 10._rtype*bfb_log10(ze_rain(k)*1.e18_rtype)
      else
         qv(k) = qv(k)+qr(k)
         th_atm(k) = th_atm(k)-exner(k)*qr(k)*latent_heat_vapor(k)*inv_cp
         vap_liq_exchange(k) = vap_liq_exchange(k) - qr(k)
         qr(k) = 0._rtype
         nr(k) = 0._rtype
      endif

      ! ice:

      qi_not_small:  if (qi(k).ge.qsmall) then

         !impose lower limits to prevent taking log of # < 0
         ni(k) = max(ni(k),nsmall)

         qi_incld=qi(k)/cld_frac_i(k)
         ni_incld=ni(k)/cld_frac_i(k)
         qm_incld=qm(k)/cld_frac_i(k)
         bm_incld=bm(k)/cld_frac_i(k)

         call calc_bulkRhoRime(qi_incld,qm_incld,bm_incld,rhop)
         qm(k)=qm_incld*cld_frac_i(k)
         bm(k)=bm_incld*cld_frac_i(k)

         call impose_max_total_ni(ni_incld,max_total_Ni,inv_rho(k))

         call find_lookupTable_indices_1a(dumi,dumjj,dumii,dumzz,dum1,dum4,          &
              dum5,dum6,isize,rimsize,densize,     &
              qi_incld,ni_incld,           &
              qm_incld,rhop)
         !qm(k),zitot(k),rhop)

         call access_lookup_table(dumjj,dumii,dumi, 2,dum1,dum4,dum5,table_val_qi_fallspd)
         call access_lookup_table(dumjj,dumii,dumi, 6,dum1,dum4,dum5,table_val_ice_eff_radius)
         call access_lookup_table(dumjj,dumii,dumi, 7,dum1,dum4,dum5,table_val_ni_lammax)
         call access_lookup_table(dumjj,dumii,dumi, 8,dum1,dum4,dum5,table_val_ni_lammin)
         call access_lookup_table(dumjj,dumii,dumi, 9,dum1,dum4,dum5,table_val_ice_reflectivity)
         call access_lookup_table(dumjj,dumii,dumi,11,dum1,dum4,dum5,table_val_ice_mean_diam)
         call access_lookup_table(dumjj,dumii,dumi,12,dum1,dum4,dum5,table_val_ice_bulk_dens)

         ! impose mean ice size bounds (i.e. apply lambda limiters)
         ! note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
         ni_incld = min(ni_incld,table_val_ni_lammax*ni_incld)
         ni_incld = max(ni_incld,table_val_ni_lammin*ni_incld)
         ni(k) = ni_incld*cld_frac_i(k)

         !--this should already be done in s/r 'calc_bulkRhoRime'
         if (qm(k).lt.qsmall) then
            qm(k) = 0._rtype
            bm(k) = 0._rtype
         endif
         !==

         ! note that reflectivity from lookup table is normalized, so we need to multiply by N
         diag_vm_qi(k)   = table_val_qi_fallspd*rhofaci(k)  
         diag_eff_radius_qi(k)  = table_val_ice_eff_radius          ! units are in m
         diag_diam_qi(k)        = table_val_ice_mean_diam
         rho_qi(k)  = table_val_ice_bulk_dens
         ! note factor of air density below is to convert from m^6/kg to m^6/m^3
         ze_ice(k) = ze_ice(k) + 0.1892_rtype*table_val_ice_reflectivity*ni_incld*rho(k)   ! sum contribution from each ice category (note: 0.1892 = 0.176/0.93)
         ze_ice(k) = max(ze_ice(k),1.e-22_rtype)

         !above formula for ze only makes sense for in-cloud vals, but users expect cell-ave output.
         ze_ice(k) = ze_ice(k)*cld_frac_i(k)
         diag_ze_ice(k) = 10._rtype*bfb_log10(ze_ice(k)*1.e18_rtype)
      else

         qv(k) = qv(k) + qi(k)
         th_atm(k) = th_atm(k) - exner(k)*qi(k)*latent_heat_sublim(k)*inv_cp
         qi(k) = 0._rtype
         ni(k) = 0._rtype
         qm(k) = 0._rtype
         bm(k) = 0._rtype
         diag_diam_qi(k) = 0.0_rtype
         diag_eff_radius_qi(k) = 0.0_rtype
         diag_vm_qi(k) = 0.0_rtype

      endif qi_not_small

      ! sum ze components and convert to dBZ
      diag_equiv_reflectivity(k) = 10._rtype*bfb_log10((ze_rain(k) + ze_ice(k))*1.e18_rtype)

      ! if qr is very small then set Nr to 0 (needs to be done here after call
      ! to ice lookup table because a minimum Nr of nsmall will be set otherwise even if qr=0)
      if (qr(k).lt.qsmall) then
         nr(k) = 0._rtype
      endif

   enddo k_loop_final_diagnostics

 end subroutine p3_main_part3

  !==========================================================================================!

  SUBROUTINE p3_main(qc,nc,qr,nr,th_atm,qv,dt,qi,qm,ni,bm,                                                                                                               &
       pres,dz,nc_nuceat_tend,nccn_prescribed,ni_activated,frzimm,frzcnt,frzdep,inv_qc_relvar,it,precip_liq_surf,precip_ice_surf,its,ite,kts,kte,diag_eff_radius_qc,     &
       diag_eff_radius_qi,rho_qi,do_predict_nc, do_prescribed_CCN,p3_autocon_coeff,p3_accret_coeff,p3_qc_autocon_expon,p3_nc_autocon_expon,p3_qc_accret_expon,           &
       p3_wbf_coeff,p3_mincdnc,p3_max_mean_rain_size,p3_embryonic_rain_size,                                                                                             &
       dpres,exner,qv2qi_depos_tend,precip_total_tend,nevapr,qr_evap_tend,precip_liq_flux,precip_ice_flux,rflx,sflx,cflx,cld_frac_r,cld_frac_l,cld_frac_i,               &
       p3_tend_out,mu_c,lamc,liq_ice_exchange,vap_liq_exchange,                                                                                                          &
       vap_ice_exchange,qv_prev,t_prev,col_location,do_precip_off,nccnst,diag_equiv_reflectivity,diag_ze_rain,diag_ze_ice                                                                     &
#ifdef SCREAM_CONFIG_IS_CMAKE
       ,elapsed_s &
#endif
      )

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

    implicit none

    !----- Input/ouput arguments:  ----------------------------------------------------------!

    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qc         ! cloud, mass mixing ratio         kg kg-1
    ! note: Nc may be specified or predicted (set by do_predict_nc)
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: nc         ! cloud, number mixing ratio       #  kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qr         ! rain, mass mixing ratio          kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: nr         ! rain, number mixing ratio        #  kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qi      ! ice, total mass mixing ratio     kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qm      ! ice, rime mass mixing ratio      kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: ni      ! ice, total number mixing ratio   #  kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: bm      ! ice, rime volume mixing ratio    m3 kg-1

    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qv         ! water vapor mixing ratio         kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: th_atm         ! potential temperature            K
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: pres       ! pressure                         Pa
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: dz        ! vertical grid spacing            m
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: nc_nuceat_tend      ! IN ccn activated number tendency kg-1 s-1
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: nccn_prescribed
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: ni_activated       ! IN actived ice nuclei concentration  1/kg
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: frzimm,frzcnt,frzdep ! From macrophysics aerop (CNT scheme) [#/cm3]
    real(rtype), intent(in)                                     :: dt         ! model time step                  s

    real(rtype), intent(out),   dimension(its:ite)              :: precip_liq_surf    ! precipitation rate, liquid       m s-1
    real(rtype), intent(out),   dimension(its:ite)              :: precip_ice_surf    ! precipitation rate, solid        m s-1
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_eff_radius_qc  ! effective radius, cloud          m
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_eff_radius_qi  ! effective radius, ice            m
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: rho_qi  ! bulk density of ice              kg m-3
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: mu_c       ! Size distribution shape parameter for radiation
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: lamc       ! Size distribution slope parameter for radiation

    integer, intent(in)                                  :: its,ite    ! array bounds (horizontal)
    integer, intent(in)                                  :: kts,kte    ! array bounds (vertical)
    integer, intent(in)                                  :: it         ! time step counter NOTE: starts at 1 for first time step

    logical(btype), intent(in)                           :: do_predict_nc ! .T. (.F.) for prediction (specification) of Nc

    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: dpres       ! pressure thickness               Pa
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: exner      ! Exner expression

    ! OUTPUT for PBUF variables used by other parameterizations
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: qv2qi_depos_tend    ! qitend due to deposition/sublimation
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: precip_total_tend      ! Total precipitation (rain + snow)
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: nevapr     ! evaporation of total precipitation (rain + snow)
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: qr_evap_tend  ! evaporation of rain
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: precip_liq_flux       ! grid-box average rain flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: precip_ice_flux       ! grid-box average ice/snow flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: rflx       ! grid-box average rain flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: sflx       ! grid-box average ice/snow flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: cflx       ! grid-box average cloud droplets flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: liq_ice_exchange ! sum of liq-ice phase change tendenices
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: vap_liq_exchange ! sum of vap-liq phase change tendenices
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: vap_ice_exchange ! sum of vap-ice phase change tendenices

    ! INPUT for prescribed CCN option
    logical(btype), intent(in)                                  :: do_prescribed_CCN

    ! INPUT for idealization options
    logical(btype), intent(in)                                  :: do_precip_off
    real(rtype), intent(in)                                     :: nccnst

    ! INPUT for p3 tuning parameters
    real(rtype), intent(in)                                     :: p3_autocon_coeff         ! autconversion coefficient
    real(rtype), intent(in)                                     :: p3_accret_coeff          ! accretion coefficient
    real(rtype), intent(in)                                     :: p3_qc_autocon_expon      ! autconversion qc exponent
    real(rtype), intent(in)                                     :: p3_nc_autocon_expon      ! autconversion nc exponent
    real(rtype), intent(in)                                     :: p3_qc_accret_expon       ! accretion qc and qr exponent
    real(rtype), intent(in)                                     :: p3_wbf_coeff             ! WBF coefficient
    real(rtype), intent(in)                                     :: p3_mincdnc               ! Lower bound of Nc
    real(rtype), intent(in)                                     :: p3_max_mean_rain_size    ! max mean rain size allowed
    real(rtype), intent(in)                                     :: p3_embryonic_rain_size   ! embryonic rain size from autoconversion

    ! INPUT needed for PBUF variables used by other parameterizations

    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: cld_frac_i, cld_frac_l, cld_frac_r ! Ice, Liquid and Rain cloud fraction
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: qv_prev, t_prev                    ! qv and t from previous p3_main call
    ! AaronDonahue, the following variable (p3_tend_out) is a catch-all for passing P3-specific variables outside of p3_main
    ! so that they can be written as ouput.  NOTE TO C++ PORT: This variable is entirely optional and doesn't need to be
    ! included in the port to C++, or can be changed if desired.
    real(rtype), intent(out),   dimension(its:ite,kts:kte,49)   :: p3_tend_out ! micro physics tendencies
    real(rtype), intent(in),    dimension(its:ite,3)            :: col_location
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: inv_qc_relvar
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_equiv_reflectivity,diag_ze_rain,diag_ze_ice  ! equivalent reflectivity [dBZ]

#ifdef SCREAM_CONFIG_IS_CMAKE
    real(rtype), optional, intent(out) :: elapsed_s ! duration of main loop in seconds
#endif

    !
    !----- Local variables and parameters:  -------------------------------------------------!
    !

    ! These outputs are no longer provided by p3_main.
    real(rtype), dimension(its:ite,kts:kte) :: diag_vm_qi ! mass-weighted fall speed of ice  m s-1
    real(rtype), dimension(its:ite,kts:kte) :: diag_diam_qi  ! mean diameter of ice             m
    real(rtype), dimension(its:ite,kts:kte) :: pratot   ! accretion of cloud by rain
    real(rtype), dimension(its:ite,kts:kte) :: prctot   ! autoconversion of cloud to rain

    real(rtype), dimension(its:ite,kts:kte) :: mu_r  ! shape parameter of rain
    real(rtype), dimension(its:ite,kts:kte) :: t_atm     ! temperature at the beginning of the microhpysics step [K]

    ! 2D size distribution and fallspeed parameters:

    real(rtype), dimension(its:ite,kts:kte) :: lamr
    real(rtype), dimension(its:ite,kts:kte) :: logn0r

    real(rtype), dimension(its:ite,kts:kte) :: nu
    real(rtype), dimension(its:ite,kts:kte) :: cdist
    real(rtype), dimension(its:ite,kts:kte) :: cdist1
    real(rtype), dimension(its:ite,kts:kte) :: cdistr

    ! Variables needed for in-cloud calculations
    real(rtype), dimension(its:ite,kts:kte) :: inv_cld_frac_i, inv_cld_frac_l, inv_cld_frac_r ! Inverse cloud fractions (1/cld)
    real(rtype), dimension(its:ite,kts:kte) :: qc_incld, qr_incld, qi_incld, qm_incld ! In cloud mass-mixing ratios
    real(rtype), dimension(its:ite,kts:kte) :: nc_incld, nr_incld, ni_incld, bm_incld ! In cloud number concentrations

    real(rtype), dimension(its:ite,kts:kte)      :: inv_dz,inv_rho,ze_ice,ze_rain,prec,rho,       &
         rhofacr,rhofaci,acn,latent_heat_sublim,latent_heat_vapor,latent_heat_fusion,qv_sat_l,qv_sat_i,qv_supersat_i,       &
         tmparr1,inv_exner

    ! -- scalar locals -- !

    real(rtype) :: inv_dt, timeScaleFactor

    integer :: ktop,kbot,kdir,i

    logical(btype) :: is_nucleat_possible, is_hydromet_present

    !--These will be added as namelist parameters in the future
    logical(btype), parameter :: debug_ON     = .true.  !.true. to switch on debugging checks/traps throughout code  TODO: Turn this back off as default once the tlay error is found.
    logical(btype), parameter :: debug_ABORT  = .false.  !.true. will result in forced abort in s/r 'check_values'

    real(rtype),dimension(its:ite,kts:kte) :: qc_old, nc_old, qr_old, nr_old, qi_old, ni_old, qv_old, th_atm_old
    integer :: knc
    real(rtype) :: mincdnc

#ifdef SCREAM_CONFIG_IS_CMAKE
    integer :: clock_count1, clock_count_rate, clock_count_max, clock_count2, clock_count_diff
#endif

    !-----------------------------------------------------------------------------------!
    !  End of variables/parameters declarations
    !-----------------------------------------------------------------------------------!

    ! direction of vertical leveling:
    !PMC got rid of 'model' option so we could just replace ktop with kts everywhere...
    ktop = kts        !k of top level
    kbot = kte        !k of bottom level
    kdir = -1         !(k: 1=top, nk=bottom)

    !PMC deleted 'threshold size difference' calculation for multicategory here

    inv_dz    = 1._rtype/dz  ! inverse of thickness of layers
    inv_dt        = 1._rtype/dt   ! inverse model time step

    ! Compute time scale factor over which to apply soft rain lambda limiter
    ! note: '1./max(30.,dt)' = '1.*min(1./30., 1./dt)'
    timeScaleFactor = min(1._rtype/120._rtype, inv_dt)

    precip_liq_surf   = 0._rtype
    precip_ice_surf   = 0._rtype
    pratot    = 0._rtype
    prctot    = 0._rtype
    prec      = 0._rtype
    mu_r      = 0._rtype
    diag_equiv_reflectivity   = -99._rtype
    diag_ze_rain = -99.0_rtype 
    diag_ze_ice  = -99.0_rtype

    ze_ice    = 1.e-22_rtype
    ze_rain   = 1.e-22_rtype
    diag_eff_radius_qc = 10.e-6_rtype ! default value
    diag_eff_radius_qi = 25.e-6_rtype ! default value
    diag_vm_qi  = 0._rtype
    diag_diam_qi   = 0._rtype
    rho_qi = 0._rtype

    qv2qi_depos_tend = 0._rtype
    precip_total_tend   = 0._rtype
    nevapr  = 0._rtype
    precip_liq_flux    = 0._rtype
    precip_ice_flux    = 0._rtype
    rflx    = 0._rtype
    sflx    = 0._rtype
    cflx    = 0._rtype
    p3_tend_out = 0._rtype

    inv_cld_frac_i = 1.0_rtype/cld_frac_i
    inv_cld_frac_l = 1.0_rtype/cld_frac_l
    inv_cld_frac_r = 1.0_rtype/cld_frac_r

    qr_evap_tend = 0._rtype
    liq_ice_exchange = 0._rtype
    vap_liq_exchange = 0._rtype
    vap_ice_exchange = 0._rtype

    mu_c = 0.0_rtype
    lamc = 0.0_rtype
    ! AaronDonahue added exner term to replace all instances of th_atm(i,k)/t(i,k), since th_atm(i,k) is updated but t_atm(i,k) is not, and this was
    ! causing energy conservation errors.
    inv_exner = 1._rtype/exner        !inverse of Exner expression, used when converting potential temp to temp
    t_atm       = th_atm    *inv_exner    !compute temperature from theta (value at beginning of microphysics step)
    qv      = max(qv,0._rtype)        !clip water vapor to prevent negative values passed in (beginning of microphysics)
    ! AaronDonahue added this load of latent heat to be consistent with E3SM, since the inconsistentcy was causing water conservation errors.
    call get_latent_heat(its,ite,kts,kte,latent_heat_vapor,latent_heat_sublim,latent_heat_fusion)

   ! initialize microphysics processes tendency output
    qc_old = qc         ! Liq. microphysics tendency, initialize
    nc_old = nc         ! Liq. # microphysics tendency, initialize
    qr_old = qr         ! Rain microphysics tendency, initialize
    nr_old = nr         ! Rain # microphysics tendency, initialize
    qi_old = qi   ! Ice  microphysics tendency, initialize
    ni_old = ni   ! Ice  # microphysics tendency, initialize
    qv_old = qv         ! Vapor  microphysics tendency, initialize
    th_atm_old = th_atm         ! Pot. Temp. microphysics tendency, initialize
    mincdnc = p3_mincdnc
#ifdef SCREAM_CONFIG_IS_CMAKE
    call system_clock(clock_count1, clock_count_rate, clock_count_max)
#endif

    !==
    !-----------------------------------------------------------------------------------!
    i_loop_main: do i = its,ite  ! main i-loop (around the entire scheme)

      ! ... update column in the debug_info module
      ! call get_debug_column_id(i)
      ! if (debug_ON) call check_values(qv,T,i,it,debug_ABORT,100,col_location)

      if (debug_ON) then
         tmparr1(i,:) = th_atm(i,:)*inv_exner(i,:)!(pres(i,:)*1.e-5)**(rd*inv_cp)
         call check_values(qv(i,:),tmparr1(i,:),kts,kte,it,debug_ABORT,100,col_location(i,:))
      endif


       call p3_main_part1(kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt, &
            pres(i,:), dpres(i,:), dz(i,:), nc_nuceat_tend(i,:), nccn_prescribed(i,:), exner(i,:), inv_exner(i,:), &
            inv_cld_frac_l(i,:), inv_cld_frac_i(i,:), inv_cld_frac_r(i,:), latent_heat_vapor(i,:), latent_heat_sublim(i,:), latent_heat_fusion(i,:), &
            t_atm(i,:), rho(i,:), inv_rho(i,:), qv_sat_l(i,:), qv_sat_i(i,:), qv_supersat_i(i,:), rhofacr(i,:), &
            rhofaci(i,:), acn(i,:), qv(i,:), th_atm(i,:), qc(i,:), nc(i,:), qr(i,:), nr(i,:), &
            qi(i,:), ni(i,:), qm(i,:), bm(i,:), qc_incld(i,:), qr_incld(i,:), &
            qi_incld(i,:), qm_incld(i,:), nc_incld(i,:), nr_incld(i,:), &
            ni_incld(i,:), bm_incld(i,:), is_nucleat_possible, is_hydromet_present, nccnst)

      if (debug_ON) then
         tmparr1(i,:) = th_atm(i,:)*inv_exner(i,:)!(pres(i,:)*1.e-5)**(rd*inv_cp)
         call check_values(qv(i,:),tmparr1(i,:),kts,kte,it,debug_ABORT,200,col_location(i,:))
      endif

       !jump to end of i-loop if is_nucleat_possible=.false.  (i.e. skip everything)
       if (.not. (is_nucleat_possible .or. is_hydromet_present)) goto 333

       call p3_main_part2(kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt, inv_dt, &
            p3_autocon_coeff,p3_accret_coeff,p3_qc_autocon_expon,p3_nc_autocon_expon,p3_qc_accret_expon, &
            p3_wbf_coeff,p3_embryonic_rain_size, p3_max_mean_rain_size, &
            pres(i,:), dpres(i,:), dz(i,:), nc_nuceat_tend(i,:), exner(i,:), inv_exner(i,:), &
            inv_cld_frac_l(i,:), inv_cld_frac_i(i,:), inv_cld_frac_r(i,:), ni_activated(i,:), inv_qc_relvar(i,:), &
            cld_frac_i(i,:), cld_frac_l(i,:), cld_frac_r(i,:), qv_prev(i,:), t_prev(i,:), &
            t_atm(i,:), rho(i,:), inv_rho(i,:), qv_sat_l(i,:), &
            qv_sat_i(i,:), qv_supersat_i(i,:), rhofacr(i,:), rhofaci(i,:), acn(i,:), qv(i,:), th_atm(i,:), &
            qc(i,:), nc(i,:), qr(i,:), nr(i,:), qi(i,:), ni(i,:), qm(i,:), &
            bm(i,:), latent_heat_vapor(i,:), latent_heat_sublim(i,:), latent_heat_fusion(i,:), qc_incld(i,:), qr_incld(i,:), &
            qi_incld(i,:), qm_incld(i,:), nc_incld(i,:), nr_incld(i,:), ni_incld(i,:), &
            bm_incld(i,:), mu_c(i,:), nu(i,:), lamc(i,:), cdist(i,:), cdist1(i,:), &
            cdistr(i,:), mu_r(i,:), lamr(i,:), logn0r(i,:), qv2qi_depos_tend(i,:), precip_total_tend(i,:), &
            nevapr(i,:), qr_evap_tend(i,:), vap_liq_exchange(i,:), vap_ice_exchange(i,:), &
            liq_ice_exchange(i,:), pratot(i,:), prctot(i,:), frzimm(i,:), frzcnt(i,:), frzdep(i,:), p3_tend_out(i,:,:), is_hydromet_present, &
	    do_precip_off, nccnst)


       ! measure microphysics processes tendency output
       p3_tend_out(i,:,42) = ( qc(i,:)      - qc_old(i,:) ) * inv_dt       ! Liq. microphysics tendency, measure
       p3_tend_out(i,:,43) = ( nc(i,:)      - nc_old(i,:) ) * inv_dt       ! Liq. # microphysics tendency, measure
       p3_tend_out(i,:,44) = ( qr(i,:)      - qr_old(i,:) ) * inv_dt       ! Rain microphysics tendency, measure
       p3_tend_out(i,:,45) = ( nr(i,:)      - nr_old(i,:) ) * inv_dt       ! Rain # microphysics tendency, measure
       p3_tend_out(i,:,46) = ( qi(i,:)      - qi_old(i,:) ) * inv_dt       ! Ice  microphysics tendency, measure
       p3_tend_out(i,:,47) = ( ni(i,:)      - ni_old(i,:) ) * inv_dt       ! Ice  # microphysics tendency, measure
       p3_tend_out(i,:,48) = ( qv(i,:)      - qv_old(i,:) ) * inv_dt       ! Vapor  microphysics tendency, measure
       p3_tend_out(i,:,49) = ( th_atm(i,:)  - th_atm_old(i,:) ) * inv_dt   ! Pot. Temp. microphysics tendency, measure
       !NOTE: At this point, it is possible to have negative (but small) nc, nr, ni.  This is not
       !      a problem; those values get clipped to zero in the sedimentation section (if necessary).
       !      (This is not done above simply for efficiency purposes.)

       !      if (debug_ON) then
       !         tmparr1(i,:) = th(i,:)*inv_exner(i,:)!(pres(i,:)*1.e-5)**(rd*inv_cp)
       !         call check_values(qv,tmparr1,i,it,debug_ABORT,300,col_location)
       !      endif

       if (debug_ON) then
         tmparr1(i,:) = th_atm(i,:)*inv_exner(i,:)!(pres(i,:)*1.e-5)**(rd*inv_cp)
         call check_values(qv(i,:),tmparr1(i,:),kts,kte,it,debug_ABORT,300,col_location(i,:))
       endif

       if (.not. is_hydromet_present) goto 333

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
         qc_incld(i,:),rho(i,:),inv_rho(i,:),cld_frac_l(i,:),acn(i,:),inv_dz(i,:), &
         dt,inv_dt,dnu,do_predict_nc, &
         qc(i,:),nc(i,:),nc_incld(i,:),mu_c(i,:),lamc(i,:),precip_liq_surf(i),cflx(i,:),p3_tend_out(i,:,36),p3_tend_out(i,:,37))

       !------------------------------------------------------------------------------------------!
       ! Rain sedimentation:  (adaptive substepping)
       p3_tend_out(i,:,38) = qr(i,:) ! Rain sedimentation tendency, initialize
       p3_tend_out(i,:,39) = nr(i,:) ! Rain # sedimentation tendency, initialize

       call rain_sedimentation(kts,kte,ktop,kbot,kdir, &
         qr_incld(i,:),rho(i,:),inv_rho(i,:),rhofacr(i,:),cld_frac_r(i,:),inv_dz(i,:),dt,inv_dt,p3_max_mean_rain_size, &
         qr(i,:),nr(i,:),nr_incld(i,:),mu_r(i,:),lamr(i,:),precip_liq_surf(i),precip_liq_flux(i,:),rflx(i,:),p3_tend_out(i,:,38), &
         p3_tend_out(i,:,39))

       !------------------------------------------------------------------------------------------!
       ! Ice sedimentation:  (adaptive substepping)
       p3_tend_out(i,:,40) = qi(i,:) ! Ice sedimentation tendency, initialize
       p3_tend_out(i,:,41) = ni(i,:) ! Ice # sedimentation tendency, initialize

       call ice_sedimentation(kts,kte,ktop,kbot,kdir,    &
         rho(i,:),inv_rho(i,:),rhofaci(i,:),cld_frac_i(i,:),inv_dz(i,:),dt,inv_dt, &
         qi(i,:),qi_incld(i,:),ni(i,:),qm(i,:),qm_incld(i,:),bm(i,:),bm_incld(i,:),ni_incld(i,:), &
         precip_ice_surf(i),precip_ice_flux(i,:),sflx(i,:),p3_tend_out(i,:,40),p3_tend_out(i,:,41))


       !*** Add tendencies for the next processes ***  

       !p3_tend_out(i,:,36) = qc(i,:)
       !p3_tend_out(i,:,37) = nc(i,:)
       !p3_tend_out(i,:,38) = qr(i,:)
       !p3_tend_out(i,:,39) = nr(i,:)
       !p3_tend_out(i,:,40) = qi(i,:)
       !p3_tend_out(i,:,41) = ni(i,:)

       !.......................................
       ! homogeneous freezing of cloud and rain
       call homogeneous_freezing(kts,kte,ktop,kbot,kdir,t_atm(i,:),exner(i,:),latent_heat_fusion(i,:),  &
         qc(i,:),nc(i,:),qr(i,:),nr(i,:),qi(i,:),ni(i,:),qm(i,:),bm(i,:),th_atm(i,:))

       !.........................................................
       ! Instantenous melting of ice/snow at T = t_snow_melt = 2c    
       call ice_complete_melting(kts,kte,ktop,kbot,kdir,qi(i,:),ni(i,:),qm(i,:),latent_heat_fusion(i,:),exner(i,:),th_atm(i,:), & 
            qr(i,:),nr(i,:),qc(i,:),nc(i,:))

         do knc=kbot,ktop,kdir
          if ((mincdnc.gt.0._rtype).and.(qc(i,knc).ge.qsmall)) then
           nc(i,knc) = max(nc(i,knc),mincdnc*cld_frac_l(i,knc)/rho(i,knc))
           nc_incld(i,knc) = max(nc_incld(i,knc),mincdnc/rho(i,knc))
          end if
         end do 

       !...................................................
       ! final checks to ensure consistency of mass/number
       ! and compute diagnostic fields for output
       call p3_main_part3(kts, kte, kbot, ktop, kdir, p3_max_mean_rain_size,                                                  &
            exner(i,:), cld_frac_l(i,:), cld_frac_r(i,:), cld_frac_i(i,:),                                                    &
            rho(i,:), inv_rho(i,:), rhofaci(i,:), qv(i,:), th_atm(i,:), qc(i,:), nc(i,:), qr(i,:), nr(i,:), qi(i,:), ni(i,:), &
            qm(i,:), bm(i,:), latent_heat_vapor(i,:), latent_heat_sublim(i,:),                                                &
            mu_c(i,:), nu(i,:), lamc(i,:), mu_r(i,:), lamr(i,:), vap_liq_exchange(i,:),                                       &
            ze_rain(i,:), ze_ice(i,:), diag_vm_qi(i,:), diag_eff_radius_qi(i,:), diag_diam_qi(i,:), rho_qi(i,:), diag_equiv_reflectivity(i,:), diag_eff_radius_qc(i,:), &
            diag_ze_rain(i,:),diag_ze_ice(i,:),mincdnc)

       !   if (debug_ON) call check_values(qv,Ti,it,debug_ABORT,800,col_location)

       !..............................................
       ! merge ice categories with similar properties

       !   note:  this should be relocated to above, such that the diagnostic
       !          ice properties are computed after merging

       !PMC nCat deleted nCat>1 stuff

       !.....................................................

333    continue

       if (debug_ON) then
          tmparr1(i,:) = th_atm(i,:)*inv_exner(i,:)!(pres(i,:)*1.e-5)**(rd*inv_cp)
          call check_values(qv(i,:),tmparr1(i,:),kts,kte,it,debug_ABORT,900,col_location(i,:))
       endif

       !.....................................................

    enddo i_loop_main

#ifdef SCREAM_CONFIG_IS_CMAKE
    call system_clock(clock_count2, clock_count_rate, clock_count_max)
    clock_count_diff = clock_count2 - clock_count1
    if (present(elapsed_s)) then
      elapsed_s = real(clock_count_diff) / real(clock_count_rate)
    endif
#endif

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
    iproc1 = ice_table_vals(dumjj,dumii,dumi,index)+(dum1-real(dumi))*(ice_table_vals(dumjj,dumii,       &
         dumi+1,index)-ice_table_vals(dumjj,dumii,dumi,index))

    ! linearly interpolate to get process rates for rimed fraction index + 1

    gproc1 = ice_table_vals(dumjj,dumii+1,dumi,index)+(dum1-real(dumi))*(ice_table_vals(dumjj,dumii+1,   &
         dumi+1,index)-ice_table_vals(dumjj,dumii+1,dumi,index))

    tmp1   = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

    ! get value at density index + 1

    ! first interpolate for current rimed fraction index

    iproc1 = ice_table_vals(dumjj+1,dumii,dumi,index)+(dum1-real(dumi))*(ice_table_vals(dumjj+1,dumii,   &
         dumi+1,index)-ice_table_vals(dumjj+1,dumii,dumi,index))

    ! linearly interpolate to get process rates for rimed fraction index + 1

    gproc1 = ice_table_vals(dumjj+1,dumii+1,dumi,index)+(dum1-real(dumi))*(ice_table_vals(dumjj+1,       &
         dumii+1,dumi+1,index)-ice_table_vals(dumjj+1,dumii+1,dumi,index))

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
    dproc1  = collect_table_vals(dumjj,dumii,dumi,dumj,index)+(dum1-real(dumi))*                &
         (collect_table_vals(dumjj,dumii,dumi+1,dumj,index)-collect_table_vals(dumjj,dumii,dumi,    &
         dumj,index))

    dproc2  = collect_table_vals(dumjj,dumii,dumi,dumj+1,index)+(dum1-real(dumi))*             &
         (collect_table_vals(dumjj,dumii,dumi+1,dumj+1,index)-collect_table_vals(dumjj,dumii,dumi,  &
         dumj+1,index))

    iproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)

    ! rime fraction index + 1

    dproc1  = collect_table_vals(dumjj,dumii+1,dumi,dumj,index)+(dum1-real(dumi))*             &
         (collect_table_vals(dumjj,dumii+1,dumi+1,dumj,index)-collect_table_vals(dumjj,dumii+1,     &
         dumi,dumj,index))

    dproc2  = collect_table_vals(dumjj,dumii+1,dumi,dumj+1,index)+(dum1-real(dumi))*           &
         (collect_table_vals(dumjj,dumii+1,dumi+1,dumj+1,index)-collect_table_vals(dumjj,dumii+1,   &
         dumi,dumj+1,index))

    gproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)
    tmp1    = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

    ! density index + 1

    ! current rime fraction index

    dproc1  = collect_table_vals(dumjj+1,dumii,dumi,dumj,index)+(dum1-real(dumi))*             &
         (collect_table_vals(dumjj+1,dumii,dumi+1,dumj,index)-collect_table_vals(dumjj+1,dumii,     &
         dumi,dumj,index))

    dproc2  = collect_table_vals(dumjj+1,dumii,dumi,dumj+1,index)+(dum1-real(dumi))*           &
         (collect_table_vals(dumjj+1,dumii,dumi+1,dumj+1,index)-collect_table_vals(dumjj+1,dumii,   &
         dumi,dumj+1,index))

    iproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)

    ! rime fraction index + 1

    dproc1  = collect_table_vals(dumjj+1,dumii+1,dumi,dumj,index)+(dum1-real(dumi))*           &
         (collect_table_vals(dumjj+1,dumii+1,dumi+1,dumj,index)-collect_table_vals(dumjj+1,dumii+1, &
         dumi,dumj,index))

    dproc2  = collect_table_vals(dumjj+1,dumii+1,dumi,dumj+1,index)+(dum1-real(dumi))*         &
         (collect_table_vals(dumjj+1,dumii+1,dumi+1,dumj+1,index)-collect_table_vals(dumjj+1,       &
         dumii+1,dumi,dumj+1,index))

    gproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)
    tmp2    = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

    ! interpolate over density to get final values
    proc    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

    return
  END SUBROUTINE access_lookup_table_coll


  !======================================================================================!

  subroutine find_lookupTable_indices_1a(dumi,dumjj,dumii,dumzz,dum1,dum4,dum5,dum6,      &
       isize,rimsize,densize,qi,ni,qm,   &
       rhop)

    !------------------------------------------------------------------------------------------!
    ! Finds indices in 3D ice (only) lookup table
    !------------------------------------------------------------------------------------------!

    implicit none

    ! arguments:
    integer, intent(out) :: dumi,dumjj,dumii,dumzz
    real(rtype),    intent(out) :: dum1,dum4,dum5,dum6
    integer, intent(in)  :: isize,rimsize,densize
    real(rtype),    intent(in)  :: qi,ni,qm,rhop

    !------------------------------------------------------------------------------------------!
    ! find index for qi (normalized ice mass mixing ratio = qi/ni)
    !             dum1 = (log10(qi)+16.)/0.70757  !orig
    !             dum1 = (log10(qi)+16.)*1.41328
    ! we are inverting this equation from the lookup table to solve for i:
    ! qi/ni=261.7**((i+10)*0.1)*1.e-18
    dum1 = (bfb_log10(qi/ni)+18._rtype)*lookup_table_1a_dum1_c-10._rtype ! For computational efficiency
    dumi = int(dum1)
    ! set limits (to make sure the calculated index doesn't exceed range of lookup table)
    dum1 = min(dum1,real(isize,rtype))
    dum1 = max(dum1,1._rtype)
    dumi = max(1,dumi)
    dumi = min(isize-1,dumi)

    ! find index for rime mass fraction
    dum4  = (qm/qi)*3._rtype + 1._rtype
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
       dumlr = bfb_cbrt(qr/(pi*rho_h2o*nr))
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
  subroutine get_cloud_dsd2(qc,nc,mu_c,rho,nu,dnu,lamc,cdist,cdist1)

    implicit none

    !arguments:
    real(rtype), dimension(:), intent(in)  :: dnu
    real(rtype),     intent(in)            :: qc,rho
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
          nc   = 6._rtype*(lamc*lamc*lamc)*qc/(pi*rho_h2o*(mu_c+3._rtype)*(mu_c+2._rtype)*(mu_c+1._rtype))
       elseif (lamc.gt.lammax) then
          lamc = lammax
          nc   = 6._rtype*(lamc*lamc*lamc)*qc/(pi*rho_h2o*(mu_c+3._rtype)*(mu_c+2._rtype)*(mu_c+1._rtype))
       endif

       cdist  = nc*(mu_c+1._rtype)/lamc
       cdist1 = nc/bfb_gamma(mu_c+1._rtype)

    else

       lamc   = 0._rtype
       cdist  = 0._rtype
       cdist1 = 0._rtype

    endif

   return

  end subroutine get_cloud_dsd2


  !===========================================================================================
  subroutine get_rain_dsd2(qr,p3_max_mean_rain_size,nr,mu_r,lamr,cdistr,logn0r)

    ! Computes and returns rain size distribution parameters

    implicit none

    !arguments:
    real(rtype),     intent(in)            :: qr
    real(rtype),     intent(in)            :: p3_max_mean_rain_size
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
       lammin = (mu_r+1._rtype)*(1._rtype)/p3_max_mean_rain_size   ! breakup is explicitly included (mean size set by namelist parameter)

       ! apply lambda limiters for rain
       if (lamr.lt.lammin) then
          lamr = lammin
          nr   = bfb_exp(3._rtype*bfb_log(lamr)+bfb_log(qr)+bfb_log(bfb_gamma(mu_r+1._rtype))-bfb_log(bfb_gamma(mu_r+4._rtype)))/(cons1)
       elseif (lamr.gt.lammax) then
          lamr = lammax
          nr   = bfb_exp(3._rtype*bfb_log(lamr)+bfb_log(qr)+bfb_log(bfb_gamma(mu_r+1._rtype))-bfb_log(bfb_gamma(mu_r+4._rtype)))/(cons1)
       endif

       cdistr  = nr/bfb_gamma(mu_r+1._rtype)
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
    !  and adjusts qm and bm appropriately.
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
  subroutine impose_max_total_ni(ni_local,max_total_ni,inv_rho_local)

    !--------------------------------------------------------------------------------
    ! Impose maximum total ice number concentration (total of all ice categories).
    ! If the sum of all ni(:) exceeds maximum allowable, each category to preserve
    ! ratio of number between categories.
    !--------------------------------------------------------------------------------

    implicit none

    !arguments:
    real(rtype), intent(inout)               :: ni_local      !PMC - scalar now that nCat deleted.
    real(rtype), intent(in)                  :: max_total_ni,inv_rho_local

    !local variables:
    real(rtype)                              :: dum

    if (ni_local.ge.1.e-20_rtype) then
       dum = max_total_ni*inv_rho_local/ni_local
       ni_local = ni_local*min(dum,1._rtype)
    endif

  end subroutine impose_max_total_ni


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
               '** WARNING IN P3_MAIN -- src, gcol, lon, lat, lvl, tstep, T:',source_ind,', ', &
               int(col_loc(1)),', ',col_loc(2),', ',col_loc(3),', ',k,', ',timestepcount,', ',T(k)
          trap = .true.
       endif
       if (.not.(Qv(k)>=0. .and. Qv(k)<Q_high)) then
          write(iulog,'(a60,i5,a2,i8,a2,f8.4,a2,f8.4,a2,i4,a2,i8,a2,e16.8)') &
               '** WARNING IN P3_MAIN -- src, gcol, lon, lat, lvl, tstep, Qv:',source_ind,', ', &
               int(col_loc(1)),', ',col_loc(2),', ',col_loc(3),', ',k,', ',timestepcount,', ',Qv(k)
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
               source_ind,' in file:',&
               __FILE__,&
               ' at line:',__LINE__
          call endscreamrun(err_msg)
       endif
    endif

   return

  end subroutine check_values

  subroutine ice_cldliq_collection(rho,t_atm,rhofaci,    &
  table_val_qc2qi_collect,qi_incld,qc_incld,ni_incld,nc_incld,    &
             qccol,nc_collect_tend,qc2qr_ice_shed_tend,ncshdc)

   !.......................
   ! collection of droplets

   ! here we multiply rates by air density, air density fallspeed correction
   ! factor, and collection efficiency since these parameters are not
   ! included in lookup table calculations
   ! for T < 273.15, assume collected cloud water is instantly frozen
   ! note 'f1pr' values are normalized, so we need to multiply by N

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: table_val_qc2qi_collect  ! collection of cloud water by ice
   real(rtype), intent(in) :: qi_incld
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: ni_incld
   real(rtype), intent(in) :: nc_incld


   real(rtype), intent(out) :: qccol
   real(rtype), intent(out) :: nc_collect_tend
   real(rtype), intent(out) :: qc2qr_ice_shed_tend
   real(rtype), intent(out) :: ncshdc

   if (qi_incld .ge.qsmall .and. qc_incld .ge.qsmall) then
      if  (t_atm .le.T_zerodegc) then
         qccol = rhofaci*table_val_qc2qi_collect*qc_incld*eci*rho*ni_incld
         nc_collect_tend = rhofaci*table_val_qc2qi_collect*nc_incld*eci*rho*ni_incld
      else if (t_atm .gt. T_zerodegc) then
         ! for T > 273.15, assume cloud water is collected and shed as rain drops
         ! sink for cloud water mass and number, note qcshed is source for rain mass
         qc2qr_ice_shed_tend = rhofaci*table_val_qc2qi_collect*qc_incld*eci*rho*ni_incld
         nc_collect_tend = rhofaci*table_val_qc2qi_collect*nc_incld*eci*rho*ni_incld
         ! source for rain number, assume 1 mm drops are shed
         ncshdc = qc2qr_ice_shed_tend*inv_dropmass
      end if
   end if

   return

  end subroutine ice_cldliq_collection


  subroutine ice_rain_collection(rho,t_atm,rhofaci,    &
  logn0r,table_val_nr_collect,table_val_qr2qi_collect,qi_incld,ni_incld,qr_incld,    &
  qrcol, nr_collect_tend)

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
   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: logn0r
   real(rtype), intent(in) :: table_val_nr_collect !collection of rain number by ice
   real(rtype), intent(in) :: table_val_qr2qi_collect !collection of rain mass by ice
   real(rtype), intent(in) :: qi_incld
   real(rtype), intent(in) :: ni_incld
   real(rtype), intent(in) :: qr_incld

   real(rtype), intent(out) :: qrcol
   real(rtype), intent(out) :: nr_collect_tend

   if (qi_incld.ge.qsmall .and. qr_incld.ge.qsmall) then
      if (t_atm.le.T_zerodegc) then
         ! note: table_val_qr2qi_collect and logn0r are already calculated as log_10
         qrcol = bfb_pow(10._rtype,(table_val_qr2qi_collect+logn0r))*rho*rhofaci*eri*ni_incld
         nr_collect_tend = bfb_pow(10._rtype,(table_val_nr_collect+logn0r))*rho*rhofaci*eri*ni_incld
      else if (t_atm .gt. T_zerodegc) then
         ! rain number sink due to collection
         ! for T > 273.15, assume collected rain number is shed as
         ! 1 mm drops
         ! note that melting of ice number is scaled to the loss
         ! rate of ice mass due to melting
         ! collection of rain above freezing does not impact total rain mass
         nr_collect_tend  = bfb_pow(10._rtype,(table_val_nr_collect + logn0r))*rho*rhofaci*eri*ni_incld
         ! for now neglect shedding of ice collecting rain above freezing, since snow is
         ! not expected to shed in these conditions (though more hevaily rimed ice would be
         ! expected to lead to shedding)
      endif
   endif

   return

  end subroutine ice_rain_collection

  subroutine ice_self_collection(rho,rhofaci,    &
  table_val_ni_self_collect,eii,qm_incld,qi_incld,ni_incld,    &
             ni_selfcollect_tend)

   ! self-collection of ice

   ! here we multiply rates by collection efficiency, air density,
   ! and air density correction factor since these are not included
   ! in the lookup table calculations
   ! note 'f1pr' values are normalized, so we need to multiply by N

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: table_val_ni_self_collect ! ice collection within a category
   real(rtype), intent(in) :: eii
   real(rtype), intent(in) :: qm_incld
   real(rtype), intent(in) :: qi_incld
   real(rtype), intent(in) :: ni_incld

   real(rtype), intent(out) :: ni_selfcollect_tend

   real(rtype) :: tmp1, Eii_fact

   if (qi_incld.ge.qsmall) then
      ! Determine additional collection efficiency factor to be applied to ice-ice collection.
      ! The computed values of qicol and nicol are multipiled by Eii_fact to gradually shut off collection
      ! if ice is highly rimed.
      if (qm_incld>0._rtype) then
         tmp1 = qm_incld/qi_incld   !rime mass fraction
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

      ni_selfcollect_tend = table_val_ni_self_collect*rho*eii*Eii_fact*rhofaci*ni_incld*ni_incld

   endif

   return

end subroutine ice_self_collection

!PMC note - indentation pattern changes here.

subroutine ice_melting(rho,t_atm,pres,rhofaci,    &
table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion,dv,sc,mu,kap,qv,qi_incld,ni_incld,    &
           qi2qr_melt_tend,ni2nr_melt_tend)
   ! melting
   ! need to add back accelerated melting due to collection of ice mass by rain (pracsw1)
   ! note 'f1pr' values are normalized, so we need to multiply by N
   ! currently enhanced melting from collision is neglected
   ! include RH dependence

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: pres
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: table_val_qi2qr_melting ! melting
   real(rtype), intent(in) :: table_val_qi2qr_vent_melt ! melting (ventilation term)
   real(rtype), intent(in) :: latent_heat_vapor
   real(rtype), intent(in) :: latent_heat_fusion
   real(rtype), intent(in) :: dv
   real(rtype), intent(in) :: sc
   real(rtype), intent(in) :: mu
   real(rtype), intent(in) :: kap
   real(rtype), intent(in) :: qv
   real(rtype), intent(in) :: qi_incld
   real(rtype), intent(in) :: ni_incld

   real(rtype), intent(out) :: qi2qr_melt_tend
   real(rtype), intent(out) :: ni2nr_melt_tend

   real(rtype) :: qsat0

   if (qi_incld .ge.qsmall .and. t_atm.gt.T_zerodegc) then
      qsat0 = qv_sat( T_zerodegc,pres,0 )

      qi2qr_melt_tend = ((table_val_qi2qr_melting+table_val_qi2qr_vent_melt*bfb_cbrt(sc)*bfb_sqrt(rhofaci*rho/mu))*((t_atm-   &
      T_zerodegc)*kap-rho*latent_heat_vapor*dv*(qsat0-qv))*2._rtype*pi/latent_heat_fusion)*ni_incld

      qi2qr_melt_tend = max(qi2qr_melt_tend,0._rtype)
      ni2nr_melt_tend = qi2qr_melt_tend*(ni_incld/qi_incld)

   endif

   return

end subroutine ice_melting


subroutine ice_cldliq_wet_growth(rho,t_atm,pres,rhofaci,    &
table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion,dv,kap,mu,sc,    &
qv,qc_incld,qi_incld,ni_incld,qr_incld,    &
           log_wetgrowth,qrcol,qccol,qwgrth,nr_ice_shed_tend,qc2qr_ice_shed_tend)

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: pres
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: table_val_qi2qr_melting ! melting
   real(rtype), intent(in) :: table_val_qi2qr_vent_melt ! melting (ventilation term)
   real(rtype), intent(in) :: latent_heat_vapor
   real(rtype), intent(in) :: latent_heat_fusion
   real(rtype), intent(in) :: dv
   real(rtype), intent(in) :: kap
   real(rtype), intent(in) :: mu
   real(rtype), intent(in) :: sc
   real(rtype), intent(in) :: qv
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: qi_incld
   real(rtype), intent(in) :: ni_incld
   real(rtype), intent(in) :: qr_incld

   logical(btype), intent(inout) :: log_wetgrowth
   real(rtype), intent(inout) :: qrcol
   real(rtype), intent(inout) :: qccol
   real(rtype), intent(inout) :: qwgrth
   real(rtype), intent(inout) :: nr_ice_shed_tend
   real(rtype), intent(inout) :: qc2qr_ice_shed_tend

   real(rtype) :: qsat0, dum, dum1

   if (qi_incld.ge.qsmall .and. qc_incld+qr_incld.ge.1.e-6_rtype .and. t_atm.lt.T_zerodegc) then
      qsat0=qv_sat( T_zerodegc,pres,0 )

      qwgrth = ((table_val_qi2qr_melting + table_val_qi2qr_vent_melt*bfb_cbrt(sc)*bfb_sqrt(rhofaci*rho/mu))*       &
      2._rtype*pi*(rho*latent_heat_vapor*dv*(qsat0-qv)-(t_atm-T_zerodegc)*           &
      kap)/(latent_heat_fusion+cpw*(t_atm-T_zerodegc)))*ni_incld

      qwgrth = max(qwgrth,0._rtype)
      dum    = max(0._rtype,(qccol+qrcol)-qwgrth)
      if (dum.ge.1.e-10_rtype) then
         nr_ice_shed_tend = nr_ice_shed_tend + dum*1.923e+6_rtype   ! 1/5.2e-7, 5.2e-7 is the mass of a 1 mm raindrop
         if ((qccol+qrcol).ge.1.e-10_rtype) then
            dum1  = 1._rtype/(qccol+qrcol)
            qc2qr_ice_shed_tend = qc2qr_ice_shed_tend + dum*qccol*dum1
            qccol = qccol - dum*qccol*dum1
            qrcol = qrcol - dum*qrcol*dum1
         endif
         ! densify due to wet growth
         log_wetgrowth = .true.
      endif


   end if

   return

end subroutine ice_cldliq_wet_growth


subroutine calc_ice_relaxation_timescale(rho,t_atm,rhofaci,     &
table_val_qi2qr_melting,table_val_qi2qr_vent_melt,dv,mu,sc,qi_incld,ni_incld,    &
epsi,epsi_tot)

   !-----------------------------
   ! calcualte total inverse ice relaxation timescale combined for all ice categories
   ! note 'f1pr' values are normalized, so we need to multiply by N
   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: table_val_qi2qr_melting ! melting
   real(rtype), intent(in) :: table_val_qi2qr_vent_melt ! melting (ventilation term)
   real(rtype), intent(in) :: dv
   real(rtype), intent(in) :: mu
   real(rtype), intent(in) :: sc
   real(rtype), intent(in) :: qi_incld
   real(rtype), intent(in) :: ni_incld

   real(rtype), intent(out) :: epsi
   real(rtype), intent(inout) :: epsi_tot

   if (qi_incld.ge.qsmall .and. t_atm.lt.T_zerodegc) then
      epsi = ((table_val_qi2qr_melting+table_val_qi2qr_vent_melt*bfb_cbrt(sc)*bfb_sqrt(rhofaci*rho/mu))*2._rtype*pi* &
      rho*dv)*ni_incld
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
      dum1 = revap_table_vals(dumii,dumjj)+(rdumii-real(dumii))*                            &
             (revap_table_vals(dumii+1,dumjj)-revap_table_vals(dumii,dumjj))

      !interoplate value at mu_r+1
      dum2 = revap_table_vals(dumii,dumjj+1)+(rdumii-real(dumii))*                          &
             (revap_table_vals(dumii+1,dumjj+1)-revap_table_vals(dumii,dumjj+1))
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


subroutine calc_rime_density(t_atm,rhofaci,    &
table_val_qi_fallspd,acn,lamc, mu_c,qc_incld,qccol,    &
           vtrmi1,rho_qm_cloud)

   !.........................
   ! calculate rime density

   !     FUTURE:  Add source term for bm (=qccol/rho_qm_cloud) so that all process rates calculations
   !              are done together, before conservation.

   ! NOTE: Tc (ambient) is assumed for the surface temperature.  Technically,
   ! we should diagose graupel surface temperature from heat balance equation.
   ! (but the ambient temperature is a reasonable approximation; tests show
   ! very little sensitivity to different assumed values, Milbrandt and Morrison 2012).

   ! Compute rime density: (based on parameterization of Cober and List, 1993 [JAS])
   ! for simplicty use mass-weighted ice and droplet/rain fallspeeds

   implicit none

   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: table_val_qi_fallspd !mass-weighted fallspeed
   real(rtype), intent(in) :: acn
   real(rtype), intent(in) :: lamc
   real(rtype), intent(in) :: mu_c
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: qccol

   real(rtype), intent(out) :: vtrmi1
   real(rtype), intent(out) :: rho_qm_cloud

   real(rtype) :: iTc, Vt_qc,  D_c, V_impact, Ri

   iTc = 0.0_rtype
   Vt_qc = 0.0_rtype
   D_c  = 0.0_rtype
   V_impact = 0.0_rtype
   Ri = 0.0_rtype

   ! if (qi_incld(i,k).ge.qsmall .and. t_atm(i,k).lt.T_zerodegc) then
   !  NOTE:  condition applicable for cloud only; modify when rain is added back
   if (qccol.ge.qsmall .and. t_atm.lt.T_zerodegc) then
      ! get mass-weighted mean ice fallspeed
      vtrmi1 = table_val_qi_fallspd*rhofaci
      iTc   = 1._rtype/min(-0.001_rtype,t_atm-T_zerodegc)

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
            rho_qm_cloud  = (0.051_rtype + 0.114_rtype*Ri - 0.0055_rtype*bfb_square(Ri))*1000._rtype
         else
            ! for Ri > 8 assume a linear fit between 8 and 12,
            ! rhorime = 900 kg m-3 at Ri = 12
            ! this is somewhat ad-hoc but allows a smoother transition
            ! in rime density up to wet growth
            rho_qm_cloud  = 611._rtype+72.25_rtype*(Ri-8._rtype)
         endif
      else
         rho_qm_cloud = 400._rtype
      endif    !if qc>qsmall
   else
      vtrmi1 = 0._rtype ! no velocity if no ice
      rho_qm_cloud = 400._rtype
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

subroutine cldliq_immersion_freezing(t_atm,lamc,mu_c,cdist1,qc_incld,inv_qc_relvar,    &
           qc2qi_hetero_freeze_tend,nc2ni_immers_freeze_tend)

   !............................................................
   ! contact and immersion freezing droplets

   implicit none

   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: lamc
   real(rtype), intent(in) :: mu_c
   real(rtype), intent(in) :: cdist1
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: inv_qc_relvar

   real(rtype), intent(out) :: qc2qi_hetero_freeze_tend
   real(rtype), intent(out) :: nc2ni_immers_freeze_tend

   real(rtype) :: dum1, dum2, Q_nuc, N_nuc, sbgrd_var_coef

   if (qc_incld.ge.qsmall .and. t_atm.le.T_rainfrz) then
      ! for future: calculate gamma(mu_c+4) in one place since its used multiple times  !AaronDonahue, TODO
      dum1 = bfb_exp(aimm*(T_zerodegc-t_atm))
      dum2 = bfb_cube(1._rtype/lamc)
      sbgrd_var_coef = subgrid_variance_scaling(inv_qc_relvar, 2._rtype)
      Q_nuc = sbgrd_var_coef*cons6*cdist1*bfb_gamma(7._rtype+mu_c)*dum1*bfb_square(dum2)
      N_nuc = cons5*cdist1*bfb_gamma(mu_c+4._rtype)*dum1*dum2
      qc2qi_hetero_freeze_tend = Q_nuc
      nc2ni_immers_freeze_tend = N_nuc
   endif

   return

end subroutine cldliq_immersion_freezing

subroutine rain_immersion_freezing(t_atm,    &
lamr, mu_r, cdistr, qr_incld,    &
qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend)

   !............................................................
   ! immersion freezing of rain
   ! for future: get rid of log statements below for rain freezing

   implicit none

   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: mu_r
   real(rtype), intent(in) :: lamr
   real(rtype), intent(in) :: cdistr
   real(rtype), intent(in) :: qr_incld

   real(rtype), intent(out) :: qr2qi_immers_freeze_tend
   real(rtype), intent(out) :: nr2ni_immers_freeze_tend

   real(rtype) :: Q_nuc, N_nuc

   if (qr_incld.ge.qsmall .and. t_atm.le.T_rainfrz) then

      Q_nuc = cons6*bfb_exp(bfb_log(cdistr) + bfb_log(bfb_gamma(7._rtype+mu_r)) - 6._rtype*bfb_log(lamr))*bfb_exp(aimm*(T_zerodegc-t_atm))
      N_nuc = cons5*bfb_exp(bfb_log(cdistr) + bfb_log(bfb_gamma(mu_r+4._rtype)) - 3._rtype*bfb_log(lamr))*bfb_exp(aimm*(T_zerodegc-t_atm))

      qr2qi_immers_freeze_tend = Q_nuc
      nr2ni_immers_freeze_tend = N_nuc

   endif

   return

end subroutine rain_immersion_freezing


subroutine ice_nucleation(t_atm,inv_rho,ni,ni_activated,qv_supersat_i,inv_dt,do_predict_nc, do_prescribed_CCN,   &
   qinuc,ni_nucleat_tend)

   !................................................................
   ! deposition/condensation-freezing nucleation
   ! allow ice nucleation if < -15 C and > 5% ice supersaturation
   ! use CELL-AVERAGE values, freezing of vapor

   implicit none

   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: inv_rho
   real(rtype), intent(in) :: ni
   real(rtype), intent(in) :: ni_activated
   real(rtype), intent(in) :: qv_supersat_i
   real(rtype), intent(in) :: inv_dt
   logical(btype), intent(in) :: do_predict_nc, do_prescribed_CCN

   real(rtype), intent(inout) :: qinuc
   real(rtype), intent(inout) :: ni_nucleat_tend


   real(rtype) :: dum, N_nuc, Q_nuc

   if ( t_atm .lt.T_icenuc .and. qv_supersat_i.ge.0.05_rtype ) then
      if(.not. do_predict_nc .or. do_prescribed_CCN) then
!         ! dum = exp(-0.639+0.1296*100.*qv_supersat_i(i,k))*1000.*inv_rho(i,k)  !Meyers et al. (1992)
         dum = 0.005_rtype*bfb_exp(0.304_rtype*(T_zerodegc-t_atm))*1000._rtype*inv_rho   !Cooper (1986)
         dum = min(dum,100.e3_rtype*inv_rho)
         N_nuc = max(0._rtype,(dum-ni)*inv_dt)
         if (N_nuc.ge.1.e-20_rtype) then
            Q_nuc = max(0._rtype,(dum-ni)*mi0*inv_dt)
            qinuc = Q_nuc
            ni_nucleat_tend = N_nuc
         endif
      else
      ! Ice nucleation predicted by the aerosol scheme
      ! +++ [JS]: Here we add the cirrus/aerosol contribution to the Cooper scheme under dedicated switch +++
            if(do_Cooper_inP3)then
               if(t_atm > 236.18_rtype) then
                  dum = 0.005_rtype*bfb_exp(0.304_rtype*(T_zerodegc-t_atm))*1000._rtype*inv_rho   !Cooper (1986)
                  dum = min(dum,100.e3_rtype*inv_rho)
                  N_nuc = max(0._rtype,(dum-ni)*inv_dt)
                  if (N_nuc.ge.1.e-20_rtype) then
                     Q_nuc = max(0._rtype,(dum-ni)*mi0*inv_dt)
                     qinuc = Q_nuc
                     ni_nucleat_tend = N_nuc
                  endif
               endif   
               ni_nucleat_tend = ni_nucleat_tend + max(0._rtype, (ni_activated - ni)*inv_dt)
               qinuc = ni_nucleat_tend * mi0
            else          
               ni_nucleat_tend = max(0._rtype, (ni_activated - ni)*inv_dt)
               qinuc = ni_nucleat_tend * mi0
            endif   
      endif
   endif

end subroutine ice_nucleation

subroutine droplet_self_collection(rho,inv_rho,qc_incld,mu_c,nu,nc2nr_autoconv_tend,    &
   nc_selfcollect_tend)

   !............................
   ! self-collection of droplets

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: inv_rho
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: mu_c
   real(rtype), intent(in) :: nu
   real(rtype), intent(in) :: nc2nr_autoconv_tend

   real(rtype), intent(out) :: nc_selfcollect_tend

   if (qc_incld.ge.qsmall) then

      if (iparam.eq.1) then
         !Seifert and Beheng (2001)
         nc_selfcollect_tend = -kc*(1.e-3_rtype*rho*qc_incld)**2*(nu+2._rtype)/(nu+1._rtype)*         &
              1.e+6_rtype*inv_rho+nc2nr_autoconv_tend
      elseif (iparam.eq.2) then
         !Beheng (994)
         nc_selfcollect_tend = -5.5e+16_rtype*inv_rho*mu_c**(-0.63_rtype)*(1.e-3_rtype*rho*qc_incld)**2
      elseif (iparam.eq.3) then
         !Khroutdinov and Kogan (2000)
         nc_selfcollect_tend = 0._rtype
      endif

   endif

end subroutine droplet_self_collection

subroutine cloud_rain_accretion(rho,inv_rho,qc_incld,nc_incld,qr_incld,inv_qc_relvar,    &
  p3_accret_coeff,p3_qc_accret_expon,     &
  qc2qr_accret_tend,nc_accret_tend)

  !............................
  ! accretion of cloud by rain

  implicit none

  real(rtype), intent(in) :: rho
  real(rtype), intent(in) :: inv_rho
  real(rtype), intent(in) :: qc_incld
  real(rtype), intent(in) :: nc_incld
  real(rtype), intent(in) :: qr_incld
  real(rtype), intent(in) :: inv_qc_relvar
  real(rtype), intent(in) :: p3_accret_coeff
  real(rtype), intent(in) :: p3_qc_accret_expon

  real(rtype), intent(out) :: qc2qr_accret_tend
  real(rtype), intent(out) :: nc_accret_tend

  real(rtype) :: dum, dum1, sbgrd_var_coef

  if (qr_incld.ge.qsmall .and. qc_incld.ge.qsmall) then

     if (iparam.eq.1) then
        !Seifert and Beheng (2001)
        dum   = 1._rtype-qc_incld/(qc_incld+qr_incld)
        dum1  = (dum/(dum+5.e-4_rtype))**4
        qc2qr_accret_tend = kr*rho*0.001_rtype*qc_incld*qr_incld*dum1
        nc_accret_tend = qc2qr_accret_tend*rho*0.001_rtype*(nc_incld*rho*1.e-6_rtype)/(qc_incld*rho*   &
             0.001_rtype)*1.e+6_rtype*inv_rho
     elseif (iparam.eq.2) then
        !Beheng (994)
        qc2qr_accret_tend = 6._rtype*rho*(qc_incld*qr_incld)
        nc_accret_tend = qc2qr_accret_tend*rho*1.e-3_rtype*(nc_incld*rho*1.e-6_rtype)/(qc_incld*rho*1.e-3_rtype)* &
             1.e+6_rtype*inv_rho
     elseif (iparam.eq.3) then
        !Khroutdinov and Kogan (2000)
        !print*,'p3_qc_accret_expon = ',p3_qc_accret_expon
        sbgrd_var_coef = subgrid_variance_scaling(inv_qc_relvar, p3_qc_accret_expon)
        qc2qr_accret_tend = sbgrd_var_coef*p3_accret_coeff*bfb_pow(qc_incld*qr_incld, p3_qc_accret_expon)
        nc_accret_tend = qc2qr_accret_tend*nc_incld/qc_incld
     endif

     if (qc2qr_accret_tend.eq.0._rtype) nc_accret_tend = 0._rtype
     if (nc_accret_tend.eq.0._rtype) qc2qr_accret_tend = 0._rtype

  endif

end subroutine cloud_rain_accretion

subroutine rain_self_collection(rho,qr_incld,nr_incld,    &
   nr_selfcollect_tend)

   !.....................................
   ! self-collection and breakup of rain
   ! (breakup following modified Verlinde and Cotton scheme)

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: qr_incld
   real(rtype), intent(in) :: nr_incld
   real(rtype), intent(out) :: nr_selfcollect_tend

   real(rtype) :: dum, dum1, dum2

   if (qr_incld.ge.qsmall) then

      ! include breakup
      dum1 = 280.e-6_rtype

      ! use mass-mean diameter (do this by using
      ! the old version of lambda w/o mu dependence)
      ! note there should be a factor of 6^(1/3), but we
      ! want to keep breakup threshold consistent so 'dum'
      ! is expressed in terms of lambda rather than mass-mean D

      dum2 = bfb_cbrt(qr_incld/(pi*rho_h2o*nr_incld))
      if (dum2.lt.dum1) then
         dum = 1._rtype
      else if (dum2.ge.dum1) then
         dum = 2._rtype - bfb_exp(2300._rtype*(dum2-dum1))
      endif

      if (iparam.eq.1) then
         nr_selfcollect_tend = dum*kr*1.e-3_rtype*qr_incld*nr_incld*rho
      elseif (iparam.eq.2 .or. iparam.eq.3) then
         nr_selfcollect_tend = dum*5.78_rtype*nr_incld*qr_incld*rho
      endif

   endif

end subroutine rain_self_collection


subroutine cloud_water_autoconversion(rho,qc_incld,nc_incld,inv_qc_relvar,          &
   p3_autocon_coeff,p3_qc_autocon_expon,p3_nc_autocon_expon,p3_embryonic_rain_size, &
   do_precip_off,qc2qr_autoconv_tend,nc2nr_autoconv_tend,ncautr)

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: qc_incld
   real(rtype), intent(in) :: nc_incld
   real(rtype), intent(in) :: inv_qc_relvar
   real(rtype), intent(in) :: p3_autocon_coeff
   real(rtype), intent(in) :: p3_qc_autocon_expon
   real(rtype), intent(in) :: p3_nc_autocon_expon
   real(rtype), intent(in) :: p3_embryonic_rain_size

   logical(btype), intent(in) :: do_precip_off

   real(rtype), intent(out) :: qc2qr_autoconv_tend
   real(rtype), intent(out) :: nc2nr_autoconv_tend
   real(rtype), intent(out) :: ncautr

   real(rtype) :: sbgrd_var_coef

   qc_not_small: if (qc_incld.ge.1.e-8_rtype) then

      !Khroutdinov and Kogan (2000)
      !print*,'p3_qc_autocon_expon = ',p3_qc_autocon_expon

      sbgrd_var_coef = subgrid_variance_scaling(inv_qc_relvar, p3_qc_autocon_expon)
      qc2qr_autoconv_tend = sbgrd_var_coef*p3_autocon_coeff*bfb_pow(qc_incld,p3_qc_autocon_expon)*bfb_pow(nc_incld*1.e-6_rtype*rho,p3_nc_autocon_expon)
      
      ! Note: ncautr is change in Nr; nc2nr_autoconv_tend is change in Nc
      ncautr = qc2qr_autoconv_tend*cons3*(1._rtype/bfb_pow(p3_embryonic_rain_size,3._rtype))
      nc2nr_autoconv_tend = qc2qr_autoconv_tend*nc_incld/qc_incld

      if (qc2qr_autoconv_tend .eq.0._rtype .or. do_precip_off) nc2nr_autoconv_tend = 0._rtype
      if (nc2nr_autoconv_tend.eq.0._rtype .or. do_precip_off) qc2qr_autoconv_tend  = 0._rtype
      if (do_precip_off) ncautr = 0._rtype

   endif qc_not_small

end subroutine cloud_water_autoconversion

subroutine back_to_cell_average(cld_frac_l,cld_frac_r,cld_frac_i,                                                       &
   qc2qr_accret_tend,qr2qv_evap_tend,qc2qr_autoconv_tend,                                                               &
   nc_accret_tend,nc_selfcollect_tend,nc2nr_autoconv_tend,nr_selfcollect_tend,nr_evap_tend,ncautr,qi2qv_sublim_tend,    &
   nr_ice_shed_tend,qc2qi_hetero_freeze_tend, qrcol,qc2qr_ice_shed_tend,qi2qr_melt_tend,qccol,qr2qi_immers_freeze_tend, &
   ni2nr_melt_tend,nc_collect_tend,ncshdc,nc2ni_immers_freeze_tend,nr_collect_tend,ni_selfcollect_tend,                 &
   qidep,nr2ni_immers_freeze_tend,ni_sublim_tend,qinuc,ni_nucleat_tend,qiberg,                                          &
   ncheti_cnt, qcheti_cnt, nicnt, qicnt, ninuc_cnt, qinuc_cnt)

   ! Here we map the microphysics tendency rates back to CELL-AVERAGE quantities for updating
   ! cell-average quantities.

   implicit none

   ! Intersection of cloud fractions for combination of ice (i), rain (r) and liquid (l)
   real(rtype), intent(in) :: cld_frac_l
   real(rtype), intent(in) :: cld_frac_r
   real(rtype), intent(in) :: cld_frac_i

   real(rtype), intent(inout) :: qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend, nc_accret_tend, &
        nc_selfcollect_tend, nc2nr_autoconv_tend, nr_selfcollect_tend, nr_evap_tend, ncautr
   real(rtype), intent(inout) :: qi2qv_sublim_tend, nr_ice_shed_tend, qc2qi_hetero_freeze_tend, qrcol, qc2qr_ice_shed_tend, &
        qi2qr_melt_tend, qccol, qr2qi_immers_freeze_tend, ni2nr_melt_tend, &
        nc_collect_tend, ncshdc, nc2ni_immers_freeze_tend, nr_collect_tend, ni_selfcollect_tend, qidep
   real(rtype), intent(inout) :: nr2ni_immers_freeze_tend, ni_sublim_tend, qinuc, ni_nucleat_tend, qiberg
   real(rtype), intent(inout) :: ncheti_cnt,qcheti_cnt,nicnt,qicnt,ninuc_cnt,qinuc_cnt

   real(rtype) :: ir_cldm, il_cldm, lr_cldm

   ir_cldm = min(cld_frac_i,cld_frac_r)  ! Intersection of ICE and RAIN cloud
   il_cldm = min(cld_frac_i,cld_frac_l)  ! Intersection of ICE and LIQUID cloud
   lr_cldm = min(cld_frac_l,cld_frac_r)  ! Intersection of LIQUID and RAIN cloud

   ! Some process rates take place within the intersection of liquid, rain and ice cloud fractions.
   ! We calculate the intersection as the minimum between combinations of cloud fractions and use
   ! these values to map back to cell-average quantities where applicable.

          ! map warm-phase process rates to cell-avg
   qc2qr_accret_tend   = qc2qr_accret_tend*lr_cldm     ! Accretion of liquid to rain
   qr2qv_evap_tend   = qr2qv_evap_tend*cld_frac_r       ! Evaporation of rain
   qc2qr_autoconv_tend   = qc2qr_autoconv_tend*cld_frac_l       ! Autoconversion of liquid
   nc_accret_tend   = nc_accret_tend*lr_cldm     ! Number change due to accretion
   nc_selfcollect_tend   = nc_selfcollect_tend*cld_frac_l       ! Self collection occurs locally in liq. cloud
   nc2nr_autoconv_tend  = nc2nr_autoconv_tend*cld_frac_l      ! Impact of autoconversion on number
   nr_selfcollect_tend   = nr_selfcollect_tend*cld_frac_r       ! Self collection occurs locally in rain cloud
   nr_evap_tend   = nr_evap_tend*cld_frac_r       ! Change in rain number due to evaporation
   ncautr  = ncautr*lr_cldm    ! Autoconversion of rain drops within rain/liq cloud

   ! map ice-phase  process rates to cell-avg
   qi2qv_sublim_tend   = qi2qv_sublim_tend*cld_frac_i       ! Sublimation of ice in ice cloud
   nr_ice_shed_tend  = nr_ice_shed_tend*il_cldm    ! Rain # increase due to shedding from rain-ice collisions, occurs when ice and liquid interact
   qc2qi_hetero_freeze_tend  = qc2qi_hetero_freeze_tend*il_cldm    ! Immersion freezing of cloud drops
   qrcol   = qrcol*ir_cldm     ! Collection of rain mass by ice
   qc2qr_ice_shed_tend   = qc2qr_ice_shed_tend*il_cldm     ! Rain mass growth due to shedding of fain drops after collisions with ice, occurs when ice and liquid interact
   qi2qr_melt_tend   = qi2qr_melt_tend*cld_frac_i       ! Melting of ice
   qccol   = qccol*il_cldm     ! Collection of water by ice
   qr2qi_immers_freeze_tend  = qr2qi_immers_freeze_tend*cld_frac_r      ! Immersion freezing of rain
   ni2nr_melt_tend   = ni2nr_melt_tend*cld_frac_i       ! Change in number due to melting
   nc_collect_tend   = nc_collect_tend*il_cldm     ! Cloud # change due to collection of cld water by ice
   ncshdc  = ncshdc*il_cldm    ! Number change due to shedding, occurs when ice and liquid interact
   nc2ni_immers_freeze_tend  = nc2ni_immers_freeze_tend*cld_frac_l      ! Number change associated with freexzing of cld drops
   nr_collect_tend   = nr_collect_tend*ir_cldm     ! Rain number change due to collection from ice
   ni_selfcollect_tend   = ni_selfcollect_tend*cld_frac_i       ! Ice self collection
   qidep   = qidep*cld_frac_i       ! Vapor deposition to ice phase
   nr2ni_immers_freeze_tend  = nr2ni_immers_freeze_tend*cld_frac_r      ! Change in number due to immersion freezing of rain
   ni_sublim_tend   = ni_sublim_tend*cld_frac_i       ! Number change due to sublimation of ice
   qiberg  = qiberg*il_cldm    ! Bergeron process
     ! AaronDonahue: These variables are related to aerosol activation and their usage will be changed in a later PR.
   qinuc   = qinuc             ! Deposition and condensation-freezing nucleation, already cell-averaged
   ni_nucleat_tend   = ni_nucleat_tend             ! Number change due to deposition and condensation-freezing, already cell-averaged

   ncheti_cnt = ncheti_cnt*cld_frac_l
   qcheti_cnt = qcheti_cnt*cld_frac_l
   nicnt = nicnt*cld_frac_l
   qicnt = qicnt*cld_frac_l
   ninuc_cnt = ninuc_cnt*cld_frac_l
   qinuc_cnt = qinuc_cnt*cld_frac_l

end subroutine back_to_cell_average

subroutine prevent_ice_overdepletion(pres,t_atm,qv,latent_heat_vapor,latent_heat_sublim,inv_dt,dt,qidep,qinuc,qinuc_cnt,    &
   qi2qv_sublim_tend,qr2qv_evap_tend)

   !-- Limit ice process rates to prevent overdepletion of sources such that
   !   the subsequent adjustments are done with maximum possible rates for the
   !   time step.  (note: most ice rates are adjusted here since they must be done
   !   simultaneously (outside of iice-loops) to distribute reduction proportionally
   !   amongst categories.
   !PMC - might need to rethink above statement since only one category now.

   implicit none

   real(rtype), intent(in) :: pres
   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: qv
   real(rtype), intent(in) :: latent_heat_sublim,latent_heat_vapor
   real(rtype), intent(in) :: inv_dt,dt
   real(rtype), intent(in) :: qidep,qinuc,qinuc_cnt

   real(rtype), intent(inout) :: qi2qv_sublim_tend,qr2qv_evap_tend

   real(rtype) :: dumqv_sat_i, qdep_satadj, qrevp_satadj
   real(rtype) :: qtmp_all,ttmp_all,qv_sat_l,q_sink,t_sink,dumqv_sat_l,qv_source_evp,qi2qv_sublim_satadj


   qtmp_all = qv - (qidep + qinuc + qinuc_cnt)*dt + (qi2qv_sublim_tend + qr2qv_evap_tend)*dt 
   ttmp_all = t_atm + ((qidep-qi2qv_sublim_tend+qinuc+qinuc_cnt)*latent_heat_sublim*inv_cp + (-qr2qv_evap_tend*latent_heat_vapor*inv_cp))*dt

   qv_sat_l = qv_sat(ttmp_all,pres,0)

   ! If water vapor mass exceeds ice saturation value, we limit only source terms (e.g., sublimation, rain evp) 
   if(qtmp_all > qv_sat_l)then

      ! ... First, rain evaporation is limited to the sub-saturation defined by the water vapor sink terms (deposition, ice nucleation)
      q_sink = qv - (qidep + qinuc + qinuc_cnt)*dt
      t_sink = t_atm + ((qidep+qinuc+qinuc_cnt)*latent_heat_sublim*inv_cp)*dt
      dumqv_sat_l = qv_sat(t_sink,pres,0)
      qv_source_evp = qr2qv_evap_tend + qi2qv_sublim_tend
      qrevp_satadj = (q_sink-dumqv_sat_l)/(1._rtype + bfb_square(latent_heat_vapor)*dumqv_sat_l/(cp*rv* bfb_square(t_sink) ))*inv_dt
      qr2qv_evap_tend    = qr2qv_evap_tend*min(1._rtype,max(0._rtype,-qrevp_satadj)/max(qv_source_evp, 1.e-20_rtype))

      ! ... Next, ice-sublimation is limited in the same way but with sublim LH 
      dumqv_sat_i = qv_sat(t_sink,pres,1)
      qi2qv_sublim_satadj = (q_sink-dumqv_sat_i)/(1._rtype + bfb_square(latent_heat_sublim)*dumqv_sat_i/(cp*rv* bfb_square(t_sink) ))*inv_dt
      qi2qv_sublim_tend  = qi2qv_sublim_tend*min(1._rtype,max(0._rtype,-qi2qv_sublim_satadj)/max(qv_source_evp, 1.e-20_rtype))
   
   endif


end subroutine prevent_ice_overdepletion

subroutine water_vapor_conservation(qv,qidep,qinuc,qi2qv_sublim_tend,qr2qv_evap_tend,qinuc_cnt,dt)
  ! If water vapor sinks cause qv<0 by the end of the step, rescale them such that qv=0 at the end of the step

  implicit none

  real(rtype), intent(in)     :: qv,qi2qv_sublim_tend,qr2qv_evap_tend,dt
  real(rtype), intent(inout) :: qidep,qinuc,qinuc_cnt
  real(rtype)                :: qv_avail, qv_sink, ratio

  qv_avail = qv + (qi2qv_sublim_tend+qr2qv_evap_tend)*dt
  if(use_hetfrz_classnuc)then
   qv_sink = (qidep + qinuc + qinuc_cnt)*dt
  else
   qv_sink = (qidep + qinuc)*dt
  endif

  if (qv_sink > qv_avail .and. qv_sink>1.e-20_rtype) then
     ratio = qv_avail/qv_sink
     qidep = qidep*ratio
     qinuc = qinuc*ratio
     if(use_hetfrz_classnuc)then
      qinuc_cnt = qinuc_cnt*ratio
     endif 
  endif

  return
end subroutine water_vapor_conservation

subroutine ice_supersat_conservation(qidep,qinuc,qi2qv_sublim_tend,qr2qv_evap_tend,qinuc_cnt,cld_frac_i,qv,qv_sat_i,latent_heat_sublim,T_atm,dt)
  !Make sure ice processes don't drag qv below ice supersaturation

  implicit none

  real(rtype), intent(in) :: cld_frac_i,qv,qv_sat_i,latent_heat_sublim,T_atm,dt,qi2qv_sublim_tend,qr2qv_evap_tend
  real(rtype), intent(inout) :: qidep,qinuc,qinuc_cnt

  real(rtype) :: qv_sink, qv_avail, fract

   if(use_hetfrz_classnuc)then
      qv_sink = qidep + qinuc + qinuc_cnt   ! in [kg/kg] cell-avg values
   else
      qv_sink = qidep + qinuc               ! in [kg/kg] cell-avg values
  endif

  if (qv_sink > qsmall .and. cld_frac_i > 1.0e-20_rtype) then
     ! --- Available water vapor for deposition/nucleation
   qv_avail = (qv + (qi2qv_sublim_tend + qr2qv_evap_tend)*dt - qv_sat_i) / &
          (1.0_rtype + bfb_square(latent_heat_sublim)*qv_sat_i / (cp*rv*bfb_square(T_atm)) ) / dt

     ! --- Only excess water vapor can be limited
     qv_avail = max(qv_avail,0.0_rtype)

     if (qv_sink >  qv_avail) then
        fract = qv_avail / qv_sink
        if(use_hetfrz_classnuc)then
         qinuc_cnt = qinuc_cnt * fract
        endif
        qinuc = qinuc * fract
        qidep = qidep * fract
     endif
  endif

  return
end subroutine ice_supersat_conservation

subroutine nc_conservation(nc, nc_selfcollect_tend, dt, nc_collect_tend, nc2ni_immers_freeze_tend, &
     nc_accret_tend, nc2nr_autoconv_tend, ncheti_cnt, nicnt)
  !Make sure sinks of nc don't force end-of-step nc below 0. Rescale them if they do.

  implicit none

  real(rtype), intent(in) :: nc,nc_selfcollect_tend,dt
  real(rtype), intent(inout) :: nc_collect_tend,nc2ni_immers_freeze_tend,&
                                nc_accret_tend,nc2nr_autoconv_tend,ncheti_cnt,nicnt
  real(rtype) :: sink_nc, source_nc, ratio

  if(use_hetfrz_classnuc)then
      sink_nc = (nc_collect_tend + ncheti_cnt + nc_accret_tend + nc2nr_autoconv_tend + nicnt)*dt
  else
      sink_nc = (nc_collect_tend + nc2ni_immers_freeze_tend + nc_accret_tend + nc2nr_autoconv_tend)*dt
  endif 
  
  source_nc = nc + nc_selfcollect_tend*dt
  if(sink_nc > source_nc) then
     ratio = source_nc/sink_nc
     nc_collect_tend  = nc_collect_tend*ratio
     nc_accret_tend  = nc_accret_tend*ratio
     nc2nr_autoconv_tend = nc2nr_autoconv_tend*ratio
     if(use_hetfrz_classnuc)then
      ncheti_cnt = ncheti_cnt*ratio
      nicnt = nicnt*ratio
     else
       nc2ni_immers_freeze_tend = nc2ni_immers_freeze_tend*ratio
     endif
  endif

  return
end subroutine nc_conservation

subroutine nr_conservation(nr,ni2nr_melt_tend,nr_ice_shed_tend,ncshdc,nc2nr_autoconv_tend,dt,nr_collect_tend,nmltratio,&
     nr2ni_immers_freeze_tend,nr_selfcollect_tend,nr_evap_tend)
  !Make sure sinks of nr don't force end-of-step nr below 0. Rescale them if they do.

  implicit none

  real(rtype), intent(in) :: nr,ni2nr_melt_tend,nr_ice_shed_tend,ncshdc,nc2nr_autoconv_tend,dt,nmltratio
  real(rtype), intent(inout) :: nr_collect_tend,nr2ni_immers_freeze_tend,nr_selfcollect_tend,nr_evap_tend
  real(rtype) :: sink_nr, source_nr, ratio

  sink_nr = (nr_collect_tend + nr2ni_immers_freeze_tend + nr_selfcollect_tend + nr_evap_tend)*dt
  source_nr = nr + (ni2nr_melt_tend*nmltratio + nr_ice_shed_tend + ncshdc + nc2nr_autoconv_tend)*dt
  if(sink_nr > source_nr) then
     ratio = source_nr/sink_nr
     nr_collect_tend  = nr_collect_tend*ratio
     nr2ni_immers_freeze_tend = nr2ni_immers_freeze_tend*ratio
     nr_selfcollect_tend  = nr_selfcollect_tend*ratio
     nr_evap_tend  = nr_evap_tend*ratio
  endif

  return
end subroutine nr_conservation

subroutine ni_conservation(ni, ni_nucleat_tend, nr2ni_immers_freeze_tend, nc2ni_immers_freeze_tend, ncheti_cnt, nicnt, ninuc_cnt, dt, &
     ni2nr_melt_tend,ni_sublim_tend,ni_selfcollect_tend)
  !Make sure ni doesn't go below zero

  implicit none

  real(rtype), intent(in) :: ni,ni_nucleat_tend,nr2ni_immers_freeze_tend,nc2ni_immers_freeze_tend,dt, &
                             ncheti_cnt,nicnt,ninuc_cnt
  real(rtype), intent(inout) :: ni2nr_melt_tend,ni_sublim_tend,ni_selfcollect_tend
  real(rtype) :: sink_ni, source_ni, ratio

  sink_ni = (ni2nr_melt_tend + ni_sublim_tend + ni_selfcollect_tend)*dt
  if(use_hetfrz_classnuc)then
      source_ni = ni + (ni_nucleat_tend+nr2ni_immers_freeze_tend+ncheti_cnt+nicnt+ninuc_cnt)*dt
   else
      source_ni = ni + (ni_nucleat_tend+nr2ni_immers_freeze_tend+nc2ni_immers_freeze_tend)*dt
  endif
  if(sink_ni > source_ni) then
     ratio = source_ni/sink_ni
     ni2nr_melt_tend  = ni2nr_melt_tend*ratio
     ni_sublim_tend = ni_sublim_tend*ratio
     ni_selfcollect_tend = ni_selfcollect_tend*ratio
  endif

  return
end subroutine ni_conservation

subroutine cloud_water_conservation(qc,dt,    &
   qc2qr_autoconv_tend,qc2qr_accret_tend,qccol,qc2qi_hetero_freeze_tend,qc2qr_ice_shed_tend,qiberg,qi2qv_sublim_tend,qidep, &
   qcheti_cnt,qicnt)

   implicit none

   real(rtype), intent(in) :: qc, dt
   real(rtype), intent(inout) :: qc2qr_autoconv_tend, qc2qr_accret_tend, qccol, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, &
                                 qiberg, qi2qv_sublim_tend, qidep, qcheti_cnt, qicnt
   real(rtype) :: sinks, ratio

   if(use_hetfrz_classnuc)then
      sinks = (qc2qr_autoconv_tend+qc2qr_accret_tend+qccol+qcheti_cnt+qc2qr_ice_shed_tend+qiberg+qicnt)*dt
   else
      sinks = (qc2qr_autoconv_tend+qc2qr_accret_tend+qccol+qc2qi_hetero_freeze_tend+qc2qr_ice_shed_tend+qiberg)*dt
   endif

   if (sinks .gt. qc .and. sinks.ge.1.e-20_rtype) then
      ratio  = qc/sinks
      qc2qr_autoconv_tend  = qc2qr_autoconv_tend*ratio
      qc2qr_accret_tend  = qc2qr_accret_tend*ratio
      qccol  = qccol*ratio
      qc2qr_ice_shed_tend  = qc2qr_ice_shed_tend*ratio
      qiberg = qiberg*ratio
      if(use_hetfrz_classnuc)then
         qcheti_cnt = qcheti_cnt*ratio
         qicnt = qicnt*ratio
      else
         qc2qi_hetero_freeze_tend = qc2qi_hetero_freeze_tend*ratio
      endif   
   else
      ratio = 1.0 ! If not limiting sinks on qc then most likely did not run out of qc
   endif

   !PMC: ratio is also frac of step w/ liq. thus we apply qiberg for
   !"ratio" of timestep and vapor deposition and sublimation  for the
   !remaining frac of the timestep.  Only limit if there will be cloud
   !water to begin with.
   if (qc .gt. 1.e-20_rtype) then
      qidep  = qidep*(1._rtype-ratio)
      qi2qv_sublim_tend  = qi2qv_sublim_tend*(1._rtype-ratio)
   end if


end subroutine cloud_water_conservation

subroutine rain_water_conservation(qr,qc2qr_autoconv_tend,qc2qr_accret_tend,qi2qr_melt_tend,qc2qr_ice_shed_tend,dt,    &
   qr2qv_evap_tend,qrcol,qr2qi_immers_freeze_tend)

   implicit none

   real(rtype), intent(in) :: qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, dt
   real(rtype), intent(inout) :: qr2qv_evap_tend, qrcol, qr2qi_immers_freeze_tend

   real(rtype) :: sinks, sources, ratio

   sinks   = (qr2qv_evap_tend+qrcol+qr2qi_immers_freeze_tend)*dt
   sources = qr + (qc2qr_autoconv_tend+qc2qr_accret_tend+qi2qr_melt_tend+qc2qr_ice_shed_tend)*dt
   if (sinks.gt.sources .and. sinks.ge.1.e-20_rtype) then
      ratio  = sources/sinks
      qr2qv_evap_tend  = qr2qv_evap_tend*ratio
      qrcol  = qrcol*ratio
      qr2qi_immers_freeze_tend = qr2qi_immers_freeze_tend*ratio
   endif

end subroutine rain_water_conservation

subroutine ice_water_conservation(qi,qidep,qinuc,qiberg,qrcol,qccol,qr2qi_immers_freeze_tend,qc2qi_hetero_freeze_tend,dt,qinuc_cnt,qcheti_cnt,qicnt,    &
                                  qi2qv_sublim_tend,qi2qr_melt_tend)

   implicit none

   real(rtype), intent(in) :: qi, qidep, qinuc, qrcol, qccol, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, qiberg, dt, &
                              qinuc_cnt, qcheti_cnt,qicnt
   real(rtype), intent(inout) :: qi2qv_sublim_tend, qi2qr_melt_tend
   real(rtype) :: sinks, sources, ratio

   sinks   = (qi2qv_sublim_tend+qi2qr_melt_tend)*dt

   if(use_hetfrz_classnuc)then
      sources = qi + (qidep+qinuc+qrcol+qccol+qr2qi_immers_freeze_tend+qiberg+qinuc_cnt+qcheti_cnt+qicnt)*dt
   else
      sources = qi + (qidep+qinuc+qrcol+qccol+qr2qi_immers_freeze_tend+qc2qi_hetero_freeze_tend+qiberg)*dt
   endif

   if (sinks.gt.sources .and. sinks.ge.1.e-20_rtype) then
      ratio = sources/sinks
      qi2qv_sublim_tend = qi2qv_sublim_tend*ratio
      qi2qr_melt_tend = qi2qr_melt_tend*ratio
   endif

end subroutine ice_water_conservation


subroutine update_prognostic_ice(qc2qi_hetero_freeze_tend,qccol,qc2qr_ice_shed_tend,   &
   nc_collect_tend,nc2ni_immers_freeze_tend,ncshdc,                                    &
   qrcol,nr_collect_tend,qr2qi_immers_freeze_tend,nr2ni_immers_freeze_tend,nr_ice_shed_tend,                                   &
   qi2qr_melt_tend,ni2nr_melt_tend,qi2qv_sublim_tend,qidep,qinuc,ni_nucleat_tend,ni_selfcollect_tend,ni_sublim_tend,qiberg,    &
   exner,latent_heat_sublim,latent_heat_fusion,                                        &
   do_predict_nc,log_wetgrowth,dt,nmltratio,rho_qm_cloud,                              &
   ncheti_cnt, nicnt, ninuc_cnt, qcheti_cnt, qicnt, qinuc_cnt,                         &
   th_atm,qv,qi,ni,qm,bm,qc,nc,qr,nr,qi_wetDepos)

   !-- ice-phase dependent processes:
   implicit none

   real(rtype), intent(in) :: qc2qi_hetero_freeze_tend
   real(rtype), intent(in) :: qccol
   real(rtype), intent(in) :: qc2qr_ice_shed_tend
   real(rtype), intent(in) :: nc_collect_tend
   real(rtype), intent(in) :: nc2ni_immers_freeze_tend
   real(rtype), intent(in) :: ncshdc

   real(rtype), intent(in) :: qrcol
   real(rtype), intent(in) :: nr_collect_tend
   real(rtype), intent(in) :: qr2qi_immers_freeze_tend
   real(rtype), intent(in) :: nr2ni_immers_freeze_tend
   real(rtype), intent(in) :: nr_ice_shed_tend

   real(rtype), intent(in) :: qi2qr_melt_tend
   real(rtype), intent(in) :: ni2nr_melt_tend
   real(rtype), intent(in) :: qi2qv_sublim_tend
   real(rtype), intent(in) :: qidep
   real(rtype), intent(in) :: qinuc
   real(rtype), intent(in) :: ni_nucleat_tend
   real(rtype), intent(in) :: ni_selfcollect_tend
   real(rtype), intent(in) :: ni_sublim_tend
   real(rtype), intent(in) :: qiberg
   real(rtype), intent(in) :: exner
   real(rtype), intent(in) :: latent_heat_fusion
   real(rtype), intent(in) :: latent_heat_sublim

   logical(btype), intent(in) :: do_predict_nc
   logical(btype), intent(in) :: log_wetgrowth
   real(rtype), intent(in) :: dt
   real(rtype), intent(in) :: nmltratio
   real(rtype), intent(in) :: rho_qm_cloud
   real(rtype), intent(in) :: ncheti_cnt, nicnt, ninuc_cnt, qcheti_cnt, qicnt, qinuc_cnt

   real(rtype), intent(inout) :: th_atm,qi_wetDepos
   real(rtype), intent(inout) :: qv
   real(rtype), intent(inout) :: qc
   real(rtype), intent(inout) :: nc
   real(rtype), intent(inout) :: qr
   real(rtype), intent(inout) :: nr
   real(rtype), intent(inout) :: qi
   real(rtype), intent(inout) :: ni
   real(rtype), intent(inout) :: qm
   real(rtype), intent(inout) :: bm

   real(rtype) :: dum

   if(.not. use_hetfrz_classnuc)then
      qc = qc + (-qc2qi_hetero_freeze_tend-qccol-qc2qr_ice_shed_tend-qiberg)*dt
   else
      qc = qc + (-qcheti_cnt-qicnt-qccol-qc2qr_ice_shed_tend-qiberg)*dt
   endif

   if (do_predict_nc) then
      if(.not. use_hetfrz_classnuc)then
         nc = nc + (-nc_collect_tend-nc2ni_immers_freeze_tend)*dt
      else
         nc = nc + (-nc_collect_tend-ncheti_cnt-nicnt)*dt
      endif   
   endif

   qr = qr + (-qrcol+qi2qr_melt_tend-qr2qi_immers_freeze_tend+qc2qr_ice_shed_tend)*dt

   ! apply factor to source for rain number from melting of ice, (ad-hoc
   ! but accounts for rapid evaporation of small melting ice particles)
   nr = nr + (-nr_collect_tend-nr2ni_immers_freeze_tend+nmltratio*ni2nr_melt_tend+nr_ice_shed_tend+ncshdc)*dt

   if (qi.ge.qsmall) then
      ! add sink terms, assume density stays constant for sink terms
      bm = bm - ((qi2qv_sublim_tend+qi2qr_melt_tend)/qi)*dt*bm
      qm = qm - ((qi2qv_sublim_tend+qi2qr_melt_tend)*qm/qi)*dt
      qi = qi - (qi2qv_sublim_tend+qi2qr_melt_tend)*dt
   endif

   if(.not. use_hetfrz_classnuc)then
      dum = (qrcol+qccol+qr2qi_immers_freeze_tend+qc2qi_hetero_freeze_tend)*dt
      qi = qi + (qidep+qinuc+qiberg)*dt + dum
      qi_wetDepos = (qidep+qinuc+qiberg) + (qrcol+qccol+qr2qi_immers_freeze_tend+qc2qi_hetero_freeze_tend) ! total ice process rate for wet-deopsition
      qm = qm + dum
      bm = bm + (qrcol*inv_rho_rimeMax+qccol/rho_qm_cloud+(qr2qi_immers_freeze_tend+ &
                 qc2qi_hetero_freeze_tend)*inv_rho_rimeMax)*dt
      ni = ni + (ni_nucleat_tend-ni2nr_melt_tend-ni_sublim_tend-ni_selfcollect_tend+ &
                 nr2ni_immers_freeze_tend+nc2ni_immers_freeze_tend)*dt
   else
      dum = (qrcol+qccol+qr2qi_immers_freeze_tend+qcheti_cnt+qicnt)*dt
      qi = qi + (qidep+qinuc+qiberg+qinuc_cnt)*dt + dum
      qi_wetDepos = (qidep+qinuc+qiberg+qinuc_cnt) + (qrcol+qccol+qr2qi_immers_freeze_tend+qcheti_cnt+qicnt) ! total ice process rate for wet-deopsition 
      qm = qm + dum
      bm = bm + (qrcol*inv_rho_rimeMax+qccol/rho_qm_cloud+(qr2qi_immers_freeze_tend+ &
                 qcheti_cnt+qicnt)*inv_rho_rimeMax)*dt
      ni = ni + (ni_nucleat_tend-ni2nr_melt_tend-ni_sublim_tend-ni_selfcollect_tend+ &
                nr2ni_immers_freeze_tend+ncheti_cnt+nicnt+ninuc_cnt)*dt
   endif     

   !PMC nCat deleted interactions_loop

   if (qm.lt.0._rtype) then
      qm = 0._rtype
      bm = 0._rtype
   endif

   ! densify under wet growth
   ! -- to be removed post-v2.1.  Densification automatically happens
   !    during wet growth due to parameterized rime density --
   if (log_wetgrowth) then
      qm = qi
      bm = qm*inv_rho_rimeMax
   endif

   ! densify in above freezing conditions and melting
   ! -- future work --
   !   Ideally, this will be treated with the predicted liquid fraction in ice.
   !   Alternatively, it can be simplified by tending qm -- qi
   !   and bm such that rho_rim (qm/bm) --> rho_liq during melting.
   ! ==

   if(.not. use_hetfrz_classnuc)then
      qv = qv + (-qidep+qi2qv_sublim_tend-qinuc)*dt
      th_atm = th_atm + exner*((qidep-qi2qv_sublim_tend+qinuc)*latent_heat_sublim*inv_cp +(qrcol+qccol+   &
                        qc2qi_hetero_freeze_tend+qr2qi_immers_freeze_tend-qi2qr_melt_tend+qiberg)* latent_heat_fusion*inv_cp)*dt
   else
      qv = qv + (-qidep+qi2qv_sublim_tend-qinuc-qinuc_cnt)*dt
      th_atm = th_atm + exner*((qidep-qi2qv_sublim_tend+qinuc+qinuc_cnt)*latent_heat_sublim*inv_cp +(qrcol+qccol+   &
                        qcheti_cnt+qicnt+qr2qi_immers_freeze_tend-qi2qr_melt_tend+qiberg)* latent_heat_fusion*inv_cp)*dt
   endif                     

end subroutine update_prognostic_ice

subroutine update_prognostic_liquid(qc2qr_accret_tend,nc_accret_tend,qc2qr_autoconv_tend,nc2nr_autoconv_tend, &
     ncautr,nc_selfcollect_tend, qr2qv_evap_tend,nr_evap_tend,nr_selfcollect_tend,         &
    do_predict_nc, nccnst, do_prescribed_CCN, inv_rho,exner,latent_heat_vapor,dt,          &
    th_atm,qv,qc,nc,qr,nr)

   !-- warm-phase only processes:
   implicit none

   real(rtype), intent(in) :: qc2qr_accret_tend
   real(rtype), intent(in) :: nc_accret_tend
   real(rtype), intent(in) :: qc2qr_autoconv_tend
   real(rtype), intent(in) :: nc2nr_autoconv_tend
   real(rtype), intent(in) :: ncautr
   real(rtype), intent(in) :: nc_selfcollect_tend
   real(rtype), intent(in) :: qr2qv_evap_tend
   real(rtype), intent(in) :: nr_evap_tend
   real(rtype), intent(in) :: nr_selfcollect_tend


   logical(btype), intent(in) :: do_predict_nc, do_prescribed_CCN
   real(rtype), intent(in) :: nccnst
   real(rtype), intent(in) :: inv_rho
   real(rtype), intent(in) :: exner
   real(rtype), intent(in) :: latent_heat_vapor
   real(rtype), intent(in) :: dt

   real(rtype), intent(inout) :: th_atm
   real(rtype), intent(inout) :: qv
   real(rtype), intent(inout) :: qc
   real(rtype), intent(inout) :: nc
   real(rtype), intent(inout) :: qr
   real(rtype), intent(inout) :: nr

   qc = qc + (-qc2qr_accret_tend-qc2qr_autoconv_tend)*dt
   qr = qr + (qc2qr_accret_tend+qc2qr_autoconv_tend-qr2qv_evap_tend)*dt

   if (do_predict_nc .or. do_prescribed_CCN) then
      nc = nc + (-nc_accret_tend-nc2nr_autoconv_tend+nc_selfcollect_tend)*dt
   else
      nc = nccnst*inv_rho
   endif
   if (iparam.eq.1 .or. iparam.eq.2) then
      nr = nr + (0.5_rtype*nc2nr_autoconv_tend-nr_selfcollect_tend-nr_evap_tend)*dt
   else
      nr = nr + (ncautr-nr_selfcollect_tend-nr_evap_tend)*dt
   endif

   qv = qv + qr2qv_evap_tend*dt
   th_atm = th_atm + exner*(-qr2qv_evap_tend*latent_heat_vapor*    &
        inv_cp)*dt

end subroutine update_prognostic_liquid



subroutine ice_deposition_sublimation(qi_incld,ni_incld,t_atm,    &
qv_sat_l,qv_sat_i,epsi,abi,qv,p3_wbf_coeff,    &
qidep,qi2qv_sublim_tend,ni_sublim_tend,qiberg)

   implicit none

   real(rtype), intent(in)  :: qi_incld
   real(rtype), intent(in)  :: ni_incld
   real(rtype), intent(in)  :: t_atm
   real(rtype), intent(in)  :: qv_sat_l
   real(rtype), intent(in)  :: qv_sat_i
   real(rtype), intent(in)  :: epsi
   real(rtype), intent(in)  :: abi
   real(rtype), intent(in)  :: qv
   real(rtype), intent(in)  :: p3_wbf_coeff
   real(rtype), intent(out) :: qidep
   real(rtype), intent(out) :: qi2qv_sublim_tend
   real(rtype), intent(out) :: ni_sublim_tend
   real(rtype), intent(out) :: qiberg

   real(rtype) :: oabi

   oabi = 1._rtype/abi
   if (qi_incld>=qsmall) then
      !Compute deposition/sublimation
      qidep = epsi * oabi * (qv - qv_sat_i)
      !Split into deposition or sublimation.
      if (t_atm < T_zerodegc .and. qidep>0._rtype) then
         qi2qv_sublim_tend=0._rtype
      else
      ! make qi2qv_sublim_tend positive for consistency with other evap/sub processes
         qi2qv_sublim_tend=-min(qidep,0._rtype)
         qidep=0._rtype
      end if
      !sublimation occurs @ any T. Not so for berg.
      if (t_atm < T_zerodegc) then
         !Compute bergeron rate assuming cloud for whole step.
         qiberg = max(epsi*oabi*(qv_sat_l - qv_sat_i), 0._rtype)
         qiberg = qiberg*p3_wbf_coeff
      else !T>frz
         qiberg=0._rtype
      end if !T<frz
      ni_sublim_tend = qi2qv_sublim_tend*(ni_incld/qi_incld)
   else
      qiberg = 0._rtype
      qidep  = 0._rtype
      qi2qv_sublim_tend  = 0._rtype
      ni_sublim_tend  = 0._rtype
   end if

   return

end subroutine ice_deposition_sublimation

subroutine rain_evap_tscale_weight(dt_over_tau,weight)
  !Returns weighting between 0 and 1 for how much of the instantaneous
  !evaporation rate and how much of the equilibrium evaporation rate to
  !blend to get the timestep-average rain evaporation rate

  real(rtype), intent(in) :: dt_over_tau  !microphysics timestep divided by effective evap timescale
  real(rtype), intent(out) :: weight

  !expm1 is exp(x)-1. This impl is more accurate than exp near x=0.
  weight= -bfb_expm1(-dt_over_tau)/dt_over_tau

  return
end subroutine rain_evap_tscale_weight

subroutine rain_evap_equilib_tend(A_c,ab,tau_eff,tau_r, tend)
  !In equilibrium, the total evaporation must balance the tendency A_c from
  !all other processes. The rain evaporation is the fraction (1/tau_r)/(1/tau_eff)
  !of the total tendency and ab corrects for saturation changes due to evaporative
  !cooling.

  real(rtype), intent(in) :: A_c, ab, tau_eff, tau_r
  real(rtype), intent(out) :: tend

  !sign convention: Negative A_c causes a supersaturation deficit which needs to be removed
  !by evaporation (which is signed positive when active) to maintain equilibrium. Thus need
  !a negative sign here since other terms are always positive.
  tend = -A_c/ab*tau_eff/tau_r

  return
end subroutine rain_evap_equilib_tend

subroutine rain_evap_instant_tend(ssat_r, ab, tau_r,tend)
  !The instantaneous rain evap tendency is just the absolute supersaturation
  !ssat_r divided by the supersaturation removal timescale for rain tau_r
  !corrected for the effect of evaporative cooling on saturation ab.

  real(rtype), intent(in) :: ssat_r, ab, tau_r
  real(rtype), intent(out) :: tend

  !sign convention: ssat_r must be <0 for evap, other terms are always positive,
  !and we want evap rate positive... so put a minus sign in front.

  tend = -ssat_r/(ab*tau_r)

  return
end subroutine rain_evap_instant_tend


subroutine evaporate_rain(qr_incld,qc_incld,nr_incld,qi_incld, &
cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i, &
ab,abi,epsr,epsi_tot,t,t_prev,latent_heat_sublim,dqsdt,dt, &
qr2qv_evap_tend,nr_evap_tend)

  !Evaporation is basically (qv - sv_sat)/(tau_eff*ab) where tau_eff
  !is the total effective supersaturation removal timescale
  !and ab is the psychrometric correction for condensational heating
  !changing qv_sat. This formulation depends sensitively on ssat_r, which
  !can change rapidly within a timestep because liquid saturation
  !adjustment has a relaxation timescale of seconds. For accuracy and
  !stability, we analytically integrate ssat_r over the timestep under
  !the simplifying assumption that all processes other than saturation
  !relaxation are a constant source/sink term A_c. See Morrison+Milbrandt 2015
  !https://doi.org/10.1175/JAS-D-14-0065.1 and Morrison+Grabowski 2008
  !https://doi.org/10.1175/2007JAS2374.1 for details.

   implicit none

   real(rtype), intent(in)  :: qr_incld
   real(rtype), intent(in)  :: qc_incld
   real(rtype), intent(in)  :: nr_incld
   real(rtype), intent(in)  :: qi_incld
   real(rtype), intent(in)  :: cld_frac_l
   real(rtype), intent(in)  :: cld_frac_r
   real(rtype), intent(in)  :: qv_sat_l,qv_sat_i
   real(rtype), intent(in)  :: ab,abi
   real(rtype), intent(in)  :: epsr,epsi_tot
   real(rtype), intent(in)  :: qv,qv_prev
   real(rtype), intent(in)  :: t,t_prev,latent_heat_sublim,dqsdt,dt
   real(rtype), intent(out) :: qr2qv_evap_tend
   real(rtype), intent(out) :: nr_evap_tend
   real(rtype) :: cld_frac, eps_eff, tau_eff, tau_r, ssat_r, A_c, sup_r,inv_dt
   real(rtype) :: equilib_evap_tend, tscale_weight, instant_evap_tend

   !Initialize variables
   qr2qv_evap_tend = 0.0_rtype
   nr_evap_tend = 0.0_rtype
   inv_dt=1._rtype/dt

   !Compute absolute supersaturation.
   !Ignore the difference between clear-sky and cell-ave qv and T
   !because micro lacks the info to reliably reconstruct macrophys
   !subgrid variability
   ssat_r = qv - qv_sat_l !absolute supersaturation

   !Cloud fraction in clear-sky conditions has been set to mincld
   !to avoid divide-by-zero problems. Because rain evap only happens
   !in rainy portions outside cloud, setting clear-sky cloud fraction
   !to mincld reduces evaporating area. We fix that here by computing
   !a temporary cloud fraction which is zero if cloud condensate is small.
   if (qc_incld + qi_incld < 1.e-6_rtype) then
         cld_frac = 0._rtype
   else
         cld_frac = cld_frac_l
   end if

   !Only evaporate in the rainy area outside cloud when subsaturated
   !Note: ignoring case where cell initially supersaturated but other processes would make
   !it subsaturated within 1 timestep.
   if (cld_frac_r > cld_frac .and. ssat_r<0._rtype .and. qr_incld >= qsmall ) then

      !Compute total effective inverse saturation removal timescale eps_eff
      !qc saturation is handled by macrophysics so the qc saturation removal timescale is
      !not included here. Below freezing, eps_eff is the sum of the inverse saturation
      !removal timescales for liquid and ice. The ice term has extra scaling terms to convert
      !it from being relative to ice to liquid instead. See Eq C3 of Morrison+Milbrandt 2015
      !https://doi.org/10.1175/JAS-D-14-0065.1
      if (t < 273.15_rtype) then
         eps_eff   = epsr + epsi_tot*(1.0_rtype + latent_heat_sublim*inv_cp*dqsdt)/abi
      else
         eps_eff   = epsr
      endif

      !Set lower bound on eps_eff to prevent division by zero
      eps_eff  = max(1.e-20_rtype,eps_eff)
      tau_eff = 1.0_rtype/eps_eff

      !Compute the constant source/sink term A_c for analytic integration. See Eq C4 in
      !Morrison+Milbrandt 2015 https://doi.org/10.1175/JAS-D-14-0065.1
      if (t < 273.15_rtype) then
         A_c = (qv - qv_prev)*inv_dt - dqsdt*(t-t_prev)*inv_dt - (qv_sat_l - qv_sat_i)* &
               (1.0_rtype + latent_heat_sublim*inv_cp*dqsdt)/abi*epsi_tot
      else
         A_c = (qv - qv_prev)*inv_dt - dqsdt*(t-t_prev)*inv_dt
      endif

      !Now compute evap rate

      !note that qr_incld<qsmall => qr2qv_evap_tend already initialized to zero above.

      !If there's negligible qr and conditions are subsaturated, evaporate all qr
      if (qr_incld < 1e-12_rtype .and. qv/qv_sat_l < 0.999_rtype) then
         qr2qv_evap_tend = qr_incld*inv_dt

      !If sizable qr, compute tend.
      else
         !Timestep-averaged evap can be written as the weighted average of instantaneous
         !and equilibrium evap rates with weighting tscale_weight. L'Hospital's rule
         !shows tscale_weight is 1 in the limit of small dt. It approaches 0 as dt
         !gets big.
         call rain_evap_tscale_weight(dt/tau_eff,tscale_weight)

         !tau_r is only used in this "else" branch so only define it here.
         !outside this branch qr_incld could be < qsmall, in which case epsr=0.
         tau_r = 1._rtype/epsr

         !in limit of very long timescales, evap must balance A_c.
         !(1/tau_r)/(1/tau_eff) is the fraction of this total tendency assigned to rain
         !Will be >0 if A_c>0: increased supersat from other procs must be balanced by
         !evaporation to stay in equilibrium.
         call  rain_evap_equilib_tend(A_c,ab,tau_eff,tau_r,equilib_evap_tend)

         !in limit of short timesteps, evap can be calculated from ssat_r at the
         !beginning of the timestep
         !ssat_r<0 when evap occurs and evap_tend is positive when evaporating, so added
         !neg in front
         call  rain_evap_instant_tend(ssat_r, ab, tau_r, instant_evap_tend)

         qr2qv_evap_tend = instant_evap_tend*tscale_weight &
              + equilib_evap_tend*(1._rtype-tscale_weight)

      end if


      !Limit evap from exceeding saturation deficit. Analytic integration
      !would prevent this from happening if A_c was part of microphysics
      !timestepping, but it isn't.
      qr2qv_evap_tend = min(qr2qv_evap_tend,-ssat_r*inv_dt/ab)

      !To maintain equilibrium, the equilibrium evaporation tendency must be
      !negative (adding mass) if A_c (other processes) are losing mass. We don't
      !allow rain evap to also condense by forcing qr2qv_evap_tend to be positive
      qr2qv_evap_tend = max(0._rtype, qr2qv_evap_tend)

      !We can't evaporate more rain mass than we had to start with
      !Sanity check: We're applying to rainy region outside cloud here because
      !qr inside cloud should be protected from evap. Conversion to rainy-area
      !ave just below scales by (cldfrac_r - cld_frac)/cldfrac_r < 1 so
      !total qr isn't pushed negative.
      qr2qv_evap_tend = min(qr2qv_evap_tend,qr_incld*inv_dt)

      !Evap rate so far is an average over the rainy area outside clouds.
      !Turn this into an average over the entire raining area
      qr2qv_evap_tend = qr2qv_evap_tend*(cld_frac_r-cld_frac)/cld_frac_r

      !Let nr remove drops proportionally to mass change
      nr_evap_tend = qr2qv_evap_tend*(nr_incld/qr_incld)


   end if !cld_frac_r>cldfrac and ssat_r<0

   return

end subroutine evaporate_rain

subroutine get_time_space_phys_variables( &
t_atm,pres,rho,latent_heat_vapor,latent_heat_sublim,qv_sat_l,qv_sat_i, &
mu,dv,sc,dqsdt,dqsidt,ab,abi,kap,eii)

   implicit none

   real(rtype), intent(in)  :: t_atm
   real(rtype), intent(in)  :: pres
   real(rtype), intent(in)  :: rho
   real(rtype), intent(in)  :: latent_heat_vapor
   real(rtype), intent(in)  :: latent_heat_sublim
   real(rtype), intent(in)  :: qv_sat_l
   real(rtype), intent(in)  :: qv_sat_i
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
   mu     = 1.496e-6_rtype*bfb_pow(t_atm,1.5_rtype)/(t_atm+120._rtype)
   dv     = 8.794e-5_rtype*bfb_pow(t_atm,1.81_rtype)/pres
   sc     = mu/(rho*dv)
   dum    = 1._rtype/(rv*bfb_square(t_atm))
   dqsdt  = latent_heat_vapor*qv_sat_l*dum
   dqsidt = latent_heat_sublim*qv_sat_i*dum
   ab     = 1._rtype+dqsdt*latent_heat_vapor*inv_cp
   abi    = 1._rtype+dqsidt*latent_heat_sublim*inv_cp
   kap    = 1.414e+3_rtype*mu

   ! very simple temperature dependent aggregation efficiency
   if (t_atm.lt.253.15_rtype) then
      eii = 0.001_rtype
   else if (t_atm.ge.253.15_rtype.and.t_atm.lt.273.15_rtype) then
      eii = 0.001_rtype + (t_atm-253.15_rtype)*(0.3_rtype - 0.001_rtype)/20._rtype  ! linear ramp from 0.1 to 1 between 253.15 and 268.15 K
   else
      eii = 0.3_rtype
   end if

   return

end subroutine get_time_space_phys_variables

subroutine cloud_sedimentation(kts,kte,ktop,kbot,kdir,   &
   qc_incld,rho,inv_rho,cld_frac_l,acn,inv_dz,&
   dt,inv_dt,dnu,do_predict_nc, &
   qc, nc, nc_incld,mu_c,lamc,precip_liq_surf,cflx,qc_tend,nc_tend)

   implicit none
   integer, intent(in) :: kts, kte
   integer, intent(in) :: ktop, kbot, kdir

   real(rtype), intent(in), dimension(kts:kte) :: rho
   real(rtype), intent(in), dimension(kts:kte) :: inv_rho
   real(rtype), intent(in), dimension(kts:kte) :: cld_frac_l
   real(rtype), intent(in), dimension(kts:kte) :: acn
   real(rtype), intent(in), dimension(kts:kte) :: inv_dz
   real(rtype), intent(in) :: dt
   real(rtype), intent(in) :: inv_dt
   real(rtype), dimension(:), intent(in) :: dnu
   logical(btype), intent(in) :: do_predict_nc

   real(rtype), intent(inout), dimension(kts:kte), target :: qc
   real(rtype), intent(inout), dimension(kts:kte), target :: nc
   real(rtype), intent(inout), dimension(kts:kte) :: qc_incld
   real(rtype), intent(inout), dimension(kts:kte) :: nc_incld
   real(rtype), intent(inout), dimension(kts:kte) :: mu_c
   real(rtype), intent(inout), dimension(kts:kte) :: lamc
   real(rtype), intent(inout) :: precip_liq_surf
   real(rtype), intent(inout), dimension(kts:kte+1) :: cflx
   real(rtype), intent(inout), dimension(kts:kte) :: qc_tend
   real(rtype), intent(inout), dimension(kts:kte) :: nc_tend

   logical(btype) :: log_qxpresent
   integer :: k
   integer :: k_qxtop, k_qxbot
   integer, parameter :: num_arrays = 2
   type(realptr), dimension(num_arrays) :: vs, fluxes, qnr

   real(rtype) :: dt_left
   real(rtype) :: dt_sub
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

      two_moment: if (do_predict_nc) then  !2-moment cloud:
         substep_sedi_c2: do while (dt_left.gt.1.e-4_rtype)

            Co_max = 0._rtype
            V_qc = 0._rtype
            V_nc = 0._rtype
            dt_sub = 0._rtype
            kloop_sedi_c2: do k = k_qxtop,k_qxbot,-kdir

               qc_notsmall_c2: if (qc_incld(k)>qsmall) then
                  !-- compute Vq, Vn
                  call get_cloud_dsd2(qc_incld(k),nc_incld(k),mu_c(k),rho(k),nu,dnu,   &
                       lamc(k),tmp1,tmp2)

                  !get_cloud_dsd2 keeps the drop-size distribution within reasonable
                  !bounds by modifying nc_incld. The next line maintains consistency
                  !between nc_incld and nc
                  nc(k) = nc_incld(k)*cld_frac_l(k)

                  dum = 1._rtype / bfb_pow(lamc(k), bcn)
                  V_qc(k) = acn(k)*bfb_gamma(4._rtype+bcn+mu_c(k))*dum/(bfb_gamma(mu_c(k)+4._rtype))
                  V_nc(k) = acn(k)*bfb_gamma(1._rtype+bcn+mu_c(k))*dum/(bfb_gamma(mu_c(k)+1._rtype))

               endif qc_notsmall_c2
               Co_max = max(Co_max, V_qc(k)*dt_left*inv_dz(k))

            enddo kloop_sedi_c2

            call generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, &
                 prt_accum, inv_dz, inv_rho, rho, num_arrays, vs, fluxes, qnr, dt_sub)
            do k = k_qxbot,k_qxtop,kdir
                  cflx(k+1) = cflx(k+1) + flux_qx(k)*dt_sub
            enddo
                 
            !Update _incld values with end-of-step cell-ave values
            !Note that cld_frac_l is set in interface to have min of mincld=1e-4
            !so dividing by it is fine.
            qc_incld(:) = qc(:)/cld_frac_l(:)
            nc_incld(:) = nc(:)/cld_frac_l(:)

         enddo substep_sedi_c2
      else
         substep_sedi_c1: do while (dt_left.gt.1.e-4_rtype)

            Co_max  = 0._rtype
            V_qc = 0._rtype

            kloop_sedi_c1: do k = k_qxtop,k_qxbot,-kdir
               qc_notsmall_c1: if (qc_incld(k)>qsmall) then
                  call get_cloud_dsd2(qc_incld(k),nc_incld(k),mu_c(k),rho(k),nu,dnu,   &
                       lamc(k),tmp1,tmp2)

                  !get_cloud_dsd2 keeps the drop-size distribution within reasonable
                  !bounds by modifying nc_incld. The next line maintains consistency
                  !between nc_incld and nc
                  nc(k) = nc_incld(k)*cld_frac_l(k)

                  dum = 1._rtype / bfb_pow(lamc(k), bcn)
                  V_qc(k) = acn(k)*bfb_gamma(4._rtype+bcn+mu_c(k))*dum/(bfb_gamma(mu_c(k)+4._rtype))
               endif qc_notsmall_c1

               Co_max = max(Co_max, V_qc(k)*dt_left*inv_dz(k))
            enddo kloop_sedi_c1

            call generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, &
                 prt_accum, inv_dz, inv_rho, rho, 1, vs, fluxes, qnr, dt_sub)
            do k = k_qxbot,k_qxtop,kdir
                  cflx(k+1) = cflx(k+1) + flux_qx(k)*dt_sub
            enddo

            !Update _incld values with end-of-step cell-ave values
            !Note that cld_frac_l is set in interface to have min of mincld=1e-4
            !so dividing by it is fine.
            qc_incld(:) = qc(:)/cld_frac_l(:)
            nc_incld(:) = nc(:)/cld_frac_l(:)

         enddo substep_sedi_c1

      endif two_moment
      cflx(:) = cflx(:)*inv_dt
      precip_liq_surf = precip_liq_surf + prt_accum*inv_rho_h2o*inv_dt

   endif qc_present

   qc_tend(:) = ( qc(:) - qc_tend(:) ) * inv_dt ! Liq. sedimentation tendency, measure
   nc_tend(:) = ( nc(:) - nc_tend(:) ) * inv_dt ! Liq. # sedimentation tendency, measure

end subroutine cloud_sedimentation

subroutine rain_sedimentation(kts,kte,ktop,kbot,kdir,                             &
   qr_incld,rho,inv_rho,rhofacr,cld_frac_r,inv_dz,dt,inv_dt,p3_max_mean_rain_size,&
   qr,nr,nr_incld,mu_r,lamr,precip_liq_surf,precip_liq_flux,rflx,qr_tend,nr_tend)

   implicit none
   integer, intent(in) :: kts, kte
   integer, intent(in) :: ktop, kbot, kdir

   real(rtype), intent(in), dimension(kts:kte) :: rho
   real(rtype), intent(in), dimension(kts:kte) :: inv_rho
   real(rtype), intent(in), dimension(kts:kte) :: rhofacr
   real(rtype), intent(in), dimension(kts:kte) :: cld_frac_r
   real(rtype), intent(in), dimension(kts:kte) :: inv_dz
   real(rtype), intent(in) :: dt
   real(rtype), intent(in) :: inv_dt
   real(rtype), intent(in) :: p3_max_mean_rain_size

   real(rtype), intent(inout), target, dimension(kts:kte) :: qr
   real(rtype), intent(inout), target, dimension(kts:kte) :: nr
   real(rtype), intent(inout), dimension(kts:kte) :: qr_incld
   real(rtype), intent(inout), dimension(kts:kte) :: nr_incld
   real(rtype), intent(inout), dimension(kts:kte) :: mu_r
   real(rtype), intent(inout), dimension(kts:kte) :: lamr
   real(rtype), intent(inout) :: precip_liq_surf
   real(rtype), intent(inout), dimension(kts:kte+1) :: precip_liq_flux
   real(rtype), intent(inout), dimension(kts:kte+1) :: rflx
   real(rtype), intent(inout), dimension(kts:kte) :: qr_tend
   real(rtype), intent(inout), dimension(kts:kte) :: nr_tend

   logical(btype) :: log_qxpresent
   integer :: k
   integer :: k_qxtop, k_qxbot
   integer, parameter :: num_arrays = 2
   type(realptr), dimension(num_arrays) :: vs, fluxes, qnr

   real(rtype) :: dt_left
   real(rtype) :: dt_sub
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
         dt_sub = 0._rtype

         kloop_sedi_r1: do k = k_qxtop,k_qxbot,-kdir

            qr_notsmall_r1: if (qr_incld(k)>qsmall) then

               call compute_rain_fall_velocity(qr_incld(k), rhofacr(k), p3_max_mean_rain_size, nr_incld(k), &
                    mu_r(k), lamr(k), V_qr(k), V_nr(k))

               !in compute_rain_fall_velocity, get_rain_dsd2 keeps the drop-size
               !distribution within reasonable bounds by modifying nr_incld.
               !The next line maintains consistency between nr_incld and nr.
               nr(k) = nr_incld(k)*cld_frac_r(k)

            endif qr_notsmall_r1

            Co_max = max(Co_max, V_qr(k)*dt_left*inv_dz(k))
            !            Co_max = max(Co_max, max(V_nr(k),V_qr(k))*dt_left*inv_dz(i,k))

         enddo kloop_sedi_r1

         call generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, &
              prt_accum, inv_dz, inv_rho, rho, num_arrays, vs, fluxes, qnr, dt_sub)

         !-- AaronDonahue, precip_liq_flux output
         do k = k_qxbot,k_qxtop,kdir
            precip_liq_flux(k+1) = precip_liq_flux(k+1) + flux_qx(k)*dt_sub ! AaronDonahue
            rflx(k+1) = rflx(k+1) + flux_qx(k)*dt_sub
         enddo

         !Update _incld values with end-of-step cell-ave values
         !Note that cld_frac_r is set in interface to have min of mincld=1e-4
         !so dividing by it is fine.
         qr_incld(:) = qr(:)/cld_frac_r(:)
         nr_incld(:) = nr(:)/cld_frac_r(:)

      enddo substep_sedi_r
      precip_liq_flux(:) = precip_liq_flux(:)*inv_dt
      rflx(:) = rflx(:)*inv_dt
      precip_liq_surf = precip_liq_surf + prt_accum*inv_rho_h2o*inv_dt

   endif qr_present

   qr_tend(:) = ( qr(:) - qr_tend(:) ) * inv_dt ! Rain sedimentation tendency, measure
   nr_tend(:) = ( nr(:) - nr_tend(:) ) * inv_dt ! Rain # sedimentation tendency, measure

end subroutine rain_sedimentation

subroutine compute_rain_fall_velocity(qr_incld, rhofacr, p3_max_mean_rain_size, nr_incld, mu_r, lamr, V_qr, V_nr)

   real(rtype), intent(in) :: qr_incld
   real(rtype), intent(in) :: rhofacr
   real(rtype), intent(in) :: p3_max_mean_rain_size
   real(rtype), intent(inout) :: nr_incld
   real(rtype), intent(out) :: mu_r
   real(rtype), intent(out) :: lamr
   real(rtype), intent(out) :: V_qr
   real(rtype), intent(out) :: V_nr

   real(rtype) :: tmp1, tmp2, dum1, dum2, inv_dum3, rdumii, rdumjj
   integer :: dumii, dumjj

   !Compute Vq, Vn:

   call get_rain_dsd2(qr_incld,p3_max_mean_rain_size,nr_incld,mu_r,lamr,tmp1,tmp2)

   call find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3,mu_r,lamr)

   !mass-weighted fall speed:

   dum1 = vm_table_vals(dumii,dumjj)+(rdumii-real(dumii))*                       &
      (vm_table_vals(dumii+1,dumjj)-vm_table_vals(dumii,dumjj))       !at mu_r
   dum2 = vm_table_vals(dumii,dumjj+1)+(rdumii-real(dumii))*                     &
      (vm_table_vals(dumii+1,dumjj+1)-vm_table_vals(dumii,dumjj+1))   !at mu_r+1

   V_qr = dum1 + (rdumjj-real(dumjj))*(dum2-dum1)         !interpolated
   V_qr = V_qr*rhofacr                                    !corrected for air density

   ! number-weighted fall speed:
   dum1 = vn_table_vals(dumii,dumjj)+(rdumii-real(dumii))*                       &
      (vn_table_vals(dumii+1,dumjj)-vn_table_vals(dumii,dumjj))       !at mu_r
   dum2 = vn_table_vals(dumii,dumjj+1)+(rdumii-real(dumii))*                     &
      (vn_table_vals(dumii+1,dumjj+1)-vn_table_vals(dumii,dumjj+1))   !at mu_r+1

   V_nr = dum1+(rdumjj-real(dumjj))*(dum2-dum1)            !interpolated
   V_nr = V_nr*rhofacr               !corrected for air density
end subroutine compute_rain_fall_velocity

subroutine ice_sedimentation(kts,kte,ktop,kbot,kdir,    &
   rho,inv_rho,rhofaci,cld_frac_i,inv_dz,dt,inv_dt,  &
   qi,qi_incld,ni,qm,qm_incld,bm,bm_incld,ni_incld,precip_ice_surf,precip_ice_flux,sflx,qi_tend,ni_tend)

   implicit none
   integer, intent(in) :: kts, kte
   integer, intent(in) :: ktop, kbot, kdir

   real(rtype), intent(in), dimension(kts:kte) :: rho
   real(rtype), intent(in), dimension(kts:kte) :: inv_rho
   real(rtype), intent(in), dimension(kts:kte) :: rhofaci
   real(rtype), intent(in), dimension(kts:kte) :: cld_frac_i
   real(rtype), intent(in), dimension(kts:kte) :: inv_dz
   real(rtype), intent(in) :: dt
   real(rtype), intent(in) :: inv_dt

   real(rtype), intent(inout), dimension(kts:kte), target :: qi
   real(rtype), intent(inout), dimension(kts:kte) :: qi_incld
   real(rtype), intent(inout), dimension(kts:kte), target :: ni
   real(rtype), intent(inout), dimension(kts:kte) :: ni_incld
   real(rtype), intent(inout), dimension(kts:kte), target :: qm
   real(rtype), intent(inout), dimension(kts:kte) :: qm_incld
   real(rtype), intent(inout), dimension(kts:kte), target :: bm
   real(rtype), intent(inout), dimension(kts:kte) :: bm_incld

   real(rtype), intent(inout) :: precip_ice_surf
   real(rtype), intent(inout), dimension(kts:kte+1) :: precip_ice_flux
   real(rtype), intent(inout), dimension(kts:kte+1) :: sflx
   real(rtype), intent(inout), dimension(kts:kte) :: qi_tend
   real(rtype), intent(inout), dimension(kts:kte) :: ni_tend

   real(rtype) :: dt_sub
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
   real(rtype) :: table_val_ni_fallspd ! number-weighted fallspeed            See lines  731 -  808  uns
   real(rtype) :: table_val_qi_fallspd ! mass-weighted fallspeed              See lines  731 -  808  ums
   real(rtype) :: table_val_ni_lammax ! minimum ice number (lambda limiter)  See lines  704 -  705  nlarge
   real(rtype) :: table_val_ni_lammin ! maximum ice number (lambda limiter)  See lines  704 -  705  nsmall

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
   qnr(1)%p => qi
   qnr(2)%p => ni
   qnr(3)%p => qm
   qnr(4)%p => bm

   !find top, determine qxpresent
   do k = ktop,kbot,-kdir
      if (qi(k).ge.qsmall) then
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
         if (qi(k).ge.qsmall) then
            k_qxbot = k
            exit
         endif
      enddo

      substep_sedi_i: do while (dt_left.gt.1.e-4_rtype)

         Co_max = 0._rtype
         V_qit = 0._rtype
         V_nit = 0._rtype
         dt_sub = 0._rtype

         kloop_sedi_i1: do k = k_qxtop,k_qxbot,-kdir

            !-- compute Vq, Vn (get values from lookup table)
            qi_notsmall_i1: if (qi_incld(k)>qsmall) then

               !--Compute Vq, Vn:
               ni_incld(k) = max(ni_incld(k),nsmall) !impose lower limits to prevent log(<0)
               call calc_bulkRhoRime(qi_incld(k),qm_incld(k),bm_incld(k),rhop)
               qm(k)=qm_incld(k)*cld_frac_i(k)
               bm(k)=bm_incld(k)*cld_frac_i(k)

               !if (.not. tripleMoment_on) zitot(i,k) = diag_mom6(qi(i,k),ni(i,k),rho(i,k))
               call find_lookupTable_indices_1a(dumi,dumjj,dumii,dumzz,dum1,dum4,      &
                    dum5,dum6,isize,rimsize,densize,                                   &
                    qi_incld(k),ni_incld(k),qm_incld(k),                               &
                    rhop)
               call access_lookup_table(dumjj,dumii,dumi, 1,dum1,dum4,dum5,table_val_ni_fallspd)
               call access_lookup_table(dumjj,dumii,dumi, 2,dum1,dum4,dum5,table_val_qi_fallspd)
               call access_lookup_table(dumjj,dumii,dumi, 7,dum1,dum4,dum5,table_val_ni_lammax)
               call access_lookup_table(dumjj,dumii,dumi, 8,dum1,dum4,dum5,table_val_ni_lammin)
               !-impose mean ice size bounds (i.e. apply lambda limiters)
               ! note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
               ni_incld(k) = min(ni_incld(k),table_val_ni_lammax*ni_incld(k))
               ni_incld(k) = max(ni_incld(k),table_val_ni_lammin*ni_incld(k))
               ni(k) = ni_incld(k)*cld_frac_i(k)
               !zitot(i,k) = min(zitot(i,k),table_val_qi_fallspd0)  !adjust Zi if needed to make sure mu_i is in bounds
               !zitot(i,k) = max(zitot(i,k),table_val_qi_fallspd1)
               V_qit(k) = table_val_qi_fallspd*rhofaci(k)     !mass-weighted  fall speed (with density factor)
               V_nit(k) = table_val_ni_fallspd*rhofaci(k)     !number-weighted    fall speed (with density factor)
               !==

            endif qi_notsmall_i1

            Co_max = max(Co_max, V_qit(k)*dt_left*inv_dz(k))

         enddo kloop_sedi_i1

         call generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, &
              dt_left, prt_accum, inv_dz, inv_rho, rho, num_arrays, vs, fluxes, qnr, dt_sub)

         do k = k_qxbot,k_qxtop,kdir
            sflx(k+1) = sflx(k+1) + flux_qit(k)*dt_sub
         enddo     

         !update _incld variables
         !Note that cld_frac_i is set in interface to have min of mincld=1e-4
         !so dividing by it is fine.
         qi_incld(:) = qi(:)/cld_frac_i(:)
         ni_incld(:) = ni(:)/cld_frac_i(:)
         qm_incld(:) = qm(:)/cld_frac_i(:)
         bm_incld(:) = bm(:)/cld_frac_i(:)

      enddo substep_sedi_i
      sflx(:)=sflx(:)*inv_dt
      precip_ice_surf = precip_ice_surf + prt_accum*inv_rho_h2o*inv_dt

   endif qi_present

   qi_tend(:) = ( qi(:) - qi_tend(:) ) * inv_dt ! Ice sedimentation tendency, measure
   ni_tend(:) = ( ni(:) - ni_tend(:) ) * inv_dt ! Ice # sedimentation tendency, measure

end subroutine ice_sedimentation

subroutine generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, &
     prt_accum, inv_dz, inv_rho, rho, num_arrays, vs, fluxes, qnx, dt_sub)

   implicit none

   integer, intent(in) :: kts, kte, kdir, k_qxtop, kbot, num_arrays
   integer, intent(inout) :: k_qxbot
   real(rtype), intent(in) :: Co_max
   real(rtype), intent(inout) :: dt_left, prt_accum
   real(rtype), dimension(kts:kte), intent(in) :: inv_dz
   real(rtype), dimension(kts:kte), intent(in) :: inv_rho
   real(rtype), dimension(kts:kte), intent(in) :: rho

   type(realptr), intent(in), dimension(num_arrays), target :: vs, fluxes, qnx

   integer :: tmpint1, k_temp, i
   real(rtype),intent(inout) :: dt_sub

   !-- compute dt_sub
   tmpint1 = int(Co_max+1._rtype)    !number of substeps remaining if dt_sub were constant
   dt_sub  = min(dt_left, dt_left/float(tmpint1))

   ! -- Move bottom cell down by 1 if not at ground already
   if (k_qxbot.eq.kbot) then
      k_temp = k_qxbot
   else
      k_temp = k_qxbot-kdir
   endif

   call calc_first_order_upwind_step(kts, kte, kdir, k_temp, k_qxtop, dt_sub, rho, inv_rho, inv_dz, num_arrays, fluxes, vs, qnx)

   !accumulated precip during time step
   if (k_qxbot.eq.kbot) prt_accum = prt_accum + fluxes(1)%p(kbot)*dt_sub

   dt_left = dt_left - dt_sub  !update time remaining for sedimentation
   if (k_qxbot.ne.kbot) k_qxbot = k_qxbot - kdir

end subroutine generalized_sedimentation

subroutine calc_first_order_upwind_step(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dz, num_arrays, fluxes, vs, qnx)

  implicit none

  integer, intent(in) :: kts, kte, kdir, kbot, k_qxtop, num_arrays
  real(rtype), intent(in) :: dt_sub
  real(rtype), dimension(kts:kte), intent(in) :: rho, inv_rho, inv_dz
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
     fluxdiv = -fluxes(i)%p(k) * inv_dz(k)
     !- update prognostic variables
     qnx(i)%p(k) = qnx(i)%p(k) + fluxdiv*dt_sub*inv_rho(k)

     do k = k_qxtop-kdir,kbot,-kdir
        !-- compute flux divergence
        fluxdiv = (fluxes(i)%p(k+kdir) - fluxes(i)%p(k))*inv_dz(k)
        !-- update prognostic variables
        qnx(i)%p(k) = qnx(i)%p(k) + fluxdiv*dt_sub*inv_rho(k)
     end do
  end do

end subroutine calc_first_order_upwind_step

subroutine homogeneous_freezing(kts,kte,ktop,kbot,kdir,t_atm,exner,latent_heat_fusion,    &
   qc,nc,qr,nr,qi,ni,qm,bm,th_atm)

   !.......................................
   ! homogeneous freezing of cloud and rain

   implicit none
   integer, intent(in) :: kts, kte
   integer, intent(in) :: ktop, kbot, kdir
   real(rtype), intent(in), dimension(kts:kte) :: t_atm
   real(rtype), intent(in), dimension(kts:kte) :: exner
   real(rtype), intent(in), dimension(kts:kte) :: latent_heat_fusion

   real(rtype), intent(inout), dimension(kts:kte) :: qc
   real(rtype), intent(inout), dimension(kts:kte) :: nc
   real(rtype), intent(inout), dimension(kts:kte) :: qr
   real(rtype), intent(inout), dimension(kts:kte) :: nr

   real(rtype), intent(inout), dimension(kts:kte) :: qi
   real(rtype), intent(inout), dimension(kts:kte) :: ni
   real(rtype), intent(inout), dimension(kts:kte) :: qm
   real(rtype), intent(inout), dimension(kts:kte) :: bm
   real(rtype), intent(inout), dimension(kts:kte) :: th_atm

   real(rtype) :: Q_nuc
   real(rtype) :: N_nuc
   integer :: k

   k_loop_fz:  do k = kbot,ktop,kdir
      if (qc(k).ge.qsmall .and. (th_atm(k)/exner(k)).lt.T_homogfrz) then
         Q_nuc = qc(k)
         N_nuc = max(nc(k),nsmall)

         qm(k) = qm(k) + Q_nuc
         qi(k) = qi(k) + Q_nuc
         bm(k) = bm(k) + Q_nuc*inv_rho_rimeMax
         ni(k) = ni(k) + N_nuc
         th_atm(k) = th_atm(k) + exner(k)*Q_nuc*latent_heat_fusion(k)*inv_cp
         qc(k) = 0._rtype
         nc(k) = 0._rtype

      endif

      if (qr(k).ge.qsmall .and. (th_atm(k)/exner(k)) .lt.T_homogfrz) then
         Q_nuc = qr(k)
         N_nuc = max(nr(k),nsmall)

         qm(k) = qm(k) + Q_nuc
         qi(k) = qi(k) + Q_nuc
         bm(k) = bm(k) + Q_nuc*inv_rho_rimeMax
         ni(k) = ni(k) + N_nuc
         th_atm(k) = th_atm(k) + exner(k)*Q_nuc*latent_heat_fusion(k)*inv_cp
         qr(k) = 0._rtype
         nr(k) = 0._rtype
      endif

   enddo k_loop_fz

end subroutine homogeneous_freezing

subroutine ice_complete_melting(kts,kte,ktop,kbot,kdir,qi,ni,qm,latent_heat_fusion,exner,th_atm, & 
            qr,nr,qc,nc)

   implicit none
   
   integer, intent(in) :: kts, kte
   integer, intent(in) :: ktop, kbot, kdir
   real(rtype), intent(in), dimension(kts:kte) :: exner,latent_heat_fusion
   real(rtype), intent(inout), dimension(kts:kte) :: qi, ni, qc, nc, qr, nr, qm, th_atm
   
   real(rtype) :: t_snow_melt,del_mass,del_num,rv_tmp,rv,frac_mlt,equiv_mass 
   integer :: k

   ! ... For now we take the uper limit for complete melting to be like in MG2 2c
   t_snow_melt = 273.15 + 2.0_rtype         
   k_loop_mlt:  do k = kbot,ktop,kdir
      if (qi(k).ge.qsmall .and. (th_atm(k)/exner(k)) > t_snow_melt) then
         del_mass = qi(k)
         del_num = max(ni(k),nsmall)
         if(th_atm(k)/exner(k) - del_mass*latent_heat_fusion(k)/cp <  t_snow_melt)then
            equiv_mass = (th_atm(k)/exner(k) - t_snow_melt)*cp/latent_heat_fusion(k)
            frac_mlt = equiv_mass/del_mass
            frac_mlt = max(0.0_rtype,frac_mlt)
            frac_mlt = min(1.0_rtype,frac_mlt)
         else
            frac_mlt = 1.0_rtype
         endif
         rv_tmp = 3.0_rtype*qi(k)/ni(k)/4.0_rtype/pi/900.0_rtype  ! density of pure ice [kg/m3]
         rv = 1.0e6_rtype*bfb_cbrt(rv_tmp)                        ! in [um]
         if((qm(k)/qi(k)) < 0.1_rtype)then
            if(rv < 100.0_rtype)then
               ! ... Ice crystas melt into cld droplets
               qi(k) = max((1.0_rtype - frac_mlt)*del_mass,0.0_rtype)
               qm(k) = max((1.0_rtype - frac_mlt)*del_mass,0.0_rtype)
               ni(k) = max((1.0_rtype - frac_mlt)*del_num,0.0_rtype)
               qc(k) = qc(k) + frac_mlt*del_mass
               nc(k) = nc(k) + frac_mlt*del_num
            else
               ! ... Unrimed snow melt to rain
               qi(k) = max((1.0_rtype - frac_mlt)*del_mass,0.0_rtype)
               qm(k) = max((1.0_rtype - frac_mlt)*del_mass,0.0_rtype)
               ni(k) = max((1.0_rtype - frac_mlt)*del_num,0.0_rtype)
               qr(k) = qr(k) + frac_mlt*del_mass
               nr(k) = nr(k) + frac_mlt*del_num
            endif
         else
            ! ... Medium rimed snow melt to rain
            qi(k) = max((1.0_rtype - frac_mlt)*del_mass,0.0_rtype)
            qm(k) = max((1.0_rtype - frac_mlt)*del_mass,0.0_rtype)
            ni(k) = max((1.0_rtype - frac_mlt)*del_num,0.0_rtype)
            qr(k) = qr(k) + frac_mlt*del_mass
            nr(k) = nr(k) + frac_mlt*del_num
         endif
         th_atm(k) = th_atm(k) - exner(k)*frac_mlt*del_mass*latent_heat_fusion(k)/cp
      endif
   enddo k_loop_mlt
   
return
end subroutine ice_complete_melting   

subroutine CNT_couple (frzimm,frzcnt,frzdep,rho,qc_incld,nc_incld,Iflag, &
                       ncheti_cnt,qcheti_cnt,nicnt,qicnt,ninuc_cnt,qinuc_cnt,odt)                   

implicit none

integer, intent(in) :: Iflag
real(rtype), intent(in) :: frzimm,frzcnt,frzdep
real(rtype), intent(in) :: rho,qc_incld,nc_incld,odt
real(rtype), intent(inout) :: ncheti_cnt,qcheti_cnt,nicnt,qicnt,ninuc_cnt,qinuc_cnt

real(rtype) :: mi0l

! minimum mass of new crystal due to freezing of cloud droplets done
! externally (kg)
real(rtype) :: mi0l_min

mi0l_min = 4.0_rtype/3.0_rtype*pi*rho_h2o*(4.e-6_rtype)**3.0_rtype
mi0l = qc_incld/max(nc_incld,1.0e6_rtype/rho)
mi0l = max(mi0l_min, mi0l)

select case (Iflag)
case (1)  ! cloud droplet immersion freezing         
   if(qc_incld > qsmall) then
      ncheti_cnt = frzimm*1.0e6_rtype/rho ! frzimm input is in [#/cm3]
      qcheti_cnt = ncheti_cnt*mi0l
   else
      ncheti_cnt = 0.0_rtype
      qcheti_cnt = 0.0_rtype 
   endif
case (2)  ! deposition freezing / contact freezing 
   if(qc_incld > qsmall) then
      nicnt = frzcnt*1.0e6_rtype/rho
      qicnt = nicnt*mi0l
      ninuc_cnt = frzdep*1.0e6_rtype/rho
      qinuc_cnt = ninuc_cnt*mi0
   else
      nicnt = 0.0_rtype
      qicnt = 0.0_rtype
      ninuc_cnt = 0.0_rtype
      qinuc_cnt = 0.0_rtype 
   endif
end select

return
end subroutine CNT_couple

end module micro_p3
