module mam_soaexch_vbs

implicit none

contains
  subroutine mam_soaexch_vbs_1subarea(                          &
     nstep,             lchnk,                                  &
     i,                 k,                jsub,                 &
     latndx,            lonndx,           lund,                 &
     dtsubstep,                                                 &
     temp,              pmid,             aircon,               &
     n_mode,                                                    &
     qgas_cur,          qgas_avg,                               &
     qaer_cur,                                                  &
     qnum_cur,                                                  &
     qwtr_cur,                                                  &
     uptkaer                                                    )
  !
  ! calculate soa condensation/evaporation for i,k,jsub over time dtsubstep
  ! for the mam vbs soa mechanism
  !
  ! currently this does a modified version of the nvsoa mechanism
  !    of shrivastava et al. (2015) (doi 10.1002/2014jd022563)
  ! the main difference is that soa and condensible organic vapors from
  !    biomass-burning (and biofuel), biogenic, and fossil-fuel sources
  !    are all lumped together rather than treated separately
  !

  ! uses
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use modal_aero_data, only: lptr2_soa_a_amode
  use cam_abortutils,  only: endrun
  use modal_aero_data, only: ntot_amode, nsoag, nsoa, npoa
  use modal_aero_amicphys, only: igas_soag, igas_soag_end, nufi, mode_aging_optaa, &
                                 max_aer, max_gas, npca, mw_gas, iaer_pom, iaer_soa, max_mode
  use cam_logfile,     only: iulog
  use cam_history,     only: outfld
  use phys_debug_util, only: phys_debug_col
  
      implicit none

! arguments
      integer,  intent(in) :: nstep                 ! model time-step number
      integer,  intent(in) :: lchnk                 ! chunk identifier
      integer,  intent(in) :: i, k                  ! column and level indices
      integer,  intent(in) :: jsub                  ! sub-area index
      integer,  intent(in) :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in) :: lund                  ! logical unit for diagnostic output
      integer,  intent(in) :: n_mode                ! current number of modes (including temporary)

      real(r8), intent(in) :: dtsubstep        ! current integration timestep (s)
      real(r8), intent(in) :: temp             ! temperature (K)
      real(r8), intent(in) :: pmid             ! pressure at model levels (Pa)
      real(r8), intent(in) :: aircon           ! air molar concentration (kmol/m3)

      real(r8), intent(inout), dimension( 1:max_gas ) :: &
         qgas_cur, qgas_avg
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qnum_cur
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qwtr_cur
      real(r8), intent(in   ), dimension( 1:max_gas, 1:max_mode ) :: &
         uptkaer
      !real(r8), intent(inout ), dimension( 1:max_gas, 1:max_mode ) :: &
       !  uptkaer ! QZR test uptkaer =100 changes soa budget
! local
      integer, parameter :: maxd_poaspec = npoa
      integer, parameter :: maxd_soaspec = 7

      integer :: iaer, igas, ip
      integer :: isoa_bgn, isoa_end
      integer :: ll
      integer :: n, niter, niter_max
      integer :: ntot_poaspec, ntot_soaspec
      integer :: ntot_soamode

      integer :: icol

      logical, parameter :: flag_pcarbon_opoa_frac_zero   = .true.

      logical :: skip_soamode(max_mode)   ! true if this mode does not have soa


      real(r8), parameter :: a_min1 = 1.0e-20
      real(r8), parameter :: g_min1 = 1.0e-20
      real(r8), parameter :: alpha_astem = 0.05_r8 ! parameter used in calc of time step
      real(r8), parameter :: dtsub_fixed = -1.0    ! fixed sub-step for time integration (s)
!     real(r8), parameter :: dtsub_fixed = 10.0    ! fixed sub-step for time integration (s)
      real(r8), parameter :: rgas = 8.3144_r8      ! gas constant in J/K/mol

      real(r8) :: a_ooa_sum_tmp(max_mode)          ! total ooa (=soa+opoa) in a mode
      real(r8) :: a_opoa(max_mode)                 ! oxidized-poa aerosol mixrat (mol/mol at actual mw)
      real(r8) :: a_soa(maxd_soaspec,max_mode)     ! soa aerosol mixrat (mol/mol at actual mw)
      real(r8) :: a_soa_prv(maxd_soaspec,max_mode) ! soa aerosol mixrat at beginning of current time sub-step
      real(r8) :: a_soa_tmp(maxd_soaspec,max_mode) ! temporary soa aerosol mixrat (mol/mol)
      real(r8) :: beta(maxd_soaspec,max_mode)      ! dtcur*xferrate
      real(r8) :: c0_soa_298(maxd_soaspec)         ! c0_soa_298 = soa gas equilib vapor conc (ug/m3) at 298 k
      real(r8) :: delh_vap_soa(maxd_soaspec)       ! delh_vap_soa = heat of vaporization for gas soa (J/mol)
      real(r8) :: del_g_soa_tmp(maxd_soaspec)
      real(r8) :: dtcur                            ! current time step (s)
      real(r8) :: dtfull                           ! full time step (s)
      real(r8) :: dtmax                            ! = (dtfull-tcur)
      real(r8) :: dtsum_g_soa_avg
      real(r8) :: g0_soa(maxd_soaspec)             ! ambient soa gas equilib mixrat (mol/mol at actual mw)
      real(r8) :: g_soa(maxd_soaspec)              ! soa gas mixrat (mol/mol at actual mw)
      real(r8) :: g_soa_prv(maxd_soaspec)          ! soa gas mixrat at beginning of current time sub-step
      real(r8) :: g_soa_avg(maxd_soaspec)          ! time-averaged soa gas mixrat
      real(r8) :: g_star(maxd_soaspec,max_mode)    ! soa gas mixrat that is in equilib
                                                   ! with each aerosol mode (mol/mol)
      real(r8) :: mw_poa(maxd_poaspec)             ! actual molec wght of poa
      real(r8) :: mw_soa(maxd_soaspec)             ! actual molec wght of soa
      real(r8) :: opoa_frac(maxd_poaspec,max_mode) ! fraction of poa that is opoa
      real(r8) :: phi(maxd_soaspec,max_mode)       ! "relative driving force"
      real(r8) :: p0_soa(maxd_soaspec)             ! soa gas equilib vapor presssure (atm)
      real(r8) :: p0_soa_298(maxd_soaspec)         ! p0_soa_298 = soa gas equilib vapor presssure (atm) at 298 k
      real(r8) :: sat(maxd_soaspec,max_mode)       ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
                                                   !    used by the numerical integration scheme -- it is not a saturation rato!
      real(r8) :: tcur                             ! current integration time (from 0 s)

      real(r8) :: tmpa, tmpb, tmpc

      real(r8) :: tot_soa(maxd_soaspec),tempvar(maxd_soaspec)             ! g_soa + sum( a_soa(:) )
      real(r8) :: uptkaer_soag(maxd_soaspec,max_mode)
      real(r8) :: uptkaer_soag_tmp(maxd_soaspec,max_mode)

      character(len=128) :: msg

!check if it enters this subroutine
        write(msg,'(a)') 'Within subroutine mam_soaexch_vbs_1subarea'
!        call endrun( msg )


      if (nsoag /= 7 .or. nsoa /= 1) then
         write(msg,'(a,2(1x,i5))') 'mam_soaexch_vbs_1subarea - bad nsoag, nsoa =', nsoag, nsoa
         call endrun( msg )
      end if

      ntot_poaspec = npoa
!
! in the nvsoa treatment, the condensible organic vapors partition to the particle phase,
!    and the resulting aerosol species are initially semi-volatile (svsoa)
! once in the particle phase, they quickly age (oligomerize) to non-volatile soa (nvsoa)
!
! ntot_soaspec here is the number of condensible organic vapors and corresponding svsoa species
!
! most of the calculations in this routine involve the dynamic partitioning (i.e., mass transfer),
!    and the nvsoa soa does not participate (or affect) the partitioning
! after partitioning has been calculated for the time step (dtsubstep),
!    all of the svsoa is immediately converted to nvsoa
! as a result, the final (i.e., end of time step) svsoa mixing ratios are always zero, 
!    so the svsoa species do not need to be transported in the model,
!    and are just "temporary variables" within this routine
!
      ntot_soaspec = nsoag

      isoa_bgn = igas_soag
      isoa_end = igas_soag_end

! calc ntot_soamode = "last" mode on which soa is allowed to condense
      ntot_soamode = 0
      do n = 1, ntot_amode
         if (n == nufi) cycle
         if (mode_aging_optaa(n) > 0) ntot_soamode = n
         if (lptr2_soa_a_amode(n,1) > 0) ntot_soamode = n
      end do
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if ( i*k == top_lev .and. ldiagd1 ) write(lund,'(/a,5i5)') &
         'ntot_amode, ntot_amode_extd, n_mode, ntot_soamode', &
          ntot_amode, ntot_amode_extd, n_mode, ntot_soamode
#endif

! in the shrivastava et al. (2015) nvsoa treatment, poa is assumed to not contribute
!    to the raoults law calculation of effective saturation vapor pressures, so opoa_frac = 0.0
      opoa_frac = 0.0_r8
! for primary carbon mode, set opoa_frac=0 for consistency with older code
! (this could be changed)
      !if ( flag_pcarbon_opoa_frac_zero ) then
       !  if (npca > 0) opoa_frac(:,npca) = 0.0_r8
      !end if   !QZR commented out 11March2022

      c0_soa_298 = 0.1_r8
      delh_vap_soa = 100.0e3
      ! c0 an delh values from shrivastava et al. (2015)
      !                      SOAG15     SOAG24    SOAG35     SOAG34     SOAG33   SOAG32    SOAG31
      delh_vap_soa(1:7) = (/   82.0_r8,  88.0_r8,   82.0_r8,  88.0_r8,  94.0_r8, 100.0_r8, 106.0_r8 /) * 1.0e3_r8
      c0_soa_298(  1:7) = (/ 1000.0_r8, 100.0_r8, 1000.0_r8, 100.0_r8,  10.0_r8,   1.0_r8,   0.1_r8 /)

      do ll = 1, ntot_soaspec
         igas = igas_soag + ll - 1
         ! convert sat vapor conc from ug/m^3 to mol/m^3 then to mol/liter
         tmpa = (c0_soa_298(ll)*1.0e-6_r8/mw_gas(igas)) * 1.0e-3_r8  
         ! calc sat vapor pressure (atm) from molar-conc and temp [ 0.082056 = gas constant in (atm/deg-K/(mol/liter)) ]
         p0_soa_298(ll) = 0.082056_r8*tmpa*temp  
      end do

! calc soa gas saturation molar-mixing-ratio at local temp and air-pressure
      do ll = 1, ntot_soaspec
         p0_soa(ll) = p0_soa_298(ll) * &
                  exp( -(delh_vap_soa(ll)/rgas)*((1.0/temp)-(1.0/298.0)) )
         g0_soa(ll) = 1.01325e5*p0_soa(ll)/pmid
      end do

      niter_max = 1000
      niter = 0
      dtfull = dtsubstep
      tcur = 0.0_r8
      dtcur = 0.0_r8
      phi(:,:) = 0.0_r8
      g_star(:,:) = 0.0_r8
      g_soa(:) = 0.0_r8
      g_soa_prv(:) = 0.0_r8
      a_opoa(:) = 0.0_r8
      a_soa(:,:) = 0.0_r8
      a_soa_prv(:,:) = 0.0_r8
      g_soa_avg(:) = 0.0_r8
      dtsum_g_soa_avg = 0.0_r8
      !uptkaer(:,:) = 0.0_r8 !QZR
      uptkaer_soag(:,:) = 0.0_r8

! initialize g_soa_prv, a_soa_prv, uptkaer_soag
! the a_soa, which involves only the svsoa species, and is initially zero
      do ll = 1, ntot_soaspec
         igas = igas_soag + ll - 1
         g_soa_prv(ll) = qgas_cur(igas)
         do n = 1, ntot_soamode
            !uptkaer(igas,n) = 100.0_r8 ! check if soa changes QZR test
            uptkaer_soag(ll,n) = uptkaer(igas,n) !100.0_r8 !uptkaer(igas,n)
         end do

         iaer = iaer_soa + ll - 1
         do n = 1, ntot_soamode
            a_soa_prv(ll,n) = 0.0_r8
         end do
      end do

!Find out if LAT-LON column exists in this lchunk or not (if icol==0, column
!with LAT-LON doesn’t exist in this chunk)
      !icol = phys_debug_col(lchnk)!(state%lchnk) !QZR commented out

      !if (icol > 0) then !MS
       ! write(iulog,*)'Incoming qgas_cur,uptkaer_soag',qgas_cur(igas_soag:igas_soagzz),uptkaer_soag(igas_soag:igas_soagzz,1), &
        !  uptkaer_soag(igas_soag:igas_soagzz,2),uptkaer_soag(igas_soag:igas_soagzz,3)
      !endif  !QZR commented out diagnostic write

! calc oxygenated poa (which does not change during the soa uptake integration)
      do n = 1, ntot_soamode
         a_opoa(n) = 0.0_r8
         do ll = 1, ntot_poaspec
            a_opoa(n) = a_opoa(n) + opoa_frac(ll,n) * max( qaer_cur(iaer_pom+ll-1,n), 0.0_r8 )
         end do
      end do


!
! main integration loop -- does multiple substeps to reach dtfull
!

time_loop: &
      do while (tcur < dtfull-1.0e-3_r8 )

      niter = niter + 1
      if (niter > niter_max) exit


! set qxxx_prv to be current value
      if (niter > 1) then
         g_soa_prv = g_soa
         a_soa_prv = a_soa
      end if


! determine which modes have non-zero transfer rates
!    and are involved in the soa gas-aerosol transfer
! for diameter = 1 nm and number = 1 #/cm3, xferrate ~= 1e-9 s-1
      do n = 1, ntot_soamode
         skip_soamode(n) = .true.
         do ll = 1, ntot_soaspec
            if (uptkaer_soag(ll,n) > 1.0e-15_r8) then
               uptkaer_soag_tmp(ll,n) = uptkaer_soag(ll,n)
               skip_soamode(n) = .false.
            else
               uptkaer_soag_tmp(ll,n) = 0.0_r8
            end if
         end do
      end do

! load incoming soag and soaa into temporary arrays
! force things to be non-negative
! calc tot_soa(ll)
! calc a_opoa (always slightly >0)
!
! *** questions ***
! > why not use qgas and qaer instead of g_soa and a_soa
! > why not calc the following on every substep because 
!      nuc and coag may change things:
!      skip)soamode, uptkaer_soag_tmp, tot_soa, a_opoa
! > include gasprod for soa ??
! > create qxxx_bgn = qxxx_cur at the very beginning (is it needed)
!
      do ll = 1, ntot_soaspec
         g_soa(ll) = max( g_soa_prv(ll), 0.0_r8 )
         tot_soa(ll) = g_soa(ll)
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            a_soa(ll,n) = max( a_soa_prv(ll,n), 0.0_r8 )
            tot_soa(ll) = tot_soa(ll) + a_soa(ll,n)
         end do
      end do


! determine time step
      tmpa = 0.0_r8  ! time integration parameter for all soa species
      do n = 1, ntot_soamode
         if ( skip_soamode(n) ) cycle
         a_ooa_sum_tmp(n) = a_opoa(n) + sum( a_soa(1:ntot_soaspec,n) )
      end do
      do ll = 1, ntot_soaspec
         tmpb = 0.0_r8  ! time integration parameter for a single soa species
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            sat(ll,n) = g0_soa(ll)/max( a_ooa_sum_tmp(n), a_min1 )
            g_star(ll,n) = sat(ll,n)*a_soa(ll,n)
            phi(ll,n) = (g_soa(ll) - g_star(ll,n))/max( g_soa(ll), g_star(ll,n), g_min1 )
            tmpb = tmpb + uptkaer_soag_tmp(ll,n)*abs(phi(ll,n))
         end do
         tmpa = max( tmpa, tmpb )
      end do

      if (dtsub_fixed > 0.0_r8) then
         dtcur = dtsub_fixed
         tcur = tcur + dtcur
      else
         dtmax = dtfull-tcur
         if (dtmax*tmpa <= alpha_astem) then
! here alpha_astem/tmpa >= dtmax, so this is final substep
            dtcur = dtmax
            tcur = dtfull
         else
            dtcur = alpha_astem/tmpa
            tcur = tcur + dtcur
         end if
      end if

      !if (icol > 0) then !MS
       !     write(iulog,*)'within timestep',dtfull,dtcur,tcur,sat(:,:),g_star(:,:),phi(:,:),tmpb,tmpa,a_soa(:,:),tot_soa(:)
      !endif !QZR commented out

! step 1 - for modes where soa is condensing, estimate "new" a_soa(ll,n)
!    using an explicit calculation with "old" g_soa
!    and g_star(ll,n) calculated using "old" a_soa(ll,n)
! do this to get better estimate of "new" a_soa(ll,n) and sat(ll,n)
      do n = 1, ntot_soamode
         if ( skip_soamode(n) ) cycle
         do ll = 1, ntot_soaspec
            ! first ll loop calcs a_soa_tmp(ll,n) & a_ooa_sum_tmp
            a_soa_tmp(ll,n) = a_soa(ll,n)
            beta(ll,n) = dtcur*uptkaer_soag_tmp(ll,n)
            del_g_soa_tmp(ll) = g_soa(ll) - g_star(ll,n)
            if (del_g_soa_tmp(ll) > 0.0_r8) then
               a_soa_tmp(ll,n) = a_soa(ll,n) + beta(ll,n)*del_g_soa_tmp(ll)
            end if
         end do
         a_ooa_sum_tmp(n) = a_opoa(n) + sum( a_soa_tmp(1:ntot_soaspec,n) )
         do ll = 1, ntot_soaspec
            ! second ll loop calcs sat & g_star
            if (del_g_soa_tmp(ll) > 0.0_r8) then
               sat(ll,n) = g0_soa(ll)/max( a_ooa_sum_tmp(n), a_min1 )
               g_star(ll,n) = sat(ll,n)*a_soa_tmp(ll,n)   ! this just needed for diagnostics
            end if
         end do
      end do

      !if (icol > 0) then !MS
       !        write(iulog,*)'at the end of  step 1',beta(:,:),uptkaer_soag_tmp(:,:),a_soa_tmp(:,:),sat(:,:),g_star(:,:)
      !endif !QZR commented out

! step 2 - implicit in g_soa and semi-implicit in a_soa,
!    with g_star(ll,n) calculated semi-implicitly
      do ll = 1, ntot_soaspec
         tmpa = 0.0_r8
         tmpb = 0.0_r8
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            tmpa = tmpa + a_soa(ll,n)/(1.0_r8 + beta(ll,n)*sat(ll,n))
            tmpb = tmpb + beta(ll,n)/(1.0_r8 + beta(ll,n)*sat(ll,n))
         end do

         g_soa(ll) = (tot_soa(ll) - tmpa)/(1.0_r8 + tmpb)
         g_soa(ll) = max( 0.0_r8, g_soa(ll) )
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            a_soa(ll,n) = (a_soa(ll,n) + beta(ll,n)*g_soa(ll))/   &
                       (1.0_r8 + beta(ll,n)*sat(ll,n))
         end do
      end do

      !if (icol > 0) then !MS
       !         write(iulog,*)'at the end of  step 2',beta(:,:),sat(:,:),a_soa(:,:),g_soa(:)
      !endif !QZR commented out

! increment g_soa_avg
      do ll = 1, ntot_soaspec
         tmpc = g_soa(ll) - g_soa_prv(ll)
         g_soa_avg(ll) = g_soa_avg(ll) + dtcur*(g_soa_prv(ll) + 0.5_r8*tmpc)
      end do


      dtsum_g_soa_avg = dtsum_g_soa_avg + dtcur

      end do time_loop


! update mix ratios for soa gas species,
! and convert g_soa_avg from sum_over[ qgas*dt_cut ] to an time averaged value
      do ll = 1, ntot_soaspec
         igas = igas_soag + ll - 1
         qgas_cur(igas) = g_soa(ll)
         qgas_avg(igas) = max( 0.0_r8, g_soa_avg(ll)/dtsum_g_soa_avg )
      end do

      !if (icol > 0) then !MS
       !  write(iulog,*)'after time loop',qgas_cur(:),qgas_avg(:),dtsum_g_soa_avg
      !endif !QZR commented out

! Print qaer_cur but for each species and sum of all mode QZR MS
     do ll = 1, ntot_soaspec
      tempvar(ll)=0.0_r8
     enddo


     do ll = 1, ntot_soaspec
        do n = 1, ntot_soamode
          if ( skip_soamode(n) ) cycle
          tempvar(ll)=tempvar(ll)+a_soa(ll,n)
        enddo
     enddo

! for each mode, sum up all of the newly condensed svsoa
! and "oligomerize" it to the single nvsoa species
      do n = 1, ntot_soamode
         if ( skip_soamode(n) ) cycle
         tmpa = 0.0_r8
         do ll = 1, ntot_soaspec
            tmpa = tmpa + a_soa(ll,n)
         end do
         iaer = iaer_soa
         qaer_cur(iaer,n) = qaer_cur(iaer,n) + tmpa
      end do

!Check for constant qgas_cur(igas_soag:igas_soagzz) &
! qaer_cur(iaer_soa,1), qaer_cur(iaer_soa,2), qaer_cur(iaer_soa,3)

     ! qgas_cur(igas_soag) = 1.0_r8
     ! qgas_cur(igas_soag+1) = 2.0_r8
     ! qgas_cur(igas_soag+2) = 3.0_r8
     ! qgas_cur(igas_soag+3) = 4.0_r8
     ! qgas_cur(igas_soag+4) = 5.0_r8
     ! qgas_cur(igas_soag+5) = 6.0_r8
     ! qgas_cur(igas_soag+6) = 7.0_r8


     ! qaer_cur(iaer_soa,1) = 5.0_r8
     ! qaer_cur(iaer_soa,2) = 10.0_r8
     ! qaer_cur(iaer_soa,3) = 15.0_r8



!Find out if LAT-LON column exists in this lchunk or not (if icol==0, column
!with LAT-LON doesn’t exist in this chunk)
     ! icol = phys_debug_col(lchnk)!(state%lchnk)
 
!if icol>0; column with LAT-LON exists in this chunk
!To print qgas_cur for 1, ntot_soaspec and qaer_cur for 1, ntot_soamode

      !if (icol > 0) then 
        !write(iulog,*)’q:’,state%q(icol,1,1) 
        ! write(iulog,*)'MAM_soaexch_vbs qgas_cur, ntot_soaspec, qaer_cur=',qgas_cur,ntot_soaspec,qaer_cur
        !write(iulog,*)'MAM_soaexch_vbs ntot_soamode =',ntot_soamode
        !write(iulog,*)'MAM_soaexch_vbs max_aer, max_gas, npca, mw_gas, iaer_pom, iaer_soa, max_mode:', max_aer, max_gas, npca, mw_gas, iaer_pom, iaer_soa, max_mode
       ! write(iulog,*)'MAM_soaexchvbs: qaer_cur_sumofallmodes or tempvar(igas_soag:igas_soagzz)',tempvar(igas_soag:igas_soagzz)
       ! write(iulog,*)'MAM_saoexch_vbs: lchnk,lev,col, pmid, temp ', lchnk, k, i, pmid, temp 
       ! write(iulog,*)'MAM_soaexch_vbs: mw_gas(igas_soag:igas_soagzz), opoa_frac(1:1) ', mw_gas(igas_soag:igas_soagzz), opoa_frac(1,1) 
       ! write(iulog,*)'MAM_soaexch_vbs qgas_cur(igas_soagigas_soagzz), iaer_pom, qaer_cur(iaer_pom,1) =',qgas_cur(igas_soag:igas_soagzz), iaer_pom, qaer_cur(iaer_pom,1)
       !  write(iulog,*)'MAM_soaexch_vbs qaer_cur(iaer_soa,1) mode 1, iaer_soa=',qaer_cur(iaer_soa,1), iaer_soa
       !  write(iulog,*)'MAM_soaexch_vbs qaer_cur(iaer_soa,2) mode 2, iaer_soa=',qaer_cur(iaer_soa,2), iaer_soa
       !  write(iulog,*)'MAM_soaexch_vbs qaer_cur(iaer_soa,2) mode 3, iaer_soa=',qaer_cur(iaer_soa,3), iaer_soa
       !  write(iulog,*)'g0_soa(1:7)=',g0_soa(1:7)
       ! write(iulog,*)'igas_soag, igas_soagzz,uptkaer(igas_soag:igas_soagzz,nmode) MAM_soaexch_vbs=',igas_soag, igas_soagzz, uptkaer(igas_soag:igas_soagzz,:)
      !endif !QZR commented out


      return
      end subroutine mam_soaexch_vbs_1subarea





end module mam_soaexch_vbs
