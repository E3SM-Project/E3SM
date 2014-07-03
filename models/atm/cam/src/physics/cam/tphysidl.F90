#define MODHS 1
subroutine tphysidl(ztodt, state, tend)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  algorithm 1: Held/Suarez IDEALIZED physics
!  algorithm 2: Held/Suarez IDEALIZED physics (Williamson modified stratosphere
!  algorithm 3: Held/Suarez IDEALIZED physics (Lin/Williamson modified strato/meso-sphere
!  algorithm 4: Boer/Denis  IDEALIZED physics
!
! Author: J. Olson
! 
!-----------------------------------------------------------------------
   use shr_kind_mod,       only: r8 => shr_kind_r8
   use ppgrid,             only: pcols, pver
   use phys_grid,          only: get_rlat_all_p
   use physics_types,      only: physics_state, physics_tend, physics_ptend, &
                                 physics_ptend_init, physics_update
   use physconst,          only: gravit, cappa, rair, cpair
   use abortutils,         only: endrun
   use ref_pres,           only: pref_mid_norm, psurf_ref
   use cam_history,        only: outfld
   use cam_logfile,        only: iulog
   use time_manager,       only: get_nstep
   use check_energy,       only: check_energy_chng

   implicit none

!
! Input arguments
!
   real(r8), intent(in) :: ztodt                   ! Two times model timestep (2 delta-t)

   type(physics_state),  intent(inout) :: state
!
! Output arguments
!
   type(physics_tend ), intent(inout) :: tend
!
!---------------------------Local workspace-----------------------------
!
   type(physics_ptend) :: ptend

   integer :: lchnk                                ! chunk identifier
   integer :: ncol                                 ! number of atmospheric columns
   integer :: nstep                                ! timestep number

   real(r8) clat(pcols)                        ! latitudes(radians) for columns
   real(r8) pmid(pcols,pver)                   ! mid-point pressure
   integer  i,k                                ! Longitude, level indices

   real(r8) tmp                                ! temporary
   real(r8) kf                                 ! 1./efolding_time for wind dissipation
   real(r8) ka                                 ! 1./efolding_time for temperature diss.
   real(r8) kaa                                ! 1./efolding_time for temperature diss.
   real(r8) ks                                 ! 1./efolding_time for temperature diss.
   real(r8) kv                                 ! 1./efolding_time (normalized) for wind
   real(r8) kt                                 ! 1./efolding_time for temperature diss.
   real(r8) trefa                              ! "radiative equilibrium" T
   real(r8) trefc                              ! used in calc of "radiative equilibrium" T
   real(r8) cossq(pcols)                       ! coslat**2
   real(r8) cossqsq(pcols)                     ! coslat**4
   real(r8) sinsq(pcols)                       ! sinlat**2
   real(r8) onemsig                            ! 1. - sigma_reference
   real(r8) efoldf                             ! efolding time for wind dissipation
   real(r8) efolda                             ! efolding time for T dissipation
   real(r8) efoldaa                            ! efolding time for T dissipation
   real(r8) efolds                             ! efolding time for T dissipation
   real(r8) efold_strat                        ! efolding time for T dissipation in Strat
   real(r8) efold_meso                         ! efolding time for T dissipation in Meso
   real(r8) efoldv                             ! efolding time for wind dissipation
   real(r8) p_infint                           ! effective top of model
   real(r8) constw                             ! constant
   real(r8) lapsew                             ! lapse rate
   real(r8) p0strat                            ! threshold pressure
   real(r8) phi0                               ! threshold latitude
   real(r8) dphi0                              ! del-latitude
   real(r8) a0                                 ! coefficient
   real(r8) aeq                                ! 100 mb
   real(r8) apole                              ! 2   mb
   real(r8) pi                                 ! 3.14159...
   real(r8) coslat(pcols)                      ! cosine(latitude)
   real(r8) acoslat                            ! abs(acos(coslat))
   real(r8) constc                             ! constant
   real(r8) lapsec                             ! lapse rate
   real(r8) lapse                              ! lapse rate
   real(r8) h0                                 ! scale height (7 km)
   real(r8) sigmab                             ! threshold sigma level
   real(r8) pressmb                            ! model pressure in mb
   real(r8) t00                                ! minimum reference temperature
   integer  idlflag                            ! Flag to choose which idealized physics
   real(r8) :: zero(pcols)
!
!-----------------------------------------------------------------------
!
   idlflag = 1

   lchnk = state%lchnk
   ncol  = state%ncol

   nstep   = get_nstep()
   zero    = 0._r8

   call get_rlat_all_p(lchnk, ncol, clat)
   do i=1,ncol
      coslat (i) = cos(clat(i))
      sinsq  (i) = sin(clat(i))*sin(clat(i))
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq (i)
   end do
   do k = 1, pver
      do i = 1, ncol
         pmid(i,k) = state%pmid(i,k)
      end do
   end do

   ! initialize individual parameterization tendencies
   call physics_ptend_init(ptend, state%psetcols, 'tphysidl', ls=.true., lu=.true., lv=.true.)

   if (idlflag == 1) then
!
!-----------------------------------------------------------------------
!
! Held/Suarez IDEALIZED physics algorithm:
!
!   Held, I. M., and M. J. Suarez, 1994: A proposal for the
!   intercomparison of the dynamical cores of atmospheric general
!   circulation models.
!   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
!
!-----------------------------------------------------------------------
!
! Compute idealized radiative heating rates (as dry static energy)
!
      efoldf =  1._r8
      efolda = 40._r8
      efolds =  4._r8
      sigmab =  0.7_r8
      t00    = 200._r8
!
      onemsig = 1._r8 - sigmab
!
      ka = 1._r8/(86400._r8*efolda)
      ks = 1._r8/(86400._r8*efolds)
!
      do k=1,pver
         if (pref_mid_norm(k) > sigmab) then
            do i=1,ncol
               kt = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
#ifdef MODHS
               tmp   = kt/(1._r8+ ztodt*kt)
#else
	       tmp = kt
#endif
               trefc   = 315._r8 - 60._r8*sinsq(i)
               trefa = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
               trefa    = max(t00,trefa)
               ptend%s(i,k) = (trefa - state%t(i,k))*tmp*cpair
            end do
         else
#ifdef MODHS
            tmp   = ka/(1._r8+ ztodt*ka)
#else
	    tmp = ka
#endif
            do i=1,ncol
               trefc   = 315._r8 - 60._r8*sinsq(i)
               trefa = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
               trefa    = max(t00,trefa)
               ptend%s(i,k) = (trefa - state%t(i,k))*tmp*cpair
            end do
         endif
      end do
!
! Add diffusion near the surface for the wind fields
!
      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0._r8
            ptend%v(i,k) = 0._r8
         end do
      end do

      kf = 1._r8/(86400._r8*efoldf)
!
      do k=1,pver
         if (pref_mid_norm(k) > sigmab) then
            kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
#ifdef MODHS
            tmp = -kv/(1._r8+ ztodt*kv)
#else
            tmp = -kv
#endif
            do i=1,ncol
               ptend%u(i,k) = tmp*state%u(i,k)
               ptend%v(i,k) = tmp*state%v(i,k)
            end do
         endif
      end do

   elseif (idlflag == 2) then
!
!-----------------------------------------------------------------------
!
! Modified Held/Suarez IDEALIZED physics algorithm
! (modified with Williamson stratosphere):
!
!   Williamson, D. L., J. G. Olson and B. A. Boville, 1998: A comparison
!   of semi--Lagrangian and Eulerian tropical climate simulations.
!   Mon. Wea. Rev., vol 126, pp. 1001-1012.
!
!-----------------------------------------------------------------------
!
! Compute idealized radiative heating rates (as dry static energy)
!
      efoldf  =  1._r8
      efolda  = 40._r8
      efoldaa = 40._r8
      efolds  =  4._r8
      sigmab  =  0.7_r8
      t00     = 200._r8
!
      onemsig = 1._r8 - sigmab
!
      ka  = 1._r8/(86400._r8*efolda)
      kaa = 1._r8/(86400._r8*efoldaa)
      ks  = 1._r8/(86400._r8*efolds)
!
      pi     = 4._r8*atan(1._r8)
      phi0   = 60._r8*pi/180._r8
      dphi0  = 15._r8*pi/180._r8
      a0     = 2.65_r8/dphi0
      aeq    = 10000._r8
      apole  = 200._r8
      lapsew = -3.345e-03_r8
      constw = rair*lapsew/gravit
      lapsec =  2.00e-03_r8
      constc = rair*lapsec/gravit
      do k=1,pver
         if (pref_mid_norm(k) > sigmab) then
            do i=1,ncol
               kt = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
               acoslat = abs(acos(coslat(i)))
               p0strat = aeq - (aeq - apole)*0.5_r8*(1._r8 + tanh(a0*(acoslat - phi0)))
               tmp     = kt/(1._r8+ ztodt*kt)
               trefc   = 315._r8 - 60._r8*sinsq(i)
               trefa   = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)** &
                                                                                     cappa
               trefa   = max(t00,trefa)
               if (pmid(i,k) < 10000._r8) then
                  trefa = t00*((pmid(i,k)/10000._r8))**constc
                  tmp   = kaa/(1._r8+ ztodt*kaa)
               endif
               if (pmid(i,k) < p0strat) then
                  trefa = trefa + t00*( ((pmid(i,k)/p0strat))**constw - 1._r8 )
                  tmp   = kaa/(1._r8+ ztodt*kaa)
               endif
               ptend%s(i,k) = (trefa - state%t(i,k))*tmp*cpair
            end do
         else
            do i=1,ncol
               acoslat = abs(acos(coslat(i)))
               p0strat = aeq - (aeq - apole)*0.5_r8*(1._r8 + tanh(a0*(acoslat - phi0)))
               tmp     = ka/(1._r8+ ztodt*ka)
               trefc   = 315._r8 - 60._r8*sinsq(i)
               trefa   = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)** &
                                                                                     cappa
               trefa   = max(t00,trefa)
               if (pmid(i,k) < 10000._r8) then
                  trefa = t00*((pmid(i,k)/10000._r8))**constc
                  tmp   = kaa/(1._r8+ ztodt*kaa)
               endif
               if (pmid(i,k) < p0strat) then
                  trefa = trefa + t00*( ((pmid(i,k)/p0strat))**constw - 1._r8 )
                  tmp   = kaa/(1._r8+ ztodt*kaa)
               endif
               ptend%s(i,k) = (trefa - state%t(i,k))*tmp*cpair
            end do
         endif
      end do
!
! Add diffusion near the surface for the wind fields
!
      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0._r8
            ptend%v(i,k) = 0._r8
         end do
      end do

      kf = 1._r8/(86400._r8*efoldf)
!
      do k=1,pver
         if (pref_mid_norm(k) > sigmab) then
            kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
            tmp = -kv/(1._r8+ ztodt*kv)
            do i=1,ncol
               ptend%u(i,k) = tmp*state%u(i,k)
               ptend%v(i,k) = tmp*state%v(i,k)
            end do
         endif
      end do

   elseif (idlflag == 3) then
!
!-----------------------------------------------------------------------
!
! Held/Suarez IDEALIZED physics algorithm:
! (modified with Lin/Williamson stratosphere/mesosphere):
!
!   Held, I. M., and M. J. Suarez, 1994: A proposal for the
!   intercomparison of the dynamical cores of atmospheric general
!   circulation models.
!   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
!
!-----------------------------------------------------------------------
!
! Compute idealized radiative heating rates (as dry static energy)
!
      efoldf      =  1._r8
      efolda      = 40._r8
      efolds      =  4._r8
      efold_strat = 40._r8
      efold_meso  = 10._r8
      efoldv      = 0.5_r8
      sigmab      = 0.7_r8
      lapse       = 0.00225_r8
      h0          = 7000._r8
      t00         = 200._r8
      p_infint    = 0.01_r8
!
      onemsig = 1._r8 - sigmab
!
      ka = 1._r8/(86400._r8*efolda)
      ks = 1._r8/(86400._r8*efolds)
!
      do k=1,pver
         if (pref_mid_norm(k) > sigmab) then
            do i=1,ncol
               kt    = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
               tmp   = kt/(1._r8+ ztodt*kt)
               trefc = 315._r8 - 60._r8*sinsq(i)
               trefa = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
               trefa = max(t00,trefa)
               ptend%s(i,k) = (trefa - state%t(i,k))*tmp*cpair
            end do
         else
            do i=1,ncol
               tmp     = ka/(1._r8+ ztodt*ka)
               pressmb = pmid(i,k)*0.01_r8
               trefc   = 315._r8 - 60._r8*sinsq(i)
               trefa   = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)** &
                                                                                     cappa
               trefa   = max(t00,trefa)
               if (pressmb <= 100._r8 .and. pressmb > 1._r8) then
                  trefa = t00 + lapse*h0*coslat(i)*log(100._r8/pressmb)
                  tmp   = (efold_strat-efold_meso)*log(pressmb)/log(100._r8)
                  tmp   = efold_meso + tmp
                  tmp   = 1._r8/(86400._r8*tmp)
                  tmp   = tmp/(1._r8+ ztodt*tmp)
               endif
               if (pressmb <= 1._r8 .and. pressmb > 0.01_r8) then
                  trefa = t00 + lapse*h0*coslat(i)*log(100._r8*pressmb)
                  tmp   = 1._r8/(86400._r8*efold_meso)
                  tmp   = tmp/(1._r8+ ztodt*tmp)
               endif
               if (pressmb <= 0.01_r8) then
                  tmp   = 1._r8/(86400._r8*efold_meso)
                  tmp   = tmp/(1._r8+ ztodt*tmp)
               endif
               ptend%s(i,k) = (trefa - state%t(i,k))*tmp*cpair
            end do
         endif
      end do
!
! Add diffusion near the surface for the wind fields
!
      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0._r8
            ptend%v(i,k) = 0._r8
         end do
      end do

      kf = 1._r8/(86400._r8*efoldf)
!
      do k=1,pver
         if (pref_mid_norm(k) > sigmab) then
            kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
            tmp = -kv/(1._r8+ ztodt*kv)
            do i=1,ncol
               ptend%u(i,k) = tmp*state%u(i,k)
               ptend%v(i,k) = tmp*state%v(i,k)
            end do
         else
            do i=1,ncol
               pressmb  = pmid(i,k)*0.01_r8
               if (pressmb <= 100._r8) then
                  kv       = 1._r8/(86400._r8*efoldv)
                  tmp      = 1._r8 + tanh(1.5_r8*log10(p_infint/pressmb))
                  kv       = kv*tmp
                  tmp      = -kv/(1._r8+ ztodt*kv)
                  ptend%u(i,k) = tmp*state%u(i,k)
                  ptend%v(i,k) = tmp*state%v(i,k)
               endif
            end do
         endif
      end do

   else
      write(iulog,*) 'TPHYSIDL: flag for choosing desired type of idealized ', &
                 'physics ("idlflag") is set incorrectly.'
      write(iulog,*) 'The valid options are 1, 2, or 3.'
      write(iulog,*) 'idlflag is currently set to: ',idlflag
      call endrun('tphysidl: invalid option')
   endif

   ! update the state and total physics tendency
   call physics_update(state, ptend, ztodt, tend)

   ! Can't turn on conservation error messages unless the appropriate heat
   ! surface flux is computed and supplied as an argument to
   ! check_energy_chng to account for how the ideal physics forcings are
   ! changing the total energy.
   call check_energy_chng(state, tend, "tphysidl", nstep, ztodt, zero, zero, zero, zero)

   call outfld('QRS', tend%dtdt, pcols, lchnk)

end subroutine tphysidl

