#define MODHS 1
subroutine tphysidl(ztodt, state, tend, ideal_phys_option)
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
                                 physics_ptend_init
   use physics_update_mod, only: physics_update
   use physconst,          only: gravit, cappa, rair, cpair
   use cam_abortutils,         only: endrun
   use ref_pres,           only: pref_mid_norm, psurf_ref
   use cam_history,        only: outfld
   use cam_logfile,        only: iulog
   use time_manager,       only: get_nstep
   use check_energy,       only: check_energy_chng
   use ref_pres,           only: ptop_ref  !--JH--

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

   ! Flag to choose which idealized physics
   !   1: Held-Suarez
   !   2: Held-Suarez with Williamson modification for stratosphere
   !   3: Held-Suarez with Lin-Williamson modification for stratosphere/mesosphere
   character(len=*), intent(in) :: ideal_phys_option
!
!---------------------------Local workspace-----------------------------
!
   type(physics_ptend) :: ptend

   integer :: lchnk                            ! chunk identifier
   integer :: ncol                             ! number of atmospheric columns
   integer :: nstep                            ! timestep number

   real(r8) clat(pcols)                        ! latitudes(radians) for columns
   real(r8) pmid(pcols,pver)                   ! mid-point pressure
   integer  i,k                                ! Longitude, level indices
   
   real(r8) tmp                                ! temporary
   real(r8) :: zero(pcols)
   
   real(r8) cossq(pcols)                       ! coslat**2
   real(r8) cossqsq(pcols)                     ! coslat**4
   real(r8) sinsq(pcols)                       ! sinlat**2
   real(r8) coslat(pcols)                      ! cosine(latitude)
   real(r8) teq(pcols, pver)                   ! for writing out eq temp.

   ! JH
   ! ---------- Forcing parameters for original Held-Suarez ----------
   !
   real(r8) efoldf                             ! efolding time for wind dissipation
   real(r8) efolda                             ! efolding time for T dissipation
   real(r8) efolds                             ! efolding time for T dissipation
   real(r8) sigmab                             ! threshold sigma level
   real(r8) t00                                ! minimum reference temperature
   real(r8) kf                                 ! 1./efolding_time for wind dissipation
   real(r8) ka                                 ! 1./efolding_time for temperature diss.
   real(r8) ks                                 ! 1./efolding_time for temperature diss. 
   real(r8) kv                                 ! 1./efolding_time (normalized) for wind
   real(r8) kt                                 ! 1./efolding_time for temperature diss.
   real(r8) trefa                              ! "radiative equilibrium" T
   real(r8) trefc                              ! used in calc of "radiative equilibrium" T
   
   ! JH
   ! ---------- Forcing parameters for Williamson modification ----------
   !
   real(r8) pi                                 ! 3.14159...
   real(r8) efoldaa                            ! efolding time for T dissipation
   real(r8) kaa                                ! 1./efolding_time for temperature diss.
   real(r8) p0strat                            ! threshold pressure
   real(r8) phi0                               ! threshold latitude
   real(r8) dphi0                              ! del-latitude
   real(r8) a0                                 ! coefficient
   real(r8) aeq                                ! 100 mb
   real(r8) apole                              ! 2   mb
   real(r8) lapsew                             ! lapse rate
   real(r8) constw                             ! constant
   real(r8) lapsec                             ! lapse rate
   real(r8) constc                             ! constant
   real(r8) acoslat                            ! abs(acos(coslat))
   real(r8) pmax                               ! pressure given by max(pmid, apole)
   
   ! JH
   ! ---------- Forcing parameters for Lin/Williamson modification ----------
   ! 
   real(r8) onemsig                            ! 1. - sigma_reference
   real(r8) efold_strat                        ! efolding time for T dissipation in Strat
   real(r8) efold_meso                         ! efolding time for T dissipation in Meso
   real(r8) efoldv                             ! efolding time for wind dissipation
   real(r8) p_infint                           ! effective top of model
   real(r8) lapse                              ! lapse rate
   real(r8) h0                                 ! scale height (7 km)
   real(r8) pressmb                            ! model pressure in mb
   
   ! JH
   ! ---------- Parameters for FV3-type sponge Rayleigh friction ----------
   ! 
   real(r8)   :: pih                                        ! 0.5*pi
   real(r8)   :: kr                                         ! RF friction coefficient
   real(r8)   :: num                                        ! kr definition numerator
   real(r8)   :: den                                        ! kr definition denominator
   real(r8)   :: fv3_rf_cutoff = 0.001_r8                   ! 100 Pa in sigma coordinate
   real(r8)   :: fv3_tau_rev   = 1.0_r8/(86400.0_r8*3.0_r8) ! 1/(3 days)

!
!-----------------------------------------------------------------------
!
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


   if (trim(ideal_phys_option) == 'held-suarez') then
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

   elseif (trim(ideal_phys_option) == 'held-suarez-williamson') then
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
      pi      = 4._r8*atan(1._r8)
      pih     = 2._r8*atan(1._r8)
      phi0   = 60._r8*pi/180._r8
      dphi0  = 15._r8*pi/180._r8
      efoldaa = 40._r8               ! T relaxation e-folding in days
      kaa     = 1._r8/(86400._r8*efoldaa)
      a0     = 2.65_r8/dphi0
      aeq    = 10000._r8             ! JH: 100 mb; apply vanilla HS below this pressure
      apole  = 200._r8              ! JH: original value from Williamson+ modificaiton
      !apole  = (ptop_ref * 0.9_r8)   ! JH: place apole slightly above model top
      lapsew = -3.345e-03_r8         ! lapse rate at the poles
      constw = rair*lapsew/gravit    
      lapsec =  2.00e-03_r8          ! lapse rate in tropics
      constc = rair*lapsec/gravit
     
      ! JH
      ! --- inherited from standard HS --- 
      efoldf  =  1._r8               ! e-folding time for wind dissipation
      efolda  = 40._r8               ! e-folding time for temp
      efolds  =  4._r8               ! e-folding time for temp
      sigmab  =  0.7_r8              ! threshold for surface wind damping, latitudinal temp damping
      t00     = 200._r8              ! min temp eq
      onemsig = 1._r8 - sigmab        
      ka  = 1._r8/(86400._r8*efolda) 
      ks  = 1._r8/(86400._r8*efolds)  
      
      
      do k=1,pver
         if (pref_mid_norm(k) > sigmab) then
            ! below HS threshold of sigma=0.7; give T damping timescale latitude dependence
            do i=1,ncol  
               ! ---- standard HS
               kt = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
               tmp     = kt/(1._r8 + ztodt*kt)
               trefc   = 315._r8 - 60._r8*sinsq(i)
               trefa   = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
               trefa   = max(t00,trefa)
               
               ! ---- HSW modifications
               acoslat = abs(acos(coslat(i)))
               p0strat = aeq - (aeq - apole)*0.5_r8*(1._r8 + tanh(a0*(acoslat - phi0)))
               if (pmid(i,k) < 10000._r8) then
                  ! above HSW p_d threshold; temp decreases with const. lapse rate
                  trefa = t00*((pmid(i,k)/10000._r8))**constc
                  tmp   = kaa/(1._r8+ ztodt*kaa)
               endif
               if (pmid(i,k) < p0strat) then
                  ! above HSW p_i threshold; temp increases with const. lapse rate
                  trefa = trefa + t00*( ((pmid(i,k)/p0strat))**constw - 1._r8 )
                  tmp   = kaa/(1._r8+ ztodt*kaa)
               endif
               ptend%s(i,k) = (trefa - state%t(i,k))*tmp*cpair
               teq(i, k) = trefa
            end do
         else
            ! above HS threshold of sigma=0.7; T damping timescale is constant
            do i=1,ncol
               
               ! instread of pmid(i,k), use max(pmid(i, k), apole) to force the equillibrium
               ! profile to a constant everywhere above apole (~the Williamson model top)
               pmax = MAX(pmid(i, k), apole)

               tmp     = ka/(1._r8+ ztodt*ka)
               trefc   = 315._r8 - 60._r8*sinsq(i)
               trefa   = (trefc - 10._r8*cossq(i)*log((pmax/psurf_ref)))*(pmax/psurf_ref)**cappa
               trefa   = max(t00,trefa)
               
               ! ---- HSW modifications
               acoslat = abs(acos(coslat(i)))
               p0strat = aeq - (aeq - apole)*0.5_r8*(1._r8 + tanh(a0*(acoslat - phi0)))
               if (pmax < 10000._r8) then
                  ! above HSW p_d threshold; temp decreases with const. lapse rate
                  trefa = t00*((pmax/10000._r8))**constc
                  tmp   = kaa/(1._r8+ ztodt*kaa)
               endif
               if (pmax < p0strat) then
                  ! above HSW p_i threshold; temp increases with const. lapse rate
                  trefa = trefa + t00*( ((pmax/p0strat))**constw - 1._r8 )
                  tmp   = kaa/(1._r8+ ztodt*kaa)
               endif
               ptend%s(i,k) = (trefa - state%t(i,k))*tmp*cpair
               teq(i, k) = trefa
            end do
         endif
      end do
      
      !
      ! Add diffusion
      !
      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0._r8
            ptend%v(i,k) = 0._r8
         end do
      end do

      kf = 1._r8/(86400._r8*efoldf)
      
      do k=1,pver
      
         ! add rayleigh friction near the surface for the wind fields
         if (pref_mid_norm(k) > sigmab) then
            kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
            tmp = -kv/(1._r8+ ztodt*kv)
            do i=1,ncol
               ptend%u(i,k) = tmp*state%u(i,k)
               ptend%v(i,k) = tmp*state%v(i,k)
            end do
         endif

         ! JH
         ! add rayleigh friction in sponge layer
         ! pref_mid_norm(1) serves as the position of the model top
         ! (full level, actual half level has lower pressure)
         ! apply RF at the cutoff sigma level
         if(pref_mid_norm(k) < fv3_rf_cutoff) then
             num = pih*log(fv3_rf_cutoff/pref_mid_norm(k))
             den = log(fv3_rf_cutoff/pref_mid_norm(1))
                kr = fv3_tau_rev * (sin(num/den))**2._r8    ! FV3 RF coefficient
             do i = 1, ncol
                 ptend%u(i,k) = -kr*state%u(i,k)   ! tendency via explicit time stepping
                 ptend%v(i,k) = -kr*state%v(i,k)   ! tendency via explicit time stepping
             end do
         endif
      end do

   elseif (trim(ideal_phys_option) == 'held-suarez-lin-williamson') then
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
      write(iulog,*) 'ideal_phys_option: string for choosing desired type of idealized ', &
                     'physics is set incorrectly.'
      write(iulog,*) 'The valid options are held-suarez, held-suarez-williamson, or held-suarez-lin-williamson.'
      write(iulog,*) 'idlflag is currently set to: ', trim(ideal_phys_option)
      call endrun('ideal_phys_option: invalid option')
   endif

   ! write heating rate, equillibrium temperature to history file
   !call outfld('HS_HEAT', ptend%s(:,:), ncol, lchnk)
   call outfld('HS_HEAT', teq(:,:), ncol, lchnk)


   ! update the state and total physics tendency
   call physics_update(state, ptend, ztodt, tend)

   ! Can't turn on conservation error messages unless the appropriate heat
   ! surface flux is computed and supplied as an argument to
   ! check_energy_chng to account for how the ideal physics forcings are
   ! changing the total energy.
   call check_energy_chng(state, tend, "tphysidl", nstep, ztodt, zero, zero, zero, zero)

   call outfld('QRS', tend%dtdt, pcols, lchnk)

end subroutine tphysidl

