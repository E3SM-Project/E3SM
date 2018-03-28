 !  SVN:$Id: ice_therm_bl99.F90 1182 2017-03-16 19:29:26Z njeffery $
!=========================================================================
!
! Update ice and snow internal temperatures
! using Bitz and Lipscomb 1999 thermodynamics
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW
!          Elizabeth C. Hunke, LANL
!
! 2012: Split from ice_therm_vertical.F90

      module ice_therm_bl99

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c2, p01, p1, p5, puny, &
          rhoi, rhos, hs_min, cp_ice, cp_ocn, depressT, Lfresh, ksno, kice
      use ice_colpkg_shared, only: conduct, calc_Tsfc, solve_zsal
      use ice_therm_shared, only: ferrmax, l_brine, hfrazilmin
      use ice_warnings, only: add_warning

      implicit none
      save

      private
      public :: surface_fluxes, temperature_changes

      real (kind=dbl_kind), parameter :: &
         betak   = 0.13_dbl_kind, & ! constant in formula for k (W m-1 ppt-1)
         kimin   = 0.10_dbl_kind    ! min conductivity of saline ice (W m-1 deg-1)

!=======================================================================

      contains

!=======================================================================
!
! Compute new surface temperature and internal ice and snow
! temperatures.  Include effects of salinity on sea ice heat
! capacity in a way that conserves energy (Bitz and Lipscomb, 1999).
!
! New temperatures are computed iteratively by solving a tridiagonal
! system of equations; heat capacity is updated with each iteration.
! Finite differencing is backward implicit.
!
! See Bitz, C.M., and W.H. Lipscomb, 1999:
! An energy-conserving thermodynamic model of sea ice,
! J. Geophys. Res., 104, 15,669-15,677.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine temperature_changes (dt,                 & 
                                      nilyr,    nslyr,    &
                                      rhoa,     flw,      &
                                      potT,     Qa,       &
                                      shcoef,   lhcoef,   &
                                      fswsfc,   fswint,   &
                                      Sswabs,   Iswabs,   &
                                      hilyr,    hslyr,    &
                                      zqin,     zTin,     &
                                      zqsn,     zTsn,     &
                                      zSin,               & 
                                      Tsf,      Tbot,     &
                                      fsensn,   flatn,    &
                                      flwoutn,  fsurfn,   &
                                      fcondtopn,fcondbot, &
                                      einit,    l_stop,   &
                                      stop_label)

      use ice_therm_shared, only: surface_heat_flux, dsurface_heat_flux_dTsf

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), &
         intent(in) :: &
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         Tbot            ! ice bottom surface temperature (deg C)

      real (kind=dbl_kind), &
         intent(inout) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         fswint          ! SW absorbed in ice interior below surface (W m-2)

      real (kind=dbl_kind), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         einit           ! initial energy of melting (J m-2)

      real (kind=dbl_kind), dimension (nslyr), &
         intent(inout) :: &
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (nilyr), &
         intent(inout) :: &
         Iswabs          ! SW radiation absorbed in ice layers (W m-2)

      real (kind=dbl_kind), intent(inout):: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         flwoutn         ! upward LW at surface (W m-2)

      real (kind=dbl_kind), intent(out):: &
         fcondbot        ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), &
         intent(inout):: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), dimension (nilyr), &
         intent(inout) :: &
         zqin        , & ! ice layer enthalpy (J m-3)
         zTin            ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (nilyr), &
         intent(in) :: &
         zSin            ! internal ice layer salinities

      real (kind=dbl_kind), dimension (nslyr), &
         intent(inout) :: &
         zqsn        , & ! snow layer enthalpy (J m-3)
         zTsn            ! internal snow layer temperatures

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      character (len=*), intent(out) :: &
         stop_label   ! abort error message
 
     ! local variables

      integer (kind=int_kind), parameter :: &
         nitermax = 100  ! max number of iterations in temperature solver

      real (kind=dbl_kind), parameter :: &
         Tsf_errmax = 5.e-4_dbl_kind ! max allowed error in Tsf
                                     ! recommend Tsf_errmax < 0.01 K

      integer (kind=int_kind) :: &
         k           , & ! ice layer index
         niter       , & ! iteration counter in temperature solver
         nmat            ! matrix dimension

      logical (kind=log_kind) :: &
         l_snow      , & ! true if snow temperatures are computed
         l_cold          ! true if surface temperature is computed

      real (kind=dbl_kind) :: &
         Tsf_start   , & ! Tsf at start of iteration
         dTsf        , & ! Tsf - Tsf_start
         dTi1        , & ! Ti1(1) - Tin_start(1)
         dfsurf_dT   , & ! derivative of fsurf wrt Tsf
         avg_Tsi     , & ! = 1. if new snow/ice temps avg'd w/starting temps
         enew            ! new energy of melting after temp change (J m-2)

      real (kind=dbl_kind) :: &
         dTsf_prev   , & ! dTsf from previous iteration
         dTi1_prev   , & ! dTi1 from previous iteration
         dfsens_dT   , & ! deriv of fsens wrt Tsf (W m-2 deg-1)
         dflat_dT    , & ! deriv of flat wrt Tsf (W m-2 deg-1)
         dflwout_dT  , & ! deriv of flwout wrt Tsf (W m-2 deg-1)
         dt_rhoi_hlyr, & ! dt/(rhoi*hilyr)
         einex       , & ! excess energy from dqmat to ocean
         ferr            ! energy conservation error (W m-2)

      real (kind=dbl_kind), dimension (nilyr) :: &
         Tin_init    , & ! zTin at beginning of time step
         Tin_start   , & ! zTin at start of iteration
         dTmat       , & ! zTin - matrix solution before limiting
         dqmat       , & ! associated enthalpy difference
         Tmlts           ! melting temp, -depressT * salinity

      real (kind=dbl_kind), dimension (nslyr) :: &
         Tsn_init    , & ! zTsn at beginning of time step
         Tsn_start   , & ! zTsn at start of iteration
         etas            ! dt / (rho * cp * h) for snow layers

      real (kind=dbl_kind), dimension (nilyr+nslyr+1) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs         , & ! rhs of tri-diagonal matrix equation
         Tmat            ! matrix output temperatures

      real (kind=dbl_kind), dimension (nilyr) :: &
         etai            ! dt / (rho * cp * h) for ice layers

      real (kind=dbl_kind), dimension(nilyr+nslyr+1):: &
         kh              ! effective conductivity at interfaces (W m-2 deg-1)

      real (kind=dbl_kind) :: &
         ci          , & ! specific heat of sea ice (J kg-1 deg-1)
         avg_Tsf     , & ! = 1. if Tsf averaged w/Tsf_start, else = 0.
         Iswabs_tmp  , & ! energy to melt through fraction frac of layer
         Sswabs_tmp  , & ! same for snow
         dswabs      , & ! difference in swabs and swabs_tmp
         frac        , & ! fraction of layer that can be melted through
         dTemp           ! minimum temperature difference for absorption

      logical (kind=log_kind) :: &
         converged       ! = true when local solution has converged

      logical (kind=log_kind) , dimension (nilyr) :: &
         reduce_kh       ! reduce conductivity when T exceeds Tmlt

      character(len=char_len_long) :: &
         warning ! warning message
      
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      converged  = .false.
      l_snow     = .false.
      l_cold     = .true.
      fcondbot   = c0
      dTsf_prev  = c0
      dTi1_prev  = c0
      dfsens_dT  = c0
      dflat_dT   = c0
      dflwout_dT = c0  
      einex      = c0
      dt_rhoi_hlyr = dt / (rhoi*hilyr)  ! hilyr > 0
      if (hslyr > hs_min/real(nslyr,kind=dbl_kind)) &
           l_snow = .true.

      do k = 1, nslyr
         Tsn_init (k) = zTsn(k) ! beginning of time step
         Tsn_start(k) = zTsn(k) ! beginning of iteration
         if (l_snow) then
            etas(k) = dt/(rhos*cp_ice*hslyr)
         else
            etas(k) = c0
         endif
      enddo                     ! k

      do k = 1, nilyr
         Tin_init (k) =  zTin(k)   ! beginning of time step
         Tin_start(k) =  zTin(k)   ! beginning of iteration
         Tmlts    (k) = -zSin(k) * depressT
      enddo

      !-----------------------------------------------------------------
      ! Compute thermal conductivity at interfaces (held fixed during
      !  subsequent iterations).
      ! Ice and snow interfaces are combined into one array (kh) to
      !  simplify the logic.
      !-----------------------------------------------------------------

      call conductivity (l_snow,                    &
                         nilyr,    nslyr,           &
                         hilyr,    hslyr,           &
                         zTin,     kh,      zSin)    

      !-----------------------------------------------------------------
      ! Check for excessive absorbed solar radiation that may result in
      ! temperature overshoots. Convergence is particularly difficult
      ! if the starting temperature is already very close to the melting 
      ! temperature and extra energy is added.   In that case, or if the
      ! amount of energy absorbed is greater than the amount needed to
      ! melt through a given fraction of a layer, we put the extra 
      ! energy into the surface.
      ! NOTE: This option is not available if the atmosphere model
      !       has already computed fsurf.  (Unless we adjust fsurf here)
      !-----------------------------------------------------------------
!mclaren: Should there be an if calc_Tsfc statement here then?? 

#ifdef CCSMCOUPLED
      frac = c1
      dTemp = p01
#else
      frac = 0.9_dbl_kind
      dTemp = 0.02_dbl_kind
#endif
      if (solve_zsal) dTemp = p1  ! lower tolerance with dynamic salinity
      do k = 1, nilyr

         Iswabs_tmp = c0 ! all Iswabs is moved into fswsfc
         if (Tin_init(k) <= Tmlts(k) - dTemp) then
            if (l_brine) then
               ci = cp_ice - Lfresh * Tmlts(k) / (Tin_init(k)**2)
               Iswabs_tmp = min(Iswabs(k), &
                                frac*(Tmlts(k)-Tin_init(k))*ci/dt_rhoi_hlyr)
            else
               ci = cp_ice
               Iswabs_tmp = min(Iswabs(k), &
                                frac*(-Tin_init(k))*ci/dt_rhoi_hlyr)
            endif
         endif
         if (Iswabs_tmp < puny) Iswabs_tmp = c0

         dswabs = min(Iswabs(k) - Iswabs_tmp, fswint)

         fswsfc   = fswsfc + dswabs
         fswint   = fswint - dswabs
         Iswabs(k) = Iswabs_tmp

      enddo

#ifdef CCSMCOUPLED
      frac = 0.9_dbl_kind
#endif
      do k = 1, nslyr
         if (l_snow) then

            Sswabs_tmp = c0
            if (Tsn_init(k) <= -dTemp) then
               Sswabs_tmp = min(Sswabs(k), &
                                -frac*Tsn_init(k)/etas(k))
            endif
            if (Sswabs_tmp < puny) Sswabs_tmp = c0

            dswabs = min(Sswabs(k) - Sswabs_tmp, fswint)

            fswsfc   = fswsfc + dswabs
            fswint   = fswint - dswabs
            Sswabs(k) = Sswabs_tmp

         endif
      enddo

      !-----------------------------------------------------------------
      ! Solve for new temperatures.
      ! Iterate until temperatures converge with minimal energy error.
      !-----------------------------------------------------------------
      converged  = .false.

      do niter = 1, nitermax

      !-----------------------------------------------------------------
      ! Identify cells, if any, where calculation has not converged.
      !-----------------------------------------------------------------

         if (.not.converged) then

      !-----------------------------------------------------------------
      ! Allocate and initialize
      !-----------------------------------------------------------------

            converged = .true.
            dfsurf_dT = c0
            avg_Tsi   = c0
            enew      = c0
            einex     = c0

      !-----------------------------------------------------------------
      ! Update specific heat of ice layers.
      ! To ensure energy conservation, the specific heat is a function of
      ! both the starting temperature and the (latest guess for) the
      ! final temperature.
      !-----------------------------------------------------------------

            do k = 1, nilyr

               if (l_brine) then
                  ci = cp_ice - Lfresh*Tmlts(k) /  &
                                (zTin(k)*Tin_init(k))
               else
                  ci = cp_ice
               endif
               etai(k) = dt_rhoi_hlyr / ci

            enddo

            if (calc_Tsfc) then

      !-----------------------------------------------------------------
      ! Update radiative and turbulent fluxes and their derivatives
      ! with respect to Tsf.
      !-----------------------------------------------------------------

               ! surface heat flux
               call surface_heat_flux(Tsf    , fswsfc, &
                                      rhoa   , flw   , &
                                      potT   , Qa    , &
                                      shcoef , lhcoef, &
                                      flwoutn, fsensn, &
                                      flatn  , fsurfn)

               ! derivative of heat flux with respect to surface temperature
               call dsurface_heat_flux_dTsf(Tsf      , fswsfc    , &
                                            rhoa     , flw       , &
                                            potT     , Qa        , &
                                            shcoef   , lhcoef    , &
                                            dfsurf_dT, dflwout_dT, &
                                            dfsens_dT, dflat_dT  )

      !-----------------------------------------------------------------
      ! Compute conductive flux at top surface, fcondtopn.
      ! If fsurfn < fcondtopn and Tsf = 0, then reset Tsf to slightly less
      !  than zero (but not less than -puny).
      !-----------------------------------------------------------------
               
               if (l_snow) then
                  fcondtopn = kh(1) * (Tsf - zTsn(1))
               else
                  fcondtopn = kh(1+nslyr) * (Tsf - zTin(1))
               endif

               if (Tsf >= c0 .and. fsurfn < fcondtopn) &
                    Tsf = -puny

      !-----------------------------------------------------------------
      ! Save surface temperature at start of iteration
      !-----------------------------------------------------------------
               
               Tsf_start = Tsf

               if (Tsf < c0) then
                  l_cold = .true.
               else
                  l_cold = .false.
               endif

      !-----------------------------------------------------------------
      ! Compute elements of tridiagonal matrix.
      !-----------------------------------------------------------------
               
               call get_matrix_elements_calc_Tsfc (nilyr, nslyr, &
                                   l_snow,      l_cold,      &
                                   Tsf,         Tbot,        &
                                   fsurfn,      dfsurf_dT,   &
                                   Tin_init,    Tsn_init,    &
                                   kh,          Sswabs,      &
                                   Iswabs,                   &
                                   etai,        etas,        &
                                   sbdiag,      diag,        &
                                   spdiag,      rhs)   

            else
               
               call get_matrix_elements_know_Tsfc (nilyr, nslyr, &
                                   l_snow,      Tbot,        &
                                   Tin_init,    Tsn_init,    &
                                   kh,          Sswabs,      &
                                   Iswabs,                   &
                                   etai,        etas,        &
                                   sbdiag,      diag,        &
                                   spdiag,      rhs,         &
                                   fcondtopn)

            endif  ! calc_Tsfc

      !-----------------------------------------------------------------
      ! Solve tridiagonal matrix to obtain the new temperatures.
      !-----------------------------------------------------------------

            nmat = nslyr + nilyr + 1  ! matrix dimension

            call tridiag_solver (nmat,    sbdiag(:),   &
                                 diag(:), spdiag(:),   &
                                 rhs(:),  Tmat(:))

      !-----------------------------------------------------------------
      ! Determine whether the computation has converged to an acceptable
      ! solution.  Five conditions must be satisfied:
      !
      !    (1) Tsf <= 0 C.
      !    (2) Tsf is not oscillating; i.e., if both dTsf(niter) and
      !        dTsf(niter-1) have magnitudes greater than puny, then
      !        dTsf(niter)/dTsf(niter-1) cannot be a negative number
      !        with magnitude greater than 0.5.  
      !    (3) abs(dTsf) < Tsf_errmax
      !    (4) If Tsf = 0 C, then the downward turbulent/radiative 
      !        flux, fsurfn, must be greater than or equal to the downward
      !        conductive flux, fcondtopn.
      !    (5) The net energy added to the ice per unit time must equal 
      !        the net change in internal ice energy per unit time,
      !        withinic the prescribed error ferrmax.
      !
      ! For briny ice (the standard case), zTsn and zTin are limited
      !  to prevent them from exceeding their melting temperatures.
      !  (Note that the specific heat formula for briny ice assumes
      !  that T < Tmlt.)  
      ! For fresh ice there is no limiting, since there are cases
      !  when the only convergent solution has zTsn > 0 and/or zTin > 0.
      !  Above-zero temperatures are then reset to zero (with melting 
      !  to conserve energy) in the thickness_changes subroutine.
      !-----------------------------------------------------------------

            if (calc_Tsfc) then

      !-----------------------------------------------------------------
      ! Reload Tsf from matrix solution
      !-----------------------------------------------------------------

               if (l_cold) then
                  if (l_snow) then
                     Tsf = Tmat(1)
                  else
                     Tsf = Tmat(1+nslyr)
                  endif
               else                ! melting surface
                  Tsf = c0
               endif

      !-----------------------------------------------------------------
      ! Initialize convergence flag (true until proven false), dTsf,
      !  and temperature-averaging coefficients.
      ! Average only if test 1 or 2 fails.
      ! Initialize energy.
      !-----------------------------------------------------------------
               
               dTsf = Tsf - Tsf_start
               avg_Tsf  = c0

      !-----------------------------------------------------------------
      ! Condition 1: check for Tsf > 0
      ! If Tsf > 0, set Tsf = 0, then average zTsn and zTin to force
      ! internal temps below their melting temps.
      !-----------------------------------------------------------------

               if (Tsf > puny) then
                  Tsf = c0
                  dTsf = -Tsf_start
                  if (l_brine) avg_Tsi = c1   ! avg with starting temp
                  converged = .false.

      !-----------------------------------------------------------------
      ! Condition 2: check for oscillating Tsf
      ! If oscillating, average all temps to increase rate of convergence.
      !-----------------------------------------------------------------

               elseif (niter > 1 &                ! condition (2)
                    .and. Tsf_start <= -puny &
                    .and. abs(dTsf) > puny &
                    .and. abs(dTsf_prev) > puny &
                    .and. -dTsf/(dTsf_prev+puny*puny) > p5) then

                  if (l_brine) then ! average with starting temp
                     avg_Tsf  = c1    
                     avg_Tsi = c1
                  endif
                  dTsf = p5 * dTsf
                  converged = .false.
               endif

!!!            dTsf_prev = dTsf

      !-----------------------------------------------------------------
      ! If condition 2 failed, average new surface temperature with
      !  starting value.
      !-----------------------------------------------------------------
               Tsf  = Tsf &
                    + avg_Tsf * p5 * (Tsf_start - Tsf)

            endif   ! calc_Tsfc

            do k = 1, nslyr

      !-----------------------------------------------------------------
      ! Reload zTsn from matrix solution
      !-----------------------------------------------------------------

               if (l_snow) then
                  zTsn(k) = Tmat(k+1)
               else
                  zTsn(k) = c0
               endif
               if (l_brine) zTsn(k) = min(zTsn(k), c0)

      !-----------------------------------------------------------------
      ! If condition 1 or 2 failed, average new snow layer
      !  temperatures with their starting values.
      !-----------------------------------------------------------------
               zTsn(k) = zTsn(k) &
                       + avg_Tsi*p5*(Tsn_start(k)-zTsn(k))

      !-----------------------------------------------------------------
      ! Compute zqsn and increment new energy.
      !-----------------------------------------------------------------
               zqsn(k) = -rhos * (Lfresh - cp_ice*zTsn(k))
               enew  = enew + hslyr * zqsn(k)

               Tsn_start(k) = zTsn(k) ! for next iteration

            enddo                  ! nslyr

            dTmat(:) = c0
            dqmat(:) = c0
            reduce_kh(:) = .false.
            do k = 1, nilyr

      !-----------------------------------------------------------------
      ! Reload zTin from matrix solution
      !-----------------------------------------------------------------

               zTin(k) = Tmat(k+1+nslyr)

               if (l_brine .and. zTin(k) > Tmlts(k) - puny) then
                  dTmat(k) = zTin(k) - Tmlts(k)
                  dqmat(k) = rhoi * dTmat(k) &
                           * (cp_ice - Lfresh * Tmlts(k)/zTin(k)**2)
! use this for the case that Tmlt changes by an amount dTmlt=Tmltnew-Tmlt(k)
!                             + rhoi * dTmlt &
!                             * (cp_ocn - cp_ice + Lfresh/zTin(k))
                  zTin(k) = Tmlts(k)
                  reduce_kh(k) = .true.
               endif

      !-----------------------------------------------------------------
      ! Condition 2b: check for oscillating zTin(1)
      ! If oscillating, average all ice temps to increase rate of convergence.
      !-----------------------------------------------------------------

               if (k==1 .and. .not.calc_Tsfc) then
                  dTi1 = zTin(k) - Tin_start(k)

                  if (niter > 1 &                    ! condition 2b    
                      .and. abs(dTi1) > puny &
                      .and. abs(dTi1_prev) > puny &
                      .and. -dTi1/(dTi1_prev+puny*puny) > p5) then

                     if (l_brine) avg_Tsi = c1
                     dTi1 = p5 * dTi1
                     converged = .false.
                  endif
                  dTi1_prev = dTi1
               endif   ! k = 1 .and. calc_Tsfc = F

      !-----------------------------------------------------------------
      ! If condition 1 or 2 failed, average new ice layer
      !  temperatures with their starting values.
      !-----------------------------------------------------------------
               zTin(k) = zTin(k) &
                       + avg_Tsi*p5*(Tin_start(k)-zTin(k))

      !-----------------------------------------------------------------
      ! Compute zqin and increment new energy.
      !-----------------------------------------------------------------
               if (l_brine) then
                  zqin(k) = -rhoi * (cp_ice*(Tmlts(k)-zTin(k)) &
                                   + Lfresh*(c1-Tmlts(k)/zTin(k)) &
                                   - cp_ocn*Tmlts(k))
               else
                  zqin(k) = -rhoi * (-cp_ice*zTin(k) + Lfresh)
               endif
               enew = enew + hilyr * zqin(k)
               einex = einex + hilyr * dqmat(k)

               Tin_start(k) = zTin(k) ! for next iteration

            enddo                  ! nilyr

            if (calc_Tsfc) then

      !-----------------------------------------------------------------
      ! Condition 3: check for large change in Tsf
      !-----------------------------------------------------------------
               
               if (abs(dTsf) > Tsf_errmax) then
                  converged = .false.
               endif

      !-----------------------------------------------------------------
      ! Condition 4: check for fsurfn < fcondtopn with Tsf >= 0
      !-----------------------------------------------------------------
               
               fsurfn = fsurfn + dTsf*dfsurf_dT
               if (l_snow) then
                  fcondtopn = kh(1) * (Tsf-zTsn(1))
               else
                  fcondtopn = kh(1+nslyr) * (Tsf-zTin(1))
               endif

               if (Tsf >= c0 .and. fsurfn < fcondtopn) then
                  converged = .false.
               endif

               dTsf_prev = dTsf

            endif                   ! calc_Tsfc

      !-----------------------------------------------------------------
      ! Condition 5: check for energy conservation error
      ! Change in internal ice energy should equal net energy input.
      !-----------------------------------------------------------------

            fcondbot = kh(1+nslyr+nilyr) * &
                       (zTin(nilyr) - Tbot)

            ! Flux extra energy out of the ice
            fcondbot = fcondbot + einex/dt 

            ferr = abs( (enew-einit)/dt &
                 - (fcondtopn - fcondbot + fswint) )

            ! factor of 0.9 allows for roundoff errors later
            if (ferr > 0.9_dbl_kind*ferrmax) then         ! condition (5)

               converged = .false.

               ! reduce conductivity for next iteration
               do k = 1, nilyr
                  if (reduce_kh(k) .and. dqmat(k) > c0) then
                     frac = max(0.5*(c1-ferr/abs(fcondtopn-fcondbot)),p1)
!                     frac = p1
                     kh(k+nslyr+1) = kh(k+nslyr+1) * frac
                     kh(k+nslyr)   = kh(k+nslyr+1)
                  endif
               enddo

            endif               ! ferr 

         endif ! convergence

      enddo                     ! temperature iteration niter

      !-----------------------------------------------------------------
      ! Check for convergence failures.
      !-----------------------------------------------------------------
      if (.not.converged) then
         write(warning,*) 'Thermo iteration does not converge,'
         call add_warning(warning)
         write(warning,*) 'Ice thickness:',  hilyr*nilyr
         call add_warning(warning)
         write(warning,*) 'Snow thickness:', hslyr*nslyr
         call add_warning(warning)
         write(warning,*) 'dTsf, Tsf_errmax:',dTsf_prev, &
              Tsf_errmax
         call add_warning(warning)
         write(warning,*) 'Tsf:', Tsf
         call add_warning(warning)
         write(warning,*) 'fsurf:', fsurfn
         call add_warning(warning)
         write(warning,*) 'fcondtop, fcondbot, fswint', &
              fcondtopn, fcondbot, fswint
         call add_warning(warning)
         write(warning,*) 'fswsfc', fswsfc
         call add_warning(warning)
         write(warning,*) 'Iswabs',(Iswabs(k),k=1,nilyr)
         call add_warning(warning)
         write(warning,*) 'Flux conservation error =', ferr
         call add_warning(warning)
         write(warning,*) 'Initial snow temperatures:'
         call add_warning(warning)
         write(warning,*) (Tsn_init(k),k=1,nslyr)
         call add_warning(warning)
         write(warning,*) 'Initial ice temperatures:'
         call add_warning(warning)
         write(warning,*) (Tin_init(k),k=1,nilyr)
         call add_warning(warning)
         write(warning,*) 'Matrix ice temperature diff:'
         call add_warning(warning)
         write(warning,*) (dTmat(k),k=1,nilyr)
         call add_warning(warning)
         write(warning,*) 'dqmat*hilyr/dt:'
         call add_warning(warning)
         write(warning,*) (hilyr*dqmat(k)/dt,k=1,nilyr)
         call add_warning(warning)
         write(warning,*) 'Final snow temperatures:'
         call add_warning(warning)
         write(warning,*) (zTsn(k),k=1,nslyr)
         call add_warning(warning)
         write(warning,*) 'Matrix ice temperature diff:'
         call add_warning(warning)
         write(warning,*) (dTmat(k),k=1,nilyr)
         call add_warning(warning)
         write(warning,*) 'dqmat*hilyr/dt:'
         call add_warning(warning)
         write(warning,*) (hilyr*dqmat(k)/dt,k=1,nilyr)
         call add_warning(warning)
         write(warning,*) 'Final ice temperatures:'
         call add_warning(warning)
         write(warning,*) (zTin(k),k=1,nilyr)
         call add_warning(warning)
         write(warning,*) 'Ice melting temperatures:'
         call add_warning(warning)
         write(warning,*) (Tmlts(k),k=1,nilyr)
         call add_warning(warning)
         write(warning,*) 'Ice bottom temperature:', Tbot
         call add_warning(warning)
         write(warning,*) 'dT initial:'
         call add_warning(warning)
         write(warning,*) (Tmlts(k)-Tin_init(k),k=1,nilyr)
         call add_warning(warning)
         write(warning,*) 'dT final:'
         call add_warning(warning)
         write(warning,*) (Tmlts(k)-zTin(k),k=1,nilyr)
         call add_warning(warning)
         write(warning,*) 'zSin'
         call add_warning(warning)
         write(warning,*) (zSin(k),k=1,nilyr)
         call add_warning(warning)
         l_stop = .true.
         stop_label = "temperature_changes: Thermo iteration does not converge"
         return
      endif

      if (calc_Tsfc) then

         ! update fluxes that depend on Tsf
         flwoutn = flwoutn + dTsf_prev * dflwout_dT
         fsensn  = fsensn  + dTsf_prev * dfsens_dT
         flatn   = flatn   + dTsf_prev * dflat_dT

      endif                        ! calc_Tsfc

      end subroutine temperature_changes

!=======================================================================
!
! Compute thermal conductivity at interfaces (held fixed during
!  the subsequent iteration).
!
! NOTE: Ice conductivity must be >= kimin
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine conductivity (l_snow,                  &
                               nilyr,    nslyr,         &
                               hilyr,    hslyr,         &
                               zTin,     kh,       zSin)

      logical (kind=log_kind), intent(in) :: &
         l_snow          ! true if snow temperatures are computed

      integer (kind=int_kind), intent(in) :: & 
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         hilyr       , & ! ice layer thickness (same for all ice layers)
         hslyr           ! snow layer thickness (same for all snow layers)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         zTin         , & ! internal ice layer temperatures
         zSin             ! internal ice layer salinities

      real (kind=dbl_kind), dimension (nilyr+nslyr+1), &
         intent(out) :: &
         kh              ! effective conductivity at interfaces (W m-2 deg-1)

      ! local variables

      integer (kind=int_kind) :: &
         k               ! vertical index

      real (kind=dbl_kind), dimension (nilyr) :: &
         kilyr           ! thermal cond at ice layer midpoints (W m-1 deg-1)

      real (kind=dbl_kind), dimension (nslyr) :: &
         kslyr           ! thermal cond at snow layer midpoints (W m-1 deg-1)

      ! interior snow layers (simple for now, but may be fancier later)
      do k = 1, nslyr
         kslyr(k) = ksno
      enddo                     ! nslyr

      ! interior ice layers
      if (conduct == 'MU71') then
         ! Maykut and Untersteiner 1971 form (with Wettlaufer 1991 constants)
         do k = 1, nilyr
            kilyr(k) = kice + betak*zSin(k)/min(-puny,zTin(k))
            kilyr(k) = max (kilyr(k), kimin)
         enddo                     ! nilyr
      else
         ! Pringle et al JGR 2007 'bubbly brine'
         do k = 1, nilyr
            kilyr(k) = (2.11_dbl_kind - 0.011_dbl_kind*zTin(k) &
                     + 0.09_dbl_kind*zSin(k)/min(-puny,zTin(k))) &
                     * rhoi / 917._dbl_kind
            kilyr(k) = max (kilyr(k), kimin)
         enddo                     ! nilyr
      endif ! conductivity

      ! top snow interface, top and bottom ice interfaces
         ! top of snow layer; top surface of top ice layer
      if (l_snow) then
         kh(1)       = c2 * kslyr(1) / hslyr
         kh(1+nslyr) = c2 * kslyr(nslyr) * kilyr(1) / &
                       ( kslyr(nslyr)*hilyr +  &
                         kilyr(1    )*hslyr )
      else
         kh(1)       = c0
         kh(1+nslyr) = c2 * kilyr(1) / hilyr
      endif

      ! bottom surface of bottom ice layer
      kh(1+nslyr+nilyr) = c2 * kilyr(nilyr) / hilyr

      ! interior snow interfaces

      if (nslyr > 1) then
         do k = 2, nslyr
            if (l_snow) then
               kh(k) = c2 * kslyr(k-1) * kslyr(k) / &
                         ((kslyr(k-1) + kslyr(k))*hslyr)
            else
               kh(k) = c0
            endif
         enddo                     ! nilyr
      endif ! nslyr > 1

      ! interior ice interfaces
      do k = 2, nilyr
         kh(k+nslyr) = c2 * kilyr(k-1) * kilyr(k) / &
                         ((kilyr(k-1) + kilyr(k))*hilyr)
      enddo                     ! nilyr

      end subroutine conductivity

!=======================================================================
!
! Compute radiative and turbulent fluxes and their derivatives
! with respect to Tsf.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine surface_fluxes (Tsf,        fswsfc,            &
                                 rhoa,       flw,               &
                                 potT,       Qa,                &
                                 shcoef,     lhcoef,            &
                                 flwoutn,    fsensn,            &
                                 flatn,      fsurfn,            &
                                 dflwout_dT, dfsens_dT,         &
                                 dflat_dT,   dfsurf_dT)

      use ice_therm_shared, only: surface_heat_flux, dsurface_heat_flux_dTsf

      real (kind=dbl_kind), intent(in) :: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), intent(in) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

      real (kind=dbl_kind), &
         intent(inout) :: &
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsurfn          ! net flux to top surface, excluding fcondtopn

      real (kind=dbl_kind), &
         intent(inout) :: &
         dfsens_dT   , & ! deriv of fsens wrt Tsf (W m-2 deg-1)
         dflat_dT    , & ! deriv of flat wrt Tsf (W m-2 deg-1)
         dflwout_dT      ! deriv of flwout wrt Tsf (W m-2 deg-1)

      real (kind=dbl_kind), &
         intent(inout) :: &
         dfsurf_dT       ! derivative of fsurfn wrt Tsf

      ! surface heat flux
      call surface_heat_flux(Tsf,     fswsfc, &
                             rhoa,    flw,    &
                             potT,    Qa,     &
                             shcoef,  lhcoef, &
                             flwoutn, fsensn, &
                             flatn,   fsurfn)

      ! derivative of heat flux with respect to surface temperature
      call dsurface_heat_flux_dTsf(Tsf,       fswsfc,     &
                                   rhoa,      flw,        &
                                   potT,      Qa,         &
                                   shcoef,    lhcoef,     &
                                   dfsurf_dT, dflwout_dT, &
                                   dfsens_dT, dflat_dT)

      end subroutine surface_fluxes

!=======================================================================
!
! Compute terms in tridiagonal matrix that will be solved to find
!  the new vertical temperature profile
! This routine is for the case in which Tsfc is being computed.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! March 2004 by William H. Lipscomb for multiple snow layers
! April 2008 by E. C. Hunke, divided into two routines based on calc_Tsfc 

      subroutine get_matrix_elements_calc_Tsfc (nilyr, nslyr, &
                                      l_snow,   l_cold,           &
                                      Tsf,      Tbot,             &
                                      fsurfn,   dfsurf_dT,        &
                                      Tin_init, Tsn_init,         &
                                      kh,       Sswabs,           &
                                      Iswabs,                     &
                                      etai,     etas,             &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs)

      integer (kind=int_kind), intent(in) :: & 
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      logical (kind=log_kind), &
         intent(in) :: &
         l_snow      , & ! true if snow temperatures are computed
         l_cold          ! true if surface temperature is computed

      real (kind=dbl_kind), intent(in) :: &
         Tsf             ! ice/snow top surface temp (deg C)

      real (kind=dbl_kind), intent(in) :: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn (W/m^2)
         Tbot            ! ice bottom surface temperature (deg C)

      real (kind=dbl_kind), intent(in) :: &
         dfsurf_dT       ! derivative of fsurf wrt Tsf

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         etai        , & ! dt / (rho*cp*h) for ice layers
         Tin_init    , & ! ice temp at beginning of time step
         Sswabs      , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabs      , & ! absorbed SW flux in ice layers
         etas        , & ! dt / (rho*cp*h) for snow layers
         Tsn_init        ! snow temp at beginning of time step
                         ! Note: no absorbed SW in snow layers

      real (kind=dbl_kind), dimension (nslyr+nilyr+1), &
         intent(in) :: &
         kh              ! effective conductivity at layer interfaces

      real (kind=dbl_kind), dimension (nslyr+nilyr+1), &
         intent(inout) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      ! local variables

      integer (kind=int_kind) :: &
         k, ki, kr       ! vertical indices and row counters

      !-----------------------------------------------------------------
      ! Initialize matrix elements.
      ! Note: When we do not need to solve for the surface or snow
      !       temperature, we solve dummy equations with solution T = 0.
      !       Ice layers are fully initialized below.
      !-----------------------------------------------------------------

      do k = 1, nslyr+1
         sbdiag(k) = c0
         diag  (k) = c1
         spdiag(k) = c0
         rhs   (k) = c0
      enddo
            
      !-----------------------------------------------------------------
      ! Compute matrix elements
      !
      ! Four possible cases to solve:
      !   (1) Cold surface (Tsf < 0), snow present
      !   (2) Melting surface (Tsf = 0), snow present
      !   (3) Cold surface (Tsf < 0), no snow
      !   (4) Melting surface (Tsf = 0), no snow
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Tsf equation for case of cold surface (with or without snow)
      !-----------------------------------------------------------------

      if (l_cold) then
         if (l_snow) then
            k = 1
         else                ! no snow
            k = 1 + nslyr
         endif
         kr = k

         sbdiag(kr) = c0
         diag  (kr) = dfsurf_dT - kh(k)
         spdiag(kr) = kh(k)
         rhs   (kr) = dfsurf_dT*Tsf - fsurfn
      endif                  ! l_cold

      !-----------------------------------------------------------------
      ! top snow layer
      !-----------------------------------------------------------------
!           k = 1
!           kr = 2

      if (l_snow) then
         if (l_cold) then
            sbdiag(2) = -etas(1) * kh(1)
            spdiag(2) = -etas(1) * kh(2)
            diag  (2) = c1 &
                      + etas(1) * (kh(1) + kh(2))
            rhs   (2) = Tsn_init(1) &
                      + etas(1) * Sswabs(1)
         else                ! melting surface
            sbdiag(2) = c0
            spdiag(2) = -etas(1) * kh(2)
            diag  (2) = c1 &
                      + etas(1) * (kh(1) + kh(2))
            rhs   (2) = Tsn_init(1) &
                      + etas(1)*kh(1)*Tsf &
                      + etas(1) * Sswabs(1)
         endif               ! l_cold
      endif                  ! l_snow

      !-----------------------------------------------------------------
      ! remaining snow layers
      !-----------------------------------------------------------------

      if (nslyr > 1) then

         do k = 2, nslyr
            kr = k + 1

            if (l_snow) then
               sbdiag(kr) = -etas(k) * kh(k)
               spdiag(kr) = -etas(k) * kh(k+1)
               diag  (kr) = c1 &
                          + etas(k) * (kh(k) + kh(k+1))
               rhs   (kr) = Tsn_init(k) &
                          + etas(k) * Sswabs(k)
            endif
         enddo                  ! nslyr

      endif                     ! nslyr > 1

      if (nilyr > 1) then

      !-----------------------------------------------------------------
      ! top ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

         if (l_snow .or. l_cold) then
            sbdiag(kr) = -etai(ki) * kh(k)
            spdiag(kr) = -etai(ki) * kh(k+1)
            diag  (kr) = c1 &
                       + etai(ki) * (kh(k) + kh(k+1))
            rhs   (kr) = Tin_init(ki) &
                       + etai(ki)*Iswabs(ki)
         else    ! no snow, warm surface
            sbdiag(kr) = c0
            spdiag(kr) = -etai(ki) * kh(k+1)
            diag  (kr) = c1 &
                       + etai(ki) * (kh(k) + kh(k+1))
            rhs   (kr) = Tin_init(ki) &
                       + etai(ki)*Iswabs(ki) &
                       + etai(ki)*kh(k)*Tsf
         endif

      !-----------------------------------------------------------------
      ! bottom ice layer
      !-----------------------------------------------------------------

         ki = nilyr
         k  = ki + nslyr
         kr = k + 1
 
         sbdiag(kr) = -etai(ki) * kh(k)
         spdiag(kr) = c0
         diag  (kr) = c1  &
                    + etai(ki) * (kh(k) + kh(k+1))
         rhs   (kr) = Tin_init(ki) &
                    + etai(ki)*Iswabs(ki) &
                    + etai(ki)*kh(k+1)*Tbot
      
      else         ! nilyr = 1

      !-----------------------------------------------------------------
      ! single ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

         if (l_snow .or. l_cold) then
            sbdiag(kr) = -etai(ki) * kh(k)
            spdiag(kr) = c0
            diag  (kr) = c1 &
                       + etai(ki) * (kh(k) + kh(k+1))
            rhs   (kr) = Tin_init(ki) &
                       + etai(ki) * Iswabs(ki) &
                       + etai(ki) * kh(k+1)*Tbot
         else   ! no snow, warm surface
            sbdiag(kr) = c0
            spdiag(kr) = c0
            diag  (kr) = c1 &
                       + etai(ki) * (kh(k) + kh(k+1))
            rhs   (kr) = Tin_init(ki) &
                       + etai(ki) * Iswabs(ki) &
                       + etai(ki) * kh(k)*Tsf &
                       + etai(ki) * kh(k+1)*Tbot
         endif
 
      endif        ! nilyr > 1

      !-----------------------------------------------------------------
      ! interior ice layers
      !-----------------------------------------------------------------

      do ki = 2, nilyr-1
           
         k  = ki + nslyr
         kr = k + 1

            sbdiag(kr) = -etai(ki) * kh(k)
            spdiag(kr) = -etai(ki) * kh(k+1)
            diag  (kr) = c1 &
                       + etai(ki) * (kh(k) + kh(k+1))
            rhs   (kr) = Tin_init(ki) &
                       + etai(ki)*Iswabs(ki)
      enddo                     ! nilyr

      end subroutine get_matrix_elements_calc_Tsfc

!=======================================================================
!
! Compute terms in tridiagonal matrix that will be solved to find
!  the new vertical temperature profile
! This routine is for the case in which Tsfc is already known.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! March 2004 by William H. Lipscomb for multiple snow layers
! April 2008 by E. C. Hunke, divided into two routines based on calc_Tsfc 

      subroutine get_matrix_elements_know_Tsfc (nilyr, nslyr, &
                                      l_snow,   Tbot,             &
                                      Tin_init, Tsn_init,         &
                                      kh,       Sswabs,           &
                                      Iswabs,                     &
                                      etai,     etas,             &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs,              &
                                      fcondtopn)

      integer (kind=int_kind), intent(in) :: & 
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      logical (kind=log_kind), &
         intent(in) :: &
         l_snow          ! true if snow temperatures are computed

      real (kind=dbl_kind), intent(in) :: &
         Tbot            ! ice bottom surface temperature (deg C)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         etai        , & ! dt / (rho*cp*h) for ice layers
         Tin_init    , & ! ice temp at beginning of time step
         Sswabs      , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabs      , & ! absorbed SW flux in ice layers
         etas        , & ! dt / (rho*cp*h) for snow layers
         Tsn_init        ! snow temp at beginning of time step
                         ! Note: no absorbed SW in snow layers

      real (kind=dbl_kind), dimension (nslyr+nilyr+1), &
         intent(in) :: &
         kh              ! effective conductivity at layer interfaces

      real (kind=dbl_kind), dimension (nslyr+nilyr+1), &
         intent(inout) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      real (kind=dbl_kind), intent(in),  &
         optional :: &
         fcondtopn       ! conductive flux at top sfc, positive down (W/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         k, ki, kr       ! vertical indices and row counters

      !-----------------------------------------------------------------
      ! Initialize matrix elements.
      ! Note: When we do not need to solve for the surface or snow
      !       temperature, we solve dummy equations with solution T = 0.
      !       Ice layers are fully initialized below.
      !-----------------------------------------------------------------

      do k = 1, nslyr+1
         sbdiag(k) = c0
         diag  (k) = c1
         spdiag(k) = c0
         rhs   (k) = c0
      enddo
            
      !-----------------------------------------------------------------
      ! Compute matrix elements
      !
      ! Four possible cases to solve:
      !   (1) Cold surface (Tsf < 0), snow present
      !   (2) Melting surface (Tsf = 0), snow present
      !   (3) Cold surface (Tsf < 0), no snow
      !   (4) Melting surface (Tsf = 0), no snow
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! top snow layer
      !-----------------------------------------------------------------
!        k = 1
!        kr = 2

      if (l_snow) then
         sbdiag(2) = c0
         spdiag(2) = -etas(1) * kh(2)
         diag  (2) = c1 &
                   + etas(1) * kh(2)
         rhs   (2) = Tsn_init(1) &
                   + etas(1) * Sswabs(1) &
                   + etas(1) * fcondtopn
      endif   ! l_snow

      !-----------------------------------------------------------------
      ! remaining snow layers
      !-----------------------------------------------------------------

      if (nslyr > 1) then

         do k = 2, nslyr
            kr = k + 1

            if (l_snow) then
               sbdiag(kr) = -etas(k) * kh(k)
               spdiag(kr) = -etas(k) * kh(k+1)
               diag  (kr) = c1 &
                          + etas(k) * (kh(k) + kh(k+1))
               rhs   (kr) = Tsn_init(k) &
                          + etas(k) * Sswabs(k)
            endif

         enddo                  ! nslyr

      endif                     ! nslyr > 1

      if (nilyr > 1) then

      !-----------------------------------------------------------------
      ! top ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

         if (l_snow) then

            sbdiag(kr) = -etai(ki) * kh(k)
            spdiag(kr) = -etai(ki) * kh(k+1)
            diag  (kr) = c1 &
                       + etai(ki) * (kh(k) + kh(k+1))
            rhs   (kr) = Tin_init(ki) &
                       + etai(ki) * Iswabs(ki)
         else                  
            sbdiag(kr) = c0
            spdiag(kr) = -etai(ki) * kh(k+1)
            diag  (kr) = c1 &
                       + etai(ki) * kh(k+1)
            rhs   (kr) = Tin_init(ki) &
                       + etai(ki) * Iswabs(ki) &
                       + etai(ki) * fcondtopn
         endif  ! l_snow

      !-----------------------------------------------------------------
      ! bottom ice layer
      !-----------------------------------------------------------------

         ki = nilyr
         k  = ki + nslyr
         kr = k + 1
      
         sbdiag(kr) = -etai(ki) * kh(k)
         spdiag(kr) = c0
         diag  (kr) = c1  &
                    + etai(ki) * (kh(k) + kh(k+1))
         rhs   (kr) = Tin_init(ki) &
                    + etai(ki)*Iswabs(ki) &
                    + etai(ki)*kh(k+1)*Tbot
      
      else         ! nilyr = 1

      !-----------------------------------------------------------------
      ! single ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

         if (l_snow) then
            sbdiag(kr) = -etai(ki) * kh(k)
            spdiag(kr) = c0
            diag  (kr) = c1 &
                       + etai(ki) * (kh(k) + kh(k+1))
            rhs   (kr) = Tin_init(ki) &
                       + etai(ki) * Iswabs(ki) &
                       + etai(ki) * kh(k+1)*Tbot
         else
            sbdiag(kr) = c0
            spdiag(kr) = c0
            diag  (kr) = c1 &
                       + etai(ki) * kh(k+1)
            rhs   (kr) = Tin_init(ki) &
                       + etai(ki) * Iswabs(ki) &
                       + etai(ki) * fcondtopn &
                       + etai(ki) * kh(k+1)*Tbot
         endif

      endif        ! nilyr > 1

      !-----------------------------------------------------------------
      ! interior ice layers
      !-----------------------------------------------------------------

      do ki = 2, nilyr-1
           
         k  = ki + nslyr
         kr = k + 1

         sbdiag(kr) = -etai(ki) * kh(k)
         spdiag(kr) = -etai(ki) * kh(k+1)
         diag  (kr) = c1 &
                    + etai(ki) * (kh(k) + kh(k+1))
         rhs   (kr) = Tin_init(ki) &
                    + etai(ki)*Iswabs(ki)

      enddo                     ! nilyr

      end subroutine get_matrix_elements_know_Tsfc

!=======================================================================
!
! Tridiagonal matrix solver--used to solve the implicit vertical heat
! equation in ice and snow
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine tridiag_solver (nmat,     sbdiag,   &
                                 diag,     spdiag,   &
                                 rhs,      xout)

      integer (kind=int_kind), intent(in) :: &
         nmat            ! matrix dimension

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         xout            ! solution vector

      ! local variables

      integer (kind=int_kind) :: &
         k               ! row counter

      real (kind=dbl_kind) :: &
         wbeta           ! temporary matrix variable

      real (kind=dbl_kind), dimension(nmat) :: &
         wgamma          ! temporary matrix variable

      wbeta = diag(1)
      xout(1) = rhs(1) / wbeta

      do k = 2, nmat
         wgamma(k) = spdiag(k-1) / wbeta
         wbeta = diag(k) - sbdiag(k)*wgamma(k)
         xout(k) = (rhs(k) - sbdiag(k)*xout(k-1)) &
                    / wbeta
      enddo                     ! k

      do k = nmat-1, 1, -1
         xout(k) = xout(k) - wgamma(k+1)*xout(k+1)
      enddo                     ! k

      end subroutine tridiag_solver

!=======================================================================

      end module ice_therm_bl99

!=======================================================================
