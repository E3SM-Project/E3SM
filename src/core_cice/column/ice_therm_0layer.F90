!  SVN:$Id: ice_therm_0layer.F90 1178 2017-03-08 19:24:07Z eclare $
!=========================================================================
!
! Update ice and snow internal temperatures
! using zero-layer thermodynamics
!
! authors: Alison McLaren, UK MetOffice
!          Elizabeth C. Hunke, LANL
!
! 2012: Split from ice_therm_vertical.F90

      module ice_therm_0layer

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, p5, puny, &
          kseaice, ksno
      use ice_therm_bl99, only: surface_fluxes
      use ice_warnings, only: add_warning

      implicit none

      private
      public :: zerolayer_temperature

!=======================================================================

      contains

!=======================================================================
!
! Compute new surface temperature using zero layer model of Semtner
! (1976).
!
! New temperatures are computed iteratively by solving a
! surface flux balance equation (i.e. net surface flux from atmos
! equals conductive flux from the top to the bottom surface).
!
! author:  Alison McLaren, Met Office
!         (but largely taken from temperature_changes)

      subroutine zerolayer_temperature(dt,                 & 
                                       nilyr,    nslyr,    &
                                       rhoa,     flw,      &
                                       potT,     Qa,       &
                                       shcoef,   lhcoef,   &
                                       fswsfc,             &
                                       hilyr,    hslyr,    &
                                       Tsf,      Tbot,     &
                                       fsensn,   flatn,    &
                                       flwoutn,  fsurfn,   &
                                       fcondtopn,fcondbot, &
                                       l_stop,   stop_label)

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), intent(in) :: &
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         Tbot        , & ! ice bottom surface temperature (deg C)
         fswsfc          ! SW absorbed at ice/snow surface (W m-2)

      real (kind=dbl_kind), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr           ! snow layer thickness (m)

      real (kind=dbl_kind), intent(inout):: &
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn       ! downward cond flux at top surface (W m-2)

      real (kind=dbl_kind), intent(out):: &
         fcondbot        ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), &
         intent(inout) :: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      character (len=*), intent(out) :: &
         stop_label   ! abort error message

      ! local variables

      logical (kind=log_kind), parameter :: &
         l_zerolayerchecks = .true.

      integer (kind=int_kind), parameter :: &
         nitermax = 50   ! max number of iterations in temperature solver

      real (kind=dbl_kind), parameter :: &
         Tsf_errmax = 5.e-4_dbl_kind ! max allowed error in Tsf
                                     ! recommend Tsf_errmax < 0.01 K

      integer (kind=int_kind) :: &
         niter           ! iteration counter in temperature solver

      real (kind=dbl_kind) :: &
         Tsf_start   , & ! Tsf at start of iteration
         dTsf        , & ! Tsf - Tsf_start
         dfsurf_dT       ! derivative of fsurfn wrt Tsf

      real (kind=dbl_kind) :: &
         dTsf_prev   , & ! dTsf from previous iteration
         dfsens_dT   , & ! deriv of fsens wrt Tsf (W m-2 deg-1)
         dflat_dT    , & ! deriv of flat wrt Tsf (W m-2 deg-1)
         dflwout_dT      ! deriv of flwout wrt Tsf (W m-2 deg-1)

      real (kind=dbl_kind) :: &
         kh          , & ! effective conductivity
         diag        , & ! diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix equation

      real (kind=dbl_kind) :: &
         heff        , & ! effective ice thickness (m)
                         ! ( hice + hsno*kseaice/ksnow)
         kratio      , & ! ratio of ice and snow conductivies
         avg_Tsf         ! = 1. if Tsf averaged w/Tsf_start, else = 0.

      logical (kind=log_kind) :: &
         converged      ! = true when local solution has converged

      character(len=char_len_long) :: &
         warning ! warning message
      
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      fcondbot = c0

      converged = .false.
      
      dTsf_prev = c0

      !-----------------------------------------------------------------
      ! Solve for new temperatures.
      ! Iterate until temperatures converge with minimal temperature
      ! change.
      !-----------------------------------------------------------------

      do niter = 1, nitermax

         if (.not. converged) then
            
      !-----------------------------------------------------------------
      ! Update radiative and turbulent fluxes and their derivatives
      ! with respect to Tsf.
      !-----------------------------------------------------------------
            
            call surface_fluxes (Tsf,        fswsfc,            &
                                 rhoa,       flw,               &
                                 potT,       Qa,                &
                                 shcoef,     lhcoef,            &
                                 flwoutn,    fsensn,            &
                                 flatn,      fsurfn,            &
                                 dflwout_dT, dfsens_dT,         &
                                 dflat_dT,   dfsurf_dT)

      !-----------------------------------------------------------------
      ! Compute effective ice thickness (includes snow) and thermal 
      ! conductivity 
      !-----------------------------------------------------------------

            kratio = kseaice/ksno
   
            heff = hilyr + kratio * hslyr
            kh = kseaice / heff

      !-----------------------------------------------------------------
      ! Compute conductive flux at top surface, fcondtopn.
      ! If fsurfn < fcondtopn and Tsf = 0, then reset Tsf to slightly less
      !  than zero (but not less than -puny).
      !-----------------------------------------------------------------

            fcondtopn = kh * (Tsf - Tbot)
            
            if (fsurfn < fcondtopn) &
                 Tsf = min (Tsf, -puny)
            
      !-----------------------------------------------------------------
      ! Save surface temperature at start of iteration
      !-----------------------------------------------------------------

            Tsf_start = Tsf

      !-----------------------------------------------------------------
      ! Solve surface balance equation to obtain the new temperatures.
      !-----------------------------------------------------------------

            diag = dfsurf_dT - kh
            rhs  = dfsurf_dT*Tsf - fsurfn   &
                 - kh*Tbot
            Tsf = rhs / diag

      !-----------------------------------------------------------------
      ! Determine whether the computation has converged to an acceptable
      ! solution.  Four conditions must be satisfied:
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
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Initialize convergence flag (true until proven false), dTsf,
      !  and temperature-averaging coefficients.
      ! Average only if test 1 or 2 fails.
      ! Initialize energy.
      !-----------------------------------------------------------------

            converged = .true.
            dTsf = Tsf - Tsf_start
            avg_Tsf = c0

      !-----------------------------------------------------------------
      ! Condition 1: check for Tsf > 0
      ! If Tsf > 0, set Tsf = 0 and leave converged=.true.
      !-----------------------------------------------------------------

            if (Tsf > puny) then
               Tsf = c0
               dTsf = -Tsf_start

      !-----------------------------------------------------------------
      ! Condition 2: check for oscillating Tsf
      ! If oscillating, average all temps to increase rate of convergence.
      ! It is possible that this may never occur.
      !-----------------------------------------------------------------

            elseif (niter > 1 &                ! condition (2)
              .and. Tsf_start <= -puny &
              .and. abs(dTsf) > puny &
              .and. abs(dTsf_prev) > puny &
              .and. -dTsf/(dTsf_prev+puny*puny) > p5) then

               avg_Tsf  = c1  ! average with starting temp  
               dTsf = p5 * dTsf
               converged = .false.
            endif

      !-----------------------------------------------------------------
      ! If condition 2 failed, average new surface temperature with
      !  starting value.
      !-----------------------------------------------------------------
            Tsf = Tsf &
                + avg_Tsf * p5 * (Tsf_start - Tsf)

      !-----------------------------------------------------------------
      ! Condition 3: check for large change in Tsf
      !-----------------------------------------------------------------

            if (abs(dTsf) > Tsf_errmax) then
               converged = .false.
            endif

      !-----------------------------------------------------------------
      ! Condition 4: check for fsurfn < fcondtopn with Tsf > 0
      !-----------------------------------------------------------------

            fsurfn = fsurfn + dTsf*dfsurf_dT
            fcondtopn = kh * (Tsf-Tbot)

            if (Tsf > -puny .and. fsurfn < fcondtopn) then
               converged = .false.
            endif

            fcondbot = fcondtopn

            dTsf_prev = dTsf

         endif ! converged

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
         write(warning,*) 'fsurfn:', fsurfn
         call add_warning(warning)
         write(warning,*) 'fcondtopn, fcondbot', &
                          fcondtopn, fcondbot
         call add_warning(warning)
         l_stop = .true.
         stop_label = "zerolayer_temperature: Thermo iteration does not converge"
         return
      endif

      !-----------------------------------------------------------------
      ! Check that if Tsfc < 0, then fcondtopn = fsurfn
      !-----------------------------------------------------------------

      if (l_zerolayerchecks) then
         if (Tsf < c0 .and. & 
              abs(fcondtopn-fsurfn) > puny) then

            write(warning,*) 'fcondtopn does not equal fsurfn,'
            call add_warning(warning)
            write(warning,*) 'Tsf=',Tsf
            call add_warning(warning)
            write(warning,*) 'fcondtopn=',fcondtopn
            call add_warning(warning)
            write(warning,*) 'fsurfn=',fsurfn
            call add_warning(warning)
            l_stop = .true.
            stop_label = "zerolayer_temperature: fcondtopn /= fsurfn"
            return
         endif
      endif                     ! l_zerolayerchecks

      ! update fluxes that depend on Tsf
      flwoutn = flwoutn + dTsf_prev * dflwout_dT
      fsensn  = fsensn  + dTsf_prev * dfsens_dT
      flatn   = flatn   + dTsf_prev * dflat_dT

      end subroutine zerolayer_temperature

!=======================================================================

      end module ice_therm_0layer

!=======================================================================
