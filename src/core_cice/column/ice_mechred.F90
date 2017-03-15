!  SVN:$Id: ice_mechred.F90 1175 2017-03-02 19:53:26Z akt $
!=======================================================================

! Driver for ice mechanical redistribution (ridging)
!
! See these references:
!
! Flato, G. M., and W. D. Hibler III, 1995: Ridging and strength
!  in modeling the thickness distribution of Arctic sea ice,
!  J. Geophys. Res., 100, 18,611-18,626.
!
! Hibler, W. D. III, 1980: Modeling a variable thickness sea ice
!  cover, Mon. Wea. Rev., 108, 1943-1973, 1980.
!
! Lipscomb, W. H., E. C. Hunke, W. Maslowski, and J. Jakacki, 2007: 
!  Improving ridging schemes for high-resolution sea ice models.
!  J. Geophys. Res. 112, C03S91, doi:10.1029/2005JC003355.
! 
! Rothrock, D. A., 1975: The energetics of the plastic deformation of
!  pack ice by ridging, J. Geophys. Res., 80, 4514-4519.
!
! Thorndike, A. S., D. A. Rothrock, G. A. Maykut, and R. Colony, 
!  1975: The thickness distribution of sea ice, J. Geophys. Res., 
!  80, 4501-4513. 
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: New options for participation and redistribution (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)

      module ice_mechred

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c2, c10, c20, c25, &
          p05, p15, p25, p333, p5, &
          puny, Lfresh, rhoi, rhos, rhow, gravit
      use ice_itd, only: column_sum, &
                         column_conservation_check
      use ice_warnings, only: add_warning

      implicit none
      save

      private
      public :: ridge_ice, asum_ridging, ridge_itd

      real (kind=dbl_kind), parameter :: & 
         Cs = p25         , & ! fraction of shear energy contrbtng to ridging 
         fsnowrdg = p5    , & ! snow fraction that survives in ridging 
         Gstar  = p15     , & ! max value of G(h) that participates 
                              ! (krdg_partic = 0) 
         astar  = p05     , & ! e-folding scale for G(h) participation 
!echmod         astar  = p1        , & ! e-folding scale for G(h) participation 
                              ! (krdg_partic = 1) 
         maxraft= c1      , & ! max value of hrmin - hi = max thickness 
                              ! of ice that rafts (m) 
         Hstar  = c25         ! determines mean thickness of ridged ice (m) 
                              ! (krdg_redist = 0) 
                              ! Flato & Hibler (1995) have Hstar = 100 

      logical (kind=log_kind), parameter :: &
         l_conservation_check = .false.  ! if true, check conservation
                                         ! (useful for debugging)

!=======================================================================

      contains

!=======================================================================

! Compute changes in the ice thickness distribution due to divergence
! and shear.
!
! author: William H. Lipscomb, LANL

      subroutine ridge_ice (dt,          ndtd,       &
                            ncat,        n_aero,     &
                            nilyr,       nslyr,      &
                            ntrcr,       hin_max,    &
                            rdg_conv,    rdg_shear,  &
                            aicen,       trcrn,      &
                            vicen,       vsnon,      &
                            aice0,                   &
                            trcr_depend, trcr_base,  &
                            n_trcr_strata,           &
                            nt_strata,   l_stop,     &
                            stop_label,              &
                            krdg_partic, krdg_redist,&
                            mu_rdg,                  &
                            dardg1dt,    dardg2dt,   &
                            dvirdgdt,    opening,    &
                            fpond,                   &
                            fresh,       fhocn,      &
                            tr_brine,    faero_ocn,  &
                            aparticn,    krdgn,      &
                            aredistn,    vredistn,   &
                            dardg1ndt,   dardg2ndt,  &
                            dvirdgndt,               &
                            araftn,      vraftn)

      use ice_colpkg_tracers, only: nt_qice, nt_qsno, nt_fbri, nt_sice

      integer (kind=int_kind), intent(in) :: &
         ndtd       , & ! number of dynamics subcycles
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nslyr , & ! number of snow layers
         n_aero, & ! number of aerosol tracers
         ntrcr     ! number of tracers in use

      real (kind=dbl_kind), intent(in) :: &
         mu_rdg , & ! gives e-folding scale of ridged ice (m^.5) 
         dt             ! time step

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max   ! category limits (m)

      real (kind=dbl_kind), intent(in) :: &
         rdg_conv   , & ! normalized energy dissipation due to convergence (1/s)
         rdg_shear      ! normalized energy dissipation due to shear (1/s)
 
      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen      , & ! concentration of ice
         vicen      , & ! volume per unit area of ice          (m)
         vsnon          ! volume per unit area of snow         (m)
 
      real (kind=dbl_kind), dimension (:,:), intent(inout) :: & 
         trcrn          ! ice tracers 

      real (kind=dbl_kind), intent(inout) :: & 
         aice0          ! concentration of open water

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      logical (kind=log_kind), intent(out) :: &
         l_stop         ! if true, abort on return

      character (len=*), intent(out) :: &
         stop_label   ! diagnostic information for abort

      integer (kind=int_kind), intent(in) :: &
         krdg_partic  , & ! selects participation function
         krdg_redist      ! selects redistribution function

      logical (kind=log_kind), intent(in) :: &
         tr_brine       ! if .true., brine height differs from ice thickness

      ! optional history fields
      real (kind=dbl_kind), intent(inout), optional :: &
         dardg1dt   , & ! rate of fractional area loss by ridging ice (1/s)
         dardg2dt   , & ! rate of fractional area gain by new ridges (1/s)
         dvirdgdt   , & ! rate of ice volume ridged (m/s)
         opening    , & ! rate of opening due to divergence/shear (1/s)
         fpond      , & ! fresh water flux to ponds (kg/m^2/s)
         fresh      , & ! fresh water flux to ocean (kg/m^2/s)
         fhocn          ! net heat flux to ocean (W/m^2)

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         dardg1ndt  , & ! rate of fractional area loss by ridging ice (1/s)
         dardg2ndt  , & ! rate of fractional area gain by new ridges (1/s)
         dvirdgndt  , & ! rate of ice volume ridged (m/s)
         aparticn   , & ! participation function
         krdgn      , & ! mean ridge thickness/thickness of ridging ice
         araftn     , & ! rafting ice area
         vraftn     , & ! rafting ice volume 
         aredistn   , & ! redistribution function: fraction of new ridge area
         vredistn       ! redistribution function: fraction of new ridge volume

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         faero_ocn      ! aerosol flux to ocean (kg/m^2/s)

      ! local variables

      real (kind=dbl_kind), dimension (ncat) :: &
         eicen          ! energy of melting for each ice layer (J/m^2)
 
      real (kind=dbl_kind), dimension (ncat) :: &
         esnon, &       ! energy of melting for each snow layer (J/m^2)
         vbrin, &       ! ice volume with defined by brine height (m)
         sicen          ! Bulk salt in h ice (ppt*m)

      real (kind=dbl_kind) :: &
         asum       , & ! sum of ice and open water area
         aksum      , & ! ratio of area removed to area ridged
         msnow_mlt  , & ! mass of snow added to ocean (kg m-2)
         esnow_mlt  , & ! energy needed to melt snow in ocean (J m-2)
         mpond      , & ! mass of pond added to ocean (kg m-2)
         closing_net, & ! net rate at which area is removed    (1/s)
                        ! (ridging ice area - area of new ridges) / dt
         divu_adv   , & ! divu as implied by transport scheme  (1/s)
         opning     , & ! rate of opening due to divergence/shear
                        ! opning is a local variable;
                        ! opening is the history diagnostic variable
         ardg1      , & ! fractional area loss by ridging ice
         ardg2      , & ! fractional area gain by new ridges
         virdg      , & ! ice volume ridged
         aopen          ! area opening due to divergence/shear

      real (kind=dbl_kind), dimension (n_aero) :: &
         maero          ! aerosol mass added to ocean (kg m-2)

      real (kind=dbl_kind), dimension (0:ncat) :: &
         apartic        ! participation function; fraction of ridging
                        ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (ncat) :: &
         hrmin      , & ! minimum ridge thickness
         hrmax      , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp      , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg       , & ! mean ridge thickness/thickness of ridging ice
         ardg1n     , & ! area of ice ridged
         ardg2n     , & ! area of new ridges
         virdgn     , & ! ridging ice volume
         mraftn         ! rafting ice mask 

      real (kind=dbl_kind) :: &
         vice_init, vice_final, & ! ice volume summed over categories
         vsno_init, vsno_final, & ! snow volume summed over categories
         eice_init, eice_final, & ! ice energy summed over layers
         vbri_init, vbri_final, & ! ice volume in fbri*vicen summed over categories
         sice_init ,sice_final, & ! ice bulk salinity summed over categories
         esno_init, esno_final    ! snow energy summed over layers

      integer (kind=int_kind), parameter :: &
         nitermax = 20  ! max number of ridging iterations

      integer (kind=int_kind) :: &
         n          , & ! thickness category index
         niter      , & ! iteration counter
         k          , & ! vertical index
         it             ! tracer index

      real (kind=dbl_kind) :: &
         dti            ! 1 / dt

      logical (kind=log_kind) :: &
         iterate_ridging ! if true, repeat the ridging

      character (len=char_len) :: &
         fieldid        ! field identifier

      character(len=char_len_long) :: &
         warning ! warning message

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.

      msnow_mlt = c0
      esnow_mlt = c0
      maero (:) = c0
      mpond     = c0
      ardg1     = c0
      ardg2     = c0
      virdg     = c0
      ardg1n(:) = c0
      ardg2n(:) = c0
      virdgn(:) = c0
      mraftn(:) = c0
      aopen     = c0

      !-----------------------------------------------------------------
      ! Compute area of ice plus open water before ridging.
      !-----------------------------------------------------------------

      call asum_ridging (ncat, aicen, aice0, asum)

      !-----------------------------------------------------------------
      ! Compute the area opening and closing.
      !-----------------------------------------------------------------

      call ridge_prep (dt,                      &
                       ncat,      hin_max,      &
                       rdg_conv,  rdg_shear,    &
                       asum,      closing_net,  &
                       divu_adv,  opning)

      !-----------------------------------------------------------------
      ! Compute initial values of conserved quantities. 
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         do n = 1, ncat
         eicen(n) = c0
         esnon(n) = c0
         sicen(n) = c0
         vbrin(n) = c0

         do k = 1, nilyr
            eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
            sicen(n) = sicen(n) + trcrn(nt_sice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo
         do k = 1, nslyr
            esnon(n) = esnon(n) + trcrn(nt_qsno+k-1,n) &
                     * vsnon(n)/real(nslyr,kind=dbl_kind)
         enddo
         vbrin(n) = vicen(n)
         if (tr_brine) vbrin(n) =  trcrn(nt_fbri,n) * vicen(n)
         enddo ! n

         call column_sum (ncat,                     &
                          vicen, vice_init)
         call column_sum (ncat,                     &
                          vsnon, vsno_init)
         call column_sum (ncat,                     &
                          eicen, eice_init)
         call column_sum (ncat,                     &
                          esnon, esno_init)
         call column_sum (ncat,                     &
                          sicen, sice_init)
         call column_sum (ncat,                     &
                          vbrin, vbri_init)

      endif            

      rdg_iteration: do niter = 1, nitermax

      !-----------------------------------------------------------------
      ! Compute the thickness distribution of ridging ice
      ! and various quantities associated with the new ridged ice.
      !-----------------------------------------------------------------

         call ridge_itd (ncat,        aice0,      &
                         aicen,       vicen,      &
                         krdg_partic, krdg_redist, &
                         mu_rdg,                   &
                         aksum,       apartic,    &
                         hrmin,       hrmax,      &
                         hrexp,       krdg,       &
                         aparticn,    krdgn,      &
                         mraftn)    

      !-----------------------------------------------------------------
      ! Redistribute area, volume, and energy.
      !-----------------------------------------------------------------

         call ridge_shift (ntrcr,       dt,          &
                           ncat,        hin_max,     &
                           aicen,       trcrn,       &
                           vicen,       vsnon,       &
                           aice0,       trcr_depend, &
                           trcr_base,   n_trcr_strata,&
                           nt_strata,   krdg_redist, &
                           aksum,       apartic,     &
                           hrmin,       hrmax,       &
                           hrexp,       krdg,        &
                           closing_net, opning,      &
                           ardg1,       ardg2,       &
                           virdg,       aopen,       &
                           ardg1n,      ardg2n,      &
                           virdgn,                   &
                           nslyr,       n_aero,      &
                           msnow_mlt,   esnow_mlt,   &
                           maero,       mpond,       &
                           l_stop,      stop_label,  &
                           aredistn,    vredistn)    

         if (l_stop) return

      !-----------------------------------------------------------------
      ! Make sure the new area = 1.  If not (because the closing
      ! and opening rates were reduced above), prepare to ridge again
      ! with new rates.
      !-----------------------------------------------------------------

         call asum_ridging (ncat, aicen, aice0, asum)

         if (abs(asum - c1) < puny) then
            iterate_ridging = .false.
            closing_net = c0
            opning      = c0
         else
            iterate_ridging = .true.
            divu_adv = (c1 - asum) / dt
            closing_net = max(c0, -divu_adv)
            opning = max(c0, divu_adv)
         endif

      !-----------------------------------------------------------------
      ! If done, exit.  If not, prepare to ridge again.
      !-----------------------------------------------------------------

         if (iterate_ridging) then
            write(warning,*) 'Repeat ridging, niter =', niter
            call add_warning(warning)
         else
            exit rdg_iteration
         endif

         if (niter == nitermax) then
            write(warning,*) ' '
            call add_warning(warning)
            write(warning,*) 'Exceeded max number of ridging iterations'
            call add_warning(warning)
            write(warning,*) 'max =',nitermax
            call add_warning(warning)
            l_stop = .true.
            stop_label = "ridge_ice: Exceeded max number of ridging iterations"
            return
         endif

      enddo rdg_iteration                    ! niter

      !-----------------------------------------------------------------
      ! Compute final values of conserved quantities. 
      ! Check for conservation (allowing for snow thrown into ocean).
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         do n = 1, ncat
         eicen(n) = c0
         esnon(n) = c0
         sicen(n) = c0
         vbrin(n) = c0

         do k = 1, nilyr
            eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
            sicen(n) = sicen(n) + trcrn(nt_sice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo
         do k = 1, nslyr
            esnon(n) = esnon(n) + trcrn(nt_qsno+k-1,n) &
                     * vsnon(n)/real(nslyr,kind=dbl_kind)
         enddo
         vbrin(n) =  vicen(n)
         if (tr_brine)  vbrin(n) =  trcrn(nt_fbri,n) * vbrin(n)
         enddo ! n

         call column_sum (ncat,                     &
                          vicen, vice_final)
         call column_sum (ncat,                     &
                          vsnon, vsno_final)
         call column_sum (ncat,                     &
                          eicen, eice_final)
         call column_sum (ncat,                     &
                          esnon, esno_final)
         call column_sum (ncat,                     &
                          sicen, sice_final)
         call column_sum (ncat,                     &
                          vbrin, vbri_final)

         vsno_final = vsno_final + msnow_mlt/rhos
         esno_final = esno_final + esnow_mlt

         fieldid = 'vice, ridging'
         call column_conservation_check (fieldid,               &
                                         vice_init, vice_final, &
                                         puny,                  &
                                         l_stop)
         fieldid = 'vsno, ridging'
         call column_conservation_check (fieldid,               &
                                         vsno_init, vsno_final, &
                                         puny,                  &
                                         l_stop)
         fieldid = 'eice, ridging'
         call column_conservation_check (fieldid,               &
                                         eice_init, eice_final, &
                                         puny*Lfresh*rhoi,      &
                                         l_stop)
         fieldid = 'esno, ridging'
         call column_conservation_check (fieldid,               &
                                         esno_init, esno_final, &
                                         puny*Lfresh*rhos,      &
                                         l_stop)
         fieldid = 'sice, ridging'
         call column_conservation_check (fieldid,               &
                                         sice_init, sice_final, &
                                         puny,                  &
                                         l_stop)
         fieldid = 'vbrin, ridging'
         call column_conservation_check (fieldid,               &
                                         vbri_init, vbri_final, &
                                         puny*c10,              &
                                         l_stop)
         if (l_stop) then
            stop_label = 'ridge_ice: Column conservation error'
            return
         endif

      endif                     ! l_conservation_check            

      !-----------------------------------------------------------------
      ! Compute ridging diagnostics.
      !-----------------------------------------------------------------

      dti = c1/dt

      if (present(dardg1dt)) then
         dardg1dt = ardg1*dti
      endif
      if (present(dardg2dt)) then
         dardg2dt = ardg2*dti
      endif
      if (present(dvirdgdt)) then
         dvirdgdt = virdg*dti
      endif
      if (present(opening)) then
         opening = aopen*dti
      endif
      if (present(dardg1ndt)) then
         do n = 1, ncat
            dardg1ndt(n) = ardg1n(n)*dti
         enddo
      endif
      if (present(dardg2ndt)) then
         do n = 1, ncat
            dardg2ndt(n) = ardg2n(n)*dti
         enddo
      endif
      if (present(dvirdgndt)) then
         do n = 1, ncat
            dvirdgndt(n) = virdgn(n)*dti
         enddo
      endif
      if (present(araftn)) then
         do n = 1, ncat
            araftn(n) = mraftn(n)*ardg2n(n)
!            araftn(n) = mraftn(n)*ardg1n(n)*p5
         enddo
      endif
      if (present(vraftn)) then
         do n = 1, ncat
            vraftn(n) = mraftn(n)*virdgn(n)
         enddo
      endif

      !-----------------------------------------------------------------
      ! Update fresh water and heat fluxes due to snow melt.
      !-----------------------------------------------------------------

      ! use thermodynamic time step (ndtd*dt here) to average properly
      dti = c1/(ndtd*dt)

      if (present(fresh)) then
         fresh = fresh + msnow_mlt*dti
      endif
      if (present(fhocn)) then
         fhocn = fhocn + esnow_mlt*dti
      endif
      if (present(faero_ocn)) then
         do it = 1, n_aero
            faero_ocn(it) = faero_ocn(it) + maero(it)*dti
         enddo
      endif
      if (present(fpond)) then
         fpond = fpond - mpond ! units change later
      endif

      !-----------------------------------------------------------------
      ! Check for fractional ice area > 1.
      !-----------------------------------------------------------------

      if (abs(asum - c1) > puny) then
         l_stop = .true.
         stop_label = "ridge_ice: total area > 1"

         write(warning,*) ' '
         call add_warning(warning)
         write(warning,*) 'Ridging error: total area > 1'
         call add_warning(warning)
         write(warning,*) 'area:', asum
         call add_warning(warning)
         write(warning,*) 'n, aicen:'
         call add_warning(warning)
         write(warning,*)  0, aice0
         call add_warning(warning)
         do n = 1, ncat
            write(warning,*) n, aicen(n)
            call add_warning(warning)
         enddo
         return
      endif

      end subroutine ridge_ice

!=======================================================================

! Find the total area of ice plus open water in each grid cell.
!
! This is similar to the aggregate_area subroutine except that the
! total area can be greater than 1, so the open water area is
! included in the sum instead of being computed as a residual.
!
! author: William H. Lipscomb, LANL

      subroutine asum_ridging (ncat, aicen, aice0, asum)

      integer (kind=int_kind), intent(in) :: & 
         ncat        ! number of thickness categories

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen          ! concentration of ice in each category

      real (kind=dbl_kind), intent(in) :: &
         aice0          ! concentration of open water

      real (kind=dbl_kind), intent(out):: &
         asum           ! sum of ice and open water area

      ! local variables

      integer (kind=int_kind) :: n

      asum = aice0
      do n = 1, ncat
            asum = asum + aicen(n)
      enddo

      end subroutine asum_ridging

!=======================================================================

! Initialize arrays, compute area of closing and opening
!
! author: William H. Lipscomb, LANL

      subroutine ridge_prep (dt,                          &
                             ncat,       hin_max,         &
                             rdg_conv,   rdg_shear,       &
                             asum,       closing_net,     &
                             divu_adv,   opning)

      integer (kind=int_kind), intent(in) :: & 
         ncat        ! number of thickness categories

      real (kind=dbl_kind), intent(in) :: &
         dt             ! time step (s)

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max   ! category limits (m)

      real (kind=dbl_kind), intent(in) :: &
         rdg_conv   , & ! normalized energy dissipation due to convergence (1/s)
         rdg_shear      ! normalized energy dissipation due to shear (1/s)

      real (kind=dbl_kind), intent(inout):: &
         asum           ! sum of ice and open water area

      real (kind=dbl_kind), &
         intent(out):: &
         closing_net, & ! net rate at which area is removed    (1/s)
         divu_adv   , & ! divu as implied by transport scheme  (1/s)
         opning         ! rate of opening due to divergence/shear

      ! local variables

      real (kind=dbl_kind), parameter :: &
         big = 1.0e+8_dbl_kind

      ! Set hin_max(ncat) to a big value to ensure that all ridged ice
      ! is thinner than hin_max(ncat).
      hin_max(ncat) = big

      !-----------------------------------------------------------------
      ! Compute the net rate of closing due to convergence
      ! and shear, based on Flato and Hibler (1995).
      !
      ! For the elliptical yield curve:
      !    rdg_conv  = -min (divu, 0)
      !    rdg_shear = (1/2) * (Delta - abs(divu))
      ! Note that the shear term also accounts for divergence.
      !
      ! The energy dissipation rate is equal to the net closing rate
      ! times the ice strength.
      !
      ! NOTE: The NET closing rate is equal to the rate that open water
      !  area is removed, plus the rate at which ice area is removed by
      !  ridging, minus the rate at which area is added in new ridges.
      !  The GROSS closing rate is equal to the first two terms (open
      !  water closing and thin ice ridging) without the third term
      !  (thick, newly ridged ice).
      !
      ! rdg_conv is calculated differently in EAP (update_ice_rdg) and 
      ! represents closing_net directly.  In that case, rdg_shear=0.
      !-----------------------------------------------------------------

      closing_net = Cs*rdg_shear + rdg_conv

      !-----------------------------------------------------------------
      ! Compute divu_adv, the divergence rate given by the transport/
      ! advection scheme, which may not be equal to divu as computed
      ! from the velocity field.
      !
      ! If divu_adv < 0, make sure the closing rate is large enough
      ! to give asum = 1.0 after ridging.
      !-----------------------------------------------------------------

      divu_adv = (c1-asum) / dt

      if (divu_adv < c0) closing_net = max(closing_net, -divu_adv)

      !-----------------------------------------------------------------
      ! Compute the (non-negative) opening rate that will give
      ! asum = 1.0 after ridging.
      !-----------------------------------------------------------------

      opning = closing_net + divu_adv

      end subroutine ridge_prep

!=======================================================================

! Compute the thickness distribution of the ice and open water
! participating in ridging and of the resulting ridges.
!
! This version includes new options for ridging participation and
!  redistribution.
! The new participation scheme (krdg_partic = 1) improves stability
!  by increasing the time scale for large changes in ice strength.
! The new exponential redistribution function (krdg_redist = 1) improves 
!  agreement between ITDs of modeled and observed ridges.   
!
! author: William H. Lipscomb, LANL
!
! 2006: Changed subroutine name to ridge_itd
!       Added new options for ridging participation and redistribution.  

      subroutine ridge_itd (ncat,        aice0,           &
                            aicen,       vicen,           &
                            krdg_partic, krdg_redist,     &
                            mu_rdg,                       &
                            aksum,       apartic,         &
                            hrmin,       hrmax,           &
                            hrexp,       krdg,            &
                            aparticn,    krdgn,           &
                            mraft)

      integer (kind=int_kind), intent(in) :: & 
         ncat        ! number of thickness categories

      real (kind=dbl_kind), intent(in) :: &
         mu_rdg , & ! gives e-folding scale of ridged ice (m^.5) 
         aice0       ! concentration of open water

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen   , & ! concentration of ice
         vicen       ! volume per unit area of ice (m)

      integer (kind=int_kind), intent(in) :: &
         krdg_partic  , & ! selects participation function
         krdg_redist      ! selects redistribution function

      real (kind=dbl_kind), intent(out):: &
         aksum       ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (0:ncat), intent(out) :: &
         apartic     ! participation function; fraction of ridging
                     ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         hrmin   , & ! minimum ridge thickness
         hrmax   , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp   , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg        ! mean ridge thickness/thickness of ridging ice

      ! diagnostic, category values
      real (kind=dbl_kind), dimension(:), intent(out), optional :: &
         aparticn, & ! participation function
         krdgn       ! mean ridge thickness/thickness of ridging ice

      real (kind=dbl_kind), dimension (:), intent(out), optional :: &
         mraft       ! rafting ice mask 

      ! local variables

      integer (kind=int_kind) :: &
         n           ! thickness category index

      real (kind=dbl_kind), parameter :: &
         Gstari   = c1/Gstar, &
         astari   = c1/astar

      real (kind=dbl_kind), dimension(-1:ncat) :: &
         Gsum        ! Gsum(n) = sum of areas in categories 0 to n

      real (kind=dbl_kind) :: &
         work        ! temporary work array

      real (kind=dbl_kind) :: &
         hi      , & ! ice thickness for each cat (m)
         hrmean  , & ! mean ridge thickness (m)
         xtmp        ! temporary variable

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      Gsum   (-1) = c0        ! by definition
!      Gsum   (0)  = c1         ! to avoid divzero below

      if (aice0 > puny) then
         Gsum(0) = aice0
      else
         Gsum(0) = Gsum(-1)
      endif
      apartic(0)  = c0

      do n = 1, ncat
         Gsum   (n) = c1    ! to avoid divzero below
         apartic(n) = c0
         hrmin  (n) = c0
         hrmax  (n) = c0
         hrexp  (n) = c0
         krdg   (n) = c1

      !-----------------------------------------------------------------
      ! Compute the thickness distribution of ice participating in ridging.
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! First compute the cumulative thickness distribution function Gsum,
      !  where Gsum(n) is the fractional area in categories 0 to n.
      ! Ignore categories with very small areas.
      !-----------------------------------------------------------------

         if (aicen(n) > puny) then
            Gsum(n) = Gsum(n-1) + aicen(n)
         else
            Gsum(n) = Gsum(n-1)
         endif
      enddo

      ! normalize

      work = c1 / Gsum(ncat)
      do n = 0, ncat
         Gsum(n) = Gsum(n) * work
      enddo

      !-----------------------------------------------------------------
      ! Compute the participation function apartic; this is analogous to
      ! a(h) = b(h)g(h) as defined in Thorndike et al. (1975).
      !
      !                area lost from category n due to ridging/closing
      !  apartic(n) = --------------------------------------------------
      !                  total area lost due to ridging/closing
      !
      !-----------------------------------------------------------------

      if (krdg_partic == 0) then  ! Thornike et al. 1975 formulation

      !-----------------------------------------------------------------
      ! Assume b(h) = (2/Gstar) * (1 - G(h)/Gstar).
      ! The expressions for apartic are found by integrating b(h)g(h) between
      ! the category boundaries.
      !-----------------------------------------------------------------

         do n = 0, ncat
            if (Gsum(n) < Gstar) then
               apartic(n) = Gstari*(Gsum(n  ) - Gsum(n-1)) * &
                      (c2 - Gstari*(Gsum(n-1) + Gsum(n  )))
            elseif (Gsum(n-1) < Gstar) then
               apartic(n) = Gstari*(Gstar - Gsum(n-1)) * &
                      (c2 - Gstari*(Gstar + Gsum(n-1)))
            endif
         enddo                  ! n

      elseif (krdg_partic==1) then   ! exponential dependence on G(h)

      !-----------------------------------------------------------------
      ! b(h) = exp(-G(h)/astar)
      ! apartic(n) = [exp(-G(n-1)/astar - exp(-G(n)/astar] / [1-exp(-1/astar)]. 
      ! The expression for apartic is found by integrating b(h)g(h)
      ! between the category boundaries.
      !-----------------------------------------------------------------

         ! precompute exponential terms using Gsum as work array
         xtmp = c1 / (c1 - exp(-astari))
         Gsum(-1) = exp(-Gsum(-1)*astari) * xtmp
         do n = 0, ncat
            Gsum(n) = exp(-Gsum(n)*astari) * xtmp
            apartic(n) = Gsum(n-1) - Gsum(n)
         enddo                  ! n

      endif                     ! krdg_partic

      !-----------------------------------------------------------------
      ! Compute variables related to ITD of ridged ice:
      ! 
      ! krdg   = mean ridge thickness / thickness of ridging ice
      ! hrmin  = min ridge thickness
      ! hrmax  = max ridge thickness (krdg_redist = 0)
      ! hrexp  = ridge e-folding scale (krdg_redist = 1)
      !----------------------------------------------------------------

      if (krdg_redist == 0) then  ! Hibler 1980 formulation

      !-----------------------------------------------------------------
      ! Assume ridged ice is uniformly distributed between hrmin and hrmax. 
      ! 
      ! This parameterization is a modified version of Hibler (1980). 
      ! In the original paper the min ridging thickness is hrmin = 2*hi, 
      !  and the max thickness is hrmax = 2*sqrt(hi*Hstar). 
      ! 
      ! Here the min thickness is hrmin = min(2*hi, hi+maxraft), 
      !  so thick ridging ice is not required to raft. 
      !
      !-----------------------------------------------------------------

         do n = 1, ncat
            if (aicen(n) > puny) then 
               hi = vicen(n) / aicen(n) 
               hrmin(n) = min(c2*hi, hi + maxraft) 
               hrmax(n) = c2*sqrt(Hstar*hi) 
               hrmax(n) = max(hrmax(n), hrmin(n)+puny) 
               hrmean = p5 * (hrmin(n) + hrmax(n)) 
               krdg(n) = hrmean / hi 

               ! diagnostic rafting mask not implemented
            endif 
         enddo                  ! n

      else               ! krdg_redist = 1; exponential redistribution
 
      !----------------------------------------------------------------- 
      ! The ridge ITD is a negative exponential: 
      ! 
      !  g(h) ~ exp[-(h-hrmin)/hrexp], h >= hrmin 
      ! 
      ! where hrmin is the minimum thickness of ridging ice and 
      ! hrexp is the e-folding thickness.
      ! 
      ! Here, assume as above that hrmin = min(2*hi, hi+maxraft).
      ! That is, the minimum ridge thickness results from rafting,
      !  unless the ice is thicker than maxraft.
      !
      ! Also, assume that hrexp = mu_rdg*sqrt(hi).
      ! The parameter mu_rdg is tuned to give e-folding scales mostly
      !  in the range 2-4 m as observed by upward-looking sonar.
      !
      ! Values of mu_rdg in the right column give ice strengths
      !  roughly equal to values of Hstar in the left column
      !  (within ~10 kN/m for typical ITDs):
      !
      !   Hstar     mu_rdg
      !
      !     25        3.0
      !     50        4.0
      !     75        5.0
      !    100        6.0
      !----------------------------------------------------------------- 

         do n = 1, ncat
            if (aicen(n) > puny) then
               hi = vicen(n) / aicen(n)
               hi = max(hi,puny)
               hrmin(n) = min(c2*hi, hi + maxraft)
               hrexp(n) = mu_rdg * sqrt(hi)
               krdg(n) = (hrmin(n) + hrexp(n)) / hi

   !echmod:  check computational efficiency
               ! diagnostic rafting mask
               if (present(mraft)) then
                  mraft(n) = max(c0, sign(c1, hi+maxraft-hrmin(n)))
                  xtmp = mraft(n)*((c2*hi+hrexp(n))/hi - krdg(n))
                  mraft(n) = max(c0, sign(c1, puny-abs(xtmp)))
               endif
            endif
         enddo

      endif                     ! krdg_redist

      !----------------------------------------------------------------
      ! Compute aksum = net ice area removed / total area participating.
      ! For instance, if a unit area of ice with h = 1 participates in
      !  ridging to form a ridge with a = 1/3 and h = 3, then
      !  aksum = 1 - 1/3 = 2/3.
      !---------------------------------------------------------------- 

      aksum = apartic(0) ! area participating = area removed

      do n = 1, ncat
         ! area participating > area removed
         aksum = aksum + apartic(n) * (c1 - c1/krdg(n)) 
      enddo

      ! diagnostics
      if (present(aparticn)) then
         do n = 1, ncat
            aparticn(n) = apartic(n)
         enddo
      endif
      if (present(krdgn)) then
         do n = 1, ncat
            krdgn(n) = krdg(n)
         enddo
      endif

      end subroutine ridge_itd

!=======================================================================

! Remove area, volume, and energy from each ridging category
! and add to thicker ice categories.
!
! Tracers:  Ridging conserves ice volume and therefore conserves volume
! tracers. It does not conserve ice area, and therefore a portion of area 
! tracers are lost (corresponding to the net closing).  Area tracers on 
! ice that participates in ridging are carried onto the resulting ridged
! ice (except the portion that are lost due to closing).  Therefore, 
! tracers must be decremented if they are lost to the ocean during ridging
! (e.g. snow, ponds) or if they are being carried only on the level ice 
! area. 
!
! author: William H. Lipscomb, LANL

      subroutine ridge_shift (ntrcr,       dt,              &
                              ncat,        hin_max,         &
                              aicen,       trcrn,           &
                              vicen,       vsnon,           &
                              aice0,       trcr_depend,     &   
                              trcr_base,   n_trcr_strata,   &
                              nt_strata,   krdg_redist,     &
                              aksum,       apartic,         &
                              hrmin,       hrmax,           &
                              hrexp,       krdg,            &
                              closing_net, opning,          &
                              ardg1,       ardg2,           &
                              virdg,       aopen,           &
                              ardg1nn,     ardg2nn,         &
                              virdgnn,                      &
                              nslyr,       n_aero,          &
                              msnow_mlt,   esnow_mlt,       &
                              maero,       mpond,           &
                              l_stop,      stop_label,      &
                              aredistn,    vredistn)

      use ice_colpkg_tracers, only: nt_qsno, nt_fbri, &
                             nt_alvl, nt_vlvl, nt_aero, tr_aero, &
                             nt_apnd, nt_hpnd, tr_pond_topo, &
                             colpkg_compute_tracers
                           
      integer (kind=int_kind), intent(in) :: & 
         ncat  , & ! number of thickness categories
         nslyr , & ! number of snow layers
         ntrcr , & ! number of tracers in use
         n_aero, & ! number of aerosol tracers
         krdg_redist      ! selects redistribution function

      real (kind=dbl_kind), intent(in) :: &
         dt             ! time step (s)

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension(0:ncat), intent(in) :: &
         hin_max   ! category limits (m)

      real (kind=dbl_kind), intent(inout) :: &
         aice0          ! concentration of open water

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen      , & ! concentration of ice
         vicen      , & ! volume per unit area of ice          (m)
         vsnon          ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn          ! ice tracers

      real (kind=dbl_kind), intent(in) :: &
         aksum          ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (0:ncat), intent(in) :: &
         apartic        ! participation function; fraction of ridging
                        ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         hrmin      , & ! minimum ridge thickness
         hrmax      , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp      , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg           ! mean ridge thickness/thickness of ridging ice

      real (kind=dbl_kind), intent(inout) :: &
         closing_net, & ! net rate at which area is removed    (1/s)
         opning     , & ! rate of opening due to divergence/shear (1/s)
         ardg1      , & ! fractional area loss by ridging ice
         ardg2      , & ! fractional area gain by new ridges
         virdg      , & ! ice volume ridged (m)
         aopen          ! area opened due to divergence/shear

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         ardg1nn    , & ! area of ice ridged
         ardg2nn    , & ! area of new ridges
         virdgnn        ! ridging ice volume

      real (kind=dbl_kind), intent(inout) :: &
         msnow_mlt  , & ! mass of snow added to ocean (kg m-2)
         esnow_mlt  , & ! energy needed to melt snow in ocean (J m-2)
         mpond          ! mass of pond added to ocean (kg m-2)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         maero          ! aerosol mass added to ocean (kg m-2)

      logical (kind=log_kind), intent(inout) :: &
         l_stop         ! if true, abort on return

      character (len=*), intent(out) :: &
         stop_label   ! diagnostic information for abort

      real (kind=dbl_kind), dimension (:), intent(inout), optional :: &
         aredistn   , & ! redistribution function: fraction of new ridge area
         vredistn       ! redistribution function: fraction of new ridge volume

      ! local variables

      integer (kind=int_kind) :: &
         n, nr      , & ! thickness category indices
         k          , & ! ice layer index
         it         , & ! tracer index
         ntr        , & ! tracer index
         itl            ! loop index

      real (kind=dbl_kind), dimension (ncat) :: &
         aicen_init , & ! ice area before ridging
         vicen_init , & ! ice volume before ridging
         vsnon_init     ! snow volume before ridging

      real (kind=dbl_kind), dimension(ntrcr,ncat) :: &
         atrcrn         ! aicen*trcrn

      real (kind=dbl_kind), dimension(3) :: &
         trfactor       ! base quantity on which tracers are carried

      real (kind=dbl_kind) :: &
         work       , & ! temporary variable
         closing_gross  ! rate at which area removed, not counting
                        ! area of new ridges

! ECH note:  the following arrays only need be defined on iridge cells
      real (kind=dbl_kind) :: &
         afrac      , & ! fraction of category area ridged
         ardg1n     , & ! area of ice ridged
         ardg2n     , & ! area of new ridges
         virdgn     , & ! ridging ice volume
         vsrdgn     , & ! ridging snow volume 
         dhr        , & ! hrmax - hrmin
         dhr2       , & ! hrmax^2 - hrmin^2
         farea      , & ! fraction of new ridge area going to nr
         fvol           ! fraction of new ridge volume going to nr

      real (kind=dbl_kind) :: &
         esrdgn         ! ridging snow energy

      real (kind=dbl_kind) :: &
         hi1        , & ! thickness of ridging ice
         hexp       , & ! ridge e-folding thickness
         hL, hR     , & ! left and right limits of integration
         expL, expR , & ! exponentials involving hL, hR
         tmpfac     , & ! factor by which opening/closing rates are cut
         wk1            ! work variable

      character(len=char_len_long) :: &
         warning ! warning message

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Save initial state variables
      !-----------------------------------------------------------------

         aicen_init(n) = aicen(n)
         vicen_init(n) = vicen(n)
         vsnon_init(n) = vsnon(n)

      !-----------------------------------------------------------------
      ! Define variables equal to aicen*trcrn, vicen*trcrn, vsnon*trcrn
      !-----------------------------------------------------------------

         do it = 1, ntrcr
            atrcrn(it,n) = trcrn(it,n)*(trcr_base(it,1) * aicen(n) &
                                      + trcr_base(it,2) * vicen(n) &
                                      + trcr_base(it,3) * vsnon(n))
            if (n_trcr_strata(it) > 0) then    ! additional tracer layers
               do itl = 1, n_trcr_strata(it)
                  ntr = nt_strata(it,itl)
                  atrcrn(it,n) = atrcrn(it,n) * trcrn(ntr,n)
               enddo
            endif
         enddo

      enddo ! ncat

      !-----------------------------------------------------------------
      ! Based on the ITD of ridging and ridged ice, convert the net
      !  closing rate to a gross closing rate.
      ! NOTE: 0 < aksum <= 1
      !-----------------------------------------------------------------

      closing_gross = closing_net / aksum

      !-----------------------------------------------------------------
      ! Reduce the closing rate if more than 100% of the open water
      ! would be removed.  Reduce the opening rate proportionately.
      !-----------------------------------------------------------------

      if (apartic(0) > c0) then
         wk1 = apartic(0) * closing_gross * dt
         if (wk1 > aice0) then
            tmpfac = aice0 / wk1
            closing_gross = closing_gross * tmpfac
            opning = opning * tmpfac
         endif
      endif

      !-----------------------------------------------------------------
      ! Reduce the closing rate if more than 100% of any ice category
      ! would be removed.  Reduce the opening rate proportionately.
      !-----------------------------------------------------------------
      do n = 1, ncat
         if (aicen(n) > puny .and. apartic(n) > c0) then
            wk1 = apartic(n) * closing_gross * dt
            if (wk1 > aicen(n)) then
               tmpfac = aicen(n) / wk1
               closing_gross = closing_gross * tmpfac
               opning = opning * tmpfac
            endif
         endif
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Compute change in open water area due to closing and opening.
      !-----------------------------------------------------------------

      aice0 = aice0 - apartic(0)*closing_gross*dt + opning*dt

      if (aice0 < -puny) then
         l_stop = .true.
         stop_label = 'Ridging error: aice0 < 0'
         write(warning,*) stop_label
         call add_warning(warning)
         write(warning,*) 'aice0:', aice0
         call add_warning(warning)
         return

      elseif (aice0 < c0) then    ! roundoff error
         aice0 = c0
      endif

      aopen = opning*dt  ! optional diagnostic

      !-----------------------------------------------------------------
      ! Compute the area, volume, and energy of ice ridging in each
      !  category, along with the area of the resulting ridge.
      !-----------------------------------------------------------------

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify grid cells with nonzero ridging
      !-----------------------------------------------------------------

         if (aicen_init(n) > puny .and. apartic(n) > c0 &
                                  .and. closing_gross > c0) then

      !-----------------------------------------------------------------
      ! Compute area of ridging ice (ardg1n) and of new ridge (ardg2n).
      ! Make sure ridging fraction <=1.  (Roundoff errors can give
      !  ardg1 slightly greater than aicen.) 
      !-----------------------------------------------------------------

            ardg1n = apartic(n)*closing_gross*dt

            if (ardg1n > aicen_init(n) + puny) then
               l_stop = .true.
               stop_label = 'Ridging error: ardg > aicen'
               write(warning,*) stop_label
               call add_warning(warning)
               write(warning,*) 'n, ardg, aicen:', &
                    n, ardg1n, aicen_init(n)
               call add_warning(warning)
               return
            else
               ardg1n = min(aicen_init(n), ardg1n)
            endif

            ardg2n = ardg1n / krdg(n)
            afrac = ardg1n / aicen_init(n)

      !-----------------------------------------------------------------
      ! Subtract area, volume, and energy from ridging category n.
      ! Note: Tracer values are unchanged.
      !-----------------------------------------------------------------

            virdgn = vicen_init(n) * afrac
            vsrdgn = vsnon_init(n) * afrac

            aicen(n) = aicen(n) - ardg1n
            vicen(n) = vicen(n) - virdgn           
            vsnon(n) = vsnon(n) - vsrdgn

      !-----------------------------------------------------------------
      ! Increment ridging diagnostics
      !-----------------------------------------------------------------

            ardg1 = ardg1 + ardg1n
            ardg2 = ardg2 + ardg2n
            virdg = virdg + virdgn

            ardg1nn(n) = ardg1n
            ardg2nn(n) = ardg2n
            virdgnn(n) = virdgn

      !-----------------------------------------------------------------
      !  Place part of the snow and tracer lost by ridging into the ocean.
      !-----------------------------------------------------------------

            msnow_mlt = msnow_mlt + rhos*vsrdgn*(c1-fsnowrdg)

            if (tr_aero) then
               do it = 1, n_aero
                  maero(it) = maero(it) &
                            + vsrdgn*(c1-fsnowrdg) &
                            *(trcrn(nt_aero  +4*(it-1),n)   &
                            + trcrn(nt_aero+1+4*(it-1),n))
               enddo
            endif

            if (tr_pond_topo) then
               mpond = mpond + ardg1n * trcrn(nt_apnd,n) &
                                      * trcrn(nt_hpnd,n)
            endif

      !-----------------------------------------------------------------
      ! Compute quantities used to apportion ice among categories
      ! in the nr loop below
      !-----------------------------------------------------------------

            dhr  = hrmax(n) - hrmin(n)
            dhr2 = hrmax(n) * hrmax(n) - hrmin(n) * hrmin(n)

      !-----------------------------------------------------------------
      ! Increment energy needed to melt snow in ocean.
      ! Note that esnow_mlt < 0; the ocean must cool to melt snow.
      !-----------------------------------------------------------------

            do k = 1, nslyr
               esrdgn = vsrdgn * trcrn(nt_qsno+k-1,n) &
                               / real(nslyr,kind=dbl_kind)
               esnow_mlt = esnow_mlt + esrdgn*(c1-fsnowrdg)
            enddo

      !-----------------------------------------------------------------
      ! Subtract area- and volume-weighted tracers from category n.
      !-----------------------------------------------------------------

            do it = 1, ntrcr

               trfactor(1) = trcr_base(it,1)*ardg1n
               trfactor(2) = trcr_base(it,2)*virdgn
               trfactor(3) = trcr_base(it,3)*vsrdgn

               work = c0
               do k = 1, 3
                  work = work + trfactor(k)*trcrn(it,n)
               enddo
               if (n_trcr_strata(it) > 0) then    ! additional tracer layers
                  do itl = 1, n_trcr_strata(it)
                     ntr = nt_strata(it,itl)
                     work = work * trcrn(ntr,n)
                  enddo
               endif
               atrcrn(it,n) = atrcrn(it,n) - work

            enddo                  ! ntrcr

      !-----------------------------------------------------------------
      ! Add area, volume, and energy of new ridge to each category nr.
      !-----------------------------------------------------------------

            do nr = 1, ncat

               if (krdg_redist == 0) then ! Hibler 1980 formulation

      !-----------------------------------------------------------------
      ! Compute the fraction of ridged ice area and volume going to
      !  thickness category nr.
      !-----------------------------------------------------------------

                  if (hrmin(n) >= hin_max(nr) .or. &
                      hrmax(n) <= hin_max(nr-1)) then
                     hL = c0
                     hR = c0
                  else
                     hL = max (hrmin(n), hin_max(nr-1))
                     hR = min (hrmax(n), hin_max(nr))
                  endif

                  farea = (hR-hL) / dhr
                  fvol  = (hR*hR - hL*hL) / dhr2

               else         ! krdg_redist = 1; 2005 exponential formulation

      !-----------------------------------------------------------------
      ! Compute the fraction of ridged ice area and volume going to
      !  thickness category nr.
      !-----------------------------------------------------------------

                  if (nr < ncat) then

                     hi1  = hrmin(n)
                     hexp = hrexp(n)

                     if (hi1 >= hin_max(nr)) then
                        farea = c0
                        fvol  = c0
                     else
                        hL = max (hi1, hin_max(nr-1))
                        hR = hin_max(nr)
                        expL = exp(-(hL-hi1)/hexp)
                        expR = exp(-(hR-hi1)/hexp)
                        farea = expL - expR
                        fvol  = ((hL + hexp)*expL  &
                                    - (hR + hexp)*expR) / (hi1 + hexp)
                     endif

                  else             ! nr = ncat

                     hi1  = hrmin(n)
                     hexp = hrexp(n)

                     hL = max (hi1, hin_max(nr-1))
                     expL = exp(-(hL-hi1)/hexp)
                     farea = expL
                     fvol  = (hL + hexp)*expL / (hi1 + hexp)

                  endif            ! nr < ncat

                  ! diagnostics
                  if (n ==1) then  ! only for thinnest ridging ice
                     if (present(aredistn)) then
                        aredistn(nr) = farea*ardg2n
                     endif
                     if (present(vredistn)) then
                        vredistn(nr) = fvol*virdgn
                     endif
                  endif

               endif               ! krdg_redist

      !-----------------------------------------------------------------
      ! Transfer ice area, ice volume, and snow volume to category nr.
      !-----------------------------------------------------------------

               aicen(nr) = aicen(nr) + farea*ardg2n
               vicen(nr) = vicen(nr) + fvol *virdgn
               vsnon(nr) = vsnon(nr) + fvol *vsrdgn*fsnowrdg

      !-----------------------------------------------------------------
      ! Transfer area-weighted and volume-weighted tracers to category nr.
      ! Note: The global sum aicen*trcrn of ice area tracers 
      !       (trcr_depend = 0) is not conserved by ridging.
      !       However, ridging conserves the global sum of volume
      !       tracers (trcr_depend = 1 or 2).
      ! Tracers associated with level ice, or that are otherwise lost
      ! from ridging ice, are not transferred.
      ! We assume that all pond water is lost from ridging ice.
      !-----------------------------------------------------------------

               do it = 1, ntrcr

                  if (it /= nt_alvl .and. it /= nt_vlvl) then
                     trfactor(1) = trcr_base(it,1)*ardg2n*farea
                     trfactor(2) = trcr_base(it,2)*virdgn*fvol
                     trfactor(3) = trcr_base(it,3)*vsrdgn*fvol*fsnowrdg
                  else
                     trfactor(1) = c0
                     trfactor(2) = c0
                     trfactor(3) = c0
                  endif

                  work = c0
                  do k = 1, 3
                     work = work + trfactor(k)*trcrn(it,n)
                  enddo
                  if (n_trcr_strata(it) > 0) then    ! additional tracer layers
                     do itl = 1, n_trcr_strata(it)
                        ntr = nt_strata(it,itl)
                        if (ntr == nt_fbri) then  ! brine fraction only
                           work = work * trcrn(ntr,n)
                        else
                           work = c0
                        endif
                     enddo
                  endif
                  atrcrn(it,nr) = atrcrn(it,nr) + work

               enddo               ! ntrcr

            enddo                  ! nr (new ridges)

         endif                     ! nonzero ridging

      enddo                        ! n (ridging categories)

      !-----------------------------------------------------------------
      ! Compute new tracers
      !-----------------------------------------------------------------

      do n = 1, ncat
            call colpkg_compute_tracers (ntrcr,       trcr_depend,   &
                                         atrcrn(:,n), aicen(n),      &
                                         vicen(n),    vsnon(n),      &
                                         trcr_base,   n_trcr_strata, &
                                         nt_strata,   trcrn(:,n))
      enddo

      end subroutine ridge_shift

!=======================================================================

      end module ice_mechred

!=======================================================================
