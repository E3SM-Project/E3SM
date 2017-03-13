!  SVN:$Id: ice_therm_itd.F90 1178 2017-03-08 19:24:07Z eclare $
!=======================================================================
!
! Thermo calculations after call to coupler, related to ITD:
! ice thickness redistribution, lateral growth and melting.
!
! NOTE: The thermodynamic calculation is split in two for load balancing.
!       First ice_therm_vertical computes vertical growth rates and coupler
!       fluxes.  Then ice_therm_itd does thermodynamic calculations not
!       needed for coupling.
!       
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb  
! 2006: Streamlined for efficiency by Elizabeth Hunke
! 2014: Column package created by Elizabeth Hunke
!
      module ice_therm_itd

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c2, c3, c4, c6, c10, &
          p001, p1, p333, p5, p666, puny, bignum, &
          rhos, rhoi, Lfresh, ice_ref_salinity
      use ice_warnings, only: add_warning
      
      
      implicit none
      save
      
      private
      public :: linear_itd, add_new_ice, lateral_melt

      logical (kind=log_kind), parameter :: &
         l_conservation_check = .false.   ! if true, check conservation
                                          ! (useful for debugging)

!=======================================================================

      contains

!=======================================================================
!
! ITD scheme that shifts ice among categories
!
! See Lipscomb, W. H.  Remapping the thickness distribution in sea
!     ice models. 2001, J. Geophys. Res., Vol 106, 13989--14000.
!
! Using the thermodynamic "velocities", interpolate to find the
! velocities in thickness space at the category boundaries, and
! compute the new locations of the boundaries.  Then for each
! category, compute the thickness distribution function,  g(h),
! between hL and hR, the left and right boundaries of the category.
! Assume g(h) is a linear polynomial that satisfies two conditions:
!
! (1) The ice area implied by g(h) equals aicen(n).
! (2) The ice volume implied by g(h) equals aicen(n)*hicen(n).
!
! Given g(h), at each boundary compute the ice area and volume lying
! between the original and new boundary locations.  Transfer area
! and volume across each boundary in the appropriate direction, thus
! restoring the original boundaries.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine linear_itd (ncat,        hin_max,     &
                             nilyr,       nslyr,       &
                             ntrcr,       trcr_depend, & 
                             trcr_base,   n_trcr_strata,&
                             nt_strata,                &
                             aicen_init,  vicen_init,  & 
                             aicen,       trcrn,       & 
                             vicen,       vsnon,       & 
                             aice,        aice0,       & 
                             fpond,       l_stop,      &
                             stop_label)

      use ice_itd, only: aggregate_area, shift_ice, & 
                         column_sum, column_conservation_check
      use ice_colpkg_tracers, only: nt_qice, nt_qsno, nt_fbri, nt_sice, &
                             tr_pond_topo, nt_apnd, nt_hpnd, tr_brine
      use ice_therm_shared, only: hi_min

      integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of thickness categories
         nilyr   , & ! number of ice layers
         nslyr   , & ! number of snow layers
         ntrcr       ! number of tracers in use

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max      ! category boundaries (m)

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen_init, & ! initial ice concentration (before vertical thermo)
         vicen_init    ! initial ice volume               (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen  , & ! ice concentration
         vicen  , & ! volume per unit area of ice      (m)
         vsnon      ! volume per unit area of snow     (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), intent(inout) :: &
         aice  , & ! concentration of ice
         aice0 , & ! concentration of open water
         fpond     ! fresh water flux to ponds (kg/m^2/s)

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

     ! character (char_len), intent(out) :: stop_label
      character (len=*), intent(out) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         n, nd        , & ! category indices
         k                ! ice layer index

      real (kind=dbl_kind) :: &
         slope        , & ! rate of change of dhice with hice
         dh0          , & ! change in ice thickness at h = 0
         da0          , & ! area melting from category 1
         damax        , & ! max allowed reduction in category 1 area
         etamin, etamax,& ! left and right limits of integration
         x1           , & ! etamax - etamin
         x2           , & ! (etamax^2 - etamin^2) / 2
         x3           , & ! (etamax^3 - etamin^3) / 3
         wk1, wk2         ! temporary variables

      real (kind=dbl_kind), dimension(0:ncat) :: &
         hbnew            ! new boundary locations

      real (kind=dbl_kind), dimension(ncat) :: &
         g0           , & ! constant coefficient in g(h)
         g1           , & ! linear coefficient in g(h)
         hL           , & ! left end of range over which g(h) > 0
         hR               ! right end of range over which g(h) > 0

      real (kind=dbl_kind), dimension(ncat) :: &
         hicen        , & ! ice thickness for each cat     (m)
         hicen_init   , & ! initial ice thickness for each cat (pre-thermo)
         dhicen       , & ! thickness change for remapping (m)
         daice        , & ! ice area transferred across boundary
         dvice            ! ice volume transferred across boundary

      real (kind=dbl_kind), dimension (ncat) :: &
         eicen, &     ! energy of melting for each ice layer (J/m^2)
         esnon, &     ! energy of melting for each snow layer (J/m^2)
         vbrin, &     ! ice volume defined by brine height (m)
         sicen        ! Bulk salt in h ice (ppt*m)

      real (kind=dbl_kind) :: &
         vice_init, vice_final, & ! ice volume summed over categories
         vsno_init, vsno_final, & ! snow volume summed over categories
         eice_init, eice_final, & ! ice energy summed over categories
         esno_init, esno_final, & ! snow energy summed over categories
         sice_init, sice_final, & ! ice bulk salinity summed over categories
         vbri_init, vbri_final    ! briny ice volume summed over categories

      ! NOTE: Third index of donor, daice, dvice should be ncat-1,
      !       except that compilers would have trouble when ncat = 1 
      integer (kind=int_kind), dimension(ncat) :: &
         donor            ! donor category index

      logical (kind=log_kind) :: &
         remap_flag       ! remap ITD if remap_flag is true

      character (len=char_len) :: &
         fieldid           ! field identifier

      logical (kind=log_kind), parameter :: &
         print_diags = .false.    ! if true, prints when remap_flag=F

      character(len=char_len_long) :: &
         warning ! warning message
      
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.

      hin_max(ncat) = 999.9_dbl_kind ! arbitrary big number

      do n = 1, ncat
         donor(n) = 0
         daice(n) = c0
         dvice(n) = c0
      enddo

      !-----------------------------------------------------------------
      ! Compute volume and energy sums that linear remapping should
      !  conserve.
      !-----------------------------------------------------------------

      if (l_conservation_check) then

      do n = 1, ncat

         eicen(n) = c0
         esnon(n) = c0
         vbrin(n) = c0
         sicen(n) = c0

         do k = 1, nilyr
            eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo
         do k = 1, nslyr
            esnon(n) = esnon(n) + trcrn(nt_qsno+k-1,n) &
                     * vsnon(n)/real(nslyr,kind=dbl_kind)
         enddo

         if (tr_brine) then
            vbrin(n) = vbrin(n) + trcrn(nt_fbri,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         endif

         do k = 1, nilyr
            sicen(n) = sicen(n) + trcrn(nt_sice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo

      enddo ! n

      call column_sum (ncat, vicen, vice_init)
      call column_sum (ncat, vsnon, vsno_init)
      call column_sum (ncat, eicen, eice_init)
      call column_sum (ncat, esnon, esno_init)
      call column_sum (ncat, sicen, sice_init)
      call column_sum (ncat, vbrin, vbri_init)

      endif ! l_conservation_check

      !-----------------------------------------------------------------
      ! Initialize remapping flag.
      ! Remapping is done wherever remap_flag = .true.
      ! In rare cases the category boundaries may shift too far for the
      !  remapping algorithm to work, and remap_flag is set to .false.
      ! In these cases the simpler 'rebin' subroutine will shift ice
      !  between categories if needed.
      !-----------------------------------------------------------------

      remap_flag = .true.

      !-----------------------------------------------------------------
      ! Compute thickness change in each category.
      !-----------------------------------------------------------------

      do n = 1, ncat

         if (aicen_init(n) > puny) then
             hicen_init(n) = vicen_init(n) / aicen_init(n)
         else
             hicen_init(n) = c0
         endif               ! aicen_init > puny

         if (aicen (n) > puny) then
             hicen (n) = vicen(n) / aicen(n) 
             dhicen(n) = hicen(n) - hicen_init(n)
         else
             hicen (n) = c0
             dhicen(n) = c0
         endif               ! aicen > puny

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Compute new category boundaries, hbnew, based on changes in
      ! ice thickness from vertical thermodynamics.
      !-----------------------------------------------------------------

      hbnew(0) = hin_max(0)

      do n = 1, ncat-1

         if (hicen_init(n)   > puny .and. &
             hicen_init(n+1) > puny) then
             ! interpolate between adjacent category growth rates
             slope = (dhicen(n+1) - dhicen(n)) / &
                 (hicen_init(n+1) - hicen_init(n))
             hbnew(n) = hin_max(n) + dhicen(n) &
                      + slope * (hin_max(n) - hicen_init(n))
         elseif (hicen_init(n) > puny) then ! hicen_init(n+1)=0
             hbnew(n) = hin_max(n) + dhicen(n)
         elseif (hicen_init(n+1) > puny) then ! hicen_init(n)=0
             hbnew(n) = hin_max(n) + dhicen(n+1)
         else
             hbnew(n) = hin_max(n)
         endif

      !-----------------------------------------------------------------
      ! Check that each boundary lies between adjacent values of hicen.
      ! If not, set remap_flag = .false.
      ! Write diagnosis outputs if remap_flag was changed to false
      !-----------------------------------------------------------------

         if (aicen(n) > puny .and. hicen(n) >= hbnew(n)) then
            remap_flag = .false.

            if (print_diags) then
               write(warning,*) 'ITD: hicen(n) > hbnew(n)'
               call add_warning(warning)
               write(warning,*) 'cat ',n
               call add_warning(warning)
               write(warning,*) 'hicen(n) =', hicen(n)
               call add_warning(warning)
               write(warning,*) 'hbnew(n) =', hbnew(n)
               call add_warning(warning)
            endif

         elseif (aicen(n+1) > puny .and. hicen(n+1) <= hbnew(n)) then
            remap_flag = .false.

            if (print_diags) then
               write(warning,*) 'ITD: hicen(n+1) < hbnew(n)'
               call add_warning(warning)
               write(warning,*) 'cat ',n
               call add_warning(warning)
               write(warning,*) 'hicen(n+1) =', hicen(n+1)
               call add_warning(warning)
               write(warning,*) 'hbnew(n) =', hbnew(n)
               call add_warning(warning)
            endif
         endif

      !-----------------------------------------------------------------
      ! Check that hbnew(n) lies between hin_max(n-1) and hin_max(n+1).
      ! If not, set remap_flag = .false.
      ! (In principle we could allow this, but it would make the code
      ! more complicated.)
      ! Write diagnosis outputs if remap_flag was changed to false
      !-----------------------------------------------------------------

         if (hbnew(n) > hin_max(n+1)) then
            remap_flag = .false.

            if (print_diags) then
               write(warning,*) 'ITD hbnew(n) > hin_max(n+1)'
               call add_warning(warning)
               write(warning,*) 'cat ',n
               call add_warning(warning)
               write(warning,*) 'hbnew(n) =', hbnew(n)
               call add_warning(warning)
               write(warning,*) 'hin_max(n+1) =', hin_max(n+1)
               call add_warning(warning)
            endif
         endif

         if (hbnew(n) < hin_max(n-1)) then
            remap_flag = .false.

            if (print_diags) then
               write(warning,*) 'ITD: hbnew(n) < hin_max(n-1)'
               call add_warning(warning)
               write(warning,*) 'cat ',n
               call add_warning(warning)
               write(warning,*) 'hbnew(n) =', hbnew(n)
               call add_warning(warning)
               write(warning,*) 'hin_max(n-1) =', hin_max(n-1)
               call add_warning(warning)
            endif
         endif

      enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Compute hbnew(ncat)
      !-----------------------------------------------------------------

      if (aicen(ncat) > puny) then
         hbnew(ncat) = c3*hicen(ncat) - c2*hbnew(ncat-1)
      else
         hbnew(ncat) = hin_max(ncat)
      endif
      hbnew(ncat) = max(hbnew(ncat),hin_max(ncat-1))

      !-----------------------------------------------------------------
      ! Identify cells where the ITD is to be remapped
      !-----------------------------------------------------------------

      if (remap_flag) then

      !-----------------------------------------------------------------
      ! Compute g(h) for category 1 at start of time step
      ! (hicen = hicen_init)
      !-----------------------------------------------------------------

         call fit_line(aicen(1),   hicen_init(1), &
                       hbnew(0),   hin_max   (1), &
                       g0   (1),   g1        (1), &
                       hL   (1),   hR        (1))

      !-----------------------------------------------------------------
      ! Find area lost due to melting of thin (category 1) ice
      !-----------------------------------------------------------------

         if (aicen(1) > puny) then

            dh0 = dhicen(1)
            if (dh0 < c0) then   ! remove area from category 1
               dh0 = min(-dh0,hin_max(1))   ! dh0 --> |dh0|

      !-----------------------------------------------------------------
      ! Integrate g(1) from 0 to dh0 to estimate area melted
      !-----------------------------------------------------------------

               ! right integration limit (left limit = 0)
               etamax = min(dh0,hR(1)) - hL(1)

               if (etamax > c0) then
                  x1 = etamax
                  x2 = p5 * etamax*etamax
                  da0 = g1(1)*x2 + g0(1)*x1 ! ice area removed

               ! constrain new thickness <= hicen_init
                  damax = aicen(1) * (c1-hicen(1)/hicen_init(1)) ! damax > 0
                  da0 = min (da0, damax)

               ! remove area, conserving volume
                  hicen(1) = hicen(1) * aicen(1) / (aicen(1)-da0)
                  aicen(1) = aicen(1) - da0

                  if (tr_pond_topo) &
                     fpond = fpond - (da0 * trcrn(nt_apnd,1) & 
                                          * trcrn(nt_hpnd,1))

               endif            ! etamax > 0

            else                ! dh0 >= 0
               hbnew(0) = min(dh0,hin_max(1))  ! shift hbnew(0) to right
            endif

         endif                  ! aicen(1) > puny

      !-----------------------------------------------------------------
      ! Compute g(h) for each ice thickness category.
      !-----------------------------------------------------------------

         do n = 1, ncat

            call fit_line(aicen(n),   hicen(n), &
                          hbnew(n-1), hbnew(n), &
                          g0   (n),   g1   (n), &
                          hL   (n),   hR   (n))

      !-----------------------------------------------------------------
      ! Compute area and volume to be shifted across each boundary.
      !-----------------------------------------------------------------

            donor(n) = 0
            daice(n) = c0
            dvice(n) = c0
         enddo

         do n = 1, ncat-1

            if (hbnew(n) > hin_max(n)) then ! transfer from n to n+1

               ! left and right integration limits in eta space
               etamin = max(hin_max(n), hL(n)) - hL(n)
               etamax = min(hbnew(n),   hR(n)) - hL(n)
               donor(n) = n

            else             ! hbnew(n) <= hin_max(n); transfer from n+1 to n

               ! left and right integration limits in eta space
               etamin = c0
               etamax = min(hin_max(n), hR(n+1)) - hL(n+1)
               donor(n) = n+1

            endif            ! hbnew(n) > hin_max(n)

            if (etamax > etamin) then
               x1  = etamax - etamin
               wk1 = etamin*etamin
               wk2 = etamax*etamax
               x2  = p5 * (wk2 - wk1)
               wk1 = wk1*etamin
               wk2 = wk2*etamax
               x3  = p333 * (wk2 - wk1)
               nd  = donor(n)
               daice(n) = g1(nd)*x2 + g0(nd)*x1
               dvice(n) = g1(nd)*x3 + g0(nd)*x2 + daice(n)*hL(nd)
            endif               ! etamax > etamin

            ! If daice or dvice is very small, shift no ice.

            nd = donor(n)

            if (daice(n) < aicen(nd)*puny) then
               daice(n) = c0
               dvice(n) = c0
               donor(n) = 0
            endif 

            if (dvice(n) < vicen(nd)*puny) then
               daice(n) = c0
               dvice(n) = c0
               donor(n) = 0
            endif

            ! If daice is close to aicen or dvice is close to vicen,
            ! shift entire category

            if (daice(n) > aicen(nd)*(c1-puny)) then
               daice(n) = aicen(nd)
               dvice(n) = vicen(nd)
            endif

            if (dvice(n) > vicen(nd)*(c1-puny)) then
               daice(n) = aicen(nd)
               dvice(n) = vicen(nd)
            endif

         enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Shift ice between categories as necessary
      !-----------------------------------------------------------------

         ! maintain qsno negative definiteness
         do n = 1, ncat
            do k = nt_qsno, nt_qsno+nslyr-1
               trcrn(k,n) = trcrn(k,n) + rhos*Lfresh
            enddo
         enddo
 
         call shift_ice (ntrcr,    ncat,        &
                         trcr_depend,           &
                         trcr_base,             &
                         n_trcr_strata,         &
                         nt_strata,             &
                         aicen,    trcrn,       &
                         vicen,    vsnon,       &
                         hicen,    donor,       &
                         daice,    dvice,       &
                         l_stop,   stop_label)
         if (l_stop) return

         ! maintain qsno negative definiteness
         do n = 1, ncat
            do k = nt_qsno, nt_qsno+nslyr-1
               trcrn(k,n) = trcrn(k,n) - rhos*Lfresh
            enddo
         enddo

      !-----------------------------------------------------------------
      ! Make sure hice(1) >= minimum ice thickness hi_min.
      !-----------------------------------------------------------------

         if (hi_min > c0 .and. aicen(1) > puny .and. hicen(1) < hi_min) then

            da0 = aicen(1) * (c1 - hicen(1)/hi_min)
            aicen(1) = aicen(1) - da0
            hicen(1) = hi_min

            if (tr_pond_topo) &
               fpond = fpond - (da0 * trcrn(nt_apnd,1) & 
                                    * trcrn(nt_hpnd,1))
         endif

      endif ! remap_flag

      !-----------------------------------------------------------------
      ! Update fractional ice area in each grid cell.
      !-----------------------------------------------------------------

      call aggregate_area (ncat, aicen, aice, aice0)

      !-----------------------------------------------------------------
      ! Check volume and energy conservation.
      !-----------------------------------------------------------------

      if (l_conservation_check) then

      do n = 1, ncat

         eicen(n) = c0
         esnon(n) = c0
         vbrin(n) = c0
         sicen(n) = c0

         do k = 1, nilyr
            eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo
         do k = 1, nslyr
            esnon(n) = esnon(n) + trcrn(nt_qsno+k-1,n) &
                     * vsnon(n)/real(nslyr,kind=dbl_kind)
         enddo

         if (tr_brine) then
            vbrin(n) = vbrin(n) + trcrn(nt_fbri,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         endif

         do k = 1, nilyr
            sicen(n) = sicen(n) + trcrn(nt_sice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo

      enddo ! n

      call column_sum (ncat, vicen, vice_final)
      call column_sum (ncat, vsnon, vsno_final)
      call column_sum (ncat, eicen, eice_final)
      call column_sum (ncat, esnon, esno_final)
      call column_sum (ncat, sicen, sice_final)
      call column_sum (ncat, vbrin, vbri_final)

      fieldid = 'vice, ITD remap'
      call column_conservation_check (fieldid,               &
                                      vice_init, vice_final, &
                                      puny,                  &
                                      l_stop)
      fieldid = 'vsno, ITD remap'
      call column_conservation_check (fieldid,               &
                                      vsno_init, vsno_final, &
                                      puny,                  &
                                      l_stop)
      fieldid = 'eice, ITD remap'
      call column_conservation_check (fieldid,               &
                                      eice_init, eice_final, &
                                      puny*Lfresh*rhoi,      &
                                      l_stop)
      fieldid = 'esno, ITD remap'
      call column_conservation_check (fieldid,               &
                                      esno_init, esno_final, &
                                      puny*Lfresh*rhos,      &
                                      l_stop)
      fieldid = 'sicen, ITD remap'
      call column_conservation_check (fieldid,               &
                                      sice_init, sice_final, &
                                      puny,                  &
                                      l_stop)
      fieldid = 'vbrin, ITD remap'
      call column_conservation_check (fieldid,               &
                                      vbri_init, vbri_final, &
                                      puny*c10,              &
                                      l_stop)
         if (l_stop) then
            stop_label = 'linear_itd: Column conservation error'
            return
         endif

      endif                     ! conservation check

      end subroutine linear_itd

!=======================================================================
!
! Fit g(h) with a line, satisfying area and volume constraints.
! To reduce roundoff errors caused by large values of g0 and g1,
! we actually compute g(eta), where eta = h - hL, and hL is the
! left boundary.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine fit_line (aicen,    hice,            &
                           hbL,      hbR,             &
                           g0,       g1,              &
                           hL,       hR)

      real (kind=dbl_kind), intent(in) :: &
         aicen           ! concentration of ice

      real (kind=dbl_kind), intent(in) :: &
         hbL, hbR    , & ! left and right category boundaries
         hice            ! ice thickness

      real (kind=dbl_kind), intent(out):: &
         g0, g1      , & ! coefficients in linear equation for g(eta)
         hL          , & ! min value of range over which g(h) > 0
         hR              ! max value of range over which g(h) > 0

      ! local variables

      real  (kind=dbl_kind) :: &
         h13         , & ! hbL + 1/3 * (hbR - hbL)
         h23         , & ! hbL + 2/3 * (hbR - hbL)
         dhr         , & ! 1 / (hR - hL)
         wk1, wk2        ! temporary variables

      !-----------------------------------------------------------------
      ! Compute g0, g1, hL, and hR for each category to be remapped.
      !-----------------------------------------------------------------

         if (aicen > puny .and. hbR - hbL > puny) then

         ! Initialize hL and hR

            hL = hbL
            hR = hbR

         ! Change hL or hR if hicen(n) falls outside central third of range

            h13 = p333 * (c2*hL + hR)
            h23 = p333 * (hL + c2*hR)
            if (hice < h13) then
               hR = c3*hice - c2*hL
            elseif (hice > h23) then
               hL = c3*hice - c2*hR
            endif

         ! Compute coefficients of g(eta) = g0 + g1*eta

            dhr = c1 / (hR - hL)
            wk1 = c6 * aicen * dhr
            wk2 = (hice - hL) * dhr
            g0 = wk1 * (p666 - wk2)
            g1 = c2*dhr * wk1 * (wk2 - p5)

         else

            g0 = c0
            g1 = c0
            hL = c0
            hR = c0

         endif                  ! aicen > puny

      end subroutine fit_line

!=======================================================================
!
! Given some added new ice to the base of the existing ice, recalculate 
! vertical tracer so that new grid cells are all the same size. 
!
! author: A. K. Turner, LANL
!
      subroutine update_vertical_tracers(nilyr, trc, h1, h2, trc0)

      integer (kind=int_kind), intent(in) :: &
         nilyr ! number of ice layers

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
           trc ! vertical tracer

      real (kind=dbl_kind), intent(in) :: &
         h1, & ! old thickness
         h2, & ! new thickness
         trc0  ! tracer value of added ice on ice bottom
           
      ! local variables

      real(kind=dbl_kind), dimension(nilyr) :: trc2 ! updated tracer temporary

      ! vertical indices for old and new grid
      integer :: k1, k2

      real (kind=dbl_kind) :: &
         z1a, z1b, & ! upper, lower boundary of old cell/added new ice at bottom
         z2a, z2b, & ! upper, lower boundary of new cell
         overlap , & ! overlap between old and new cell
         rnilyr

        rnilyr = real(nilyr,dbl_kind)

        ! loop over new grid cells
        do k2 = 1, nilyr

           ! initialize new tracer
           trc2(k2) = c0

           ! calculate upper and lower boundary of new cell
           z2a = ((k2 - 1) * h2) / rnilyr
           z2b = (k2       * h2) / rnilyr

           ! loop over old grid cells
           do k1 = 1, nilyr

              ! calculate upper and lower boundary of old cell
              z1a = ((k1 - 1) * h1) / rnilyr
              z1b = (k1       * h1) / rnilyr
              
              ! calculate overlap between old and new cell
              overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)

              ! aggregate old grid cell contribution to new cell
              trc2(k2) = trc2(k2) + overlap * trc(k1)

           enddo ! k1

           ! calculate upper and lower boundary of added new ice at bottom
           z1a = h1
           z1b = h2
           
           ! calculate overlap between added ice and new cell
           overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)
           ! aggregate added ice contribution to new cell
           trc2(k2) = trc2(k2) + overlap * trc0
           ! renormalize new grid cell
           trc2(k2) = (rnilyr * trc2(k2)) / h2

        enddo ! k2

        ! update vertical tracer array with the adjusted tracer
        trc = trc2

      end subroutine update_vertical_tracers

!=======================================================================
!
! Given the fraction of ice melting laterally in each grid cell
!  (computed in subroutine frzmlt_bottom_lateral), melt ice.
!
! author: C. M. Bitz, UW
! 2003:   Modified by William H. Lipscomb and Elizabeth C. Hunke, LANL
!
      subroutine lateral_melt (dt,         ncat,       &
                               nilyr,      nslyr,      &
                               n_aero,     fpond,      &
                               fresh,      fsalt,      &
                               fhocn,      faero_ocn,  &
                               rside,      meltl,      &
                               aicen,      vicen,      &
                               vsnon,      trcrn,      &
                               fzsal,      flux_bio,   &
                               nbtrcr,     nblyr)

      use ice_colpkg_tracers, only: nt_qice, nt_qsno, nt_aero, tr_aero, &
                             tr_pond_topo, nt_apnd, nt_hpnd, bio_index
      use ice_colpkg_shared, only: z_tracers , hs_ssl, solve_zsal
      use ice_zbgc, only: lateral_melt_bgc               

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of thickness categories
         nilyr   , & ! number of ice layers
         nblyr   , & ! number of bio layers
         nslyr   , & ! number of snow layers
         n_aero  , & ! number of aerosol tracers
         nbtrcr      ! number of bio tracers

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume per unit area of ice          (m)
         vsnon       ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcrn     ! tracer array

      real (kind=dbl_kind), intent(in) :: &
         rside     ! fraction of ice that melts laterally

      real (kind=dbl_kind), intent(inout) :: &
         fpond     , & ! fresh water flux to ponds (kg/m^2/s)
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt     , & ! salt flux to ocean (kg/m^2/s)
         fhocn     , & ! net heat flux to ocean (W/m^2)
         meltl     , & ! lateral ice melt         (m/step-->cm/day)
         fzsal         ! salt flux from zsalinity (kg/m2/s)
  
      real (kind=dbl_kind), dimension(nbtrcr), &
         intent(inout) :: &
         flux_bio  ! biology tracer flux from layer bgc (mmol/m^2/s)  

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         faero_ocn     ! aerosol flux to ocean (kg/m^2/s)

      ! local variables

      integer (kind=int_kind) :: &
         n           , & ! thickness category index
         k               ! layer index

      real (kind=dbl_kind) :: &
         dfhocn  , & ! change in fhocn
         dfpond  , & ! change in fpond
         dfresh  , & ! change in fresh
         dfsalt  , & ! change in fsalt
         dvssl   , & ! snow surface layer volume
         dvint       ! snow interior layer

      real (kind=dbl_kind), dimension (ncat) :: &
         vicen_init   ! volume per unit area of ice (m)

      if (rside > c0) then ! grid cells with lateral melting.

         do n = 1, ncat

      !-----------------------------------------------------------------
      ! Melt the ice and increment fluxes.
      !-----------------------------------------------------------------

            ! fluxes to coupler
            ! dfresh > 0, dfsalt > 0, dfpond > 0

            dfresh = (rhos*vsnon(n) + rhoi*vicen(n)) * rside / dt
            dfsalt = rhoi*vicen(n)*ice_ref_salinity*p001 * rside / dt
            fresh  = fresh + dfresh
            fsalt  = fsalt + dfsalt

            if (tr_pond_topo) then
               dfpond = aicen(n) &
                      * trcrn(nt_apnd,n) & 
                      * trcrn(nt_hpnd,n) &
                      * rside
               fpond = fpond - dfpond
            endif

            ! history diagnostics
            meltl = meltl + vicen(n)*rside

            ! state variables
            vicen_init(n) = vicen(n)
            aicen(n) = aicen(n) * (c1 - rside)
            vicen(n) = vicen(n) * (c1 - rside)
            vsnon(n) = vsnon(n) * (c1 - rside)

            do k = 1, nilyr
               ! enthalpy tracers do not change (e/v constant)
               ! heat flux to coupler for ice melt (dfhocn < 0)
               dfhocn = trcrn(nt_qice+k-1,n)*rside / dt &
                      * vicen(n)/real(nilyr,kind=dbl_kind)
               fhocn  = fhocn + dfhocn
            enddo                  ! nilyr

            do k = 1, nslyr
               ! heat flux to coupler for snow melt (dfhocn < 0)
               dfhocn = trcrn(nt_qsno+k-1,n)*rside / dt &
                      * vsnon(n)/real(nslyr,kind=dbl_kind)
               fhocn  = fhocn + dfhocn
            enddo                  ! nslyr

            if (tr_aero) then
               do k = 1, n_aero
                  faero_ocn(k) = faero_ocn(k) + (vsnon(n) &
                               *(trcrn(nt_aero  +4*(k-1),n)   &
                               + trcrn(nt_aero+1+4*(k-1),n))  &
                                              +  vicen(n) &
                               *(trcrn(nt_aero+2+4*(k-1),n)   &
                               + trcrn(nt_aero+3+4*(k-1),n))) &
                               * rside / dt
               enddo ! k
            endif    ! tr_aero

      !-----------------------------------------------------------------
      ! Biogeochemistry
      !-----------------------------------------------------------------     

            if (z_tracers) then   ! snow tracers
               dvssl  = min(p5*vsnon(n), hs_ssl*aicen(n))       !snow surface layer
               dvint  = vsnon(n)- dvssl                         !snow interior
               do k = 1, nbtrcr
                  flux_bio(k) = flux_bio(k) &
                              + (trcrn(bio_index(k)+nblyr+1,n)*dvssl  &
                              +  trcrn(bio_index(k)+nblyr+2,n)*dvint) &
                              * rside / dt
               enddo
            endif

         enddo       ! n

         if (solve_zsal .or. z_tracers) &
            call lateral_melt_bgc(dt,                         &
                                  ncat,        nblyr,         &
                                  rside,       vicen_init,    &
                                  trcrn,       fzsal,         &
                                  flux_bio,    nbtrcr)

      endif          ! rside

      end subroutine lateral_melt

!=======================================================================
!
! Given the volume of new ice grown in open water, compute its area
! and thickness and add it to the appropriate category or categories.
!
! NOTE: Usually all the new ice is added to category 1.  An exception is
!       made if there is no open water or if the new ice is too thick
!       for category 1, in which case ice is distributed evenly over the
!       entire cell.  Subroutine rebin should be called in case the ice
!       thickness lies outside category bounds after new ice formation.
!
! When ice must be added to categories above category 1, the mushy
! formulation (ktherm=2) adds it only to the bottom of the ice.  When
! added to only category 1, all formulations combine the new ice and
! existing ice tracers as bulk quantities.
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!         Adrian Turner, LANL
!
      subroutine add_new_ice (ncat,      nilyr,    nblyr,  &
                              n_aero,    dt,         &
                              ntrcr,     nltrcr,            &
                              hin_max,   ktherm,     &
                              aicen,     trcrn,      &
                              vicen,     vsnon1,     &
                              aice0,     aice,       &
                              frzmlt,    frazil,     &
                              frz_onset, yday,       &
                              update_ocn_f,          &
                              fresh,     fsalt,      &
                              Tf,        sss,        &
                              salinz,    phi_init,   &
                              dSin0_frazil,          &
                              bgrid,      cgrid,      igrid,    &
                              nbtrcr,    flux_bio,   &
                              ocean_bio, fzsal,      &
                              frazil_diag,           &
                              l_stop,    stop_label)

      use ice_itd, only: column_sum, &
                         column_conservation_check 
      use ice_colpkg_tracers, only: nt_Tsfc, nt_iage, nt_FY, nt_sice, nt_qice, &
                             nt_alvl, nt_vlvl, nt_aero, nt_apnd, &
                             tr_pond_cesm, tr_pond_lvl, tr_pond_topo, &
                             tr_iage, tr_FY, tr_lvl, tr_aero, tr_brine
      use ice_colpkg_shared, only: solve_zsal, skl_bgc, initbio_frac, salt_loss, rhosi
      use ice_mushy_physics, only: liquidus_temperature_mush, enthalpy_mush
      use ice_therm_shared, only: hfrazilmin
      use ice_zbgc, only: add_new_ice_bgc
      use ice_zbgc_shared, only: bgc_tracer_type

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr , & ! number of bio layers
         ntrcr , & ! number of tracers
         nltrcr, & ! number of zbgc tracers
         n_aero, & ! number of aerosol tracers
         ktherm    ! type of thermodynamics (0 0-layer, 1 BL99, 2 mushy)

      real (kind=dbl_kind), dimension(0:ncat), intent(in) :: &
         hin_max      ! category boundaries (m)

      real (kind=dbl_kind), intent(in) :: &
         dt    , & ! time step (s)
         aice  , & ! total concentration of ice
         frzmlt, & ! freezing/melting potential (W/m^2)
         Tf    , & ! freezing temperature (C)
         sss   , & ! sea surface salinity (ppt)
         vsnon1    ! category 1 snow volume per ice area (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen , & ! concentration of ice
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature

      real (kind=dbl_kind), intent(inout) :: &
         aice0     , & ! concentration of open water
         frazil    , & ! frazil ice growth        (m/step-->cm/day)
         frazil_diag,& ! frazil ice growth diagnostic (m/step-->cm/day)
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt         ! salt flux to ocean (kg/m^2/s)

      real (kind=dbl_kind), intent(inout), optional :: &
         frz_onset ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(in), optional :: &
         yday      ! day of year

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         salinz     ! initial salinity profile

      real (kind=dbl_kind), intent(in) :: &
         phi_init     , & ! initial frazil liquid fraction
         dSin0_frazil     ! initial frazil bulk salinity reduction from sss

      logical (kind=log_kind), intent(in) :: &
         update_ocn_f ! if true, update fresh water and salt fluxes

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

     ! character (char_len), intent(out) :: stop_label
      character (len=*), intent(out) :: stop_label

      ! BGC
      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      integer (kind=int_kind), intent(in) :: &
         nbtrcr          ! number of biology tracers

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         flux_bio   ! tracer flux to ocean from biology (mmol/m^2/s) 
        
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         ocean_bio   ! ocean concentration of biological tracer

      ! zsalinity
      real (kind=dbl_kind),  intent(inout) :: &  
         fzsal      ! salt flux to ocean from zsalinity (kg/m^2/s)

      ! local variables

      integer (kind=int_kind) :: &
         i, j         , & ! horizontal indices
         n            , & ! ice category index
         k            , & ! ice layer index
         it               ! aerosol tracer index

      real (kind=dbl_kind) :: &
         ai0new       , & ! area of new ice added to cat 1
         vi0new       , & ! volume of new ice added to cat 1
         hsurp        , & ! thickness of new ice added to each cat
         fnew         , & ! heat flx to open water for new ice (W/m^2)
         hi0new       , & ! thickness of new ice
         hi0max       , & ! max allowed thickness of new ice
         vsurp        , & ! volume of new ice added to each cat
         vtmp         , & ! total volume of new and old ice
         area1        , & ! starting fractional area of existing ice
         alvl         , & ! starting level ice area
         rnilyr       , & ! real(nilyr)
         dfresh       , & ! change in fresh
         dfsalt       , & ! change in fsalt
         vi0tmp       , & ! frzmlt part of frazil
         Ti           , & ! frazil temperature
         qi0new       , & ! frazil ice enthalpy
         Si0new       , & ! frazil ice bulk salinity
         vi0_init     , & ! volume of new ice
         vice1        , & ! starting volume of existing ice
         vice_init, vice_final, & ! ice volume summed over categories
         eice_init, eice_final    ! ice energy summed over categories

      real (kind=dbl_kind), dimension (nilyr) :: &
         Sprofile         ! salinity profile used for new ice additions

      character (len=char_len) :: &
         fieldid           ! field identifier

      real (kind=dbl_kind), dimension (ncat) :: &
         eicen, &     ! energy of melting for each ice layer (J/m^2)
         aicen_init, &    ! fractional area of ice
         vicen_init       ! volume per unit area of ice (m)

      ! BGC
      real (kind=dbl_kind) :: &
         vbri1       , & ! starting volume of existing brine
         vbri_init   , & ! brine volume summed over categories
         vbri_final      ! brine volume summed over categories

      real (kind=dbl_kind), dimension (ncat) :: &
         vbrin           ! trcrn(nt_fbri,n)*vicen(n) 

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------

      l_stop = .false.

      rnilyr = real(nilyr,kind=dbl_kind)

      if (ncat > 1) then
         hi0max = hin_max(1)*0.9_dbl_kind  ! not too close to boundary
      else
         hi0max = bignum                   ! big number
      endif

      do n = 1, ncat
         aicen_init(n) = aicen(n)
         vicen_init(n) = vicen(n)
         eicen(n) = c0
      enddo

      if (l_conservation_check) then

         do n = 1, ncat
         do k = 1, nilyr
            eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo
         enddo
         call column_sum (ncat, vicen, vice_init)
         call column_sum (ncat, eicen, eice_init)

      endif ! l_conservation_check

      !-----------------------------------------------------------------
      ! Compute average enthalpy of new ice.
      ! Sprofile is the salinity profile used when adding new ice to
      ! all categories, for ktherm/=2, and to category 1 for all ktherm.
      !
      ! NOTE:  POP assumes new ice is fresh!
      !-----------------------------------------------------------------

      if (ktherm == 2) then  ! mushy
         if (sss > c2 * dSin0_frazil) then
            Si0new = sss - dSin0_frazil
         else
            Si0new = sss**2 / (c4*dSin0_frazil)
         endif
         do k = 1, nilyr
            Sprofile(k) = Si0new
         enddo
         Ti = min(liquidus_temperature_mush(Si0new/phi_init), -p1)
         qi0new = enthalpy_mush(Ti, Si0new)
      else
         do k = 1, nilyr
            Sprofile(k) = salinz(k)
         enddo
         qi0new = -rhoi*Lfresh
      endif    ! ktherm

      !-----------------------------------------------------------------
      ! Compute the volume, area, and thickness of new ice.
      !-----------------------------------------------------------------

      fnew = max (frzmlt, c0)    ! fnew > 0 iff frzmlt > 0
      vi0new = -fnew*dt / qi0new ! note sign convention, qi < 0
      vi0_init = vi0new          ! for bgc

      ! increment ice volume and energy
      if (l_conservation_check) then
         vice_init = vice_init + vi0new
         eice_init = eice_init + vi0new*qi0new
      endif

      ! history diagnostics
      frazil = vi0new
      if (solve_zsal) &  
             fzsal = fzsal  &
                  - rhosi*vi0new/dt*p001*sss*salt_loss

      if (present(frz_onset) .and. present(yday)) then
         if (frazil > puny .and. frz_onset < puny) frz_onset = yday
      endif

      !-----------------------------------------------------------------
      ! Update fresh water and salt fluxes.
      !
      ! NOTE: POP assumes fresh water and salt flux due to frzmlt > 0
      !       is NOT included in fluxes fresh and fsalt.
      !-----------------------------------------------------------------

      if (update_ocn_f) then
         if (ktherm <= 1) then
            dfresh = -rhoi*vi0new/dt 
            dfsalt = ice_ref_salinity*p001*dfresh
            fresh  = fresh  + dfresh
            fsalt  = fsalt  + dfsalt
         ! elseif (ktherm == 2) the fluxes are added elsewhere
         endif
      else ! update_ocn_f = false
         if (ktherm == 2) then ! return mushy-layer frazil to ocean (POP)
            vi0tmp = fnew*dt / (rhoi*Lfresh)
            dfresh = -rhoi*(vi0new - vi0tmp)/dt
            dfsalt = ice_ref_salinity*p001*dfresh
            fresh  = fresh + dfresh
            fsalt  = fsalt + dfsalt
            frazil_diag = frazil - vi0tmp
         ! elseif ktherm==1 do nothing
         endif
      endif

      !-----------------------------------------------------------------
      ! Decide how to distribute the new ice.
      !-----------------------------------------------------------------

      hsurp  = c0
      ai0new = c0

      if (vi0new > c0) then

         ! new ice area and thickness
         ! hin_max(0) < new ice thickness < hin_max(1)
         if (aice0 > puny) then
            hi0new = max(vi0new/aice0, hfrazilmin)
            if (hi0new > hi0max .and. aice0+puny < c1) then
               ! distribute excess volume over all categories (below)
               hi0new = hi0max
               ai0new = aice0
               vsurp      = vi0new - ai0new*hi0new
               hsurp  = vsurp / aice
               vi0new = ai0new*hi0new
            else
               ! put ice in a single category, with hsurp = 0
               ai0new = vi0new/hi0new
            endif
         else                ! aice0 < puny
            hsurp = vi0new/aice ! new thickness in each cat
            vi0new = c0
         endif               ! aice0 > puny
      endif                  ! vi0new > puny

      !-----------------------------------------------------------------
      ! Distribute excess ice volume among ice categories by increasing
      ! ice thickness, leaving ice area unchanged.
      !
      ! NOTE: If new ice contains globally conserved tracers
      !       (e.g., isotopes from seawater), code must be added here.
      !
      ! The mushy formulation (ktherm=2) puts the new ice only at the
      ! bottom of existing ice and adjusts the layers accordingly.
      ! The other formulations distribute the new ice throughout the 
      ! existing ice column.
      !-----------------------------------------------------------------

      if (hsurp > c0) then   ! add ice to all categories

         do n = 1, ncat

            vsurp = hsurp * aicen(n)

            ! update ice age due to freezing (new ice age = dt)
            vtmp = vicen(n) + vsurp
            if (tr_iage .and. vtmp > puny) &
                trcrn(nt_iage,n) = &
               (trcrn(nt_iage,n)*vicen(n) + dt*vsurp) / vtmp

            if (tr_lvl .and. vicen(n) > puny) then
                trcrn(nt_vlvl,n) = &
               (trcrn(nt_vlvl,n)*vicen(n) + &
                trcrn(nt_alvl,n)*vsurp) / vtmp
            endif

            if (tr_aero .and. vtmp > puny) then
               do it = 1, n_aero
                  trcrn(nt_aero+2+4*(it-1),n) = &
                  trcrn(nt_aero+2+4*(it-1),n)*vicen(n) / vtmp
                  trcrn(nt_aero+3+4*(it-1),n) = &
                  trcrn(nt_aero+3+4*(it-1),n)*vicen(n) / vtmp
               enddo
            endif

            ! update category volumes
            vicen(n) = vtmp

            if (ktherm == 2) then
               vsurp = hsurp * aicen(n)  ! note - save this above?
               vtmp = vicen(n) - vsurp   ! vicen is the new volume
               if (vicen(n) > c0) then
                  call update_vertical_tracers(nilyr, &
                              trcrn(nt_qice:nt_qice+nilyr-1,n), &
                              vtmp, vicen(n), qi0new)
                  call update_vertical_tracers(nilyr, &
                              trcrn(nt_sice:nt_sice+nilyr-1,n), &
                              vtmp, vicen(n), Si0new)
               endif
            else
               do k = 1, nilyr
                  ! factor of nilyr cancels out
                  vsurp = hsurp * aicen(n)  ! note - save this above?
                  vtmp = vicen(n) - vsurp      ! vicen is the new volume
                  if (vicen(n) > c0) then
                     ! enthalpy
                     trcrn(nt_qice+k-1,n) = &
                    (trcrn(nt_qice+k-1,n)*vtmp + qi0new*vsurp) / vicen(n)
                     ! salinity
                     if (.not. solve_zsal) & 
                     trcrn(nt_sice+k-1,n) = &
                    (trcrn(nt_sice+k-1,n)*vtmp + Sprofile(k)*vsurp) / vicen(n) 
                  endif
               enddo               ! k
            endif                  ! ktherm

         enddo                     ! n

      endif ! hsurp > 0

      !-----------------------------------------------------------------
      ! Combine new ice grown in open water with category 1 ice.
      ! Assume that vsnon and esnon are unchanged.
      ! The mushy formulation assumes salt from frazil is added uniformly
      ! to category 1, while the others use a salinity profile.
      !-----------------------------------------------------------------

      if (vi0new > c0) then  ! add ice to category 1

         area1    = aicen(1)   ! save
         vice1    = vicen(1)   ! save
         aicen(1) = aicen(1) + ai0new
         aice0    = aice0    - ai0new
         vicen(1) = vicen(1) + vi0new

         trcrn(nt_Tsfc,1) = &
            (trcrn(nt_Tsfc,1)*area1 + Tf*ai0new)/aicen(1) 
         trcrn(nt_Tsfc,1) = min (trcrn(nt_Tsfc,1), c0)

         if (tr_FY) then
            trcrn(nt_FY,1) = &
           (trcrn(nt_FY,1)*area1 + ai0new)/aicen(1)
            trcrn(nt_FY,1) = min(trcrn(nt_FY,1), c1)
         endif

         if (vicen(1) > puny) then
            if (tr_iage) &
               trcrn(nt_iage,1) = &
              (trcrn(nt_iage,1)*vice1 + dt*vi0new)/vicen(1)

            if (tr_aero) then
               do it = 1, n_aero
                  trcrn(nt_aero+2+4*(it-1),1) = &
                  trcrn(nt_aero+2+4*(it-1),1)*vice1/vicen(1)
                  trcrn(nt_aero+3+4*(it-1),1) = &
                  trcrn(nt_aero+3+4*(it-1),1)*vice1/vicen(1)
               enddo
            endif

            if (tr_lvl) then
                alvl = trcrn(nt_alvl,1)
                trcrn(nt_alvl,1) = &
               (trcrn(nt_alvl,1)*area1 + ai0new)/aicen(1)
                trcrn(nt_vlvl,1) = &
               (trcrn(nt_vlvl,1)*vice1 + vi0new)/vicen(1)
            endif

            if (tr_pond_cesm .or. tr_pond_topo) then
               trcrn(nt_apnd,1) = &
               trcrn(nt_apnd,1)*area1/aicen(1)
            elseif (tr_pond_lvl) then
               if (trcrn(nt_alvl,1) > puny) then
                  trcrn(nt_apnd,1) = &
                  trcrn(nt_apnd,1) * alvl*area1 / (trcrn(nt_alvl,1)*aicen(1))
               endif
            endif
         endif

         do k = 1, nilyr
               
            if (vicen(1) > c0) then
               ! factor of nilyr cancels out
               ! enthalpy
               trcrn(nt_qice+k-1,1) = &
              (trcrn(nt_qice+k-1,1)*vice1 + qi0new*vi0new)/vicen(1)
               ! salinity
               if (.NOT. solve_zsal)&
               trcrn(nt_sice+k-1,1) = &
              (trcrn(nt_sice+k-1,1)*vice1 + Sprofile(k)*vi0new)/vicen(1)
            endif
         enddo

      endif ! vi0new > 0

      if (l_conservation_check) then

         do n = 1, ncat
            eicen(n) = c0
            do k = 1, nilyr
               eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                        * vicen(n)/real(nilyr,kind=dbl_kind)
            enddo
         enddo
         call column_sum (ncat, vicen, vice_final)
         call column_sum (ncat, eicen, eice_final)

         fieldid = 'vice, add_new_ice'
         call column_conservation_check (fieldid,               &
                                         vice_init, vice_final, &
                                         puny,                  &
                                         l_stop)
         fieldid = 'eice, add_new_ice'
         call column_conservation_check (fieldid,               &
                                         eice_init, eice_final, &
                                         puny*Lfresh*rhoi,      &
                                         l_stop)
         if (l_stop) then
            stop_label = 'add_new_ice: Column conservation error'
            return
         endif

      endif ! l_conservation_check

      !-----------------------------------------------------------------
      ! Biogeochemistry
      !-----------------------------------------------------------------     
     if (tr_brine .or. nbtrcr > 0) &
        call add_new_ice_bgc(dt,         nblyr,                &
                             ncat, nilyr, nltrcr, &
                             bgrid,      cgrid,      igrid,    &
                             aicen_init, vicen_init, vi0_init, &
                             aicen,      vicen,      vsnon1,   &
                             vi0new,     ntrcr,      trcrn,    &
                             nbtrcr,     sss,        ocean_bio,&
                             flux_bio,   hsurp,                &
                             l_stop,     stop_label,           &
                             l_conservation_check)

      end subroutine add_new_ice

!=======================================================================

      end module ice_therm_itd

!=======================================================================
