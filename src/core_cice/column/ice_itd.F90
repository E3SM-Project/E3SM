!  SVN:$Id: ice_itd.F90 1178 2017-03-08 19:24:07Z eclare $
!=======================================================================

! Routines to initialize the ice thickness distribution and
! utilities to redistribute ice among categories. These routines
! are not specific to a particular numerical implementation.
!
! See Bitz, C.M., and W.H. Lipscomb, 1999:
! An energy-conserving thermodynamic model of sea ice,
! J. Geophys. Res., 104, 15,669--15,677.
!
! See Bitz, C.M., M.M. Holland, A.J. Weaver, M. Eby, 2001:
! Simulating the ice-thickness distribution in a climate model,
! J. Geophys. Res., 106, 2441--2464.
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb and Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
!
! 2004 WHL: Added multiple snow layers, block structure, cleanup_itd
! 2006 ECH: Added WMO standard ice thickness categories as kcatbound=2
!           Streamlined for efficiency 
!           Converted to free source form (F90)
! 2014 ECH: Converted to column package

      module ice_itd

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c2, p001, puny, p5, &
          Lfresh, rhos, ice_ref_salinity, hs_min, cp_ice, Tocnfrz, rhoi
      use ice_warnings, only: &
          add_warning

      implicit none
      save

      private
      public :: aggregate_area, shift_ice, column_sum, &
                column_conservation_check, cleanup_itd, reduce_area

!=======================================================================

      contains

!=======================================================================

! Aggregate ice area (but not other state variables) over thickness 
! categories.
!
! authors: William H. Lipscomb, LANL

      subroutine aggregate_area (ncat, aicen, aice, aice0)

      integer (kind=int_kind), intent(in) :: &
         ncat      ! number of thickness categories

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen     ! concentration of ice

      real (kind=dbl_kind), intent(inout) :: &
         aice, &   ! concentration of ice
         aice0     ! concentration of open water

      ! local variables

      integer (kind=int_kind) :: n

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      aice = c0
      do n = 1, ncat
         aice = aice + aicen(n)
      enddo                     ! n

      ! open water fraction
      aice0 = max (c1 - aice, c0)

      end subroutine aggregate_area

!=======================================================================

! Rebins thicknesses into defined categories
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL

      subroutine rebin (ntrcr,    trcr_depend,     &
                        trcr_base,                 &
                        n_trcr_strata,             &
                        nt_strata,                 &
                        aicen,    trcrn,           &
                        vicen,    vsnon,           &
                        ncat,     hin_max,         &
                        l_stop,   stop_label)

      integer (kind=int_kind), intent(in) :: &
         ntrcr , & ! number of tracers in use
         ncat      ! number of thickness categories

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (ncat), intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice           (m)
         vsnon     ! volume per unit area of snow          (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension(0:ncat), intent(in) :: &
         hin_max   ! category limits (m)

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      character (char_len), intent(out) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         n         ! category index

      logical (kind=log_kind) :: &
         shiftflag          ! = .true. if ice must be shifted

      integer (kind=int_kind), dimension (ncat) :: &
         donor              ! donor category index

      real (kind=dbl_kind), dimension (ncat) :: &
         daice          , & ! ice area transferred
         dvice          , & ! ice volume transferred
         hicen              ! ice thickness for each cat (m)

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.

      do n = 1, ncat
         donor(n) = 0
         daice(n) = c0
         dvice(n) = c0

      !-----------------------------------------------------------------
      ! Compute ice thickness.
      !-----------------------------------------------------------------
         if (aicen(n) > puny) then
            hicen(n) = vicen(n) / aicen(n)
         else
            hicen(n) = c0
         endif
      enddo                     ! n

      !-----------------------------------------------------------------
      ! make sure thickness of cat 1 is at least hin_max(0)
      !-----------------------------------------------------------------

      if (aicen(1) > puny) then
         if (hicen(1) <= hin_max(0) .and. hin_max(0) > c0 ) then
            aicen(1) = vicen(1) / hin_max(0)
            hicen(1) = hin_max(0)
         endif
      endif

      !-----------------------------------------------------------------
      ! If a category thickness is not in bounds, shift the
      ! entire area, volume, and energy to the neighboring category
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Move thin categories up
      !-----------------------------------------------------------------

      do n = 1, ncat-1          ! loop over category boundaries

      !-----------------------------------------------------------------
      ! identify thicknesses that are too big
      !-----------------------------------------------------------------
         shiftflag = .false.
         if (aicen(n) > puny .and. &
             hicen(n) > hin_max(n)) then
            shiftflag = .true.
            donor(n) = n
            daice(n) = aicen(n)
            dvice(n) = vicen(n)
         endif

         if (shiftflag) then

      !-----------------------------------------------------------------
      ! shift ice between categories
      !-----------------------------------------------------------------

            call shift_ice (ntrcr,    ncat,       &
                            trcr_depend,          &
                            trcr_base,            &
                            n_trcr_strata,        &
                            nt_strata,            &
                            aicen,    trcrn,      &
                            vicen,    vsnon,      &
                            hicen,    donor,      &
                            daice,    dvice,      &
                            l_stop,   stop_label)

      !-----------------------------------------------------------------
      ! reset shift parameters
      !-----------------------------------------------------------------

            donor(n) = 0
            daice(n) = c0
            dvice(n) = c0

         endif                  ! shiftflag

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Move thick categories down
      !-----------------------------------------------------------------

      do n = ncat-1, 1, -1      ! loop over category boundaries

      !-----------------------------------------------------------------
      ! identify thicknesses that are too small
      !-----------------------------------------------------------------

         shiftflag = .false.
         if (aicen(n+1) > puny .and. &
             hicen(n+1) <= hin_max(n)) then
            shiftflag = .true.
            donor(n) = n+1
            daice(n) = aicen(n+1)
            dvice(n) = vicen(n+1)
         endif

         if (shiftflag) then

      !-----------------------------------------------------------------
      ! shift ice between categories
      !-----------------------------------------------------------------

            call shift_ice (ntrcr,    ncat,       &
                            trcr_depend,          &
                            trcr_base,            &
                            n_trcr_strata,        &
                            nt_strata,            &
                            aicen,    trcrn,      &
                            vicen,    vsnon,      &
                            hicen,    donor,      &
                            daice,    dvice,      &
                            l_stop,   stop_label)

      !-----------------------------------------------------------------
      ! reset shift parameters
      !-----------------------------------------------------------------

            donor(n) = 0
            daice(n) = c0
            dvice(n) = c0

         endif                  ! shiftflag

      enddo                     ! n

      end subroutine rebin

!=======================================================================

! Reduce area when ice melts for special case of ncat=1
!
! Use CSM 1.0-like method of reducing ice area
! when melting occurs: assume only half the ice volume
! change goes to thickness decrease, the other half
! to reduction in ice fraction
!
! authors: C. M. Bitz, UW
! modified by: Elizabeth C. Hunke, LANL

      subroutine reduce_area (hin_max,            &
                              aicen,     vicen,   &
                              aicen_init,vicen_init)

      real (kind=dbl_kind), intent(in) :: &
         hin_max       ! lowest category boundary

      real (kind=dbl_kind), intent(inout) :: &
         aicen     , & ! concentration of ice
         vicen         ! volume per unit area of ice          (m)

      real (kind=dbl_kind), intent(in) :: &
         aicen_init, & ! old ice area for category 1 (m)
         vicen_init    ! old ice volume for category 1 (m)

      ! local variables

      real (kind=dbl_kind) :: &
         hi0       , & ! initial hi
         hi1       , & ! current hi
         dhi           ! hi1 - hi0

            hi0 = c0
            if (aicen_init > c0) &
                hi0 = vicen_init / aicen_init

            hi1 = c0
            if (aicen > c0) &
                hi1 = vicen / aicen

            ! make sure thickness of cat 1 is at least hin_max(0)
            if (hi1 <= hin_max .and. hin_max > c0 ) then
               aicen = vicen / hin_max
               hi1 = hin_max
            endif

            if (aicen > c0) then
               dhi = hi1 - hi0
               if (dhi < c0) then
                  hi1  = vicen / aicen
                  aicen = c2 * vicen / (hi1 + hi0)
               endif
            endif

      end subroutine reduce_area

!=======================================================================

! Shift ice across category boundaries, conserving area, volume, and
! energy.
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL

      subroutine shift_ice (ntrcr,    ncat,        &
                            trcr_depend,           &
                            trcr_base,             &
                            n_trcr_strata,         &
                            nt_strata,             &
                            aicen,    trcrn,       &
                            vicen,    vsnon,       &
                            hicen,    donor,       &
                            daice,    dvice,       &
                            l_stop,   stop_label)

      use ice_colpkg_tracers, only: colpkg_compute_tracers

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         ntrcr     ! number of tracers in use

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers

      ! NOTE: Third index of donor, daice, dvice should be ncat-1,
      !       except that compilers would have trouble when ncat = 1 
      integer (kind=int_kind), dimension(:), intent(in) :: &
         donor             ! donor category index

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         daice         , & ! ice area transferred across boundary
         dvice         , & ! ice volume transferred across boundary
         hicen             ! ice thickness for each cat        (m)

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      character (char_len), intent(out) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         n             , & ! thickness category index
         nr            , & ! receiver category
         nd            , & ! donor category
         it            , & ! tracer index
         ntr           , & ! tracer index
         itl               ! loop index

      real (kind=dbl_kind), dimension(ntrcr,ncat) :: &
         atrcrn            ! aicen*trcrn   

      real (kind=dbl_kind) :: &
         dvsnow        , & ! snow volume transferred
         datrcr            ! aicen*train transferred

      logical (kind=log_kind) :: &
        daice_negative     , & ! true if daice < -puny
        dvice_negative     , & ! true if dvice < -puny
        daice_greater_aicen, & ! true if daice > aicen
        dvice_greater_vicen    ! true if dvice > vicen

      real (kind=dbl_kind) :: &
        worka, workb

      character(len=char_len_long) :: &
        warning ! warning message

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.

      !-----------------------------------------------------------------
      ! Define variables equal to aicen*trcrn, vicen*trcrn, vsnon*trcrn
      !-----------------------------------------------------------------

      do n = 1, ncat
         do it = 1, ntrcr
            atrcrn(it,n) = trcrn(it,n)*(trcr_base(it,1) * aicen(n) &
                                      + trcr_base(it,2) * vicen(n) &
                                      + trcr_base(it,3) * vsnon(n))
            if (n_trcr_strata(it) > 0) then
               do itl = 1, n_trcr_strata(it)
                  ntr = nt_strata(it,itl)
                  atrcrn(it,n) = atrcrn(it,n) * trcrn(ntr,n)
               enddo
            endif
         enddo ! it
      enddo    ! n

      !-----------------------------------------------------------------
      ! Check for daice or dvice out of range, allowing for roundoff error
      !-----------------------------------------------------------------

      do n = 1, ncat-1

         daice_negative      = .false.
         dvice_negative      = .false.
         daice_greater_aicen = .false.
         dvice_greater_vicen = .false.

            if (donor(n) > 0) then
               nd = donor(n)

               if (daice(n) < c0) then
                  if (daice(n) > -puny*aicen(nd)) then
                     daice(n) = c0 ! shift no ice
                     dvice(n) = c0
                  else
                     daice_negative = .true.
                  endif
               endif
         
               if (dvice(n) < c0) then
                  if (dvice(n) > -puny*vicen(nd)) then   
                     daice(n) = c0 ! shift no ice
                     dvice(n) = c0
                  else
                     dvice_negative = .true.
                  endif
               endif

               if (daice(n) > aicen(nd)*(c1-puny)) then
                  if (daice(n) < aicen(nd)*(c1+puny)) then
                     daice(n) = aicen(nd)
                     dvice(n) = vicen(nd)
                  else
                     daice_greater_aicen = .true.
                  endif
               endif    

               if (dvice(n) > vicen(nd)*(c1-puny)) then
                  if (dvice(n) < vicen(nd)*(c1+puny)) then
                     daice(n) = aicen(nd)
                     dvice(n) = vicen(nd)
                  else
                     dvice_greater_vicen = .true.
                  endif
               endif

            endif               ! donor > 0

      !-----------------------------------------------------------------
      ! error messages
      !-----------------------------------------------------------------

         if (daice_negative) then
               if (donor(n) > 0 .and.  &
                   daice(n) <= -puny*aicen(nd)) then
                  write(warning,*) ' '
                  call add_warning(warning)
                  write(warning,*) 'shift_ice: negative daice'
                  call add_warning(warning)
                  write(warning,*) 'boundary, donor cat:', n, nd
                  call add_warning(warning)
                  write(warning,*) 'daice =', daice(n)
                  call add_warning(warning)
                  write(warning,*) 'dvice =', dvice(n)
                  call add_warning(warning)
                  l_stop = .true.
                  stop_label = 'shift_ice: negative daice'
               endif
         endif
         if (l_stop) return

         if (dvice_negative) then
               if (donor(n) > 0 .and.  &
                   dvice(n) <= -puny*vicen(nd)) then
                  write(warning,*) ' '
                  call add_warning(warning)
                  write(warning,*) 'shift_ice: negative dvice'
                  call add_warning(warning)
                  write(warning,*) 'boundary, donor cat:', n, nd
                  call add_warning(warning)
                  write(warning,*) 'daice =', daice(n)
                  call add_warning(warning)
                  write(warning,*) 'dvice =', dvice(n)
                  call add_warning(warning)
                  l_stop = .true.
                  stop_label = 'shift_ice: negative dvice'
               endif
         endif
         if (l_stop) return

         if (daice_greater_aicen) then
               if (donor(n) > 0) then
                  nd = donor(n)
                  if (daice(n) >= aicen(nd)*(c1+puny)) then
                     write(warning,*) ' '
                     call add_warning(warning)
                     write(warning,*) 'shift_ice: daice > aicen'
                     call add_warning(warning)
                     write(warning,*) 'boundary, donor cat:', n, nd
                     call add_warning(warning)
                     write(warning,*) 'daice =', daice(n)
                     call add_warning(warning)
                     write(warning,*) 'aicen =', aicen(nd)
                     call add_warning(warning)
                     l_stop = .true.
                     stop_label = 'shift_ice: daice > aicen'
                  endif
               endif
         endif
         if (l_stop) return

         if (dvice_greater_vicen) then
               if (donor(n) > 0) then
                  nd = donor(n)
                  if (dvice(n) >= vicen(nd)*(c1+puny)) then
                     write(warning,*) ' '
                     call add_warning(warning)
                     write(warning,*) 'shift_ice: dvice > vicen'
                     call add_warning(warning)
                     write(warning,*) 'boundary, donor cat:', n, nd
                     call add_warning(warning)
                     write(warning,*) 'dvice =', dvice(n)
                     call add_warning(warning)
                     write(warning,*) 'vicen =', vicen(nd)
                     call add_warning(warning)
                     l_stop = .true.
                     stop_label = 'shift_ice: dvice > vicen'
                  endif
               endif
         endif
         if (l_stop) return

      !-----------------------------------------------------------------
      ! transfer volume and energy between categories
      !-----------------------------------------------------------------

         if (daice(n) > c0) then ! daice(n) can be < puny

            nd = donor(n)
            worka = daice(n) / aicen(nd)
            if (nd  ==  n) then
               nr = nd+1
            else                ! nd = n+1
               nr = n
            endif

            aicen(nd) = aicen(nd) - daice(n)
            aicen(nr) = aicen(nr) + daice(n)

            vicen(nd) = vicen(nd) - dvice(n)
            vicen(nr) = vicen(nr) + dvice(n)

            dvsnow = vsnon(nd) * worka
            vsnon(nd) = vsnon(nd) - dvsnow
            vsnon(nr) = vsnon(nr) + dvsnow
            workb = dvsnow

            do it = 1, ntrcr
               nd = donor(n)
               if (nd == n) then
                  nr = nd+1
               else             ! nd = n+1
                  nr = n
               endif

               datrcr = trcrn(it,nd)*(trcr_base(it,1) * daice(n) &
                                    + trcr_base(it,2) * dvice(n) &
                                    + trcr_base(it,3) * workb)
               if (n_trcr_strata(it) > 0) then
                  do itl = 1, n_trcr_strata(it)
                     ntr = nt_strata(it,itl)
                     datrcr = datrcr * trcrn(ntr,nd)
                  enddo
               endif

               atrcrn(it,nd) = atrcrn(it,nd) - datrcr
               atrcrn(it,nr) = atrcrn(it,nr) + datrcr
            
            enddo ! ntrcr
         endif    ! daice
      enddo       ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Update ice thickness and tracers
      !-----------------------------------------------------------------

      do n = 1, ncat

         if (aicen(n) > puny) then
            hicen(n) = vicen (n) / aicen(n)                
         else
            hicen(n) = c0
         endif

      !-----------------------------------------------------------------
      ! Compute new tracers
      !-----------------------------------------------------------------

         call colpkg_compute_tracers (ntrcr,       trcr_depend, &
                                      atrcrn(:,n), aicen(n),    &
                                      vicen(n),    vsnon(n),    &
                                      trcr_base,   n_trcr_strata,  &
                                      nt_strata,   trcrn(:,n))

      enddo                     ! ncat

      end subroutine shift_ice

!=======================================================================

! For each grid cell, sum field over all ice categories.
!
! author: William H. Lipscomb, LANL

      subroutine column_sum (nsum, xin, xout)

      integer (kind=int_kind), intent(in) :: &
         nsum             ! number of categories/layers

      real (kind=dbl_kind), dimension (nsum), &
         intent(in) :: &
         xin              ! input field

      real (kind=dbl_kind), intent(out) :: &
         xout             ! output field

      ! local variables

      integer (kind=int_kind) :: &
         n                ! category/layer index

      xout = c0
      do n = 1, nsum
         xout = xout + xin(n)
      enddo                 ! n

      end subroutine column_sum

!=======================================================================

! For each physical grid cell, check that initial and final values
! of a conserved field are equal to within a small value.
!
! author: William H. Lipscomb, LANL

      subroutine column_conservation_check (fieldid,          &
                                            x1,       x2,     &
                                            max_err,          &
                                            l_stop)

      real (kind=dbl_kind), intent(in) :: &
         x1            , & ! initial field
         x2                ! final field

      real (kind=dbl_kind), intent(in) :: &
         max_err           ! max allowed error

      character (len=char_len), intent(in) :: &
         fieldid           ! field identifier

      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, abort on return

      character(len=char_len_long) :: &
         warning ! warning message
      
      ! local variables

      if (abs (x2-x1) > max_err) then
         l_stop = .true.
         write(warning,*) ' '
         call add_warning(warning)
         write(warning,*) 'Conservation error: ', trim(fieldid)
         call add_warning(warning)
         write(warning,*) 'Initial value =', x1
         call add_warning(warning)
         write(warning,*) 'Final value =',   x2
         call add_warning(warning)
         write(warning,*) 'Difference =', x2 - x1
         call add_warning(warning)
      endif

      end subroutine column_conservation_check

!=======================================================================

! Cleanup subroutine that rebins thickness categories if necessary,
!  eliminates very small ice areas while conserving mass and energy, 
!  aggregates state variables, and does a boundary call.  
! It is a good idea to call this subroutine after the thermodynamics
!  (thermo_vertical/thermo_itd) and again after the dynamics 
!  (evp/transport/ridging).
!
! author: William H. Lipscomb, LANL

      subroutine cleanup_itd (dt,          ntrcr,      &
                              nilyr,       nslyr,      &
                              ncat,        hin_max,    &
                              aicen,       trcrn,      &
                              vicen,       vsnon,      &
                              aice0,       aice,       &   
                              n_aero,                  &
                              nbtrcr,      nblyr,      &
                              l_stop,      stop_label, &
                              tr_aero,                 &
                              tr_pond_topo,            &
                              heat_capacity,           & 
                              first_ice,               &
                              trcr_depend, trcr_base,  &
                              n_trcr_strata,nt_strata, &
                              fpond,       fresh,      &
                              fsalt,       fhocn,      &
                              faero_ocn,   fzsal,      &
                              flux_bio,    limit_aice_in)

      integer (kind=int_kind), intent(in) :: & 
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr , & ! number of bio layers
         nslyr , & ! number of snow layers
         ntrcr , & ! number of tracers in use
         nbtrcr, & ! number of bio tracers in use
         n_aero    ! number of aerosol tracers
 
      real (kind=dbl_kind), intent(in) :: & 
         dt        ! time step 
 
      real (kind=dbl_kind), dimension(0:ncat), intent(in) :: &
         hin_max   ! category boundaries (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: & 
         aicen , & ! concentration of ice 
         vicen , & ! volume per unit area of ice          (m) 
         vsnon     ! volume per unit area of snow         (m) 
 
      real (kind=dbl_kind), dimension (:,:), intent(inout) :: & 
         trcrn     ! ice tracers 

      real (kind=dbl_kind), intent(inout) :: & 
         aice  , & ! total ice concentration
         aice0     ! concentration of open water 
     
      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      logical (kind=log_kind), intent(in) :: &
         tr_aero,      & ! aerosol flag
         tr_pond_topo, & ! topo pond flag
         heat_capacity   ! if false, ice and snow have zero heat capacity

      logical (kind=log_kind), dimension(ncat),intent(inout) :: &
         first_ice   ! For bgc and S tracers. set to true if zapping ice.

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      character (char_len), intent(out) :: stop_label

      ! ice-ocean fluxes (required for strict conservation)

      real (kind=dbl_kind), intent(inout), optional :: &
         fpond    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean        (kg/m^2/s)
         fhocn    , & ! net heat flux to ocean     (W/m^2)
         fzsal        ! net salt flux to ocean from zsalinity (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout), optional :: &
         flux_bio     ! net tracer flux to ocean from biology (mmol/m^2/s)

      real (kind=dbl_kind), dimension (:), &
         intent(inout), optional :: &
         faero_ocn    ! aerosol flux to ocean     (kg/m^2/s)

      logical (kind=log_kind), intent(in), optional ::   &
         limit_aice_in      ! if false, allow aice to be out of bounds
                            ! may want to allow this for unit tests

      ! local variables

      integer (kind=int_kind) :: &
         n        , & ! category index
         it           ! tracer index

      real (kind=dbl_kind) &
         dfpond   , & ! zapped pond water flux (kg/m^2/s)
         dfresh   , & ! zapped fresh water flux (kg/m^2/s)
         dfsalt   , & ! zapped salt flux   (kg/m^2/s)
         dfhocn   , & ! zapped energy flux ( W/m^2)
         dfzsal       ! zapped salt flux for zsalinity (kg/m^2/s)

      real (kind=dbl_kind), dimension (n_aero) :: &
         dfaero_ocn   ! zapped aerosol flux   (kg/m^2/s)

      real (kind=dbl_kind), dimension (ntrcr) :: &
         dflux_bio    ! zapped biology flux  (mmol/m^2/s)

      logical (kind=log_kind) ::   &
         limit_aice         ! if true, check for aice out of bounds

      character(len=char_len_long) :: &
         warning ! warning message
      
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      if (present(limit_aice_in)) then
         limit_aice = limit_aice_in
      else
         limit_aice = .true.
      endif

      l_stop = .false.

      dfpond = c0
      dfresh = c0
      dfsalt = c0
      dfhocn = c0
      dfaero_ocn(:) = c0
      dflux_bio(:) = c0
      dfzsal = c0

      !-----------------------------------------------------------------
      ! Compute total ice area.
      !-----------------------------------------------------------------

      call aggregate_area (ncat, aicen, aice, aice0)

      if (limit_aice) then  ! check for aice out of bounds
         if (aice > c1+puny .or. aice < -puny) then
            l_stop = .true.
            stop_label = 'aggregate ice area out of bounds'
            write(warning,*) 'aice:', aice
            call add_warning(warning)
            do n = 1, ncat
               write(warning,*) 'n, aicen:', n, aicen(n)
               call add_warning(warning)
            enddo
            return
         endif
      endif                     ! limit_aice

      !-----------------------------------------------------------------
      ! Identify grid cells with ice.
      !-----------------------------------------------------------------

      if (aice > puny) then

      !-----------------------------------------------------------------
      ! Make sure ice in each category is within its thickness bounds.
      ! NOTE: The rebin subroutine is needed only in the rare cases
      !       when the linear_itd subroutine cannot transfer ice
      !       correctly (e.g., very fast ice growth).
      !-----------------------------------------------------------------

         call rebin (ntrcr,      trcr_depend, &
                     trcr_base,               &
                     n_trcr_strata,           &
                     nt_strata,               &
                     aicen,      trcrn,       &
                     vicen,      vsnon,       &
                     ncat,       hin_max,     &
                     l_stop,     stop_label)

      endif ! aice > puny

      !-----------------------------------------------------------------
      ! Zero out ice categories with very small areas.
      !-----------------------------------------------------------------

      if (limit_aice) then
         call zap_small_areas (dt,           ntrcr,         &
                               ncat,         n_aero,        &
                               nblyr,                       &
                               nilyr,        nslyr,         &
                               aice,         aice0,         &
                               aicen,        trcrn,         &
                               vicen,        vsnon,         &
                               dfpond,                      &
                               dfresh,       dfsalt,        &
                               dfhocn,       dfaero_ocn,    &
                               tr_aero,      tr_pond_topo,  &
                               first_ice,    nbtrcr,        &
                               dfzsal,       dflux_bio,     & 
                               l_stop,       stop_label)

         if (l_stop) then
            write(warning,*) 'aice:', aice
            call add_warning(warning)
            do n = 1, ncat
               write(warning,*) 'n, aicen:', n, aicen(n)
               call add_warning(warning)
            enddo
            return
         endif

      endif   ! l_limit_aice

    !-------------------------------------------------------------------
    ! Zap snow that has out of bounds temperatures
    !-------------------------------------------------------------------

      call zap_snow_temperature(dt,            ncat,     &
                                heat_capacity, nblyr,    &
                                nslyr,         aicen,    &
                                trcrn,         vsnon,    &
                                dfresh,        dfhocn,   &
                                dfaero_ocn,    tr_aero,  &
                                dflux_bio,     nbtrcr,   &
                                n_aero,        ntrcr)

    !-------------------------------------------------------------------
    ! Update ice-ocean fluxes for strict conservation
    !-------------------------------------------------------------------

      if (present(fpond)) &
           fpond        = fpond        + dfpond 
      if (present(fresh)) &
           fresh        = fresh        + dfresh 
      if (present(fsalt)) &
           fsalt        = fsalt        + dfsalt
      if (present(fhocn)) &
           fhocn        = fhocn        + dfhocn
      if (present(faero_ocn)) then
         do it = 1, n_aero
           faero_ocn(it) = faero_ocn(it) + dfaero_ocn(it)
         enddo
      endif
      if (present(flux_bio)) then
         do it = 1, nbtrcr
           flux_bio (it) = flux_bio(it) + dflux_bio(it)
         enddo
      endif
      if (present(fzsal)) &
           fzsal        = fzsal         + dfzsal     

      !----------------------------------------------------------------
      ! If using zero-layer model (no heat capacity), check that the 
      ! energy of snow and ice is correct. 
      !----------------------------------------------------------------

      if ((.not. heat_capacity) .and. aice > puny) then
         call zerolayer_check (ncat,       nilyr,    &
                               nslyr,      aicen,    &
                               vicen,      vsnon,    &
                               trcrn,      l_stop,   &
                               stop_label)
      endif

      end subroutine cleanup_itd

!=======================================================================

! For each ice category in each grid cell, remove ice if the fractional
! area is less than puny.
!
! author: William H. Lipscomb, LANL

      subroutine zap_small_areas (dt,        ntrcr,        &
                                  ncat,      n_aero,       &
                                  nblyr,                   &
                                  nilyr,     nslyr,        &
                                  aice,      aice0,        &
                                  aicen,     trcrn,        &
                                  vicen,     vsnon,        &
                                  dfpond,                  &
                                  dfresh,    dfsalt,       &
                                  dfhocn,    dfaero_ocn,   &
                                  tr_aero,   tr_pond_topo, &
                                  first_ice, nbtrcr,       &
                                  dfzsal,    dflux_bio,    &
                                  l_stop,    stop_label)

      use ice_colpkg_tracers, only: nt_Tsfc, nt_qice, nt_qsno, nt_aero, &
                             nt_apnd, nt_hpnd, nt_fbri, tr_brine, nt_bgc_S, &
                             bio_index
      use ice_colpkg_shared, only:  solve_zsal, skl_bgc, z_tracers, min_salin, &
                             sk_l, rhosi
      use ice_zbgc_shared, only: zap_small_bgc

      integer (kind=int_kind), intent(in) :: &
         ncat     , & ! number of thickness categories
         nilyr    , & ! number of ice layers
         nblyr    , & ! number of bio layers
         nslyr    , & ! number of snow layers
         ntrcr    , & ! number of tracers in use
         n_aero   , & ! number of aerosol tracers
         nbtrcr       ! number of biology tracers

      real (kind=dbl_kind), intent(in) :: &
         dt           ! time step

      real (kind=dbl_kind), intent(inout) :: &
         aice     , & ! total ice concentration
         aice0        ! concentration of open water

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice          (m)
         vsnon        ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn        ! ice tracers

      real (kind=dbl_kind), intent(out) :: &
         dfpond   , & ! zapped pond water flux (kg/m^2/s)
         dfresh   , & ! zapped fresh water flux (kg/m^2/s)
         dfsalt   , & ! zapped salt flux   (kg/m^2/s)
         dfhocn   , & ! zapped energy flux ( W/m^2)
         dfzsal       ! zapped salt flux from zsalinity(kg/m^2/s) 

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         dfaero_ocn   ! zapped aerosol flux   (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dflux_bio     ! zapped bio tracer flux from biology (mmol/m^2/s)

      logical (kind=log_kind), intent(in) :: &
         tr_aero, &   ! aerosol flag
         tr_pond_topo ! pond flag

      logical (kind=log_kind), dimension (:),intent(inout) :: &
         first_ice    ! For bgc tracers.  Set to true if zapping ice 

      logical (kind=log_kind), intent(out) :: &
         l_stop       ! if true, abort on return

      character (char_len), intent(out) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         n, k, it, & !counting indices
         blevels

      real (kind=dbl_kind) :: xtmp      ! temporary variables
      real (kind=dbl_kind) , dimension (1):: trcr_skl    
      real (kind=dbl_kind) , dimension (nblyr+1):: bvol     

      l_stop = .false.

      !-----------------------------------------------------------------
      ! I. Zap categories with very small areas.
      !-----------------------------------------------------------------
      dfzsal = c0
      
      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Count categories to be zapped.
      !-----------------------------------------------------------------

         if (aicen(n) < -puny) then
            l_stop = .true.
            stop_label = 'Zap ice: negative ice area'
            return
         elseif (abs(aicen(n)) /= c0 .and. &
                 abs(aicen(n)) <= puny) then

      !-----------------------------------------------------------------
      ! Account for tracers important for conservation
      !-----------------------------------------------------------------

            if (tr_pond_topo) then
               xtmp = aicen(n) &
                    * trcrn(nt_apnd,n) * trcrn(nt_hpnd,n)
               dfpond = dfpond - xtmp
            endif

            if (tr_aero) then
               do it = 1, n_aero
                  xtmp = (vicen(n)*(trcrn(nt_aero+2+4*(it-1),n)     &
                                  + trcrn(nt_aero+3+4*(it-1),n)))/dt
                  dfaero_ocn(it) = dfaero_ocn(it) + xtmp
               enddo
            endif

           if (solve_zsal) then
               do it = 1, nblyr
                  xtmp = rhosi*trcrn(nt_fbri,n)*vicen(n)*p001&
                        *trcrn(nt_bgc_S+it-1,n)/ &
                         real(nblyr,kind=dbl_kind)/dt
                  dfzsal = dfzsal + xtmp
               enddo                 ! n
           endif

         if (skl_bgc .and. nbtrcr > 0) then
               blevels = 1
               bvol(1) =  aicen(n)*sk_l
               it = 1
               do it = 1, nbtrcr
                  trcr_skl(1) = trcrn(bio_index(it),n)
                  call zap_small_bgc(blevels, dflux_bio(it), &
                       dt, bvol(1:blevels), trcr_skl(blevels))
               enddo
         elseif (z_tracers .and. nbtrcr > 0) then
               blevels = nblyr + 1
               bvol(:) = vicen(n)/real(nblyr,kind=dbl_kind)*trcrn(nt_fbri,n)
               bvol(1) = p5*bvol(1)
               bvol(blevels) = p5*bvol(blevels)
               do it = 1, nbtrcr
                  call zap_small_bgc(blevels, dflux_bio(it), &
                       dt, bvol(1:blevels),trcrn(bio_index(it):bio_index(it)+blevels-1,n))
               enddo
         endif

      !-----------------------------------------------------------------
      ! Zap ice energy and use ocean heat to melt ice
      !-----------------------------------------------------------------

            do k = 1, nilyr
               xtmp = trcrn(nt_qice+k-1,n) / dt &
                    * vicen(n)/real(nilyr,kind=dbl_kind) ! < 0
               dfhocn = dfhocn + xtmp
               trcrn(nt_qice+k-1,n) = c0
            enddo                  ! k

      !-----------------------------------------------------------------
      ! Zap ice and snow volume, add water and salt to ocean
      !-----------------------------------------------------------------

            xtmp = (rhoi*vicen(n)) / dt
            dfresh = dfresh + xtmp

            xtmp = rhoi*vicen(n)*ice_ref_salinity*p001 / dt
            dfsalt = dfsalt + xtmp

            aice0 = aice0 + aicen(n)
            aicen(n) = c0
            vicen(n) = c0
            trcrn(nt_Tsfc,n) = Tocnfrz

      !-----------------------------------------------------------------
      ! Zap snow
      !-----------------------------------------------------------------
            call zap_snow(dt,            nslyr,    &
                          trcrn(:,n),    vsnon(n), &
                          dfresh,        dfhocn,   &
                          dfaero_ocn,    tr_aero,  &
                          dflux_bio,     nbtrcr,   &
                          n_aero,        ntrcr,    &
                          aicen(n),      nblyr)

      !-----------------------------------------------------------------
      ! Zap tracers
      !-----------------------------------------------------------------
         
            if (ntrcr >= 2) then
               do it = 2, ntrcr
                  if (tr_brine .and. it == nt_fbri) then
                     trcrn(it,n) = c1
                  else
                     trcrn(it,n) = c0
                  endif
               enddo
            endif
            first_ice(n) = .true.

         endif ! aicen
      enddo                     ! n

      !-----------------------------------------------------------------
      ! II. Count cells with excess ice (aice > c1) due to roundoff errors.
      !     Zap a little ice in each category so that aice = c1.
      !-----------------------------------------------------------------

      if (aice > (c1+puny)) then
         l_stop = .true.
         stop_label = 'Zap ice: excess ice area'
         return
      elseif (aice > c1 .and. aice < (c1+puny)) then

         do n = 1, ncat

      !-----------------------------------------------------------------
      ! Account for tracers important for conservation
      !-----------------------------------------------------------------

            if (tr_pond_topo) then
               xtmp = aicen(n) &
                    * trcrn(nt_apnd,n) * trcrn(nt_hpnd,n) &
                    * (aice-c1)/aice
               dfpond = dfpond - xtmp
            endif

            if (tr_aero) then
               do it = 1, n_aero
                  xtmp = (vsnon(n)*(trcrn(nt_aero  +4*(it-1),n)     &
                                  + trcrn(nt_aero+1+4*(it-1),n))    &
                       +  vicen(n)*(trcrn(nt_aero+2+4*(it-1),n)     &
                                  + trcrn(nt_aero+3+4*(it-1),n)))   &
                       * (aice-c1)/aice / dt
                  dfaero_ocn(it) = dfaero_ocn(it) + xtmp
               enddo               ! it
            endif

      !----------------------------------------------------------------- 
      ! Zap ice energy and use ocean heat to melt ice 
      !----------------------------------------------------------------- 
       
            do k = 1, nilyr 
               xtmp = trcrn(nt_qice+k-1,n) &
                    * vicen(n)/real(nilyr,kind=dbl_kind) &
                    * (aice-c1)/aice / dt ! < 0 
               dfhocn = dfhocn + xtmp 
            enddo                  ! k 
 
      !----------------------------------------------------------------- 
      ! Zap snow energy and use ocean heat to melt snow 
      !----------------------------------------------------------------- 
 
            do k = 1, nslyr 
               xtmp = trcrn(nt_qsno+k-1,n) &
                    * vsnon(n)/real(nslyr,kind=dbl_kind) &
                    * (aice-c1)/aice / dt ! < 0 
               dfhocn = dfhocn + xtmp 
            enddo                  ! k
 
      !-----------------------------------------------------------------
      ! Zap ice and snow volume, add water and salt to ocean
      !-----------------------------------------------------------------

            xtmp = (rhoi*vicen(n) + rhos*vsnon(n)) &
                 * (aice-c1)/aice / dt 
            dfresh = dfresh + xtmp 
 
            xtmp = rhoi*vicen(n)*ice_ref_salinity*p001 &
                 * (aice-c1)/aice / dt
            dfsalt = dfsalt + xtmp 
 
            if (solve_zsal) then
            do k = 1,nblyr
               xtmp = rhosi*trcrn(nt_fbri,n)*vicen(n)*p001&
                    /real(nblyr,kind=dbl_kind)*trcrn(nt_bgc_S+k-1,n) &
                    * (aice-c1)/aice /dt
               dfzsal = dfzsal + xtmp
            enddo

            if (vicen(n) > vicen(n)*trcrn(nt_fbri,n)) then
               xtmp = (vicen(n)-vicen(n)*trcrn(nt_fbri,n))*(aice-c1)/&
                      aice*p001*rhosi*min_salin/dt
               dfzsal = dfzsal + xtmp
            endif
            endif ! solve_zsal

            aicen(n) = aicen(n) * (c1/aice) 
            vicen(n) = vicen(n) * (c1/aice) 
            vsnon(n) = vsnon(n) * (c1/aice)
 
      ! Note: Tracers are unchanged.

         enddo                     ! n

      !-----------------------------------------------------------------
      ! Correct aice
      !-----------------------------------------------------------------

         aice = c1
         aice0 = c0

      endif ! aice

      end subroutine zap_small_areas

!=======================================================================

      subroutine zap_snow(dt,         nslyr,    &
                          trcrn,      vsnon,    &
                          dfresh,     dfhocn,   &
                          dfaero_ocn, tr_aero,  &
                          dflux_bio,  nbtrcr,   &
                          n_aero,     ntrcr,    &
                          aicen,      nblyr)

      use ice_colpkg_tracers, only: nt_qsno, nt_aero, bio_index
      use ice_colpkg_shared, only: hs_ssl, z_tracers

      integer (kind=int_kind), intent(in) :: &
         nslyr    , & ! number of snow layers
         n_aero   , & ! number of aerosol tracers
         ntrcr    , & ! number of tracers in use
         nblyr    , & ! number of bio  layers
         nbtrcr
 
      real (kind=dbl_kind), intent(in) :: &
         dt           ! time step

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         trcrn        ! ice tracers

      real (kind=dbl_kind), intent(in) :: &
         aicen        ! ice area fraction

      real (kind=dbl_kind), intent(inout) :: &
         vsnon    , & ! volume per unit area of snow         (m)
         dfresh   , & ! zapped fresh water flux (kg/m^2/s)
         dfhocn       ! zapped energy flux ( W/m^2)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dfaero_ocn   ! zapped aerosol flux   (kg/m^2/s)

     real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dflux_bio     ! zapped bio tracer flux from biology (mmol/m^2/s)

      logical (kind=log_kind), intent(in) :: &
         tr_aero      ! aerosol flag

      ! local variables

      integer (kind=int_kind) :: &
         k, it        ! counting indices

      real (kind=dbl_kind) :: xtmp, dvssl, dvint

      ! aerosols
      if (tr_aero) then
         do it = 1, n_aero
            xtmp = (vsnon*(trcrn(nt_aero  +4*(it-1))     &
                         + trcrn(nt_aero+1+4*(it-1))))/dt
            dfaero_ocn(it) = dfaero_ocn(it) + xtmp
         enddo                 ! it
      endif ! tr_aero

      if (z_tracers) then
            dvssl  = min(p5*vsnon, hs_ssl*aicen)   !snow surface layer
            dvint  = vsnon- dvssl                  !snow interior

            do it = 1, nbtrcr
               xtmp = (trcrn(bio_index(it)+nblyr+1)*dvssl + &
                       trcrn(bio_index(it)+nblyr+2)*dvint)/dt
               dflux_bio(it) = dflux_bio(it) + xtmp
            enddo                 ! it

      endif ! z_tracers

      ! snow enthalpy tracer
      do k = 1, nslyr 
         xtmp = trcrn(nt_qsno+k-1) / dt &
              * vsnon/real(nslyr,kind=dbl_kind) ! < 0
         dfhocn = dfhocn + xtmp
         trcrn(nt_qsno+k-1) = c0
      enddo                  ! k

      ! snow volume
      xtmp = (rhos*vsnon) / dt
      dfresh = dfresh + xtmp
      vsnon = c0

      end subroutine zap_snow

!=======================================================================
   
      subroutine zap_snow_temperature(dt,         ncat,     &
                                      heat_capacity,        &
                                      nblyr,                &
                                      nslyr,      aicen,    &
                                      trcrn,      vsnon,    &
                                      dfresh,     dfhocn,   &
                                      dfaero_ocn, tr_aero,  &
                                      dflux_bio,  nbtrcr,   &
                                      n_aero,     ntrcr)

      use ice_colpkg_tracers, only: nt_qsno 
      use ice_therm_shared, only: Tmin

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nslyr , & ! number of snow layers
         n_aero, & ! number of aerosol tracers
         nbtrcr, & ! number of z-tracers in use
         nblyr , & ! number of bio  layers in ice
         ntrcr     ! number of tracers in use

      real (kind=dbl_kind), intent(in) :: &
         dt           ! time step

      logical (kind=log_kind), intent(in) :: &
         heat_capacity   ! if false, ice and snow have zero heat capacity

      real (kind=dbl_kind), dimension (:), intent(in) :: & 
         aicen        ! concentration of ice 

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         vsnon        ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn        ! ice tracers

      real (kind=dbl_kind), intent(inout) :: &
         dfresh   , & ! zapped fresh water flux (kg/m^2/s)
         dfhocn       ! zapped energy flux ( W/m^2)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dfaero_ocn   ! zapped aerosol flux   (kg/m^2/s)

      real (kind=dbl_kind), dimension (:),intent(inout) :: &
         dflux_bio    ! zapped biology flux  (mmol/m^2/s)

      logical (kind=log_kind), intent(in) :: &
         tr_aero      ! aerosol flag

      ! local variables

      integer (kind=int_kind) :: &
         n, k, it     ! counting indices

      real (kind=dbl_kind) :: &
         rnslyr   , & ! real(nslyr)
         hsn      , & ! snow thickness (m)
         zqsn     , & ! snow layer enthalpy (J m-2)
         zTsn     , & ! snow layer temperature (C)
         Tmax         ! maximum allowed snow temperature

      logical :: &
         l_zap        ! logical whether zap snow

      character(len=char_len_long) :: &
         warning ! warning message
      
      rnslyr = real(nslyr,kind=dbl_kind)
      
      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Determine cells to zap
      !-----------------------------------------------------------------

         l_zap = .false.

         if (aicen(n) > puny) then

         ! snow thickness
         hsn = vsnon(n) / aicen(n)

         ! check each snow layer - zap all if one is bad
         do k = 1, nslyr

            ! snow enthalpy and max temperature
            if (hsn > hs_min .and. heat_capacity) then
               ! zqsn < 0              
               zqsn = trcrn(nt_qsno+k-1,n)
               Tmax = -zqsn*puny*rnslyr / (rhos*cp_ice*vsnon(n))
            else
               zqsn = -rhos * Lfresh
               Tmax = puny
            endif
                     
            ! snow temperature
            zTsn = (Lfresh + zqsn/rhos)/cp_ice

            ! check for zapping
            if (zTsn < Tmin .or. zTsn > Tmax) then
               l_zap = .true.
               write(warning,*) "zap_snow_temperature: temperature out of bounds!"
               call add_warning(warning)
               write(warning,*) "k:"   , k
               call add_warning(warning)
               write(warning,*) "zTsn:", zTsn
               call add_warning(warning)
               write(warning,*) "Tmin:", Tmin
               call add_warning(warning)
               write(warning,*) "Tmax:", Tmax
               call add_warning(warning)
               write(warning,*) "zqsn:", zqsn
               call add_warning(warning)
            endif

         enddo ! k

         endif ! aicen > puny

      !-----------------------------------------------------------------
      ! Zap the cells
      !-----------------------------------------------------------------
         if (l_zap) &
            call zap_snow(dt,            nslyr,    &
                          trcrn(:,n),    vsnon(n), &
                          dfresh,        dfhocn,   &
                          dfaero_ocn,    tr_aero,  &
                          dflux_bio,     nbtrcr,   &
                          n_aero,        ntrcr,    &
                          aicen(n),      nblyr)

      enddo ! n

      end subroutine zap_snow_temperature

!=======================================================================
! Checks that the snow and ice energy in the zero layer thermodynamics
! model still agrees with the snow and ice volume.
! If there is an error, the model will abort.
! This subroutine is only called if heat_capacity = .false.
!
! author: Alison McLaren, Met Office
!         May 2010:  ECH replaced eicen, esnon with trcrn but did not test 
! the changes.  The loop below runs over n=1,ncat and I added loops 
! over k, making the test more stringent.

      subroutine zerolayer_check (ncat,        nilyr,      &
                                  nslyr,       aicen,      &
                                  vicen,       vsnon,      &
                                  trcrn,       l_stop,     &
                                  stop_label)

      use ice_colpkg_tracers, only: nt_qice, nt_qsno

      integer (kind=int_kind), intent(in) :: & 
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), dimension (:), intent(inout) :: & 
         aicen , & ! concentration of ice 
         vicen , & ! volume per unit area of ice          (m) 
         vsnon     ! volume per unit area of snow         (m) 

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers
      
      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      character (char_len), intent(out) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         k     , & ! vertical index
         n         ! category index

      real (kind=dbl_kind), parameter :: &
         max_error = puny*Lfresh*rhos ! max error in zero layer energy check
                                      ! (so max volume error = puny)

      real (kind=dbl_kind), dimension (ncat) :: &
         eicen     ! energy of melting for each ice layer (J/m^2) 
 
      real (kind=dbl_kind), dimension (ncat) :: &
         esnon     ! energy of melting for each snow layer (J/m^2) 

      logical (kind=log_kind) :: &
         ice_energy_correct  , & ! zero layer ice energy check
         snow_energy_correct     ! zero layer snow energy check

      real (kind=dbl_kind) :: &
         worka, workb

      character(len=char_len_long) :: &
         warning ! warning message
      
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.

      !----------------------------------------------------------------
      ! Calculate difference between ice and snow energies and the
      ! energy values derived from the ice and snow volumes
      !----------------------------------------------------------------

      ice_energy_correct  = .true.
      snow_energy_correct = .true.

      worka = c0
      workb = c0

      do n = 1, ncat

         eicen(n) = c0
         do k = 1, nilyr
            eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                     * vicen(n) / real(nilyr,kind=dbl_kind)
         enddo
         worka = eicen(n) + rhoi * Lfresh * vicen(n)
         esnon(n) = c0
         do k = 1, nslyr
            esnon(n) = esnon(n) + trcrn(nt_qsno+k-1,n) &
                     * vsnon(n) / real(nslyr,kind=dbl_kind)
         enddo
         workb = esnon(n) + rhos * Lfresh * vsnon(n)

         if(abs(worka) > max_error) ice_energy_correct = .false.
         if(abs(workb) > max_error) snow_energy_correct = .false.

      !----------------------------------------------------------------
      ! If there is a problem, abort with error message
      !----------------------------------------------------------------

         if (.not. ice_energy_correct) then

            if (abs(worka) > max_error) then
               l_stop = .true.
               stop_label = 'zerolayer check - wrong ice energy'
               write(warning,*) stop_label
               call add_warning(warning)
               write(warning,*) 'n:', n
               call add_warning(warning)
               write(warning,*) 'eicen =', eicen(n)
               call add_warning(warning)
               write(warning,*) 'error=',  worka
               call add_warning(warning)
               write(warning,*) 'vicen =', vicen(n)
               call add_warning(warning)
               write(warning,*) 'aicen =', aicen(n)
               call add_warning(warning)
            endif

         endif
         if (l_stop) return

         if (.not. snow_energy_correct) then

            if (abs(workb) > max_error) then
               l_stop = .true.
               stop_label = 'zerolayer check - wrong snow energy'
               write(warning,*) stop_label
               call add_warning(warning)
               write(warning,*) 'n:', n
               call add_warning(warning)
               write(warning,*) 'esnon =', esnon(n)
               call add_warning(warning)
               write(warning,*) 'error=',  workb
               call add_warning(warning)
               write(warning,*) 'vsnon =', vsnon(n)
               call add_warning(warning)
               write(warning,*) 'aicen =', aicen(n)
               call add_warning(warning)
               return
            endif

         endif

      enddo  ! ncat

      end subroutine zerolayer_check

!=======================================================================

      end module ice_itd

!=======================================================================









