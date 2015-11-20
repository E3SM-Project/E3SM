!=======================================================================
!BOP
!
! !MODULE: ice_itd - initialize and redistribute ice in the ITD
!
! !DESCRIPTION:
!
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
! !REVISION HISTORY:
!  SVN:$Id: ice_itd.F90 138 2008-07-08 20:39:37Z eclare $
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
!
! !INTERFACE:
!
      module ice_itd
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_domain_size
      use ice_constants
      use ice_fileunits
      use ice_exit
!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) :: &
         kitd        , & ! type of itd conversions
                         !   0 = delta function
                         !   1 = linear remap
         kcatbound   , & !   0 = old category boundary formula
                         !   1 = new formula giving round numbers
                         !   2 = WMO standard
         ilyr1 (ncat), & ! array position of top ice layer in each cat
         ilyrn (ncat), & ! array position of bottom ice layer in each cat
         slyr1 (ncat), & ! array position of top snow layer in each cat
         slyrn (ncat)    ! array position of bottom snow layer in each cat

      real (kind=dbl_kind), parameter :: &
         hi_min = p01    ! minimum ice thickness allowed (m)

      real (kind=dbl_kind) :: &
         hin_max(0:ncat) ! category limits (m)

      character (len=35) :: c_hi_range(ncat)

!-------------------------------------------------------------------
! a note regarding hi_min and hin_max(0):
! both represent a minimum ice thickness.  hin_max(0) is
! intended to be used for particular numerical implementations
! of category conversions in the ice thickness distribution.
! hi_min is a more general purpose parameter, but is specifically
! for maintaining stability in the thermodynamics.
! hin_max(0) = 0.1 m for the delta function itd
! hin_max(0) = 0.0 m for linear remapping
!
! Also note that the upper limit on the thickest category
! is only used for the linear remapping scheme
! and it is not a true upper limit on the thickness
!-------------------------------------------------------------------

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_itd - initalize area fraction and thickness boundaries for ITD
!
! !INTERFACE:
!
      subroutine init_itd
!
! !DESCRIPTION:
!
! Initialize area fraction and thickness boundaries for the itd model
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!          C. M. Bitz, UW
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
           n    ! thickness category index

      real (kind=dbl_kind) :: &
           cc1, cc2, cc3, & ! parameters for kcatbound = 0
           x1           , &
           rn           , & ! real(n)
           rncat        , & ! real(ncat)
           d1           , & ! parameters for kcatbound = 1 (m)
           d2

      real (kind=dbl_kind), dimension(5) :: wmo5 ! data for wmo itd
      real (kind=dbl_kind), dimension(6) :: wmo6 ! data for wmo itd
      real (kind=dbl_kind), dimension(7) :: wmo7 ! data for wmo itd

      character(len=8) :: c_hinmax1,c_hinmax2
      character(len=2) :: c_nc

      rncat = real(ncat, kind=dbl_kind)
      d1 = 3.0_dbl_kind / rncat
      d2 = 0.5_dbl_kind / rncat

      !-----------------------------------------------------------------
      ! Choose category boundaries based on one of three options.
      !
      ! The first formula (kcatbound = 0) was used in Lipscomb (2001) 
      !  and in CICE versions 3.0 and 3.1.
      !
      ! The second formula is more user-friendly in the sense that it
      !  is easy to obtain round numbers for category boundaries:
      !
      !    H(n) = n * [d1 + d2*(n-1)] 
      ! 
      ! Default values are d1 = 300/ncat, d2 = 50/ncat.
      ! For ncat = 5, boundaries in cm are 60, 140, 240, 360, which are 
      !  close to the standard values given by the first formula.
      ! For ncat = 10, boundaries in cm are 30, 70, 120, 180, 250, 330,
      !  420, 520, 630.    
      !
      ! The third option provides support for World Meteorological
      !  Organization classification based on thickness.  The full
      !  WMO thickness distribution is used if ncat = 7;  if ncat=5 
      !  or ncat = 6, some of the thinner categories are combined.
      ! For ncat = 5,  boundaries are         30, 70, 120, 200, >200 cm.
      ! For ncat = 6,  boundaries are     15, 30, 70, 120, 200, >200 cm.
      ! For ncat = 7,  boundaries are 10, 15, 30, 70, 120, 200, >200 cm.
      !-----------------------------------------------------------------

      if (kcatbound == 0) then   ! original scheme

         if (kitd == 1) then
            ! linear remapping itd category limits
            cc1 = c3/rncat
            cc2 = c15*cc1
            cc3 = c3

            hin_max(0) = c0     ! minimum ice thickness, m
         else
            ! delta function itd category limits
            cc1 = max(1.1_dbl_kind/rncat,c1*hi_min)
            cc2 = c25*cc1
            cc3 = 2.25_dbl_kind

            ! hin_max(0) should not be zero
            ! use some caution in making it less than 0.10
            hin_max(0) = hi_min ! minimum ice thickness, m
         endif                  ! kitd

         do n = 1, ncat
            x1 = real(n-1,kind=dbl_kind) / rncat
            hin_max(n) = hin_max(n-1) &
                       + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
         enddo

      elseif (kcatbound == 1) then  ! new scheme

         hin_max(0) = c0
         do n = 1, ncat
            rn = real(n, kind=dbl_kind)
            hin_max(n) = rn * (d1 + (rn-c1)*d2)
         enddo

      elseif (kcatbound == 2) then  ! WMO standard

        if (ncat == 5) then
         ! thinnest 3 categories combined
         data wmo5 / 0.30_dbl_kind, 0.70_dbl_kind, &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind  /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo5(n)
         enddo
       elseif (ncat == 6) then
         ! thinnest 2 categories combined
         data wmo6 / 0.15_dbl_kind, &
                    0.30_dbl_kind, 0.70_dbl_kind,  &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo6(n)
         enddo
       elseif (ncat == 7) then
         ! all thickness categories 
         data wmo7 / 0.10_dbl_kind, 0.15_dbl_kind, &
                    0.30_dbl_kind, 0.70_dbl_kind,  &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind  /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo7(n)
         enddo
       else
         write (nu_diag,*) 'kcatbound=3 (WMO) must have ncat=5, 6 or 7'
         stop
       endif

      endif ! kcatbound

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'hin_max(n-1) < Cat n < hin_max(n)'
         do n = 1, ncat
            write (nu_diag,*) hin_max(n-1),' < Cat ',n, ' < ',hin_max(n)
            ! Write integer n to character string
            write (c_nc, '(i2)') n    

            ! Write hin_max to character string
            write (c_hinmax1, '(f5.3)') hin_max(n-1)
            write (c_hinmax2, '(f5.3)') hin_max(n)

            ! Save character string to write to history file
            c_hi_range(n)=c_hinmax1//'m < hi Cat '//c_nc//' < '// &
                          c_hinmax2//'m'
         enddo
         write (nu_diag,*) ' '
      endif

      !-----------------------------------------------------------------
      ! vectors identifying first and last layer in each category
      !-----------------------------------------------------------------
      ilyr1(1) = 1                       ! if nilyr  = 4
      ilyrn(1) = nilyr                   !   ilyr1 = { 1,5,9 }
      do n = 2, ncat                     !   ilyrn = { 4,8,12} etc
         ilyr1(n) = ilyrn(n-1) + 1
         ilyrn(n) = ilyrn(n-1) + nilyr
      enddo

      slyr1(1) = 1
      slyrn(1) = nslyr
      do n = 2, ncat
         slyr1(n) = slyrn(n-1) + 1
         slyrn(n) = slyrn(n-1) + nslyr
      enddo

      end subroutine init_itd

!=======================================================================
!BOP
!
! !IROUTINE: aggregate - aggregate ice state variables
!
! !INTERFACE:
!
      subroutine aggregate (nx_block, ny_block, &
                            aicen,    trcrn,    &
                            vicen,    vsnon,    &
                            eicen,    esnon,    &
                            aice,     trcr,     &
                            vice,     vsno,     &
                            eice,     esno,     &
                            aice0,    tmask,    &
                            ntrcr,    trcr_depend)
!
! !DESCRIPTION:
!
! Aggregate ice state variables over thickness categories.
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW
!          W. H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ntrcr                 ! number of tracers in use

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(in) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(in) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(in) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         tmask     ! land/boundary mask, thickness (T-cell)

      integer (kind=int_kind), dimension (max_ntrcr), intent(in) :: &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block),  &
         intent(out) :: &
         aice  , & ! concentration of ice
         vice  , & ! volume per unit area of ice          (m)
         vsno  , & ! volume per unit area of snow         (m)
         eice  , & ! energy of melt. of ice           (J/m^2)
         esno  , & ! energy of melt. of snow layer    (J/m^2)
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr),  &
         intent(out) :: &
         trcr      ! ice tracers
!
!EOP
!
      integer (kind=int_kind) :: &
        icells                ! number of ocean/ice cells

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
        indxi, &              ! compressed indices in i/j directions
        indxj

      integer (kind=int_kind) :: &
        i, j, k, n, it, &
        ij                    ! combined i/j horizontal index

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
        atrcr      ! sum of aicen*trcrn or vicen*trcrn or vsnon*trcrn

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      icells = 0
      do j = 1, ny_block
      do i = 1, nx_block
         if (tmask(i,j)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif                  ! tmask

         aice0(i,j) = c1
         aice (i,j) = c0
         vice (i,j) = c0
         vsno (i,j) = c0
         eice (i,j) = c0
         esno (i,j) = c0
      enddo
      enddo


      allocate (atrcr(icells,ntrcr))

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      atrcr(:,:) = c0

      do n = 1, ncat

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            aice(i,j) = aice(i,j) + aicen(i,j,n)
            vice(i,j) = vice(i,j) + vicen(i,j,n)
            vsno(i,j) = vsno(i,j) + vsnon(i,j,n)
         enddo                  ! ij

         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then  ! ice area tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcr(ij,it) = atrcr(ij,it)  &
                                + trcrn(i,j,it,n)*aicen(i,j,n)
               enddo            ! ij

            elseif (trcr_depend(it) == 1) then  ! ice volume tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcr(ij,it) = atrcr(ij,it)  &
                                + trcrn(i,j,it,n)*vicen(i,j,n)
               enddo            ! ij

            elseif (trcr_depend(it) ==2) then ! snow volume tracer

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcr(ij,it) = atrcr(ij,it)  &
                                + trcrn(i,j,it,n)*vsnon(i,j,n)
               enddo            ! ij

            endif               ! trcr_depend
         enddo                  ! ntrcr

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               eice(i,j) = eice(i,j) + eicen(i,j,ilyr1(n)+k-1)
            enddo
         enddo                  ! nilyr

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               esno(i,j) = esno(i,j) + esnon(i,j,slyr1(n)+k-1)
            enddo
         enddo                  ! nslyr

      enddo                     ! ncat

      ! Open water fraction

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         aice0(i,j) = max (c1 - aice(i,j), c0)
      enddo                     ! ij

      ! Tracers

      call compute_tracers (nx_block,     ny_block,   &
                            icells,   indxi,   indxj, &
                            ntrcr,    trcr_depend,    &
                            atrcr, aice,    &
                            vice ,   vsno,  &
                            trcr)

      deallocate (atrcr)

      end subroutine aggregate

!=======================================================================
!BOP
!
! !IROUTINE: aggregate_area - aggregate ice area
!
! !INTERFACE:
!
      subroutine aggregate_area (nx_block, ny_block,        &
                                 aicen,    aice,     aice0)
!
! !DESCRIPTION:
!
! Aggregate ice area (but not other state variables) over thickness 
! categories.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          modified Jan 2004 by Clifford Chen, Fujitsu
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block  ! block dimensions

      real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
         aicen     ! concentration of ice

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         aice, &   ! concentration of ice
         aice0     ! concentration of open water
!
!EOP
!
      integer (kind=int_kind) :: i, j, n

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      aice(:,:) = c0

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aice(i,j) = aice(i,j) + aicen(i,j,n)
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      do j = 1, ny_block
      do i = 1, nx_block

         ! open water fraction
         aice0(i,j) = max (c1 - aice(i,j), c0)

      enddo                     ! i
      enddo                     ! j

      end subroutine aggregate_area

!=======================================================================
!BOP
!
! !IROUTINE: rebin - rebins thicknesses into defined categories
!
! !INTERFACE:
!
      subroutine rebin (nx_block, ny_block,        &
                        icells,   indxi,    indxj, &
                        ntrcr,    trcr_depend,     &
                        aicen,    trcrn,           &
                        vicen,    vsnon,           &
                        eicen,    esnon,           &
                        l_stop,                    &
                        istop,    jstop)
!
! !DESCRIPTION:
!
! Rebins thicknesses into defined categories
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of grid cells with ice
         ntrcr                 ! number of tracers in use

       integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj      ! compressed i/j indices

      integer (kind=int_kind), dimension (max_ntrcr), intent(in) :: &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice           (m)
         vsnon     ! volume per unit area of snow          (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(inout) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(inout) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop    ! indices of grid cell where model aborts
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j          , & ! horizontal indices
         n            , & ! category index
         ij                ! combined horizontal index

      logical (kind=log_kind) :: &
         shiftflag          ! = .true. if ice must be shifted

      integer (kind=int_kind), dimension (icells,ncat) :: &
         donor              ! donor category index

      real (kind=dbl_kind), dimension (icells,ncat) :: &
         daice          , & ! ice area transferred
         dvice          , & ! ice volume transferred
         hicen              ! ice thickness for each cat (m)

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

      do n = 1, ncat
         do ij = 1, icells       ! aice(i,j) > puny
            i = indxi(ij)
            j = indxj(ij)

            donor(ij,n) = 0
            daice(ij,n) = c0
            dvice(ij,n) = c0

      !-----------------------------------------------------------------
      ! Compute ice thickness.
      !-----------------------------------------------------------------
            if (aicen(i,j,n) > puny) then
               hicen(ij,n) = vicen(i,j,n) / aicen(i,j,n)
            else
               hicen(ij,n) = c0
            endif
         enddo                  ! ij
      enddo                     ! n

      !-----------------------------------------------------------------
      ! make sure thickness of cat 1 is at least hin_max(0)
      !-----------------------------------------------------------------
      do ij = 1, icells       ! aice(i,j) > puny
         i = indxi(ij)
         j = indxj(ij)

         if (aicen(i,j,1) > puny) then
            if (hicen(ij,1) <= hin_max(0) .and. hin_max(0) > c0 ) then
               aicen(i,j,1) = vicen(i,j,1) / hin_max(0)
               hicen(ij,1) = hin_max(0)
            endif
         endif
      enddo                     ! ij

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
         do ij = 1, icells       ! aice(i,j) > puny
            i = indxi(ij)
            j = indxj(ij)

            if (aicen(i,j,n) > puny .and. &
                hicen(ij,n) > hin_max(n)) then
               shiftflag = .true.
               donor(ij,n) = n
               daice(ij,n) = aicen(i,j,n)
               dvice(ij,n) = vicen(i,j,n)
            endif
         enddo                  ! ij

         if (shiftflag) then

      !-----------------------------------------------------------------
      ! shift ice between categories
      !-----------------------------------------------------------------
            call shift_ice (nx_block, ny_block,    &
                            indxi,    indxj,       &
                            icells,                &
                            ntrcr,    trcr_depend, &
                            aicen,    trcrn,       &
                            vicen,    vsnon,       &
                            eicen,    esnon,       &
                            hicen,    donor,       &
                            daice,    dvice,       &
                            l_stop,                &
                            istop,    jstop)

            if (l_stop) return

      !-----------------------------------------------------------------
      ! reset shift parameters
      !-----------------------------------------------------------------

         do ij = 1, icells       ! aice(i,j) > puny
            donor(ij,n) = 0
            daice(ij,n) = c0
            dvice(ij,n) = c0
         enddo

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
         do ij = 1, icells       ! aice(i,j) > puny
            i = indxi(ij)
            j = indxj(ij)

            if (aicen(i,j,n+1) > puny .and. &
                hicen(ij,n+1) <= hin_max(n)) then
               shiftflag = .true.
               donor(ij,n) = n+1
               daice(ij,n) = aicen(i,j,n+1)
               dvice(ij,n) = vicen(i,j,n+1)
            endif
         enddo                  ! ij

         if (shiftflag) then

      !-----------------------------------------------------------------
      ! shift ice between categories
      !-----------------------------------------------------------------
            call shift_ice (nx_block, ny_block,    &
                            indxi,    indxj,       &
                            icells,                &
                            ntrcr,    trcr_depend, &
                            aicen,    trcrn,       &
                            vicen,    vsnon,       &
                            eicen,    esnon,       &
                            hicen,    donor,       &
                            daice,    dvice,       &
                            l_stop,                &
                            istop,    jstop)

            if (l_stop) return

      !-----------------------------------------------------------------
      ! reset shift parameters
      !-----------------------------------------------------------------

         do ij = 1, icells       ! aice(i,j) > puny
            donor(ij,n) = 0
            daice(ij,n) = c0
            dvice(ij,n) = c0
         enddo

         endif                  ! shiftflag

      enddo                     ! n


      end subroutine rebin

!=======================================================================
!BOP
!
! !IROUTINE: reduce_area - reduce area when ice melts for special case ncat=1
!
! !INTERFACE:
!
      subroutine reduce_area (nx_block, ny_block, &
                              ilo, ihi, jlo, jhi, &
                              tmask,              &
                              aicen,     vicen,   &
                              aicen_init,vicen_init)
!
! !DESCRIPTION:
!
! Reduce area when ice melts for special case of ncat=1
!
! Use CSM 1.0-like method of reducing ice area
! when melting occurs: assume only half the ice volume
! change goes to thickness decrease, the other half
! to reduction in ice fraction
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW
! modified by: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         tmask     ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         aicen , & ! concentration of ice
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         aicen_init, & ! old ice area for category 1 (m)
         vicen_init    ! old ice volume for category 1 (m)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        ! horizontal indices

      real (kind=dbl_kind) :: &
         hi0     , & ! initial hi
         hi1     , & ! current hi
         dhi         ! hi1 - hi0

      do j = jlo, jhi
      do i = ilo, ihi
         if (tmask(i,j)) then

            hi0 = c0
            if (aicen_init(i,j) > c0) &
                hi0 = vicen_init(i,j) / aicen_init(i,j)

            hi1 = c0
            if (aicen(i,j) > c0) &
                hi1 = vicen(i,j) / aicen(i,j)

            ! make sure thickness of cat 1 is at least hin_max(0)
            if (hi1 <= hin_max(0) .and. hin_max(0) > c0 ) then
               aicen(i,j) = vicen(i,j) / hin_max(0)
               hi1 = hin_max(0)
            endif

            if (aicen(i,j) > c0) then
               dhi = hi1 - hi0
               if (dhi < c0) then
                  hi1  = vicen(i,j) / aicen(i,j)
                  aicen(i,j) = c2 * vicen(i,j) / (hi1 + hi0)
               endif
            endif

         endif                  ! tmask
      enddo                     ! i
      enddo                     ! j

      end subroutine reduce_area

!=======================================================================
!BOP
!
! !IROUTINE: shift_ice - shift ice across category boundaries
!
! !INTERFACE:
!
      subroutine shift_ice (nx_block, ny_block,    &
                            indxi,    indxj,       &
                            icells,                &
                            ntrcr,    trcr_depend, &
                            aicen,    trcrn,       &
                            vicen,    vsnon,       &
                            eicen,    esnon,       &
                            hicen,    donor,       &
                            daice,    dvice,       &
                            l_stop,                &
                            istop,    jstop)
!
! !DESCRIPTION:
!
! Shift ice across category boundaries, conserving area, volume, and
! energy.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of ocean/ice cells
         ntrcr                 ! number of tracers in use

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi             , & ! compressed indices in i/j directions
         indxj

      integer (kind=int_kind), dimension (max_ntrcr), intent(in) :: &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(inout) :: &
         eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(inout) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      ! NOTE: Third index of donor, daice, dvice should be ncat-1,
      !       except that compilers would have trouble when ncat = 1 
      integer (kind=int_kind), dimension(icells,ncat), &
         intent(in) :: &
         donor             ! donor category index

      real (kind=dbl_kind), dimension(icells,ncat), &
           intent(inout) :: &
         daice         , & ! ice area transferred across boundary
         dvice         , & ! ice volume transferred across boundary
         hicen             ! ice thickness for each cat        (m)

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop    ! indices of grid cell where model aborts
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, m       , & ! horizontal indices
         n             , & ! thickness category index
         nr            , & ! receiver category
         nd            , & ! donor category
         k             , & ! ice layer index
         it            , & ! tracer index
         ilo,ihi,jlo,jhi   ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(icells,max_ntrcr,ncat) :: &
         atrcrn            ! aicen*trcrn

      real (kind=dbl_kind) :: &
         dvsnow        , & ! snow volume transferred
         desnow        , & ! snow energy transferred
         deice         , & ! ice energy transferred
         datrcr            ! aicen*train transferred

      integer (kind=int_kind), dimension (icells) :: &
        indxii       , & ! compressed indices for i/j directions
        indxjj       , &
        indxij

      integer (kind=int_kind) :: &
        ishift      , & ! number of cells with ice to transfer
        ij              ! combined i/j horizontal index

      logical (kind=log_kind) :: &
        daice_negative     , & ! true if daice < -puny
        dvice_negative     , & ! true if dvice < -puny
        daice_greater_aicen, & ! true if daice > aicen
        dvice_greater_vicen    ! true if dvice > vicen

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka, &
         workb

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

      worka(:,:) = c0
      workb(:,:) = c0

      !-----------------------------------------------------------------
      ! Define variables equal to aicen*trcrn, vicen*trcrn, vsnon*trcrn
      !-----------------------------------------------------------------

      do n = 1, ncat
         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then ! ice area tracer
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcrn(ij,it,n) = aicen(i,j,n)*trcrn(i,j,it,n)
               enddo
            elseif (trcr_depend(it) ==1) then  ! ice volume tracer
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcrn(ij,it,n) = vicen(i,j,n)*trcrn(i,j,it,n)
               enddo
            elseif (trcr_depend(it) ==2) then  ! snow volume tracer
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcrn(ij,it,n) = vsnon(i,j,n)*trcrn(i,j,it,n)
               enddo
            endif
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Check for daice or dvice out of range, allowing for roundoff error
      !-----------------------------------------------------------------

      do n = 1, ncat-1

         daice_negative = .false.
         dvice_negative = .false.
         daice_greater_aicen = .false.
         dvice_greater_vicen = .false.


         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (donor(ij,n) > 0) then
               nd = donor(ij,n)

               if (daice(ij,n) < c0) then
                  if (daice(ij,n) > -puny*aicen(i,j,nd)) then
                     daice(ij,n) = c0 ! shift no ice
                     dvice(ij,n) = c0
                  else
                     daice_negative = .true.
                  endif
               endif
         
               if (dvice(ij,n) < c0) then
                  if (dvice(ij,n) > -puny*vicen(i,j,nd)) then   
                     daice(ij,n) = c0 ! shift no ice
                     dvice(ij,n) = c0
                  else
                     dvice_negative = .true.
                  endif
               endif

               if (daice(ij,n) > aicen(i,j,nd)*(c1-puny)) then
                  if (daice(ij,n) < aicen(i,j,nd)*(c1+puny)) then
                     daice(ij,n) = aicen(i,j,nd)
                     dvice(ij,n) = vicen(i,j,nd)
                  else
                     daice_greater_aicen = .true.
                  endif
               endif    

               if (dvice(ij,n) > vicen(i,j,nd)*(c1-puny)) then
                  if (dvice(ij,n) < vicen(i,j,nd)*(c1+puny)) then
                     daice(ij,n) = aicen(i,j,nd)
                     dvice(ij,n) = vicen(i,j,nd)
                  else
                     dvice_greater_vicen = .true.
                  endif
               endif

            endif               ! donor > 0
         enddo                  ! ij

      !-----------------------------------------------------------------
      ! error messages
      !-----------------------------------------------------------------

         if (daice_negative) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

               if (donor(ij,n) > 0 .and.  &
                   daice(ij,n) <= -puny*aicen(i,j,nd)) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'shift_ice: negative daice'
                  write(nu_diag,*) 'i, j:', i, j
                  write(nu_diag,*) 'boundary, donor cat:', n, nd
                  write(nu_diag,*) 'daice =', daice(ij,n)
                  write(nu_diag,*) 'dvice =', dvice(ij,n)
                  l_stop = .true.
                  istop = i
                  jstop = j
               endif
            enddo
         endif
         if (l_stop) return

         if (dvice_negative) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

               if (donor(ij,n) > 0 .and.  &
                   dvice(ij,n) <= -puny*vicen(i,j,nd)) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'shift_ice: negative dvice'
                  write(nu_diag,*) 'i, j:', i, j
                  write(nu_diag,*) 'boundary, donor cat:', n, nd
                  write(nu_diag,*) 'daice =', daice(ij,n)
                  write(nu_diag,*) 'dvice =', dvice(ij,n)
                  l_stop = .true.
                  istop = i
                  jstop = j
               endif
            enddo
         endif
         if (l_stop) return

         if (daice_greater_aicen) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

               if (donor(ij,n) > 0) then
                  nd = donor(ij,n)
                  if (daice(ij,n) >= aicen(i,j,nd)*(c1+puny)) then
                     write(nu_diag,*) ' '
                     write(nu_diag,*) 'shift_ice: daice > aicen'
                     write(nu_diag,*) 'i, j:', i, j
                     write(nu_diag,*) 'boundary, donor cat:', n, nd
                     write(nu_diag,*) 'daice =', daice(ij,n)
                     write(nu_diag,*) 'aicen =', aicen(i,j,nd)
                     l_stop = .true.
                     istop = i
                     jstop = j
                  endif
               endif
            enddo
         endif
         if (l_stop) return

         if (dvice_greater_vicen) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

               if (donor(ij,n) > 0) then
                  nd = donor(ij,n)
                  if (dvice(ij,n) >= vicen(i,j,nd)*(c1+puny)) then
                     write(nu_diag,*) ' '
                     write(nu_diag,*) 'shift_ice: dvice > vicen'
                     write(nu_diag,*) 'i, j:', i, j
                     write(nu_diag,*) 'boundary, donor cat:', n, nd
                     write(nu_diag,*) 'dvice =', dvice(ij,n)
                     write(nu_diag,*) 'vicen =', vicen(i,j,nd)
                     l_stop = .true.
                     istop = i
                     jstop = j
                  endif
               endif
            enddo
         endif
         if (l_stop) return

      !-----------------------------------------------------------------
      ! transfer volume and energy between categories
      !-----------------------------------------------------------------

         ishift = 0
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

           if (daice(ij,n) > c0) then ! daice(n) can be < puny
             ishift = ishift + 1
             indxii(ishift) = i
             indxjj(ishift) = j
             indxij(ishift) = ij
           endif   ! tmask
         enddo

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, ishift
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

            nd = donor(m,n)
            worka(i,j) = dvice(m,n) / vicen(i,j,nd)
            if (nd  ==  n) then
               nr = nd+1
            else                ! nd = n+1
               nr = n
            endif

            aicen(i,j,nd) = aicen(i,j,nd) - daice(m,n)
            aicen(i,j,nr) = aicen(i,j,nr) + daice(m,n)

            vicen(i,j,nd) = vicen(i,j,nd) - dvice(m,n)
            vicen(i,j,nr) = vicen(i,j,nr) + dvice(m,n)

            dvsnow = vsnon(i,j,nd) * worka(i,j)
            vsnon(i,j,nd) = vsnon(i,j,nd) - dvsnow
            vsnon(i,j,nr) = vsnon(i,j,nr) + dvsnow
            workb(i,j) = dvsnow

         enddo                  ! ij

         do it = 1, ntrcr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, ishift
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               nd = donor(m,n)
               if (nd == n) then
                  nr = nd+1
               else             ! nd = n+1
                  nr = n
               endif

               if (trcr_depend(it) == 0) then
                  datrcr = daice(m,n)*trcrn(i,j,it,nd)
               elseif (trcr_depend(it) == 1) then
                  datrcr = dvice(m,n)*trcrn(i,j,it,nd)
               elseif (trcr_depend(it) == 2) then
                  datrcr = workb(i,j)  *trcrn(i,j,it,nd)
               endif

               atrcrn(m,it,nd) = atrcrn(m,it,nd) - datrcr
               atrcrn(m,it,nr) = atrcrn(m,it,nr) + datrcr
            enddo               ! ij
         enddo                  ! ntrcr

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, ishift
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               nd = donor(m,n)
               if (nd == n) then
                  nr = nd+1
               else             ! nd = n+1
                  nr = n
               endif

               deice = eicen(i,j,ilyr1(nd)+k-1) * worka(i,j)
               eicen(i,j,ilyr1(nd)+k-1) = &
                    eicen(i,j,ilyr1(nd)+k-1) - deice
               eicen(i,j,ilyr1(nr)+k-1) = &
                    eicen(i,j,ilyr1(nr)+k-1) + deice
            enddo               ! ij
         enddo                  ! nilyr

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, ishift
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               nd = donor(m,n)
               if (nd == n) then
                  nr = nd+1
               else             ! nd = n+1
                  nr = n
               endif

               desnow = esnon(i,j,slyr1(nd)+k-1) * worka(i,j)
               esnon(i,j,slyr1(nd)+k-1) = &
                    esnon(i,j,slyr1(nd)+k-1) - desnow
               esnon(i,j,slyr1(nr)+k-1) = &
                    esnon(i,j,slyr1(nr)+k-1) + desnow
            enddo               ! ij
         enddo                  ! nslyr

      enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Update ice thickness and tracers
      !-----------------------------------------------------------------

      do n = 1, ncat

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (aicen(i,j,n) > puny) then
               hicen(ij,n)   = vicen (i,j,n)   / aicen(i,j,n)
            else
               hicen(ij,n)   = c0
            endif
         enddo

         call compute_tracers (nx_block,        ny_block,       &
                               icells,          indxi,   indxj, &
                               ntrcr,           trcr_depend,    &
                               atrcrn(:,:,n),   aicen(:,:,  n), &
                               vicen (:,:,  n), vsnon(:,:,  n), &
                               trcrn(:,:,:,n))

      enddo                     ! ncat

      end subroutine shift_ice

!=======================================================================
!BOP
!
! !IROUTINE: column_sum - sum field over all ice categories
!
! !INTERFACE:
!
      subroutine column_sum (nx_block, ny_block,       &
                             icells,   indxi,   indxj, &
                             nsum,                     &
                             xin,      xout)
!
! !DESCRIPTION:
!
! For each grid cell, sum field over all ice categories.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         nsum              , & ! number of categories/layers
         icells                ! number of ice/ocean grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj          ! compressed i/j indices

      real (kind=dbl_kind), dimension (nx_block,ny_block,nsum), &
           intent(in) :: &
           xin              ! input field

      real (kind=dbl_kind), dimension (icells), intent(out) :: &
           xout             ! output field
!
!EOP
!
      integer (kind=int_kind) :: &
           i, j, ij     , & ! horizontal indices
           n                ! category/layer index

      do ij = 1, icells
         xout(ij) = c0
      enddo

      do n = 1, nsum
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            xout(ij) = xout(ij) + xin(i,j,n)
         enddo                  ! ij
      enddo                     ! n

      end subroutine column_sum

!=======================================================================
!BOP
!
! !IROUTINE: column_conservation_check
!
! !INTERFACE:
!
      subroutine column_conservation_check (nx_block, ny_block,       &
                                            icells,   indxi,   indxj, &
                                            fieldid,                  &
                                            x1,       x2,             &
                                            max_err,  l_stop,         &
                                            istop,    jstop)
!
! !DESCRIPTION:
!
! For each physical grid cell, check that initial and final values
! of a conserved field are equal to within a small value.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of ice/ocean grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj     ! compressed i/j indices

      real (kind=dbl_kind), dimension(icells), intent(in) :: &
         x1            , & ! initial field
         x2                ! final field

      real (kind=dbl_kind), intent(in) :: &
         max_err           ! max allowed error

      character (len=char_len), intent(in) :: &
         fieldid           ! field identifier

      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, abort on return

      integer (kind=int_kind), intent(inout) :: &
         istop, jstop      ! indices of grid cell where model aborts
!
!EOP
!
      integer (kind=int_kind) :: &
         ij                    ! horizontal indices

      do ij = 1, icells
         if (abs (x2(ij)-x1(ij)) > max_err) then
            l_stop = .true.
            istop = indxi(ij)
            jstop = indxj(ij)

            write (nu_diag,*) ' '
            write (nu_diag,*) 'Conservation error: ', trim(fieldid)
            write (nu_diag,*) 'i, j =', istop, jstop
            write (nu_diag,*) 'Initial value =', x1(ij)
            write (nu_diag,*) 'Final value =',   x2(ij)
            write (nu_diag,*) 'Difference =', x2(ij) - x1(ij)
         endif
      enddo

      end subroutine column_conservation_check

!=======================================================================
!BOP
!
! !IROUTINE: compute_tracers - compute tracer fields
!
! !INTERFACE:
!
      subroutine compute_tracers (nx_block, ny_block,       &
                                  icells,   indxi,   indxj, &
                                  ntrcr,    trcr_depend,    &
                                  atrcrn,   aicen,          &
                                  vicen,    vsnon,          &
                                  trcrn)
!
! !DESCRIPTION:
!
! Compute tracer fields.
! Given atrcrn = aicen*trcrn (or vicen*trcrn, vsnon*trcrn), compute trcrn.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!          
! !USES:
!
      use ice_state, only: nt_Tsfc
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of ice/ocean grid cells
         ntrcr                 ! number of tracers in use

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj       ! compressed i/j indices

      integer (kind=int_kind), dimension (max_ntrcr), intent(in) :: &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (icells,ntrcr), &
         intent(in) :: &
         atrcrn    ! aicen*trcrn or vicen*trcrn or vsnon*trcrn

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr), &
         intent(out) :: &
         trcrn     ! ice tracers
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, it, ij       ! counting indices


      trcrn(:,:,:) = c0

      !-----------------------------------------------------------------
      ! Compute new tracers
      !-----------------------------------------------------------------

      do it = 1, ntrcr
         if (it == nt_Tsfc) then      ! surface temperature
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (aicen(i,j) > puny) then
                  trcrn(i,j,it) = atrcrn(ij,it) / aicen(i,j)
               else
                  trcrn(i,j,it) = Tocnfrz
               endif
            enddo

         elseif (trcr_depend(it) == 0) then ! ice area tracers
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (aicen(i,j) > puny) then
                  trcrn(i,j,it) = atrcrn(ij,it) / aicen(i,j)
               else
                  trcrn(i,j,it) = c0
               endif
            enddo

         elseif (trcr_depend(it) == 1) then ! ice volume tracers
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (vicen(i,j) > puny) then
                  trcrn(i,j,it) = atrcrn(ij,it) / vicen(i,j)
               else
                  trcrn(i,j,it) = c0
               endif
            enddo

         elseif (trcr_depend(it) == 2) then ! snow volume tracers
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (vsnon(i,j) > puny) then
                  trcrn(i,j,it) = atrcrn(ij,it) / vsnon(i,j)
               else
                  trcrn(i,j,it) = c0
               endif
            enddo

         endif                  ! trcr_depend
      enddo                     ! ntrcr

      end subroutine compute_tracers

!=======================================================================
!BOP
!
! !IROUTINE: cleanup_itd - rebin if needed, eliminate small ice areas,
!                          and aggregate over categories
!
! !INTERFACE:
!
      subroutine cleanup_itd (nx_block,    ny_block,   &
                              ilo, ihi,    jlo, jhi,   &
                              dt,          ntrcr,      &
                              aicen,       trcrn,      &
                              vicen,       vsnon,      &
                              eicen,       esnon,      &
                              aice0,       aice,       &
                              trcr_depend, fresh,      &
                              fsalt,       fhocn,      &
                              fsoot,       tr_aero,    &
                              heat_capacity, l_stop,   &
                              istop,         jstop,    &
                              limit_aice_in)
!
! !DESCRIPTION:
!
! Cleanup subroutine that rebins thickness categories if necessary,
!  eliminates very small ice areas while conserving mass and energy, 
!  aggregates state variables, and does a boundary call.  
! It is a good idea to call this subroutine after the thermodynamics
!  (thermo_vertical/thermo_itd) and again after the dynamics 
!  (evp/transport/ridging).
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions 
         ilo,ihi,jlo,jhi   , & ! beginning and end of physical domain
         ntrcr                 ! number of tracers in use
 
      real (kind=dbl_kind), intent(in) :: & 
         dt        ! time step 
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),  &
         intent(inout) :: & 
         aicen , & ! concentration of ice 
         vicen , & ! volume per unit area of ice          (m) 
         vsnon     ! volume per unit area of snow         (m) 
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat),  &
         intent(inout) :: & 
         trcrn     ! ice tracers 
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),  &
         intent(inout) :: & 
         eicen     ! energy of melting for each ice layer (J/m^2) 
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),  &
         intent(inout) :: & 
         esnon     ! energy of melting for each snow layer (J/m^2) 
 
      real (kind=dbl_kind), dimension (nx_block,ny_block),  &
         intent(inout) :: & 
         aice  , & ! total ice concentration
         aice0     ! concentration of open water 
     
      integer (kind=int_kind), dimension(max_ntrcr), intent(in) :: & 
         trcr_depend  ! tracer dependency information

      logical (kind=log_kind), intent(in) :: &
         tr_aero,      &
         heat_capacity   ! if false, ice and snow have zero heat capacity

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop ! indices of grid cell where model aborts

      ! ice-ocean fluxes (required for strict conservation)
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout), optional :: &
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean        (kg/m^2/s)
         fhocn        ! net heat flux to ocean     (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,n_aeromx), &
         intent(inout), optional :: &
         fsoot        ! soot flux to ocean        (kg/m^2/s)

      logical (kind=log_kind), intent(in), optional ::   &
         limit_aice_in      ! if false, allow aice to be out of bounds
                            ! may want to allow this for unit tests
!    
!EOP
!
      integer (kind=int_kind) :: &
         i, j             , & ! horizontal indices
         n                , & ! category index
         icells               ! number of grid cells with ice

       integer (kind=int_kind), dimension (nx_block*ny_block) :: &
         indxi, indxj      ! compressed i/j indices

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         dfresh   , & ! zapped fresh water flux (kg/m^2/s)
         dfsalt   , & ! zapped salt flux   (kg/m^2/s)
         dfhocn       ! zapped energy flux ( W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,n_aeromx) :: &
         dfsoot    ! zapped soot flux   (kg/m^2/s)

      logical (kind=log_kind) ::   &
         limit_aice         ! if true, check for aice out of bounds

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      if (present(limit_aice_in)) then
         limit_aice = limit_aice_in
      else
         limit_aice = .true.
      endif

      l_stop = .false.
      istop = 0
      jstop = 0

      !-----------------------------------------------------------------
      ! Compute total ice area.
      !-----------------------------------------------------------------

      call aggregate_area (nx_block, ny_block, &
                           aicen(:,:,:), &
                           aice,     aice0)


      if (limit_aice) then  ! check for aice out of bounds
      
         do j = jlo,jhi
         do i = ilo,ihi
            if (aice(i,j) > c1+puny .or. aice(i,j) < -puny) then
               l_stop = .true.
               istop = i
               jstop = j
            endif
         enddo
         enddo

         if (l_stop) then      ! area out of bounds
            i = istop
            j = jstop
            write(nu_diag,*) ' '
            write(nu_diag,*) 'aggregate ice area out of bounds'
            write(nu_diag,*) 'my_task, i, j, aice:', &
                              my_task, i, j, aice(i,j)
            do n = 1, ncat
               write(nu_diag,*) 'n, aicen:', n, aicen(i,j,n)
            enddo
            return
         endif                  ! l_stop
      endif                     ! limit_aice

      !-----------------------------------------------------------------
      ! Identify grid cells with ice.
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo,jhi
      do i = ilo,ihi
         if (aice(i,j) > puny) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Make sure ice in each category is within its thickness bounds.
      ! NOTE: The rebin subroutine is needed only in the rare cases
      !       when the linear_itd subroutine cannot transfer ice
      !       correctly (e.g., very fast ice growth).
      !-----------------------------------------------------------------

      call rebin (nx_block,     ny_block,       &
                  icells,       indxi, indxj,   &
                  ntrcr,        trcr_depend,    &
                  aicen(:,:,:), trcrn(:,:,:,:), &
                  vicen(:,:,:), vsnon(:,:,:),   &
                  eicen(:,:,:), esnon(:,:,:),   &
                  l_stop,                       &
                  istop,      jstop)

      if (l_stop) return

      !-----------------------------------------------------------------
      ! Zero out ice categories with very small areas.
      !-----------------------------------------------------------------

      if (limit_aice) then
         call zap_small_areas (nx_block,     ny_block,       &
                               ilo, ihi,     jlo, jhi,       &
                               dt,           ntrcr,          &
                               aice,         aice0,          &
                               aicen(:,:,:), trcrn(:,:,:,:), &
                               vicen(:,:,:), vsnon(:,:,:),   &
                               eicen(:,:,:), esnon(:,:,:),   &
                               dfresh,       dfsalt,         &
                               dfhocn,       dfsoot,         &
                               tr_aero,                      &
                               l_stop,                       &
                               istop,        jstop)
         if (l_stop) return
      endif   ! l_limit_aice

    !-------------------------------------------------------------------
    ! Update ice-ocean fluxes for strict conservation
    !-------------------------------------------------------------------

      if (present(fresh)) &
           fresh     (:,:) = fresh(:,:)      + dfresh(:,:) 
      if (present(fsalt)) &
           fsalt     (:,:) = fsalt(:,:)      + dfsalt(:,:)
      if (present(fhocn)) &
           fhocn     (:,:) = fhocn(:,:)      + dfhocn(:,:)
      if (present(fsoot)) &
           fsoot   (:,:,:) = fsoot(:,:,:)    + dfsoot(:,:,:)

      !----------------------------------------------------------------
      ! If using zero-layer model (no heat capacity), check that the 
      ! energy of snow and ice is correct. 
      !----------------------------------------------------------------

      if (.not. heat_capacity) then

         call zerolayer_check(nx_block,    ny_block,   &
                              icells,  indxi,   indxj, &
                              aicen,                   &
                              vicen,       vsnon,      &
                              eicen,       esnon,      &
                              l_stop,                  &
                              istop,       jstop)

      endif

      end subroutine cleanup_itd

!=======================================================================
!BOP
!
! !IROUTINE: zap_small_areas - eliminate very small ice areas
!
! !INTERFACE:
!
      subroutine zap_small_areas (nx_block, ny_block, &
                                  ilo, ihi, jlo, jhi, &
                                  dt,       ntrcr,    &
                                  aice,     aice0,    &
                                  aicen,    trcrn,    &
                                  vicen,    vsnon,    &
                                  eicen,    esnon,    &
                                  dfresh,   dfsalt,   &
                                  dfhocn,   dfsoot,   &
                                  tr_aero,            &
                                  l_stop,             &
                                  istop,    jstop)
!
! !DESCRIPTION:
!
! For each ice category in each grid cell, remove ice if the fractional
! area is less than puny.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
      use ice_state, only: nt_Tsfc, nt_aero
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi   , & ! beginning and end of physical domain
         ntrcr                 ! number of tracers in use

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         aice     , & ! total ice concentration
         aice0        ! concentration of open water

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
         intent(inout) :: &
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice          (m)
         vsnon        ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(inout) :: &
         eicen        ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(inout) :: &
         esnon        ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: &
         trcrn        ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         dfresh   , & ! zapped fresh water flux (kg/m^2/s)
         dfsalt   , & ! zapped salt flux   (kg/m^2/s)
         dfhocn       ! zapped energy flux ( W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,n_aeromx), &
         intent(out) :: &
         dfsoot    ! zapped soot flux   (kg/m^2/s)

      logical (kind=log_kind), intent(in) :: &
         tr_aero

      logical (kind=log_kind), intent(out) :: &
         l_stop   ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop ! indices of grid cell where model aborts
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j, n, k, it  , & ! counting indices
         icells         , & ! number of cells with ice to zap
         ij                 ! combined i/j horizontal index

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
        indxi       , & ! compressed indices for i/j directions
        indxj

      real (kind=dbl_kind) :: xtmp      ! temporary variable

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

      dfresh(:,:) = c0
      dfsalt(:,:) = c0
      dfhocn(:,:) = c0
      dfsoot(:,:,:) = c0

      !-----------------------------------------------------------------
      ! Zap categories with very small areas.
      !-----------------------------------------------------------------

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Count categories to be zapped.
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n) < -puny) then
               write (nu_diag,*) 'Zap ice: negative ice area'
               write (nu_diag,*) 'i, j, n, aicen =', &
                                  i, j, n, aicen(i,j,n)
               l_stop = .true.
               istop = i
               jstop = j
               return
            elseif ((aicen(i,j,n) >= -puny .and. aicen(i,j,n) < c0) .or. &
                    (aicen(i,j,n) > c0 .and. aicen(i,j,n) <= puny)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo
         enddo

      !-----------------------------------------------------------------
      ! Zap ice energy and use ocean heat to melt ice
      !-----------------------------------------------------------------

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               xtmp = eicen(i,j,ilyr1(n)+k-1) / dt ! < 0
               dfhocn(i,j) = dfhocn(i,j) + xtmp
               eicen(i,j,ilyr1(n)+k-1) = c0

            enddo               ! ij
         enddo                  ! k

      !-----------------------------------------------------------------
      ! Zap snow energy and use ocean heat to melt snow
      !-----------------------------------------------------------------

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               xtmp = esnon(i,j,slyr1(n)+k-1) / dt ! < 0
               dfhocn(i,j) = dfhocn(i,j) + xtmp
               esnon(i,j,slyr1(n)+k-1) = c0

            enddo               ! ij
         enddo                  ! k

      !-----------------------------------------------------------------
      ! Zap ice and snow volume, add water and salt to ocean
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            xtmp = (rhoi*vicen(i,j,n) + rhos*vsnon(i,j,n)) / dt
            dfresh(i,j) = dfresh(i,j) + xtmp

            xtmp = rhoi*vicen(i,j,n)*ice_ref_salinity*p001 / dt
            dfsalt(i,j) = dfsalt(i,j) + xtmp

            aice0(i,j) = aice0(i,j) + aicen(i,j,n)
            aicen(i,j,n) = c0
            vicen(i,j,n) = c0
            vsnon(i,j,n) = c0
            trcrn(i,j,nt_Tsfc,n) = Tocnfrz

         enddo                  ! ij

         if (tr_aero) then
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          do ij = 1, icells
           i = indxi(ij)
           j = indxj(ij)
           do it=1,n_aero
            xtmp &
              = (vsnon(i,j,n)*(trcrn(i,j,nt_aero  +4*(it-1),n)   &
                              +trcrn(i,j,nt_aero+1+4*(it-1),n))  &
              +  vicen(i,j,n)*(trcrn(i,j,nt_aero+2+4*(it-1),n)   &
                              +trcrn(i,j,nt_aero+3+4*(it-1),n))) &
              / dt
            dfsoot(i,j,it) = dfsoot(i,j,it) + xtmp
           enddo                 ! n
          enddo                  ! ij
         endif

      !-----------------------------------------------------------------
      ! Zap tracers
      !-----------------------------------------------------------------
         
         if (ntrcr >= 2) then
            do it = 1, ntrcr   ! this assumes nt_Tsfc = 1
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  trcrn(i,j,it,n) = c0
               enddo
            enddo
         endif

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Count cells with excess ice (aice > c1) due to roundoff errors.
      ! Zap a little ice in each category so that aice = c1.
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (aice(i,j) > (c1+puny)) then
            write (nu_diag,*) 'Zap ice: excess ice area'
            write (nu_diag,*) 'i, j, aice =', &
                               i, j, aice(i,j)
            l_stop = .true.
            istop = i
            jstop = j
            return
         elseif (aice(i,j) > c1 .and. aice(i,j) < (c1+puny)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo
      enddo

      do n = 1, ncat

      !----------------------------------------------------------------- 
      ! Zap ice energy and use ocean heat to melt ice 
      !----------------------------------------------------------------- 
       
         do k = 1, nilyr 
!DIR$ CONCURRENT !Cray 
!cdir nodep      !NEC 
!ocl novrec      !Fujitsu 
            do ij = 1, icells 
               i = indxi(ij) 
               j = indxj(ij) 
 
               xtmp = eicen(i,j,ilyr1(n)+k-1)  &
                    * (aice(i,j)-c1)/aice(i,j) / dt ! < 0 
               dfhocn(i,j) = dfhocn(i,j) + xtmp 
               eicen(i,j,ilyr1(n)+k-1) = eicen(i,j,ilyr1(n)+k-1) &
                                        * (c1/aice(i,j))
 
            enddo               ! ij 
         enddo                  ! k 
 
      !----------------------------------------------------------------- 
      ! Zap snow energy and use ocean heat to melt snow 
      !----------------------------------------------------------------- 

         do k = 1, nslyr 
!DIR$ CONCURRENT !Cray 
!cdir nodep      !NEC 
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij) 
               j = indxj(ij) 
 
               xtmp = esnon(i,j,slyr1(n)+k-1)  &
                    * (aice(i,j)-c1)/aice(i,j) / dt ! < 0 
               dfhocn(i,j) = dfhocn(i,j) + xtmp 
               esnon(i,j,slyr1(n)+k-1) = esnon(i,j,slyr1(n)+k-1) &
                                        *(c1/aice(i,j))
 
            enddo               ! ij
         enddo                  ! k
 
      !-----------------------------------------------------------------
      ! Zap ice and snow volume, add water and salt to ocean
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray 
!cdir nodep      !NEC 
!ocl novrec      !Fujitsu 
         do ij = 1, icells 
            i = indxi(ij) 
            j = indxj(ij) 
 
            xtmp = (rhoi*vicen(i,j,n) + rhos*vsnon(i,j,n)) &
                 * (aice(i,j)-c1)/aice(i,j) / dt 
            dfresh(i,j) = dfresh(i,j) + xtmp 
 
            xtmp = rhoi*vicen(i,j,n)*ice_ref_salinity*p001 &
                 * (aice(i,j)-c1)/aice(i,j) / dt
            dfsalt(i,j) = dfsalt(i,j) + xtmp 
 
            aicen(i,j,n) = aicen(i,j,n) * (c1/aice(i,j)) 
            vicen(i,j,n) = vicen(i,j,n) * (c1/aice(i,j)) 
            vsnon(i,j,n) = vsnon(i,j,n) * (c1/aice(i,j))
 
         enddo                  ! ij

      ! Note: Tracers are unchanged.

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         if (tr_aero) then
          do ij = 1, icells
           i = indxi(ij)
           j = indxj(ij)
           do it=1,n_aero
            xtmp &
              = (vsnon(i,j,n)*(trcrn(i,j,nt_aero  +4*(it-1),n)   &
                              +trcrn(i,j,nt_aero+1+4*(it-1),n))  &
              +  vicen(i,j,n)*(trcrn(i,j,nt_aero+2+4*(it-1),n)   &
                              +trcrn(i,j,nt_aero+3+4*(it-1),n))) &
              * (aice(i,j)-c1)/aice(i,j) / dt
            dfsoot(i,j,it) = dfsoot(i,j,it) + xtmp
           enddo                 ! n
          enddo                  ! ij
         endif

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Correct aice
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray 
!cdir nodep      !NEC 
!ocl novrec      !Fujitsu 
      do ij = 1, icells 
         i = indxi(ij) 
         j = indxj(ij) 
         aice(i,j) = c1
         aice0(i,j) = c0
      enddo

      end subroutine zap_small_areas

!=======================================================================
!BOP
!
! !IROUTINE: zerolayer_check - check that snow and ice energy is
!                         correct when using zero layer thermodynamics
!
! !INTERFACE:
!
      subroutine zerolayer_check (nx_block,    ny_block,   &
                                  icells,  indxi,   indxj, &
                                  aicen,                   &
                                  vicen,       vsnon,      &
                                  eicen,       esnon,      &
                                  l_stop,                  &
                                  istop,       jstop)
!
! !DESCRIPTION:
!
! Checks that the snow and ice energy in the zero layer thermodynamics
! model still agrees with the snow and ice volume.
! If there is an error, the model will abort.
! This subroutine is only called if heat_capacity = .false.
!
! !REVISION HISTORY:
!
! author: Alison McLaren, Met Office
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions 
         icells                ! number of grid cells with ice

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj      ! compressed i/j indices
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),  &
         intent(inout) :: & 
         aicen , & ! concentration of ice 
         vicen , & ! volume per unit area of ice          (m) 
         vsnon     ! volume per unit area of snow         (m) 
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),  &
         intent(in) :: & 
         eicen     ! energy of melting for each ice layer (J/m^2) 
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),  &
         intent(in) :: & 
         esnon     ! energy of melting for each snow layer (J/m^2) 
      
      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop ! indices of grid cell where model aborts
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j             , & ! horizontal indices
         n                , & ! category index
         ij                   ! combined horizontal index

      real (kind=dbl_kind), parameter :: &
         max_error = puny*Lfresh*rhos ! max error in zero layer energy check
                                      ! (so max volume error = puny)

      logical (kind=log_kind) :: &
         ice_energy_correct  , & ! zero layer ice energy check
         snow_energy_correct     ! zero layer snow energy check

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka, &
         workb

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

      worka(:,:) = c0
      workb(:,:) = c0

      !----------------------------------------------------------------
      ! Calculate difference between ice and snow energies and the
      ! energy values derived from the ice and snow volumes
      !----------------------------------------------------------------

      ice_energy_correct  = .true.
      snow_energy_correct = .true.

      do n=1,ncat

         do ij=1,icells
            i=indxi(ij)
            j=indxj(ij)

            worka(i,j) = eicen(i,j,n) + rhoi * Lfresh * vicen(i,j,n)
            workb(i,j) = esnon(i,j,n) + rhos * Lfresh * vsnon(i,j,n)

            if(abs(worka(i,j)) > max_error) then
               ice_energy_correct = .false.
            endif

            if(abs(workb(i,j)) > max_error) then
               snow_energy_correct = .false.
            endif
         enddo

      !----------------------------------------------------------------
      ! If there is a problem, abort with error message
      !----------------------------------------------------------------

         if (.not. ice_energy_correct) then

            do ij=1,icells
               i=indxi(ij)
               j=indxj(ij)

               if(abs(worka(i,j)) > max_error) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) &
                    'zerolayer check - wrong ice energy'
                  write(nu_diag,*) 'i, j, n:', i,j,n
                  write(nu_diag,*) 'eicen =', eicen(i,j,n)
                  write(nu_diag,*) 'error=',  worka(i,j)
                  write(nu_diag,*) 'vicen =', vicen(i,j,n)
                  write(nu_diag,*) 'aicen =', aicen(i,j,n)
                  l_stop = .true.
                  istop = i
                  jstop = j
               endif
            enddo

         endif
         if (l_stop) return

         if (.not. snow_energy_correct) then

            do ij=1,icells
               i=indxi(ij)
               j=indxj(ij)

               if(abs(workb(i,j)) > max_error) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) &
                    'zerolayer_check - wrong snow energy'
                  write(nu_diag,*) 'i, j, n:', i,j,n
                  write(nu_diag,*) 'esnon =', esnon(i,j,n)
                  write(nu_diag,*) 'error=',  workb(i,j)
                  write(nu_diag,*) 'vsnon =', vsnon(i,j,n)
                  write(nu_diag,*) 'aicen =', aicen(i,j,n)
                  l_stop = .true.
                  istop = i
                  jstop = j
                  return
               endif
            enddo

         endif

      enddo  ! ncat

      end subroutine zerolayer_check

!=======================================================================

      end module ice_itd

!=======================================================================









