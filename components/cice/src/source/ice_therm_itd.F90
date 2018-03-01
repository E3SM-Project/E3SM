!=======================================================================
!BOP
!
! !MODULE: ice_therm_itd - thermo calculations after call to coupler
!
! !DESCRIPTION:
!
! Thermo calculations after call to coupler, related to ITD:
! ice thickness redistribution, lateral growth and melting.
!
! NOTE: The thermodynamic calculation is split in two for load balancing.
!       First ice_therm_vertical computes vertical growth rates and coupler
!       fluxes.  Then ice_therm_itd does thermodynamic calculations not
!       needed for coupling.
!       
! !REVISION HISTORY:
!  SVN:$Id: ice_therm_itd.F90 48 2007-01-09 23:57:33Z eclare $
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb.  
! 2006: Streamlined for efficiency by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_therm_itd
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_domain_size
      use ice_constants
      use ice_fileunits
!
!EOP
!
      implicit none
      save

      real (kind=dbl_kind), parameter, private :: &
         hfrazilmin = 0.05_dbl_kind ! min thickness of new frazil ice (m)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: linear_itd - ITD scheme that shifts ice among categories
!
! !INTERFACE:
!
      subroutine linear_itd (nx_block,    ny_block,    & 
                             icells, indxi, indxj,     & 
                             ntrcr,       trcr_depend, & 
                             aicen_init,  vicen_init,  & 
                             aicen,       trcrn,       & 
                             vicen,       vsnon,       & 
                             eicen,       esnon,       & 
                             aice,        aice0,       & 
                             l_stop,                   & 
                             istop,       jstop)

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
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_itd, only: hin_max, hi_min, aggregate_area, shift_ice, & 
                         column_sum, column_conservation_check 
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

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen_init, & ! initial ice concentration (before vertical thermo)
         vicen_init    ! initial ice volume               (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         aicen  , & ! ice concentration
         vicen  , & ! volume per unit area of ice      (m)
         vsnon      ! volume per unit area of snow     (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(inout) :: &
         eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(inout) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         aice  , & ! concentration of ice
         aice0     ! concentration of open water

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop    ! indices of grid cell where model aborts 
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j         , & ! horizontal indices
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

      real (kind=dbl_kind), dimension(icells,0:ncat) :: &
         hbnew            ! new boundary locations

      real (kind=dbl_kind), dimension(icells) :: &
         work             ! temporary work array (for hbnew)

      integer (kind=int_kind), dimension(icells) :: &
         indxii, indxjj,& ! compressed i/j indices
         indxij           ! compressed 1D i/j indices

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         g0           , & ! constant coefficient in g(h)
         g1           , & ! linear coefficient in g(h)
         hL           , & ! left end of range over which g(h) > 0
         hR               ! right end of range over which g(h) > 0

      real (kind=dbl_kind), dimension(icells,ncat) :: &
         hicen        , & ! ice thickness for each cat     (m)
         hicen_init   , & ! initial ice thickness for each cat (pre-thermo)
         dhicen       , & ! thickness change for remapping (m)
         daice        , & ! ice area transferred across boundary
         dvice            ! ice volume transferred across boundary

      real (kind=dbl_kind), dimension(icells) :: &
         vice_init, vice_final, & ! ice volume summed over categories
         vsno_init, vsno_final, & ! snow volume summed over categories
         eice_init, eice_final, & ! ice energy summed over categories
         esno_init, esno_final    ! snow energy summed over categories

      ! NOTE: Third index of donor, daice, dvice should be ncat-1,
      !       except that compilers would have trouble when ncat = 1 
      integer (kind=int_kind), dimension(icells,ncat) :: &
         donor            ! donor category index

      logical (kind=log_kind), dimension(icells) :: &
         remap_flag       ! remap ITD if remap_flag(ij) is true

      character (len=char_len) :: &
         fieldid           ! field identifier

      logical (kind=log_kind), parameter :: &
         l_conservation_check = .true.   ! if true, check conservation
!         l_conservation_check = .false.   ! if true, check conservation
                                         ! (useful for debugging)

       integer (kind=int_kind) :: &
         iflag         , & ! number of grid cells with remap_flag = .true.
         ij, m             ! combined horizontal indices

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

      hin_max(ncat) = 999.9_dbl_kind ! arbitrary big number

      !-----------------------------------------------------------------
      ! Compute volume and energy sums that linear remapping should
      !  conserve.
      !-----------------------------------------------------------------

      if (l_conservation_check) then
         call column_sum (nx_block, ny_block,       &
                          icells,   indxi,   indxj, &
                          ncat,                     &
                          vicen,    vice_init)
         call column_sum (nx_block, ny_block,       &
                          icells,   indxi,   indxj, &
                          ncat,                     &
                          vsnon,    vsno_init)
         call column_sum (nx_block, ny_block,       &
                          icells,   indxi,   indxj, &
                          ntilyr,                   &
                          eicen,    eice_init)
         call column_sum (nx_block, ny_block,       &
                          icells,   indxi,   indxj, &
                          ntslyr,                   &
                          esnon,    esno_init)
      endif

      !-----------------------------------------------------------------
      ! Initialize remapping flag.
      ! Remapping is done wherever remap_flag = .true.
      ! In rare cases the category boundaries may shift too far for the
      !  remapping algorithm to work, and remap_flag is set to .false.
      ! In these cases the simpler 'rebin' subroutine will shift ice
      !  between categories if needed.
      !-----------------------------------------------------------------

      do ij = 1, icells
         remap_flag(ij) = .true.
      enddo

      !-----------------------------------------------------------------
      ! Compute thickness change in each category.
      !-----------------------------------------------------------------

      do n = 1, ncat
         do ij = 1, icells       ! aice(i,j) > puny
            i = indxi(ij)
            j = indxj(ij)

            if (aicen_init(i,j,n) > puny) then
               hicen_init (ij,n) = vicen_init(i,j,n) /  &
                                    aicen_init(i,j,n) 
            else
               hicen_init(ij,n) = c0
            endif               ! aicen_init > puny

            if (aicen(i,j,n) > puny) then
               hicen (ij,n) = vicen(i,j,n) / aicen(i,j,n) 
               dhicen(ij,n) = hicen(ij,n) - hicen_init(ij,n)
            else
               hicen (ij,n) = c0
               dhicen(ij,n) = c0
            endif               ! aicen > puny

         enddo                  ! ij
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Compute new category boundaries, hbnew, based on changes in
      ! ice thickness from vertical thermodynamics.
      !-----------------------------------------------------------------

      do ij = 1, icells       ! aice(i,j) > puny
         hbnew(ij,0) = hin_max(0)
      enddo

      do n = 1, ncat-1

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells       ! aice(i,j) > puny
            i = indxi(ij)
            j = indxj(ij)

            if (hicen_init(ij,n)   > puny .and. &
                hicen_init(ij,n+1) > puny) then
                 ! interpolate between adjacent category growth rates
               slope = (dhicen(ij,n+1) - dhicen(ij,n)) / &
                       (hicen_init(ij,n+1) - hicen_init(ij,n))
               hbnew(ij,n) = hin_max(n) + dhicen(ij,n) &
                            + slope * (hin_max(n) - hicen_init(ij,n))
            elseif (hicen_init(ij,n) > puny) then ! hicen_init(n+1)=0
               hbnew(ij,n) = hin_max(n) + dhicen(ij,n)
            elseif (hicen_init(ij,n+1) > puny) then ! hicen_init(n)=0
               hbnew(ij,n) = hin_max(n) + dhicen(ij,n+1)
            else
               hbnew(ij,n) = hin_max(n)
            endif

      !-----------------------------------------------------------------
      ! Check that each boundary lies between adjacent values of hicen.
      ! If not, set remap_flag = .false.
      ! Write diagnosis outputs if remap_flag was changed to false
      !-----------------------------------------------------------------

            if (aicen(i,j,n) > puny .and. &
                hicen(ij,n) >= hbnew(ij,n)) then
               remap_flag(ij) = .false.

                  write(nu_diag,*) my_task,':',i,j, &
                       'ITD: hicen(n) > hbnew(n)'
                  write(nu_diag,*) 'cat ',n
                  write(nu_diag,*) my_task,':',i,j, &
                       'hicen(n) =', hicen(ij,n)
                  write(nu_diag,*) my_task,':',i,j, &
                       'hbnew(n) =', hbnew(ij,n)

            elseif (aicen(i,j,n+1) > puny .and. &
                    hicen(ij,n+1) <= hbnew(ij,n)) then
               remap_flag(ij) = .false.

                  write(nu_diag,*) my_task,':',i,j, &
                       'ITD: hicen(n+1) < hbnew(n)'
                  write(nu_diag,*) 'cat ',n
                  write(nu_diag,*) my_task,':',i,j, &
                       'hicen(n+1) =', hicen(ij,n+1)
                  write(nu_diag,*) my_task,':',i,j, &
                       'hbnew(n) =', hbnew(ij,n)
            endif

      !-----------------------------------------------------------------
      ! Check that hbnew(n) lies between hin_max(n-1) and hin_max(n+1).
      ! If not, set remap_flag = .false.
      ! (In principle we could allow this, but it would make the code
      ! more complicated.)
      ! Write diagnosis outputs if remap_flag was changed to false
      !-----------------------------------------------------------------

            if (hbnew(ij,n) > hin_max(n+1)) then
               remap_flag(ij) = .false.

                  write(nu_diag,*) my_task,':',i,j, &
                       'ITD hbnew(n) > hin_max(n+1)'
                  write(nu_diag,*) 'cat ',n
                  write(nu_diag,*) my_task,':',i,j, &
                       'hbnew(n) =', hbnew(ij,n)
                  write(nu_diag,*) my_task,':',i,j, &
                       'hin_max(n+1) =', hin_max(n+1)
            endif

            if (hbnew(ij,n) < hin_max(n-1)) then
               remap_flag(ij) = .false.

                  write(nu_diag,*) my_task,':',i,j, &
                       'ITD: hbnew(n) < hin_max(n-1)'
                  write(nu_diag,*) 'cat ',n
                  write(nu_diag,*) my_task,':',i,j, &
                       'hbnew(n) =', hbnew(ij,n)
                  write(nu_diag,*) my_task,':',i,j, &
                       'hin_max(n-1) =', hin_max(n-1)
            endif

         enddo                  ! ij

      enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Compute hbnew(ncat)
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells         ! aice(i,j) > puny
         i = indxi(ij)
         j = indxj(ij)
         if (aicen(i,j,ncat) > puny) then
            hbnew(ij,ncat) = c3*hicen(ij,ncat) - c2*hbnew(ij,ncat-1)
         else
            hbnew(ij,ncat) = hin_max(ncat)
         endif
         hbnew(ij,ncat) = max(hbnew(ij,ncat),hin_max(ncat-1))
      enddo

      !-----------------------------------------------------------------
      ! Identify cells where the ITD is to be remapped
      !-----------------------------------------------------------------

      iflag = 0
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (remap_flag(ij)) then
            iflag = iflag + 1
            indxii(iflag) = i
            indxjj(iflag) = j
            indxij(iflag) = ij
         endif
      enddo

      allocate(g0(iflag,ncat))
      allocate(g1(iflag,ncat))
      allocate(hL(iflag,ncat))
      allocate(hR(iflag,ncat))

      !-----------------------------------------------------------------
      ! Compute g(h) for category 1 at start of time step
      ! (hicen = hicen_init)
      !-----------------------------------------------------------------

      do ij = 1, icells       ! aice(i,j) > puny
         work(ij) = hin_max(1)
      enddo

      call fit_line(nx_block,       ny_block,          &
                    iflag,          icells,            &
                    indxii,         indxjj,    indxij, &
                    aicen(:,:,1),   hicen_init(:,1),   &
                    hbnew(:,0),   work     (:),        &
                    g0   (:,1),   g1        (:,1),     &
                    hL   (:,1),   hR        (:,1))

      !-----------------------------------------------------------------
      ! Find area lost due to melting of thin (category 1) ice
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, iflag    ! remap_flag = .true.
         i = indxii(ij)
         j = indxjj(ij)
         m = indxij(ij)

         if (aicen(i,j,1) > puny) then

            dh0 = dhicen(m,1)

            if (dh0 < c0) then   ! remove area from category 1
               dh0 = min(-dh0,hin_max(1))   ! dh0 --> |dh0|

      !-----------------------------------------------------------------
      ! Integrate g(1) from 0 to dh0 to estimate area melted
      !-----------------------------------------------------------------

               ! right integration limit (left limit = 0)
               etamax = min(dh0,hR(ij,1)) - hL(ij,1)

               if (etamax > c0) then
                  x1 = etamax
                  x2 = p5 * etamax*etamax
                  da0 = g1(ij,1)*x2 + g0(ij,1)*x1 ! ice area removed

               ! constrain new thickness <= hicen_init
                  damax = aicen(i,j,1) &
                        * (c1-hicen(m,1)/hicen_init(m,1)) ! damax > 0
                  da0 = min (da0, damax)

               ! remove area, conserving volume
                  hicen(m,1) = hicen(m,1) &
                               * aicen(i,j,1) / (aicen(i,j,1)-da0)
                  aicen(i,j,1) = aicen(i,j,1) - da0
               endif            ! etamax > 0

            else                ! dh0 >= 0
               hbnew(m,0) = min(dh0,hin_max(1))  ! shift hbnew(0) to right
            endif

         endif                  ! aicen(i,j,1) > puny
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Compute g(h) for each ice thickness category.
      !-----------------------------------------------------------------

      do n = 1, ncat
         call fit_line(nx_block,       ny_block,          &
                       iflag,          icells,            &
                       indxii,         indxjj,    indxij, &
                       aicen(:,:,n),   hicen(:,n),        &
                       hbnew(:,n-1), hbnew(:,n),          &
                       g0   (:,n),   g1   (:,n),          &
                       hL   (:,n),   hR   (:,n))

      enddo

      !-----------------------------------------------------------------
      ! Compute area and volume to be shifted across each boundary.
      !-----------------------------------------------------------------

      do n = 1, ncat
         do ij = 1, icells
            donor(ij,n) = 0
            daice(ij,n) = c0
            dvice(ij,n) = c0
         enddo
      enddo

      do n = 1, ncat-1

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, iflag   ! remap_flag = .true.
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

            if (hbnew(m,n) > hin_max(n)) then ! transfer from n to n+1

               ! left and right integration limits in eta space
               etamin = max(hin_max(n), hL(ij,n)) - hL(ij,n)
               etamax = min(hbnew(m,n), hR(ij,n)) - hL(ij,n)
               donor(m,n) = n

            else             ! hbnew(n) <= hin_max(n); transfer from n+1 to n

               ! left and right integration limits in eta space
               etamin = c0
               etamax = min(hin_max(n), hR(ij,n+1)) - hL(ij,n+1)
               donor(m,n) = n+1

            endif            ! hbnew(n) > hin_max(n)

            if (etamax > etamin) then
               x1  = etamax - etamin
               wk1 = etamin*etamin
               wk2 = etamax*etamax
               x2  = p5 * (wk2 - wk1)
               wk1 = wk1*etamin
               wk2 = wk2*etamax
               x3  = p333 * (wk2 - wk1)
               nd  = donor(m,n)
               daice(m,n) = g1(ij,nd)*x2 + g0(ij,nd)*x1
               dvice(m,n) = g1(ij,nd)*x3 + g0(ij,nd)*x2 &
                            + daice(m,n)*hL(ij,nd)
            endif               ! etamax > etamin

            ! If daice or dvice is very small, shift no ice.

            nd = donor(m,n)

            if (daice(m,n) < aicen(i,j,nd)*puny) then
               daice(m,n) = c0
               dvice(m,n) = c0
               donor(m,n) = 0
            endif 

            if (dvice(m,n) < vicen(i,j,nd)*puny) then
               daice(m,n) = c0
               dvice(m,n) = c0
               donor(m,n) = 0
            endif

            ! If daice is close to aicen or dvice is close to vicen,
            ! shift entire category

            if (daice(m,n) > aicen(i,j,nd)*(c1-puny)) then
               daice(m,n) = aicen(i,j,nd)
               dvice(m,n) = vicen(i,j,nd)
            endif

            if (dvice(m,n) > vicen(i,j,nd)*(c1-puny)) then
               daice(m,n) = aicen(i,j,nd)
               dvice(m,n) = vicen(i,j,nd)
            endif

         enddo                  ! ij
      enddo                     ! boundaries, 1 to ncat-1

      deallocate(g0)
      deallocate(g1)
      deallocate(hL)
      deallocate(hR)

      !-----------------------------------------------------------------
      ! Shift ice between categories as necessary
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
      ! Make sure hice(i,j,1) >= minimum ice thickness hi_min.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, iflag          ! remap_flag = .true.
         i = indxii(ij)
         j = indxjj(ij)
         m = indxij(ij)
         if (hi_min > c0 .and. &
              aicen(i,j,1) > puny .and. hicen(m,1) < hi_min) then
            aicen(i,j,1) = aicen(i,j,1) * hicen(m,1)/hi_min
            hicen(m,1) = hi_min
         endif
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Update fractional ice area in each grid cell.
      !-----------------------------------------------------------------
      call aggregate_area (nx_block, ny_block, &
                           aicen,              &
                           aice,     aice0)

      if (l_stop) return

      !-----------------------------------------------------------------
      ! Check volume and energy conservation.
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         call column_sum (nx_block, ny_block,       &
                          icells,   indxi,   indxj, &
                          ncat,                     &
                          vicen,    vice_final)
         fieldid = 'vice, ITD remap'
         call column_conservation_check (nx_block,  ny_block,      &
                                         icells,   indxi,   indxj, &
                                         fieldid,                  &
                                         vice_init, vice_final,    &
                                         puny,      l_stop,        &
                                         istop,     jstop)
         if (l_stop) return

         call column_sum (nx_block, ny_block,       &
                          icells,   indxi,   indxj, &
                          ncat,                     &
                          vsnon,    vsno_final)
         fieldid = 'vsno, ITD remap'
         call column_conservation_check (nx_block,  ny_block,      &
                                         icells,   indxi,   indxj, &
                                         fieldid,                  &
                                         vsno_init, vsno_final,    &
                                         puny,      l_stop,        &
                                         istop,     jstop)
         if (l_stop) return

         call column_sum (nx_block, ny_block,       &
                          icells,   indxi,   indxj, &
                          ntilyr,                   &
                          eicen,    eice_final)
         fieldid = 'eice, ITD remap'
         call column_conservation_check (nx_block,   ny_block,     &
                                         icells,   indxi,   indxj, &
                                         fieldid,                  &
                                         eice_init,  eice_final,   &
                                         puny*Lfresh*rhoi,         &
                                         l_stop,                   &
                                         istop,     jstop)
         if (l_stop) return

         call column_sum (nx_block, ny_block,       &
                          icells,   indxi,   indxj, &
                          ntslyr,                   &
                          esnon,    esno_final)
         fieldid = 'esno, ITD remap'
         call column_conservation_check (nx_block,   ny_block,     &
                                         icells,   indxi,   indxj, &
                                         fieldid,                  &
                                         esno_init,  esno_final,   &
                                         puny*Lfresh*rhos,         &
                                         l_stop,                   &
                                         istop,     jstop)
         if (l_stop) return

      endif                     ! conservation check

      end subroutine linear_itd

!=======================================================================
!BOP
!
! !IROUTINE: fit_line - fit g(h) with a line using area, volume constraints
!
! !INTERFACE:
!

      subroutine fit_line (nx_block, ny_block,        &
                           iflag,   icells,           &
                           indxii,   indxjj,  indxij, &
                           aicen,    hice,            &
                           hbL,      hbR,             &
                           g0,       g1,              &
                           hL,       hR)
!
! !DESCRIPTION:
!
! Fit g(h) with a line, satisfying area and volume constraints.
! To reduce roundoff errors caused by large values of g0 and g1,
! we actually compute g(eta), where eta = h - hL, and hL is the
! left boundary.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of grid cells with ice
         iflag                 ! number of grid cells with remap_flag = .true.

       integer (kind=int_kind), dimension (icells), &
         intent(in) :: &
         indxii, indxjj, & ! compressed i/j indices (from iflag)
         indxij            ! compressed i/j indices (from icells)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aicen           ! concentration of ice

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hbL, hbR    , & ! left and right category boundaries
         hice            ! ice thickness

      real (kind=dbl_kind), dimension (iflag), intent(out):: &
         g0, g1      , & ! coefficients in linear equation for g(eta)
         hL          , & ! min value of range over which g(h) > 0
         hR              ! max value of range over which g(h) > 0
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j         , & ! horizontal indices
         ij, m           ! combined horizontal indices

      real  (kind=dbl_kind) :: &
         h13         , & ! hbL + 1/3 * (hbR - hbL)
         h23         , & ! hbL + 2/3 * (hbR - hbL)
         dhr         , & ! 1 / (hR - hL)
         wk1, wk2        ! temporary variables

      !-----------------------------------------------------------------
      ! Compute g0, g1, hL, and hR for each category to be remapped.
      !-----------------------------------------------------------------

       do ij = 1, iflag
         i = indxii(ij)
         j = indxjj(ij)
         m = indxij(ij)

         if (aicen(i,j) > puny .and. hbR(m) - hbL(m) > puny) then

         ! Initialize hL and hR

            hL(ij) = hbL(m)
            hR(ij) = hbR(m)

         ! Change hL or hR if hicen(n) falls outside central third of range

            h13 = p333 * (c2*hL(ij) + hR(ij))
            h23 = p333 * (hL(ij) + c2*hR(ij))
            if (hice(m) < h13) then
               hR(ij) = c3*hice(m) - c2*hL(ij)
            elseif (hice(m) > h23) then
               hL(ij) = c3*hice(m) - c2*hR(ij)
            endif

         ! Compute coefficients of g(eta) = g0 + g1*eta

            dhr = c1 / (hR(ij) - hL(ij))
            wk1 = c6 * aicen(i,j) * dhr
            wk2 = (hice(m) - hL(ij)) * dhr
            g0(ij) = wk1 * (p666 - wk2)
            g1(ij) = c2*dhr * wk1 * (wk2 - p5)

         else

            g0(ij) = c0
            g1(ij) = c0
            hL(ij) = c0
            hR(ij) = c0

         endif                  ! aicen > puny

      enddo                     ! ij

      end subroutine fit_line

!=======================================================================
!BOP
!
! !ROUTINE: add_new_ice - add frazil ice to ice thickness distribution
!
! !DESCRIPTION:
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
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine add_new_ice (nx_block,  ny_block, &
                              ntrcr,     icells,   &
                              indxi,     indxj,    &
                              tmask,     dt,       &
                              aicen,     trcrn,    &
                              vicen,     eicen,    &
                              aice0,     aice,     &
                              frzmlt,    frazil,   &
                              frz_onset, yday,     &
                              fresh,     fsalt,    &
                              Tf,        l_stop,   &
                              istop,     jstop)
!
! !USES:
!
      use ice_itd, only: hin_max, ilyr1, column_sum, &
                         column_conservation_check
      use ice_state, only: nt_Tsfc, nt_iage, nt_FY, nt_aero, &
                           tr_iage, tr_FY, tr_aero, &
                           nt_alvl, nt_vlvl, tr_lvl
      use ice_flux, only: update_ocn_f

! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         icells                ! number of ice/ocean grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj          ! compressed i/j indices

      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
         tmask     ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aice  , & ! total concentration of ice
         frzmlt, & ! freezing/melting potential (W/m^2)
         Tf        ! freezing temperature (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         aicen , & ! concentration of ice
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(inout) :: &
         eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         aice0     , & ! concentration of open water
         frazil    , & ! frazil ice growth        (m/step-->cm/day)
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt         ! salt flux to ocean (kg/m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout), optional :: &
         frz_onset ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(in), optional :: &
         yday      ! day of year

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop    ! indices of grid cell where model aborts
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j         , & ! horizontal indices
         n            , & ! ice category index
         k            , & ! ice layer index
         it               ! aerosol tracer index

      real (kind=dbl_kind), dimension (icells) :: &
         ai0new       , & ! area of new ice added to cat 1
         vi0new       , & ! volume of new ice added to cat 1
         hsurp        , & ! thickness of new ice added to each cat
         vlyr             ! ice layer volume

      real (kind=dbl_kind), dimension (icells) :: &
         vice_init, vice_final  ! ice volume summed over categories

      real (kind=dbl_kind) :: &
         fnew         , & ! heat flx to open water for new ice (W/m^2)
         hi0new       , & ! thickness of new ice
         hi0max       , & ! max allowed thickness of new ice
         qi0(nilyr)   , & ! frazil ice enthalpy
         qi0av        , & ! mean value of qi0 for new ice (J kg-1)
         vsurp        , & ! volume of new ice added to each cat
         area1        , & ! starting fractional area of existing ice
         vice1        , & ! starting volume of existing ice
         rnilyr       , & ! real(nilyr)
         dfresh       , & ! change in fresh
         dfsalt       , & ! change in fsalt
         vtmp

      integer (kind=int_kind) :: &
         jcells, kcells     , & ! grid cell counters
         ij, m                  ! combined i/j horizontal indices

      integer (kind=int_kind), dimension (icells) :: &
         indxij2,  indxij3  , & ! compressed i/j indices
         indxi2, indxj2     , &
         indxi3, indxj3

      character (len=char_len) :: &
         fieldid           ! field identifier

      l_stop = .false.
      istop = 0
      jstop = 0

      jcells = 0
      kcells = 0

      if (ncat > 1) then
         hi0max = hin_max(1)*0.9_dbl_kind  ! not too close to boundary
      else
         hi0max = bignum                   ! big number
      endif

      ! initial ice volume in each grid cell
      call column_sum (nx_block, ny_block,       &
                       icells,   indxi,   indxj, &
                       ncat,                     &
                       vicen,    vice_init)

      !-----------------------------------------------------------------
      ! Compute average enthalpy of new ice.
      !
      ! POP assumes new ice is fresh.  Otherwise, it would be better
      ! to do something like this:
      !  qi0(i,j,k) = -rhoi * (cp_ice*(Tmlt(k)-Tf(i,j))
      !             + Lfresh*(1.-Tmlt(k)/Tf(i,j)) - cp_ocn*Tmlt(k))
      !-----------------------------------------------------------------

      rnilyr = real(nilyr,kind=dbl_kind)
      qi0av = c0
      do k = 1, nilyr
         qi0(k) = -rhoi*Lfresh  ! note sign convention, qi < 0
         qi0av  = qi0av + qi0(k)
      enddo
      qi0av = qi0av/rnilyr

      !-----------------------------------------------------------------
      ! Compute the volume, area, and thickness of new ice.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         fnew = max (frzmlt(i,j), c0)   ! fnew > 0 iff frzmlt > 0
         vi0new(ij) = -fnew*dt / qi0av ! note sign convention, qi < 0

         ! increment ice volume
         vice_init(ij) = vice_init(ij) + vi0new(ij)

         ! history diagnostics
         frazil(i,j) = vi0new(ij)

         if (present(frz_onset) .and. present(yday)) then
            if (frazil(i,j) > puny .and. frz_onset(i,j) < puny) &
                 frz_onset(i,j) = yday
         endif

      !-----------------------------------------------------------------
      ! Update fresh water and salt fluxes.
      !
      ! NOTE: POP assumes fresh water and salt flux due to frzmlt > 0
      !       is NOT included in fluxes fresh and fsalt.
      !-----------------------------------------------------------------

         if (update_ocn_f) then
            dfresh = -rhoi*vi0new(ij)/dt
            dfsalt = ice_ref_salinity*p001*dfresh
            fresh(i,j)      = fresh(i,j)      + dfresh
            fsalt(i,j)      = fsalt(i,j)      + dfsalt
         endif

      !-----------------------------------------------------------------
      ! Decide how to distribute the new ice.
      !-----------------------------------------------------------------

         hsurp(ij)  = c0
         ai0new(ij) = c0

         if (vi0new(ij) > c0) then

            ! new ice area and thickness
            ! hin_max(0) < new ice thickness < hin_max(1)
            if (aice0(i,j) > puny) then
               hi0new = max(vi0new(ij)/aice0(i,j), hfrazilmin)
               if (hi0new > hi0max .and. aice0(i,j)+puny < c1) then
                  ! distribute excess volume over all categories (below)
                  hi0new = hi0max
                  ai0new(ij) = aice0(i,j)
                  vsurp       = vi0new(ij) - ai0new(ij)*hi0new
                  hsurp(ij)  = vsurp / aice(i,j)
                  vi0new(ij) = ai0new(ij)*hi0new
               else
                  ! put ice in a single category, with hsurp = 0
                  ai0new(ij) = vi0new(ij)/hi0new
               endif
            else                ! aice0 < puny
               hsurp(ij) = vi0new(ij)/aice(i,j) ! new thickness in each cat
               vi0new(ij) = c0
            endif               ! aice0 > puny
         endif                  ! vi0new > puny

      !-----------------------------------------------------------------
      ! Identify grid cells receiving new ice.
      !-----------------------------------------------------------------

         i = indxi(ij)
         j = indxj(ij)

         if (vi0new(ij) > c0) then  ! add ice to category 1
            jcells = jcells + 1
            indxi2(jcells) = i
            indxj2(jcells) = j
            indxij2(jcells) = ij
         endif

         if (hsurp(ij) > c0) then   ! add ice to all categories
            kcells = kcells + 1
            indxi3(kcells) = i
            indxj3(kcells) = j
            indxij3(kcells) = ij
         endif

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Distribute excess ice volume among ice categories by increasing
      ! ice thickness, leaving ice area unchanged.
      !
      ! NOTE: If new ice contains globally conserved tracers
      !       (e.g., isotopes from seawater), code must be added here.
      !-----------------------------------------------------------------

      do n = 1, ncat

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, kcells
            i = indxi3(ij)
            j = indxj3(ij)
            m = indxij3(ij)

            vsurp = hsurp(m) * aicen(i,j,n)

            ! update ice age due to freezing (new ice age = dt)
            vtmp = vicen(i,j,n) + vsurp
            if (vtmp > puny) then

            if (tr_iage) trcrn(i,j,nt_iage,n)  &
                   = (trcrn(i,j,nt_iage,n)*vicen(i,j,n) + dt*vsurp) &
                   / vtmp

            if (tr_lvl) &
                trcrn(i,j,nt_vlvl,n) = &
               (trcrn(i,j,nt_vlvl,n)*vicen(i,j,n) + &
                trcrn(i,j,nt_alvl,n)*vsurp) / vtmp

            if (tr_aero) then
             do it=1,n_aero
               trcrn(i,j,nt_aero+2+4*(it-1),n)  &
                   = trcrn(i,j,nt_aero+2+4*(it-1),n)*vicen(i,j,n) &
                   / vtmp
               trcrn(i,j,nt_aero+3+4*(it-1),n)  &
                   = trcrn(i,j,nt_aero+3+4*(it-1),n)*vicen(i,j,n) &
                   / vtmp
             enddo
            endif

            endif

            ! update category volumes
            vicen(i,j,n) = vicen(i,j,n) + vsurp
            vlyr(m) = vsurp/rnilyr

         enddo                  ! ij

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, kcells
               i = indxi3(ij)
               j = indxj3(ij)
               m = indxij3(ij)

               eicen(i,j,ilyr1(n)+k-1) = &
                    eicen(i,j,ilyr1(n)+k-1) + qi0(k)*vlyr(m)
            enddo               ! ij
         enddo                  ! k

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Combine new ice grown in open water with category 1 ice.
      ! Assume that vsnon and esnon are unchanged.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, jcells
         i = indxi2(ij)
         j = indxj2(ij)
         m = indxij2(ij)

         area1 = aicen(i,j,1)   ! save
         vice1 = vicen(i,j,1)   ! save
         aicen(i,j,1) = aicen(i,j,1) + ai0new(m)
         aice0(i,j)   = aice0(i,j)   - ai0new(m)
         vicen(i,j,1) = vicen(i,j,1) + vi0new(m)
         trcrn(i,j,nt_Tsfc,1) = (Tf(i,j)*ai0new(m) + trcrn(i,j,nt_Tsfc,1)*area1) &
                      / aicen(i,j,1)
         trcrn(i,j,nt_Tsfc,1) = min (trcrn(i,j,nt_Tsfc,1), c0)

         if (vicen(i,j,1) > puny) then

           if (tr_iage) trcrn(i,j,nt_iage,1) = &
              (trcrn(i,j,nt_iage,1)*vice1 + dt*vi0new(m))/vicen(i,j,1)

           if (tr_aero) then
              do it=1,n_aero
                trcrn(i,j,nt_aero+2+4*(it-1),1) = &
                  trcrn(i,j,nt_aero+2+4*(it-1),1)*vice1/vicen(i,j,1)
                trcrn(i,j,nt_aero+3+4*(it-1),1) = &
                  trcrn(i,j,nt_aero+3+4*(it-1),1)*vice1/vicen(i,j,1)
              enddo
           endif

         endif

         if (tr_lvl .and. aicen(i,j,1) > puny) then
             trcrn(i,j,nt_alvl,1) = &
            (trcrn(i,j,nt_alvl,1)*area1 + ai0new(m))/aicen(i,j,1)
             trcrn(i,j,nt_vlvl,1) = &
            (trcrn(i,j,nt_vlvl,1)*vice1 + vi0new(m))/vicen(i,j,1)
         endif

         if (tr_FY) &
             trcrn(i,j,nt_FY,1) = (ai0new(m) + trcrn(i,j,nt_FY,1)*area1) &
                                / aicen(i,j,1)

         vlyr(m)    = vi0new(m) / rnilyr
      enddo                     ! ij

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, jcells
            i = indxi2(ij)
            j = indxj2(ij)
            m = indxij2(ij)
            eicen(i,j,k) = eicen(i,j,k) + qi0(k)*vlyr(m)
         enddo
      enddo

      call column_sum (nx_block, ny_block,       &
                       icells,   indxi,   indxj, &
                       ncat,                     &
                       vicen,    vice_final)

      fieldid = 'vice, add_new_ice'
      call column_conservation_check (nx_block,  ny_block,      &
                                      icells,   indxi,   indxj, &
                                      fieldid,                  &
                                      vice_init, vice_final,    &
                                      puny,      l_stop,        &
                                      istop,     jstop)
      if (l_stop) return

      end subroutine add_new_ice

!=======================================================================
!BOP
!
! !ROUTINE: lateral_melt - melt ice laterally
!
! !DESCRIPTION:
!
! Given the fraction of ice melting laterally in each grid cell
!  (computed in subroutine frzmlt_bottom_lateral), melt ice.
!
! !REVISION HISTORY:
!
! author: C. M. Bitz, UW
! 2003:   Modified by William H. Lipscomb and Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine lateral_melt (nx_block,   ny_block,   &
                               ilo, ihi,   jlo, jhi,   &
                               dt,                     &
                               fresh,      fsalt,      &
                               fhocn,      fsoot,      &
                               rside,      meltl,      &
                               aicen,      vicen,      &
                               vsnon,      eicen,      &
                               esnon,      trcrn)
!
! !USES:
!
      use ice_itd, only: ilyr1, slyr1
      use ice_state, only: nt_aero, tr_aero
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume per unit area of ice          (m)
         vsnon       ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(inout) :: &
         eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(inout) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(in) :: &
         trcrn     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         rside     ! fraction of ice that melts laterally

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt     , & ! salt flux to ocean (kg/m^2/s)
         fhocn     , & ! net heat flux to ocean (W/m^2)
         meltl         ! lateral ice melt         (m/step-->cm/day)

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_aeromx), &
         intent(inout) :: &
         fsoot      ! 
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         n           , & ! thickness category index
         k           , & ! layer index
         ij          , & ! horizontal index, combines i and j loops
         icells      , & ! number of cells with aice > puny
         it              ! tracer index for aerosols

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with aice > puny

      real (kind=dbl_kind) :: &
         dfhocn  , & ! change in fhocn
         dfresh  , & ! change in fresh
         dfsalt      ! change in fsalt

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify grid cells with lateral melting.
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (rside(i,j) > c0) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo                  ! i
         enddo                  ! j

      !-----------------------------------------------------------------
      ! Melt the ice and increment fluxes.
      !-----------------------------------------------------------------

         if (tr_aero) then
          do k=1,n_aero
           do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
             fsoot(i,j,k)      = fsoot(i,j,k) &
               + (vsnon(i,j,n) &
               *(trcrn(i,j,nt_aero  +4*(k-1),n)   &
                +trcrn(i,j,nt_aero+1+4*(k-1),n))  &
               +  vicen(i,j,n) &
               *(trcrn(i,j,nt_aero+2+4*(k-1),n)   &
                +trcrn(i,j,nt_aero+3+4*(k-1),n))) &
               * rside(i,j) / dt
            enddo
          enddo
         endif

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            ! fluxes to coupler
            ! dfresh > 0, dfsalt > 0

            dfresh = (rhos*vsnon(i,j,n) + rhoi*vicen(i,j,n)) &
                   * rside(i,j) / dt
            dfsalt = rhoi*vicen(i,j,n)*ice_ref_salinity*p001 &
                   * rside(i,j) / dt

            fresh(i,j)      = fresh(i,j)      + dfresh
            fsalt(i,j)      = fsalt(i,j)      + dfsalt

            ! history diagnostics
            meltl(i,j) = meltl(i,j) + vicen(i,j,n)*rside(i,j)

            ! state variables
            aicen(i,j,n) = aicen(i,j,n) * (c1 - rside(i,j))
            vicen(i,j,n) = vicen(i,j,n) * (c1 - rside(i,j))
            vsnon(i,j,n) = vsnon(i,j,n) * (c1 - rside(i,j))

         enddo                  ! ij

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               ! heat flux to coupler for ice melt (dfhocn < 0)

               dfhocn = eicen(i,j,ilyr1(n)+k-1)*rside(i,j) / dt
               fhocn(i,j)      = fhocn(i,j)      + dfhocn

               ! ice energy
               eicen(i,j,ilyr1(n)+k-1) = eicen(i,j,ilyr1(n)+k-1) &
                                       * (c1 - rside(i,j))
            enddo               ! ij
         enddo                  ! nilyr

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               ! heat flux to coupler for snow melt (dfhocn < 0)

               dfhocn = esnon(i,j,slyr1(n)+k-1)*rside(i,j) / dt
               fhocn(i,j)      = fhocn(i,j)      + dfhocn

               ! snow energy
               esnon(i,j,slyr1(n)+k-1) = esnon(i,j,slyr1(n)+k-1) &
                                       * (c1 - rside(i,j))
            enddo               ! ij
         enddo                  ! nslyr

      enddo  ! n

      end subroutine lateral_melt

!=======================================================================

      end module ice_therm_itd

!=======================================================================
