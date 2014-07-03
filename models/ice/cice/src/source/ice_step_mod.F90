!=======================================================================
!
!BOP
!
! !MODULE: ice_step_mod
!
! !DESCRIPTION:
!
!  Contains CICE component driver routines common to all drivers.
!
! !REVISION HISTORY:
!  SVN:$Id: $
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2008 ECH: created module by moving subroutines from drivers/cice4/
!
! !INTERFACE:
!
      module ice_step_mod
!
! !USES:
!
      use ice_atmo
      use ice_calendar
      use ice_communicate
      use ice_diagnostics
      use ice_domain
      use ice_dyn_evp
      use ice_fileunits
      use ice_flux
      use ice_grid
      use ice_history
      use ice_restart
      use ice_itd
      use ice_kinds_mod
      use ice_mechred
      use ice_ocean
      use ice_orbital
      use ice_shortwave
      use ice_state
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use perf_mod, only: t_startf, t_stopf, t_barrierf

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: step_therm2, step_dynamics, &
                prep_radiation, step_radiation
!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: prep_radiation - step pre-thermo radiation
!
! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
! authors: Mariana Vertenstein, NCAR
!
! !INTERFACE:

      subroutine prep_radiation(dt)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j,n,iblk    ! block index

      if (calc_Tsfc) then

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call prep_radiation_iblk(dt, iblk)
         end do
         !$OMP END PARALLEL DO

      else    ! .not. calc_Tsfc

         ! Initialize for safety
         do iblk = 1, nblocks
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            fswsfcn(i,j,n,iblk) = c0
            fswintn(i,j,n,iblk) = c0
            fswthrun(i,j,n,iblk) = c0
         enddo   ! i
         enddo   ! j
         enddo   ! ncat
            Iswabsn(:,:,:,iblk) = c0
            Sswabsn(:,:,:,iblk) = c0
         enddo   ! iblk

      endif    ! calc_Tsfc

      end subroutine prep_radiation

!=======================================================================
!BOP
!
! !ROUTINE: prep_radiation_iblk - step pre-thermo radiation
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: David A. Bailey, NCAR
!
! !INTERFACE:

      subroutine prep_radiation_iblk (dt, iblk)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n           , & ! thickness category index
         il1, il2    , & ! ice layer indices for eice
         sl1, sl2        ! snow layer indices for esno

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      ! snow variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsn             ! snow horizontal fraction
      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) :: &
         rhosnwn     , & ! snow density (kg/m3)
         rsnwn           ! snow grain radius (micro-meters)

      ! pond variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fpn         , & ! pond fraction
         hpn             ! pond depth (m)

      real (kind=dbl_kind) :: netsw, netsw_old, ar

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      l_stop = .false.

      fswfac(:,:,iblk) = c1

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Compute netsw scaling factor (new netsw / old netsw)
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         if (aice(i,j,iblk) > c0 .and. scale_factor(i,j,iblk) > puny) then
            netsw = swvdr(i,j,iblk)*(c1 - alvdr_gbm(i,j,iblk)) &
                  + swvdf(i,j,iblk)*(c1 - alvdf_gbm(i,j,iblk)) &
                  + swidr(i,j,iblk)*(c1 - alidr_gbm(i,j,iblk)) &
                  + swidf(i,j,iblk)*(c1 - alidf_gbm(i,j,iblk))
            scale_factor(i,j,iblk) = netsw / scale_factor(i,j,iblk)
         else
            scale_factor(i,j,iblk) = c1
         endif
         fswfac(i,j,iblk) = scale_factor(i,j,iblk) ! for history
      enddo               ! i
      enddo               ! j

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify cells with nonzero ice area
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n,iblk) > puny) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

      !-----------------------------------------------------------------
      ! Scale absorbed solar radiation for change in net shortwave
      !-----------------------------------------------------------------

         il1 = ilyr1(n)
         il2 = ilyrn(n)
         sl1 = slyr1(n)
         sl2 = slyrn(n)

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            fswsfcn(i,j,n,iblk) = scale_factor(i,j,iblk)*fswsfcn (i,j,n,iblk)
            fswintn(i,j,n,iblk) = scale_factor(i,j,iblk)*fswintn (i,j,n,iblk)
            fswthrun(i,j,n,iblk)= scale_factor(i,j,iblk)*fswthrun(i,j,n,iblk)
            Sswabsn(i,j,sl1:sl2,iblk) = &
                    scale_factor(i,j,iblk)*Sswabsn(i,j,sl1:sl2,iblk)
            Iswabsn(i,j,il1:il2,iblk) = &
                    scale_factor(i,j,iblk)*Iswabsn(i,j,il1:il2,iblk)
         enddo
      enddo                  ! ncat

      end subroutine prep_radiation_iblk

!=======================================================================
!BOP
!
! !ROUTINE: step_therm2 - step post-coupler thermodynamics
!
! !DESCRIPTION:
!
!-----------------------------------------------------------------------
! Wrapper for driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting. Needed for 
! introducing OpenMP threading more simply.
!
! !REVISION HISTORY:
!
! author: Mariana Vertenstein, NCAR
!
! !INTERFACE:
!
      subroutine step_therm2 (dt)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
         iblk,  &  ! block index
         i, j

!      call t_barrierf('cice_step2_therm_BARRIER',MPI_COMM_ICE)
      call t_startf('cice_step2_therm')
      call ice_timer_start(timer_column)  ! column physics
      call ice_timer_start(timer_thermo)  ! thermodynamics
!      call ice_timer_start(timer_tmp)  ! temporary timer

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call step_therm2_iblk(dt, iblk)
      end do
      !$OMP END PARALLEL DO

      !-------------------------------------------------------------------
      ! Ghost cell updates for state variables.
      !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)
      call bound_state (aicen, trcrn, &
                        vicen, vsnon, &
                        eicen, esnon)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Aggregate the updated state variables. 
      !----------------------------------------------------------------- 
 
         call aggregate (nx_block,          ny_block,             &
                         aicen(:,:,:,iblk), trcrn(:,:,:,:,iblk),  &
                         vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
                         eicen(:,:,:,iblk), esnon(:,:,  :,iblk),  &
                         aice (:,:,  iblk), trcr (:,:,:,  iblk),  &
                         vice (:,:,  iblk), vsno (:,:,    iblk),  &
                         eice (:,:,  iblk), esno (:,:,    iblk),  &
                         aice0(:,:,  iblk), tmask(:,:,    iblk),  &
                         ntrcr,             trcr_depend) 


      !-----------------------------------------------------------------
      ! Compute thermodynamic area and volume tendencies.
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            daidtt(i,j,iblk) = (aice(i,j,iblk) - daidtt(i,j,iblk)) / dt
            dvidtt(i,j,iblk) = (vice(i,j,iblk) - dvidtt(i,j,iblk)) / dt
         enddo
         enddo

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call t_stopf('cice_step2_therm')
!      call ice_timer_stop(timer_tmp)  ! temporary timer
      call ice_timer_stop(timer_thermo)  ! column physics
      call ice_timer_stop(timer_column)  ! column physics

      end subroutine step_therm2

!=======================================================================
!BOP
!
! !ROUTINE: step_therm2_iblk - step post-coupler thermodynamics
!
! !DESCRIPTION:
!
!-----------------------------------------------------------------------
! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! NOTE: Ocean fluxes are initialized here.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:

      subroutine step_therm2_iblk (dt, iblk)
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk ! block index
!
!EOP
!
!lipscomb - delete hicen later?
!      real (kind=dbl_kind), &
!         dimension (nx_block,ny_block,ncat,max_blocks) :: &
!         hicen           ! ice thickness (m)

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, n

      integer (kind=int_kind) :: &
         icells          ! number of ice/ocean cells 

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for ice/ocean cells

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts

      l_stop = .false.

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Let rain drain through to the ocean.
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         fresh     (i,j,iblk) = fresh(i,j,iblk)       &
              + frain(i,j,iblk)*aice(i,j,iblk)
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Given thermodynamic growth rates, transport ice between
      ! thickness categories.
      !-----------------------------------------------------------------

      call ice_timer_start(timer_catconv, iblk) ! category conversions

      if (kitd == 1) then
      !-----------------------------------------------------------------
      ! Compute fractional ice area in each grid cell.
      !-----------------------------------------------------------------
         call aggregate_area (nx_block,          ny_block, &
                              aicen(:,:,:,iblk),           &
                              aice (:,:,  iblk), aice0(:,:,iblk))

      !-----------------------------------------------------------------
      ! Identify grid cells with ice.
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo,jhi
         do i = ilo,ihi
            if (aice(i,j,iblk) > puny) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo
         enddo

         if (icells > 0) then

            call linear_itd (nx_block, ny_block,       &
                             icells, indxi, indxj,     &
                             ntrcr,  trcr_depend,      &
                             aicen_init(:,:,:,iblk),   &
                             vicen_init(:,:,:,iblk),   &
                             aicen     (:,:,:,iblk),   &
                             trcrn     (:,:,:,:,iblk), & 
                             vicen     (:,:,:,iblk),   &
                             vsnon     (:,:,:,iblk),   &
                             eicen     (:,:,:,iblk),   &
                             esnon     (:,:,:,iblk),   &
                             aice      (:,:,  iblk),   &
                             aice0     (:,:,  iblk),   &
                             l_stop,                   &
                             istop,    jstop)

            if (l_stop) then
               write (nu_diag,*) 'istep1, my_task, iblk =', &
                                  istep1, my_task, iblk
               write (nu_diag,*) 'Global block:', this_block%block_id
               if (istop > 0 .and. jstop > 0) &
                    write(nu_diag,*) 'Global i and j:', &
                                      this_block%i_glob(istop), &
                                      this_block%j_glob(jstop) 
               call abort_ice ('ice: Linear ITD error')
            endif

         endif

      endif  ! kitd

      call ice_timer_stop(timer_catconv, iblk)  ! category conversions

      !-----------------------------------------------------------------
      ! Add frazil ice growing in leads.
      !-----------------------------------------------------------------

      ! identify ice-ocean cells
      icells = 0
      do j = 1, ny_block
      do i = 1, nx_block
         if (tmask(i,j,iblk)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo               ! i
      enddo               ! j

      call add_new_ice (nx_block,              ny_block, &
                        ntrcr,                 icells,   &
                        indxi,                 indxj,    &
                        tmask    (:,:,  iblk), dt,       &
                        aicen    (:,:,:,iblk),           &
                        trcrn    (:,:,:,:,iblk),         &
                        vicen    (:,:,:,iblk),           &
                        eicen    (:,:,:,iblk),           &
                        aice0    (:,:,  iblk),           &
                        aice     (:,:,  iblk),           &
                        frzmlt   (:,:,  iblk),           &
                        frazil   (:,:,  iblk),           &
                        frz_onset(:,:,  iblk), yday,     &
                        fresh    (:,:,  iblk),           &
                        fsalt    (:,:,  iblk),           &
                        Tf       (:,:,  iblk), l_stop,   &
                        istop, jstop)

      if (l_stop) then
         write (nu_diag,*) 'istep1, my_task, iblk =', &
                            istep1, my_task, iblk
         write (nu_diag,*) 'Global block:', this_block%block_id
         if (istop > 0 .and. jstop > 0) &
              write(nu_diag,*) 'Global i and j:', &
                                this_block%i_glob(istop), &
                                this_block%j_glob(jstop) 
         call abort_ice ('ice: add_new_ice error')
      endif

      !-----------------------------------------------------------------
      ! Melt ice laterally.
      !-----------------------------------------------------------------
      call lateral_melt (nx_block, ny_block,     &
                         ilo, ihi, jlo, jhi,     &
                         dt,                     &
                         fresh     (:,:,  iblk), &
                         fsalt     (:,:,  iblk), &    
                         fhocn     (:,:,  iblk), &
                         fsoot     (:,:,:,iblk), &
                         rside     (:,:,  iblk), &
                         meltl     (:,:,  iblk), &
                         aicen     (:,:,:,iblk), &
                         vicen     (:,:,:,iblk), &
                         vsnon     (:,:,:,iblk), &
                         eicen     (:,:,:,iblk), &
                         esnon     (:,:,:,iblk), &
                         trcrn     (:,:,:,:,iblk) )

      !-----------------------------------------------------------------
      ! For the special case of a single category, adjust the area and
      ! volume (assuming that half the volume change decreases the
      ! thickness, and the other half decreases the area).  
      !-----------------------------------------------------------------

!NOTE - this does not work - hicen_init is not defined - ECH

!         if (ncat==1) &
!              call reduce_area (nx_block, ny_block,     &
!                                ilo, ihi, jlo, jhi,     &
!                                tmask     (:,:,  iblk), &
!                                aicen     (:,:,:,iblk), &
!                                vicen     (:,:,:,iblk), &
!                                hicen_init(:,:,1,iblk), &
!                                hicen     (:,:,1,iblk)) 

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

      call cleanup_itd (nx_block,             ny_block,             &
                        ilo, ihi,             jlo, jhi,             &
                        dt,                   ntrcr,                &
                        aicen   (:,:,:,iblk), trcrn (:,:,:,:,iblk), &
                        vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
                        eicen   (:,:,:,iblk), esnon (:,:,  :,iblk), &
                        aice0   (:,:,  iblk), aice      (:,:,iblk), &
                        trcr_depend,                                &
                        fresh   (:,:,  iblk), fsalt   (:,:,  iblk), &
                        fhocn   (:,:,  iblk), fsoot   (:,:,:,iblk), &
                        tr_aero,                                    &
                        heat_capacity,        l_stop,               &
                        istop,                jstop)

      if (l_stop) then
         write (nu_diag,*) 'istep1, my_task, iblk =', &
                            istep1, my_task, iblk
         write (nu_diag,*) 'Global block:', this_block%block_id
         if (istop > 0 .and. jstop > 0) &
              write(nu_diag,*) 'Global i and j:', &
                                this_block%i_glob(istop), &
                                this_block%j_glob(jstop) 
         call abort_ice ('ice: ITD cleanup error')
      endif

      end subroutine step_therm2_iblk

!=======================================================================
!BOP
!
! !ROUTINE: step_dynamics - step ice dynamics, transport, and ridging
!
! !DESCRIPTION:
!
! Run one time step of dynamics, horizontal transport, and ridging.
! NOTE: The evp and transport modules include boundary updates, so
!       they cannot be done inside a single block loop.  Ridging
!       and cleanup, on the other hand, are single-column operations. 
!       They are called with argument lists inside block loops
!       to increase modularity.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!
! !INTERFACE:

      subroutine step_dynamics(dt_dyn,dt_thm)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt_dyn  , & ! dynamic time step
         dt_thm      ! thermodynamic time step for diagnostics
!
!EOP
!
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: & 
         iblk        , & ! block index 
         i,j         , & ! horizontal indices
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts

      call init_history_dyn     ! initialize dynamic history variables

      dynCnt = dynCnt + 1

      !-----------------------------------------------------------------
      ! Elastic-viscous-plastic ice dynamics
      !-----------------------------------------------------------------

!      call t_barrierf ('cice_step_evp_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_step_evp')
      if (kdyn == 1) call evp (dt_dyn)
      call t_stopf ('cice_step_evp')

      !-----------------------------------------------------------------
      ! Horizontal ice transport
      !-----------------------------------------------------------------

!      call t_barrierf ('cice_step_horz_transport_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_step_horz_transport')
      if (advection == 'upwind') then
         call transport_upwind (dt_dyn)    ! upwind
      else
         call transport_remap (dt_dyn)     ! incremental remapping
      endif
      call t_stopf ('cice_step_horz_transport')

      !-----------------------------------------------------------------
      ! Ridging
      !-----------------------------------------------------------------

      call ice_timer_start(timer_column)
      call ice_timer_start(timer_ridge)
!      call t_barrierf ('cice_step_ridge_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_step_ridge')

      l_stop = .false.

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,&
      !$OMP	                icells,indxi,indxj,l_stop,istop,jstop)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk), iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Identify ice-ocean cells.
      ! Note:  We can not define icells here using aice>puny because
      !        aice has not yet been updated since the transport (and
      !        it may be out of whack, which the ridging helps fix).-ECH
      !-----------------------------------------------------------------
           
         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
!PERF  There is a performance advantage here to only perform the ridging over a limited set of points
!PERF            if (tmask(i,j,iblk) .and. aice(i,j,iblk) > c0) then
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

         if (icells > 0) then

         call ridge_ice (nx_block,             ny_block,                 &
                         dt_dyn,               dt_thm,                   &
                         ntrcr,                icells,                   &
                         indxi,                indxj,                    &
!!                         Delt    (:,:,  iblk), divu      (:,:,  iblk), &
                         rdg_conv(:,:,  iblk), rdg_shear (:,:,  iblk),   &
                         aicen   (:,:,:,iblk), trcrn     (:,:,:,:,iblk), &
                         vicen   (:,:,:,iblk), vsnon     (:,:,:,iblk),   &
                         eicen   (:,:,:,iblk), esnon     (:,:,:,iblk),   &
                         aice0   (:,:,  iblk),                           &
                         trcr_depend,          l_stop,                   &
                         istop,                jstop,                    &   
                         dardg1dt(:,:,iblk),   dardg2dt  (:,:,iblk),     &
                         dvirdgdt(:,:,iblk),   opening   (:,:,iblk),     &
                         fresh   (:,:,iblk),   fhocn     (:,:,iblk),     &
                         fsoot   (:,:,:,iblk))

         if (l_stop) then
            write (nu_diag,*) 'istep1, my_task, iblk =', &
                               istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0) &
                 write(nu_diag,*) 'Global i and j:', &
                                  this_block%i_glob(istop), &
                                  this_block%j_glob(jstop) 
            call abort_ice ('ice: Ridging error')
         endif

         endif

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_stop(timer_ridge)
      call t_stopf ('cice_step_ridge')

!      call t_barrierf ('cice_step_column_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_step_column')

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,&
      !$OMP	                icells,indxi,indxj,l_stop,istop,jstop)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk), iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

         call cleanup_itd (nx_block,             ny_block,             &
                           ilo, ihi,             jlo, jhi,             &
                           dt_thm,               ntrcr,                &
                           aicen   (:,:,:,iblk), trcrn (:,:,:,:,iblk), &
                           vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
                           eicen   (:,:,:,iblk), esnon (:,:,  :,iblk), &
                           aice0   (:,:,  iblk), aice      (:,:,iblk), &
                           trcr_depend,                                &
                           fresh   (:,:,  iblk), fsalt   (:,:,  iblk), &
                           fhocn   (:,:,  iblk), fsoot   (:,:,:,iblk), &
                           tr_aero,                                    &
                           heat_capacity,        l_stop,               &
                           istop,                jstop)

         if (l_stop) then
            write (nu_diag,*) 'istep1, my_task, iblk =', &
                               istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0) &
                 write(nu_diag,*) 'Global i and j:', &
                                  this_block%i_glob(istop), &
                                  this_block%j_glob(jstop) 
            call abort_ice ('ice: ITD cleanup error')
         endif

      enddo              ! iblk
      !$OMP END PARALLEL DO

      call t_stopf ('cice_step_column')

      !-------------------------------------------------------------------
      ! Ghost cell updates for state variables.
      !-------------------------------------------------------------------

!      call t_barrierf ('cice_step_bound_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_step_bound')
      call ice_timer_start(timer_bound)
      call bound_state (aicen, trcrn, &
                        vicen, vsnon, &
                        eicen, esnon)
      call ice_timer_stop(timer_bound)
      call t_stopf ('cice_step_bound')

!      call t_barrierf ('cice_step_agg_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_step_agg')

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Aggregate the updated state variables. 
      !----------------------------------------------------------------- 
 
         call aggregate (nx_block,          ny_block,             &
                         aicen(:,:,:,iblk), trcrn(:,:,:,:,iblk),  &
                         vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
                         eicen(:,:,:,iblk), esnon(:,:,  :,iblk),  &
                         aice (:,:,  iblk), trcr (:,:,:,  iblk),  &
                         vice (:,:,  iblk), vsno (:,:,    iblk),  &
                         eice (:,:,  iblk), esno (:,:,    iblk),  &
                         aice0(:,:,  iblk), tmask(:,:,    iblk),  &
                         ntrcr,             trcr_depend) 

      !-----------------------------------------------------------------
      ! Compute dynamic area and volume tendencies.
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo,jhi
         do i = ilo,ihi
            dvidtd(i,j,iblk) = (vice(i,j,iblk) - dvidtd(i,j,iblk)) /dt_dyn
            daidtd(i,j,iblk) = (aice(i,j,iblk) - daidtd(i,j,iblk)) /dt_dyn
         enddo
         enddo

      enddo              ! iblk
      !$OMP END PARALLEL DO

      call t_stopf ('cice_step_agg')
      call ice_timer_stop(timer_column)

      end subroutine step_dynamics

!=======================================================================
!BOP
!
! !ROUTINE: step_radiation - step pre-coupler radiation
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: Mariana Vertenstein, NCAR
!
! !INTERFACE:

      subroutine step_radiation(dt)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j,n,iblk    ! block index


      alvdr(:,:,:) = c0
      alvdf(:,:,:) = c0
      alidr(:,:,:) = c0
      alidf(:,:,:) = c0
      Sswabsn(:,:,:,:) = c0
#ifdef AEROFRC
      dalvdr_noaero(:,:,:) = c0
      dalvdf_noaero(:,:,:) = c0
      dalidr_noaero(:,:,:) = c0
      dalidf_noaero(:,:,:) = c0
      dSswabsn_noaero(:,:,:,:) = c0
#endif
#ifdef PONDFRC
      dalvdr_nopond(:,:,:) = c0
      dalvdf_nopond(:,:,:) = c0
      dalidr_nopond(:,:,:) = c0
      dalidf_nopond(:,:,:) = c0
      dSswabsn_nopond(:,:,:,:) = c0
#endif
!      do iblk = 1,nblocks
!         alvdr(:,:,iblk) = c0
!         alvdf(:,:,iblk) = c0
!         alidr(:,:,iblk) = c0
!         alidf(:,:,iblk) = c0
!         Sswabsn(:,:,:,iblk) = c0
!      enddo

      if (calc_Tsfc) then

         !$OMP PARALLEL DO PRIVATE(iblk)
!         call t_barrierf('cice_step_radiationib_BARRIER',MPI_COMM_ICE)
         do iblk = 1, nblocks
            call t_startf('cice_step_radiationib')
            call step_radiation_iblk(dt, iblk)
            call t_stopf('cice_step_radiationib')
         end do
         !$OMP END PARALLEL DO

      else    ! .not. calc_Tsfc

         ! Initialize for safety
         do iblk = 1, nblocks
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            alvdrn(i,j,n,iblk) = c0
            alidrn(i,j,n,iblk) = c0
            alvdfn(i,j,n,iblk) = c0
            alidfn(i,j,n,iblk) = c0
            fswsfcn(i,j,n,iblk) = c0
            fswintn(i,j,n,iblk) = c0
            fswthrun(i,j,n,iblk) = c0
         enddo   ! i
         enddo   ! j
         enddo   ! ncat
            Iswabsn(:,:,:,iblk) = c0
            Sswabsn(:,:,:,iblk) = c0
         enddo   ! iblk

      endif    ! calc_Tsfc

      end subroutine step_radiation

!=======================================================================
!BOP
!
! !ROUTINE: step_radiation_iblk - step pre-coupler radiation
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: David A. Bailey, NCAR
!
! !INTERFACE:

      subroutine step_radiation_iblk (dt, iblk)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n           , & ! thickness category index
         il1, il2    , & ! ice layer indices for eice
         sl1, sl2        ! snow layer indices for esno

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      ! snow variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsn             ! snow horizontal fraction
      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) :: &
         rhosnwn     , & ! snow density (kg/m3)
         rsnwn           ! snow grain radius (micro-meters)

      ! pond variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fpn         , & ! pond fraction
         hpn             ! pond depth (m)

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      l_stop = .false.

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Compute cosine of solar zenith angle.
      ! This is used by the delta-Eddington shortwave module.
      ! Albedos are aggregated in merge_fluxes only for cells w/ coszen > 0.
      ! For basic shortwave, simply set coszen to a constant between 0 and 1.
      !-----------------------------------------------------------------

      if (trim(shortwave) == 'dEdd') then ! delta Eddington

!         call t_barrierf('cice_step_rib_swdedd1_BARRIER',MPI_COMM_ICE)
         call t_startf('cice_step_rib_swdedd1')
         ! identify ice-ocean cells
         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

         call compute_coszen (nx_block,         ny_block,       &
                              icells,                           &
                              indxi,            indxj,          &
                              tlat  (:,:,iblk), tlon(:,:,iblk), &
                              coszen(:,:,iblk), dt)

         call t_stopf('cice_step_rib_swdedd1')
      else                     ! basic (ccsm3) shortwave
         coszen(:,:,iblk) = p5 ! sun above the horizon
      endif

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify cells with nonzero ice area
      !-----------------------------------------------------------------
           
         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n,iblk) > puny) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

      !-----------------------------------------------------------------
      ! Solar radiation: albedo and absorbed shortwave
      !-----------------------------------------------------------------

         il1 = ilyr1(n)
         il2 = ilyrn(n)
         sl1 = slyr1(n)
         sl2 = slyrn(n)

         if (trim(shortwave) == 'dEdd') then   ! delta Eddington

      ! note that rhoswn, rsnw, fp, hp and Sswabs ARE NOT dimensioned with ncat
      ! BPB 19 Dec 2006

            ! set snow properties
!            call t_barrierf('cice_step_rib_swdedd2_BARRIER',MPI_COMM_ICE)
            call t_startf('cice_step_rib_swdedd2')
            call shortwave_dEdd_set_snow(nx_block, ny_block,           &
                              icells,                                  &
                              indxi,               indxj,              &
                              aicen(:,:,n,iblk),   vsnon(:,:,n,iblk),  &
                              trcrn(:,:,nt_Tsfc,n,iblk), fsn,          &
                              rhosnwn,             rsnwn)
            call t_stopf('cice_step_rib_swdedd2')


            if (.not. tr_pond) then

               ! set pond properties
!               call t_barrierf('cice_step_rib_swdedd3_BARRIER',MPI_COMM_ICE)
               call t_startf('cice_step_rib_swdedd3')
               call shortwave_dEdd_set_pond(nx_block, ny_block,            &
                                 icells,                                   &
                                 indxi,               indxj,               &
                                 aicen(:,:,n,iblk),                        &
                                 trcrn(:,:,nt_Tsfc,n,iblk),                &
                                 fsn,                 fpn,                 &
                                 hpn)
               call t_stopf('cice_step_rib_swdedd3')
            else

               fpn(:,:) = apondn(:,:,n,iblk)
               hpn(:,:) = hpondn(:,:,n,iblk)

            endif

#ifdef AEROFRC
            ! Diagnose aerosol forcing

            if (tr_aero) then

            tr_aero = .false.

!            call t_barrierf('cice_step_rib_swdedd4_BARRIER',MPI_COMM_ICE)
            call t_startf('cice_step_rib_swdedd4')
            call shortwave_dEdd(nx_block,        ny_block,            &
                                icells,                                 &
                                indxi,             indxj,               &
                                coszen(:,:, iblk),                      &
                                aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                                vsnon(:,:,n,iblk), fsn,                 &
                                rhosnwn,           rsnwn,               &
                                fpn,               hpn,                 &
                                trcrn(:,:,:,n,iblk), tarea(:,:,iblk),   &
                                swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                                swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                                dalvdrn_noaero(:,:,n,iblk),              &
                                dalvdfn_noaero(:,:,n,iblk),              &
                                dalidrn_noaero(:,:,n,iblk),              &
                                dalidfn_noaero(:,:,n,iblk),              &
                                dfswsfcn_noaero(:,:,n,iblk),             &
                                dfswintn_noaero(:,:,n,iblk),             &
                                dfswthrun_noaero(:,:,n,iblk),            &
                                dSswabsn_noaero(:,:,sl1:sl2,iblk),       &
                                dIswabsn_noaero(:,:,il1:il2,iblk),       &
                                dalbicen_noaero(:,:,n,iblk),             &
                                dalbsnon_noaero(:,:,n,iblk),             &
                                dalbpndn_noaero(:,:,n,iblk))
            call t_stopf('cice_step_rib_swdedd4')
            tr_aero = .true.

            endif
#endif

#ifdef PONDFRC
            ! Diagnose pond forcing

            if (tr_pond) then

            tr_aero = .false.

            fpn(:,:) = c0
            hpn(:,:) = c0

!            call t_barrierf('cice_step_rib_swdedd5_BARRIER',MPI_COMM_ICE)
            call t_startf('cice_step_rib_swdedd5')
            call shortwave_dEdd(nx_block,        ny_block,            &
                                icells,                                 &
                                indxi,             indxj,               &
                                coszen(:,:, iblk),                      &
                                aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                                vsnon(:,:,n,iblk), fsn,                 &
                                rhosnwn,           rsnwn,               &
                                fpn,               hpn,                 &
                                trcrn(:,:,:,n,iblk), tarea(:,:,iblk),   &
                                swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                                swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                                dalvdrn_nopond(:,:,n,iblk),              &
                                dalvdfn_nopond(:,:,n,iblk),              &
                                dalidrn_nopond(:,:,n,iblk),              &
                                dalidfn_nopond(:,:,n,iblk),              &
                                dfswsfcn_nopond(:,:,n,iblk),             &
                                dfswintn_nopond(:,:,n,iblk),             &
                                dfswthrun_nopond(:,:,n,iblk),            &
                                dSswabsn_nopond(:,:,sl1:sl2,iblk),       &
                                dIswabsn_nopond(:,:,il1:il2,iblk),       &
                                dalbicen_nopond(:,:,n,iblk),             &
                                dalbsnon_nopond(:,:,n,iblk),             &
                                dalbpndn_nopond(:,:,n,iblk))
            call t_stopf('cice_step_rib_swdedd5')

            fpn(:,:) = apondn(:,:,n,iblk)
            hpn(:,:) = hpondn(:,:,n,iblk)

            tr_aero = .true.

            endif
#endif
#ifdef CCSM3FRC
!            call t_barrierf('cice_step_rib_swdedd6_BARRIER',MPI_COMM_ICE)
            call t_startf('cice_step_rib_swdedd6')
            call shortwave_ccsm3(nx_block,     ny_block,           &
                           icells,                                 &
                           indxi,             indxj,               &
                           aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                           vsnon(:,:,n,iblk),                      &
                           trcrn(:,:,nt_Tsfc,n,iblk),              &
                           swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                           swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                           dalvdrn_ccsm3(:,:,n,iblk),               &
                           dalidrn_ccsm3(:,:,n,iblk),               &
                           dalvdfn_ccsm3(:,:,n,iblk),               &
                           dalidfn_ccsm3(:,:,n,iblk),               &
                           dfswsfcn_ccsm3(:,:,n,iblk),              &
                           dfswintn_ccsm3(:,:,n,iblk),              &
                           dfswthrun_ccsm3(:,:,n,iblk),             &
                           dIswabsn_ccsm3(:,:,il1:il2,iblk),        &
                           dalbicen_ccsm3(:,:,n,iblk),              &
                           dalbsnon_ccsm3(:,:,n,iblk))
            call t_stopf('cice_step_rib_swdedd6')

#endif
!            call t_barrierf('cice_step_rib_swdedd7_BARRIER',MPI_COMM_ICE)
            call t_startf('cice_step_rib_swdedd7')
            call shortwave_dEdd(nx_block,        ny_block,            &
                                icells,                                 &
                                indxi,             indxj,               &
                                coszen(:,:, iblk),                      &
                                aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                                vsnon(:,:,n,iblk), fsn,                 &
                                rhosnwn,           rsnwn,               &
                                fpn,               hpn,                 &
                                trcrn(:,:,:,n,iblk), tarea(:,:,iblk),   &
                                swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                                swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                                alvdrn(:,:,n,iblk),alvdfn(:,:,n,iblk),  &
                                alidrn(:,:,n,iblk),alidfn(:,:,n,iblk),  &
                                fswsfcn(:,:,n,iblk),fswintn(:,:,n,iblk),&
                                fswthrun(:,:,n,iblk),                   &
                                Sswabsn(:,:,sl1:sl2,iblk),              &
                                Iswabsn(:,:,il1:il2,iblk),              &
                                albicen(:,:,n,iblk),                    &
                                albsnon(:,:,n,iblk),albpndn(:,:,n,iblk))
            call t_stopf('cice_step_rib_swdedd7')

         else

            Sswabsn(:,:,sl1:sl2,iblk) = c0

!            call t_barrierf('cice_step_rib_ccsm3_BARRIER',MPI_COMM_ICE)
            call t_startf('cice_step_rib_ccsm3')
            call shortwave_ccsm3(nx_block,     ny_block,           &
                           icells,                                 &
                           indxi,             indxj,               &
                           aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                           vsnon(:,:,n,iblk),                      &
                           trcrn(:,:,nt_Tsfc,n,iblk),              &
                           swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                           swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                           alvdrn(:,:,n,iblk),alidrn(:,:,n,iblk),  &
                           alvdfn(:,:,n,iblk),alidfn(:,:,n,iblk),  &
                           fswsfcn(:,:,n,iblk),fswintn(:,:,n,iblk),&
                           fswthrun(:,:,n,iblk),                   &
                           Iswabsn(:,:,il1:il2,iblk),              &
                           albicen(:,:,n,iblk),albsnon(:,:,n,iblk))

            call t_stopf('cice_step_rib_ccsm3')
         endif

#ifdef AEROFRC
            dalvdrn_noaero(:,:,n,iblk) = dalvdrn_noaero(:,:,n,iblk)-alvdrn(:,:,n,iblk)
            dalvdfn_noaero(:,:,n,iblk) = dalvdfn_noaero(:,:,n,iblk)-alvdfn(:,:,n,iblk)
            dalidrn_noaero(:,:,n,iblk) = dalidrn_noaero(:,:,n,iblk)-alidrn(:,:,n,iblk)
            dalidfn_noaero(:,:,n,iblk) = dalidfn_noaero(:,:,n,iblk)-alidfn(:,:,n,iblk)
            dfswsfcn_noaero(:,:,n,iblk) = dfswsfcn_noaero(:,:,n,iblk)-fswsfcn(:,:,n,iblk)
            dfswintn_noaero(:,:,n,iblk) = dfswintn_noaero(:,:,n,iblk)-fswintn(:,:,n,iblk)
            dfswthrun_noaero(:,:,n,iblk) = dfswthrun_noaero(:,:,n,iblk)-fswthrun(:,:,n,iblk)
            dfswabsn_noaero(:,:,n,iblk) = dfswsfcn_noaero(:,:,n,iblk)+dfswintn_noaero(:,:,n,iblk)+dfswthrun_noaero(:,:,n,iblk)
            dalbicen_noaero(:,:,n,iblk) = dalbicen_noaero(:,:,n,iblk)-albicen(:,:,n,iblk)
            dalbsnon_noaero(:,:,n,iblk) = dalbsnon_noaero(:,:,n,iblk)-albsnon(:,:,n,iblk)
            dalbpndn_noaero(:,:,n,iblk) = dalbpndn_noaero(:,:,n,iblk)-albpndn(:,:,n,iblk)
            dSswabsn_noaero(:,:,sl1:sl2,iblk) = dSswabsn_noaero(:,:,sl1:sl2,iblk)-Sswabsn(:,:,sl1:sl2,iblk)
            dIswabsn_noaero(:,:,il1:il2,iblk) = dIswabsn_noaero(:,:,il1:il2,iblk)-Iswabsn(:,:,il1:il2,iblk)
#endif
#ifdef CCSM3FRC
            dalvdrn_ccsm3(:,:,n,iblk) = dalvdrn_ccsm3(:,:,n,iblk)-alvdrn(:,:,n,iblk)
            dalvdfn_ccsm3(:,:,n,iblk) = dalvdfn_ccsm3(:,:,n,iblk)-alvdfn(:,:,n,iblk)
            dalidrn_ccsm3(:,:,n,iblk) = dalidrn_ccsm3(:,:,n,iblk)-alidrn(:,:,n,iblk)
            dalidfn_ccsm3(:,:,n,iblk) = dalidfn_ccsm3(:,:,n,iblk)-alidfn(:,:,n,iblk)
            dfswsfcn_ccsm3(:,:,n,iblk) = dfswsfcn_ccsm3(:,:,n,iblk)-fswsfcn(:,:,n,iblk)
            dfswintn_ccsm3(:,:,n,iblk) = dfswintn_ccsm3(:,:,n,iblk)-fswintn(:,:,n,iblk)
            dfswthrun_ccsm3(:,:,n,iblk) = dfswthrun_ccsm3(:,:,n,iblk)-fswthrun(:,:,n,iblk)
            dfswabsn_ccsm3(:,:,n,iblk) = dfswsfcn_ccsm3(:,:,n,iblk)+dfswintn_ccsm3(:,:,n,iblk)+dfswthrun_ccsm3(:,:,n,iblk)
            dalbicen_ccsm3(:,:,n,iblk) = dalbicen_ccsm3(:,:,n,iblk)-albicen(:,:,n,iblk)
            dalbsnon_ccsm3(:,:,n,iblk) = dalbsnon_ccsm3(:,:,n,iblk)-albsnon(:,:,n,iblk)
            dIswabsn_ccsm3(:,:,il1:il2,iblk) = dIswabsn_ccsm3(:,:,il1:il2,iblk)-Iswabsn(:,:,il1:il2,iblk)
#endif
#ifdef PONDFRC
            dalvdrn_nopond(:,:,n,iblk) = dalvdrn_nopond(:,:,n,iblk)-alvdrn(:,:,n,iblk)
            dalvdfn_nopond(:,:,n,iblk) = dalvdfn_nopond(:,:,n,iblk)-alvdfn(:,:,n,iblk)
            dalidrn_nopond(:,:,n,iblk) = dalidrn_nopond(:,:,n,iblk)-alidrn(:,:,n,iblk)
            dalidfn_nopond(:,:,n,iblk) = dalidfn_nopond(:,:,n,iblk)-alidfn(:,:,n,iblk)
            dfswsfcn_nopond(:,:,n,iblk) = dfswsfcn_nopond(:,:,n,iblk)-fswsfcn(:,:,n,iblk)
            dfswintn_nopond(:,:,n,iblk) = dfswintn_nopond(:,:,n,iblk)-fswintn(:,:,n,iblk)
            dfswthrun_nopond(:,:,n,iblk) = dfswthrun_nopond(:,:,n,iblk)-fswthrun(:,:,n,iblk)
            dfswabsn_nopond(:,:,n,iblk) = dfswsfcn_nopond(:,:,n,iblk)+dfswintn_nopond(:,:,n,iblk)+dfswthrun_nopond(:,:,n,iblk)
            dalbicen_nopond(:,:,n,iblk) = dalbicen_nopond(:,:,n,iblk)-albicen(:,:,n,iblk)
            dalbsnon_nopond(:,:,n,iblk) = dalbsnon_nopond(:,:,n,iblk)-albsnon(:,:,n,iblk)
            dalbpndn_nopond(:,:,n,iblk) = dalbpndn_nopond(:,:,n,iblk)-albpndn(:,:,n,iblk)
            dSswabsn_nopond(:,:,sl1:sl2,iblk) = dSswabsn_nopond(:,:,sl1:sl2,iblk)-Sswabsn(:,:,sl1:sl2,iblk)
            dIswabsn_nopond(:,:,il1:il2,iblk) = dIswabsn_nopond(:,:,il1:il2,iblk)-Iswabsn(:,:,il1:il2,iblk)
#endif
      enddo                  ! ncat

      end subroutine step_radiation_iblk

!=======================================================================

      end module ice_step_mod

!=======================================================================
