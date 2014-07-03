!=======================================================================
!
!BOP
!
! !MODULE: ice_restoring
!
! !DESCRIPTION:
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! !REVISION HISTORY:
!  SVN:$Id: $
!
! authors: Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      module ice_restoring
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks
      use ice_calendar, only: dt
      use ice_domain
      use ice_domain_size
      use ice_communicate, only: my_task, master_task
      use ice_constants
      use ice_exit
      use ice_fileunits
      use ice_forcing, only: trestore, trest
      use ice_state
      use ice_timers
!
!EOP
!
      implicit none
      save

      logical (kind=log_kind) :: &
         restore_ice                 ! restore ice state if true

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable :: &
         aicen_rest , & ! concentration of ice
         vicen_rest , & ! volume per unit area of ice          (m)
         vsnon_rest , & ! volume per unit area of snow         (m)
         eicen_rest , & ! energy of melting for each ice layer  (J/m^2)
         esnon_rest     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable :: &
         trcrn_rest     ! tracers

!=======================================================================

      contains

!=======================================================================
!BOP
! !IROUTINE: ice_HaloRestore_init
! !INTERFACE:

 subroutine ice_HaloRestore_init

! !DESCRIPTION:
!  Allocates and initializes arrays needed for restoring the ice state 
!  in cells surrounding the grid.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   allocate (aicen_rest(nx_block,ny_block,ncat,max_blocks), &
             vicen_rest(nx_block,ny_block,ncat,max_blocks), &
             vsnon_rest(nx_block,ny_block,ncat,max_blocks), &
             eicen_rest(nx_block,ny_block,ntilyr,max_blocks), &
             esnon_rest(nx_block,ny_block,ntslyr,max_blocks), &
             trcrn_rest(nx_block,ny_block,ntrcr,ncat,max_blocks))

   ! initialize to the default initial ice state
   ! these could be set to values read from a file...

      call ice_timer_start(timer_bound)
      call bound_state (aicen, trcrn, &
                        vicen, vsnon, &
                        eicen, esnon)
      call ice_timer_stop(timer_bound)

   aicen_rest(:,:,:,:) = aicen(:,:,:,:)
   vicen_rest(:,:,:,:) = vicen(:,:,:,:)
   vsnon_rest(:,:,:,:) = vsnon(:,:,:,:)
   eicen_rest(:,:,:,:) = eicen(:,:,:,:)
   esnon_rest(:,:,:,:) = esnon(:,:,:,:)
   trcrn_rest(:,:,:,:,:) = trcrn(:,:,1:ntrcr,:,:)

   if (my_task == master_task) &
      write (nu_diag,*) 'ice restoring timescale = ',trestore,' days' 

 end subroutine ice_HaloRestore_init

!=======================================================================

!BOP
! !IROUTINE: ice_HaloRestore
! !INTERFACE:

 subroutine ice_HaloRestore

! !DESCRIPTION:
!  This subroutine is intended for restoring the ice state to desired
!  values in cells surrounding the grid.
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
! !REVISION HISTORY:
!  same as module
!
! !USES:

   use ice_distribution

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,nt,n,      &! dummy loop indices
     ilo,ihi,jlo,jhi,    &! beginning and end of physical domain
     ibc,                &! ghost cell column or row
     npad                 ! padding column/row counter

   type (block) :: &
     this_block  ! block info for current block

   real (dbl_kind) :: &
     ctime                ! dt/trest

   call ice_timer_start(timer_bound)

!-----------------------------------------------------------------------
!
!  Initialize
!
!-----------------------------------------------------------------------

      ! for now, use same restoring constant as for SST
      if (trestore == 0) then
         trest = dt          ! use data instantaneously
      else
         trest = real(trestore,kind=dbl_kind) * secday ! seconds
      endif
      ctime = dt/trest

!-----------------------------------------------------------------------
!
!  Restore values in cells surrounding the grid
!
!-----------------------------------------------------------------------

   do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, ny_block
            do i = 1, ilo
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
            do n = 1, ntilyr
            do j = 1, ny_block
            do i = 1, ilo
               eicen(i,j,n,iblk) = eicen(i,j,n,iblk) &
                  + (eicen_rest(i,j,n,iblk)-eicen(i,j,n,iblk))*ctime
            enddo
            enddo
            enddo
            do n = 1, ntslyr
            do j = 1, ny_block
            do i = 1, ilo
               esnon(i,j,n,iblk) = esnon(i,j,n,iblk) &
                  + (esnon_rest(i,j,n,iblk)-esnon(i,j,n,iblk))*ctime
            enddo
            enddo
            enddo
         endif

      elseif (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block + 1
            npad = 0
            do i = nx_block, 1, - 1
               if (this_block%i_glob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + this_block%j_glob(j)
                  enddo
               endif
               if (npad == 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = 1, ny_block
            do i = ihi, ibc
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
            do n = 1, ntilyr
            do j = 1, ny_block
            do i = ihi, ibc
               eicen(i,j,n,iblk) = eicen(i,j,n,iblk) &
                  + (eicen_rest(i,j,n,iblk)-eicen(i,j,n,iblk))*ctime
            enddo
            enddo
            enddo
            do n = 1, ntslyr
            do j = 1, ny_block
            do i = ihi, ibc
               esnon(i,j,n,iblk) = esnon(i,j,n,iblk) &
                  + (esnon_rest(i,j,n,iblk)-esnon(i,j,n,iblk))*ctime
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, jlo
            do i = 1, nx_block
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
            do n = 1, ntilyr
            do j = 1, jlo
            do i = 1, nx_block
               eicen(i,j,n,iblk) = eicen(i,j,n,iblk) &
                  + (eicen_rest(i,j,n,iblk)-eicen(i,j,n,iblk))*ctime
            enddo
            enddo
            enddo
            do n = 1, ntslyr
            do j = 1, jlo
            do i = 1, nx_block
               esnon(i,j,n,iblk) = esnon(i,j,n,iblk) &
                  + (esnon_rest(i,j,n,iblk)-esnon(i,j,n,iblk))*ctime
            enddo
            enddo
            enddo
         endif

      elseif (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_boundary_type) /= 'cyclic' .and. &
             trim(ns_boundary_type) /= 'tripole' .and. &
             trim(ns_boundary_type) /= 'tripoleT') then
            ! locate ghost cell row (avoid padding)
            ibc = ny_block + 1
            npad = 0
            do j = ny_block, 1, - 1
               if (this_block%j_glob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + this_block%i_glob(i)
                  enddo
               endif
               if (npad == 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = jhi, ibc
            do i = 1, nx_block
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
            do n = 1, ntilyr
            do j = jhi, ibc
            do i = 1, nx_block
               eicen(i,j,n,iblk) = eicen(i,j,n,iblk) &
                  + (eicen_rest(i,j,n,iblk)-eicen(i,j,n,iblk))*ctime
            enddo
            enddo
            enddo
            do n = 1, ntslyr
            do j = jhi, ibc
            do i = 1, nx_block
               esnon(i,j,n,iblk) = esnon(i,j,n,iblk) &
                  + (esnon_rest(i,j,n,iblk)-esnon(i,j,n,iblk))*ctime
            enddo
            enddo
            enddo
         endif
      endif

   enddo ! iblk

   call ice_timer_stop(timer_bound)

 end subroutine ice_HaloRestore

!=======================================================================

      end module ice_restoring

!=======================================================================
