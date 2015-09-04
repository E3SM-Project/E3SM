#define _SINGLEMSG 1
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ice_gather_scatter

 module ice_gather_scatter

! !DESCRIPTION:
!  This module contains routines for gathering data to a single
!  processor from a distributed array, and scattering data from a
!  single processor to a distributed array.
!
!  NOTE: The arrays gathered and scattered are assumed to have
!        horizontal dimensions (nx_block, ny_block).
!
! !REVISION HISTORY:
!  SVN:$Id: ice_gather_scatter.F90 131 2008-05-30 16:53:40Z eclare $
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
! Jan. 2008: Elizabeth Hunke replaced old routines with new POP
!              infrastructure, added specialized routine scatter_global_stress

! !USES:

   use ice_kinds_mod
   use ice_communicate
   use ice_constants
   use ice_blocks
   use ice_distribution
   use ice_domain
   use ice_domain_size
   use ice_exit

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: gather_global,      &
             scatter_global,     &
             gatherArray,        & 
             scatter_global_stress

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  overload module functions
!
!-----------------------------------------------------------------------

   interface gather_global
     module procedure gather_global_dbl,  &
                      gather_global_real, &
                      gather_global_int
   end interface 

   interface scatter_global
     module procedure scatter_global_dbl,  &
                      scatter_global_real, &
                      scatter_global_int
   end interface 

   interface gatherArray
     module procedure gatherArray_dbl
   end interface

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains


 subroutine gatherArray_dbl(array_g,array,length,root)

   include 'mpif.h'

   real(dbl_kind) :: array_g(:)  ! The concatonated array
   real(dbl_kind) :: array(:)    ! the local piece of the array
   integer(int_kind) :: length   ! number of elements in the array
   integer(int_kind) :: root     ! root to which to collect the array

   integer(int_kind) :: ierr

#ifdef _USE_FLOW_CONTROL
   integer (int_kind) :: &
     i                  ,&! loop index
     nprocs             ,&! number of processes in communicator
     itask              ,&! task loop index
     mtag               ,&! MPI message tag
     displs             ,&! MPI message length
     rcv_request        ,&! request id
     signal               ! MPI handshaking variable

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: status

   signal = 1
   nprocs = get_num_procs()

   if (root .eq. my_task) then
      do itask=0, nprocs-1
         if (itask .ne. root) then
            mtag = 3*mpitag_gs+itask
            displs = itask*length
            call mpi_irecv( array_g(displs+1), length, &
                            MPI_REAL8, itask, mtag, MPI_COMM_ICE, &
                            rcv_request, ierr )
            call mpi_send ( signal, 1, mpi_integer, itask, mtag, &
                            MPI_COMM_ICE, ierr )
            call mpi_wait ( rcv_request, status, ierr )
         end if
      end do

! copy local data
      displs = root*length
      do i=1,length
         array_g(displs+i) = array(i)
      enddo

   else
      mtag = 3*mpitag_gs+my_task

      call mpi_recv ( signal, 1, mpi_integer, root, mtag, MPI_COMM_ICE, &
                      status, ierr )
      call mpi_rsend ( array, length, MPI_REAL8, root, mtag, MPI_COMM_ICE, &
                      ierr )
   endif
#else
   call MPI_Gather(array,length,MPI_REAL8,array_g, length,MPI_REAL8,root, &
                   MPI_COMM_ICE,ierr)
#endif

 end subroutine gatherArray_dbl


!***********************************************************************
!BOP
! !IROUTINE: gather_global
! !INTERFACE:

 subroutine gather_global_dbl(ARRAY_G, ARRAY, dst_task, src_dist)

! !DESCRIPTION:
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific inteface for double precision arrays 
!  corresponding to the generic interface gather_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific inteface based
!  on the data type of the input argument).


! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist   ! distribution of blocks in the source array

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

! !OUTPUT PARAMETERS:

   real (dbl_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     src_task       ,&! source of message
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

#ifdef _SINGLEMSG
   real (dbl_kind), dimension(:,:,:), allocatable :: &
     msg_buffer
#else
   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer
#endif

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind) :: ib, ig, itask, nprocs, maxBlocks
   integer (int_kind) :: iig,ijg,it
   integer (int_kind) :: msgLen, msgTag

#ifdef _USE_FLOW_CONTROL
   integer (int_kind) :: &
     rcv_request    ,&! request id
     signal           ! MPI handshaking variable

   signal = 1
#endif

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   nprocs = get_num_procs()

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%blockLocation(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%blockLocation(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = c0
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest
#ifdef _SINGLEMSG
     allocate (msg_buffer(nx_block,ny_block,max_blocks))
     do itask = 0,nprocs-1
        if(itask /= dst_task) then 
           maxBlocks = src_dist%BlockCnt(itask+1)
           msgLen = nx_block*ny_block*maxBlocks
           msgTag = 3*mpitag_gs+itask
           if(maxBlocks>0) then 

#ifdef _USE_FLOW_CONTROL
             call MPI_IRECV(msg_buffer(:,:,1:maxBlocks),msgLen, &
               mpiR8, itask, msgTag, MPI_COMM_ICE, rcv_request, &
               ierr)

             call MPI_SEND(signal, 1, mpi_integer, itask, &
                           msgTag, MPI_COMM_ICE, ierr)

             call MPI_WAIT(rcv_request, status, ierr)
#else
             call MPI_RECV(msg_buffer(:,:,1:maxBlocks),msgLen, &
               mpiR8, itask, msgTag, MPI_COMM_ICE, status, ierr)
#endif

             do ib=1,maxBlocks
               ig = src_dist%blockIndex(itask+1,ib)
      	       this_block = get_block(ig,ig)
               do j=this_block%jlo,this_block%jhi
                  do i=this_block%ilo,this_block%ihi
                     iig = this_block%i_glob(i) 
                     ijg = this_block%j_glob(j) 
                     ARRAY_G(iig,ijg) = msg_buffer(i,j,ib)
                  end do
               end do
              enddo
            endif
        endif
     enddo

     deallocate(msg_buffer)
#else
     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) > 0 .and. &
           src_dist%blockLocation(n) /= my_task+1) then

         this_block = get_block(n,n)

#ifdef _USE_FLOW_CONTROL
         call MPI_IRECV(msg_buffer, size(msg_buffer), &
                        mpiR8, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                        MPI_COMM_ICE, rcv_request, ierr)

         call MPI_SEND(signal, 1, mpi_integer, &
                       src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

         call MPI_WAIT(rcv_request, status, ierr)
#else
         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpiR8, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)
#endif

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)
#endif

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

#ifdef _SINGLEMSG

     if(nblocks>0) then 
        msgLen = nx_block*ny_block*nblocks
        msgTag = 3*mpitag_gs+my_task
#ifdef _USE_FLOW_CONTROL
        call MPI_RECV(signal, 1, mpi_integer, dst_task, msgTag, &
                      MPI_COMM_ICE, status, ierr)
	call MPI_RSEND(ARRAY(:,:,1:nblocks),msgLen,mpiR8,dst_task, &
           msgTag,MPI_COMM_ICE,ierr)
#else
	call MPI_SEND(ARRAY(:,:,1:nblocks),msgLen,mpiR8,dst_task, &
           msgTag,MPI_COMM_ICE,ierr)
#endif
     endif

#else
     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then

         nsends = nsends + 1
         src_block = src_dist%blockLocalID(n)
#ifdef _USE_FLOW_CONTROL
         call MPI_RECV(signal, 1, mpi_integer, dst_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)
         call MPI_IRSEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpiR8, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
#else
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpiR8, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
#endif
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)
#endif

   endif

!-----------------------------------------------------------------------

 end subroutine gather_global_dbl

!***********************************************************************

 subroutine gather_global_real(ARRAY_G, ARRAY, dst_task, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist       ! distribution of blocks in the source array

   real (real_kind), dimension(:,:,:), intent(in) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (real_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

#ifdef _SINGLEMSG
   real (real_kind), dimension(:,:,:), allocatable :: &
     msg_buffer
#else
   real (real_kind), dimension(:,:), allocatable :: &
     msg_buffer
#endif

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind) :: ib, ig, itask, nprocs, maxBlocks
   integer (int_kind) :: iig,ijg
   integer (int_kind) :: msgLen, msgTag

#ifdef _USE_FLOW_CONTROL
   integer (int_kind) :: &
     rcv_request    ,&! request id
     signal           ! MPI handshaking variable

   signal = 1
#endif

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   nprocs = get_num_procs()

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%blockLocation(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

       !*** fill land blocks with special values

       else if (src_dist%blockLocation(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = spval
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

#ifdef _SINGLEMSG

     allocate (msg_buffer(nx_block,ny_block,max_blocks))
     do itask = 0,nprocs-1
        if( itask /= dst_task) then
           maxBlocks = src_dist%BlockCnt(itask+1)
           msgLen = nx_block*ny_block*maxBlocks
           msgTag = 3*mpitag_gs+itask
           if( maxBlocks>0) then 

#ifdef _USE_FLOW_CONTROL
              call MPI_IRECV(msg_buffer(:,:,1:maxBlocks), msgLen, &
                   mpiR4, itask, msgTag ,MPI_COMM_ICE, &
                   rcv_request, ierr)

              call MPI_SEND(signal, 1, mpi_integer, itask, &
                            msgTag, MPI_COMM_ICE, ierr)

              call MPI_WAIT(rcv_request, status, ierr)
#else
              call MPI_RECV(msg_buffer(:,:,1:maxBlocks), msgLen, &
                   mpiR4, itask, msgTag ,MPI_COMM_ICE, status, ierr)
#endif

              do ib=1,maxBlocks
                 ig = src_dist%blockIndex(itask+1,ib)
                 this_block = get_block(ig,ig)
                 do j=this_block%jlo,this_block%jhi
                   do i=this_block%ilo,this_block%ihi
                     iig = this_block%i_glob(i) 
                     ijg = this_block%j_glob(j) 
                     ARRAY_G(iig,ijg) = msg_buffer(i,j,ib)
                   end do
                 end do
              enddo
           endif
        endif
     enddo

     deallocate(msg_buffer)

#else
     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) > 0 .and. &
           src_dist%blockLocation(n) /= my_task+1) then

         this_block = get_block(n,n)

#ifdef _USE_FLOW_CONTROL
         call MPI_IRECV(msg_buffer, size(msg_buffer), &
                        mpiR4, src_dist%blockLocation(n)-1, &
                        3*mpitag_gs+n, MPI_COMM_ICE, rcv_request, ierr)

         call MPI_SEND(signal, 1, mpi_integer, &
                       src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

         call MPI_WAIT(rcv_request, status, ierr)
#else
         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpiR4, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)
#endif

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do
     deallocate(msg_buffer)
#endif


!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else
#ifdef _SINGLEMSG

     if(nblocks>0) then
        msgLen = nx_block*ny_block*nblocks
        msgTag = 3*mpitag_gs+my_task 
#ifdef _USE_FLOW_CONTROL
        call MPI_RECV(signal, 1, mpi_integer, dst_task, msgTag, &
                      MPI_COMM_ICE, status, ierr)
	call MPI_RSEND(ARRAY(:,:,1:nblocks),msgLen,mpiR4,dst_task, &
           msgTag,MPI_COMM_ICE,ierr)
#else
	call MPI_SEND(ARRAY(:,:,1:nblocks),msgLen,mpiR4,dst_task, &
           msgTag,MPI_COMM_ICE,ierr)
#endif
     endif

#else

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then

         nsends = nsends + 1
         src_block = src_dist%blockLocalID(n)
#ifdef _USE_FLOW_CONTROL
         call MPI_RECV(signal, 1, mpi_integer, dst_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)
         call MPI_IRSEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpiR4, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
#else
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpiR4, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
#endif
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)
#endif

   endif

!-----------------------------------------------------------------------

 end subroutine gather_global_real

!***********************************************************************

 subroutine gather_global_int(ARRAY_G, ARRAY, dst_task, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist       ! distribution of blocks in the source array

   integer (int_kind), dimension(:,:,:), intent(in) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

#ifdef _SINGLEMSG
   integer (int_kind), dimension(:,:,:), allocatable :: &
     msg_buffer
#else
   integer (int_kind), dimension(:,:), allocatable :: &
     msg_buffer
#endif

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind) :: ib, ig, len, itask, nprocs, maxBlocks
   integer (int_kind) :: iig,ijg
   integer (int_kind) :: msgLen, msgTag

#ifdef _USE_FLOW_CONTROL
   integer (int_kind) :: &
     rcv_request    ,&! request id
     signal           ! MPI handshaking variable

   signal = 1
#endif

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   nprocs = get_num_procs()

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%blockLocation(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%blockLocation(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = 0
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

#ifdef _SINGLEMSG

     allocate (msg_buffer(nx_block,ny_block,max_blocks))
     do itask = 0,nprocs-1
        if( itask /= dst_task) then
           maxBlocks = src_dist%BlockCnt(itask+1)
           msgLen = nx_block*ny_block*max_blocks
           msgTag = 3*mpitag_gs+itask
           if(maxBLocks>0) then 
#ifdef _USE_FLOW_CONTROL
             call MPI_IRECV(msg_buffer(:,:,1:maxBlocks),msgLen, &
                  mpi_integer, itask, msgTag, MPI_COMM_ICE, &
                  rcv_request, ierr)

             call MPI_SEND(signal, 1, mpi_integer, itask, &
                           msgTag, MPI_COMM_ICE, ierr)

             call MPI_WAIT(rcv_request, status, ierr)
#else
             call MPI_RECV(msg_buffer(:,:,1:maxBlocks),msgLen, &
                   mpi_integer, itask, msgTag, MPI_COMM_ICE, status, ierr)
#endif
             do ib=1,maxBlocks
                ig = src_dist%blockIndex(itask+1,ib)
                this_block = get_block(ig,ig)
                do j=this_block%jlo,this_block%jhi
                   do i=this_block%ilo,this_block%ihi
                      iig = this_block%i_glob(i) 
                      ijg = this_block%j_glob(j) 
                      ARRAY_G(iig,ijg) = msg_buffer(i,j,ib)
                   end do
                end do
              enddo
           endif
        endif
     enddo
     deallocate(msg_buffer)

#else
     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) > 0 .and. &
           src_dist%blockLocation(n) /= my_task+1) then

         this_block = get_block(n,n)

#ifdef _USE_FLOW_CONTROL
         call MPI_IRECV(msg_buffer, size(msg_buffer), &
                        mpi_integer, src_dist%blockLocation(n)-1, &
                        3*mpitag_gs+n, MPI_COMM_ICE, rcv_request, ierr)

         call MPI_SEND(signal, 1, mpi_integer, &
                       src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

         call MPI_WAIT(rcv_request, status, ierr)
#else
         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpi_integer, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)
#endif

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do
     deallocate(msg_buffer)
#endif


!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

#ifdef _SINGLEMSG

     if(nblocks>0) then
        msgLen = nx_block*ny_block*nblocks
        msgTag = 3*mpitag_gs+my_task
#ifdef _USE_FLOW_CONTROL
        call MPI_RECV(signal, 1, mpi_integer, dst_task, msgTag, &
                      MPI_COMM_ICE, status, ierr)
        call MPI_RSEND(ARRAY(:,:,1:nblocks),msgLen, &
             mpi_integer, dst_task, msgTag, MPI_COMM_ICE, ierr)
#else
        call MPI_SEND(ARRAY(:,:,1:nblocks),msgLen, &
             mpi_integer, dst_task, msgTag, MPI_COMM_ICE, ierr)
#endif
     endif

#else

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then

         nsends = nsends + 1
         src_block = src_dist%blockLocalID(n)
#ifdef _USE_FLOW_CONTROL
         call MPI_RECV(signal, 1, mpi_integer, dst_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)
         call MPI_IRSEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpi_integer, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
#else
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpi_integer, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
#endif
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)
#endif

   endif

!-----------------------------------------------------------------------

 end subroutine gather_global_int

!EOC
!***********************************************************************
!BOP
! !IROUTINE: scatter_global
! !INTERFACE:

 subroutine scatter_global_dbl(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

! !DESCRIPTION:
!  This subroutine scatters a global-sized array to a distributed array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface scatter_global.

! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (dbl_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !OUTPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
   if (this_block%tripoleTFlag) then
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 2
        yoffset = 0
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 1
        yoffset = -1
     case (field_loc_Eface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_Nface)    ! cell face location
        xoffset = 2
        yoffset = -1
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   else
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 1
        yoffset = 1
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 0
        yoffset = 0
     case (field_loc_Eface)    ! cell face location
        xoffset = 0
        yoffset = 1
     case (field_loc_Nface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   endif

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call abort_ice('Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) > 0 .and. &
           dst_dist%blockLocation(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  ! for yoffset=0 or 1, yoffset2=0,0
                  ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
                  do yoffset2=0,max(yoffset,0)-yoffset
                    jsrc = ny_global + yoffset + yoffset2 + &
                         (this_block%j_glob(j) + ny_global)
                    do i=1,nx_block
                      if (this_block%i_glob(i) /= 0) then
                         isrc = nx_global + xoffset - this_block%i_glob(i)
                         if (isrc < 1) isrc = isrc + nx_global
                         if (isrc > nx_global) isrc = isrc - nx_global
                         msg_buffer(i,j-yoffset2) = isign * ARRAY_G(isrc,jsrc)
                      endif
                    end do
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpiR8, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  ! for yoffset=0 or 1, yoffset2=0,0
                  ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
                  do yoffset2=0,max(yoffset,0)-yoffset
                    jsrc = ny_global + yoffset + yoffset2 + &
                         (this_block%j_glob(j) + ny_global)
                    do i=1,nx_block
                      if (this_block%i_glob(i) /= 0) then
                         isrc = nx_global + xoffset - this_block%i_glob(i)
                         if (isrc < 1) isrc = isrc + nx_global
                         if (isrc > nx_global) isrc = isrc - nx_global
                         ARRAY(i,j-yoffset2,dst_block) &
                           = isign * ARRAY_G(isrc,jsrc)
                      endif
                    end do
                  end do

               endif
            end do

         endif
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%blockLocalID(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpiR8, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

   !-----------------------------------------------------------------
   ! Ensure unused ghost cell values are 0
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
      do n=1,nblocks_tot
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         if (dst_block > 0) then

         ! north edge
         do j = this_block%jhi+1,ny_block
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! east edge
         do j = 1, ny_block
         do i = this_block%ihi+1,nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! south edge
         do j = 1, this_block%jlo-1
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! west edge
         do j = 1, ny_block
         do i = 1, this_block%ilo-1
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo

         endif
      enddo
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_dbl

!***********************************************************************

 subroutine scatter_global_real(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a global-sized array to a distributed array.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (real_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (real_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (real_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = 0._real_kind

   this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
   if (this_block%tripoleTFlag) then
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 2
        yoffset = 0
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 1
        yoffset = 1
     case (field_loc_Eface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_Nface)    ! cell face location
        xoffset = 2
        yoffset = 1
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   else
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 1
        yoffset = 1
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 0
        yoffset = 0
     case (field_loc_Eface)    ! cell face location
        xoffset = 0
        yoffset = 1
     case (field_loc_Nface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   endif

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call abort_ice('Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) > 0 .and. &
           dst_dist%blockLocation(n)-1 /= my_task) then

         msg_buffer = 0._real_kind
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  ! for yoffset=0 or 1, yoffset2=0,0
                  ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
                  do yoffset2=0,max(yoffset,0)-yoffset
                    jsrc = ny_global + yoffset + yoffset2 + &
                         (this_block%j_glob(j) + ny_global)
                    do i=1,nx_block
                      if (this_block%i_glob(i) /= 0) then
                         isrc = nx_global + xoffset - this_block%i_glob(i)
                         if (isrc < 1) isrc = isrc + nx_global
                         if (isrc > nx_global) isrc = isrc - nx_global
                         msg_buffer(i,j-yoffset2) = isign * ARRAY_G(isrc,jsrc)
                      endif
                    end do
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpiR4, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  ! for yoffset=0 or 1, yoffset2=0,0
                  ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
                  do yoffset2=0,max(yoffset,0)-yoffset
                    jsrc = ny_global + yoffset + yoffset2 + &
                         (this_block%j_glob(j) + ny_global)
                    do i=1,nx_block
                      if (this_block%i_glob(i) /= 0) then
                         isrc = nx_global + xoffset - this_block%i_glob(i)
                         if (isrc < 1) isrc = isrc + nx_global
                         if (isrc > nx_global) isrc = isrc - nx_global
                         ARRAY(i,j-yoffset2,dst_block) &
                           = isign * ARRAY_G(isrc,jsrc)
                      endif
                    end do
                  end do

               endif
            end do

         endif
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%blockLocalID(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpiR4, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

   !-----------------------------------------------------------------
   ! Ensure unused ghost cell values are 0
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
      do n=1,nblocks_tot
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         if (dst_block > 0) then

         ! north edge
         do j = this_block%jhi+1,ny_block
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = 0._real_kind
         enddo
         enddo
         ! east edge
         do j = 1, ny_block
         do i = this_block%ihi+1,nx_block
            ARRAY (i,j,dst_block) = 0._real_kind
         enddo
         enddo
         ! south edge
         do j = 1, this_block%jlo-1
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = 0._real_kind
         enddo
         enddo
         ! west edge
         do j = 1, ny_block
         do i = 1, this_block%ilo-1
            ARRAY (i,j,dst_block) = 0._real_kind
         enddo
         enddo

         endif
      enddo
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_real

!***********************************************************************

 subroutine scatter_global_int(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a global-sized array to a distributed array.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   integer (int_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = 0

   this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
   if (this_block%tripoleTFlag) then
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 2
        yoffset = 0
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 1
        yoffset = 1
     case (field_loc_Eface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_Nface)    ! cell face location
        xoffset = 2
        yoffset = 1
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   else
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 1
        yoffset = 1
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 0
        yoffset = 0
     case (field_loc_Eface)    ! cell face location
        xoffset = 0
        yoffset = 1
     case (field_loc_Nface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   endif

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call abort_ice('Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) > 0 .and. &
           dst_dist%blockLocation(n)-1 /= my_task) then

         msg_buffer = 0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  ! for yoffset=0 or 1, yoffset2=0,0
                  ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
                  do yoffset2=0,max(yoffset,0)-yoffset
                    jsrc = ny_global + yoffset + yoffset2 + &
                         (this_block%j_glob(j) + ny_global)
                    do i=1,nx_block
                      if (this_block%i_glob(i) /= 0) then
                         isrc = nx_global + xoffset - this_block%i_glob(i)
                         if (isrc < 1) isrc = isrc + nx_global
                         if (isrc > nx_global) isrc = isrc - nx_global
                         msg_buffer(i,j-yoffset2) = isign * ARRAY_G(isrc,jsrc)
                      endif
                    end do
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpi_integer, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  ! for yoffset=0 or 1, yoffset2=0,0
                  ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
                  do yoffset2=0,max(yoffset,0)-yoffset
                    jsrc = ny_global + yoffset + yoffset2 + &
                         (this_block%j_glob(j) + ny_global)
                    do i=1,nx_block
                      if (this_block%i_glob(i) /= 0) then
                         isrc = nx_global + xoffset - this_block%i_glob(i)
                         if (isrc < 1) isrc = isrc + nx_global
                         if (isrc > nx_global) isrc = isrc - nx_global
                         ARRAY(i,j-yoffset2,dst_block) &
                           = isign * ARRAY_G(isrc,jsrc)
                      endif
                    end do
                  end do

               endif
            end do

         endif
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%blockLocalID(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpi_integer, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

   !-----------------------------------------------------------------
   ! Ensure unused ghost cell values are 0
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
      do n=1,nblocks_tot
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         if (dst_block > 0) then

         ! north edge
         do j = this_block%jhi+1,ny_block
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = 0
         enddo
         enddo
         ! east edge
         do j = 1, ny_block
         do i = this_block%ihi+1,nx_block
            ARRAY (i,j,dst_block) = 0
         enddo
         enddo
         ! south edge
         do j = 1, this_block%jlo-1
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = 0
         enddo
         enddo
         ! west edge
         do j = 1, ny_block
         do i = 1, this_block%ilo-1
            ARRAY (i,j,dst_block) = 0
         enddo
         enddo

         endif
      enddo
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_int

!EOC
!***********************************************************************
!BOP
! !IROUTINE: scatter_global_stress
! !INTERFACE:

 subroutine scatter_global_stress(ARRAY, ARRAY_G1, ARRAY_G2, &
                                  src_task, dst_dist)

! !DESCRIPTION:
!  This subroutine scatters global stresses to a distributed array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  Ghost cells in the stress tensor must be handled separately on tripole
!  grids, because matching the corner values requires 2 different arrays.

! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (dbl_kind), dimension(:,:), intent(in) :: &
     ARRAY_G1,     &! array containing global field on src_task
     ARRAY_G2       ! array containing global field on src_task

! !OUTPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
   if (this_block%tripoleTFlag) then
     xoffset = 2  ! treat stresses as cell-centered scalars (they are not 
     yoffset = 0  ! shared with neighboring grid cells)
   else
     xoffset = 1  ! treat stresses as cell-centered scalars (they are not 
     yoffset = 1  ! shared with neighboring grid cells)
   endif
   isign   = 1

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) > 0 .and. &
           dst_dist%blockLocation(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G1(this_block%i_glob(i),&
                                          this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G1(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G1(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = isign * ARRAY_G2(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpiR8, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G1(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G1(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G1(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  ! for yoffset=0 or 1, yoffset2=0,0
                  ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
                  do yoffset2=0,max(yoffset,0)-yoffset
                    jsrc = ny_global + yoffset + yoffset2 + &
                         (this_block%j_glob(j) + ny_global)
                    do i=1,nx_block
                      if (this_block%i_glob(i) /= 0) then
                         isrc = nx_global + xoffset - this_block%i_glob(i)
                         if (isrc < 1) isrc = isrc + nx_global
                         if (isrc > nx_global) isrc = isrc - nx_global
                         ARRAY(i,j-yoffset2,dst_block) &
                           = isign * ARRAY_G2(isrc,jsrc)
                      endif
                    end do
                  end do

               endif
            end do

         endif
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%blockLocalID(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpiR8, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_stress

!***********************************************************************

 end module ice_gather_scatter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
