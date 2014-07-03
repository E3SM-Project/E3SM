!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ice_distribution

 module ice_distribution

!
! !DESCRIPTION:
!  This module provides data types and routines for distributing
!  blocks across processors.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_distribution.F90 118 2008-04-08 20:57:17Z eclare $
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP by William H. Lipscomb, LANL
! Jan. 2008: Elizabeth Hunke updated to new POP infrastructure
!
! !USES:

   use ice_kinds_mod
   use ice_domain_size
   use ice_communicate
   use ice_blocks
   use ice_exit
   use ice_fileunits, only: nu_diag, nu_timing, ice_stdout, flush_fileunit
   use ice_spacecurve
   use ice_broadcast
!   use ice_probability, only:  EstimateCost, BuildProbabilityStats, WoriteProbabilityStats
   use ice_probability_tools

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: distrb  ! distribution data type
      integer (int_kind) :: &
         nprocs            ,&! number of processors in this dist
         communicator      ,&! communicator to use in this dist
         numLocalBlocks      ! number of blocks distributed to this
                             !   local processor

      integer (int_kind), dimension(:), pointer :: &
         blockLocation     ,&! processor location for all blocks
         blockLocalID      ,&! local  block id for all blocks
         blockGlobalID       ! global block id for each local block

      integer (int_kind), dimension(:), pointer ::  blockCnt
      integer (int_kind), dimension(:,:), pointer :: blockIndex

   end type

   integer (int_kind), public, parameter ::  &    ! types of blocks:
         lndType     = 0,                    &    ! 	Land
         icefreeType = 1,                    &    !     ice free (ocean only)
         iceType     = 2                          !     sea ice 

!    integer, parameter :: numCoeff = 5



! !PUBLIC MEMBER FUNCTIONS:

   public :: create_distribution, &
             ice_distributionGet,         &
             ice_distributionGetBlockLoc, &
             ice_distributionGetBlockID, &
             create_local_block_ids

! !PUBLIC DATA MEMBERS:

   character (char_len), public :: &
       processor_shape       ! 'square-pop' (approx) POP default config
                             ! 'square-ice' like square-pop but better for ice
                             ! 'slenderX1' (NPX x 1)
                             ! 'slenderX2' (NPX x 2)

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: create_distribution
! !INTERFACE:

 function create_distribution(dist_type, nprocs, minBlock, maxBlock, work_per_block, prob_per_block, blockType, bStats, &
      FixMaxBlock, maxDil)

! !DESCRIPTION:
!  This routine determines the distribution of blocks across processors
!  by call the appropriate subroutine based on distribution type
!  requested.  Currently three distributions are supported:
!  2-d Cartesian distribution (cartesian), a load-balanced
!  distribution using a rake algorithm based on an input amount of work 
!  per block, and a space-filling-curve algorithm.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      dist_type             ! method for distributing blocks
                            !  either cartesian or rake

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), intent(in) :: minBlock ! minimum number of blocks to use
   integer (int_kind), intent(in) :: maxBlock ! maximum number of blocks to use

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

   real (dbl_kind), dimension(:), intent(in) :: &
      prob_per_block        ! probability that block contains sea-ice

   integer (int_kind), dimension(:), intent(in) :: &
      blockType             ! type of block

   real (dbl_kind), dimension(:,:), intent(in)  :: bStats  ! block statistics 

   logical, intent(in) :: FixMaxBlock
   real (real_kind), intent(in) :: maxDil

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distribution   ! resulting structure describing
                            !  distribution of blocks

!EOP

!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, nb, nc        ! dummy counters

!----------------------------------------------------------------------
!
!  select the appropriate distribution type
!
!----------------------------------------------------------------------

   select case (trim(dist_type))

   case('roundrobin')

      create_distribution = create_distrb_roundrobin(nprocs, work_per_block)

   case('blkrobin')

      create_distribution = create_distrb_blkrobin(nprocs, work_per_block)

   case('blkcart')

      create_distribution = create_distrb_blkcart(nprocs, work_per_block)

   case('cartesian')

      create_distribution = create_distrb_cart(nprocs, work_per_block)

   case('rake')

      create_distribution = create_distrb_rake(nprocs, work_per_block)

   case('spacecurve')

!DBG      write(nu_diag,*) 'before call to create_distrb_spacecurve'
      create_distribution = create_distrb_spacecurve(nprocs, minBlock, maxBlock, &
                            work_per_block,prob_per_block,blockType,bStats, FixMaxBlock, maxDil )
!DBG    stop 'create_distribution: after call to create_distrb_spacecurve'

   case default

      call abort_ice('ice distribution: unknown distribution type')

   end select

   if(my_task == master_task) then 
      write(nu_diag,*) ' '
      nc = 0
      do n = 1,nprocs
         nb = create_distribution%blockcnt(n)
         if (nb > 0) then
            nc = nc + 1
!            write(nu_diag,*) ' Blocks on proc : ',n,nb,':', &
!               create_distribution%blockindex(n,1:nb)
         else
!            write(nu_diag,*) ' Blocks on proc : ',n,nb
         endif
      enddo
      write(nu_diag,*) ' Active processors: ',nc,MAXVAL(create_distribution%blockLocation)
      write(nu_diag,*) ' '
   endif
!DBG write(nu_diag,*) 'end of create_distribution'
!-----------------------------------------------------------------------
!EOC

 end function create_distribution

!***********************************************************************
!BOP
! !IROUTINE: create_local_block_ids
! !INTERFACE:

 subroutine create_local_block_ids(block_ids, distribution)

! !DESCRIPTION:
!  This routine determines which blocks in an input distribution are
!  located on the local processor and creates an array of block ids
!  for all local blocks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which local
                             !  blocks required

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), pointer :: &
      block_ids              ! array of block ids for every block
                             ! that resides on the local processor
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, bid, bcount        ! dummy counters

   logical (log_kind) :: dbug

!-----------------------------------------------------------------------
!
!  first determine number of local blocks to allocate array
!
!-----------------------------------------------------------------------

   bcount = 0
   do n=1,size(distribution%blockLocation)
      if (distribution%blockLocation(n) == my_task+1) bcount = bcount + 1
   end do


   if (bcount > 0) allocate(block_ids(bcount))

!-----------------------------------------------------------------------
!
!  now fill array with proper block ids
!
!-----------------------------------------------------------------------

!   dbug = .true.
   dbug = .false.
   if (bcount > 0) then
      do n=1,size(distribution%blockLocation)
         if (distribution%blockLocation(n) == my_task+1) then
            block_ids(distribution%blockLocalID(n)) = n

            if (dbug) then
            write(nu_diag,*) 'block id, proc, local_block: ', &
                             block_ids(distribution%blockLocalID(n)), &
                             distribution%blockLocation(n), &
                             distribution%blockLocalID(n)
            endif
         endif
      end do
   endif

!EOC

 end subroutine create_local_block_ids

!***********************************************************************
!BOP
! !IROUTINE: proc_decomposition
! !INTERFACE:

 subroutine proc_decomposition(nprocs, nprocs_x, nprocs_y)

! !DESCRIPTION:
!  This subroutine attempts to find an optimal (nearly square)
!  2d processor decomposition for a given number of processors.
!
! !REVISION HISTORY:
!  same as module
!
! !USES:

   use ice_domain_size

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs                       ! total number or processors

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      nprocs_x, nprocs_y           ! number of procs in each dimension

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      iguess, jguess               ! guesses for nproc_x,y

   real (real_kind) :: &
      square                       ! square root of nprocs

!----------------------------------------------------------------------
!
!  start with an initial guess
!
!----------------------------------------------------------------------

   square = sqrt(real(nprocs,kind=real_kind))
   nprocs_x = 0
   nprocs_y = 0

   if (processor_shape == 'square-pop') then ! make as square as possible
      iguess = nint(square)
      jguess = nprocs/iguess
   elseif (processor_shape == 'square-ice') then ! better for bipolar ice
      jguess = nint(square)
      iguess = nprocs/jguess
   elseif (processor_shape == 'slenderX1') then ! 1 proc in y direction
      jguess = 1
      iguess = nprocs/jguess
   elseif (processor_shape == 'slenderX2') then ! 2 proc in y direction
      jguess = min(2, nprocs)
      iguess = nprocs/jguess
   else                                  ! abort
      call abort_ice('ice: processor_shape not supported, '//trim(processor_shape))
   endif

!----------------------------------------------------------------------
!
!  try various decompositions to find the best
!
!----------------------------------------------------------------------

   proc_loop: do
   if (processor_shape == 'square-pop') then
      jguess = nprocs/iguess
   else
      iguess = nprocs/jguess
   endif

      if (iguess*jguess == nprocs) then ! valid decomp

         !*** if the blocks can be evenly distributed, it is a
         !*** good decomposition
         if (mod(nblocks_x,iguess) == 0 .and. &
             mod(nblocks_y,jguess) == 0) then
            nprocs_x = iguess
            nprocs_y = jguess
            exit proc_loop

         !*** if the blocks can be evenly distributed in a
         !*** transposed direction, it is a good decomposition
         else if (mod(nblocks_x,jguess) == 0 .and. &
                mod(nblocks_y,iguess) == 0) then
            nprocs_x = jguess
            nprocs_y = iguess
            exit proc_loop

         !*** A valid decomposition, but keep searching for
         !***  a better one
         else
            if (nprocs_x == 0) then
               nprocs_x = iguess
               nprocs_y = jguess
            endif
            if (processor_shape == 'square-pop') then
               iguess = iguess - 1
               if (iguess == 0) then
                  exit proc_loop
               else
                  cycle proc_loop
               endif
            else
               jguess = jguess - 1
               if (jguess == 0) then
                  exit proc_loop
               else
                  cycle proc_loop
               endif
            endif
         endif

      else ! invalid decomp - keep trying

         if (processor_shape == 'square-pop') then
            iguess = iguess - 1
            if (iguess == 0) then
               exit proc_loop
            else
               cycle proc_loop
            endif
         else
            jguess = jguess - 1
            if (jguess == 0) then
               exit proc_loop
            else
               cycle proc_loop
            endif
         endif
      endif

   end do proc_loop

   if (nprocs_x == 0) then
      call abort_ice('ice: Unable to find 2d processor config')
   endif

   if (my_task == master_task) then
     write(nu_diag,'(a23,i4,a3,i4)') '  Processors (X x Y) = ', &
                                        nprocs_x,' x ',nprocs_y
   endif

!----------------------------------------------------------------------
!EOC

 end subroutine proc_decomposition

!**********************************************************************
!BOP
! !IROUTINE: ice_distributionDestroy
! !INTERFACE:

 subroutine ice_distributionDestroy(distribution)

! !DESCRIPTION:
!  This routine destroys a defined distribution by deallocating
!  all memory associated with the distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (distrb), intent(inout) :: &
      distribution          ! distribution to destroy

! !OUTPUT PARAMETERS:

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: istat  ! status flag for deallocate

!----------------------------------------------------------------------
!
!  reset scalars
!
!----------------------------------------------------------------------

   distribution%nprocs       = 0
   distribution%communicator   = 0
   distribution%numLocalBlocks = 0

!----------------------------------------------------------------------
!
!  deallocate arrays
!
!----------------------------------------------------------------------

   deallocate(distribution%blockLocation, stat=istat)
   deallocate(distribution%blockLocalID , stat=istat)
   deallocate(distribution%blockGlobalID, stat=istat)
   deallocate(distribution%blockCnt, stat=istat)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_distributionDestroy

!***********************************************************************
!BOP
! !IROUTINE: ice_distributionGet
! !INTERFACE:

 subroutine ice_distributionGet(distribution,&
                            nprocs, communicator, numLocalBlocks, &
                            blockLocation, blockLocalID, blockGlobalID)


! !DESCRIPTION:
!  This routine extracts information from a distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

! !OUTPUT PARAMETERS:

      integer (int_kind), intent(out), optional ::   &
         nprocs          ,&! number of processors in this dist
         communicator      ,&! communicator to use in this dist
         numLocalBlocks      ! number of blocks distributed to this
                             !   local processor

      integer (int_kind), dimension(:), pointer, optional :: &
         blockLocation     ,&! processor location for all blocks
         blockLocalID      ,&! local  block id for all blocks
         blockGlobalID       ! global block id for each local block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  depending on which optional arguments are present, extract the
!  requested info
!
!-----------------------------------------------------------------------

   if (present(nprocs))       nprocs       = distribution%nprocs
   if (present(communicator))   communicator   = distribution%communicator
   if (present(numLocalBlocks)) numLocalBlocks = distribution%numLocalBlocks

   if (present(blockLocation)) then
      if (associated(distribution%blockLocation)) then
         blockLocation => distribution%blockLocation
      else
        call abort_ice( &
            'ice_distributionGet: blockLocation not allocated')
         return
      endif
   endif

   if (present(blockLocalID)) then
      if (associated(distribution%blockLocalID)) then
         blockLocalID = distribution%blockLocalID
      else
        call abort_ice( &
            'ice_distributionGet: blockLocalID not allocated')
         return
      endif
   endif

   if (present(blockGlobalID)) then
      if (associated(distribution%blockGlobalID)) then
         blockGlobalID = distribution%blockGlobalID
      else
        call abort_ice( &
            'ice_distributionGet: blockGlobalID not allocated')
         return
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_distributionGet

!***********************************************************************
!BOP
! !IROUTINE: ice_distributionGetBlockLoc
! !INTERFACE:

 subroutine ice_distributionGetBlockLoc(distribution, blockID, &
                                        processor, localID)


! !DESCRIPTION:
!  Given a distribution of blocks and a global block ID, return
!  the processor and local index for the block.  A zero for both
!  is returned in the case that the block has been eliminated from
!  the distribution (i.e. has no active points).
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

   integer (int_kind), intent(in) :: &
      blockID                ! global block id for which location requested

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) ::  &
      processor,            &! processor on which block resides
      localID                ! local index for this block on this proc

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for valid blockID
!
!-----------------------------------------------------------------------

   if (blockID < 0 .or. blockID > nblocks_tot) then
     call abort_ice( &
         'ice_distributionGetBlockLoc: invalid block id')
      return
   endif

!-----------------------------------------------------------------------
!
!  extract the location from the distribution data structure
!
!-----------------------------------------------------------------------

   processor = distribution%blockLocation(blockID)
   localID   = distribution%blockLocalID (blockID)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_distributionGetBlockLoc

!***********************************************************************
!BOP
! !IROUTINE: ice_distributionGetBlockID
! !INTERFACE:

 subroutine ice_distributionGetBlockID(distribution, localID, &
                                       blockID)


! !DESCRIPTION:
!  Given a distribution of blocks and a local block index, return
!  the global block id for the block.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

   integer (int_kind), intent(in) ::  &
      localID                ! local index for this block on this proc

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      blockID                ! global block id for this local block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for valid localID
!
!-----------------------------------------------------------------------

   if (localID < 0 .or. localID > distribution%numLocalBlocks) then
     call abort_ice( &
         'ice_distributionGetBlockID: invalid local id')
      return
   endif

!-----------------------------------------------------------------------
!
!  extract the global ID from the distribution data structure
!
!-----------------------------------------------------------------------

   blockID   = distribution%blockGlobalID (localID)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_distributionGetBlockID

!***********************************************************************
!BOP
! !IROUTINE: create_distrb_roundrobin
! !INTERFACE:

 function create_distrb_roundrobin(nprocs, workPerBlock) result(newDistrb)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a simple roundrobin algorithm. Mean for prescribed ice or
!  standalone CAM mode.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs            ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j,                  &! dummy loop indices
      istat,                 &! status flag for allocation
      iblock, jblock,        &!
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID                 ! block location on this processor

   integer (int_kind), dimension(:), allocatable :: &
      proc_tmp           ! temp processor id
   
   integer (int_kind) :: pid,n

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (newDistrb%blockLocation(nblocks_tot), &
             newDistrb%blockLocalID (nblocks_tot), stat=istat)

   allocate (newDistrb%blockCnt(nprocs))
!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   allocate(proc_tmp(nprocs))
   processor = 0
   globalID = 0
   proc_tmp = 0

   allocate(newDistrb%blockIndex(nprocs,max_blocks))
   newDistrb%blockIndex(:,:) = 0

   do j=1,nblocks_y
   do i=1,nblocks_x
      
      globalID = globalID + 1

      if (workPerBlock(globalID) /= 0) then
         processor = mod(processor,nprocs) + 1
         proc_tmp(processor) = proc_tmp(processor) + 1
         localID = proc_tmp(processor)
         newDistrb%blockLocation(globalID) = processor
         newDistrb%blockLocalID (globalID) = localID
         newDistrb%blockIndex(processor,localID) = globalID
      else  ! no work - eliminate block from distribution
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif

   end do
   end do

   newDistrb%numLocalBlocks = proc_tmp(my_task+1)
   newDistrb%blockCnt(:) = proc_tmp(:)
   deallocate(proc_tmp)

!   write(nu_diag,*) 'my_task,newDistrb%numLocalBlocks',&
!      my_task,newDistrb%numLocalBlocks

!----------------------------------------------------------------------
!
!  now store the local info
!
!----------------------------------------------------------------------

   globalID = 0

   if (newDistrb%numLocalBlocks > 0) then
      allocate (newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
                stat=istat)

      processor = my_task + 1
      do localID = 1,newDistrb%numLocalBlocks
         newDistrb%blockGlobalID (localID) = newDistrb%blockIndex(processor,&
                                             localID)
      enddo
   endif

!----------------------------------------------------------------------
!EOC
 end function create_distrb_roundrobin
 
!***********************************************************************
!BOP
! !IROUTINE: create_distrb_blkcart
! !INTERFACE:

 function create_distrb_blkcart(nprocs, workPerBlock) result(newDistrb)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a simple blkcart algorithm. Mean for prescribed ice or
!  standalone CAM mode.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs            ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, i2, j2,          &! dummy loop indices
      istat,                 &! status flag for allocation
      iblock, jblock,        &!
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID,               &! block location on this processor
      blktogether,           &! number of blocks together
      cnt                     ! counter

   integer (int_kind), dimension(:), allocatable :: &
      proc_tmp           ! temp processor id
   
   integer (int_kind) :: pid,n

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (newDistrb%blockLocation(nblocks_tot), &
             newDistrb%blockLocalID (nblocks_tot), stat=istat)

   allocate (newDistrb%blockCnt(nprocs))
!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   allocate(proc_tmp(nprocs))
   proc_tmp = 0

   allocate(newDistrb%blockIndex(nprocs,max_blocks))
   newDistrb%blockIndex(:,:) = 0

   blktogether = max(1,nint(float(nblocks_x*nblocks_y)/float(4*nprocs)))

   ! --- two phases, resetnt processor and cnt for each phase
   ! --- phase 1 is south to north, east to west on the left half of the domain
   ! --- phase 2 is north to south, east to west on the right half of the domain

   if (mod(nblocks_x,2) /= 0) then
      call abort_ice( &
         'create_distrb_blkcart: nblocks_x not divisible by 2')
      return
   endif

   do n=1,2
   processor = 1
   cnt = 0
   do j2=1,nblocks_y
   do i2=1,nblocks_x/2
      
      if (n == 1) then
         i = i2
         j = j2
      else
         i = nblocks_x/2 + i2
         j = nblocks_y - j2 + 1
      endif

      globalID = (j-1)*nblocks_x + i
      if (cnt >= blktogether) then
         processor = mod(processor,nprocs) + 1
         cnt = 0
      endif
      cnt = cnt + 1

      if (workPerBlock(globalID) /= 0) then
         proc_tmp(processor) = proc_tmp(processor) + 1
         localID = proc_tmp(processor)
         newDistrb%blockLocation(globalID) = processor
         newDistrb%blockLocalID (globalID) = localID
         newDistrb%blockIndex(processor,localID) = globalID
      else  ! no work - eliminate block from distribution
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif

   end do
   end do
   end do

   newDistrb%numLocalBlocks = proc_tmp(my_task+1)
   newDistrb%blockCnt(:) = proc_tmp(:)
   deallocate(proc_tmp)

!   write(nu_diag,*) 'my_task,newDistrb%numLocalBlocks',&
!      my_task,newDistrb%numLocalBlocks

!----------------------------------------------------------------------
!
!  now store the local info
!
!----------------------------------------------------------------------

   globalID = 0

   if (newDistrb%numLocalBlocks > 0) then
      allocate (newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
                stat=istat)

      processor = my_task + 1
      do localID = 1,newDistrb%numLocalBlocks
         newDistrb%blockGlobalID (localID) = newDistrb%blockIndex(processor,&
                                             localID)
      enddo
   endif

!----------------------------------------------------------------------
!EOC
 end function create_distrb_blkcart
 
!***********************************************************************
!BOP
! !IROUTINE: create_distrb_blkrobin
! !INTERFACE:

 function create_distrb_blkrobin(nprocs, workPerBlock) result(newDistrb)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a simple blkrobin algorithm. Mean for prescribed ice or
!  standalone CAM mode.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs            ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j,                  &! dummy loop indices
      istat,                 &! status flag for allocation
      iblock, jblock,        &!
      mblocks,               &! estimate of max blocks per pe
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID                 ! block location on this processor

   integer (int_kind), dimension(:), allocatable :: &
      proc_tmp           ! temp processor id

   logical (log_kind), dimension(:), allocatable :: &
      bfree              ! map of assigned blocks
   
   integer (int_kind) :: pid,n,cnt, blktogether, i2, j2
   integer (int_kind) :: totblocks, nchunks
   logical (log_kind) :: keepgoing

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (newDistrb%blockLocation(nblocks_tot), &
             newDistrb%blockLocalID (nblocks_tot), stat=istat)

   allocate (newDistrb%blockCnt(nprocs))
!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   allocate(proc_tmp(nprocs))
   processor = 0
   globalID = 0
   proc_tmp = 0

   allocate(newDistrb%blockIndex(nprocs,max_blocks))
   newDistrb%blockIndex(:,:) = 0

   allocate(bfree(nblocks_x*nblocks_y))
   bfree=.true.

   totblocks = 0
   do j=1,nblocks_y
   do i=1,nblocks_x
      globalID = (j-1)*nblocks_x + i
      if (workPerBlock(globalID) /= 0) then
         totblocks=totblocks+1
      else  ! no work - eliminate block from distribution
         bfree(globalID) = .false.
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif
   enddo
   enddo

   mblocks = totblocks/nprocs
   if (mod(totblocks,nprocs) > 0) mblocks=mblocks+1

   blktogether = max(1,nint(float(totblocks)/float(6*nprocs)))

!   write(nu_diag,*) 'ice_distrb_blkrobin totblocks = ',totblocks,nblocks_y*nblocks_x
 
   !------------------------------
   ! southern group of blocks
   !   weave back and forth in i vs j
   !   go south to north, low - high pes
   !------------------------------

   processor=1
   cnt = 0
   keepgoing = .true.
   do j=1,nblocks_y
   do i=1,nblocks_x
      if (mod(j,2) == 0) then
         i2 = nblocks_x - i + 1
      else
         i2 = i
      endif
      globalID = (j-1)*nblocks_x + i2
      if (cnt >= blktogether) then
         processor = mod(processor,nprocs) + 1
         cnt = 0
         if (processor == 1) keepgoing = .false.
      endif
!      write(nu_diag,'(a,6i7,l2)') 'tcx ',i,j,globalID,cnt,blktogether,processor,keepgoing

      if (keepgoing) then
         if (bfree(globalID)) then
         if (workPerBlock(globalID) /= 0) then
            proc_tmp(processor) = proc_tmp(processor) + 1
            localID = proc_tmp(processor)
            newDistrb%blockLocation(globalID) = processor
            newDistrb%blockLocalID (globalID) = localID
            newDistrb%blockIndex(processor,localID) = globalID
            cnt = cnt + 1
            totblocks = totblocks-1
            bfree(globalID) = .false.

         else  ! no work - eliminate block from distribution
            bfree(globalID) = .false.
            newDistrb%blockLocation(globalID) = 0
            newDistrb%blockLocalID (globalID) = 0
         endif
         endif  ! bfree
      endif
   end do
   end do

!   write(nu_diag,*) 'ice_distrb_blkrobin totblocks left after southern = ',totblocks

   !------------------------------
   ! northern group of blocks
   !   weave back and forth in i vs j
   !   go north to south, high - low pes
   !------------------------------

   processor=nprocs
   cnt = 0
   keepgoing = .true.
   do j=nblocks_y,1,-1
   do i=1,nblocks_x
      if (mod(j,2) == 1) then
         i2 = nblocks_x - i + 1
      else
         i2 = i
      endif
      globalID = (j-1)*nblocks_x + i2
      if (cnt >= blktogether) then
         processor = mod(processor+nprocs-2,nprocs) + 1
         cnt = 0
         if (processor == nprocs) keepgoing = .false.
      endif

      if (keepgoing) then
         if (bfree(globalID)) then
         if (workPerBlock(globalID) /= 0) then
            proc_tmp(processor) = proc_tmp(processor) + 1
            localID = proc_tmp(processor)
            newDistrb%blockLocation(globalID) = processor
            newDistrb%blockLocalID (globalID) = localID
            newDistrb%blockIndex(processor,localID) = globalID
            cnt = cnt + 1
            totblocks = totblocks - 1
            bfree(globalID) = .false.

         else  ! no work - eliminate block from distribution
            bfree(globalID) = .false.
            newDistrb%blockLocation(globalID) = 0
            newDistrb%blockLocalID (globalID) = 0
         endif
         endif  ! bfree
      endif
   end do
   end do

!   write(nu_diag,*) 'ice_distrb_blkrobin totblocks left after northern = ',totblocks

   !------------------------------
   ! central group of blocks
   !   weave back and forth in i vs j
   !   go north to south, low - high / low - high pes
   !------------------------------

   nchunks = 2*nprocs
   blktogether = max(1,nint(float(totblocks)/float(nchunks)))
   processor=1
   cnt = 0
   do j=nblocks_y,1,-1
   do i=1,nblocks_x
      if (mod(j,2) == 1) then
         i2 = nblocks_x - i + 1
      else
         i2 = i
      endif
      globalID = (j-1)*nblocks_x + i2
      if (totblocks > 0) then
      do while (proc_tmp(processor) >= mblocks .or. cnt >= blktogether)
         nchunks = nchunks - 1
         if (nchunks == 0) then
            blktogether = 1
         else
            blktogether = max(1,nint(float(totblocks)/float(nchunks)))
         endif
         cnt = 0
         processor = mod(processor,nprocs) + 1
      enddo
      endif

!      write(nu_diag,*) 'ice_distrb_blkrobin central ',i,j,totblocks,cnt,nchunks,blktogether,processor

      if (bfree(globalID)) then
      if (workPerBlock(globalID) /= 0) then
         proc_tmp(processor) = proc_tmp(processor) + 1
         localID = proc_tmp(processor)
         newDistrb%blockLocation(globalID) = processor
         newDistrb%blockLocalID (globalID) = localID
         newDistrb%blockIndex(processor,localID) = globalID
         cnt = cnt + 1
         totblocks = totblocks-1
         bfree(globalID) = .false.

      else  ! no work - eliminate block from distribution
         bfree(globalID) = .false.
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif
      endif  ! bfree
   end do
   end do

   newDistrb%numLocalBlocks = proc_tmp(my_task+1)
   newDistrb%blockCnt(:) = proc_tmp(:)
   deallocate(proc_tmp)
   deallocate(bfree)

!----------------------------------------------------------------------
!
!  now store the local info
!
!----------------------------------------------------------------------

   globalID = 0

   if (newDistrb%numLocalBlocks > 0) then
      allocate (newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
                stat=istat)

      processor = my_task + 1
      do localID = 1,newDistrb%numLocalBlocks
         newDistrb%blockGlobalID (localID) = newDistrb%blockIndex(processor,&
                                             localID)
      enddo
   endif

!----------------------------------------------------------------------
!EOC
 end function create_distrb_blkrobin
 
!***********************************************************************
!BOP
! !IROUTINE: create_distrb_cart
! !INTERFACE:

 function create_distrb_cart(nprocs, workPerBlock) result(newDistrb)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a 2-d Cartesian distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs            ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j,                  &! dummy loop indices
      istat,                 &! status flag for allocation
      iblock, jblock,        &!
      is, ie, js, je,        &! start, end block indices for each proc
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID,               &! block location on this processor
      nprocsX,             &! num of procs in x for global domain
      nprocsY,             &! num of procs in y for global domain
      numBlocksXPerProc,     &! num of blocks per processor in x
      numBlocksYPerProc       ! num of blocks per processor in y

   integer (int_kind), dimension(:), allocatable :: &
      proc_tmp           ! temp processor id
   
   integer (int_kind) :: pid,n

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

   ! This assumes that we really do know what we are doing.
   if (processor_shape == 'blocks') then
      nprocsX = nblocks_x / max_blocks
      nprocsY = nblocks_y / max_blocks
   else
      call proc_decomposition(nprocs, nprocsX, nprocsY)
   endif
                                  

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (newDistrb%blockLocation(nblocks_tot), &
             newDistrb%blockLocalID (nblocks_tot), stat=istat)

   allocate (newDistrb%blockCnt(nprocs))
!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   numBlocksXPerProc = (nblocks_x-1)/nprocsX + 1
   numBlocksYPerProc = (nblocks_y-1)/nprocsY + 1

   do j=1,nprocsY
   do i=1,nprocsX
      processor = (j-1)*nprocsX + i    ! number the processors 
                                         ! left to right, bot to top

      is = (i-1)*numBlocksXPerProc + 1   ! starting block in i
      ie =  i   *numBlocksXPerProc       ! ending   block in i
      if (ie > nblocks_x) ie = nblocks_x
      js = (j-1)*numBlocksYPerProc + 1   ! starting block in j
      je =  j   *numBlocksYPerProc       ! ending   block in j
      if (je > nblocks_y) je = nblocks_y

      localID        = 0  ! initialize counter for local index
      do jblock = js,je
      do iblock = is,ie
         globalID = (jblock - 1)*nblocks_x + iblock
         if (workPerBlock(globalID) /= 0) then
            localID = localID + 1
            newDistrb%blockLocation(globalID) = processor
            newDistrb%blockLocalID (globalID) = localID
         else  ! no work - eliminate block from distribution
            newDistrb%blockLocation(globalID) = 0
            newDistrb%blockLocalID (globalID) = 0
         endif
      end do
      end do

      ! if this is the local processor, set number of local blocks
      if (my_task == processor - 1) then
         newDistrb%numLocalBlocks = localID
      endif

   end do
   end do

!----------------------------------------------------------------------
!
!  now store the local info
!
!----------------------------------------------------------------------

   if (newDistrb%numLocalBlocks > 0) then
      allocate (newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
                stat=istat)

      do j=1,nprocsY
      do i=1,nprocsX
         processor = (j-1)*nprocsX + i

         if (processor == my_task + 1) then
            is = (i-1)*numBlocksXPerProc + 1   ! starting block in i
            ie =  i   *numBlocksXPerProc       ! ending   block in i
            if (ie > nblocks_x) ie = nblocks_x
            js = (j-1)*numBlocksYPerProc + 1   ! starting block in j
            je =  j   *numBlocksYPerProc       ! ending   block in j
            if (je > nblocks_y) je = nblocks_y

            localID        = 0  ! initialize counter for local index
            do jblock = js,je
            do iblock = is,ie
               globalID = (jblock - 1)*nblocks_x + iblock
               if (workPerBlock(globalID) /= 0) then
                  localID = localID + 1
                  newDistrb%blockGlobalID (localID) = globalID
               endif
            end do
            end do

         endif

      end do
      end do

   endif

   allocate(proc_tmp(nprocs))
   proc_tmp = 0

   allocate(newDistrb%blockIndex(nprocs,max_blocks))
   newDistrb%blockIndex(:,:) = 0

   do n=1,nblocks_tot
      pid = newDistrb%blockLocation(n)
      if(pid>0) then
        proc_tmp(pid) = proc_tmp(pid) + 1
        if(proc_tmp(pid) <= max_blocks) then 
            newDistrb%blockIndex(pid,proc_tmp(pid)) = n
        endif
      endif
   enddo

   newDistrb%blockCnt(:) = proc_tmp(:)
   deallocate(proc_tmp)


!----------------------------------------------------------------------
!EOC

 end function create_distrb_cart

!**********************************************************************
!BOP
! !IROUTINE: create_distrb_rake
! !INTERFACE:

 function create_distrb_rake(nprocs, workPerBlock) result(newDistrb)

! !DESCRIPTION:
!  This  function distributes blocks across processors in a
!  load-balanced manner based on the amount of work per block.
!  A rake algorithm is used in which the blocks are first distributed
!  in a Cartesian distribution and then a rake is applied in each
!  Cartesian direction.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      newDistrb           ! resulting structure describing
                          ! load-balanced distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) ::    &
      i,j,n              ,&! dummy loop indices
      pid                ,&! dummy for processor id
      istat              ,&! status flag for allocates
      localBlock         ,&! local block position on processor
      numOcnBlocks       ,&! number of ocean blocks
      maxWork            ,&! max amount of work in any block
      nprocsX          ,&! num of procs in x for global domain
      nprocsY            ! num of procs in y for global domain

   integer (int_kind), dimension(:), allocatable :: &
      priority           ,&! priority for moving blocks
      workTmp            ,&! work per row or column for rake algrthm
      procTmp              ! temp processor id for rake algrthm

   type (distrb) :: dist  ! temp hold distribution

!----------------------------------------------------------------------
!
!  first set up as Cartesian distribution
!
!----------------------------------------------------------------------

   dist = create_distrb_cart(nprocs, workPerBlock)
                                    
!----------------------------------------------------------------------
!
!  if the number of blocks is close to the number of processors,
!  only do a 1-d rake on the entire distribution
!
!----------------------------------------------------------------------

   numOcnBlocks = count(workPerBlock /= 0)

   if (numOcnBlocks <= 2*nprocs) then

      allocate(priority(nblocks_tot), stat=istat)

      !*** initialize priority array

      do j=1,nblocks_y
      do i=1,nblocks_x
         n=(j-1)*nblocks_x + i
         if (workPerBlock(n) > 0) then
            priority(n) = maxWork + n - workPerBlock(n)
         else
            priority(n) = 0
         endif
      end do
      end do

      allocate(workTmp(nblocks_tot), procTmp(nblocks_tot), stat=istat)

      workTmp(:) = 0
      do i=1,nprocs
         procTmp(i) = i
         do n=1,nblocks_tot
            if (dist%blockLocation(n) == i) then
               workTmp(i) = workTmp(i) + workPerBlock(n)
            endif
         end do
      end do

      call ice_distributionRake (workTmp, procTmp, workPerBlock, &
                                 priority, dist)

      deallocate(workTmp, procTmp, stat=istat)

!----------------------------------------------------------------------
!
!  otherwise re-distribute blocks using a rake in each direction
!
!----------------------------------------------------------------------

   else

      maxWork = maxval(workPerBlock)

      call proc_decomposition(dist%nprocs, nprocsX, nprocsY)

!----------------------------------------------------------------------
!
!     load-balance using a rake algorithm in the x-direction first
!
!----------------------------------------------------------------------

      allocate(priority(nblocks_tot), stat=istat)

      !*** set highest priority such that eastern-most blocks
      !*** and blocks with the least amount of work are
      !*** moved first

      do j=1,nblocks_y
      do i=1,nblocks_x
         n=(j-1)*nblocks_x + i
         if (workPerBlock(n) > 0) then
            priority(n) = (maxWork + 1)*(nblocks_x + i) - &
                          workPerBlock(n)
         else
            priority(n) = 0
         endif
      end do
      end do

      allocate(workTmp(nprocsX), procTmp(nprocsX), stat=istat)

      do j=1,nprocsY

         workTmp(:) = 0
         do i=1,nprocsX
            pid = (j-1)*nprocsX + i
            procTmp(i) = pid
            do n=1,nblocks_tot
               if (dist%blockLocation(n) == pid) then
                  workTmp(i) = workTmp(i) + workPerBlock(n)
               endif
            end do
         end do

         call ice_distributionRake (workTmp, procTmp, workPerBlock, &
                                    priority, dist)
      end do
   
      deallocate(workTmp, procTmp, stat=istat)


!----------------------------------------------------------------------
!
!     use a rake algorithm in the y-direction now
!
!----------------------------------------------------------------------

      !*** set highest priority for northern-most blocks

      do j=1,nblocks_y
      do i=1,nblocks_x
         n=(j-1)*nblocks_x + i
         if (workPerBlock(n) > 0) then
            priority(n) = (maxWork + 1)*(nblocks_y + j) - &
                          workPerBlock(n)
         else
            priority(n) = 0
         endif
      end do
      end do

      allocate(workTmp(nprocsY), procTmp(nprocsY), stat=istat)

      do i=1,nprocsX

         workTmp(:) = 0
         do j=1,nprocsY
            pid = (j-1)*nprocsX + i
            procTmp(j) = pid
            do n=1,nblocks_tot
               if (dist%blockLocation(n) == pid) then
                  workTmp(j) = workTmp(j) + workPerBlock(n)
               endif
            end do
         end do

         call ice_distributionRake (workTmp, procTmp, workPerBlock, &
                                    priority, dist)

      end do

      deallocate(workTmp, procTmp, priority, stat=istat)

   endif  ! 1d or 2d rake

!----------------------------------------------------------------------
!
!  create new distribution with info extracted from the temporary
!  distribution
!
!----------------------------------------------------------------------

   newDistrb%nprocs     = nprocs
   newDistrb%communicator = dist%communicator

   allocate(newDistrb%blockLocation(nblocks_tot), &
            newDistrb%blockLocalID(nblocks_tot), stat=istat)
   allocate (newDistrb%blockCnt(nprocs))

   allocate(procTmp(nprocs), stat=istat)

   allocate(newDistrb%blockIndex(nprocs,max_blocks))
   newDistrb%blockIndex(:,:) = 0

   procTmp = 0
   do n=1,nblocks_tot
      pid = dist%blockLocation(n)  ! processor id
      newDistrb%blockLocation(n) = pid

      if (pid > 0) then
         procTmp(pid) = procTmp(pid) + 1
         newDistrb%blockLocalID (n) = procTmp(pid)
	 if(procTmp(pid) <= max_blocks) then 
            newDistrb%blockIndex(pid,procTmp(pid)) = n
         endif
      else
         newDistrb%blockLocalID (n) = 0
      endif
   end do

   newDistrb%numLocalBlocks = procTmp(my_task+1)

   if (minval(procTmp) < 1) then
      call abort_ice( &
         'create_distrb_rake: processors left with no blocks')
      return
   endif

   newDistrb%blockCnt(:) = procTmp(:) 

   deallocate(procTmp, stat=istat)

   if (istat > 0) then
      call abort_ice( &
         'create_distrb_rake: error allocating last procTmp')
      return
   endif

   allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
            stat=istat)

   if (istat > 0) then
      call abort_ice( &
         'create_distrb_rake: error allocating blockGlobalID')
      return
   endif

   localBlock = 0
   do n=1,nblocks_tot
      if (newDistrb%blockLocation(n) == my_task+1) then
         localBlock = localBlock + 1
         newDistrb%blockGlobalID(localBlock) = n
      endif
   end do



!----------------------------------------------------------------------

   call ice_distributionDestroy(dist)

!----------------------------------------------------------------------
!EOC

 end function create_distrb_rake

!**********************************************************************
!BOP
! !IROUTINE: create_distrb_spacecurve
! !INTERFACE:

 function create_distrb_spacecurve(nprocs, minBlock, maxBlock, work_per_block,prob_per_block,blockType, bStats, &
      FixMaxBlock, maxDil )

! !Description:
!  This function distributes blocks across processors in a
!  load-balanced manner using space-filling curves
!
! !REVISION HISTORY:
!  added by J. Dennis 3/10/06

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), intent(in) :: minBlock, maxBlock

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

   real(dbl_kind), dimension(:), intent(in) :: &
      prob_per_block        ! probability sea-ice within block 

   integer (int_kind), dimension(:), intent(in)  :: &
      blockType            ! type of block

    real (dbl_kind),  dimension(:,:), intent(in) :: bStats

    logical, intent(in) :: FixMaxBlock

     real (real_kind) :: maxDil

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_spacecurve  ! resulting structure describing
                                ! load-balanced distribution of blocks
!EOP
!BOC
!----------------------------------------------------------------------

!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,k,n              ,&! dummy loop indices
      pid                  ,&! dummy for processor id
      localID              ,&! local block position on processor
      max_work             ,&! max amount of work in any block
      nprocs_x             ,&! num of procs in x for global domain
      nprocs_y               ! num of procs in y for global domain

   character(char_len) :: fname

   integer (int_kind) :: maxB

   integer (int_kind), dimension(:),allocatable :: &
        idxT_i,idxT_j,Lindx,Lindx2       ! Temporary indices for SFC

   integer (int_kind), dimension(:,:),allocatable :: &
        Mesh            ,&!   !arrays to hold Space-filling curve
        Mesh2             !

   integer (int_kind) :: &
        nblocksL,nblocks, &! Number of blocks local and total
        ii,extra,tmp1,    &! loop tempories used for
        s1,ig              ! partitioning curve

   logical, parameter :: Debug = .FALSE.

   integer (int_kind), dimension(:), allocatable :: &
      priority           ,&! priority for moving blocks
      work_tmp           ,&! work per row or column for rake algrthm
      proc_tmp           ,&! temp processor id for rake algrthm
      block_count          ! counter to determine local block indx
 
   integer (int_kind), allocatable, dimension(:) ::  &
          blockLocation,   &! block to processor mapping
          distance,        &! location in uncompressed SFC
          type_on_curve,   &! type of blocks 
          work_on_curve

   type (distrb) :: dist  ! temp hold distribution

   integer (int_kind) :: numIce, minblocks
   type (factor_t) :: xdim, ydim
   
   integer (int_kind) :: it,jj,i2,j2
   integer (int_kind) :: curveSize, sb_x, sb_y, itmp,numfac
   integer (int_kind) :: subNum,sfcNum
   logical (log_kind) :: foundX 
   
   integer (int_kind) :: ns,ntmp
   integer (int_kind) :: numLocalBlocks

   real (dbl_kind), allocatable, dimension(:) :: &
       Cost_per_proc , &
       Cost_per_block

   real (dbl_kind), allocatable, dimension(:,:) :: cStats  ! block statistics on SFC curve 

   integer (int_kind), allocatable, dimension(:) :: work_per_proc,work_per_block2
   integer (int_kind) :: ierr  ! error return code 

   character(len=char_len) :: partitioning_type 
   logical, parameter :: verbose = .FALSE.
!------------------------------------------------------
! Space filling curves only work if:
!
!    nblocks_x = nblocks_y
!       nblocks_x = 2^m 3^n 5^p where m,n,p are integers
!------------------------------------------------------
   
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1: nblocks_x,nblocks_y: ',nblocks_x,nblocks_y
   if((.not. IsFactorable(nblocks_y)) .or. (.not. IsFactorable(nblocks_x))) then
     if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1.1'
     create_distrb_spacecurve = create_distrb_cart(nprocs, work_per_block)
     if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1.2'
     return
   endif
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1.3'

   !-----------------------------------------------
   ! Factor the numbers of blocks in each dimension
   !-----------------------------------------------
   xdim = Factor(nblocks_x)
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1.4'
   ydim = Factor(nblocks_y)
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #1.5'
   numfac = xdim%numfact
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #2'

   !---------------------------------------------
   ! Match the common factors to create SFC curve
   !---------------------------------------------
   curveSize=1
   do it=1,numfac
      call MatchFactor(xdim,ydim,itmp,foundX)
      curveSize = itmp*curveSize
   enddo
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #3'

   !--------------------------------------
   ! determine the size of the sub-blocks 
   ! within the space-filling curve 
   !--------------------------------------
   sb_x = ProdFactor(xdim)
   sb_y = ProdFactor(ydim)

   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #4'
   call create_communicator(dist%communicator, nprocs)

   dist%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #5'
   allocate(blockLocation(nblocks_tot),type_on_curve(nblocks_tot))
   allocate(distance(nblocks_tot))
   allocate(work_on_curve(nblocks_tot))
   allocate (dist%blockLocation(nblocks_tot), &
             dist%blockLocalID (nblocks_tot))

   allocate (dist%blockCnt(nprocs))

   blockLocation = 0
   dist%blockLocation=0
   dist%blockLocalID =0
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #6'

!----------------------------------------------------------------------
!  Create the array to hold the SFC and indices into it
!----------------------------------------------------------------------
   allocate(Mesh(curveSize,curveSize))
   allocate(Mesh2(nblocks_x,nblocks_y))
   allocate(idxT_i(nblocks_tot),idxT_j(nblocks_tot),Lindx(nblocks_tot))
   allocate(Lindx2(nblocks_tot))
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #7'

   Mesh  = 0
   Mesh2 = 0
!----------------------------------------------------------------------
!  Cost function estimation... this is a potential replace for work_per_block 
!----------------------------------------------------------------------
   allocate(Cost_per_proc(nblocks_tot))
   allocate(Cost_per_block(nblocks_tot))
   allocate(work_per_proc(nprocs))
   allocate(work_per_block2(nblocks_tot))
   call EstimateCost(bStats,nblocks_tot,Cost_per_block)
   blockLocation=0

   do i=1,nblocks_tot
      work_per_block2(i) = NINT(10.0*ABS(Cost_per_block(i)),kind=int_kind)
   enddo

   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #8'

!----------------------------------------------------------------------
!  Generate the space-filling curve
!----------------------------------------------------------------------
   call GenSpaceCurve(Mesh)
   Mesh = Mesh + 1    ! make it 1-based indexing
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #8.1'
   if(Debug) then
     if(my_task ==0) call PrintCurve(Mesh)
   endif
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #8.2'
   !-----------------------------------------------
   ! Reindex the SFC to address internal sub-blocks  
   !-----------------------------------------------
   do j=1,curveSize
   do i=1,curveSize
      sfcNum = (Mesh(i,j) - 1)*(sb_x*sb_y) + 1
      do jj=1,sb_y
      do ii=1,sb_x
         subNum = (jj-1)*sb_x + (ii-1)
         i2 = (i-1)*sb_x + ii
         j2 = (j-1)*sb_y + jj
         Mesh2(i2,j2) = sfcNum + subNum
      enddo
      enddo
   enddo
   enddo
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #9'
   !------------------------------------------------
   ! create a linear array of i,j coordinates of SFC
   !------------------------------------------------
   idxT_i=0;idxT_j=0;Lindx=0;Lindx2=0
   do j=1,nblocks_y
     do i=1,nblocks_x
        n = (j-1)*nblocks_x + i
        ig = Mesh2(i,j)
        if(work_per_block(n) /= 0) then
            idxT_i(ig)=i;idxT_j(ig)=j
        endif
        Lindx(n) = ig
     enddo
   enddo
   do i=1,nblocks_tot
      type_on_curve(Lindx(i)) = blockType(i)
   enddo
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #10'

   ! ------------------------------
   ! compress out the land blocks
   ! ------------------------------
   ii=0
   do i=1,nblocks_tot
      if(IdxT_i(i) .gt. 0) then
         ii=ii+1
!         Mesh3(idxT_i(i),idxT_j(i)) = ii
         n = (idxT_j(i)-1)*nblocks_x + idxT_i(i)
         Lindx2(n) = ii 
      endif
   enddo
   nblocks=ii
!DBG   if(my_task == 0) then 
!DBG     write(nu_diag,*) 'work_per_block2:',work_per_block2
!DBG   endif
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #11'
   allocate(cStats(numCoeff,nblocks)) 
   do i=1,nblocks_tot
      if(Lindx2(i)>0) then 
         work_on_curve(Lindx2(i)) = work_per_block2(i) 
         type_on_curve(Lindx2(i))  = blockType(i)
         cStats(:,Lindx2(i)) = bStats(:,i)
         distance(Lindx2(i)) = Lindx(i)
      endif
   enddo
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #12'

   if(lprint_stats) then 
      fname = 'WorkPerBlock.bin'
      call WriteIntegerArray(fname,nblocks,work_on_curve) 
!DBG     write(nu_diag,*) 'work_per_block2:',work_per_block2
!DBG     write(nu_diag,*) 'work_on_curve:', work_on_curve(1:nblocks)
!     open(nu_timing,file='WorkPerBlock.bin',recl=4*nblocks, &
!          form = 'unformatted', access='direct',status='unknown')
!     write(nu_timing,rec=1) work_on_curve(1:nblocks)
!     close(nu_timing)
   endif 

   if(lprint_stats) then 
      call WriteProbabilityStats(bStats,nblocks_tot)
   endif
     
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #13'

   maxB=MIN(max_blocks,maxBlock)
   partitioning_type = 'weight'
   select case(partitioning_type)
      case ('type') 
         ! KLUDGE this is just for testing need to come up with a general solution
         numIce     = COUNT(blockType .eq. iceType)
         minblocks  = CEILING(REAL(numIce)/REAL(nprocs),kind=int_kind)
         !   write(nu_diag,*) 'before TypePartition: {minblocks,maxblocks}: ',minblocks,maxB 
         call  TypePartition(type_on_curve,minblocks,maxB,blockLocation)
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #14'
         if(MAXVAL(blockLocation) > nprocs) then 
              write(nu_diag,*) 'ERROR: problem with partitioning: insufficient processors'
         endif
         ! re-index blockLocation from curve to physical ordering
         do i=1,nblocks_tot
             dist%blockLocation(i) = blockLocation(Lindx(i))
         enddo
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #15'
      case ('weight')
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #16'
         call PartitionCurve(work_on_curve(1:nblocks),work_per_proc, &
		blockLocation(1:nblocks),distance(1:nblocks), nprocs,minBlock, maxB,cStats,FixMaxBlock, maxDil, ierr)
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #17'
!DBG         write(nu_diag,*) 'After PartitionCurve:'
         if(ierr < 0) then 
             call abort_ice('create_distrb_spacecurve: PartitionCurve failed')
         endif 
!DBG         write(nu_diag,*) 'before broadcast_array:'
         call broadcast_array(blockLocation,master_task)
!DBG         write(nu_diag,*) 'After broadcast_array'
         ! re-index blockLocation from curve to physical ordering
         numLocalBlocks=0
         do i=1,nblocks_tot
             if(Lindx2(i)>0) then  
                dist%blockLocation(i) = blockLocation(Lindx2(i))
                numLocalBlocks=numLocalBlocks+1
             endif
         enddo
    end select
!   call qsort(dist%blockLocation(1:numLocalBlocks))
!   write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: dist%blockLocation(:)', dist%blockLocation(1:nblocks_tot)
!   write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: dist%blockLocation(1:numLocalBlocks)', &
!        dist%blockLocation(1:numLocalBlocks)
  
!   call ConvertStatsBlock2Proc(dist%blockLocation,bStats,cStats)
!DBG   write(nu_diag,*) 'Before call to BuildProbabilityStats2'
   call BuildProbabilityStats2(dist%blockLocation,cStats)
!DBG   write(nu_diag,*) 'After call to BuildProbabilityStats2'

   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #19'
   call EstimateCost(cStats,nprocs,Cost_per_proc)
!DBG   write(nu_diag,*) 'before WriteProbabilityStats'
!DBG   write(nu_diag,*) 'after WriteProbabilityStats'

   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #20'
   if(lprint_stats) then 
      fname = 'Q.bin'
      call WriteIntegerArray(fname,nblocks_tot,dist%blockLocation) 

      fname = 'Cost.bin'
      call WriteDblArray(fname,nprocs,Cost_per_proc)

      fname = 'perfm.bin'
      call WriteDblArray(fname,numCoeff,perfmodel)
   endif
   
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #21'

!----------------------------------------------------------------------
!  Reset the dist data structure
!----------------------------------------------------------------------
   allocate(proc_tmp(nprocs))
   proc_tmp = 0

   allocate(dist%blockIndex(nprocs,max_blocks))
   dist%blockIndex(:,:) = 0
   do n=1,nblocks_tot
      pid = dist%blockLocation(n)
      if(pid>0) then
        proc_tmp(pid) = proc_tmp(pid) + 1
        dist%blockLocalID(n) = proc_tmp(pid)
        if(proc_tmp(pid) <= max_blocks) then 
            dist%blockIndex(pid,proc_tmp(pid)) = n
        endif
      endif
   enddo
   dist%blockCnt(:) = proc_tmp(:) 

   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #22'
   !---------------------------------------
   ! Set the number of active local blocks
   !---------------------------------------
   dist%numLocalBlocks = proc_tmp(my_task+1)
   if (dist%numLocalBlocks>0) then 
      allocate(dist%blockGlobalID(dist%numLocalBlocks))
      localID=1
      do n=1,nblocks_tot
         pid = dist%blockLocation(n)
         if(pid == my_task+1) then 
              dist%blockGlobalID(localID) = n
              localID=localID+1
         endif
      enddo
   endif
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #23'
   if(Debug) then
      if(my_task==0) write(nu_diag,*) 'dist%blockLocation:= ',dist%blockLocation
      write(nu_diag,*) 'IAM: ',my_task,' SpaceCurve: Number of blocks {total,local} :=', &
                nblocks_tot,nblocks,proc_tmp(my_task+1)
   endif
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate block'  
   !---------------------------------
   ! Deallocate temporary arrays
   !---------------------------------
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(blockLocation)'  
   deallocate(blockLocation)
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(type_on_curve)'  
   deallocate(type_on_curve)
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(work_on_curve)'  
   deallocate(work_on_curve)
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(proc_tmp)'  
   deallocate(proc_tmp)
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(Mesh,Mesh2)'  
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #24'
   deallocate(Mesh,Mesh2)
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(idxT_i,idxT_j,Lindx,Lindx2)'  
   deallocate(idxT_i,idxT_j,Lindx,Lindx2)
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(Cost_per_proc)'  
   deallocate(Cost_per_proc) 
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(Cost_per_block)'  
   deallocate(Cost_per_block)
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(work_per_proc)'  
   deallocate(work_per_proc)
!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before deallocate(work_per_block2)'  
   deallocate(work_per_block2)  
    
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #25'
   if(verbose .and. my_task == 0) then 
      write(nu_diag,*) 'create_distrb_spacecurve: blockCnt ',dist%blockCnt
      write(nu_diag,*) 'create_distrb_spacecurve: blockIndex ',dist%blockIndex
   endif

!DBG   write(nu_diag,*) 'create_distrb_spacecurve: before assignment of result'  
!----------------------------------------------------------------------
   create_distrb_spacecurve = dist  ! return the result
!----------------------------------------------------------------------

   if (verbose .and. my_task==0) then 
      write(nu_diag,*) 'create_distrb_spacecurve%blockGlobalID:  ',create_distrb_spacecurve%blockGlobalID
      write(nu_diag,*) 'create_distrb_spacecurve%blockCnt: ',create_distrb_spacecurve%blockCnt
   endif
!DBG   write(nu_diag,*) 'At the end of create_distrb_spacecurve'  
   if(verbose) write(nu_diag,*) 'IAM: ',my_task,'create_distrb_spacecurve: point #26'
!EOC

 end function create_distrb_spacecurve

 subroutine TypePartition(blockType,minblocks,maxblocks,blockLocation)

   integer(kind=int_kind), intent(in)     :: blockType(nblocks_tot)
   integer, intent(in)  :: minblocks, & ! Minimum number of blocks per processor 
                           maxblocks    ! Maximum number of blocks per processor 
   integer(kind=int_kind), intent(inout)  :: blockLocation(nblocks_tot)

   integer :: i,ip,cur,next 
   integer :: cntLnd,tcntLnd, &
              cntIcefree,tcntIcefree, &
              cntIce,tcntIce


   tcntIce     = 0
   tcntLnd     = 0
   tcntIcefree = 0
   cntIce      = 0
   cntLnd      = 0
   cntIcefree  = 0
   ip = 1
   do i=1,nblocks_tot
      cur  = blockType(i)
      if(i<nblocks_tot) then 
         next = blockType(i+1)
      else
         next = cur
      endif

      !------------ 
      ! Land point 
      !------------ 
      if(cur == lndType) then 
          tcntLnd = tcntLnd + 1
      endif

     !----------------------
     ! ice Free  point
     !----------------------
     if (cur == icefreeType) then 
        tcntIcefree = tcntIcefree + 1
        cntIcefree = cntIcefree + 1
        blockLocation(i) = ip
        if ((cntIcefree == maxblocks) .or. (next /= cur)) then
           ip = ip + 1
           cntIcefree=0
        endif
     endif


     !----------------------
     ! ice Free  point
     !----------------------
     if (cur == iceType) then 
        tcntIce = tcntIce + 1
        cntIce  = cntIce + 1
        blockLocation(i) = ip
        if ((cntIce == minblocks) .or. (next /= cur)) then
           ip = ip + 1
           cntIce=0
        endif
     endif
   enddo

   if(my_task == 0) then 
      write(*,23) tcntLnd+tcntIce+tcntIcefree, tcntIce, tcntIcefree, tcntLnd
!      write(nu_diag,*) 'TypePartition: land blks: ',tcntLnd,' Ice blks: ', &
!                tcntIce,' IceFree blks: ',tcntIcefree
      write(*,24) MAXVAL(blockLocation) 
!      write(nu_diag,*) 'TypePartition: Partitioned across ',MAXVAL(blockLocation),' processors'
   endif
  
23   format('Total blocks: ',i5,' Ice blocks: ',i5,' IceFree blocks: ',i5,' Land blocks: ',i5)
24   format('Partitioned across ',i5,' processors')
  

 end subroutine TypePartition

  subroutine PartitionCurve(work_per_block, work_per_proc, blockLocation, distance, &
             nproc, min_blocks, max_blocks, Stats, FixMaxBlock, maxDil, ierr)

    integer (int_kind), intent(inout) :: work_per_block(:)
    integer (int_kind), intent(inout) :: work_per_proc(:)
    integer (int_kind), intent(inout) :: blockLocation(:)
    integer (int_kind), intent(inout) :: distance(:)
    integer (int_kind), intent(in)    :: nproc
    integer (int_kind), intent(in)    :: min_blocks
    integer (int_kind), intent(in)    :: max_blocks
    real    (dbl_kind), intent(in)    :: Stats(:,:) 
    logical,            intent(in)    :: FixMaxBlock
    real    (real_kind), intent(in)   :: maxDil
    integer (int_kind), intent(inout) :: ierr

    integer :: cnt

    integer :: nb,anProc
    integer :: n,ip,i,imax
    real :: totalCost, avgCost, maxCost,dtCost
    integer, allocatable :: saveblockLocation(:)
    integer :: maxBlocks,minBlocks
    integer :: ivalue
    real (real_kind) :: maxValue
    real :: maxValue_save
    real :: minValue
    logical :: contLoop

    integer :: maxB,minB

    logical, parameter :: verbose = .TRUE.
    integer :: it,maxiter
    integer :: save_maxB 
    real :: save_maxCost
    
    integer :: aminBlocks, amaxBlocks 

    real :: aWork,minWork,maxWork,maxWorkBlock
    real :: minCostBlock, maxCostBlock
    real :: maxCost_old
    real (real_kind) :: amaxDil
    real (dbl_kind), allocatable, dimension(:) :: cost_per_proc
    real (dbl_kind), allocatable, dimension(:) :: cost_per_block
    real (dbl_kind), allocatable, dimension(:,:) :: pStats 


    allocate(pStats(numCoeff,nblocks_tot))

    maxB=max_blocks
!    maxDil = 5.0

    if ( my_task .eq. master_task) then  

       nb = SIZE(blockLocation)  ! number of blocks
       allocate(saveblockLocation(nb))
       allocate(cost_per_proc(nb))
       allocate(cost_per_block(nb))

       ! --------------------------------------------
       ! Estimate the computational cost of each block  
       ! --------------------------------------------
       call EstimateCost(Stats,nb,cost_per_block)

       totalCost = SUM(cost_per_block)
       save_maxCost = totalCost
       maxCostBlock = MAXVAL(cost_per_block)
       minCostBlock = MINVAL(cost_per_block)
       save_maxB = maxB
       avgCost = (totalCost/nproc)

       write(nu_diag,*) 'PartitionCurve: nblocks,nproc ',nb,nproc
       write(nu_diag,213) totalCost, avgCost, minCostBlock, maxCostBlock
!DBG       write(nu_diag,*) 'PartitionCurve: totalCost,avgCost, maxCostBlock: ',totalCost,avgCost,maxCostBlock
!DBG       write(nu_diag,*) distance


       minB = CEILING(real(nb)/real(nproc),kind=int_kind)
       if(maxB < minB ) then
          write(nu_diag,*) 'ERROR: unable to partition ',nb,' blocks across ',nproc,' processors'
          write(nu_diag,*) 'ERROR: Either increase max_blocks := ',maxB
          write(nu_diag,*) 'ERROR: Either increase number of processors'
          ierr = -2
          return
       endif
       if(FixMaxBlock) then 
         minB = maxB
       endif
       do while (maxB  >= minB )  

          dtCost = 1.0
          maxValue = maxCostBlock*real(maxB)
          minValue = maxCostBlock
          maxCost_old = maxValue

          contLoop = .true.
          maxiter = 20
          it = 1
          do while(contLoop )

            cost_per_proc=0.0
            call wPartition(cost_per_block,blockLocation, distance, nproc,min_blocks, maxB,maxValue,maxDil,aminBlocks, &
                 amaxBlocks, amaxDil)
            anProc = MAXVAL(blockLocation)
            call ConvertStatsBlock2Proc(blockLocation,Stats,pStats)
            call EstimateCost(pStats,anProc,cost_per_proc)
            maxCost = MAXVAL(cost_per_proc)
   

            if(lprint_stats) then 
		write(nu_diag,211) it,anProc,aminBlocks,amaxBlocks,maxB,minValue, maxValue,maxCost, amaxDil
            endif

            if(maxCost > maxValue) then
               minValue =  maxValue
               dtCost = (maxCost_old-minValue)/2.0
               maxValue = maxCost_old - dtCost
            else
               dtCost = (maxCost-minValue)/2.0
               maxValue = maxCost - dtCost
               maxCost_old = maxCost
            endif

            if(maxCost == maxCostBlock) contLoop = .false.
            if(dtCost < 1.0e-5) contLoop = .false.
            if( anProc == nproc .and. it >= maxiter)  contLoop = .false.
            if ((save_maxCost > maxCost .and. amaxBlocks <= maxB) .or. &
                (save_maxCost == maxCost .and. save_maxB > maxB)  ) then
                save_maxCost = maxCost
                save_maxB = maxB  
                saveblockLocation = blockLocation
            endif
            it=it+1
          enddo
          maxB = maxB - 1
       enddo
       blockLocation = saveblockLocation
       write(nu_diag,214) perfmodel_name

       write(nu_diag,*) '-------------------------wSFC-----------------------'
       call PrintPartitionLB(blockLocation,nproc,Stats) 

 211 format('Partition loop: it: ',i4,' anProc: ',i4,' a{min,max}Blocks ', &
              2(i4),' maxB: ',i4,' [min,max]Value: ',2(f14.4),' maxCost: ',f14.4,' maxDilation: ',f14.4 )  
 213 format('PartitionCurve: TotalCost: ',f12.4,' avgCost: ',f12.4, &
	    ' minBlockCost: ',f12.4,' maxBlockCost: ',f12.4 ) 
 214 format('Using performance model: ',a20)

       deallocate(cost_per_proc)
       deallocate(cost_per_block)
       deallocate(saveblockLocation)
    endif

    deallocate(pStats)

    ierr = 0

  end subroutine PartitionCurve

  subroutine wPartition(cost_per_block, blockLocation, distance, nproc, min_blocks, max_blocks, maxValue, maxDil, aminBlocks, &
       amaxBlocks,amaxDil)

    real (dbl_kind), intent(in) :: cost_per_block(:)
    integer (int_kind), intent(inout) :: blockLocation(:)
    integer (int_kind), intent(inout) :: distance(:)
    integer (int_kind), intent(in) :: nproc
    integer (int_kind), intent(in) :: min_blocks
    integer (int_kind), intent(in) :: max_blocks
    real (real_kind), intent(in) :: maxvalue
    real (real_kind), intent(in)  :: maxDil    ! maximum alloable dilation of domains
    integer (int_kind), intent(inout) :: aminBlocks, amaxBlocks 
    real (real_kind), intent(out) :: amaxDil

    integer (int_kind)  :: n,ip,i,numB
    real :: totalCost,avgCost,sumTMP,sumTMP2
    
!   integer (int_kind), allocatable :: work_per_proc(:)
    
    integer :: minDist,Dist
    integer :: sloc
    real    :: dilation, dilation2

    logical, parameter :: Info = .FALSE.
    logical  :: break_loop

    n = SIZE(blockLocation)  ! number of blocks

    totalCost = SUM(cost_per_block)

    aminBlocks = n
    amaxBlocks = 0
    amaxDil  = 1.0
    avgCost = (totalCost/nproc)

!DBG    write(nu_diag,*) 'cost_per_block: ',cost_per_block

    ip = 1
    i=1
!    cost_per_block=0
    do while (i<=n)
        sumTMP = 0
        numB = 0
        break_loop = .FALSE.
        sloc = distance(i)
        do while ( ((sumTMP<maxValue) .or. (ip ==  nproc))  &
		.and. (i<=n) &
		.and. (.not. break_loop))
          sumTMP2  = sumTMP + cost_per_block(i)
          minDist  = numB + 1
          Dist     = distance(i) - sloc + 1
          dilation2 = real(Dist)/real(minDist)
          if(((sumTMP2 <= maxValue) .and. (numB < max_blocks) .and. dilation2 <= maxDil) &
!          if(((sumTMP2 <= maxValue) .and. (numB < max_blocks)) &
                 .or. (ip == nproc) .or. (numB < min_blocks) ) then
              blockLocation(i) = ip
              i=i+1
              sumTMP = sumTMP2
              dilation = dilation2
              numB = numB+1
          else
              break_loop = .TRUE.
          endif
        enddo
        if(aminBlocks > numB)  aminBlocks = numB
        if(amaxBlocks < numB)  amaxBlocks = numB
        if(amaxDil < dilation) amaxDil = dilation
        ip = ip+1
    enddo
    
!    call abort_ice('wPartition: at the end of the subroutine')

  end subroutine wPartition

!**********************************************************************
!BOP
! !IROUTINE: ice_distributionRake
! !INTERFACE:

 subroutine ice_distributionRake (procWork, procID, blockWork, &
                                  priority, distribution)

! !DESCRIPTION:
!  This subroutine performs a rake algorithm to distribute the work
!  along a vector of processors.  In the rake algorithm, a work
!  threshold is first set.  Then, moving from left to right, work
!  above that threshold is raked to the next processor in line.
!  The process continues until the end of the vector is reached
!  and then the threshold is reduced by one for a second rake pass.
!  In this implementation, a priority for moving blocks is defined
!  such that the rake algorithm chooses the highest priority
!  block to be moved to the next processor.  This can be used
!  for example to always choose the eastern-most block or to
!  ensure a block does not stray too far from its neighbors.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in), dimension(:) :: &
      blockWork          ,&! amount of work per block
      procID               ! global processor number

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), intent(inout), dimension(:) :: &
      procWork           ,&! amount of work per processor
      priority             ! priority for moving a given block

   type (distrb), intent(inout) :: &
      distribution         ! distribution to change

! !OUTPUT PARAMETERS:

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, n,                  &! dummy loop indices
      np1,                   &! n+1 corrected for cyclical wrap
      iproc, inext,          &! processor ids for current and next 
      nprocs, numBlocks,   &! number of blocks, processors
      lastPriority,          &! priority for most recent block
      minPriority,           &! minimum priority
      lastLoc,               &! location for most recent block
      meanWork, maxWork,     &! mean,max work per processor
      diffWork, residual,    &! work differences and residual work
      numTransfers            ! counter for number of block transfers

!----------------------------------------------------------------------
!
!  initialization
!
!----------------------------------------------------------------------

   nprocs  = size(procWork)
   numBlocks = size(blockWork)

   !*** compute mean,max work per processor

   meanWork = sum(procWork)/nprocs + 1
   maxWork  = maxval(procWork)
   residual = mod(meanWork,nprocs)

   minPriority = 1000000
   do n=1,nprocs
      iproc = procID(n)
      do i=1,numBlocks
         if (distribution%blockLocation(i) == iproc) then
            minPriority = min(minPriority,priority(i))
         endif
      end do
   end do

!----------------------------------------------------------------------
!
!  do two sets of transfers
!
!----------------------------------------------------------------------

   transferLoop: do

!----------------------------------------------------------------------
!
!     do rake across the processors
!
!----------------------------------------------------------------------

      numTransfers = 0
      do n=1,nprocs
         if (n < nprocs) then
            np1   = n+1
         else
            np1   = 1
         endif
         iproc = procID(n)
         inext = procID(np1)

         if (procWork(n) > meanWork) then !*** pass work to next

            diffWork = procWork(n) - meanWork

            rake1: do while (diffWork > 1)

               !*** attempt to find a block with the required
               !*** amount of work and with the highest priority
               !*** for moving (eg boundary blocks first)

               lastPriority = 0
               lastLoc = 0

               do i=1,numBlocks
                  if (distribution%blockLocation(i) == iproc) then
                     if (priority(i) > lastPriority ) then
                        lastPriority = priority(i)
                        lastLoc = i
                     endif
                  endif
               end do
               if (lastLoc == 0) exit rake1 ! could not shift work

               numTransfers = numTransfers + 1
               distribution%blockLocation(lastLoc) = inext
               if (np1 == 1) priority(lastLoc) = minPriority
               diffWork = diffWork - blockWork(lastLoc)

               procWork(n  ) = procWork(n  )-blockWork(lastLoc)
               procWork(np1) = procWork(np1)+blockWork(lastLoc)
            end do rake1
         endif

      end do

!----------------------------------------------------------------------
!
!     increment meanWork by one and repeat
!
!----------------------------------------------------------------------

      meanWork = meanWork + 1
      if (numTransfers == 0 .or. meanWork > maxWork) exit transferLoop

   end do transferLoop

!----------------------------------------------------------------------
!EOC

end subroutine ice_distributionRake

!***********************************************************************
    subroutine PrintPartitionLB(Location,n,bStats)

      integer :: Location(:)
      integer :: n
      real (dbl_kind), dimension(:,:) :: bStats 

      real :: maxCost,minCost
      integer :: anProc
      integer :: minB, maxB 
      real :: aCost
      real (dbl_kind), allocatable, dimension(:) :: cost_per_proc(:)
      integer :: i,ncnt

      real (dbl_kind), allocatable, dimension(:,:) :: pStats

      allocate(cost_per_proc(n))
      allocate(pStats(numCoeff,n))
      
       
      call ConvertStatsBlock2Proc(Location,bStats,pStats)
      call EstimateCost(pStats,n,cost_per_proc)

!DBG      write(nu_diag,*) 'PrintPartitinoLB: Location:',Location
!DBG      write(nu_diag,*) 'PrintPartitinoLB: cost_per_proc:',cost_per_proc

      maxCost = MAXVAL(cost_per_proc)
      minCost = MINVAL(cost_per_proc)
      aCost = SUM(cost_per_proc)/real(n)
      anProc=MAXVAL(Location)
      maxB=0
      minB=n
      do i=1,n
           ncnt = COUNT(Location==i)  
           maxB = MAX(ncnt,maxB)
           minB = MIN(ncnt,minB)
      enddo

!DBG      write(nu_diag,*) 'maxCost: ',maxCost
!DBG      write(nu_diag,*) 'minCost: ',minCost
!DBG      write(nu_diag,*) 'aCost: ',aCost

!      write(nu_diag,*) 'PrintPartitionLB: on ',anProc,' processors Avg,Min,Max work/proc, imbalance  ', &
!                aWork,minWork,maxWork,ABS(aWork-maxWork)/aWork
      write(nu_diag,212)  anProc, minB, maxB,aCost,minCost,maxCost,ABS(aCost-maxCost)/aCost
      deallocate(cost_per_proc)
      deallocate(pStats)

 212 format('PrintPartitionLB: on ',i4,' procs a{min,max}Blocks: ',(2(i4)),' Avg,Min,Max Cost/proc: ',(3(f10.4)), &
          ' imbalance: ',f8.2)

    end subroutine PrintPartitionLB

   subroutine EstimateCost(coeffMatrix,n,Cost)

     real (dbl_kind) :: coeffMatrix(:,:)
     integer (int_kind) :: n
     real (dbl_kind):: Cost(:)

     real (dbl_kind) :: tmp

     integer (int_kind) :: i,j

     Cost=0.0_dbl_kind
     do i=1,n
        tmp = 0.0_dbl_kind
        do j=1,numCoeff
           tmp = tmp + coeffMatrix(j,i) *perfmodel(j)
        enddo
        Cost(i) = tmp
     enddo

   end subroutine EstimateCost

   subroutine ConvertStatsBlock2Proc(Location,bStats,pStats)

      integer (int_kind) :: Location(:)
      real (dbl_kind), intent(in)  :: bStats(:,:)
      real (dbl_kind), intent(out) :: pStats(:,:)

      integer (int_kind) :: i,ip,n

      n = size(Location)
      pStats = 0.0d0
      do i=1,n
         ip = Location(i)
         if(ip > 0) then 
            pStats(:,ip) = pStats(:,ip) + bStats(:,i) 
         endif
      enddo

   end subroutine ConvertStatsBlock2Proc

   subroutine WriteProbabilityStats(coeffMatrix,n)

    real(dbl_kind)  :: coeffMatrix(:,:)
    integer (int_kind) :: n


    if(my_task == master_task) then
       open(nu_timing,file='probStats.bin',recl=8*numCoeff*n, &
            form = 'unformatted', access = 'direct', status = 'unknown')
       write(nu_timing,rec=1) coeffMatrix(:,1:n)
       close(nu_timing)
    endif


   end subroutine WriteProbabilityStats

   subroutine WriteIntegerArray(fname,n,array)
     character(char_len) :: fname
     integer (int_kind) :: n
     integer (int_kind) :: array(:)

     if(my_task == master_task) then 
        open(nu_timing,file=TRIM(fname),recl=4*n, &
          form = 'unformatted', access = 'direct', status = 'unknown')
        write(nu_timing,rec=1) array(1:n)
        close(nu_timing)
     endif

   end subroutine WriteIntegerArray

   subroutine WriteDblArray(fname,n,array)
     character(char_len) :: fname
     real (dbl_kind) :: array(:)
     integer (int_kind) :: n

     if(my_task == master_task) then 
        open(nu_timing,file=TRIM(fname),recl=8*n, &
          form = 'unformatted', access = 'direct', status = 'unknown')
        write(nu_timing,rec=1) array(1:n)
        close(nu_timing)
     endif

   end subroutine WriteDblArray

end module ice_distribution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
