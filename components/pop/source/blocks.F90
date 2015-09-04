!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module blocks

!BOP
! !MODULE: blocks
!
! !DESCRIPTION: 
!  This module contains data types and tools for decomposing a global
!  horizontal domain into a set of blocks.  It contains a data type 
!  for describing each block and contains routines for creating and 
!  querying the block decomposition for a global domain.
!
! !REVISION HISTORY:
!  SVN:$Id: blocks.F90 44198 2013-02-25 22:43:22Z mlevy@ucar.edu $
!
! !USES:

   use kinds_mod
   use exit_mod
   use domain_size, only: nx_global, ny_global, block_size_x, block_size_y
   use communicate, only: my_task

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: block   ! block data type
      integer (int_kind) :: &
         block_id          ,&! global block number
         local_id          ,&! local address of block in current distrib
         ib, ie, jb, je    ,&! begin,end indices for physical domain
         iblock, jblock      ! cartesian i,j position for bloc

      integer (int_kind), dimension(:), pointer :: &
         i_glob, j_glob     ! global domain location for each point
   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: create_blocks       ,&
             destroy_blocks      ,&
             get_block           ,&
             get_block_parameter ,&
             get_block_ids_from_coords

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public :: &
      nghost = 2       ! number of ghost cells around each block

   integer (int_kind), parameter, public :: &! size of block domain in
      nx_block = block_size_x + 2*nghost,   &!  x,y dir including ghost
      ny_block = block_size_y + 2*nghost     !  cells 

! !PUBLIC DATA MEMBERS:

   integer (int_kind), public :: &
      nblocks_tot      ,&! total number of blocks in decomposition
      nblocks_x        ,&! tot num blocks in i direction
      nblocks_y          ! tot num blocks in j direction

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module private data
!
!-----------------------------------------------------------------------

   type (block), dimension(:), allocatable :: &
      all_blocks         ! block information for all blocks in domain

   integer (int_kind), dimension(:,:), allocatable, target :: &
      i_global,         &! global i index for each point in each block
      j_global           ! global j index for each point in each block

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: create_blocks
! !INTERFACE:

 subroutine create_blocks(nx_global, ny_global, ew_boundary_type, &
                                                ns_boundary_type)

! !DESCRIPTION:
!  This subroutine decomposes the global domain into blocks and
!  fills the data structures with all the necessary block information.
!
! !REVISION HISTORY: 
!  same as module!
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nx_global, ny_global           ! global domain size in x,y

   character (*), intent(in) :: &
      ew_boundary_type,  &! type of boundary in logical east-west dir
      ns_boundary_type    ! type of boundary in logical north-south dir

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, ip1, j, jp1, n    ,&! loop indices
      iblock, jblock       ,&! block loop indices
      is, ie, js, je         ! temp start, end indices

!----------------------------------------------------------------------
!
!  compute number of blocks and cartesian decomposition
!  if the requested block size does not divide the global domain
!  size evenly, add additional block space to accomodate padding
!
!----------------------------------------------------------------------

   nblocks_x   = (nx_global-1)/block_size_x + 1
   nblocks_y   = (ny_global-1)/block_size_y + 1
   nblocks_tot = nblocks_x*nblocks_y

!----------------------------------------------------------------------
!
!  allocate block arrays
!
!----------------------------------------------------------------------

   allocate(all_blocks(nblocks_tot))
   allocate(i_global(nx_block,nblocks_tot), &
            j_global(ny_block,nblocks_tot))

!----------------------------------------------------------------------
!
!  fill block data structures for all blocks in domain
!
!----------------------------------------------------------------------

   n = 0
   do jblock=1,nblocks_y
      js = (jblock-1)*block_size_y + 1
      je = js + block_size_y - 1
      if (js > ny_global) call exit_POP(sigAbort, &
               'create_blocks: Bad block decomp: ny_block too large?')
      if (je > ny_global) je = ny_global ! pad array

      do iblock=1,nblocks_x
         n = n + 1  ! global block id

         is = (iblock-1)*block_size_x + 1
         ie = is + block_size_x - 1
         if (is > nx_global) call exit_POP(sigAbort, &
               'create_blocks: Bad block decomp: nx_block too large?')
         if (ie > nx_global) ie = nx_global

         all_blocks(n)%block_id = n
         all_blocks(n)%iblock   = iblock
         all_blocks(n)%jblock   = jblock
         all_blocks(n)%ib       = nghost + 1
         all_blocks(n)%jb       = nghost + 1
         all_blocks(n)%ie       = nx_block - nghost ! default value
         all_blocks(n)%je       = ny_block - nghost ! default value

         do j=1,ny_block
            j_global(j,n) = js - nghost + j - 1


            !*** southern ghost cells

            if (j_global(j,n) < 1) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) + ny_global
               case ('closed')
                  j_global(j,n) = 0
               case ('tripole')
                  j_global(j,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown n-s bndy type')
               end select
            endif

            !*** padding required

            if (j_global(j,n) > ny_global + nghost) then
               j_global(j,n) = 0   ! padding

            !*** northern ghost cells

            else if (j_global(j,n) > ny_global) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) - ny_global
               case ('closed')
                  j_global(j,n) = 0
               case ('tripole')
                  j_global(j,n) = -j_global(j,n)
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown n-s bndy type')
               end select

            !*** set last physical point if padded domain

            else if (j_global(j,n) == ny_global .and. &
                     j > all_blocks(n)%jb) then
               all_blocks(n)%je = j   ! last physical point in padded domain
            endif
         end do

         all_blocks(n)%j_glob => j_global(:,n)

         do i=1,nx_block
            i_global(i,n) = is - nghost + i - 1

            !*** western ghost cells

            if (i_global(i,n) < 1) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global(i,n) = i_global(i,n) + nx_global
               case ('closed')
                  i_global(i,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown e-w bndy type')
               end select
            endif

            !*** padded domain - fill padded region with zero

            if (i_global(i,n) > nx_global + nghost) then
               i_global(i,n) = 0

            !*** eastern ghost cells

            else if (i_global(i,n) > nx_global) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global(i,n) = i_global(i,n) - nx_global
               case ('closed')
                  i_global(i,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown e-w bndy type')
               end select

            !*** last physical point in padded domain

            else if (i_global(i,n) == nx_global .and. &
                     i > all_blocks(n)%ib) then
               all_blocks(n)%ie = i
            endif
         end do

         all_blocks(n)%i_glob => i_global(:,n)

      end do
   end do

!EOC
!----------------------------------------------------------------------

end subroutine create_blocks

!***********************************************************************
!BOP
! !IROUTINE: get_block
! !INTERFACE:

 function get_block(block_id,local_id)

! !DESCRIPTION:
!  This function returns the block data structure for the block
!  associated with the input block id.
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      block_id,   &! global block id for requested block info
      local_id     ! local  block id to assign to this block

! !OUTPUT PARAMETERS:

   type (block) :: &
      get_block    ! block information returned for requested block

!EOP
!BOC
!----------------------------------------------------------------------
!
!  check for valid id.  if valid, return block info for requested block
!
!----------------------------------------------------------------------

   if (block_id < 1 .or. block_id > nblocks_tot) then
      call exit_POP(sigAbort,'get_block: invalid block_id')
   endif

   get_block = all_blocks(block_id)
   get_block%local_id = local_id

!----------------------------------------------------------------------
!EOC

 end function get_block

!**********************************************************************
!BOP
! !IROUTINE: get_block_parameter
! !INTERFACE:

 subroutine get_block_parameter(block_id, local_id, ib, ie, jb, je, &
                                iblock, jblock, i_glob, j_glob)

! !DESCRIPTION:
!  This routine returns requested parts of the block data type
!  for the block associated with the input block id
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      block_id   ! global block id for which parameters are requested

! !OUTPUT PARAMETERS:

   !(optional) parts of block data type to extract if requested

   integer (int_kind), intent(out), optional :: &
      local_id         ,&! local id assigned to block in current distrb
      ib, ie, jb, je   ,&! begin,end indices for physical domain
      iblock, jblock     ! cartesian i,j position for bloc

   integer (int_kind), dimension(:), pointer, optional :: &
      i_glob, j_glob     ! global domain location for each point

!EOP
!BOC
!----------------------------------------------------------------------
!
!  extract each component of data type if requested
!
!----------------------------------------------------------------------

   if (block_id < 1 .or. block_id > nblocks_tot) then
      call exit_POP(sigAbort,'get_block_parameter: invalid block_id')
   endif

   if (present(local_id)) local_id = all_blocks(block_id)%local_id
   if (present(ib      )) ib       = all_blocks(block_id)%ib
   if (present(ie      )) ie       = all_blocks(block_id)%ie
   if (present(jb      )) jb       = all_blocks(block_id)%jb
   if (present(je      )) je       = all_blocks(block_id)%je
   if (present(iblock  )) iblock   = all_blocks(block_id)%iblock
   if (present(jblock  )) jblock   = all_blocks(block_id)%jblock
   if (present(i_glob  )) i_glob   = all_blocks(block_id)%i_glob
   if (present(j_glob  )) j_glob   = all_blocks(block_id)%j_glob

!----------------------------------------------------------------------
!EOC

 end subroutine get_block_parameter

!**********************************************************************
!BOP
! !IROUTINE: destroy_blocks
! !INTERFACE:

 subroutine destroy_blocks

! !DESCRIPTION:
!  This subroutine deallocates the array with block information.
!
! !REVISION HISTORY:
!  same as module
!EOP
!----------------------------------------------------------------------
!BOC

   deallocate(all_blocks)

!EOC
!----------------------------------------------------------------------

 end subroutine destroy_blocks


!***********************************************************************
!BOP
! !IROUTINE: get_block_ids_from_coords
! !INTERFACE:

 subroutine get_block_ids_from_coords(num_blocks, block_ids, &
      g_imin, g_jmin, g_imax, g_jmax)

! !DESCRIPTION:
!  This function returns the global block id(s) of the block(s) 
!  that contain the given global coords (if min and max are given,
!  then we find all ids within the box defined by (g_imin, g_jmin) and
!   (g_imax, g_jmax
!
!
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      g_imin, g_jmin !global coordinates in i and j

   integer (int_kind), intent(in), optional :: &
      g_imax, g_jmax !global coordinates in i and j


! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      num_blocks    ! number of block ids returned
   
   integer (int_kind), intent(out), dimension(:), pointer :: &
        block_ids


! ! local variables 

   integer (int_kind) :: imin_pos, &
        jmin_pos, imax_pos, jmax_pos, id, ni, nj, i, j,count
   real (r4) :: tmp_r

   !make sure these are valid coords, i.e. in the global domain

   if ( (nx_global < g_imin)  .or. &
        (ny_global < g_jmin)) then
      call exit_POP(sigAbort, &
           'get_block_from_coords: invalid coords - too large')
   end if

   if (present(g_imax)) then
      if (present(g_jmax)) then
         if ( (nx_global < g_imax)  .or. &
              (ny_global < g_jmax)) then
            call exit_POP(sigAbort, &
                 'get_block_ids_from_coords: invalid coords - too large')
         end if
      else
         call exit_POP(sigAbort, &
              'get_block_ids_fromPcoords:  must supply BOTH g_imax AND g_jmax.')
      end if
   end if


   !lets's find out what the block position containing the min coords are is

   imin_pos = ceiling(real(g_imin, r4)/real(block_size_x, r4))
   jmin_pos = ceiling(real(g_jmin, r4)/real(block_size_y, r4))
   

   !now, which global block id is this min?
   id = (jmin_pos-1)*nblocks_x + imin_pos

   ! set return 

   if (.not. present(g_imax)) then
      !just one to return
      num_blocks = 1
      allocate(block_ids(num_blocks))
      block_ids(1) = id
   else
      !need to figure out all the block ids in the described box

      !get the block position for the max coords
      imax_pos = ceiling(real(g_imax, r4)/real(block_size_x, r4))
      jmax_pos = ceiling(real(g_jmax, r4)/real(block_size_y, r4))

      !how many blocks?
      ni = imax_pos - imin_pos + 1
      nj = jmax_pos - jmin_pos + 1

      num_blocks = ni*nj

      allocate(block_ids(num_blocks))
      count = 0

      do j = jmin_pos, jmax_pos
         do i = imin_pos, imax_pos            
            id = (j-1)*nblocks_x + i
            count =  count + 1
            block_ids(count) = id
         end do
      end do

   end if

!   print *, 'IAM (in blocs.F90): ', my_task, 'num blocks: ', num_blocks, &
!        'blocks= ', block_ids(1:num_blocks)



 end subroutine get_block_ids_from_coords


!***********************************************************************

 end module blocks

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
