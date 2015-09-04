!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ice_blocks

 module ice_blocks

!
! !DESCRIPTION: 
!  This module contains data types and tools for decomposing a global
!  horizontal domain into a set of blocks.  It contains a data type 
!  for describing each block and contains routines for creating and 
!  querying the block decomposition for a global domain.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_blocks.F90 100 2008-01-29 00:25:32Z eclare $
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP by William H. Lipscomb, LANL
!
! !USES:

   use ice_kinds_mod
   use ice_exit
   use ice_domain_size

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: block   ! block data type
      integer (int_kind) :: &
         block_id           ,&! global block number
         local_id           ,&! local address of block in current distrib
         ilo, ihi, jlo, jhi ,&! begin, end indices for physical domain
         iblock, jblock       ! cartesian i,j position for block

      logical (log_kind) :: &
         tripole,           & ! flag is true if block is at tripole bndy
         tripoleTFlag         ! tripole boundary is a T-fold

      integer (int_kind), dimension(:), pointer :: &
         i_glob, j_glob     ! global domain location for each point
   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: create_blocks       ,&
             get_block           ,&
             get_block_parameter ,&
             ice_blocksGetNbrID

   public :: LocateTrouble,PrintTrouble
   interface PrintTrouble
       module procedure PrintTrouble_dbl, PrintTrouble_int
   end interface

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public :: &
      nghost = 1       ! number of ghost cells around each block

   integer (int_kind), parameter, public :: &! size of block domain in
      nx_block = block_size_x + 2*nghost,   &!  x,y dir including ghost
      ny_block = block_size_y + 2*nghost     !  cells 

   ! predefined directions for neighbor id routine
   ! Note: the directions that are commented out are implemented in 
   !       POP but not in CICE.  If the tripole cut were in the south
   !       instead of the north, these would need to be used (and also
   !       implemented in ice_boundary.F90).
   integer (int_kind), parameter, public :: &
      ice_blocksNorth          =  1,      & ! (i  ,j+1)
      ice_blocksSouth          =  2,      & ! (i  ,j-1)
      ice_blocksEast           =  3,      & ! (i+1,j  )
      ice_blocksWest           =  4,      & ! (i-1,j  )
      ice_blocksNorthEast      =  5,      & ! (i+1,j+1)
      ice_blocksNorthWest      =  6,      & ! (i-1,j+1)
      ice_blocksSouthEast      =  7,      & ! (i+1,j-1)
      ice_blocksSouthWest      =  8         ! (i-1,j-1)
   integer (int_kind), parameter, public :: &
!      ice_blocksNorth2         =  9,      & ! (i  ,j+2)
!      ice_blocksSouth2         = 10,      & ! (i  ,j-2)
      ice_blocksEast2          = 11,      & ! (i+2,j  )
      ice_blocksWest2          = 12         ! (i-2,j  )
!      ice_blocksNorthEast2     = 13,      & ! (i+2,j+2)
!      ice_blocksNorthWest2     = 14,      & ! (i-2,j+2)
!      ice_blocksSouthEast2     = 15,      & ! (i+2,j-2)
!      ice_blocksSouthWest2     = 16         ! (i-2,j-2)
   integer (int_kind), parameter, public :: &
      ice_blocksEastNorthEast  = 17,      & ! (i+2,j+1)
!      ice_blocksEastSouthEast  = 18,      & ! (i+2,j-1)
      ice_blocksWestNorthWest  = 19         ! (i-2,j+1)
!      ice_blocksWestSouthWest  = 20,      & ! (i-2,j-1)
!      ice_blocksNorthNorthEast = 21,      & ! (i+1,j-2)
!      ice_blocksSouthSouthEast = 22,      & ! (i+1,j-2)
!      ice_blocksNorthNorthWest = 23,      & ! (i-1,j+2)
!      ice_blocksSouthSouthWest = 24         ! (i-1,j-2)

! !PUBLIC DATA MEMBERS:

   integer (int_kind), public :: &
      nblocks_tot      ,&! total number of blocks in decomposition
      nblocks_x        ,&! tot num blocks in i direction
      nblocks_y          ! tot num blocks in j direction

   integer (int_kind), public, parameter ::  &
		trouble_ig = 2, &
		trouble_jg = 2
				
   integer (int_kind), public :: trouble_il,trouble_jl,trouble_ibl

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module private data
!
!-----------------------------------------------------------------------

   type (block), dimension(:), allocatable :: &
      all_blocks         ! block information for all blocks in domain

   integer (int_kind), dimension(:,:),allocatable :: &
      all_blocks_ij      ! block index stored in Cartesian order
                         !   useful for determining block index
                         !   of neighbor blocks

   integer (int_kind), dimension(:,:), allocatable, target :: &
      i_global,         &! global i index for each point in each block
      j_global           ! global j index for each point in each block

!EOC
!***********************************************************************

contains

   subroutine PrintTrouble_dbl(name,array)

     use ice_fileunits, only: nu_diag

     character(len=*), intent(in) :: name
     real(kind=dbl_kind), intent(in)         :: array(nx_block,ny_block) 
       
     write(nu_diag,*) TRIM(name),'(i,j) =',array(trouble_il,trouble_jl)

   end subroutine PrintTrouble_dbl

   subroutine PrintTrouble_int(name,array)

     use ice_fileunits, only: nu_diag

     character(len=*), intent(in) :: name
     integer(int_kind), intent(in)         :: array(nx_block,ny_block) 
        
     write(nu_diag,*) TRIM(name),'(i,j) =',array(trouble_il,trouble_jl)

   end subroutine PrintTrouble_int

   subroutine LocateTrouble

     type (block) :: this_block
     integer(int_kind) :: n,gbid,ib,ie,jb,je,ig,jg,i,j
     integer(int_kind) :: j_save,i_save,ib_save
     logical  :: foundi, foundj
    

      do n=1,nblocks_tot
	  foundi=.false.;foundj=.false.
	  gbid = all_blocks(n)%block_id
	  ib   = all_blocks(n)%ilo
	  ie   = all_blocks(n)%ihi
	  jb   = all_blocks(n)%jlo
          je   = all_blocks(n)%jhi
	  do i=ib,ie
	     ig=all_blocks(n)%i_glob(i) 
	     if(ig == trouble_ig) then 
		i_save=i
		foundi=.true.
	     endif
          enddo
	  do j=jb,je
	     jg=all_blocks(n)%j_glob(j) 
	     if(jg == trouble_jg) then
		j_save=j;foundj=.true.
	     endif
          enddo
          if(foundi .and. foundj) then 
		ib_save = n
	  endif
      enddo
      trouble_il  = i_save
      trouble_jl  = j_save
      trouble_ibl = ib_save
!      print *,'Trouble located: ib, i, j: ',ib_save,i_save,j_save

   end subroutine LocateTrouble




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
!  same as module
!
! !USES:

   use ice_fileunits, only: nu_diag
   use ice_communicate, only: my_task, master_task

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

   logical (log_kind) :: dbug

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
   allocate(all_blocks_ij(nblocks_x,nblocks_y))

!----------------------------------------------------------------------
!
!  fill block data structures for all blocks in domain
!
!----------------------------------------------------------------------

   n = 0
   do jblock=1,nblocks_y
      js = (jblock-1)*block_size_y + 1
      if (js > ny_global) call abort_ice(&
            'ice: create_blocks: Bad block decomp: ny_block too large?')
      je = js + block_size_y - 1
      if (je > ny_global) je = ny_global ! pad array

      do iblock=1,nblocks_x
         n = n + 1  ! global block id

         is = (iblock-1)*block_size_x + 1
         if (is > nx_global) call abort_ice(&
            'ice: create_blocks: Bad block decomp: nx_block too large?')
         ie = is + block_size_x - 1
         if (ie > nx_global) ie = nx_global

         all_blocks(n)%block_id = n
         all_blocks(n)%iblock   = iblock
         all_blocks(n)%jblock   = jblock
         all_blocks(n)%ilo      = nghost + 1
         all_blocks(n)%jlo      = nghost + 1
         all_blocks(n)%ihi      = nx_block - nghost ! default value
         all_blocks(n)%jhi      = ny_block - nghost ! default value

         if (jblock == nblocks_y .and. &
             (ns_boundary_type == 'tripole' .or. &
             ns_boundary_type == 'tripoleT')) then
             all_blocks(n)%tripole = .true.
         else
             all_blocks(n)%tripole = .false.
         endif
         all_blocks(n)%tripoleTFlag = (ns_boundary_type == 'tripoleT')

         all_blocks_ij(iblock,jblock) = n

         do j=1,ny_block
            j_global(j,n) = js - nghost + j - 1

            !*** southern ghost cells

            if (j_global(j,n) < 1) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) + ny_global
               case ('open')
                  j_global(j,n) = nghost - j + 1
               case ('closed')
                  j_global(j,n) = 0
               case ('tripole')
                  j_global(j,n) = nghost - j + 1 ! open
               case ('tripoleT')
                  j_global(j,n) = -j_global(j,n) + 1 ! open
               case default
                  call abort_ice(&
                           'ice: create_blocks: unknown n-s bndy type')
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
               case ('open')
                  j_global(j,n) = 2*ny_global - j_global(j,n) + 1
               case ('closed')
                  j_global(j,n) = 0
               case ('tripole')
                  j_global(j,n) = -j_global(j,n)
               case ('tripoleT')
                  j_global(j,n) = -j_global(j,n)
               case default
                  call abort_ice(&
                           'ice: create_blocks: unknown n-s bndy type')
               end select

            !*** set last physical point if padded domain

            else if (j_global(j,n) == ny_global .and. &
                     j > all_blocks(n)%jlo      .and. &
                     j < all_blocks(n)%jhi) then
               all_blocks(n)%jhi = j   ! last physical point in padded domain
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
               case ('open')
                  i_global(i,n) = nghost - i + 1
               case ('closed')
                  i_global(i,n) = 0
               case default
                  call abort_ice(&
                           'ice: create_blocks: unknown e-w bndy type')
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
               case ('open')
                  i_global(i,n) = 2*nx_global - i_global(i,n) + 1
               case ('closed')
                  i_global(i,n) = 0
               case default
                  call abort_ice(&
                           'ice: create_blocks: unknown e-w bndy type')
               end select

            !*** last physical point in padded domain

            else if (i_global(i,n) == nx_global .and. &
                     i > all_blocks(n)%ilo      .and. &
                     i < all_blocks(n)%ihi) then
               all_blocks(n)%ihi = i
            endif
         end do

         all_blocks(n)%i_glob => i_global(:,n)

      end do
   end do

!   dbug = .true.
   dbug = .false.
   if (dbug) then
      if (my_task == master_task) then
         write(nu_diag,*) 'block i,j locations'
         do n = 1, nblocks_tot
            write(nu_diag,*) 'block id, iblock, jblock:', &
            all_blocks(n)%block_id, &
            all_blocks(n)%iblock,   & 
            all_blocks(n)%jblock
         enddo
      endif
   endif

!EOC
!----------------------------------------------------------------------

end subroutine create_blocks

!***********************************************************************
!BOP
! !IROUTINE: ice_blocksGetNbrID
! !INTERFACE:

 function ice_blocksGetNbrID(blockID, direction, iBoundary, jBoundary) &
                             result (nbrID)

! !DESCRIPTION:
!  This function returns the block id of a neighboring block in a
!  requested direction.  Directions:
!      ice\_blocksNorth             (i  ,j+1)
!      ice\_blocksSouth             (i  ,j-1)
!      ice\_blocksEast              (i+1,j  )
!      ice\_blocksWest              (i-1,j  )
!      ice\_blocksNorthEast         (i+1,j+1)
!      ice\_blocksNorthWest         (i-1,j+1)
!      ice\_blocksSouthEast         (i  ,j-1)
!      ice\_blocksSouthWest         (i-1,j-1)
!      ice\_blocksNorth2            (i  ,j+2)
!      ice\_blocksSouth2            (i  ,j-2)
!      ice\_blocksEast2             (i+2,j  )
!      ice\_blocksWest2             (i-2,j  )
!      ice\_blocksNorthEast2        (i+2,j+2)
!      ice\_blocksNorthWest2        (i-2,j+2)
!      ice\_blocksSouthEast2        (i+2,j-2)
!      ice\_blocksSouthWest2        (i-2,j-2)
!      ice\_blocksEastNorthEast     (i+2,j+1)
!      ice\_blocksEastSouthEast     (i+2,j-1)
!      ice\_blocksWestNorthWest     (i-2,j+1)
!      ice\_blocksWestSouthWest     (i-2,j-1)
!      ice\_blocksNorthNorthEast    (i+1,j-2)
!      ice\_blocksSouthSouthEast    (i+1,j-2)
!      ice\_blocksNorthNorthWest    (i-1,j+2)
!      ice\_blocksSouthSouthWest    (i-1,j-2)
!

! !INPUT PARAMETERS:

   integer (int_kind), intent(in)  :: &
      blockID,       &! id of block for which neighbor id requested
      direction       ! direction for which to look for neighbor -
                      !   must be one of the predefined module
                      !   variables for block direction

   character (*), intent(in) :: &
      iBoundary,     &! determines what to do at edges of domain
      jBoundary       !  options are - open, closed, cyclic, tripole, tripoleT

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      nbrID           ! block ID of neighbor in requested dir

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
    
   integer (int_kind) :: &
      iBlock, jBlock,  &! i,j block location of current block
      inbr,   jnbr      ! i,j block location of neighboring block

!----------------------------------------------------------------------
!
!  retrieve info for current block
!
!----------------------------------------------------------------------

   call get_block_parameter(blockID, iblock=iBlock, jblock=jBlock)

!----------------------------------------------------------------------
!
!  compute i,j block location of neighbor
!
!----------------------------------------------------------------------

   select case(direction)

   case (ice_blocksNorth)

      inbr = iBlock
      jnbr = jBlock + 1
      if (jnbr > nblocks_y) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = 1
         case ('tripole':'tripoleT')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  nblocks_x - iBlock + 1 
            jnbr = -jBlock
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown north boundary')
         end select
      endif

   case (ice_blocksSouth)

      inbr = iBlock
      jnbr = jBlock - 1
      if (jnbr < 1) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = nblocks_y
         case ('tripole')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('tripoleT')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown south boundary')
         end select
      endif

   case (ice_blocksEast )

      inbr = iBlock + 1
      jnbr = jBlock
      if (inbr > nblocks_x) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = 1
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown east boundary')
         end select
      endif

   case (ice_blocksWest )

      jnbr = jBlock
      inbr = iBlock - 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = nblocks_x
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown west boundary')
         end select
      endif

   case (ice_blocksNorthEast)

      inbr = iBlock + 1
      jnbr = jBlock + 1
      if (inbr > nblocks_x) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = 1
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr > nblocks_y) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = 1
         case ('tripole':'tripoleT')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  nblocks_x - iBlock 
            if (inbr == 0) inbr = nblocks_x
            jnbr = -jBlock
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown north boundary')
         end select
      endif

   case (ice_blocksNorthWest)

      inbr = iBlock - 1
      jnbr = jBlock + 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = nblocks_x
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr > nblocks_y) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = 1
         case ('tripole':'tripoleT')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  nblocks_x - iBlock + 2 
            if (inbr > nblocks_x) inbr = 1
            jnbr = -jBlock
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown north boundary')
         end select
      endif

   case (ice_blocksSouthEast )

      inbr = iBlock + 1
      jnbr = jBlock - 1
      if (inbr > nblocks_x) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = 1
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = nblocks_y
         case ('tripole')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('tripoleT')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown south boundary')
         end select
      endif

   case (ice_blocksSouthWest )
      inbr = iBlock - 1
      jnbr = jBlock - 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = nblocks_x
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = nblocks_y
         case ('tripole')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('tripoleT')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown south boundary')
         end select
      endif

   case (ice_blocksEast2)

      inbr = iBlock + 2
      jnbr = jBlock
      if (inbr > nblocks_x) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - nblocks_x
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown east boundary')
         end select
      endif

   case (ice_blocksWest2)
      jnbr = jBlock
      inbr = iBlock - 2
      if (inbr < 1) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = nblocks_x + inbr
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown west boundary')
         end select
      endif

   case (ice_blocksEastNorthEast)

      inbr = iBlock + 2
      jnbr = jBlock + 1
      if (inbr > nblocks_x) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - nblocks_x
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr > nblocks_y) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = jnbr - nblocks_y
         case ('tripole':'tripoleT')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  nblocks_x - iBlock - 1 
            if (inbr <= 0) inbr = inbr + nblocks_x
            jnbr = -jBlock
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown north boundary')
         end select
      endif

   case (ice_blocksWestNorthWest)

      inbr = iBlock - 2
      jnbr = jBlock + 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = nblocks_x + inbr
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr > nblocks_y) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = jnbr + nblocks_y
         case ('tripole':'tripoleT')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  nblocks_x - iBlock + 3
            if (inbr > nblocks_x) inbr = inbr - nblocks_x
            jnbr = -jBlock
         case default
            call abort_ice( &
               'ice_blocksGetNbrID: unknown north boundary')
         end select
      endif

   case default

      call abort_ice( &
          'ice_blocksGetNbrID: unknown direction')
      return

   end select

!----------------------------------------------------------------------
!
!  now get block id for this neighbor block
!
!----------------------------------------------------------------------

   if (inbr > 0 .and. jnbr > 0) then
      nbrID = all_blocks_ij(inbr,jnbr)
   else if (inbr > 0 .and. jnbr < 0) then  ! tripole upper boundary
      !*** return negative value to flag tripole
      nbrID = -all_blocks_ij(inbr,abs(jnbr))
   else
      nbrID = 0   ! neighbor outside domain
   endif

!----------------------------------------------------------------------
!EOC

 end function ice_blocksGetNbrID

!**********************************************************************
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
      call abort_ice('ice: get_block: invalid block_id')
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

 subroutine get_block_parameter(block_id, local_id,           & 
                                ilo, ihi, jlo, jhi,           &
                                iblock, jblock, tripole,      &
                                i_glob, j_glob)

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
      local_id           ,&! local id assigned to block in current distrb
      ilo, ihi, jlo, jhi ,&! begin,end indices for physical domain
      iblock, jblock       ! cartesian i,j position for bloc

   logical (log_kind), intent(out), optional :: &
      tripole              ! flag is true if block on tripole bndy

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
      call abort_ice('ice: get_block_parameter: invalid block_id')
   endif

   if (present(local_id)) local_id = all_blocks(block_id)%local_id
   if (present(ilo     )) ilo      = all_blocks(block_id)%ilo
   if (present(ihi     )) ihi      = all_blocks(block_id)%ihi
   if (present(jlo     )) jlo      = all_blocks(block_id)%jlo
   if (present(jhi     )) jhi      = all_blocks(block_id)%jhi
   if (present(iblock  )) iblock   = all_blocks(block_id)%iblock
   if (present(jblock  )) jblock   = all_blocks(block_id)%jblock
   if (present(i_glob  )) i_glob   => all_blocks(block_id)%i_glob
   if (present(j_glob  )) j_glob   => all_blocks(block_id)%j_glob
   if (present(tripole )) tripole  = all_blocks(block_id)%tripole

!----------------------------------------------------------------------
!EOC

 end subroutine get_block_parameter

!**********************************************************************

 end module ice_blocks

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
