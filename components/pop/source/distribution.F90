!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 module distribution

!BOP
! !MODULE: distribution
!
! !DESCRIPTION:
!  This module provides data types and routines for distributing
!  blocks across processors.
!
! !REVISION HISTORY:
!  CVS:$Id: distribution.F90 43213 2013-01-11 18:24:48Z mlevy@ucar.edu $
!  CVS:$Name:  $

! !USES:

   use kinds_mod
   use communicate
   use blocks
   use POP_SpaceCurveMod
   use exit_mod
   use io_types, only: stdout

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: distrb  ! distribution data type
      integer (int_kind) :: &
         nprocs            ,&! number of processors in this dist
         communicator        ! communicator to use in this dist

      integer (int_kind), dimension(:), pointer :: &
         proc              ,&! processor location for this block
         local_block         ! block position in local array on proc
   end type

   logical, public :: sameDistribution

! !PUBLIC MEMBER FUNCTIONS:

   public :: create_distribution, &
             create_local_block_ids

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: create_distribution
! !INTERFACE:

 function create_distribution(dist_type, nprocs, workPerBlock)

! !DESCRIPTION:
!  This routine determines the distribution of blocks across processors
!  by call the appropriate subroutine based on distribution type
!  requested.  Currently only two distributions are supported:
!  2-d Cartesian distribution (cartesian) and a load-balanced
!  distribution (balanced) based on an input amount of work per
!  block.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      dist_type             ! method for distributing blocks
                            !  either cartesian or balanced

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distribution   ! resulting structure describing
                            !  distribution of blocks

!EOP
!BOC

   integer (int_kind) :: maxWork
   integer (int_kind), allocatable ::  work_per_block(:)

!----------------------------------------------------------------------
!
!  select the appropriate distribution type
!
!----------------------------------------------------------------------

   select case (trim(dist_type))

   case('cartesian')

      !------------------------------------------------
      ! The following comments and code were contributed 
      ! by John Dennis, CISL
      !
      ! This particular partitioning algorithm does not 
      ! handle land block elimination anyway 
      ! KLUDGE: probably should do something better here 
      !------------------------------------------------
      allocate(work_per_block(size(workPerBlock)))
      maxWork = MAXVAL(workPerBlock)
      work_per_block = maxWork
      create_distribution = create_distrb_cart(nprocs, work_per_block)
      deallocate(work_per_block)

   case('balanced')

      !------------------------------------------------
      ! The following comments and code were contributed 
      ! by John Dennis, CISL
      !
      ! This particular partitioning algorithm does not 
      ! handle land block elimination anyway 
      ! KLUDGE: probably should do something better here 
      !------------------------------------------------
      allocate(work_per_block(size(workPerBlock)))
      maxWork = MAXVAL(workPerBlock)
      work_per_block = maxWork
      create_distribution = create_distrb_balanced(nprocs, &
						   work_per_block)
      deallocate(work_per_block) 

   case('spacecurve')

      create_distribution = create_distrb_spacecurve(nprocs, &
						   workPerBlock)
   case('blockone')

      create_distribution = create_distrb_blockone(nprocs, &
						   workPerBlock)
   case default

      call exit_POP(sigAbort,'distribution: unknown distribution type')

   end select

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

!-----------------------------------------------------------------------
!
!  first determine number of local blocks to allocate array
!
!-----------------------------------------------------------------------

   bcount = 0
   do n=1,size(distribution%proc)
      if (distribution%proc(n) == my_task+1) bcount = bcount + 1
   end do

   if (bcount > 0) allocate(block_ids(bcount))

!-----------------------------------------------------------------------
!
!  now fill array with proper block ids
!
!-----------------------------------------------------------------------

   if (bcount > 0) then
      do n=1,size(distribution%proc)
         if (distribution%proc(n) == my_task+1) then
            block_ids(distribution%local_block(n)) = n
         endif
      end do
   endif

!EOC

 end subroutine create_local_block_ids

!***********************************************************************
!BOP
! !IROUTINE: create_distrb_cart
! !INTERFACE:

 function create_distrb_cart(nprocs, work_per_block)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a 2-d Cartesian distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs     ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_cart  ! resulting structure describing Cartesian
                          !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n               ,&! dummy loop indices
      iblock, jblock, nblck ,&!
      is, ie, js, je        ,&! start, end block indices for each proc
      local_block           ,&! block location on this processor
      nprocs_x              ,&! num of procs in x for global domain
      nprocs_y              ,&! num of procs in y for global domain
      nblocks_x_loc         ,&! num of blocks per processor in x
      nblocks_y_loc           ! num of blocks per processor in y

   type (distrb) :: dist  ! temp hold distribution

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(dist%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   dist%nprocs = nprocs

   call proc_decomposition(dist%nprocs, nprocs_x, nprocs_y)

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (dist%proc       (nblocks_tot), &
             dist%local_block(nblocks_tot))

!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   nblocks_x_loc = (nblocks_x-1)/nprocs_x + 1
   nblocks_y_loc = (nblocks_y-1)/nprocs_y + 1

   do j=1,nprocs_y
   do i=1,nprocs_x
      n = (j-1)*nprocs_x + i

      is = (i-1)*nblocks_x_loc + 1
      ie =  i   *nblocks_x_loc
      if (ie > nblocks_x) ie = nblocks_x
      js = (j-1)*nblocks_y_loc + 1
      je =  j   *nblocks_y_loc
      if (je > nblocks_y) je = nblocks_y

      local_block = 0
      do jblock = js,je
      do iblock = is,ie
         nblck = (jblock - 1)*nblocks_x + iblock
         if (work_per_block(nblck) /= 0) then
            local_block = local_block + 1
            dist%proc(nblck) = n
            dist%local_block(nblck) = local_block
         else
            dist%proc(nblck) = 0
            dist%local_block(nblck) = 0
         endif
      end do
      end do
   end do
   end do

!----------------------------------------------------------------------

   create_distrb_cart = dist  ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_cart
!***********************************************************************
!BOP
! !IROUTINE: create_distrb_blockone
! !INTERFACE:

 function create_distrb_blockone(nprocs, work_per_block)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a simple blocked distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs     ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_blockone  ! resulting structure describing blocked
                              !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j                  ,&! dummy loop indices
      iblock, jblock, nblck ,&!
      pe                    ,&! task number
      pecnt                   ! count blocks per pe
   real(r8)           :: &
      totwork               ,&! total work to do
      pework                ,&! work assigned to current pe
      petargetwork            ! target work for pe

   type (distrb) :: dist  ! temp hold distribution

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(dist%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   dist%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (dist%proc       (nblocks_tot), &
             dist%local_block(nblocks_tot))

!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   totwork = 0.0_r8
   do nblck = 1,nblocks_tot
      totwork = totwork + float(work_per_block(nblck))
   enddo

   pe = 1
   pecnt = 0
   pework = 0.0_r8
   petargetwork = totwork/float(nprocs-pe+1)
   do nblck = 1,nblocks_tot

! tcraig reminder of block ordering
!      nblck = (jblock - 1)*nblocks_x + iblock

      if (pecnt > 0 .and. pe < nprocs) then
         if (abs(pework-petargetwork) < abs(pework+work_per_block(nblck)-petargetwork)) then
            write(stdout,*) 'debug: cre_dist_rr ',pe,pecnt,pework,petargetwork
            pe = min(pe + 1,nprocs)
            pecnt = 0
            pework = 0.0
            petargetwork = totwork/float(nprocs-pe+1)
         endif
      endif
      pecnt = pecnt + 1
      pework = pework + float(work_per_block(nblck))
      totwork = totwork - float(work_per_block(nblck))

      if (work_per_block(nblck) /= 0) then
         dist%proc(nblck) = pe
         dist%local_block(nblck) = pecnt
      else
         dist%proc(nblck) = 0
         dist%local_block(nblck) = 0
      endif
   end do
   write(stdout,*) 'debug: cre_dist_rr ',pe,pecnt,pework,petargetwork

!----------------------------------------------------------------------

   create_distrb_blockone = dist  ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_blockone
!**********************************************************************
!BOP
! !IROUTINE: create_distrb_spacecurve
! !INTERFACE:

 function create_distrb_spacecurve(nprocs,work_per_block)

! !Description:
!  This function distributes blocks across processors in a 
!  load-balanced manner using space-filling curves
!
! !REVISION HISTORY:
!  added by J. Dennis 3/10/06 

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

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
      local_block          ,&! local block position on processor
      max_work             ,&! max amount of work in any block
      nprocs_x             ,&! num of procs in x for global domain
      nprocs_y               ! num of procs in y for global domain

   integer (int_kind), dimension(:),allocatable :: &
	idxT_i,idxT_j

   integer (int_kind), dimension(:,:),allocatable :: Mesh, Mesh2, Mesh3
   integer (int_kind) :: nblocksL,nblocks,ii,extra,i2,j2,tmp1,s1,ig

   integer (int_kind) :: ierr
   logical, parameter :: Debug = .FALSE.

   integer (int_kind), dimension(:), allocatable :: &
      priority           ,&! priority for moving blocks
      work_tmp           ,&! work per row or column for rake algrthm
      proc_tmp           ,&! temp processor id for rake algrthm
      block_count          ! counter to determine local block indx

   type (distrb) :: dist  ! temp hold distribution

   type (factor_t) :: xdim,ydim
   integer (int_kind) :: it,jj
   integer (int_kind) :: curveSize,sb_x,sb_y,itmp,numfac
   integer (int_kind) :: subNum, sfcNum 
   logical            :: foundx

!----------------------------------------------------------------------
!
!  first set up as Cartesian distribution
!  retain the Cartesian distribution if nblocks_tot = nprocs
!  to avoid processors with no work
!
!----------------------------------------------------------------------
   !------------------------------------------------------
   ! Space filling curves only work if:
   ! 
   ! 	nblocks_x = 2^m1 3^n1 5^o1 where m1,n1,o1 are integers
   ! 	nblocks_y = 2^m2 3^n2 5^o2 where m2,n2,o2 are integers
   !------------------------------------------------------
   if((.not. IsFactorable(nblocks_y)) .or. (.not. IsFactorable(nblocks_x))) then 
     create_distrb_spacecurve = create_distrb_cart(nprocs, work_per_block)
     return
   endif

   !-----------------------------------------------
   ! Factor the numbers of blocks in each dimension
   !-----------------------------------------------
   xdim = Factor(nblocks_x)
   ydim = Factor(nblocks_y)
   numfac = xdim%numfact

   !---------------------------------------------
   ! Match the common factors to create SFC curve
   !---------------------------------------------
   curveSize=1
   do it=1,numfac
      call MatchFactor(xdim,ydim,itmp,foundX)
      curveSize = itmp*curveSize
   enddo
   !--------------------------------------
   ! determine the size of the sub-blocks 
   ! within the space-filling curve 
   !--------------------------------------
   sb_x = ProdFactor(xdim)
   sb_y = ProdFactor(ydim)

   call create_communicator(dist%communicator, nprocs)

   dist%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (dist%proc       (nblocks_tot), &
             dist%local_block(nblocks_tot))
   dist%proc=0
   dist%local_block=0


!----------------------------------------------------------------------
!  Create the array to hold the SFC
!----------------------------------------------------------------------
   allocate(Mesh(curveSize,curveSize))
   allocate(Mesh2(nblocks_x,nblocks_y),Mesh3(nblocks_x,nblocks_y))
   Mesh  = 0
   Mesh2 = 0
   Mesh3 = 0

   allocate(idxT_i(nblocks_tot),idxT_j(nblocks_tot))


!----------------------------------------------------------------------
!  Generate the space-filling curve
!----------------------------------------------------------------------
   call GenSpaceCurve(Mesh)
   Mesh = Mesh + 1 ! make it 1-based indexing
   if(Debug) then 
     if(my_task ==0) call PrintCurve(Mesh)
   endif 
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
   !------------------------------------------------
   ! create a linear array of i,j coordinates of SFC
   !------------------------------------------------
   idxT_i=0;idxT_j=0
   do j=1,nblocks_y
     do i=1,nblocks_x
	n = (j-1)*nblocks_x + i
	ig = Mesh2(i,j)
	if(work_per_block(n) /= 0) then 
	    idxT_i(ig)=i;idxT_j(ig)=j
	endif
     enddo
   enddo
   !-----------------------------
   ! Compress out the land blocks 
   !-----------------------------
   ii=0
   do i=1,nblocks_tot
      if(IdxT_i(i) .gt. 0) then 
	 ii=ii+1
	 Mesh3(idxT_i(i),idxT_j(i)) = ii
      endif
   enddo
   if(Debug) then 
     if(my_task==0) call PrintCurve(Mesh3)
   endif
   
   nblocks=ii  
   nblocksL = nblocks/nprocs
   ! every cpu gets nblocksL blocks, but the first 'extra' get nblocksL+1
   extra = mod(nblocks,nprocs)
   s1 = extra*(nblocksL+1)
   ! split curve into two curves:
   ! 1 ... s1  s2 ... nblocks
   !
   !  s1 = extra*(nblocksL+1)         (count be 0)
   !  s2 = s1+1
   !
   ! First region gets nblocksL+1 blocks per partition
   ! Second region gets nblocksL blocks per partition
   if(Debug) print *,'nprocs,extra,nblocks,nblocksL,s1: ', &
		nprocs,extra,nblocks,nblocksL,s1

   do j=1,nblocks_y
   do i=1,nblocks_x
      n = (j-1)*nblocks_x + i
!      i2 = idxT_i(n)
!      j2 = idxT_j(n)
      ii = Mesh3(i,j)
      if(ii>0) then 
        !DBG if(my_task ==0) print *,'i,j,ii:= ',i,j,ii
        if(ii<=s1) then 
           ii=ii-1
           tmp1 = ii/(nblocksL+1)
	   dist%proc(n) = tmp1+1
        else
	   ii=ii-s1-1
	   tmp1 = ii/nblocksL
	   dist%proc(n) = extra + tmp1 + 1
        endif
      endif
   enddo
   enddo

!----------------------------------------------------------------------
!  Reset the dist data structure 
!----------------------------------------------------------------------

   allocate(proc_tmp(nprocs))
   proc_tmp = 0

   do n=1,nblocks_tot
      pid = dist%proc(n)
      if(pid>0) then 
        proc_tmp(pid) = proc_tmp(pid) + 1 	
        dist%local_block(n) = proc_tmp(pid)
      endif
   enddo

   if(Debug) then 
      if(my_task==0) print *,'dist%proc:= ',dist%proc
      print *,'IAM: ',my_task,' SpaceCurve: Number of blocks {total,local} :=', &
		nblocks_tot,nblocks,proc_tmp(my_task+1)
   endif

   deallocate(proc_tmp)
   ierr=1

   deallocate(Mesh,Mesh2,Mesh3)
   deallocate(idxT_i,idxT_j)
!----------------------------------------------------------------------
   create_distrb_spacecurve = dist  ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_spacecurve

!**********************************************************************
!BOP
! !IROUTINE: create_distrb_balanced
! !INTERFACE:

 function create_distrb_balanced(nprocs, work_per_block)

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
      work_per_block        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_balanced  ! resulting structure describing
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
      local_block          ,&! local block position on processor
      max_work             ,&! max amount of work in any block
      nprocs_x             ,&! num of procs in x for global domain
      nprocs_y               ! num of procs in y for global domain

   integer (int_kind), dimension(:), allocatable :: &
      priority           ,&! priority for moving blocks
      work_tmp           ,&! work per row or column for rake algrthm
      proc_tmp           ,&! temp processor id for rake algrthm
      block_count          ! counter to determine local block indx

   type (distrb) :: dist  ! temp hold distribution

!----------------------------------------------------------------------
!
!  first set up as Cartesian distribution
!  retain the Cartesian distribution if nblocks_tot = nprocs
!  to avoid processors with no work
!
!----------------------------------------------------------------------

   dist = create_distrb_cart(nprocs, work_per_block)
   if (nblocks_tot == nprocs) then
      create_distrb_balanced = dist  ! return the result
      return
   endif

!----------------------------------------------------------------------
!
!  now re-distribute blocks using a rake in each direction
!
!----------------------------------------------------------------------

   max_work = maxval(work_per_block)

   call proc_decomposition(dist%nprocs, nprocs_x, nprocs_y)

!----------------------------------------------------------------------
!
!  load-balance using a rake algorithm in the x-direction first
!
!----------------------------------------------------------------------

   allocate(priority(nblocks_tot))

   !*** set highest priority such that eastern-most blocks
   !*** and blocks with the least amount of work are
   !*** moved first

   do j=1,nblocks_y
   do i=1,nblocks_x
      n=(j-1)*nblocks_x + i
      if (work_per_block(n) > 0) then
         priority(n) = (max_work + 1)*(nblocks_x + i) - &
                       work_per_block(n)
      else
         priority(n) = 0
      endif
   end do
   end do

   allocate(work_tmp(nprocs_x), &
            proc_tmp(nprocs_x))

   do j=1,nprocs_y

      work_tmp(:) = 0
      do i=1,nprocs_x
         pid = (j-1)*nprocs_x + i
         proc_tmp(i) = pid
         do n=1,nblocks_tot
            if (dist%proc(n) == pid) then
               work_tmp(i) = work_tmp(i) + work_per_block(n)
            endif
         end do
      end do

      call rake (work_tmp, proc_tmp, work_per_block, priority, dist)

   end do

   deallocate(work_tmp, proc_tmp)

!----------------------------------------------------------------------
!
!  use a rake algorithm in the y-direction now
!
!----------------------------------------------------------------------

   !*** set highest priority for northern-most blocks

   do j=1,nblocks_y
   do i=1,nblocks_x
      n=(j-1)*nblocks_x + i
      if (work_per_block(n) > 0) then
         priority(n) = (max_work + 1)*(nblocks_y + j) - &
                       work_per_block(n)
      else
         priority(n) = 0
      endif
   end do
   end do

   allocate(work_tmp(nprocs_y), &
            proc_tmp(nprocs_y))

   do i=1,nprocs_x

      work_tmp(:) = 0
      do j=1,nprocs_y
         pid = (j-1)*nprocs_x + i
         proc_tmp(j) = pid
         do n=1,nblocks_tot
            if (dist%proc(n) == pid) then
               work_tmp(j) = work_tmp(j) + work_per_block(n)
            endif
         end do
      end do

      call rake (work_tmp, proc_tmp, work_per_block, priority, dist)

   end do

   deallocate(work_tmp, proc_tmp)
   deallocate(priority)

!----------------------------------------------------------------------
!
!  reset local_block info based on new distribution
!
!----------------------------------------------------------------------

   allocate(proc_tmp(nprocs))
   proc_tmp = 0

   do pid=1,nprocs
      local_block = 0
      do n=1,nblocks_tot
         if (dist%proc(n) == pid) then
            local_block = local_block + 1
            dist%local_block(n) = local_block
            proc_tmp(pid) = proc_tmp(pid) + 1
         endif
      end do
   end do

   if (minval(proc_tmp) < 1) then
      call exit_POP(sigAbort,'Load-balanced distribution failed')
   endif

   deallocate(proc_tmp)

!----------------------------------------------------------------------

   create_distrb_balanced = dist  ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_balanced

!**********************************************************************
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

   real (r4) :: &
      square                       ! square root of nprocs

!----------------------------------------------------------------------
!
!  start with an initial guess that is closest to square decomp
!
!----------------------------------------------------------------------

   square = sqrt(real(nprocs))
   nprocs_x = 0
   nprocs_y = 0

   iguess = nint(square)

!----------------------------------------------------------------------
!
!  try various decompositions to find the best
!
!----------------------------------------------------------------------

   proc_loop: do
      jguess = nprocs/iguess

      if (iguess*jguess == nprocs) then ! valid decomp

         !***
         !*** if the blocks can be evenly distributed, it is a
         !*** good decomposition
         !***

         if (mod(nblocks_x,iguess) == 0 .and. &
             mod(nblocks_y,jguess) == 0) then
            nprocs_x = iguess
            nprocs_y = jguess
            exit proc_loop

         !***
         !*** if the blocks can be evenly distributed in a
         !*** transposed direction, it is a good decomposition
         !***

         else if (mod(nblocks_x,jguess) == 0 .and. &
                mod(nblocks_y,iguess) == 0) then
            nprocs_x = jguess
            nprocs_y = iguess
            exit proc_loop

         !***
         !*** A valid decomposition, but keep searching for
         !***  a better one
         !***

         else
            if (nprocs_x == 0) then
               nprocs_x = iguess
               nprocs_y = jguess
            endif
            iguess = iguess - 1
            if (iguess == 0) then
               exit proc_loop
            else
               cycle proc_loop
            endif
         endif

      else ! invalid decomp - keep trying

         iguess = iguess - 1
         if (iguess == 0) then
            exit proc_loop
         else
            cycle proc_loop
         endif
      endif
   end do proc_loop

   if (nprocs_x == 0) then
      call exit_POP(sigAbort,'Unable to find 2d processor config')
   endif

!----------------------------------------------------------------------
!EOC

 end subroutine proc_decomposition

!**********************************************************************
!BOP
! !IROUTINE: rake
! !INTERFACE:

 subroutine rake (proc_work, proc_id, block_work, priority, dist)

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

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), intent(inout), dimension(:) :: &
      proc_work           ,&! amount of work per processor
      priority              ! priority for moving a given block

   integer (int_kind), intent(in), dimension(:) :: &
      block_work          ,&! amount of work per block
      proc_id               ! global processor number

   type (distrb), intent(inout) :: &
      dist                  ! distribution to change

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n, m, np1, &
      iproc, inext, &
      nprocs, nblocks, &
      last_priority, last_loc, &
      residual, &
      work_mean, work_max, work_diff, &
      iter, niters, itransfer, ntransfers, &
      min_priority

!----------------------------------------------------------------------
!
!  initialization
!
!----------------------------------------------------------------------

   nprocs  = size(proc_work)
   nblocks = size(block_work)

   !*** mean work per processor

   work_mean = sum(proc_work)/nprocs + 1
   work_max  = maxval(proc_work)
   residual = mod(work_mean,nprocs)

   min_priority = 1000000
   do n=1,nprocs
      iproc = proc_id(n)
      do i=1,nblocks
         if (dist%proc(i) == iproc) then
            min_priority = min(min_priority,priority(i))
         endif
      end do
   end do

!----------------------------------------------------------------------
!
!  do two sets of transfers
!
!----------------------------------------------------------------------

   transfer_loop: do

!----------------------------------------------------------------------
!
!     do rake across the processors
!
!----------------------------------------------------------------------

      ntransfers = 0
       do n=1,nprocs
          if (n < nprocs) then
             np1   = n+1
          else
             np1   = 1
          endif
         iproc = proc_id(n)
         inext = proc_id(np1)

         if (proc_work(n) > work_mean) then !*** pass work to next
            work_diff = proc_work(n) - work_mean

            rake1: do while (work_diff > 1)

               !*** attempt to find a block with the required
               !*** amount of work and with the highest priority
               !*** for moving (eg boundary blocks first)

               last_priority = 0
               last_loc = 0
               do i=1,nblocks
                  if (dist%proc(i) == iproc) then
                     if (priority(i) > last_priority ) then
                        last_priority = priority(i)
                        last_loc = i
                     endif
                  endif
               end do
               if (last_loc == 0) exit rake1 ! could not shift work

               ntransfers = ntransfers + 1
               dist%proc(last_loc) = inext
               if (np1 == 1) priority(last_loc) = min_priority
               work_diff = work_diff - block_work(last_loc)

               proc_work(n  ) = proc_work(n  )-block_work(last_loc)
               proc_work(np1) = proc_work(np1)+block_work(last_loc)
            end do rake1
         endif

      end do

!----------------------------------------------------------------------
!
!     increment work_mean by one and repeat
!
!----------------------------------------------------------------------

      work_mean = work_mean + 1
      if (ntransfers == 0 .or. work_mean > work_max) exit transfer_loop

   end do transfer_loop

!----------------------------------------------------------------------
!EOC

end subroutine rake

!***********************************************************************

end module distribution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
