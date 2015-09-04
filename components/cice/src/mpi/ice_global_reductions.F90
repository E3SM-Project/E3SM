!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!BOP
! !MODULE: ice_global_reductions

 module ice_global_reductions

! !DESCRIPTION:
!  This module contains all the routines for performing global
!  reductions like global sums, minvals, maxvals, etc.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_global_reductions.F90 112 2008-03-13 21:06:56Z eclare $
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
! Feb. 2008: Updated from POP version by Elizabeth C. Hunke, LANL
!
! !USES:

   use ice_kinds_mod
   use ice_communicate
   use ice_constants
   use ice_blocks
   use ice_distribution
   use ice_domain_size

   implicit none

   private
   include 'mpif.h'
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: global_sum,      &
             global_sum_prod, &
             global_maxval,   &
             global_minval,   &
             init_global_reductions

   public :: sum_vector_dbl
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  generic interfaces for module procedures
!
!-----------------------------------------------------------------------

   interface global_sum
     module procedure global_sum_dbl,              &
                      global_sum_real,             &
                      global_sum_int,              &
                      global_sum_scalar_dbl,       &
                      global_sum_scalar_real,      &
                      global_sum_scalar_int
   end interface

   interface global_sum_prod
     module procedure global_sum_prod_dbl,         &
                      global_sum_prod_real,        &
                      global_sum_prod_int
   end interface

   interface global_maxval
     module procedure global_maxval_dbl,           &
                      global_maxval_real,          &
                      global_maxval_int,           &
                      global_maxval_scalar_dbl,    &
                      global_maxval_scalar_real,   &
                      global_maxval_scalar_int
   end interface

   interface global_minval
     module procedure global_minval_dbl,           &
                      global_minval_real,          &
                      global_minval_int,           &
                      global_minval_scalar_dbl,    &
                      global_minval_scalar_real,   &
                      global_minval_scalar_int
   end interface

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   logical(log_kind) :: ltripole_grid  ! in lieu of use domain

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_global_reductions
! !INTERFACE:

 subroutine init_global_reductions(tripole_flag)

! !DESCRIPTION:
!  Initializes necessary buffers for global reductions.
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:
!
   logical(log_kind), intent(in) :: tripole_flag
!
!EOP
!BOC

! This flag is apparently never used; if it were used, it might need
! a corresponding tripoleTFlag to be defined.
   ltripole_grid = tripole_flag

!EOC

 end subroutine init_global_reductions

 subroutine sum_vector_dbl(local_vector,global_vector, dist)

!-----------------------------------------------------------------------
!
!  this function returns the sum of vector value across processors
!
!-----------------------------------------------------------------------

   type (distrb), intent(in) :: &
      dist                 ! distribution from which this is called

   real (dbl_kind), intent(inout) :: &
      local_vector(:)                ! local vector to be compared

   real (dbl_kind) :: global_vector(:)   ! resulting global sum

   integer (int_kind) :: ierr ! MPI error flag

   integer (int_kind) :: len
!-----------------------------------------------------------------------

   len = size(local_vector)
   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_vector, global_vector, len, &
                            mpiR8, MPI_SUM, dist%communicator, ierr)
      else
         global_vector = c0
      endif
   else
      global_vector = local_vector
   endif

!-----------------------------------------------------------------------

 end subroutine sum_vector_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

 function global_sum_dbl(array, dist, field_loc, mMask, lMask) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_sum
!  function corresponding to double precision arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

! !USES:

! !INPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      array                ! array to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind), dimension(:), allocatable :: &
      blockSum,     &! sum of local block domain
      localSum,     &! sum of all local block domains
      globalSumTmp   ! higher precision global sum

   integer (int_kind) :: &
      i,j,iblock,n, &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      ierr,         &! mpi error flag
      blockID,      &! block location
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator, &! communicator for this distribution
      nreduce,      &! mpi count
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block     ! holds local block information

!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   nreduce = nblocks_tot
#else
   nreduce = 1
#endif
   allocate(blockSum(nreduce), &
            globalSumTmp(nreduce))
   blockSum     = 0.0_dbl_kind
   globalSumTmp = 0.0_dbl_kind
   globalSum    = 0.0_dbl_kind

   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,       &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

#ifdef REPRODUCIBLE
      n = blockID
#else
      n = 1
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum(n) = &
            blockSum(n) + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum(n) = &
               blockSum(n) + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum(n) = blockSum(n) + array(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
            field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
 
         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum(n) = &
                     blockSum(n) - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                     blockSum(n) = blockSum(n) - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum(n) = blockSum(n) - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif
   end do

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(blockSum, globalSumTmp, nreduce, &
                         mpiR8, MPI_SUM, communicator, ierr)
   endif

   do n=1,nreduce
      globalSum = globalSum + globalSumTmp(n)
   enddo

   deallocate(blockSum, globalSumTmp)

!-----------------------------------------------------------------------
!EOC

 end function global_sum_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

 function global_sum_real(array, dist, field_loc, mMask, lMask) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_sum
!  function corresponding to real arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

! !USES:

! !INPUT PARAMETERS:

   real (real_kind), dimension(:,:,:), intent(in) :: &
      array                ! array to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (real_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (real_kind) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   real (dbl_kind) :: &
      blockSum,     &! sum of local block domain
      localSum,     &! sum of all local block domains
      globalSumTmp   ! higher precision global sum
#else
   real (real_kind) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains
#endif

   integer (int_kind) :: &
      i,j,iblock,   &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      ierr,         &! mpi error flag
      blockID,      &! block location
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator, &! communicator for this distribution
      maxiglob       ! maximum non-redundant value of i_global
 
   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block     ! holds local block information

!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   localSum  = 0.0_dbl_kind
#else
   localSum  = 0.0_real_kind
#endif
   globalSum = 0.0_real_kind

   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,       &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

#ifdef REPRODUCIBLE
      blockSum = 0.0_dbl_kind
#else
      blockSum = 0.0_real_kind
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum = &
            blockSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum = &
               blockSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum = blockSum + array(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
            field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
 
         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum = &
                     blockSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                     blockSum = blockSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum = blockSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSumTmp, 1, &
                         mpiR8, MPI_SUM, communicator, ierr)
      globalSum = globalSumTmp
   endif
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         mpiR4, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------
!EOC

 end function global_sum_real

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

 function global_sum_int(array, dist, field_loc, mMask, lMask) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_sum
!  function corresponding to integer arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

! !USES:

! !INPUT PARAMETERS:

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      array                ! array to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   integer (int_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,   &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      ierr,         &! mpi error flag
      blockID,      &! block location
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator, &! communicator for this distribution
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block     ! holds local block information

!-----------------------------------------------------------------------

   localSum  = 0_int_kind
   globalSum = 0_int_kind

   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,       &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockSum = 0

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum = &
            blockSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum = &
               blockSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum = blockSum + array(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
            field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
 
         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum = &
                     blockSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                     blockSum = blockSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum = blockSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function global_sum_int

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

 function global_sum_scalar_dbl(scalar, dist) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of a set of scalars distributed across
!  a parallel machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_sum
!  function corresponding to double precision scalars.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

! !USES:

! !INPUT PARAMETERS:

   real (dbl_kind), intent(in) :: &
      scalar               ! scalar to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,         &! mpi error flag
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator   ! communicator for this distribution

!#ifdef REPRODUCIBLE
!   real (r16_kind) :: &
!      scalarTmp, globalSumTmp  ! higher precision for reproducibility
!#endif

!-----------------------------------------------------------------------
!
!  get communicator for MPI calls
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!  REPRODUCIBLE option is commented out because MPI does not handle 
!  REAL16 correctly.
!
!-----------------------------------------------------------------------

!#ifdef REPRODUCIBLE
!   if (my_task < numProcs) then
!      scalarTmp = scalar
!      call MPI_ALLREDUCE(scalarTmp, globalSumTmp, 1, &
!                         mpiR16, MPI_SUM, communicator, ierr)
!      globalSum = globalSumTmp
!   endif
!#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalSum, 1, &
                         mpiR8, MPI_SUM, communicator, ierr)
   endif
!#endif

!-----------------------------------------------------------------------
!EOC

 end function global_sum_scalar_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

 function global_sum_scalar_real(scalar, dist) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of a set of scalars distributed across
!  a parallel machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_sum
!  function corresponding to real scalars.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

! !USES:

! !INPUT PARAMETERS:

   real (real_kind), intent(in) :: &
      scalar               ! scalar to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

! !OUTPUT PARAMETERS:

   real (real_kind) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,         &! mpi error flag
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator   ! communicator for this distribution

#ifdef REPRODUCIBLE
   real (dbl_kind) :: &
      scalarTmp, globalSumTmp  ! higher precision for reproducibility
#endif

!-----------------------------------------------------------------------
!
!  get communicator for MPI calls
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   if (my_task < numProcs) then
      scalarTmp = scalar
      call MPI_ALLREDUCE(scalarTmp, globalSumTmp, 1, &
                         mpiR8, MPI_SUM, communicator, ierr)
      globalSum = globalSumTmp
   endif
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalSum, 1, &
                         mpiR4, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------
!EOC

 end function global_sum_scalar_real

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

 function global_sum_scalar_int(scalar, dist) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of a set of scalars distributed across
!  a parallel machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_sum
!  function corresponding to integer scalars.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

! !USES:

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      scalar               ! scalar to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,         &! mpi error flag
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator   ! communicator for this distribution

!-----------------------------------------------------------------------
!
!  get communicator for MPI calls
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalSum, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function global_sum_scalar_int

!***********************************************************************
!BOP
! !IROUTINE: global_sum_prod
! !INTERFACE:

 function global_sum_prod_dbl (array1, array2, dist, field_loc, &
                               mMask, lMask) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a product of
!  two 2-d arrays.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  global_sum_prod function corresponding to double precision arrays.
!  The generic interface is identical but will handle real and integer 
!  2-d slabs.

! !USES:

! !INPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      array1, array2       ! arrays whose product is to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for arrays

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind), dimension(:), allocatable :: &
      blockSum,     &! sum of local block domain
      localSum,     &! sum of all local block domains
      globalSumTmp   ! higher precision global sum

   integer (int_kind) :: &
      i,j,iblock,n,    &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      nreduce,         &! mpi count
      maxiglob          ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block     ! holds local block information

!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   nreduce = nblocks_tot
#else
   nreduce = 1
#endif
   allocate(blockSum(nreduce), &
            globalSumTmp(nreduce))
   blockSum     = 0.0_dbl_kind
   globalSumTmp = 0.0_dbl_kind
   globalSum    = 0.0_dbl_kind

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

#ifdef REPRODUCIBLE
      n = blockID
#else
      n = 1
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum(n) = &
            blockSum(n) + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum(n) = &
               blockSum(n) + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum(n) = blockSum(n) + array1(i,j,iblock)*array2(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
            field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
 
         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum(n) = &
                     blockSum(n) - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                        blockSum(n) = blockSum(n) - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum(n) = blockSum(n) - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

   end do

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(blockSum, globalSumTmp, nreduce, &
                         mpiR8, MPI_SUM, communicator, ierr)
   endif

   do n=1,nreduce
      globalSum = globalSum + globalSumTmp(n)
   enddo

   deallocate(blockSum, globalSumTmp)

!-----------------------------------------------------------------------
!EOC

 end function global_sum_prod_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_sum_prod
! !INTERFACE:

 function global_sum_prod_real (array1, array2, dist, field_loc, &
                                mMask, lMask) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a product of
!  two 2-d arrays.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  global_sum_prod function corresponding to single precision arrays.
!  The generic interface is identical but will handle real and integer 
!  2-d slabs.

! !USES:

! !INPUT PARAMETERS:

   real (real_kind), dimension(:,:,:), intent(in) :: &
      array1, array2       ! arrays whose product is to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for arrays

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (real_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (real_kind) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   real (dbl_kind) :: &
      blockSum,      &! sum of local block domain
      localSum,      &! sum of all local block domains
      globalSumTmp    ! higher precision for reproducibility
#else
   real (real_kind) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains
#endif

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      maxiglob          ! maximum non-redundant value of i_global
 
   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   localSum  = 0.0_dbl_kind
#else
   localSum  = 0.0_real_kind
#endif
   globalSum = 0.0_real_kind

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

#ifdef REPRODUCIBLE
      blockSum = 0.0_dbl_kind
#else
      blockSum = 0.0_real_kind
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum = &
            blockSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum = &
               blockSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum = blockSum + array1(i,j,iblock)*array2(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
            field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
 
         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum = &
                     blockSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                        blockSum = blockSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum = blockSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSumTmp, 1, &
                         mpiR8, MPI_SUM, communicator, ierr)
      globalSum = globalSumTmp
   endif
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         mpiR4, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------
!EOC

 end function global_sum_prod_real

!***********************************************************************
!BOP
! !IROUTINE: global_sum_prod
! !INTERFACE:

 function global_sum_prod_int (array1, array2, dist, field_loc, &
                               mMask, lMask) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a product of
!  two 2-d arrays.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  global_sum_prod function corresponding to integer arrays.
!  The generic interface is identical but will handle real and integer 
!  2-d slabs.

! !USES:

! !INPUT PARAMETERS:

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      array1, array2       ! arrays whose product is to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for arrays

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   integer (int_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      maxiglob          ! maximum non-redundant value of i_global
 
   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------

   localSum  = 0_int_kind
   globalSum = 0_int_kind

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockSum = 0

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum = &
            blockSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum = &
               blockSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum = blockSum + array1(i,j,iblock)*array2(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
            field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
 
         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum = &
                     blockSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                        blockSum = blockSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     blockSum = blockSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function global_sum_prod_int

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

 function global_maxval_dbl (array, dist, lMask) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of the physical domain of a 2-d field
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to double precision arrays.  

! !INPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
      globalMaxval         ! resulting maximum value of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind) ::    &
      blockMaxval,     &! sum of local block domain
      localMaxval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------

   localMaxval  = -HUGE(0.0_dbl_kind)
   globalMaxval = -HUGE(0.0_dbl_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMaxval = -HUGE(0.0_dbl_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMaxval = max(blockMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMaxval = max(blockMaxval,array(i,j,iblock))
         end do
         end do
      endif

      localMaxval = max(localMaxval,blockMaxval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
                         mpiR8, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_maxval_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

 function global_maxval_real (array, dist, lMask) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of the physical domain of a 2-d field
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to single precision arrays.  

! !INPUT PARAMETERS:

   real (real_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (real_kind) :: &
      globalMaxval         ! resulting maximum value of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (real_kind) ::    &
      blockMaxval,     &! sum of local block domain
      localMaxval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------

   localMaxval  = -HUGE(0.0_real_kind)
   globalMaxval = -HUGE(0.0_real_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMaxval = -HUGE(0.0_real_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMaxval = max(blockMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMaxval = max(blockMaxval,array(i,j,iblock))
         end do
         end do
      endif

      localMaxval = max(localMaxval,blockMaxval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
                         mpiR4, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_maxval_real

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

 function global_maxval_int (array, dist, lMask) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of the physical domain of a 2-d field
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to integer arrays.  

! !INPUT PARAMETERS:

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      globalMaxval         ! resulting maximum value of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::    &
      blockMaxval,     &! sum of local block domain
      localMaxval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------

   localMaxval  = -HUGE(0_int_kind)
   globalMaxval = -HUGE(0_int_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMaxval = -HUGE(0_int_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMaxval = max(blockMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMaxval = max(blockMaxval,array(i,j,iblock))
         end do
         end do
      endif

      localMaxval = max(localMaxval,blockMaxval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_maxval_int

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

 function global_maxval_scalar_dbl (scalar, dist) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to double precision scalars.  

! !INPUT PARAMETERS:

   real (dbl_kind), intent(in) :: &
      scalar               ! scalar for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
      globalMaxval         ! resulting maximum value

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
                         mpiR8, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_maxval_scalar_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

 function global_maxval_scalar_real (scalar, dist) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to single precision scalars.  

! !INPUT PARAMETERS:

   real (real_kind), intent(in) :: &
      scalar               ! scalar for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

! !OUTPUT PARAMETERS:

   real (real_kind) :: &
      globalMaxval         ! resulting maximum value

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
                         mpiR4, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_maxval_scalar_real

!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

 function global_maxval_scalar_int (scalar, dist) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to single precision scalars.  

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      scalar               ! scalar for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      globalMaxval         ! resulting maximum value

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_maxval_scalar_int

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

 function global_minval_dbl (array, dist, lMask) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of the physical domain of a 2-d field
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_minval
!  function corresponding to double precision arrays.  

! !INPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
      globalMinval         ! resulting minimum value of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind) ::    &
      blockMinval,     &! sum of local block domain
      localMinval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------

   localMinval  = HUGE(0.0_dbl_kind)
   globalMinval = HUGE(0.0_dbl_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMinval = HUGE(0.0_dbl_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMinval = min(blockMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMinval = min(blockMinval,array(i,j,iblock))
         end do
         end do
      endif

      localMinval = min(localMinval,blockMinval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
                         mpiR8, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_minval_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

 function global_minval_real (array, dist, lMask) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of the physical domain of a 2-d field
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_minval
!  function corresponding to single precision arrays.  

! !INPUT PARAMETERS:

   real (real_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (real_kind) :: &
      globalMinval         ! resulting minimum value of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (real_kind) ::    &
      blockMinval,     &! sum of local block domain
      localMinval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------

   localMinval  = HUGE(0.0_real_kind)
   globalMinval = HUGE(0.0_real_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMinval = HUGE(0.0_real_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMinval = min(blockMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMinval = min(blockMinval,array(i,j,iblock))
         end do
         end do
      endif

      localMinval = min(localMinval,blockMinval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
                         mpiR4, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_minval_real

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

 function global_minval_int (array, dist, lMask) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of the physical domain of a 2-d field
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_minval
!  function corresponding to integer arrays.  

! !INPUT PARAMETERS:

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      globalMinval         ! resulting minimum value of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::    &
      blockMinval,     &! sum of local block domain
      localMinval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------

   localMinval  = HUGE(0_int_kind)
   globalMinval = HUGE(0_int_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMinval = HUGE(0_int_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMinval = min(blockMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMinval = min(blockMinval,array(i,j,iblock))
         end do
         end do
      endif

      localMinval = min(localMinval,blockMinval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
                         MPI_INTEGER, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_minval_int

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

 function global_minval_scalar_dbl (scalar, dist) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_minval
!  function corresponding to double precision scalars.  

! !INPUT PARAMETERS:

   real (dbl_kind), intent(in) :: &
      scalar               ! scalar for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
      globalMinval         ! resulting minimum value

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMinval, 1, &
                         mpiR8, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_minval_scalar_dbl

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

 function global_minval_scalar_real (scalar, dist) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_minval
!  function corresponding to single precision scalars.  

! !INPUT PARAMETERS:

   real (real_kind), intent(in) :: &
      scalar               ! scalar for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

! !OUTPUT PARAMETERS:

   real (real_kind) :: &
      globalMinval         ! resulting minimum value

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMinval, 1, &
                         mpiR4, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_minval_scalar_real

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

 function global_minval_scalar_int (scalar, dist) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_minval
!  function corresponding to single precision scalars.  

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      scalar               ! scalar for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      globalMinval         ! resulting minimum value

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMinval, 1, &
                         MPI_INTEGER, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------

 end function global_minval_scalar_int

!***********************************************************************

 end module ice_global_reductions

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
