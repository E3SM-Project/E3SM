!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ice_gather_scatter

 module ice_gather_scatter

! !DESCRIPTION:
!  This module contains routines that mimic the behavior of the mpi
!  version in the case of a single processor:  gathering data to a single
!  processor from a distributed array, and scattering data from a
!  single processor to a distributed array.
!
! !REVISION HISTORY:
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

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

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
     src_block        ! local block index in source distribution

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      this_block = get_block(n,n)

      !*** copy local blocks

      if (src_dist%blockLocation(n) /= 0) then

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

      else !*** fill land blocks with special values

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = spval_dbl
         end do
         end do

      endif

   end do

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
      i,j,n             ,&! dummy loop counters
      src_block           ! local block index in source distribution

   type (block) :: &
      this_block   ! block info for current block

!-----------------------------------------------------------------------
!
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      this_block = get_block(n,n)

      !*** copy local blocks

      if (src_dist%blockLocation(n) /= 0) then

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

      else !*** fill land blocks with zeroes

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = 0._real_kind
         end do
         end do

      endif

   end do

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
      i,j,n,             &! dummy loop counters
      src_block           ! local block index in source distribution

   type (block) :: &
      this_block    ! block info for current block

!-----------------------------------------------------------------------
!
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      this_block = get_block(n,n)

      !*** copy local blocks

      if (src_dist%blockLocation(n) /= 0) then

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

      else !*** fill land blocks with zeroes

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = 0
         end do
         end do

      endif

   end do

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
!  This subroutine scatters a global-sized array on the processor
!  src\_task to a distribution of blocks given by dst\_dist.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays
!  corresponding to the generic interface scatter_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific interface based
!  on the data type of the input argument).

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
     isrc, jsrc,         &! source addresses
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     dst_block            ! local block index in dest distribution

   type (block) :: &
     this_block  ! block info for current block

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
   case (field_type_noupdate) ! ghost cells not needed - use cell center
      isign =  1
   case default
      call abort_ice('Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  copy blocks of global array into local block distribution
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      if (dst_dist%blockLocation(n) /= 0) then

         this_block = get_block(n,n)
         dst_block  = dst_dist%blockLocalID(n)

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
      endif ! dst block not land
   end do  ! block loop

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
!  This subroutine scatters a global-sized array on the processor
!  src\_task to a distribution of blocks given by dst\_dist.
!
!-----------------------------------------------------------------------
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
     isrc, jsrc,         &! source addresses
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     dst_block            ! local block index in dest distribution

   type (block) :: &
     this_block  ! block info for current block

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
   case (field_type_noupdate) ! ghost cells not needed - use cell center
      isign =  1
   case default
      call abort_ice('Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  copy blocks of global array into local block distribution
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      if (dst_dist%blockLocation(n) /= 0) then

         this_block = get_block(n,n)
         dst_block  = dst_dist%blockLocalID(n)

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
      endif ! dst block not land
   end do  ! block loop

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
!  This subroutine scatters a global-sized array on the processor
!  src\_task to a distribution of blocks given by dst\_dist.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   integer (int_kind), dimension(:,:), intent(in) :: &
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

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,bid,         &! dummy loop indices
      isrc, jsrc,        &! source addresses
      xoffset, yoffset,  &! offsets for tripole boundary conditions
      isign,             &! sign factor for tripole boundary conditions
      yoffset2,          &!
      dst_block           ! local block index in dest distribution

   type (block) :: &
      this_block  ! block info for current block

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
   case (field_type_noupdate) ! ghost cells not needed - use cell center
      isign =  1
   case default
      call abort_ice('Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  copy blocks of global array into local block distribution
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      if (dst_dist%blockLocation(n) /= 0) then

         this_block = get_block(n,n)
         dst_block  = dst_dist%blockLocalID(n)

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
      endif ! dst block not land
   end do  ! block loop

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
     isrc, jsrc,         &! source addresses
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     dst_block            ! local block index in dest distribution

   type (block) :: &
     this_block  ! block info for current block

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
!  copy blocks of global array into local block distribution
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      if (dst_dist%blockLocation(n) /= 0) then

         this_block = get_block(n,n)
         dst_block  = dst_dist%blockLocalID(n)

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
      endif ! dst block not land
   end do  ! block loop

!-----------------------------------------------------------------------

 end subroutine scatter_global_stress

!***********************************************************************

 end module ice_gather_scatter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
