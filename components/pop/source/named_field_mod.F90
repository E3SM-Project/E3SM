!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module named_field_mod

!BOP
! !MODULE: named_field_mod
! !DESCRIPTION:
!  Provide an interface to declare, set, and get named fields.
!
!  Fields are declared with a string name. The set and get interfaces
!  require an integer index which is returned by the register
!  subroutine. This index is also retrievable with the get_index
!  subroutine.
!
!  Generic interfaces are available for named_field_set and
!  named_field_get so that they may be
!  1) called within threaded regions and
!  2) called with a single call in non-threaded regions.
!
! !REVISION HISTORY:
!  SVN:$Id:$
!

! !USES:

   use kinds_mod
   use exit_mod, only: exit_POP, sigAbort
   use blocks, only: nx_block, ny_block
   use domain, only: nblocks_clinic

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: named_field_register   ,&
             named_field_get_index  ,&
             named_field_set        ,&
             named_field_get

!EOP
!BOC
!-----------------------------------------------------------------------
!  module private types and data
!-----------------------------------------------------------------------

   type :: named_field
      character (char_len)                 :: name
      real (r8), dimension(:,:,:), pointer :: field
   end type

   type (named_field), dimension(:), allocatable :: &
      named_field_array   ! array of registered named fields

   integer (int_kind) :: &
      named_field_cnt = 0 ! number of registered named fields

!-----------------------------------------------------------------------
!  generic interface definitions
!-----------------------------------------------------------------------

   interface named_field_set
      module procedure named_field_set_1_block, &
                       named_field_set_all_blocks
   end interface

   interface named_field_get
      module procedure named_field_get_1_block, &
                       named_field_get_all_blocks
   end interface

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: named_field_register
! !INTERFACE:

   subroutine named_field_register(name, index)

! !DESCRIPTION:
!  Register a named field.
!  It is a fatal error to register a previously registered name.
!
!  This should only be called once per task.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in)       :: &
      name              ! name of field to be registered

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      index             ! returned index of registered field

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   type (named_field), dimension(:), allocatable :: &
      named_field_array_tmp ! Tempory array to store current named_field
                            ! array contents. Used when growing the array.

!-----------------------------------------------------------------------
!  check to see if name is already registered
!-----------------------------------------------------------------------

   call named_field_get_index(name, index, exit_on_err=.false.)
   if (index > 0) then
      call exit_POP(sigAbort, &
         'named_field_register: name already registered, ' /&
                                                            &/ name)
   end if

!-----------------------------------------------------------------------
!  increase the length of the array named_field_array, preserving
!  the current contents if named_field_cnt > 0
!-----------------------------------------------------------------------

   if (named_field_cnt > 0) then
      allocate(named_field_array_tmp(named_field_cnt))
      named_field_array_tmp = named_field_array
      deallocate(named_field_array)
   endif

   allocate(named_field_array(named_field_cnt+1))

   if (named_field_cnt > 0) then
      named_field_array(1:named_field_cnt) = named_field_array_tmp
      deallocate(named_field_array_tmp)
   endif

   named_field_cnt = named_field_cnt + 1
   index = named_field_cnt
   named_field_array(named_field_cnt)%name = name
   allocate(named_field_array(named_field_cnt)%field(nx_block,ny_block,nblocks_clinic))

!-----------------------------------------------------------------------
!EOC

   end subroutine named_field_register

!***********************************************************************
!BOP
! !IROUTINE: named_field_get_index
! !INTERFACE:

   subroutine named_field_get_index(name, index, exit_on_err)

! !DESCRIPTION:
!  Search for a named field and return its index.
!  If the name is not found then exit_POP is called, unless
!     the optional argument exit_on_err is set to .false., in which
!     case index is set to 0.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in)                :: &
      name              ! name of field to be looked up

   logical (log_kind), intent(in), optional :: &
      exit_on_err       ! Is exit_POP called if name not found?

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out)          :: index

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      loc_exit_on_err   ! local copy of exit_on_err

   integer (int_kind) :: &
      i                 ! loop index

!-----------------------------------------------------------------------

   if (.not. present(exit_on_err)) then
      loc_exit_on_err = .true.
   else
      loc_exit_on_err = exit_on_err
   endif

   index = 0
   do i=1,named_field_cnt
      if (named_field_array(i)%name == name) then
         index = i
         exit
      endif
   enddo

   if (index == 0 .and. loc_exit_on_err) then
      call exit_POP(sigAbort, &
         'named_field_get_index: name not found, ' /&
                                                    &/ name)
   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine named_field_get_index

!***********************************************************************
!BOP
! !IROUTINE: named_field_set_1_block
! !INTERFACE:

   subroutine named_field_set_1_block(index, block, field)

! !DESCRIPTION:
!  Set the field component of a named field to field. This is an instance
!  of named_field_set, for use with a single block of field.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      index,           &! index into array of available named fields
      block             ! local block address (in baroclinic distribution)

   real (r8), dimension(:,:), intent(in) :: &
      field             ! array of data to be stored, single block

!EOP
!BOC
!-----------------------------------------------------------------------

   if (index < 1 .or. index > named_field_cnt) then
      call exit_POP(sigAbort, 'name_field_set: index out of range')
   end if

   named_field_array(index)%field(:,:,block) = field

!-----------------------------------------------------------------------
!EOC

   end subroutine named_field_set_1_block

!***********************************************************************
!BOP
! !IROUTINE: named_field_set_all_blocks
! !INTERFACE:

   subroutine named_field_set_all_blocks(index, field)

! !DESCRIPTION:
!  Set the field component of a named field to field. This is an instance
!  of named_field_set, for use with all blocks of field.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      index             ! index into array of available named fields

   real (r8), dimension(:,:,:), intent(in) :: &
      field             ! array of data to be stored, all blocks

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! local block number

!-----------------------------------------------------------------------

   if (index < 1 .or. index > named_field_cnt) then
      call exit_POP(sigAbort, 'name_field_set: index out of range')
   end if

   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock=1,nblocks_clinic
      named_field_array(index)%field(:,:,iblock) = field(:,:,iblock)
   end do
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

   end subroutine named_field_set_all_blocks

!***********************************************************************
!BOP
! !IROUTINE: named_field_get_1_block
! !INTERFACE:

   subroutine named_field_get_1_block(index, block, field)

! !DESCRIPTION:
!  Get the contents of a named field and store the results in field. This
!  is an instance of named_field_get, for use with a single block of field.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      index,           &! index into array of available named fields
      block             ! local block address (in baroclinic distribution)

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:), intent(out) :: &
      field             ! array where data is to be stored, single block

!EOP
!BOC
!-----------------------------------------------------------------------

   if (index < 1 .or. index > named_field_cnt) then
      call exit_POP(sigAbort, 'name_field_get: index out of range')
   end if

   field = named_field_array(index)%field(:,:,block)

!-----------------------------------------------------------------------
!EOC

   end subroutine named_field_get_1_block

!***********************************************************************
!BOP
! !IROUTINE: named_field_get_all_blocks
! !INTERFACE:

   subroutine named_field_get_all_blocks(index, field)

! !DESCRIPTION:
!  Get the contents of a named field and store the results in field. This
!  is an instance of named_field_get, for use with all blocks of field.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      index             ! index into array of available named fields

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(out) :: &
      field             ! array where data is to be stored, all blocks

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! local block number

!-----------------------------------------------------------------------

   if (index < 1 .or. index > named_field_cnt) then
      call exit_POP(sigAbort, 'name_field_get: index out of range')
   end if

   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock=1,nblocks_clinic
      field(:,:,iblock) = named_field_array(index)%field(:,:,iblock)
   end do
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

   end subroutine named_field_get_all_blocks

 end module named_field_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
