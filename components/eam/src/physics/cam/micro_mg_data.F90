module micro_mg_data

!
! Packing and time averaging for the MG interface.
!
! Use is as follows:
!
! 1) Figure out which columns will do averaging (mgncol) and the number of
!    levels where the microphysics will run (nlev).
!
! 2) Create an MGPacker object and assign it as follows:
!
!      packer = MGPacker(pcols, pver, mgcols, top_lev)
!
!    Where [pcols, pver] is the shape of the ultimate input/output arrays
!    that are defined at level midpoints.
!
! 3) Create a post-processing array of type MGPostProc:
!
!      post_proc = MGPostProc(packer)
!
! 4) Add pairs of pointers for packed and unpacked representations, already
!    associated with buffers of the correct dimensions:
!
!      call post_proc%add_field(unpacked_pointer, packed_pointer, &
!             fillvalue, accum_mean)
!
!    The third value is the default value used to "unpack" for points with
!    no "packed" part, and the fourth value is the method used to
!    accumulate values over time steps. These two arguments can be omitted,
!    in which case the default value will be 0 and the accumulation method
!    will take the mean.
!
! 5) Use the packed fields in MG, and for each MG iteration, do:
!
!      call post_proc%accumulate()
!
! 6) Perform final accumulation and scatter values into the unpacked arrays:
!
!      call post_proc%process_and_unpack()
!
! 7) Destroy the object when complete:
!
!      call post_proc%finalize()
!
! Caveat: MGFieldPostProc will hit a divide-by-zero error if you try to
!         take the mean over 0 steps.
!

! This include header defines CPP macros that only have an effect for debug
! builds.
#include "shr_assert.h"

use shr_kind_mod, only: r8 => shr_kind_r8
use shr_log_mod, only: &
     errMsg => shr_log_errMsg, &
     OOBMsg => shr_log_OOBMsg
use shr_sys_mod, only: shr_sys_abort

implicit none
private

public :: MGPacker
public :: MGFieldPostProc
public :: accum_null
public :: accum_mean
public :: MGPostProc

type :: MGPacker
   ! Unpacked array dimensions.
   integer :: pcols
   integer :: pver
   ! Calculated packed dimensions, stored for convenience.
   integer :: mgncol
   integer :: nlev
   ! Which columns are packed.
   integer, allocatable :: mgcols(:)
   ! Topmost level to copy into the packed array.
   integer :: top_lev
 contains
   procedure, private :: pack_1D
   procedure, private :: pack_2D
   procedure, private :: pack_3D
   generic :: pack => pack_1D, pack_2D, pack_3D
   procedure :: pack_interface
   procedure, private :: unpack_1D
   procedure, private :: unpack_1D_array_fill
   procedure, private :: unpack_2D
   procedure, private :: unpack_2D_array_fill
   procedure, private :: unpack_3D
   procedure, private :: unpack_3D_array_fill
   generic :: unpack => unpack_1D, unpack_1D_array_fill, &
        unpack_2D, unpack_2D_array_fill, unpack_3D, unpack_3D_array_fill
   procedure :: finalize => MGPacker_finalize
end type MGPacker

interface MGPacker
   module procedure new_MGPacker
end interface

! Enum for time accumulation/averaging methods.
integer, parameter :: accum_null = 0
integer, parameter :: accum_mean = 1

type :: MGFieldPostProc
   integer :: accum_method = -1
   integer :: rank = -1
   integer :: num_steps = 0
   real(r8) :: fillvalue = 0._r8
   real(r8), pointer :: unpacked_1D(:) => null()
   real(r8), pointer :: packed_1D(:) => null()
   real(r8), allocatable :: buffer_1D(:)
   real(r8), pointer :: unpacked_2D(:,:) => null()
   real(r8), pointer :: packed_2D(:,:) => null()
   real(r8), allocatable :: buffer_2D(:,:)
 contains
   procedure :: accumulate => MGFieldPostProc_accumulate
   procedure :: process_and_unpack => MGFieldPostProc_process_and_unpack
   procedure :: unpack_only => MGFieldPostProc_unpack_only
   procedure :: finalize => MGFieldPostProc_finalize
end type MGFieldPostProc

interface MGFieldPostProc
   module procedure MGFieldPostProc_1D
   module procedure MGFieldPostProc_2D
end interface MGFieldPostProc

#define VECTOR_NAME MGFieldPostProcVec
#define TYPE_NAME type(MGFieldPostProc)
#define THROW(string) call shr_sys_abort(string)

public :: VECTOR_NAME

#include "dynamic_vector_typedef.inc"

type MGPostProc
   type(MGPacker) :: packer
   type(MGFieldPostProcVec) :: field_procs
 contains
   procedure, private :: add_field_1D
   procedure, private :: add_field_2D
   generic :: add_field => add_field_1D, add_field_2D
   procedure :: accumulate => MGPostProc_accumulate
   procedure :: process_and_unpack => MGPostProc_process_and_unpack
   procedure :: unpack_only => MGPostProc_unpack_only
   procedure :: finalize => MGPostProc_finalize
   procedure, private :: MGPostProc_copy
! Sarat: Following results in a segfault with PGI(17.7) compilers on Summitdev.
! This workaround is still needed for Intel compilers (v16 etc.) as of 12/19/2017.
#if !defined(SUMMITDEV_PGI)
   generic :: assignment(=) => MGPostProc_copy
#endif
end type MGPostProc

interface MGPostProc
   module procedure new_MGPostProc
end interface MGPostProc

contains

function new_MGPacker(pcols, pver, mgcols, top_lev)
  integer, intent(in) :: pcols, pver
  integer, intent(in) :: mgcols(:)
  integer, intent(in) :: top_lev

  type(MGPacker) :: new_MGPacker

  new_MGPacker%pcols = pcols
  new_MGPacker%pver = pver
  new_MGPacker%mgncol = size(mgcols)
  new_MGPacker%nlev = pver - top_lev + 1

  allocate(new_MGPacker%mgcols(new_MGPacker%mgncol))
  new_MGPacker%mgcols = mgcols
  new_MGPacker%top_lev = top_lev

end function new_MGPacker

! Rely on the fact that intent(out) forces the compiler to deallocate all
! allocatable components and restart the type from scratch. Although
! compiler support for finalization varies, this seems to be one of the few
! cases where all major compilers are reliable, and humans are not.
subroutine MGPacker_finalize(self)
  class(MGPacker), intent(out) :: self
end subroutine MGPacker_finalize

function pack_1D(self, unpacked) result(packed)
  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: unpacked(:)

  real(r8) :: packed(self%mgncol)

  SHR_ASSERT_FL(size(unpacked) == self%pcols, __FILE__, __LINE__)

  packed = unpacked(self%mgcols)

end function pack_1D

! Separation of pack and pack_interface is to workaround a PGI bug.
function pack_2D(self, unpacked) result(packed)
  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: unpacked(:,:)

  real(r8) :: packed(self%mgncol,self%nlev)

  SHR_ASSERT_FL(size(unpacked, 1) == self%pcols, __FILE__, __LINE__)

  packed = unpacked(self%mgcols,self%top_lev:)

end function pack_2D

function pack_interface(self, unpacked) result(packed)
  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: unpacked(:,:)

  real(r8) :: packed(self%mgncol,self%nlev+1)

  packed = unpacked(self%mgcols,self%top_lev:)

end function pack_interface

function pack_3D(self, unpacked) result(packed)
  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: unpacked(:,:,:)

  real(r8) :: packed(self%mgncol,self%nlev,size(unpacked, 3))

  SHR_ASSERT_FL(size(unpacked,1) == self%pcols, __FILE__, __LINE__)

  packed = unpacked(self%mgcols,self%top_lev:,:)

end function pack_3D

function unpack_1D(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:)
  real(r8), intent(in) :: fill

  real(r8) :: unpacked(self%pcols)

  SHR_ASSERT_FL(size(packed) == self%mgncol, __FILE__, __LINE__)

  unpacked = fill
  unpacked(self%mgcols) = packed

end function unpack_1D

function unpack_1D_array_fill(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:)
  real(r8), intent(in) :: fill(:)

  real(r8) :: unpacked(self%pcols)

  SHR_ASSERT_FL(size(packed) == self%mgncol, __FILE__, __LINE__)

  unpacked = fill
  unpacked(self%mgcols) = packed

end function unpack_1D_array_fill

function unpack_2D(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:,:)
  real(r8), intent(in) :: fill

  real(r8) :: unpacked(self%pcols,self%pver+size(packed, 2)-self%nlev)

  SHR_ASSERT_FL(size(packed, 1) == self%mgncol, __FILE__, __LINE__)

  unpacked = fill
  unpacked(self%mgcols,self%top_lev:) = packed

end function unpack_2D

function unpack_2D_array_fill(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:,:)
  real(r8), intent(in) :: fill(:,:)

  real(r8) :: unpacked(self%pcols,self%pver+size(packed, 2)-self%nlev)

  SHR_ASSERT_FL(size(packed, 1) == self%mgncol, __FILE__, __LINE__)

  unpacked = fill
  unpacked(self%mgcols,self%top_lev:) = packed

end function unpack_2D_array_fill

function unpack_3D(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:,:,:)
  real(r8), intent(in) :: fill

  real(r8) :: unpacked(self%pcols,self%pver,size(packed, 3))

  SHR_ASSERT_FL(size(packed, 1) == self%mgncol, __FILE__, __LINE__)

  unpacked = fill
  unpacked(self%mgcols,self%top_lev:,:) = packed

end function unpack_3D

function unpack_3D_array_fill(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:,:,:)
  real(r8), intent(in) :: fill(:,:,:)

  real(r8) :: unpacked(self%pcols,self%pver,size(packed, 3))

  SHR_ASSERT_FL(size(packed, 1) == self%mgncol, __FILE__, __LINE__)

  unpacked = fill
  unpacked(self%mgcols,self%top_lev:,:) = packed

end function unpack_3D_array_fill

function MGFieldPostProc_1D(unpacked_ptr, packed_ptr, fillvalue, &
     accum_method) result(field_proc)
  real(r8), pointer, intent(in) :: unpacked_ptr(:)
  real(r8), pointer, intent(in) :: packed_ptr(:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  type(MGFieldPostProc) :: field_proc

  field_proc%rank = 1
  field_proc%unpacked_1D => unpacked_ptr
  field_proc%packed_1D => packed_ptr
  if (present(fillvalue)) then
     field_proc%fillvalue = fillvalue
  else
     field_proc%fillvalue = 0._r8
  end if
  if (present(accum_method)) then
     field_proc%accum_method = accum_method
  else
     field_proc%accum_method = accum_mean
  end if

end function MGFieldPostProc_1D

function MGFieldPostProc_2D(unpacked_ptr, packed_ptr, fillvalue, &
     accum_method) result(field_proc)
  real(r8), pointer, intent(in) :: unpacked_ptr(:,:)
  real(r8), pointer, intent(in) :: packed_ptr(:,:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  type(MGFieldPostProc) :: field_proc

  field_proc%rank = 2
  field_proc%unpacked_2D => unpacked_ptr
  field_proc%packed_2D => packed_ptr
  if (present(fillvalue)) then
     field_proc%fillvalue = fillvalue
  else
     field_proc%fillvalue = 0._r8
  end if
  if (present(accum_method)) then
     field_proc%accum_method = accum_method
  else
     field_proc%accum_method = accum_mean
  end if

end function MGFieldPostProc_2D

! Use the same intent(out) trick as for MGPacker, which is actually more
! useful here.
subroutine MGFieldPostProc_finalize(self)
  class(MGFieldPostProc), intent(out) :: self
end subroutine MGFieldPostProc_finalize

subroutine MGFieldPostProc_accumulate(self)
  class(MGFieldPostProc), intent(inout) :: self

  select case (self%accum_method)
  case (accum_null)
     ! "Null" method does nothing.
  case (accum_mean)
     ! Allocation is done on the first accumulation step to allow the
     ! MGFieldPostProc to be copied after construction without copying the
     ! allocated array (until this function is first called).
     self%num_steps = self%num_steps + 1
     select case (self%rank)
     case (1)
        SHR_ASSERT_FL(associated(self%packed_1D), __FILE__, __LINE__)
        if (.not. allocated(self%buffer_1D)) then
           allocate(self%buffer_1D(size(self%packed_1D)))
           self%buffer_1D = 0._r8
        end if
        self%buffer_1D = self%buffer_1D + self%packed_1D
     case (2)
        SHR_ASSERT_FL(associated(self%packed_2D), __FILE__, __LINE__)
        if (.not. allocated(self%buffer_2D)) then
           ! Awkward; in F2008 can be replaced by source/mold.
           allocate(self%buffer_2D(&
                size(self%packed_2D, 1),size(self%packed_2D, 2)))
           self%buffer_2D = 0._r8
        end if
        self%buffer_2D = self%buffer_2D + self%packed_2D
     case default
        call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
             " Unsupported rank for MGFieldPostProc accumulation.")
     end select
  case default
     call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
          " Unrecognized MGFieldPostProc accumulation method.")
  end select

end subroutine MGFieldPostProc_accumulate

subroutine MGFieldPostProc_process_and_unpack(self, packer)
  class(MGFieldPostProc), intent(inout) :: self
  class(MGPacker), intent(in) :: packer

  select case (self%accum_method)
  case (accum_null)
     ! "Null" method just leaves the value as the last time step, so don't
     ! actually need to do anything.
  case (accum_mean)
     select case (self%rank)
     case (1)
        SHR_ASSERT_FL(associated(self%packed_1D), __FILE__, __LINE__)
        self%packed_1D = self%buffer_1D/self%num_steps
     case (2)
        SHR_ASSERT_FL(associated(self%packed_2D), __FILE__, __LINE__)
        self%packed_2D = self%buffer_2D/self%num_steps
     case default
        call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
             " Unsupported rank for MGFieldPostProc accumulation.")
     end select
  case default
     call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
          " Unrecognized MGFieldPostProc accumulation method.")
  end select

  call self%unpack_only(packer)

end subroutine MGFieldPostProc_process_and_unpack

subroutine MGFieldPostProc_unpack_only(self, packer)
  class(MGFieldPostProc), intent(inout) :: self
  class(MGPacker), intent(in) :: packer

  select case (self%rank)
  case (1)
     SHR_ASSERT_FL(associated(self%unpacked_1D), __FILE__, __LINE__)
     self%unpacked_1D = packer%unpack(self%packed_1D, self%fillvalue)
  case (2)
     SHR_ASSERT_FL(associated(self%unpacked_2D), __FILE__, __LINE__)
     self%unpacked_2D = packer%unpack(self%packed_2D, self%fillvalue)
  case default
     call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
          " Unsupported rank for MGFieldPostProc unpacking.")
  end select

end subroutine MGFieldPostProc_unpack_only

#include "dynamic_vector_procdef.inc"

function new_MGPostProc(packer) result(post_proc)
  type(MGPacker), intent(in) :: packer

  type(MGPostProc) :: post_proc

  post_proc%packer = packer
  call post_proc%field_procs%clear()

end function new_MGPostProc

! Can't use the same intent(out) trick, because PGI doesn't get the
! recursive deallocation right.
subroutine MGPostProc_finalize(self)
  class(MGPostProc), intent(inout) :: self

  integer :: i

  call self%packer%finalize()
  do i = 1, self%field_procs%vsize()
     call self%field_procs%data(i)%finalize()
  end do
  call self%field_procs%clear()
  call self%field_procs%shrink_to_fit()

end subroutine MGPostProc_finalize

subroutine add_field_1D(self, unpacked_ptr, packed_ptr, fillvalue, &
     accum_method)
  class(MGPostProc), intent(inout) :: self
  real(r8), pointer, intent(in) :: unpacked_ptr(:)
  real(r8), pointer, intent(in) :: packed_ptr(:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method

  call self%field_procs%push_back(MGFieldPostProc(unpacked_ptr, &
       packed_ptr, fillvalue, accum_method))

end subroutine add_field_1D

subroutine add_field_2D(self, unpacked_ptr, packed_ptr, fillvalue, &
     accum_method)
  class(MGPostProc), intent(inout) :: self
  real(r8), pointer, intent(in) :: unpacked_ptr(:,:)
  real(r8), pointer, intent(in) :: packed_ptr(:,:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method

  call self%field_procs%push_back(MGFieldPostProc(unpacked_ptr, &
       packed_ptr, fillvalue, accum_method))

end subroutine add_field_2D

subroutine MGPostProc_accumulate(self)
  class(MGPostProc), intent(inout) :: self

  integer :: i

  do i = 1, self%field_procs%vsize()
     call self%field_procs%data(i)%accumulate()
  end do

end subroutine MGPostProc_accumulate

subroutine MGPostProc_process_and_unpack(self)
  class(MGPostProc), intent(inout) :: self

  integer :: i

  do i = 1, self%field_procs%vsize()
     call self%field_procs%data(i)%process_and_unpack(self%packer)
  end do

end subroutine MGPostProc_process_and_unpack

subroutine MGPostProc_unpack_only(self)
  class(MGPostProc), intent(inout) :: self

  integer :: i

  do i = 1, self%field_procs%vsize()
     call self%field_procs%data(i)%unpack_only(self%packer)
  end do

end subroutine MGPostProc_unpack_only

! This is necessary only to work around Intel/PGI bugs.
subroutine MGPostProc_copy(lhs, rhs)
  class(MGPostProc), intent(out) :: lhs
  type(MGPostProc), intent(in) :: rhs

  lhs%packer = rhs%packer
  lhs%field_procs = rhs%field_procs
end subroutine MGPostProc_copy

end module micro_mg_data
