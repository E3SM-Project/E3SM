#ifdef INTEL_MKL
include 'mkl_vsl.f90'
#endif

module shr_RandNum_mod

! this module contains the driver for the available versions of random number generator

use mersennetwister_mod, only : new_RandomNumberSequence, getRandomReal, &
                                randomNumberSequence, finalize_RandomNumberSequence
use dSFMT_interface,     only : dSFMT_init, dSFMT_end, get_rand_arr_close_open, dSFMT_t
use kissvec_mod,         only : kissvec

#ifdef INTEL_MKL
use mkl_vsl_type
use mkl_vsl
#endif
implicit none

private

public :: ShrRandGen

public :: ShrIntrinsicRandGen
public :: ShrKissRandGen
public :: ShrF95MtRandGen
public :: ShrDsfmtRandGen
#ifdef INTEL_MKL
public :: ShrMklMtRandGen
#endif

integer, parameter :: r8 = selected_real_kind(12)
integer, parameter :: i8 = selected_int_kind(12)

! Abstract base class for random number generators
type, abstract :: ShrRandGen
 contains
   procedure(gen_random), deferred :: random
   procedure(gen_finalize), deferred :: finalize
end type ShrRandGen

abstract interface
   subroutine gen_random(self, array)
     import
     class(ShrRandGen), intent(inout) :: self
     real(r8), dimension(:,:), intent(out) :: array
   end subroutine gen_random

   subroutine gen_finalize(self)
     import
     class(ShrRandGen), intent(inout) :: self
   end subroutine gen_finalize
end interface

! Fortran 90 "random_number" intrinsic
type, extends(ShrRandGen) :: ShrIntrinsicRandGen
   integer, allocatable, private :: seed(:,:)
 contains
   procedure, non_overridable :: random => intrinsic_random
   procedure, non_overridable :: finalize => intrinsic_finalize
end type ShrIntrinsicRandGen

interface ShrIntrinsicRandGen
   module procedure ShrIntrinsicRandGen_constructor
end interface ShrIntrinsicRandGen

! KISS (Keep It Simple Stupid) pseudorandom number generator
type, extends(ShrRandGen) :: ShrKissRandGen
   integer, allocatable, private :: seed(:,:)
 contains
   procedure, non_overridable :: random => kiss_random
   procedure, non_overridable :: finalize => kiss_finalize
end type ShrKissRandGen

interface ShrKissRandGen
   module procedure ShrKissRandGen_constructor
end interface ShrKissRandGen

! Fortran 95 implementation of Mersene Twister 19937
type, extends(ShrRandGen) :: ShrF95MtRandGen
   type(randomNumberSequence), allocatable, private :: wrapped(:)
 contains
   procedure, non_overridable :: random => f95_mt_random
   procedure, non_overridable :: finalize => f95_mt_finalize
end type ShrF95MtRandGen

interface ShrF95MtRandGen
   module procedure ShrF95MtRandGen_constructor
end interface ShrF95MtRandGen

! Double precision SIMD Fast Mersene Twister 19937
type, extends(ShrRandGen) :: ShrDsfmtRandGen
   type(dSFMT_t), allocatable, private :: wrapped(:)
 contains
   procedure, non_overridable :: random => dsfmt_random
   procedure, non_overridable :: finalize => dsfmt_finalize
end type ShrDsfmtRandGen

interface ShrDsfmtRandGen
   module procedure ShrDsfmtRandGen_constructor
end interface ShrDsfmtRandGen

! Intel Math Kernel Library Mersenne Twister
#ifdef INTEL_MKL
type, extends(ShrRandGen) :: ShrMklMtRandGen
  type (VSL_STREAM_STATE), allocatable :: wrapped(:)
 contains
   procedure, non_overridable :: random => mkl_mt_random
   procedure, non_overridable :: finalize => mkl_mt_finalize
end type ShrMklMtRandGen

interface ShrMklMtRandGen
   module procedure ShrMklMtRandGen_constructor
end interface ShrMklMtRandGen
#endif

contains

function ShrIntrinsicRandGen_constructor(seed) result(rand_gen)
  integer, intent(in) :: seed(:,:)
  type(ShrIntrinsicRandGen) :: rand_gen

  allocate(rand_gen%seed(size(seed, 1),size(seed, 2)))
  rand_gen%seed = seed

end function ShrIntrinsicRandGen_constructor

subroutine intrinsic_random(self, array)
  class(ShrIntrinsicRandGen), intent(inout) :: self
  real(r8), dimension(:,:), intent(out) :: array

  integer :: i

  do i = 1, size(self%seed, 1)
     call random_seed(put=self%seed(i,:))
     call random_number(array(i,:))
     call random_seed(get=self%seed(i,:))
  end do

end subroutine intrinsic_random

subroutine intrinsic_finalize(self)
  class(ShrIntrinsicRandGen), intent(inout) :: self

  if ( allocated(self%seed) ) deallocate(self%seed)

end subroutine intrinsic_finalize

function ShrKissRandGen_constructor(seed) result(rand_gen)
  integer, intent(in) :: seed(:,:)
  type(ShrKissRandGen) :: rand_gen

  allocate(rand_gen%seed(size(seed, 1),4))
  rand_gen%seed = seed

end function ShrKissRandGen_constructor

subroutine kiss_random(self, array)
  class(ShrKissRandGen), intent(inout) :: self
  real(r8), dimension(:,:), intent(out) :: array

  integer :: nstream, i

  nstream = size(self%seed, 1)

  do i = 1, size(array, 2)
     call kissvec(self%seed(:,1), self%seed(:,2), self%seed(:,3), &
          self%seed(:,4), array(:,i), nstream)
  end do

end subroutine kiss_random

subroutine kiss_finalize(self)
  class(ShrKissRandGen), intent(inout) :: self

  if ( allocated(self%seed) ) deallocate(self%seed)

end subroutine kiss_finalize

function ShrF95MtRandGen_constructor(seed) result(rand_gen)
  integer, intent(in) :: seed(:)
  type(ShrF95MtRandGen) :: rand_gen

  integer :: i

  allocate(rand_gen%wrapped(size(seed)))
  do i = 1, size(seed)
     rand_gen%wrapped(i) = new_RandomNumberSequence(seed=seed(i))
  end do

end function ShrF95MtRandGen_constructor

subroutine f95_mt_random(self, array)
  class(ShrF95MtRandGen), intent(inout) :: self
  real(r8), dimension(:,:), intent(out) :: array

  integer :: nstream, i, j

  nstream = size(self%wrapped)

  do i = 1, size(array, 2)
     do j = 1, nstream
        array(j,i) = getRandomReal(self%wrapped(j))
     end do
  end do

end subroutine f95_mt_random

subroutine f95_mt_finalize(self)
  class(ShrF95MtRandGen), intent(inout) :: self

  integer :: i

  if (allocated(self%wrapped)) then
     do i = 1, size(self%wrapped)
        call finalize_RandomNumberSequence(self%wrapped(i))
     end do
     deallocate(self%wrapped)
  end if

end subroutine f95_mt_finalize

function ShrDsfmtRandGen_constructor(seed, cache_size) result(rand_gen)
  integer, intent(in) :: seed(:)
  integer, intent(in) :: cache_size
  type(ShrDsfmtRandGen) :: rand_gen

  integer :: i

  allocate(rand_gen%wrapped(size(seed)))
  do i = 1, size(seed)
     call dSFMT_init(seed(i), cache_size, rand_gen%wrapped(i))
  end do

end function ShrDsfmtRandGen_constructor

subroutine dsfmt_random(self, array)
  class(ShrDsfmtRandGen), intent(inout) :: self
  real(r8), dimension(:,:), intent(out) :: array

  integer :: nrandom, i

  nrandom = size(array, 2)

  do i = 1, size(self%wrapped)
     call get_rand_arr_close_open(self%wrapped(i), array(i,:), nrandom)
  end do

end subroutine dsfmt_random

subroutine dsfmt_finalize(self)
  class(ShrDsfmtRandGen), intent(inout) :: self

  integer :: i

  if (allocated(self%wrapped)) then
     do i = 1, size(self%wrapped)
        call dSFMT_end(self%wrapped(i))
     end do
     deallocate(self%wrapped)
  end if

end subroutine dsfmt_finalize

#ifdef INTEL_MKL
function ShrMklMtRandGen_constructor(seed) result(rand_gen)
  integer, intent(in) :: seed(:)
  type(ShrMklMtRandGen) :: rand_gen

  integer :: err

  integer :: i

  allocate(rand_gen%wrapped(size(seed)))
  do i = 1, size(seed)
     err = vslnewstream(rand_gen%wrapped(i), VSL_BRNG_SFMT19937, seed=seed(i))
  end do

end function ShrMklMtRandGen_constructor

subroutine mkl_mt_random(self, array)
  class(ShrMklMtRandGen), intent(inout) :: self
  real(r8), dimension(:,:), intent(out) :: array

  real(r8), parameter :: range_start = 0.0_r8, range_end = 1.0_r8

  integer :: err

  integer :: nrandom, i

  nrandom = size(array, 2)

  do i = 1, size(self%wrapped)
     err = vdrnguniform(VSL_RNG_METHOD_UNIFORM_STD, self%wrapped(i), &
          nrandom, array(i,:), range_start, range_end)
  end do

end subroutine mkl_mt_random

subroutine mkl_mt_finalize(self)
  class(ShrMklMtRandGen), intent(inout) :: self

  integer :: err

  integer :: i

  if (allocated(self%wrapped)) then
     do i = 1, size(self%wrapped)
        err = vslDeleteStream(self%wrapped(i))
     end do
     deallocate(self%wrapped)
  end if

end subroutine mkl_mt_finalize
#endif

end module shr_RandNum_mod
