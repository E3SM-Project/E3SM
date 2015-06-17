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

integer, parameter :: r8 = selected_real_kind(12)
integer, parameter :: i8 = selected_int_kind(12)

integer :: errcode

real(r8) :: a=0.0_r8, b=1.0_r8
integer  :: brng, method

type shr_rand_t
  character(len=32) :: type
  integer           :: nstream
  integer           :: nrandom
  integer, allocatable :: seed(:,:)
#ifdef INTEL_MKL
  type (VSL_STREAM_STATE), allocatable :: vsl_stream(:)
#endif
  type (dSFMT_t),              allocatable :: rng(:)
  type (randomNumberSequence), allocatable :: randomseq(:)
end type

public shr_genRandNum, shr_RandNum_init, shr_RandNum_term, shr_rand_t

contains

subroutine shr_genRandNum(randStream, array, nrandom_reset)

  type (shr_rand_t),        intent(inout) :: randStream
  real(r8), dimension(:,:), intent(inout) :: array
  integer, optional,         intent(in)   :: nrandom_reset

  integer :: i, n, nstream, nrandom
  integer, dimension(1) :: seed

  nstream = randStream%nstream
  nrandom = randStream%nrandom

  select case (randStream%type)

#ifdef INTEL_MKL
! intel math kernel library SIMD fast merseene twister 19937
    case ("SFMT_MKL")
    do n=1,nstream
      errcode = vdrnguniform( method, randStream%vsl_stream(n), nrandom, array(n,:), a, b)
    enddo
#endif

! keep it simple stupid
    case("KISSVEC")
    if (present(nrandom_reset)) then
      nrandom = nrandom_reset
    endif
    do i=1,nrandom
       call kissvec( randStream%seed(:,1), randStream%seed(:,2), &
            randStream%seed(:,3), randStream%seed(:,4), array(:,i), nstream )
    enddo

! fortran-95 implementation of merseene twister 19937
    case("MT19937")
    do i=1,nrandom
      do n=1,nstream
        array(n,i) = getRandomReal( randStream%randomseq(n) )
      enddo
    enddo

! fortran-90 intrinsic pseudorandom number generator
    case("F90_INTRINSIC")
    do n=1,nstream
      call random_seed( put=randStream%seed(n,:) )
      call random_number( array(n,:) )
    enddo

! SIMD-oriented fast mersenne twister 19937
  case("DSFMT_F03")
  do n=1,nstream
    call get_rand_arr_close_open( randStream%rng(n), array(n,:), nrandom )
  enddo

end select

end subroutine shr_genRandNum

integer function shr_RandNum_seed_size(type)
  character(len=*),  intent(in)    :: type

  select case (to_upper(type))
#ifdef INTEL_MKL
  case ("SFMT_MKL")
     shr_RandNum_seed_size = 1
#endif
  case("KISSVEC")
     shr_RandNum_seed_size = 4
  case("MT19937", "DSFMT_F03")
     shr_RandNum_seed_size = 1
  case("F90_INTRINSIC")
     call random_seed(size=shr_RandNum_seed_size)
  end select
end function shr_RandNum_seed_size

subroutine shr_RandNum_init( randStream, nstream, nrandom, type, seed )

  type (shr_rand_t), intent(inout) :: randStream
  integer,           intent(in)    :: nstream  ! number of streams of random numbers
  integer,           intent(in)    :: nrandom  ! number of random numbers in streams
  character(len=*),  intent(in)    :: type
  integer,           intent(in)    :: seed(:,:)

  integer :: i, n

  randStream%nstream = nstream
  randStream%nrandom = nrandom
  randStream%type    = to_upper(type)

  ! It might be useful to assert that the "seed" argument is of shape
  ! [nstream, shr_Randnum_seed_size(type)]

  select case (randStream%type)

#ifdef INTEL_MKL
! intel math kernel library SIMD fast merseene twister 19937
  case ("SFMT_MKL")
     method = VSL_RNG_METHOD_UNIFORM_STD
     brng   = VSL_BRNG_SFMT19937

     if( .NOT. allocated(randStream%vsl_stream) ) allocate(randStream%vsl_stream(nstream))
     do n=1,nstream
        errcode = vslnewstream( randStream%vsl_stream(n), brng, seed=seed(n,1) )
     enddo
#endif

! keep it simple stupid
  case("KISSVEC")

     if (allocated(randStream%seed)) deallocate(randStream%seed)
     allocate(randStream%seed(nstream,4))

     randStream%seed = seed

! fortran-95 implementation of mersenne twister 19937
  case("MT19937")
     if( .NOT. allocated(randStream%randomseq) ) allocate(randStream%randomseq(nstream))

     do n=1,nstream
        randStream%randomseq(n) = new_RandomNumberSequence( seed=seed(n,1) )
     enddo

! fortran-90 intrinsic pseudorandom number generator
  case("F90_INTRINSIC")
     if (allocated(randStream%seed)) deallocate(randStream%seed)
     allocate(randStream%seed(nstream,size(seed, 2)))

     randStream%seed = seed

! SIMD-oriented fast merseene twister 19937
  case("DSFMT_F03")
     if( .NOT. allocated(randStream%rng) ) allocate(randStream%rng(nstream))
     do n=1,nstream
        call dSFMT_init(seed(n,1), nrandom, randStream%rng(n) )
     enddo

  end select

end subroutine shr_RandNum_init

subroutine shr_RandNum_term( randStream )

type (shr_rand_t), intent(inout) :: randStream

integer :: n, nstream,ierr

  nstream = randStream%nstream

  select case (randStream%type)

#ifdef INTEL_MKL
! intel math kernel library SIMD fast merseene twister 19937
  case ("SFMT_MKL")
     if( allocated (randStream%vsl_stream) ) then
        do n=1,nstream
           ierr = vslDeleteStream(randStream%vsl_stream(n))
        enddo
        deallocate (randStream%vsl_stream)
     endif
#endif

! keep it simple stupid
  case("KISSVEC")
     if ( allocated (randStream%seed) ) deallocate (randStream%seed)

! fortran-95 implementation of merseene twister 19937
  case("MT19937")
     do n=1,nstream
        call finalize_RandomNumberSequence(randStream%randomseq(n))
     enddo
     if ( allocated(randStream%randomseq) ) deallocate (randStream%randomseq)

! fortran-90 intrinsic pseudorandom number generator
  case("F90_INTRINSIC")
     if ( allocated(randStream%seed) ) deallocate(randStream%seed)

! SIMD-oriented fast merseene twister 19937
  case("DSFMT_F03")
     do n=1,nstream
        call dSFMT_end(randStream%rng(n))
     enddo
     if ( allocated (randStream%rng) ) deallocate (randStream%rng)

  end select

end subroutine shr_RandNum_term

function to_upper(strIn) result(strOut)

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i, j

     do i = 1,len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function to_upper

end module shr_RandNum_mod
