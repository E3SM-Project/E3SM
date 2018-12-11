! Fortran-95 implementation of the Mersenne Twister 19937, following
!   the C implementation described below (code mt19937ar-cok.c, dated 2002/2/10),
!   adapted cosmetically by making the names more general.
! Users must declare one or more variables of type randomNumberSequence in the calling
!   procedure which are then initialized using a required seed. If the
!   variable is not initialized the random numbers will all be 0.
! For example:
! program testRandoms
!   use RandomNumbers
!   type(randomNumberSequence) :: randomNumbers
!   integer                    :: i
!
!   randomNumbers = new_RandomNumberSequence(seed = 100)
!   do i = 1, 10
!     print ('(f12.10, 2x)'), getRandomReal(randomNumbers)
!   end do
! end program testRandoms
!
! Fortran-95 implementation by
!   Robert Pincus
!   NOAA-CIRES Climate Diagnostics Center
!   Boulder, CO 80305
!   email: Robert.Pincus@colorado.edu
!
! This documentation in the original C program reads:
! -------------------------------------------------------------
!    A C-program for MT19937, with initialization improved 2002/2/10.
!    Coded by Takuji Nishimura and Makoto Matsumoto.
!    This is a faster version by taking Shawn Cokus's optimization,
!    Matthe Bellew's simplification, Isaku Wada's real version.
!
!    Before using, initialize the state by using init_genrand(seed)
!    or init_by_array(init_key, key_length).
!
!    Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
!    All rights reserved.
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!
!      3. The names of its contributors may not be used to endorse or promote
!         products derived from this software without specific prior written
!         permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!    A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
!    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!
!    Any feedback is very welcome.
!    http://www.math.keio.ac.jp/matumoto/emt.html
!    email: matumoto@math.keio.ac.jp
! -------------------------------------------------------------

module MersenneTwister_mod
! -------------------------------------------------------------
  implicit none
  private

  ! Algorithm parameters
  ! -------
  ! Period parameters
  integer, parameter :: i8 = selected_int_kind(13)
  integer, parameter :: blockSize = 624,         &
                        M         = 397,         &
                        MATRIX_A  = -1727483681, & ! constant vector a         (0x9908b0dfUL)
                        UMASK     = -2147483648_i8,& ! most significant w-r bits (0x80000000UL)
                        LMASK     =  2147483647    ! least significant r bits  (0x7fffffffUL)
  ! Tempering parameters
  integer, parameter :: TMASKB= -1658038656, & ! (0x9d2c5680UL)
                        TMASKC= -272236544     ! (0xefc60000UL)
  ! -------

  ! The type containing the state variable
  type randomNumberSequence
    integer                            :: currentElement ! = blockSize
    integer, dimension(0:blockSize -1) :: state ! = 0
  end type randomNumberSequence

  interface new_RandomNumberSequence
    module procedure initialize_scalar, initialize_vector
  end interface new_RandomNumberSequence

  public :: randomNumberSequence
  public :: new_RandomNumberSequence, finalize_RandomNumberSequence, &
            getRandomInt, getRandomPositiveInt, getRandomReal
! -------------------------------------------------------------
contains
  ! -------------------------------------------------------------
  ! Private functions
  ! ---------------------------
  function mixbits(u, v)
    integer, intent( in) :: u, v
    integer              :: mixbits

    mixbits = ior(iand(u, UMASK), iand(v, LMASK))
  end function mixbits
  ! ---------------------------
  function twist(u, v)
    integer, intent( in) :: u, v
    integer              :: twist

    ! Local variable
    integer, parameter, dimension(0:1) :: t_matrix = (/ 0, MATRIX_A /)

    twist = ieor(ishft(mixbits(u, v), -1), t_matrix(iand(v, 1)))
    twist = ieor(ishft(mixbits(u, v), -1), t_matrix(iand(v, 1)))
  end function twist
  ! ---------------------------
  subroutine nextState(twister)
    type(randomNumberSequence), intent(inout) :: twister

    ! Local variables
    integer :: k

    do k = 0, blockSize - M - 1
      twister%state(k) = ieor(twister%state(k + M), &
                              twist(twister%state(k), twister%state(k + 1)))
    end do
    do k = blockSize - M, blockSize - 2
      twister%state(k) = ieor(twister%state(k + M - blockSize), &
                              twist(twister%state(k), twister%state(k + 1)))
    end do
    twister%state(blockSize - 1) = ieor(twister%state(M - 1), &
                                        twist(twister%state(blockSize - 1), twister%state(0)))
    twister%currentElement = 0

  end subroutine nextState
  ! ---------------------------
  elemental function temper(y)
    integer, intent(in) :: y
    integer             :: temper

    integer :: x

    ! Tempering
    x      = ieor(y, ishft(y, -11))
    x      = ieor(x, iand(ishft(x,  7), TMASKB))
    x      = ieor(x, iand(ishft(x, 15), TMASKC))
    temper = ieor(x, ishft(x, -18))
  end function temper
  ! -------------------------------------------------------------
  ! Public (but hidden) functions
  ! --------------------
  function initialize_scalar(seed) result(twister)
    integer,       intent(in   ) :: seed
    type(randomNumberSequence)                :: twister

    integer :: i
    ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. In the previous versions,
    !   MSBs of the seed affect only MSBs of the array state[].
    !   2002/01/09 modified by Makoto Matsumoto

    twister%state(0) = iand(seed, -1)
    do i = 1,  blockSize - 1 ! ubound(twister%state)
       twister%state(i) = 1812433253 * ieor(twister%state(i-1), &
                                            ishft(twister%state(i-1), -30)) + i
       twister%state(i) = iand(twister%state(i), -1) ! for >32 bit machines
    end do
    twister%currentElement = blockSize
  end function initialize_scalar
  ! -------------------------------------------------------------
  function initialize_vector(seed) result(twister)
    integer, dimension(0:), intent(in) :: seed
    type(randomNumberSequence)                      :: twister

    integer :: i, j, k, nFirstLoop, nWraps

    nWraps  = 0
    twister = initialize_scalar(19650218)

    nFirstLoop = max(blockSize, size(seed))
    do k = 1, nFirstLoop
       i = mod(k + nWraps, blockSize)
       j = mod(k - 1,      size(seed))
       if(i == 0) then
         twister%state(i) = twister%state(blockSize - 1)
         twister%state(1) = ieor(twister%state(1),                                 &
                                 ieor(twister%state(1-1),                          &
                                      ishft(twister%state(1-1), -30)) * 1664525) + &
                            seed(j) + j ! Non-linear
         twister%state(i) = iand(twister%state(i), -1) ! for >32 bit machines
         nWraps = nWraps + 1
       else
         twister%state(i) = ieor(twister%state(i),                                 &
                                 ieor(twister%state(i-1),                          &
                                      ishft(twister%state(i-1), -30)) * 1664525) + &
                            seed(j) + j ! Non-linear
         twister%state(i) = iand(twister%state(i), -1) ! for >32 bit machines
      end if
    end do

    !
    ! Walk through the state array, beginning where we left off in the block above
    !
    do i = mod(nFirstLoop, blockSize) + nWraps + 1, blockSize - 1
      twister%state(i) = ieor(twister%state(i),                                 &
                              ieor(twister%state(i-1),                          &
                                   ishft(twister%state(i-1), -30)) * 1566083941) - i ! Non-linear
      twister%state(i) = iand(twister%state(i), -1) ! for >32 bit machines
    end do

    twister%state(0) = twister%state(blockSize - 1)

    do i = 1, mod(nFirstLoop, blockSize) + nWraps
      twister%state(i) = ieor(twister%state(i),                                 &
                              ieor(twister%state(i-1),                          &
                                   ishft(twister%state(i-1), -30)) * 1566083941) - i ! Non-linear
      twister%state(i) = iand(twister%state(i), -1) ! for >32 bit machines
    end do

    twister%state(0) = UMASK
    twister%currentElement = blockSize

  end function initialize_vector
  ! -------------------------------------------------------------
  ! Public functions
  ! --------------------
  function getRandomInt(twister)
    type(randomNumberSequence), intent(inout) :: twister
    integer                      :: getRandomInt
    ! Generate a random integer on the interval [0,0xffffffff]
    !   Equivalent to genrand_int32 in the C code.
    !   Fortran doesn't have a type that's unsigned like C does,
    !   so this is integers in the range -2**31 - 2**31
    ! All functions for getting random numbers call this one,
    !   then manipulate the result

    if(twister%currentElement >= blockSize) call nextState(twister)

    getRandomInt = temper(twister%state(twister%currentElement))
    twister%currentElement = twister%currentElement + 1

  end function getRandomInt
  ! --------------------
  function getRandomPositiveInt(twister)
    type(randomNumberSequence), intent(inout) :: twister
    integer                      :: getRandomPositiveInt
    ! Generate a random integer on the interval [0,0x7fffffff]
    !   or [0,2**31]
    !   Equivalent to genrand_int31 in the C code.

    ! Local integers
    integer :: localInt

    localInt = getRandomInt(twister)
    getRandomPositiveInt = ishft(localInt, -1)

  end function getRandomPositiveInt
  ! --------------------
  function getRandomReal(twister)

    type(randomNumberSequence), intent(inout) :: twister
    double precision             :: getRandomReal
    ! Generate a random number on [0,1]
    !   Equivalent to genrand_real1 in the C code
    !   The result is stored as double precision but has 32 bit resolution

    integer :: localInt

    localInt = getRandomInt(twister)
    if(localInt < 0) then
      getRandomReal = dble(localInt + 2.0d0**32)/(2.0d0**32 - 1.0d0)
    else
      getRandomReal = dble(localInt            )/(2.0d0**32 - 1.0d0)
    end if
  end function getRandomReal
  ! --------------------
  subroutine finalize_RandomNumberSequence(twister)
    type(randomNumberSequence), intent(inout) :: twister

      twister%currentElement = blockSize
      twister%state(:) = 0
  end subroutine finalize_RandomNumberSequence
  ! --------------------
end module MersenneTwister_mod

