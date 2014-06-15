! A C-program for MT19937, with initialization improved 2002/1/26.
! Coded by Takuji Nishimura and Makoto Matsumoto.                 

! Code converted to Fortran 95 by José Rui Faustino de Sousa
! Date: 2002-02-01

! Enhanced version by José Rui Faustino de Sousa
! Date: 2003-04-30

! Interface:
!
! Kinds:
!   genrand_intg
!     Integer kind used must be at least 32 bits.
!   genrand_real
!     Real kind used 
!
! Types:
!   genrand_state
!     Internal representation of the RNG state.
!   genrand_srepr
!     Public representation of the RNG state. Should be used to save the RNG state.
!
! Procedures:
!   assignment(=)
!     Converts from type genrand_state to genrand_srepr and vice versa.
!   genrand_init
!     Internal RNG state initialization subroutine accepts either an genrand_intg integer
!     or a vector as seed or a new state using "put=" returns the present state using
!     "get=". If it is called with "get=" before being seeded with "put=" returns a state
!     initialized with a default seed.
!   genrand_int32
!     Subroutine returns an array or scalar whose elements are random integer on the
!     [0,0xffffffff] interval.
!   genrand_int31
!     Subroutine returns an array or scalar whose elements are random integer on the
!     [0,0x7fffffff] interval.
!   genrand_real1
!     Subroutine returns an array or scalar whose elements are random real on the
!     [0,1] interval.
!   genrand_real2
!     Subroutine returns an array or scalar whose elements are random real on the
!     [0,1[ interval.
!   genrand_real3
!     Subroutine returns an array or scalar whose elements are random real on the
!     ]0,1[ interval.
!   genrand_res53
!     Subroutine returns an array or scalar whose elements are random real on the
!     [0,1[ interval with 53-bit resolution.

! Before using, initialize the state by using genrand_init( put=seed )  

! This library is free software.                                  
! This library is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of  
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            

! Copyright (C) 1997, 2002 Makoto Matsumoto and Takuji Nishimura. 
! Any feedback is very welcome.                                   
! http://www.math.keio.ac.jp/matumoto/emt.html                    
! email: matumoto@math.keio.ac.jp                                 
module mt95
  
  implicit none

  public  :: genrand_init, assignment(=)
  public  :: genrand_int32, genrand_int31, genrand_real1
  public  :: genrand_real2, genrand_real3, genrand_res53
  private :: uiadd, uisub, uimlt, uidiv, uimod
  private :: init_by_type, init_by_scalar, init_by_array, next_state
  private :: genrand_encode, genrand_decode, genrand_load_state, genrand_dump_state
  private :: genrand_int32_0d, genrand_int32_1d, genrand_int32_2d, genrand_int32_3d
  private :: genrand_int32_4d, genrand_int32_5d, genrand_int32_6d, genrand_int32_7d
  private :: genrand_int31_0d, genrand_int31_1d, genrand_int31_2d, genrand_int31_3d
  private :: genrand_int31_4d, genrand_int31_5d, genrand_int31_6d, genrand_int31_7d
  private :: genrand_real1_0d, genrand_real1_1d, genrand_real1_2d, genrand_real1_3d
  private :: genrand_real1_4d, genrand_real1_5d, genrand_real1_6d, genrand_real1_7d
  private :: genrand_real2_0d, genrand_real2_1d, genrand_real2_2d, genrand_real2_3d
  private :: genrand_real2_4d, genrand_real2_5d, genrand_real2_6d, genrand_real2_7d
  private :: genrand_real3_0d, genrand_real3_1d, genrand_real3_2d, genrand_real3_3d
  private :: genrand_real3_4d, genrand_real3_5d, genrand_real3_6d, genrand_real3_7d
  private :: genrand_res53_0d, genrand_res53_1d, genrand_res53_2d, genrand_res53_3d
  private :: genrand_res53_4d, genrand_res53_5d, genrand_res53_6d, genrand_res53_7d

  intrinsic :: selected_int_kind, selected_real_kind

  integer, public, parameter  :: genrand_intg = selected_int_kind( 9 )
  integer, public, parameter  :: genrand_real = selected_real_kind( 15 )

  integer, private, parameter :: wi = genrand_intg
  integer, private, parameter :: wr = genrand_real

  ! Period parameters   
  integer(kind=wi), private, parameter :: n = 624_wi
  integer(kind=wi), private, parameter :: m = 397_wi

  integer(kind=wi), private, parameter :: default_seed = 5489_wi

  integer(kind=wi), private, parameter :: fbs = 32_wi
  integer(kind=wi), private, parameter :: hbs = fbs / 2_wi
  integer(kind=wi), private, parameter :: qbs = hbs / 2_wi
  integer(kind=wi), private, parameter :: tbs = 3_wi * qbs

  real(kind=wr), private, parameter :: p231       = 2147483648.0_wr
  real(kind=wr), private, parameter :: p232       = 4294967296.0_wr
  real(kind=wr), private, parameter :: p232_1     = p232 - 1.0_wr
  real(kind=wr), private, parameter :: pi232      = 1.0_wr / p232
  real(kind=wr), private, parameter :: pi232_1    = 1.0_wr / p232_1
  real(kind=wr), private, parameter :: pi227      = 1.0_wr / 134217728.0_wr
  real(kind=wr), private, parameter :: pi253      = 1.0_wr / 9007199254740992.0_wr
  real(kind=wr), private, parameter :: p231d232_1 = p231 / p232_1
  real(kind=wr), private, parameter :: p231_5d232 = ( p231 + 0.5_wr ) / p232

  character(len=*), private, parameter  :: alph = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
  character(len=*), private, parameter  :: sepr = "&"
  integer(kind=wi), private, parameter  :: alps = 62_wi
  integer(kind=wi), private, parameter  :: clen = ( n + 1_wi ) * 7_wi !n * ( ceiling( fbs * log( 2.0_core_rknd ) / log( alps ) ) + 1 )

  type, public :: genrand_state
    private
    logical(kind=wi)                :: ini = .false._wi
    integer(kind=wi)                :: cnt = n+1_wi
    integer(kind=wi), dimension(n)  :: val = 0_wi
  end type genrand_state

  type, public :: genrand_srepr
    character(len=clen) :: repr
  end type genrand_srepr

  type(genrand_state), private, save  :: state

  interface assignment( = )
    module procedure genrand_load_state
    module procedure genrand_dump_state
  end interface assignment( = )

  interface genrand_init
    module procedure init_by_type
    module procedure init_by_scalar
    module procedure init_by_array
  end interface genrand_init

  interface genrand_int32
    module procedure genrand_int32_0d
    module procedure genrand_int32_1d
    module procedure genrand_int32_2d
    module procedure genrand_int32_3d
    module procedure genrand_int32_4d
    module procedure genrand_int32_5d
    module procedure genrand_int32_6d
    module procedure genrand_int32_7d
  end interface genrand_int32

  interface genrand_int31
    module procedure genrand_int31_0d
    module procedure genrand_int31_1d
    module procedure genrand_int31_2d
    module procedure genrand_int31_3d
    module procedure genrand_int31_4d
    module procedure genrand_int31_5d
    module procedure genrand_int31_6d
    module procedure genrand_int31_7d
  end interface genrand_int31

  interface genrand_real1
    module procedure genrand_real1_0d
    module procedure genrand_real1_1d
    module procedure genrand_real1_2d
    module procedure genrand_real1_3d
    module procedure genrand_real1_4d
    module procedure genrand_real1_5d
    module procedure genrand_real1_6d
    module procedure genrand_real1_7d
  end interface genrand_real1

  interface genrand_real2
    module procedure genrand_real2_0d
    module procedure genrand_real2_1d
    module procedure genrand_real2_2d
    module procedure genrand_real2_3d
    module procedure genrand_real2_4d
    module procedure genrand_real2_5d
    module procedure genrand_real2_6d
    module procedure genrand_real2_7d
  end interface genrand_real2

  interface genrand_real3
    module procedure genrand_real3_0d
    module procedure genrand_real3_1d
    module procedure genrand_real3_2d
    module procedure genrand_real3_3d
    module procedure genrand_real3_4d
    module procedure genrand_real3_5d
    module procedure genrand_real3_6d
    module procedure genrand_real3_7d
  end interface genrand_real3

  interface genrand_res53
    module procedure genrand_res53_0d
    module procedure genrand_res53_1d
    module procedure genrand_res53_2d
    module procedure genrand_res53_3d
    module procedure genrand_res53_4d
    module procedure genrand_res53_5d
    module procedure genrand_res53_6d
    module procedure genrand_res53_7d
  end interface genrand_res53

  contains

  elemental function uiadd( a, b ) result( c )

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
    s1 = a1 + b1
    s2 = a2 + b2 + ibits( s1, hbs, hbs )
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )
    return

  end function uiadd
  
  elemental function uisub( a, b ) result( c )

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
    s1 = a1 - b1
    s2 = a2 - b2 + ibits( s1, hbs, hbs )
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )
    return

  end function uisub
  
  elemental function uimlt( a, b ) result( c )

    intrinsic :: ibits, ior, ishft

    integer(kind=wi), intent(in)  :: a, b

    integer(kind=wi)  :: c

    integer(kind=wi)  :: a0, a1, a2, a3
    integer(kind=wi)  :: b0, b1, b2, b3
    integer(kind=wi)  :: p0, p1, p2, p3

    a0 = ibits( a, 0, qbs )
    a1 = ibits( a, qbs, qbs )
    a2 = ibits( a, hbs, qbs )
    a3 = ibits( a, tbs, qbs )
    b0 = ibits( b, 0, qbs )
    b1 = ibits( b, qbs, qbs )
    b2 = ibits( b, hbs, qbs )
    b3 = ibits( b, tbs, qbs )
    p0 = a0 * b0
    p1 = a1 * b0 + a0 * b1 + ibits( p0, qbs, tbs )
    p2 = a2 * b0 + a1 * b1 + a0 * b2 + ibits( p1, qbs, tbs )
    p3 = a3 * b0 + a2 * b1 + a1 * b2 + a0 * b3 + ibits( p2, qbs, tbs )
    c  = ior( ishft( p1, qbs ), ibits( p0, 0, qbs ) )
    c  = ior( ishft( p2, hbs ), ibits( c, 0, hbs ) )
    c  = ior( ishft( p3, tbs ), ibits( c, 0, tbs ) )
    return

  end function uimlt

  elemental function uidiv( a, b ) result( c )
    
    intrinsic :: btest, ishft

    integer(kind=wi), intent(in)  :: a, b

    integer(kind=wi)  :: c

    integer(kind=wi)  :: dl, rl

    if ( btest( a, fbs-1 ) ) then
      if ( btest( b, fbs-1 ) ) then
        if ( a < b ) then
          c = 0
        else
          c = 1
        end if
      else
        dl = ishft( ishft( a, -1 ) / b, 1 )
        rl = uisub( a, uimlt( b, dl ) )
        if ( rl < b ) then
          c = dl
        else
          c = uiadd( dl, 1 )
        end if
      end if
    else
      if ( btest( b, fbs-1 ) ) then
        c = 0
      else
        c = a / b
      end if
    end if
    return

  end function uidiv

  elemental function uimod( a, b ) result( c )
    
    intrinsic :: modulo, btest, ishft

    integer(kind=wi), intent(in)  :: a, b

    integer(kind=wi)  :: c

    integer(kind=wi)  :: dl, rl

    if ( btest( a, fbs-1 ) ) then
      if ( btest( b, fbs-1 ) ) then
        if ( a < b ) then
          c = a
        else
          c = uisub( a, b )
        end if
      else
        dl = ishft( ishft( a, -1 ) / b, 1 )
        rl = uisub( a, uimlt( b, dl ) )
        if ( rl < b ) then
          c = rl
        else
          c = uisub( rl, b )
        end if
      end if
    else
      if ( btest( b, fbs-1 ) ) then
        c = a
      else
        c = modulo( a, b )
      end if
    end if
    return

  end function uimod

  subroutine init_by_type( put, get )

    intrinsic :: present

    type(genrand_state), optional, intent(in ) :: put
    type(genrand_state), optional, intent(out) :: get

    if ( present( put ) ) then
      if ( put%ini ) state = put
    else if ( present( get ) ) then
      if ( .not. state%ini ) call init_by_scalar( default_seed )
      get = state
    else
      call init_by_scalar( default_seed )
    end if
    return

  end subroutine init_by_type

  ! initializes mt[N] with a seed
  subroutine init_by_scalar( put )

    intrinsic :: ishft, ieor, ibits

    integer(kind=wi), parameter :: mult_a = 1812433253_wi !z'6C078965'

    integer(kind=wi), intent(in)  :: put

    integer(kind=wi)  :: i

    state%ini = .true._wi
    state%val(1) = ibits( put, 0, fbs )
    do i = 2, n, 1
      state%val(i) = ieor( state%val(i-1), ishft( state%val(i-1), -30 ) )
      state%val(i) = uimlt( state%val(i), mult_a )
      state%val(i) = uiadd( state%val(i), i-1_wi )
      ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. 
      ! In the previous versions, MSBs of the seed affect   
      ! only MSBs of the array mt[].                        
      ! 2002/01/09 modified by Makoto Matsumoto             
      state%val(i) = ibits( state%val(i), 0, fbs )
      ! for >32 bit machines 
    end do
    state%cnt = n + 1_wi
    return

  end subroutine init_by_scalar

  ! initialize by an array with array-length 
  ! init_key is the array for initializing keys 
  ! key_length is its length
  subroutine init_by_array( put )

    intrinsic :: size, max, ishft, ieor, ibits
    
    integer(kind=wi), parameter :: seed_d =    19650218_wi !z'12BD6AA'
    integer(kind=wi), parameter :: mult_a =     1664525_wi !z'19660D'
    integer(kind=wi), parameter :: mult_b =  1566083941_wi !z'5D588B65'
    integer(kind=wi), parameter :: msb1_d = ishft( 1_wi, fbs-1 ) !z'80000000'

    integer(kind=wi), dimension(:), intent(in)  :: put

    integer(kind=wi)  :: i, j, k, tp, key_length

    call init_by_scalar( seed_d )
    key_length = size( put, dim=1 )
    i = 2_wi
    j = 1_wi
    do k = max( n, key_length ), 1, -1
      tp = ieor( state%val(i-1), ishft( state%val(i-1), -30 ) )
      tp = uimlt( tp, mult_a )
      state%val(i) = ieor( state%val(i), tp )
      state%val(i) = uiadd( state%val(i), uiadd( put(j), j-1_wi ) ) ! non linear 
      state%val(i) = ibits( state%val(i), 0, fbs ) ! for WORDSIZE > 32 machines 
      i = i + 1_wi
      j = j + 1_wi
      if ( i > n ) then
        state%val(1) = state%val(n)
        i = 2_wi
      end if
      if ( j > key_length) j = 1_wi
    end do
    do k = n-1, 1, -1
      tp = ieor( state%val(i-1), ishft( state%val(i-1), -30 ) )
      tp = uimlt( tp, mult_b )
      state%val(i) = ieor( state%val(i), tp )
      state%val(i) = uisub( state%val(i), i-1_wi ) ! non linear 
      state%val(i) = ibits( state%val(i), 0, fbs ) ! for WORDSIZE > 32 machines 
      i = i + 1_wi
      if ( i > n ) then
        state%val(1) = state%val(n)
        i = 2_wi
      end if
    end do
    state%val(1) = msb1_d ! MSB is 1; assuring non-zero initial array
    return

  end subroutine init_by_array

  subroutine next_state( )

    intrinsic :: ishft, ieor, btest, ibits, mvbits

    integer(kind=wi), parameter :: matrix_a = -1727483681_wi !z'9908b0df'

    integer(kind=wi)  :: i, mld

    if ( .not. state%ini ) call init_by_scalar( default_seed )
    do i = 1, n-m, 1
      mld = ibits( state%val(i+1), 0, 31 )
      call mvbits( state%val(i), 31, 1, mld, 31 )
      state%val(i) = ieor( state%val(i+m), ishft( mld, -1 ) )
      if ( btest( state%val(i+1), 0 ) ) state%val(i) = ieor( state%val(i), matrix_a )
    end do
    do i = n-m+1, n-1, 1
      mld = ibits( state%val(i+1), 0, 31 )
      call mvbits( state%val(i), 31, 1, mld, 31 )
      state%val(i) = ieor( state%val(i+m-n), ishft( mld, -1 ) )
      if ( btest( state%val(i+1), 0 ) ) state%val(i) = ieor( state%val(i), matrix_a )
    end do
    mld = ibits( state%val(1), 0, 31 )
    call mvbits( state%val(n), 31, 1, mld, 31 )
    state%val(n) = ieor( state%val(m), ishft( mld, -1 ) )
    if ( btest( state%val(1), 0 ) ) state%val(n) = ieor( state%val(n), matrix_a )
    state%cnt = 1_wi
    return

  end subroutine next_state

  elemental subroutine genrand_encode( chr, val )
    
    intrinsic :: len

    character(len=*), intent(out) :: chr
    integer(kind=wi), intent(in ) :: val
    
    integer(kind=wi)  :: i, m, d

    d = val
    chr = ""
    do i = 1, len( chr ), 1
      m = uimod( d, alps ) + 1
      chr(i:i) = alph(m:m)
      d = uidiv( d, alps )
      if ( d == 0 ) exit
    end do
    return

  end subroutine genrand_encode

  elemental subroutine genrand_decode( val, chr )
    
    intrinsic :: len, len_trim, trim, adjustl, scan

    integer(kind=wi), intent(out) :: val
    character(len=*), intent(in ) :: chr
    
    integer(kind=wi)        :: i, e, p
    character(len=len(chr)) :: c

    e = 1
    c = trim( adjustl( chr ) )
    val = 0
    do i = 1, len_trim( c ), 1
      p = scan( alph, c(i:i) ) - 1
      if( p >= 0 ) then
        val = uiadd( val, uimlt( p, e ) )
        e = uimlt( e, alps )
      end if
    end do
    return

  end subroutine genrand_decode

  elemental subroutine genrand_load_state( stt, rpr )

    intrinsic :: scan

    type(genrand_state), intent(out)  :: stt
    type(genrand_srepr), intent(in )  :: rpr

    integer(kind=wi)    :: i, j
    character(len=clen) :: c

    i = 1
    c = rpr%repr
    do
      j = scan( c, sepr )
      if ( j /= 0 ) then
        call genrand_decode( stt%val(i), c(:j-1) )
        i = i + 1
        c = c(j+1:)
      else
        exit
      end if
    end do
    call genrand_decode( stt%cnt, c )
    stt%ini = .true._wi
    return

  end subroutine genrand_load_state

  elemental subroutine genrand_dump_state( rpr, stt )

    intrinsic :: len_trim

    type(genrand_srepr), intent(out) :: rpr
    type(genrand_state), intent(in ) :: stt

    integer(kind=wi)  :: i, j

    j = 1
    rpr%repr = ""
    do i = 1, n, 1
      call genrand_encode( rpr%repr(j:), stt%val(i) )
      j = len_trim( rpr%repr ) + 1
      rpr%repr(j:j) = sepr
      j = j + 1
    end do
    call genrand_encode( rpr%repr(j:), stt%cnt )
    return

  end subroutine genrand_dump_state

  ! generates a random number on [0,0xffffffff]-interval
  subroutine genrand_int32_0d( y )

    intrinsic :: ieor, iand, ishft

    integer(kind=wi), parameter :: temper_a = -1658038656_wi !z'9D2C5680'
    integer(kind=wi), parameter :: temper_b =  -272236544_wi !z'EFC60000'

    integer(kind=wi), intent(out) :: y
    
    if ( state%cnt > n ) call next_state( )
    y = state%val(state%cnt)
    state%cnt = state%cnt + 1_wi
    ! Tempering 
    y = ieor( y, ishft( y, -11 ) )
    y = ieor( y, iand( ishft( y,  7 ), temper_a ) )
    y = ieor( y, iand( ishft( y, 15 ), temper_b ) )
    y = ieor( y, ishft( y, -18 ) )
    return

  end subroutine genrand_int32_0d

  subroutine genrand_int32_1d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:), intent(out) :: y

    integer(kind=wi)  :: i

    do i = 1, size( y, 1 ), 1
      call genrand_int32_0d( y(i) )
    end do
    return

  end subroutine genrand_int32_1d

  subroutine genrand_int32_2d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:), intent(out) :: y
    
    integer(kind=wi)  :: i

    do i = 1, size( y, 2 ), 1
      call genrand_int32_1d( y(:,i) )
    end do
    return

  end subroutine genrand_int32_2d

  subroutine genrand_int32_3d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:,:), intent(out) :: y
    
    integer(kind=wi)  :: i

    do i = 1, size( y, 3 ), 1
      call genrand_int32_2d( y(:,:,i) )
    end do
    return

  end subroutine genrand_int32_3d

  subroutine genrand_int32_4d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:,:,:), intent(out) :: y
    
    integer(kind=wi)  :: i

    do i = 1, size( y, 4 ), 1
      call genrand_int32_3d( y(:,:,:,i) )
    end do
    return

  end subroutine genrand_int32_4d

  subroutine genrand_int32_5d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:,:,:,:), intent(out) :: y
    
    integer(kind=wi)  :: i

    do i = 1, size( y, 5 ), 1
      call genrand_int32_4d( y(:,:,:,:,i) )
    end do
    return

  end subroutine genrand_int32_5d

  subroutine genrand_int32_6d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:,:,:,:,:), intent(out) :: y
    
    integer(kind=wi)  :: i

    do i = 1, size( y, 6 ), 1
      call genrand_int32_5d( y(:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_int32_6d

  subroutine genrand_int32_7d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:,:,:,:,:,:), intent(out) :: y
    
    integer(kind=wi)  :: i

    do i = 1, size( y, 7 ), 1
      call genrand_int32_6d( y(:,:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_int32_7d

  ! generates a random number on [0,0x7fffffff]-interval
  subroutine genrand_int31_0d( y )

    intrinsic :: ishft

    integer(kind=wi), intent(out) :: y

    call genrand_int32_0d( y )
    y = ishft( y, -1 )
    return

  end subroutine genrand_int31_0d

  subroutine genrand_int31_1d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:), intent(out) :: y

    integer(kind=wi)  :: i

    do i = 1, size( y, 1 ), 1
      call genrand_int31_0d( y(i) )
    end do
    return

  end subroutine genrand_int31_1d

  subroutine genrand_int31_2d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:), intent(out) :: y

    integer(kind=wi)  :: i

    do i = 1, size( y, 2 ), 1
      call genrand_int31_1d( y(:,i) )
    end do
    return

  end subroutine genrand_int31_2d

  subroutine genrand_int31_3d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:,:), intent(out) :: y

    integer(kind=wi)  :: i

    do i = 1, size( y, 3 ), 1
      call genrand_int31_2d( y(:,:,i) )
    end do
    return

  end subroutine genrand_int31_3d

  subroutine genrand_int31_4d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:,:,:), intent(out) :: y

    integer(kind=wi)  :: i

    do i = 1, size( y, 4 ), 1
      call genrand_int31_3d( y(:,:,:,i) )
    end do
    return

  end subroutine genrand_int31_4d

  subroutine genrand_int31_5d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:,:,:,:), intent(out) :: y

    integer(kind=wi)  :: i

    do i = 1, size( y, 5 ), 1
      call genrand_int31_4d( y(:,:,:,:,i) )
    end do
    return

  end subroutine genrand_int31_5d

  subroutine genrand_int31_6d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:,:,:,:,:), intent(out) :: y

    integer(kind=wi)  :: i

    do i = 1, size( y, 6 ), 1
      call genrand_int31_5d( y(:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_int31_6d

  subroutine genrand_int31_7d( y )

    intrinsic :: size

    integer(kind=wi), dimension(:,:,:,:,:,:,:), intent(out) :: y

    integer(kind=wi)  :: i

    do i = 1, size( y, 7 ), 1
      call genrand_int31_6d( y(:,:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_int31_7d

  ! generates a random number on [0,1]-real-interval
  subroutine genrand_real1_0d( r )

    intrinsic :: real

    real(kind=wr), intent(out)  :: r

    integer(kind=wi)  :: a

    call genrand_int32_0d( a )
    r = real( a, kind=wr ) * pi232_1 + p231d232_1
    ! divided by 2^32-1
    return

  end subroutine genrand_real1_0d

  subroutine genrand_real1_1d( r )

    intrinsic :: size

    real(kind=wr), dimension(:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 1 ), 1
      call genrand_real1_0d( r(i) )
    end do
    return

  end subroutine genrand_real1_1d

  subroutine genrand_real1_2d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 2 ), 1
      call genrand_real1_1d( r(:,i) )
    end do
    return

  end subroutine genrand_real1_2d

  subroutine genrand_real1_3d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 3 ), 1
      call genrand_real1_2d( r(:,:,i) )
    end do
    return

  end subroutine genrand_real1_3d

  subroutine genrand_real1_4d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 4 ), 1
      call genrand_real1_3d( r(:,:,:,i) )
    end do
    return

  end subroutine genrand_real1_4d

  subroutine genrand_real1_5d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 5 ), 1
      call genrand_real1_4d( r(:,:,:,:,i) )
    end do
    return

  end subroutine genrand_real1_5d

  subroutine genrand_real1_6d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 6 ), 1
      call genrand_real1_5d( r(:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_real1_6d

  subroutine genrand_real1_7d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 7 ), 1
      call genrand_real1_6d( r(:,:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_real1_7d

  ! generates a random number on [0,1)-real-interval
  subroutine genrand_real2_0d( r )

    intrinsic :: real

    real(kind=wr), intent(out)  :: r

    integer(kind=wi)  :: a

    call genrand_int32_0d( a )
    r = real( a, kind=wr ) * pi232 + 0.5_wr
    ! divided by 2^32
    return

  end subroutine genrand_real2_0d

  subroutine genrand_real2_1d( r )

    intrinsic :: size

    real(kind=wr), dimension(:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 1 ), 1
      call genrand_real2_0d( r(i) )
    end do
    return

  end subroutine genrand_real2_1d

  subroutine genrand_real2_2d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 2 ), 1
      call genrand_real2_1d( r(:,i) )
    end do
    return

  end subroutine genrand_real2_2d

  subroutine genrand_real2_3d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 3 ), 1
      call genrand_real2_2d( r(:,:,i) )
    end do
    return

  end subroutine genrand_real2_3d

  subroutine genrand_real2_4d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 4 ), 1
      call genrand_real2_3d( r(:,:,:,i) )
    end do
    return

  end subroutine genrand_real2_4d

  subroutine genrand_real2_5d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 5 ), 1
      call genrand_real2_4d( r(:,:,:,:,i) )
    end do
    return

  end subroutine genrand_real2_5d

  subroutine genrand_real2_6d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 6 ), 1
      call genrand_real2_5d( r(:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_real2_6d

  subroutine genrand_real2_7d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 7 ), 1
      call genrand_real2_6d( r(:,:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_real2_7d

  ! generates a random number on (0,1)-real-interval
  subroutine genrand_real3_0d( r )

    intrinsic :: real

    real(kind=wr), intent(out)  :: r

    integer(kind=wi)  :: a

    call genrand_int32_0d( a )
    r = real( a, kind=wr ) * pi232 + p231_5d232
    ! divided by 2^32 
    return

  end subroutine genrand_real3_0d

  subroutine genrand_real3_1d( r )

    intrinsic :: size

    real(kind=wr), dimension(:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 1 ), 1
      call genrand_real3_0d( r(i) )
    end do
    return

  end subroutine genrand_real3_1d

  subroutine genrand_real3_2d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 2 ), 1
      call genrand_real3_1d( r(:,i) )
    end do
    return

  end subroutine genrand_real3_2d

  subroutine genrand_real3_3d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 3 ), 1
      call genrand_real3_2d( r(:,:,i) )
    end do
    return

  end subroutine genrand_real3_3d

  subroutine genrand_real3_4d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 4 ), 1
      call genrand_real3_3d( r(:,:,:,i) )
    end do
    return

  end subroutine genrand_real3_4d

  subroutine genrand_real3_5d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 5 ), 1
      call genrand_real3_4d( r(:,:,:,:,i) )
    end do
    return

  end subroutine genrand_real3_5d

  subroutine genrand_real3_6d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 6 ), 1
      call genrand_real3_5d( r(:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_real3_6d

  subroutine genrand_real3_7d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 7 ), 1
      call genrand_real3_6d( r(:,:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_real3_7d

  ! generates a random number on [0,1) with 53-bit resolution
  subroutine genrand_res53_0d( r )

    intrinsic :: ishft, real

    real(kind=wr), intent(out)  :: r

    integer(kind=wi)  :: a, b

    call genrand_int32_0d( a )
    call genrand_int32_0d( b )
    a = ishft( a, -5 )
    b = ishft( b, -6 )
    r = real( a, kind=wr ) * pi227 + real( b, kind=wr ) * pi253
    return

  end subroutine genrand_res53_0d

  subroutine genrand_res53_1d( r )

    intrinsic :: size

    real(kind=wr), dimension(:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 1 ), 1
      call genrand_res53_0d( r(i) )
    end do
    return

  end subroutine genrand_res53_1d

  subroutine genrand_res53_2d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 2 ), 1
      call genrand_res53_1d( r(:,i) )
    end do
    return

  end subroutine genrand_res53_2d

  subroutine genrand_res53_3d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 3 ), 1
      call genrand_res53_2d( r(:,:,i) )
    end do
    return

  end subroutine genrand_res53_3d

  subroutine genrand_res53_4d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 4 ), 1
      call genrand_res53_3d( r(:,:,:,i) )
    end do
    return

  end subroutine genrand_res53_4d

  subroutine genrand_res53_5d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 5 ), 1
      call genrand_res53_4d( r(:,:,:,:,i) )
    end do
    return

  end subroutine genrand_res53_5d

  subroutine genrand_res53_6d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 6 ), 1
      call genrand_res53_5d( r(:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_res53_6d

  subroutine genrand_res53_7d( r )

    intrinsic :: size

    real(kind=wr), dimension(:,:,:,:,:,:,:), intent(out)  :: r

    integer(kind=wi)  :: i

    do i = 1, size( r, 7 ), 1
      call genrand_res53_6d( r(:,:,:,:,:,:,i) )
    end do
    return

  end subroutine genrand_res53_7d
  ! These real versions are due to Isaku Wada, 2002/01/09 added 
  ! Altered by José Sousa genrand_real[1-3] will not return exactely
  ! the same values but should have the same properties and are faster

end module mt95

