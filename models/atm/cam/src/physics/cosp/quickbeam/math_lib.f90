! MATH_LIB: Mathematics procedures for F90
! Compiled/Modified:
!   07/01/06  John Haynes (haynes@atmos.colostate.edu)
! 
! gamma (function)
! path_integral (function)
! avint (subroutine)
  
  module math_lib
  implicit none

  contains

! ----------------------------------------------------------------------------
! function GAMMA
! ----------------------------------------------------------------------------
  function gamma(x)
  implicit none
!
! Purpose:
!   Returns the gamma function
!
! Input:
!   [x]   value to compute gamma function of
!
! Returns:
!   gamma(x)
!
! Coded:
!   02/02/06  John Haynes (haynes@atmos.colostate.edu)
!   (original code of unknown origin)

! ----- INPUTS -----
  real*8, intent(in) :: x
  
! ----- OUTPUTS -----
  real*8 :: gamma

! ----- INTERNAL -----  
  real*8 :: pi,ga,z,r,gr
  real*8 :: g(26)
  integer :: k,m1,m
       
  pi = acos(-1.)    
  if (x ==int(x)) then
    if (x > 0.0) then
      ga=1.0
      m1=x-1
      do k=2,m1
        ga=ga*k
      enddo
    else
      ga=1.0+300
    endif
  else
    if (abs(x) > 1.0) then
      z=abs(x)
      m=int(z)
      r=1.0
      do k=1,m
        r=r*(z-k)
      enddo
      z=z-m
    else
      z=x
    endif
    data g/1.0,0.5772156649015329, &
           -0.6558780715202538, -0.420026350340952d-1, &
           0.1665386113822915,-.421977345555443d-1, &
           -.96219715278770d-2, .72189432466630d-2, &
           -.11651675918591d-2, -.2152416741149d-3, &
           .1280502823882d-3, -.201348547807d-4, &
           -.12504934821d-5, .11330272320d-5, &
           -.2056338417d-6, .61160950d-8, &
           .50020075d-8, -.11812746d-8, &
           .1043427d-9, .77823d-11, &
          -.36968d-11, .51d-12, &
          -.206d-13, -.54d-14, .14d-14, .1d-15/
    gr=g(26)
    do k=25,1,-1
      gr=gr*z+g(k)
    enddo 
    ga=1.0/(gr*z)
    if (abs(x) > 1.0) then
      ga=ga*r
      if (x < 0.0) ga=-pi/(x*ga*sin(pi*x))
    endif
  endif
  gamma = ga
  return
  end function gamma
  
! ----------------------------------------------------------------------------
! function PATH_INTEGRAL 
! ----------------------------------------------------------------------------
  function path_integral(f,s,i1,i2)
  use m_mrgrnk
  use array_lib
  implicit none
!
! Purpose:
!   evalues the integral (f ds) between f(index=i1) and f(index=i2)
!   using the AVINT procedure
!
! Inputs:
!   [f]    functional values
!   [s]    abscissa values
!   [i1]   index of lower limit
!   [i2]   index of upper limit
!
! Returns:
!   result of path integral
!
! Notes:
!   [s] may be in forward or reverse numerical order
!
! Requires:
!   mrgrnk package
!
! Created:
!   02/02/06  John Haynes (haynes@atmos.colostate.edu)

! ----- INPUTS -----  
  real*8, intent(in), dimension(:) :: f,s  
  integer, intent(in) :: i1, i2

! ---- OUTPUTS -----
  real*8 :: path_integral  
  
! ----- INTERNAL -----    
  real*8 :: sumo, deltah, val
  integer*4 :: nelm, j
  integer*4, dimension(i2-i1+1) :: idx
  real*8, dimension(i2-i1+1) :: f_rev, s_rev

  nelm = i2-i1+1
  if (nelm > 3) then
    call mrgrnk(s(i1:i2),idx)
    s_rev = s(idx)
    f_rev = f(idx)
    call avint(f_rev(i1:i2),s_rev(i1:i2),(i2-i1+1), &
      s_rev(i1),s_rev(i2), val)
    path_integral = val
  else
     sumo = 0.
     do j=i1,i2
       deltah = abs(s(i1+1)-s(i1))
       sumo = sumo + f(j)*deltah
    enddo
    path_integral = sumo
  endif 
  ! print *, sumo

  return
  end function path_integral
  
! ----------------------------------------------------------------------------
! subroutine AVINT
! ----------------------------------------------------------------------------
  subroutine avint ( ftab, xtab, ntab, a_in, b_in, result )
  implicit none
!
! Purpose:
!   estimate the integral of unevenly spaced data
!
! Inputs:
!   [ftab]     functional values
!   [xtab]     abscissa values
!   [ntab]     number of elements of [ftab] and [xtab]
!   [a]        lower limit of integration
!   [b]        upper limit of integration
!
! Outputs:
!   [result]   approximate value of integral
!
! Reference:
!   From SLATEC libraries, in public domain
!
!***********************************************************************
!
!  AVINT estimates the integral of unevenly spaced data.
!
!  Discussion:
!
!    The method uses overlapping parabolas and smoothing.
!
!  Modified:
!
!    30 October 2000
!    4 January 2008, A. Bodas-Salcedo. Error control for XTAB taken out of
!                    loop to allow vectorization.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    P E Hennion,
!    Algorithm 77,
!    Interpolation, Differentiation and Integration,
!    Communications of the Association for Computing Machinery,
!    Volume 5, page 96, 1962.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FTAB(NTAB), the function values,
!    FTAB(I) = F(XTAB(I)).
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the abscissas at which the
!    function values are given.  The XTAB's must be distinct
!    and in ascending order.
!
!    Input, integer NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.  A should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.  B should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.

  integer, intent(in) :: ntab

  integer,parameter :: KR8 = selected_real_kind(15,300)
  real ( kind = KR8 ), intent(in) :: a_in
  real ( kind = KR8 ) a
  real ( kind = KR8 ) atemp
  real ( kind = KR8 ), intent(in) :: b_in
  real ( kind = KR8 ) b
  real ( kind = KR8 ) btemp
  real ( kind = KR8 ) ca
  real ( kind = KR8 ) cb
  real ( kind = KR8 ) cc
  real ( kind = KR8 ) ctemp
  real ( kind = KR8 ), intent(in) :: ftab(ntab)
  integer i
  integer ihi
  integer ilo
  integer ind
  real ( kind = KR8 ), intent(out) :: result
  real ( kind = KR8 ) sum1
  real ( kind = KR8 ) syl
  real ( kind = KR8 ) term1
  real ( kind = KR8 ) term2
  real ( kind = KR8 ) term3
  real ( kind = KR8 ) x1
  real ( kind = KR8 ) x2
  real ( kind = KR8 ) x3
  real ( kind = KR8 ), intent(in) :: xtab(ntab)
  logical lerror
  
  lerror = .false.
  a = a_in
  b = b_in  
  
  if ( ntab < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Fatal error!'
    write ( *, '(a,i6)' ) '  NTAB is less than 3.  NTAB = ', ntab
    stop
  end if
 
  do i = 2, ntab
    if ( xtab(i) <= xtab(i-1) ) then
       lerror = .true.
       exit
    end if
  end do
  
  if (lerror) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AVINT - Fatal error!'
      write ( *, '(a)' ) '  XTAB(I) is not greater than XTAB(I-1).'
      write ( *, '(a,i6)' ) '  Here, I = ', i
      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ', xtab(i-1)
      write ( *, '(a,g14.6)' ) '  XTAB(I) =   ', xtab(i)
      stop  
  end if
 
  result = 0.0D+00
 
  if ( a == b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Warning!'
    write ( *, '(a)' ) '  A = B, integral=0.'
    return
  end if
!
!  If B < A, temporarily switch A and B, and store sign.
!
  if ( b < a ) then
    syl = b
    b = a
    a = syl
    ind = -1
  else
    syl = a
    ind = 1
  end if
!
!  Bracket A and B between XTAB(ILO) and XTAB(IHI).
!
  ilo = 1
  ihi = ntab

  do i = 1, ntab
    if ( a <= xtab(i) ) then
      exit
    end if
    ilo = ilo + 1
  end do

  ilo = max ( 2, ilo )
  ilo = min ( ilo, ntab - 1 )

  do i = 1, ntab
    if ( xtab(i) <= b ) then
      exit
    end if
    ihi = ihi - 1
  end do
  
  ihi = min ( ihi, ntab - 1 )
  ihi = max ( ilo, ihi - 1 )
!
!  Carry out approximate integration from XTAB(ILO) to XTAB(IHI).
!
  sum1 = 0.0D+00
 
  do i = ilo, ihi
 
    x1 = xtab(i-1)
    x2 = xtab(i)
    x3 = xtab(i+1)
    
    term1 = ftab(i-1) / ( ( x1 - x2 ) * ( x1 - x3 ) )
    term2 = ftab(i)   / ( ( x2 - x1 ) * ( x2 - x3 ) )
    term3 = ftab(i+1) / ( ( x3 - x1 ) * ( x3 - x2 ) )
 
    atemp = term1 + term2 + term3

    btemp = - ( x2 + x3 ) * term1 &
            - ( x1 + x3 ) * term2 &
            - ( x1 + x2 ) * term3

    ctemp = x2 * x3 * term1 + x1 * x3 * term2 + x1 * x2 * term3
 
    if ( i <= ilo ) then
      ca = atemp
      cb = btemp
      cc = ctemp
    else
      ca = 0.5D+00 * ( atemp + ca )
      cb = 0.5D+00 * ( btemp + cb )
      cc = 0.5D+00 * ( ctemp + cc )
    end if
 
    sum1 = sum1 &
          + ca * ( x2**3 - syl**3 ) / 3.0D+00 &
          + cb * 0.5D+00 * ( x2**2 - syl**2 ) &
          + cc * ( x2 - syl )
 
    ca = atemp
    cb = btemp
    cc = ctemp
 
    syl = x2
 
  end do
 
  result = sum1 &
        + ca * ( b**3 - syl**3 ) / 3.0D+00 &
        + cb * 0.5D+00 * ( b**2 - syl**2 ) &
        + cc * ( b - syl )
!
!  Restore original values of A and B, reverse sign of integral
!  because of earlier switch.
!
  if ( ind /= 1 ) then
    ind = 1
    syl = b
    b = a
    a = syl
    result = -result
  end if
 
  return
  end subroutine avint
  
  end module math_lib