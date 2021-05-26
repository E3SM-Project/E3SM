module fftpack51D
!-----------------------------------------------------------------------
! Purpose: 
! 
! selected Fast Fourier Transform routines taken from FFTPACK5.1
! Used for Explicit Scalar Momentum Tendency (ESMT) code 
!
! Revision history: 
! Jan, 2018 - Walter Hannah
!              initial version - adapted from https://people.sc.fsu.edu/~jburkardt/f_src/fftpack5.1d/fftpack5.1d.f90
!
!---------------------------------------------------------------------------
  
  use params, only: crm_rknd

  implicit none

  public rfft1i
  public rfft1f
  public rfft1b

  contains

!========================================================================================
!========================================================================================
subroutine rfft1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! RFFT1I: initialization for RFFT1B and RFFT1F.
!
!  Discussion:
!
!    RFFT1I initializes array WSAVE for use in its companion routines
!    RFFT1B and RFFT1F.  The prime factorization of N together with a
!    tabulation of the trigonometric functions are computed and stored
!    in array WSAVE.  Separate WSAVE arrays are required for different
!    values of N.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    25 March 2005
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be
!    transformed.  The transform is most efficient when N is a product of
!    small primes.
!
!    Output, real ( kind = crm_rknd ) WSAVE(LENSAV), containing the prime factors of
!    N and also containing certain trigonometric values which will be used in
!    routines RFFT1B or RFFT1F.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough.
!
  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) n
  real ( kind = crm_rknd ) wsave(lensav)

  ier = 0

  if ( lensav < n + int ( log ( real ( n, kind = 4 ) ) ) + 4 ) then
    ier = 2
    call xerfft ( 'rfft1i ', 3 )
    return
  end if

  if ( n == 1 ) then
    return
  end if

  call rffti1 ( n, wsave(1), wsave(n+1) )

  return
end
!========================================================================================
!========================================================================================
subroutine rffti1 ( n, wa, fac )

!*****************************************************************************80
!
!! RFFTI1 is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number for which factorization
!    and other information is needed.
!
!    Output, real ( kind = crm_rknd ) WA(N), trigonometric information.
!
!    Output, real ( kind = crm_rknd ) FAC(15), factorization information.
!    FAC(1) is N, FAC(2) is NF, the number of factors, and FAC(3:NF+2) are the
!    factors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) argh
  real ( kind = 8 ) argld
  real ( kind = crm_rknd ) fac(15)
  real ( kind = crm_rknd ) fi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipm
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) ld
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nfm1
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry
  real ( kind = 8 ) tpi
  real ( kind = crm_rknd ) wa(n)

  nl = n
  nf = 0
  j = 0

  do while ( 1 < nl )

    j = j + 1

    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if

    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nf = nf + 1
      fac(nf+2) = real ( ntry, kind = 4 )
      nl = nq
!
!  If 2 is a factor, make sure it appears first in the list of factors.
!
      if ( ntry == 2 ) then
        if ( nf /= 1 ) then
          do i = 2, nf
            ib = nf - i + 2
            fac(ib+2) = fac(ib+1)
          end do
          fac(3) = 2.0E+00
        end if
      end if

    end do

  end do

  fac(1) = real ( n, kind = 4 )
  fac(2) = real ( nf, kind = 4 )
  tpi = 8.0D+00 * atan ( 1.0D+00 )
  argh = tpi / real ( n, kind = 4 )
  is = 0
  nfm1 = nf-1
  l1 = 1

  do k1 = 1, nfm1
    ip = int ( fac(k1+2) )
    ld = 0
    l2 = l1 * ip
    ido = n / l2
    ipm = ip - 1
    do j = 1, ipm
      ld = ld + l1
      i = is
      argld = real ( ld, kind = 4 ) * argh
      fi = 0.0E+00
      do ii = 3, ido, 2
        i = i + 2
        fi = fi + 1.0E+00
        arg = fi * argld
        wa(i-1) = real ( cos ( arg ), kind = 4 )
        wa(i) = real ( sin ( arg ), kind = 4 )
      end do
      is = is + ido
    end do
    l1 = l2
  end do

  return
end
!========================================================================================
!========================================================================================
subroutine rfft1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! RFFT1F: real single precision forward fast Fourier transform, 1D.
!
!  Discussion:
!
!    RFFT1F computes the one-dimensional Fourier transform of a periodic
!    sequence within a real array.  This is referred to as the forward
!    transform or Fourier analysis, transforming the sequence from physical
!    to spectral space.  This transform is normalized since a call to
!    RFFT1F followed by a call to RFFT1B (or vice-versa) reproduces the
!    original array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    25 March 2005
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be
!    transformed.  The transform is most efficient when N is a product of
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in
!    array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = crm_rknd ) R(LENR), on input, contains the sequence
!    to be transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.
!    LENR must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = crm_rknd ) WSAVE(LENSAV).  WSAVE's contents must be
!    initialized with a call to RFFT1I before the first call to routine RFFT1F
!    or RFFT1B for a given transform length N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = crm_rknd ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.
!    LENWRK must be at least N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough:
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough.
!
  implicit none

  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) n
  real ( kind = crm_rknd ) work(lenwrk)
  real ( kind = crm_rknd ) wsave(lensav)
  real ( kind = crm_rknd ) r(lenr)

  ier = 0

  if ( lenr < inc * ( n - 1 ) + 1 ) then
    ier = 1
    call xerfft ( 'rfft1f ', 6 )
    return
  end if

  if ( lensav < n + int ( log ( real ( n, kind = 4 ) ) ) + 4 ) then
    ier = 2
    call xerfft ( 'rfft1f ', 8 )
    return
  end if

  if ( lenwrk < n ) then
    ier = 3
    call xerfft ( 'rfft1f ', 10 )
    return
  end if

  if ( n == 1 ) then
    return
  end if

  call rfftf1 ( n, inc, r, work, wsave, wsave(n+1) )

  return
end
!========================================================================================
!========================================================================================
subroutine rfftf1 ( n, in, c, ch, wa, fac )

!*****************************************************************************80
!
!! RFFTF1 is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) in
  integer ( kind = 4 ) n

  real ( kind = crm_rknd ) c(in,*)
  real ( kind = crm_rknd ) ch(*)
  real ( kind = crm_rknd ) fac(15)
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) kh
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = crm_rknd ) sn
  real ( kind = crm_rknd ) tsn
  real ( kind = crm_rknd ) tsnm
  real ( kind = crm_rknd ) wa(n)

  nf = int ( fac(2) )
  na = 1
  l2 = n
  iw = n

  do k1 = 1, nf

    kh = nf - k1
    ip = int ( fac(kh+3) )
    l1 = l2 / ip
    ido = n / l2
    idl1 = ido * l1
    iw = iw - ( ip - 1 ) * ido
    na = 1 - na

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call r1f4kf ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2), wa(ix3) )
      else
        call r1f4kf ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2), wa(ix3) )
      end if

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call r1f2kf ( ido, l1, c, in, ch, 1, wa(iw) )
      else
        call r1f2kf ( ido, l1, ch, 1, c, in, wa(iw) )
      end if

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call r1f3kf ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2) )
      else
        call r1f3kf ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2) )
      end if

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call r1f5kf ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call r1f5kf ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

    else

      if ( ido == 1 ) then
        na = 1 - na
      end if

      if ( na == 0 ) then
        call r1fgkf ( ido, ip, l1, idl1, c, c, c, in, ch, ch, 1, wa(iw) )
        na = 1
      else
        call r1fgkf ( ido, ip, l1, idl1, ch, ch, ch, 1, c, c, in, wa(iw) )
        na = 0
      end if

    end if

    l2 = l1

  end do

  sn = 1.0E+00 / real ( n, kind = 4 )
  tsn = 2.0E+00 / real ( n, kind = 4 )
  tsnm = -tsn
  modn = mod ( n, 2 )
  nl = n - 2

  if ( modn /= 0 ) then
    nl = n - 1
  end if

  if ( na == 0 ) then

    c(1,1) = sn * ch(1)
    do j = 2, nl, 2
      c(1,j) = tsn * ch(j)
      c(1,j+1) = tsnm * ch(j+1)
    end do

    if ( modn == 0 ) then
      c(1,n) = sn * ch(n)
    end if

  else

    c(1,1) = sn * c(1,1)

    do j = 2, nl, 2
      c(1,j) = tsn * c(1,j)
      c(1,j+1) = tsnm * c(1,j+1)
    end do

    if ( modn == 0 ) then
      c(1,n) = sn * c(1,n)
    end if

  end if

  return
end
!========================================================================================
!========================================================================================
subroutine rfft1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! RFFT1B: real single precision backward fast Fourier transform, 1D.
!
!  Discussion:
!
!    RFFT1B computes the one-dimensional Fourier transform of a periodic
!    sequence within a real array.  This is referred to as the backward
!    transform or Fourier synthesis, transforming the sequence from
!    spectral to physical space.  This transform is normalized since a
!    call to RFFT1B followed by a call to RFFT1F (or vice-versa) reproduces
!    the original array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    25 March 2005
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be
!    transformed.  The transform is most efficient when N is a product of
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations,
!    in array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = crm_rknd ) R(LENR), on input, the data to be
!    transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.
!    LENR must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = crm_rknd ) WSAVE(LENSAV).  WSAVE's contents must be
!    initialized with a call to RFFT1I before the first call to routine
!    RFFT1F or RFFT1B for a given transform length N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = crm_rknd ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.
!    LENWRK must be at least N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough.
!
  implicit none

  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) n
  real ( kind = crm_rknd ) r(lenr)
  real ( kind = crm_rknd ) work(lenwrk)
  real ( kind = crm_rknd ) wsave(lensav)

  ier = 0

  if ( lenr < inc * ( n - 1 ) + 1 ) then
    ier = 1
    call xerfft ( 'rfft1b ', 6 )
    return
  end if

  if ( lensav < n + int ( log ( real ( n, kind = 4 ) ) ) + 4 ) then
    ier = 2
    call xerfft ( 'rfft1b ', 8 )
    return
  end if

  if ( lenwrk < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RFFT1B - Fatal error!'
    write ( *, '(a)' ) '  LENWRK < N:'
    write ( *, '(a,i6)' ) '  LENWRK = ', lenwrk
    write ( *, '(a,i6)' ) '  N = ', n
    ier = 3
    call xerfft ( 'rfft1b ', 10 )
    return
  end if

  if ( n == 1 ) then
    return
  end if

  call rfftb1 ( n, inc, r, work, wsave, wsave(n+1) )

  return
end
!========================================================================================
!========================================================================================
subroutine rfftb1 ( n, in, c, ch, wa, fac )

!*****************************************************************************80
!
!! RFFTB1 is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) in
  integer ( kind = 4 ) n

  real ( kind = crm_rknd ) c(in,*)
  real ( kind = crm_rknd ) ch(*)
  real ( kind = crm_rknd ) fac(15)
  real ( kind = crm_rknd ) half
  real ( kind = crm_rknd ) halfm
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = crm_rknd ) wa(n)

  nf = int ( fac(2) )
  na = 0

  do k1 = 1, nf

    ip = int ( fac(k1+2) )
    na = 1 - na

    if ( 5 < ip ) then
      if ( k1 /= nf ) then
        na = 1 - na
      end if
    end if

  end do

  half = 0.5E+00
  halfm = -0.5E+00
  modn = mod ( n, 2 )
  nl = n - 2
  if ( modn /= 0 ) then
    nl = n - 1
  end if

  if ( na == 0 ) then

    do j = 2, nl, 2
      c(1,j) = half * c(1,j)
      c(1,j+1) = halfm * c(1,j+1)
    end do

  else

    ch(1) = c(1,1)
    ch(n) = c(1,n)

    do j = 2, nl, 2
      ch(j) = half*c(1,j)
      ch(j+1) = halfm*c(1,j+1)
    end do

  end if

  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = int ( fac(k1+2) )
    l2 = ip * l1
    ido = n / l2
    idl1 = ido * l1

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call r1f4kb ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2), wa(ix3) )
      else
        call r1f4kb ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call r1f2kb ( ido, l1, c, in, ch, 1, wa(iw) )
      else
        call r1f2kb ( ido, l1, ch, 1, c, in, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call r1f3kb ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2) )
      else
        call r1f3kb ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call r1f5kb ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call r1f5kb ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call r1fgkb ( ido, ip, l1, idl1, c, c, c, in, ch, ch, 1, wa(iw) )
      else
        call r1fgkb ( ido, ip, l1, idl1, ch, ch, ch, 1, c, c, in, wa(iw) )
      end if

      if ( ido == 1 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * ido

  end do

  return
end
!========================================================================================
!========================================================================================
subroutine r1f2kf ( ido, l1, cc, in1, ch, in2, wa1 )

!*****************************************************************************80
!
!! R1F2KF is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) ch(in2,ido,2,l1)
  real ( kind = crm_rknd ) cc(in1,ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = crm_rknd ) wa1(ido)

  do k = 1, l1
    ch(1,1,1,k)   = cc(1,1,k,1) + cc(1,1,k,2)
    ch(1,ido,2,k) = cc(1,1,k,1) - cc(1,1,k,2)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2

    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i
        ch(1,i,1,k) = cc(1,i,k,1)   + ( wa1(i-2) * cc(1,i,k,2) &
                                      - wa1(i-1) * cc(1,i-1,k,2) )
        ch(1,ic,2,k) = -cc(1,i,k,1) + ( wa1(i-2) * cc(1,i,k,2) &
                                      - wa1(i-1) * cc(1,i-1,k,2) )
        ch(1,i-1,1,k) = cc(1,i-1,k,1)  + ( wa1(i-2) * cc(1,i-1,k,2) &
                                         + wa1(i-1) * cc(1,i,k,2))
        ch(1,ic-1,2,k) = cc(1,i-1,k,1) - ( wa1(i-2) * cc(1,i-1,k,2) &
                                         + wa1(i-1) * cc(1,i,k,2))
      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    ch(1,1,2,k) = -cc(1,ido,k,2)
    ch(1,ido,1,k) = cc(1,ido,k,1)
  end do

  return
end
!========================================================================================
!========================================================================================
subroutine r1f3kf ( ido, l1, cc, in1, ch, in2, wa1, wa2 )

!*****************************************************************************80
!
!! R1F3KF is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) arg
  real ( kind = crm_rknd ) cc(in1,ido,l1,3)
  real ( kind = crm_rknd ) ch(in2,ido,3,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = crm_rknd ) taui
  real ( kind = crm_rknd ) taur
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)

  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 3.0E+00
  taur = cos ( arg )
  taui = sin ( arg )

  do k = 1, l1
    ch(1,1,1,k) = cc(1,1,k,1)          + ( cc(1,1,k,2) + cc(1,1,k,3) )
    ch(1,1,3,k) =                 taui * ( cc(1,1,k,3) - cc(1,1,k,2) )
    ch(1,ido,2,k) = cc(1,1,k,1) + taur * ( cc(1,1,k,2) + cc(1,1,k,3) )
  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      ch(1,i-1,1,k) = cc(1,i-1,k,1)+((wa1(i-2)*cc(1,i-1,k,2)+ &
        wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3)))
      ch(1,i,1,k) = cc(1,i,k,1)+((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3)))
      ch(1,i-1,3,k) = (cc(1,i-1,k,1)+taur*((wa1(i-2)* &
        cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)* &
        cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))))+(taui*((wa1(i-2)* &
        cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa2(i-2)* &
        cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))))
      ch(1,ic-1,2,k) = (cc(1,i-1,k,1)+taur*((wa1(i-2)* &
        cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)* &
        cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))))-(taui*((wa1(i-2)* &
        cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa2(i-2)* &
        cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))))
      ch(1,i,3,k) = (cc(1,i,k,1)+taur*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))))+(taui*((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2))))
      ch(1,ic,2,k) = (taui*((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2))))-(cc(1,i,k,1)+taur*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))))
    end do
  end do

  return
end
!========================================================================================
!========================================================================================
subroutine r1f4kf ( ido, l1, cc, in1, ch, in2, wa1, wa2, wa3 )

!*****************************************************************************80
!
!! R1F4KF is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) cc(in1,ido,l1,4)
  real ( kind = crm_rknd ) ch(in2,ido,4,l1)
  real ( kind = crm_rknd ) hsqt2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)
  real ( kind = crm_rknd ) wa3(ido)

  hsqt2 = sqrt ( 2.0E+00 ) / 2.0E+00

  do k = 1, l1
    ch(1,1,1,k)   = ( cc(1,1,k,2) + cc(1,1,k,4) ) &
                  + ( cc(1,1,k,1) + cc(1,1,k,3) )
    ch(1,ido,4,k) = ( cc(1,1,k,1) + cc(1,1,k,3) ) &
                  - ( cc(1,1,k,2) + cc(1,1,k,4) )
    ch(1,ido,2,k) = cc(1,1,k,1) - cc(1,1,k,3)
    ch(1,1,3,k)   = cc(1,1,k,4) - cc(1,1,k,2)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2

    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i
        ch(1,i-1,1,k) = ((wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
          cc(1,i,k,2))+(wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
          cc(1,i,k,4)))+(cc(1,i-1,k,1)+(wa2(i-2)*cc(1,i-1,k,3)+ &
          wa2(i-1)*cc(1,i,k,3)))
        ch(1,ic-1,4,k) = (cc(1,i-1,k,1)+(wa2(i-2)*cc(1,i-1,k,3)+ &
          wa2(i-1)*cc(1,i,k,3)))-((wa1(i-2)*cc(1,i-1,k,2)+ &
          wa1(i-1)*cc(1,i,k,2))+(wa3(i-2)*cc(1,i-1,k,4)+ &
          wa3(i-1)*cc(1,i,k,4)))
        ch(1,i,1,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
          cc(1,i-1,k,2))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
          cc(1,i-1,k,4)))+(cc(1,i,k,1)+(wa2(i-2)*cc(1,i,k,3)- &
          wa2(i-1)*cc(1,i-1,k,3)))
        ch(1,ic,4,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
          cc(1,i-1,k,2))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
          cc(1,i-1,k,4)))-(cc(1,i,k,1)+(wa2(i-2)*cc(1,i,k,3)- &
          wa2(i-1)*cc(1,i-1,k,3)))
        ch(1,i-1,3,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
          cc(1,i-1,k,2))-(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
          cc(1,i-1,k,4)))+(cc(1,i-1,k,1)-(wa2(i-2)*cc(1,i-1,k,3)+ &
          wa2(i-1)*cc(1,i,k,3)))
        ch(1,ic-1,2,k) = (cc(1,i-1,k,1)-(wa2(i-2)*cc(1,i-1,k,3)+ &
          wa2(i-1)*cc(1,i,k,3)))-((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
          cc(1,i-1,k,2))-(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
          cc(1,i-1,k,4)))
        ch(1,i,3,k) = ((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
          cc(1,i,k,4))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
          cc(1,i,k,2)))+(cc(1,i,k,1)-(wa2(i-2)*cc(1,i,k,3)- &
          wa2(i-1)*cc(1,i-1,k,3)))
        ch(1,ic,2,k) = ((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
          cc(1,i,k,4))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
          cc(1,i,k,2)))-(cc(1,i,k,1)-(wa2(i-2)*cc(1,i,k,3)- &
          wa2(i-1)*cc(1,i-1,k,3)))
       end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    ch(1,ido,1,k) = (hsqt2*(cc(1,ido,k,2)-cc(1,ido,k,4)))+ cc(1,ido,k,1)
    ch(1,ido,3,k) = cc(1,ido,k,1)-(hsqt2*(cc(1,ido,k,2)- cc(1,ido,k,4)))
    ch(1,1,2,k) = (-hsqt2*(cc(1,ido,k,2)+cc(1,ido,k,4)))- cc(1,ido,k,3)
    ch(1,1,4,k) = (-hsqt2*(cc(1,ido,k,2)+cc(1,ido,k,4)))+ cc(1,ido,k,3)
  end do

  return
end
!========================================================================================
!========================================================================================
subroutine r1f5kf ( ido, l1, cc, in1, ch, in2, wa1, wa2, wa3, wa4 )

!*****************************************************************************80
!
!! R1F5KF is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) arg
  real ( kind = crm_rknd ) cc(in1,ido,l1,5)
  real ( kind = crm_rknd ) ch(in2,ido,5,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = crm_rknd ) ti11
  real ( kind = crm_rknd ) ti12
  real ( kind = crm_rknd ) tr11
  real ( kind = crm_rknd ) tr12
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)
  real ( kind = crm_rknd ) wa3(ido)
  real ( kind = crm_rknd ) wa4(ido)

  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 5.0E+00
  tr11 = cos ( arg )
  ti11 = sin ( arg )
  tr12 = cos ( 2.0E+00 * arg )
  ti12 = sin ( 2.0E+00 * arg )

  do k = 1, l1

    ch(1,1,1,k) = cc(1,1,k,1) + ( cc(1,1,k,5) + cc(1,1,k,2) ) &
                              + ( cc(1,1,k,4) + cc(1,1,k,3) )

    ch(1,ido,2,k) = cc(1,1,k,1) + tr11 * ( cc(1,1,k,5) + cc(1,1,k,2) ) &
                                + tr12 * ( cc(1,1,k,4) + cc(1,1,k,3) )

    ch(1,1,3,k) =                 ti11 * ( cc(1,1,k,5) - cc(1,1,k,2) ) &
                                + ti12 * ( cc(1,1,k,4) - cc(1,1,k,3) )

    ch(1,ido,4,k) = cc(1,1,k,1) + tr12 * ( cc(1,1,k,5) + cc(1,1,k,2) ) &
                                + tr11 * ( cc(1,1,k,4) + cc(1,1,k,3) )

    ch(1,1,5,k) =                 ti12 * ( cc(1,1,k,5) - cc(1,1,k,2) ) &
                                - ti11 * ( cc(1,1,k,4) - cc(1,1,k,3) )
  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      ch(1,i-1,1,k) = cc(1,i-1,k,1)+((wa1(i-2)*cc(1,i-1,k,2)+ &
        wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
        cc(1,i,k,5)))+((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))+(wa3(i-2)*cc(1,i-1,k,4)+ &
        wa3(i-1)*cc(1,i,k,4)))
      ch(1,i,1,k) = cc(1,i,k,1)+((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
        cc(1,i-1,k,5)))+((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
        cc(1,i-1,k,4)))
      ch(1,i-1,3,k) = cc(1,i-1,k,1)+tr11* &
        ( wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2) &
        +wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5))+tr12* &
        ( wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3) &
        +wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))+ti11* &
        ( wa1(i-2)*cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2) &
        -(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))+ti12* &
        ( wa2(i-2)*cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3) &
        -(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4)))
      ch(1,ic-1,2,k) = cc(1,i-1,k,1)+tr11* &
        ( wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2) &
        +wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5))+tr12* &
        ( wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3) &
        +wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))-(ti11* &
        ( wa1(i-2)*cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2) &
        -(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))+ti12* &
        ( wa2(i-2)*cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3) &
        -(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))
      ch(1,i,3,k) = (cc(1,i,k,1)+tr11*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
        cc(1,i-1,k,5)))+tr12*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
        cc(1,i-1,k,4))))+(ti11*((wa4(i-2)*cc(1,i-1,k,5)+ &
        wa4(i-1)*cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2)))+ti12*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
        cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))))
      ch(1,ic,2,k) = (ti11*((wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
        cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2)))+ti12*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
        cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))))-(cc(1,i,k,1)+tr11*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
        cc(1,i-1,k,5)))+tr12*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
        cc(1,i-1,k,4))))
      ch(1,i-1,5,k) = (cc(1,i-1,k,1)+tr12*((wa1(i-2)* &
        cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)* &
        cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5)))+tr11*((wa2(i-2)* &
        cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))+(wa3(i-2)* &
        cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))))+(ti12*((wa1(i-2)* &
        cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa4(i-2)* &
        cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))-ti11*((wa2(i-2)* &
        cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))-(wa3(i-2)* &
        cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))
      ch(1,ic-1,4,k) = (cc(1,i-1,k,1)+tr12*((wa1(i-2)* &
        cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)* &
        cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5)))+tr11*((wa2(i-2)* &
        cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))+(wa3(i-2)* &
        cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))))-(ti12*((wa1(i-2)* &
        cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa4(i-2)* &
        cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))-ti11*((wa2(i-2)* &
        cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))-(wa3(i-2)* &
        cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))
      ch(1,i,5,k) = (cc(1,i,k,1)+tr12*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
        cc(1,i-1,k,5)))+tr11*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
        cc(1,i-1,k,4))))+(ti12*((wa4(i-2)*cc(1,i-1,k,5)+ &
        wa4(i-1)*cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2)))-ti11*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
        cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))))
      ch(1,ic,4,k) = (ti12*((wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
        cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2)))-ti11*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
        cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))))-(cc(1,i,k,1)+tr12*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
        cc(1,i-1,k,5)))+tr11*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
        cc(1,i-1,k,4))))
     end do
  end do

  return
end
!========================================================================================
!========================================================================================
subroutine r1fgkf ( ido, ip, l1, idl1, cc, c1, c2, in1, ch, ch2, in2, wa )

!*****************************************************************************80
!
!! R1FGKF is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) ai1
  real ( kind = crm_rknd ) ai2
  real ( kind = crm_rknd ) ar1
  real ( kind = crm_rknd ) ar1h
  real ( kind = crm_rknd ) ar2
  real ( kind = crm_rknd ) ar2h
  real ( kind = crm_rknd ) arg
  real ( kind = crm_rknd ) c1(in1,ido,l1,ip)
  real ( kind = crm_rknd ) c2(in1,idl1,ip)
  real ( kind = crm_rknd ) cc(in1,ido,ip,l1)
  real ( kind = crm_rknd ) ch(in2,ido,l1,ip)
  real ( kind = crm_rknd ) ch2(in2,idl1,ip)
  real ( kind = crm_rknd ) dc2
  real ( kind = crm_rknd ) dcp
  real ( kind = crm_rknd ) ds2
  real ( kind = crm_rknd ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nbd
  real ( kind = crm_rknd ) tpi
  real ( kind = crm_rknd ) wa(ido)

  tpi = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 )
  arg = tpi / real ( ip, kind = 4 )
  dcp = cos ( arg )
  dsp = sin ( arg )
  ipph = ( ip + 1 ) / 2
  ipp2 = ip + 2
  idp2 = ido + 2
  nbd = ( ido - 1 ) / 2

  if ( ido == 1 ) then

    do ik = 1, idl1
      c2(1,ik,1) = ch2(1,ik,1)
    end do

  else

    do ik = 1, idl1
      ch2(1,ik,1) = c2(1,ik,1)
    end do

    do j = 2, ip
      do k = 1, l1
        ch(1,1,k,j) = c1(1,1,k,j)
      end do
    end do

    if ( l1 < nbd ) then

      is = -ido

      do j = 2, ip
        is = is + ido
        do k = 1, l1
          idij = is
          do i = 3, ido, 2
            idij = idij + 2
            ch(1,i-1,k,j) = wa(idij-1)*c1(1,i-1,k,j)+wa(idij) *c1(1,i,k,j)
            ch(1,i,k,j) = wa(idij-1)*c1(1,i,k,j)-wa(idij) *c1(1,i-1,k,j)
          end do
        end do
      end do

    else

      is = -ido

      do j = 2, ip
        is = is + ido
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          do k = 1, l1
            ch(1,i-1,k,j) = wa(idij-1)*c1(1,i-1,k,j)+wa(idij) *c1(1,i,k,j)
            ch(1,i,k,j) = wa(idij-1)*c1(1,i,k,j)-wa(idij) *c1(1,i-1,k,j)
          end do
        end do
      end do

    end if

    if ( nbd < l1 ) then

      do j = 2, ipph
        jc = ipp2 - j
        do i = 3, ido, 2
          do k = 1, l1
            c1(1,i-1,k,j) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
            c1(1,i-1,k,jc) = ch(1,i,k,j)-ch(1,i,k,jc)
            c1(1,i,k,j) = ch(1,i,k,j)+ch(1,i,k,jc)
            c1(1,i,k,jc) = ch(1,i-1,k,jc)-ch(1,i-1,k,j)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ipp2 - j
        do k = 1, l1
          do i = 3, ido, 2
            c1(1,i-1,k,j) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
            c1(1,i-1,k,jc) = ch(1,i,k,j)-ch(1,i,k,jc)
            c1(1,i,k,j) = ch(1,i,k,j)+ch(1,i,k,jc)
            c1(1,i,k,jc) = ch(1,i-1,k,jc)-ch(1,i-1,k,j)
          end do
        end do
      end do

    end if

  end if

  do j = 2, ipph
    jc = ipp2 - j
    do k = 1, l1
      c1(1,1,k,j) = ch(1,1,k,j)+ch(1,1,k,jc)
      c1(1,1,k,jc) = ch(1,1,k,jc)-ch(1,1,k,j)
    end do
  end do

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ipp2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 = dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      ch2(1,ik,l) = c2(1,ik,1)+ar1*c2(1,ik,2)
      ch2(1,ik,lc) = ai1*c2(1,ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph
      jc = ipp2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2 = dc2 * ai2 + ds2 * ar2
      ar2 = ar2h
      do ik = 1, idl1
        ch2(1,ik,l) = ch2(1,ik,l)+ar2*c2(1,ik,j)
        ch2(1,ik,lc) = ch2(1,ik,lc)+ai2*c2(1,ik,jc)
      end do
    end do

  end do

  do j = 2, ipph
    do ik = 1, idl1
      ch2(1,ik,1) = ch2(1,ik,1)+c2(1,ik,j)
    end do
  end do

  if ( ido < l1 ) then

    do i = 1, ido
      do k = 1, l1
        cc(1,i,1,k) = ch(1,i,k,1)
      end do
    end do

  else

    do k = 1, l1
      do i = 1, ido
        cc(1,i,1,k) = ch(1,i,k,1)
      end do
    end do

  end if

  do j = 2, ipph
    jc = ipp2 - j
    j2 = j+j
    do k = 1, l1
      cc(1,ido,j2-2,k) = ch(1,1,k,j)
      cc(1,1,j2-1,k) = ch(1,1,k,jc)
    end do
  end do

  if ( ido == 1 ) then
    return
  end if

  if ( nbd < l1 ) then

    do j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      do i = 3, ido, 2
        ic = idp2 - i
        do k = 1, l1
          cc(1,i-1,j2-1,k) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
          cc(1,ic-1,j2-2,k) = ch(1,i-1,k,j)-ch(1,i-1,k,jc)
          cc(1,i,j2-1,k) = ch(1,i,k,j)+ch(1,i,k,jc)
          cc(1,ic,j2-2,k) = ch(1,i,k,jc)-ch(1,i,k,j)
        end do
      end do
   end do

  else

    do j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      do k = 1, l1
        do i = 3, ido, 2
          ic = idp2 - i
          cc(1,i-1,j2-1,k)  = ch(1,i-1,k,j) + ch(1,i-1,k,jc)
          cc(1,ic-1,j2-2,k) = ch(1,i-1,k,j) - ch(1,i-1,k,jc)
          cc(1,i,j2-1,k)    = ch(1,i,k,j)   + ch(1,i,k,jc)
          cc(1,ic,j2-2,k)   = ch(1,i,k,jc)  - ch(1,i,k,j)
        end do
      end do
    end do

  end if

  return
end
!========================================================================================
!========================================================================================
subroutine r1f2kb ( ido, l1, cc, in1, ch, in2, wa1 )

!*****************************************************************************80
!
!! R1F2KB is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) cc(in1,ido,2,l1)
  real ( kind = crm_rknd ) ch(in2,ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = crm_rknd ) wa1(ido)

  do k = 1, l1
    ch(1,1,k,1) = cc(1,1,1,k) + cc(1,ido,2,k)
    ch(1,1,k,2) = cc(1,1,1,k) - cc(1,ido,2,k)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2
    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i

        ch(1,i-1,k,1) = cc(1,i-1,1,k) + cc(1,ic-1,2,k)
        ch(1,i,k,1)   = cc(1,i,1,k)   - cc(1,ic,2,k)

        ch(1,i-1,k,2) = wa1(i-2) * ( cc(1,i-1,1,k) - cc(1,ic-1,2,k) ) &
                      - wa1(i-1) * ( cc(1,i,1,k)   + cc(1,ic,2,k) )
        ch(1,i,k,2)   = wa1(i-2) * ( cc(1,i,1,k)   + cc(1,ic,2,k) ) &
                      + wa1(i-1) * ( cc(1,i-1,1,k) - cc(1,ic-1,2,k) )

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    ch(1,ido,k,1) =     cc(1,ido,1,k) + cc(1,ido,1,k)
    ch(1,ido,k,2) = - ( cc(1,1,2,k)   + cc(1,1,2,k) )
  end do

  return
end
!========================================================================================
!========================================================================================
subroutine r1f3kb ( ido, l1, cc, in1, ch, in2, wa1, wa2 )

!*****************************************************************************80
!
!! R1F3KB is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) arg
  real ( kind = crm_rknd ) cc(in1,ido,3,l1)
  real ( kind = crm_rknd ) ch(in2,ido,l1,3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = crm_rknd ) taui
  real ( kind = crm_rknd ) taur
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)

  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 3.0E+00
  taur = cos ( arg )
  taui = sin ( arg )

  do k = 1, l1
    ch(1,1,k,1) = cc(1,1,1,k) + 2.0E+00 * cc(1,ido,2,k)
    ch(1,1,k,2) = cc(1,1,1,k) + 2.0E+00 * taur * cc(1,ido,2,k) &
                              - 2.0E+00 * taui * cc(1,1,3,k)
    ch(1,1,k,3) = cc(1,1,1,k) + 2.0E+00 * taur * cc(1,ido,2,k) &
                              + 2.0E+00 * taui * cc(1,1,3,k)
  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      ch(1,i-1,k,1) = cc(1,i-1,1,k)+(cc(1,i-1,3,k)+cc(1,ic-1,2,k))
      ch(1,i,k,1) = cc(1,i,1,k)+(cc(1,i,3,k)-cc(1,ic,2,k))
      ch(1,i-1,k,2) = wa1(i-2)* &
        ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))- &
        (taui*(cc(1,i,3,k)+cc(1,ic,2,k)))) -wa1(i-1)* &
        ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))+ &
        (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))))
      ch(1,i,k,2) = wa1(i-2)* &
        ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))+ &
        (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))) +wa1(i-1)* &
        ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))- &
        (taui*(cc(1,i,3,k)+cc(1,ic,2,k))))
      ch(1,i-1,k,3) = wa2(i-2)* &
        ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))+ &
        (taui*(cc(1,i,3,k)+cc(1,ic,2,k)))) -wa2(i-1)* &
        ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))- &
        (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))))
      ch(1,i,k,3) = wa2(i-2)* &
        ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))- &
        (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))) +wa2(i-1)* &
        ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))+ &
        (taui*(cc(1,i,3,k)+cc(1,ic,2,k))))
    end do
  end do

  return
end
!========================================================================================
!========================================================================================
subroutine r1f4kb ( ido, l1, cc, in1, ch, in2, wa1, wa2, wa3 )

!*****************************************************************************80
!
!! R1F4KB is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) cc(in1,ido,4,l1)
  real ( kind = crm_rknd ) ch(in2,ido,l1,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = crm_rknd ) sqrt2
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)
  real ( kind = crm_rknd ) wa3(ido)

  sqrt2 = sqrt ( 2.0E+00 )

  do k = 1, l1
    ch(1,1,k,3) = ( cc(1,1,1,k)   + cc(1,ido,4,k) ) &
                - ( cc(1,ido,2,k) + cc(1,ido,2,k) )
    ch(1,1,k,1) = ( cc(1,1,1,k)   + cc(1,ido,4,k) ) &
                + ( cc(1,ido,2,k) + cc(1,ido,2,k) )
    ch(1,1,k,4) = ( cc(1,1,1,k)   - cc(1,ido,4,k) ) &
                + ( cc(1,1,3,k)   + cc(1,1,3,k) )
    ch(1,1,k,2) = ( cc(1,1,1,k)   - cc(1,ido,4,k) ) &
                - ( cc(1,1,3,k)   + cc(1,1,3,k) )
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2

    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i
        ch(1,i-1,k,1) = (cc(1,i-1,1,k)+cc(1,ic-1,4,k)) &
          +(cc(1,i-1,3,k)+cc(1,ic-1,2,k))
        ch(1,i,k,1) = (cc(1,i,1,k)-cc(1,ic,4,k)) &
          +(cc(1,i,3,k)-cc(1,ic,2,k))
        ch(1,i-1,k,2) = wa1(i-2)*((cc(1,i-1,1,k)-cc(1,ic-1,4,k)) &
          -(cc(1,i,3,k)+cc(1,ic,2,k)))-wa1(i-1) &
          *((cc(1,i,1,k)+cc(1,ic,4,k))+(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))
        ch(1,i,k,2) = wa1(i-2)*((cc(1,i,1,k)+cc(1,ic,4,k)) &
          +(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))+wa1(i-1) &
          *((cc(1,i-1,1,k)-cc(1,ic-1,4,k))-(cc(1,i,3,k)+cc(1,ic,2,k)))
        ch(1,i-1,k,3) = wa2(i-2)*((cc(1,i-1,1,k)+cc(1,ic-1,4,k)) &
          -(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))-wa2(i-1) &
          *((cc(1,i,1,k)-cc(1,ic,4,k))-(cc(1,i,3,k)-cc(1,ic,2,k)))
        ch(1,i,k,3) = wa2(i-2)*((cc(1,i,1,k)-cc(1,ic,4,k)) &
          -(cc(1,i,3,k)-cc(1,ic,2,k)))+wa2(i-1) &
          *((cc(1,i-1,1,k)+cc(1,ic-1,4,k))-(cc(1,i-1,3,k) &
          +cc(1,ic-1,2,k)))
        ch(1,i-1,k,4) = wa3(i-2)*((cc(1,i-1,1,k)-cc(1,ic-1,4,k)) &
          +(cc(1,i,3,k)+cc(1,ic,2,k)))-wa3(i-1) &
          *((cc(1,i,1,k)+cc(1,ic,4,k))-(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))
        ch(1,i,k,4) = wa3(i-2)*((cc(1,i,1,k)+cc(1,ic,4,k)) &
          -(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))+wa3(i-1) &
          *((cc(1,i-1,1,k)-cc(1,ic-1,4,k))+(cc(1,i,3,k)+cc(1,ic,2,k)))
      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    ch(1,ido,k,1) = ( cc(1,ido,1,k) + cc(1,ido,3,k) ) &
                  + ( cc(1,ido,1,k) + cc(1,ido,3,k))
    ch(1,ido,k,2) = sqrt2 * ( ( cc(1,ido,1,k) - cc(1,ido,3,k) ) &
                            - ( cc(1,1,2,k)   + cc(1,1,4,k) ) )
    ch(1,ido,k,3) = ( cc(1,1,4,k) - cc(1,1,2,k) ) &
                  + ( cc(1,1,4,k) - cc(1,1,2,k) )
    ch(1,ido,k,4) = -sqrt2 * ( ( cc(1,ido,1,k) - cc(1,ido,3,k) ) &
                             + ( cc(1,1,2,k) + cc(1,1,4,k) ) )
  end do

  return
end
!========================================================================================
!========================================================================================
subroutine r1f5kb ( ido, l1, cc, in1, ch, in2, wa1, wa2, wa3, wa4 )

!*****************************************************************************80
!
!! R1F5KB is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) arg
  real ( kind = crm_rknd ) cc(in1,ido,5,l1)
  real ( kind = crm_rknd ) ch(in2,ido,l1,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = crm_rknd ) ti11
  real ( kind = crm_rknd ) ti12
  real ( kind = crm_rknd ) tr11
  real ( kind = crm_rknd ) tr12
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)
  real ( kind = crm_rknd ) wa3(ido)
  real ( kind = crm_rknd ) wa4(ido)

  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 5.0E+00
  tr11 = cos ( arg )
  ti11 = sin ( arg )
  tr12 = cos ( 2.0E+00 * arg )
  ti12 = sin ( 2.0E+00 * arg )

  do k = 1, l1

    ch(1,1,k,1) = cc(1,1,1,k) + 2.0E+00 * cc(1,ido,2,k) &
                              + 2.0E+00 * cc(1,ido,4,k)

    ch(1,1,k,2) = ( cc(1,1,1,k) &
      +   tr11 * 2.0E+00 * cc(1,ido,2,k) + tr12 * 2.0E+00 * cc(1,ido,4,k) ) &
      - ( ti11 * 2.0E+00 * cc(1,1,3,k)   + ti12 * 2.0E+00 * cc(1,1,5,k))

    ch(1,1,k,3) = ( cc(1,1,1,k) &
      +   tr12 * 2.0E+00 * cc(1,ido,2,k) + tr11 * 2.0E+00 * cc(1,ido,4,k) ) &
      - ( ti12 * 2.0E+00 * cc(1,1,3,k)   - ti11 * 2.0E+00 * cc(1,1,5,k))

    ch(1,1,k,4) = ( cc(1,1,1,k) &
      +   tr12 * 2.0E+00 * cc(1,ido,2,k) + tr11 * 2.0E+00 * cc(1,ido,4,k) ) &
      + ( ti12 * 2.0E+00 * cc(1,1,3,k)   - ti11 * 2.0E+00 * cc(1,1,5,k))

    ch(1,1,k,5) = ( cc(1,1,1,k) &
      +   tr11 * 2.0E+00 * cc(1,ido,2,k) + tr12 * 2.0E+00 * cc(1,ido,4,k) ) &
      + ( ti11 * 2.0E+00 * cc(1,1,3,k)   + ti12 * 2.0E+00 * cc(1,1,5,k) )

  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      ch(1,i-1,k,1) = cc(1,i-1,1,k)+(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +(cc(1,i-1,5,k)+cc(1,ic-1,4,k))
      ch(1,i,k,1) = cc(1,i,1,k)+(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +(cc(1,i,5,k)-cc(1,ic,4,k))
      ch(1,i-1,k,2) = wa1(i-2)*((cc(1,i-1,1,k)+tr11* &
        (cc(1,i-1,3,k)+cc(1,ic-1,2,k))+tr12 &
        *(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa1(i-1)*((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))+(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))
      ch(1,i,k,2) = wa1(i-2)*((cc(1,i,1,k)+tr11*(cc(1,i,3,k) &
        -cc(1,ic,2,k))+tr12*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti11*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))+ti12 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))+wa1(i-1) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k) &
        +cc(1,ic-1,2,k))+tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k))) &
        -(ti11*(cc(1,i,3,k)+cc(1,ic,2,k))+ti12 &
        *(cc(1,i,5,k)+cc(1,ic,4,k))))
      ch(1,i-1,k,3) = wa2(i-2) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa2(i-1) &
        *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
      cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))
      ch(1,i,k,3) = wa2(i-2) &
        *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa2(i-1) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k))))
      ch(1,i-1,k,4) = wa3(i-2) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa3(i-1) &
        *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
      cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        -(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))
      ch(1,i,k,4) = wa3(i-2) &
        *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        -(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa3(i-1) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k))))
      ch(1,i-1,k,5) = wa4(i-2) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa4(i-1) &
        *((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))-(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))
      ch(1,i,k,5) = wa4(i-2) &
        *((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))-(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa4(i-1) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k))))
    end do
  end do

  return
end
!========================================================================================
!========================================================================================
subroutine r1fgkb ( ido, ip, l1, idl1, cc, c1, c2, in1, ch, ch2, in2, wa )

!*****************************************************************************80
!
!! R1FGKB is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    27 March 2009
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) ai1
  real ( kind = crm_rknd ) ai2
  real ( kind = crm_rknd ) ar1
  real ( kind = crm_rknd ) ar1h
  real ( kind = crm_rknd ) ar2
  real ( kind = crm_rknd ) ar2h
  real ( kind = crm_rknd ) arg
  real ( kind = crm_rknd ) c1(in1,ido,l1,ip)
  real ( kind = crm_rknd ) c2(in1,idl1,ip)
  real ( kind = crm_rknd ) cc(in1,ido,ip,l1)
  real ( kind = crm_rknd ) ch(in2,ido,l1,ip)
  real ( kind = crm_rknd ) ch2(in2,idl1,ip)
  real ( kind = crm_rknd ) dc2
  real ( kind = crm_rknd ) dcp
  real ( kind = crm_rknd ) ds2
  real ( kind = crm_rknd ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nbd
  real ( kind = crm_rknd ) tpi
  real ( kind = crm_rknd ) wa(ido)

  tpi = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 )
  arg = tpi / real ( ip, kind = 4 )
  dcp = cos ( arg )
  dsp = sin ( arg )
  idp2 = ido + 2
  nbd = ( ido - 1 ) / 2
  ipp2 = ip + 2
  ipph = ( ip + 1 ) / 2

  if ( ido < l1 ) then
    do i = 1, ido
      do k = 1, l1
        ch(1,i,k,1) = cc(1,i,1,k)
      end do
    end do
  else
    do k = 1, l1
      do i = 1, ido
        ch(1,i,k,1) = cc(1,i,1,k)
      end do
    end do
  end if

  do j = 2, ipph
    jc = ipp2 - j
    j2 = j + j
    do k = 1, l1
      ch(1,1,k,j) = cc(1,ido,j2-2,k)+cc(1,ido,j2-2,k)
      ch(1,1,k,jc) = cc(1,1,j2-1,k)+cc(1,1,j2-1,k)
    end do
  end do

  if ( ido == 1 ) then

  else if ( nbd < l1 ) then

    do j = 2, ipph
      jc = ipp2 - j
      do i = 3, ido, 2
        ic = idp2 - i
        do k = 1, l1
          ch(1,i-1,k,j) = cc(1,i-1,2*j-1,k)+cc(1,ic-1,2*j-2,k)
          ch(1,i-1,k,jc) = cc(1,i-1,2*j-1,k)-cc(1,ic-1,2*j-2,k)
          ch(1,i,k,j) = cc(1,i,2*j-1,k)-cc(1,ic,2*j-2,k)
          ch(1,i,k,jc) = cc(1,i,2*j-1,k)+cc(1,ic,2*j-2,k)
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ipp2 - j
      do k = 1, l1
        do i = 3, ido, 2
          ic = idp2 - i
          ch(1,i-1,k,j) = cc(1,i-1,2*j-1,k)+cc(1,ic-1,2*j-2,k)
          ch(1,i-1,k,jc) = cc(1,i-1,2*j-1,k)-cc(1,ic-1,2*j-2,k)
          ch(1,i,k,j) = cc(1,i,2*j-1,k)-cc(1,ic,2*j-2,k)
          ch(1,i,k,jc) = cc(1,i,2*j-1,k)+cc(1,ic,2*j-2,k)
        end do
      end do
    end do

  end if

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ipp2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 = dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      c2(1,ik,l) = ch2(1,ik,1)+ar1*ch2(1,ik,2)
      c2(1,ik,lc) = ai1*ch2(1,ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ipp2 - j
      ar2h = dc2*ar2-ds2*ai2
      ai2 = dc2*ai2+ds2*ar2
      ar2 = ar2h

      do ik = 1, idl1
        c2(1,ik,l) = c2(1,ik,l)+ar2*ch2(1,ik,j)
        c2(1,ik,lc) = c2(1,ik,lc)+ai2*ch2(1,ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    do ik = 1, idl1
      ch2(1,ik,1) = ch2(1,ik,1)+ch2(1,ik,j)
    end do
  end do

  do j = 2, ipph
    jc = ipp2 - j
    do k = 1, l1
      ch(1,1,k,j) = c1(1,1,k,j)-c1(1,1,k,jc)
      ch(1,1,k,jc) = c1(1,1,k,j)+c1(1,1,k,jc)
    end do
  end do

  if ( ido == 1 ) then

  else if ( nbd < l1 ) then

    do j = 2, ipph
      jc = ipp2 - j
      do i = 3, ido, 2
        do k = 1, l1
          ch(1,i-1,k,j)  = c1(1,i-1,k,j) - c1(1,i,k,jc)
          ch(1,i-1,k,jc) = c1(1,i-1,k,j) + c1(1,i,k,jc)
          ch(1,i,k,j)    = c1(1,i,k,j)   + c1(1,i-1,k,jc)
          ch(1,i,k,jc)   = c1(1,i,k,j)   - c1(1,i-1,k,jc)
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ipp2 - j
      do k = 1, l1
        do i = 3, ido, 2
          ch(1,i-1,k,j) = c1(1,i-1,k,j)-c1(1,i,k,jc)
          ch(1,i-1,k,jc) = c1(1,i-1,k,j)+c1(1,i,k,jc)
          ch(1,i,k,j) = c1(1,i,k,j)+c1(1,i-1,k,jc)
          ch(1,i,k,jc) = c1(1,i,k,j)-c1(1,i-1,k,jc)
        end do
      end do
    end do

  end if

  if ( ido == 1 ) then
    return
  end if

  do ik = 1, idl1
    c2(1,ik,1) = ch2(1,ik,1)
  end do

  do j = 2, ip
    do k = 1, l1
      c1(1,1,k,j) = ch(1,1,k,j)
    end do
  end do

  if ( l1 < nbd ) then

    is = -ido
    do j = 2, ip
       is = is + ido
       do k = 1, l1
         idij = is
         do i = 3, ido, 2
           idij = idij + 2
           c1(1,i-1,k,j) = wa(idij-1)*ch(1,i-1,k,j)-wa(idij)* ch(1,i,k,j)
           c1(1,i,k,j) = wa(idij-1)*ch(1,i,k,j)+wa(idij)* ch(1,i-1,k,j)
         end do
       end do
    end do

  else

    is = -ido

    do j = 2, ip
      is = is + ido
      idij = is
      do i = 3, ido, 2
        idij = idij + 2
        do k = 1, l1
           c1(1,i-1,k,j) = wa(idij-1) * ch(1,i-1,k,j) - wa(idij) * ch(1,i,k,j)
           c1(1,i,k,j)   = wa(idij-1) * ch(1,i,k,j)   + wa(idij) * ch(1,i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
!========================================================================================
!========================================================================================

!========================================================================================
!========================================================================================

!========================================================================================
!========================================================================================

!========================================================================================
! XERFFT is an error handler for the FFTPACK routines.
!========================================================================================
subroutine xerfft ( srname, info )

  !*****************************************************************************80
  !
  !  Discussion:
  !
  !    XERFFT is an error handler for FFTPACK version 5.1 routines.
  !    It is called by an FFTPACK 5.1 routine if an input parameter has an
  !    invalid value.  A message is printed and execution stops.
  !
  !    Installers may consider modifying the stop statement in order to
  !    call system-specific exception-handling facilities.
  !
  !  License:
  !
  !    Licensed under the GNU General Public License (GPL).
  !    Copyright (C) 1995-2004, Scientific Computing Division,
  !    University Corporation for Atmospheric Research
  !
  !  Modified:
  !
  !    31 July 2011
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Paul Swarztrauber,
  !    Vectorizing the Fast Fourier Transforms,
  !    in Parallel Computations,
  !    edited by G. Rodrigue,
  !    Academic Press, 1982.
  !
  !    Paul Swarztrauber,
  !    Fast Fourier Transform Algorithms for Vector Computers,
  !    Parallel Computing, pages 45-63, 1984.
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) SRNAME, the name of the calling routine.
  !
  !    Input, integer ( kind = 4 ) INFO, an error code.  When a single invalid 
  !    parameter in the parameter list of the calling routine has been detected, 
  !    INFO is the position of that parameter.  In the case when an illegal 
  !    combination of LOT, JUMP, N, and INC has been detected, the calling 
  !    subprogram calls XERFFT with INFO = -1.
  !
  implicit none

  integer ( kind = 4 ) info
  character ( len = * ) srname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XERFFT - Fatal error!'

  if ( 1 <= info ) then
    write ( *, '(a,a,a,i3,a)') '  On entry to ', trim ( srname ), &
      ' parameter number ', info, ' had an illegal value.'
  else if ( info == -1 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameters LOT, JUMP, N and INC are inconsistent.'
  else if ( info == -2 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter L is greater than LDIM.'
  else if ( info == -3 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter M is greater than MDIM.'
  else if ( info == -5 ) then
    write( *, '(a,a,a,a)') '  Within ', trim ( srname ), &
      ' input error returned by lower level routine.'
  else if ( info == -6 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter LDIM is less than 2*(L/2+1).'
  end if

  stop
end
!========================================================================================
!========================================================================================
function xercon ( inc, jump, n, lot )

!*****************************************************************************80
!
!! XERCON checks INC, JUMP, N and LOT for consistency.
!
!  Discussion:
!
!    Positive integers INC, JUMP, N and LOT are "consistent" if,
!    for any values I1 and I2 < N, and J1 and J2 < LOT,
!
!      I1 * INC + J1 * JUMP = I2 * INC + J2 * JUMP
!
!    can only occur if I1 = I2 and J1 = J2.
!
!    For multiple FFT's to execute correctly, INC, JUMP, N and LOT must
!    be consistent, or else at least one array element will be
!    transformed more than once.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    25 March 2005
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INC, JUMP, N, LOT, the parameters to check.
!
!    Output, logical XERCON, is TRUE if the parameters are consistent.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnew
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lcm
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) n
  logical              xercon

  i = inc
  j = jump

  do while ( j /= 0 )
    jnew = mod ( i, j )
    i = j
    j = jnew
  end do
!
!  LCM = least common multiple of INC and JUMP.
!
  lcm = ( inc * jump ) / i

  if ( lcm <= ( n - 1 ) * inc .and. lcm <= ( lot - 1 ) * jump ) then
    xercon = .false.
  else
    xercon = .true.
  end if

  return
end
!========================================================================================
!========================================================================================
end module fftpack51D