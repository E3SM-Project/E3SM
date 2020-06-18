module fftpack5
!-----------------------------------------------------------------------
! Purpose: 
! 
! selected Fast Fourier Transform routines taken from FFTPACK5.1
! Used for Explicit Scalar Momentum Tendency (ESMT) code 
!
! Revision history: 
! Jan, 2018 - Walter Hannah
!              initial version - adapted from https://michalakes.svn.cloudforge.com/rrtmmic/HWRF/external/fftpack/fftpack5
!
!---------------------------------------------------------------------------
  
  use params, only: crm_rknd

  implicit none

  public rfftmi
  public rfftmf
  public rfftmb

  contains
  
!========================================================================================
! RFFTMI: initialization for RFFTMB and RFFTMF.
!========================================================================================
subroutine rfftmi ( n, wsave, lensav, ier )

  !*****************************************************************************80
  !
  !  Discussion:
  !
  !    RFFTMI initializes array WSAVE for use in its companion routines 
  !    RFFTMB and RFFTMF.  The prime factorization of N together with a 
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
  !    Input, integer ( kind = 4 ) N, the length of each sequence to be 
  !    transformed.  The transform is most efficient when N is a product of 
  !    small primes.
  !
  !    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
  !    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
  !
  !    Output, real ( kind = crm_rknd ) WSAVE(LENSAV), work array containing the prime 
  !    factors of N and also containing certain trigonometric 
  !    values which will be used in routines RFFTMB or RFFTMF.
  !
  !    Output, integer ( kind = 4 ) IER, error flag.
  !    0, successful exit;
  !    2, input parameter LENSAV not big enough.
  !
  implicit none

  integer ( kind = 4 )     ,intent(in   ) :: lensav
  integer ( kind = 4 )     ,intent(in   ) :: n
  integer ( kind = 4 )     ,intent(  out) :: ier
  real ( kind = crm_rknd ) ,intent(inout) :: wsave(lensav)

  ier = 0

  if ( lensav < n + int ( log ( real ( n, kind = crm_rknd ) ) ) + 4 ) then
    ier = 2
    call xerfft ( 'rfftmi ', 3 )
    return
  end if

  if ( n == 1 ) then
    return
  end if

  call mrfti1 ( n, wsave(1), wsave(n+1) )

  return
end
!========================================================================================
! RFFTMF: real single precision forward FFT, 1D, multiple vectors.
!========================================================================================
subroutine rfftmf ( lot, jump, n, inc, r, lenr, wsave, lensav, &
              work, lenwrk, ier )

  !*****************************************************************************80
  !
  !  Discussion:
  !
  !    RFFTMF computes the one-dimensional Fourier transform of multiple 
  !    periodic sequences within a real array.  This transform is referred 
  !    to as the forward transform or Fourier analysis, transforming the 
  !    sequences from physical to spectral space.
  !
  !    This transform is normalized since a call to RFFTMF followed
  !    by a call to RFFTMB (or vice-versa) reproduces the original array
  !    within roundoff error.
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
  !    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed 
  !    within array R.
  !
  !    Input, integer ( kind = 4 ) JUMP, the increment between the locations, in 
  !    array R, of the first elements of two consecutive sequences to be 
  !    transformed.
  !
  !    Input, integer ( kind = 4 ) N, the length of each sequence to be 
  !    transformed.  The transform is most efficient when N is a product of 
  !    small primes.
  !
  !    Input, integer ( kind = 4 ) INC, the increment between the locations, 
  !    in array R, of two consecutive elements within the same sequence.
  !
  !    Input/output, real ( kind = crm_rknd ) R(LENR), real array containing LOT 
  !    sequences, each having length N.  R can have any number of dimensions, but 
  !    the total number of locations must be at least LENR.  On input, the
  !    physical data to be transformed, on output the spectral data.
  !
  !    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
  !    LENR must be at least (LOT-1)*JUMP + INC*(N-1) + 1.
  !
  !    Input, real ( kind = crm_rknd ) WSAVE(LENSAV).  WSAVE's contents must be 
  !    initialized with a call to RFFTMI before the first call to routine RFFTMF 
  !    or RFFTMB for a given transform length N.  
  !
  !    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
  !    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
  !
  !    Workspace, real ( kind = crm_rknd ) WORK(LENWRK).
  !
  !    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
  !    LENWRK must be at least LOT*N.
  !
  !    Output, integer ( kind = 4 ) IER, error flag.
  !    0, successful exit;
  !    1, input parameter LENR not big enough;
  !    2, input parameter LENSAV not big enough;
  !    3, input parameter LENWRK not big enough;
  !    4, input parameters INC, JUMP, N, LOT are not consistent.
  !
  implicit none
  
  integer ( kind = 4 )        ,intent(in   ) :: lot
  integer ( kind = 4 )        ,intent(in   ) :: jump
  integer ( kind = 4 )        ,intent(in   ) :: n
  integer ( kind = 4 )        ,intent(in   ) :: inc
  integer ( kind = 4 )        ,intent(in   ) :: lenr
  integer ( kind = 4 )        ,intent(in   ) :: lensav
  integer ( kind = 4 )        ,intent(in   ) :: lenwrk
  integer ( kind = 4 )        ,intent(  out) :: ier
  real    ( kind = crm_rknd ) ,intent(inout) :: r(lenr)
  real    ( kind = crm_rknd ) ,intent(inout) :: work(lenwrk)
  real    ( kind = crm_rknd ) ,intent(inout) :: wsave(lensav)
  logical xercon

  ier = 0

  if ( lenr < ( lot - 1 ) * jump + inc * ( n - 1 ) + 1 ) then
    ier = 1
    call xerfft ( 'rfftmf ', 6 )
    return
  end if

  if ( lensav < n + int ( log ( real ( n, kind = crm_rknd ) ) ) + 4 ) then
    ier = 2
    call xerfft ( 'rfftmf ', 8 )
    return
  end if

  if ( lenwrk < lot * n ) then
    ier = 3
    call xerfft ( 'rfftmf ', 10 )
    return
  end if

  if ( .not. xercon ( inc, jump, n, lot ) ) then
    ier = 4
    call xerfft ( 'rfftmf ', -1 )
    return
  end if

  if ( n == 1 ) then
    return
  end if

  call mrftf1 ( lot, jump, n, inc, r, work, wsave, wsave(n+1) )

  return
end
!========================================================================================
! RFFTMB: real single precision backward FFT, 1D, multiple vectors.
!========================================================================================
subroutine rfftmb ( lot, jump, n, inc, r, lenr, wsave, lensav, &
              work, lenwrk, ier )

  !*****************************************************************************80
  !
  !  Discussion:
  !
  !    RFFTMB computes the one-dimensional Fourier transform of multiple 
  !    periodic sequences within a real array.  This transform is referred 
  !    to as the backward transform or Fourier synthesis, transforming the
  !    sequences from spectral to physical space.
  !
  !    This transform is normalized since a call to RFFTMB followed
  !    by a call to RFFTMF (or vice-versa) reproduces the original
  !    array  within roundoff error.
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
  !    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed 
  !    within array R.
  !
  !    Input, integer ( kind = 4 ) JUMP, the increment between the locations, in
  !    array R, of the first elements of two consecutive sequences to be 
  !    transformed.
  !
  !    Input, integer ( kind = 4 ) N, the length of each sequence to be 
  !    transformed.  The transform is most efficient when N is a product of 
  !    small primes.
  !
  !    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
  !    array R, of two consecutive elements within the same sequence.
  !
  !    Input/output, real ( kind = crm_rknd ) R(LENR), real array containing LOT 
  !    sequences, each having length N.  R can have any number of dimensions, 
  !    but the total number of locations must be at least LENR.  On input, the
  !    spectral data to be transformed, on output the physical data.
  !
  !    Input, integer ( kind = 4 ) LENR, the dimension of the R array. 
  !    LENR must be at least (LOT-1)*JUMP + INC*(N-1) + 1.
  !
  !    Input, real ( kind = crm_rknd ) WSAVE(LENSAV).  WSAVE's contents must be 
  !    initialized with a call to RFFTMI before the first call to routine RFFTMF 
  !    or RFFTMB for a given transform length N.  
  !
  !    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array. 
  !    LENSAV must  be at least N + INT(LOG(REAL(N))) + 4.
  !
  !    Workspace, real ( kind = crm_rknd ) WORK(LENWRK).
  !
  !    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
  !    LENWRK must be at least LOT*N.
  !
  !    Output, integer ( kind = 4 ) IER, error flag.
  !    0, successful exit;
  !    1, input parameter LENR not big enough;
  !    2, input parameter LENSAV not big enough;
  !    3, input parameter LENWRK not big enough;
  !    4, input parameters INC, JUMP, N, LOT are not consistent.
  !
  implicit none

  integer ( kind = 4 )        ,intent(in   ) :: lenr
  integer ( kind = 4 )        ,intent(in   ) :: lensav
  integer ( kind = 4 )        ,intent(in   ) :: lenwrk
  integer ( kind = 4 )        ,intent(in   ) :: inc
  integer ( kind = 4 )        ,intent(in   ) :: jump
  integer ( kind = 4 )        ,intent(in   ) :: lot
  integer ( kind = 4 )        ,intent(in   ) :: n
  integer ( kind = 4 )        ,intent(  out) :: ier
  real    ( kind = crm_rknd ) ,intent(inout) :: r(lenr)
  real    ( kind = crm_rknd ) ,intent(inout) :: work(lenwrk)
  real    ( kind = crm_rknd ) ,intent(inout) :: wsave(lensav)
  logical xercon

  ier = 0

  if ( lenr < ( lot - 1 ) * jump + inc * ( n - 1 ) + 1 ) then
    ier = 1
    call xerfft ( 'rfftmb ', 6 )
    return
  end if

  if ( lensav < n + int ( log ( real ( n, kind = crm_rknd ) ) ) + 4 ) then
    ier = 2
    call xerfft ( 'rfftmb ', 8 )
    return
  end if

  if ( lenwrk < lot * n ) then
    ier = 3
    call xerfft ( 'rfftmb ', 10 )
    return
  end if

  if ( .not. xercon ( inc, jump, n, lot ) ) then
    ier = 4
    call xerfft ( 'rfftmb ', -1 )
    return
  end if

  if ( n == 1 ) then
    return
  end if

  call mrftb1 ( lot, jump, n, inc, r, work, wsave, wsave(n+1) )

  return
end

!========================================================================================
! MRADB2 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mradb2 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1)

  !*****************************************************************************80
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
  implicit none

  integer ( kind = 4 )     ,intent(in   ) :: ido
  integer ( kind = 4 )     ,intent(in   ) :: in1
  integer ( kind = 4 )     ,intent(in   ) :: in2
  integer ( kind = 4 )     ,intent(in   ) :: l1
  integer ( kind = 4 )     ,intent(in   ) :: im1
  integer ( kind = 4 )     ,intent(in   ) :: im2
  integer ( kind = 4 )     ,intent(in   ) :: m
  real ( kind = crm_rknd ) ,intent(in   ) :: cc(in1,ido,2,l1)
  real ( kind = crm_rknd ) ,intent(in   ) :: wa1(ido)
  real ( kind = crm_rknd ) ,intent(inout) :: ch(in2,ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,1,k,1) = cc(m1,1,1,k) + cc(m1,ido,2,k)
      ch(m2,1,k,2) = cc(m1,1,1,k) - cc(m1,ido,2,k)
    end do
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2

    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch(m2,i-1,k,1) = cc(m1,i-1,1,k) + cc(m1,ic-1,2,k)
          ch(m2,i,k,1)   = cc(m1,i,1,k) - cc(m1,ic,2,k)
          ch(m2,i-1,k,2) = wa1(i-2) * ( cc(m1,i-1,1,k) - cc(m1,ic-1,2,k) ) &
                         - wa1(i-1) * ( cc(m1,i,1,k)   + cc(m1,ic,2,k) )
          ch(m2,i,k,2)   = wa1(i-2) * ( cc(m1,i,1,k)   + cc(m1,ic,2,k) ) &
                         + wa1(i-1) * ( cc(m1,i-1,1,k) - cc(m1,ic-1,2,k) )
        end do
      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,ido,k,1) = cc(m1,ido,1,k) + cc(m1,ido,1,k)
      ch(m2,ido,k,2) = -( cc(m1,1,2,k) + cc(m1,1,2,k) )
    end do
  end do

  return
end
!========================================================================================
! MRADB3 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mradb3 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2)

  !*****************************************************************************80
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
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) arg
  real ( kind = crm_rknd ) ,intent(in   ) :: cc(in1,ido,3,l1)
  real ( kind = crm_rknd ) ,intent(inout) :: ch(in2,ido,l1,3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = crm_rknd ) taui
  real ( kind = crm_rknd ) taur
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2
  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 3.0E+00
  taur = cos ( arg )
  taui = sin ( arg )

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,1,k,1) = cc(m1,1,1,k) + 2.0E+00 * cc(m1,ido,2,k)
      ch(m2,1,k,2) = cc(m1,1,1,k) + ( 2.0E+00 * taur ) * cc(m1,ido,2,k) &
        - ( 2.0E+00 * taui ) * cc(m1,1,3,k)
      ch(m2,1,k,3) = cc(m1,1,1,k) + ( 2.0E+00 * taur ) * cc(m1,ido,2,k) &
        + 2.0E+00 * taui * cc(m1,1,3,k)
    end do
  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      m2 = m2s
      do m1 = 1, m1d, im1

        m2 = m2 + im2

        ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))

        ch(m2,i,k,1) = cc(m1,i,1,k)+(cc(m1,i,3,k)-cc(m1,ic,2,k))

        ch(m2,i-1,k,2) = wa1(i-2)* &
          ((cc(m1,i-1,1,k)+taur*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))- &
          (taui*(cc(m1,i,3,k)+cc(m1,ic,2,k)))) - wa1(i-1)* &
          ((cc(m1,i,1,k)+taur*(cc(m1,i,3,k)-cc(m1,ic,2,k)))+ &
          (taui*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))))

        ch(m2,i,k,2) = wa1(i-2)* &
          ((cc(m1,i,1,k)+taur*(cc(m1,i,3,k)-cc(m1,ic,2,k)))+ &
          (taui*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))) + wa1(i-1)* &
          ((cc(m1,i-1,1,k)+taur*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))- &
          (taui*(cc(m1,i,3,k)+cc(m1,ic,2,k))))

        ch(m2,i-1,k,3) = wa2(i-2)* &
          ((cc(m1,i-1,1,k)+taur*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))+ &
          (taui*(cc(m1,i,3,k)+cc(m1,ic,2,k)))) - wa2(i-1)* &
          ((cc(m1,i,1,k)+taur*(cc(m1,i,3,k)-cc(m1,ic,2,k)))- &
          (taui*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))))

        ch(m2,i,k,3) = wa2(i-2)* &
          ((cc(m1,i,1,k)+taur*(cc(m1,i,3,k)-cc(m1,ic,2,k)))- &
          (taui*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))) + wa2(i-1)* &
          ((cc(m1,i-1,1,k)+taur*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))+ &
          (taui*(cc(m1,i,3,k)+cc(m1,ic,2,k))))

      end do
    end do
  end do

  return
end
!========================================================================================
! MRADB4 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mradb4 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3)

  !*****************************************************************************80
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
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) ,intent(in   ) :: cc(in1,ido,4,l1)
  real ( kind = crm_rknd ) ,intent(inout) :: ch(in2,ido,l1,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = crm_rknd ) sqrt2
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)
  real ( kind = crm_rknd ) wa3(ido)

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2
  sqrt2 = sqrt ( 2.0E+00 )

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,1,k,3) = (cc(m1,1,1,k)+cc(m1,ido,4,k)) &
        -(cc(m1,ido,2,k)+cc(m1,ido,2,k))
      ch(m2,1,k,1) = (cc(m1,1,1,k)+cc(m1,ido,4,k)) &
        +(cc(m1,ido,2,k)+cc(m1,ido,2,k))
      ch(m2,1,k,4) = (cc(m1,1,1,k)-cc(m1,ido,4,k)) &
        +(cc(m1,1,3,k)+cc(m1,1,3,k))
      ch(m2,1,k,2) = (cc(m1,1,1,k)-cc(m1,ido,4,k)) &
        -(cc(m1,1,3,k)+cc(m1,1,3,k))
    end do
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2

    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch(m2,i-1,k,1) = (cc(m1,i-1,1,k)+cc(m1,ic-1,4,k)) &
            +(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))
          ch(m2,i,k,1) = (cc(m1,i,1,k)-cc(m1,ic,4,k)) &
            +(cc(m1,i,3,k)-cc(m1,ic,2,k))
          ch(m2,i-1,k,2) = wa1(i-2)*((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k)) &
            -(cc(m1,i,3,k)+cc(m1,ic,2,k)))-wa1(i-1) &
            *((cc(m1,i,1,k)+cc(m1,ic,4,k))+(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))
          ch(m2,i,k,2) = wa1(i-2)*((cc(m1,i,1,k)+cc(m1,ic,4,k)) &
            +(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))) + wa1(i-1) &
            *((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k))-(cc(m1,i,3,k)+cc(m1,ic,2,k)))
          ch(m2,i-1,k,3) = wa2(i-2)*((cc(m1,i-1,1,k)+cc(m1,ic-1,4,k)) &
            -(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))) - wa2(i-1) &
            *((cc(m1,i,1,k)-cc(m1,ic,4,k))-(cc(m1,i,3,k)-cc(m1,ic,2,k)))
          ch(m2,i,k,3) = wa2(i-2)*((cc(m1,i,1,k)-cc(m1,ic,4,k)) &
            -(cc(m1,i,3,k)-cc(m1,ic,2,k))) + wa2(i-1) &
            *((cc(m1,i-1,1,k)+cc(m1,ic-1,4,k))-(cc(m1,i-1,3,k) &
            +cc(m1,ic-1,2,k)))
          ch(m2,i-1,k,4) = wa3(i-2)*((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k)) &
            +(cc(m1,i,3,k)+cc(m1,ic,2,k))) - wa3(i-1) &
            *((cc(m1,i,1,k)+cc(m1,ic,4,k))-(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))
          ch(m2,i,k,4) = wa3(i-2)*((cc(m1,i,1,k)+cc(m1,ic,4,k)) &
            -(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))) + wa3(i-1) &
            *((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k))+(cc(m1,i,3,k)+cc(m1,ic,2,k)))
        end do
      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,ido,k,1) = (cc(m1,ido,1,k)+cc(m1,ido,3,k)) &
        +(cc(m1,ido,1,k)+cc(m1,ido,3,k))
      ch(m2,ido,k,2) = sqrt2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
        -(cc(m1,1,2,k)+cc(m1,1,4,k)))
      ch(m2,ido,k,3) = (cc(m1,1,4,k)-cc(m1,1,2,k)) &
        +(cc(m1,1,4,k)-cc(m1,1,2,k))
      ch(m2,ido,k,4) = -sqrt2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
        +(cc(m1,1,2,k)+cc(m1,1,4,k)))
    end do
  end do

  return
end
!========================================================================================
! MRADB5 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mradb5 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3,wa4)

  !*****************************************************************************80
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
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) arg
  real ( kind = crm_rknd ) ,intent(in   ) :: cc(in1,ido,5,l1)
  real ( kind = crm_rknd ) ,intent(inout) :: ch(in2,ido,l1,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = crm_rknd ) ti11
  real ( kind = crm_rknd ) ti12
  real ( kind = crm_rknd ) tr11
  real ( kind = crm_rknd ) tr12
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)
  real ( kind = crm_rknd ) wa3(ido)
  real ( kind = crm_rknd ) wa4(ido)

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2
  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 5.E+00
  tr11 = cos ( arg )
  ti11 = sin ( arg )
  tr12 = cos ( 2.0E+00 * arg )
  ti12 = sin ( 2.0E+00 * arg )

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,1,k,1) = cc(m1,1,1,k) + 2.0E+00 * cc(m1,ido,2,k) &
        + 2.0E+00 * cc(m1,ido,4,k)
      ch(m2,1,k,2) = ( cc(m1,1,1,k) + tr11 * 2.0E+00 * cc(m1,ido,2,k) &
        + tr12 * 2.0E+00 * cc(m1,ido,4,k) ) - ( ti11 * 2.0E+00 * cc(m1,1,3,k) &
        + ti12 * 2.0E+00 * cc(m1,1,5,k) )
      ch(m2,1,k,3) = ( cc(m1,1,1,k) + tr12 * 2.0E+00 * cc(m1,ido,2,k) &
        + tr11 * 2.0E+00 * cc(m1,ido,4,k) ) - ( ti12 * 2.0E+00 * cc(m1,1,3,k) &
        - ti11 * 2.0E+00 * cc(m1,1,5,k) )
      ch(m2,1,k,4) = ( cc(m1,1,1,k) + tr12 * 2.0E+00 * cc(m1,ido,2,k) &
        + tr11 * 2.0E+00 * cc(m1,ido,4,k) ) + ( ti12 * 2.0E+00 * cc(m1,1,3,k) &
        - ti11 * 2.0E+00 * cc(m1,1,5,k) )
      ch(m2,1,k,5) = ( cc(m1,1,1,k) + tr11 * 2.0E+00 * cc(m1,ido,2,k) &
        + tr12 * 2.0E+00 * cc(m1,ido,4,k) ) + ( ti11 * 2.0E+00 * cc(m1,1,3,k) &
        + ti12 * 2.0E+00 * cc(m1,1,5,k) )
    end do
  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2
  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
          +(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k))
        ch(m2,i,k,1) = cc(m1,i,1,k)+(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
          +(cc(m1,i,5,k)-cc(m1,ic,4,k))
        ch(m2,i-1,k,2) = wa1(i-2)*((cc(m1,i-1,1,k)+tr11* &
          (cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))+tr12 &
          *(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(ti11*(cc(m1,i,3,k) &
          +cc(m1,ic,2,k))+ti12*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
          -wa1(i-1)*((cc(m1,i,1,k)+tr11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
          +tr12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))+(ti11*(cc(m1,i-1,3,k) &
          -cc(m1,ic-1,2,k))+ti12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))
        ch(m2,i,k,2) = wa1(i-2)*((cc(m1,i,1,k)+tr11*(cc(m1,i,3,k) &
          -cc(m1,ic,2,k))+tr12*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
          +(ti11*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))+ti12 &
          *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) + wa1(i-1) &
          *((cc(m1,i-1,1,k)+tr11*(cc(m1,i-1,3,k) &
          +cc(m1,ic-1,2,k))+tr12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k))) &
          -(ti11*(cc(m1,i,3,k)+cc(m1,ic,2,k))+ti12 &
          *(cc(m1,i,5,k)+cc(m1,ic,4,k))))
        ch(m2,i-1,k,3) = wa2(i-2) &
          *((cc(m1,i-1,1,k)+tr12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
          +tr11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(ti12*(cc(m1,i,3,k) &
          +cc(m1,ic,2,k))-ti11*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
          -wa2(i-1) &
          *((cc(m1,i,1,k)+tr12*(cc(m1,i,3,k)- &
        cc(m1,ic,2,k))+tr11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
          +(ti12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-ti11 &
          *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))
        ch(m2,i,k,3) = wa2(i-2) &
          *((cc(m1,i,1,k)+tr12*(cc(m1,i,3,k)- &
          cc(m1,ic,2,k))+tr11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
          +(ti12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-ti11 &
          *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
          + wa2(i-1) &
          *((cc(m1,i-1,1,k)+tr12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
          +tr11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(ti12*(cc(m1,i,3,k) &
          +cc(m1,ic,2,k))-ti11*(cc(m1,i,5,k)+cc(m1,ic,4,k))))
        ch(m2,i-1,k,4) = wa3(i-2) &
          *((cc(m1,i-1,1,k)+tr12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
          +tr11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(ti12*(cc(m1,i,3,k) &
          +cc(m1,ic,2,k))-ti11*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
          -wa3(i-1) &
          *((cc(m1,i,1,k)+tr12*(cc(m1,i,3,k)- &
          cc(m1,ic,2,k))+tr11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
          -(ti12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-ti11 &
          *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))
        ch(m2,i,k,4) = wa3(i-2) &
          *((cc(m1,i,1,k)+tr12*(cc(m1,i,3,k)- &
          cc(m1,ic,2,k))+tr11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
          -(ti12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-ti11 &
          *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
          + wa3(i-1) &
          *((cc(m1,i-1,1,k)+tr12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
          +tr11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(ti12*(cc(m1,i,3,k) &
          +cc(m1,ic,2,k))-ti11*(cc(m1,i,5,k)+cc(m1,ic,4,k))))
        ch(m2,i-1,k,5) = wa4(i-2) &
          *((cc(m1,i-1,1,k)+tr11*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
          +tr12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(ti11*(cc(m1,i,3,k) &
          +cc(m1,ic,2,k))+ti12*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
          -wa4(i-1) &
          *((cc(m1,i,1,k)+tr11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
          +tr12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))-(ti11*(cc(m1,i-1,3,k) &
          -cc(m1,ic-1,2,k))+ti12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))
        ch(m2,i,k,5) = wa4(i-2) &
          *((cc(m1,i,1,k)+tr11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
          +tr12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))-(ti11*(cc(m1,i-1,3,k) &
          -cc(m1,ic-1,2,k))+ti12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
          + wa4(i-1) &
          *((cc(m1,i-1,1,k)+tr11*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
          +tr12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(ti11*(cc(m1,i,3,k) &
          +cc(m1,ic,2,k))+ti12*(cc(m1,i,5,k)+cc(m1,ic,4,k))))
      end do
    end do
  end do

  return
end
!========================================================================================
! MRADBG is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mradbg (m,ido,ip,l1,idl1,cc,c1,c2,im1,in1,ch,ch2,im2,in2,wa)

  !*****************************************************************************80
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
  real ( kind = crm_rknd ) ,intent(inout) :: ch(in2,ido,l1,ip)
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
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) nbd
  real ( kind = crm_rknd ) tpi
  real ( kind = crm_rknd ) wa(ido)

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2
  tpi = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 )
  arg = tpi / real ( ip, kind = crm_rknd )
  dcp = cos ( arg )
  dsp = sin ( arg )
  idp2 = ido + 2
  nbd = ( ido - 1 ) / 2
  ipp2 = ip + 2
  ipph = ( ip + 1 ) / 2

  if ( ido < l1 ) then

    do i = 1, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch(m2,i,k,1) = cc(m1,i,1,k)
        end do
      end do
    end do

  else

    do k = 1, l1
      do i = 1, ido
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch(m2,i,k,1) = cc(m1,i,1,k)
        end do
      end do
    end do

  end if

  do j = 2, ipph
    jc = ipp2 - j
    j2 = j + j
    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch(m2,1,k,j) = cc(m1,ido,j2-2,k) + cc(m1,ido,j2-2,k)
        ch(m2,1,k,jc) = cc(m1,1,j2-1,k) + cc(m1,1,j2-1,k)
      end do
    end do
  end do

  if ( ido == 1 ) then

  else if ( nbd < l1 ) then

    do j = 2, ipph
      jc = ipp2 - j
      do i = 3, ido, 2
        ic = idp2 - i
        do k = 1, l1
          m2 = m2s
          do m1 = 1, m1d, im1
            m2 = m2 + im2
            ch(m2,i-1,k,j)  = cc(m1,i-1,2*j-1,k) + cc(m1,ic-1,2*j-2,k)
            ch(m2,i-1,k,jc) = cc(m1,i-1,2*j-1,k) - cc(m1,ic-1,2*j-2,k)
            ch(m2,i,k,j)    = cc(m1,i,2*j-1,k)   - cc(m1,ic,2*j-2,k)
            ch(m2,i,k,jc)   = cc(m1,i,2*j-1,k)   + cc(m1,ic,2*j-2,k)
          end do
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ipp2 - j
      do k = 1, l1
        do i = 3, ido, 2
          ic = idp2 - i
          m2 = m2s
          do m1 = 1, m1d, im1
            m2 = m2 + im2
            ch(m2,i-1,k,j)  = cc(m1,i-1,2*j-1,k) + cc(m1,ic-1,2*j-2,k)
            ch(m2,i-1,k,jc) = cc(m1,i-1,2*j-1,k) - cc(m1,ic-1,2*j-2,k)
            ch(m2,i,k,j)    = cc(m1,i,2*j-1,k)   - cc(m1,ic,2*j-2,k)
            ch(m2,i,k,jc)   = cc(m1,i,2*j-1,k)   + cc(m1,ic,2*j-2,k)
          end do
        end do
      end do
    end do

  end if

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ipp2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        c2(m1,ik,l)  = ch2(m2,ik,1) + ar1 * ch2(m2,ik,2)
        c2(m1,ik,lc) =                ai1 * ch2(m2,ik,ip)
      end do
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph
      jc = ipp2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2  = dc2 * ai2 + ds2 * ar2
      ar2 = ar2h
      do ik = 1, idl1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          c2(m1,ik,l)  = c2(m1,ik,l)  + ar2 * ch2(m2,ik,j)
          c2(m1,ik,lc) = c2(m1,ik,lc) + ai2 * ch2(m2,ik,jc)
        end do
      end do
    end do

  end do

  do j = 2, ipph
    do ik = 1, idl1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch2(m2,ik,1) = ch2(m2,ik,1) + ch2(m2,ik,j)
      end do
    end do
  end do

  do j = 2, ipph
    jc = ipp2 - j
    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch(m2,1,k,j)  = c1(m1,1,k,j) - c1(m1,1,k,jc)
        ch(m2,1,k,jc) = c1(m1,1,k,j) + c1(m1,1,k,jc)
      end do
    end do
  end do

  if ( ido == 1 ) then

  else if ( nbd < l1 ) then

    do j = 2, ipph
      jc = ipp2 - j
      do i = 3, ido, 2
        do k = 1, l1
          m2 = m2s
          do m1 = 1, m1d, im1
            m2 = m2 + im2
            ch(m2,i-1,k,j)  = c1(m1,i-1,k,j) - c1(m1,i,k,jc)
            ch(m2,i-1,k,jc) = c1(m1,i-1,k,j) + c1(m1,i,k,jc)
            ch(m2,i,k,j)    = c1(m1,i,k,j)   + c1(m1,i-1,k,jc)
            ch(m2,i,k,jc)   = c1(m1,i,k,j)   - c1(m1,i-1,k,jc)
          end do
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ipp2 - j
      do k = 1, l1
        do i = 3, ido, 2
          m2 = m2s
          do m1 = 1, m1d, im1
            m2 = m2 + im2
            ch(m2,i-1,k,j)  = c1(m1,i-1,k,j) - c1(m1,i,k,jc)
            ch(m2,i-1,k,jc) = c1(m1,i-1,k,j) + c1(m1,i,k,jc)
            ch(m2,i,k,j)    = c1(m1,i,k,j)   + c1(m1,i-1,k,jc)
            ch(m2,i,k,jc)   = c1(m1,i,k,j)   - c1(m1,i-1,k,jc)
          end do
        end do
      end do
    end do

  end if

  if ( ido == 1 ) then
    return
  end if

  do ik = 1, idl1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      c2(m1,ik,1) = ch2(m2,ik,1)
    end do
  end do

  do j = 2, ip
    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        c1(m1,1,k,j) = ch(m2,1,k,j)
      end do
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
          m2 = m2s
          do m1 = 1, m1d, im1
            m2 = m2 + im2
            c1(m1,i-1,k,j) = wa(idij-1) * ch(m2,i-1,k,j) &
                           - wa(idij)   * ch(m2,i,k,j)
            c1(m1,i,k,j) =   wa(idij-1) * ch(m2,i,k,j) &
                           + wa(idij)   * ch(m2,i-1,k,j)
          end do
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
          m2 = m2s
          do m1 = 1, m1d, im1
            m2 = m2 + im2
            c1(m1,i-1,k,j) = wa(idij-1) * ch(m2,i-1,k,j) &
                           - wa(idij)   * ch(m2,i,k,j)
            c1(m1,i,k,j) =   wa(idij-1) * ch(m2,i,k,j) &
                           + wa(idij)   * ch(m2,i-1,k,j)
          end do
        end do
      end do
    end do

  end if

  return
end
!========================================================================================
! MRADF2 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mradf2 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1)

  !*****************************************************************************80
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
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) cc(in1,ido,l1,2)
  real ( kind = crm_rknd ) ,intent(inout) :: ch(in2,ido,2,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = crm_rknd ) wa1(ido)

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,1,1,k)   = cc(m1,1,k,1) + cc(m1,1,k,2)
      ch(m2,ido,2,k) = cc(m1,1,k,1) - cc(m1,1,k,2)
    end do
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2

    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch(m2,i,1,k) =    cc(m1,i,k,1)   + ( wa1(i-2) * cc(m1,i,k,2) &
                                             - wa1(i-1) * cc(m1,i-1,k,2) )
          ch(m2,ic,2,k)  = -cc(m1,i,k,1)   + ( wa1(i-2) * cc(m1,i,k,2) &
                                             - wa1(i-1) * cc(m1,i-1,k,2) )
          ch(m2,i-1,1,k)  = cc(m1,i-1,k,1) + ( wa1(i-2) * cc(m1,i-1,k,2) &
                                             + wa1(i-1) * cc(m1,i,k,2))
          ch(m2,ic-1,2,k) = cc(m1,i-1,k,1) - ( wa1(i-2) * cc(m1,i-1,k,2) &
                                             + wa1(i-1) * cc(m1,i,k,2) )
        end do
      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,1,2,k) = -cc(m1,ido,k,2)
      ch(m2,ido,1,k) = cc(m1,ido,k,1)
    end do
  end do

  return
end
!========================================================================================
! MRADF3 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mradf3 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2)

  !*****************************************************************************80
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
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) arg
  real ( kind = crm_rknd ) cc(in1,ido,l1,3)
  real ( kind = crm_rknd ) ,intent(inout) :: ch(in2,ido,3,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = crm_rknd ) taui
  real ( kind = crm_rknd ) taur
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2
  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 3.0E+00
  taur = cos ( arg )
  taui = sin ( arg )

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,1,1,k)   = cc(m1,1,k,1)        + ( cc(m1,1,k,2) + cc(m1,1,k,3) )
      ch(m2,1,3,k)   =                taui * ( cc(m1,1,k,3) - cc(m1,1,k,2) )
      ch(m2,ido,2,k) = cc(m1,1,k,1) + taur * ( cc(m1,1,k,2) + cc(m1,1,k,3) )
     end do
  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+((wa1(i-2)*cc(m1,i-1,k,2)+ &
          wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)*cc(m1,i-1,k,3) + wa2(i-1)* &
          cc(m1,i,k,3)))
        ch(m2,i,1,k) = cc(m1,i,k,1)+((wa1(i-2)*cc(m1,i,k,2)- &
          wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
          cc(m1,i-1,k,3)))
        ch(m2,i-1,3,k) = (cc(m1,i-1,k,1)+taur*((wa1(i-2)* &
          cc(m1,i-1,k,2) + wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)* &
          cc(m1,i-1,k,3) + wa2(i-1)*cc(m1,i,k,3))))+(taui*((wa1(i-2)* &
          cc(m1,i,k,2) - wa1(i-1)*cc(m1,i-1,k,2))-(wa2(i-2)* &
          cc(m1,i,k,3) - wa2(i-1)*cc(m1,i-1,k,3))))
        ch(m2,ic-1,2,k) = (cc(m1,i-1,k,1)+taur*((wa1(i-2)* &
          cc(m1,i-1,k,2) + wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)* &
          cc(m1,i-1,k,3) + wa2(i-1)*cc(m1,i,k,3))))-(taui*((wa1(i-2)* &
          cc(m1,i,k,2) - wa1(i-1)*cc(m1,i-1,k,2)) - ( wa2(i-2)* &
          cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))))
        ch(m2,i,3,k) = (cc(m1,i,k,1)+taur*((wa1(i-2)*cc(m1,i,k,2)- &
          wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3) - wa2(i-1)* &
          cc(m1,i-1,k,3))))+(taui*((wa2(i-2)*cc(m1,i-1,k,3) + wa2(i-1)* &
          cc(m1,i,k,3))-(wa1(i-2)*cc(m1,i-1,k,2) + wa1(i-1)* &
          cc(m1,i,k,2))))
        ch(m2,ic,2,k) = (taui*((wa2(i-2)*cc(m1,i-1,k,3) + wa2(i-1)* &
          cc(m1,i,k,3))-(wa1(i-2)*cc(m1,i-1,k,2) + wa1(i-1)* &
          cc(m1,i,k,2))))-(cc(m1,i,k,1)+taur*((wa1(i-2)*cc(m1,i,k,2)- &
          wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
          cc(m1,i-1,k,3))))
      end do
    end do
  end do

  return
end
!========================================================================================
! MRADF4 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mradf4 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3)

  !*****************************************************************************80
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
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real    ( kind = crm_rknd ) ,intent(in   ) :: cc(in1,ido,l1,4)
  real    ( kind = crm_rknd ) ,intent(inout) :: ch(in2,ido,4,l1)
  real    ( kind = crm_rknd ) ,intent(inout) :: wa1(ido)
  real    ( kind = crm_rknd ) ,intent(inout) :: wa2(ido)
  real    ( kind = crm_rknd ) ,intent(inout) :: wa3(ido)
  real    ( kind = crm_rknd ) hsqt2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  

! m   = lot  = nz
! im1 = (jump = nx) or (1)
! im2 = (1)         or (jump = nx)
! in1 = (inc = 1)   or (lot = nz)
! in2 = (lot = nz)  or (inc = 1)

! m1d = { (nz-1)*nx+1 } or { (nz-1)*1+1 = nz }

  hsqt2 = sqrt ( 2.0E+00 ) / 2.0E+00
  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,1,1,k) =   ( cc(m1,1,k,2) + cc(m1,1,k,4) ) &
                     + ( cc(m1,1,k,1) + cc(m1,1,k,3) )
      ch(m2,ido,4,k) = ( cc(m1,1,k,1) + cc(m1,1,k,3) ) &
                      -( cc(m1,1,k,2) + cc(m1,1,k,4) )
      ch(m2,ido,2,k) =   cc(m1,1,k,1) - cc(m1,1,k,3)
      ch(m2,1,3,k) =     cc(m1,1,k,4) - cc(m1,1,k,2)
    end do
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2

    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch(m2,i-1,1,k) = ((wa1(i-2)*cc(m1,i-1,k,2) + wa1(i-1)* &
            cc(m1,i,k,2))+(wa3(i-2)*cc(m1,i-1,k,4) + wa3(i-1)* &
            cc(m1,i,k,4)))+(cc(m1,i-1,k,1)+(wa2(i-2)*cc(m1,i-1,k,3)+ &
            wa2(i-1)*cc(m1,i,k,3)))
          ch(m2,ic-1,4,k) = (cc(m1,i-1,k,1)+(wa2(i-2)*cc(m1,i-1,k,3)+ &
            wa2(i-1)*cc(m1,i,k,3)))-((wa1(i-2)*cc(m1,i-1,k,2)+ &
            wa1(i-1)*cc(m1,i,k,2))+(wa3(i-2)*cc(m1,i-1,k,4)+ &
            wa3(i-1)*cc(m1,i,k,4)))
          ch(m2,i,1,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
            cc(m1,i-1,k,2))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
            cc(m1,i-1,k,4)))+(cc(m1,i,k,1)+(wa2(i-2)*cc(m1,i,k,3)- &
            wa2(i-1)*cc(m1,i-1,k,3)))
          ch(m2,ic,4,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
            cc(m1,i-1,k,2))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
            cc(m1,i-1,k,4)))-(cc(m1,i,k,1)+(wa2(i-2)*cc(m1,i,k,3)- &
            wa2(i-1)*cc(m1,i-1,k,3)))
          ch(m2,i-1,3,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
            cc(m1,i-1,k,2))-(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
            cc(m1,i-1,k,4)))+(cc(m1,i-1,k,1)-(wa2(i-2)*cc(m1,i-1,k,3)+ &
            wa2(i-1)*cc(m1,i,k,3)))
          ch(m2,ic-1,2,k) = (cc(m1,i-1,k,1)-(wa2(i-2)*cc(m1,i-1,k,3)+ &
            wa2(i-1)*cc(m1,i,k,3)))-((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
            cc(m1,i-1,k,2))-(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
            cc(m1,i-1,k,4)))
          ch(m2,i,3,k) = ((wa3(i-2)*cc(m1,i-1,k,4) + wa3(i-1)* &
            cc(m1,i,k,4))-(wa1(i-2)*cc(m1,i-1,k,2) + wa1(i-1)* &
            cc(m1,i,k,2)))+(cc(m1,i,k,1)-(wa2(i-2)*cc(m1,i,k,3)- &
            wa2(i-1)*cc(m1,i-1,k,3)))
          ch(m2,ic,2,k) = ((wa3(i-2)*cc(m1,i-1,k,4) + wa3(i-1)* &
            cc(m1,i,k,4))-(wa1(i-2)*cc(m1,i-1,k,2) + wa1(i-1)* &
            cc(m1,i,k,2)))-(cc(m1,i,k,1)-(wa2(i-2)*cc(m1,i,k,3)- &
            wa2(i-1)*cc(m1,i-1,k,3)))
        end do
      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,ido,1,k) = cc(m1,ido,k,1) &
        + (  hsqt2 * ( cc(m1,ido,k,2) - cc(m1,ido,k,4) ) )
      ch(m2,ido,3,k) = cc(m1,ido,k,1) &
        - (  hsqt2 * ( cc(m1,ido,k,2) - cc(m1,ido,k,4) ) )
      ch(m2,1,2,k) =  -cc(m1,ido,k,3) &
        + ( -hsqt2 * ( cc(m1,ido,k,2) + cc(m1,ido,k,4) ) )
      ch(m2,1,4,k) =   cc(m1,ido,k,3) &
        + ( -hsqt2 * ( cc(m1,ido,k,2) + cc(m1,ido,k,4) ) )
    end do
  end do

  return
end
!========================================================================================
! MRADF5 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mradf5 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3,wa4)

  !*****************************************************************************80
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
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = crm_rknd ) arg
  real ( kind = crm_rknd ) cc(in1,ido,l1,5)
  real ( kind = crm_rknd ) ,intent(inout) :: ch(in2,ido,5,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = crm_rknd ) ti11
  real ( kind = crm_rknd ) ti12
  real ( kind = crm_rknd ) tr11
  real ( kind = crm_rknd ) tr12
  real ( kind = crm_rknd ) wa1(ido)
  real ( kind = crm_rknd ) wa2(ido)
  real ( kind = crm_rknd ) wa3(ido)
  real ( kind = crm_rknd ) wa4(ido)

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2
  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 5.0E+00
  tr11 = cos ( arg )
  ti11 = sin ( arg )
  tr12 = cos ( 2.0E+00 * arg )
  ti12 = sin ( 2.0E+00 * arg )

  do k = 1, l1
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch(m2,1,1,k) = cc(m1,1,k,1) + ( cc(m1,1,k,5) + cc(m1,1,k,2) ) &
                                  + ( cc(m1,1,k,4) + cc(m1,1,k,3) )
      ch(m2,ido,2,k) = cc(m1,1,k,1) + tr11 * ( cc(m1,1,k,5) + cc(m1,1,k,2) ) &
                                    + tr12 * ( cc(m1,1,k,4) + cc(m1,1,k,3) )
      ch(m2,1,3,k) = ti11 * ( cc(m1,1,k,5) - cc(m1,1,k,2) ) &
                   + ti12 * ( cc(m1,1,k,4) - cc(m1,1,k,3) )
      ch(m2,ido,4,k) = cc(m1,1,k,1) + tr12 * ( cc(m1,1,k,5) + cc(m1,1,k,2) ) &
                                    + tr11 * ( cc(m1,1,k,4) + cc(m1,1,k,3) )
      ch(m2,1,5,k) = ti12 * ( cc(m1,1,k,5) - cc(m1,1,k,2) ) &
                   - ti11 * ( cc(m1,1,k,4) - cc(m1,1,k,3) )
    end do
  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+((wa1(i-2)*cc(m1,i-1,k,2)+ &
          wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
          cc(m1,i,k,5)))+((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
          cc(m1,i,k,3))+(wa3(i-2)*cc(m1,i-1,k,4)+ &
          wa3(i-1)*cc(m1,i,k,4)))
        ch(m2,i,1,k) = cc(m1,i,k,1)+((wa1(i-2)*cc(m1,i,k,2)- &
          wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
          cc(m1,i-1,k,5)))+((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
          cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
          cc(m1,i-1,k,4)))
        ch(m2,i-1,3,k) = cc(m1,i-1,k,1)+tr11* &
          ( wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2) &
          +wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5))+tr12* &
          ( wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3) &
          +wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))+ti11* &
          ( wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2) &
          -(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))+ti12* &
          ( wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3) &
          -(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4)))
        ch(m2,ic-1,2,k) = cc(m1,i-1,k,1)+tr11* &
          ( wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2) &
          +wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5))+tr12* &
          ( wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3) &
          +wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))-(ti11* &
          ( wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2) &
          -(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))+ti12* &
          ( wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3) &
          -(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))
        ch(m2,i,3,k) = (cc(m1,i,k,1)+tr11*((wa1(i-2)*cc(m1,i,k,2)- &
          wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
          cc(m1,i-1,k,5)))+tr12*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
          cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
          cc(m1,i-1,k,4))))+(ti11*((wa4(i-2)*cc(m1,i-1,k,5)+ &
          wa4(i-1)*cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
          cc(m1,i,k,2)))+ti12*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
          cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
          cc(m1,i,k,3))))
        ch(m2,ic,2,k) = (ti11*((wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
          cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
          cc(m1,i,k,2)))+ti12*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
          cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
          cc(m1,i,k,3))))-(cc(m1,i,k,1)+tr11*((wa1(i-2)*cc(m1,i,k,2)- &
          wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
          cc(m1,i-1,k,5)))+tr12*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
          cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
          cc(m1,i-1,k,4))))
        ch(m2,i-1,5,k) = (cc(m1,i-1,k,1)+tr12*((wa1(i-2)* &
          cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)* &
          cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5)))+tr11*((wa2(i-2)* &
          cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))+(wa3(i-2)* &
          cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))))+(ti12*((wa1(i-2)* &
          cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa4(i-2)* &
          cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))-ti11*((wa2(i-2)* &
          cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))-(wa3(i-2)* &
          cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))
        ch(m2,ic-1,4,k) = (cc(m1,i-1,k,1)+tr12*((wa1(i-2)* &
          cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)* &
          cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5)))+tr11*((wa2(i-2)* &
          cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))+(wa3(i-2)* &
          cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))))-(ti12*((wa1(i-2)* &
          cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa4(i-2)* &
          cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))-ti11*((wa2(i-2)* &
          cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))-(wa3(i-2)* &
          cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))
        ch(m2,i,5,k) = (cc(m1,i,k,1)+tr12*((wa1(i-2)*cc(m1,i,k,2)- &
          wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
          cc(m1,i-1,k,5)))+tr11*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
          cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
          cc(m1,i-1,k,4))))+(ti12*((wa4(i-2)*cc(m1,i-1,k,5)+ &
          wa4(i-1)*cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
          cc(m1,i,k,2)))-ti11*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
          cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
          cc(m1,i,k,3))))
        ch(m2,ic,4,k) = (ti12*((wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
          cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
          cc(m1,i,k,2)))-ti11*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
          cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
          cc(m1,i,k,3))))-(cc(m1,i,k,1)+tr12*((wa1(i-2)*cc(m1,i,k,2)- &
          wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
          cc(m1,i-1,k,5)))+tr11*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
          cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
          cc(m1,i-1,k,4))))
      end do
    end do
  end do

  return
end
!========================================================================================
! MRADFG is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mradfg (m,ido,ip,l1,idl1,cc,c1,c2,im1,in1,ch,ch2,im2,in2,wa)

  !*****************************************************************************80
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
  real ( kind = crm_rknd ) ,intent(inout) :: ch(in2,ido,l1,ip)
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
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) nbd
  real ( kind = crm_rknd ) tpi
  real ( kind = crm_rknd ) wa(ido)

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2
  tpi = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 )
  arg = tpi / real ( ip, kind = crm_rknd )
  dcp = cos ( arg )
  dsp = sin ( arg )
  ipph = ( ip + 1 ) / 2
  ipp2 = ip + 2
  idp2 = ido + 2
  nbd = ( ido - 1 ) / 2

  if ( ido == 1 ) then
    do ik = 1, idl1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        c2(m1,ik,1) = ch2(m2,ik,1)
      end do
    end do

  else

    do ik = 1, idl1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch2(m2,ik,1) = c2(m1,ik,1)
      end do
    end do

    do j = 2, ip
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch(m2,1,k,j) = c1(m1,1,k,j)
        end do
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
            m2 = m2s
            do m1 = 1, m1d, im1
              m2 = m2 + im2
              ch(m2,i-1,k,j) = wa(idij-1) * c1(m1,i-1,k,j) &
                             + wa(idij)   * c1(m1,i,k,j)
              ch(m2,i,k,j) =   wa(idij-1) * c1(m1,i,k,j)   &
                             - wa(idij)   * c1(m1,i-1,k,j)
            end do
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
            m2 = m2s
            do m1 = 1, m1d, im1
              m2 = m2 + im2
              ch(m2,i-1,k,j) = wa(idij-1) * c1(m1,i-1,k,j) &
                             + wa(idij)   * c1(m1,i,k,j)
              ch(m2,i,k,j) =   wa(idij-1) * c1(m1,i,k,j)   &
                             - wa(idij)   * c1(m1,i-1,k,j)
            end do
          end do
        end do
      end do

    end if

    if ( nbd < l1 ) then

      do j = 2, ipph
        jc = ipp2 - j
        do i = 3, ido, 2
          do k = 1, l1
            m2 = m2s
            do m1 = 1, m1d, im1
              m2 = m2 + im2
              c1(m1,i-1,k,j)  = ch(m2,i-1,k,j)  + ch(m2,i-1,k,jc)
              c1(m1,i-1,k,jc) = ch(m2,i,k,j)    - ch(m2,i,k,jc)
              c1(m1,i,k,j)    = ch(m2,i,k,j)    + ch(m2,i,k,jc)
              c1(m1,i,k,jc)   = ch(m2,i-1,k,jc) - ch(m2,i-1,k,j)
            end do
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ipp2 - j
        do k = 1, l1
          do i = 3, ido, 2
            m2 = m2s
            do m1 = 1, m1d, im1
              m2 = m2 + im2
              c1(m1,i-1,k,j)  = ch(m2,i-1,k,j)  + ch(m2,i-1,k,jc)
              c1(m1,i-1,k,jc) = ch(m2,i,k,j)    - ch(m2,i,k,jc)
              c1(m1,i,k,j)    = ch(m2,i,k,j)    + ch(m2,i,k,jc)
              c1(m1,i,k,jc)   = ch(m2,i-1,k,jc) - ch(m2,i-1,k,j)
            end do
          end do
        end do
      end do

    end if

  end if

  do j = 2, ipph
    jc = ipp2 - j
    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        c1(m1,1,k,j)  = ch(m2,1,k,j)  + ch(m2,1,k,jc)
        c1(m1,1,k,jc) = ch(m2,1,k,jc) - ch(m2,1,k,j)
      end do
    end do
  end do

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ipp2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch2(m2,ik,l)  = c2(m1,ik,1) + ar1 * c2(m1,ik,2)
        ch2(m2,ik,lc) =               ai1 * c2(m1,ik,ip)
      end do
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph
      jc = ipp2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2  = dc2 * ai2 + ds2 * ar2
      ar2 = ar2h
      do ik = 1, idl1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch2(m2,ik,l)  = ch2(m2,ik,l)  + ar2 * c2(m1,ik,j)
          ch2(m2,ik,lc) = ch2(m2,ik,lc) + ai2 * c2(m1,ik,jc)
        end do
      end do
    end do

  end do

  do j = 2, ipph
    do ik = 1, idl1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch2(m2,ik,1) = ch2(m2,ik,1) + c2(m1,ik,j)
      end do
    end do
  end do

  if ( ido < l1 ) then

    do i = 1, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          cc(m1,i,1,k) = ch(m2,i,k,1)
        end do
      end do
    end do

  else

    do k = 1, l1
      do i = 1, ido
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          cc(m1,i,1,k) = ch(m2,i,k,1)
        end do
      end do
    end do

  end if

  do j = 2, ipph
    jc = ipp2 - j
    j2 = j + j
    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        cc(m1,ido,j2-2,k) = ch(m2,1,k,j)
        cc(m1,1,j2-1,k) = ch(m2,1,k,jc)
      end do
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
          m2 = m2s
          do m1 = 1, m1d, im1
            m2 = m2 + im2
            cc(m1,i-1,j2-1,k)  = ch(m2,i-1,k,j) + ch(m2,i-1,k,jc)
            cc(m1,ic-1,j2-2,k) = ch(m2,i-1,k,j) - ch(m2,i-1,k,jc)
            cc(m1,i,j2-1,k)    = ch(m2,i,k,j)   + ch(m2,i,k,jc)
            cc(m1,ic,j2-2,k)   = ch(m2,i,k,jc)  - ch(m2,i,k,j)
          end do
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
          m2 = m2s
          do m1 = 1, m1d, im1
            m2 = m2 + im2
            cc(m1,i-1,j2-1,k)  = ch(m2,i-1,k,j) + ch(m2,i-1,k,jc)
            cc(m1,ic-1,j2-2,k) = ch(m2,i-1,k,j) - ch(m2,i-1,k,jc)
            cc(m1,i,j2-1,k)    = ch(m2,i,k,j)   + ch(m2,i,k,jc)
            cc(m1,ic,j2-2,k)   = ch(m2,i,k,jc)  - ch(m2,i,k,j)
          end do
        end do
      end do
    end do

  end if

  return
end
!========================================================================================
! MRFTB1 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mrftb1 (m,im,n,inc,c,ch,wa,fac)

  !*****************************************************************************80
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
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = crm_rknd ) ,intent(inout) :: c(inc,*)
  real ( kind = crm_rknd ) ,intent(inout) :: ch(m,*)
  real ( kind = crm_rknd ) fac(15)
  real ( kind = crm_rknd ) half
  real ( kind = crm_rknd ) halfm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) im
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m2
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
      m2 = 1 - im
      do i = 1, m
        m2 = m2 + im
        c(m2,j) = half * c(m2,j)
        c(m2,j+1) = halfm * c(m2,j+1)
      end do
    end do

  else

    m2 = 1 - im

    do i = 1, m
      m2 = m2 + im
      ch(i,1) = c(m2,1)
      ch(i,n) = c(m2,n)
    end do

    do j = 2, nl, 2
      m2 = 1 - im
      do i = 1, m
        m2 = m2 + im
        ch(i,j) = half * c(m2,j)
        ch(i,j+1) = halfm * c(m2,j+1)
      end do
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
        call mradb4 ( m, ido, l1, c, im, inc, ch, 1, m, wa(iw), wa(ix2), &
          wa(ix3) )
      else
        call mradb4 ( m, ido, l1, ch, 1, m, c, im, inc, wa(iw), wa(ix2), &
          wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call mradb2 ( m, ido, l1, c, im, inc, ch, 1, m, wa(iw) )
      else
        call mradb2 ( m, ido, l1, ch, 1, m, c, im, inc, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call mradb3 ( m, ido, l1, c, im, inc, ch, 1, m, wa(iw), wa(ix2) )
      else
        call mradb3 ( m, ido, l1, ch, 1, m, c, im, inc, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call mradb5 ( m, ido, l1, c, im, inc, ch, 1, m, wa(iw), wa(ix2), &
          wa(ix3), wa(ix4) )
      else
        call mradb5 ( m, ido, l1, ch, 1, m, c, im, inc, wa(iw), wa(ix2), &
          wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call mradbg ( m, ido, ip, l1, idl1, c, c, c, im, inc, ch, ch, 1, &
          m, wa(iw) )
      else
        call mradbg ( m, ido, ip, l1, idl1, ch, ch, ch, 1, m, c, c, im, &
          inc, wa(iw) )
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
! MRFTF1 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mrftf1 (m,im,n,inc,c,ch,wa,fac)

  !*****************************************************************************80
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
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) im
  integer ( kind = 4 ) n
  integer ( kind = 4 ) inc
  real ( kind = crm_rknd ) ,intent(inout) :: c(inc,*)
  real ( kind = crm_rknd ) ,intent(inout) :: ch(m,*)
  real ( kind = crm_rknd ) ,intent(inout) :: wa(n)
  real ( kind = crm_rknd ) ,intent(inout) :: fac(15)
  integer ( kind = 4 ) i
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
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = crm_rknd ) sn
  real ( kind = crm_rknd ) tsn
  real ( kind = crm_rknd ) tsnm

  nf = int ( fac(2) )   ! number of factors
  na = 1
  l2 = n
  iw = n

! n = nx
! m = lot  = nz

! do i = 1,15
! write(*,*)'whannah - fac(',i,') : '+fac(i)
! end do

  do k1 = 1, nf

    kh = nf - k1
    ip = int ( fac(kh+3) )
! ip = fac(nf-k1+3)
! nf = fac(2)
! fac = wsave(n+1)
! ip = fac( fac(2) - k1 + 3 )

! write(*,222) k1,ip,nf,na,(kh+3),fac(kh+3)
! 222 format('whannah - k1=',i2,' ip=',i2,' nf=',i2,' na=',i2,' fac(',i2,')=',f4.1)

    l1 = l2 / ip
    ido = n / l2
    idl1 = ido * l1
    iw = iw - ( ip - 1 ) * ido
    na = 1 - na

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call mradf4 ( m, ido, l1, c, im, inc, ch, 1,m, wa(iw), wa(ix2), &
          wa(ix3) )
      else
        call mradf4 ( m, ido, l1, ch, 1, m, c, im, inc, wa(iw), wa(ix2), &
          wa(ix3) )
      end if

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call mradf2 ( m, ido, l1, c, im, inc, ch, 1, m, wa(iw) )
      else
        call mradf2 ( m, ido, l1, ch, 1, m, c, im, inc, wa(iw) )
      end if

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call mradf3 ( m, ido, l1, c, im, inc, ch, 1, m, wa(iw), wa(ix2) )
      else
        call mradf3 ( m, ido, l1, ch, 1, m, c, im, inc, wa(iw), wa(ix2) )
      end if

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call mradf5 ( m, ido, l1, c, im, inc, ch, 1, m, wa(iw), wa(ix2), &
          wa(ix3), wa(ix4) )
      else
        call mradf5 ( m, ido, l1, ch, 1, m, c, im, inc, wa(iw), wa(ix2), &
          wa(ix3), wa(ix4) )
      end if

    else

      if ( ido == 1 ) then
        na = 1 - na
      end if

      if ( na == 0 ) then
        call mradfg ( m, ido, ip, l1, idl1, c, c, c, im, inc, ch, ch, 1, &
          m, wa(iw) )
        na = 1
      else
        call mradfg ( m, ido, ip, l1, idl1, ch, ch, ch, 1, m, c, c, im, &
          inc, wa(iw) )
        na = 0
      end if

    end if

    l2 = l1

  end do

  sn  = 1.0E+00 / real ( n, kind = crm_rknd )
  tsn = 2.0E+00 / real ( n, kind = crm_rknd )
  tsnm = -tsn
  modn = mod ( n, 2 )

  if ( modn /= 0 ) then
    nl = n - 1
  else
    nl = n - 2
  end if

  if ( na == 0 ) then

    m2 = 1-im
    do i = 1, m
      m2 = m2 + im
      c(m2,1) = sn * ch(i,1)
    end do

    do j = 2, nl, 2
      m2 = 1 - im
      do i = 1, m
        m2 = m2 + im
        c(m2,j) = tsn * ch(i,j)
        c(m2,j+1) = tsnm * ch(i,j+1)
      end do
    end do

    if ( modn == 0 ) then
      m2 = 1 - im
      do i = 1, m
        m2 = m2 + im
        c(m2,n) = sn * ch(i,n)
      end do
    end if

  else

    m2 = 1-im
    do i = 1, m
      m2 = m2 + im
      c(m2,1) = sn * c(m2,1)
    end do

    do j = 2, nl, 2
      m2 = 1 - im
      do i = 1, m
        m2 = m2 + im
        c(m2,j) = tsn * c(m2,j)
        c(m2,j+1) = tsnm * c(m2,j+1)
      end do
    end do

    if ( modn == 0 ) then

      m2 = 1 - im

      do i = 1, m
        m2 = m2 + im
        c(m2,n) = sn * c(m2,n)
      end do

    end if

  end if

  return
end
!========================================================================================
! MRFTI1 is an FFTPACK5.1 auxilliary function.
!========================================================================================
subroutine mrfti1 (n,wa,fac)

  !*****************************************************************************80
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
  !    Input, integer ( kind = 4 ) N, the number for which factorization and 
  !    other information is needed.
  !
  !    Output, real ( kind = crm_rknd ) WA(N), trigonometric information.
  !
  !    Output, real ( kind = crm_rknd ) FAC(15), factorization information.  FAC(1) is 
  !    N, FAC(2) is NF, the number of factors, and FAC(3:NF+2) are the factors.
  !
  implicit none

  integer ( kind = 4 )     ,intent(in   ) :: n
  real ( kind = crm_rknd ) ,intent(inout) :: wa(n)
  real ( kind = crm_rknd ) ,intent(inout) :: fac(15)

  real ( kind = 8 ) arg
  real ( kind = 8 ) argh
  real ( kind = 8 ) argld
  
  real ( kind = crm_rknd ) fi
  real ( kind = 8 ) tpi
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
  integer ( kind = 4 ) ntryh(4)
  
  

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
      fac(nf+2) = real ( ntry, kind = crm_rknd )
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

  fac(1) = real ( n , kind = crm_rknd )
  fac(2) = real ( nf, kind = crm_rknd )
  tpi = 8.0D+00 * atan ( 1.0D+00 )
  argh = tpi / real ( n, kind = crm_rknd )
  is = 0
  nfm1 = nf - 1
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
      argld = real ( ld, kind = crm_rknd ) * argh
      fi = 0.0E+00
      do ii = 3, ido, 2
        i = i + 2
        fi = fi + 1.0E+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do
      is = is + ido
    end do
    l1 = l2
  end do

  return
end
!========================================================================================
! XERCON checks INC, JUMP, N and LOT for consistency.
!========================================================================================
function xercon ( inc, jump, n, lot )

  !*****************************************************************************80
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
  logical xercon

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
end module fftpack5