! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History:
! July 2006: John Haynes      - Initial version
! May 2015:  Dustin Swales    - Modified for COSPv2.0
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
module math_lib
  USE COSP_KINDS,     ONLY: wp
  use mod_cosp_error, ONLY: errorMessage
  implicit none

contains
  ! ##########################################################################  
  !                           function PATH_INTEGRAL 
  ! ##########################################################################  
  function path_integral(f,s,i1,i2)
    use m_mrgrnk
    use array_lib
    implicit none
    ! ########################################################################
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
    ! ########################################################################

    ! INPUTS
    real(wp),intent(in), dimension(:) :: &
         f,  & ! Functional values
         s     ! Abscissa values  
    integer, intent(in) :: &
         i1, & ! Index of lower limit
         i2    ! Index of upper limit

    ! OUTPUTS
    real(wp) :: path_integral  
  
    ! Internal variables
    real(wp) :: sumo, deltah, val
    integer :: nelm, j
    integer, dimension(i2-i1+1) :: idx
    real(wp), dimension(i2-i1+1) :: f_rev, s_rev

    nelm = i2-i1+1
    if (nelm > 3) then
       call mrgrnk(s(i1:i2),idx)
       s_rev = s(idx)
       f_rev = f(idx)
       call avint(f_rev(i1:i2),s_rev(i1:i2),(i2-i1+1), &
            s_rev(i1),s_rev(i2), val)
       path_integral = val
    else
       sumo = 0._wp
       do j=i1,i2
          deltah = abs(s(i1+1)-s(i1))
          sumo = sumo + f(j)*deltah
       enddo
       path_integral = sumo
    endif
    
    return
  end function path_integral
  
  ! ##########################################################################
  !                            subroutine AVINT
  ! ##########################################################################
  subroutine avint ( ftab, xtab, ntab, a_in, b_in, result )
    implicit none
    ! ########################################################################  
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
    ! ##########################################################################  

    ! INPUTS
    integer,intent(in) :: &
         ntab    ! Number of elements of [ftab] and [xtab]
    real(wp),intent(in) :: &
         a_in, & ! Lower limit of integration
         b_in    ! Upper limit of integration
    real(wp),intent(in),dimension(ntab) :: &
         ftab, & ! Functional values
         xtab    ! Abscissa value
    
    ! OUTPUTS
    real(wp),intent(out) :: result  ! Approximate value of integral

    ! Internal varaibles
    real(wp) :: a, atemp, b, btemp,ca,cb,cc,ctemp,sum1,syl,term1,term2,term3,x1,x2,x3
    integer  :: i,ihi,ilo,ind
    logical  :: lerror
  
    lerror = .false.
    a = a_in
    b = b_in  
  
    if ( ntab < 3 ) then
       call errorMessage('FATAL ERROR(optics/math_lib.f90:AVINT): Ntab is less than 3.')
       return
    end if
    
    do i = 2, ntab
       if ( xtab(i) <= xtab(i-1) ) then
          lerror = .true.
          exit
       end if
    end do
    
    if (lerror) then
       call errorMessage('FATAL ERROR(optics/math_lib.f90:AVINT): Xtab(i) is not greater than Xtab(i-1).')
       return
    end if
    
!ds    result = 0.0D+00
    result = 0._wp
    
    if ( a == b ) then
       call errorMessage('WARNING(optics/math_lib.f90:AVINT): A=B => integral=0')
       return
    end if
    
    !  If B < A, temporarily switch A and B, and store sign.
    if ( b < a ) then
       syl = b
       b = a
       a = syl
       ind = -1
    else
       syl = a
       ind = 1
    end if
    
    !  Bracket A and B between XTAB(ILO) and XTAB(IHI).
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
    
    !  Carry out approximate integration from XTAB(ILO) to XTAB(IHI).
    sum1 = 0._wp
!ds    sum1 = 0.0D+00
    
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
          ca = 0.5_wp * ( atemp + ca )
          cb = 0.5_wp * ( btemp + cb )
          cc = 0.5_wp * ( ctemp + cc )
!ds          ca = 0.5D+00 * ( atemp + ca )
!ds          cb = 0.5D+00 * ( btemp + cb )
!ds          cc = 0.5D+00 * ( ctemp + cc )
       end if
       
       sum1 = sum1 + ca * ( x2**3 - syl**3 ) / 3._wp &
            + cb * 0.5_wp * ( x2**2 - syl**2 ) + cc * ( x2 - syl )
!ds       sum1 = sum1 + ca * ( x2**3 - syl**3 ) / 3.0D+00 &
!ds            + cb * 0.5D+00 * ( x2**2 - syl**2 ) + cc * ( x2 - syl )
       
       ca = atemp
       cb = btemp
       cc = ctemp
       
       syl = x2
       
    end do

    result = sum1 + ca * ( b**3 - syl**3 ) / 3._wp &
         + cb * 0.5_wp * ( b**2 - syl**2 ) + cc * ( b - syl )
!ds    result = sum1 + ca * ( b**3 - syl**3 ) / 3.0D+00 &
!ds         + cb * 0.5D+00 * ( b**2 - syl**2 ) + cc * ( b - syl )

    !  Restore original values of A and B, reverse sign of integral
    !  because of earlier switch.
    if ( ind /= 1 ) then
       ind = 1
       syl = b
       b = a
       a = syl
       result = -result
    end if
    
    return
  end subroutine avint
  ! ######################################################################################
  ! SUBROUTINE gamma
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
  ! ######################################################################################
  function gamma(x)
    ! Inputs 
    real(wp), intent(in) :: x

    ! Outputs
    real(wp) :: gamma
    
    ! Local variables
    real(wp) :: pi,ga,z,r,gr
    integer  :: k,m1,m  
    
    ! Parameters
    real(wp),dimension(26),parameter :: &
         g = (/1.0,0.5772156649015329, -0.6558780715202538, -0.420026350340952e-1,     &
               0.1665386113822915,-0.421977345555443e-1,-0.96219715278770e-2,              &
               0.72189432466630e-2,-0.11651675918591e-2, -0.2152416741149e-3,                &
               0.1280502823882e-3, -0.201348547807e-4, -0.12504934821e-5, 0.11330272320e-5,  &
               -0.2056338417e-6, 0.61160950e-8,0.50020075e-8, -0.11812746e-8, 0.1043427e-9,   &
               0.77823e-11, -0.36968e-11, 0.51e-12, -0.206e-13, -0.54e-14, 0.14e-14, 0.1e-15/)  
!ds    real(wp),dimension(26),parameter :: &
!ds         g = (/1.0d0,0.5772156649015329d0, -0.6558780715202538d0, -0.420026350340952d-1,     &
!ds               0.1665386113822915d0,-0.421977345555443d-1,-0.96219715278770d-2,              &
!ds               0.72189432466630d-2,-0.11651675918591d-2, -0.2152416741149d-3,                &
!ds               0.1280502823882d-3, -0.201348547807d-4, -0.12504934821d-5, 0.11330272320d-5,  &
!ds               -0.2056338417d-6, 0.61160950d-8,0.50020075d-8, -0.11812746d-8, 0.1043427d-9,   &
!ds               0.77823d-11, -0.36968d-11, 0.51d-12, -0.206d-13, -0.54d-14, 0.14d-14, 0.1d-15/)  
    
    pi = acos(-1._wp)    
    if (x ==int(x)) then
       if (x > 0.0) then
          ga=1._wp
          m1=x-1
          do k=2,m1
             ga=ga*k
          enddo
       else
          ga=1._wp+300
       endif
    else
       if (abs(x) > 1.0) then
          z=abs(x)
          m=int(z)
          r=1._wp
          do k=1,m
             r=r*(z-k)
          enddo
          z=z-m
       else
          z=x
       endif
       gr=g(26)
       do k=25,1,-1
          gr=gr*z+g(k)
       enddo
       ga=1._wp/(gr*z)
       if (abs(x) > 1.0) then
          ga=ga*r
          if (x < 0.0) ga=-pi/(x*ga*sin(pi*x))
       endif
    endif
    gamma = ga
    return
  end function gamma
end module math_lib
