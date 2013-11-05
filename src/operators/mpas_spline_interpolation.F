! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_spline_interpolation
!
!> \brief   MPAS Vector reconstruction module
!> \author  Mark Petersen
!> \date    04/02/13
!> \details 
!> This module provides routines for performing spline interpolation.
!
!-----------------------------------------------------------------------
module mpas_spline_interpolation

  use mpas_kind_types

  implicit none

  private

  public ::   mpas_cubic_spline_coefficients, &
              mpas_interpolate_cubic_spline, &
              mpas_integrate_cubic_spline, &
              mpas_integrate_column_cubic_spline, &
              mpas_interpolate_linear, &
              mpas_test_interpolate

  contains

!***********************************************************************
!
!  routine mpas_cubic_spline_coefficients
!
!> \brief   MPAS Cubic spline coefficients routine
!> \author  Mark Petersen
!> \date    04/02/13
!> \details 
!>  This routine computes second derivatives at nodes.  
!>  This must be run before any of the other cubic spine functions.
!>
!>  Given arrays x(1:n) and y(1:n) containing a function,
!>  i.e., y(i) = f(x(i)), with x monotonically increasing
!>  this routine returns an array y2ndDer(1:n) that contains 
!>  the second derivatives of the interpolating function at x(1:n). 
!>  This routine uses boundary conditions for a natural spline, 
!>  with zero second derivative on that boundary.
!
!-----------------------------------------------------------------------
 subroutine mpas_cubic_spline_coefficients(x,y,n,y2ndDer)  !{{{

! INPUT PARAMETERS:

  integer, intent(in) :: &
    n     !< Input: number of nodes
  real(kind=RKIND), intent(in), dimension(n) :: &
    x,   &!< Input: location of nodes
    y     !< Input: value at nodes

! OUTPUT PARAMETERS:

  real(kind=RKIND), intent(out), dimension(n) :: &
    y2ndDer    !< Output: dy^2/dx^2 at each node

!  local variables:

  integer :: i
  real(kind=RKIND) :: &
    temp,xRatio,a(n)  

   y2ndDer(1)=0.0
   y2ndDer(n)=0.0
   a(1)=0.0

   do i=2,n-1  
      xRatio=(x(i)-x(i-1))/(x(i+1)-x(i-1))  
      temp=1.0/(2.0+xRatio*y2ndDer(i-1))
      y2ndDer(i)=temp*(xRatio-1.0)
      a(i) = temp*(6.0*((y(i+1)-y(i))/(x(i+1)-x(i)) &
          -(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1)) &
          -xRatio*a(i-1)) 
   enddo

   do i=n-1,1,-1  
      y2ndDer(i)=y2ndDer(i)*y2ndDer(i+1)+a(i)  
   enddo

  end subroutine mpas_cubic_spline_coefficients!}}}

!***********************************************************************
!
!  routine mpas_interpolate_cubic_spline
!
!> \brief   MPAS Cubic spline interpolation routine
!> \author  Mark Petersen
!> \date    04/02/13
!> \details 
!>  Given the arrays x(1:n) and y(1:n), which tabulate a function,
!>  and given the array y2ndDer(1:n), which is the output from 
!>  CubicSplineCoefficients above, this routine returns the 
!>  cubic-spline interpolated values of yOut(1:nOut) at xOut(1:nOut).
!>  This subroutine assumes that both x and xOut are monotonically
!>  increasing, and that all values of xOut are within the first and
!>  last values of x.
!
!-----------------------------------------------------------------------
  subroutine mpas_interpolate_cubic_spline( &!{{{
                x,y,y2ndDer,n, &
                xOut,yOut,nOut)  

! INPUT PARAMETERS:

  real (kind=RKIND), dimension(n), intent(in) :: &
    x,         &!< Input: node location, input grid
    y,       &!< Input: interpolation variable, input grid
    y2ndDer     !< Input: 2nd derivative of y at nodes

  real (kind=RKIND), dimension(nOut), intent(in) :: &
    xOut          !< Input: node location, output grid

  integer, intent(in) :: &
    n,      &!< Input: number of nodes, input grid
    nOut       !< Input: number of nodes, output grid

! OUTPUT PARAMETERS:

  real (kind=RKIND), dimension(nOut), intent(out) :: &
    yOut        !< Output: interpolation variable, output grid

!  local variables:

  integer :: &
    kIn, kOut ! counters

  real (kind=RKIND) :: &
    a, b, h

  kOut = 1

  kInLoop: do kIn = 1,n-1

    h = x(kIn+1)-x(kIn)

    do while(xOut(kOut) < x(kIn+1)) 

      a = (x(kIn+1)-xOut(kOut))/h  
      b = (xOut(kOut)-x (kIn) )/h  
      yOut(kOut) = a*y(kIn) + b*y(kIn+1) &
        + ((a**3-a)*y2ndDer(kIn) + (b**3-b)*y2ndDer(kIn+1)) &
         *(h**2)/6.0

      kOut = kOut + 1

      if (kOut>nOut) exit kInLoop

    enddo
  
  enddo kInLoop

end subroutine mpas_interpolate_cubic_spline!}}}

!***********************************************************************
!
!  routine mpas_integrate_cubic_spline
!
!> \brief   MPAS Cubic spline integration routine
!> \author  Mark Petersen
!> \date    04/02/13
!> \details 
!>  Given the arrays x(1:n) and y(1:n), which tabulate a function,
!>  and given the array y2ndDer(1:n), which is the output from 
!>  CubicSplineCoefficients above, this routine returns y_integral,
!>  the integral of y from x1 to x2.  The integration formula was 
!>  created by analytically integrating a cubic spline between each node.
!>  This subroutine assumes that x is monotonically increasing, and
!>  that x1 < x2.
!
!-----------------------------------------------------------------------
subroutine mpas_integrate_cubic_spline(x,y,y2ndDer,n,x1,x2,y_integral)  !{{{

! INPUT PARAMETERS:

  integer, intent(in) :: &
    n     !< Input: number of nodes
  real(kind=RKIND), intent(in), dimension(n) :: &
    x,   &!< Input: location of nodes
    y,   &!< Input: value at nodes
    y2ndDer    !< Input: dy^2/dx^2 at each node
  real(kind=RKIND), intent(in) :: &
    x1,x2 !< Input: limits of integration

! OUTPUT PARAMETERS:

  real(kind=RKIND), intent(out) :: &
    y_integral  !< Output: integral of y

!  local variables:
  
  integer :: i,j,k
  real(kind=RKIND) :: h,h2, A2,B2, F1,F2, eps1

  if (x1<x(1).or.x2>x(n).or.x1>x2) then
    print *, 'error on integration bounds'
  endif

  y_integral = 0.0
  eps1 = 1e-14*x2

  do j=1,n-1  ! loop through sections
    ! section x(j) ... x(j+1)

    if (x2<=x(j)  +eps1) exit
    if (x1>=x(j+1)-eps1) cycle

      h = x(j+1) - x(j)
      h2 = h**2

      ! left side:
      if (x1<x(j)) then
        F1 = -y(j)*h*0.5 + y2ndDer(j)*h**3/24.0
      else
        A2 = (x(j+1)-x1  )**2/h2
        B2 = (x1    -x(j))**2/h2
        F1 = 0.5*h*( -y(j)*A2 + y(j+1)*B2 &
             + y2ndDer(j)  *h2*(-0.5*A2**2 + A2)/6.0 &
             + y2ndDer(j+1)*h2*( 0.5*B2**2 - B2)/6.0 )
      endif

      ! right side:
      if (x2>x(j+1)) then
        F2 = y(j+1)*h*0.5 - y2ndDer(j+1)*h**3/24.0
      else
        A2 = (x(j+1)-x2  )**2/h2
        B2 = (x2    -x(j))**2/h2
        F2 = 0.5*h*( -y(j)*A2 + y(j+1)*B2 &
             + y2ndDer(j)  *h2*(-0.5*A2**2 + A2)/6.0 &
             + y2ndDer(j+1)*h2*( 0.5*B2**2 - B2)/6.0 )
      endif

      y_integral = y_integral + F2 - F1

  enddo ! j

  end subroutine mpas_integrate_cubic_spline!}}}

!***********************************************************************
!
!  routine mpas_integrate_column_cubic_spline
!
!> \brief   MPAS Cubic spline column integration routine
!> \author  Mark Petersen
!> \date    04/02/13
!> \details 
!>  Given the arrays x(1:n) and y(1:n), which tabulate a function,
!>  and given the array y2ndDer(1:n), which is the output from 
!>  CubicSplineCoefficients above, this routine returns 
!>  y_integral(1:nOut), the integral of y.
!>  This is a cumulative integration, so that
!>  y_integral(j) holds the integral of y from x(1) to xOut(j).
!>  The integration formula was created by analytically integrating a 
!>  cubic spline between each node.
!>  This subroutine assumes that both x and xOut are monotonically
!>  increasing, and that all values of xOut are within the first and
!
!-----------------------------------------------------------------------
  subroutine mpas_integrate_column_cubic_spline( &!{{{
               x,y,y2ndDer,n, &
               xOut,y_integral, nOut)  

! INPUT PARAMETERS:

  integer, intent(in) :: &
    n,   &!< Input: number of nodes
    nOut  !< Input: number of output locations to compute integral
  real(kind=RKIND), intent(in), dimension(n) :: &
    x,   &!< Input: location of nodes
    y,   &!< Input: value at nodes
    y2ndDer    !< Input: dy^2/dx^2 at each node
  real(kind=RKIND), dimension(nOut), intent(in) :: &
    xOut  !< Input: output locations to compute integral

! OUTPUT PARAMETERS:

  real(kind=RKIND), dimension(nOut), intent(out) :: &
    y_integral  !< Output: integral from 0 to xOut

!  local variables:

  integer :: i,j,k
  real(kind=RKIND) :: h,h2, A2,B2, F1,F2, eps1

  y_integral = 0.0
  j = 1
  h = x(j+1) - x(j)
  h2 = h**2
  F1 = -y(j)*h*0.5 + y2ndDer(j)*h**3/24.0
  eps1 = 0.0 ! note: could use 1e-12*xOut(nOut)

  k_loop: do k = 1,nOut

    if (k>1) y_integral(k) = y_integral(k-1)

    do while(xOut(k) > x(j+1)-eps1) 
      F2 = y(j+1)*h*0.5 - y2ndDer(j+1)*h**3/24.0
      
      y_integral(k) = y_integral(k) + F2 - F1
      j = j+1
      h = x(j+1) - x(j)
      h2 = h**2
      F1 = -y(j)*h*0.5 + y2ndDer(j)*h**3/24.0
      if (abs(xOut(k) - x(j+1))<eps1) cycle k_loop
    enddo

    A2 = (x(j+1)  - xOut(k))**2/h2
    B2 = (xOut(k) - x(j)   )**2/h2
    F2 = 0.5*h*( -y(j)*A2 + y(j+1)*B2 &
             + y2ndDer(j)  *h2*(-0.5*A2**2 + A2)/6.0 &
             + y2ndDer(j+1)*h2*( 0.5*B2**2 - B2)/6.0 )

    y_integral(k) = y_integral(k) + F2 - F1

    if (k < nOut) then
      A2 = (x(j+1)  -xOut(k))**2/h2
      B2 = (xOut(k) -x(j)   )**2/h2
      F1 = 0.5*h*( -y(j)*A2 + y(j+1)*B2 &
             + y2ndDer(j)  *h2*(-0.5*A2**2 + A2)/6.0 &
             + y2ndDer(j+1)*h2*( 0.5*B2**2 - B2)/6.0 )
    endif

  enddo k_loop

 end subroutine mpas_integrate_column_cubic_spline!}}}

!***********************************************************************
!
!  routine mpas_interpolate_linear
!
!> \brief   MPAS Linear interpolation routine
!> \author  Mark Petersen
!> \date    04/02/13
!> \details 
!>  Given the arrays x(1:n) and y(1:n), which tabulate a function,
!>  this routine returns the linear interpolated values of yOut(1:nOut)
!>  at xOut(1:nOut).
!>  This subroutine assumes that both x and xOut are monotonically
!>  increasing, and that all values of xOut are within the first and
!>  last values of x.
!
!-----------------------------------------------------------------------
 subroutine mpas_interpolate_linear( &!{{{
                x,y,n, &
                xOut,yOut,nOut)  


! !INPUT PARAMETERS:

  real (kind=RKIND), dimension(n), intent(in) :: &
    x,         &!< Input: node location, input grid
    y         !< Input: interpolation variable, input grid

  real (kind=RKIND), dimension(nOut), intent(in) :: &
    xOut          !< Input: node location, output grid

  integer, intent(in) :: &
    N,      &!< Input: number of nodes, input grid
    NOut       !< Input: number of nodes, output grid

! !OUTPUT PARAMETERS:

  real (kind=RKIND), dimension(nOut), intent(out) :: &
    yOut        !< Output: interpolation variable, output grid

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

  integer :: &
    kIn, kOut ! counters

  kOut = 1

  kInLoop: do kIn = 1,n-1

    do while(xOut(kOut) < x(kIn+1)) 

      yOut(kOut) = y(kIn)  &
        + (y(kIn+1)-y(kIn)) &
         /(x(kIn+1)  -x(kIn)  ) &
         *(xOut(kOut)  -x(kIn)  )

      kOut = kOut + 1

      if (kOut>nOut) exit kInLoop

    enddo
  
  enddo kInLoop

  end subroutine mpas_interpolate_linear!}}}

!***********************************************************************
!
!  routine mpas_test_interpolate
!
!> \brief   MPAS Interpolation test routine
!> \author  Mark Petersen
!> \date    04/02/13
!> \details 
!>  Test routine to show how to operate the cubic spline subroutines
!
!-----------------------------------------------------------------------
  subroutine mpas_test_interpolate!{{{

  integer, parameter :: &
    n = 10
  real (kind=RKIND), dimension(n) :: &
    y, x, y2ndDer

  integer, parameter :: &
    nOut = 100
  real (kind=RKIND), dimension(nOut) :: &
    yOut, xOut

  integer :: &
    k

!-----------------------------------------------------------------------
!
!  Create x, y, xOut
!
!-----------------------------------------------------------------------

   do k=1,n
      x(k) = k-4
      ! trig function:
      y(k) = sin(x(k)/2)
   enddo

   do k=1,nOut
      xOut(k) = x(1) + k/(nOut+1.0)*(x(n)-x(1))
   enddo

!-----------------------------------------------------------------------
!
!  Interpolate
!
!-----------------------------------------------------------------------

   ! First, compute second derivative values at each node, y2ndDer.
   call mpas_cubic_spline_coefficients(x,y,n,y2ndDer)

   ! Compute interpolated values yOut.
   call mpas_interpolate_cubic_spline( &
      x,y,y2ndDer,n, &
      xOut,yOut,nOut)

   ! The following output can be copied directly into Matlab
   print *, 'subplot(2,1,1)'
   print '(a,10f8.4,a)', 'x = [',x,'];'
   print '(a,10f8.4,a)', 'y = [',y,'];'
   print '(a,100f8.4,a)', 'xOut = [',xOut,'];'
   print '(a,100f8.4,a)', 'yOut = [',yOut,'];'
   print *, "plot(x,y,'-*r',xOut,yOut,'x')"

   ! Compute interpolated values yOut.
   call mpas_integrate_column_cubic_spline( &
      x,y,y2ndDer,n, &
      xOut,yOut,nOut)  

   ! The following output can be copied directly into Matlab
   print *, 'subplot(2,1,2)'
   print '(a,10f8.4,a)', 'x = [',x,'];'
   print '(a,10f8.4,a)', 'y = 2*cos(-3/2) -2*cos(x/2);'
   print '(a,100f8.4,a)', 'xOut = [',xOut,'];'
   print '(a,100f8.4,a)', 'yOut = [',yOut,'];'
   print *, "plot(x,y,'-*r',xOut,yOut,'x')"

  end subroutine mpas_test_interpolate!}}}

end module mpas_spline_interpolation

