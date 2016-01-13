!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_integrate.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> integration of functions using Romberg integration
module glint_integrate

  use glimmer_global, only : dp

  implicit none

  public :: romberg_int, romberg_int_prms

!-------------------------------------------------------------

contains

!-------------------------------------------------------------

  !> double precision function to perform Romberg Integration on function. 
  !!
  !! This routine is an implementation of ACM algorithm 60, by F. L. Bauer.
  !! (Comm. ACM, vol. 4, issue 6, June 1961).

  recursive real(dp) function romberg_int(fct,lgr,rgr)

    implicit none

    interface 
       !> function to be integrated
       function fct(x)
         use glimmer_global, only : dp
         implicit none
         real(dp), intent(in) :: x !< the argument
         real(dp) :: fct
       end function fct
    end interface
    
    real(dp),intent(in) :: lgr    !< Lower bound
    real(dp),intent(in) :: rgr    !< Upper bound
    integer,parameter :: ord = 6

    real(dp),dimension(ord+1) :: t
    real(dp) :: l,u,m
    integer :: f,h,j,n

    l=rgr-lgr
    t(1)=(fct(lgr)+fct(rgr))/2.0
    n=1

    do h=1,ord
       u=0
       m=l/(2*n)

       do j=1,2*n-1,2
          u=u+fct(lgr+j*m)
       end do

       t(h+1)=((u/n)+t(h))/2.0
       f=1
       
       do j=h,1,-1
          f=f*4
          t(j)=t(j+1)+(t(j+1)-t(j))/(f-1)
       end do

       n=2*n

    end do

    romberg_int=t(1)*l

  end function romberg_int

!-------------------------------------------------------------

  !> double precision function to perform Romberg Integration on function. 
  !!
  !! This routine is an implementation of ACM algorithm 60, by F. L. Bauer.
  !! (Comm. ACM, vol. 4, issue 6, June 1961).

  recursive real(dp) function romberg_int_prms(fct,lgr,rgr,params)

    implicit none

    interface 
       !> the function to be integrated
       function fct(x,params)
         use glimmer_global, only : dp
         implicit none
         real(dp), intent(in) :: x !< the argument
         real(dp), intent(in), dimension(:) :: params !< an array of function parameters
         real(dp) :: fct
       end function fct
    end interface
    
    real(dp),intent(in) :: lgr    !< Lower bound
    real(dp),intent(in) :: rgr    !< Upper bound
    real(dp),intent(in),dimension(:) :: params !< parameters for function
    integer,parameter :: ord = 6

    real(dp),dimension(ord+1) :: t
    real(dp) :: l,u,m
    integer :: f,h,j,n

    l=rgr-lgr
    t(1)=(fct(lgr,params)+fct(rgr,params))/2.0
    n=1

    do h=1,ord
       u=0
       m=l/(2*n)

       do j=1,2*n-1,2
          u=u+fct(lgr+j*m,params)
       end do

       t(h+1)=((u/n)+t(h))/2.0
       f=1
       
       do j=h,1,-1
          f=f*4
          t(j)=t(j+1)+(t(j+1)-t(j))/(f-1)
       end do

       n=2*n

    end do

    romberg_int_prms=t(1)*l

  end function romberg_int_prms

!-------------------------------------------------------------

end module glint_integrate

!-------------------------------------------------------------
