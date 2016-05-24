!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_mpinterp.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glint_mpinterp

  !*FD Uses SLAP to calculate the field needed to perform
  !*FD mean-preserving interpolation on a sphere

  use glimmer_global, only: dp
  use glimmer_physcon, only: pi
  implicit none

  type mpinterp
     integer :: nx,ny !*FD Grid sizes
     integer :: lenw,leniw !*FD Lengths of work arrays
     ! Column vectors
     real(dp),dimension(:),  pointer :: rhs   => null() !*FD Right-hand side
     real(dp),dimension(:),  pointer :: answ  => null() !*FD Answer
     ! Sparse matrix storage
     integer, dimension(:),  pointer :: row   => null() !*FD Row indicies
     integer, dimension(:),  pointer :: col   => null() !*FD Column indices
     real(dp),dimension(:),  pointer :: arr   => null() !*FD Array elements
     ! Work arrays
     real(dp),dimension(:),  pointer :: rwork => null() !*FD Real work array
     integer, dimension(:),  pointer :: iwork => null() !*FD Int work array
     ! Grid-box areas
     real(dp),dimension(:,:),pointer :: areas => null() !*FD Grid-box areas
  end type mpinterp

  private
  public mpinterp,new_mpinterp,mean_preserve_interp

contains

  subroutine new_mpinterp(params,grid)

    use glint_global_grid

    type(mpinterp),   intent(inout) :: params
    type(global_grid),intent(in)    :: grid

    real(dp),dimension(:),allocatable :: dx1,dx2,dx3,dx4
    real(dp),dimension(:),allocatable :: dy1,dy2,dy3,dy4

    real(dp),dimension(:,:),allocatable :: fa,ga,ha,ja
    real(dp),dimension(:,:),allocatable :: fb,gb,hb,jb
    real(dp),dimension(:,:),allocatable :: fc,gc,hc,jc
    real(dp),dimension(:,:),allocatable :: fd,gd,hd,jd
    real(dp),dimension(:,:),allocatable :: xta,xtb,xtc,xtd
    real(dp),dimension(:,:),allocatable :: yta,ytb,ytc,ytd
    real(dp),dimension(:,:),allocatable :: xyta,xytb,xytc,xytd
    real(dp),dimension(:,:),allocatable :: ph1,ph2,ph3,ph4
    real(dp),dimension(:,:),allocatable :: ph5,ph6,ph7,ph8,ph9

    real(dp),dimension(:),allocatable :: lons,lonb,lats,latb

    integer :: i,j,ii

    params%nx = grid%nx
    params%ny = grid%ny

    ! Allocate derived type elements
    if (associated(params%rhs))   deallocate(params%rhs)
    if (associated(params%answ))  deallocate(params%answ)
    if (associated(params%row))   deallocate(params%row)
    if (associated(params%col))   deallocate(params%col)
    if (associated(params%arr))   deallocate(params%arr)
    if (associated(params%rwork)) deallocate(params%rwork)
    if (associated(params%iwork)) deallocate(params%iwork)
    if (associated(params%areas)) deallocate(params%areas)

    allocate(params%rhs  (params%nx*params%ny))
    allocate(params%answ (params%nx*params%ny))
    allocate(params%row  (params%nx*params%ny*9))
    allocate(params%col  (params%nx*params%ny*9))
    allocate(params%arr  (params%nx*params%ny*9))
    allocate(params%areas(params%nx,params%ny))

    ! Allocate temporary work arrays
    ! 1D arrays
    allocate(dx1(params%nx),dx2(params%nx),dx3(params%nx),dx4(params%nx))
    allocate(dy1(params%ny),dy2(params%ny),dy3(params%ny),dy4(params%ny))

    ! 2D arrays
    allocate(fa(params%nx,params%ny),ga(params%nx,params%ny))
    allocate(ha(params%nx,params%ny),ja(params%nx,params%ny))
    allocate(fb(params%nx,params%ny),gb(params%nx,params%ny))
    allocate(hb(params%nx,params%ny),jb(params%nx,params%ny))
    allocate(fc(params%nx,params%ny),gc(params%nx,params%ny))
    allocate(hc(params%nx,params%ny),jc(params%nx,params%ny))
    allocate(fd(params%nx,params%ny),gd(params%nx,params%ny))
    allocate(hd(params%nx,params%ny),jd(params%nx,params%ny))

    allocate(xta(params%nx,params%ny),xtb(params%nx,params%ny))
    allocate(xtc(params%nx,params%ny),xtd(params%nx,params%ny))
    allocate(yta(params%nx,params%ny),ytb(params%nx,params%ny))
    allocate(ytc(params%nx,params%ny),ytd(params%nx,params%ny))
    allocate(xyta(params%nx,params%ny),xytb(params%nx,params%ny))
    allocate(xytc(params%nx,params%ny),xytd(params%nx,params%ny))

    allocate(ph1(params%nx,params%ny),ph2(params%nx,params%ny))
    allocate(ph3(params%nx,params%ny),ph4(params%nx,params%ny))
    allocate(ph5(params%nx,params%ny),ph6(params%nx,params%ny))
    allocate(ph7(params%nx,params%ny),ph8(params%nx,params%ny))
    allocate(ph9(params%nx,params%ny))

    ! Local copy of grid data 
    ! (necessary because we deal in radians)
    allocate(lons(params%nx),lonb(0:params%nx))
    allocate(lats(params%ny),latb(0:params%ny))

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  xls.f90 - part of the Glimmer-CISM ice model             + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
    dy2(params%ny) = latb(params%ny-1)-lats(params%ny)
    dy3(params%ny) = lats(params%ny)  -latb(params%ny)
    dy4(params%ny) = latb(params%ny)  +pi+lats(params%ny)

    ! Quadrant constants
    do i=1,params%nx
       do j=1,params%ny
          ! F constants
          fa(i,j) = dx2(i) * (sin(latb(j)+dy3(j)) - sin(latb(j)))
          fb(i,j) = dx3(i) * (sin(latb(j)+dy3(j)) - sin(latb(j)))
          fc(i,j) = dx3(i) * (sin(lats(j)+dy2(j)) - sin(lats(j)))
          fd(i,j) = dx2(i) * (sin(lats(j)+dy2(j)) - sin(lats(j)))
          ! G constants
          ga(i,j) = (dx2(i)**2/2.0) * (sin(latb(j)+dy3(j)) - sin(latb(j)))
          gb(i,j) = (dx3(i)**2/2.0) * (sin(latb(j)+dy3(j)) - sin(latb(j)))
          gc(i,j) = (dx3(i)**2/2.0) * (sin(lats(j)+dy2(j)) - sin(lats(j)))
          gd(i,j) = (dx2(i)**2/2.0) * (sin(lats(j)+dy2(j)) - sin(lats(j)))
          ! H constants
          ha(i,j) = dx2(i) * (dy3(j)*sin(latb(j)+dy3(j)) + cos(latb(j)+dy3(j)) - cos(latb(j)))
          hb(i,j) = dx3(i) * (dy3(j)*sin(latb(j)+dy3(j)) + cos(latb(j)+dy3(j)) - cos(latb(j)))
          hc(i,j) = dx3(i) * (dy2(j)*sin(lats(j)+dy2(j)) + cos(lats(j)+dy2(j)) - cos(lats(j)))
          hd(i,j) = dx2(i) * (dy2(j)*sin(lats(j)+dy2(j)) + cos(lats(j)+dy2(j)) - cos(lats(j)))
          ! J constants
          ja(i,j) = (dx2(i)**2/2.0) * (dy3(j)*sin(latb(j)+dy3(j)) + cos(latb(j)+dy3(j)) - cos(latb(j)))
          jb(i,j) = (dx3(i)**2/2.0) * (dy3(j)*sin(latb(j)+dy3(j)) + cos(latb(j)+dy3(j)) - cos(latb(j)))
          jc(i,j) = (dx3(i)**2/2.0) * (dy2(j)*sin(lats(j)+dy2(j)) + cos(lats(j)+dy2(j)) - cos(lats(j)))
          jd(i,j) = (dx2(i)**2/2.0) * (dy2(j)*sin(lats(j)+dy2(j)) + cos(lats(j)+dy2(j)) - cos(lats(j)))
       end do
    end do
    params%areas = fa+fb+fc+fd

    do i=1,params%nx
       do j=1,params%ny
          xta(i,j) = (ga(i,j)+fa(i,j)*dx1(i))/(dx1(i)+dx2(i))
          xtb(i,j) = (gb(i,j))               /(dx3(i)+dx4(i))
          xtc(i,j) = (gc(i,j))               /(dx3(i)+dx4(i))
          xtd(i,j) = (gd(i,j)+fd(i,j)*dx1(i))/(dx1(i)+dx2(i))
          
          yta(i,j) = (ha(i,j)+fa(i,j)*dy4(j))/(dy3(j)+dy4(j))
          ytb(i,j) = (hb(i,j)+fb(i,j)*dy4(j))/(dy3(j)+dy4(j))
          ytc(i,j) = (hc(i,j))               /(dy1(j)+dy2(j))
          ytd(i,j) = (hd(i,j))               /(dy1(j)+dy2(j))

          xyta(i,j) = (ja(i,j)+ha(i,j)*dx1(i)+ga(i,j)*dy4(j)+fa(i,j)*dx1(i)*dy4(j))/ &
               ((dx1(i)+dx2(i))*(dy3(j)+dy4(j)))
          xytb(i,j) = (jb(i,j)+gb(i,j)*dy4(j))/((dx3(i)+dx4(i))*(dy3(j)+dy4(j)))
          xytc(i,j) = jc(i,j)                 /((dx3(i)+dx4(i))*(dy1(j)+dy2(j)))
          xytd(i,j) = (jd(i,j)+hd(i,j)*dx1(i))/((dx1(i)+dx2(i))*(dy1(j)+dy2(j)))
       end do
    end do

    ! Calculate PHIs (main part)
    ph1 = fa - xta - yta + xyta
    ph2 = xta - xyta + fb -  xtb - ytb + xytb
    ph3 = xtb - xytb
    ph4 = yta - xyta + fd -  xtd - ytd + xytd
    ph5 = xyta + ytb - xytb + xtd - xytd + fc - xtc - ytc + xytc
    ph6 = xytb + xtc - xytc
    ph7 = ytd - xytd
    ph8 = ytc - xytc + xytd
    ph9 = xytc

    params%row=0
    params%col=0
    params%arr=0.0

    ! Matrix main body
    ii=1
    do i=1,params%nx
       do j=1,params%ny
          params%row(ii) = linloc(i,j,params%nx,params%ny)
          params%col(ii) = linloc(i-1,j+1,params%nx,params%ny)
          params%arr(ii) = ph1(i,j)
          ii=ii+1
          params%row(ii) = linloc(i,j,params%nx,params%ny)
          params%col(ii) = linloc(i,j+1,params%nx,params%ny)
          params%arr(ii) = ph2(i,j)
          ii=ii+1
          params%row(ii) = linloc(i,j,params%nx,params%ny)
          params%col(ii) = linloc(i+1,j+1,params%nx,params%ny)
          params%arr(ii) = ph3(i,j)
          ii=ii+1
          params%row(ii) = linloc(i,j,params%nx,params%ny)
          params%col(ii) = linloc(i-1,j,params%nx,params%ny)
          params%arr(ii) = ph4(i,j)
          ii=ii+1
          params%row(ii) = linloc(i,j,params%nx,params%ny)
          params%col(ii) = linloc(i,j,params%nx,params%ny)
          params%arr(ii) = ph5(i,j)
          ii=ii+1
          params%row(ii) = linloc(i,j,params%nx,params%ny)
          params%col(ii) = linloc(i+1,j,params%nx,params%ny)
          params%arr(ii) = ph6(i,j)
          ii=ii+1
          params%row(ii) = linloc(i,j,params%nx,params%ny)
          params%col(ii) = linloc(i-1,j-1,params%nx,params%ny)
          params%arr(ii) = ph7(i,j)
          ii=ii+1
          params%row(ii) = linloc(i,j,params%nx,params%ny)
          params%col(ii) = linloc(i,j-1,params%nx,params%ny)
          params%arr(ii) = ph8(i,j)
          ii=ii+1
          params%row(ii) = linloc(i,j,params%nx,params%ny)
          params%col(ii) = linloc(i+1,j-1,params%nx,params%ny)
          params%arr(ii) = ph9(i,j)
          ii=ii+1
       end do
    end do

    params%lenw  = 10*params%nx*params%ny+8*params%nx*params%ny
    params%leniw = 10*params%nx*params%ny+4*params%nx*params%ny+12
    allocate(params%rwork(params%lenw),params%iwork(params%leniw))

    deallocate(dx1,dx2,dx3,dx4)
    deallocate(dy1,dy2,dy3,dy4)
    deallocate(fa,ga,ha,ja)
    deallocate(fb,gb,hb,jb)
    deallocate(fc,gc,hc,jc)
    deallocate(fd,gd,hd,jd)
    deallocate(xta,xtb,xtc,xtd)
    deallocate(yta,ytb,ytc,ytd)
    deallocate(xyta,xytb,xytc,xytd)
    deallocate(ph1,ph2,ph3,ph4)
    deallocate(ph5,ph6,ph7,ph8,ph9)
    deallocate(lons,lonb,lats,latb)

  end subroutine new_mpinterp

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_preserve_interp(params,in,out,zeros)

    type(mpinterp),         intent(inout) :: params
    real(dp),dimension(:,:),intent(in)    :: in
    real(dp),dimension(:,:),intent(out)   :: out
    logical, dimension(:,:),intent(out)   :: zeros
    
    integer :: i,j,ii

    integer :: iter,ierr
    real(dp) :: err

    ! Copy right-hand side
    do i=1,params%nx
       do j=1,params%ny
          params%rhs(linloc(i,j,params%nx,params%ny)) = in(i,j)*params%areas(i,j)
       end do
    end do

    params%answ = params%rhs

    call dslucs(params%nx*params%ny,   &  ! n  ... order of matrix a (in)
                params%rhs,            &  ! b  ... right hand side vector (in)
                params%answ,           &  ! x  ... initial quess/final solution vector (in/out)
                params%nx*params%ny*9, &  ! nelt ... number of non-zeroes in A (in)
                params%row,            &  ! ia  ... sparse matrix format of A (in)
                params%col,            &  ! ja  ... sparse matrix format of A (in)
                params%arr,            &  ! a   ... matrix (in)
                0,                     &  ! isym ... storage method (0 is complete) (in)
                2,                     &  ! itol ... convergence criteria (2 recommended) (in)
                1.0d-12,               &  ! tol  ... criteria for convergence (in)
                101,                   &  ! itmax ... maximum number of iterations (in)
                iter,                  &  ! iter  ... returned number of iterations (out)
                err,                   &  ! err   ... error estimate of solution (out)
                ierr,                  &  ! ierr  ... returned error message (0 is ok) (out)
                0,                     &  ! iunit ... unit for error writes during iteration (0 no write) (in)
                params%rwork,          &  ! rwork ... workspace for SLAP routines (in)
                params%lenw,           &  ! lenw
                params%iwork,          &  ! iwork ... workspace for SLAP routines (in)
                params%leniw)             ! leniw

    ! Rejig answer back into 2d array 
    do i=1,params%nx
       do j=1,params%ny
          out(i,j)=params%answ(linloc(i,j,params%nx,params%ny))
       end do
    end do

    ! Check for zeros
    zeros = (in==0.0)

  end subroutine mean_preserve_interp

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function linloc(i,j,nx,ny)

    integer :: i,j
    integer :: nx,ny
    integer :: ii,jj

    ii=i
    jj=j

    do
       if (jj>=1.and.jj<=ny) exit
       if (jj==0) then
          jj=1
          ii=ii+nx/2
       end if
       if (jj==ny+1) then
          jj=ny
          ii=ii+nx/2
       end if
    end do

    do
       if (ii>=1.and.ii<=nx) exit
       if (ii<1) ii=ii+nx
       if (ii>nx) ii=ii-nx
    end do

    linloc = (ii-1)*ny+jj

  end function linloc

end module glint_mpinterp
