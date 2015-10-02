!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   advect_2ndmo.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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
! implementation of the conservation of 2nd order moments algorithm
! M. Prather 1987, JGR Vol 91, D6, p 6671-6681

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

#define _S0(i,j)  adv%moments(i,j)%s0
#define _SX(i,j)  adv%moments(i,j)%sx
#define _SXX(i,j)  adv%moments(i,j)%sxx
#define _SY(i,j)  adv%moments(i,j)%sy
#define _SYY(i,j)  adv%moments(i,j)%syy
#define _SXY(i,j)  adv%moments(i,j)%syx
#define _SM(i,j)  adv%moments(i,j)%sm

#define _XF0(i)  adv%mx(i)%s0
#define _XFX(i)  adv%mx(i)%sx
#define _XFXX(i)  adv%mx(i)%sxx
#define _XFY(i)  adv%mx(i)%sy
#define _XFYY(i)  adv%mx(i)%syy
#define _XFXY(i)  adv%mx(i)%syx
#define _XFM(i)  adv%mx(i)%sm

#define _YF0(i)  adv%my(i)%s0
#define _YFX(i)  adv%my(i)%sx
#define _YFXX(i)  adv%my(i)%sxx
#define _YFY(i)  adv%my(i)%sy
#define _YFYY(i)  adv%my(i)%syy
#define _YFXY(i)  adv%my(i)%syx
#define _YFM(i)  adv%my(i)%sm

module advect_2ndmo

  use glimmer_global, only : rk

  private
  public :: advect_type, init_advect, advect, get_species, set_species, finalise_advect

  type moment_type
     real(kind=rk) :: s0 = 0.
     real(kind=rk) :: sx = 0.
     real(kind=rk) :: sxx = 0.
     real(kind=rk) :: sy = 0.
     real(kind=rk) :: syy = 0.
     real(kind=rk) :: syx = 0.
     real(kind=rk) :: sm = 0.
  end type moment_type

  type advect_type
     private
     integer :: nx, ny   !*FD size of grid
     real(kind=rk) :: deltax, deltay
     type(moment_type), dimension(:,:), pointer :: moments => NULL()
     type(moment_type), dimension(:), pointer :: mx => NULL()
     type(moment_type), dimension(:), pointer :: my => NULL()
  end type advect_type

contains
  
  subroutine init_advect(adv, species, deltax, deltay)
    implicit none
    type(advect_type) :: adv
    real(kind=rk), dimension(:,:), intent(in) :: species
    real(kind=rk), intent(in) :: deltax, deltay

    ! set grid spacing
    adv%deltax = deltax
    adv%deltay = deltay
    ! get grid size from size of species array
    adv%nx = size(species,1)
    adv%ny = size(species,2)
    ! allocate memory
    allocate(adv%moments(adv%nx,adv%ny))
    allocate(adv%mx(adv%nx))
    allocate(adv%my(adv%ny))

    ! initialise zeroth order moments (cell average concentration)
    call set_species(adv,species)
    ! initialise total mass contained in each cell
    adv%moments(:,:)%sm = adv%deltax * adv%deltay
  end subroutine init_advect

  subroutine advect(adv, velox, veloy, deltat)
    implicit none
    type(advect_type) :: adv
    real(kind=rk), dimension(:,:), intent(in) :: velox
    real(kind=rk), dimension(:,:), intent(in) :: veloy
    real(kind=rk), intent(in) :: deltat

    ! advect moments in x direction
    call advect_x(adv,velox, deltat)
    ! advect moments in y direction
    call advect_y(adv,veloy, deltat)
  end subroutine advect

  subroutine set_species(adv, species)
    implicit none
    type(advect_type) :: adv
    real(kind=rk), dimension(:,:), intent(in) :: species
    integer i,j

    do j = 2,adv%ny - 1
       do i = 2,adv%nx - 1
          adv%moments(i,j)%s0 = species(i,j)
       end do
    end do
  end subroutine set_species

  subroutine get_species(adv, species)
    implicit none
    type(advect_type) :: adv
    real(kind=rk), dimension(:,:) :: species
    
    species(:,:) = adv%moments(:,:)%s0
  end subroutine get_species

  subroutine finalise_advect(adv)
    implicit none
    type(advect_type) :: adv

    deallocate(adv%moments)
    deallocate(adv%mx)
    deallocate(adv%my)
  end subroutine finalise_advect

  !*****************************************************************************
  ! private procedures
  !*****************************************************************************

  subroutine advect_x(adv,velo, deltat)
    implicit none
    type(advect_type) :: adv
    real(kind=rk), dimension(:,:), intent(in) :: velo
    real(kind=rk), intent(in) :: deltat
    
    ! local variables
    integer i,il,j
    real(kind=rk) :: alf, alfq, alf1, alf1q
    real(kind=rk) :: slpmax, s1max, s1new, s2new
    real(kind=rk) :: temptm

    do j = 2,adv%ny - 1
       ! place limits on appropriate moments to limit flux
       do i = 2,adv%nx - 1
          slpmax = max(real(0.,kind=rk),_S0(i,j))
          s1max = 1.5*slpmax
          s1new = min( s1max, max(-s1max,_SX(i,j)) )
          s2new = min( (slpmax+slpmax-0.3334*abs(s1new)), max(abs(s1new)-slpmax,_SXX(i,j)) )
          _SXY(i,j) = min(slpmax,max(-slpmax,_SXY(i,j)))
          _SX(i,j) = s1new
          _SXX(i,j) = s2new
       end do

       i = 1
       do il = 2,adv%nx - 1
          if (velo(i,j) .ge. 0.) then
             ! flux from (i) to (i+1) when velo greater than 0
             ! create temporary moments/masses for partial cell in transit
             _XFM(i)  = velo(i,j)*deltat*adv%deltay
             alf      = _XFM(i) / _SM(i,j)
             alfq     = alf*alf
             alf1     = 1. - alf
             alf1q    = alf1 * alf1
             _XF0(i)  = alf*(_S0(i,j) +alf1*(_SX(i,j) + (alf1-alf) *_SXX(i,j) ))
             _XFX(i)  = alfq*(_SX(i,j) + 3.*alf1*_SXX(i,j))
             _XFXX(i) = alf*alfq*_SXX(i,j)
             _XFY(i)  = alf*(_SY(i,j) + alf1*_SXY(i,j))
             _XFXY(i) = alfq*_SXY(i,j)
             _XFYY(i) = alf*_SYY(i,j)
             
             ! readjust moments remaining in cell
             _SM(i,j)  = _SM(i,j) - _XFM(i)
             _S0(i,j)  = _S0(i,j) - _XF0(i)
             _SX(i,j)  = alf1q*(_SX(i,j) - 3.*alf*_SXX(i,j))
             _SXX(i,j) = alf1*alf1q*_SXX(i,j)
             _SY(i,j)  = _SY(i,j) - _XFY(i)
             _SYY(i,j) = _SYY(i,j) - _XFYY(i)
             _SXY(i,j) = alf1q*_SXY(i,j)
          else 
             ! flux from (i+1) to (i) when velo less than 0
             _XFM(i) = -velo(i,j)*deltat*adv%deltay
             alf = _XFM(i)/_SM(il,j)
             alfq = alf*alf
             alf1 = 1.-alf
             alf1q = alf1*alf1
             _XF0(i)  = alf*(_S0(il,j) +alf1*(_SX(il,j) + (alf1-alf) *_SXX(il,j) ))
             _XFX(i)  = alfq*(_SX(il,j) + 3.*alf1*_SXX(il,j))
             _XFXX(i) = alf*alfq*_SXX(il,j)
             _XFY(i)  = alf*(_SY(il,j) + alf1*_SXY(il,j))
             _XFXY(i) = alfq*_SXY(il,j)
             _XFYY(i) = alf*_SYY(il,j)
             
             ! readjust moments remaining in cell
             _SM(il,j)  = _SM(il,j) - _XFM(i)
             _S0(il,j)  = _S0(il,j) - _XF0(i)
             _SX(il,j)  = alf1q*(_SX(il,j) - 3.*alf*_SXX(il,j))
             _SXX(il,j) = alf1*alf1q*_SXX(il,j)
             _SY(il,j)  = _SY(il,j) - _XFY(i)
             _SYY(il,j) = _SYY(il,j) - _XFYY(i)
             _SXY(il,j) = alf1q*_SXY(il,j)
          end if
          i = il
       end do

       ! put temporary moments into appropriate neighbouring cells
       i = 1
       do il = 2,adv%nx - 1
          if (velo(i,j) .ge. 0.) then
             ! flux from (i) to (i+1) when velo greater than 0
             _SM(il,j) = _SM(il,j) + _XFM(i)
             alf = _XFM(i)/_SM(il,j)
             alf1 = 1. - alf
             temptm = alf*_S0(il,j) - alf1*_XF0(i)
             _S0(il,j) = _S0(il,j) + _XF0(i)
             _SXX(il,j) = alf*alf*_XFXX(i) + alf1*alf1*_SXX(il,j) &
                  + 5.*(alf*alf1* (_SX(il,j)-_XFX(i)) - (alf1-alf)*temptm)
             _SX(il,j) = alf*_XFX(i) + alf1*_SX(il,j) + 3.*temptm
             _SXY(il,j) = alf*_XFXY(i) + alf1*_SXY(il,j) &
                  + 3.*(-alf1*_XFY(i) + alf*_SY(il,j))
             _SY(il,j) = _SY(il,j) + _XFY(i)
             _SYY(il,j) = _SYY(il,j) + _XFYY(i)
          else
              ! flux from (i+1) to (i) when velo less than 0
             _SM(i,j) = _SM(i,j) + _XFM(i)
             alf = _XFM(i)/_SM(i,j)
             alf1 = 1. - alf
             temptm = alf*_S0(i,j) - alf1*_XF0(i)
             _S0(i,j) = _S0(i,j) + _XF0(i)
             _SXX(i,j) = alf*alf*_XFXX(i) + alf1*alf1*_SXX(i,j) &
                  + 5.*(alf*alf1* (_SX(i,j)-_XFX(i)) - (alf1-alf)*temptm)
             _SX(i,j) = alf*_XFX(i) + alf1*_SX(i,j) + 3.*temptm
             _SXY(i,j) = alf*_XFXY(i) + alf1*_SXY(i,j) &
                  + 3.*(-alf1*_XFY(i) + alf*_SY(i,j))
             _SY(i,j) = _SY(i,j) + _XFY(i)
             _SYY(i,j) = _SYY(i,j) + _XFYY(i)

          end if
          i = il
       end do

    end do
    adv%moments(:,:)%sm = adv%deltax * adv%deltay
  end subroutine advect_x

  subroutine advect_y(adv,velo, deltat)
    implicit none
    type(advect_type) :: adv
    real(kind=rk), dimension(:,:), intent(in) :: velo
    real(kind=rk), intent(in) :: deltat
    
    ! local variables
    integer i,j,jl
    real(kind=rk) :: alf, alfq, alf1, alf1q
    real(kind=rk) :: slpmax, s1max, s1new, s2new
    real(kind=rk) :: temptm

    do i = 2,adv%nx - 1
       ! place limits on appropriate moments to limit flux
       do j = 2,adv%ny - 1
          slpmax = max(real(0.,kind=rk),_S0(i,j))
          s1max = 1.5*slpmax
          s1new = min( s1max, max(-s1max,_SY(i,j)) )
          s2new = min( (slpmax+slpmax-0.3334*abs(s1new)), max(abs(s1new)-slpmax,_SYY(i,j)) )
          _SXY(i,j) = min(slpmax,max(-slpmax,_SXY(i,j)))
          _SY(i,j) = s1new
          _SYY(i,j) = s2new
       end do

       j = 1
       do jl = 2,adv%ny - 1
          if (velo(i,j) .ge. 0.) then
             ! flux from (j) to (j+1) when velo greater than 0
             ! create temporary moments/masses for partial cell in transit
             _YFM(j)  = velo(i,j)*deltat*adv%deltax
             alf      = _YFM(j) / _SM(i,j)
             alfq     = alf*alf
             alf1     = 1. - alf
             alf1q    = alf1 * alf1
             _YF0(j)  = alf*(_S0(i,j) +alf1*(_SY(i,j) + (alf1-alf) *_SYY(i,j) ))
             _YFY(j)  = alfq*(_SY(i,j) + 3.*alf1*_SYY(i,j))
             _YFYY(j) = alf*alfq*_SYY(i,j)
             _YFX(j)  = alf*(_SX(i,j) + alf1*_SXY(i,j))
             _YFXY(j) = alfq*_SXY(i,j)
             _YFXX(j) = alf*_SXX(i,j)
             
             ! readjust moments remaining in cell
             _SM(i,j)  = _SM(i,j) - _YFM(j)
             _S0(i,j)  = _S0(i,j) - _YF0(j)
             _SY(i,j)  = alf1q*(_SY(i,j) - 3.*alf*_SYY(i,j))
             _SYY(i,j) = alf1*alf1q*_SYY(i,j)
             _SX(i,j)  = _SX(i,j) - _YFX(j)
             _SXX(i,j) = _SXX(i,j) - _YFXX(j)
             _SXY(i,j) = alf1q*_SXY(i,j)
          else
             ! flux from (j+1) to (j) when velo less than 0
             _YFM(j) = -velo(i,j)*deltat*adv%deltax
             alf = _YFM(j)/_SM(i,jl)
             alfq = alf*alf
             alf1 = 1.-alf
             alf1q = alf1*alf1
             _YF0(j)  = alf*(_S0(i,jl) +alf1*(_SY(i,jl) + (alf1-alf) *_SYY(i,jl) ))
             _YFY(j)  = alfq*(_SY(i,jl) + 3.*alf1*_SYY(i,jl))
             _YFYY(j) = alf*alfq*_SYY(i,jl)
             _YFX(j)  = alf*(_SX(i,jl) + alf1*_SXY(i,jl))
             _YFXY(j) = alfq*_SXY(i,jl)
             _YFXX(j) = alf*_SXX(i,jl)
             
             ! readjust moments remaining in cell
             _SM(i,jl)  = _SM(i,jl) - _YFM(j)
             _S0(i,jl)  = _S0(i,jl) - _YF0(j)
             _SY(i,jl)  = alf1q*(_SY(i,jl) - 3.*alf*_SYY(i,jl))
             _SYY(i,jl) = alf1*alf1q*_SYY(i,jl)
             _SX(i,jl)  = _SX(i,jl) - _YFX(j)
             _SXX(i,jl) = _SXX(i,jl) - _YFXX(j)
             _SXY(i,jl) = alf1q*_SXY(i,jl)
          end if
          j = jl
       end do          

       ! put temporary moments into appropriate neighbouring cells
       j = 1
       do jl = 2,adv%ny - 1
          if (velo(i,j) .ge. 0.) then
             ! flux from (j) to (j+1) when velo greater than 0
             _SM(i,jl) = _SM(i,jl) + _YFM(j)
             alf = _YFM(j)/_SM(i,jl)
             alf1 = 1. - alf
             temptm = alf*_S0(i,jl) - alf1*_YF0(j)
             _S0(i,jl) = _S0(i,jl) + _YF0(j)
             _SYY(i,jl) = alf*alf*_YFYY(j) + alf1*alf1*_SYY(i,jl) &
                  + 5.*(alf*alf1* (_SY(i,jl)-_YFY(j)) - (alf1-alf)*temptm)
             _SY(i,jl) = alf*_YFY(j) + alf1*_SY(i,jl) + 3.*temptm
             _SXY(i,jl) = alf*_YFXY(j) + alf1*_SXY(i,jl) &
                  + 3.*(-alf1*_YFX(j) + alf*_SX(i,jl))
             _SX(i,jl) = _SX(i,jl) + _YFX(j)
             _SXX(i,jl) = _SXX(i,jl) + _YFXX(j)
          else
              ! flux from (j+1) to (j) when velo less than 0
             _SM(i,j) = _SM(i,j) + _YFM(j)
             alf = _YFM(j)/_SM(i,j)
             alf1 = 1. - alf
             temptm = alf*_S0(i,j) - alf1*_YF0(j)
             _S0(i,j) = _S0(i,j) + _YF0(j)
             _SYY(i,j) = alf*alf*_YFYY(j) + alf1*alf1*_SYY(i,j) &
                  + 5.*(alf*alf1* (_SY(i,j)-_YFY(j)) - (alf1-alf)*temptm)
             _SY(i,j) = alf*_YFY(j) + alf1*_SY(i,j) + 3.*temptm
             _SXY(i,j) = alf*_YFXY(j) + alf1*_SXY(i,j) &
                  + 3.*(-alf1*_YFX(j) + alf*_SX(i,j))
             _SX(i,j) = _SX(i,j) + _YFX(j)
             _SXX(i,j) = _SXX(i,j) + _YFXX(j)

          end if
          j = jl
       end do

    end do
    adv%moments(:,:)%sm = adv%deltax * adv%deltay
  end subroutine advect_y
end module advect_2ndmo
