!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_grid_operators.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

! Various grid operators for the Glide dycore, including routines for computing gradients
! and switching between staggered and unstaggered grids

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_nan.inc"
#include "glide_mask.inc"

module glide_grid_operators

    use glimmer_global, only : dp
    implicit none

contains

!----------------------------------------------------------------------------

  subroutine glide_geometry_derivs(model)

! Compute geometric quantities needed by the glide dycore:
! stagthck (given thck), along with the gradients
! dusrfdew/dns (given usrf) and dthckdew/dns (given thck).

    use glide_types, only: glide_global_type

    implicit none

    type(glide_global_type), intent(inout) :: model

    ! Interpolate ice thickness to velocity points
 
    call stagvarb(model%geometry% thck, &
                  model%geomderv% stagthck,      &
                  model%general%  ewn,           &
                  model%general%  nsn)

    ! Compute EW and NS gradients in usrf and thck

    call geomders(model%numerics, &
                  model%geometry% usrf, &
                  model%geomderv% stagthck,&
                  model%geomderv% dusrfdew, &
                  model%geomderv% dusrfdns)

    call geomders(model%numerics, &
                  model%geometry% thck, &
                  model%geomderv% stagthck,&
                  model%geomderv% dthckdew, &
                  model%geomderv% dthckdns)

!TODO - Should the following code, which is in stagthickness, be added?

    ! Make sure the derivatives are 0 where stagthck = 0

!     where (model%geomderv%stagthck == 0.d0)
!            model%geomderv%dusrfdew = 0.d0
!            model%geomderv%dusrfdns = 0.d0
!            model%geomderv%dthckdew = 0.d0
!            model%geomderv%dthckdns = 0.d0
!     endwhere

  end subroutine glide_geometry_derivs

!---------------------------------------------------------------

  subroutine stagvarb(ipvr,opvr,ewn,nsn)

  ! Interpolate a scalar variable such as ice thickness from cell centers to cell corners.

  !NOTE: This subroutine, used by the glide SIA dycore, is different from 
  !      stagthickness, which is used by the glam HO dycore.  In stagthickness, zero-thickness 
  !      values are ignored when thickness is averaged over four adjacent grid cells.
  !      In stagvarb, zero-thickness values are included in the average.
  !      The glam approach works better for calving. For now we have kept the old glide
  !      approach in stagvarb for backward compatibility (and because calving is less
  !      important to treat realistically in a shallow-ice model).
  !TODO - Switch to the stagthickness approach?

    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:)  :: ipvr
    
    integer, intent(in) :: ewn,nsn

    opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                             ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0

  end subroutine stagvarb

!----------------------------------------------------------------------------

  subroutine stagvarb_3d(ipvr, opvr, ewn, nsn, upn)
    real(dp), intent(in), dimension(:,:,:) :: ipvr
    real(dp), intent(out), dimension(:,:,:) :: opvr
    integer, intent(in) :: ewn, nsn, upn
    integer :: k

    do k = 1, upn
        call stagvarb(ipvr(k,:,:), opvr(k,:,:), ewn, nsn)
    end do

  end subroutine stagvarb_3d

!----------------------------------------------------------------------------

  subroutine stagvarb_mask(ipvr,opvr,ewn,nsn,geometry_mask)

    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:)  :: ipvr
    
    integer, intent(in) :: ewn,nsn
    integer, intent(in), dimension(:,:) :: geometry_mask
    integer :: ew,ns,n
    real(dp) :: tot

        opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                                 ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0

        do ns = 1,nsn-1
            do ew = 1,ewn-1

                !If any of our staggering points are shelf front, ignore zeros when staggering
                if (any(GLIDE_NO_ICE(geometry_mask(ew:ew+1, ns:ns+1)))) then
                    n = 0
                    tot = 0
    
                    if (GLIDE_HAS_ICE(geometry_mask(ew,ns))) then
                        tot = tot + ipvr(ew,ns)
                        n   = n   + 1
                    end if
                    if (GLIDE_HAS_ICE(geometry_mask(ew+1,ns))) then
                        tot = tot + ipvr(ew+1,ns)
                        n   = n   + 1
                    end if
                    if (GLIDE_HAS_ICE(geometry_mask(ew,ns+1))) then
                        tot = tot + ipvr(ew,ns+1)
                        n   = n   + 1
                    end if
                    if (GLIDE_HAS_ICE(geometry_mask(ew+1,ns+1))) then
                        tot = tot + ipvr(ew+1,ns+1)
                        n   = n   + 1
                    end if
                    if (n > 0) then
                        opvr(ew,ns) = tot/n
                    else
                        opvr(ew,ns) = 0
                    end if
                
                !Standard Staggering
                else
                        opvr(ew,ns) = (ipvr(ew+1,ns) + ipvr(ew,ns+1) + &
                                      ipvr(ew+1,ns+1) + ipvr(ew,ns)) / 4.0d0
                end if
  
        end do
    end do

  end subroutine stagvarb_mask

!----------------------------------------------------------------------------

  subroutine stagvarb_3d_mask(ipvr, opvr, ewn, nsn, upn, geometry_mask)
    real(dp), intent(in), dimension(:,:,:) :: ipvr
    real(dp), intent(out), dimension(:,:,:) :: opvr
    integer, intent(in) :: ewn, nsn, upn
    integer, intent(in), dimension(:,:) :: geometry_mask
    integer :: k

    do k = 1, upn
        call stagvarb_mask(ipvr(k,:,:), opvr(k,:,:), ewn, nsn, geometry_mask)
    end do

  end subroutine stagvarb_3d_mask

!----------------------------------------------------------------------------

  subroutine geomders(numerics,ipvr,stagthck,opvrew,opvrns)

    use glimmer_global, only : dp
    use glide_types, only: glide_numerics

    implicit none 

    type(glide_numerics) :: numerics
    real(dp), intent(out), dimension(:,:) :: opvrew, opvrns
    real(dp), intent(in), dimension(:,:) :: ipvr, stagthck

    real(dp) :: dew2, dns2 
    integer :: ew,ns,ewn,nsn

    ! Obviously we don't need to do this every time,
    ! but will do so for the moment.
    dew2 = 1.d0/(2.0d0 * numerics%dew)
    dns2 = 1.d0/(2.0d0 * numerics%dns)
    ewn=size(ipvr,1)
    nsn=size(ipvr,2)

    do ns=1,nsn-1
       do ew = 1,ewn-1
          if (stagthck(ew,ns) /= 0.0d0) then
             opvrew(ew,ns) = (ipvr(ew+1,ns+1)+ipvr(ew+1,ns)-ipvr(ew,ns)-ipvr(ew,ns+1)) * dew2
             opvrns(ew,ns) = (ipvr(ew+1,ns+1)+ipvr(ew,ns+1)-ipvr(ew,ns)-ipvr(ew+1,ns)) * dns2
          else
             opvrew(ew,ns) = 0.
             opvrns(ew,ns) = 0.
          end if
       end do
    end do
    
  end subroutine geomders

!----------------------------------------------------------------------------

    !TODO - This subroutine is not currently used.  Remove it?

    subroutine unstagger_field_2d(f_stag, f, periodic_x, periodic_y)

    ! Copies a staggered grid onto a nonstaggered grid.  This version
    ! assumes periodic boundary conditions.

        use nan_mod, only : NaN

        real(dp), dimension(:,:), intent(in) :: f_stag
        real(dp), dimension(:,:), intent(out) :: f
        logical, intent(in) :: periodic_x, periodic_y

        real(dp), dimension(4) :: pts

        real(dp) :: s,n

        integer :: i,j, k,i1, i2, j1, j2, ni, nj
        
        ni = size(f, 1)
        nj = size(f, 2)

        do i = 1, size(f, 1)
            do j = 1, size(f, 2)
                s = 0
                n = 0
                
                i1 = i-1
                i2 = i
                
                !If we're unstaggering with periodic boundaries, we cross over to the
                !other side of the domain when we "de-average".  Otherwise, we just ignore
                !the point that's off the domain.
                if (i1 == 0) then
                    if (periodic_y) then
                        i1 = ni - 1
                    else
                        i1 = 1
                    end if
                end if
    
                if (i2 == ni) then
                    if (periodic_y) then
                        i2 = 1
                    else
                        i2 = ni - 1
                    end if
                end if
    
                j1 = j-1
                j2 = j
    
                if (j1 == 0) then
                    if (periodic_y) then
                        j1 = nj - 1
                    else
                        j1 = 1
                    end if
                end if
    
                if (j2 == nj) then
                    if (periodic_x) then
                        j2 = 1
                    else
                        j2 = nj - 1
                    end if
                end if
                
                !Place the points into an array, loop over them, and average
                !all the points that AREN'T NaN.
                pts = (/f_stag(i1, j1), f_stag(i2, j1), f_stag(i1, j2), f_stag(i2, j2)/)
            
                do k=1,4
                    if (.not. (IS_NAN(pts(k)))) then
                        s = s + pts(k)
                        n = n + 1
                    end if
                end do
                
                if (n /= 0) then
                    f(i,j) = s/n
                else
                    f(i,j) = NaN
                end if

                !If the upper left of this location is not a number, then the staggered point must be
                !not a number.  This is to prevent dirichlet boundary conditions from being duplicated.
                if (IS_NAN(f_stag(i2,j2))) then
                    f(i,j) = NaN
                end if
            end do
        end do
    
    end subroutine unstagger_field_2d

!----------------------------------------------------------------------------

    !TODO - This subroutine is not used.  Remove it?

    subroutine unstagger_field_3d(f, f_stag, periodic_x, periodic_y)

        real(dp), dimension(:,:,:) :: f, f_stag
        logical, intent(in) :: periodic_x, periodic_y

        integer :: i

        do i = 1,size(f,1)
            call unstagger_field_2d(f(i,:,:), f_stag(i,:,:), periodic_x, periodic_y)
        end do
        
    end subroutine unstagger_field_3d

!----------------------------------------------------------------------------

    !TODO - This is not used.  Remove it?

    subroutine periodic_boundaries(m, apply_to_x, apply_to_y, nlayers_arg)

      use parallel
        !*FD Applies periodic boundary conditions to a 2D array
        real(dp), dimension(:,:), intent(inout) :: m
        integer :: maxx, maxy
        logical :: apply_to_x, apply_to_y
        integer, optional :: nlayers_arg

        integer :: nlayers 

        if (present(nlayers_arg)) then
            nlayers = nlayers_arg
        else
            nlayers = 1
        end if

        maxx = size(m, 1)
        maxy = size(m, 2)
       
        if (apply_to_x) then
            m( 1 : nlayers, : ) = m( maxx-nlayers*2 + 1 : maxx - nlayers, :)
            m( maxx-nlayers+1 : maxx, : ) = m(nlayers + 1 : nlayers*2, : )
        end if

        if (apply_to_y) then
            m( :, 1 : nlayers ) = m( :, maxy-nlayers*2 + 1 : maxy - nlayers )
            m( :, maxy-nlayers+1 : maxy ) = m( :, nlayers + 1 : nlayers*2 )
        end if
        
        !If both directions are periodic, treat the corners specially.
        if(apply_to_x .and. apply_to_y) then
            m(1:nlayers, 1:nlayers) = m(maxx-nlayers*2+1:maxx-nlayers, maxy-nlayers*2+1:maxy-nlayers)
            m(nlayers+1:2*nlayers, nlayers+1:2*nlayers) = m(maxx-nlayers+1:maxx, maxy-nlayers+1:maxy)
            
            m(1:nlayers, maxy-nlayers+1:maxy) = m(maxx-nlayers*2+1:maxx-nlayers, nlayers+1:2*nlayers)
            m(nlayers+1:2*nlayers, maxy-nlayers*2+1:maxy-nlayers) = m(maxx-nlayers+1:maxx, 1:nlayers)
        end if

        call parallel_velo_halo(m)

    end subroutine periodic_boundaries

!----------------------------------------------------------------------------
    
    !TODO - This is not used.  Remove it?

    subroutine periodic_boundaries_3d(m, apply_to_x, apply_to_y, nlayers_arg)

        !*FD Applies periodic boundary conditions to a 3D array
        real(dp), dimension(:,:,:), intent(inout) :: m
        logical :: apply_to_x, apply_to_y
        integer, optional :: nlayers_arg
    
        integer :: i
        
        do i = 1, size(m,1)
            call periodic_boundaries(m(i,:,:), apply_to_x, apply_to_y, nlayers_arg)
        end do

    end subroutine periodic_boundaries_3d
 
!---------------------------------------------------------------------------------

end module glide_grid_operators

!----------------------------------------------------------------------------
