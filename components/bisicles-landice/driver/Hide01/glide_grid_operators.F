!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_grid_operators.F90 - part of the Community Ice Sheet Model (CISM)  
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

    !NOTE: The following commented-out code is included in stagthickness.
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
  !      The glam approach works better for calving. 
  !TODO: Add a flag that allows zero-thickness values to be omitted from the gradient (e.g., for flwa and temp).

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

end module glide_grid_operators

!----------------------------------------------------------------------------
