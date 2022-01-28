!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_ground.F90 - part of the Community Ice Sheet Model (CISM)  
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

!TODO - Change module and file names to something more appropriate (glide_calving?)

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"
module glide_ground

  use glide_types
  use glimmer_global, only: dp
  use parallel

  implicit none

contains
!-------------------------------------------------------------------------------  

  subroutine glide_calve_ice(whichcalving,                   &
                             thck,         relx,             &  
                             topg,         mask,             &
                             marine_limit, calving_fraction, &    
                             eus,          calving_thck)

    ! Calve ice according to one of several alternative methods
 
    use glimmer_paramets, only: thk0

    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    !TODO: Change mask to thkmask?  The argument passed in is model%geometry%thkmask.

    integer,                intent(in)    :: whichcalving   !> option for calving law
    real(dp),dimension(:,:),intent(inout) :: thck           !> ice thickness
    real(dp),dimension(:,:),intent(in)    :: relx           !> relaxed bedrock topography
    real(dp),dimension(:,:),intent(in)    :: topg           !> present bedrock topography
    integer, dimension(:,:), intent(in)   :: mask           !> grid type mask
    real(dp), intent(in)                  :: marine_limit   !> lower limit on topography elevation for ice to be present
    real(dp), intent(in) :: calving_fraction                !> fraction of ice lost when calving; used with whichcalving = 2
    real(dp), intent(in) :: eus                             !> eustatic sea level
    real(dp),dimension(:,:),intent(out) :: calving_thck     !> thickness lost due to calving

    integer :: ew,ns

    !---------------------------------------------------------------------
   
    calving_thck(:,:) = 0.d0

    select case (whichcalving)

    case(CALVING_NONE)    ! do nothing

        
    case(CALVING_FLOAT_ZERO) ! set thickness to zero if ice is floating

      where (GLIDE_IS_FLOAT(mask))
         calving_thck = thck
         thck = 0.0d0
      end where

    case(CALVING_FLOAT_FRACTION) ! remove fraction of ice when floating

       !WHL - Changed definition of calving_fraction; now it is the fraction lost
      do ns = 2,size(thck,2)-1
         do ew = 2,size(thck,1)-1
            if (GLIDE_IS_CALVING(mask(ew,ns))) then
!!!               calving_thck(ew,ns) = (1.d0-calving_fraction)*thck(ew,ns)
!!!               thck(ew,ns) =  calving_fraction*thck(ew,ns)
               calving_thck(ew,ns) = calving_fraction * thck(ew,ns)
               thck(ew,ns) =  thck(ew,ns) - calving_thck(ew,ns)
               !mask(ew,ns) = ior(mask(ew,ns), GLIDE_MASK_OCEAN)
            end if
         end do
      end do

      ! if uncomment above mask update, then call parallel_halo(mask)

    case(CALVING_RELX_THRESHOLD) ! Set thickness to zero if relaxed bedrock is below a given depth

       where (relx <= marine_limit+eus)
          calving_thck = thck
          thck = 0.0d0
       end where
       
    case(CALVING_TOPG_THRESHOLD) ! Set thickness to zero at marine edge if present bedrock is below a given level

       where (GLIDE_IS_MARINE_ICE_EDGE(mask) .and. topg < marine_limit+eus)
          calving_thck = thck
          thck = 0.0d0
       end where

    ! Huybrechts grounding line scheme for Greenland initialization

    case(CALVING_HUYBRECHTS)   ! used to be case(7)

       !WHL - Previously, this code assumed that eus and relx have units of meters.
       !      Changed to be consistent with dimensionless thickness units.  
!       if(eus > -80.d0) then
!          where (relx <= 2.d0*eus)
!             calving_thck = thck
!             thck = 0.0d0
!          end where
!       elseif (eus <= -80.d0) then
!          where (relx <= (2.d0*eus - 0.25d0*(eus + 80.d0)**2.d0))
!             calving_thck = thck
!             thck = 0.0d0
!          end where
!       end if
       if (eus*thk0 > -80.d0) then
          where (relx*thk0 <= 2.d0*eus*thk0)
             calving_thck = thck
             thck = 0.0d0
          end where
       elseif (eus*thk0 <= -80.d0) then
          where (relx*thk0 <= (2.d0*eus*thk0 - 0.25d0*(eus*thk0 + 80.d0)**2.d0))
             calving_thck = thck
             thck = 0.0d0
          end where
       end if
       
    end select
    
  end subroutine glide_calve_ice
  
!-------------------------------------------------------------------------
!WHL - Removed subroutine calc_gline_flux
!WHL - Removed functions get_ground_thck, get_ground_line
!WHL - Removed subroutines update_ground_line, set_ground_line, lin_reg_xg
!      (Associated with unsupporting calving cases)
!-------------------------------------------------------------------------

end module glide_ground

!---------------------------------------------------------------------------
