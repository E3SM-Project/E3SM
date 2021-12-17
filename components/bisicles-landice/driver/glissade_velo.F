!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo.F90 - part of the Community Ice Sheet Model (CISM)  
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

!TODO - Are all these includes needed?
#ifdef HAVE_CONFIG_H 
#include "config.inc" 
#endif
#include "glide_nan.inc"
#include "glide_mask.inc"

module glissade_velo

    use parallel

    ! Driver for Glissade velocity solvers

    implicit none
    
contains
        
    subroutine glissade_velo_driver(model)

      ! Glissade higher-order velocity driver

      use glimmer_global, only : dp
      use glimmer_physcon, only: gn, scyr
      use glimmer_paramets, only: thk0, len0, vel0, vis0, tau0, evs0
      use glimmer_log
      use glide_types
      use glissade_velo_higher, only: glissade_velo_higher_solve
      use glissade_velo_sia, only: glissade_velo_sia_solve
      use glide_mask

      type(glide_global_type),intent(inout) :: model

      integer :: i, j        

      !-------------------------------------------------------------------
      ! Call the velocity solver.
      ! The standard glissade higher-order solver is glissade_velo_higher_solve.
      ! There is an additional local shallow-ice solver, glissade_velo_sia_solve.
      !-------------------------------------------------------------------
      
      if (model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) then
 
         call glissade_velo_sia_solve (model,                                       &
                                       model%general%ewn,      model%general%nsn,   &
                                       model%general%upn)
                       
      else   ! standard higher-order solve
             ! can be BP, L1L2, SSA or SIA, depending on model%options%which_ho_approx

         !-------------------------------------------------------------------
         ! Compute mask for staggered grid. This is needed as an input to calcbeta
         ! (which used to be called here but now is called from glissade_velo_higher_solve).
         ! TODO - Remove the use of stagmask in the Glissade solver?
         !-------------------------------------------------------------------

         call glide_set_mask(model%numerics,                                     &
                             model%geomderv%stagthck, model%geomderv%stagtopg,   &
                             model%general%ewn-1,     model%general%nsn-1,       &
                             model%climate%eus,       model%geometry%stagmask)

         if (model%options%which_ho_nonlinear == HO_NONLIN_PICARD ) then ! Picard (standard solver)

            ! Note: The geometry fields (thck, topg, and usrf) must be updated in halos
            !        before calling glissade_velo_higher_solve.
            !       These updates are done in subroutine glissade_diagnostic_variable_solve
            !        in module glissade.F90.

            call t_startf('glissade_velo_higher_solver')
            call glissade_velo_higher_solve(model,                                             &
                                            model%general%ewn,      model%general%nsn,         &
                                            model%general%upn)
            call t_stopf('glissade_velo_higher_solver')

         else if (model%options%which_ho_nonlinear == HO_NONLIN_JFNK) then

            !TODO - Create a JFNK solver?

            call write_log('JFNK not supported for Glissade velocity solver', GM_FATAL)

         else   

            call write_log('Invalid which_ho_nonlinear option.', GM_FATAL)

         end if  ! which_ho_nonlinear

      endif   ! which_ho_approx

    end subroutine glissade_velo_driver

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module glissade_velo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
