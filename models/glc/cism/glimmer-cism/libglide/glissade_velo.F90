!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!
!TODO - Are all these includes needed?
#ifdef HAVE_CONFIG_H 
#include "config.inc" 
#endif
#include "glide_nan.inc"
#include "glide_mask.inc"

!TODO - What is this?
#define shapedbg(x) write(*,*) "x", shape(x)

module glissade_velo

    use parallel

    ! Driver for glam and glissade higher-order velocity solvers

    use glimmer_global, only : dp
    use glimmer_physcon, only: gn, scyr
    use glimmer_paramets, only: thk0, len0, vel0, vis0
    use glimmer_log
    use glide_types
    use glide_grid_operators, only: stagvarb
    use glide_mask

    use glam_strs2, only: glam_velo_solver, JFNK_velo_solver
    use glissade_velo_higher, only: glissade_velo_higher_solve

    implicit none
    
contains
        
    subroutine glissade_velo_driver(model)

        ! Glissade higher-order velocity driver

        use glide_mask

!!        use glimmer_horiz_bcs, only: horiz_bcs_stag_scalar
        
        type(glide_global_type),intent(inout) :: model

        !For HO masking
        logical :: empty
        integer :: totpts
        real(sp), dimension(model%general%ewn-1, model%general%nsn-1) :: stagmassb

        integer, dimension(model%general%ewn-1, model%general%nsn-1)  :: geom_mask_stag
        real(dp), dimension(model%general%ewn-1, model%general%nsn-1) :: latbc_norms_stag

!WHL - temporary velocity arrays to remove scaling
        real(dp), dimension(model%general%upn, model%general%ewn-1, model%general%nsn-1) ::  &
           uvel, vvel    ! uvel and vvel with scale factor (vel0) removed

!WHL - debug
        integer :: i, j

        !-------------------------------------------------------------------
        ! Velocity prep that is independent of the solver
        !-------------------------------------------------------------------

!TODO - Would this be a good place to compute geometry derivatives?

!TODO - Verify that glide_set_mask works correctly when the input field is on the velo grid.
!       Would be safer to call a set_mask_staggered subroutine?

        !Compute the "geometry mask" (type of square) for the staggered grid

        call glide_set_mask(model%numerics,                                     &
                            model%geomderv%stagthck, model%geomderv%stagtopg,   &
                            model%general%ewn-1,     model%general%nsn-1,       &
                            model%climate%eus,       geom_mask_stag)

!        call stag_parallel_halo ( geom_mask_stag )
!        call horiz_bcs_stag_scalar(geom_mask_stag)


!TODO - What exactly does this do?  Is it solver-dependent?
        !Augment masks with kinematic boundary condition info
!TODO Adding the kinematic bc to thkmask is not needed.  Can be removed.
        call augment_kinbc_mask(model%geometry%thkmask, model%velocity%kinbcmask)
        call augment_kinbc_mask(geom_mask_stag, model%velocity%kinbcmask)

!TODO - Remove this call?  Don't think it is ever used.
        !Compute the normal vectors to the marine margin for the staggered grid
!!        call glide_marine_margin_normal(model%geomderv%stagthck, geom_mask_stag, latbc_norms_stag)

        ! save the final mask to 'dynbcmask' for exporting to netCDF output file
        model%velocity%dynbcmask = geom_mask_stag

        !-------------------------------------------------------------------
        ! Compute the velocity field
        !-------------------------------------------------------------------

        if (model%options%whichdycore == DYCORE_GLAM) then    ! glam finite-difference dycore

           if (model%options%which_ho_nonlinear == HO_NONLIN_PICARD ) then ! Picard (standard solver)

             !TODO - Are all these options still supported?
             !WHL - Removed periodic_ew, periodic_ns

             call t_startf('glam_velo_solver')
              call glam_velo_solver( model%general%ewn,       model%general%nsn,                 &
                                     model%general%upn,                                          &
                                     model%numerics%dew,      model%numerics%dns,                &
                                     model%numerics%sigma,    model%numerics%stagsigma,          &
                                     model%geometry%thck,     model%geometry%usrf,               &
                                     model%geometry%lsrf,     model%geometry%topg,               &
                                     model%geomderv%dthckdew, model%geomderv%dthckdns,           &
                                     model%geomderv%dusrfdew, model%geomderv%dusrfdns,           &
                                     model%geomderv%dlsrfdew, model%geomderv%dlsrfdns,           & 
                                     model%geomderv%stagthck, model%temper%flwa,                 &
                                     model%basalproc%minTauf,                                    & 
                                     model%velocity%btraction,                                   & 
                                     geom_mask_stag,                                             &
                                     model%options%which_ho_babc,                                &
                                     model%options%which_ho_efvs,                                &
                                     model%options%which_ho_resid,                               &
                                     model%options%which_ho_nonlinear,                           &
                                     model%options%which_ho_sparse,                              &
                                     model%velocity%beta,                                        & 
                                     model%velocity%uvel, model%velocity%vvel,                   &
                                     model%velocity%uflx, model%velocity%vflx,                   &
                                     model%stress%efvs )
             call t_stopf('glam_velo_solver')

           else if ( model%options%which_ho_nonlinear == HO_NONLIN_JFNK ) then ! JFNK

!TODO - Create a JFNK solver that can work with an arbitrary calcF routine
!       (e.g., variational as well as Payne-Price)
! noxsolve could eventually go here 

             call t_startf('JFNK_velo_solver')
              call JFNK_velo_solver (model, geom_mask_stag) 
             call t_stopf('JFNK_velo_solver')

           else   
              call write_log('Invalid which_ho_nonlinear option.',GM_FATAL)
           end if

        else  ! glissade finite-element dycore

          !----------------------------------------------------------------
          ! Note: The glissade solver uses SI units.
          ! Thus we have grid cell dimensions and ice thickness in meters,
          !  velocity in m/s, and the rate factor in Pa^(-n) s(-1).
          !----------------------------------------------------------------

           if (model%options%which_ho_nonlinear == HO_NONLIN_PICARD ) then ! Picard (standard solver)

              ! Note: The geometry fields (thck, topg, and usrf) must be updated in halos
              !        before calling glissade_velo_higher-solver.
              !       These updates are done in subroutine glissade_diagnostic_variable_solve
              !        in module glissade.F90.

              print*, ' '
              print*, 'Call glissade_velo_higher_solve'
              print*, 'thk0, thklim =', thk0, model%numerics%thklim
              print*, 'scaled flwa(1,26,19) =', model%temper%flwa(1,26,19)
              print*, 'unscaled flwa(1,26,19) =', model%temper%flwa(1,26,19) *vis0*scyr
 
              ! compute the unscaled velocity (m/yr)
              uvel(:,:,:) = vel0 * model%velocity%uvel(:,:,:)
              vvel(:,:,:) = vel0 * model%velocity%vvel(:,:,:)

              call glissade_velo_higher_solve(model%general%ewn,      model%general%nsn,         &
                                              model%general%upn,                                 &  
                                              model%numerics%sigma,                              &
                                              nhalo,                                             &  
                                              len0 * model%numerics%dew,                         &
                                              len0 * model%numerics%dns,                         &
                                              thk0 * model%geometry%thck,                        &
                                              thk0 * model%geometry%usrf,                        &
                                              thk0 * model%geometry%topg,                        &
                                              real(model%climate%eus,dp),                        &
                                              thk0 * model%numerics%thklim,                      &
                                              vis0 * model%temper%flwa,                   &
!!                                              model%velocity%uvel,  model%velocity%vvel,       &
                                              uvel,                   vvel,                      &
!                                              model%velocity%beta,               &  ! add this one later 
!                                              model%options%which_ho_babc,       &  ! add this one later
                                              model%options%which_ho_efvs,       &
                                              model%options%which_ho_resid,      &
                                              model%options%which_ho_nonlinear,  &
                                              model%options%which_ho_sparse)

              ! rescale the velocity since the rest of the code expects it
              model%velocity%uvel(:,:,:) = uvel(:,:,:) / vel0
              model%velocity%vvel(:,:,:) = vvel(:,:,:) / vel0

           else if ( model%options%which_ho_nonlinear == HO_NONLIN_JFNK ) then ! JFNK

!TODO - Create a JFNK solver that can work with an arbitrary calcF routine
!       (i.e., glissade as well as glam)
! noxsolve could eventually go here 

             call t_startf('JFNK_velo_solver')
!!              call JFNK_velo_solver (model, geom_mask_stag)   ! needs to be generalized to glissade solver
             call t_stopf('JFNK_velo_solver')

           else   
              call write_log('Invalid which_ho_nonlinear option.',GM_FATAL)
           end if

        end if   ! whichdycore

        !-------------------------------------------------------------------
        ! Velocity-related computations that are independent of the solver
        !-------------------------------------------------------------------

        ! compute the velocity norm (for diagnostic output)

        model%velocity%velnorm = sqrt(model%velocity%uvel**2 + model%velocity%vvel**2)

        ! WHL - Copy uvel and vvel to arrays uvel_icegrid and vvel_icegrid.
        !       These arrays have horizontal dimensions (nx,ny) instead of (nx-1,ny-1).
        !       Thus they are better suited for I/O if we have periodic BC,
        !        where the velocity field we are solving for has global dimensions (nx,ny).
        !       Since uvel and vvel are not defined for i = nx or for j = ny, the
        !        uvel_icegrid and vvel_icegrid arrays will have values of zero at these points.
        !       But these are halo points, so when we write netCDF I/O it shouldn't matter;
        !        we should have the correct values at physical points.

        model%velocity%uvel_icegrid(:,:,:) = 0.d0
        model%velocity%vvel_icegrid(:,:,:) = 0.d0

        do j = 1, model%general%nsn-1
           do i = 1, model%general%ewn-1
              model%velocity%uvel_icegrid(:,i,j) = model%velocity%uvel(:,i,j)
              model%velocity%vvel_icegrid(:,i,j) = model%velocity%vvel(:,i,j)             
           enddo
        enddo
        
    end subroutine glissade_velo_driver

end module glissade_velo
