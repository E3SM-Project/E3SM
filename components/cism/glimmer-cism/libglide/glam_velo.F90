!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glam_velo.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!WHL - This code used to be part of glissade_velo.F90, which contained calls to both the glam
!      and glissade velo solvers. I moved the glam part to this module, leaving just the
!      glissade part in glissade_velo.F90.

module glam_velo

  use parallel

  ! Driver for glam higher-order velocity solver

  implicit none
    
contains
        
      subroutine glam_velo_driver(model)

        ! Glissade higher-order velocity driver
        !TODO - Determine how much of the following is needed by glam solver only.

        use glimmer_global, only : dp
        use glimmer_physcon, only: gn, scyr
        use glimmer_paramets, only: thk0, len0, vel0, vis0
        use glimmer_log
        use glide_types
        use glam_strs2, only: glam_velo_solver, JFNK_velo_solver
        use glissade_velo_higher, only: glissade_velo_higher_solve
        
        use glam_grid_operators,  only: glam_geometry_derivs, df_field_2d_staggered
        use glide_grid_operators, only: stagvarb    !TODO - Is this needed?  Seems redundant with df_field_2d_staggered
        use glide_mask
!!        use glimmer_horiz_bcs, only: horiz_bcs_stag_scalar
        
        type(glide_global_type),intent(inout) :: model

        integer, dimension(model%general%ewn-1, model%general%nsn-1)  :: geom_mask_stag
        real(dp), dimension(model%general%ewn-1, model%general%nsn-1) :: latbc_norms_stag

!WHL - temporary velocity arrays to remove scaling
!!        real(dp), dimension(model%general%upn, model%general%ewn-1, model%general%nsn-1) ::  &
!!             uvel, vvel    ! uvel and vvel with scale factor (vel0) removed

        !WHL - debug
!!        integer :: i, j

        !-------------------------------------------------------------------
        ! Velocity prep; compute geometry info.
        !-------------------------------------------------------------------


        !TODO - The next chunk of code needs work.  Several calls are repeated.
        !       We should work out which calls are actually needed.

        ! ------------------------------------------------------------------------ 
        ! Now that geometry (thck, topg, lsrf, usrf) is finalized for the time step, 
        ! calculate derivatives that may be needed for the velocity solve.
        ! ------------------------------------------------------------------------     

        !HALO TODO - Make sure these geometry derivs are computed everywhere they are needed
        !       (all locally owned velocity points?)


        !TODO - The subroutine glam_geometry_derivs calls subroutine stagthickness to compute stagthck.
        !       Similarly for dthckdew/ns and dusrfdew/ns
        !       I don't know why we need to call the next three subroutines as well as glam_geometry_derivs.
        !       This calculation of stagthck differs from that in glam_geometry_derivs which calls stagthickness() 
        !        in glide_grids.F90  Which do we want to use?  
        !        stagthickness() seems to be noisier but there are notes in there about some issue related to margins.

        ! SFP: not sure if these are all needed here or not. Halo updates for usrf and thck are needed in order 
        ! for periodic bcs to work. Otherwise, global halos do not contain correct values and, presumably, the gradients
        ! calculated below are incorrect in and near the global halos.
        ! Calls were added here for other staggered variables (stagusrf, stagtopg, and staglsrf), first providing halo
        ! updates to the non-stag vars, then calc. their stag values. This was done because debug lines show that these
        ! stag fields did not have the correct values in their global halos. This may be ok if they are not used at all 
        ! by the dycores called here, but I added them for consistency. More testing needed to determine if they are
        ! essential or not.

        ! SFP: for consistency, I added these calls, so that all scalars interpolated to the stag mesh
        ! first have had their global halos updated. As w/ above calls to halo updates, these may be better 
        ! placed elsewhere. The only call originally here was the one to calc stagthck.

        !TODO - Should we replace these with calls to df_field_2d_staggered?

        call stagvarb(model%geometry%usrf, model%geomderv%stagusrf,&
                      model%general%ewn,   model%general%nsn)

        call stagvarb(model%geometry%lsrf, model%geomderv%staglsrf,&
                      model%general%ewn,   model%general%nsn)

        call stagvarb(model%geometry%topg, model%geomderv%stagtopg,&
                      model%general%ewn,   model%general%nsn)

        call stagvarb(model%geometry%thck, model%geomderv%stagthck,&    ! SFP: this call was already here. Calls to calc 
                      model%general%ewn,   model%general%nsn)           ! stagusrf, staglsrf, and stagtopg were added


        call df_field_2d_staggered(model%geometry%usrf, &
                                   model%numerics%dew,      model%numerics%dns, &
                                   model%geomderv%dusrfdew, model%geomderv%dusrfdns, &
                                   model%geometry%thck,     model%numerics%thklim )

        call df_field_2d_staggered(model%geometry%thck, &
                                   model%numerics%dew,      model%numerics%dns, &
                                   model%geomderv%dthckdew, model%geomderv%dthckdns, &
                                   model%geometry%thck,     model%numerics%thklim )

        !SFP: W.r.t WHL comment below, I went the other route above - that is, did halo updates for the non-stag
        !fields first, then did the subroutine calls to calc. fields on the unstag mesh. I think this makes sure
        !you are not populating the stag field global halos with bad information that may have been sitting in the 
        !associated non-stag field halos in the case that you forgot to update them. Maybe?

        !TODO - Not sure these are needed.
        !Halo updates required for inputs to glide_stress?
        call staggered_parallel_halo (model%geomderv%dusrfdew)
        !       call horiz_bcs_stag_vector_ew(model%geomderv%dusrfdew)

        call staggered_parallel_halo (model%geomderv%dusrfdns)
        !       call horiz_bcs_stag_vector_ns(model%geomderv%dusrfdns)

        call staggered_parallel_halo (model%geomderv%dthckdew)
        !       call horiz_bcs_stag_vector_ew(model%geomderv%dthckdew)

        call staggered_parallel_halo (model%geomderv%dthckdns)
        !       call horiz_bcs_stag_vector_ns(model%geomderv%dthckdns)

        ! call parallel_halo(model%geometry%thkmask) in earlier glide_set_mask call

        ! Compute the new geometry derivatives for this time step
        ! TODO Merge glam_geometry_derivs with the above calculation.

        !SFP: For some reason, this next call IS needed. It does not affect the results of the periodic ismip-hom test case either
        ! way (that is, if it is active or commented out), or the dome test case. But for some reason, if it is not active, it
        ! messes up both shelf test cases. There must be some important derivs being calculated within this call that are NOT
        ! being explicitly calculated above. 

        ! Compute stagthck, staglsrf, stagtopg, dusrfdew/dns, dthckdew/dns, dlsrfdew/dns, d2thckdew2/dns2, d2usrfdew2/dns2

        call glam_geometry_derivs(model)

        !WHL - This is the end of the geometry calculations that need to be streamlined.

        !TODO - Verify that glide_set_mask works correctly when the input field is on the velo grid.
        !       Would be safer to call a set_mask_staggered subroutine?

        !Compute the "geometry mask" (type of square) for the staggered grid

        call glide_set_mask(model%numerics,                                     &
                            model%geomderv%stagthck, model%geomderv%stagtopg,   &
                            model%general%ewn-1,     model%general%nsn-1,       &
                            model%climate%eus,       geom_mask_stag)

        !        call stag_parallel_halo ( geom_mask_stag )
        !        call horiz_bcs_stag_scalar(geom_mask_stag)

        !Augment masks with kinematic boundary condition info
        call augment_kinbc_mask(model%geometry%thkmask, model%velocity%kinbcmask)
        call augment_kinbc_mask(geom_mask_stag, model%velocity%kinbcmask)

        ! save the final mask to 'dynbcmask' for exporting to netCDF output file
        model%velocity%dynbcmask = geom_mask_stag

        !-------------------------------------------------------------------
        ! Compute the velocity field
        !-------------------------------------------------------------------


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
                                  model%temper%bwat,       model%basalproc%mintauf,           & 
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

      end subroutine glam_velo_driver

end module glam_velo
