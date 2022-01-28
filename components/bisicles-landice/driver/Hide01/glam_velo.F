!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glam_velo.F90 - part of the Community Ice Sheet Model (CISM)  
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

!
#ifdef HAVE_CONFIG_H 
#include "config.inc" 
#endif
#include "glide_nan.inc"
#include "glide_mask.inc"

!NOTE - What is shapedbg?
#define shapedbg(x) write(*,*) "x", shape(x)

module glam_velo

  use parallel
  use glimmer_global, only : dp

  ! Driver for glam higher-order velocity solver

  implicit none
  
  private
  public :: glam_velo_driver, glam_basal_friction

contains
        
  subroutine glam_velo_driver(model)

        ! Glissade higher-order velocity driver

        use glimmer_log
        use glide_types
        use glam_strs2, only: glam_velo_solver, JFNK_velo_solver
        use glam_grid_operators,  only: glam_geometry_derivs, df_field_2d_staggered
        use glide_grid_operators, only: stagvarb
        use glide_mask
        use glide_stress
        use glimmer_paramets, only: tau0, vel0
        use glimmer_physcon, only: scyr

        type(glide_global_type),intent(inout) :: model

        logical, parameter :: verbose_glam_velo = .false.
        integer :: i, j, k

        !-------------------------------------------------------------------
        ! Velocity prep; compute geometry info.
        !-------------------------------------------------------------------

        !NOTE - The next chunk of code needs work.  Several calls are repeated.
        !       We should work out which calls are actually needed.

        ! ------------------------------------------------------------------------ 
        ! Now that geometry (thck, topg, lsrf, usrf) is finalized for the time step, 
        ! calculate derivatives that may be needed for the velocity solve.
        ! ------------------------------------------------------------------------     

        !NOTE - Make sure these geometry derivs are computed everywhere they are needed
        !       (all locally owned velocity points?)


        !NOTE - The subroutine glam_geometry_derivs calls subroutine stagthickness to compute stagthck.
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

        !NOTE - Should we replace these with calls to df_field_2d_staggered?

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

        !NOTE - Not sure halo updates are needed for dusrfdew, etc.
        !Halo updates required for inputs to glide_stress?
        call staggered_parallel_halo (model%geomderv%dusrfdew)
        call staggered_parallel_halo (model%geomderv%dusrfdns)
        call staggered_parallel_halo (model%geomderv%dthckdew)
        call staggered_parallel_halo (model%geomderv%dthckdns)
        ! call parallel_halo(model%geometry%thkmask) in earlier glide_set_mask call

        ! Compute the new geometry derivatives for this time step
        ! NOTE Merge glam_geometry_derivs with the above calculation.

        !SFP: For some reason, this next call IS needed. It does not affect the results of the periodic ismip-hom test case either
        ! way (that is, if it is active or commented out), or the dome test case. But for some reason, if it is not active, it
        ! messes up both shelf test cases. There must be some important derivs being calculated within this call that are NOT
        ! being explicitly calculated above. 

        ! Compute stagthck, staglsrf, stagtopg, dusrfdew/dns, dthckdew/dns, dlsrfdew/dns, d2thckdew2/dns2, d2usrfdew2/dns2

        call glam_geometry_derivs(model)

        !WHL - This is the end of the geometry calculations that need to be streamlined.

        !NOTE - Verify that glide_set_mask works correctly when the input field is on the velo grid.
        !       Would be safer to call a set_mask_staggered subroutine?

        !Compute the "geometry mask" (type of square) for the staggered grid

        call glide_set_mask(model%numerics,                                     &
                            model%geomderv%stagthck, model%geomderv%stagtopg,   &
                            model%general%ewn-1,     model%general%nsn-1,       &
                            model%climate%eus,       model%geometry%stagmask)

        !        call stag_parallel_halo (model%geometry%stagmask)

        !Augment masks with kinematic boundary condition info
        call augment_kinbc_mask(model%geometry%thkmask, model%velocity%kinbcmask)
        call augment_kinbc_mask(model%geometry%stagmask, model%velocity%kinbcmask)

        ! save the final mask to 'dynbcmask' for exporting to netCDF output file
        model%velocity%dynbcmask = model%geometry%stagmask

        !-------------------------------------------------------------------
        ! Compute the velocity field
        !-------------------------------------------------------------------

        if (model%options%which_ho_nonlinear == HO_NONLIN_PICARD ) then ! Picard (standard solver)

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
                                  model%velocity%btraction,                                   & 
                                  model%geometry%stagmask,                                    &
                                  model%options%which_ho_babc,                                &
                                  model%options%which_ho_efvs,                                &
                                  model%options%which_ho_resid,                               &
                                  model%options%which_ho_nonlinear,                           &
                                  model%options%which_ho_sparse,                              &
                                  model%velocity%beta_internal,                               &   ! beta weighted by f_ground
                                  model%velocity%beta,                                        &   ! fixed, external beta 
                                  model%velocity%ho_beta_const,                               &
                                  model%basalproc%mintauf,                                    &
                                  model%temper%bwat,                                          &
                                  model%basal_physics,                                        &
                                  model%velocity%uvel, model%velocity%vvel,                   &
                                  model%velocity%uflx, model%velocity%vflx,                   &
                                  model%stress%efvs )
           call t_stopf('glam_velo_solver')

        else if ( model%options%which_ho_nonlinear == HO_NONLIN_JFNK ) then ! JFNK

           ! noxsolve could eventually go here 
           !NOTE - Remove model%geometry%stagmask from argument list; just pass in model
           !       (model%geometry%stagmask used to be called geom_mask_stag, which was not part of model derived type)

           call t_startf('JFNK_velo_solver')
           call JFNK_velo_solver (model, model%geometry%stagmask) 
           call t_stopf('JFNK_velo_solver')

        else   
           call write_log('Invalid which_ho_nonlinear option.',GM_FATAL)
        end if    ! which_ho_nonlinear

        ! Compute internal stresses
        call glide_calcstrsstr(model)

        !WHL - debug - output internal stresses and velocity at a diagnostic point
        if (verbose_glam_velo .and. this_rank==model%numerics%rdiag_local) then
           i = model%numerics%idiag_local
           j = model%numerics%jdiag_local
           print*, ' '
           print*, ' '
           print*, 'i, j =', i, j
           print*, 'k, tau_xz, tau_yz, tau_xx, tau_yy, tau_xy, tau_eff:'
           do k = 1, model%general%upn-1
              print*, k, tau0*model%stress%tau%xz(k,i,j), tau0*model%stress%tau%yz(k,i,j), &
                         tau0*model%stress%tau%xx(k,i,j), tau0*model%stress%tau%yy(k,i,j), &
                         tau0*model%stress%tau%xy(k,i,j), tau0*model%stress%tau%scalar(k,i,j) 
           enddo
           print*, 'New velocity: rank, i, j =', this_rank, i, j
           print*, 'k, uvel, vvel:'
           do k = 1, model%general%upn
              print*, k, vel0*scyr*model%velocity%uvel(k,i,j), vel0*scyr*model%velocity%vvel(k,i,j)
           enddo
        endif

      end subroutine glam_velo_driver

!=======================================================================

  subroutine glam_basal_friction (ewn,        nsn,             &
                                  ice_mask,   floating_mask,   &
                                  ubas,       vbas,            &
                                  btraction,  bfricflx)

    ! Compute frictional heat source due to sliding at the bed
    ! Based on a subroutine that used to be in glissade_temp.F90
    !  but now is used only by Glam

    use glimmer_paramets, only: vel0, vel_scale

    !-----------------------------------------------------------------
    ! Input/output arguments
    !-----------------------------------------------------------------

    integer, intent(in) :: ewn, nsn         ! grid dimensions
    integer, dimension(:,:), intent(in) ::   &
         ice_mask,      & ! = 1 if thck > thklim, else = 0
         floating_mask    ! = 1 if ice is floating, else = 0
    real(dp), dimension(:,:), intent(in) :: ubas, vbas   ! basal velocity
    real(dp), dimension(:,:,:), intent(in) :: btraction  ! basal traction
    real(dp), dimension(:,:), intent(out) :: bfricflx    ! basal friction heat flux (W m-2)

    !-----------------------------------------------------------------
    ! Local arguments
    !-----------------------------------------------------------------

    real(dp) :: slterm       ! sliding friction
    integer :: ew, ns, i, j
    integer :: slide_count   ! number of neighbor cells with nonzero sliding

    bfricflx(:,:) = 0.d0

    ! compute heat source due to basal friction
    ! Note: slterm and bfricflx are defined to be >= 0

    do ns = 2, nsn-1
       do ew = 2, ewn-1

          slterm = 0.d0
          slide_count = 0

          ! Note: btraction is computed in glam_strs2.F90

          !WHL - Using thklim instead of thklim_temp because ice thinner than thklim
          !      is assumed to be at rest.

          if (ice_mask(ew,ns)==1 .and. floating_mask(ew,ns)==0) then
             do j = ns-1,ns
                do i = ew-1,ew

                   !SCALING - WHL: Multiplied ubas by vel0/vel_scale so we get the same result in these two cases:
                   !           (1) With scaling:     vel0 = vel_scale = 500/scyr, and ubas is non-dimensional
                   !           (2) Without scaling:  vel0 = 1, vel_scale = 500/scyr, and ubas is in m/s.

!!!                   if (abs(ubas(i,j)) > 1.0d-6 .or.   &
!!!                       abs(vbas(i,j)) > 1.0d-6) then
                   if ( abs(ubas(i,j))*(vel0/vel_scale) > 1.0d-6 .or.   &
                        abs(vbas(i,j))*(vel0/vel_scale) > 1.0d-6 ) then
                      slide_count = slide_count + 1
                      slterm = slterm + btraction(1,i,j)*ubas(i,j) + btraction(2,i,j)*vbas(i,j) 
                   end if
                end do
             end do

          endif  ! ice_mask = 1, floating_mask = 0

          ! include sliding contrib only if temperature node is surrounded by sliding velo nodes
          !NOTE - This may result in non-conservation of energy. 

          if (slide_count == 4) then
             slterm = 0.25d0 * slterm
          else
             slterm = 0.0d0
          end if

          bfricflx(ew,ns) = slterm

       enddo    ! ns
    enddo       ! ew

  end subroutine glam_basal_friction

!===============================================================================

end module glam_velo

!===============================================================================
