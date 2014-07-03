! WJS (1-30-12): The following (turning optimization off) is included as a workaround for
! LONG (infinite???) compile times with xlf, at least in IBM XL Fortran for AIX, V12.1 on bluefire
#ifdef CPRIBM
@PROCESS OPT(0)
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_thck.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#include "glide_nan.inc"

module glide_thck

  use glimmer_global, only : dp
  use glide_types
  use glimmer_sparse
  use glimmer_sparse_type

  !DEBUG only
!!  use xls

  !TODO - Remove oldglide when code comparisons are complete
  use glimmer_paramets, only: oldglide

  implicit none

  private

  public :: init_thck, thck_nonlin_evolve, thck_lin_evolve, stagleapthck, &
            glide_thck_index, glide_calclsrf

  ! debugging Picard iteration
  integer, private, parameter :: picard_unit=101
  real(dp),private, parameter :: picard_interval=500.d0
  integer, private            :: picard_max=0

contains

  subroutine init_thck(model)

    !*FD initialise work data for ice thickness evolution
    use glimmer_log
    use glimmer_paramets, only: GLC_DEBUG
    implicit none
    type(glide_global_type) :: model

      ! Removed this messy array
!!    model%solver_data%fc2 = (/ model%numerics%alpha * model%numerics%dt / (2.0d0 * model%numerics%dew * model%numerics%dew), &
!!                               model%numerics%dt,                                   &
!!                               (1.0d0 - model%numerics%alpha) / model%numerics%alpha, &
!!                               1.0d0 / model%numerics%alpha,                        &
!!                               model%numerics%alpha * model%numerics%dt / (2.0d0 * model%numerics%dns * model%numerics%dns), &
!!                               0.0d0 /) 

    ! WJS: The following code has syntax errors; simply commenting it out for now
    ! if (GLC_DEBUG) then
    !    call write_log('Logging Picard iterations')
    !    if (main_task) then
    !       open(picard_unit,name='picard_info.data',status='unknown')
    !       write(picard_unit,*) '#time    max_iter'
    !    end if
    ! end if

!TODO -  Make sure the arrays allocated here are deallocated at the end of the run.
!        Might want to move allocation/deallocation to subroutines in glide_types.

    ! allocate memory for ADI scheme

    if (model%options%whichevol == EVOL_ADI) then
       allocate(model%thckwk%alpha(max(model%general%ewn, model%general%nsn)))
       allocate(model%thckwk%beta (max(model%general%ewn, model%general%nsn)))
       allocate(model%thckwk%gamma(max(model%general%ewn, model%general%nsn)))
       allocate(model%thckwk%delta(max(model%general%ewn, model%general%nsn)))
    end if

  end subroutine init_thck

!---------------------------------------------------------------------------------

  subroutine thck_lin_evolve(model,newtemps)

    !*FD this subroutine solves the linearised ice thickness equation by computing the
    !*FD diffusivity from quantities of the previous time step

    use glide_velo
    use glimmer_paramets, only: GLC_DEBUG
    use glide_grid_operators, only: glide_geometry_derivs

    implicit none

    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A

    if (model%geometry%empty) then

       model%geometry%thck = dmax1(0.0d0, model%geometry%thck + model%climate%acab * model%numerics%dt)
       if (GLC_DEBUG) then
          print *, "* thck empty - net accumulation added", model%numerics%time
       end if
    else

       !Note: glide_geometry_derivs is called at the beginning of glide_tstep_p1,
       !      and the geometry has not changed, so stagthck and the geometry
       !      derivatives are still up to date.  A call might be needed here
       !      if glide_tstep_p2 were called out of order.

!!       call glide_geometry_derivs(model)

       ! calculate basal velos
       if (newtemps) then

          call slipvelo(model,                &
                        1,                             &
                        model%velocity% btrc,          &
                        model%velocity% ubas,          &
                        model%velocity% vbas)

          ! calculate Glen's A if necessary
          call velo_integrate_flwa(model%velowk,  &
                                   model%geomderv%stagthck,   &
                                   model%temper%flwa)

       end if

       call slipvelo(model,                &
                     2,                             &
                     model%velocity% btrc,          &
                     model%velocity% ubas,          &
                     model%velocity% vbas)

       ! calculate diffusivity

       call velo_calc_diffu(model%velowk,            model%geomderv%stagthck,  &
                            model%geomderv%dusrfdew, model%geomderv%dusrfdns,  &
                            model%velocity%diffu)

        ! get new thicknesses

        call thck_evolve(model,    &
                         model%velocity%diffu, model%velocity%diffu, &
                         .true.,   &
                         model%geometry%thck,  model%geometry%thck)

!--- MJH: Since the linear evolution uses a diffusivity based on the old geometry, the
!    velocity calculated here will also be based on the old geometry.  If it is
!    desired to calculate a velocity for the new evolved geometry, then the diffusivity
!    and other things need to be updated before calculating velocity (commented out with !* ).
!    If using this block starting with !* , delete the call to slipvelo with option 3 below.
!*       ! Update geometry information for new thickness before calculating velocity
!*       model%geometry%usrf = model%geometry%lsrf + model%geometry%thck  ! usrf needed for slope calculations in geometry_derivs
!*       call geometry_derivs(model)  !this updates stagthck and the slopes
!*       call velo_calc_diffu(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
!*            model%geomderv%dusrfdns,model%velocity%diffu)
!*       call slipvelo(model,                &
!*            0,                             &
!*            model%velocity% btrc,          &
!*            model%velocity% ubas,          &
!*            model%velocity% vbas)
!----

       ! calculate horizontal velocity field
       ! (These calls must appear after thck_evolve, as thck_evolve uses ubas,
       ! which slipvelo mutates)

       call slipvelo(model,                         &
                     3,                             &
                     model%velocity%btrc,           &
                     model%velocity%ubas,           &
                     model%velocity%vbas)

       call velo_calc_velo(model%velowk,            model%geomderv%stagthck,   &
                           model%geomderv%dusrfdew, model%geomderv%dusrfdns,   &
                           model%temper%flwa,       model%velocity%diffu,      &
                           model%velocity%ubas,     model%velocity%vbas,       &
                           model%velocity%uvel,     model%velocity%vvel,       &
                           model%velocity%uflx,     model%velocity%vflx,&
                           model%velocity%velnorm)

    end if

  end subroutine thck_lin_evolve

!---------------------------------------------------------------------------------

  subroutine thck_nonlin_evolve(model,newtemps)

    !*FD this subroutine solves the ice thickness equation by doing an outer, 
    !*FD non-linear iteration to update the diffusivities and in inner, linear
    !*FD iteration to calculate the new ice thickness distrib

    use glide_velo
    use glide_setup
    use glide_nonlin !For unstable manifold correction
    use glimmer_paramets, only: thk0, thk_scale, GLC_DEBUG
    use glide_grid_operators, only: glide_geometry_derivs

    !EIB! use glide_deriv, only : df_field_2d_staggered 
    implicit none
    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A

    ! local variables
    integer, parameter :: pmax=50                       !*FD maximum Picard iterations

    real(dp), parameter :: tol=1.0d-6
    real(dp) :: residual
    integer p
    logical first_p

#ifdef USE_UNSTABLE_MANIFOLD
    ! local variables used by unstable manifold correction
    real(dp), dimension(model%general%ewn*model%general%nsn) :: umc_new_vec   
    real(dp), dimension(model%general%ewn*model%general%nsn) :: umc_old_vec 
    real(dp), dimension(model%general%ewn*model%general%nsn) :: umc_correction_vec
    logical :: umc_continue_iteration
    integer :: linearize_start

    umc_correction_vec = 0
    umc_new_vec = 0
    umc_old_vec = 0
#endif

    if (model%geometry%empty) then

       model%geometry%thck = dmax1(0.0d0, model%geometry%thck + model%climate%acab * model%numerics%dt)
       if (GLC_DEBUG) then
          print *, "* thck empty - net accumulation added", model%numerics%time
       end if
    else

       !Note: glide_geometry_derivs is called at the beginning of glide_tstep_p1,
       !      and the geometry has not changed, so stagthck and the geometry
       !      derivatives are still up to date.  A call might be needed here
       !      if glide_tstep_p2 were called out of order.
       !      This subroutine must be called during each Picard iteration below.

!!          call glide_geometry_derivs(model)

       ! calculate basal velos
       if (newtemps) then

          call slipvelo(model,                         &
                        1,                             &
                        model%velocity% btrc,          &
                        model%velocity% ubas,          &
                        model%velocity% vbas)

          ! calculate Glen's A if necessary
          call velo_integrate_flwa(model%velowk,             &
                                   model%geomderv%stagthck,  &
                                   model%temper%flwa)

       end if

       first_p = .true.
       model%thckwk%oldthck = model%geometry%thck

       ! do Picard iteration

       model%thckwk%oldthck2 = model%geometry%thck

       do p = 1, pmax

          ! update stagthck, dusrfdew/dns, dthckdew/dns

          call glide_geometry_derivs(model)   

          ! flag = 2: compute basal contribution to diffusivity
          call slipvelo(model,                         &
                        2,                             &
                        model%velocity% btrc,          &
                        model%velocity% ubas,          &
                        model%velocity% vbas)

          ! calculate diffusivity
          call velo_calc_diffu(model%velowk,            model%geomderv%stagthck,  &
                               model%geomderv%dusrfdew, model%geomderv%dusrfdns,  &
                               model%velocity%diffu)

          ! get new thickness
          call thck_evolve(model, model%velocity%diffu, model%velocity%diffu, &
                           first_p, model%geometry%thck, model%geometry%thck)

          first_p = .false.

!TODO - Is this option ever used?  If so, then replace the ifdef with a logical option?
#ifdef USE_UNSTABLE_MANIFOLD
          linearize_start = 1
          call linearize_2d(umc_new_vec, linearize_start, model%geometry%thck)
          linearize_start = 1
          call linearize_2d(umc_old_vec, linearize_start, model%thckwk%oldthck2)
          umc_continue_iteration = unstable_manifold_correction(umc_new_vec, umc_old_vec, &
                                                                umc_correction_vec, size(umc_correction_vec),&
                                                                tol)
          !Only the old thickness might change as a result of this call
           linearize_start = 1
           call delinearize_2d(umc_old_vec, linearize_start, model%thckwk%oldthck2)
          
          if (umc_continue_iteration) then
            exit
          end if
#else
!SCALING - Multiply thickness residual by thk0/thk_scale so we get the same result in these two cases:
!           (1) Old Glimmer with scaling:         thk0 = thk_scale = 2000 m, and thck is non-dimensional
!           (2) New Glimmer-CISM without scaling: thk0 = 1, thk_scale = 2000 m, and thck is in true meters.

!!!          residual = maxval(abs(model%geometry%thck-model%thckwk%oldthck2))
          residual = maxval( abs(model%geometry%thck-model%thckwk%oldthck2) * (thk0/thk_scale) )

          if (residual <= tol) then
             exit
          end if

          model%thckwk%oldthck2 = model%geometry%thck
#endif

       end do   ! Picard iteration

       if (GLC_DEBUG) then
          picard_max = max(picard_max,p)
          if (model%numerics%tinc > mod(model%numerics%time,picard_interval)) then
             write(picard_unit,*) model%numerics%time, p
             picard_max = 0
          end if
       end if

       ! Note: the values for stagthck, slopes, diffu, and ubas (from option 2 call to slipvelo)
       ! will be outdated from the previous Picard iteration.  
       ! To ensure exact restarts are possible, calculate these one more time so that
       ! they can be reconstructed with the restart values of thk and flwa
       ! This will change answers very slightly (to within the Picard convergence tolerance)
       ! relative to older versions of the code.   --MJH 1/9/13

!WHL - oldglide does not update the diffusivity here
!      By skipping this update, I get the same velocities as oldglide on the 
!       first timestep of the dome test case (within double-precision roundoff).  
!      Including this update, velocites agree only to ~4 sig digits.
!       

     if (.not. oldglide) then  ! update the diffusivity before computing the final velocity

       call glide_geometry_derivs(model)   ! stagvarb, geomders as in old Glide code

       ! flag = 2: basal contribution to diffusivity
       call slipvelo(model,                         &
                     2,                             &
                     model%velocity%btrc,           &
                     model%velocity%ubas,           &
                     model%velocity%vbas)

       ! calculate diffusivity
       call velo_calc_diffu(model%velowk,            model%geomderv%stagthck,  &
                            model%geomderv%dusrfdew, model%geomderv%dusrfdns,  &
                            model%velocity%diffu)

     endif   ! oldglide = F

       ! calculate horizontal velocity field

       ! flag = 3: Calculate the basal velocity from the diffusivities
       call slipvelo(model,                         &
                     3,                             &
                     model%velocity%btrc,           &
                     model%velocity%ubas,           &
                     model%velocity%vbas)

       call velo_calc_velo(model%velowk,            model%geomderv%stagthck,   &
                           model%geomderv%dusrfdew, model%geomderv%dusrfdns,   &
                           model%temper%flwa,       model%velocity%diffu,      &
                           model%velocity%ubas,     model%velocity%vbas,       &
                           model%velocity%uvel,     model%velocity%vvel,       &
                           model%velocity%uflx,     model%velocity%vflx,&
                           model%velocity%velnorm)

    end if    ! model%geometry%empty

  end subroutine thck_nonlin_evolve

!---------------------------------------------------------------------------------

  subroutine thck_evolve(model, diffu_x, diffu_y, calc_rhs, old_thck, new_thck)

    !*FD set up sparse matrix and solve matrix equation to find new ice thickness distribution
    !*FD this routine does not override the old thickness distribution

    use glimmer_log
    use glimmer_paramets, only: vel0, thk0, GLC_DEBUG

    implicit none

    ! subroutine arguments -------------------------------------------------------------

    type(glide_global_type) :: model
    logical,intent(in) :: calc_rhs                      !*FD set to true when rhs should be calculated 
                                                        !*FD i.e. when doing lin solution or first picard iteration
    real(dp), intent(in), dimension(:,:) :: diffu_x
    real(dp), intent(in), dimension(:,:) :: diffu_y
    real(dp), intent(in), dimension(:,:) :: old_thck    !*FD contains ice thicknesses from previous time step
    real(dp), intent(inout), dimension(:,:) :: new_thck !*FD on entry contains first guess for new ice thicknesses
                                                        !*FD on exit contains ice thicknesses of new time step

    ! local variables ------------------------------------------------------------------

    real(dp), dimension(5) :: sumd 
    real(dp) :: err
    integer :: linit
    integer :: ew,ns

    real(dp) :: alpha_dt_ew, alpha_dt_ns  ! factors used repeatedly in matrix elements

    alpha_dt_ew = model%numerics%alpha * model%numerics%dt / (2.0d0 * model%numerics%dew * model%numerics%dew)
    alpha_dt_ns = model%numerics%alpha * model%numerics%dt / (2.0d0 * model%numerics%dns * model%numerics%dns)

    ! Zero the arrays holding the sparse matrix
    call sparse_clear(model%solver_data%matrix)

    ! Set the order of the matrix
    model%solver_data%matrix%order = model%geometry%totpts

    !EIB! old way
    ! the number of grid points
    !model%solver_data%pcgsize(1) = model%geometry%totpts
    ! Zero the arrays holding the sparse matrix
    !model%solver_data%pcgval = 0.0
    !model%solver_data%pcgcol = 0 
    !model%solver_data%pcgrow = 0
    !model%solver_data%ct = 1

    ! Boundary Conditions ---------------------------------------------------------------

    ! BCs are for scalar points in outer layer of cells

    ! north and south BC

    do ew = 1,model%general%ewn

       ns=1
       if (model%geometry%thck_index(ew,ns) /= 0) then
          call sparse_insert_val(model%solver_data%matrix, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,ns), 1d0)
          !EIB! old way
          !call putpcgc(model%solver_data,1.0d0, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,ns))
          if (calc_rhs) then
             model%solver_data%rhsd(model%geometry%thck_index(ew,ns)) = old_thck(ew,ns) 
          end if
          model%solver_data%answ(model%geometry%thck_index(ew,ns)) = new_thck(ew,ns)
       end if

       ns=model%general%nsn
       if (model%geometry%thck_index(ew,ns) /= 0) then
          call sparse_insert_val(model%solver_data%matrix, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,ns), 1d0)
          !EIB! old way
          !call putpcgc(model%solver_data,1.0d0, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,ns))
          if (calc_rhs) then
             model%solver_data%rhsd(model%geometry%thck_index(ew,ns)) = old_thck(ew,ns) 
          end if
          model%solver_data%answ(model%geometry%thck_index(ew,ns)) = new_thck(ew,ns)
       end if

    end do

    ! east and west BC

    if (model%options%periodic_ew) then

       do ns=2,model%general%nsn-1
          ew = 1
          if (model%geometry%thck_index(ew,ns) /= 0) then
             call findsums(model%general%ewn-2,model%general%ewn-1,ns-1,ns)
             call generate_row(model%general%ewn-2,ew,ew+1,ns-1,ns,ns+1)
          end if

          ew=model%general%ewn
          if (model%geometry%thck_index(ew,ns) /= 0) then
             call findsums(1,2,ns-1,ns)
             call generate_row(ew-1,ew,3,ns-1,ns,ns+1)
          end if
       end do

    else

       do ns=2,model%general%nsn-1

          ew=1
          if (model%geometry%thck_index(ew,ns) /= 0) then
             call sparse_insert_val(model%solver_data%matrix, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,ns), 1d0)
             !EIB! old way
             !call putpcgc(model%solver_data,1.0d0, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,ns))
             if (calc_rhs) then
                model%solver_data%rhsd(model%geometry%thck_index(ew,ns)) = old_thck(ew,ns) 
             end if
             model%solver_data%answ(model%geometry%thck_index(ew,ns)) = new_thck(ew,ns)
          end if

          ew=model%general%ewn
          if (model%geometry%thck_index(ew,ns) /= 0) then
             call sparse_insert_val(model%solver_data%matrix, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,ns), 1d0)
             !EIB! old way
             !call putpcgc(model%solver_data,1.0d0, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,ns))
             if (calc_rhs) then
                model%solver_data%rhsd(model%geometry%thck_index(ew,ns)) = old_thck(ew,ns) 
             end if
             model%solver_data%answ(model%geometry%thck_index(ew,ns)) = new_thck(ew,ns)
          end if

       end do

    end if   ! periodic_ew

    ! ice interior -------------------------------------------------------------------------

    do ns = 2,model%general%nsn-1
       do ew = 2,model%general%ewn-1

          if (model%geometry%thck_index(ew,ns) /= 0) then
                
             call findsums(ew-1, ew, ns-1, ns)
             call generate_row(ew-1, ew, ew+1, ns-1, ns, ns+1)

          end if
       end do
    end do

    !TODO - EIB - not needed?
    ! Calculate the total number of points
    !model%solver_data%pcgsize(2) = model%solver_data%ct - 1 

    ! Solve the system using SLAP
    !EIB! call slapsolv(model,linit,err)   

    call sparse_easy_solve(model%solver_data%matrix,                       &
                           model%solver_data%rhsd, model%solver_data%answ, &
                           err, linit)

    ! Rejig the solution onto a 2D array

    do ns = 1,model%general%nsn
       do ew = 1,model%general%ewn 
          if (model%geometry%thck_index(ew,ns) /= 0) then
             new_thck(ew,ns) = model%solver_data%answ(model%geometry%thck_index(ew,ns))
          end if
       end do
    end do

    new_thck = max(0.0d0, new_thck)

    if (GLC_DEBUG) then
       print *, "* thck ", model%numerics%time, linit, model%geometry%totpts, &
            real(thk0 * new_thck(model%general%ewn/2+1,model%general%nsn/2+1)), &
            real(vel0 * maxval(abs(model%velocity%ubas))), real(vel0*maxval(abs(model%velocity%vbas))) 
    end if

    !TODO Why are lsrf and usrf calculated here?  This is confusing because model%geometry%thck has only been updated 
    ! because new_thck points to it, but that was only the case because of the way this subroutine is called, and would
    ! not generally be true.  This calculation should be made with new_thck, if it's going to be made here at all!

    ! calculate upper and lower surface

    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

  contains

    subroutine generate_row(ewm, ew, ewp, &
                            nsm, ns, nsp)

      ! calculate row of sparse matrix equation

      implicit none

      integer, intent(in) :: ewm,ew,ewp  ! ew index to left, central, right node
      integer, intent(in) :: nsm,ns,nsp  ! ns index to lower, central, upper node

      !fill matrix using the new API
      call sparse_insert_val(model%solver_data%matrix, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ewm,ns), sumd(1)) ! point (ew-1,ns)
      call sparse_insert_val(model%solver_data%matrix, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ewp,ns), sumd(2)) ! point (ew+1,ns)
      call sparse_insert_val(model%solver_data%matrix, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,nsm), sumd(3)) ! point (ew,ns-1)
      call sparse_insert_val(model%solver_data%matrix, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,nsp), sumd(4)) ! point (ew,ns+1)
      call sparse_insert_val(model%solver_data%matrix, model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,ns),  1d0 + sumd(5))! point (ew,ns)

    !EIB! old way
      ! fill sparse matrix
    !  call putpcgc(model%solver_data,sumd(1), model%geometry%thck_index(ewm,ns), model%geometry%thck_index(ew,ns))       ! point (ew-1,ns)
    !  call putpcgc(model%solver_data,sumd(2), model%geometry%thck_index(ewp,ns), model%geometry%thck_index(ew,ns))       ! point (ew+1,ns)
    !  call putpcgc(model%solver_data,sumd(3), model%geometry%thck_index(ew,nsm), model%geometry%thck_index(ew,ns))       ! point (ew,ns-1)
    !  call putpcgc(model%solver_data,sumd(4), model%geometry%thck_index(ew,nsp), model%geometry%thck_index(ew,ns))       ! point (ew,ns+1)
    !  call putpcgc(model%solver_data,1.0d0 + sumd(5), model%geometry%thck_index(ew,ns), model%geometry%thck_index(ew,ns))! point (ew,ns)

      ! calculate RHS
      if (calc_rhs) then

         model%solver_data%rhsd(model%geometry%thck_index(ew,ns)) =                   &
              old_thck(ew,ns) * (1.0d0 - ((1.0d0-model%numerics%alpha)/model%numerics%alpha) * sumd(5))          &
            - ((1.0d0 - model%numerics%alpha) / model%numerics%alpha) *               &
                                              (old_thck(ewm,ns) * sumd(1)             &
                                             + old_thck(ewp,ns) * sumd(2)             &
                                             + old_thck(ew,nsm) * sumd(3)             &
                                             + old_thck(ew,nsp) * sumd(4))            &
            - (1.d0 / model%numerics%alpha) * (model%geometry%lsrf(ew,ns)  * sumd(5)  &
                                             + model%geometry%lsrf(ewm,ns) * sumd(1)  &
                                             + model%geometry%lsrf(ewp,ns) * sumd(2)  &
                                             + model%geometry%lsrf(ew,nsm) * sumd(3)  &
                                             + model%geometry%lsrf(ew,nsp) * sumd(4)) &
            + model%climate%acab(ew,ns) * model%numerics%dt

         if (model%options%basal_mbal==1) then   ! basal melt rate included in continuity equation
             model%solver_data%rhsd(model%geometry%thck_index(ew,ns)) =                     &
                   model%solver_data%rhsd(model%geometry%thck_index(ew,ns))                 &
                 - model%temper%bmlt(ew,ns) * model%numerics%dt  ! basal melt is positive for mass loss
         end if

      end if   ! calc_rhs

      model%solver_data%answ(model%geometry%thck_index(ew,ns)) = new_thck(ew,ns)      

    end subroutine generate_row

!---------------------------------------------------------------

    subroutine findsums(ewm, ew, nsm, ns)

      ! calculate diffusivities

      implicit none
      integer, intent(in) :: ewm,ew  ! ew index to left, right
      integer, intent(in) :: nsm,ns  ! ns index to lower, upper

      ! calculate sparse matrix elements
      sumd(1) = alpha_dt_ew * (&
               (diffu_x(ewm,nsm) + diffu_x(ewm,ns)) + &
               (model%velocity%ubas (ewm,nsm) + model%velocity%ubas (ewm,ns)))
      sumd(2) = alpha_dt_ew * (&
               (diffu_x(ew,nsm) + diffu_x(ew,ns)) + &
               (model%velocity%ubas (ew,nsm) + model%velocity%ubas (ew,ns)))
      sumd(3) = alpha_dt_ns * (&
               (diffu_y(ewm,nsm) + diffu_y(ew,nsm)) + &
               (model%velocity%ubas (ewm,nsm) + model%velocity%ubas (ew,nsm)))
      sumd(4) = alpha_dt_ns * (&
               (diffu_y(ewm,ns) + diffu_y(ew,ns)) + &
               (model%velocity%ubas (ewm,ns) + model%velocity%ubas (ew,ns)))
      sumd(5) = - (sumd(1) + sumd(2) + sumd(3) + sumd(4))

      !EIB! old way
      !sumd(1) = alpha_dt_ew * (&
      !     (model%velocity%diffu(ewm,nsm) + model%velocity%diffu(ewm,ns)) + &
      !     (model%velocity%ubas (ewm,nsm) + model%velocity%ubas (ewm,ns)))
      !sumd(2) = alpha_dt_ew * (&
      !     (model%velocity%diffu(ew,nsm) + model%velocity%diffu(ew,ns)) + &
      !     (model%velocity%ubas (ew,nsm) + model%velocity%ubas (ew,ns)))
      !sumd(3) = alpha_dt_ns * (&
      !     (model%velocity%diffu(ewm,nsm) + model%velocity%diffu(ew,nsm)) + &
      !     (model%velocity%ubas (ewm,nsm) + model%velocity%ubas (ew,nsm)))
      !sumd(4) = alpha_dt_ns * (&
      !     (model%velocity%diffu(ewm,ns) + model%velocity%diffu(ew,ns)) + &
      !     (model%velocity%ubas (ewm,ns) + model%velocity%ubas (ew,ns)))
      !sumd(5) = - (sumd(1) + sumd(2) + sumd(3) + sumd(4))

    end subroutine findsums

  end subroutine thck_evolve

!---------------------------------------------------------------

!WHL - This subroutine used to be called glide_maskthck and located in its own module, 
!       but I put it in glide_thck.F90 since it is used only for the glide thickness calculation. 

  subroutine glide_thck_index(thck,             acab,  &
                              thck_index,       totpts,  &
                              include_adjacent, empty)
    
    ! Compute an integer mask for the glide thickness calculation.
    ! The mask generally includes ice-covered cells (thck > 0), cells adjacent to 
    !  ice-covered cells, and cells with a positive mass balance (acab > 0).

    !-------------------------------------------------------------------------
    ! Subroutine arguments
    !-------------------------------------------------------------------------

    real(dp),dimension(:,:),intent(in)  :: thck       !*FD Ice thickness
    real(dp),dimension(:,:),intent(in)  :: acab       !*FD Mass balance
    integer, dimension(:,:),intent(out) :: thck_index !*FD integer index (1, 2, 3, ..., totpts)
    integer,                intent(out) :: totpts     !*FD Total number of points in mask
    logical, intent(in)       :: include_adjacent     ! If true, points with no ice but that are adjacent
                                                      ! to points with ice are included in the mask
    logical,                intent(out) :: empty      !*FD true if no points in mask

    !-------------------------------------------------------------------------
    ! Internal variables
    !-------------------------------------------------------------------------

    logical,dimension(size(thck,2))   :: full
    integer :: covtot
    integer :: ew,ns,ewn,nsn

!!    integer,dimension(size(thck,2),2) :: band   ! no longer used
!!    integer, dimension(4) :: dom  ! used to be an output argument, but no longer used

    !-------------------------------------------------------------------------

    ewn=size(thck,1) ; nsn=size(thck,2)

    thck_index = 0
    covtot  = 0 

    !-------------------------------------------------------------------------

    do ns = 1,nsn

      full(ns) = .false.

      do ew = 1,ewn

        if ( thckcrit(thck(max(1,ew-1):min(ewn,ew+1),max(1,ns-1):min(nsn,ns+1)), acab(ew,ns)) ) then

          covtot = covtot + 1
          thck_index(ew,ns) = covtot 

          if ( .not. full(ns) ) then
!!            band(ns,1) = ew
            full(ns)   = .true.
          else
!!            band(ns,2) = ew
          end if
               
        end if
      end do
    end do
  
    totpts = covtot

!!    dom(1:2) = (/ewn,1/)
    empty = .true.

    do ns = 1,nsn
           
      if (full(ns)) then

        if (empty) then
          empty  = .false.
!!          dom(3) = ns
        end if

!!        dom(4) = ns
!!        dom(1) = min0(dom(1),band(ns,1))
!!        dom(2) = max0(dom(2),band(ns,2))
      end if

    end do

  contains

    logical function thckcrit(ca,cb)

      implicit none

      real(dp),dimension(:,:),intent(in) :: ca 
      real(dp),               intent(in) :: cb


!TODO - Is there any case in which we would not want to include adjacent cells
!       in the mask for the thickness calculation?

      if (.not. include_adjacent) then

        ! Include only points with ice in the mask
        ! ca(2,2) corresponds to the current (ew,ns)
 
        if ( ca(2,2) > 0.d0 .or. cb > 0.d0) then
          thckcrit = .true.
        else
          thckcrit = .false.
        end if

      else

        ! If the thickness in the region under consideration
        ! or the mass balance is positive, thckcrit is .true.
        ! This means that the mask includes points that have no
        ! ice but are adjacent to points that do have ice

        if ( any((ca(:,:) > 0.d0)) .or. cb > 0.d0 ) then
          thckcrit = .true.
        else
          thckcrit = .false.
        end if

      end if

    end function thckcrit

  end subroutine glide_thck_index

  !-----------------------------------------------------------------------------
  ! ADI routines
  !-----------------------------------------------------------------------------

  subroutine stagleapthck(model,newtemps)
    
    !*FD this subroutine solves the ice sheet thickness equation using the ADI scheme
    !*FD diffusivities are updated for each half time step

    !TODO The ADI scheme has not been checked for consistency with the new time-stepping convention.  

    use glide_velo
    use glimmer_utils, only: tridiag
    use glimmer_paramets, only: GLC_DEBUG
    use glide_grid_operators, only: glide_geometry_derivs
    implicit none

    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A

    ! local variables
    integer ew,ns, n

    if (model%geometry%empty) then

       model%geometry%thck = dmax1(0.0d0, model%geometry%thck + model%climate%acab * model%numerics%dt)
       if (GLC_DEBUG) then
          print *, "* thck empty - net accumulation added", model%numerics%time
       end if

    else

       !Note: glide_geometry_derivs is called at the beginning of glide_tstep_p1,
       !      and the geometry has not changed, so stagthck and the geometry
       !      derivatives are still up to date.  A call might be needed here
       !      if glide_tstep_p2 were called out of order.

!!       call glide_geometry_derivs(model)

       ! calculate basal velos

       if (newtemps) then
          call slipvelo(model,                &
               1,                             &
               model%velocity% btrc,          &
               model%velocity% ubas,          &
               model%velocity% vbas)

          ! calculate Glen's A if necessary
          call velo_integrate_flwa(model%velowk,model%geomderv%stagthck,model%temper%flwa)
       end if

       call slipvelo(model,                &
            2,                             &
            model%velocity% btrc,          &
            model%velocity% ubas,          &
            model%velocity% vbas)

       ! calculate diffusivity

       call velo_calc_diffu(model%velowk,            model%geomderv%stagthck,  &
                            model%geomderv%dusrfdew, model%geomderv%dusrfdns,  &
                            model%velocity%diffu)

       model%velocity%total_diffu(:,:) = model%velocity%diffu(:,:) + model%velocity%ubas(:,:)

       ! first ADI step, solve thickness equation along rows j

       n = model%general%ewn
       do ns=2,model%general%nsn-1

          call adi_tri ( model%thckwk%alpha,                 &
                         model%thckwk%beta,                  &
                         model%thckwk%gamma,                 &
                         model%thckwk%delta,                 &
                         model%geometry%thck(:,ns),          &
                         model%geometry%lsrf(:,ns),          &
                         model%climate%acab(:,ns),           &
                         model%velocity%vflx(:,ns),          &
                         model%velocity%vflx(:,ns-1),        &
                         model%velocity%total_diffu(:,ns),   &
                         model%velocity%total_diffu(:,ns-1), &
                         model%numerics%dt,                  &
                         model%numerics%dew,                 &
                         model%numerics%dns )
                    !EIB! gc2 acab input, not sure why the difference
                    !model%climate%acab(:,ns)-real(model%options%basal_mbal)*real(model%temper%bmlt(:,ns),sp),           &

          call tridiag(model%thckwk%alpha(1:n),    &
                       model%thckwk%beta(1:n),     &
                       model%thckwk%gamma(1:n),    &
                       model%thckwk%oldthck(:,ns), &
                       model%thckwk%delta(1:n))
       end do

       model%thckwk%oldthck(:,:) = max(model%thckwk%oldthck(:,:), 0.d0)

       ! second ADI step, solve thickness equation along columns i
       n = model%general%nsn
       do ew=2,model%general%ewn-1
          call adi_tri ( model%thckwk%alpha,                 &
                         model%thckwk%beta,                  &
                         model%thckwk%gamma,                 &
                         model%thckwk%delta,                 &
                         model%thckwk%oldthck(ew,:),         &
                         model%geometry%lsrf(ew, :),         &
                         model%climate%acab(ew, :),          &
                         model%velocity%uflx(ew,:),          &
                         model%velocity%uflx(ew-1,:),        &
                         model%velocity%total_diffu(ew,:),   &
                         model%velocity%total_diffu(ew-1,:), &
                         model%numerics%dt,                  &
                         model%numerics%dns,                 &
                         model%numerics%dew )
                      !EIB! again, input difference
                      !model%climate%acab(ew, :)-real(model%options%basal_mbal)*real(model%temper%bmlt(ew, :),sp),          &
                      
          call tridiag(model%thckwk%alpha(1:n),    &
                       model%thckwk%beta(1:n),     &
                       model%thckwk%gamma(1:n),    &
                       model%geometry%thck(ew, :), &
                       model%thckwk%delta(1:n))
       end do

       model%geometry%thck(:,:) = max(model%geometry%thck(:,:), 0.d0)

       ! Apply boundary conditions
       model%geometry%thck(1,:) = 0.d0
       model%geometry%thck(model%general%ewn,:) = 0.d0
       model%geometry%thck(:,1) = 0.d0
       model%geometry%thck(:,model%general%nsn) = 0.d0

       ! calculate horizontal velocity field

       call slipvelo(model,                &
            3,                             &
            model%velocity%btrc,           &
            model%velocity%ubas,           &
            model%velocity%vbas)

       call velo_calc_velo(model%velowk,            model%geomderv%stagthck,  &
                           model%geomderv%dusrfdew, model%geomderv%dusrfdns,  &
                           model%temper%flwa,       model%velocity%diffu,     &
                           model%velocity%ubas,     model%velocity%vbas,      &
                           model%velocity%uvel,     model%velocity%vvel,      &
                           model%velocity%uflx,     model%velocity%vflx,      &
                           model%velocity%velnorm)

    end if   ! empty

    !------------------------------------------------------------
    ! calculate upper and lower surface
    !------------------------------------------------------------

    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

  end subroutine stagleapthck

!---------------------------------------------------------------------------------

  subroutine adi_tri(a,b,c,d,thk,tpg,mb,flx_p,flx_m,dif_p,dif_m,dt,ds1, ds2)

    !*FD construct tri-diagonal matrix system for a column/row

    implicit none
    
    real(dp), dimension(:), intent(out) :: a !*FD alpha (subdiagonal)
    real(dp), dimension(:), intent(out) :: b !*FD alpha (diagonal)
    real(dp), dimension(:), intent(out) :: c !*FD alpha (superdiagonal)
    real(dp), dimension(:), intent(out) :: d !*FD right-hand side
    
    real(dp), dimension(:), intent(in) :: thk   !*FD ice thickness
    real(dp), dimension(:), intent(in) :: tpg   !*FD lower surface of ice
    real(dp), dimension(:), intent(in) :: mb    !*FD mass balance
    real(dp), dimension(:), intent(in) :: flx_p !*FD flux +1/2
    real(dp), dimension(:), intent(in) :: flx_m !*FD flux -1/2
    real(dp), dimension(:), intent(in) :: dif_p !*FD diffusivity +1/2
    real(dp), dimension(:), intent(in) :: dif_m !*FD diffusivity -1/2
    
    real(dp), intent(in) :: dt !*FD time step
    real(dp), intent(in) :: ds1, ds2 !*FD spatial steps inline and transversal

    ! local variables
    real(dp) :: f1, f2, f3
    integer :: i,n
    
    n = size(thk)

    f1 = dt/(4*ds1*ds1)
    f2 = dt/(4*ds2)
    f3 = dt/2.

    a(:) = 0.
    b(:) = 0.
    c(:) = 0.
    d(:) = 0.

    a(1) = 0.
    do i=2,n
       a(i) = f1*(dif_m(i-1)+dif_p(i-1))
    end do
    do i=1,n-1
       c(i) = f1*(dif_m(i)+dif_p(i))
    end do
    c(n) = 0.
    b(:) = -(a(:)+c(:))

    ! calculate RHS
    do i=2,n-1
       d(i) = thk(i) - &
            f2 * (flx_p(i-1) + flx_p(i) - flx_m(i-1) - flx_m(i)) + &
            f3 * mb(i) - &
            a(i)*tpg(i-1) - b(i)*tpg(i) - c(i)*tpg(i+1)
    end do

    b(:) = 1.+b(:)

  end subroutine adi_tri

!-------------------------------------------------------------------------

  subroutine glide_calclsrf(thck,topg,eus,lsrf)

    ! Calculates the elevation of the lower surface of the ice, 
    ! by considering whether it is floating or not.
    !
    ! NOTE: This subroutine computes over all grid cells, not just locally owned.
    !       Halos should be updated before it is called.

    use glimmer_physcon, only : rhoi, rhoo

    implicit none

    real(dp), intent(in),  dimension(:,:) :: thck !*FD Ice thickness
    real(dp), intent(in),  dimension(:,:) :: topg !*FD Bedrock topography elevation
    real(dp), intent(in)                  :: eus  !*FD global sea level
    real(dp), intent(out), dimension(:,:) :: lsrf !*FD Lower ice surface elevation

    real(dp), parameter :: con = - rhoi / rhoo

    where (topg - eus < con * thck)
       lsrf = con * thck
    elsewhere
       lsrf = topg
    end where

  end subroutine glide_calclsrf

!---------------------------------------------------------------------------------

!TODO - This subroutine is not used.  Remove it?

  subroutine filterthck(thck,ewn,nsn)

    implicit none

    real(dp), dimension(:,:), intent(inout) :: thck
    real(dp), dimension(:,:), allocatable :: smth
    integer :: ewn,nsn

    real(dp), parameter :: f = 0.1d0 / 16.0d0
    integer :: count
    integer :: ns,ew

    allocate(smth(ewn,nsn))
    count = 1

    do ns = 3,nsn-2
      do ew = 3,ewn-2

        if (all((thck(ew-2:ew+2,ns) > 0.0d0)) .and. all((thck(ew,ns-2:ns+2) > 0.0d0))) then
          smth(ew,ns) =  thck(ew,ns) + f * &
                        (thck(ew-2,ns) - 4.0d0 * thck(ew-1,ns) + 12.0d0 * thck(ew,ns) - &
                         4.0d0 * thck(ew+1,ns) + thck(ew+2,ns) + &
                         thck(ew,ns-2) - 4.0d0 * thck(ew,ns-1) - &
                         4.0d0 * thck(ew,ns+1) + thck(ew,ns+2))
          count = count + 1
        else
          smth(ew,ns) = thck(ew,ns)
        end if

      end do
    end do

    thck(3:ewn-2,3:nsn-2) = smth(3:ewn-2,3:nsn-2)
    print *, count

    deallocate(smth)            

  end subroutine filterthck

!----------------------------------------------------------------------

!TODO - This subroutine is not used.  Remove it?

  subroutine swapbndh(bc,a,b,c,d)

    implicit none

    real(dp), intent(out), dimension(:) :: a, c
    real(dp), intent(in), dimension(:) :: b, d
    integer, intent(in) :: bc

    if (bc == 0) then
      a = b
      c = d
    end if

  end subroutine swapbndh

!---------------------------------------------------------------------------------

end module glide_thck

!---------------------------------------------------------------------------------
