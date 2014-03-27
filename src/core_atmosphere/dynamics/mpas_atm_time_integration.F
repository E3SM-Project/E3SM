! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module atm_time_integration

   use mpas_grid_types
   use mpas_kind_types
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   ! Added only clause to keep xlf90 from getting confused from the overloaded abs intrinsic in mpas_timekeeping
   use mpas_timekeeping, only: MPAS_Time_type, MPAS_TimeInterval_type, &
                               mpas_set_time, mpas_set_timeInterval, mpas_get_time, operator(+), add_t_ti

#ifdef DO_PHYSICS
   use mpas_atmphys_driver_microphysics
   use mpas_atmphys_todynamics
   use mpas_atmphys_utilities
#endif

   contains


   subroutine atm_timestep(domain, dt, timeStamp, itimestep)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp
      integer, intent(in) :: itimestep


      type (block_type), pointer :: block
      type (MPAS_Time_type) :: currTime
      type (MPAS_TimeInterval_type) :: dtInterval
      character (len=StrKIND), pointer :: xtime
      character (len=StrKIND) :: xtime_new
      type (mpas_pool_type), pointer :: state
      character (len=StrKIND), pointer :: config_time_integration


      call mpas_pool_get_config(domain % blocklist % configs, 'config_time_integration', config_time_integration)

      if (trim(config_time_integration) == 'SRK3') then
         call atm_srk3(domain, dt, itimestep)
      else
         write(0,*) 'Unknown time integration option '//trim(config_time_integration)
         write(0,*) 'Currently, only ''SRK3'' is supported.'
         call mpas_dmpar_abort(domain % dminfo)
      end if

      call mpas_set_time(currTime, dateTimeString=timeStamp)
      call mpas_set_timeInterval(dtInterval, dt=dt)
      currTime = currTime + dtInterval
      call mpas_get_time(currTime, dateTimeString=xtime_new)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_array(state, 'xtime', xtime, 2)
         xtime = xtime_new
         block => block % next
      end do

   end subroutine atm_timestep


   subroutine atm_srk3(domain, dt, itimestep)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step using 
   !   time-split RK3 scheme
   !
   ! Nonhydrostatic atmospheric solver
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data and some diagnostics of state.
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds,
   !                  and some diagnostics in diag 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: itimestep


      integer :: iCell, k, iEdge
      type (block_type), pointer :: block

      integer :: rk_step, number_of_sub_steps
      integer :: iScalar

      real (kind=RKIND), dimension(3) :: rk_timestep, rk_sub_timestep
      integer, dimension(3) :: number_sub_steps
      integer :: small_step
      logical, parameter :: debug = .false.
!      logical, parameter :: debug = .true.
      logical, parameter :: debug_mass_conservation = .true.

      real (kind=RKIND) :: scalar_min, scalar_max
      real (kind=RKIND) :: global_scalar_min, global_scalar_max

      integer, pointer :: config_number_of_sub_steps
      logical, pointer :: config_scalar_advection
      logical, pointer :: config_positive_definite
      logical, pointer :: config_monotonic
      real (kind=RKIND), pointer :: config_dt
      character (len=StrKIND), pointer :: config_microp_scheme

      integer, pointer :: num_scalars, index_qv, nCells, nCellsSolve, nEdges, nEdgesSolve, nVertLevels

      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: diag_physics
      type (mpas_pool_type), pointer :: mesh
      type (mpas_pool_type), pointer :: tend
      type (mpas_pool_type), pointer :: tend_physics

      type (field2DReal), pointer :: theta_m_field
      type (field3DReal), pointer :: scalars_field
      type (field2DReal), pointer :: pressure_p_field
      type (field2DReal), pointer :: rtheta_p_field
      type (field2DReal), pointer :: rtheta_pp_field
      type (field2DReal), pointer :: tend_u_field
      type (field2DReal), pointer :: u_field
      type (field2DReal), pointer :: w_field
      type (field2DReal), pointer :: rw_p_field
      type (field2DReal), pointer :: ru_p_field
      type (field2DReal), pointer :: rho_pp_field
      type (field2DReal), pointer :: pv_edge_field
      type (field2DReal), pointer :: rho_edge_field

      real (kind=RKIND), dimension(:,:), pointer :: w
      real (kind=RKIND), dimension(:,:), pointer :: u, uReconstructZonal, uReconstructMeridional, uReconstructX, uReconstructY, uReconstructZ
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, scalars_1, scalars_2

      real (kind=RKIND), dimension(:,:), pointer :: rqvdynten


      !
      ! Retrieve configuration options
      !
      call mpas_pool_get_config(domain % blocklist % configs, 'config_number_of_sub_steps', config_number_of_sub_steps)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_scalar_advection', config_scalar_advection)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_positive_definite', config_positive_definite)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_monotonic', config_monotonic)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', config_dt)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_microp_scheme', config_microp_scheme)

      !
      ! Retrieve field structures
      !
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)

      !
      ! Retrieve fields
      !
      call mpas_pool_get_field(state, 'theta_m', theta_m_field, 1)
      call mpas_pool_get_field(state, 'scalars', scalars_field, 1)
      call mpas_pool_get_field(diag, 'pressure_p', pressure_p_field)
      call mpas_pool_get_field(diag, 'rtheta_p', rtheta_p_field)


      !
      ! Initialize RK weights
      !

      number_of_sub_steps = config_number_of_sub_steps
      rk_timestep(1) = dt/3.
      rk_timestep(2) = dt/2.
      rk_timestep(3) = dt

      rk_sub_timestep(1) = dt/3.
      rk_sub_timestep(2) = dt/real(number_of_sub_steps)
      rk_sub_timestep(3) = dt/real(number_of_sub_steps)

      number_sub_steps(1) = 1
      number_sub_steps(2) = number_of_sub_steps/2
      number_sub_steps(3) = number_of_sub_steps

      if(debug) write(0,*) ' copy step in rk solver '

! theta_m
      call mpas_dmpar_exch_halo_field(theta_m_field)
 
! scalars
      call mpas_dmpar_exch_halo_field(scalars_field)

! pressure_p
      call mpas_dmpar_exch_halo_field(pressure_p_field)

! rtheta_p
      call mpas_dmpar_exch_halo_field(rtheta_p_field)


      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'diag', diag)
         call atm_rk_integration_setup(state, diag)
         block => block % next
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! BEGIN Runge-Kutta loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      do rk_step = 1, 3  ! Runge-Kutta loop

         if(debug) write(0,*) ' rk substep ', rk_step

         block => domain % blocklist
         do while (associated(block))
            ! The coefficients are set for owned cells (cqw) and for all edges of owned cells, 
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_get_subpool(block % structs, 'diag', diag)
            call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
            call atm_compute_moist_coefficients( block % dimensions, state, diag, mesh )    !MGD could do away with dimensions arg
            block => block % next
         end do


         if (debug) write(0,*) ' compute_dyn_tend '
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_get_subpool(block % structs, 'diag', diag)
            call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block % structs, 'tend', tend)

            call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

            call atm_compute_dyn_tend( tend, state, diag, mesh, block % configs, nVertLevels, rk_step, dt )

            block => block % next
         end do
         if (debug) write(0,*) ' finished compute_dyn_tend '

#ifdef DO_PHYSICS
         if (debug) write(0,*) ' add physics tendencies '
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_get_subpool(block % structs, 'diag', diag)
            call mpas_pool_get_subpool(block % structs, 'tend', tend)
            call mpas_pool_get_subpool(block % structs, 'tend_physics', tend_physics)
            call physics_addtend( block, &
                         mesh, &
                         state, &     
                         diag, &
                         tend, &
                         tend_physics, &
                         block % configs, &
                         rk_step )
            block => block % next
         end do
         if (debug) write(0,*) ' finished add physics tendencies '
#endif

         !***********************************
         !  need tendencies at all edges of owned cells -
         !  we are solving for all edges of owned cells to minimize communications
         !  during the acoustic substeps
         !***********************************

! tend_u
         call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tend)
         call mpas_pool_get_field(tend, 'u', tend_u_field)
         call mpas_dmpar_exch_halo_field(tend_u_field, (/ 1 /))

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_get_subpool(block % structs, 'diag', diag)
            call mpas_pool_get_subpool(block % structs, 'tend', tend)

            call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

            call atm_set_smlstep_pert_variables( tend, diag, mesh, block % configs )
            call atm_compute_vert_imp_coefs( state, mesh, diag, block % configs, nVertLevels, rk_sub_timestep(rk_step) )

            block => block % next
         end do

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! begin acoustic steps loop
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do small_step = 1, number_sub_steps(rk_step)

            if(debug) write(0,*) ' acoustic step ',small_step
      
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
               call mpas_pool_get_subpool(block % structs, 'state', state)
               call mpas_pool_get_subpool(block % structs, 'diag', diag)
               call mpas_pool_get_subpool(block % structs, 'tend', tend)

               call mpas_pool_get_dimension(mesh, 'nCells', nCells)
               call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

               call atm_advance_acoustic_step( state, diag, tend,  mesh, block % configs, nCells, nVertLevels, rk_sub_timestep(rk_step) )

               block => block % next
            end do

            if(debug) write(0,*) ' acoustic step complete '
  
! rtheta_pp
! This is the only communications needed during the acoustic steps because we solve for u on all edges of owned cells

            call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
            call mpas_pool_get_field(diag, 'rtheta_pp', rtheta_pp_field)
            call mpas_dmpar_exch_halo_field(rtheta_pp_field, (/ 1 /))
 
         end do  ! end of acoustic steps loop

         !CR: SMALLER STENCIL?: call mpas_dmpar_exch_halo_field(block % diag % rw_p, (/ 1 /))
         call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
         call mpas_pool_get_field(diag, 'rw_p', rw_p_field)
         call mpas_dmpar_exch_halo_field(rw_p_field)

         !CR: SMALLER STENCIL?: call mpas_dmpar_exch_halo_field(block % diag % ru_p, (/ 2 /))
         call mpas_pool_get_field(diag, 'ru_p', ru_p_field)
         call mpas_dmpar_exch_halo_field(ru_p_field)

         call mpas_pool_get_field(diag, 'rho_pp', rho_pp_field)
         call mpas_dmpar_exch_halo_field(rho_pp_field)

         ! the second layer of halo cells must be exchanged before calling atm_recover_large_step_variables
         call mpas_pool_get_field(diag, 'rtheta_pp', rtheta_pp_field)
         call mpas_dmpar_exch_halo_field(rtheta_pp_field, (/ 2 /))

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_get_subpool(block % structs, 'diag', diag)
            call mpas_pool_get_subpool(block % structs, 'tend', tend)

            call atm_recover_large_step_variables( state, diag, tend, mesh, block % configs, rk_timestep(rk_step), number_sub_steps(rk_step), rk_step  )

            block => block % next
         end do

! u
         !CR: SMALLER STENCIL?: call mpas_dmpar_exch_halo_field(block % state % time_levs(2) % state % u, (/ 3 /))
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
         call mpas_pool_get_field(state, 'u', u_field, 2)
         call mpas_dmpar_exch_halo_field(u_field)

         ! scalar advection: RK3 scheme of Skamarock and Gassmann (2011). 
         ! PD or monotonicity constraints applied only on the final Runge-Kutta substep.

         if (config_scalar_advection) then

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'tend', tend)
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_get_subpool(block % structs, 'diag', diag)
            call mpas_pool_get_subpool(block % structs, 'mesh', mesh)

            call mpas_pool_get_dimension(mesh, 'nCells', nCells)
            call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
            call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
            call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)

            !
            ! Note: The advance_scalars_mono routine can be used without limiting, and thus, encompasses 
            !       the functionality of the advance_scalars routine; however, it is noticeably slower, 
            !       so we use the advance_scalars routine for the first two RK substeps.
            !
            if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
               call atm_advance_scalars( tend, state, diag, mesh, block % configs, num_scalars, nVertLevels, rk_timestep(rk_step) )
            else
               block % domain = domain 
               call atm_advance_scalars_mono( block, tend, state, diag, mesh, block % configs, nCells, nEdges, nVertLevels, rk_timestep(rk_step) )
            end if
            block => block % next
         end do

         else
 
            write(0,*) ' no scalar advection '

         end if

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_get_subpool(block % structs, 'diag', diag)
            call mpas_pool_get_subpool(block % structs, 'mesh', mesh)

            call atm_compute_solve_diagnostics( dt, state, 2, diag, mesh, block % configs )

            block => block % next
         end do

         if(debug) write(0,*) ' diagnostics complete '

! w
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
         call mpas_pool_get_field(state, 'w', w_field, 2)
         call mpas_dmpar_exch_halo_field(w_field)

! pv_edge
         call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
         call mpas_pool_get_field(diag, 'pv_edge', pv_edge_field)
         call mpas_dmpar_exch_halo_field(pv_edge_field)

! rho_edge
         call mpas_pool_get_field(diag, 'rho_edge', rho_edge_field)
         call mpas_dmpar_exch_halo_field(rho_edge_field)

! scalars
         if (rk_step < 3) then
            call mpas_pool_get_field(state, 'scalars', scalars_field, 2)
            call mpas_dmpar_exch_halo_field(scalars_field)
         end if

      end do ! rk_step loop

!...  compute full velocity vectors at cell centers:
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'diag', diag)
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)

         call mpas_pool_get_array(state, 'u', u, 2)
         call mpas_pool_get_array(diag, 'uReconstructX', uReconstructX)
         call mpas_pool_get_array(diag, 'uReconstructY', uReconstructY)
         call mpas_pool_get_array(diag, 'uReconstructZ', uReconstructZ)
         call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
         call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)

         call mpas_reconstruct(mesh, u,                &
                               uReconstructX,          &
                               uReconstructY,          &
                               uReconstructZ,          &
                               uReconstructZonal,      &
                               uReconstructMeridional  &
                              )

         block => block % next
      end do

!... call to parameterizations of cloud microphysics. calculation of the tendency of water vapor to horizontal and
!... vertical advection needed for the Tiedtke parameterization of convection.

#ifdef DO_PHYSICS
      block => domain % blocklist
      do while(associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'diag', diag)
         call mpas_pool_get_subpool(block % structs, 'diag_physics', diag_physics)
         call mpas_pool_get_subpool(block % structs, 'tend_physics', tend_physics)
         call mpas_pool_get_subpool(block % structs, 'tend', tend)
         call mpas_pool_get_array(state, 'scalars', scalars_1, 1)
         call mpas_pool_get_array(state, 'scalars', scalars_2, 2)
         call mpas_pool_get_dimension(state, 'index_qv', index_qv)
         call mpas_pool_get_array(tend_physics, 'rqvdynten', rqvdynten)

         !NOTE: The calculation of the tendency due to horizontal and vertical advection for the water vapor mixing ratio
         !requires that the subroutine atm_advance_scalars_mono was called on the third Runge Kutta step, so that a halo
         !update for the scalars at time_levs(1) is applied. A halo update for the scalars at time_levs(2) is done above. 
         if (config_monotonic) then
            rqvdynten(:,:) = ( scalars_2(index_qv,:,:) - scalars_1(index_qv,:,:) ) / config_dt
         else
            rqvdynten(:,:) = 0._RKIND
         end if

         !simply set to zero negative mixing ratios of different water species (for now):
         where ( scalars_2(:,:,:) < 0.0)  &
            scalars_2(:,:,:) = 0.0

         !call microphysics schemes:
         if (config_microp_scheme .ne. 'off')  &
            call microphysics_driver ( block % configs, mesh, state, 2, diag, diag_physics, tend, itimestep )

         block => block % next
      end do
#endif

      102 format(' global min, max scalar',i4,2(1x,e17.10))
      write(0,*)
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)

         call mpas_pool_get_array(state, 'w', w, 2)
         call mpas_pool_get_array(state, 'u', u, 2)
         call mpas_pool_get_array(state, 'scalars', scalars, 2)
         call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
         call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

         scalar_min = 0.
         scalar_max = 0.
         do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            scalar_min = min(scalar_min, w(k,iCell))
            scalar_max = max(scalar_max, w(k,iCell))
         end do
         end do
         call mpas_dmpar_min_real(domain % dminfo, scalar_min, global_scalar_min)
         call mpas_dmpar_max_real(domain % dminfo, scalar_max, global_scalar_max)
         write(0,*) 'global min, max w ',global_scalar_min, global_scalar_max

         scalar_min = 0.
         scalar_max = 0.
         do iEdge = 1, nEdgesSolve
         do k = 1, nVertLevels
            scalar_min = min(scalar_min, u(k,iEdge))
            scalar_max = max(scalar_max, u(k,iEdge))
         end do
         end do
         call mpas_dmpar_min_real(domain % dminfo, scalar_min, global_scalar_min)
         call mpas_dmpar_max_real(domain % dminfo, scalar_max, global_scalar_max)
         write(0,*) 'global min, max u ',global_scalar_min, global_scalar_max

         do iScalar = 1, num_scalars
            scalar_min = 0.
            scalar_max = 0.
            do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               scalar_min = min(scalar_min, scalars(iScalar,k,iCell))
               scalar_max = max(scalar_max, scalars(iScalar,k,iCell))
            end do
            end do
            call mpas_dmpar_min_real(domain % dminfo, scalar_min, global_scalar_min)
            call mpas_dmpar_max_real(domain % dminfo, scalar_max, global_scalar_max)
            write(0,102) iScalar,global_scalar_min,global_scalar_max
         end do

         block => block % next
      end do

   end subroutine atm_srk3

!---

   subroutine atm_rk_integration_setup( state, diag )

      implicit none

      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag

      real (kind=RKIND), dimension(:,:), pointer :: ru
      real (kind=RKIND), dimension(:,:), pointer :: ru_save
      real (kind=RKIND), dimension(:,:), pointer :: rw
      real (kind=RKIND), dimension(:,:), pointer :: rw_save
      real (kind=RKIND), dimension(:,:), pointer :: rtheta_p
      real (kind=RKIND), dimension(:,:), pointer :: rtheta_p_save
      real (kind=RKIND), dimension(:,:), pointer :: rho_p
      real (kind=RKIND), dimension(:,:), pointer :: rho_p_save

      real (kind=RKIND), dimension(:,:), pointer :: u_1, u_2
      real (kind=RKIND), dimension(:,:), pointer :: w_1, w_2
      real (kind=RKIND), dimension(:,:), pointer :: theta_m_1, theta_m_2
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz_1, rho_zz_2
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars_1, scalars_2

      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'ru_save', ru_save)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'rw_save', rw_save)
      call mpas_pool_get_array(diag, 'rtheta_p', rtheta_p)
      call mpas_pool_get_array(diag, 'rtheta_p_save', rtheta_p_save)
      call mpas_pool_get_array(diag, 'rho_p', rho_p)
      call mpas_pool_get_array(diag, 'rho_p_save', rho_p_save)

      call mpas_pool_get_array(state, 'u', u_1, 1)
      call mpas_pool_get_array(state, 'u', u_2, 2)
      call mpas_pool_get_array(state, 'w', w_1, 1)
      call mpas_pool_get_array(state, 'w', w_2, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m_1, 1)
      call mpas_pool_get_array(state, 'theta_m', theta_m_2, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_1, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_2, 2)
      call mpas_pool_get_array(state, 'scalars', scalars_1, 1)
      call mpas_pool_get_array(state, 'scalars', scalars_2, 2)

      ru_save(:,:) = ru(:,:)
      rw_save(:,:) = rw(:,:)
      rtheta_p_save(:,:) = rtheta_p(:,:)
      rho_p_save(:,:) = rho_p(:,:)

      u_2(:,:) = u_1(:,:)
      w_2(:,:) = w_1(:,:)
      theta_m_2(:,:) = theta_m_1(:,:)
      rho_zz_2(:,:) = rho_zz_1(:,:)
      scalars_2(:,:,:) = scalars_1(:,:,:)

   end subroutine atm_rk_integration_setup

!-----

   subroutine atm_compute_moist_coefficients( dims, state, diag, mesh )

      ! the moist coefficients cqu and cqw serve to transform the inverse dry density (1/rho_d) 
      ! into the inverse full (moist) density (1/rho_m).

      implicit none

      type (mpas_pool_type), intent(in) :: dims
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: mesh


      integer :: iEdge, iCell, k, cell1, cell2, iq
      integer, pointer :: nCells, nEdges, nVertLevels, nCellsSolve
      real (kind=RKIND) :: qtot
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: moist_start, moist_end
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: cqw
      real (kind=RKIND), dimension(:,:), pointer :: cqu

      call mpas_pool_get_dimension(dims, 'nCells', nCells)
      call mpas_pool_get_dimension(dims, 'nEdges', nEdges)
      call mpas_pool_get_dimension(dims, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(state, 'moist_start', moist_start)
      call mpas_pool_get_dimension(state, 'moist_end', moist_end)

      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_array(diag, 'cqw', cqw)
      call mpas_pool_get_array(diag, 'cqu', cqu)


      do iCell = 1, nCellsSolve
         do k = 2, nVertLevels
            qtot = 0.
            do iq = moist_start, moist_end
               qtot = qtot + 0.5 * (scalars(iq, k, iCell) + scalars(iq, k-1, iCell))
            end do
            cqw(k,iCell) = 1./(1.+qtot)
         end do
      end do

      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k = 1, nVertLevels
               qtot = 0.
               do iq = moist_start, moist_end
                  qtot = qtot + 0.5 * ( scalars(iq, k, cell1) + scalars(iq, k, cell2) )
               end do
               cqu(k,iEdge) = 1./( 1. + qtot)
            end do
         end if
      end do

   end subroutine atm_compute_moist_coefficients

!---

   subroutine atm_compute_vert_imp_coefs(state, mesh, diag, configs, nVertLevels, dts)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute coefficients for vertically implicit gravity-wave/acoustic computations
   !
   ! Input: state - current model state
   !        mesh - grid metadata
   !
   ! Output: diag - cofrz, cofwr, cofwz, coftz, cofwt, a, alpha and gamma
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(in)    :: state
      type (mpas_pool_type), intent(in)    :: mesh
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in)    :: configs
      integer, intent(in)                  :: nVertLevels          ! for allocating stack variables
      real (kind=RKIND), intent(in)        :: dts


      integer :: iCell, k, iq

      integer, pointer :: nCells, nCellsSolve
      real (kind=RKIND), dimension(:,:), pointer :: zz, cqw, p, t, rb, rtb, pb, rt
      real (kind=RKIND), dimension(:,:), pointer :: cofwr, cofwz, coftz, cofwt, a_tri, alpha_tri, gamma_tri
      real (kind=RKIND), dimension(:), pointer :: cofrz, rdzw, fzm, fzp, rdzu
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      real (kind=RKIND), dimension( nVertLevels ) :: b_tri,c_tri
      real (kind=RKIND), pointer :: epssm
      real (kind=RKIND) :: dtseps, c2, qtot, rcv

      integer, pointer :: moist_start, moist_end

!  set coefficients

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_config(configs, 'config_epssm', epssm)

      call mpas_pool_get_array(mesh, 'rdzu', rdzu)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zz', zz)

      call mpas_pool_get_array(diag, 'cqw', cqw)
      call mpas_pool_get_array(diag, 'exner', p)
      call mpas_pool_get_array(diag, 'exner_base', pb)
      call mpas_pool_get_array(diag, 'rtheta_p', rt)
      call mpas_pool_get_array(diag, 'rtheta_base', rtb)
      call mpas_pool_get_array(diag, 'rho_base', rb)

      call mpas_pool_get_array(diag, 'alpha_tri', alpha_tri)
      call mpas_pool_get_array(diag, 'gamma_tri', gamma_tri)
      call mpas_pool_get_array(diag, 'a_tri', a_tri)
      call mpas_pool_get_array(diag, 'cofwr', cofwr)
      call mpas_pool_get_array(diag, 'cofwz', cofwz)
      call mpas_pool_get_array(diag, 'coftz', coftz)
      call mpas_pool_get_array(diag, 'cofwt', cofwt)
      call mpas_pool_get_array(diag, 'cofrz', cofrz)

      call mpas_pool_get_array(state, 'theta_m', t, 2)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_dimension(state, 'moist_start', moist_start)
      call mpas_pool_get_dimension(state, 'moist_end', moist_end)


      dtseps = .5*dts*(1.+epssm)
      rcv = rgas/(cp-rgas)
      c2 = cp*rcv

      do k=1,nVertLevels
         cofrz(k) = dtseps*rdzw(k)
      end do

      do iCell = 1, nCellsSolve  !  we only need to do cells we are solving for, not halo cells

         do k=2,nVertLevels
            cofwr(k,iCell) =.5*dtseps*gravity*(fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))
         end do
         coftz(1,iCell) = 0.0
         do k=2,nVertLevels
            cofwz(k,iCell) = dtseps*c2*(fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))  &
                 *rdzu(k)*cqw(k,iCell)*(fzm(k)*p (k,iCell)+fzp(k)*p (k-1,iCell))
            coftz(k,iCell) = dtseps*   (fzm(k)*t (k,iCell)+fzp(k)*t (k-1,iCell))
         end do
         coftz(nVertLevels+1,iCell) = 0.0
         do k=1,nVertLevels

            qtot = 0.
            do iq = moist_start, moist_end
               qtot = qtot + scalars(iq, k, iCell)
            end do

            cofwt(k,iCell) = .5*dtseps*rcv*zz(k,iCell)*gravity*rb(k,iCell)/(1.+qtot)  &
                                *p(k,iCell)/((rtb(k,iCell)+rt(k,iCell))*pb(k,iCell))
         end do

         a_tri(1,iCell) = 0.  ! note, this value is never used
         b_tri(1) = 1.    ! note, this value is never used
         c_tri(1) = 0.    ! note, this value is never used
         gamma_tri(1,iCell) = 0.
         alpha_tri(1,iCell) = 0.  ! note, this value is never used

         do k=2,nVertLevels
            a_tri(k,iCell) = -cofwz(k  ,iCell)* coftz(k-1,iCell)*rdzw(k-1)*zz(k-1,iCell)   &
                         +cofwr(k  ,iCell)* cofrz(k-1  )                       &
                         -cofwt(k-1,iCell)* coftz(k-1,iCell)*rdzw(k-1)
            b_tri(k) = 1.                                                  &
                         +cofwz(k  ,iCell)*(coftz(k  ,iCell)*rdzw(k  )*zz(k  ,iCell)   &
                                      +coftz(k  ,iCell)*rdzw(k-1)*zz(k-1,iCell))   &
                         -coftz(k  ,iCell)*(cofwt(k  ,iCell)*rdzw(k  )             &
                                       -cofwt(k-1,iCell)*rdzw(k-1))            &
                         +cofwr(k,  iCell)*(cofrz(k    )-cofrz(k-1))
            c_tri(k) =   -cofwz(k  ,iCell)* coftz(k+1,iCell)*rdzw(k  )*zz(k  ,iCell)   &
                         -cofwr(k  ,iCell)* cofrz(k    )                       &
                         +cofwt(k  ,iCell)* coftz(k+1,iCell)*rdzw(k  )
         end do
         do k=2,nVertLevels
            alpha_tri(k,iCell) = 1./(b_tri(k)-a_tri(k,iCell)*gamma_tri(k-1,iCell))
            gamma_tri(k,iCell) = c_tri(k)*alpha_tri(k,iCell)
         end do

      end do ! loop over cells

   end subroutine atm_compute_vert_imp_coefs

!------------------------

   subroutine atm_set_smlstep_pert_variables( tend, diag, mesh, configs )

      ! following Klemp et al MWR 2007, we use preturbation variables
      ! in the acoustic-step integration.  This routine computes those 
      ! perturbation variables.  state variables are reconstituted after 
      ! the acousstic steps in subroutine atm_recover_large_step_variables


      implicit none

      type (mpas_pool_type), intent(inout) :: tend
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in) :: configs

      integer :: iCell, iEdge, k, cell1, cell2
      real (kind=RKIND), pointer :: coef_3rd_order
      integer, pointer :: config_theta_adv_order
      integer, pointer :: nCellsSolve, nCells, nVertLevels, nEdges
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp, dvEdge, areaCell
      real (kind=RKIND) :: flux
      real (kind=RKIND), dimension(:,:), pointer :: ruAvg, wwAvg
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
      real (kind=RKIND), dimension(:,:), pointer :: zz
      real (kind=RKIND), dimension(:,:), pointer :: w_tend, u_tend
      real (kind=RKIND), dimension(:,:), pointer :: rho_pp, rho_p_save, rho_p
      real (kind=RKIND), dimension(:,:), pointer :: ru_p, ru, ru_save
      real (kind=RKIND), dimension(:,:), pointer :: rtheta_pp, rtheta_p_save, rtheta_p, rtheta_pp_old
      real (kind=RKIND), dimension(:,:), pointer :: rw_p, rw_save, rw


      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(tend, 'w', w_tend)
      call mpas_pool_get_array(tend, 'u', u_tend)

      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)

      call mpas_pool_get_array(diag, 'rho_pp', rho_pp)
      call mpas_pool_get_array(diag, 'rho_p_save', rho_p_save)
      call mpas_pool_get_array(diag, 'rho_p', rho_p)

      call mpas_pool_get_array(diag, 'ru_p', ru_p)
      call mpas_pool_get_array(diag, 'ru_save', ru_save)
      call mpas_pool_get_array(diag, 'ru', ru)

      call mpas_pool_get_array(diag, 'rtheta_pp', rtheta_pp)
      call mpas_pool_get_array(diag, 'rtheta_p_save', rtheta_p_save)
      call mpas_pool_get_array(diag, 'rtheta_p', rtheta_p)
      call mpas_pool_get_array(diag, 'rtheta_pp_old', rtheta_pp_old)

      call mpas_pool_get_array(diag, 'rw_p', rw_p)
      call mpas_pool_get_array(diag, 'rw_save', rw_save)
      call mpas_pool_get_array(diag, 'rw', rw)

      if (config_theta_adv_order /= 3) coef_3rd_order = 0.0

      ! set the acoustic step perturbation variables by subtracting the RK timestep variables
      ! from their at the previous RK substep.

      rho_pp = rho_p_save - rho_p
      ru_p = ru_save - ru
      rtheta_pp = rtheta_p_save - rtheta_p
      rtheta_pp_old = rtheta_pp
      rw_p = rw_save - rw

      ! we solve for omega instead of w (see Klemp et al MWR 2007),
      ! so here we change the w_p tendency to an omega_p tendency

      do iCell = 1, nCellsSolve
         do k = 2, nVertLevels
            w_tend(k,iCell) = ( fzm(k) * zz(k,iCell) + fzp(k) * zz(k-1,iCell)   ) * w_tend(k,iCell)
         end do
      end do

      ! here we need to compute the omega tendency in a manner consistent with our diagnosis of omega.
      ! this requires us to use the same flux divergence as is used in the theta eqn - see Klemp et al MWR 2003.

      do iEdge = 1, nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k = 2, nVertLevels
            flux = fzm(k) * u_tend(k,iEdge) + fzp(k) * u_tend(k-1,iEdge)
            w_tend(k,cell2) = w_tend(k,cell2)   &
                     + (zb(k,2,iEdge) + coef_3rd_order * sign(1.0_RKIND, u_tend(k,iEdge)) * zb3(k,2,iEdge)) * flux   &
                     * (fzm(k) * zz(k,cell2) + fzp(k) * zz(k-1,cell2)) 
            w_tend(k,cell1) = w_tend(k,cell1)   &
                     - (zb(k,1,iEdge) + coef_3rd_order * sign(1.0_RKIND, u_tend(k,iEdge)) * zb3(k,1,iEdge)) * flux   &
                     * (fzm(k) * zz(k,cell1) + fzp(k) * zz(k-1,cell1)) 
         end do

      end do

      !  ruAvg and wwAvg will store the mass fluxes averaged over the acoustic steps for the subsequent scalar transport.

      ruAvg(:,:) = 0.0
      wwAvg(:,:) = 0.0

   end subroutine atm_set_smlstep_pert_variables

!-------------------------------

   subroutine atm_advance_acoustic_step( state, diag, tend, mesh, configs, nCells, nVertLevels, dts )

      !  This subroutine performs the entire acoustic step update, following Klemp et al MWR 2007,
      !  using forward-backward vertically implicit integration.  
      !  The gravity-waves are included in the acoustic-step integration.
      !  The input state variables that are updated are ru_p, rw_p (note that this is (rho*omega)_p here),
      !  rtheta_p, and rho_pp.  The time averaged mass flux is accumulated in ruAvg and wwAvg

      implicit none

      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: tend
      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in)    :: configs
      integer, intent(in) :: nCells                  ! for allocating stack variables
      integer, intent(in) :: nVertLevels             ! for allocating stack variables
      real (kind=RKIND), intent(in) :: dts


      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, theta_m, ru_p, rw_p, rtheta_pp,  &
                                                    rtheta_pp_old, zz, exner, cqu, ruAvg,    &
                                                    wwAvg, rho_pp, cofwt, coftz, zx,         &
                                                    a_tri, alpha_tri, gamma_tri, dss,        &
                                                    tend_ru, tend_rho, tend_rt, tend_rw,     &
                                                    zgrid, cofwr, cofwz, w, h_divergence
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp, rdzw, dcEdge, AreaCell, cofrz, dvEdge

      real (kind=RKIND), dimension(:,:), pointer :: cpr, cpl, pzp, pzm
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) :: c2, rcv
      real (kind=RKIND), dimension( nVertLevels ) :: du
      real (kind=RKIND), dimension( nVertLevels + 1 ) :: dpzx
      real (kind=RKIND), dimension( nVertLevels, nCells+1 ) :: ts, rs

      integer :: cell1, cell2, iEdge, iCell, k
      real (kind=RKIND) :: pgrad, flux, resm
      real (kind=RKIND), pointer :: epssm, smdiv

      real (kind=RKIND), pointer :: cf1, cf2, cf3
      real (kind=RKIND) :: pr, pl
      integer :: kr, kl

      integer, pointer :: nEdges, nCellsSolve

      logical, parameter :: debug = .false.
      logical, parameter :: debug1 = .false.
      logical, pointer :: newpx

!--
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(state, 'w', w, 2)

      call mpas_pool_get_array(diag, 'rtheta_pp', rtheta_pp)
      call mpas_pool_get_array(diag, 'rtheta_pp_old', rtheta_pp_old)
      call mpas_pool_get_array(diag, 'h_divergence', h_divergence)
      call mpas_pool_get_array(diag, 'ru_p', ru_p)
      call mpas_pool_get_array(diag, 'rw_p', rw_p)
      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'cqu', cqu)
      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
      call mpas_pool_get_array(diag, 'rho_pp', rho_pp)
      call mpas_pool_get_array(diag, 'cofwt', cofwt)
      call mpas_pool_get_array(diag, 'coftz', coftz)
      call mpas_pool_get_array(diag, 'cofrz', cofrz)
      call mpas_pool_get_array(diag, 'cofwr', cofwr)
      call mpas_pool_get_array(diag, 'cofwz', cofwz)
      call mpas_pool_get_array(diag, 'a_tri', a_tri)
      call mpas_pool_get_array(diag, 'alpha_tri', alpha_tri)
      call mpas_pool_get_array(diag, 'gamma_tri', gamma_tri)

      call mpas_pool_get_array(mesh, 'dss', dss)
      call mpas_pool_get_array(mesh, 'pzp', pzp)
      call mpas_pool_get_array(mesh, 'pzm', pzm)

      call mpas_pool_get_array(tend, 'u', tend_ru)
      call mpas_pool_get_array(tend, 'rho_zz', tend_rho)
      call mpas_pool_get_array(tend, 'theta_m', tend_rt)
      call mpas_pool_get_array(tend, 'w', tend_rw)

      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zx', zx)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)

      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)

      call mpas_pool_get_array(mesh, 'cpr', cpr)
      call mpas_pool_get_array(mesh, 'cpl', cpl)

      call mpas_pool_get_config(configs, 'config_newpx', newpx) 

      ! epssm is the offcentering coefficient for the vertically implicit integration.
      ! smdiv is the 3D divergence-damping coefficient.
      call mpas_pool_get_config(configs, 'config_epssm', epssm) 
      call mpas_pool_get_config(configs, 'config_smdiv', smdiv) 


      rcv = rgas/(cp-rgas)
      c2 = cp*rcv
      resm   = (1.-epssm)/(1.+epssm)

      ts = 0.
      rs = 0.

      ! acoustic step divergence damping - forward weight rtheta_pp - see Klemp et al MWR 2007
      rtheta_pp_old = rtheta_pp + smdiv*(rtheta_pp - rtheta_pp_old)

      if (debug) write(0,*) ' updating ru_p '

      ! forward-backward acoustic step integration.
      ! begin by updating the horizontal velocity u, 
      ! and accumulating the contribution from the updated u to the other tendencies.

      ! we are looping over all edges, but only computing on edges of owned cells. This will include updates of
      ! all owned edges plus some edges that are owned by other blocks.  We perform these redundant computations
      ! so that we do not have to communicate updates of u to update the cell variables (rho, w, and theta). 

      do iEdge = 1, nEdges
 
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! update edges for block-owned cells
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

            if (newpx) then

               k = 1
               pr  =   cpr(k  ,iEdge)*zz(k  ,cell2)*rtheta_pp_old(k  ,cell2)   &
                     + cpr(k+1,iEdge)*zz(k+1,cell2)*rtheta_pp_old(k+1,cell2)   &
                     + cpr(k+2,iEdge)*zz(k+2,cell2)*rtheta_pp_old(k+2,cell2)

               pl  =   cpl(k  ,iEdge)*zz(k  ,cell1)*rtheta_pp_old(k  ,cell1)   &
                     + cpl(k+1,iEdge)*zz(k+1,cell1)*rtheta_pp_old(k+1,cell1)   &
                     + cpl(k+2,iEdge)*zz(k+2,cell1)*rtheta_pp_old(k+2,cell1)
               pgrad = 2./(zz(k,cell1)+zz(k,cell2))*(pr-pl)/dcEdge(iEdge)
               pgrad = 0.5*c2*(exner(k,cell1)+exner(k,cell2))*pgrad
               du(k) = dts*(tend_ru(k,iEdge) - cqu(k,iEdge) * pgrad) 

               do k=2,nVertLevels

                  kr = min(nVertLevels,k+ nint(.5-sign(0.5_RKIND,zx(k,iEdge)+zx(k+1,iEdge))))
                  kl = min(nVertLevels,2*k+1-kr)
                  pr = zz(k,cell2)*rtheta_pp_old(k ,cell2)+.5*(zgrid(k   ,cell1) +zgrid(k +1,cell1)   &
                                                              -zgrid(k   ,cell2) -zgrid(k +1,cell2))  &
                                                             /(zgrid(kr+1,cell2) -zgrid(kr-1,cell2))  &
                    *(zz(kr,cell2)*rtheta_pp_old(kr,cell2)-zz(kr-1,cell2)*rtheta_pp_old(kr-1,cell2))
                  pl = zz(k,cell1)*rtheta_pp_old(k ,cell1)+.5*(zgrid(k   ,cell2) +zgrid(k +1,cell2)   &
                                                              -zgrid(k   ,cell1) -zgrid(k +1,cell1))  &
                                                             /(zgrid(kl+1,cell1) -zgrid(kl-1,cell1))  &
                    *(zz(kl,cell1)*rtheta_pp_old(kl,cell1)-zz(kl-1,cell1)*rtheta_pp_old(kl-1,cell1))
                  pgrad = 2./(zz(k,cell1)+zz(k,cell2))*(pr-pl)/dcEdge(iEdge)
                  pgrad = 0.5*c2*(exner(k,cell1)+exner(k,cell2))*pgrad
                  du(k) = dts*(tend_ru(k,iEdge) - cqu(k,iEdge) * pgrad) 
               end do

            else

               k = 1
               dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                               &
                         *(pzm(k,cell2)*(zz(k+1,cell2)*rtheta_pp_old(k+1,cell2)        &
                                        -zz(k  ,cell2)*rtheta_pp_old(k  ,cell2))       &
                          +pzm(k,cell1)*(zz(k+1,cell1)*rtheta_pp_old(k+1,cell1)        &
                                        -zz(k  ,cell1)*rtheta_pp_old(k  ,cell1))       &
                          +pzp(k,cell2)*(zz(k+2,cell2)*rtheta_pp_old(k+2,cell2)        &
                                        -zz(k  ,cell2)*rtheta_pp_old(k  ,cell2))       &
                          +pzp(k,cell1)*(zz(k+2,cell1)*rtheta_pp_old(k+2,cell1)        &
                                        -zz(k  ,cell1)*rtheta_pp_old(k  ,cell1)))

               do k=2,nVertLevels-1
                  dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                                   &
                                   *(pzp(k,cell2)*(zz(k+1,cell2)*rtheta_pp_old(k+1,cell2)     &
                                                  -zz(k  ,cell2)*rtheta_pp_old(k  ,cell2))    &
                                    +pzm(k,cell2)*(zz(k  ,cell2)*rtheta_pp_old(k  ,cell2)     &
                                                  -zz(k-1,cell2)*rtheta_pp_old(k-1,cell2))    &
                                    +pzp(k,cell1)*(zz(k+1,cell1)*rtheta_pp_old(k+1,cell1)     &
                                                  -zz(k  ,cell1)*rtheta_pp_old(k  ,cell1))    &
                                    +pzm(k,cell1)*(zz(k  ,cell1)*rtheta_pp_old(k  ,cell1)     &
                                                  -zz(k-1,cell1)*rtheta_pp_old(k-1,cell1)))
               end do

               k = nVertLevels
               dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                                   &
                                *(pzm(k,cell2)*(zz(k  ,cell2)*rtheta_pp_old(k  ,cell2)     &
                                               -zz(k-1,cell2)*rtheta_pp_old(k-1,cell2))    &
                                 +pzm(k,cell1)*(zz(k  ,cell1)*rtheta_pp_old(k  ,cell1)     &
                                               -zz(k-1,cell1)*rtheta_pp_old(k-1,cell1)))

               do k=1,nVertLevels
                  pgrad =     ((rtheta_pp_old(k,cell2)*zz(k,cell2)                    &
                               -rtheta_pp_old(k,cell1)*zz(k,cell1))/dcEdge(iEdge)     &
                            -dpzx(k))/(.5*(zz(k,cell2)+zz(k,cell1)))
                  pgrad = 0.5*c2*(exner(k,cell1)+exner(k,cell2))*pgrad
                  du(k) = dts*(tend_ru(k,iEdge) - cqu(k,iEdge) * pgrad) 
               end do
            end if

            do k=1,nVertLevels

               ! full update of ru_p

               ru_p(k,iEdge) = ru_p(k,iEdge) + du(k)

               ! add horizontal fluxes using updated ru_p into density update, rtheta update and w update

               flux = dts*dvEdge(iEdge)*ru_p(k,iEdge)
               rs(k,cell1) = rs(k,cell1)-flux/AreaCell(cell1)
               rs(k,cell2) = rs(k,cell2)+flux/AreaCell(cell2)
   
               flux = flux*0.5*(theta_m(k,cell2)+theta_m(k,cell1))
               ts(k,cell1) = ts(k,cell1)-flux/AreaCell(cell1)
               ts(k,cell2) = ts(k,cell2)+flux/AreaCell(cell2)

               ! accumulate ru_p for use later in scalar transport

               ruAvg(k,iEdge) = ruAvg(k,iEdge) + ru_p(k,iEdge)

            end do

         end if ! end test for block-owned cells

      end do ! end loop over edges

      ! saving rtheta_pp before update for use in divergence damping in next acoustic step

      rtheta_pp_old(:,:) = rtheta_pp(:,:)

      ! vertically implicit acoustic and gravity wave integration.
      ! this follows Klemp et al MWR 2007, with the addition of an implicit Rayleigh damping of w
      ! serves as a gravity-wave absorbing layer, from Klemp et al 2008.

      do iCell = 1, nCellsSolve

         do k=1, nVertLevels
            rs(k,iCell) = rho_pp(k,iCell) + dts*tend_rho(k,iCell) + rs(k,iCell)      &
                            - cofrz(k)*resm*(rw_p(k+1,iCell)-rw_p(k,iCell))
            ts(k,iCell) = rtheta_pp(k,iCell) + dts*tend_rt(k,iCell) + ts(k,iCell)    &
                               - resm*rdzw(k)*(coftz(k+1,iCell)*rw_p(k+1,iCell)      &
                               -coftz(k,iCell)*rw_p(k,iCell))
         end do

         do k=2, nVertLevels

            wwavg(k,iCell) = wwavg(k,iCell) + 0.5*(1.-epssm)*rw_p(k,iCell)

            rw_p(k,iCell) = rw_p(k,iCell) +  dts*tend_rw(k,iCell)                       &
                       - cofwz(k,iCell)*((zz(k  ,iCell)*ts (k  ,iCell)                  &
                                     -zz(k-1,iCell)*ts (k-1,iCell))                     &
                               +resm*(zz(k  ,iCell)*rtheta_pp(k  ,iCell)                &
                                     -zz(k-1,iCell)*rtheta_pp(k-1,iCell)))              &
                       - cofwr(k,iCell)*((rs (k,iCell)+rs (k-1,iCell))                  &
                               +resm*(rho_pp(k,iCell)+rho_pp(k-1,iCell)))               &
                       + cofwt(k  ,iCell)*(ts (k  ,iCell)+resm*rtheta_pp(k  ,iCell))    &
                       + cofwt(k-1,iCell)*(ts (k-1,iCell)+resm*rtheta_pp(k-1,iCell))
         end do

         ! tridiagonal solve sweeping up and then down the column

         do k=2,nVertLevels
            rw_p(k,iCell) = (rw_p(k,iCell)-a_tri(k,iCell)*rw_p(k-1,iCell))*alpha_tri(k,iCell)
         end do

         do k=nVertLevels,1,-1
            rw_p(k,iCell) = rw_p(k,iCell) - gamma_tri(k,iCell)*rw_p(k+1,iCell)     
         end do

         ! the implicit Rayleigh damping on w (gravity-wave absorbing) 

         do k=2,nVertLevels
            rw_p(k,iCell) = (rw_p(k,iCell)-dts*dss(k,iCell)*               &
                        (fzm(k)*zz (k,iCell)+fzp(k)*zz (k-1,iCell))        &
                        *(fzm(k)*rho_zz(k,iCell)+fzp(k)*rho_zz(k-1,iCell))       &
                                 *w(k,iCell)    )/(1.+dts*dss(k,iCell))
 
            ! accumulate (rho*omega)' for use later in scalar transport

            wwAvg(k,iCell) = wwAvg(k,iCell) + 0.5*(1.+epssm)*rw_p(k,iCell)
 
         end do

         ! update rho_pp and theta_pp given updated rw_p

         do k=1,nVertLevels
            rho_pp(k,iCell) = rs(k,iCell) - cofrz(k) *(rw_p(k+1,iCell)-rw_p(k  ,iCell))
            rtheta_pp(k,iCell) = ts(k,iCell) - rdzw(k)*(coftz(k+1,iCell)*rw_p(k+1,iCell)  &
                               -coftz(k  ,iCell)*rw_p(k  ,iCell))
         end do

      end do !  end of loop over cells

   end subroutine atm_advance_acoustic_step

!------------------------

   subroutine atm_recover_large_step_variables( state, diag, tend, mesh, configs, dt, ns, rk_step )

      ! reconstitute state variables from acoustic-step perturbation variables 
      ! after the acoustic steps.  The perturbation variables were originally set in
      ! subroutine atm_set_smlstep_pert_variables prior to their acoustic-steps update.
      ! we are also computing a few other state-derived variables here.

      implicit none

      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: tend
      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: ns, rk_step
      real (kind=RKIND), intent(in) :: dt


      real (kind=RKIND), dimension(:,:), pointer :: wwAvg, rw_save, w, rw, rw_p, rtheta_p, rtheta_pp,   &
                                                    rtheta_p_save, rt_diabatic_tend, rho_p, rho_p_save, &
                                                    rho_pp, rho_zz, rho_base, ruAvg, ru_save, ru_p, u, ru, &
                                                    exner, exner_base, rtheta_base, pressure_p,         &
                                                    zz, theta_m, pressure_b, qvapor
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp, dvEdge, areaCell
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3 
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer :: iCell, iEdge, k, cell1, cell2
      integer, pointer :: nVertLevels, nCells, nCellsSolve, nEdges, nEdgesSolve
      real (kind=RKIND) :: rcv, p0, flux
      real (kind=RKIND), pointer :: cf1, cf2, cf3, coef_3rd_order
      integer, pointer :: config_theta_adv_order
      integer, pointer :: index_qv

      logical, parameter :: debug=.false.


      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
      call mpas_pool_get_array(diag, 'rw_save', rw_save)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'rw_p', rw_p)
      call mpas_pool_get_array(state, 'w', w, 2)

      call mpas_pool_get_array(diag, 'rtheta_p', rtheta_p)
      call mpas_pool_get_array(diag, 'rtheta_p_save', rtheta_p_save)
      call mpas_pool_get_array(diag, 'rtheta_pp', rtheta_pp)
      call mpas_pool_get_array(diag, 'rtheta_base', rtheta_base)
      call mpas_pool_get_array(tend, 'rt_diabatic_tend', rt_diabatic_tend)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)

      call mpas_pool_get_dimension(state, 'index_qv', index_qv)

      qvapor => scalars(index_qv,:,:)   ! MGD does this actually work?

      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(diag, 'rho_p', rho_p)
      call mpas_pool_get_array(diag, 'rho_p_save', rho_p_save)
      call mpas_pool_get_array(diag, 'rho_pp', rho_pp)
      call mpas_pool_get_array(diag, 'rho_base', rho_base)

      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(diag, 'ru_save', ru_save)
      call mpas_pool_get_array(diag, 'ru_p', ru_p)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(state, 'u', u, 2)

      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'exner_base', exner_base)

      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_b)

      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)

      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)

      rcv = rgas/(cp-rgas)
      p0 = 1.e+05  ! this should come from somewhere else...

      if (config_theta_adv_order /=3) coef_3rd_order = 0.0

      ! compute new density everywhere so we can compute u from ru.
      ! we will also need it to compute theta_m below

      do iCell = 1, nCells

         do k = 1, nVertLevels

            rho_p(k,iCell) = rho_p(k,iCell) + rho_pp(k,iCell)

            rho_zz(k,iCell) = rho_p(k,iCell) + rho_base(k,iCell)
         end do

         w(1,iCell) = 0.
         do k = 2, nVertLevels
            wwAvg(k,iCell) = rw(k,iCell) + (wwAvg(k,iCell) / float(ns))

            rw(k,iCell) = rw(k,iCell) + rw_p(k,iCell)


          ! pick up part of diagnosed w from omega
            w(k,iCell) = rw(k,iCell)/( (fzm(k)*zz (k,iCell)+fzp(k)*zz (k-1,iCell))   &
                                      *(fzm(k)*rho_zz(k,iCell)+fzp(k)*rho_zz(k-1,iCell)) )
         end do
         w(nVertLevels+1,iCell) = 0.

         if (rk_step == 3) then
            do k = 1, nVertLevels
               rtheta_p(k,iCell) = rtheta_p(k,iCell) + rtheta_pp(k,iCell) &
                                 - dt * rho_zz(k,iCell) * rt_diabatic_tend(k,iCell)
            end do
         else
            do k = 1, nVertLevels
               rtheta_p(k,iCell) = rtheta_p(k,iCell) + rtheta_pp(k,iCell)
            end do
         end if

         do k = 1, nVertLevels
            theta_m(k,iCell) = (rtheta_p(k,iCell) + rtheta_base(k,iCell))/rho_zz(k,iCell)
            exner(k,iCell) = (zz(k,iCell)*(rgas/p0)*(rtheta_p(k,iCell)+rtheta_base(k,iCell)))**rcv
            ! pressure_p is perturbation pressure
            pressure_p(k,iCell) = zz(k,iCell) * rgas * (exner(k,iCell)*rtheta_p(k,iCell)+rtheta_base(k,iCell)  &
                                                          * (exner(k,iCell)-exner_base(k,iCell)))
         end do

      end do

      ! recover time-averaged ruAvg on all edges of owned cells (for upcoming scalar transport).  
      ! we solved for these in the acoustic-step loop.  
      ! we will compute ru and u here also, given we are here, even though we only need them on nEdgesSolve

      ! Avoid FP errors caused by a potential division by zero below by 
      ! initializing the "garbage cell" of rho_zz to a non-zero value
      rho_zz(:,nCells+1) = 1.0

      do iEdge = 1, nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k = 1, nVertLevels
            ruAvg(k,iEdge) = ru(k,iEdge) + (ruAvg(k,iEdge) / float(ns))
            ru(k,iEdge) = ru(k,iEdge) + ru_p(k,iEdge)
            u(k,iEdge) = 2.*ru(k,iEdge)/(rho_zz(k,cell1)+rho_zz(k,cell2))
         end do

         !  finish recovering w from (rho*omega)_p.  as when we formed (rho*omega)_p from u and w, we need
         !  to use the same flux-divergence operator as is used for the horizontal theta transport
         !  (See Klemp et al 2003).

         flux = cf1*ru(1,iEdge) + cf2*ru(2,iEdge) + cf3*ru(3,iEdge)
         w(1,cell2) = w(1,cell2) - (zb(1,2,iEdge) + sign(1.0_RKIND,flux)*coef_3rd_order*zb3(1,2,iEdge))  &
                                *flux/(cf1*rho_zz(1,cell2)+cf2*rho_zz(2,cell2)+cf3*rho_zz(3,cell2))
         w(1,cell1) = w(1,cell1) + (zb(1,1,iEdge) + sign(1.0_RKIND,flux)*coef_3rd_order*zb3(1,1,iEdge))  &
                                *flux/(cf1*rho_zz(1,cell1)+cf2*rho_zz(2,cell1)+cf3*rho_zz(3,cell1))

         do k = 2, nVertLevels
            flux = (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
            w(k,cell2) = w(k,cell2) - (zb(k,2,iEdge)+sign(1.0_RKIND,flux)*coef_3rd_order*zb3(k,2,iEdge)) &
                                 *flux/(fzm(k)*rho_zz(k,cell2)+fzp(k)*rho_zz(k-1,cell2))
            w(k,cell1) = w(k,cell1) + (zb(k,1,iEdge)+sign(1.0_RKIND,flux)*coef_3rd_order*zb3(k,1,iEdge)) &
                                 *flux/(fzm(k)*rho_zz(k,cell1)+fzp(k)*rho_zz(k-1,cell1))
         end do

      end do

   end subroutine atm_recover_large_step_variables

!---------------------------------------------------------------------------------------

   subroutine atm_advance_scalars( tend, state, diag, mesh, configs, num_scalars, nVertLevels, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Integrate scalar equations - explicit transport plus other tendencies
   !
   ! Input: s - current model state, 
   !            including tendencies from sources other than resolved transport.
   !        grid - grid metadata
   !
   ! input scalars in state are uncoupled (i.e. not mulitplied by density)
   ! 
   ! Output: updated uncoupled scalars (scalars in state).
   ! Note: scalar tendencies are also modified by this routine.
   !
   ! This routine DOES NOT apply any positive definite or monotonic renormalizations.
   !
   ! The transport scheme is from Skamarock and Gassmann MWR 2011.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (mpas_pool_type), intent(in) :: tend
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(in) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: num_scalars      ! for allocating stack variables
      integer, intent(in) :: nVertLevels      ! for allocating stack variables
      real (kind=RKIND) :: dt

      integer :: i, iCell, iEdge, k, iScalar, cell1, cell2
      real (kind=RKIND) :: scalar_weight

      real (kind=RKIND), dimension(:,:,:), pointer :: scalar_old, scalar_new, scalar_tend
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: uhAvg, h_old, h_new, wwAvg, rho_edge, rho_zz, zgrid, kdiff
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell, qv_init
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer, dimension(:,:), pointer :: advCellsForEdge
      integer, dimension(:), pointer :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      real (kind=RKIND), dimension( num_scalars, nVertLevels ) :: flux_arr

      real (kind=RKIND), dimension( num_scalars, nVertLevels + 1 ) :: wdtn
      integer, pointer :: nCellsSolve, nEdges

      real (kind=RKIND), dimension(:), pointer :: fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4
      real (kind=RKIND), pointer :: coef_3rd_order

      real (kind=RKIND), pointer :: h_theta_eddy_visc2, v_theta_eddy_visc2

      real (kind=RKIND) :: flux3, flux4
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3

      integer, pointer :: config_scalar_vadv_order

      integer, parameter :: hadv_opt = 2

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
          ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0

      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_scalar_vadv_order', config_scalar_vadv_order)

      call mpas_pool_get_array(state, 'scalars', scalar_old, 1)
      call mpas_pool_get_array(state, 'scalars', scalar_new, 2)
      call mpas_pool_get_array(state, 'rho_zz', h_old, 1)
      call mpas_pool_get_array(state, 'rho_zz', h_new, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)

      call mpas_pool_get_array(diag, 'kdiff', kdiff)
      call mpas_pool_get_array(diag, 'ruAvg', uhAvg)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)

      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(tend, 'scalars_tend', scalar_tend)

      call mpas_pool_get_array(mesh, 'fzm', fnm)
      call mpas_pool_get_array(mesh, 'fzp', fnp)
      call mpas_pool_get_array(mesh, 'rdzw', rdnw)
      call mpas_pool_get_array(mesh, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(mesh, 'meshScalingDel4', meshScalingDel4)

      call mpas_pool_get_array(mesh, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(mesh, 'advCellsForEdge', advCellsForEdge)
      call mpas_pool_get_array(mesh, 'adv_coefs', adv_coefs)
      call mpas_pool_get_array(mesh, 'adv_coefs_3rd', adv_coefs_3rd)

      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)
      call mpas_pool_get_array(mesh, 'qv_init', qv_init)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)

      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)

      call mpas_pool_get_config(configs, 'config_h_theta_eddy_visc2', h_theta_eddy_visc2)
      call mpas_pool_get_config(configs, 'config_v_theta_eddy_visc2', v_theta_eddy_visc2)


#ifndef DO_PHYSICS
      scalar_tend = 0.  !  testing purposes - we have no sources or sinks
#endif

      !
      ! Runge Kutta integration, so we compute fluxes from scalar_new values, update starts from scalar_old
      !
      !  horizontal flux divergence, accumulate in scalar_tend

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then  ! only for owned cells
  
            !  flux_arr stores the value of the scalar at the edge.
            !  a better name perhaps would be scalarEdge

            flux_arr(:,:) = 0.
            do i=1,nAdvCellsForEdge(iEdge)
               iCell = advCellsForEdge(i,iEdge)
               do k=1,nVertLevels
                  scalar_weight = adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(i,iEdge)
                  do iScalar=1,num_scalars
                     flux_arr(iScalar,k) = flux_arr(iScalar,k) + scalar_weight* scalar_new(iScalar,k,iCell)
                  end do
               end do
            end do

            ! here we add the horizontal flux divergence into the scalar tendency.
            ! note that the scalar tendency is modified.

            do k=1,nVertLevels
            do iScalar=1,num_scalars
                  scalar_tend(iScalar,k,cell1) = scalar_tend(iScalar,k,cell1) &
                         - uhAvg(k,iEdge)*flux_arr(iScalar,k)/areaCell(cell1)
                  scalar_tend(iScalar,k,cell2) = scalar_tend(iScalar,k,cell2) &
                         + uhAvg(k,iEdge)*flux_arr(iScalar,k)/areaCell(cell2)
            end do
            end do

         end if
      end do

      !
      !  vertical flux divergence and update of the scalars
      !

      ! zero fluxes at top and bottom

      wdtn(:,1) = 0.
      wdtn(:,nVertLevels+1) = 0.

      if (config_scalar_vadv_order == 2) then

         do iCell=1,nCellsSolve
            do k = 2, nVertLevels
               do iScalar=1,num_scalars
                  wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
               end do
            end do
            do k=1,nVertLevels
               do iScalar=1,num_scalars
                  scalar_new(iScalar,k,iCell) = (   scalar_old(iScalar,k,iCell)*h_old(k,iCell) &
                        + dt*( scalar_tend(iScalar,k,iCell) -rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k)) ) )/h_new(k,iCell)
               end do
            end do
         end do

      else if ( config_scalar_vadv_order == 3 ) then

         do iCell=1,nCellsSolve

            k = 2
            do iScalar=1,num_scalars
               wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
            end do
             
            do k=3,nVertLevels-1
               do iScalar=1,num_scalars
                  wdtn(iScalar,k) = flux3( scalar_new(iScalar,k-2,iCell),scalar_new(iScalar,k-1,iCell),  &
                                           scalar_new(iScalar,k  ,iCell),scalar_new(iScalar,k+1,iCell),  &
                                           wwAvg(k,iCell), coef_3rd_order )
               end do
            end do
            k = nVertLevels
            do iScalar=1,num_scalars
               wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
            end do

            do k=1,nVertLevels
               do iScalar=1,num_scalars
                  scalar_new(iScalar,k,iCell) = (   scalar_old(iScalar,k,iCell)*h_old(k,iCell) &
                        + dt*( scalar_tend(iScalar,k,iCell) -rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k)) ) )/h_new(k,iCell)
               end do
            end do

         end do

      else if ( config_scalar_vadv_order == 4 ) then

         do iCell=1,nCellsSolve

            k = 2
            do iScalar=1,num_scalars
               wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
            end do
            do k=3,nVertLevels-1
               do iScalar=1,num_scalars
                  wdtn(iScalar,k) = flux4( scalar_new(iScalar,k-2,iCell),scalar_new(iScalar,k-1,iCell),  &
                                           scalar_new(iScalar,k  ,iCell),scalar_new(iScalar,k+1,iCell), wwAvg(k,iCell) )
               end do
            end do
            k = nVertLevels
            do iScalar=1,num_scalars
               wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
            end do

            do k=1,nVertLevels
               do iScalar=1,num_scalars
                  scalar_new(iScalar,k,iCell) = (   scalar_old(iScalar,k,iCell)*h_old(k,iCell) &
                        + dt*( scalar_tend(iScalar,k,iCell) -rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k)) ) )/h_new(k,iCell)
               end do
            end do

         end do
                                                                                        
      else 

         write(0,*) ' bad value for config_scalar_vadv_order - ', config_scalar_vadv_order

      end if

   end subroutine atm_advance_scalars

!---------------------------

   subroutine atm_advance_scalars_mono(block, tend, state, diag, mesh, configs, nCells, nEdges, nVertLevels, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Integrate scalar equations - transport plus other tendencies
   !
   ! Input: s - current model state, 
   !            including tendencies from sources other than resolved transport.
   !        grid - grid metadata
   !
   ! input scalars in state are uncoupled (i.e. not mulitplied by density)
   ! 
   ! Output: updated uncoupled scalars (scalars in s_new).
   ! Note: scalar tendencies are also modified by this routine.
   !
   ! This routine DOES apply positive definite or monotonic renormalizations.
   !
   ! The transport scheme is from Skamarock and Gassmann MWR 2011.
   !
   ! The positive-definite or monotonic renormalization is from Zalesak JCP 1979
   !   as used in the RK3 scheme as described in Wang et al MWR 2009
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (block_type), intent(inout), target :: block
      type (mpas_pool_type), intent(in)    :: tend
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(in)    :: diag
      type (mpas_pool_type), intent(in)    :: mesh
      type (mpas_pool_type), intent(in)    :: configs
      integer, intent(in)                  :: nCells           ! for allocating stack variables
      integer, intent(in)                  :: nEdges           ! for allocating stack variables
      integer, intent(in)                  :: nVertLevels      ! for allocating stack variables
      real (kind=RKIND), intent(in)        :: dt

      integer :: i, iCell, iEdge, k, iScalar, cell1, cell2
      real (kind=RKIND) :: flux, scalar_weight

      real (kind=RKIND), dimension(:,:,:), pointer :: scalar_tend
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: uhAvg, h_old, h_new, wwAvg, rho_edge, rho_zz, zgrid, kdiff
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell, qv_init
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell

      integer, dimension(:,:), pointer :: advCellsForEdge
      integer, dimension(:), pointer :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars_old, scalars_new
      type (field3DReal), pointer :: scalars_old_field

      type (field3DReal), pointer :: tempField
      type (field3DReal), target :: tempFieldTarget

      real (kind=RKIND), dimension( nVertLevels, nCells ) :: scalar_old, scalar_new
      real (kind=RKIND), dimension( nVertLevels, nCells ) :: s_max, s_min
      real (kind=RKIND), dimension( 2, nVertLevels, nCells ), target :: scale_arr

      integer, parameter :: SCALE_IN = 1, SCALE_OUT = 2

      real (kind=RKIND), dimension( nVertLevels, nEdges ) :: flux_arr
      real (kind=RKIND), dimension( nVertLevels + 1, nCells ) :: wdtn

      integer, pointer :: nCellsSolve, num_scalars
      integer :: icellmax, kmax

      real (kind=RKIND), dimension(:), pointer :: fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4
      integer, dimension(:), pointer :: nEdgesOnCell
      real (kind=RKIND), pointer :: coef_3rd_order

      real (kind=RKIND), pointer :: h_theta_eddy_visc2, v_theta_eddy_visc2

      real (kind=RKIND) :: flux3, flux4, flux_upwind
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3, scmin,scmax
      real (kind=RKIND) :: s_min_update, s_max_update, s_upwind, scale_factor

      integer, parameter :: hadv_opt = 2
      real (kind=RKIND), parameter :: eps=1.e-20
      logical, parameter :: debug_print = .false.

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
          ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0


      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_h_theta_eddy_visc2', h_theta_eddy_visc2)
      call mpas_pool_get_config(configs, 'config_v_theta_eddy_visc2', v_theta_eddy_visc2)

      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)

      call mpas_pool_get_array(diag, 'kdiff', kdiff)
      call mpas_pool_get_array(diag, 'ruAvg', uhAvg)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)

      call mpas_pool_get_array(tend, 'scalars_tend', scalar_tend)

      call mpas_pool_get_array(state, 'rho_zz', h_old, 1)
      call mpas_pool_get_array(state, 'rho_zz', h_new, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(state, 'scalars', scalars_old, 1)
      call mpas_pool_get_array(state, 'scalars', scalars_new, 2)
      call mpas_pool_get_field(state, 'scalars', scalars_old_field, 1)

      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'fzm', fnm)
      call mpas_pool_get_array(mesh, 'fzp', fnp)
      call mpas_pool_get_array(mesh, 'rdzw', rdnw)
      call mpas_pool_get_array(mesh, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(mesh, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(mesh, 'advCellsForEdge', advCellsForEdge)
      call mpas_pool_get_array(mesh, 'adv_coefs', adv_coefs)
      call mpas_pool_get_array(mesh, 'adv_coefs_3rd', adv_coefs_3rd)
      call mpas_pool_get_array(mesh, 'qv_init', qv_init)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)

#ifndef DO_PHYSICS
      scalar_tend = 0.  !  testing purposes - we have no sources or sinks
#endif

      !  for positive-definite or monotonic option, we first update scalars using the tendency from sources other than
      !  the resolved transport (these should constitute a positive definite update).  
      !  Note, however, that we enforce positive-definiteness in this update.
      !  The transport will maintain this positive definite solution and optionally, shape preservation (monotonicity).

      do iCell = 1, nCellsSolve
      do k = 1, nVertLevels
      do iScalar = 1,num_scalars
         scalars_old(iScalar,k,iCell) = scalars_old(iScalar,k,iCell)+dt*scalar_tend(iScalar,k,iCell) / h_old(k,iCell)
         scalar_tend(iScalar,k,iCell) = 0.
      end do
      end do
      end do

      !  halo exchange

      call mpas_dmpar_exch_halo_field(scalars_old_field)

      !
      ! Runge Kutta integration, so we compute fluxes from scalar_new values, update starts from scalar_old
      !

      ! do one scalar at a time

      do iScalar = 1, num_scalars
         write(0,*) ' mono transport for scalar ',iScalar

         do iCell = 1, nCells
         do k = 1, nVertLevels
            scalar_old(k,iCell) = scalars_old(iScalar,k,iCell)
            scalar_new(k,iCell) = scalars_new(iScalar,k,iCell)
         end do
         end do

         if (debug_print) then
            scmin = scalar_old(1,1)
            scmax = scalar_old(1,1)
            do iCell = 1, nCells
            do k=1, nVertLevels
               scmin = min(scmin,scalar_old(k,iCell))
               scmax = max(scmax,scalar_old(k,iCell))
            end do
            end do
            write(0,*) ' scmin, scmin old in ',scmin,scmax

            scmin = scalar_new(1,1)
            scmax = scalar_new(1,1)
            do iCell = 1, nCells
            do k=1, nVertLevels
               scmin = min(scmin,scalar_new(k,iCell))
               scmax = max(scmax,scalar_new(k,iCell))
            end do
            end do
            write(0,*) ' scmin, scmin new in ',scmin,scmax
         end if


      !
      !  vertical flux divergence, and min and max bounds for flux limiter
      !

 
         do iCell=1,nCellsSolve

            ! zero flux at top and bottom
            wdtn(1,iCell) = 0.
            wdtn(nVertLevels+1,iCell) = 0.

            k = 1
            s_max(k,iCell) = max(scalar_old(1,iCell),scalar_old(2,iCell))
            s_min(k,iCell) = min(scalar_old(1,iCell),scalar_old(2,iCell))

            k = 2
            wdtn(k,iCell) = wwAvg(k,iCell)*(fnm(k)*scalar_new(k,iCell)+fnp(k)*scalar_new(k-1,iCell))
            s_max(k,iCell) = max(scalar_old(k-1,iCell),scalar_old(k,iCell),scalar_old(k+1,iCell))
            s_min(k,iCell) = min(scalar_old(k-1,iCell),scalar_old(k,iCell),scalar_old(k+1,iCell))
             
            do k=3,nVertLevels-1
               wdtn(k,iCell) = flux3( scalar_new(k-2,iCell),scalar_new(k-1,iCell),  &
                                      scalar_new(k  ,iCell),scalar_new(k+1,iCell),  &
                                      wwAvg(k,iCell), coef_3rd_order )
               s_max(k,iCell) = max(scalar_old(k-1,iCell),scalar_old(k,iCell),scalar_old(k+1,iCell))
               s_min(k,iCell) = min(scalar_old(k-1,iCell),scalar_old(k,iCell),scalar_old(k+1,iCell))
            end do
 
            k = nVertLevels
            wdtn(k,iCell) = wwAvg(k,iCell)*(fnm(k)*scalar_new(k,iCell)+fnp(k)*scalar_new(k-1,iCell))
            s_max(k,iCell) = max(scalar_old(k,iCell),scalar_old(k-1,iCell))
            s_min(k,iCell) = min(scalar_old(k,iCell),scalar_old(k-1,iCell))

      ! pull s_min and s_max from the (horizontal) surrounding cells

            do i=1, nEdgesOnCell(iCell)
               do k=1, nVertLevels
                  s_max(k,iCell) = max(s_max(k,iCell),scalar_old(k, cellsOnCell(i,iCell)))
                  s_min(k,iCell) = min(s_min(k,iCell),scalar_old(k, cellsOnCell(i,iCell)))
               end do
            end do

         end do

      !
      !  horizontal flux divergence

         flux_arr(:,:) = 0.
         do iEdge=1,nEdges

            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then  ! only for owned cells
  
               do i=1,nAdvCellsForEdge(iEdge)
                  iCell = advCellsForEdge(i,iEdge)
                  do k=1,nVertLevels
                     scalar_weight = uhAvg(k,iEdge)*(adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(i,iEdge))
                     flux_arr(k,iEdge) = flux_arr(k,iEdge) + scalar_weight* scalar_new(k,iCell)
                  end do
               end do

            end if

         end do

!  vertical flux divergence for upwind update, we will put upwind update into scalar_new, and put factor of dt in fluxes

         do iCell = 1, nCellsSolve

            k = 1
            scalar_new(k,iCell) = scalar_old(k,iCell)*h_old(k,iCell)

            do k = 2, nVertLevels
               scalar_new(k,iCell) = scalar_old(k,iCell)*h_old(k,iCell)
               flux_upwind = dt*(max(0.0_RKIND,wwAvg(k,iCell))*scalar_old(k-1,iCell) + min(0.0_RKIND,wwAvg(k,iCell))*scalar_old(k,iCell))
               scalar_new(k-1,iCell) = scalar_new(k-1,iCell) - flux_upwind*rdnw(k-1)
               scalar_new(k  ,iCell) = scalar_new(k  ,iCell) + flux_upwind*rdnw(k)
               wdtn(k,iCell) = dt*wdtn(k,iCell) - flux_upwind
            end do

! scale_arr(SCALE_IN,:,:) and scale_arr(SCALE_OUT:,:) are used here to store the incoming and outgoing perturbation flux 
! contributions to the update:  first the vertical flux component, then the horizontal

            do k=1,nVertLevels
               scale_arr(SCALE_IN, k,iCell) = - rdnw(k)*(min(0.0_RKIND,wdtn(k+1,iCell))-max(0.0_RKIND,wdtn(k,iCell)))
               scale_arr(SCALE_OUT,k,iCell) = - rdnw(k)*(max(0.0_RKIND,wdtn(k+1,iCell))-min(0.0_RKIND,wdtn(k,iCell)))
            end do

         end do

!  horizontal flux divergence for upwind update

         !  upwind flux computation

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then  ! only for owned cells
               do k=1, nVertLevels
                  flux_upwind = dvEdge(iEdge) * dt *   &
                         (max(0.0_RKIND,uhAvg(k,iEdge))*scalar_old(k,cell1) + min(0.0_RKIND,uhAvg(k,iEdge))*scalar_old(k,cell2))
                  flux_arr(k,iEdge) = dt*flux_arr(k,iEdge) - flux_upwind
                  scalar_new(k,cell1) = scalar_new(k,cell1) - flux_upwind / areaCell(cell1)
                  scalar_new(k,cell2) = scalar_new(k,cell2) + flux_upwind / areaCell(cell2)
 
                  scale_arr(SCALE_OUT,k,cell1) = scale_arr(SCALE_OUT,k,cell1) - max(0.0_RKIND,flux_arr(k,iEdge)) / areaCell(cell1)
                  scale_arr(SCALE_IN, k,cell1) = scale_arr(SCALE_IN, k,cell1) - min(0.0_RKIND,flux_arr(k,iEdge)) / areaCell(cell1)
                  scale_arr(SCALE_OUT,k,cell2) = scale_arr(SCALE_OUT,k,cell2) + min(0.0_RKIND,flux_arr(k,iEdge)) / areaCell(cell2)
                  scale_arr(SCALE_IN, k,cell2) = scale_arr(SCALE_IN, k,cell2) + max(0.0_RKIND,flux_arr(k,iEdge)) / areaCell(cell2)
 
               end do
            end if
         end do

!  next, the limiter

         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               s_min_update = (scalar_new(k,iCell)+scale_arr(SCALE_OUT,k,iCell))/h_new(k,iCell)
               s_max_update = (scalar_new(k,iCell)+scale_arr(SCALE_IN,k,iCell))/h_new(k,iCell)
               s_upwind = scalar_new(k,iCell)/h_new(k,iCell)

#ifdef SINGLE_PRECISION_IFORT
               ! ifort 13.0.1 seems to have problems optimizing this particular
               ! bit of code without generating NaNs in single precision, so we write it differently
               if ((s_max_update - s_upwind) /= 0.0) then
                  scale_factor = (s_max(k,iCell)-s_upwind)/(s_max_update-s_upwind)
                  scale_arr(SCALE_IN,k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
               else
                  scale_arr(SCALE_IN,k,iCell) = 1.0
               end if

               if ((s_upwind - s_min_update) /= 0.0) then
                  scale_factor = (s_upwind-s_min(k,iCell))/(s_upwind-s_min_update)
                  scale_arr(SCALE_OUT,k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
               else
                  scale_arr(SCALE_OUT,k,iCell) = 1.0
               end if
#else
               scale_factor = (s_max(k,iCell)-s_upwind)/(s_max_update-s_upwind+eps)
               scale_arr(SCALE_IN,k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
 
               scale_factor = (s_upwind-s_min(k,iCell))/(s_upwind-s_min_update+eps)
               scale_arr(SCALE_OUT,k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
#endif
            end do
         end do

!
!  communicate scale factors here.
!  communicate only first halo row in these next two exchanges
!
         tempField => tempFieldTarget

         tempField % block => block
         tempField % dimSizes(1) = 2
         tempField % dimSizes(2) = nVertLevels
         tempField % dimSizes(3) = nCells
         tempField % sendList => block % parinfo % cellsToSend
         tempField % recvList => block % parinfo % cellsToRecv
         tempField % copyList => block % parinfo % cellsToCopy
         tempField % prev => null()
         tempField % next => null()

         tempField % array => scale_arr
         call mpas_dmpar_exch_halo_field(tempField, (/ 1 /))

!
!  rescale the fluxes
!
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
               do k = 1, nVertLevels
                  flux = flux_arr(k,iEdge)
                  flux = max(0.0_RKIND,flux) * min(scale_arr(SCALE_OUT,k,cell1), scale_arr(SCALE_IN, k,cell2)) &
                       + min(0.0_RKIND,flux) * min(scale_arr(SCALE_IN, k,cell1), scale_arr(SCALE_OUT,k,cell2))
                  flux_arr(k,iEdge) = flux
               end do
            end if
         end do
 
       ! rescale the vertical flux
 
         do iCell=1,nCellsSolve
            do k = 2, nVertLevels
               flux =  wdtn(k,iCell)
               flux = max(0.0_RKIND,flux) * min(scale_arr(SCALE_OUT,k-1,iCell), scale_arr(SCALE_IN,k  ,iCell)) &
                    + min(0.0_RKIND,flux) * min(scale_arr(SCALE_OUT,k  ,iCell), scale_arr(SCALE_IN,k-1,iCell))
               wdtn(k,iCell) = flux
            end do
         end do
!
!  do the scalar update now that we have the fluxes
!
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then  ! only for owned cells
               do k=1,nVertLevels
                  scalar_new(k,cell1) = scalar_new(k,cell1) - flux_arr(k,iEdge)/areaCell(cell1)
                  scalar_new(k,cell2) = scalar_new(k,cell2) + flux_arr(k,iEdge)/areaCell(cell2)
               end do
            end if
         end do

         do iCell=1,nCellsSolve
            do k=1,nVertLevels
               scalar_new(k,iCell) = (   scalar_new(k,iCell)  &
                   + (-rdnw(k)*(wdtn(k+1,iCell)-wdtn(k,iCell)) ) )/h_new(k,iCell)
            end do
         end do

         if(debug_print) then

            scmin = scalar_new(1,1)
            scmax = scalar_new(1,1)
            do iCell = 1, nCellsSolve
            do k=1, nVertLevels
               scmax = max(scmax,scalar_new(k,iCell))
               scmin = min(scmin,scalar_new(k,iCell))
               if (s_max(k,iCell) < scalar_new(k,iCell)) then
                  write(32,*) ' over - k,iCell,s_min,s_max,scalar_new ',k,iCell,s_min(k,iCell),s_max(k,iCell),scalar_new(k,iCell)
               end if
               if (s_min(k,iCell) > scalar_new(k,iCell)) then
                  write(32,*) ' under - k,iCell,s_min,s_max,scalar_new ',k,iCell,s_min(k,iCell),s_max(k,iCell),scalar_new(k,iCell)
               end if
            end do
            end do
            write(0,*) ' scmin, scmax new out ',scmin,scmax
            write(0,*) ' icell_min, k_min ',icellmax, kmax

         end if

         ! the update should be positive definite. but roundoff can sometimes leave small negative values
         ! hence the enforcement of PD in the copy back to the model state.

         do iCell = 1, nCells
         do k=1, nVertLevels
            scalars_new(iScalar,k,iCell) = max(0.0_RKIND,scalar_new(k,iCell))
         end do
         end do

      end do !  loop over scalars

   end subroutine atm_advance_scalars_mono

!----

   subroutine atm_compute_dyn_tend(tend, state, diag, mesh, configs, nVertLevels, rk_step, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute height and normal wind tendencies, as well as diagnostic variables
   !
   ! Input: state - current model state
   !        mesh - grid metadata
   !        diag - some grid diagnostics
   !
   ! Output: tend - tendencies: tend_u, tend_w, tend_theta and tend_rho
   !                these are all coupled-variable tendencies.
   !         various other quantities in diag: Smagorinsky eddy viscosity
   !                
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (mpas_pool_type), intent(inout) :: tend
      type (mpas_pool_type), intent(in) :: state
      type (mpas_pool_type), intent(in) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: nVertLevels              ! for allocating stack variables
      integer, intent(in) :: rk_step
      real (kind=RKIND), intent(in) :: dt


      logical, parameter :: rk_diffusion = .false.

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, iq
      real (kind=RKIND) :: flux, workpv

      integer, pointer :: nCells, nEdges, nVertices, nCellsSolve, nEdgesSolve
      integer, pointer :: moist_start, moist_end
      real (kind=RKIND), pointer :: h_mom_eddy_visc2, v_mom_eddy_visc2
      real (kind=RKIND), pointer :: h_theta_eddy_visc2, v_theta_eddy_visc2
      real (kind=RKIND) :: h_mom_eddy_visc4
      real (kind=RKIND) :: h_theta_eddy_visc4
      real (kind=RKIND) :: u_diffusion
      real (kind=RKIND), dimension(:), pointer ::  fEdge, dvEdge, dcEdge, areaCell, areaTriangle, meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge, zgrid, rho_edge, rho_zz, ru, u, v, tend_u, &
                                                    circulation, divergence, vorticity, ke, pv_edge, theta_m, rw, tend_rho, &
                                                    rt_diabatic_tend, tend_theta, tend_w, w, cqw, rb, rr, pp, pressure_b, zz, zx, cqu, & 
                                                    h_divergence, kdiff
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      real (kind=RKIND), dimension(:,:), pointer :: tend_u_euler, tend_w_euler, tend_theta_euler

      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgesOnCell, edgesOnEdge, cellsOnCell
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, latEdge, angleEdge, u_init

      real (kind=RKIND), dimension( nVertLevels + 1 ) :: wduz, wdwz, wdtz, dpzx
      real (kind=RKIND), dimension( nVertLevels ) :: u_mix, ru_edge_w, q
      real (kind=RKIND) :: theta_turb_flux, z1, z2, z3, z4, zm, z0, zp, r
      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2

      integer, dimension(:,:), pointer :: advCellsForEdge
      integer, dimension(:), pointer :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      real (kind=RKIND) :: scalar_weight

      real (kind=RKIND), dimension(:), pointer :: rdzu, rdzw, fzm, fzp, qv_init
      real (kind=RKIND), dimension(:,:), pointer :: t_init 

      real (kind=RKIND), dimension(:,:), pointer :: cpr, cpl, pzp, pzm
      integer :: kr, kl

      real (kind=RKIND), allocatable, dimension(:,:) :: divergence_ru, qtot 
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_theta, delsq_divergence
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity
      real (kind=RKIND), pointer :: cf1, cf2, cf3
      real (kind=RKIND) :: pr, pl
      real (kind=RKIND) :: prandtl_inv

      logical, parameter :: debug = .false.

      logical, parameter :: curvature = .true.
      real (kind=RKIND), pointer :: r_earth
      real (kind=RKIND), dimension(:,:), pointer :: ur_cell, vr_cell

      real (kind=RKIND), parameter :: c_s = 0.125
!      real (kind=RKIND), parameter :: c_s = 0.25
      real (kind=RKIND), dimension( nVertLevels ) :: d_diag, d_off_diag, flux_arr
      real (kind=RKIND), dimension(:,:), pointer :: defc_a, defc_b
      logical :: delsq_horiz_mixing

      real(kind=RKIND), dimension(:,:), pointer :: tend_w_pgf, tend_w_buoy

      real (kind=RKIND), pointer :: coef_3rd_order
      logical, pointer :: newpx
      logical, pointer :: config_mix_full
      integer, pointer :: config_theta_vadv_order
      character (len=StrKIND), pointer :: config_horiz_mixing
      integer, pointer :: config_theta_adv_order
      integer, pointer :: config_w_vadv_order
      integer, pointer :: config_w_adv_order
      real (kind=RKIND), pointer :: config_del4u_div_factor
      real (kind=RKIND), pointer :: config_h_theta_eddy_visc4
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc4
      real (kind=RKIND), pointer :: config_visc4_2dsmag
      real (kind=RKIND), pointer :: config_len_disp

      real (kind=RKIND) :: flux3, flux4
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
                ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0

!-----------

      call mpas_pool_get_config(mesh, 'sphere_radius', r_earth)
      call mpas_pool_get_config(configs, 'config_newpx', newpx)
      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_mix_full', config_mix_full)
      call mpas_pool_get_config(configs, 'config_theta_vadv_order', config_theta_vadv_order)
      call mpas_pool_get_config(configs, 'config_horiz_mixing', config_horiz_mixing)
      call mpas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)
      call mpas_pool_get_config(configs, 'config_w_vadv_order', config_w_vadv_order)
      call mpas_pool_get_config(configs, 'config_w_adv_order', config_w_adv_order)
      call mpas_pool_get_config(configs, 'config_del4u_div_factor', config_del4u_div_factor)
      call mpas_pool_get_config(configs, 'config_h_theta_eddy_visc4', config_h_theta_eddy_visc4)
      call mpas_pool_get_config(configs, 'config_h_mom_eddy_visc4', config_h_mom_eddy_visc4)
      call mpas_pool_get_config(configs, 'config_visc4_2dsmag', config_visc4_2dsmag)
      call mpas_pool_get_config(configs, 'config_len_disp', config_len_disp)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(state, 'u', u, 2)
      call mpas_pool_get_array(state, 'w', w, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)

      call mpas_pool_get_array(diag, 'uReconstructZonal', ur_cell)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', vr_cell)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)
      call mpas_pool_get_array(diag, 'rho_base', rb)
      call mpas_pool_get_array(diag, 'rho_p', rr)
      call mpas_pool_get_array(diag, 'v', v)
      call mpas_pool_get_array(diag, 'kdiff', kdiff)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'circulation', circulation)
      call mpas_pool_get_array(diag, 'divergence', divergence)
      call mpas_pool_get_array(diag, 'vorticity', vorticity)
      call mpas_pool_get_array(diag, 'ke', ke)
      call mpas_pool_get_array(diag, 'pv_edge', pv_edge)
      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
      call mpas_pool_get_array(diag, 'h_divergence', h_divergence)

      call mpas_pool_get_array(mesh, 'pzp', pzp)
      call mpas_pool_get_array(mesh, 'pzm', pzm)
      call mpas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'fEdge', fEdge)
      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zx', zx)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'latEdge', latEdge)
      call mpas_pool_get_array(mesh, 'angleEdge', angleEdge)
      call mpas_pool_get_array(mesh, 'defc_a', defc_a)
      call mpas_pool_get_array(mesh, 'defc_b', defc_b)
      call mpas_pool_get_array(mesh, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(mesh, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(mesh, 'u_init', u_init)
      call mpas_pool_get_array(mesh, 't_init', t_init)
      call mpas_pool_get_array(mesh, 'qv_init', qv_init)

      call mpas_pool_get_array(mesh, 'rdzu', rdzu)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'cpr', cpr)
      call mpas_pool_get_array(mesh, 'cpl', cpl)

      call mpas_pool_get_array(tend, 'u', tend_u)
      call mpas_pool_get_array(tend, 'theta_m', tend_theta)
      call mpas_pool_get_array(tend, 'w', tend_w)
      call mpas_pool_get_array(tend, 'rho_zz', tend_rho)
      call mpas_pool_get_array(tend, 'rt_diabatic_tend', rt_diabatic_tend)
      call mpas_pool_get_array(tend, 'u_euler', tend_u_euler)
      call mpas_pool_get_array(tend, 'theta_euler', tend_theta_euler)
      call mpas_pool_get_array(tend, 'w_euler', tend_w_euler)
      call mpas_pool_get_array(tend, 'w_pgf', tend_w_pgf)
      call mpas_pool_get_array(tend, 'w_buoy', tend_w_buoy)

      call mpas_pool_get_array(diag, 'cqw', cqw)
      call mpas_pool_get_array(diag, 'cqu', cqu)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_dimension(state, 'moist_start', moist_start)
      call mpas_pool_get_dimension(state, 'moist_end', moist_end)

      call mpas_pool_get_config(configs, 'config_h_mom_eddy_visc2', h_mom_eddy_visc2)
      call mpas_pool_get_config(configs, 'config_v_mom_eddy_visc2', v_mom_eddy_visc2)
      call mpas_pool_get_config(configs, 'config_h_theta_eddy_visc2', h_theta_eddy_visc2)
      call mpas_pool_get_config(configs, 'config_v_theta_eddy_visc2', v_theta_eddy_visc2)

      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(mesh, 'advCellsForEdge', advCellsForEdge)
      call mpas_pool_get_array(mesh, 'adv_coefs', adv_coefs)
      call mpas_pool_get_array(mesh, 'adv_coefs_3rd', adv_coefs_3rd)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)


      prandtl_inv = 1.0_RKIND/prandtl

      write(0,*) ' rk_step in compute_dyn_tend ',rk_step


      delsq_horiz_mixing = .false.
      if (config_horiz_mixing == "2d_smagorinsky" .and. (rk_step == 1 .or. rk_diffusion)) then

         ! Smagorinsky eddy viscosity, based on horizontal deformation (in this case on model coordinate surfaces).
         ! The integration coefficients were precomputed and stored in defc_a and defc_b

         do iCell = 1, nCells
            d_diag(:) = 0.
            d_off_diag(:) = 0.
            do iEdge = 1, nEdgesOnCell(iCell)
               do k=1, nVertLevels
                  d_diag(k)     = d_diag(k)     + defc_a(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell))  &
                                                - defc_b(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
                  d_off_diag(k) = d_off_diag(k) + defc_b(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell))  &
                                                + defc_a(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
               end do
            end do
            do k=1, nVertLevels
               ! here is the Smagorinsky formulation, 
               ! followed by imposition of an upper bound on the eddy viscosity
               kdiff(k,iCell) = (c_s * config_len_disp)**2 * sqrt(d_diag(k)**2 + d_off_diag(k)**2)
               kdiff(k,iCell) = min(kdiff(k,iCell),(0.01*config_len_disp**2)/dt)
            end do
         end do
!ldf (2012-10-10):
         h_mom_eddy_visc4   = config_visc4_2dsmag * config_len_disp**3
         h_theta_eddy_visc4 = h_mom_eddy_visc4
         delsq_horiz_mixing = .true.
         write(0,*) '... config_visc4_2dsmag = ', config_visc4_2dsmag
         write(0,*) '... h_mom_eddy_visc4    = ', h_mom_eddy_visc4
         write(0,*) '... h_theta_eddy_visc4  = ', h_theta_eddy_visc4
      else if ( config_horiz_mixing == "2d_fixed") then
         h_mom_eddy_visc4   = config_h_mom_eddy_visc4
         h_theta_eddy_visc4 = config_h_theta_eddy_visc4
         delsq_horiz_mixing = .true.
!ldf (2012-10-10):
      end if

      tend_u(:,:) = 0.0

      ! tendency for density.
      ! accumulate total water here for later use in w tendency calculation.

      allocate(divergence_ru(nVertLevels, nCells+1))
      allocate(qtot(nVertLevels, nCells+1))

      divergence_ru(:,:) = 0.0
      h_divergence(:,:) = 0.

      ! accumulate horizontal mass-flux

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
            flux = ru(k,iEdge)*dvEdge(iEdge)
            divergence_ru(k,cell1) = divergence_ru(k,cell1) + flux
            divergence_ru(k,cell2) = divergence_ru(k,cell2) - flux
         end do
      end do

      qtot(:,:)=0.

      ! compute horiontal mass-flux divergence, add vertical mass flux divergence to complete tend_rho

      do iCell = 1,nCells
         r = 1.0 / areaCell(iCell)
         do k = 1,nVertLevels
            divergence_ru(k,iCell) = divergence_ru(k,iCell) * r
            h_divergence(k,iCell) = divergence_ru(k,iCell)
            tend_rho(k,iCell) = -divergence_ru(k,iCell)-rdzw(k)*(rw(k+1,iCell)-rw(k,iCell))

            do iq = moist_start, moist_end
               qtot(k,iCell) = qtot(k,iCell) + scalars(iq, k, iCell)
            end do

         end do
      end do    

      !
      ! Compute u (normal) velocity tendency for each edge (cell face)
      !

      do iEdge=1,nEdgesSolve

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! horizontal pressure gradient 

         if (newpx)  then

            k = 1
            pr  = cpr(k,iEdge)*pp(k,cell2)+cpr(k+1,iEdge)*pp(k+1,cell2)+cpr(k+2,iEdge)*pp(k+2,cell2)
            pl  = cpl(k,iEdge)*pp(k,cell1)+cpl(k+1,iEdge)*pp(k+1,cell1)+cpl(k+2,iEdge)*pp(k+2,cell1)
            tend_u(k,iEdge) =  - cqu(k,iEdge)*2./(zz(k,cell1)+zz(k,cell2))*(pr-pl)/dcEdge(iEdge)

            do k=2,nVertLevels

               kr = min(nVertLevels,k+ nint(.5-sign(0.5_RKIND,zx(k,iEdge)+zx(k+1,iEdge))))
               kl = min(nVertLevels,2*k+1-kr)

               pr = pp(k,cell2)+.5*(zgrid(k   ,cell1)+zgrid(k +1,cell1)-zgrid(k ,cell2)-zgrid(k +1,cell2))   &
                                  /(zgrid(kr+1,cell2)-zgrid(kr-1,cell2))*( pp(kr,cell2)-pp   (kr-1,cell2))
               pl = pp(k,cell1)+.5*(zgrid(k   ,cell2)+zgrid(k +1,cell2)-zgrid(k ,cell1)-zgrid(k +1,cell1))   &
                                  /(zgrid(kl+1,cell1)-zgrid(kl-1,cell1))*( pp(kl,cell1)-pp   (kl-1,cell1))
               tend_u(k,iEdge) =  - cqu(k,iEdge)*2./(zz(k,cell1)+zz(k,cell2))*(pr-pl)/dcEdge(iEdge)

            end do

         else
            k = 1

            dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                  &
                         *(pzm(k,cell2)*(pp(k+1,cell2)-pp(k,cell2))    &
                          +pzm(k,cell1)*(pp(k+1,cell1)-pp(k,cell1))    &
                          +pzp(k,cell2)*(pp(k+2,cell2)-pp(k,cell2))    &
                          +pzp(k,cell1)*(pp(k+2,cell1)-pp(k,cell1))) 
  
            do k = 2, nVertLevels-1

               dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                  &
                          *(pzp(k,cell2)*(pp(k+1,cell2)-pp(k  ,cell2))    &
                           +pzm(k,cell2)*(pp(k  ,cell2)-pp(k-1,cell2))    &
                           +pzp(k,cell1)*(pp(k+1,cell1)-pp(k  ,cell1))    &
                           +pzm(k,cell1)*(pp(k  ,cell1)-pp(k-1,cell1)))   

            end do

            k = nVertLevels
            dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                  &
                       *(pzm(k,cell2)*(pp(k  ,cell2)-pp(k-1,cell2))    &
                        +pzm(k,cell1)*(pp(k  ,cell1)-pp(k-1,cell1)))   

            do k=1,nVertLevels

               tend_u(k,iEdge) =  - cqu(k,iEdge)*((pp(k,cell2)-pp(k,cell1))/dcEdge(iEdge)   &
                                       - dpzx(k) ) / (.5*(zz(k,cell2)+zz(k,cell1)))
            end do

         end if

         ! vertical transport of u

         wduz(1) = 0.
         k = 2
         wduz(k) =  0.5*( rw(k,cell1)+rw(k,cell2) )*(fzm(k)*u(k,iEdge)+fzp(k)*u(k-1,iEdge))  
         do k=3,nVertLevels-1
            wduz(k) = flux3( u(k-2,iEdge),u(k-1,iEdge),u(k,iEdge),u(k+1,iEdge),0.5*(rw(k,cell1)+rw(k,cell2)), 1.0_RKIND )
         end do
         k = nVertLevels
         wduz(k) =  0.5*( rw(k,cell1)+rw(k,cell2) )*(fzm(k)*u(k,iEdge)+fzp(k)*u(k-1,iEdge))  

         wduz(nVertLevels+1) = 0.

         do k=1,nVertLevels
            tend_u(k,iEdge) = tend_u(k,iEdge) - rdzw(k)*(wduz(k+1)-wduz(k)) 
         end do

         ! Next, nonlinear Coriolis term (q) following Ringler et al JCP 2009

         q(:) = 0.0
         do j = 1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(j,iEdge)
            do k=1,nVertLevels
               workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe))
               q(k) = q(k) + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv * rho_edge(k,eoe)
            end do
         end do

         do k=1,nVertLevels

            ! horizontal ke gradient and vorticity terms in the vector invariant formulation
            ! of the horizontal momentum equation
            tend_u(k,iEdge) = tend_u(k,iEdge) + rho_edge(k,iEdge)* (q(k) - (ke(k,cell2) - ke(k,cell1))       &
                                                                 / dcEdge(iEdge))                            &
                                             - u(k,iEdge)*0.5*(divergence_ru(k,cell1)+divergence_ru(k,cell2)) 
            if (curvature) then

               ! curvature terms for the sphere

               tend_u(k,iEdge) = tend_u(k,iEdge) &
                                - 2.*omega*cos(angleEdge(iEdge))*cos(latEdge(iEdge))  &
                                  *rho_edge(k,iEdge)*.25*(w(k,cell1)+w(k+1,cell1)+w(k,cell2)+w(k+1,cell2))          & 
                                - u(k,iEdge)*.25*(w(k+1,cell1)+w(k,cell1)+w(k,cell2)+w(k+1,cell2))                  &
                                  *rho_edge(k,iEdge)/r_earth
            end if
         end do
      end do

      deallocate(divergence_ru)

      !
      !  horizontal mixing for u
      !  mixing terms are integrated using forward-Euler, so this tendency is only computed in the
      !  first Runge-Kutta substep and saved for use in later RK substeps 2 and 3.
      !

      if (rk_step == 1 .or. rk_diffusion) then

      tend_u_euler = 0.

      if (delsq_horiz_mixing) then

         if ((h_mom_eddy_visc2 > 0.0) .and. (config_horiz_mixing == "2d_fixed")) then
            do iEdge=1, nEdgesSolve
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               vertex1 = verticesOnEdge(1,iEdge)
               vertex2 = verticesOnEdge(2,iEdge)
   
               do k=1,nVertLevels
  
                  !
                  ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
                  !                    only valid for h_mom_eddy_visc2 == constant
                  !
                  ! Note that we impose a lower bound on the edge length used in the derivative of the vorticity;
                  ! this is done to avoid an overly stringent stability constraint for small edge lengths that can
                  ! occur on some variable-resolution meshes.
                  !
                  u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                                 -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / max(dvEdge(iEdge),0.25*dcEdge(iEdge))
                  u_diffusion = rho_edge(k,iEdge)*h_mom_eddy_visc2 * u_diffusion
                  u_diffusion = u_diffusion * meshScalingDel2(iEdge)
  
                  tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + u_diffusion
               end do
            end do

         else if ( config_horiz_mixing == "2d_smagorinsky") then

            do iEdge=1, nEdgesSolve
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               vertex1 = verticesOnEdge(1,iEdge)
               vertex2 = verticesOnEdge(2,iEdge)
 
               do k=1,nVertLevels
                  !
                  ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
                  !                    only valid for h_mom_eddy_visc2 == constant
                  !
                  u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                                 -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / max(dvEdge(iEdge),0.25*dcEdge(iEdge))
                  u_diffusion = rho_edge(k,iEdge)* 0.5*(kdiff(k,cell1)+kdiff(k,cell2)) * u_diffusion
                  u_diffusion = u_diffusion * meshScalingDel2(iEdge)
 
                  tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + u_diffusion
               end do
            end do
         end if

      end if ! delsq_horiz_mixing for u

      if ((h_mom_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_fixed") .or. &
          (h_mom_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_smagorinsky")) then

         ! del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
         ! First, storage to hold the result from the first del^2 computation.

         allocate(delsq_divergence(nVertLevels, nCells+1))
         allocate(delsq_u(nVertLevels, nEdges+1))
         allocate(delsq_circulation(nVertLevels, nVertices+1))
         allocate(delsq_vorticity(nVertLevels, nVertices+1))

         delsq_u(:,:) = 0.0

         do iEdge=1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,nVertLevels

               !
               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc4 == constant
               !
               u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / max(dvEdge(iEdge), 0.25*dcEdge(iEdge))

               delsq_u(k,iEdge) = delsq_u(k,iEdge) + u_diffusion
            end do
         end do

         delsq_circulation(:,:) = 0.0
         do iEdge=1,nEdges
            do k=1,nVertLevels
               delsq_circulation(k,verticesOnEdge(1,iEdge)) = delsq_circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * delsq_u(k,iEdge)
               delsq_circulation(k,verticesOnEdge(2,iEdge)) = delsq_circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * delsq_u(k,iEdge)
            end do
         end do
         do iVertex=1,nVertices
            r = 1.0 / areaTriangle(iVertex)
            do k=1,nVertLevels
               delsq_vorticity(k,iVertex) = delsq_circulation(k,iVertex) * r
            end do
         end do

         delsq_divergence(:,:) = 0.0
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,nVertLevels
               delsq_divergence(k,cell1) = delsq_divergence(k,cell1) + delsq_u(k,iEdge)*dvEdge(iEdge)
               delsq_divergence(k,cell2) = delsq_divergence(k,cell2) - delsq_u(k,iEdge)*dvEdge(iEdge)
            end do
         end do
         do iCell = 1,nCells
            r = 1.0 / areaCell(iCell)
            do k = 1,nVertLevels
               delsq_divergence(k,iCell) = delsq_divergence(k,iCell) * r
            end do
         end do

         do iEdge=1,nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,nVertLevels

               !
               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc4 == constant
               !
               ! Here, we scale the diffusion on the divergence part a factor of config_del4u_div_factor 
               !    relative to the rotational part.  The stability constraint on the divergence component is much less
               !    stringent than the rotational part, and this flexibility may be useful.
               !
               u_diffusion =  rho_edge(k,iEdge) * ( config_del4u_div_factor * ( delsq_divergence(k,cell2)  - delsq_divergence(k,cell1) ) / dcEdge(iEdge)  &
                           -( delsq_vorticity(k,vertex2) - delsq_vorticity(k,vertex1) ) / max(dvEdge(iEdge), 0.25*dcEdge(iEdge)) &
                                                  )

               u_diffusion = u_diffusion * meshScalingDel4(iEdge)
               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) - h_mom_eddy_visc4 * u_diffusion
            end do
         end do

         deallocate(delsq_divergence)
         deallocate(delsq_u)
         deallocate(delsq_circulation)
         deallocate(delsq_vorticity)

      end if

      !
      !  vertical mixing for u - 2nd order filter in physical (z) space
      !
      if ( v_mom_eddy_visc2 > 0.0 ) then

         if (config_mix_full) then

            do iEdge=1,nEdgesSolve

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               do k=2,nVertLevels-1

                  z1 = 0.5*(zgrid(k-1,cell1)+zgrid(k-1,cell2))
                  z2 = 0.5*(zgrid(k  ,cell1)+zgrid(k  ,cell2))
                  z3 = 0.5*(zgrid(k+1,cell1)+zgrid(k+1,cell2))
                  z4 = 0.5*(zgrid(k+2,cell1)+zgrid(k+2,cell2))

                  zm = 0.5*(z1+z2)
                  z0 = 0.5*(z2+z3)
                  zp = 0.5*(z3+z4)

                  tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  &
                                     (u(k+1,iEdge)-u(k  ,iEdge))/(zp-z0)                      &
                                    -(u(k  ,iEdge)-u(k-1,iEdge))/(z0-zm) )/(0.5*(zp-zm))
               end do
            end do

         else  ! idealized cases where we mix on the perturbation from the initial 1-D state

            do iEdge=1,nEdgesSolve

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               do k=1,nVertLevels
#ifdef ROTATED_GRID
                  u_mix(k) = u(k,iEdge) - u_init(k) * sin( angleEdge(iEdge) )
#else
                  u_mix(k) = u(k,iEdge) - u_init(k) * cos( angleEdge(iEdge) )
#endif
               end do

               do k=2,nVertLevels-1

                  z1 = 0.5*(zgrid(k-1,cell1)+zgrid(k-1,cell2))
                  z2 = 0.5*(zgrid(k  ,cell1)+zgrid(k  ,cell2))
                  z3 = 0.5*(zgrid(k+1,cell1)+zgrid(k+1,cell2))
                  z4 = 0.5*(zgrid(k+2,cell1)+zgrid(k+2,cell2))

                  zm = 0.5*(z1+z2)
                  z0 = 0.5*(z2+z3)
                  zp = 0.5*(z3+z4)

                  tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  &
                                     (u_mix(k+1)-u_mix(k  ))/(zp-z0)                      &
                                    -(u_mix(k  )-u_mix(k-1))/(z0-zm) )/(0.5*(zp-zm))
               end do
            end do

         end if

      end if

      end if ! (rk_step 1 test for computing mixing terms)

!  add in mixing for u

      do iEdge=1,nEdgesSolve
         do k=1,nVertLevels
            tend_u(k,iEdge) = tend_u(k,iEdge) + tend_u_euler(k,iEdge)
         end do
      end do

!----------- rhs for w

      tend_w(:,:) = 0.
      if(rk_step .eq. 1) then
         tend_w_pgf(:,:)  = 0.
         tend_w_buoy(:,:) = 0.
      endif

      !
      !  horizontal advection for w
      !

      if (config_w_adv_order == 2) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
               do k=2,nVertLevels
                  flux = dvEdge(iEdge) * (fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge) ) &
                                        *(w(k,cell1) + w(k,cell2))*0.5 
                  tend_w(k,cell1) = tend_w(k,cell1) - flux
                  tend_w(k,cell2) = tend_w(k,cell2) + flux
               end do
            end if
         end do

      else if (config_w_adv_order == 3) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               do k=2,nVertLevels
                  ru_edge_w(k) = fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge)
               end do

               flux_arr(:) = 0.

               ! flux_arr stores the value of w at the cell edge used in the horizontal transport

               do i=1,nAdvCellsForEdge(iEdge)
                  iCell = advCellsForEdge(i,iEdge)
                  do k=2,nVertLevels
                     scalar_weight = adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,ru_edge_w(k))*adv_coefs_3rd(i,iEdge)
                     flux_arr(k) = flux_arr(k) + scalar_weight* w(k,iCell)
                  end do
               end do

               do k=2,nVertLevels
                  tend_w(k,cell1) = tend_w(k,cell1) - ru_edge_w(k)*flux_arr(k)
                  tend_w(k,cell2) = tend_w(k,cell2) + ru_edge_w(k)*flux_arr(k)
               end do

            end if
         end do

      else  if (config_w_adv_order == 4) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               do k=2,nVertLevels

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * w(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * w(k,cell2)
                  do i=1, nEdgesOnCell(cell1)
                     if ( cellsOnCell(i,cell1) <= nCells) &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * w(k,cellsOnCell(i,cell1))
                  end do
                  do i=1, nEdgesOnCell(cell2)
                     if ( cellsOnCell(i,cell2) <= nCells) &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * w(k,cellsOnCell(i,cell2))
                  end do

                  flux = dvEdge(iEdge) * (fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge)) * (  &
                                          0.5*(w(k,cell1) + w(k,cell2))                   &
                                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. )

                  tend_w(k,cell1) = tend_w(k,cell1) - flux
                  tend_w(k,cell2) = tend_w(k,cell2) + flux
               end do

            end if

         end do
      end if

      if (curvature) then

         do iCell = 1, nCellsSolve
            do k=2,nVertLevels
               tend_w(k,iCell) = tend_w(k,iCell) + (rho_zz(k,iCell)*fzm(k)+rho_zz(k-1,iCell)*fzp(k))*          &
                                         ( (fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))**2.             &
                                          +(fzm(k)*vr_cell(k,iCell)+fzp(k)*vr_cell(k-1,iCell))**2. )/r_earth   &
                                   + 2.*omega*cos(latCell(iCell))                                              &
                                          *(fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))                 &
                                          *(rho_zz(k,iCell)*fzm(k)+rho_zz(k-1,iCell)*fzp(k))

            end do
         end do

      end if

      !
      !  horizontal mixing for w - we could combine this with advection directly (i.e. as a turbulent flux),
      !  but here we can also code in hyperdiffusion if we wish (2nd order at present)
      !

      if (rk_step == 1 .or. rk_diffusion) then

         tend_w_euler = 0.

         if (delsq_horiz_mixing) then

            if ((h_mom_eddy_visc2 > 0.0) .and. (config_horiz_mixing == "2d_fixed")) then

               do iEdge=1,nEdges
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)
                  if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

                     ! horizontal flux divergence of the gradient (i.e. del^2)
                     ! note, for w, even though we use theta_* local scratch variables
                     do k=2,nVertLevels
                        theta_turb_flux = h_mom_eddy_visc2*(w(k,cell2) - w(k,cell1))/dcEdge(iEdge)
                        theta_turb_flux = theta_turb_flux * meshScalingDel2(iEdge)
                        flux = 0.5*dvEdge (iEdge) * (rho_edge(k,iEdge)+rho_edge(k-1,iEdge)) * theta_turb_flux
                        tend_w_euler(k,cell1) = tend_w_euler(k,cell1) + flux/areaCell(cell1)
                        tend_w_euler(k,cell2) = tend_w_euler(k,cell2) - flux/areaCell(cell2)
                     end do

                  end if
               end do

            else if (config_horiz_mixing == "2d_smagorinsky") then

               do iEdge=1,nEdges
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)
                  if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
   
                     do k=2,nVertLevels
                        theta_turb_flux = 0.25*(kdiff(k,cell1)+kdiff(k,cell2)+kdiff(k-1,cell1)+kdiff(k-1,cell2))  &
                                              *(w(k,cell2) - w(k,cell1))/dcEdge(iEdge)
                        theta_turb_flux = theta_turb_flux * meshScalingDel2(iEdge)
                        flux = 0.5*dvEdge (iEdge) * (rho_edge(k,iEdge)+rho_edge(k-1,iEdge)) * theta_turb_flux
                        tend_w_euler(k,cell1) = tend_w_euler(k,cell1) + flux/areaCell(cell1)
                        tend_w_euler(k,cell2) = tend_w_euler(k,cell2) - flux/areaCell(cell2)
                     end do

                  end if
               end do
            end if
         end if   ! delsq_horiz_mixing

         if ((h_mom_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_fixed") .or. &
             (h_mom_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_smagorinsky")) then


            ! del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
            !
            ! First, storage to hold the result from the first del^2 computation.
            !  we copied code from the theta mixing, hence the theta* names.

            allocate(delsq_theta(nVertLevels, nCells+1))

            delsq_theta(:,:) = 0.

            do iEdge=1,nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               do k=2,nVertLevels
                  delsq_theta(k,cell1) = delsq_theta(k,cell1) + dvEdge(iEdge)*0.5*(rho_edge(k,iEdge)+rho_edge(k-1,iEdge))*(w(k,cell2) - w(k,cell1))/dcEdge(iEdge)
                  delsq_theta(k,cell2) = delsq_theta(k,cell2) - dvEdge(iEdge)*0.5*(rho_edge(k,iEdge)+rho_edge(k-1,iEdge))*(w(k,cell2) - w(k,cell1))/dcEdge(iEdge)
               end do
            end do

            do iCell = 1, nCells
               r = 1.0 / areaCell(iCell)
               do k=2,nVertLevels
                  delsq_theta(k,iCell) = delsq_theta(k,iCell) * r
               end do
            end do

            do iEdge=1,nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

                  do k=2,nVertLevels
                     theta_turb_flux = h_mom_eddy_visc4*(delsq_theta(k,cell2) - delsq_theta(k,cell1))/dcEdge(iEdge)
                     theta_turb_flux = theta_turb_flux * meshScalingDel4(iEdge)
                     flux = dvEdge (iEdge) * theta_turb_flux
                     tend_w_euler(k,cell1) = tend_w_euler(k,cell1) - flux/areaCell(cell1)
                     tend_w_euler(k,cell2) = tend_w_euler(k,cell2) + flux/areaCell(cell2)
                  end do

               end if
            end do

            deallocate(delsq_theta)

         end if

      end if ! horizontal mixing for w computed in first rk_step

      !
      !  vertical advection, pressure gradient and buoyancy for w
      !

      do iCell = 1, nCells

         wdwz(1) = 0.
         if (config_w_vadv_order == 2) then

            do k=2,nVertLevels
               wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))
            end do

         else if (config_w_vadv_order == 3) then

            k = 2
            wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))
            do k=3,nVertLevels-1
               wdwz(k) = flux3( w(k-2,iCell),w(k-1,iCell),w(k,iCell),w(k+1,iCell),0.5*(rw(k,iCell)+rw(k-1,iCell)), 1.0_RKIND )
            end do
            k = nVertLevels
            wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))

         else if (config_w_vadv_order == 4) then

            k = 2
            wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))
            do k=3,nVertLevels-1
               wdwz(k) = flux4( w(k-2,iCell),w(k-1,iCell),w(k,iCell),w(k+1,iCell),0.5*(rw(k,iCell)+rw(k-1,iCell)) )
            end do
            k = nVertLevels
            wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))

         end if

         wdwz(nVertLevels+1) = 0.

      !  Note: next we are also dividing through by the cell area after the horizontal flux divergence

         do k=2,nVertLevels

            tend_w(k,iCell) = tend_w(k,iCell)/areaCell(iCell) -rdzu(k)*(wdwz(k+1)-wdwz(k))    &
                                  - cqw(k,iCell)*( rdzu(k)*(pp(k,iCell)-pp(k-1,iCell))        &
                                  + gravity*  &
                                   ( fzm(k)*(rb(k,iCell)*(qtot(k,iCell)) +         &
                                             rr(k,iCell)*(1.+qtot(k,iCell)))                  &
                                    +fzp(k)*(rb(k-1,iCell)*(qtot(k-1,iCell))  +  &
                                             rr(k-1,iCell)*(1.+qtot(k-1,iCell))) ))

            if(rk_step == 1) then
               tend_w_pgf(k,iCell)  = cqw(k,iCell)*(rdzu(k)*(pp(k,iCell)-pp(k-1,iCell))) 
               tend_w_buoy(k,iCell) = cqw(k,iCell)*gravity* &
                                   ( fzm(k)*(rb(k,iCell)*(qtot(k,iCell)) +       &
                                             rr(k,iCell)*(1.+qtot(k,iCell)))     &
                                    +fzp(k)*(rb(k-1,iCell)*(qtot(k-1,iCell))  +  &
                                             rr(k-1,iCell)*(1.+qtot(k-1,iCell))) )
            endif

         end do
      end do

      !
      !  vertical mixing for w - 2nd order 
      !

      if (rk_step == 1 .or. rk_diffusion) then

         if ( v_mom_eddy_visc2 > 0.0 ) then

            do iCell = 1, nCellsSolve
            do k=2,nVertLevels
               tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + v_mom_eddy_visc2*0.5*(rho_zz(k,iCell)+rho_zz(k-1,iCell))*(  &
                                        (w(k+1,iCell)-w(k  ,iCell))*rdzw(k)                              &
                                       -(w(k  ,iCell)-w(k-1,iCell))*rdzw(k-1) )*rdzu(k)
            end do
            end do

         end if

      end if ! mixing term computed first rk_step


! add in mixing terms for w

      do iCell = 1, nCellsSolve
         do k=2,nVertLevels
            tend_w(k,iCell) = tend_w(k,iCell) + tend_w_euler(k,iCell)
         end do
      end do

      deallocate(qtot)

!----------- rhs for theta

      tend_theta(:,:) = 0.

      !
      !  horizontal advection for theta
      !

      if (config_theta_adv_order == 2) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
               do k=1,nVertLevels
                  flux = dvEdge(iEdge) *  ru(k,iEdge) * ( 0.5*(theta_m(k,cell1) + theta_m(k,cell2)) )
                  tend_theta(k,cell1) = tend_theta(k,cell1) - flux
                  tend_theta(k,cell2) = tend_theta(k,cell2) + flux
               end do
            end if
         end do

      else if (config_theta_adv_order == 3) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               flux_arr(:) = 0.
               do i=1,nAdvCellsForEdge(iEdge)
                  iCell = advCellsForEdge(i,iEdge)
                  do k=1,nVertLevels
                     scalar_weight = adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,ru(k,iEdge))*adv_coefs_3rd(i,iEdge)
                     flux_arr(k) = flux_arr(k) + scalar_weight* theta_m(k,iCell)
                  end do
               end do

               do k=1,nVertLevels
                  tend_theta(k,cell1) = tend_theta(k,cell1) - ru(k,iEdge)*flux_arr(k)
                  tend_theta(k,cell2) = tend_theta(k,cell2) + ru(k,iEdge)*flux_arr(k)
               end do

            end if
         end do

      else  if (config_theta_adv_order == 4) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCells .and. cell2 <= nCells) then

               do k=1,nVertLevels

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * theta_m(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * theta_m(k,cell2)
                  do i=1, nEdgesOnCell(cell1)
                     if ( cellsOnCell(i,cell1) <= nCells) &
                        d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * theta_m(k,cellsOnCell(i,cell1))
                  end do
                  do i=1, nEdgesOnCell(cell2)
                     if ( cellsOnCell(i,cell2) <= nCells) &
                        d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * theta_m(k,cellsOnCell(i,cell2))
                  end do

                  flux = dvEdge(iEdge) *  ru(k,iEdge) * (                                               &
                                         0.5*(theta_m(k,cell1) + theta_m(k,cell2))                      &
                                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. )

                  tend_theta(k,cell1) = tend_theta(k,cell1) - flux
                  tend_theta(k,cell2) = tend_theta(k,cell2) + flux
               end do

            end if

         end do
      end if

      !
      !  horizontal mixing for theta_m - we could combine this with advection directly (i.e. as a turbulent flux),
      !  but here we can also code in hyperdiffusion if we wish (2nd order at present)
      !

      if (rk_step == 1 .or. rk_diffusion) then

      tend_theta_euler = 0.

      if (delsq_horiz_mixing) then
         if ( (h_theta_eddy_visc2 > 0.0) .and. (config_horiz_mixing == "2d_fixed") ) then

            do iEdge=1,nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
  
                  do k=1,nVertLevels
                     theta_turb_flux = h_theta_eddy_visc2*prandtl_inv*(theta_m(k,cell2) - theta_m(k,cell1))/dcEdge(iEdge)
                     theta_turb_flux = theta_turb_flux * meshScalingDel2(iEdge)
                     flux = dvEdge (iEdge) * rho_edge(k,iEdge) * theta_turb_flux
                     tend_theta_euler(k,cell1) = tend_theta_euler(k,cell1) + flux/areaCell(cell1)
                     tend_theta_euler(k,cell2) = tend_theta_euler(k,cell2) - flux/areaCell(cell2)
                  end do
  
               end if
            end do

         else if (  ( config_horiz_mixing == "2d_smagorinsky") ) then

            do iEdge=1,nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
 
                  do k=1,nVertLevels
                     theta_turb_flux = 0.5*(kdiff(k,cell1)+kdiff(k,cell2))*prandtl_inv  &
                                           *(theta_m(k,cell2) - theta_m(k,cell1))/dcEdge(iEdge)
                     theta_turb_flux = theta_turb_flux * meshScalingDel2(iEdge)
                     flux = dvEdge (iEdge) * rho_edge(k,iEdge) * theta_turb_flux
                     tend_theta_euler(k,cell1) = tend_theta_euler(k,cell1) + flux/areaCell(cell1)
                     tend_theta_euler(k,cell2) = tend_theta_euler(k,cell2) - flux/areaCell(cell2)
                  end do
  
               end if
            end do
         end if

      end if

      if ((h_theta_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_fixed") .or. &
          (h_theta_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_smagorinsky")) then

         allocate(delsq_theta(nVertLevels, nCells+1))

         delsq_theta(:,:) = 0.

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,nVertLevels
               delsq_theta(k,cell1) = delsq_theta(k,cell1) + dvEdge(iEdge)*rho_edge(k,iEdge)*(theta_m(k,cell2) - theta_m(k,cell1))/dcEdge(iEdge)
               delsq_theta(k,cell2) = delsq_theta(k,cell2) - dvEdge(iEdge)*rho_edge(k,iEdge)*(theta_m(k,cell2) - theta_m(k,cell1))/dcEdge(iEdge)
            end do
         end do

         do iCell = 1, nCells
            r = 1.0 / areaCell(iCell)
            do k=1,nVertLevels
               delsq_theta(k,iCell) = delsq_theta(k,iCell) * r
            end do
         end do

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               do k=1,nVertLevels
                  theta_turb_flux = h_theta_eddy_visc4*prandtl_inv*(delsq_theta(k,cell2) - delsq_theta(k,cell1))/dcEdge(iEdge)
                  theta_turb_flux = theta_turb_flux * meshScalingDel4(iEdge)
                  flux = dvEdge (iEdge) * theta_turb_flux
                  tend_theta_euler(k,cell1) = tend_theta_euler(k,cell1) - flux/areaCell(cell1)
                  tend_theta_euler(k,cell2) = tend_theta_euler(k,cell2) + flux/areaCell(cell2)
               end do

            end if
         end do

         deallocate(delsq_theta)

      end if

      end if ! theta mixing calculated first rk_step

      !
      !  vertical advection plus diabatic term
      !  Note: we are also dividing through by the cell area after the horizontal flux divergence
      !
      do iCell = 1, nCells

         wdtz(1) = 0.
         if (config_theta_vadv_order == 2) then

            do k=2,nVertLevels
               wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))
            end do

         else if (config_theta_vadv_order == 3) then

            k = 2
            wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))
            do k=3,nVertLevels-1
               wdtz(k) = flux3( theta_m(k-2,iCell),theta_m(k-1,iCell),theta_m(k,iCell),theta_m(k+1,iCell), rw(k,iCell), coef_3rd_order )
            end do
            k = nVertLevels
            wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))

         else if (config_theta_vadv_order == 4) then

            k = 2
            wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))
            do k=3,nVertLevels-1
               wdtz(k) = flux4( theta_m(k-2,iCell),theta_m(k-1,iCell),theta_m(k,iCell),theta_m(k+1,iCell), rw(k,iCell) )
            end do
            k = nVertLevels
            wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))

         end if


         wdtz(nVertLevels+1) = 0.

         do k=1,nVertLevels
            tend_theta(k,iCell) = tend_theta(k,iCell)/areaCell(iCell) -rdzw(k)*(wdtz(k+1)-wdtz(k))
            tend_theta(k,iCell) = tend_theta(k,iCell) + rho_zz(k,iCell)*rt_diabatic_tend(k,iCell)
         end do
      end do

      !
      !  vertical mixing for theta - 2nd order 
      !

      if (rk_step == 1 .or. rk_diffusion) then

      if ( v_theta_eddy_visc2 > 0.0 ) then

         if (config_mix_full) then

            do iCell = 1, nCellsSolve
               do k=2,nVertLevels-1
                  z1 = zgrid(k-1,iCell)
                  z2 = zgrid(k  ,iCell)
                  z3 = zgrid(k+1,iCell)
                  z4 = zgrid(k+2,iCell)

                  zm = 0.5*(z1+z2)
                  z0 = 0.5*(z2+z3)
                  zp = 0.5*(z3+z4)

                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
                                           (theta_m(k+1,iCell)-theta_m(k  ,iCell))/(zp-z0)                 &
                                          -(theta_m(k  ,iCell)-theta_m(k-1,iCell))/(z0-zm) )/(0.5*(zp-zm))
               end do
            end do

         else  ! idealized cases where we mix on the perturbation from the initial 1-D state

            do iCell = 1, nCellsSolve
               do k=2,nVertLevels-1
                  z1 = zgrid(k-1,iCell)
                  z2 = zgrid(k  ,iCell)
                  z3 = zgrid(k+1,iCell)
                  z4 = zgrid(k+2,iCell)

                  zm = 0.5*(z1+z2)
                  z0 = 0.5*(z2+z3)
                  zp = 0.5*(z3+z4)

                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
                                           ((theta_m(k+1,iCell)-t_init(k+1,iCell))-(theta_m(k  ,iCell)-t_init(k,iCell)))/(zp-z0)      &
                                          -((theta_m(k  ,iCell)-t_init(k,iCell))-(theta_m(k-1,iCell)-t_init(k-1,iCell)))/(z0-zm) )/(0.5*(zp-zm))
               end do
            end do

         end if

      end if

      end if ! compute theta mixing on first rk_step

      do iCell = 1, nCellsSolve
         do k=1,nVertLevels
            tend_theta(k,iCell) = tend_theta(k,iCell) + tend_theta_euler(k,iCell)
         end do
      end do

   end subroutine atm_compute_dyn_tend

!-------

   subroutine atm_compute_solve_diagnostics(dt, state, time_lev, diag, mesh, configs)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute diagnostic fields used in the tendency computations
   !
   ! Input: state (s), grid - grid metadata
   !
   ! Output: diag - computed diagnostics
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      real (kind=RKIND), intent(in) :: dt
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev                   ! which time level of state to use
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: configs


      integer :: iEdge, iCell, iVertex, k, cell1, cell2, eoe, i
      real (kind=RKIND) :: h_vertex, r

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree
      real (kind=RKIND), dimension(:), pointer :: fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, &
                                                    circulation, vorticity, ke, pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, &
                                                    divergence
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge

      logical, parameter :: hollingsworth=.true.
      real (kind=RKIND), allocatable, dimension(:,:) :: ke_vertex
      real (kind=RKIND) :: ke_fact
      real (kind=RKIND), pointer :: config_apvm_upwinding


      call mpas_pool_get_config(configs, 'config_apvm_upwinding', config_apvm_upwinding)

      call mpas_pool_get_array(state, 'rho_zz', h, time_lev)
      call mpas_pool_get_array(state, 'u', u, time_lev)

      call mpas_pool_get_array(diag, 'v', v)
      call mpas_pool_get_array(diag, 'rv', vh)
      call mpas_pool_get_array(diag, 'rho_edge', h_edge)
      call mpas_pool_get_array(diag, 'circulation', circulation)
      call mpas_pool_get_array(diag, 'vorticity', vorticity)
      call mpas_pool_get_array(diag, 'divergence', divergence)
      call mpas_pool_get_array(diag, 'ke', ke)
      call mpas_pool_get_array(diag, 'pv_edge', pv_edge)
      call mpas_pool_get_array(diag, 'pv_vertex', pv_vertex)
      call mpas_pool_get_array(diag, 'pv_cell', pv_cell)
      call mpas_pool_get_array(diag, 'gradPVn', gradPVn)
      call mpas_pool_get_array(diag, 'gradPVt', gradPVt)

      call mpas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'fVertex', fVertex)
      call mpas_pool_get_array(mesh, 'fEdge', fEdge)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)


      !
      ! Compute height on cell edges at velocity locations
      !
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
            h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
         end do
      end do


      !
      ! Compute circulation and relative vorticity at each vertex
      !
      circulation(:,:) = 0.0
      do iEdge=1,nEdges
         do k=1,nVertLevels
            circulation(k,verticesOnEdge(1,iEdge)) = circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u(k,iEdge)
            circulation(k,verticesOnEdge(2,iEdge)) = circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u(k,iEdge)
         end do
      end do
      do iVertex=1,nVertices
         do k=1,nVertLevels
            vorticity(k,iVertex) = circulation(k,iVertex) / areaTriangle(iVertex)
         end do
      end do


      !
      ! Compute the divergence at each cell center
      !
      divergence(:,:) = 0.0
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
           divergence(k,cell1) = divergence(k,cell1) + u(k,iEdge)*dvEdge(iEdge)
           divergence(k,cell2) = divergence(k,cell2) - u(k,iEdge)*dvEdge(iEdge)
         end do
      end do
      do iCell = 1,nCells
         r = 1.0 / areaCell(iCell)
         do k = 1,nVertLevels
            divergence(k,iCell) = divergence(k,iCell) * r
         end do
      end do


      !
      ! Compute kinetic energy in each cell (Ringler et al JCP 2009)
      !
      ke(:,:) = 0.0
      do iCell=1,nCells
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            do k=1,nVertLevels
               ke(k,iCell) = ke(k,iCell) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2.0
            end do
         end do
         do k=1,nVertLevels
            ke(k,iCell) = ke(k,iCell) / areaCell(iCell)
         end do
      end do

      if (hollingsworth) then

         ! Compute ke at cell vertices - AG's new KE construction, part 1
         ! *** approximation here because we don't have inner triangle areas
         !

         allocate (ke_vertex(nVertLevels,nVertices))
         do iVertex=1,nVertices
            do k=1,nVertLevels

               ke_vertex(k,iVertex) = (  dcEdge(EdgesOnVertex(1,iVertex))*dvEdge(EdgesOnVertex(1,iVertex))*u(k,EdgesOnVertex(1,iVertex))**2.0  &
                                        +dcEdge(EdgesOnVertex(2,iVertex))*dvEdge(EdgesOnVertex(2,iVertex))*u(k,EdgesOnVertex(2,iVertex))**2.0  &
                                        +dcEdge(EdgesOnVertex(3,iVertex))*dvEdge(EdgesOnVertex(3,iVertex))*u(k,EdgesOnVertex(3,iVertex))**2.0  &
                                                    ) * 0.25 / AreaTriangle(iVertex)

            end do
         end do

         ! adjust ke at cell vertices - AG's new KE construction, part 2
         !

         ke_fact = 1.0 - .375

         do iCell=1,nCells
            do k=1,nVertLevels
               ke(k,iCell) = ke_fact*ke(k,iCell)
            end do
         end do

         ! Avoid FP errors caused by a potential division by zero below by 
         ! initializing the "garbage cell" of areaCell to a non-zero value
         areaCell(nCells+1) = 1.0

         do iVertex = 1, nVertices
            do i=1,vertexDegree
               iCell = cellsOnVertex(i,iVertex)
               do k = 1,nVertLevels
                  ke(k,iCell) = ke(k,iCell) + (1.-ke_fact)*kiteAreasOnVertex(i, iVertex) * ke_vertex(k, iVertex) / areaCell(iCell)
               end do
            end do
         end do
         deallocate (ke_vertex)

      end if

      !
      ! Compute v (tangential) velocities following Thuburn et al JCP 2009
      !
      v(:,:) = 0.0
      do iEdge = 1,nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1,nVertLevels
               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
            end do
         end do
      end do


      !
      ! Compute height at vertices, pv at vertices, and average pv to edge locations
      !  ( this computes pv_vertex at all vertices bounding real cells )
      !
      do iVertex = 1,nVertices
         do k=1,nVertLevels
            h_vertex = 0.0
            do i=1,vertexDegree
               h_vertex = h_vertex + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            h_vertex = h_vertex / areaTriangle(iVertex)

            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex
         end do
      end do


      !
      ! Compute pv at the edges
      !   ( this computes pv_edge at all edges bounding real cells )
      !
      pv_edge(:,:) = 0.0
      do iVertex = 1,nVertices
         do i=1,vertexDegree
            iEdge = edgesOnVertex(i,iVertex)
            do k=1,nVertLevels
               pv_edge(k,iEdge) =  pv_edge(k,iEdge)  + 0.5 * pv_vertex(k,iVertex)
            end do
         end do
      end do

      !
      ! Compute pv at cell centers
      !    ( this computes pv_cell for all real cells )
      !
      pv_cell(:,:) = 0.0
      do iVertex = 1, nVertices
         do i=1,vertexDegree
            iCell = cellsOnVertex(i,iVertex)
            do k = 1,nVertLevels
               pv_cell(k,iCell) = pv_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * pv_vertex(k, iVertex) / areaCell(iCell)
            end do
         end do
      end do


      if (config_apvm_upwinding > 0.0) then

         !
         ! Modify PV edge with upstream bias. 
         !
         ! Compute gradient of PV in the tangent direction
         !   ( this computes gradPVt at all edges bounding real cells )
         !
         do iEdge = 1,nEdges
            do k = 1,nVertLevels
               gradPVt(k,iEdge) = (pv_vertex(k,verticesOnEdge(2,iEdge)) - pv_vertex(k,verticesOnEdge(1,iEdge))) / &
                                  dvEdge(iEdge)
            end do
         end do

         !
         ! Compute gradient of PV in normal direction
         !   (tdr: 2009-10-02: this is not correct because the pv_cell in the halo is not correct)
         !
         gradPVn(:,:) = 0.0
         do iEdge = 1,nEdges
            do k = 1,nVertLevels
               gradPVn(k,iEdge) = (pv_cell(k,cellsOnEdge(2,iEdge)) - pv_cell(k,cellsOnEdge(1,iEdge))) / &
                                    dcEdge(iEdge)
            end do
         end do

         do iEdge = 1,nEdges
            do k = 1,nVertLevels
               pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * v(k,iEdge) * dt * gradPVt(k,iEdge)
            end do
         end do

         ! Modify PV edge with upstream bias.
         !
         do iEdge = 1,nEdges
            do k = 1,nVertLevels
               pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * u(k,iEdge) *dt * gradPVn(k,iEdge)
            end do
         end do

      end if  ! apvm upwinding


   end subroutine atm_compute_solve_diagnostics

!----------

   subroutine atm_init_coupled_diagnostics( state, time_lev, diag, mesh, configs )

      implicit none
   
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev                    ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in) :: configs

      integer :: k, iCell, iEdge, iCell1, iCell2, cell1, cell2
      real (kind=RKIND), pointer :: coef_3rd_order
      integer, pointer :: config_theta_adv_order
      integer, pointer :: nCells, nEdges, nVertLevels
      integer, pointer :: index_qv
      real (kind=RKIND) :: p0, rcv, flux
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:,:), pointer :: theta_m
      real (kind=RKIND), dimension(:,:), pointer :: theta
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz
      real (kind=RKIND), dimension(:,:), pointer :: rho
      real (kind=RKIND), dimension(:,:), pointer :: rho_p
      real (kind=RKIND), dimension(:,:), pointer :: rho_base
      real (kind=RKIND), dimension(:,:), pointer :: rtheta_base
      real (kind=RKIND), dimension(:,:), pointer :: theta_base
      real (kind=RKIND), dimension(:,:), pointer :: rtheta_p
      real (kind=RKIND), dimension(:,:), pointer :: zz
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: ru
      real (kind=RKIND), dimension(:,:), pointer :: rw
      real (kind=RKIND), dimension(:,:), pointer :: u
      real (kind=RKIND), dimension(:,:), pointer :: w
      real (kind=RKIND), dimension(:,:), pointer :: pressure_p
      real (kind=RKIND), dimension(:,:), pointer :: exner
      real (kind=RKIND), dimension(:,:), pointer :: exner_base
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3 


      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)

      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)

      call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, time_lev)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'rho_p', rho_p)
      call mpas_pool_get_array(diag, 'rho_base', rho_base)
      call mpas_pool_get_array(diag, 'rtheta_base', rtheta_base)
      call mpas_pool_get_array(diag, 'theta_base', theta_base)
      call mpas_pool_get_array(diag, 'rtheta_p', rtheta_p)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(state, 'u', u, time_lev)
      call mpas_pool_get_array(state, 'w', w, time_lev)
      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'exner_base', exner_base)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)

      if (config_theta_adv_order /= 3) coef_3rd_order = 0.0

      rcv = rgas / (cp-rgas)
      p0 = 1.e5  ! this should come from somewhere else...

      do iCell=1,nCells
         do k=1,nVertLevels
            theta_m(k,iCell) = theta(k,iCell) * (1._RKIND + rvord * scalars(index_qv,k,iCell))
            rho_zz(k,iCell) = rho(k,iCell) / zz(k,iCell)
         end do
      end do

      do iEdge = 1, nEdges
         iCell1 = cellsOnEdge(1,iEdge)
         iCell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
            ru(k,iEdge) = 0.5 * u(k,iEdge) * (rho_zz(k,iCell1) + rho_zz(k,iCell2))
         end do
      end do

      ! Compute rw (i.e. rho_zz * omega) from rho_zz, w, and ru.
      ! We are reversing the procedure we use in subroutine atm_recover_large_step_variables.
      ! first, the piece that depends on w.
      do iCell=1,nCells
         rw(1,iCell) = 0.0
         rw(nVertLevels+1,iCell) = 0.0
         do k=2,nVertLevels
            rw(k,iCell) = w(k,iCell)     &
                          * (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell)) &
                          * (fzp(k) * zz(k-1,iCell) + fzm(k) * zz(k,iCell))
         end do
      end do
  
      ! next, the piece that depends on ru
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 2, nVertLevels
            flux = (fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge))
            rw(k,cell2) = rw(k,cell2)   &
                          + (zb(k,2,iEdge) + coef_3rd_order * sign(1.0_RKIND,flux) * zb3(k,2,iEdge))*flux   &
                          * (fzp(k) * zz(k-1,cell2) + fzm(k) * zz(k,cell2))
            rw(k,cell1) = rw(k,cell1)   &
                          - (zb(k,1,iEdge) + coef_3rd_order * sign(1.0_RKIND,flux) * zb3(k,1,iEdge))*flux   &
                          * (fzp(k) * zz(k-1,cell1) + fzm(k) * zz(k,cell1))
         end do
      end do

      do iCell = 1, nCells
         do k=1,nVertLevels
            rho_p(k,iCell) = rho_zz(k,iCell) - rho_base(k,iCell)
         end do
      end do

      do iCell = 1, nCells
         do k=1,nVertLevels
            rtheta_base(k,iCell) = theta_base(k,iCell) * rho_base(k,iCell)
         end do
      end do

      do iCell = 1, nCells
         do k=1,nVertLevels
            rtheta_p(k,iCell) = theta_m(k,iCell) * rho_p(k,iCell)  &
                                             + rho_base(k,iCell)   * (theta_m(k,iCell) - theta_base(k,iCell))
         end do
      end do

      do iCell=1,nCells
         do k=1,nVertLevels
            exner(k,iCell) = (zz(k,iCell) * (rgas/p0) * (rtheta_p(k,iCell) + rtheta_base(k,iCell)))**rcv
         end do
      end do

      do iCell=1,nCells
         do k=1,nVertLevels
            pressure_p(k,iCell) = zz(k,iCell) * rgas &
                                               * (  exner(k,iCell) * rtheta_p(k,iCell) &
                                                  + rtheta_base(k,iCell) * (exner(k,iCell) - exner_base(k,iCell)) &
                                                 )
         end do
      end do

   end subroutine atm_init_coupled_diagnostics


end module atm_time_integration
