! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

#ifdef _GPTL
#define MPAS_GPTL_START(T) call gptlstart(T)
#define MPAS_GPTL_STOP(T) call gptlstop(T)
#define MPAS_GPTL_PRINT(T) call gptlpr(T)
#else
#define MPAS_GPTL_START(T) !!! call gptlstart(T)
#define MPAS_GPTL_STOP(T) !!! call gptlstop(T)
#define MPAS_GPTL_PRINT(T) !!! call gptlpr(T)
#endif

module atm_time_integration

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_kind_types
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   ! Added only clause to keep xlf90 from getting confused from the overloaded abs intrinsic in mpas_timekeeping
   use mpas_timekeeping, only: MPAS_Time_type, MPAS_TimeInterval_type, &
                               mpas_set_time, mpas_set_timeInterval, mpas_get_time, operator(+), add_t_ti
   use mpas_timer

#ifdef DO_PHYSICS
   use mpas_atmphys_driver_microphysics
   use mpas_atmphys_todynamics
   use mpas_atmphys_utilities
#endif

   use mpas_atm_iau  

   integer :: timerid, secs, u_secs

   ! Used to store physics tendencies for dynamics variables
   real (kind=RKIND), allocatable, dimension(:,:) :: tend_ru_physics, tend_rtheta_physics
   
   ! Used in compute_dyn_tend
   real (kind=RKIND), allocatable, dimension(:,:) :: qtot 
   real (kind=RKIND), allocatable, dimension(:,:) :: delsq_theta, delsq_w, delsq_divergence
   real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
!   real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation    ! no longer used -> removed 
   real (kind=RKIND), allocatable, dimension(:,:) :: delsq_vorticity
   real (kind=RKIND), allocatable, dimension(:,:) :: dpdz

   ! Used in atm_advance_scalars
   real (kind=RKIND), dimension(:,:,:), allocatable :: horiz_flux_array

   ! Used in atm_advance_scalars_mono
   real (kind=RKIND), dimension(:,:), allocatable :: scalar_old_arr, scalar_new_arr
   real (kind=RKIND), dimension(:,:), allocatable :: s_max_arr, s_min_arr
   real (kind=RKIND), dimension(:,:,:), allocatable :: scale_array
   real (kind=RKIND), dimension(:,:), allocatable :: flux_array
   real (kind=RKIND), dimension(:,:), allocatable :: flux_upwind_tmp_arr
   real (kind=RKIND), dimension(:,:), allocatable :: flux_tmp_arr
   real (kind=RKIND), dimension(:,:), allocatable :: wdtn_arr
   real (kind=RKIND), dimension(:,:), allocatable :: rho_zz_int
   real (kind=RKIND), dimension(:,:,:), allocatable :: scalar_tend_array

   ! Used in compute_solve_diagnostics
   real (kind=RKIND), allocatable, dimension(:,:) :: ke_vertex
   real (kind=RKIND), allocatable, dimension(:,:) :: ke_edge

      
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
         MPAS_GPTL_START('atm_srk3')
         call atm_srk3(domain, dt, itimestep)
         MPAS_GPTL_STOP('atm_srk3')
      else
         call mpas_dmpar_global_abort('Unknown time integration option '//trim(config_time_integration), deferredAbort=.true.)
         call mpas_dmpar_global_abort('Currently, only ''SRK3'' is supported.')
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

      integer :: thread
      integer :: iCell, k, iEdge
      type (block_type), pointer :: block

      integer, pointer :: nThreads
      integer, dimension(:), pointer :: cellThreadStart, cellThreadEnd
      integer, dimension(:), pointer :: cellSolveThreadStart, cellSolveThreadEnd
      integer, dimension(:), pointer :: edgeThreadStart, edgeThreadEnd
      integer, dimension(:), pointer :: edgeSolveThreadStart, edgeSolveThreadEnd
      integer, dimension(:), pointer :: vertexThreadStart, vertexThreadEnd
      integer, dimension(:), pointer :: vertexSolveThreadStart, vertexSolveThreadEnd

      integer :: rk_step, number_of_sub_steps
      integer :: iScalar

      real (kind=RKIND), dimension(3) :: rk_timestep, rk_sub_timestep
      integer, dimension(3) :: number_sub_steps
      integer :: small_step

      !  additions for splitting scalar transport from dynamics, WCS 18 November 2014
      logical, pointer :: config_split_dynamics_transport
      integer, pointer :: config_dynamics_split
      integer :: dynamics_substep, dynamics_split
      real (kind=RKIND) :: dt_dynamics

      integer, pointer :: config_number_of_sub_steps
      integer, pointer :: config_time_integration_order
      logical, pointer :: config_scalar_advection
      logical, pointer :: config_positive_definite
      logical, pointer :: config_monotonic
      real (kind=RKIND), pointer :: config_dt
      character (len=StrKIND), pointer :: config_microp_scheme
      character (len=StrKIND), pointer :: config_convection_scheme

      integer, pointer :: num_scalars, index_qv, nCells, nCellsSolve, nEdges, nEdgesSolve, nVertices, nVerticesSolve, nVertLevels

      integer, pointer :: config_IAU_option
      real(kind=RKIND) :: weight_iau
      real(kind=RKIND),      pointer :: time_window_sec

      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: diag_physics
      type (mpas_pool_type), pointer :: mesh
      type (mpas_pool_type), pointer :: tend
      type (mpas_pool_type), pointer :: tend_physics
      type (mpas_pool_type), pointer :: tend_iau

      type (field2DReal), pointer :: theta_m_field
      type (field3DReal), pointer :: scalars_field
      type (field2DReal), pointer :: pressure_p_field
      type (field2DReal), pointer :: rtheta_p_field
      type (field2DReal), pointer :: rtheta_pp_field
      type (field2DReal), pointer :: divergence_3d_field
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

      logical, parameter :: debug = .false.


      !
      ! Retrieve configuration options
      !
      call mpas_pool_get_config(domain % blocklist % configs, 'config_number_of_sub_steps', config_number_of_sub_steps)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_time_integration_order', config_time_integration_order)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_scalar_advection', config_scalar_advection)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_positive_definite', config_positive_definite)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_monotonic', config_monotonic)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', config_dt)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_microp_scheme', config_microp_scheme)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_convection_scheme', config_convection_scheme)

      !  config variables for dynamics-transport splitting, WCS 18 November 2014
      call mpas_pool_get_config(domain % blocklist % configs, 'config_split_dynamics_transport', config_split_dynamics_transport)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dynamics_split_steps', config_dynamics_split)

      ! Retrieve weight/coefficient for IAU (Incremental Analysis Update) - HA (June-16-2016)
      weight_iau = atm_iau_coef (domain % blocklist % configs, itimestep, dt)

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
      ! allocate storage for physics tendency save
      !
      call mpas_pool_get_dimension(state, 'nCells', nCells)
      call mpas_pool_get_dimension(state, 'nEdges', nEdges)
      call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels)

      allocate(qtot(nVertLevels,nCells+1))
      allocate(tend_rtheta_physics(nVertLevels,nCells+1))
      allocate(tend_ru_physics(nVertLevels,nEdges+1))

      !
      ! Initialize RK weights
      !

      dynamics_split = config_dynamics_split
      if (config_split_dynamics_transport) then
        dt_dynamics = dt/real(dynamics_split)
        write(0,*) ' split dynamics-transport integration ',dynamics_split
      else
        dynamics_split = 1
        dt_dynamics = dt
        write(0,*) ' coupled RK3 dynamics-transport integration '
      end if
      if (.not. config_scalar_advection )  write(0,*) ' scalar advection turned off '

      number_of_sub_steps = config_number_of_sub_steps

      if(config_time_integration_order == 3) then

        rk_timestep(1) = dt_dynamics/3.
        rk_timestep(2) = dt_dynamics/2.
        rk_timestep(3) = dt_dynamics

        rk_sub_timestep(1) = dt_dynamics/3.
        rk_sub_timestep(2) = dt_dynamics/real(number_of_sub_steps)
        rk_sub_timestep(3) = dt_dynamics/real(number_of_sub_steps)

        number_sub_steps(1) = 1
        number_sub_steps(2) = max(1,number_of_sub_steps/2)
        number_sub_steps(3) = number_of_sub_steps

      else if (config_time_integration_order == 2) then

        rk_timestep(1) = dt_dynamics/2.
        rk_timestep(2) = dt_dynamics/2.
        rk_timestep(3) = dt_dynamics

        rk_sub_timestep(1) = dt_dynamics/real(number_of_sub_steps)
        rk_sub_timestep(2) = dt_dynamics/real(number_of_sub_steps)
        rk_sub_timestep(3) = dt_dynamics/real(number_of_sub_steps)

        number_sub_steps(1) = max(1,number_of_sub_steps/2)
        number_sub_steps(2) = max(1,number_of_sub_steps/2)
        number_sub_steps(3) = number_of_sub_steps

      end if
        
! theta_m
      call mpas_dmpar_exch_halo_field(theta_m_field)
 
! scalars
      call mpas_dmpar_exch_halo_field(scalars_field)

! pressure_p
      call mpas_dmpar_exch_halo_field(pressure_p_field)

! rtheta_p
      call mpas_dmpar_exch_halo_field(rtheta_p_field)


      call mpas_timer_start('atm_rk_integration_setup')

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'diag', diag)
         ! mesh is needed for atm_compute_moist_coefficients
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)

         call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

         call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
         call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

         call mpas_pool_get_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)
         call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadStart', vertexSolveThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadEnd', vertexSolveThreadEnd)

         call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
         call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

!$OMP PARALLEL DO
         do thread=1,nThreads
            call atm_rk_integration_setup(state, diag, &
                                          cellThreadStart(thread), cellThreadEnd(thread), &
                                          vertexThreadStart(thread), vertexThreadEnd(thread), &
                                          edgeThreadStart(thread), edgeThreadEnd(thread), &
                                          cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                          vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                          edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread))
         end do
!$OMP END PARALLEL DO

         block => block % next
      end do
      call mpas_timer_stop('atm_rk_integration_setup')

      call mpas_timer_start('atm_compute_moist_coefficients')
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'diag', diag)
         ! mesh is needed for atm_compute_moist_coefficients
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)

         call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

         call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
         call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

         call mpas_pool_get_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)
         call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadStart', vertexSolveThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadEnd', vertexSolveThreadEnd)

         call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
         call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

!$OMP PARALLEL DO
         do thread=1,nThreads
            call atm_compute_moist_coefficients( block % dimensions, state, diag, mesh, &     !MGD could do away with dimensions arg
                                                 cellThreadStart(thread), cellThreadEnd(thread), &
                                                 vertexThreadStart(thread), vertexThreadEnd(thread), &
                                                 edgeThreadStart(thread), edgeThreadEnd(thread), &
                                                 cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                                 vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                                 edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread))
         end do
!$OMP END PARALLEL DO

         block => block % next
      end do
      call mpas_timer_stop('atm_compute_moist_coefficients')

#ifdef DO_PHYSICS
      call mpas_timer_start('physics_addtend')
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'diag', diag)
         call mpas_pool_get_subpool(block % structs, 'tend', tend)
         call mpas_pool_get_subpool(block % structs, 'tend_physics', tend_physics)
         rk_step = 1
         dynamics_substep = 1
         call physics_get_tend( block, &
                                 mesh, &
                                state, &     
                                 diag, &
                                 tend, &
                         tend_physics, &
                      block % configs, &
                              rk_step, &
                     dynamics_substep, &
                      tend_ru_physics, &
                  tend_rtheta_physics )
         block => block % next
      end do
      call mpas_timer_stop('physics_addtend')
#endif


      DYNAMICS_SUBSTEPS : do dynamics_substep = 1, dynamics_split

         !  Compute the coefficients for the vertically implicit solve in the acoustic step.
         !  These coefficients will work for the first acoustic step in all cases.
         call mpas_timer_start('atm_compute_vert_imp_coefs')
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_get_subpool(block % structs, 'diag', diag)
            call mpas_pool_get_subpool(block % structs, 'tend', tend)

            call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

            call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

            call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
            call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
            call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
            call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

            call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
            call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
            call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
            call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

            rk_step = 1
!$OMP PARALLEL DO
            do thread=1,nThreads
               call atm_compute_vert_imp_coefs( state, mesh, diag, block % configs, nVertLevels, rk_sub_timestep(rk_step), &
                                          cellThreadStart(thread), cellThreadEnd(thread), &
                                          edgeThreadStart(thread), edgeThreadEnd(thread), &
                                          cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                          edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread))
            end do
!$OMP END PARALLEL DO
            block => block % next
         end do
         call mpas_timer_stop('atm_compute_vert_imp_coefs')


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
         ! BEGIN Runge-Kutta loop 
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

         RK3_DYNAMICS : do rk_step = 1, 3  ! Runge-Kutta loop

            ! recompute vertically implicit coefficients if necessary
            if( (config_time_integration_order == 3) .and. (rk_step == 2)) then

              !  Compute the coefficients for the vertically implicit solve in the acoustic step.
              !  These coefficients will work for the 2nd and 3rd acoustic steps (dt is the same for both).
              block => domain % blocklist
              do while (associated(block))
                 call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
                 call mpas_pool_get_subpool(block % structs, 'state', state)
                 call mpas_pool_get_subpool(block % structs, 'diag', diag)
                 call mpas_pool_get_subpool(block % structs, 'tend', tend)
   
                 call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

                 call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

                 call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
                 call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
                 call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
                 call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

                 call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
                 call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
                 call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
                 call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

!$OMP PARALLEL DO
                 do thread=1,nThreads
                    call atm_compute_vert_imp_coefs( state, mesh, diag, block % configs, nVertLevels, rk_sub_timestep(rk_step), &
                                                     cellThreadStart(thread), cellThreadEnd(thread), &
                                                     edgeThreadStart(thread), edgeThreadEnd(thread), &
                                                     cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                                     edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread))
                 end do
!$OMP END PARALLEL DO
                 block => block % next
              end do
            end if  

            call mpas_timer_start('atm_compute_dyn_tend')
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_subpool(block % structs, 'state', state)
               call mpas_pool_get_subpool(block % structs, 'diag', diag)
               call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
               call mpas_pool_get_subpool(block % structs, 'tend', tend)
               call mpas_pool_get_subpool(block % structs, 'tend_physics', tend_physics)
   
               call mpas_pool_get_dimension(mesh, 'nCells', nCells)
               call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
               call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
               call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

               call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

               call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

               call mpas_pool_get_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadStart', vertexSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadEnd', vertexSolveThreadEnd)

               call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

               allocate(delsq_theta(nVertLevels,nCells+1))
               allocate(delsq_w(nVertLevels,nCells+1))
!!               allocate(qtot(nVertLevels,nCells+1))  ! initializing this earlier in solution sequence
               allocate(delsq_divergence(nVertLevels,nCells+1))
               allocate(delsq_u(nVertLevels,nEdges+1))
!!               allocate(delsq_circulation(nVertLevels,nVertices+1))  ! no longer used -> removed 
               allocate(delsq_vorticity(nVertLevels,nVertices+1))
               allocate(dpdz(nVertLevels,nCells+1))

MPAS_GPTL_START('dyn_tend')
!$OMP PARALLEL DO
               do thread=1,nThreads
                  call atm_compute_dyn_tend( tend, tend_physics, state, diag, mesh, block % configs, nVertLevels, rk_step, dt, & 
                                             cellThreadStart(thread), cellThreadEnd(thread), &
                                             vertexThreadStart(thread), vertexThreadEnd(thread), &
                                             edgeThreadStart(thread), edgeThreadEnd(thread), &
                                             cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                             vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                             edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread))
               end do
!$OMP END PARALLEL DO
MPAS_GPTL_STOP('dyn_tend')

               deallocate(delsq_theta)
               deallocate(delsq_w)
!!               deallocate(qtot)  ! deallocation after dynamics step complete, see below
               deallocate(delsq_divergence)
               deallocate(delsq_u)
!!               deallocate(delsq_circulation)    ! no longer used -> removed 
               deallocate(delsq_vorticity)
               deallocate(dpdz)
   
               block => block % next
            end do
            call mpas_timer_stop('atm_compute_dyn_tend')

#ifdef DO_PHYSICS
!            call mpas_timer_start('physics_addtend')
!            block => domain % blocklist
!            do while (associated(block))
!               call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
!               call mpas_pool_get_subpool(block % structs, 'state', state)
!               call mpas_pool_get_subpool(block % structs, 'diag', diag)
!               call mpas_pool_get_subpool(block % structs, 'tend', tend)
!               call mpas_pool_get_subpool(block % structs, 'tend_physics', tend_physics)
!               call physics_addtend( block, &
!                                      mesh, &
!                                     state, &     
!                                      diag, &
!                                      tend, &
!                              tend_physics, &
!                           block % configs, &
!                                   rk_step, &
!                          dynamics_substep  )
!               block => block % next
!            end do
!            call mpas_timer_stop('physics_addtend')
#endif

            if (weight_iau >= 1.e-12) then      ! HA (June-16-2016)
               if (debug) write(0,*) ' add incremental analysis forcing to tendencies '
               block => domain % blocklist
 
               do while (associated(block))

                  call atm_add_tend_anal_incr ( block % structs, weight_iau, dt )

                  block => block % next
              end do
            end if

            !***********************************
            !  need tendencies at all edges of owned cells -
            !  we are solving for all edges of owned cells to minimize communications
            !  during the acoustic substeps
            !***********************************

! tend_u
            call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tend)
            call mpas_pool_get_field(tend, 'u', tend_u_field)
            call mpas_dmpar_exch_halo_field(tend_u_field, (/ 1 /))
   
            call mpas_timer_start('small_step_prep')

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
               call mpas_pool_get_subpool(block % structs, 'state', state)
               call mpas_pool_get_subpool(block % structs, 'diag', diag)
               call mpas_pool_get_subpool(block % structs, 'tend', tend)
   
               call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

               call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

               call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

               call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)
               
!$OMP PARALLEL DO
               do thread=1,nThreads
                  call atm_set_smlstep_pert_variables( tend, diag, mesh, block % configs, &
                                             cellThreadStart(thread), cellThreadEnd(thread), &
                                             edgeThreadStart(thread), edgeThreadEnd(thread), &
                                             cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                             edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread))
               end do
!$OMP END PARALLEL DO
               block => block % next
            end do
            call mpas_timer_stop('small_step_prep')

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! begin acoustic steps loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do small_step = 1, number_sub_steps(rk_step)

               call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
               call mpas_pool_get_field(diag, 'rho_pp', rho_pp_field)
               call mpas_dmpar_exch_halo_field(rho_pp_field, (/ 1 /))

               call mpas_timer_start('atm_advance_acoustic_step')
               block => domain % blocklist
               do while (associated(block))
                  call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
                  call mpas_pool_get_subpool(block % structs, 'state', state)
                  call mpas_pool_get_subpool(block % structs, 'diag', diag)
                  call mpas_pool_get_subpool(block % structs, 'tend', tend)

                  call mpas_pool_get_dimension(mesh, 'nCells', nCells)
                  call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

                  call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

                  call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
                  call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

                  call mpas_pool_get_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)
                  call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadStart', vertexSolveThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadEnd', vertexSolveThreadEnd)

                  call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
                  call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

!$OMP PARALLEL DO
                  do thread=1,nThreads
                     call atm_advance_acoustic_step( state, diag, tend,  mesh, block % configs, nCells, nVertLevels, &
                                             rk_sub_timestep(rk_step), small_step, &
                                             cellThreadStart(thread), cellThreadEnd(thread), &
                                             vertexThreadStart(thread), vertexThreadEnd(thread), &
                                             edgeThreadStart(thread), edgeThreadEnd(thread), &
                                             cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                             vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                             edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread))
                  end do
!$OMP END PARALLEL DO

                  block => block % next
               end do
               call mpas_timer_stop('atm_advance_acoustic_step')

  
! rtheta_pp
! This is the only communications needed during the acoustic steps because we solve for u on all edges of owned cells

               call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
               call mpas_pool_get_field(diag, 'rtheta_pp', rtheta_pp_field)
               call mpas_dmpar_exch_halo_field(rtheta_pp_field, (/ 1 /))

               call mpas_pool_get_field(diag, 'divergence_3d', divergence_3d_field)
               call mpas_dmpar_exch_halo_field(divergence_3d_field, (/ 1 /))

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

            call mpas_timer_start('atm_recover_large_step_variables')
            block => domain % blocklist
            do while (associated(block))

               call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
               call mpas_pool_get_subpool(block % structs, 'state', state)
               call mpas_pool_get_subpool(block % structs, 'diag', diag)
               call mpas_pool_get_subpool(block % structs, 'tend', tend)

               call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

               call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

               call mpas_pool_get_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadStart', vertexSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadEnd', vertexSolveThreadEnd)

               call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

!$OMP PARALLEL DO
               do thread=1,nThreads
                  call atm_recover_large_step_variables( state, diag, tend, mesh, block % configs, rk_timestep(rk_step), &
                                             number_sub_steps(rk_step), rk_step, &
                                             cellThreadStart(thread), cellThreadEnd(thread), &
                                             vertexThreadStart(thread), vertexThreadEnd(thread), &
                                             edgeThreadStart(thread), edgeThreadEnd(thread), &
                                             cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                             vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                             edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread))
               end do
!$OMP END PARALLEL DO

               block => block % next
            end do
            call mpas_timer_stop('atm_recover_large_step_variables')

            ! u
            !CR: SMALLER STENCIL?: call mpas_dmpar_exch_halo_field(block % state % time_levs(2) % state % u, (/ 3 /))
            call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
            call mpas_pool_get_field(state, 'u', u_field, 2)
            call mpas_dmpar_exch_halo_field(u_field)

            ! scalar advection: RK3 scheme of Skamarock and Gassmann (2011). 
            ! PD or monotonicity constraints applied only on the final Runge-Kutta substep.

            if (config_scalar_advection .and. (.not. config_split_dynamics_transport) ) then

               if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
                  call mpas_timer_start('atm_advance_scalars')
               else
                  call mpas_timer_start('atm_advance_scalars_mono')
MPAS_GPTL_START('atm_advance_scalars_mono')
               end if
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

                  call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

                  call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
                  call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

                  call mpas_pool_get_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)
                  call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadStart', vertexSolveThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadEnd', vertexSolveThreadEnd)

                  call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
                  call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
                  call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

                  allocate(scalar_old_arr(nVertLevels,nCells+1))
                  allocate(scalar_new_arr(nVertLevels,nCells+1))
                  allocate(s_max_arr(nVertLevels,nCells+1))
                  allocate(s_min_arr(nVertLevels,nCells+1))
                  allocate(scale_array(nVertLevels,2,nCells+1))
                  allocate(flux_array(nVertLevels,nEdges+1))
                  allocate(wdtn_arr(nVertLevels+1,nCells+1))
                  if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
                     allocate(horiz_flux_array(num_scalars,nVertLevels,nEdges+1))
                  else
                     allocate(flux_upwind_tmp_arr(nVertLevels,nEdges+1))
                     allocate(flux_tmp_arr(nVertLevels,nEdges+1))
                  end if

                  !
                  ! Note: The advance_scalars_mono routine can be used without limiting, and thus, encompasses 
                  !       the functionality of the advance_scalars routine; however, it is noticeably slower, 
                  !       so we use the advance_scalars routine for the first two RK substeps.
                  !
!$OMP PARALLEL DO
                  do thread=1,nThreads
                     if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
                        call atm_advance_scalars( tend, state, diag, mesh, block % configs, num_scalars, nCells, nVertLevels, rk_timestep(rk_step), &
                                                cellThreadStart(thread), cellThreadEnd(thread), &
                                                vertexThreadStart(thread), vertexThreadEnd(thread), &
                                                edgeThreadStart(thread), edgeThreadEnd(thread), &
                                                cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                                vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                                edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread), &
                                                horiz_flux_array, rk_step, config_time_integration_order, &
                                                advance_density=.false. )
                     else
         
                        block % domain = domain 
                        call atm_advance_scalars_mono( block, tend, state, diag, mesh, block % configs, nCells, nEdges, nVertLevels, rk_timestep(rk_step), &
                                                cellThreadStart(thread), cellThreadEnd(thread), &
                                                vertexThreadStart(thread), vertexThreadEnd(thread), &
                                                edgeThreadStart(thread), edgeThreadEnd(thread), &
                                                cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                                vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                                edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread), &
                                                scalar_old_arr, scalar_new_arr, s_max_arr, s_min_arr, wdtn_arr, &
                                                scale_array, flux_array, flux_upwind_tmp_arr, flux_tmp_arr, &
                                                advance_density=.false.)
                     end if
                  end do
!$OMP END PARALLEL DO

                  deallocate(scalar_old_arr)
                  deallocate(scalar_new_arr)
                  deallocate(s_max_arr)
                  deallocate(s_min_arr)
                  deallocate(scale_array)
                  deallocate(flux_array)
                  deallocate(wdtn_arr)
                  if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
                     deallocate(horiz_flux_array)
                  else
                     deallocate(flux_upwind_tmp_arr)
                     deallocate(flux_tmp_arr)
                  end if

                  block => block % next
               end do
               if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
                  call mpas_timer_stop('atm_advance_scalars')
               else
MPAS_GPTL_STOP('atm_advance_scalars_mono')
                  call mpas_timer_stop('atm_advance_scalars_mono')
               end if

            end if

            call mpas_timer_start('atm_compute_solve_diagnostics')
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_subpool(block % structs, 'state', state)
               call mpas_pool_get_subpool(block % structs, 'diag', diag)
               call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
   
               call mpas_pool_get_dimension(state, 'nEdges', nEdges)
               call mpas_pool_get_dimension(state, 'nVertices', nVertices)
               call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels)

               call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

               call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)

               call mpas_pool_get_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)

               call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)

               allocate(ke_vertex(nVertLevels,nVertices+1))
               allocate(ke_edge(nVertLevels,nEdges+1))
   
!$OMP PARALLEL DO
               do thread=1,nThreads
                  call atm_compute_solve_diagnostics(dt, state, 2, diag, mesh, block % configs, &
                                                     cellThreadStart(thread), cellThreadEnd(thread), &
                                                     vertexThreadStart(thread), vertexThreadEnd(thread), &
                                                     edgeThreadStart(thread), edgeThreadEnd(thread), rk_step)
               end do
!$OMP END PARALLEL DO

               deallocate(ke_vertex)
               deallocate(ke_edge)

               block => block % next
            end do
            call mpas_timer_stop('atm_compute_solve_diagnostics')


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
            if (config_scalar_advection .and. (.not. config_split_dynamics_transport) ) then
               call mpas_pool_get_field(state, 'scalars', scalars_field, 2)
               call mpas_dmpar_exch_halo_field(scalars_field)
            end if

         end do RK3_DYNAMICS

         if (dynamics_substep < dynamics_split) then
            call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
            call mpas_pool_get_field(state, 'theta_m', theta_m_field, 2)

            call mpas_dmpar_exch_halo_field(theta_m_field)
            call mpas_dmpar_exch_halo_field(pressure_p_field)
            call mpas_dmpar_exch_halo_field(rtheta_p_field)
         end if

         !  dynamics-transport split, WCS 18 November 2014
         !  (1) time level 1 needs to be set to time level 2
         !  (2) need to accumulate ruAvg and wwAvg over the dynamics substeps, prepare for use in transport
         !  Notes:  physics tendencies for scalars should be OK coming out of dynamics

         call mpas_timer_start('atm_rk_dynamics_substep_finish')
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_get_subpool(block % structs, 'diag', diag)

            call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

            call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
            call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
            call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
            call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

            call mpas_pool_get_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
            call mpas_pool_get_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)
            call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadStart', vertexSolveThreadStart)
            call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadEnd', vertexSolveThreadEnd)

            call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
            call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
            call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
            call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

!$OMP PARALLEL DO
            do thread=1,nThreads
               call atm_rk_dynamics_substep_finish(state, diag, dynamics_substep, dynamics_split, &
                                                cellThreadStart(thread), cellThreadEnd(thread), &
                                                vertexThreadStart(thread), vertexThreadEnd(thread), &
                                                edgeThreadStart(thread), edgeThreadEnd(thread), &
                                                cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                                vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                                edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread))
            end do
!$OMP END PARALLEL DO

            block => block % next
         end do
         call mpas_timer_stop('atm_rk_dynamics_substep_finish')

      end do DYNAMICS_SUBSTEPS


      deallocate(qtot)  !  we are finished with these now
      deallocate(tend_rtheta_physics)
      deallocate(tend_ru_physics)

      !
      !  split transport, at present RK3
      !

      if (config_scalar_advection .and. config_split_dynamics_transport) then

         rk_timestep(1) = dt/3.
         rk_timestep(2) = dt/2.
         rk_timestep(3) = dt
         !  switch for 2nd order time integration for scalar transport
         if(config_time_integration_order == 2) rk_timestep(1) = dt/2.

         RK3_SPLIT_TRANSPORT : do rk_step = 1, 3  ! Runge-Kutta loop


            if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
               call mpas_timer_start('atm_advance_scalars')
            else
               call mpas_timer_start('atm_advance_scalars_mono')
MPAS_GPTL_START('atm_advance_scalars_mono')
            end if
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

               call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

               call mpas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

               call mpas_pool_get_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadStart', vertexSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'vertexSolveThreadEnd', vertexSolveThreadEnd)

               call mpas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
               call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
               call mpas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

               allocate(scalar_old_arr(nVertLevels,nCells+1))
               allocate(scalar_new_arr(nVertLevels,nCells+1))
               allocate(s_max_arr(nVertLevels,nCells+1))
               allocate(s_min_arr(nVertLevels,nCells+1))
               allocate(scale_array(nVertLevels,2,nCells+1))
               allocate(flux_array(nVertLevels,nEdges+1))
               allocate(wdtn_arr(nVertLevels+1,nCells+1))
               allocate(rho_zz_int(nVertLevels,nCells+1))
               allocate(scalar_tend_array(num_scalars,nVertLevels,nCells+1))
               if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
                  allocate(horiz_flux_array(num_scalars,nVertLevels,nEdges+1))
               else
                  allocate(flux_upwind_tmp_arr(nVertLevels,nEdges+1))
                  allocate(flux_tmp_arr(nVertLevels,nEdges+1))
               end if

               !
               ! Note: The advance_scalars_mono routine can be used without limiting, and thus, encompasses 
               !       the functionality of the advance_scalars routine; however, it is noticeably slower, 
               !       so we use the advance_scalars routine for the first two RK substeps.
               !

               !  The latest version of atm_advance_scalars does not need the arrays scalar_tend_array or rho_zz_int
               !  We can remove scalar_tend_array????  WCS 20160921
!$OMP PARALLEL DO
               do thread=1,nThreads
                  if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
                     call atm_advance_scalars( tend, state, diag, mesh, block % configs, num_scalars, nCells, nVertLevels, rk_timestep(rk_step), &
                                             cellThreadStart(thread), cellThreadEnd(thread), &
                                             vertexThreadStart(thread), vertexThreadEnd(thread), &
                                             edgeThreadStart(thread), edgeThreadEnd(thread), &
                                             cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                             vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                             edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread), &
                                             horiz_flux_array, rk_step, config_time_integration_order, &
                                             advance_density=.true., scalar_tend=scalar_tend_array, rho_zz_int=rho_zz_int )
                  else
      
                     block % domain = domain 
                     call atm_advance_scalars_mono( block, tend, state, diag, mesh, block % configs, nCells, nEdges, nVertLevels, rk_timestep(rk_step), &
                                             cellThreadStart(thread), cellThreadEnd(thread), &
                                             vertexThreadStart(thread), vertexThreadEnd(thread), &
                                             edgeThreadStart(thread), edgeThreadEnd(thread), &
                                             cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                             vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                             edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread), &
                                             scalar_old_arr, scalar_new_arr, s_max_arr, s_min_arr, wdtn_arr, &
                                             scale_array, flux_array, flux_upwind_tmp_arr, flux_tmp_arr, &
                                             advance_density=.true., rho_zz_int=rho_zz_int)
                  end if
               end do
!$OMP END PARALLEL DO

               deallocate(scalar_old_arr)
               deallocate(scalar_new_arr)
               deallocate(s_max_arr)
               deallocate(s_min_arr)
               deallocate(scale_array)
               deallocate(flux_array)
               deallocate(wdtn_arr)
               deallocate(rho_zz_int)
               deallocate(scalar_tend_array)
               if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
                  deallocate(horiz_flux_array)
               else
                  deallocate(flux_upwind_tmp_arr)
                  deallocate(flux_tmp_arr)
               end if

               block => block % next
            end do
            if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
               call mpas_timer_stop('atm_advance_scalars')
            else
MPAS_GPTL_STOP('atm_advance_scalars_mono')
               call mpas_timer_stop('atm_advance_scalars_mono')
            end if

            if (rk_step < 3) then
               call mpas_pool_get_field(state, 'scalars', scalars_field, 2)
               call mpas_dmpar_exch_halo_field(scalars_field)
            end if

         end do RK3_SPLIT_TRANSPORT

      end if

      !
      ! reconstruct full velocity vectors at cell centers:
      !
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

      !
      ! call to parameterizations of cloud microphysics. calculation of the tendency of water vapor to horizontal and
      ! vertical advection needed for the Tiedtke parameterization of convection.
      !

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

         call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

         call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
         call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

         if(config_convection_scheme == 'cu_grell_freitas'          .or. &
            config_convection_scheme == 'cu_tiedtke'                .or. &
            config_convection_scheme == 'cu_ntiedtke') then

            call mpas_pool_get_array(tend_physics, 'rqvdynten', rqvdynten)

            !NOTE: The calculation of the tendency due to horizontal and vertical advection for the water vapor mixing ratio
            !requires that the subroutine atm_advance_scalars_mono was called on the third Runge Kutta step, so that a halo
            !update for the scalars at time_levs(1) is applied. A halo update for the scalars at time_levs(2) is done above. 
            if (config_monotonic) then
               rqvdynten(:,:) = ( scalars_2(index_qv,:,:) - scalars_1(index_qv,:,:) ) / config_dt
            else
               rqvdynten(:,:) = 0._RKIND
            end if
         end if

         !simply set to zero negative mixing ratios of different water species (for now):
         where ( scalars_2(:,:,:) < 0.0)  &
            scalars_2(:,:,:) = 0.0

         !call microphysics schemes:
         if (trim(config_microp_scheme) /= 'off')  then
            call mpas_timer_start('microphysics')
!$OMP PARALLEL DO
            do thread=1,nThreads
               call driver_microphysics ( block % configs, mesh, state, 2, diag, diag_physics, tend, itimestep, &
                                          cellSolveThreadStart(thread), cellSolveThreadEnd(thread))
            end do
!$OMP END PARALLEL DO
            call mpas_timer_stop('microphysics')
         end if
         block => block % next
      end do
#endif

      call summarize_timestep(domain)

   end subroutine atm_srk3


   subroutine atm_rk_integration_setup( state, diag, &
                                   cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)

      implicit none

      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd

      real (kind=RKIND), dimension(:,:), pointer :: ru
      real (kind=RKIND), dimension(:,:), pointer :: ru_save
      real (kind=RKIND), dimension(:,:), pointer :: rw
      real (kind=RKIND), dimension(:,:), pointer :: rw_save
      real (kind=RKIND), dimension(:,:), pointer :: rtheta_p
      real (kind=RKIND), dimension(:,:), pointer :: rtheta_p_save
      real (kind=RKIND), dimension(:,:), pointer :: rho_p
      real (kind=RKIND), dimension(:,:), pointer :: rho_p_save
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz_old_split

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
      call mpas_pool_get_array(diag, 'rho_zz_old_split', rho_zz_old_split)

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

      ru_save(:,edgeStart:edgeEnd) = ru(:,edgeStart:edgeEnd)
      rw_save(:,cellStart:cellEnd) = rw(:,cellStart:cellEnd)
      rtheta_p_save(:,cellStart:cellEnd) = rtheta_p(:,cellStart:cellEnd)
      rho_p_save(:,cellStart:cellEnd) = rho_p(:,cellStart:cellEnd)

      u_2(:,edgeStart:edgeEnd) = u_1(:,edgeStart:edgeEnd)
      w_2(:,cellStart:cellEnd) = w_1(:,cellStart:cellEnd)
      theta_m_2(:,cellStart:cellEnd) = theta_m_1(:,cellStart:cellEnd)
      rho_zz_2(:,cellStart:cellEnd) = rho_zz_1(:,cellStart:cellEnd)
      rho_zz_old_split(:,cellStart:cellEnd) =  rho_zz_1(:,cellStart:cellEnd)
      scalars_2(:,:,cellStart:cellEnd) = scalars_1(:,:,cellStart:cellEnd)

   end subroutine atm_rk_integration_setup


   subroutine atm_compute_moist_coefficients( dims, state, diag, mesh, &
                                   cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)

      ! the moist coefficients cqu and cqw serve to transform the inverse dry density (1/rho_d) 
      ! into the inverse full (moist) density (1/rho_m).

      implicit none

      type (mpas_pool_type), intent(in) :: dims
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd

      integer :: iEdge, iCell, k, cell1, cell2, iq
      integer, pointer :: nCells, nEdges, nVertLevels, nCellsSolve
      real (kind=RKIND) :: qtotal
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

!      do iCell = cellSolveStart,cellSolveEnd
      do iCell = cellStart,cellEnd
         qtot(1:nVertLevels,iCell) = 0.0
         do k = 1,nVertLevels
            do iq = moist_start, moist_end
               qtot(k,iCell) = qtot(k,iCell) + scalars(iq, k, iCell)
            end do
         end do
      end do    

!      do iCell = cellSolveStart,cellSolveEnd
      do iCell = cellStart,cellEnd
         do k = 2, nVertLevels
            qtotal = 0.5*(qtot(k,iCell)+qtot(k-1,iCell))
            cqw(k,iCell) = 1.0 / (1.0 + qtotal)
         end do
      end do

!  would need to compute qtot for all cells and an openmp barrier to use qtot below.

      do iEdge = edgeStart,edgeEnd
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k = 1, nVertLevels
               qtotal = 0.0
               do iq = moist_start, moist_end
                  qtotal = qtotal + 0.5 * ( scalars(iq, k, cell1) + scalars(iq, k, cell2) )
               end do
               cqu(k,iEdge) = 1.0 / (1.0 + qtotal)
            end do
         end if
      end do

   end subroutine atm_compute_moist_coefficients


   subroutine atm_compute_vert_imp_coefs(state, mesh, diag, configs, nVertLevels, dts, &
                                   cellStart, cellEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd)
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
      integer, intent(in) :: cellStart, cellEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd


      real (kind=RKIND), dimension(:,:), pointer :: zz, cqw, p, t, rb, rtb, pb, rt
      real (kind=RKIND), dimension(:,:), pointer :: cofwr, cofwz, coftz, cofwt, a_tri, alpha_tri, gamma_tri
      real (kind=RKIND), dimension(:), pointer :: cofrz, rdzw, fzm, fzp, rdzu
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      real (kind=RKIND), pointer :: epssm

      integer, pointer :: nCells, moist_start, moist_end


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

      call mpas_pool_get_dimension(state, 'nCells', nCells)
      call mpas_pool_get_dimension(state, 'moist_start', moist_start)
      call mpas_pool_get_dimension(state, 'moist_end', moist_end)


      call atm_compute_vert_imp_coefs_work(nCells, moist_start, moist_end, dts, epssm, &
                                   zz, cqw, p, t, rb, rtb, pb, rt, cofwr, cofwz, coftz, cofwt, &
                                   a_tri, alpha_tri, gamma_tri, cofrz, rdzw, fzm, fzp, rdzu, scalars, &
                                   cellStart, cellEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd)


   end subroutine atm_compute_vert_imp_coefs


   subroutine atm_compute_vert_imp_coefs_work(nCells, moist_start, moist_end, dts, epssm, &
                                   zz, cqw, p, t, rb, rtb, pb, rt, cofwr, cofwz, coftz, cofwt, &
                                   a_tri, alpha_tri, gamma_tri, cofrz, rdzw, fzm, fzp, rdzu, scalars, &
                                   cellStart, cellEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd)

      use mpas_atm_dimensions
 
      implicit none


      !
      ! Dummy arguments
      !
      integer, intent(in) :: nCells, moist_start, moist_end
      real (kind=RKIND), intent(in) :: dts
      real (kind=RKIND), intent(in) :: epssm

      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: zz
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: cqw
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: p
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: t
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rb
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rtb
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: pb
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rt
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: cofwr
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: cofwz
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: coftz
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: cofwt
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: a_tri
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: alpha_tri
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: gamma_tri
      real (kind=RKIND), dimension(nVertLevels) :: cofrz
      real (kind=RKIND), dimension(nVertLevels) :: rdzw
      real (kind=RKIND), dimension(nVertLevels) :: fzm
      real (kind=RKIND), dimension(nVertLevels) :: fzp
      real (kind=RKIND), dimension(nVertLevels) :: rdzu
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1) :: scalars

      integer, intent(in) :: cellStart, cellEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd


      !
      ! Local variables
      !
      integer :: iCell, k, iq
      real (kind=RKIND) :: dtseps, c2, qtotal, rcv
      real (kind=RKIND), dimension( nVertLevels ) :: b_tri, c_tri


      !  set coefficients
      dtseps = .5*dts*(1.+epssm)
      rcv = rgas/(cp-rgas)
      c2 = cp*rcv

! MGD bad to have all threads setting this variable?
      do k=1,nVertLevels
         cofrz(k) = dtseps*rdzw(k)
      end do

      do iCell = cellSolveStart,cellSolveEnd  !  we only need to do cells we are solving for, not halo cells

!DIR$ IVDEP
         do k=2,nVertLevels
            cofwr(k,iCell) =.5*dtseps*gravity*(fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))
         end do
         coftz(1,iCell) = 0.0
!DIR$ IVDEP
         do k=2,nVertLevels
            cofwz(k,iCell) = dtseps*c2*(fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))  &
                 *rdzu(k)*cqw(k,iCell)*(fzm(k)*p (k,iCell)+fzp(k)*p (k-1,iCell))
            coftz(k,iCell) = dtseps*   (fzm(k)*t (k,iCell)+fzp(k)*t (k-1,iCell))
         end do
         coftz(nVertLevels+1,iCell) = 0.0
!DIR$ IVDEP
         do k=1,nVertLevels

!            qtotal = 0.
!            do iq = moist_start, moist_end
!               qtotal = qtotal + scalars(iq, k, iCell)
!            end do
            qtotal = qtot(k,iCell)

            cofwt(k,iCell) = .5*dtseps*rcv*zz(k,iCell)*gravity*rb(k,iCell)/(1.+qtotal)  &
                                *p(k,iCell)/((rtb(k,iCell)+rt(k,iCell))*pb(k,iCell))
!            cofwt(k,iCell) = 0.
         end do

         a_tri(1,iCell) = 0.  ! note, this value is never used
         b_tri(1) = 1.    ! note, this value is never used
         c_tri(1) = 0.    ! note, this value is never used
         gamma_tri(1,iCell) = 0.
         alpha_tri(1,iCell) = 0.  ! note, this value is never used

!DIR$ IVDEP
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
!MGD VECTOR DEPENDENCE
         do k=2,nVertLevels
            alpha_tri(k,iCell) = 1./(b_tri(k)-a_tri(k,iCell)*gamma_tri(k-1,iCell))
            gamma_tri(k,iCell) = c_tri(k)*alpha_tri(k,iCell)
         end do

      end do ! loop over cells

   end subroutine atm_compute_vert_imp_coefs_work


   subroutine atm_set_smlstep_pert_variables( tend, diag, mesh, configs, &
                                   cellStart, cellEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd)

      ! following Klemp et al MWR 2007, we use preturbation variables
      ! in the acoustic-step integration.  This routine computes those 
      ! perturbation variables.  state variables are reconstituted after 
      ! the acousstic steps in subroutine atm_recover_large_step_variables


      implicit none

      type (mpas_pool_type), intent(inout) :: tend
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: cellStart, cellEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd

      integer, pointer :: nCells, nEdges, nCellsSolve
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp
      real (kind=RKIND), dimension(:,:), pointer :: ruAvg, wwAvg
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3, zb_cell, zb3_cell
      real (kind=RKIND), dimension(:,:), pointer :: zz
      real (kind=RKIND), dimension(:,:), pointer :: w_tend, u_tend
      real (kind=RKIND), dimension(:,:), pointer :: rho_pp, rho_p_save, rho_p
      real (kind=RKIND), dimension(:,:), pointer :: ru_p, ru, ru_save
      real (kind=RKIND), dimension(:,:), pointer :: rtheta_pp, rtheta_p_save, rtheta_p, zz_rtheta_pp_old
      real (kind=RKIND), dimension(:,:), pointer :: rw_p, rw_save, rw
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign


      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)

      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)
      call mpas_pool_get_array(mesh, 'zb_cell', zb_cell)
      call mpas_pool_get_array(mesh, 'zb3_cell', zb3_cell)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)

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
      call mpas_pool_get_array(diag, 'zz_rtheta_pp_old', zz_rtheta_pp_old)

      call mpas_pool_get_array(diag, 'rw_p', rw_p)
      call mpas_pool_get_array(diag, 'rw_save', rw_save)
      call mpas_pool_get_array(diag, 'rw', rw)

      call atm_set_smlstep_pert_variables_work(nCells, nEdges, nCellsSolve, &
                                   nEdgesOnCell, cellsOnEdge, edgesOnCell, fzm, fzp, ruAvg, wwAvg, zb, zb3, zb_cell, zb3_cell, &
                                   zz, w_tend, u_tend, rho_pp, rho_p_save, rho_p, ru_p, ru, ru_save, &
                                   rtheta_pp, rtheta_p_save, rtheta_p, zz_rtheta_pp_old, rw_p, rw_save, rw, edgesOnCell_sign, &
                                   cellStart, cellEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd)


   end subroutine atm_set_smlstep_pert_variables


   subroutine atm_set_smlstep_pert_variables_work(nCells, nEdges, nCellsSolve, &
                                   nEdgesOnCell, cellsOnEdge, edgesOnCell, fzm, fzp, ruAvg, wwAvg, zb, zb3, zb_cell, zb3_cell, &
                                   zz, w_tend, u_tend, rho_pp, rho_p_save, rho_p, ru_p, ru, ru_save, &
                                   rtheta_pp, rtheta_p_save, rtheta_p, zz_rtheta_pp_old, rw_p, rw_save, rw, edgesOnCell_sign, &
                                   cellStart, cellEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd)

      use mpas_atm_dimensions

      implicit none


      !
      ! Dummy arguments
      !
      integer, intent(in) :: nCells, nEdges, nCellsSolve

      integer, intent(in) :: cellStart, cellEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd

      integer, dimension(nCells+1) :: nEdgesOnCell
      integer, dimension(2,nEdges+1) :: cellsOnEdge
      integer, dimension(maxEdges,nCells+1) :: edgesOnCell
      real (kind=RKIND), dimension(nVertLevels) :: fzm
      real (kind=RKIND), dimension(nVertLevels) :: fzp
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ruAvg
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: wwAvg
      real (kind=RKIND), dimension(nVertLevels+1,2,nEdges+1) :: zb
      real (kind=RKIND), dimension(nVertLevels+1,2,nEdges+1) :: zb3
      real (kind=RKIND), dimension(nVertLevels+1,maxEdges,nCells+1) :: zb_cell
      real (kind=RKIND), dimension(nVertLevels+1,maxEdges,nCells+1) :: zb3_cell
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: zz
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: w_tend
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: u_tend
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_pp
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_p_save
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_p
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru_p
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru_save
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rtheta_pp
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rtheta_p_save
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rtheta_p
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: zz_rtheta_pp_old
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw_p
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw_save
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw
      real (kind=RKIND), dimension(maxEdges,nCells+1) :: edgesOnCell_sign

      !
      ! Local variables
      !
      integer :: iCell, iEdge, i, k
      real (kind=RKIND) :: flux

      ! we solve for omega instead of w (see Klemp et al MWR 2007),
      ! so here we change the w_p tendency to an omega_p tendency

      ! here we need to compute the omega tendency in a manner consistent with our diagnosis of omega.
      ! this requires us to use the same flux divergence as is used in the theta eqn - see Klemp et al MWR 2003.

!!      do iCell=cellStart,cellEnd
      do iCell=cellSolveStart,cellSolveEnd
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
!DIR$ IVDEP
            do k = 2, nVertLevels
               flux = edgesOnCell_sign(i,iCell) * (fzm(k) * u_tend(k,iEdge) + fzp(k) * u_tend(k-1,iEdge))
               w_tend(k,iCell) = w_tend(k,iCell)   &
                        - (zb_cell(k,i,iCell) + sign(1.0_RKIND, u_tend(k,iEdge)) * zb3_cell(k,i,iCell)) * flux
            end do
         end do
!DIR$ IVDEP
         do k = 2, nVertLevels
            w_tend(k,iCell) = ( fzm(k) * zz(k,iCell) + fzp(k) * zz(k-1,iCell)   ) * w_tend(k,iCell)
         end do
      end do

   end subroutine atm_set_smlstep_pert_variables_work


   subroutine atm_advance_acoustic_step( state, diag, tend, mesh, configs, nCells, nVertLevels, dts, small_step, &
                                   cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)

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
      integer, intent(in) :: small_step              ! acoustic step number
      integer, intent(in) :: nCells                  ! for allocating stack variables
      integer, intent(in) :: nVertLevels             ! for allocating stack variables
      real (kind=RKIND), intent(in) :: dts
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd


      real (kind=RKIND), dimension(nVertLevels) :: du
      real (kind=RKIND), dimension(nVertLevels+1) :: dpzx

      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, theta_m, ru_p, rw_p, rtheta_pp,  &
                                                    zz_rtheta_pp_old, zz, exner, cqu, ruAvg, &
                                                    wwAvg, rho_pp, cofwt, coftz, zxu,        &
                                                    a_tri, alpha_tri, gamma_tri, dss,        &
                                                    tend_ru, tend_rho, tend_rt, tend_rw,     &
                                                    zgrid, cofwr, cofwz, w, divergence_3d

! redefine ru_p to be perturbation from time t, change 3a  ! temporary
      real (kind=RKIND), dimension(:,:), pointer :: ru
      real (kind=RKIND), dimension(:,:), pointer :: ru_save
! redefine rw_p to be perturbation from time t, change 3a  ! temporary
      real (kind=RKIND), dimension(:,:), pointer :: rw
      real (kind=RKIND), dimension(:,:), pointer :: rw_save

      real (kind=RKIND), dimension(:), pointer :: fzm, fzp, rdzw, dcEdge, invDcEdge, invAreaCell, cofrz, dvEdge

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign

      real (kind=RKIND), pointer :: epssm, smdiv, smdiv_p_forward

      real (kind=RKIND), pointer :: cf1, cf2, cf3

      integer, pointer :: nEdges, nCellsSolve

      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
!      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
!  change needed for rw_p, change 6 (see rayleigh damping)
      call mpas_pool_get_array(state, 'w', w, 2)
!      call mpas_pool_get_array(state, 'w', w, 1)

      call mpas_pool_get_array(diag, 'rtheta_pp', rtheta_pp)
      call mpas_pool_get_array(diag, 'zz_rtheta_pp_old', zz_rtheta_pp_old)
      call mpas_pool_get_array(diag, 'ru_p', ru_p)
      call mpas_pool_get_array(diag, 'rw_p', rw_p)
      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'cqu', cqu)
      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
      call mpas_pool_get_array(diag, 'rho_pp', rho_pp)
      call mpas_pool_get_array(diag, 'divergence_3d', divergence_3d)
      call mpas_pool_get_array(diag, 'cofwt', cofwt)
      call mpas_pool_get_array(diag, 'coftz', coftz)
      call mpas_pool_get_array(diag, 'cofrz', cofrz)
      call mpas_pool_get_array(diag, 'cofwr', cofwr)
      call mpas_pool_get_array(diag, 'cofwz', cofwz)
      call mpas_pool_get_array(diag, 'a_tri', a_tri)
      call mpas_pool_get_array(diag, 'alpha_tri', alpha_tri)
      call mpas_pool_get_array(diag, 'gamma_tri', gamma_tri)

      call mpas_pool_get_array(mesh, 'dss', dss)

      call mpas_pool_get_array(tend, 'u', tend_ru)
      call mpas_pool_get_array(tend, 'rho_zz', tend_rho)
      call mpas_pool_get_array(tend, 'theta_m', tend_rt)
      call mpas_pool_get_array(tend, 'w', tend_rw)

      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zxu', zxu)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'invDcEdge', invDcEdge)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)

      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)

! redefine ru_p to be perturbation from time t, change 3b ! temporary
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'ru_save', ru_save)
! redefine rw_p to be perturbation from time t, change 3b ! temporary
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'rw_save', rw_save)

      ! epssm is the offcentering coefficient for the vertically implicit integration.
      ! smdiv is the 3D divergence-damping coefficients.
      call mpas_pool_get_config(configs, 'config_epssm', epssm) 
      call mpas_pool_get_config(configs, 'config_smdiv', smdiv) 
      call mpas_pool_get_config(configs, 'config_smdiv_p_forward', smdiv_p_forward) 

      call atm_advance_acoustic_step_work(nCells, nEdges, nCellsSolve, cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
                                   rho_zz, theta_m, ru_p, rw_p, rtheta_pp, zz_rtheta_pp_old, zz, exner, cqu, ruAvg, wwAvg, &
                                   rho_pp, cofwt, coftz, zxu, a_tri, alpha_tri, gamma_tri, dss, tend_ru, tend_rho, tend_rt, &
                                   tend_rw, zgrid, cofwr, cofwz, w, ru, ru_save, rw, rw_save, divergence_3d, fzm, fzp, rdzw, dcEdge, invDcEdge, &
                                   invAreaCell, cofrz, dvEdge, nEdgesOnCell, cellsOnEdge, edgesOnCell, edgesOnCell_sign, &
                                   dts, small_step, epssm, smdiv, smdiv_p_forward, cf1, cf2, cf3 &
                                   )

   end subroutine atm_advance_acoustic_step


   subroutine atm_advance_acoustic_step_work(nCells, nEdges, nCellsSolve, cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
                                   rho_zz, theta_m, ru_p, rw_p, rtheta_pp, zz_rtheta_pp_old, zz, exner, cqu, ruAvg, wwAvg, &
                                   rho_pp, cofwt, coftz, zxu, a_tri, alpha_tri, gamma_tri, dss, tend_ru, tend_rho, tend_rt, &
                                   tend_rw, zgrid, cofwr, cofwz, w, ru, ru_save, rw, rw_save, divergence_3d, fzm, fzp, rdzw, dcEdge, invDcEdge, &
                                   invAreaCell, cofrz, dvEdge, nEdgesOnCell, cellsOnEdge, edgesOnCell, edgesOnCell_sign, &
                                   dts, small_step, epssm, smdiv, smdiv_p_forward, cf1, cf2, cf3 &
                                   )

      use mpas_atm_dimensions

      implicit none


      !
      ! Dummy arguments
      !
      integer, intent(in) :: nCells, nEdges, nCellsSolve
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd

      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_zz
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: theta_m
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru_p
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw_p
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rtheta_pp

      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: zz_rtheta_pp_old
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: divergence_3d
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: zz
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: exner
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: cqu
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ruAvg
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: wwAvg
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_pp
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: cofwt
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: coftz
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: zxu

      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: a_tri
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: alpha_tri
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: gamma_tri
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: dss
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: tend_ru
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: tend_rho
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: tend_rt
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: tend_rw
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: zgrid
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: cofwr
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: cofwz

      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: w
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru_save
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw_save

      real (kind=RKIND), dimension(nVertLevels) :: fzm
      real (kind=RKIND), dimension(nVertLevels) :: fzp
      real (kind=RKIND), dimension(nVertLevels) :: rdzw
      real (kind=RKIND), dimension(nEdges+1) :: dcEdge
      real (kind=RKIND), dimension(nEdges+1) :: invDcEdge
      real (kind=RKIND), dimension(nCells+1) :: invAreaCell
      real (kind=RKIND), dimension(nVertLevels) :: cofrz
      real (kind=RKIND), dimension(nEdges+1) :: dvEdge

      integer, dimension(nCells+1) :: nEdgesOnCell
      integer, dimension(2,nEdges+1) :: cellsOnEdge
      integer, dimension(maxEdges,nCells+1) :: edgesOnCell
      real (kind=RKIND), dimension(maxEdges,nCells+1) :: edgesOnCell_sign

      integer, intent(in) :: small_step
      real (kind=RKIND), intent(in) :: dts, epssm, smdiv, smdiv_p_forward, cf1, cf2, cf3
      real (kind=RKIND), dimension(nVertLevels) :: ts, rs

  
      !
      ! Local variables
      !
      integer :: cell1, cell2, iEdge, iCell, i, k
      real (kind=RKIND) :: c2, rcv, rtheta_pp_tmp
      real (kind=RKIND) :: pgrad, flux, resm, rdts


      rcv = rgas / (cp - rgas)
      c2 = cp * rcv
      resm = (1.0 - epssm) / (1.0 + epssm)
      rdts = 1./dts

      if(small_step /= 1) then  !  not needed on first small step 

         do iCell = cellStart,cellEnd
           ! acoustic step divergence damping - forward weight rtheta_pp - see Klemp et al MWR 2007
           do k = 1,nVertLevels
              rtheta_pp_tmp = rtheta_pp(k,iCell)
              rtheta_pp(k,iCell) = (rtheta_pp(k,iCell) + smdiv_p_forward * (rtheta_pp(k,iCell)-zz_rtheta_pp_old(k,iCell)))*zz(k,iCell)
              zz_rtheta_pp_old(k,iCell) = rtheta_pp_tmp
           end do
        end do
        
!$OMP BARRIER

        ! forward-backward acoustic step integration.
        ! begin by updating the horizontal velocity u, 
        ! and accumulating the contribution from the updated u to the other tendencies.

        ! we are looping over all edges, but only computing on edges of owned cells. This will include updates of
        ! all owned edges plus some edges that are owned by other blocks.  We perform these redundant computations
        ! so that we do not have to communicate updates of u to update the cell variables (rho, w, and theta). 

        !MGD this loop will not be very load balanced with if-test below

        do iEdge=edgeStart,edgeEnd ! MGD do we really just need edges touching owned cells?
 
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)

           ! update edges for block-owned cells
           if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

!DIR$ IVDEP
              do k=1,nVertLevels
!!                 pgrad = ((zz_rtheta_pp(k,cell2)-zz_rtheta_pp_old(k,cell1))*invDcEdge(iEdge) )/(.5*(zz(k,cell2)+zz(k,cell1)))
                 pgrad = ((rtheta_pp(k,cell2)-rtheta_pp(k,cell1))*invDcEdge(iEdge) )/(.5*(zz(k,cell2)+zz(k,cell1)))
                 pgrad = cqu(k,iEdge)*0.5*c2*(exner(k,cell1)+exner(k,cell2))*pgrad
                 pgrad = pgrad + 0.5*zxu(k,iEdge)*gravity*(rho_pp(k,cell1)+rho_pp(k,cell2))
                 ru_p(k,iEdge) = ru_p(k,iEdge) + dts*(tend_ru(k,iEdge) - pgrad) &
                                  - smdiv*dcEdge(iEdge)*(divergence_3d(k,cell2)-divergence_3d(k,cell1))
              end do

              ! accumulate ru_p for use later in scalar transport
!DIR$ IVDEP
              do k=1,nVertLevels
                 ruAvg(k,iEdge) = ruAvg(k,iEdge) + ru_p(k,iEdge)
              end do

           end if ! end test for block-owned cells

        end do ! end loop over edges

      else !  this is all that us needed for ru_p update for first acoustic step in RK substep

        do iEdge=edgeStart,edgeEnd ! MGD do we really just need edges touching owned cells?

           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)

           ! update edges for block-owned cells
           if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

!DIR$ IVDEP
              do k=1,nVertLevels
                 ru_p(k,iEdge) = dts*tend_ru(k,iEdge) - smdiv*dcEdge(iEdge)*(tend_rho(k,cell2)-tend_rho(k,cell1))
              end do
!DIR$ IVDEP
              do k=1,nVertLevels
!!                 ruAvg(k,iEdge) = ruAvg(k,iEdge) + ru_p(k,iEdge)
                 ruAvg(k,iEdge) = ru_p(k,iEdge)                 
              end do

           end if ! end test for block-owned cells
        end do ! end loop over edges

      end if ! test for first acoustic step

!$OMP BARRIER

      if (small_step == 1) then  ! initialize here on first small timestep.
         do iCell=cellStart,cellEnd
            zz_rtheta_pp_old(1:nVertLevels,iCell) = 0.0
         end do
      end if

!!!OMP BARRIER -- not needed, since zz_rtheta_pp_old not used below when small_step == 1

      do iCell=cellSolveStart,cellSolveEnd  ! loop over all owned cells to solve

         ts(:) = 0.0
         rs(:) = 0.0

         if(small_step == 1) then  ! initialize here on first small timestep.
            wwAvg(1:nVertLevels+1,iCell) = 0.0            
            rho_pp(1:nVertLevels,iCell) = 0.0            
            rtheta_pp(1:nVertLevels,iCell) = 0.0            
!MGD moved to loop above over all cells
!            zz_rtheta_pp_old(1:nVertLevels,iCell) = 0.0
            rw_p(:,iCell) = 0.0
            divergence_3d(1:nVertLevels,iCell) = 0.
         else  ! reset rtheta_pp to input value;
               !  rtheta_pp_old stores input value for use in div damping on next acoustic step.
               !  Save rho_pp to compute d_rho_pp/dt to get divergence for next acoustic filter application.
            rtheta_pp(1:nVertLevels,iCell) = zz_rtheta_pp_old(1:nVertLevels,iCell)  
            divergence_3d(1:nVertLevels,iCell) = rho_pp(1:nVertLevels,iCell)
         end if
            
         do i=1,nEdgesOnCell(iCell) 
            iEdge = edgesOnCell(i,iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
            do k=1,nVertLevels
               flux = edgesOnCell_sign(i,iCell)*dts*dvEdge(iEdge)*ru_p(k,iEdge) * invAreaCell(iCell)
               rs(k) = rs(k)-flux
               ts(k) = ts(k)-flux*0.5*(theta_m(k,cell2)+theta_m(k,cell1))
            end do
         end do

      ! vertically implicit acoustic and gravity wave integration.
      ! this follows Klemp et al MWR 2007, with the addition of an implicit Rayleigh damping of w
      ! serves as a gravity-wave absorbing layer, from Klemp et al 2008.

!DIR$ IVDEP
         do k=1, nVertLevels
            rs(k) = rho_pp(k,iCell) + dts*tend_rho(k,iCell) + rs(k)                  &
                            - cofrz(k)*resm*(rw_p(k+1,iCell)-rw_p(k,iCell)) 
            ts(k) = rtheta_pp(k,iCell) + dts*tend_rt(k,iCell) + ts(k)                &
                               - resm*rdzw(k)*( coftz(k+1,iCell)*rw_p(k+1,iCell)     &
                                               -coftz(k,iCell)*rw_p(k,iCell))
         end do

!DIR$ IVDEP
         do k=2, nVertLevels
            wwavg(k,iCell) = wwavg(k,iCell) + 0.5*(1.0-epssm)*rw_p(k,iCell)
         end do

!DIR$ IVDEP
         do k=2, nVertLevels
            rw_p(k,iCell) = rw_p(k,iCell) +  dts*tend_rw(k,iCell)                       &
                       - cofwz(k,iCell)*((zz(k  ,iCell)*ts(k)                           &
                                     -zz(k-1,iCell)*ts(k-1))                            &
                               +resm*(zz(k  ,iCell)*rtheta_pp(k  ,iCell)                &
                                     -zz(k-1,iCell)*rtheta_pp(k-1,iCell)))              &
                       - cofwr(k,iCell)*((rs(k)+rs(k-1))                                &
                               +resm*(rho_pp(k,iCell)+rho_pp(k-1,iCell)))               &
                       + cofwt(k  ,iCell)*(ts(k  )+resm*rtheta_pp(k  ,iCell))           &
                       + cofwt(k-1,iCell)*(ts(k-1)+resm*rtheta_pp(k-1,iCell))
         end do

         ! tridiagonal solve sweeping up and then down the column

!MGD VECTOR DEPENDENCE
         do k=2,nVertLevels
            rw_p(k,iCell) = (rw_p(k,iCell)-a_tri(k,iCell)*rw_p(k-1,iCell))*alpha_tri(k,iCell)
         end do

!MGD VECTOR DEPENDENCE
         do k=nVertLevels,1,-1
            rw_p(k,iCell) = rw_p(k,iCell) - gamma_tri(k,iCell)*rw_p(k+1,iCell)     
         end do

         ! the implicit Rayleigh damping on w (gravity-wave absorbing) 

!DIR$ IVDEP
         do k=2,nVertLevels
            rw_p(k,iCell) = (rw_p(k,iCell) + (rw_save(k  ,iCell) - rw(k  ,iCell)) -dts*dss(k,iCell)*               &
                        (fzm(k)*zz (k,iCell)+fzp(k)*zz (k-1,iCell))        &
                        *(fzm(k)*rho_zz(k,iCell)+fzp(k)*rho_zz(k-1,iCell))       &
                                 *w(k,iCell)    )/(1.0+dts*dss(k,iCell)) &
                         - (rw_save(k  ,iCell) - rw(k  ,iCell))
         end do

         ! accumulate (rho*omega)' for use later in scalar transport
!DIR$ IVDEP
         do k=2,nVertLevels
            wwAvg(k,iCell) = wwAvg(k,iCell) + 0.5*(1.0+epssm)*rw_p(k,iCell)
         end do

         ! update rho_pp and theta_pp given updated rw_p

!DIR$ IVDEP
         do k=1,nVertLevels
            rho_pp(k,iCell) = rs(k) - cofrz(k) *(rw_p(k+1,iCell)-rw_p(k  ,iCell))
            rtheta_pp(k,iCell) = ts(k) - rdzw(k)*(coftz(k+1,iCell)*rw_p(k+1,iCell)  &
                               -coftz(k  ,iCell)*rw_p(k  ,iCell))
            divergence_3d(k,iCell) = (rho_pp(k,iCell) - divergence_3d(k,iCell))*rdts
         end do

      end do !  end of loop over cells

   end subroutine atm_advance_acoustic_step_work


   subroutine atm_recover_large_step_variables( state, diag, tend, mesh, configs, dt, ns, rk_step, &
                                       cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                       cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)

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
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd


      real (kind=RKIND), dimension(:,:), pointer :: wwAvg, rw_save, w, rw, rw_p, rtheta_p, rtheta_pp,   &
                                                    rtheta_p_save, rt_diabatic_tend, rho_p, rho_p_save, &
                                                    rho_pp, rho_zz, rho_base, ruAvg, ru_save, ru_p, u, ru, &
                                                    exner, exner_base, rtheta_base, pressure_p,         &
                                                    zz, theta_m, pressure_b
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3, zb_cell, zb3_cell
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
      integer, dimension(:), pointer :: nEdgesOnCell

      integer :: i, iCell, iEdge, k, cell1, cell2
      integer, pointer :: nVertLevels, nCells, nCellsSolve, nEdges, nEdgesSolve
      real (kind=RKIND) :: invNs, rcv, p0, flux
      real (kind=RKIND), pointer :: cf1, cf2, cf3


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
      call mpas_pool_get_array(mesh, 'zb_cell', zb_cell)
      call mpas_pool_get_array(mesh, 'zb3_cell', zb3_cell)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)

      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)


      call atm_recover_large_step_variables_work(nCells, nEdges, nCellsSolve, nEdgesSolve, dt, ns, rk_step, &
                             wwAvg, rw_save, w, rw, rw_p, rtheta_p, rtheta_pp, rtheta_p_save, rt_diabatic_tend, rho_p, &
                             rho_p_save, rho_pp, rho_zz, rho_base, ruAvg, ru_save, ru_p, u, ru, exner, exner_base, &
                             rtheta_base, pressure_p, zz, theta_m, pressure_b, scalars, fzm, fzp, &
                             zb, zb3, zb_cell, zb3_cell, edgesOnCell_sign, cellsOnEdge, edgesOnCell, nEdgesOnCell, &
                             cf1, cf2, cf3, &
                             cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                             cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)

   end subroutine atm_recover_large_step_variables


   subroutine atm_recover_large_step_variables_work(nCells, nEdges, nCellsSolve, nEdgesSolve, dt, ns, rk_step, &
                             wwAvg, rw_save, w, rw, rw_p, rtheta_p, rtheta_pp, rtheta_p_save, rt_diabatic_tend, rho_p, &
                             rho_p_save, rho_pp, rho_zz, rho_base, ruAvg, ru_save, ru_p, u, ru, exner, exner_base, &
                             rtheta_base, pressure_p, zz, theta_m, pressure_b, scalars, fzm, fzp, &
                             zb, zb3, zb_cell, zb3_cell, edgesOnCell_sign, cellsOnEdge, edgesOnCell, nEdgesOnCell, &
                             cf1, cf2, cf3, &
                             cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                             cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)

      use mpas_atm_dimensions

      implicit none


      !
      ! Dummy arguments
      !
      integer, intent(in) :: nCells, nEdges, nCellsSolve, nEdgesSolve
      integer, intent(in) :: ns, rk_step
      real (kind=RKIND), intent(in) :: dt

      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: wwAvg
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw_save
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: w
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw_p
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rtheta_p
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rtheta_pp
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rtheta_p_save
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rt_diabatic_tend
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_p
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_p_save
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_pp
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_zz
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_base
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ruAvg
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru_save
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru_p
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: u
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: exner
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: exner_base
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rtheta_base
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: pressure_p
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: zz
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: theta_m
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: pressure_b
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1) :: scalars
      real (kind=RKIND), dimension(nVertLevels) :: fzm
      real (kind=RKIND), dimension(nVertLevels) :: fzp
      real (kind=RKIND), dimension(nVertLevels+1,2,nEdges+1) :: zb
      real (kind=RKIND), dimension(nVertLevels+1,2,nEdges+1) :: zb3
      real (kind=RKIND), dimension(nVertLevels+1,maxEdges,nCells+1) :: zb_cell
      real (kind=RKIND), dimension(nVertLevels+1,maxEdges,nCells+1) :: zb3_cell
      real (kind=RKIND), dimension(maxEdges,nCells+1) :: edgesOnCell_sign
      integer, dimension(2,nEdges+1) :: cellsOnEdge
      integer, dimension(maxEdges,nCells+1) :: edgesOnCell
      integer, dimension(nCells+1) :: nEdgesOnCell

      real (kind=RKIND) :: cf1, cf2, cf3

      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd


      !
      ! Local variables
      !
      integer :: i, iCell, iEdge, k, cell1, cell2
      real (kind=RKIND) :: invNs, rcv, p0, flux


      rcv = rgas/(cp-rgas)
      p0 = 1.0e+05  ! this should come from somewhere else...

      ! Avoid FP errors caused by a potential division by zero below by 
      ! initializing the "garbage cell" of rho_zz to a non-zero value
      do k=1,nVertLevels
         rho_zz(k,nCells+1) = 1.0
      end do

      ! compute new density everywhere so we can compute u from ru.
      ! we will also need it to compute theta_m below

      invNs = 1 / real(ns,RKIND)

      do iCell=cellStart,cellEnd

!DIR$ IVDEP
         do k = 1, nVertLevels
            rho_p(k,iCell) = rho_p_save(k,iCell) + rho_pp(k,iCell)

            rho_zz(k,iCell) = rho_p(k,iCell) + rho_base(k,iCell)
         end do

         w(1,iCell) = 0.0

!DIR$ IVDEP
         do k = 2, nVertLevels
            wwAvg(k,iCell) = rw_save(k,iCell) + (wwAvg(k,iCell) * invNs)
            rw(k,iCell) = rw_save(k,iCell) + rw_p(k,iCell)

          ! pick up part of diagnosed w from omega - divide by density later
            w(k,iCell) = rw(k,iCell)/(fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))
                                      
         end do

         w(nVertLevels+1,iCell) = 0.0

         if (rk_step == 3) then
!DIR$ IVDEP
            do k = 1, nVertLevels
               rtheta_p(k,iCell) = rtheta_p_save(k,iCell) + rtheta_pp(k,iCell) &
                                 - dt * rho_zz(k,iCell) * rt_diabatic_tend(k,iCell)
               theta_m(k,iCell) = (rtheta_p(k,iCell) + rtheta_base(k,iCell))/rho_zz(k,iCell)
               exner(k,iCell) = (zz(k,iCell)*(rgas/p0)*(rtheta_p(k,iCell)+rtheta_base(k,iCell)))**rcv
               ! pressure_p is perturbation pressure
               pressure_p(k,iCell) = zz(k,iCell) * rgas * (exner(k,iCell)*rtheta_p(k,iCell)+rtheta_base(k,iCell)  &
                                                          * (exner(k,iCell)-exner_base(k,iCell)))
            end do
         else
!DIR$ IVDEP
            do k = 1, nVertLevels
               rtheta_p(k,iCell) = rtheta_p_save(k,iCell) + rtheta_pp(k,iCell)
               theta_m(k,iCell) = (rtheta_p(k,iCell) + rtheta_base(k,iCell))/rho_zz(k,iCell)
            end do
         end if

      end do

      ! recover time-averaged ruAvg on all edges of owned cells (for upcoming scalar transport).  
      ! we solved for these in the acoustic-step loop.  
      ! we will compute ru and u here also, given we are here, even though we only need them on nEdgesSolve

!$OMP BARRIER

      do iEdge=edgeStart,edgeEnd

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

!DIR$ IVDEP
         do k = 1, nVertLevels
            ruAvg(k,iEdge) = ru_save(k,iEdge) + (ruAvg(k,iEdge) * invNs)
            ru(k,iEdge) = ru_save(k,iEdge) + ru_p(k,iEdge)
            u(k,iEdge) = 2.*ru(k,iEdge)/(rho_zz(k,cell1)+rho_zz(k,cell2))
         end do
      end do

!$OMP BARRIER

      do iCell=cellStart,cellEnd

         !  finish recovering w from (rho*omega)_p.  as when we formed (rho*omega)_p from u and w, we need
         !  to use the same flux-divergence operator as is used for the horizontal theta transport
         !  (See Klemp et al 2003).

         do i=1,nEdgesOnCell(iCell)
            iEdge=edgesOnCell(i,iCell)

            flux = (cf1*ru(1,iEdge) + cf2*ru(2,iEdge) + cf3*ru(3,iEdge))
            w(1,iCell) = w(1,iCell) + edgesOnCell_sign(i,iCell) * &
                                   (zb_cell(1,i,iCell) + sign(1.0_RKIND,flux)*zb3_cell(1,i,iCell))*flux

!DIR$ IVDEP
            do k = 2, nVertLevels
               flux = (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
               w(k,iCell) = w(k,iCell) + edgesOnCell_sign(i,iCell) * &
                                    (zb_cell(k,i,iCell)+sign(1.0_RKIND,flux)*zb3_cell(k,i,iCell))*flux
            end do

         end do

         w(1,iCell) = w(1,iCell)/(cf1*rho_zz(1,iCell)+cf2*rho_zz(2,iCell)+cf3*rho_zz(3,iCell))
!DIR$ IVDEP
         do k = 2, nVertLevels
           w(k,iCell) = w(k,iCell)/(fzm(k)*rho_zz(k,iCell)+fzp(k)*rho_zz(k-1,iCell))
         end do

      end do

   end subroutine atm_recover_large_step_variables_work


   subroutine atm_advance_scalars( tend, state, diag, mesh, configs, num_scalars, nCells, nVertLevels, dt, &
                                   cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
                                   horiz_flux_arr, rk_step, config_time_integration_order, advance_density, scalar_tend, rho_zz_int)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Integrate scalar equations - explicit transport plus other tendencies
   !
   !  Wrapper for atm_advance_scalars_work() to de-reference pointers
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (mpas_pool_type), intent(in) :: tend
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(in) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: num_scalars      ! for allocating stack variables
      integer, intent(in) :: nCells           ! for allocating stack variables
      integer, intent(in) :: nVertLevels      ! for allocating stack variables
      integer, intent(in) :: rk_step    !  rk substep we are integrating
      integer, intent(in) :: config_time_integration_order  ! time integration order
      real (kind=RKIND) :: dt
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd
      logical, intent(in), optional :: advance_density
      real (kind=RKIND), dimension(:,:,:), intent(inout), optional :: scalar_tend
      real (kind=RKIND), dimension(:,:), intent(inout), optional :: rho_zz_int

      integer :: i, j, iCell, iAdvCell, iEdge, k, iScalar, cell1, cell2
      real (kind=RKIND), dimension(:), pointer :: invAreaCell
      real (kind=RKIND) :: rho_zz_new_inv

      real (kind=RKIND) :: scalar_weight

      real (kind=RKIND), dimension(:,:,:), pointer :: scalar_old, scalar_new, scalar_tend_save
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: uhAvg, rho_zz_old, rho_zz_new, wwAvg, rho_edge, zgrid, kdiff
      real (kind=RKIND), dimension(:), pointer :: dvEdge, qv_init
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: horiz_flux_arr

      integer, dimension(:,:), pointer :: advCellsForEdge, edgesOnCell
      integer, dimension(:), pointer :: nAdvCellsForEdge, nEdgesOnCell
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd, edgesOnCell_sign

      real (kind=RKIND), dimension( num_scalars, nVertLevels + 1 ) :: wdtn
      integer, pointer :: nCellsSolve, nEdges

      real (kind=RKIND), dimension(:), pointer :: fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4
      real (kind=RKIND), pointer :: coef_3rd_order

      logical :: local_advance_density

      if (present(advance_density)) then
         local_advance_density = advance_density
      else
         local_advance_density = .true.
      end if

      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)

      call mpas_pool_get_array(state, 'scalars', scalar_old, 1)
      call mpas_pool_get_array(state, 'scalars', scalar_new, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_old, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_new, 2)

      call mpas_pool_get_array(diag, 'kdiff', kdiff)
      call mpas_pool_get_array(diag, 'ruAvg', uhAvg)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)

      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      call mpas_pool_get_array(tend, 'scalars_tend', scalar_tend_save)

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

      if (local_advance_density) then
!      call atm_advance_scalars_work(num_scalars, nCells, nVertLevels, dt, &
!             cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
!             cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
!             coef_3rd_order, scalar_old, scalar_new, rho_zz_old, rho_zz_new, kdiff, &
!             uhAvg, wwAvg, deriv_two, dvEdge, &
!             cellsOnEdge, edgesOnCell, nEdgesOnCell, edgesOnCell_sign, invAreaCell, &
!             scalar_tend_save, fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4, &
!             nAdvCellsForEdge, advCellsForEdge, adv_coefs, adv_coefs_3rd, rho_edge, qv_init, zgrid, &
!             nCellsSolve, nEdges, horiz_flux_arr, &
!             local_advance_density, scalar_tend, rho_zz_int)
      call atm_advance_scalars_work_new(num_scalars, nCells, nVertLevels, dt, &
             cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
             cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
             coef_3rd_order, scalar_old, scalar_new, rho_zz_old, rho_zz_new, kdiff, &
             uhAvg, wwAvg, deriv_two, dvEdge, &
             cellsOnEdge, edgesOnCell, nEdgesOnCell, edgesOnCell_sign, invAreaCell, &
             scalar_tend_save, fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4, &
             nAdvCellsForEdge, advCellsForEdge, adv_coefs, adv_coefs_3rd, rho_edge, qv_init, zgrid, &
             nCellsSolve, nEdges, horiz_flux_arr, rk_step, config_time_integration_order, &
             local_advance_density, scalar_tend, rho_zz_int)
      else
!      call atm_advance_scalars_work(num_scalars, nCells, nVertLevels, dt, &
!             cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
!             cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
!             coef_3rd_order, scalar_old, scalar_new, rho_zz_old, rho_zz_new, kdiff, &
!             uhAvg, wwAvg, deriv_two, dvEdge, &
!             cellsOnEdge, edgesOnCell, nEdgesOnCell, edgesOnCell_sign, invAreaCell, &
!             scalar_tend_save, fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4, &
!             nAdvCellsForEdge, advCellsForEdge, adv_coefs, adv_coefs_3rd, rho_edge, qv_init, zgrid, &
!             nCellsSolve, nEdges, horiz_flux_arr, &
!             local_advance_density)
      call atm_advance_scalars_work_new(num_scalars, nCells, nVertLevels, dt, &
             cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
             cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
             coef_3rd_order, scalar_old, scalar_new, rho_zz_old, rho_zz_new, kdiff, &
             uhAvg, wwAvg, deriv_two, dvEdge, &
             cellsOnEdge, edgesOnCell, nEdgesOnCell, edgesOnCell_sign, invAreaCell, &
             scalar_tend_save, fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4, &
             nAdvCellsForEdge, advCellsForEdge, adv_coefs, adv_coefs_3rd, rho_edge, qv_init, zgrid, &
             nCellsSolve, nEdges, horiz_flux_arr, rk_step, config_time_integration_order, &
             local_advance_density)
      end if

   end subroutine atm_advance_scalars


   subroutine atm_advance_scalars_work( num_scalars_dummy, nCells, nVertLevels_dummy, dt, &
             cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
             cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
             coef_3rd_order, scalar_old, scalar_new, rho_zz_old, rho_zz_new, kdiff, &
             uhAvg, wwAvg, deriv_two, dvEdge, &
             cellsOnEdge, edgesOnCell, nEdgesOnCell, edgesOnCell_sign, invAreaCell, &
             scalar_tend_save, fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4, &
             nAdvCellsForEdge, advCellsForEdge, adv_coefs, adv_coefs_3rd, rho_edge, qv_init, zgrid, &
             nCellsSolve, nEdges, horiz_flux_arr, &
             advance_density, scalar_tend, rho_zz_int)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Integrate scalar equations - explicit transport plus other tendencies
   !
   !  this transport routine is similar to the original atm_advance_scalars, except it also advances
   !  (re-integrates) the density.  This re-integration allows the scalar transport routine to use a different 
   !  timestep than the dry dynamics, and also makes possible a spatial splitting of the scalar transport integration
   !  (and density integration).  The current integration is, however, not spatially split.
   !
   !  WCS 18 November 2014
   !-----------------------
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

      use mpas_atm_dimensions

      implicit none

      integer, intent(in) :: num_scalars_dummy ! for allocating stack variables
      integer, intent(in) :: nCells           ! for allocating stack variables
      integer, intent(in) :: nVertLevels_dummy ! for allocating stack variables
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd
      logical, intent(in) :: advance_density
      real (kind=RKIND), dimension(:,:,:), intent(in) :: scalar_old
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: scalar_new
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: scalar_tend_save
      real (kind=RKIND), dimension(:,:,:), intent(in) :: deriv_two
      real (kind=RKIND), dimension(:,:), intent(in) :: rho_zz_old
      real (kind=RKIND), dimension(:,:), intent(in) :: uhAvg, wwAvg, rho_edge, zgrid, rho_zz_new, kdiff
      real (kind=RKIND), dimension(:), intent(in) :: dvEdge, qv_init
      integer, dimension(:,:), intent(in) :: cellsOnEdge
      integer, dimension(:,:), intent(in) :: advCellsForEdge, edgesOnCell
      integer, dimension(:), intent(in) :: nAdvCellsForEdge, nEdgesOnCell
      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs, adv_coefs_3rd, edgesOnCell_sign
      real (kind=RKIND), dimension(:), intent(in) :: fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4
      real (kind=RKIND), intent(in) :: coef_3rd_order
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nEdges+1), intent(inout) :: horiz_flux_arr
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(inout), optional :: scalar_tend
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout), optional :: rho_zz_int
      real (kind=RKIND), dimension(:), intent(in) :: invAreaCell
      integer, intent(in) :: nCellsSolve, nEdges

      integer :: i, j, iCell, iAdvCell, iEdge, k, iScalar, cell1, cell2
      real (kind=RKIND) :: rho_zz_new_inv

      real (kind=RKIND) :: scalar_weight

      real (kind=RKIND), dimension( num_scalars, nVertLevels + 1 ) :: wdtn

      real (kind=RKIND), dimension(nVertLevels,10) :: scalar_weight2
      integer, dimension(10) :: ica

      real (kind=RKIND) :: flux3, flux4
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3

      logical :: local_advance_density

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
          ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0

      local_advance_density = advance_density

      !
      ! Runge Kutta integration, so we compute fluxes from scalar_new values, update starts from scalar_old
      !
      !  horizontal flux divergence, accumulate in scalar_tend


      !  horiz_flux_arr stores the value of the scalar at the edge.
      !  a better name perhaps would be scalarEdge

      do iEdge=edgeStart,edgeEnd

         select case(nAdvCellsForEdge(iEdge))

         case(10)

            do j=1,10
!DIR$ IVDEP
               do k=1,nVertLevels
                  scalar_weight2(k,j) = adv_coefs(j,iEdge) + sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(j,iEdge)
               end do
            end do
            do j=1,10
               ica(j) = advCellsForEdge(j,iEdge)
            end do
!DIR$ IVDEP
            do k = 1,nVertLevels
!DIR$ IVDEP
               do iScalar = 1,num_scalars
                  horiz_flux_arr(iScalar,k,iEdge) = &
                       scalar_weight2(k,1)  * scalar_new(iScalar,k,ica(1)) + &
                       scalar_weight2(k,2)  * scalar_new(iScalar,k,ica(2)) + &
                       scalar_weight2(k,3)  * scalar_new(iScalar,k,ica(3)) + &
                       scalar_weight2(k,4)  * scalar_new(iScalar,k,ica(4)) + &
                       scalar_weight2(k,5)  * scalar_new(iScalar,k,ica(5)) + &
                       scalar_weight2(k,6)  * scalar_new(iScalar,k,ica(6)) + &
                       scalar_weight2(k,7)  * scalar_new(iScalar,k,ica(7)) + &
                       scalar_weight2(k,8)  * scalar_new(iScalar,k,ica(8)) + &
                       scalar_weight2(k,9)  * scalar_new(iScalar,k,ica(9)) + &
                       scalar_weight2(k,10) * scalar_new(iScalar,k,ica(10))
               end do
            end do

         case default

            horiz_flux_arr(:,:,iEdge) = 0.0
            do j=1,nAdvCellsForEdge(iEdge)
               iAdvCell = advCellsForEdge(j,iEdge)
!DIR$ IVDEP
               do k=1,nVertLevels
                  scalar_weight = adv_coefs(j,iEdge) + sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(j,iEdge)
!DIR$ IVDEP
                  do iScalar=1,num_scalars
                     horiz_flux_arr(iScalar,k,iEdge) = horiz_flux_arr(iScalar,k,iEdge) + scalar_weight * scalar_new(iScalar,k,iAdvCell)
                  end do
               end do
            end do

         end select
      end do

!$OMP BARRIER

      if (local_advance_density) then
         if ((.not.present(scalar_tend)) .or. (.not.present(rho_zz_int))) then
            call mpas_dmpar_global_abort('Error: rho_zz_int or scalar_tend not supplied to atm_advance_scalars( ) when advance_density=.true.')
         end if

         do iCell=cellSolveStart,cellSolveEnd
            scalar_tend(:,:,iCell) = scalar_tend_save(:,:,iCell)
#ifndef DO_PHYSICS
            scalar_tend(:,:,iCell) = 0.0  !  testing purposes - we have no sources or sinks
#endif

            rho_zz_int(:,iCell) = 0.0

            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
   
               ! here we add the horizontal flux divergence into the scalar tendency.
               ! note that the scalar tendency is modified.
!DIR$ IVDEP
               do k=1,nVertLevels
                  rho_zz_int(k,iCell) = rho_zz_int(k,iCell) - edgesOnCell_sign(i,iCell) * uhAvg(k,iEdge)*dvEdge(iEdge) * invAreaCell(iCell)
!DIR$ IVDEP
                  do iScalar=1,num_scalars
                        scalar_tend(iScalar,k,iCell) = scalar_tend(iScalar,k,iCell) &
                               - edgesOnCell_sign(i,iCell) * uhAvg(k,iEdge)*horiz_flux_arr(iScalar,k,iEdge) * invAreaCell(iCell)
                  end do
               end do
               
            end do

!DIR$ IVDEP
            do k=1,nVertLevels
               rho_zz_int(k,iCell) = rho_zz_old(k,iCell) + dt*( rho_zz_int(k,iCell) - rdnw(k)*(wwAvg(k+1,iCell)-wwAvg(k,iCell)) )
            end do
         end do

      else

         do iCell=cellSolveStart,cellSolveEnd
#ifndef DO_PHYSICS
            scalar_tend_save(:,:,iCell) = 0.0  !  testing purposes - we have no sources or sinks
#endif

            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)

               ! here we add the horizontal flux divergence into the scalar tendency.
               ! note that the scalar tendency is modified.
!DIR$ IVDEP
               do k=1,nVertLevels
!DIR$ IVDEP
                  do iScalar=1,num_scalars
                        scalar_tend_save(iScalar,k,iCell) = scalar_tend_save(iScalar,k,iCell) &
                               - edgesOnCell_sign(i,iCell) * uhAvg(k,iEdge)*horiz_flux_arr(iScalar,k,iEdge) * invAreaCell(iCell)
                  end do
               end do
               
            end do
         end do

      end if

      !
      !  vertical flux divergence and update of the scalars
      !

      ! zero fluxes at top and bottom
      wdtn(:,1) = 0.0
      wdtn(:,nVertLevels+1) = 0.0


      do iCell=cellSolveStart,cellSolveEnd

         k = 2
         do iScalar=1,num_scalars
            wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
         end do
          
!DIR$ IVDEP
         do k=3,nVertLevels-1
!DIR$ IVDEP
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

      if (local_advance_density) then
!DIR$ IVDEP
         do k=1,nVertLevels
            rho_zz_new_inv = 1.0_RKIND / rho_zz_int(k,iCell)
!DIR$ IVDEP
            do iScalar=1,num_scalars
               scalar_new(iScalar,k,iCell) = (   scalar_old(iScalar,k,iCell)*rho_zz_old(k,iCell) &
                     + dt*( scalar_tend(iScalar,k,iCell) -rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k)) ) ) * rho_zz_new_inv
            end do
         end do
      else
!DIR$ IVDEP
         do k=1,nVertLevels
            rho_zz_new_inv = 1.0_RKIND / rho_zz_new(k,iCell)
!DIR$ IVDEP
            do iScalar=1,num_scalars
               scalar_new(iScalar,k,iCell) = (   scalar_old(iScalar,k,iCell)*rho_zz_old(k,iCell) &
                     + dt*( scalar_tend_save(iScalar,k,iCell) -rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k)) ) ) * rho_zz_new_inv
            end do
         end do
      end if

      end do

   end subroutine atm_advance_scalars_work


   subroutine atm_advance_scalars_work_new( num_scalars_dummy, nCells, nVertLevels_dummy, dt, &
             cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
             cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
             coef_3rd_order, scalar_old, scalar_new, rho_zz_old, rho_zz_new, kdiff, &
             uhAvg, wwAvg, deriv_two, dvEdge, &
             cellsOnEdge, edgesOnCell, nEdgesOnCell, edgesOnCell_sign, invAreaCell, &
             scalar_tend_save, fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4, &
             nAdvCellsForEdge, advCellsForEdge, adv_coefs, adv_coefs_3rd, rho_edge, qv_init, zgrid, &
             nCellsSolve, nEdges, horiz_flux_arr, rk_step, config_time_integration_order, &
             advance_density, scalar_tend, rho_zz_int)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Integrate scalar equations - explicit transport plus other tendencies
   !
   !  this transport routine is similar to the original atm_advance_scalars, except it also advances
   !  (re-integrates) the density.  This re-integration allows the scalar transport routine to use a different 
   !  timestep than the dry dynamics, and also makes possible a spatial splitting of the scalar transport integration
   !  (and density integration).  The current integration is, however, not spatially split.
   !
   !  WCS 18 November 2014
   !-----------------------
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

      use mpas_atm_dimensions

      implicit none

      integer, intent(in) :: num_scalars_dummy ! for allocating stack variables
      integer, intent(in) :: nCells           ! for allocating stack variables
      integer, intent(in) :: nVertLevels_dummy ! for allocating stack variables
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd
      integer, intent(in) :: rk_step, config_time_integration_order
      logical, intent(in) :: advance_density
      real (kind=RKIND), dimension(:,:,:), intent(in) :: scalar_old
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: scalar_new
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: scalar_tend_save
      real (kind=RKIND), dimension(:,:,:), intent(in) :: deriv_two
      real (kind=RKIND), dimension(:,:), intent(in) :: rho_zz_old
      real (kind=RKIND), dimension(:,:), intent(in) :: uhAvg, wwAvg, rho_edge, zgrid, rho_zz_new, kdiff
      real (kind=RKIND), dimension(:), intent(in) :: dvEdge, qv_init
      integer, dimension(:,:), intent(in) :: cellsOnEdge
      integer, dimension(:,:), intent(in) :: advCellsForEdge, edgesOnCell
      integer, dimension(:), intent(in) :: nAdvCellsForEdge, nEdgesOnCell
      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs, adv_coefs_3rd, edgesOnCell_sign
      real (kind=RKIND), dimension(:), intent(in) :: fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4
      real (kind=RKIND), intent(in) :: coef_3rd_order
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nEdges+1), intent(inout) :: horiz_flux_arr
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(inout), optional :: scalar_tend
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout), optional :: rho_zz_int
      real (kind=RKIND), dimension(:), intent(in) :: invAreaCell
      integer, intent(in) :: nCellsSolve, nEdges

      integer :: i, j, iCell, iAdvCell, iEdge, k, iScalar, cell1, cell2
      real (kind=RKIND) :: rho_zz_new_inv

      real (kind=RKIND) :: scalar_weight

      real (kind=RKIND), dimension( num_scalars, nVertLevels + 1 ) :: wdtn

      real (kind=RKIND), dimension(nVertLevels,10) :: scalar_weight2
      integer, dimension(10) :: ica

      real (kind=RKIND) :: flux3, flux4
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3

      logical :: local_advance_density

      real (kind=RKIND) :: weight_time_old, weight_time_new
      real (kind=RKIND), dimension(num_scalars,nVertLevels) :: scalar_tend_column  ! local storage to accumulate tendency

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
          ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0

      local_advance_density = advance_density

      !
      ! Runge Kutta integration, so we compute fluxes from scalar_new values, update starts from scalar_old
      !
      !  horizontal flux divergence, accumulate in scalar_tend


      !  horiz_flux_arr stores the value of the scalar at the edge.
      !  a better name perhaps would be scalarEdge

      !  weights for the time interpolation of the input density
      !
      if (.not. advance_density ) then
         weight_time_new = 1.
      else
         if((rk_step == 1) .and. config_time_integration_order == 3) weight_time_new = 1./3
         if((rk_step == 1) .and. config_time_integration_order == 2) weight_time_new = 1./2
         if(rk_step == 2) weight_time_new = 1./2
         if(rk_step == 3) weight_time_new = 1.
      end if
      weight_time_old = 1. - weight_time_new


      do iEdge=edgeStart,edgeEnd

         select case(nAdvCellsForEdge(iEdge))

         case(10)

            do j=1,10
!DIR$ IVDEP
               do k=1,nVertLevels
                  scalar_weight2(k,j) = adv_coefs(j,iEdge) + sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(j,iEdge)
               end do
            end do
            do j=1,10
               ica(j) = advCellsForEdge(j,iEdge)
            end do
!DIR$ IVDEP
            do k = 1,nVertLevels
!DIR$ IVDEP
               do iScalar = 1,num_scalars
                  horiz_flux_arr(iScalar,k,iEdge) = &
                       scalar_weight2(k,1)  * scalar_new(iScalar,k,ica(1)) + &
                       scalar_weight2(k,2)  * scalar_new(iScalar,k,ica(2)) + &
                       scalar_weight2(k,3)  * scalar_new(iScalar,k,ica(3)) + &
                       scalar_weight2(k,4)  * scalar_new(iScalar,k,ica(4)) + &
                       scalar_weight2(k,5)  * scalar_new(iScalar,k,ica(5)) + &
                       scalar_weight2(k,6)  * scalar_new(iScalar,k,ica(6)) + &
                       scalar_weight2(k,7)  * scalar_new(iScalar,k,ica(7)) + &
                       scalar_weight2(k,8)  * scalar_new(iScalar,k,ica(8)) + &
                       scalar_weight2(k,9)  * scalar_new(iScalar,k,ica(9)) + &
                       scalar_weight2(k,10) * scalar_new(iScalar,k,ica(10))
               end do
            end do

         case default

            horiz_flux_arr(:,:,iEdge) = 0.0
            do j=1,nAdvCellsForEdge(iEdge)
               iAdvCell = advCellsForEdge(j,iEdge)
!DIR$ IVDEP
               do k=1,nVertLevels
                  scalar_weight = adv_coefs(j,iEdge) + sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(j,iEdge)
!DIR$ IVDEP
                  do iScalar=1,num_scalars
                     horiz_flux_arr(iScalar,k,iEdge) = horiz_flux_arr(iScalar,k,iEdge) + scalar_weight * scalar_new(iScalar,k,iAdvCell)
                  end do
               end do
            end do

         end select
      end do

!$OMP BARRIER

!  scalar update, for each column sum fluxes over horizontal edges, add physics tendency, and add vertical flux divergence in update.

      
      do iCell=cellSolveStart,cellSolveEnd
#ifndef DO_PHYSICS
            scalar_tend_save(:,:,iCell) = 0.0  !  testing purposes - we have no sources or sinks
#endif
            scalar_tend_column(1:num_scalars,1:nVertlevels) = 0.

            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)

               ! here we add the horizontal flux divergence into the scalar tendency.
               ! note that the scalar tendency is modified.
!DIR$ IVDEP
               do k=1,nVertLevels
!DIR$ IVDEP
                  do iScalar=1,num_scalars
                        scalar_tend_column(iScalar,k) = scalar_tend_column(iScalar,k) &
                               - edgesOnCell_sign(i,iCell) * uhAvg(k,iEdge)*horiz_flux_arr(iScalar,k,iEdge)
                  end do
               end do
               
            end do

!DIR$ IVDEP
            do k=1,nVertLevels
!DIR$ IVDEP
               do iScalar=1,num_scalars
                     scalar_tend_column(iScalar,k) = scalar_tend_column(iScalar,k) * invAreaCell(iCell) + scalar_tend_save(iScalar,k,iCell)
               end do
            end do
 

      !
      !  vertical flux divergence and update of the scalars
      !
         wdtn(:,1) = 0.0
         wdtn(:,nVertLevels+1) = 0.0

         k = 2
         do iScalar=1,num_scalars
            wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
         end do
          
!DIR$ IVDEP
         do k=3,nVertLevels-1
!DIR$ IVDEP
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

!DIR$ IVDEP
         do k=1,nVertLevels
            rho_zz_new_inv = 1.0_RKIND / (weight_time_old*rho_zz_old(k,iCell) + weight_time_new*rho_zz_new(k,iCell))
!DIR$ IVDEP
            do iScalar=1,num_scalars
               scalar_new(iScalar,k,iCell) = (   scalar_old(iScalar,k,iCell)*rho_zz_old(k,iCell) &
                     + dt*( scalar_tend_column(iScalar,k) -rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k)) ) ) * rho_zz_new_inv
            end do
         end do

      end do

   end subroutine atm_advance_scalars_work_new


   subroutine atm_advance_scalars_mono(block, tend, state, diag, mesh, configs, nCells, nEdges, nVertLevels_dummy, dt, &
                                   cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
                                   scalar_old, scalar_new, s_max, s_min, wdtn, scale_arr, flux_arr, &
                                   flux_upwind_tmp, flux_tmp, advance_density, rho_zz_int)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Integrate scalar equations - transport plus other tendencies
   !
   !  wrapper routine for atm_advance_scalars_mono_work
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      use mpas_atm_dimensions

      implicit none

      type (block_type), intent(inout), target :: block
      type (mpas_pool_type), intent(in)    :: tend
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(in)    :: diag
      type (mpas_pool_type), intent(in)    :: mesh
      type (mpas_pool_type), intent(in)    :: configs
      integer, intent(in)                  :: nCells           ! for allocating stack variables
      integer, intent(in)                  :: nEdges           ! for allocating stack variables
      integer, intent(in)                  :: nVertLevels_dummy      ! for allocating stack variables
      real (kind=RKIND), intent(in)        :: dt
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: scalar_old, scalar_new
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: s_max, s_min
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(inout) :: wdtn
      real (kind=RKIND), dimension(nVertLevels,2,nCells+1), intent(inout) :: scale_arr
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: flux_arr
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: flux_upwind_tmp, flux_tmp
      logical, intent(in), optional :: advance_density
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout), optional :: rho_zz_int

      real (kind=RKIND), dimension(:,:,:), pointer :: scalar_tend
      real (kind=RKIND), dimension(:,:), pointer :: uhAvg, rho_zz_old, rho_zz_new, wwAvg
      real (kind=RKIND), dimension(:), pointer :: dvEdge, invAreaCell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, edgesOnCell
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign

      integer, dimension(:,:), pointer :: advCellsForEdge
      integer, dimension(:), pointer :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars_old, scalars_new

      integer, pointer :: nCellsSolve

      real (kind=RKIND), dimension(:), pointer :: fnm, fnp, rdnw
      integer, dimension(:), pointer :: nEdgesOnCell
      real (kind=RKIND), pointer :: coef_3rd_order

      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)

      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(diag, 'ruAvg', uhAvg)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)

      call mpas_pool_get_array(tend, 'scalars_tend', scalar_tend)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz_old, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_new, 2)
      call mpas_pool_get_array(state, 'scalars', scalars_old, 1)
      call mpas_pool_get_array(state, 'scalars', scalars_new, 2)

      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'fzm', fnm)
      call mpas_pool_get_array(mesh, 'fzp', fnp)
      call mpas_pool_get_array(mesh, 'rdzw', rdnw)
      call mpas_pool_get_array(mesh, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(mesh, 'advCellsForEdge', advCellsForEdge)
      call mpas_pool_get_array(mesh, 'adv_coefs', adv_coefs)
      call mpas_pool_get_array(mesh, 'adv_coefs_3rd', adv_coefs_3rd)

      call atm_advance_scalars_mono_work(block, state, nCells, nEdges, nVertLevels, dt, &
                                   cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
                                   coef_3rd_order, nCellsSolve, num_scalars, uhAvg, wwAvg, scalar_tend, rho_zz_old, &
                                   rho_zz_new, scalars_old, scalars_new, invAreaCell, dvEdge, cellsOnEdge, cellsOnCell, &
                                   edgesOnCell, edgesOnCell_sign, nEdgesOnCell, fnm, fnp, rdnw, nAdvCellsForEdge, &
                                   advCellsForEdge, adv_coefs, adv_coefs_3rd, scalar_old, scalar_new, s_max, s_min, &
                                   wdtn, scale_arr, flux_arr, flux_upwind_tmp, flux_tmp, &
                                   advance_density, rho_zz_int)

   end subroutine atm_advance_scalars_mono


   subroutine atm_advance_scalars_mono_work(block, state, nCells, nEdges, nVertLevels_dummy, dt, &
                                   cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd, &
                                   coef_3rd_order, nCellsSolve, num_scalars_dummy, uhAvg, wwAvg, scalar_tend, rho_zz_old, &
                                   rho_zz_new, scalars_old, scalars_new, invAreaCell, dvEdge, cellsOnEdge, cellsOnCell, &
                                   edgesOnCell, edgesOnCell_sign, nEdgesOnCell, fnm, fnp, rdnw, nAdvCellsForEdge, &
                                   advCellsForEdge, adv_coefs, adv_coefs_3rd, scalar_old, scalar_new, s_max, s_min, &
                                   wdtn, scale_arr, flux_arr, flux_upwind_tmp, flux_tmp, &
                                   advance_density, rho_zz_int)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Integrate scalar equations - transport plus other tendencies
   !
   !  this transport routine is similar to the original atm_advance_scalars_mono_work, except it also advances
   !  (re-integrates) the density.  This re-integration allows the scalar transport routine to use a different 
   !  timestep than the dry dynamics, and also makes possible a spatial splitting of the scalar transport integration
   !  (and density integration).  The current integration is, however, not spatially split.
   !
   !  WCS 18 November 2014
   !-----------------------
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

      use mpas_atm_dimensions

      implicit none

      type (block_type), intent(inout), target :: block
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in)                  :: nCells           ! for allocating stack variables
      integer, intent(in)                  :: nEdges           ! for allocating stack variables
      integer, intent(in)                  :: nVertLevels_dummy ! for allocating stack variables
      real (kind=RKIND), intent(in)        :: dt
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd
      logical, intent(in), optional :: advance_density
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout), optional :: rho_zz_int

      integer :: ii,jj
      integer, dimension(10) :: ica
      real (kind=RKIND), dimension(10,2) :: swa

      integer :: i, iCell, iEdge, k, iScalar, cell1, cell2
      real (kind=RKIND) :: flux, scalar_weight

      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: scalar_tend
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: uhAvg
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rho_zz_old, rho_zz_new
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(in) :: wwAvg
      real (kind=RKIND), dimension(:), intent(in) :: dvEdge, invAreaCell
      integer, dimension(:,:), intent(in) :: cellsOnEdge, cellsOnCell, edgesOnCell
      real (kind=RKIND), dimension(:,:), intent(in) :: edgesOnCell_sign

      integer, dimension(:,:), intent(in) :: advCellsForEdge
      integer, dimension(:), intent(in) :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs, adv_coefs_3rd
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: scalars_old, scalars_new
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: scalar_old, scalar_new
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: s_max, s_min
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(inout) :: wdtn
      real (kind=RKIND), dimension(nVertLevels,2,nCells+1), intent(inout), target :: scale_arr
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: flux_arr
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: flux_upwind_tmp, flux_tmp
      type (field3DReal), pointer :: scalars_old_field

      type (field3DReal), pointer :: tempField
      type (field3DReal), target :: tempFieldTarget


      integer, parameter :: SCALE_IN = 1, SCALE_OUT = 2

      integer, intent(in) :: nCellsSolve, num_scalars_dummy
      integer :: icellmax, kmax

      real (kind=RKIND), dimension(nVertLevels), intent(in) :: fnm, fnp, rdnw
      integer, dimension(:), intent(in) :: nEdgesOnCell
      real (kind=RKIND), intent(in) :: coef_3rd_order


      real (kind=RKIND), dimension(nVertLevels) :: flux_upwind_arr
      real (kind=RKIND) :: flux3, flux4, flux_upwind
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3, scmin,scmax
      real (kind=RKIND) :: scale_factor

      logical :: local_advance_density

      real (kind=RKIND), parameter :: eps=1.e-20

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
          ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0

      if (present(advance_density)) then
         local_advance_density = advance_density
      else
         local_advance_density = .true.
      end if

      call mpas_pool_get_field(state, 'scalars', scalars_old_field, 1)

      !  for positive-definite or monotonic option, we first update scalars using the tendency from sources other than
      !  the resolved transport (these should constitute a positive definite update).  
      !  Note, however, that we enforce positive-definiteness in this update.
      !  The transport will maintain this positive definite solution and optionally, shape preservation (monotonicity).

      do iCell=cellSolveStart,cellSolveEnd
!DIR$ IVDEP
      do k = 1,nVertLevels
!DIR$ IVDEP
      do iScalar = 1,num_scalars

#ifndef DO_PHYSICS
!TBH:  Michael, would you please check this change?  Our test uses -DDO_PHYSICS
!TBH:  so this code is not executed.  The change avoids redundant work.  
         scalar_tend(iScalar,k,iCell) = 0.0  !  testing purposes - we have no sources or sinks
#endif
         scalars_old(iScalar,k,iCell) = scalars_old(iScalar,k,iCell)+dt*scalar_tend(iScalar,k,iCell) / rho_zz_old(k,iCell)
         scalar_tend(iScalar,k,iCell) = 0.0
      end do
      end do
      end do


!$OMP BARRIER
!$OMP MASTER
MPAS_GPTL_START('mono_exch_1')
      call mpas_dmpar_exch_halo_field(scalars_old_field)
MPAS_GPTL_STOP('mono_exch_1')
!$OMP END MASTER
!$OMP BARRIER

      !
      ! Runge Kutta integration, so we compute fluxes from scalar_new values, update starts from scalar_old
      !

      if (local_advance_density) then
         if (.not.present(rho_zz_int)) then
            call mpas_dmpar_global_abort('Error: rho_zz_int not supplied to atm_advance_scalars_mono_work( ) when advance_density=.true.')
         end if

         !  begin with update of density
         do iCell=cellStart,cellEnd
            rho_zz_int(:,iCell) = 0.0
         end do
!$OMP BARRIER
         do iCell=cellSolveStart,cellSolveEnd
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
   
!DIR$ IVDEP
               do k=1,nVertLevels
                  rho_zz_int(k,iCell) = rho_zz_int(k,iCell) - edgesOnCell_sign(i,iCell) * uhAvg(k,iEdge) * dvEdge(iEdge) * invAreaCell(iCell)
               end do
               
            end do
         end do
         do iCell=cellSolveStart,cellSolveEnd
!DIR$ IVDEP
            do k=1,nVertLevels
               rho_zz_int(k,iCell) = rho_zz_old(k,iCell) + dt*( rho_zz_int(k,iCell) - rdnw(k)*(wwAvg(k+1,iCell)-wwAvg(k,iCell)) )
            end do
         end do
!$OMP BARRIER
      end if

      ! next, do one scalar at a time

      do iScalar = 1, num_scalars

         do iCell=cellStart,cellEnd
!DIR$ IVDEP
            do k=1,nVertLevels
               scalar_old(k,iCell) = scalars_old(iScalar,k,iCell)
               scalar_new(k,iCell) = scalars_new(iScalar,k,iCell)
            end do
         end do

!$OMP BARRIER

#ifdef DEBUG_TRANSPORT
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
#endif


         !
         !  vertical flux divergence, and min and max bounds for flux limiter
         !
         do iCell=cellSolveStart,cellSolveEnd

            ! zero flux at top and bottom
            wdtn(1,iCell) = 0.0
            wdtn(nVertLevels+1,iCell) = 0.0

            k = 1
            s_max(k,iCell) = max(scalar_old(1,iCell),scalar_old(2,iCell))
            s_min(k,iCell) = min(scalar_old(1,iCell),scalar_old(2,iCell))

            k = 2
            wdtn(k,iCell) = wwAvg(k,iCell)*(fnm(k)*scalar_new(k,iCell)+fnp(k)*scalar_new(k-1,iCell))
            s_max(k,iCell) = max(scalar_old(k-1,iCell),scalar_old(k,iCell),scalar_old(k+1,iCell))
            s_min(k,iCell) = min(scalar_old(k-1,iCell),scalar_old(k,iCell),scalar_old(k+1,iCell))
             
!DIR$ IVDEP
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

            !
            ! pull s_min and s_max from the (horizontal) surrounding cells
            !

            ! speclal treatment of calculations involving hexagonal cells
            ! original code retained in select "default" case
            select case(nEdgesOnCell(iCell))
            case(6)
!DIR$ IVDEP
               do k=1, nVertLevels
                  s_max(k,iCell) = max(s_max(k,iCell), &
                       scalar_old(k, cellsOnCell(1,iCell)), &
                       scalar_old(k, cellsOnCell(2,iCell)), &
                       scalar_old(k, cellsOnCell(3,iCell)), &
                       scalar_old(k, cellsOnCell(4,iCell)), &
                       scalar_old(k, cellsOnCell(5,iCell)), &
                       scalar_old(k, cellsOnCell(6,iCell)))
                  s_min(k,iCell) = min(s_min(k,iCell), &
                       scalar_old(k, cellsOnCell(1,iCell)), &
                       scalar_old(k, cellsOnCell(2,iCell)), &
                       scalar_old(k, cellsOnCell(3,iCell)), &
                       scalar_old(k, cellsOnCell(4,iCell)), &
                       scalar_old(k, cellsOnCell(5,iCell)), &
                       scalar_old(k, cellsOnCell(6,iCell)))
               enddo

            case default
               do i=1, nEdgesOnCell(iCell)
!DIR$ IVDEP
                  do k=1, nVertLevels
                     s_max(k,iCell) = max(s_max(k,iCell),scalar_old(k, cellsOnCell(i,iCell)))
                     s_min(k,iCell) = min(s_min(k,iCell),scalar_old(k, cellsOnCell(i,iCell)))
                  end do
               end do
            end select

         end do

!$OMP BARRIER

         !
         !  horizontal flux divergence
         !

         do iEdge=edgeStart,edgeEnd

            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then  ! only for owned cells
  
               ! speclal treatment of calculations involving edges between hexagonal cells
               ! original code retained in select "default" case
               ! be sure to see additional declarations near top of subroutine
               select case(nAdvCellsForEdge(iEdge))
               case(10)
                  do jj=1,10
                     ica(jj)    = advCellsForEdge(jj,iEdge)
                     swa(jj,1)  = adv_coefs(jj,iEdge) + adv_coefs_3rd(jj,iEdge)
                     swa(jj,2)  = adv_coefs(jj,iEdge) - adv_coefs_3rd(jj,iEdge)
                  enddo
                  do k=1,nVertLevels
                     ii = merge(1, 2, uhAvg(k,iEdge) > 0)
                     flux_arr(k,iEdge) = uhAvg(k,iEdge)*( &
                          swa(1,ii)*scalar_new(k,ica(1)) + swa(2,ii)*scalar_new(k,ica(2)) + &
                          swa(3,ii)*scalar_new(k,ica(3)) + swa(4,ii)*scalar_new(k,ica(4)) + &
                          swa(5,ii)*scalar_new(k,ica(5)) + swa(6,ii)*scalar_new(k,ica(6)) + &
                          swa(7,ii)*scalar_new(k,ica(7)) + swa(8,ii)*scalar_new(k,ica(8)) + &
                          swa(9,ii)*scalar_new(k,ica(9)) + swa(10,ii)*scalar_new(k,ica(10)))
                  enddo

               case default
                  do k=1,nVertLevels
                     flux_arr(k,iEdge) = 0.0_RKIND
                  enddo
                  do i=1,nAdvCellsForEdge(iEdge)
                     iCell = advCellsForEdge(i,iEdge)
!DIR$ IVDEP
                     do k=1,nVertLevels
                        scalar_weight = uhAvg(k,iEdge)*(adv_coefs(i,iEdge) + sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(i,iEdge))
                        flux_arr(k,iEdge) = flux_arr(k,iEdge) + scalar_weight* scalar_new(k,iCell)
                     end do
                  end do
               end select

            end if

         end do

!$OMP BARRIER

         !
         !  vertical flux divergence for upwind update, we will put upwind update into scalar_new, and put factor of dt in fluxes
         !

         do iCell=cellSolveStart,cellSolveEnd

            k = 1
            scalar_new(k,iCell) = scalar_old(k,iCell) * rho_zz_old(k,iCell)

!DIR$ IVDEP
            do k = 2, nVertLevels
               scalar_new(k,iCell) = scalar_old(k,iCell)*rho_zz_old(k,iCell)
               flux_upwind_arr(k) = dt*(max(0.0_RKIND,wwAvg(k,iCell))*scalar_old(k-1,iCell) + min(0.0_RKIND,wwAvg(k,iCell))*scalar_old(k,iCell))
            end do
            do k = 1, nVertLevels-1
               scalar_new(k,iCell) = scalar_new(k,iCell) - flux_upwind_arr(k+1)*rdnw(k)
            end do
!DIR$ IVDEP
            do k = 2, nVertLevels
               scalar_new(k  ,iCell) = scalar_new(k  ,iCell) + flux_upwind_arr(k)*rdnw(k)
               wdtn(k,iCell) = dt*wdtn(k,iCell) - flux_upwind_arr(k)
            end do


            !
            ! scale_arr(SCALE_IN,:,:) and scale_arr(SCALE_OUT:,:) are used here to store the incoming and outgoing perturbation flux 
            ! contributions to the update:  first the vertical flux component, then the horizontal
            !

!DIR$ IVDEP
            do k=1,nVertLevels
               scale_arr(k,SCALE_IN, iCell) = - rdnw(k)*(min(0.0_RKIND,wdtn(k+1,iCell))-max(0.0_RKIND,wdtn(k,iCell)))
               scale_arr(k,SCALE_OUT,iCell) = - rdnw(k)*(max(0.0_RKIND,wdtn(k+1,iCell))-min(0.0_RKIND,wdtn(k,iCell)))
            end do

         end do

         !
         !  horizontal flux divergence for upwind update
         !

         !  upwind flux computation
         do iEdge=edgeStart,edgeEnd
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
            do k=1, nVertLevels
               flux_upwind_tmp(k,iEdge) = dvEdge(iEdge) * dt *   &
                      (max(0.0_RKIND,uhAvg(k,iEdge))*scalar_old(k,cell1) + min(0.0_RKIND,uhAvg(k,iEdge))*scalar_old(k,cell2))
               flux_tmp(k,iEdge) = dt * flux_arr(k,iEdge) - flux_upwind_tmp(k,iEdge)
            end do
         end do
!$OMP BARRIER
         do iCell=cellSolveStart,cellSolveEnd
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)

!DIR$ IVDEP
               do k=1, nVertLevels
                  scalar_new(k,iCell) = scalar_new(k,iCell) - edgesOnCell_sign(i,iCell) * flux_upwind_tmp(k,iEdge) * invAreaCell(iCell)
 
                  scale_arr(k,SCALE_OUT,iCell) = scale_arr(k,SCALE_OUT,iCell) &
                                                 - max(0.0_RKIND,edgesOnCell_sign(i,iCell)*flux_tmp(k,iEdge)) * invAreaCell(iCell)
                  scale_arr(k,SCALE_IN, iCell) = scale_arr(k,SCALE_IN, iCell) &
                                                 - min(0.0_RKIND,edgesOnCell_sign(i,iCell)*flux_tmp(k,iEdge)) * invAreaCell(iCell)
               end do

            end do
         end do

         !
         !  next, the limiter
         !

         ! simplification of limiter calculations
         ! worked through algebra and found equivalent form
         ! added benefit that it should address ifort single prec overflow issue
      if (local_advance_density) then
         do iCell=cellSolveStart,cellSolveEnd
!DIR$ IVDEP
            do k = 1, nVertLevels

               scale_factor = (s_max(k,iCell)*rho_zz_int(k,iCell) - scalar_new(k,iCell)) / &
                    (scale_arr(k,SCALE_IN,iCell)  + eps)
               scale_arr(k,SCALE_IN,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )

               scale_factor = (s_min(k,iCell)*rho_zz_int(k,iCell) - scalar_new(k,iCell)) / &
                    (scale_arr(k,SCALE_OUT,iCell) - eps)
               scale_arr(k,SCALE_OUT,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
            end do
         end do
      else
         do iCell=cellSolveStart,cellSolveEnd
!DIR$ IVDEP
            do k = 1, nVertLevels

               scale_factor = (s_max(k,iCell)*rho_zz_new(k,iCell) - scalar_new(k,iCell)) / &
                    (scale_arr(k,SCALE_IN,iCell)  + eps)
               scale_arr(k,SCALE_IN,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )

               scale_factor = (s_min(k,iCell)*rho_zz_new(k,iCell) - scalar_new(k,iCell)) / &
                    (scale_arr(k,SCALE_OUT,iCell) - eps)
               scale_arr(k,SCALE_OUT,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
            end do
         end do
      end if

         !
         !  communicate scale factors here.
         !  communicate only first halo row in these next two exchanges
         !
!$OMP BARRIER
!$OMP MASTER
MPAS_GPTL_START('mono_exch_2')
         tempField => tempFieldTarget

         tempField % block => block
         tempField % dimSizes(1) = nVertLevels
         tempField % dimSizes(2) = 2
         tempField % dimSizes(3) = nCells
         tempField % sendList => block % parinfo % cellsToSend
         tempField % recvList => block % parinfo % cellsToRecv
         tempField % copyList => block % parinfo % cellsToCopy
         tempField % prev => null()
         tempField % next => null()
         tempField % isActive = .true.

         tempField % array => scale_arr
         call mpas_dmpar_exch_halo_field(tempField, (/ 1 /))
MPAS_GPTL_STOP('mono_exch_2')
!$OMP END MASTER
!$OMP BARRIER

         do iEdge=edgeStart,edgeEnd
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then  ! only for owned cells
!DIR$ IVDEP
               do k=1, nVertLevels
                  flux_upwind = dvEdge(iEdge) * dt *   &
                         (max(0.0_RKIND,uhAvg(k,iEdge))*scalar_old(k,cell1) + min(0.0_RKIND,uhAvg(k,iEdge))*scalar_old(k,cell2))
                  flux_arr(k,iEdge) = dt*flux_arr(k,iEdge) - flux_upwind
               end do
            end if
         end do

         !
         !  rescale the fluxes
         !

         ! moved assignment to scalar_new from separate loop (see commented code below)
         ! into the following loops. Avoids having to save elements of flux array
         do iEdge=edgeStart,edgeEnd
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
!DIR$ IVDEP
               do k = 1, nVertLevels
                  flux = flux_arr(k,iEdge)
                  flux = max(0.0_RKIND,flux) * min(scale_arr(k,SCALE_OUT,cell1), scale_arr(k,SCALE_IN, cell2)) &
                       + min(0.0_RKIND,flux) * min(scale_arr(k,SCALE_IN, cell1), scale_arr(k,SCALE_OUT,cell2))
                  flux_arr(k,iEdge) = flux
               end do
            end if
         end do
 
        !
        ! rescale the vertical flux
        !
!$OMP BARRIER
 
         do iCell=cellSolveStart,cellSolveEnd
!DIR$ IVDEP
            do k = 2, nVertLevels
               flux = wdtn(k,iCell)
               flux = max(0.0_RKIND,flux) * min(scale_arr(k-1,SCALE_OUT,iCell), scale_arr(k  ,SCALE_IN,iCell)) &
                    + min(0.0_RKIND,flux) * min(scale_arr(k  ,SCALE_OUT,iCell), scale_arr(k-1,SCALE_IN,iCell))
               wdtn(k,iCell) = flux
            end do
         end do


         !
         !  do the scalar update now that we have the fluxes
         !
         do iCell=cellSolveStart,cellSolveEnd
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
!DIR$ IVDEP
               do k=1,nVertLevels
                  scalar_new(k,iCell) = scalar_new(k,iCell) - edgesOnCell_sign(i,iCell)*flux_arr(k,iEdge) * invAreaCell(iCell)
               end do
            end do

      if (local_advance_density) then
!DIR$ IVDEP
            do k=1,nVertLevels
               scalar_new(k,iCell) = (   scalar_new(k,iCell) + (-rdnw(k)*(wdtn(k+1,iCell)-wdtn(k,iCell)) ) )/rho_zz_int(k,iCell)
            end do
      else
!DIR$ IVDEP
            do k=1,nVertLevels
               scalar_new(k,iCell) = (   scalar_new(k,iCell) + (-rdnw(k)*(wdtn(k+1,iCell)-wdtn(k,iCell)) ) )/rho_zz_new(k,iCell)
            end do
      end if
         end do

#ifdef DEBUG_TRANSPORT
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
#endif

         ! the update should be positive definite. but roundoff can sometimes leave small negative values
         ! hence the enforcement of PD in the copy back to the model state.
!$OMP BARRIER

         do iCell=cellStart,cellEnd
            do k=1, nVertLevels
               scalars_new(iScalar,k,iCell) = max(0.0_RKIND,scalar_new(k,iCell))
            end do
         end do

      end do !  loop over scalars

   end subroutine atm_advance_scalars_mono_work


   subroutine atm_compute_dyn_tend(tend, tend_physics, state, diag, mesh, configs, nVertLevels, rk_step, dt, &
                                   cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)
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

      !
      ! Dummy arguments
      !
      type (mpas_pool_type), intent(inout) :: tend
      type (mpas_pool_type), intent(inout) :: tend_physics
      type (mpas_pool_type), intent(in) :: state
      type (mpas_pool_type), intent(in) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: nVertLevels              ! for allocating stack variables
      integer, intent(in) :: rk_step
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd

      !
      ! Local variables
      !
      integer, pointer :: nCells, nEdges, nVertices, nCellsSolve, nEdgesSolve, vertexDegree, maxEdges, maxEdges2
      integer, pointer :: moist_start, moist_end, num_scalars

      real (kind=RKIND), dimension(:), pointer ::  fEdge, dvEdge, dcEdge, invDcEdge, invDvEdge, invAreaCell, invAreaTriangle, &
                                                   meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge, zgrid, rho_edge, rho_zz, ru, u, v, tend_u, &
                                                    divergence, vorticity, ke, pv_edge, theta_m, rw, tend_rho, &
                                                    rt_diabatic_tend, tend_theta, tend_w, w, cqw, rb, rr, pp, pressure_b, zz, zxu, cqu, & 
                                                    h_divergence, kdiff, edgesOnCell_sign, edgesOnVertex_sign, rw_save, ru_save

      real (kind=RKIND), dimension(:,:), pointer :: theta_m_save

      real (kind=RKIND), dimension(:,:), pointer :: exner

      real (kind=RKIND), dimension(:,:), pointer :: rr_save


      real (kind=RKIND), dimension(:,:), pointer :: tend_rtheta_adv  ! needed for Tiedtke convection scheme
      real (kind=RKIND), dimension(:,:), pointer :: rthdynten  ! needed for Grell-Freitas convection scheme

      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      real (kind=RKIND), dimension(:,:), pointer :: tend_u_euler, tend_w_euler, tend_theta_euler

      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgesOnCell, edgesOnEdge, cellsOnCell, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, latEdge, angleEdge, u_init

      integer, dimension(:,:), pointer :: advCellsForEdge
      integer, dimension(:), pointer :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd

      real (kind=RKIND), dimension(:), pointer :: rdzu, rdzw, fzm, fzp, qv_init
      real (kind=RKIND), dimension(:,:), pointer :: t_init 

      real (kind=RKIND), pointer :: cf1, cf2, cf3

      real (kind=RKIND), pointer :: r_earth
      real (kind=RKIND), dimension(:,:), pointer :: ur_cell, vr_cell

      real (kind=RKIND), dimension(:,:), pointer :: defc_a, defc_b

      real(kind=RKIND), dimension(:,:), pointer :: tend_w_pgf, tend_w_buoy

      real (kind=RKIND), pointer :: coef_3rd_order, c_s, smdiv
      logical, pointer :: config_mix_full
      character (len=StrKIND), pointer :: config_horiz_mixing
      real (kind=RKIND), pointer :: config_del4u_div_factor
      real (kind=RKIND), pointer :: config_h_theta_eddy_visc4
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc4
      real (kind=RKIND), pointer :: config_visc4_2dsmag
      real (kind=RKIND), pointer :: config_len_disp
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc2, config_v_mom_eddy_visc2
      real (kind=RKIND), pointer :: config_h_theta_eddy_visc2, config_v_theta_eddy_visc2

      logical :: inactive_rthdynten


      call mpas_pool_get_config(mesh, 'sphere_radius', r_earth)
      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_mix_full', config_mix_full)
      call mpas_pool_get_config(configs, 'config_horiz_mixing', config_horiz_mixing)
      call mpas_pool_get_config(configs, 'config_del4u_div_factor', config_del4u_div_factor)
      call mpas_pool_get_config(configs, 'config_h_theta_eddy_visc4', config_h_theta_eddy_visc4)
      call mpas_pool_get_config(configs, 'config_h_mom_eddy_visc4', config_h_mom_eddy_visc4)
      call mpas_pool_get_config(configs, 'config_h_theta_eddy_visc2', config_h_theta_eddy_visc2)
      call mpas_pool_get_config(configs, 'config_h_mom_eddy_visc2', config_h_mom_eddy_visc2)
      call mpas_pool_get_config(configs, 'config_v_theta_eddy_visc2', config_v_theta_eddy_visc2)
      call mpas_pool_get_config(configs, 'config_v_mom_eddy_visc2', config_v_mom_eddy_visc2)
      call mpas_pool_get_config(configs, 'config_visc4_2dsmag', config_visc4_2dsmag)
      call mpas_pool_get_config(configs, 'config_len_disp', config_len_disp)
      call mpas_pool_get_config(configs, 'config_smagorinsky_coef', c_s)
      call mpas_pool_get_config(configs, 'config_smdiv', smdiv) 

      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(state, 'u', u, 2)
      call mpas_pool_get_array(state, 'w', w, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m_save, 1)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)

      call mpas_pool_get_array(diag, 'uReconstructZonal', ur_cell)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', vr_cell)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)
      call mpas_pool_get_array(diag, 'rho_base', rb)
      call mpas_pool_get_array(diag, 'rho_p', rr)
      call mpas_pool_get_array(diag, 'rho_p_save', rr_save)
      call mpas_pool_get_array(diag, 'v', v)
      call mpas_pool_get_array(diag, 'kdiff', kdiff)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'ru_save', ru_save)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'rw_save', rw_save)
      call mpas_pool_get_array(diag, 'divergence', divergence)
      call mpas_pool_get_array(diag, 'vorticity', vorticity)
      call mpas_pool_get_array(diag, 'ke', ke)
      call mpas_pool_get_array(diag, 'pv_edge', pv_edge)
      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
      call mpas_pool_get_array(diag, 'h_divergence', h_divergence)
      call mpas_pool_get_array(diag, 'exner', exner)

      call mpas_pool_get_array(diag, 'tend_rtheta_adv', tend_rtheta_adv)
      call mpas_pool_get_array(tend_physics, 'rthdynten', rthdynten)

      call mpas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(mesh, 'edgesOnVertex_sign', edgesOnVertex_sign)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'invDcEdge', invDcEdge)
      call mpas_pool_get_array(mesh, 'invDvEdge', invDvEdge)
      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      call mpas_pool_get_array(mesh, 'invAreaTriangle', invAreaTriangle)
      call mpas_pool_get_array(mesh, 'fEdge', fEdge)
      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zxu', zxu)
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
      call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)
      call mpas_pool_get_dimension(mesh, 'maxEdges', maxEdges)
      call mpas_pool_get_dimension(mesh, 'maxEdges2', maxEdges2)

      call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
      call mpas_pool_get_dimension(state, 'moist_start', moist_start)
      call mpas_pool_get_dimension(state, 'moist_end', moist_end)

      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(mesh, 'advCellsForEdge', advCellsForEdge)
      call mpas_pool_get_array(mesh, 'adv_coefs', adv_coefs)
      call mpas_pool_get_array(mesh, 'adv_coefs_3rd', adv_coefs_3rd)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)

      !
      ! rthdynten is currently associated with packages, and if those packages
      ! are not active at run-time, we need to produce an rthdynten array for
      ! use in the atm_compute_dyn_tend_work routine
      !
      inactive_rthdynten = .false.
      if (.not. associated(rthdynten)) then
         allocate(rthdynten(nVertLevels,nCells+1))
         inactive_rthdynten = .true.
      end if

      call atm_compute_dyn_tend_work(nCells, nEdges, nVertices, nVertLevels, &
         nCellsSolve, nEdgesSolve, vertexDegree, maxEdges, maxEdges2, num_scalars, moist_start, moist_end, &
         fEdge, dvEdge, dcEdge, invDcEdge, invDvEdge, invAreaCell, invAreaTriangle, meshScalingDel2, meshScalingDel4, &
         weightsOnEdge, zgrid, rho_edge, rho_zz, ru, u, v, tend_u, &
         divergence, vorticity, ke, pv_edge, theta_m, rw, tend_rho, &
         rt_diabatic_tend, tend_theta, tend_w, w, cqw, rb, rr, pp, pressure_b, zz, zxu, cqu, & 
         h_divergence, kdiff, edgesOnCell_sign, edgesOnVertex_sign, rw_save, ru_save, &
         theta_m_save, exner, rr_save, scalars, tend_u_euler, tend_w_euler, tend_theta_euler, deriv_two, &
         cellsOnEdge, verticesOnEdge, edgesOnCell, edgesOnEdge, cellsOnCell, edgesOnVertex, nEdgesOnCell, nEdgesOnEdge, &
         latCell, latEdge, angleEdge, u_init, advCellsForEdge, nAdvCellsForEdge, adv_coefs, adv_coefs_3rd, &
         rdzu, rdzw, fzm, fzp, qv_init, t_init, cf1, cf2, cf3, r_earth, ur_cell, vr_cell, defc_a, defc_b, &
         tend_w_pgf, tend_w_buoy, coef_3rd_order, c_s, smdiv, config_mix_full, config_horiz_mixing, config_del4u_div_factor, &
         config_h_mom_eddy_visc2, config_v_mom_eddy_visc2, config_h_theta_eddy_visc2, config_v_theta_eddy_visc2, &
         config_h_theta_eddy_visc4, config_h_mom_eddy_visc4, config_visc4_2dsmag, config_len_disp, rk_step, dt, &
         tend_rtheta_adv, rthdynten, &
         cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
         cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)

      if (inactive_rthdynten) then
         deallocate(rthdynten)
      end if

   end subroutine atm_compute_dyn_tend


   subroutine atm_compute_dyn_tend_work(nCells, nEdges, nVertices, nVertLevels_dummy, &
      nCellsSolve, nEdgesSolve, vertexDegree, maxEdges_dummy, maxEdges2_dummy, num_scalars_dummy, moist_start, moist_end, &
      fEdge, dvEdge, dcEdge, invDcEdge, invDvEdge, invAreaCell, invAreaTriangle, meshScalingDel2, meshScalingDel4, &
      weightsOnEdge, zgrid, rho_edge, rho_zz, ru, u, v, tend_u, &
      divergence, vorticity, ke, pv_edge, theta_m, rw, tend_rho, &
      rt_diabatic_tend, tend_theta, tend_w, w, cqw, rb, rr, pp, pressure_b, zz, zxu, cqu, & 
      h_divergence, kdiff, edgesOnCell_sign, edgesOnVertex_sign, rw_save, ru_save, &
      theta_m_save, exner, rr_save, scalars, tend_u_euler, tend_w_euler, tend_theta_euler, deriv_two, &
      cellsOnEdge, verticesOnEdge, edgesOnCell, edgesOnEdge, cellsOnCell, edgesOnVertex, nEdgesOnCell, nEdgesOnEdge, &
      latCell, latEdge, angleEdge, u_init, advCellsForEdge, nAdvCellsForEdge, adv_coefs, adv_coefs_3rd, &
      rdzu, rdzw, fzm, fzp, qv_init, t_init, cf1, cf2, cf3, r_earth, ur_cell, vr_cell, defc_a, defc_b, &
      tend_w_pgf, tend_w_buoy, coef_3rd_order, c_s, smdiv, config_mix_full, config_horiz_mixing, config_del4u_div_factor, &
      config_h_mom_eddy_visc2, config_v_mom_eddy_visc2, config_h_theta_eddy_visc2, config_v_theta_eddy_visc2, &
      config_h_theta_eddy_visc4, config_h_mom_eddy_visc4, config_visc4_2dsmag, config_len_disp, rk_step, dt, &
      tend_rtheta_adv, rthdynten, &
      cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
      cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)


      use mpas_atm_dimensions


      implicit none


      !
      ! Dummy arguments
      !
      integer :: nCells, nEdges, nVertices, nVertLevels_dummy, nCellsSolve, nEdgesSolve, vertexDegree, &
                 maxEdges_dummy, maxEdges2_dummy, num_scalars_dummy, moist_start, moist_end

      real (kind=RKIND), dimension(nEdges+1) :: fEdge
      real (kind=RKIND), dimension(nEdges+1) :: dvEdge
      real (kind=RKIND), dimension(nEdges+1) :: dcEdge
      real (kind=RKIND), dimension(nEdges+1) :: invDcEdge
      real (kind=RKIND), dimension(nEdges+1) :: invDvEdge
      real (kind=RKIND), dimension(nCells+1) :: invAreaCell
      real (kind=RKIND), dimension(nVertices+1) :: invAreaTriangle
      real (kind=RKIND), dimension(nEdges+1) :: meshScalingDel2
      real (kind=RKIND), dimension(nEdges+1) :: meshScalingDel4
      real (kind=RKIND), dimension(maxEdges2,nEdges+1) :: weightsOnEdge
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: zgrid
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: rho_edge
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rho_zz
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: u
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: v
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: tend_u
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: divergence
      real (kind=RKIND), dimension(nVertLevels,nVertices+1) :: vorticity
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: ke
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: pv_edge
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: theta_m
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: tend_rho
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rt_diabatic_tend
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: tend_theta
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: tend_w
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: w
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: cqw
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rb
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rr
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: pp
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: pressure_b
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: zz
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: zxu
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: cqu
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: h_divergence
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: kdiff
      real (kind=RKIND), dimension(maxEdges,nCells+1) :: edgesOnCell_sign
      real (kind=RKIND), dimension(vertexDegree,nVertices+1) :: edgesOnVertex_sign
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: rw_save
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: ru_save

      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: theta_m_save
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: exner
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rr_save
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1) :: scalars
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: tend_u_euler
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: tend_w_euler
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: tend_theta_euler
      real (kind=RKIND), dimension(15,2,nEdges+1) :: deriv_two
      integer, dimension(2,nEdges+1) :: cellsOnEdge
      integer, dimension(2,nEdges+1) :: verticesOnEdge
      integer, dimension(maxEdges,nCells+1) :: edgesOnCell
      integer, dimension(maxEdges2,nEdges+1) :: edgesOnEdge
      integer, dimension(maxEdges,nCells+1) :: cellsOnCell
      integer, dimension(vertexDegree,nVertices+1) :: edgesOnVertex
      integer, dimension(nCells+1) :: nEdgesOnCell
      integer, dimension(nEdges+1) :: nEdgesOnEdge
      real (kind=RKIND), dimension(nCells+1) :: latCell
      real (kind=RKIND), dimension(nEdges+1) :: latEdge
      real (kind=RKIND), dimension(nEdges+1) :: angleEdge
      real (kind=RKIND), dimension(nVertLevels) :: u_init

      integer, dimension(15,nEdges+1) :: advCellsForEdge
      integer, dimension(nEdges+1) :: nAdvCellsForEdge
      real (kind=RKIND), dimension(15,nEdges+1) :: adv_coefs
      real (kind=RKIND), dimension(15,nEdges+1) :: adv_coefs_3rd

      real (kind=RKIND), dimension(nVertLevels) :: rdzu
      real (kind=RKIND), dimension(nVertLevels) :: rdzw
      real (kind=RKIND), dimension(nVertLevels) :: fzm
      real (kind=RKIND), dimension(nVertLevels) :: fzp
      real (kind=RKIND), dimension(nVertLevels) :: qv_init
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: t_init 

      real (kind=RKIND) :: cf1, cf2, cf3
      real (kind=RKIND) :: prandtl_inv, r_areaCell, rgas_cprcv

      real (kind=RKIND) :: r_earth
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: ur_cell
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: vr_cell

      real (kind=RKIND), dimension(maxEdges,nCells+1) :: defc_a
      real (kind=RKIND), dimension(maxEdges,nCells+1) :: defc_b

      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: tend_w_pgf
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1) :: tend_w_buoy

      real (kind=RKIND) :: coef_3rd_order, c_s, smdiv
      logical :: config_mix_full
      character (len=StrKIND) :: config_horiz_mixing
      real (kind=RKIND) :: config_del4u_div_factor
      real (kind=RKIND) :: config_h_theta_eddy_visc4
      real (kind=RKIND) :: config_h_mom_eddy_visc4
      real (kind=RKIND) :: config_visc4_2dsmag
      real (kind=RKIND) :: config_len_disp
      real (kind=RKIND) :: config_h_mom_eddy_visc2, config_v_mom_eddy_visc2, config_h_theta_eddy_visc2, config_v_theta_eddy_visc2

      integer, intent(in) :: rk_step
      real (kind=RKIND), intent(in) :: dt

      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: tend_rtheta_adv
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: rthdynten

      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd


      !
      ! Local variables
      !
      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, iq, iAdvCell

      !real (kind=RKIND), parameter :: c_s = 0.125
      real (kind=RKIND), dimension( nVertLevels+1 ) :: d_diag, d_off_diag, flux_arr
      real (kind=RKIND), dimension( nVertLevels + 1 ) :: wduz, wdwz, wdtz, dpzx
      real (kind=RKIND), dimension( nVertLevels ) :: ru_edge_w, q, u_mix
      real (kind=RKIND) :: theta_turb_flux, w_turb_flux, r
      real (kind=RKIND) :: scalar_weight
      real (kind=RKIND) :: inv_r_earth

      real (kind=RKIND) :: invDt, flux, workpv
      real (kind=RKIND) :: edge_sign, pr_scale, r_dc, r_dv, u_mix_scale
      real (kind=RKIND) :: h_mom_eddy_visc4, v_mom_eddy_visc2
      real (kind=RKIND) :: h_theta_eddy_visc4, v_theta_eddy_visc2
      real (kind=RKIND) :: u_diffusion

      real (kind=RKIND) :: kdiffu, z1, z2, z3, z4, zm, z0, zp

      

      real (kind=RKIND) :: flux3, flux4
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
                ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0


      prandtl_inv = 1.0_RKIND / prandtl
      invDt = 1.0_RKIND / dt
      inv_r_earth = 1.0_RKIND / r_earth

       v_mom_eddy_visc2   = config_v_mom_eddy_visc2
       v_theta_eddy_visc2 = config_v_theta_eddy_visc2

      if (rk_step == 1) then

!         tend_u_euler(1:nVertLevels,edgeStart:edgeEnd) = 0.0

         ! Smagorinsky eddy viscosity, based on horizontal deformation (in this case on model coordinate surfaces).
         ! The integration coefficients were precomputed and stored in defc_a and defc_b

         if(config_horiz_mixing == "2d_smagorinsky") then
            do iCell = cellStart,cellEnd
               d_diag(1:nVertLevels) = 0.0
               d_off_diag(1:nVertLevels) = 0.0
               do iEdge=1,nEdgesOnCell(iCell)
                  do k=1,nVertLevels
                     d_diag(k)     = d_diag(k)     + defc_a(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell))  &
                                                   - defc_b(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
                     d_off_diag(k) = d_off_diag(k) + defc_b(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell))  &
                                                   + defc_a(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
                  end do
               end do
!DIR$ IVDEP
               do k=1, nVertLevels
                  ! here is the Smagorinsky formulation, 
                  ! followed by imposition of an upper bound on the eddy viscosity
                  kdiff(k,iCell) = min((c_s * config_len_disp)**2 * sqrt(d_diag(k)**2 + d_off_diag(k)**2),(0.01*config_len_disp**2) * invDt)
               end do
            end do

            h_mom_eddy_visc4   = config_visc4_2dsmag * config_len_disp**3
            h_theta_eddy_visc4 = h_mom_eddy_visc4

         else if(config_horiz_mixing == "2d_fixed") then

            kdiff(1:nVertLevels,cellStart:cellEnd) = config_h_theta_eddy_visc2
            h_mom_eddy_visc4 = config_h_mom_eddy_visc4
            h_theta_eddy_visc4 = config_h_theta_eddy_visc4

         end if
            
      end if

      ! tendency for density.
      ! accumulate total water here for later use in w tendency calculation.

      ! accumulate horizontal mass-flux

      do iCell=cellStart,cellEnd
         h_divergence(1:nVertLevels,iCell) = 0.0
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            edge_sign = edgesOnCell_sign(i,iCell) * dvEdge(iEdge)
!DIR$ IVDEP
            do k=1,nVertLevels
               h_divergence(k,iCell) = h_divergence(k,iCell) + edge_sign * ru(k,iEdge)
            end do
         end do
      end do

      ! compute horiontal mass-flux divergence, add vertical mass flux divergence to complete tend_rho

      do iCell = cellStart,cellEnd
         r = invAreaCell(iCell)
         do k = 1,nVertLevels
            h_divergence(k,iCell) = h_divergence(k,iCell) * r
         end do
      end do    

      !
      ! dp / dz and tend_rho
      !
      ! only needed on first rk_step with pert variables defined a pert from time t
      !
      if(rk_step == 1) then

        rgas_cprcv = rgas*cp/cv
        do iCell = cellStart,cellEnd

!DIR$ IVDEP
          do k = 1,nVertLevels
            tend_rho(k,iCell) = -h_divergence(k,iCell)-rdzw(k)*(rw(k+1,iCell)-rw(k,iCell))
            dpdz(k,iCell) = -gravity*(rb(k,iCell)*(qtot(k,iCell)) + rr_save(k,iCell)*(1.+qtot(k,iCell)))
          end do
        end do
      end if

!$OMP BARRIER

      !
      ! Compute u (normal) velocity tendency for each edge (cell face)
      !

      do iEdge=edgeSolveStart,edgeSolveEnd

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! horizontal pressure gradient 

         if(rk_step == 1) then
!DIR$ IVDEP
            do k=1,nVertLevels
               tend_u_euler(k,iEdge) =  - cqu(k,iEdge)*( (pp(k,cell2)-pp(k,cell1))*invDcEdge(iEdge)/(.5*(zz(k,cell2)+zz(k,cell1))) &
                                              -0.5*zxu(k,iEdge)*(dpdz(k,cell1)+dpdz(k,cell2)) )
            end do

         end if

         ! vertical transport of u

         wduz(1) = 0.

         k = 2
         wduz(k) =  0.5*( rw(k,cell1)+rw(k,cell2))*(fzm(k)*u(k,iEdge)+fzp(k)*u(k-1,iEdge))
         do k=3,nVertLevels-1
            wduz(k) = flux3( u(k-2,iEdge),u(k-1,iEdge),u(k,iEdge),u(k+1,iEdge),0.5*(rw(k,cell1)+rw(k,cell2)), 1.0_RKIND )
         end do
         k = nVertLevels
         wduz(k) =  0.5*( rw(k,cell1)+rw(k,cell2))*(fzm(k)*u(k,iEdge)+fzp(k)*u(k-1,iEdge))

         wduz(nVertLevels+1) = 0.

!DIR$ IVDEP
         do k=1,nVertLevels
            tend_u(k,iEdge) = - rdzw(k)*(wduz(k+1)-wduz(k)) !  first use of tend_u
         end do

         ! Next, nonlinear Coriolis term (q) following Ringler et al JCP 2009

         q(:) = 0.0
         do j = 1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(j,iEdge)
            do k=1,nVertLevels
               workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe))
!  the original definition of pv_edge had a factor of 1/density.  We have removed that factor
!  given that it was not integral to any conservation property of the system
               q(k) = q(k) + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv
            end do
         end do

!DIR$ IVDEP
         do k=1,nVertLevels

            ! horizontal ke gradient and vorticity terms in the vector invariant formulation
            ! of the horizontal momentum equation
            tend_u(k,iEdge) = tend_u(k,iEdge) + rho_edge(k,iEdge)* (q(k) - (ke(k,cell2) - ke(k,cell1))       &
                                                                 * invDcEdge(iEdge))                            &
                                             - u(k,iEdge)*0.5*(h_divergence(k,cell1)+h_divergence(k,cell2)) 
#ifdef CURVATURE
            ! curvature terms for the sphere
            tend_u(k,iEdge) = tend_u(k,iEdge) &
                             - 2.*omega*cos(angleEdge(iEdge))*cos(latEdge(iEdge))  &
                               *rho_edge(k,iEdge)*.25*(w(k,cell1)+w(k+1,cell1)+w(k,cell2)+w(k+1,cell2))          & 
                             - u(k,iEdge)*.25*(w(k+1,cell1)+w(k,cell1)+w(k,cell2)+w(k+1,cell2))                  &
                               *rho_edge(k,iEdge) * inv_r_earth
#endif
         end do

      end do


      !
      !  horizontal mixing for u
      !  mixing terms are integrated using forward-Euler, so this tendency is only computed in the
      !  first Runge-Kutta substep and saved for use in later RK substeps 2 and 3.
      !

      if (rk_step == 1) then

!$OMP BARRIER

         ! del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
         ! First, storage to hold the result from the first del^2 computation.

         delsq_u(1:nVertLevels,edgeStart:edgeEnd) = 0.0

         do iEdge=edgeStart,edgeEnd
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)
            r_dc = invDcEdge(iEdge)
            r_dv = min(invDvEdge(iEdge), 4*invDcEdge(iEdge))

!DIR$ IVDEP
            do k=1,nVertLevels

               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc4 == constant
              u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) * r_dc  &
                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) * r_dv

               delsq_u(k,iEdge) = delsq_u(k,iEdge) + u_diffusion

               kdiffu = 0.5*(kdiff(k,cell1)+kdiff(k,cell2))

               ! include 2nd-orer diffusion here 
               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) &
                                       + rho_edge(k,iEdge)* kdiffu * u_diffusion * meshScalingDel2(iEdge)

            end do
         end do

         if (h_mom_eddy_visc4 > 0.0) then  ! 4th order mixing is active

!$OMP BARRIER

            do iVertex=vertexStart,vertexEnd
               delsq_vorticity(1:nVertLevels,iVertex) = 0.0
               do i=1,vertexDegree
                  iEdge = edgesOnVertex(i,iVertex)
                  edge_sign = invAreaTriangle(iVertex) * dcEdge(iEdge) * edgesOnVertex_sign(i,iVertex)
                  do k=1,nVertLevels
                     delsq_vorticity(k,iVertex) = delsq_vorticity(k,iVertex) + edge_sign * delsq_u(k,iEdge)
                  end do
               end do
            end do

            do iCell=cellStart,cellEnd
               delsq_divergence(1:nVertLevels,iCell) = 0.0
               r = invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i,iCell)
                  edge_sign = r * dvEdge(iEdge) * edgesOnCell_sign(i,iCell)
                  do k=1,nVertLevels
                     delsq_divergence(k,iCell) = delsq_divergence(k,iCell) + edge_sign * delsq_u(k,iEdge)
                  end do
               end do
            end do

!$OMP BARRIER

            do iEdge=edgeSolveStart,edgeSolveEnd
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               vertex1 = verticesOnEdge(1,iEdge)
               vertex2 = verticesOnEdge(2,iEdge)

               u_mix_scale = meshScalingDel4(iEdge)*h_mom_eddy_visc4
               r_dc = u_mix_scale * config_del4u_div_factor * invDcEdge(iEdge)
               r_dv = u_mix_scale * min(invDvEdge(iEdge), 4*invDcEdge(iEdge))

!DIR$ IVDEP
               do k=1,nVertLevels

                  ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
                  !                    only valid for h_mom_eddy_visc4 == constant
                  !
                  ! Here, we scale the diffusion on the divergence part a factor of config_del4u_div_factor 
                  !    relative to the rotational part.  The stability constraint on the divergence component is much less
                  !    stringent than the rotational part, and this flexibility may be useful.
                  !
                  u_diffusion =  rho_edge(k,iEdge) *  ( ( delsq_divergence(k,cell2)  - delsq_divergence(k,cell1) ) * r_dc  &
                                                       -( delsq_vorticity(k,vertex2) - delsq_vorticity(k,vertex1) ) * r_dv )
                  tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) - u_diffusion
                  
               end do
            end do
         
         end if ! 4th order mixing is active 

      !
      !  vertical mixing for u - 2nd order filter in physical (z) space
      !
         if ( v_mom_eddy_visc2 > 0.0 ) then

            if (config_mix_full) then  ! mix full state

               do iEdge=edgeSolveStart,edgeSolveEnd

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

               do iEdge=edgeSolveStart,edgeSolveEnd

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

            end if  ! mix perturbation state

         end if  ! vertical mixing of horizontal momentum

      end if ! (rk_step 1 test for computing mixing terms)

!$OMP BARRIER

!  add in mixing for u

      do iEdge=edgeSolveStart,edgeSolveEnd
!DIR$ IVDEP
         do k=1,nVertLevels
!            tend_u(k,iEdge) = tend_u(k,iEdge) + tend_u_euler(k,iEdge)
            tend_u(k,iEdge) = tend_u(k,iEdge) + tend_u_euler(k,iEdge) + tend_ru_physics(k,iEdge)
         end do
      end do


!----------- rhs for w


      !
      !  horizontal advection for w
      !

      do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
         tend_w(1:nVertLevels+1,iCell) = 0.0
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            edge_sign = edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * 0.5

            do k=2,nVertLevels
               ru_edge_w(k) = fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge)
            end do

            flux_arr(1:nVertLevels) = 0.0

            ! flux_arr stores the value of w at the cell edge used in the horizontal transport

            do j=1,nAdvCellsForEdge(iEdge)
               iAdvCell = advCellsForEdge(j,iEdge)
               do k=2,nVertLevels
                  scalar_weight = adv_coefs(j,iEdge) + sign(1.0_RKIND,ru_edge_w(k)) * adv_coefs_3rd(j,iEdge)
                  flux_arr(k) = flux_arr(k) + scalar_weight * w(k,iAdvCell)
               end do
            end do

!DIR$ IVDEP
            do k=2,nVertLevels
               tend_w(k,iCell) = tend_w(k,iCell) - edgesOnCell_sign(i,iCell) * ru_edge_w(k)*flux_arr(k)
            end do

         end do
      end do

#ifdef CURVATURE
      do iCell = cellSolveStart, cellSolveEnd
!DIR$ IVDEP
         do k=2,nVertLevels
            tend_w(k,iCell) = tend_w(k,iCell) + (rho_zz(k,iCell)*fzm(k)+rho_zz(k-1,iCell)*fzp(k))*          &
                                      ( (fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))**2.             &
                                       +(fzm(k)*vr_cell(k,iCell)+fzp(k)*vr_cell(k-1,iCell))**2. )/r_earth   &
                                + 2.*omega*cos(latCell(iCell))                                              &
                                       *(fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))                 &
                                       *(rho_zz(k,iCell)*fzm(k)+rho_zz(k-1,iCell)*fzp(k))

         end do
      end do
#endif


      !
      !  horizontal mixing for w - we could combine this with advection directly (i.e. as a turbulent flux),
      !  but here we can also code in hyperdiffusion if we wish (2nd order at present)
      !

      if (rk_step == 1) then

!  !OMP BARRIER  why is this openmp barrier here???

         ! del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
         !
         ! First, storage to hold the result from the first del^2 computation.
         !  we copied code from the theta mixing, hence the theta* names.


         delsq_w(1:nVertLevels,cellStart:cellEnd) = 0.0

         do iCell=cellStart,cellEnd
            tend_w_euler(1:nVertLevels+1,iCell) = 0.0
            r_areaCell = invAreaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)

               edge_sign = 0.5 * r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * invDcEdge(iEdge)

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

!DIR$ IVDEP
              do k=2,nVertLevels

                  w_turb_flux =  edge_sign*(rho_edge(k,iEdge)+rho_edge(k-1,iEdge))*(w(k,cell2) - w(k,cell1))
                  delsq_w(k,iCell) = delsq_w(k,iCell) + w_turb_flux
                  w_turb_flux = w_turb_flux * meshScalingDel2(iEdge) * 0.25 * &
                                  (kdiff(k,cell1)+kdiff(k,cell2)+kdiff(k-1,cell1)+kdiff(k-1,cell2))
                  tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + w_turb_flux
               end do
            end do
         end do

!$OMP BARRIER

         if (h_mom_eddy_visc4 > 0.0) then  ! 4th order mixing is active

            do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
               r_areaCell = h_mom_eddy_visc4 * invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i,iCell)
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell) * invDcEdge(iEdge)

                  do k=2,nVertLevels
                     tend_w_euler(k,iCell) = tend_w_euler(k,iCell) - edge_sign * (delsq_w(k,cell2) - delsq_w(k,cell1))
                  end do
           
               end do
            end do

         end if ! 4th order mixing is active 

      end if ! horizontal mixing for w computed in first rk_step

! Note for OpenMP parallelization: We could avoid allocating the delsq_w scratch
!   array, and just use the delsq_theta array as was previously done; however,
!   particularly when oversubscribing cores with threads, there is the risk that
!   some threads may reach code further below that re-uses the delsq_theta array, 
!   in which case we would need a barrier somewhere between here and that code 
!   below to ensure correct behavior.

      !
      !  vertical advection, pressure gradient and buoyancy for w
      !

      do iCell=cellSolveStart,cellSolveEnd

         wdwz(1) = 0.0

         k = 2
         wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))
         do k=3,nVertLevels-1
            wdwz(k) = flux3( w(k-2,iCell),w(k-1,iCell),w(k,iCell),w(k+1,iCell),0.5*(rw(k,iCell)+rw(k-1,iCell)), 1.0_RKIND )
         end do
         k = nVertLevels
         wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))

         wdwz(nVertLevels+1) = 0.0

      !  Note: next we are also dividing through by the cell area after the horizontal flux divergence

!DIR$ IVDEP
         do k=2,nVertLevels
            tend_w(k,iCell) = tend_w(k,iCell) * invAreaCell(iCell) -rdzu(k)*(wdwz(k+1)-wdwz(k))
         end do

         if(rk_step == 1) then
!DIR$ IVDEP
            do k=2,nVertLevels
              tend_w_euler(k,iCell) = tend_w_euler(k,iCell) - cqw(k,iCell)*(   &
                                           rdzu(k)*(pp(k,iCell)-pp(k-1,iCell)) &
                                         - (fzm(k)*dpdz(k,iCell) + fzp(k)*dpdz(k-1,iCell)) )  ! dpdz is the buoyancy term here.
            end do
          end if

      end do

      if (rk_step == 1) then

         if ( v_mom_eddy_visc2 > 0.0 ) then

            do iCell=cellSolveStart,cellSolveEnd
!DIR$ IVDEP
               do k=2,nVertLevels
                  tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + v_mom_eddy_visc2*0.5*(rho_zz(k,iCell)+rho_zz(k-1,iCell))*(  &
                                           (w(k+1,iCell)-w(k  ,iCell))*rdzw(k)                              &
                                          -(w(k  ,iCell)-w(k-1,iCell))*rdzw(k-1) )*rdzu(k)
               end do
            end do

         end if

      end if ! mixing term computed first rk_step

      ! add in mixing terms for w

      do iCell = cellSolveStart,cellSolveEnd
!DIR$ IVDEP
         do k=2,nVertLevels
            tend_w(k,iCell) = tend_w(k,iCell) + tend_w_euler(k,iCell)
         end do
      end do

!----------- rhs for theta

      !
      !  horizontal advection for theta
      !

      do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
         tend_theta(1:nVertLevels,iCell) = 0.0
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)

            flux_arr(1:nVertLevels) = 0.0

            do j=1,nAdvCellsForEdge(iEdge)
               iAdvCell = advCellsForEdge(j,iEdge)
               do k=1,nVertLevels
                  scalar_weight = adv_coefs(j,iEdge) + sign(1.0_RKIND,ru(k,iEdge))*adv_coefs_3rd(j,iEdge)
                  flux_arr(k) = flux_arr(k) + scalar_weight* theta_m(k,iAdvCell)
               end do
            end do

!DIR$ IVDEP
            do k=1,nVertLevels
               tend_theta(k,iCell) = tend_theta(k,iCell) - edgesOnCell_sign(i,iCell) * ru(k,iEdge) * flux_arr(k)
            end do

         end do
      end do

!  addition to pick up perturbation flux for rtheta_pp equation

      if(rk_step > 1) then
        do iCell=cellSolveStart,cellSolveEnd
          do i=1,nEdgesOnCell(iCell) 
            iEdge = edgesOnCell(i,iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
            do k=1,nVertLevels
               flux = edgesOnCell_sign(i,iCell)*dvEdge(iEdge)*(ru_save(k,iEdge)-ru(k,iEdge))*0.5*(theta_m_save(k,cell2)+theta_m_save(k,cell1))
               tend_theta(k,iCell) = tend_theta(k,iCell)-flux  ! division by areaCell picked up down below
            end do
          end do
        end do
      end if

      !
      !  horizontal mixing for theta_m - we could combine this with advection directly (i.e. as a turbulent flux),
      !  but here we can also code in hyperdiffusion if we wish (2nd order at present)
      !

      if (rk_step == 1) then

         delsq_theta(1:nVertLevels,cellStart:cellEnd) = 0.0

         do iCell=cellStart,cellEnd
            tend_theta_euler(1:nVertLevels,iCell) = 0.0
            r_areaCell = invAreaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
               edge_sign = r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * invDcEdge(iEdge)
               pr_scale = prandtl_inv * meshScalingDel2(iEdge)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
               do k=1,nVertLevels

!  we are computing the Smagorinsky filter at more points than needed here so as to pick up the delsq_theta for 4th order filter below

                  theta_turb_flux = edge_sign*(theta_m(k,cell2) - theta_m(k,cell1))*rho_edge(k,iEdge)
                  delsq_theta(k,iCell) = delsq_theta(k,iCell) + theta_turb_flux
                  theta_turb_flux = theta_turb_flux*0.5*(kdiff(k,cell1)+kdiff(k,cell2)) * pr_scale
                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + theta_turb_flux

               end do
            end do
          end do

!$OMP BARRIER
            
         if (h_theta_eddy_visc4 > 0.0) then  ! 4th order mixing is active

            do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
               r_areaCell = h_theta_eddy_visc4 * prandtl_inv * invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)

                  iEdge = edgesOnCell(i,iCell)
                  edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell)*invDcEdge(iEdge)

                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  do k=1,nVertLevels
                     tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) - edge_sign*(delsq_theta(k,cell2) - delsq_theta(k,cell1))
                  end do
               end do
            end do

         end if ! 4th order mixing is active 

      end if ! theta mixing calculated first rk_step

      !
      !  vertical advection plus diabatic term
      !  Note: we are also dividing through by the cell area after the horizontal flux divergence
      !
      do iCell = cellSolveStart,cellSolveEnd

         wdtz(1) = 0.0

         k = 2
         wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))  
         wdtz(k) =  wdtz(k)+(rw_save(k,icell)-rw(k,icell))*(fzm(k)*theta_m_save(k,iCell)+fzp(k)*theta_m_save(k-1,iCell))
         do k=3,nVertLevels-1
            wdtz(k) = flux3( theta_m(k-2,iCell),theta_m(k-1,iCell),theta_m(k,iCell),theta_m(k+1,iCell), rw(k,iCell), coef_3rd_order )
            wdtz(k) =  wdtz(k) + (rw_save(k,icell)-rw(k,iCell))*(fzm(k)*theta_m_save(k,iCell)+fzp(k)*theta_m_save(k-1,iCell))  ! rtheta_pp redefinition
         end do
         k = nVertLevels
         wdtz(k) =  rw_save(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))  ! rtheta_pp redefinition

         wdtz(nVertLevels+1) = 0.0

!DIR$ IVDEP
         do k=1,nVertLevels
            tend_theta(k,iCell) = tend_theta(k,iCell)*invAreaCell(iCell) -rdzw(k)*(wdtz(k+1)-wdtz(k))
            tend_rtheta_adv(k,iCell) = tend_theta(k,iCell)   !  this is for the Tiedke scheme
            rthdynten(k,iCell) = tend_theta(k,iCell)/rho_zz(k,iCell)  !  this is for the Grell-Freitas scheme
            tend_theta(k,iCell) = tend_theta(k,iCell) + rho_zz(k,iCell)*rt_diabatic_tend(k,iCell)
         end do
      end do

      !
      !  vertical mixing for theta - 2nd order 
      !

      if (rk_step == 1) then

         if ( v_theta_eddy_visc2 > 0.0 ) then  ! vertical mixing for theta_m

            if (config_mix_full) then

               do iCell = cellSolveStart,cellSolveEnd
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

               do iCell = cellSolveStart,cellSolveEnd
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

      end if ! compute vertical theta mixing on first rk_step

      do iCell = cellSolveStart,cellSolveEnd
!DIR$ IVDEP
         do k=1,nVertLevels
!            tend_theta(k,iCell) = tend_theta(k,iCell) + tend_theta_euler(k,iCell)
            tend_theta(k,iCell) = tend_theta(k,iCell) + tend_theta_euler(k,iCell) + tend_rtheta_physics(k,iCell)
         end do
      end do

   end subroutine atm_compute_dyn_tend_work


   subroutine atm_compute_solve_diagnostics(dt, state, time_lev, diag, mesh, configs, &
                                            cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                            rk_step )
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
      integer, intent(in), optional :: rk_step          ! which rk_step
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd


      integer, pointer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree
      real (kind=RKIND), dimension(:), pointer :: fVertex, fEdge, invAreaTriangle, invAreaCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, invDvEdge, invDcEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, &
                                                    vorticity, ke, pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, &
                                                    divergence
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, &
                                          kiteForCell, verticesOnCell
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnVertex_sign, edgesOnCell_sign
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge

      real (kind=RKIND), pointer :: config_apvm_upwinding


      call mpas_pool_get_config(configs, 'config_apvm_upwinding', config_apvm_upwinding)

      call mpas_pool_get_array(state, 'rho_zz', h, time_lev)
      call mpas_pool_get_array(state, 'u', u, time_lev)

      call mpas_pool_get_array(diag, 'v', v)
      call mpas_pool_get_array(diag, 'rho_edge', h_edge)
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
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(mesh, 'edgesOnVertex_sign', edgesOnVertex_sign)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'kiteForCell', kiteForCell)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'invDcEdge', invDcEdge)
      call mpas_pool_get_array(mesh, 'invDvEdge', invDvEdge)
      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      call mpas_pool_get_array(mesh, 'invAreaTriangle', invAreaTriangle)
      call mpas_pool_get_array(mesh, 'fVertex', fVertex)
      call mpas_pool_get_array(mesh, 'fEdge', fEdge)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)

      call atm_compute_solve_diagnostics_work(nCells, nEdges, nVertices, &
               vertexDegree, dt, config_apvm_upwinding, &
               fVertex, fEdge, invAreaTriangle, invAreaCell, dvEdge, dcEdge, invDvEdge, invDcEdge, &
               weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, vorticity, ke, pv_edge, pv_vertex, pv_cell, &
               gradPVn, gradPVt, divergence, cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, &
               edgesOnVertex, kiteForCell, verticesOnCell, edgesOnVertex_sign, edgesOnCell_sign, nEdgesOnCell, nEdgesOnEdge, &
               cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
               rk_step)

   end subroutine atm_compute_solve_diagnostics


   subroutine atm_compute_solve_diagnostics_work(nCells, nEdges, nVertices, &
            vertexDegree, dt, config_apvm_upwinding, &
            fVertex, fEdge, invAreaTriangle, invAreaCell, dvEdge, dcEdge, invDvEdge, invDcEdge, &
            weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, vorticity, ke, pv_edge, pv_vertex, pv_cell, &
            gradPVn, gradPVt, divergence, cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, &
            edgesOnVertex, kiteForCell, verticesOnCell, edgesOnVertex_sign, edgesOnCell_sign, nEdgesOnCell, nEdgesOnEdge, &
            cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
            rk_step)

      use mpas_atm_dimensions

      implicit none

      !
      ! Dummy arguments
      !
      integer, intent(in) :: nCells, nEdges, nVertices, vertexDegree
      real (kind=RKIND), intent(in) :: dt, config_apvm_upwinding
      real (kind=RKIND), dimension(nVertices+1) :: fVertex
      real (kind=RKIND), dimension(nEdges+1) :: fEdge
      real (kind=RKIND), dimension(nVertices+1) :: invAreaTriangle
      real (kind=RKIND), dimension(nCells+1) :: invAreaCell
      real (kind=RKIND), dimension(nEdges+1) :: dvEdge
      real (kind=RKIND), dimension(nEdges+1) :: dcEdge
      real (kind=RKIND), dimension(nEdges+1) :: invDvEdge
      real (kind=RKIND), dimension(nEdges+1) :: invDcEdge
      real (kind=RKIND), dimension(maxEdges2,nEdges+1) :: weightsOnEdge
      real (kind=RKIND), dimension(3,nVertices+1) :: kiteAreasOnVertex
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: h_edge
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: h
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: u
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: v
      real (kind=RKIND), dimension(nVertLevels,nVertices+1) :: vorticity
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: ke
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: pv_edge
      real (kind=RKIND), dimension(nVertLevels,nVertices+1) :: pv_vertex
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: pv_cell
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: gradPVn
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: gradPVt
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: divergence
      integer, dimension(2,nEdges+1) :: cellsOnEdge
      integer, dimension(3,nVertices+1) :: cellsOnVertex
      integer, dimension(2,nEdges+1) :: verticesOnEdge
      integer, dimension(maxEdges,nCells+1) :: edgesOnCell
      integer, dimension(maxEdges2,nEdges+1) :: edgesOnEdge
      integer, dimension(3,nVertices+1) :: edgesOnVertex
      integer, dimension(maxEdges,nCells+1) :: kiteForCell
      integer, dimension(maxEdges,nCells+1) :: verticesOnCell
      real (kind=RKIND), dimension(3,nVertices+1) :: edgesOnVertex_sign
      real (kind=RKIND), dimension(maxEdges,nCells+1) :: edgesOnCell_sign
      integer, dimension(nCells+1) :: nEdgesOnCell
      integer, dimension(nEdges+1) :: nEdgesOnEdge

      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd

      integer, intent(in), optional :: rk_step

      !
      ! Local variables
      !
      integer :: iEdge, iCell, iVertex, k, cell1, cell2, eoe, i, j
      real (kind=RKIND) :: h_vertex, r, s
      real (kind=RKIND) :: r1, r2

      logical, parameter :: hollingsworth=.true.
      real (kind=RKIND) :: ke_fact, efac
      logical :: reconstruct_v


      !
      ! Compute height on cell edges at velocity locations
      !
      do iEdge=edgeStart,edgeEnd
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
         do k=1,nVertLevels
            h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
         end do

!  the first openmp barrier below is set so that ke_edge is computed
!  it would be good to move this somewhere else?

         efac = dcEdge(iEdge)*dvEdge(iEdge)
         do k=1,nVertLevels
            ke_edge(k,iEdge) = efac*u(k,iEdge)**2
         end do

      end do

      !
      ! Compute circulation and relative vorticity at each vertex
      !
      do iVertex=vertexStart,vertexEnd
         vorticity(1:nVertLevels,iVertex) = 0.0
         do i=1,vertexDegree
            iEdge = edgesOnVertex(i,iVertex)
            s = edgesOnVertex_sign(i,iVertex) * dcEdge(iEdge)
!DIR$ IVDEP
            do k=1,nVertLevels
               vorticity(k,iVertex) = vorticity(k,iVertex) + s * u(k,iEdge)
            end do
         end do
!DIR$ IVDEP
         do k=1,nVertLevels
            vorticity(k,iVertex) = vorticity(k,iVertex) * invAreaTriangle(iVertex)
         end do
      end do


      !
      ! Compute the divergence at each cell center
      !
      do iCell=cellStart,cellEnd
         divergence(1:nVertLevels,iCell) = 0.0
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            s = edgesOnCell_sign(i,iCell) * dvEdge(iEdge)
!DIR$ IVDEP
            do k=1,nVertLevels
              divergence(k,iCell) = divergence(k,iCell) + s * u(k,iEdge)
            end do
         end do
         r = invAreaCell(iCell)
         do k = 1,nVertLevels
            divergence(k,iCell) = divergence(k,iCell) * r
         end do
      end do


!$OMP BARRIER

      !
      ! Compute kinetic energy in each cell (Ringler et al JCP 2009)
      !
      ! Replace 2.0 with 2 in exponentiation to avoid outside chance that
      ! compiler will actually allow "float raised to float" operation
      do iCell=cellStart,cellEnd
         ke(1:nVertLevels,iCell) = 0.0
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            do k=1,nVertLevels
!               ke(k,iCell) = ke(k,iCell) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2
               ke(k,iCell) = ke(k,iCell) + 0.25 * ke_edge(k,iEdge)
            end do
         end do
!DIR$ IVDEP
         do k=1,nVertLevels
            ke(k,iCell) = ke(k,iCell) * invAreaCell(iCell)
         end do
      end do


      if (hollingsworth) then

         ! Compute ke at cell vertices - AG's new KE construction, part 1
         ! *** approximation here because we don't have inner triangle areas
         !

         ! Replace 2.0 with 2 in exponentiation to avoid outside chance that
         ! compiler will actually allow "float raised to float" operation
         do iVertex=vertexStart,vertexEnd
            r = 0.25 * invAreaTriangle(iVertex) 
            do k=1,nVertLevels

!               ke_vertex(k,iVertex) = (  dcEdge(EdgesOnVertex(1,iVertex))*dvEdge(EdgesOnVertex(1,iVertex))*u(k,EdgesOnVertex(1,iVertex))**2  &
!                                        +dcEdge(EdgesOnVertex(2,iVertex))*dvEdge(EdgesOnVertex(2,iVertex))*u(k,EdgesOnVertex(2,iVertex))**2  &
!                                        +dcEdge(EdgesOnVertex(3,iVertex))*dvEdge(EdgesOnVertex(3,iVertex))*u(k,EdgesOnVertex(3,iVertex))**2  &
!                                      ) * r

               ke_vertex(k,iVertex) = (  ke_edge(k,EdgesOnVertex(1,iVertex))+ke_edge(k,EdgesOnVertex(2,iVertex))+ke_edge(k,EdgesOnVertex(3,iVertex)) )*r

            end do
         end do

!$OMP BARRIER

         ! adjust ke at cell vertices - AG's new KE construction, part 2
         !

         ke_fact = 1.0 - .375

         do iCell=cellStart,cellEnd
            do k=1,nVertLevels
               ke(k,iCell) = ke_fact * ke(k,iCell)
            end do
         end do


         do iCell=cellStart,cellEnd
            r = invAreaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iVertex = verticesOnCell(i,iCell)
               j = kiteForCell(i,iCell)
!DIR$ IVDEP
               do k = 1,nVertLevels
                  ke(k,iCell) = ke(k,iCell) + (1.-ke_fact)*kiteAreasOnVertex(j,iVertex) * ke_vertex(k,iVertex) * r
               end do
            end do
         end do

      end if

      !
      ! Compute v (tangential) velocities following Thuburn et al JCP 2009
      ! The tangential velocity is only used to compute the Smagorinsky coefficient

      reconstruct_v = .true.
      if(present(rk_step)) then
        if(rk_step /= 3) reconstruct_v = .false.
      end if

      if (reconstruct_v) then
        do iEdge = edgeStart,edgeEnd
          v(1:nVertLevels,iEdge) = 0.0
          do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
!DIR$ IVDEP
            do k = 1,nVertLevels
              v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
            end do
          end do
        end do
      end if

      !
      ! Compute height at vertices, pv at vertices, and average pv to edge locations
      !  ( this computes pv_vertex at all vertices bounding real cells )
      !
      ! Avoid dividing h_vertex by areaTriangle and move areaTriangle into
      ! numerator for the pv_vertex calculation
      do iVertex = vertexStart,vertexEnd
!DIR$ IVDEP
         do k=1,nVertLevels
!
! the following commented code is for the PV conserving shallow water solver.  
!            h_vertex = 0.0
!            do i=1,vertexDegree
!               h_vertex = h_vertex + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
!            end do
!            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) * areaTriangle(iVertex) / h_vertex
            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex))
         end do
      end do

!$OMP BARRIER

      !
      ! Compute pv at the edges
      !   ( this computes pv_edge at all edges bounding real cells )
      !
      do iEdge = edgeStart,edgeEnd
!DIR$ IVDEP
         do k=1,nVertLevels
            pv_edge(k,iEdge) =  0.5 * (pv_vertex(k,verticesOnEdge(1,iEdge)) + pv_vertex(k,verticesOnEdge(2,iEdge)))
         end do
      end do

      if (config_apvm_upwinding > 0.0) then

      !
      ! Compute pv at cell centers
      !    ( this computes pv_cell for all real cells )
      !  only needed for APVM upwinding
      !
      do iCell=cellStart,cellEnd
         pv_cell(1:nVertLevels,iCell) = 0.0
         r = invAreaCell(iCell)
         do i=1,nEdgesOnCell(iCell)
            iVertex = verticesOnCell(i,iCell)
            j = kiteForCell(i,iCell)
!DIR$ IVDEP
            do k = 1,nVertLevels
               pv_cell(k,iCell) = pv_cell(k,iCell) + kiteAreasOnVertex(j,iVertex) * pv_vertex(k,iVertex) * r
            end do
         end do
      end do


!$OMP BARRIER

         !
         ! Modify PV edge with upstream bias. 
         !
         ! Compute gradient of PV in the tangent direction
         !   ( this computes gradPVt at all edges bounding real cells )
         !
         ! Compute gradient of PV in normal direction
         !   (tdr: 2009-10-02: this is not correct because the pv_cell in the halo is not correct)
         !
         ! Modify PV edge with upstream bias.
         !
         ! Merged loops for calculating gradPVt, gradPVn and pv_edge
         ! Also precomputed inverses of dvEdge and dcEdge to avoid repeated divisions
         !
         r = config_apvm_upwinding * dt
         do iEdge = edgeStart,edgeEnd
            r1 = 1.0_RKIND * invDvEdge(iEdge)
            r2 = 1.0_RKIND * invDcEdge(iEdge)
!DIR$ IVDEP
            do k = 1,nVertLevels
               gradPVt(k,iEdge) = (pv_vertex(k,verticesOnEdge(2,iEdge)) - pv_vertex(k,verticesOnEdge(1,iEdge))) * r1
               gradPVn(k,iEdge) = (pv_cell(k,cellsOnEdge(2,iEdge)) - pv_cell(k,cellsOnEdge(1,iEdge))) * r2
               pv_edge(k,iEdge) = pv_edge(k,iEdge) - r * (v(k,iEdge) * gradPVt(k,iEdge) + u(k,iEdge) * gradPVn(k,iEdge))
            end do
         end do

      end if  ! apvm upwinding

   end subroutine atm_compute_solve_diagnostics_work


   subroutine atm_init_coupled_diagnostics(state, time_lev, diag, mesh, configs, &
                                       cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                       cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)

      implicit none
   
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev                    ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd

      integer :: i, k, iCell, iEdge, cell1, cell2
      integer, pointer :: nCells, nEdges, nVertLevels
      integer, pointer :: index_qv
      real (kind=RKIND) :: p0, rcv, flux
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign
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
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3, zb_cell, zb3_cell


      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)

      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)

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
      call mpas_pool_get_array(mesh, 'zb_cell', zb_cell)
      call mpas_pool_get_array(mesh, 'zb3_cell', zb3_cell)


      rcv = rgas / (cp-rgas)
      p0 = 1.e5  ! this should come from somewhere else...

      do iCell=cellStart,cellEnd
         do k=1,nVertLevels
            theta_m(k,iCell) = theta(k,iCell) * (1._RKIND + rvord * scalars(index_qv,k,iCell))
            rho_zz(k,iCell) = rho(k,iCell) / zz(k,iCell)
         end do
      end do

!$OMP BARRIER

      do iEdge=edgeStart,edgeEnd
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
            ru(k,iEdge) = 0.5 * u(k,iEdge) * (rho_zz(k,cell1) + rho_zz(k,cell2))
         end do
      end do

!$OMP BARRIER

      ! Compute rw (i.e. rho_zz * omega) from rho_zz, w, and ru.
      ! We are reversing the procedure we use in subroutine atm_recover_large_step_variables.
      ! first, the piece that depends on w.
      do iCell=cellStart,cellEnd
         rw(1,iCell) = 0.0
         rw(nVertLevels+1,iCell) = 0.0
         do k=2,nVertLevels
            rw(k,iCell) = w(k,iCell)     &
                          * (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell)) &
                          * (fzp(k) * zz(k-1,iCell) + fzm(k) * zz(k,iCell))
         end do
      end do
  
      ! next, the piece that depends on ru
      do iCell=cellStart,cellEnd
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            do k = 2,nVertLevels
            flux = (fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge))
            rw(k,iCell) = rw(k,iCell)   &
                          - edgesOnCell_sign(i,iCell) * (zb_cell(k,i,iCell) + sign(1.0_RKIND,flux) * zb3_cell(k,i,iCell))*flux   &
                          * (fzp(k) * zz(k-1,iCell) + fzm(k) * zz(k,iCell))
            end do
         end do
      end do

      do iCell=cellStart,cellEnd
         do k=1,nVertLevels
            rho_p(k,iCell) = rho_zz(k,iCell) - rho_base(k,iCell)
         end do
      end do

      do iCell=cellStart,cellEnd
         do k=1,nVertLevels
            rtheta_base(k,iCell) = theta_base(k,iCell) * rho_base(k,iCell)
         end do
      end do

      do iCell=cellStart,cellEnd
         do k=1,nVertLevels
            rtheta_p(k,iCell) = theta_m(k,iCell) * rho_p(k,iCell)  &
                                             + rho_base(k,iCell)   * (theta_m(k,iCell) - theta_base(k,iCell))
         end do
      end do

      do iCell=cellStart,cellEnd
         do k=1,nVertLevels
            exner(k,iCell) = (zz(k,iCell) * (rgas/p0) * (rtheta_p(k,iCell) + rtheta_base(k,iCell)))**rcv
         end do
      end do

      do iCell=cellStart,cellEnd
         do k=1,nVertLevels
            pressure_p(k,iCell) = zz(k,iCell) * rgas &
                                               * (  exner(k,iCell) * rtheta_p(k,iCell) &
                                                  + rtheta_base(k,iCell) * (exner(k,iCell) - exner_base(k,iCell)) &
                                                 )
         end do
      end do

   end subroutine atm_init_coupled_diagnostics


   subroutine atm_rk_dynamics_substep_finish( state, diag, dynamics_substep, dynamics_split, &
                                       cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                       cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)

      implicit none

      !  this routine resets the dry dynamics variables at the end of an rk3 substep for the case
      !  where the dry dynamics is split from the scalar transport (i.e. where the dry dynamics is
      !  using a different, usually smaller, timestep.
      !
      !  WCS 18 November 2014

      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      integer, intent(in) :: dynamics_substep, dynamics_split
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd

      real (kind=RKIND) :: inv_dynamics_split
      
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
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz_1, rho_zz_2, rho_zz_old_split
      real (kind=RKIND), dimension(:,:), pointer :: ruAvg, wwAvg, ruAvg_split, wwAvg_split

      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'ru_save', ru_save)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'rw_save', rw_save)
      call mpas_pool_get_array(diag, 'rtheta_p', rtheta_p)
      call mpas_pool_get_array(diag, 'rtheta_p_save', rtheta_p_save)
      call mpas_pool_get_array(diag, 'rho_p', rho_p)
      call mpas_pool_get_array(diag, 'rho_p_save', rho_p_save)
      call mpas_pool_get_array(diag, 'rho_zz_old_split', rho_zz_old_split)
      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(diag, 'ruAvg_split', ruAvg_split)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
      call mpas_pool_get_array(diag, 'wwAvg_split', wwAvg_split)

      call mpas_pool_get_array(state, 'u', u_1, 1)
      call mpas_pool_get_array(state, 'u', u_2, 2)
      call mpas_pool_get_array(state, 'w', w_1, 1)
      call mpas_pool_get_array(state, 'w', w_2, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m_1, 1)
      call mpas_pool_get_array(state, 'theta_m', theta_m_2, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_1, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_2, 2)

      inv_dynamics_split = 1.0_RKIND / real(dynamics_split)
      
      if (dynamics_substep < dynamics_split) then

         ru_save(:,edgeStart:edgeEnd) = ru(:,edgeStart:edgeEnd)
         rw_save(:,cellStart:cellEnd) = rw(:,cellStart:cellEnd)
         rtheta_p_save(:,cellStart:cellEnd) = rtheta_p(:,cellStart:cellEnd)
         rho_p_save(:,cellStart:cellEnd) = rho_p(:,cellStart:cellEnd)

         u_1(:,edgeStart:edgeEnd) = u_2(:,edgeStart:edgeEnd)
         w_1(:,cellStart:cellEnd) = w_2(:,cellStart:cellEnd)
         theta_m_1(:,cellStart:cellEnd) = theta_m_2(:,cellStart:cellEnd)
         rho_zz_1(:,cellStart:cellEnd) = rho_zz_2(:,cellStart:cellEnd)

      end if

      if (dynamics_substep == 1) then
         ruAvg_split(:,edgeStart:edgeEnd) = ruAvg(:,edgeStart:edgeEnd)
         wwAvg_split(:,cellStart:cellEnd) = wwAvg(:,cellStart:cellEnd)
      else
         ruAvg_split(:,edgeStart:edgeEnd) = ruAvg(:,edgeStart:edgeEnd)+ruAvg_split(:,edgeStart:edgeEnd)
         wwAvg_split(:,cellStart:cellEnd) = wwAvg(:,cellStart:cellEnd)+wwAvg_split(:,cellStart:cellEnd)
      end if

      if (dynamics_substep == dynamics_split) then
         ruAvg(:,edgeStart:edgeEnd) = ruAvg_split(:,edgeStart:edgeEnd) * inv_dynamics_split
         wwAvg(:,cellStart:cellEnd) = wwAvg_split(:,cellStart:cellEnd) * inv_dynamics_split
         rho_zz_1(:,cellStart:cellEnd) = rho_zz_old_split(:,cellStart:cellEnd)
      end if

   end subroutine atm_rk_dynamics_substep_finish


   subroutine summarize_timestep(domain)

       use ieee_arithmetic, only : ieee_is_nan

       implicit none

       type (domain_type), intent(inout) :: domain

       real (kind=RKIND), parameter :: pi_const = 2.0_RKIND*asin(1.0_RKIND)

       logical, pointer :: config_print_global_minmax_vel
       logical, pointer :: config_print_detailed_minmax_vel
       logical, pointer :: config_print_global_minmax_sca

       type (block_type), pointer :: block

       integer :: iCell, k, iEdge, iScalar
       integer, pointer :: num_scalars, nCellsSolve, nEdgesSolve, nVertLevels

       type (mpas_pool_type), pointer :: state
       type (mpas_pool_type), pointer :: diag
       type (mpas_pool_type), pointer :: mesh

       real (kind=RKIND) :: scalar_min, scalar_max
       real (kind=RKIND) :: global_scalar_min, global_scalar_max

       real (kind=RKIND), dimension(:), pointer :: latCell
       real (kind=RKIND), dimension(:), pointer :: lonCell
       real (kind=RKIND), dimension(:), pointer :: latEdge
       real (kind=RKIND), dimension(:), pointer :: lonEdge
       integer, dimension(:), pointer :: indexToCellID
       integer :: indexMax, indexMax_global
       integer :: kMax, kMax_global
       real (kind=RKIND) :: latMax, latMax_global
       real (kind=RKIND) :: lonMax, lonMax_global
       real (kind=RKIND), dimension(5) :: localVals, globalVals

       real (kind=RKIND) :: spd
       real (kind=RKIND), dimension(:,:), pointer :: w
       real (kind=RKIND), dimension(:,:), pointer :: u, v, uReconstructZonal, uReconstructMeridional, uReconstructX, uReconstructY, uReconstructZ
       real (kind=RKIND), dimension(:,:,:), pointer :: scalars, scalars_1, scalars_2

       call mpas_pool_get_config(domain % blocklist % configs, 'config_print_global_minmax_vel', config_print_global_minmax_vel)
       call mpas_pool_get_config(domain % blocklist % configs, 'config_print_detailed_minmax_vel', config_print_detailed_minmax_vel)
       call mpas_pool_get_config(domain % blocklist % configs, 'config_print_global_minmax_sca', config_print_global_minmax_sca)

      if (config_print_detailed_minmax_vel) then
         write(0,*) ' '

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_get_subpool(block % structs, 'diag', diag)
            call mpas_pool_get_subpool(block % structs, 'mesh', mesh)

            call mpas_pool_get_array(state, 'w', w, 2)
            call mpas_pool_get_array(state, 'u', u, 2)
            call mpas_pool_get_array(diag, 'v', v)
            call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)
            call mpas_pool_get_array(mesh, 'latCell', latCell)
            call mpas_pool_get_array(mesh, 'lonCell', lonCell)
            call mpas_pool_get_array(mesh, 'latEdge', latEdge)
            call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
            call mpas_pool_get_dimension(state, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_dimension(state, 'nEdgesSolve', nEdgesSolve)
            call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels)

            scalar_min = 1.0e20
            indexMax = -1
            kMax = -1
            latMax = 0.0
            lonMax = 0.0
            do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               if (w(k,iCell) < scalar_min) then
                  scalar_min = w(k,iCell)
                  indexMax = iCell
                  kMax = k
                  latMax = latCell(iCell)
                  lonMax = lonCell(iCell)
               end if
            end do
            end do
            localVals(1) = scalar_min
            localVals(2) = real(indexMax,kind=RKIND)
            localVals(3) = real(kMax,kind=RKIND)
            localVals(4) = latMax
            localVals(5) = lonMax
            call mpas_dmpar_minattributes_real(domain % dminfo, scalar_min, localVals, globalVals)
            global_scalar_min = globalVals(1)
            indexMax_global = int(globalVals(2))
            kMax_global = int(globalVals(3))
            latMax_global = globalVals(4)
            lonMax_global = globalVals(5)
            latMax_global = latMax_global * 180.0_RKIND / pi_const
            lonMax_global = lonMax_global * 180.0_RKIND / pi_const
            if (lonMax_global > 180.0) then
               lonMax_global = lonMax_global - 360.0
            end if
            write(0,'(a,f9.4,a,i4,a,f7.3,a,f8.3,a)') ' global min w: ', global_scalar_min, &
                                                   ' k=', kMax_global, ', ', latMax_global, ' lat ', lonMax_global, ' lon'

            scalar_max = -1.0e20
            indexMax = -1
            kMax = -1
            latMax = 0.0
            lonMax = 0.0
            do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               if (w(k,iCell) > scalar_max) then
                  scalar_max = w(k,iCell)
                  indexMax = iCell
                  kMax = k
                  latMax = latCell(iCell)
                  lonMax = lonCell(iCell)
               end if
            end do
            end do
            localVals(1) = scalar_max
            localVals(2) = real(indexMax,kind=RKIND)
            localVals(3) = real(kMax,kind=RKIND)
            localVals(4) = latMax
            localVals(5) = lonMax
            call mpas_dmpar_maxattributes_real(domain % dminfo, scalar_max, localVals, globalVals)
            global_scalar_max = globalVals(1)
            indexMax_global = int(globalVals(2))
            kMax_global = int(globalVals(3))
            latMax_global = globalVals(4)
            lonMax_global = globalVals(5)
            latMax_global = latMax_global * 180.0_RKIND / pi_const
            lonMax_global = lonMax_global * 180.0_RKIND / pi_const
            if (lonMax_global > 180.0) then
               lonMax_global = lonMax_global - 360.0
            end if
            write(0,'(a,f9.4,a,i4,a,f7.3,a,f8.3,a)') ' global max w: ', global_scalar_max, &
                                                   ' k=', kMax_global, ', ', latMax_global, ' lat ', lonMax_global, ' lon'

            scalar_min = 1.0e20
            indexMax = -1
            kMax = -1
            latMax = 0.0
            lonMax = 0.0
            do iEdge = 1, nEdgesSolve
            do k = 1, nVertLevels
               if (u(k,iEdge) < scalar_min) then
                  scalar_min = u(k,iEdge)
                  indexMax = iEdge
                  kMax = k
                  latMax = latEdge(iEdge)
                  lonMax = lonEdge(iEdge)
               end if
            end do
            end do
            localVals(1) = scalar_min
            localVals(2) = real(indexMax,kind=RKIND)
            localVals(3) = real(kMax,kind=RKIND)
            localVals(4) = latMax
            localVals(5) = lonMax
            call mpas_dmpar_minattributes_real(domain % dminfo, scalar_min, localVals, globalVals)
            global_scalar_min = globalVals(1)
            indexMax_global = int(globalVals(2))
            kMax_global = int(globalVals(3))
            latMax_global = globalVals(4)
            lonMax_global = globalVals(5)
            latMax_global = latMax_global * 180.0_RKIND / pi_const
            lonMax_global = lonMax_global * 180.0_RKIND / pi_const
            if (lonMax_global > 180.0) then
               lonMax_global = lonMax_global - 360.0
            end if
            write(0,'(a,f9.4,a,i4,a,f7.3,a,f8.3,a)') ' global min u: ', global_scalar_min, &
                                                   ' k=', kMax_global, ', ', latMax_global, ' lat ', lonMax_global, ' lon'

            scalar_max = -1.0e20
            indexMax = -1
            kMax = -1
            latMax = 0.0
            lonMax = 0.0
            do iEdge = 1, nEdgesSolve
            do k = 1, nVertLevels
               if (u(k,iEdge) > scalar_max) then
                  scalar_max = u(k,iEdge)
                  indexMax = iEdge
                  kMax = k
                  latMax = latEdge(iEdge)
                  lonMax = lonEdge(iEdge)
               end if
            end do
            end do
            localVals(1) = scalar_max
            localVals(2) = real(indexMax,kind=RKIND)
            localVals(3) = real(kMax,kind=RKIND)
            localVals(4) = latMax
            localVals(5) = lonMax
            call mpas_dmpar_maxattributes_real(domain % dminfo, scalar_max, localVals, globalVals)
            global_scalar_max = globalVals(1)
            indexMax_global = int(globalVals(2))
            kMax_global = int(globalVals(3))
            latMax_global = globalVals(4)
            lonMax_global = globalVals(5)
            latMax_global = latMax_global * 180.0_RKIND / pi_const
            lonMax_global = lonMax_global * 180.0_RKIND / pi_const
            if (lonMax_global > 180.0) then
               lonMax_global = lonMax_global - 360.0
            end if
            write(0,'(a,f9.4,a,i4,a,f7.3,a,f8.3,a)') ' global max u: ', global_scalar_max, &
                                                   ' k=', kMax_global, ', ', latMax_global, ' lat ', lonMax_global, ' lon'

            scalar_max = -1.0e20
            indexMax = -1
            kMax = -1
            latMax = 0.0
            lonMax = 0.0
            do iEdge = 1, nEdgesSolve
            do k = 1, nVertLevels
               spd = sqrt(u(k,iEdge)*u(k,iEdge) + v(k,iEdge)*v(k,iEdge))
               if (spd > scalar_max) then
                  scalar_max = spd
                  indexMax = iEdge
                  kMax = k
                  latMax = latEdge(iEdge)
                  lonMax = lonEdge(iEdge)
               end if
            end do
            end do
            localVals(1) = scalar_max
            localVals(2) = real(indexMax,kind=RKIND)
            localVals(3) = real(kMax,kind=RKIND)
            localVals(4) = latMax
            localVals(5) = lonMax
            call mpas_dmpar_maxattributes_real(domain % dminfo, scalar_max, localVals, globalVals)
            global_scalar_max = globalVals(1)
            indexMax_global = int(globalVals(2))
            kMax_global = int(globalVals(3))
            latMax_global = globalVals(4)
            lonMax_global = globalVals(5)
            latMax_global = latMax_global * 180.0_RKIND / pi_const
            lonMax_global = lonMax_global * 180.0_RKIND / pi_const
            if (lonMax_global > 180.0) then
               lonMax_global = lonMax_global - 360.0
            end if
            write(0,'(a,f9.4,a,i4,a,f7.3,a,f8.3,a)') ' global max wsp: ', global_scalar_max, &
                                                   ' k=', kMax_global, ', ', latMax_global, ' lat ', lonMax_global, ' lon'

            !
            ! Check for NaNs
            !
            do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               if (ieee_is_nan(w(k,iCell))) then
                  call mpas_dmpar_global_abort('NaN detected in ''w'' field.')
               end if
            end do
            end do

            do iEdge = 1, nEdgesSolve
            do k = 1, nVertLevels
               if (ieee_is_nan(u(k,iEdge))) then
                  call mpas_dmpar_global_abort('NaN detected in ''u'' field.')
               end if
            end do
            end do

            block => block % next
         end do

      else if (config_print_global_minmax_vel) then
         write(0,*) ' '

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', state)

            call mpas_pool_get_array(state, 'w', w, 2)
            call mpas_pool_get_array(state, 'u', u, 2)
            call mpas_pool_get_dimension(state, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_dimension(state, 'nEdgesSolve', nEdgesSolve)
            call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels)

            scalar_min = 0.0
            scalar_max = 0.0
            do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               scalar_min = min(scalar_min, w(k,iCell))
               scalar_max = max(scalar_max, w(k,iCell))
            end do
            end do
            call mpas_dmpar_min_real(domain % dminfo, scalar_min, global_scalar_min)
            call mpas_dmpar_max_real(domain % dminfo, scalar_max, global_scalar_max)
            write(0,*) 'global min, max w ', global_scalar_min, global_scalar_max

            scalar_min = 0.0
            scalar_max = 0.0
            do iEdge = 1, nEdgesSolve
            do k = 1, nVertLevels
               scalar_min = min(scalar_min, u(k,iEdge))
               scalar_max = max(scalar_max, u(k,iEdge))
            end do
            end do
            call mpas_dmpar_min_real(domain % dminfo, scalar_min, global_scalar_min)
            call mpas_dmpar_max_real(domain % dminfo, scalar_max, global_scalar_max)
            write(0,*) 'global min, max u ', global_scalar_min, global_scalar_max

            block => block % next
         end do
      end if

      if (config_print_global_minmax_sca) then
         if (.not. (config_print_global_minmax_vel .or. config_print_detailed_minmax_vel)) write(0,*) ' '

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', state)

            call mpas_pool_get_array(state, 'scalars', scalars, 2)
            call mpas_pool_get_dimension(state, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels)
            call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)

            do iScalar = 1, num_scalars
               scalar_min = 0.0
               scalar_max = 0.0
               do iCell = 1, nCellsSolve
               do k = 1, nVertLevels
                  scalar_min = min(scalar_min, scalars(iScalar,k,iCell))
                  scalar_max = max(scalar_max, scalars(iScalar,k,iCell))
               end do
               end do
               call mpas_dmpar_min_real(domain % dminfo, scalar_min, global_scalar_min)
               call mpas_dmpar_max_real(domain % dminfo, scalar_max, global_scalar_max)
               write(0,'(a,i4,2(1x,e17.10))') ' global min, max scalar ', iScalar, global_scalar_min, global_scalar_max
            end do

            block => block % next
         end do
      end if

   end subroutine summarize_timestep

end module atm_time_integration
