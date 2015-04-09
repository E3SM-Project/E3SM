! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_time_integration_fe
!
!> \brief MPAS land ice Forward Euler time integration scheme
!> \author Matt Hoffman
!> \date   17 April 2011
!> \details
!>  This module contains the Forward Euler time integration scheme
!
!-----------------------------------------------------------------------

module li_time_integration_fe

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use mpas_timer
   use mpas_vector_reconstruction
   use mpas_configure
   use li_velocity, only: li_velocity_solve
   use li_tendency
   use li_diagnostic_vars
   use li_setup

   implicit none
   private

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: li_time_integrator_forwardeuler

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************
   contains
!***********************************************************************


!***********************************************************************
!
!  routine li_time_integrator_forwardeuler
!
!> \brief   Forward Euler time integration scheme 
!> \author  Matthew Hoffman
!> \date    10 January 2012
!> \details 
!>  This routine performs Forward Euler time integration.
!
!-----------------------------------------------------------------------
   subroutine li_time_integrator_forwardeuler(domain, deltat, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: deltat  !< Input: time step

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object 

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      integer :: err_tmp


      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      ! (time level 1 should not be modified.)

!!!      procVertexMaskChanged = 0

! === Implicit column physics (vertical temperature diffusion) ===========
      !call ()

! === Calculate Tendencies ========================
      call mpas_timer_start("calculate tendencies")
      call calculate_tendencies(domain, deltat, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("calculate tendencies")


! === Compute new state for prognostic variables ==================================
! (once implicit column physics are added (i.e. temp diffusion), these calculations will need to be adjusted to apply to the new values as needed)
      call mpas_timer_start("calc. new prognostic vars")
      call update_prognostics(domain, deltat, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("calc. new prognostic vars")


! === Calculate diagnostic variables for new state =====================

      ! Some velocity solvers require an initial guess.
      ! If needed, that guess should be inserted into normalVelocity
      ! in time level 2 before calling li_calculate_diagnostic_vars.

      call li_calculate_diagnostic_vars(domain, timeLevel=2, solveVelo=.true., err=err_tmp)
      err = ior(err, err_tmp)


! === Cleanup & Misc. =============================

      ! === error check
      if (err == 1) then
          write (0,*) "An error has occurred in li_time_integrator_forwardeuler."
      endif

   !--------------------------------------------------------------------
   end subroutine li_time_integrator_forwardeuler



!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************


!***********************************************************************
!
!  routine calculate_tendencies
!
!> \brief   Calculates tendencies 
!> \author  Matthew Hoffman
!> \date    12 September 2013
!> \details
!>  This routine calculates tendencies of thickness (and eventually tracers)
!>  to use with the Forward Euler time integrator.
!
!-----------------------------------------------------------------------

   subroutine calculate_tendencies(domain, deltat, err)

      use mpas_timekeeping

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND) :: deltat

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (dm_info), pointer :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: tendPool

      real (kind=RKIND), dimension(:,:), pointer :: layerThickness_tend
      type (field2DReal), pointer :: layerThickness_tend_field
      integer :: allowableDtProcNumber, allowableDtMinProcNumber
      real (kind=RKIND) :: allowableDt, allowableDtMin
      logical, pointer :: config_print_thickness_advection_info
      type (MPAS_TimeInterval_type) :: allowableDtMinStringInterval
      character (len=StrKIND) :: allowableDtMinString
      integer :: err_tmp

      integer :: y, m, d, hh, mm, ss

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)

      dminfo => domain % dminfo

      ! ===
      ! === Thickness tendencies
      ! ===
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_array(tendPool, 'layerThickness', layerThickness_tend)

         ! Calculate thickness tendency using state at time n =========
         call li_tendency_thickness(meshPool, statePool, layerThickness_tend, deltat, dminfo, allowableDt, err_tmp)
         err = ior(err,err_tmp)

         block => block % next
      end do

      ! Now that we have exited the block loop, do any needed halo updates.
      ! update halos on thickness tend
      call mpas_timer_start("halo updates")
      call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)
      call mpas_pool_get_field(tendPool, 'layerThickness', layerThickness_tend_field)
      call mpas_dmpar_exch_halo_field(layerThickness_tend_field)
      call mpas_timer_stop("halo updates")

      ! If we are printing advection debug information, 
      ! then find out what the CFL limit is.  Don't do this otherwise because it
      ! is requires 2 unnecessary MPI communications.
      if (config_print_thickness_advection_info) then
          ! Determine CFL limit on all procs
          call mpas_dmpar_min_real(dminfo, allowableDt, allowableDtMin)
          ! Determine which processor has the limiting CFL
          if (allowableDt .eq. allowableDtMin) then
            allowableDtProcNumber = dminfo % my_proc_id
          else 
            allowableDtProcNumber = -1
          endif
          call mpas_dmpar_max_int(dminfo, allowableDtProcNumber, allowableDtMinProcNumber)
          call mpas_set_timeInterval(allowableDtMinStringInterval, dt=allowableDtMin, ierr=err_tmp)
          err = ior(err,err_tmp)
          call mpas_get_timeInterval(allowableDtMinStringInterval, timeString=allowableDtMinString, ierr=err_tmp)
          err = ior(err,err_tmp)
          write(6,*) '  Maximum allowable time step (yr) for all processors is (Days_hhh:mmm:sss): ' // trim(allowableDtMinString) //  '  Time step is limited by processor number ', allowableDtMinProcNumber
      endif

      if (err .gt. 0) then
           write(0,*) 'Error in calculating thickness tendency (possibly CFL violation)'
      endif


      ! ===
      ! === Tracer tendencies
      ! ===
! Commenting out tracer tendencies
!!!      ! Tracer tendencies
!!!      block => domain % blocklist
!!!      do while (associated(block))
!!!         ! Mesh information
!!!         mesh => block % mesh
!!!         ! State at time n
!!!         stateOld => block % state % time_levs(1) % state
!!!         ! Tendencies
!!!         layerThickness_tend => block % tend % layerThickness % array
!!!         tracer_tendency => block % tend % tracers % array

!!!         ! Calculate tracer tendencies ==========
!!!         ! There could be a negative layer thickness with SMB turned on!
!!!         call li_tendency_tracers(mesh, stateOld, layerThickness_tend, tracer_tendency, dt, dminfo, err)

!!!         block => block % next
!!!      end do 

!!!      if (err == 1) then
!!!          call mpas_dmpar_global_abort("An error has occurred in li_tendency_tracers. Aborting...")
!!!      endif

!!!      ! update halos on tracer tend 
!!!      select case (config_tracer_advection)
!!!      case ('None')  !===================================================
!!!          ! Do nothing - no need to waste time doing a halo update if not advecting tracers!  The tendency will be 0 everywhere
!!!      case default
!!!          call mpas_dmpar_exch_halo_field3d_real(dminfo, tracer_tendency, &
!!!                                         size(tracer_tendency,dim=1), mesh % nVertLevels, mesh % nCells, &
!!!                                         domain % blocklist % parinfo % cellsToSend, domain % blocklist % parinfo % cellsToRecv)
!!!      end select


      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in calculate_tendencies."
      endif

   !--------------------------------------------------------------------
   end subroutine calculate_tendencies



!***********************************************************************
!
!  routine update_prognostics
!
!> \brief   Updates prognostic variables 
!> \author  Matthew Hoffman
!> \date    12 September 2013
!> \details
!>  This routine updates prognostic variables using a 
!>  Forward Euler time integrator.
!
!-----------------------------------------------------------------------

   subroutine update_prognostics(domain, deltat, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND) :: deltat

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (dm_info), pointer :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool, statePool, tendPool
      integer, pointer :: nCells
      logical, pointer :: config_print_thickness_advection_info

      real (kind=RKIND), dimension(:), pointer :: thicknessNew, layerThicknessFractions
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessOld, layerThicknessNew, layerThickness_tend
!!!      real (kind=RKIND), dimension(:,:,:), pointer :: tracer_tendency, tracersNew, tracersOld
!!!      integer, dimension(:), pointer :: cellMaskOld
      integer, dimension(:), allocatable :: masktmp  ! Temporary mask for assessing new thickness field

      err = 0

      dminfo => domain % dminfo
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)

      block => domain % blocklist
      do while (associated(block))
         ! Mesh information

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)

         ! State at time n
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessOld, timeLevel=1)
!!!         tracersOld => stateOld % tracers % array
!!!         cellMaskOld => stateOld % cellMask % array

         ! State at time n+1 (advanced by dt by Forward Euler)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, timeLevel=2)
         call mpas_pool_get_array(statePool, 'thickness', thicknessNew, timeLevel=2)
!!!         tracersNew => stateNew % tracers % array

         ! Tendencies
         call mpas_pool_get_array(tendPool, 'layerThickness', layerThickness_tend)
!!!         tracer_tendency => block % tend % tracers % array


         ! Update thickness ======================
         
         ! Commented out usage for advecting thickness as a column
         !!!thicknessNew = thicknessOld + thickness_tend * dt / SecondsInYear   
         ! Commented out usage for using FCT for thickness 
         !!!stateNew % sup_thickness % array(1,:,:) = (stateOld % tracers % array(stateOld % index_temperature,:,:) * layerThicknessOld  + tracer_tendency(stateOld % index_temperature, :, :) * dt / SecondsInYear) / (layerThicknessNew+1.0e-12)


         layerThicknessNew = layerThicknessOld + layerThickness_tend * deltat
         thicknessNew = sum(layerThicknessNew, 1)


         !Optionally print some information about the new thickness
         !print *, 'thickness_tend maxval:', maxval(thickness_tend(1:mesh % nCellsSolve))       
         !print *, 'thicknessOld maxval:', maxval(thicknessOld(1:mesh % nCellsSolve))
         !print *, '  thicknessNew maxval:', maxval(thicknessNew(1:mesh % nCellsSolve))
         allocate( masktmp(nCells + 1) )
         masktmp = 0

!!!         ! if holding advance within initial extent of ice, set thickness to 0 anywhere it has expanded beyond initial extent
!!!         if (config_allow_additional_advance .eqv. .false.) then 
!!!             where ( MASK_WAS_INITIALLY_NOT_ICE(cellMaskOld) )
!!!                 thicknessNew = 0.0_RKIND
!!!             end where
!!!         endif

         ! reset negative thickness to 0.  This should not happen unless negative MB is larger than entire ice column.
         where (thicknessNew < 0.0_RKIND)
            masktmp = 1
            thicknessNew = 0.0_RKIND
!!!            stateNew % iceArea % array = 0.0_RKIND
         end where

         if (config_print_thickness_advection_info) then
            if (sum(masktmp) > 0) then
               write(6,*) '  Cells with negative thickness (set to 0):',sum(masktmp)
            endif

            ! Note how many cells have ice.
            masktmp = 0
            where (thicknessNew > 0.0_RKIND)
               masktmp = 1
            end where
            write(6,*) '  Cells with nonzero thickness:', sum(masktmp)
         endif
         deallocate(masktmp)



!!!         ! Calculate new tracer values =================
!!!         if (config_tracer_advection .ne. 'None') then
!!!           do iTracer = 1, size(tracersNew, 1)
!!!               where (layerThicknessNew > 0.0_RKIND)
!!!                   tracersNew(iTracer,:,:) = (tracersOld(iTracer,:,:) * layerThicknessOld &
!!!                       + tracer_tendency(iTracer,:,:) * dt / SecondsInYear)  /  (layerThicknessNew)
!!!               elsewhere
!!!                   ! May or may not want to assign tracer values to non-ice cells
!!!                   tracersNew(iTracer,:,:) = 0.0_RKIND
!!!               end where
!!!           end do
!!!         endif

!!!         ! Apply calving after we have updated the new state  - TODO Is this the right place?
!!!         call li_apply_calving(mesh, stateNew, err)


         block => block % next
      end do


      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in update_prognostics."
      endif


   !--------------------------------------------------------------------
   end subroutine update_prognostics



end module li_time_integration_fe


