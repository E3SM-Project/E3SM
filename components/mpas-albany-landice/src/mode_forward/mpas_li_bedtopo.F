! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_bedtopo
!
!> \MPAS land-ice bedtopo driver
!> \author Matt Hoffman and Holly Han (modified)
!> \date   20 June 2019, January 2022 (modified)
!> \details
!>  This module contains the routines for
!>  bed topography for solid earth changes
!
!-----------------------------------------------------------------------

module li_bedtopo

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log
   use mpas_timer
   use li_mask
   use li_setup
   use netcdf
#ifdef _MPI
   use mpi
#endif
#ifdef USE_FASTISOSTASY
   use fastisostasy, only: isos_class, bsl_class, isos_init, isos_init_ref, isos_init_state, &
                           isos_update, isos_end, bsl_init, bsl_update, &
                           isos_restart_write, bsl_restart_write, bsl_restart_read
   use isostasy_defs, only: wp
#endif

   implicit none
   private

   !--------------------------------------------------------------------
   ! Public parameters
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! Public member functions
   !--------------------------------------------------------------------

   public :: li_bedtopo_init, &
             li_bedtopo_finalize, &
             li_bedtopo_block_init, &
             li_bedtopo_solve, &
             li_bedtopo_write_restart

   !--------------------------------------------------------------------
   ! Private module variables
   !--------------------------------------------------------------------

   ! sea-level model timestep
   integer, save :: slmTimeStep, nglv_slm, norder_slm

   ! FastIsostasy timestep and solver dimensions
   integer, save :: fastisostasyTimeStep, nFastIsostasyX, nFastIsostasyY
#ifdef USE_FASTISOSTASY
   type (isos_class), save :: isosFastIsostasy
   type (bsl_class), save :: bslFastIsostasy
#endif

   ! Interpolation weights variables
   integer, dimension(:), allocatable :: toRowValues, toColValues
   integer, dimension(:), allocatable :: fromRowValues, fromColValues
   real, dimension(:), allocatable :: toSValues, fromSValues
   integer:: nMpas, nGrid

   ! MPI variables
   integer :: nCellsGlobal
   integer, dimension(:), allocatable :: nCellsDisplacement
   integer, dimension(:), allocatable :: indexToCellIDGathered
   integer, dimension(:), allocatable :: nCellsPerProc
   integer, pointer :: nCellsAll
   integer, pointer ::  nCellsOwned
   integer :: iCell, ilm, curProc
   real (kind=RKIND), dimension(:), allocatable :: globalArrayThickness
   real (kind=RKIND), dimension(:), allocatable :: gatheredArrayThickness
   real (kind=RKIND), dimension(:), allocatable :: globalArrayBedTopography
   real (kind=RKIND), dimension(:), allocatable :: gatheredArrayBedTopography
   real (kind=RKIND), dimension(:), allocatable :: globalArrayTopoChange
   real (kind=RKIND), dimension(:), allocatable :: gatheredArrayTopoChange
   real (kind=RKIND), dimension(:), allocatable :: globalArrayGeoidChange
   real (kind=RKIND), dimension(:), allocatable :: gatheredArrayGeoidChange

!***********************************************************************

contains

!***********************************************************************
!
!  routine li_bedtopo_init
!
!> \brief   Initializes bedtopo solver
!> \author  Matt Hoffman and Holly Han (modified)
!> \date    20 June 2019 (original), December 2021 (modified)
!> \details
!>  This routine initializes the bedtopo solver.
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_init(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

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

      character (len=StrKIND), pointer :: config_uplift_method
      integer :: err_tmp

      ! No init is needed.
      err = 0

      call mpas_pool_get_config(liConfigs, 'config_uplift_method', config_uplift_method)
      if (trim(config_uplift_method)=='sealevelmodel') then
         ! initialize the 1D sea-level model if fresh start
         call mpas_timer_start("sea-level model init")
         call slmodel_init(domain, err_tmp)
         call mpas_timer_stop("sea-level model init")
         err = ior(err, err_tmp)

         if (err /= 0) then
            call mpas_log_write('Error in li_bedtopo_init', MPAS_LOG_ERR)
         endif
      elseif (trim(config_uplift_method)=='fastisostasy') then
#ifdef USE_FASTISOSTASY
         call mpas_timer_start("FastIsostasy init")
         call fastisostasy_init(domain, err_tmp)
         call mpas_timer_stop("FastIsostasy init")
         err = ior(err, err_tmp)

         if (err /= 0) then
            call mpas_log_write('Error in li_bedtopo_init for FastIsostasy', MPAS_LOG_ERR)
         endif
#else
         call mpas_log_write("FastIsostasy needs to be included in the compilation with 'FASTISOSTASY=true'", &
              MPAS_LOG_ERR)
         err = ior(err,1)
#endif
      endif

   !--------------------------------------------------------------------

   end subroutine li_bedtopo_init



!***********************************************************************
!
!  routine li_bedtopo_block_init
!
!> \brief   Initializes blocks for bedtopo solver
!> \author  Matt Hoffman
!> \date    20 June 2019
!> \details
!>  This routine initializes each block of the bedtopo solver.
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_block_init(block, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: &
         block          !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! No block init needed.
      err = 0

   !--------------------------------------------------------------------
   end subroutine li_bedtopo_block_init



!***********************************************************************
!
!  subroutine li_bedtopo_solve
!
!> \brief   Updates bed topography
!> \author  Matt Hoffman and Holly Han (modified)
!> \date    20 June 2019 (original), December 2021 (modified)
!> \details
!>  This routine updates the bed topography.  Currently the only option
!>  is a data field passed in as input.
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_solve(domain, err)

      use mpas_timekeeping
      use li_mask
      use li_advection, only: li_update_geometry

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block
      character (len=StrKIND), pointer :: config_uplift_method
      type (mpas_pool_type), pointer :: meshPool          !< mesh information
      type (mpas_pool_type), pointer :: geometryPool      !< geometry information
      type (mpas_pool_type), pointer :: velocityPool      !< velocity information

      real (kind=RKIND), dimension(:), pointer :: bedTopography, upliftRate
      real (kind=RKIND), pointer :: deltat
      integer :: err_tmp

      err = 0
      err_tmp = 0

      ! Set needed variables and pointers
      call mpas_pool_get_config(liConfigs, 'config_uplift_method', config_uplift_method)
      if (trim(config_uplift_method)=='none') then
         ! do nothing

      elseif (trim(config_uplift_method)=='data') then

         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)

            call mpas_pool_get_array(meshPool, 'deltat', deltat)
            call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
            call mpas_pool_get_array(geometryPool, 'upliftRate', upliftRate)

            bedTopography(:) = bedTopography(:) + upliftRate(:) * deltat

            call li_update_geometry(geometryPool)
            call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)

            block => block % next
         end do

      elseif (trim(config_uplift_method)=='sealevelmodel') then
#ifdef USE_SEALEVELMODEL
         if (mpas_is_alarm_ringing(domain % clock, 'slmCouplingInterval', ierr=err_tmp)) then
            err = ior(err, err_tmp)

            slmTimeStep = slmTimeStep + 1

            call mpas_log_write("Calling the SLM. SLM timestep $i", intArgs=(/slmTimeStep/))
            call mpas_timer_start("sea-level model solve")
            call slmodel_solve(slmTimeStep, domain)
            call mpas_timer_stop("sea-level model solve")

            call mpas_reset_clock_alarm(domain % clock, 'slmCouplingInterval', ierr=err_tmp)
            err = ior(err, err_tmp)
         else
            ! do nothing for now, but could calculate uplift rate here later instead.
         endif
#else
         call mpas_log_write("The sea-level model needs to be included in the compilation with 'SLM=true'", &
              MPAS_LOG_ERR)
         err = ior(err,1)
#endif
      elseif (trim(config_uplift_method)=='fastisostasy') then
#ifdef USE_FASTISOSTASY
         if (mpas_is_alarm_ringing(domain % clock, 'fastisostasyCouplingInterval', ierr=err_tmp)) then
            err = ior(err, err_tmp)

            fastisostasyTimeStep = fastisostasyTimeStep + 1

            call mpas_log_write("Calling FastIsostasy. FastIsostasy timestep $i", &
                                intArgs=(/fastisostasyTimeStep/))
            call mpas_timer_start("FastIsostasy solve")
            call fastisostasy_solve(fastisostasyTimeStep, domain)
            call mpas_timer_stop("FastIsostasy solve")

            call mpas_reset_clock_alarm(domain % clock, 'fastisostasyCouplingInterval', ierr=err_tmp)
            err = ior(err, err_tmp)
         else
            ! do nothing for now, but could calculate uplift rate here later instead.
         endif
#else
         call mpas_log_write("FastIsostasy needs to be included in the compilation with 'FASTISOSTASY=true'", &
              MPAS_LOG_ERR)
         err = ior(err,1)
#endif
      else
         call mpas_log_write("Unknown option selected for 'config_uplift_method'", MPAS_LOG_ERR)
      endif





      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_bedtopo_solve.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------

   end subroutine li_bedtopo_solve



!***********************************************************************
!
!  routine li_bedtopo_finalize
!
!> \brief   finalizes bedtopo solver
!> \author  Matt Hoffman
!> \date    20 June 2019
!> \details
!>  This routine finalizes the bedtopo solver.
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_finalize(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      character (len=StrKIND), pointer :: config_uplift_method

      err = 0

#ifdef USE_SEALEVELMODEL
      call mpas_pool_get_config(liConfigs, 'config_uplift_method', config_uplift_method)
      if (trim(config_uplift_method)=='sealevelmodel') then
         if (curProc.eq.0) then
            deallocate(toRowValues, toColValues, toSValues)
            deallocate(fromRowValues, fromColValues, fromSValues)
         endif
         deallocate(nCellsPerProc)
         deallocate(nCellsDisplacement)
         deallocate(indexToCellIDGathered)
      endif
#endif
#ifdef USE_FASTISOSTASY
      call mpas_pool_get_config(liConfigs, 'config_uplift_method', config_uplift_method)
      if (trim(config_uplift_method)=='fastisostasy') then
         if (curProc.eq.0) then
            deallocate(toRowValues, toColValues, toSValues)
            deallocate(fromRowValues, fromColValues, fromSValues)
            call isos_end(isosFastIsostasy)
         endif
         deallocate(nCellsPerProc)
         deallocate(nCellsDisplacement)
         deallocate(indexToCellIDGathered)
      endif
#endif
   !--------------------------------------------------------------------
   end subroutine li_bedtopo_finalize



   ! private subroutines

!***********************************************************************
!
!  routine slmodel_init
!
!> \brief   Initializes the sea-level model
!> \author  Holly Kyeore Han
!> \date    January 2022
!> \details
!>  This wrapper routine initializes the sea-level solver(Han et al., 2022, GMD,
!>  https://github.com/MALI-Dev/1DSeaLevelModel_FWTW)
!
!-----------------------------------------------------------------------

   subroutine slmodel_init(domain, err)

#ifdef USE_SEALEVELMODEL
      use mpas_timekeeping
      use sl_model_mod                    !< this is part of the SLM code
      use sl_io_mod                       !< this is part of the SLM code
#endif
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

#ifdef USE_SEALEVELMODEL
      logical, pointer :: config_do_restart
      type (mpas_pool_type), pointer :: meshPool     !< mesh information
      type (mpas_pool_type), pointer :: geometryPool
      real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography
      real (kind=RKIND), dimension(:), allocatable :: meshMask
      real (kind=RKIND), dimension(:,:), allocatable :: ismIceload, ismBedtopo, ismMask
      real (kind=RKIND), dimension(:), allocatable :: thicknessSLgrid1D
      real (kind=RKIND), dimension(:), allocatable :: bedtopoSLgrid1D
      real (kind=RKIND), dimension(:), allocatable :: maskSLgrid1D
      integer :: err_tmp
      integer :: unit_num_slm                     ! SLM variable
      integer :: itersl, dtime                    ! SLM variable
      real    :: starttime                        ! SLM variable
      integer, dimension(:), pointer :: cellMask  ! integer bitmask for cells
      integer, pointer :: config_slm_coupling_interval
      type (MPAS_TimeInterval_type) :: slm_coupling_interval
      character (len=StrKIND), pointer :: simulationStartTime
      type (MPAS_Time_type) :: simulationStartTime_timeType

      ! MPI variables
      integer, dimension(:), pointer :: indexToCellID
      integer :: iProc, l, ilm, nProcs

      err = 0
      err_tmp = 0

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)

      ! Set up the alarm for the coupling time interval
      call mpas_pool_get_config(liConfigs, 'config_slm_coupling_interval', config_slm_coupling_interval)
      call mpas_set_timeInterval(slm_coupling_interval, YY=config_slm_coupling_interval, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_pool_get_array(meshPool, 'simulationStartTime', simulationStartTime)
      call mpas_set_time(simulationStartTime_timeType, dateTimeString=simulationStartTime, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_add_clock_alarm(domain%clock, 'slmCouplingInterval', alarmTime=simulationStartTime_timeType, &
           alarmTimeInterval=slm_coupling_interval, ierr=err_tmp)
      err = ior(err,err_tmp)
      if (mpas_is_alarm_ringing(domain%clock, 'slmCouplingInterval', ierr=err_tmp)) then
         err = ior(err, err_tmp)
         call mpas_reset_clock_alarm(domain%clock, 'slmCouplingInterval', ierr=err_tmp)
         err = ior(err, err_tmp)
      endif

      ! initialize interpolation
      call interpolate_init(domain, err_tmp)
      err = ior(err, err_tmp)

      ! Set needed variables for using MPI
      call mpas_pool_get_dimension(meshPool, 'nCells', nCellsAll)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsOwned)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

      ! Allocate globalArray and gatheredArray only on process 0
      call MPI_COMM_RANK( domain % dminfo % comm, curProc, err_tmp)
      err = ior(err, err_tmp)
      call MPI_COMM_SIZE( domain % dminfo % comm, nProcs, err_tmp)
      err = ior(err, err_tmp)

      ! perform the initialization on the head processor
      allocate(nCellsPerProc(nProcs))
      allocate(nCellsDisplacement(nProcs))

      ! Gather nCellsOwned
      call MPI_GATHER( nCellsOwned, 1, MPI_INTEGER, nCellsPerProc, 1, MPI_INTEGER, &
      0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)

      if (curProc.eq.0) then
         ! First, set SLM unit number to the MALI output log file unit
         unit_num_slm = domain % logInfo % outputLog % unitNum
         call sl_set_unit_num(unit_num_slm)

         call sl_drive_readnl(itersl, dtime, starttime) !SLM subroutine

         ! First, check consistency in coupling interval set up in MALI and SLM
         call check_SLM_coupling_interval(dtime, meshPool, err_tmp)
         err = ior(err, err_tmp)
         if (err /= 0) then
            call mpas_log_write("Error occurred in check_SLM_coupling_interval.", MPAS_LOG_ERR)
         endif

         ! check the number of points on the SLM Gauss-Legendre grid
         call sl_get_model_res(norder_slm, nglv_slm)
         call mpas_log_write("Solving the SLM upto $i order and degrees of Spherical Harmonics", intArgs=(/norder_slm/))
         call mpas_log_write("Number of Gauss-Legendre latitudinal nodes in the SLM grid is $i", intArgs=(/nglv_slm/))

         ! Set Displacement variable for GATHERV command
         nCellsGlobal = sum(nCellsPerProc)
         allocate(indexToCellIDGathered(nCellsGlobal))
         nCellsDisplacement(1) = 0
         if (nProcs > 1) then
            do iProc=2,nProcs
               nCellsDisplacement(iProc) = nCellsDisplacement(iProc-1) + nCellsPerProc(iProc-1)
            enddo
         endif
      else
         ! Intel requires this be allocated even though it is not meaningful on the non-destination procs
         allocate(indexToCellIDGathered(1))
      endif

      ! Gather indexToCellID
      call MPI_GATHERV( indexToCellID, nCellsOwned, MPI_INTEGER, indexToCellIDGathered, &
      nCellsPerProc, nCellsDisplacement, MPI_INTEGER, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)

      ! check if the run is a restart or not
      call mpas_pool_get_config(liConfigs, 'config_do_restart', config_do_restart)
      if (config_do_restart) then

         call find_slm_restart_timestep(meshPool, slmTimeStep, err_tmp)
         err = ior(err, err_tmp)

         ! Note: no need to call slmodel_solve on init of a restart.
         ! If this time level happens to be a coupling interval, SLM would have been solved already
         ! in the previous run that generated the restart file.
         ! While it would not hurt (other than unneeded execution time) to call SLM again if the
         ! restart time level happens to be a coupling interval, if the restart time is in between
         ! coupling intervals calling SLM here will make things out of sync.

      else

         call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         
         if (curProc.eq.0) then
            ! allocate and initialize the SLM-related variables used in this module
            allocate(globalArrayThickness(nCellsGlobal), gatheredArrayThickness(nCellsGlobal))
            allocate(globalArrayBedTopography(nCellsGlobal), gatheredArrayBedTopography(nCellsGlobal))
            allocate(meshMask(nCellsGlobal))
            allocate(ismIceload(nglv_slm,2*nglv_slm), ismBedtopo(nglv_slm,2*nglv_slm), ismMask(nglv_slm,2*nglv_slm))
            allocate(bedtopoSLgrid1D(nglv_slm*2*nglv_slm), thicknessSLgrid1D(nglv_slm*2*nglv_slm), maskSLgrid1D(nglv_slm*2*nglv_slm))
            ismIceload(:,:) = 0.0
            ismBedtopo(:,:) = 0.0
            ismMask(:,:) = 0.0
            bedtopoSLgrid1D(:) = 0.0
            thicknessSLgrid1D(:) = 0.0
            maskSLgrid1D(:) = 0.0
            call mpas_log_write("Allocating SLM variables in MALI, nglv is $i", intArgs=(/nglv_slm/))
         else
            ! Intel requires these be allocated even though they are not meaningful on the non-destination procs
            allocate(globalArrayThickness(1), gatheredArrayThickness(1))
            allocate(globalArrayBedTopography(1), gatheredArrayBedTopography(1))
            allocate(meshMask(1))
         endif

         ! Gather only the nCellsOwned from thickness and bedtopo (does not include Halos)
         call MPI_GATHERV((thickness*real(li_mask_is_grounded_ice_int(cellMask),RKIND)), &
                nCellsOwned, MPI_DOUBLE, gatheredArrayThickness, nCellsPerProc, &
                nCellsDisplacement, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
         err = ior(err, err_tmp)
         call MPI_GATHERV(bedTopography, nCellsOwned, MPI_DOUBLE, gatheredArrayBedTopography, nCellsPerProc, &
                          nCellsDisplacement, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
         err = ior(err, err_tmp)

         if (curProc.eq.0) then
            ! Rearrange data into CellID order
            do iCell = 1,nCellsGlobal
               globalArrayThickness(indexToCellIDGathered(iCell)) = gatheredArrayThickness(iCell)
               globalArrayBedTopography(indexToCellIDGathered(iCell)) = gatheredArrayBedTopography(iCell)
               meshMask(indexToCellIDGathered(iCell)) = 1
            enddo

            ! interpolate thickness, bedTopograpy, mesh mask to the Gaussian grid
            call interpolate(toColValues, toRowValues, toSvalues, globalArrayThickness, thicknessSLgrid1D)
            call interpolate(toColValues, toRowValues, toSvalues, globalArrayBedTopography, bedtopoSLgrid1D)
            call interpolate(toColValues, toRowValues, toSvalues, meshMask, maskSLgrid1D)

            ! reformat the interpolated data
            ismIceload = reshape(thicknessSLgrid1D, [nglv_slm,2*nglv_slm])
            ismBedtopo = reshape(bedtopoSLgrid1D, [nglv_slm,2*nglv_slm])
            ismMask = reshape(maskSLgrid1D, [nglv_slm,2*nglv_slm])

            ! initialize coupling time step number. initial time is 0
            slmTimeStep = 0

            ! series of calling SLM routines
            if (err == 0) then
               call sl_call_readnl
               call sl_timewindow(slmTimeStep)
               call sl_allocate_and_initialize_array
               call sl_solver_checkpoint(itersl, dtime)
               call sl_solver_init(itersl, starttime, ismIceload, ismBedtopo, ismMask)
               call sl_deallocate_array
            endif

            deallocate(ismIceload, ismBedtopo, ismMask)
            deallocate(bedtopoSLgrid1D, thicknessSLgrid1D, maskSLgrid1D)
         endif
         deallocate(globalArrayThickness)
         deallocate(gatheredArrayThickness)
         deallocate(globalArrayBedTopography)
         deallocate(gatheredArrayBedTopography)
         deallocate(meshMask)

      endif ! endif restart or not

# else
      call mpas_log_write("The sea-level model needs to be included in the compilation with 'SLM=true'", &
           MPAS_LOG_ERR)
      err = ior(err,1)
# endif

   !--------------------------------------------------------------------
   end subroutine slmodel_init



!***********************************************************************
!
!  routine slmodel_solve
!
!> \brief   Solves gravitationally consistent sea-level change
!> \author  Holly Kyeore Han
!> \date    January 2022
!> \details
!>  This wrapper routine calls the sea-level solver that takes in
!>  ice thickness and provides sea-level change (i.e., changes in the
!>  heights of the sea surface and the solid Earth surface associated
!>  with ice sheet changes. The sea-level model is taken and modified
!>  from Han et al. (2021, GMD, https://doi.org/10.5281/zenodo.5775235)
!
!-----------------------------------------------------------------------

   subroutine slmodel_solve(slmTimeStep, domain)

      use li_advection, only: li_update_geometry
#ifdef USE_SEALEVELMODEL
      use sl_model_mod                    !< this is part of the SLM code
      use sl_io_mod                       !< this is part of the SLM code
#endif
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: slmTimeStep

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

#ifdef USE_SEALEVELMODEL
      type (mpas_pool_type), pointer :: meshPool          !< mesh information
      type (mpas_pool_type), pointer :: geometryPool      !< geometry information
      type (mpas_pool_type), pointer :: velocityPool      !< velocity information

      real (kind=RKIND), dimension(:), pointer :: bedTopography, thickness
      real (kind=RKIND), dimension(:), pointer :: bedTopographyChange
      real (kind=RKIND), dimension(:), pointer :: geoidChange
      real (kind=RKIND), dimension(:), allocatable :: meshMask
      real (kind=RKIND), dimension(:,:), allocatable :: ismIceload, ismMask
      real (kind=RKIND), dimension(:,:), allocatable :: slmSLchange
      real (kind=RKIND), dimension(:,:), allocatable :: slmGeoidChange
      real (kind=RKIND), dimension(:), allocatable :: slChangeSLgrid1D
      real (kind=RKIND), dimension(:), allocatable :: geoidChangeSLgrid1D
      real (kind=RKIND), dimension(:), allocatable :: thicknessSLgrid1D
      real (kind=RKIND), dimension(:), allocatable :: maskSLgrid1D
      integer, dimension(:), pointer :: cellMask

      integer :: err, err_tmp

      integer :: itersl, dtime                 ! SLM variable
      real    :: starttime                     ! SLM variable

      err = 0
      err_tmp = 0

      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'bedTopographyChange', bedTopographyChange)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'geoidChange', geoidChange)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'velocity', velocityPool)

      ! Allocate globalArray and gatheredArray only on process 0
      call MPI_COMM_RANK(domain % dminfo % comm, curProc, err_tmp)
      err = ior(err, err_tmp)

      if (curProc.eq.0) then
         ! call the SLM modules to initialize arrays and their dimensions
         call sl_drive_readnl(itersl, dtime, starttime)
         ! allocate and initialize the SLM-related variables used in this module
         allocate(globalArrayThickness(nCellsGlobal), gatheredArrayThickness(nCellsGlobal))
         allocate(globalArrayTopoChange(nCellsGlobal), gatheredArrayTopoChange(nCellsGlobal))
         allocate(globalArrayGeoidChange(nCellsGlobal), gatheredArrayGeoidChange(nCellsGlobal))
         allocate(meshMask(nCellsGlobal))
         allocate(ismIceload(nglv_slm,2*nglv_slm), ismMask(nglv_slm,2*nglv_slm))
         allocate(slmSLchange(nglv_slm,2*nglv_slm), slmGeoidChange(nglv_slm,2*nglv_slm))
         allocate(slChangeSLgrid1D(nglv_slm*2*nglv_slm), thicknessSLgrid1D(nglv_slm*2*nglv_slm), maskSLgrid1D(nglv_slm*2*nglv_slm))
         allocate(geoidChangeSLgrid1D(nglv_slm*2*nglv_slm))
         ismIceload(:,:) = 0.0_RKIND
         ismMask(:,:) = 0.0_RKIND
         slmSLchange(:,:) = 0.0_RKIND
         slmGeoidChange(:,:) = 0.0_RKIND
         slChangeSLgrid1D(:) = 0.0_RKIND
         geoidChangeSLgrid1D(:) = 0.0_RKIND
         thicknessSLgrid1D(:) = 0.0_RKIND
         maskSLgrid1D(:) = 0.0_RKIND
      else
         ! Intel requires these be allocated even though they are not meaningful on the non-destination procs
         allocate(globalArrayThickness(1), gatheredArrayThickness(1))
         allocate(globalArrayTopoChange(1), gatheredArrayTopoChange(1))
         allocate(globalArrayGeoidChange(1), gatheredArrayGeoidChange(1))
         allocate(meshMask(1))
      endif

      ! Gather only the nCellsOwned from ice thickness (does not include Halos)
      call MPI_GATHERV((thickness*real(li_mask_is_grounded_ice_int(cellMask),RKIND)), &
                       nCellsOwned, MPI_DOUBLE, gatheredArrayThickness, nCellsPerProc, &
                       nCellsDisplacement, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)

      if (curProc.eq.0) then

         ! Rearrange thickness into CellID order
         do iCell = 1,nCellsGlobal
            globalArrayThickness(indexToCellIDGathered(iCell)) = gatheredArrayThickness(iCell)
            meshMask(indexToCellIDGathered(iCell)) = 1
         enddo

         ! interpolate thickness to Gaussian grid
         call interpolate(toColValues, toRowValues, toSvalues, globalArrayThickness, thicknessSLgrid1D)
         call interpolate(toColValues, toRowValues, toSvalues, meshMask, maskSLgrid1D)

         ! reformat the interpolated data
         ismIceload = reshape(thicknessSLgrid1D, [nglv_slm,2*nglv_slm])
         ismMask = reshape(maskSLgrid1D, [nglv_slm,2*nglv_slm])

         ! series of calling SLM routines
         call sl_call_readnl
         call sl_timewindow(slmTimeStep)
         call sl_allocate_and_initialize_array
         call sl_solver_checkpoint(itersl, dtime)
         call sl_solver(itersl, slmTimeStep, dtime, starttime, ismIceload, ismMask, slmSLchange, slmGeoidChange)
         call sl_deallocate_array

         ! reshape 2D array SLM output into 1D array
         slChangeSLgrid1D = reshape(slmSLchange, [nglv_slm*2*nglv_slm])
         geoidChangeSLgrid1D = reshape(slmGeoidchange, [nglv_slm*2*nglv_slm])

         ! interpolate sea-level change from GL grid to MALI mesh.
         ! note: in the static sea-level theory, sea level and topography are globally defined !>
         ! and negative of each other. That is, topography change is negative of sea-level change
         call interpolate(fromColValues, fromRowValues, fromSValues, -1.0_RKIND*(slChangeSLgrid1D), &
                          globalArrayTopoChange)
         call interpolate(fromColValues, fromRowValues, fromSValues, geoidChangeSLgrid1D, &
                          globalArrayGeoidChange)

         ! Rearrange back to index order
         do iCell = 1,nCellsGlobal
            gatheredArrayTopoChange(iCell) = globalArrayTopoChange(indexToCellIDGathered(iCell))
            gatheredArrayGeoidChange(iCell) = globalArrayGeoidChange(indexToCellIDGathered(iCell))
         enddo

         deallocate(ismIceload, ismMask, slmSLchange, slmGeoidChange)
         deallocate(slChangeSLgrid1D, geoidChangeSLgrid1D, thicknessSLgrid1D, maskSLgrid1D)
      endif

      ! scatter output sea-level changes to processors
      call MPI_SCATTERV(gatheredArrayTopoChange, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
                        bedTopographyChange, nCellsAll, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)
      call MPI_SCATTERV(gatheredArrayGeoidChange, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
                        geoidChange, nCellsAll, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)

      ! update bedTopography
      bedTopography(:) = bedTopography(:) + bedTopographyChange(:)

      ! Perform Halo exchange update
      ! Don't bother with geoidChange because it is only a diagnostic output
      call mpas_dmpar_field_halo_exch(domain,'bedTopography')
      call li_update_geometry(geometryPool)
      call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
      err = ior(err, err_tmp)

      ! deallocate memory
      deallocate(globalArrayThickness)
      deallocate(gatheredArrayThickness)
      deallocate(globalArrayTopoChange)
      deallocate(gatheredArrayTopoChange)
      deallocate(globalArrayGeoidChange)
      deallocate(gatheredArrayGeoidChange)
      deallocate(meshMask)
#endif

   !--------------------------------------------------------------------
   end subroutine slmodel_solve


!***********************************************************************
!
!  routine fastisostasy_init
!
!> \brief   Initializes FastIsostasy coupling
!> \author  MALI developers
!> \date    June 2026
!> \details
!>  This wrapper routine initializes FastIsostasy and the coupling
!>  infrastructure (alarm, interpolation weights, MPI gather/scatter maps).
!
!-----------------------------------------------------------------------

   subroutine fastisostasy_init(domain, err)

#ifdef USE_FASTISOSTASY
      use mpas_timekeeping
#endif
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err
      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

#ifdef USE_FASTISOSTASY
      logical, pointer :: config_do_restart
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography
      real (kind=RKIND), dimension(:), allocatable :: meshMask
      real (kind=RKIND), dimension(:), allocatable :: thicknessFastIsostasy1D
      real (kind=RKIND), dimension(:), allocatable :: bedtopoFastIsostasy1D
      real (kind=wp), dimension(:,:), allocatable :: thicknessFastIsostasy2D
      real (kind=wp), dimension(:,:), allocatable :: bedtopoFastIsostasy2D
      integer, dimension(:), pointer :: cellMask
      integer :: err_tmp
      integer, pointer :: config_fastisostasy_coupling_interval
      type (MPAS_TimeInterval_type) :: fastisostasy_coupling_interval
      character (len=StrKIND), pointer :: simulationStartTime
      type (MPAS_Time_type) :: simulationStartTime_timeType
      character (len=StrKIND), pointer :: config_fastisostasy_parameter_file
      character (len=StrKIND), pointer :: config_MALI_to_FASTISOSTASY_weights_file
      character (len=StrKIND), pointer :: config_FASTISOSTASY_to_MALI_weights_file
      real (kind=RKIND), pointer :: config_fastisostasy_grid_resolution
      character (len=StrKIND), pointer :: config_fastisostasy_restart_file
      character (len=StrKIND), pointer :: config_fastisostasy_bsl_restart_file
      character (len=StrKIND) :: restart_filename, bsl_restart_filename
      character (len=StrKIND), pointer :: xtime
      logical :: fastisostasyRestartFileExists, fastisostasyBslRestartFileExists
      real (kind=wp) :: fastisostasyInitTime

      ! MPI variables
      integer, dimension(:), pointer :: indexToCellID
      integer :: iProc, nProcs
      real (kind=wp) :: fastisostasyDx, fastisostasyDy

      err = 0
      err_tmp = 0

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)

      ! Set up the alarm for the coupling time interval
      call mpas_pool_get_config(liConfigs, 'config_fastisostasy_coupling_interval', config_fastisostasy_coupling_interval)
      call mpas_set_timeInterval(fastisostasy_coupling_interval, YY=config_fastisostasy_coupling_interval, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_pool_get_array(meshPool, 'simulationStartTime', simulationStartTime)
      call mpas_set_time(simulationStartTime_timeType, dateTimeString=simulationStartTime, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_add_clock_alarm(domain%clock, 'fastisostasyCouplingInterval', alarmTime=simulationStartTime_timeType, &
           alarmTimeInterval=fastisostasy_coupling_interval, ierr=err_tmp)
      err = ior(err,err_tmp)
      if (mpas_is_alarm_ringing(domain%clock, 'fastisostasyCouplingInterval', ierr=err_tmp)) then
         err = ior(err, err_tmp)
         call mpas_reset_clock_alarm(domain%clock, 'fastisostasyCouplingInterval', ierr=err_tmp)
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(liConfigs, 'config_MALI_to_FASTISOSTASY_weights_file', config_MALI_to_FASTISOSTASY_weights_file)
      call mpas_pool_get_config(liConfigs, 'config_FASTISOSTASY_to_MALI_weights_file', config_FASTISOSTASY_to_MALI_weights_file)
      call mpas_pool_get_config(liConfigs, 'config_fastisostasy_grid_resolution', config_fastisostasy_grid_resolution)
      call mpas_pool_get_config(liConfigs, 'config_fastisostasy_parameter_file', config_fastisostasy_parameter_file)

      if (config_fastisostasy_grid_resolution <= 0.0_RKIND) then
         call mpas_log_write("config_fastisostasy_grid_resolution must be set to the physical grid " // &
              "spacing (in meters) of the regular grid defined by config_MALI_to_FASTISOSTASY_weights_file " // &
              "and config_FASTISOSTASY_to_MALI_weights_file. It has no valid default.", MPAS_LOG_ERR)
         err = ior(err, 1)
      endif

      ! initialize interpolation
      call interpolate_init(domain, err_tmp, trim(config_MALI_to_FASTISOSTASY_weights_file), &
                            trim(config_FASTISOSTASY_to_MALI_weights_file))
      err = ior(err, err_tmp)

      ! Set needed variables for using MPI
      call mpas_pool_get_dimension(meshPool, 'nCells', nCellsAll)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsOwned)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

      call MPI_COMM_RANK( domain % dminfo % comm, curProc, err_tmp)
      err = ior(err, err_tmp)
      call MPI_COMM_SIZE( domain % dminfo % comm, nProcs, err_tmp)
      err = ior(err, err_tmp)

      allocate(nCellsPerProc(nProcs))
      allocate(nCellsDisplacement(nProcs))

      call MPI_GATHER( nCellsOwned, 1, MPI_INTEGER, nCellsPerProc, 1, MPI_INTEGER, &
      0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)

      if (curProc.eq.0) then
         nCellsGlobal = sum(nCellsPerProc)
         allocate(indexToCellIDGathered(nCellsGlobal))
         nCellsDisplacement(1) = 0
         if (nProcs > 1) then
            do iProc=2,nProcs
               nCellsDisplacement(iProc) = nCellsDisplacement(iProc-1) + nCellsPerProc(iProc-1)
            enddo
         endif
      else
         ! Intel requires this be allocated even though it is not meaningful on the non-destination procs
         allocate(indexToCellIDGathered(1))
      endif

      call MPI_GATHERV( indexToCellID, nCellsOwned, MPI_INTEGER, indexToCellIDGathered, &
      nCellsPerProc, nCellsDisplacement, MPI_INTEGER, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)

      call mpas_pool_get_config(liConfigs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(liConfigs, 'config_fastisostasy_restart_file', config_fastisostasy_restart_file)
      call mpas_pool_get_config(liConfigs, 'config_fastisostasy_bsl_restart_file', config_fastisostasy_bsl_restart_file)

      ! Determine FastIsostasy's own coupling-step counter and internal clock.
      if (config_do_restart) then
         call find_fastisostasy_restart_timestep(meshPool, fastisostasyTimeStep, err_tmp)
         err = ior(err, err_tmp)

         if (curProc.eq.0) then
            ! Construct timestamped filenames (same as what write path creates)
            call construct_fastisostasy_restart_filename(domain, config_fastisostasy_restart_file, &
                                                         restart_filename, err_tmp)
            err = ior(err, err_tmp)
            call construct_fastisostasy_restart_filename(domain, config_fastisostasy_bsl_restart_file, &
                                                         bsl_restart_filename, err_tmp)
            err = ior(err, err_tmp)

            ! Check if timestamped restart files exist
            if (err == 0) then
               inquire(file=trim(restart_filename), exist=fastisostasyRestartFileExists)
               inquire(file=trim(bsl_restart_filename), exist=fastisostasyBslRestartFileExists)

               if (.not. fastisostasyRestartFileExists) then
                  call mpas_log_write("FastIsostasy restart file not found: " // trim(restart_filename), &
                                      MPAS_LOG_ERR)
                  call mpas_pool_get_array(meshPool, 'xtime', xtime)
                  call mpas_log_write("  Current model time: " // trim(xtime), MPAS_LOG_ERR)
                  call mpas_log_write("  Filename constructed from config_fastisostasy_restart_file (" // &
                                      trim(config_fastisostasy_restart_file) // ") and current model time.", &
                                      MPAS_LOG_ERR)
                  call mpas_log_write("  This should match a MALI restart file at the same time.", MPAS_LOG_ERR)
                  call mpas_log_write("A restart with config_uplift_method='fastisostasy' requires matching " // &
                                      "FI restart files.", MPAS_LOG_ERR)
                  err = ior(err, 1)
               endif

               if (.not. fastisostasyBslRestartFileExists) then
                  call mpas_log_write("FastIsostasy BSL restart file not found: " // trim(bsl_restart_filename), &
                                      MPAS_LOG_ERR)
                  call mpas_pool_get_array(meshPool, 'xtime', xtime)
                  call mpas_log_write("  Current model time: " // trim(xtime), MPAS_LOG_ERR)
                  call mpas_log_write("  Filename constructed from config_fastisostasy_bsl_restart_file (" // &
                                      trim(config_fastisostasy_bsl_restart_file) // ") and current model time.", &
                                      MPAS_LOG_ERR)
                  call mpas_log_write("  This should match a MALI restart file at the same time.", MPAS_LOG_ERR)
                  call mpas_log_write("A restart with config_uplift_method='fastisostasy' requires matching " // &
                                      "FI restart files.", MPAS_LOG_ERR)
                  err = ior(err, 1)
               endif

               if (err == 0) then
                  call mpas_log_write("Reading FastIsostasy restart files:")
                  call mpas_log_write("  " // trim(restart_filename))
                  call mpas_log_write("  " // trim(bsl_restart_filename))
               endif
            endif
         endif

         ! Broadcast error status to all ranks
         call MPI_BCAST(err, 1, MPI_INTEGER, 0, domain % dminfo % comm, err_tmp)
      else
         fastisostasyTimeStep = 0
      endif
      fastisostasyInitTime = real(fastisostasyTimeStep, wp) * real(config_fastisostasy_coupling_interval, wp)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

      if (curProc.eq.0) then
         allocate(globalArrayThickness(nCellsGlobal), gatheredArrayThickness(nCellsGlobal))
         allocate(globalArrayBedTopography(nCellsGlobal), gatheredArrayBedTopography(nCellsGlobal))
         allocate(meshMask(nCellsGlobal))
         allocate(thicknessFastIsostasy1D(nGrid), bedtopoFastIsostasy1D(nGrid))
         meshMask(:) = 0.0_RKIND
         thicknessFastIsostasy1D(:) = 0.0_RKIND
         bedtopoFastIsostasy1D(:) = 0.0_RKIND
      else
         allocate(globalArrayThickness(1), gatheredArrayThickness(1))
         allocate(globalArrayBedTopography(1), gatheredArrayBedTopography(1))
         allocate(meshMask(1))
      endif

      call MPI_GATHERV((thickness*real(li_mask_is_grounded_ice_int(cellMask),RKIND)), &
             nCellsOwned, MPI_DOUBLE, gatheredArrayThickness, nCellsPerProc, &
             nCellsDisplacement, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)
      call MPI_GATHERV(bedTopography, nCellsOwned, MPI_DOUBLE, gatheredArrayBedTopography, nCellsPerProc, &
                       nCellsDisplacement, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)

      if (curProc.eq.0) then
         do iCell = 1,nCellsGlobal
            globalArrayThickness(indexToCellIDGathered(iCell)) = gatheredArrayThickness(iCell)
            globalArrayBedTopography(indexToCellIDGathered(iCell)) = gatheredArrayBedTopography(iCell)
            meshMask(indexToCellIDGathered(iCell)) = 1.0_RKIND
         enddo

         call interpolate(toColValues, toRowValues, toSvalues, globalArrayThickness, thicknessFastIsostasy1D)
         call interpolate(toColValues, toRowValues, toSvalues, globalArrayBedTopography, bedtopoFastIsostasy1D)

         nFastIsostasyX = int(nint(sqrt(real(nGrid, RKIND))))
         nFastIsostasyY = nFastIsostasyX
         if (nFastIsostasyX*nFastIsostasyY /= nGrid) then
            call mpas_log_write('FastIsostasy coupling currently expects a square target grid in the mapping file.', &
                                MPAS_LOG_ERR)
            err = ior(err, 1)
            return
         endif

         allocate(thicknessFastIsostasy2D(nFastIsostasyX,nFastIsostasyY))
         allocate(bedtopoFastIsostasy2D(nFastIsostasyX,nFastIsostasyY))
         thicknessFastIsostasy2D = real(reshape(thicknessFastIsostasy1D, [nFastIsostasyX,nFastIsostasyY]), wp)
         bedtopoFastIsostasy2D = real(reshape(bedtopoFastIsostasy1D, [nFastIsostasyX,nFastIsostasyY]), wp)

         ! FastIsostasy requires a square grid (isos_utils::init_domain_size aborts otherwise), so
         ! both dx and dy are set from the single config_fastisostasy_grid_resolution value.
         fastisostasyDx = real(config_fastisostasy_grid_resolution, wp)
         fastisostasyDy = real(config_fastisostasy_grid_resolution, wp)

         if (err == 0) then
            call bsl_init(bslFastIsostasy, trim(config_fastisostasy_parameter_file), fastisostasyInitTime)
            call isos_init(isosFastIsostasy, trim(config_fastisostasy_parameter_file), 'isostasy', &
                           nFastIsostasyX, nFastIsostasyY, fastisostasyDx, fastisostasyDy)

            if (config_do_restart) then
               ! Use timestamped restart filenames constructed earlier
               isosFastIsostasy%par%restart = trim(restart_filename)
               isosFastIsostasy%par%use_restart = .true.
               call bsl_restart_read(bslFastIsostasy, trim(bsl_restart_filename))
            else
               isosFastIsostasy%par%use_restart = .false.
               call isos_init_ref(isosFastIsostasy, bedtopoFastIsostasy2D, thicknessFastIsostasy2D)
            endif

            call isos_init_state(isosFastIsostasy, bedtopoFastIsostasy2D, thicknessFastIsostasy2D, &
                                 fastisostasyInitTime, bslFastIsostasy)
         endif

         deallocate(thicknessFastIsostasy2D, bedtopoFastIsostasy2D)
         deallocate(thicknessFastIsostasy1D, bedtopoFastIsostasy1D)
      endif

      deallocate(globalArrayThickness)
      deallocate(gatheredArrayThickness)
      deallocate(globalArrayBedTopography)
      deallocate(gatheredArrayBedTopography)
      deallocate(meshMask)

#else
      call mpas_log_write("FastIsostasy needs to be included in the compilation with 'FASTISOSTASY=true'", &
           MPAS_LOG_ERR)
      err = ior(err,1)
#endif

   end subroutine fastisostasy_init


!***********************************************************************
!
!  routine fastisostasy_solve
!
!> \brief   Solves FastIsostasy bedrock deformation coupling
!> \author  MALI developers
!> \date    June 2026
!> \details
!>  This wrapper routine gathers grounded ice thickness, calls FastIsostasy
!>  on rank 0, and scatters bed topography change back to MPAS ranks.
!
!-----------------------------------------------------------------------

   subroutine fastisostasy_solve(fastisostasyTimeStep, domain)

      use li_advection, only: li_update_geometry
#ifdef USE_FASTISOSTASY
      use isos_utils, only: in2out
#endif

      integer, intent(in) :: fastisostasyTimeStep
      type (domain_type), intent(inout) :: domain

#ifdef USE_FASTISOSTASY
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: velocityPool
      integer, pointer :: config_fastisostasy_coupling_interval

      real (kind=RKIND), dimension(:), pointer :: bedTopography, thickness
      real (kind=RKIND), dimension(:), pointer :: bedTopographyChange
      real (kind=RKIND), dimension(:), pointer :: geoidChange
      real (kind=RKIND), dimension(:), allocatable :: thicknessFastIsostasy1D
      real (kind=RKIND), dimension(:), allocatable :: bedtopoFastIsostasy1D
      real (kind=RKIND), dimension(:), allocatable :: topoChangeFastIsostasy1D
      real (kind=RKIND), dimension(:), allocatable :: geoidChangeFastIsostasy1D
      real (kind=RKIND), dimension(:), allocatable :: meshMask
      real (kind=wp), dimension(:,:), allocatable :: thicknessFastIsostasy2D
      real (kind=wp), dimension(:,:), allocatable :: bedtopoFastIsostasy2D
      real (kind=wp), dimension(:,:), allocatable :: bedtopoFastIsostasyUpdated2D
      real (kind=wp), dimension(:,:), allocatable :: zssRefFastIsostasy2D
      real (kind=wp), dimension(:,:), allocatable :: geoidChangeFastIsostasy2D
      real (kind=wp) :: fastisostasyTime
      integer, dimension(:), pointer :: cellMask
      integer :: err, err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'bedTopographyChange', bedTopographyChange)
      call mpas_pool_get_array(geometryPool, 'geoidChange', geoidChange)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'velocity', velocityPool)
      call mpas_pool_get_config(liConfigs, 'config_fastisostasy_coupling_interval', config_fastisostasy_coupling_interval)

      call MPI_COMM_RANK(domain % dminfo % comm, curProc, err_tmp)
      err = ior(err, err_tmp)

      if (curProc.eq.0) then
         allocate(globalArrayThickness(nCellsGlobal), gatheredArrayThickness(nCellsGlobal))
         allocate(globalArrayBedTopography(nCellsGlobal), gatheredArrayBedTopography(nCellsGlobal))
         allocate(globalArrayTopoChange(nCellsGlobal), gatheredArrayTopoChange(nCellsGlobal))
         allocate(globalArrayGeoidChange(nCellsGlobal), gatheredArrayGeoidChange(nCellsGlobal))
         allocate(meshMask(nCellsGlobal))
         allocate(thicknessFastIsostasy1D(nGrid), bedtopoFastIsostasy1D(nGrid))
         allocate(topoChangeFastIsostasy1D(nGrid))
         allocate(geoidChangeFastIsostasy1D(nGrid))
         allocate(thicknessFastIsostasy2D(nFastIsostasyX,nFastIsostasyY))
         allocate(bedtopoFastIsostasy2D(nFastIsostasyX,nFastIsostasyY))
         allocate(bedtopoFastIsostasyUpdated2D(nFastIsostasyX,nFastIsostasyY))
         allocate(zssRefFastIsostasy2D(nFastIsostasyX,nFastIsostasyY))
         allocate(geoidChangeFastIsostasy2D(nFastIsostasyX,nFastIsostasyY))
         meshMask(:) = 0.0_RKIND
         thicknessFastIsostasy1D(:) = 0.0_RKIND
         bedtopoFastIsostasy1D(:) = 0.0_RKIND
         topoChangeFastIsostasy1D(:) = 0.0_RKIND
         geoidChangeFastIsostasy1D(:) = 0.0_RKIND
      else
         allocate(globalArrayThickness(1), gatheredArrayThickness(1))
         allocate(globalArrayBedTopography(1), gatheredArrayBedTopography(1))
         allocate(globalArrayTopoChange(1), gatheredArrayTopoChange(1))
         allocate(globalArrayGeoidChange(1), gatheredArrayGeoidChange(1))
         allocate(meshMask(1))
      endif

      call MPI_GATHERV((thickness*real(li_mask_is_grounded_ice_int(cellMask),RKIND)), &
                       nCellsOwned, MPI_DOUBLE, gatheredArrayThickness, nCellsPerProc, &
                       nCellsDisplacement, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)

      call MPI_GATHERV(bedTopography, nCellsOwned, MPI_DOUBLE, gatheredArrayBedTopography, nCellsPerProc, &
                       nCellsDisplacement, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)

      if (curProc.eq.0) then
         do iCell = 1,nCellsGlobal
            globalArrayThickness(indexToCellIDGathered(iCell)) = gatheredArrayThickness(iCell)
            globalArrayBedTopography(indexToCellIDGathered(iCell)) = gatheredArrayBedTopography(iCell)
            meshMask(indexToCellIDGathered(iCell)) = 1.0_RKIND
         enddo

         call interpolate(toColValues, toRowValues, toSvalues, globalArrayThickness, thicknessFastIsostasy1D)
         call interpolate(toColValues, toRowValues, toSvalues, globalArrayBedTopography, bedtopoFastIsostasy1D)

         thicknessFastIsostasy2D = real(reshape(thicknessFastIsostasy1D, [nFastIsostasyX,nFastIsostasyY]), wp)
         bedtopoFastIsostasy2D = real(reshape(bedtopoFastIsostasy1D, [nFastIsostasyX,nFastIsostasyY]), wp)

         fastisostasyTime = real(fastisostasyTimeStep * config_fastisostasy_coupling_interval, wp)
         call bsl_update(bslFastIsostasy, fastisostasyTime)
         call isos_update(isosFastIsostasy, thicknessFastIsostasy2D, fastisostasyTime, bslFastIsostasy)
         ! Correct for sea-surface-height change
	 call in2out(zssRefFastIsostasy2D, isosFastIsostasy%ref%z_ss, isosFastIsostasy%domain)
         bedtopoFastIsostasyUpdated2D = isosFastIsostasy%out%z_bed - &
            (isosFastIsostasy%out%z_ss - zssRefFastIsostasy2D)

         ! Calculate geoid change
         geoidChangeFastIsostasy2D = isosFastIsostasy%out%z_ss - zssRefFastIsostasy2D

         topoChangeFastIsostasy1D = real(reshape(bedtopoFastIsostasyUpdated2D - bedtopoFastIsostasy2D, [nGrid]), RKIND)
         geoidChangeFastIsostasy1D = real(reshape(geoidChangeFastIsostasy2D, [nGrid]), RKIND)

         call interpolate(fromColValues, fromRowValues, fromSValues, topoChangeFastIsostasy1D, globalArrayTopoChange)
         call interpolate(fromColValues, fromRowValues, fromSValues, geoidChangeFastIsostasy1D, globalArrayGeoidChange)

         do iCell = 1,nCellsGlobal
            gatheredArrayTopoChange(iCell) = globalArrayTopoChange(indexToCellIDGathered(iCell))
            gatheredArrayGeoidChange(iCell) = globalArrayGeoidChange(indexToCellIDGathered(iCell))
         enddo

         deallocate(thicknessFastIsostasy1D, bedtopoFastIsostasy1D)
         deallocate(topoChangeFastIsostasy1D)
         deallocate(geoidChangeFastIsostasy1D)
         deallocate(thicknessFastIsostasy2D, bedtopoFastIsostasy2D, bedtopoFastIsostasyUpdated2D)
         deallocate(zssRefFastIsostasy2D, geoidChangeFastIsostasy2D)
      endif

      call MPI_SCATTERV(gatheredArrayTopoChange, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
                        bedTopographyChange, nCellsAll, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)
      call MPI_SCATTERV(gatheredArrayGeoidChange, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
                        geoidChange, nCellsAll, MPI_DOUBLE, 0, domain % dminfo % comm, err_tmp)
      err = ior(err, err_tmp)

      bedTopography(:) = bedTopography(:) + bedTopographyChange(:)

      ! Don't bother with geoidChange halo exchange because it is only a diagnostic output
      call mpas_dmpar_field_halo_exch(domain,'bedTopography')
      call li_update_geometry(geometryPool)
      call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
      err = ior(err, err_tmp)

      deallocate(globalArrayThickness)
      deallocate(gatheredArrayThickness)
      deallocate(globalArrayBedTopography)
      deallocate(gatheredArrayBedTopography)
      deallocate(globalArrayTopoChange)
      deallocate(gatheredArrayTopoChange)
      deallocate(globalArrayGeoidChange)
      deallocate(gatheredArrayGeoidChange)
      deallocate(meshMask)
#endif

   end subroutine fastisostasy_solve


!***********************************************************************
!
!  routine li_bedtopo_write_restart
!
!> \brief   Writes FastIsostasy's internal restart state
!> \author  MALI developers
!> \date    June 2026
!> \details
!>  FastIsostasy carries prognostic, history-dependent state (viscous
!>  bedrock displacement and its ODE integrator, the load reference
!>  state, and the prognostic barystatic sea level) that cannot be
!>  reconstructed from MALI's own restart fields. This routine snapshots
!>  that state to disk every time MALI writes its own restart stream, so
!>  a subsequent restart can resume the solid-Earth relaxation already in
!>  progress instead of silently re-initializing it as if already in
!>  equilibrium with the current ice/ocean load (see fastisostasy_init).
!>  This must be called only when MALI's 'restart' stream is written, and
!>  only by the task that owns isosFastIsostasy/bslFastIsostasy (rank 0).
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_write_restart(domain, err)

      use, intrinsic :: ieee_exceptions, only: ieee_all, ieee_get_halting_mode, &
                                               ieee_set_halting_mode, ieee_set_flag

      type (domain_type), intent(in) :: domain
      integer, intent(out) :: err

#ifdef USE_FASTISOSTASY
      character (len=StrKIND), pointer :: config_uplift_method
      character (len=StrKIND), pointer :: config_fastisostasy_restart_file
      character (len=StrKIND), pointer :: config_fastisostasy_bsl_restart_file
      integer, pointer :: config_fastisostasy_coupling_interval
      character (len=StrKIND) :: restart_filename, bsl_restart_filename
      real (kind=wp) :: fastisostasyTime
      integer :: err_tmp
      logical :: haltingMode(size(ieee_all))
#endif

      err = 0

#ifdef USE_FASTISOSTASY
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_uplift_method', config_uplift_method)
      if (trim(config_uplift_method) == 'fastisostasy') then

         ! This routine is called only on rank 0 (see li_core_run), which owns
         ! isosFastIsostasy and bslFastIsostasy. It must not contain collective MPI
         ! calls (e.g., no MPI_BCAST) or rank 0 would hang waiting on other ranks.
         call mpas_pool_get_config(liConfigs, 'config_fastisostasy_restart_file', config_fastisostasy_restart_file)
         call mpas_pool_get_config(liConfigs, 'config_fastisostasy_bsl_restart_file', &
                                   config_fastisostasy_bsl_restart_file)
         call mpas_pool_get_config(liConfigs, 'config_fastisostasy_coupling_interval', &
                                   config_fastisostasy_coupling_interval)

         ! Construct timestamped filenames using MPAS clock (same source as MALI restart stream)
         call construct_fastisostasy_restart_filename(domain, config_fastisostasy_restart_file, &
                                                      restart_filename, err_tmp)
         err = ior(err, err_tmp)
         call construct_fastisostasy_restart_filename(domain, config_fastisostasy_bsl_restart_file, &
                                                      bsl_restart_filename, err_tmp)
         err = ior(err, err_tmp)

         if (err == 0) then
            ! Must match how fastisostasy_solve computes the time associated with the
            ! most recent isos_update call, so that a subsequent restart resumes at
            ! exactly the time value stamped in this file (see fastisostasy_init and
            ! find_fastisostasy_restart_timestep).
            fastisostasyTime = real(fastisostasyTimeStep, wp) * real(config_fastisostasy_coupling_interval, wp)

            call mpas_log_write("Writing FastIsostasy restart files at FastIsostasy time level $i:", &
                                intArgs=(/fastisostasyTimeStep/))
            call mpas_log_write("  " // trim(restart_filename))
            call mpas_log_write("  " // trim(bsl_restart_filename))

            ! Under optimization, fesm-utils' ncio (nc_put_att) evaluates
            ! int(actual_range) for every variable type, not just NF90_INT.
            ! Volume scalars such as V_af exceed the 32-bit integer range, so
            ! this raises FE_INVALID and traps in DEBUG builds. The result is
            ! discarded for float variables, so disable halting for these writes.
            call ieee_get_halting_mode(ieee_all, haltingMode)
            call ieee_set_halting_mode(ieee_all, .false.)

            call isos_restart_write(isosFastIsostasy, trim(restart_filename), fastisostasyTime)
            call bsl_restart_write(bslFastIsostasy, trim(bsl_restart_filename), fastisostasyTime)

            call ieee_set_flag(ieee_all, .false.)
            call ieee_set_halting_mode(ieee_all, haltingMode)
         endif
      endif
#endif

   !--------------------------------------------------------------------
   end subroutine li_bedtopo_write_restart


!***********************************************************************
!
!  routine interpolate
!
!> \brief   Perform interpolation between MALI mesh and SLM grid
!> \author  Holly Han
!> \date    December 2021
!> \details
!>  This routine contains the sparse matrix multiplication
!>  algorithm to interpolate between MPAS and Gaussian Grid.
!>  Note: This routine is a copy of the inerpolation
!>  routine written by Kristin Barton in the code
!>  mpas_ocn_vel_self_attraction_loading.F in MPAS-Ocean,
!>  with addition of reformatting the interpolated data
!>
!-----------------------------------------------------------------------

   subroutine interpolate(colValues, rowValues, sValues, dataIn, dataOut)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      integer, dimension(:), intent(in) :: rowValues, colValues
      real (kind=RKIND), dimension(:), intent(in) :: sValues, dataIn

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: dataOut

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real (kind=RKIND) :: rhs
      integer :: n_S, n, nRow, nCol

      n_S = size(sValues)
      n = 1
      rhs = 0.0_RKIND

      do while (n .LE. n_S)
         nRow = rowValues(n)
         do while ( (n.LT.(n_S)) .AND. (rowValues(n).EQ.nRow) )
            nCol = colValues(n)
            rhs = rhs + dataIn(nCol) * sValues(n)
            n = n + 1
         enddo
         if ( (n.EQ.n_S) .AND. (rowValues(n).EQ.nRow) ) then
            nCol = colValues(n)
            rhs = rhs + dataIn(nCol) * sValues(n)
            n = n + 1
         endif
         dataOut(nRow) = rhs
         rhs = 0.0_RKIND
      enddo

   !--------------------------------------------------------------------
   end subroutine interpolate



!***********************************************************************
!
!  routine interpolate_init
!
!> \brief   Sets up interpolation between MALI and SLM native grids
!> \author  Holly Han
!> \date    December 2021
!> \details
!>  This routine reads in map (weight) files needed to interpolate
!>  values of ice thicknesss, bedTopography, sea-level change
!>  between native grid of MALI (unstructured) and SLM (Gaussian).
!>  It also gathers and scatters data from and to multiple processors.
!>  Note: A big portion of the routine is copied from routine
!>  'ocn_vel_self_attraction_loading_init' written by Kristin Barton
!>  in the code smpas_ocn_vel_self_attraction_loading.F in MPAS-Ocean
!
!-----------------------------------------------------------------------

   subroutine interpolate_init(domain, err, mapToFileIn, mapFromFileIn)

      use, intrinsic :: ieee_exceptions, only: ieee_all, ieee_get_halting_mode, &
                                               ieee_set_halting_mode, ieee_set_flag

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain    !< Input/output: Domain

      character (len=*), intent(in), optional :: mapToFileIn
      character (len=*), intent(in), optional :: mapFromFileIn

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: err_tmp
      character (len=StrKIND), pointer :: config_MALI_to_SLM_weights_file
      character (len=StrKIND), pointer :: config_SLM_to_MALI_weights_file

      ! NetCDF and weights file variables
      integer :: toNcId, toNsDimId, toRowId, toColId, toSId
      integer :: fromNcId, fromNsDimId, fromRowId, fromColId, fromSId
      integer:: nMpasDimId, nGridDimId, toNsLen, fromNsLen
      character (len = NF90_MAX_NAME) :: toNsName, fromNsName, nMpasName, nGridName
      integer, pointer :: n_s
      character(len=StrKIND) :: mpasToGridFile, gridToMpasFile

      ! MPI variables
      integer :: curProc

      logical :: haltingMode(size(ieee_all))

      err = 0
      err_tmp = 0

      if (present(mapToFileIn) .and. present(mapFromFileIn)) then
         mpasToGridFile = trim(mapToFileIn)
         gridToMpasFile = trim(mapFromFileIn)
      else
         call mpas_pool_get_config(liConfigs, 'config_MALI_to_SLM_weights_file', config_MALI_to_SLM_weights_file)
         call mpas_pool_get_config(liConfigs, 'config_SLM_to_MALI_weights_file', config_SLM_to_MALI_weights_file)
         mpasToGridFile = trim(config_MALI_to_SLM_weights_file)
         gridToMpasFile = trim(config_SLM_to_MALI_weights_file)
      endif

      ! Begin MPI portion
      call MPI_COMM_RANK( domain % dminfo % comm, curProc, err_tmp)
      err = ior(err, err_tmp)

      !initialize interpolation
      if (curProc.eq.0) then

         ! HDF5 (netCDF-4) raises FP exceptions internally, which trap in DEBUG builds
         call ieee_get_halting_mode(ieee_all, haltingMode)
         call ieee_set_halting_mode(ieee_all, .false.)

         ! Open netcdf weights files
         call check( nf90_open(path = mpasToGridFile, mode = nf90_nowrite, ncid = toNcId), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_open(path = gridToMpasFile, mode = nf90_nowrite, ncid = fromNcId), err_tmp)
         err = ior(err, err_tmp)

         ! Get dimension ID
         call check( nf90_inq_dimid(toNcId, "n_a", nMpasDimId), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inq_dimid(toNcId, "n_s", toNsDimId), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inq_dimid(fromNcId, "n_a", nGridDimId), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inq_dimid(fromNcId, "n_s", fromNsDimId), err_tmp)
         err = ior(err, err_tmp)

         ! Get Variable IDs
         call check( nf90_inq_varid(toNcId, "row", toRowId), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inq_varid(toNcId, "col", toColId), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inq_varid(toNcId, "S", toSId), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inq_varid(fromNcId, "row", fromRowId), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inq_varid(fromNcId, "col", fromColId), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inq_varid(fromNcId, "S", fromSId), err_tmp)
         err = ior(err, err_tmp)

         ! Get Dimension Length
         call check( nf90_inquire_dimension(toNcId, toNsDimId, toNsName, toNsLen), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inquire_dimension(fromNcId, fromNsDimId, fromNsName, fromNsLen), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inquire_dimension(toNcId, nMpasDimId, nMpasName, nMpas), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_inquire_dimension(fromNcId, nGridDimId, nGridName, nGrid), err_tmp)
         err = ior(err, err_tmp)

         ! Allocate matrices to read data into
         allocate ( toRowValues (toNsLen) )
         allocate ( toColValues (toNsLen) )
         allocate ( toSValues (toNsLen) )
         allocate ( fromRowValues (fromNsLen) )
         allocate ( fromColValues (fromNsLen) )
         allocate ( fromSValues (fromNsLen) )

         ! Retrieve data
         call check( nf90_get_var(toNcId, toColId, toColValues(:) ), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_get_var(toNcId, toRowId, toRowValues(:) ), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_get_var(toNcId, toSId, toSValues(:) ), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_get_var(fromNcId, fromColId, fromColValues(:) ), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_get_var(fromNcId, fromRowId, fromRowValues(:) ), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_get_var(fromNcId, fromSId, fromSValues(:) ), err_tmp)
         err = ior(err, err_tmp)

         call check( nf90_close(toNcId), err_tmp)
         err = ior(err, err_tmp)
         call check( nf90_close(fromNcId), err_tmp)
         err = ior(err, err_tmp)

         call ieee_set_flag(ieee_all, .false.)
         call ieee_set_halting_mode(ieee_all, haltingMode)

      endif

   !--------------------------------------------------------------------
   end subroutine interpolate_init



!***********************************************************************
!
!  routine check
!
!> \brief   Check status of netcdf operations
!> \author  Holly Han
!> \date    December 2021
!> \details
!>  This routine checks to status of the netcdf file
!
!-----------------------------------------------------------------------

   subroutine check(status, err)

      integer, intent ( in) :: status
      integer, intent(inout) :: err

      if(status /= nf90_noerr) then
         err = 1
      endif
   !--------------------------------------------------------------------
   end subroutine check


!***********************************************************************
!
!  routine check_SLM_coupling_interval
!
!> \brief   Perform various checks on the SLM coupling interval setting
!> \author  Matt Hoffman
!> \date    Feb 2024
!> \details
!>  This routine checks that the SLM coupling interval is an even multiple
!>  of the adaptive timestep force inverval and divides evenly into the
!>  restart interval.  It also checks that the coupling interval in the MALI
!>  matches the value in the SLM namelist.
!
!-----------------------------------------------------------------------

   subroutine check_SLM_coupling_interval(slm_dt1, meshPool, err)

      use mpas_timekeeping
      use mpas_stream_manager
      use mpas_derived_types, only : MPAS_STREAM_PROPERTY_RECORD_INTV

      integer, intent (in) :: slm_dt1
      type (mpas_pool_type), intent(in) :: meshPool  !< mesh information
      integer, intent(out) :: err

      ! local variables
      integer, pointer :: config_slm_coupling_interval
      logical, pointer :: config_adaptive_timestep
      character (len=StrKIND), pointer :: config_adaptive_timestep_force_interval, config_dt
      type (MPAS_TimeInterval_Type) :: coupling_interval, force_interval, dt_interval, zero_interval
      type (MPAS_Time_Type) :: start_time
      character (len=StrKIND), pointer :: simulationStartTime
      integer :: YYYY, MM, DD, H, M, S  ! time components
      integer (kind=I8KIND) :: n_intervals
      type (MPAS_TimeInterval_type) :: remainder
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_log_write("")
      call mpas_log_write("-- Checking consistency of config_slm_coupling_interval and other settings --")

      ! First, check consistency in coupling interval set up in MALI and SLM
      call mpas_pool_get_config(liConfigs, "config_slm_coupling_interval", config_slm_coupling_interval)
      if (config_slm_coupling_interval == slm_dt1) then
         call mpas_log_write("The coupling interval in MALI ($i yr) and SLM ($i yr) are consistent - check passes", &
                 intArgs=(/config_slm_coupling_interval, slm_dt1/))
      else
         call mpas_log_write("The coupling interval in MALI ($i yr) and SLM ($i yr) are inconsistent", MPAS_LOG_ERR, &
                 intArgs=(/config_slm_coupling_interval, slm_dt1/))
         err = ior(err, 1)
      endif

      ! define zero interval for comparing against below
      call mpas_set_timeInterval(zero_interval, dt = 0.0_RKIND)
      ! get start time as a reference time
      call mpas_pool_get_array(meshPool, 'simulationStartTime', simulationStartTime)
      call mpas_set_time(start_time, dateTimeString=simulationStartTime)
      ! define SLM coupling interval as a timeInterval type
      call mpas_set_timeInterval(coupling_interval, YY=config_slm_coupling_interval, MM=0, DD=0, H=0, M=0, S=0, ierr=err_tmp)
      err = ior(err, err_tmp)

      call mpas_pool_get_config(liConfigs, "config_adaptive_timestep", config_adaptive_timestep)
      if (config_adaptive_timestep) then
         ! for adaptive dt, check that config_adaptive_timestep_force_interval divides evenly into config_slm_coupling_interval
         call mpas_pool_get_config(liConfigs, "config_adaptive_timestep_force_interval", config_adaptive_timestep_force_interval)
         call mpas_set_timeInterval(force_interval, timeString=config_adaptive_timestep_force_interval, ierr=err_tmp)
         err = ior(err, err_tmp)
         call mpas_interval_division(start_time, coupling_interval, force_interval, n_intervals, remainder)
         if (remainder .EQ. zero_interval) then
            call mpas_log_write("config_adaptive_timestep_force_interval divides into config_slm_coupling_interval $i times " // &
                 "with no remainder - check passes", intArgs=(/int(n_intervals)/))
         else
            call mpas_log_write("config_adaptive_timestep_force_interval divides into config_slm_coupling_interval $i times " // &
                 "with nonzero remainder", MPAS_LOG_ERR, intArgs=(/int(n_intervals)/))
            err = ior(err, 1)
         endif
      else
        ! For fixed dt, check that dt divides evenly into config_slm_coupling_interval
         call mpas_pool_get_config(liConfigs, "config_dt", config_dt)
         call mpas_set_timeInterval(dt_interval, timeString=config_dt, ierr=err_tmp)
         err = ior(err, err_tmp)
         call mpas_interval_division(start_time, coupling_interval, dt_interval, n_intervals, remainder)
         if (remainder .EQ. zero_interval) then
            call mpas_log_write("config_dt divides into config_slm_coupling_interval $i times " // &
                 "with no remainder - check passes", intArgs=(/int(n_intervals)/))
         else
            call mpas_log_write("config_dt divides into config_slm_coupling_interval $i times " // &
                 "with nonzero remainder", MPAS_LOG_ERR, intArgs=(/int(n_intervals)/))
            err = ior(err, 1)
         endif
      endif

      ! No need to compare restart interval and coupling interval because restarts with SLM are supported for
      ! any restart interval now

      call mpas_log_write("")

   !--------------------------------------------------------------------
   end subroutine check_SLM_coupling_interval


!***********************************************************************
!
!  routine find_slm_restart_timestep
!
!> \brief   Perform various checks on the SLM coupling interval setting
!> \author  Matt Hoffman
!> \date    Feb 2024
!> \details
!>  This routine checks that the SLM coupling interval is an even multiple
!>  of the adaptive timestep force inverval and divides evenly into the
!>  restart interval.  It also checks that the coupling interval in the MALI
!>  matches the value in the SLM namelist.
!
!-----------------------------------------------------------------------

   subroutine find_slm_restart_timestep(meshPool, slmTimeStep, err)

      use mpas_timekeeping

      type (mpas_pool_type), intent(in) :: meshPool  !< mesh information
      integer, intent(out) :: slmTimeStep
      integer, intent(out) :: err

      ! local vars
      integer, pointer :: config_slm_coupling_interval
      character (len=StrKIND), pointer :: xtime, simulationStartTime
      type (MPAS_TimeInterval_Type) :: coupling_interval
      type (MPAS_Time_Type) :: start_time, curr_time
      integer (kind=I8KIND) :: n_intervals
      type (MPAS_TimeInterval_type) :: remainder
      character (len=StrKIND) :: remainder_string
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, "config_slm_coupling_interval", config_slm_coupling_interval)
      ! define SLM coupling interval as a timeInterval type
      call mpas_set_timeInterval(coupling_interval, YY=config_slm_coupling_interval, MM=0, DD=0, H=0, M=0, S=0, ierr=err_tmp)
      err = ior(err, err_tmp)

      call mpas_pool_get_array(meshPool, 'simulationStartTime', simulationStartTime)
      call mpas_pool_get_array(meshPool, 'xtime', xtime)
      call mpas_set_time(start_time, dateTimeString=simulationStartTime, ierr=err_tmp)
      err = ior(err, err_tmp)
      call mpas_set_time(curr_time, dateTimeString=xtime, ierr=err_tmp)
      err = ior(err, err_tmp)

      call mpas_interval_division(start_time, curr_time - start_time, coupling_interval, n_intervals, remainder)
      slmTimeStep = int(n_intervals)
      call mpas_get_timeInterval(remainder, start_time, timeString=remainder_string, ierr=err_tmp)
      err = ior(err, err_tmp)
      call mpas_log_write("SLM Restart: Using SLM time level $i because config_slm_coupling_interval divides into " // &
              "elapsed time that many times ", intArgs=(/int(slmTimeStep)/))
      call mpas_log_write("  That calculation implies it has been " // trim(remainder_string) // " since last SLM coupling. " // &
              "If that interval seems wrong, there may be an error in your configuration.")
      call mpas_log_write("")

   !--------------------------------------------------------------------
   end subroutine find_slm_restart_timestep


!***********************************************************************
!
!  routine construct_fastisostasy_restart_filename
!
!> \brief   Constructs timestamped FastIsostasy restart filename
!> \author  Matt Hoffman
!> \date    September 2026
!> \details
!>  Constructs a timestamped FastIsostasy restart filename by querying
!>  the MPAS TimeManager clock and appending the timestamp to the base
!>  filename. The timestamp is obtained from the same clock source used
!>  by MALI's restart stream, ensuring that FastIsostasy restart files
!>  match MALI restart file timestamps exactly.
!>  Format: <base_filename>.<YYYY-MM-DD_hh.mm.ss>.nc
!>  Example: fastisostasy_restart.2026-09-25_00.00.00.nc
!
!-----------------------------------------------------------------------

   subroutine construct_fastisostasy_restart_filename(domain, base_filename, output_filename, err)

      use mpas_timekeeping

      type (domain_type), intent(in) :: domain
      character(len=*), intent(in) :: base_filename
      character(len=StrKIND), intent(out) :: output_filename
      integer, intent(out) :: err

      ! local vars
      type (MPAS_Time_type) :: currTime
      character(len=StrKIND) :: timeStamp, timeStampFormatted, base_clean
      integer :: i, err_tmp

      err = 0

      ! Remove .nc extension if present in base filename
      base_clean = trim(base_filename)
      if (len_trim(base_clean) > 3) then
         if (base_clean(len_trim(base_clean)-2:len_trim(base_clean)) == '.nc') then
            base_clean = base_clean(1:len_trim(base_clean)-3)
         endif
      endif

      ! Query the MPAS clock for current time (same source as MALI restart stream)
      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, err_tmp)
      err = ior(err, err_tmp)
      if (err /= 0) return

      ! Extract timestamp string (format: "YYYY-MM-DD_hh:mm:ss")
      call mpas_get_time(currTime, dateTimeString=timeStamp, ierr=err_tmp)
      err = ior(err, err_tmp)
      if (err /= 0) return

      ! Format timestamp for filesystem: replace ':' with '.'
      ! Result: "YYYY-MM-DD_hh.mm.ss"
      timeStampFormatted = timeStamp
      do i = 1, len_trim(timeStampFormatted)
         if (timeStampFormatted(i:i) == ':') then
            timeStampFormatted(i:i) = '.'
         endif
      enddo

      ! Construct timestamped filename
      output_filename = trim(base_clean) // '.' // trim(timeStampFormatted) // '.nc'

   !--------------------------------------------------------------------
   end subroutine construct_fastisostasy_restart_filename


!***********************************************************************
!
!  routine find_fastisostasy_restart_timestep
!
!> \brief   Determines FastIsostasy's coupling step count on restart
!> \author  MALI developers
!> \date    June 2026
!> \details
!>  Computes how many FastIsostasy coupling intervals have elapsed since
!>  the simulation started, so that a restarted run resumes the
!>  fastisostasyTimeStep counter (and the FastIsostasy-internal time
!>  derived from it) at the correct point rather than restarting
!>  FastIsostasy's clock at year zero. Mirrors find_slm_restart_timestep.
!
!-----------------------------------------------------------------------

   subroutine find_fastisostasy_restart_timestep(meshPool, fastisostasyTimeStep, err)

      use mpas_timekeeping

      type (mpas_pool_type), intent(in) :: meshPool  !< mesh information
      integer, intent(out) :: fastisostasyTimeStep
      integer, intent(out) :: err

      ! local vars
      integer, pointer :: config_fastisostasy_coupling_interval
      character (len=StrKIND), pointer :: xtime, simulationStartTime
      type (MPAS_TimeInterval_Type) :: coupling_interval
      type (MPAS_Time_Type) :: start_time, curr_time
      integer (kind=I8KIND) :: n_intervals
      type (MPAS_TimeInterval_type) :: remainder
      character (len=StrKIND) :: remainder_string
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, "config_fastisostasy_coupling_interval", config_fastisostasy_coupling_interval)
      ! define FastIsostasy coupling interval as a timeInterval type
      call mpas_set_timeInterval(coupling_interval, YY=config_fastisostasy_coupling_interval, MM=0, DD=0, H=0, M=0, S=0, &
           ierr=err_tmp)
      err = ior(err, err_tmp)

      call mpas_pool_get_array(meshPool, 'simulationStartTime', simulationStartTime)
      call mpas_pool_get_array(meshPool, 'xtime', xtime)
      call mpas_set_time(start_time, dateTimeString=simulationStartTime, ierr=err_tmp)
      err = ior(err, err_tmp)
      call mpas_set_time(curr_time, dateTimeString=xtime, ierr=err_tmp)
      err = ior(err, err_tmp)

      call mpas_interval_division(start_time, curr_time - start_time, coupling_interval, n_intervals, remainder)
      fastisostasyTimeStep = int(n_intervals)
      call mpas_get_timeInterval(remainder, start_time, timeString=remainder_string, ierr=err_tmp)
      err = ior(err, err_tmp)
      call mpas_log_write("FastIsostasy Restart: Using FastIsostasy time level $i because " // &
              "config_fastisostasy_coupling_interval divides into elapsed time that many times ", &
              intArgs=(/int(fastisostasyTimeStep)/))
      call mpas_log_write("  That calculation implies it has been " // trim(remainder_string) // &
              " since last FastIsostasy coupling. If that interval seems wrong, there may be an error " // &
              "in your configuration.")
      call mpas_log_write("")

   !--------------------------------------------------------------------
   end subroutine find_fastisostasy_restart_timestep

!***********************************************************************

end module li_bedtopo

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
