! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration
!
!> \brief MPAS ocean time integration driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for calling
!>  the time integration scheme
!
!-----------------------------------------------------------------------

module ocn_time_integration

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer
   use mpas_io_units

   use ocn_constants
   use ocn_time_integration_rk4
   use ocn_time_integration_split

   implicit none
   private
   save

   public :: ocn_timestep, &
             ocn_timestep_init

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

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

    logical :: rk4On, splitOn

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_timestep
!
!> \brief MPAS ocean time integration driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine handles a single timestep for the ocean. It determines
!>  the time integrator that will be used for the run, and calls the
!>  appropriate one.
!
!-----------------------------------------------------------------------

   subroutine ocn_timestep(domain, dt, timeStamp)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:))
   !                 plus mesh meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp

      type (dm_info) :: dminfo
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: diagnosticsPool, statePool, meshPool

      character (len=StrKIND), pointer :: xtime
      real (kind=RKIND), pointer :: daysSinceStartOfSim
      character (len=StrKIND), pointer :: simulationStartTime
      type (MPAS_Time_type) :: xtime_timeType, simulationStartTime_timeType


      if (rk4On) then
         call ocn_time_integrator_rk4(domain, dt)
      elseif (splitOn) then
         call ocn_time_integrator_split(domain, dt)
     endif

     block => domain % blocklist
     do while (associated(block))
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

        call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)

        xtime = timeStamp

        ! compute time since start of simulation, in days
        call mpas_pool_get_array(diagnosticsPool, 'simulationStartTime', simulationStartTime)
        call mpas_pool_get_array(diagnosticsPool, 'daysSinceStartOfSim',daysSinceStartOfSim)
        call mpas_set_time(xtime_timeType, dateTimeString=xtime)
        call mpas_set_time(simulationStartTime_timeType, dateTimeString=simulationStartTime)
        call mpas_get_timeInterval(xtime_timeType - simulationStartTime_timeType,dt=daysSinceStartOfSim)
        daysSinceStartOfSim = daysSinceStartOfSim*days_per_second

        ! test ocean state for abnormal values and NaNs in the velocity field
        call test_ocean_state(dminfo, meshPool, diagnosticsPool, statePool)

        block => block % next
     end do

   contains

   subroutine test_ocean_state(dminfo, meshPool, diagnosticsPool, statePool)

      type (dm_info) :: dminfo
      type (mpas_pool_type), pointer :: statePool, meshPool, diagnosticsPool, tracersPool
      character (len=StrKIND), pointer :: xtime
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      integer, dimension(:), pointer   :: maxLevelCell
      integer, dimension(:), pointer :: indexToCellID
      integer, pointer :: indexTemperature
      integer, pointer :: indexSalinity
      integer, pointer :: nCellsSolve
      real (kind=RKIND) :: nanCheck, workValue, workLat, workLon
      integer :: workGlobalID(2), errorUnit, mpiRank, iCell, k
      logical :: errorFlag
      character(len=StrKIND) :: charMPIRank, charFilename

      !get all pointers that might be needed
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
      call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 2)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 2)

      !assume that no abnormal values exist
      errorFlag = .false.

      !now step through tests to see if abnormal values do exist
      !if so, reset errorFlag to be true

      ! test for abnormal values
      do iCell=1,nCellsSolve
         do k=1,maxLevelCell(iCell)
            if(kineticEnergyCell(k,iCell).gt.4.0_RKIND) errorFlag=.true.
            if(activeTracers(indexTemperature,k,iCell).lt.-1.9_RKIND) errorFlag=.true.
            if(activeTracers(indexTemperature,k,iCell).gt.33.0_RKIND) errorFlag=.true.
            if(activeTracers(indexSalinity,k,iCell).lt.0.0_RKIND) errorFlag=.true.
            if(layerThickness(k,iCell).lt.1.0e-2_RKIND) errorFlag=.true.
         enddo
      enddo

      !if an errorFlag exists, then
      !  1) open a file
      !  2) step through tests again and write the file
      !  3) close file

      if(errorFlag) then

        !find an open unit
        call mpas_new_unit(errorUnit)
        mpiRank = dminfo % my_proc_id
        charFilename = 'mpas_ocean_state_test_'
        if(                       mpiRank.le.     9) write(charMPIRank,'(I1)') mpiRank
        if(mpiRank.gt.    9 .and. mpiRank.le.    99) write(charMPIRank,'(I2)') mpiRank
        if(mpiRank.gt.   99 .and. mpiRank.le.   999) write(charMPIRank,'(I3)') mpiRank
        if(mpiRank.gt.  999 .and. mpiRank.le.  9999) write(charMPIRank,'(I4)') mpiRank
        if(mpiRank.gt. 9999 .and. mpiRank.le. 99999) write(charMPIRank,'(I5)') mpiRank
        if(mpiRank.gt.99999 .and. mpiRank.le.999999) write(charMPIRank,'(I6)') mpiRank
        charFilename = trim(charFilename) // trim(charMPIRank)
        open(unit=errorUnit, file=charFilename, form='formatted', status='unknown', position='append')

        !write time
        write(errorUnit,'(a80)') trim(xtime)

        !test to see if cell kinetic energy is greater than 4.0 m2/s2
        do iCell=1,nCellsSolve
           do k=1,maxLevelCell(iCell)
              if(kineticEnergyCell(k,iCell).gt.4.0_RKIND) then
                 workValue = kineticEnergyCell(k,iCell)
                 workGlobalID(1) = k
                 workGlobalID(2) = iCell
                 workLat = latCell(iCell)
                 workLon = lonCell(iCell)
                 write(errorUnit, 10) 'KE= ', workValue, 'cell= ', workGlobalID(2), &
                    'k= ',workGlobalID(1), 'lat = ', workLat, 'lon= ', workLon
                 10 format(a4,e10.3, 3x,a6,i8, 3x,a3,i4, 3x,a6,f6.2, 3x,a6,f6.2)
              endif
           enddo
        enddo

        !test to see if cell temperature is less than -1.9C
        do iCell=1,nCellsSolve
           do k=1,maxLevelCell(iCell)
              if(activeTracers(indexTemperature,k,iCell).lt.-1.9_RKIND) then
                 workValue = activeTracers(indexTemperature,k,iCell)
                 workGlobalID(1) = k
                 workGlobalID(2) = iCell
                 workLat = latCell(iCell)
                 workLon = lonCell(iCell)
                 write(errorUnit, 10) 'T= ', workValue, 'cell= ', workGlobalID(2), &
                    'k= ',workGlobalID(1), 'lat = ', workLat, 'lon= ', workLon
              endif
           enddo
        enddo

        !test to see if cell temperature is greater than 33.0
        do iCell=1,nCellsSolve
           do k=1,maxLevelCell(iCell)
              if(activeTracers(indexTemperature,k,iCell).gt.33.0_RKIND) then
                 workValue = activeTracers(indexTemperature,k,iCell)
                 workGlobalID(1) = k
                 workGlobalID(2) = iCell
                 workLat = latCell(iCell)
                 workLon = lonCell(iCell)
                 write(errorUnit, 10) 'T= ', workValue, 'cell= ', workGlobalID(2), &
                    'k= ',workGlobalID(1), 'lat = ', workLat, 'lon= ', workLon
              endif
           enddo
        enddo

        !test to see if cell salinity is less than 0
        do iCell=1,nCellsSolve
           do k=1,maxLevelCell(iCell)
              if(activeTracers(indexSalinity,k,iCell).lt.0.0_RKIND) then
                 workValue = activeTracers(indexSalinity,k,iCell)
                 workGlobalID(1) = k
                 workGlobalID(2) = iCell
                 workLat = latCell(iCell)
                 workLon = lonCell(iCell)
                 write(errorUnit, 10) 'S= ', workValue, 'cell= ', workGlobalID(2), &
                    'k= ',workGlobalID(1), 'lat = ', workLat, 'lon= ', workLon
              endif
           enddo
        enddo

        !test to see if cell thickness is less than 1.0e-2
        do iCell=1,nCellsSolve
           do k=1,maxLevelCell(iCell)
              if(layerThickness(k,iCell).lt.1.0e-2_RKIND) then
                 workValue = layerThickness(k,iCell)
                 workGlobalID(1) = k
                 workGlobalID(2) = iCell
                 workLat = latCell(iCell)
                 workLon = lonCell(iCell)
                 write(errorUnit, 10) 'S= ', workValue, 'cell= ', workGlobalID(2), &
                    'k= ',workGlobalID(1), 'lat = ', workLat, 'lon= ', workLon
              endif
           enddo
        enddo

        write(errorUnit,*) ''

        !close unit
        close(errorUnit)
        call mpas_release_unit(errorUnit)

      endif ! if(errorFlag)

      !now test for NaNs in the ocean velocity field and abort if true
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 2)
      call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
      nanCheck = sum(normalVelocity)
      if (nanCheck /= nanCheck) then
         write(stderrUnit,*) 'Ocean Abort: NaN detected'
         call mpas_dmpar_global_abort('MPAS-ocean: Ocean Abort: NaN detected')
      endif

   end subroutine test_ocean_state

   end subroutine ocn_timestep!}}}

   subroutine ocn_timestep_init(err)!{{{

      integer, intent(out) :: err

      character (len=StrKIND), pointer :: config_time_integrator

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_time_integrator', config_time_integrator)

      rk4On = .false.
      splitOn = .false.

      if (trim(config_time_integrator) == 'RK4') then
          rk4On = .true.
      elseif (trim(config_time_integrator) == 'split_explicit' &
          .or.trim(config_time_integrator) == 'unsplit_explicit') then
          splitOn = .true.
      else
          err = 1
          write (stderrUnit,*) 'Incorrect choice for config_time_integrator:', trim(config_time_integrator)
          write (stderrUnit,*) '   choices are: RK4, split_explicit, unsplit_explicit'
      endif


   end subroutine ocn_timestep_init!}}}

end module ocn_time_integration

! vim: foldmethod=marker
