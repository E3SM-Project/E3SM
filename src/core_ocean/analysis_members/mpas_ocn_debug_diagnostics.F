! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_debug_diagnostics
!
!> \brief MPAS ocean analysis mode member: debug_diagnostics
!> \author Mark Petersen
!> \date   March 2016
!> \details
!>  MPAS ocean analysis mode member: debug_diagnostics
!>  Compute diagnostics used for debugging.
!>
!-----------------------------------------------------------------------

module ocn_debug_diagnostics

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_io_units

   use ocn_constants
   use ocn_diagnostics_routines

   implicit none
   private
   save

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

   public :: ocn_init_debug_diagnostics, &
             ocn_compute_debug_diagnostics, &
             ocn_restart_debug_diagnostics, &
             ocn_finalize_debug_diagnostics

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_debug_diagnostics
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_debug_diagnostics(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_init_debug_diagnostics!}}}

!***********************************************************************
!
!  routine ocn_compute_debug_diagnostics
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_debug_diagnostics(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: debugDiagnosticsAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      !type (mpas_pool_type), pointer :: debugDiagnosticsAM

      ! Here are some example variables which may be needed for your analysis member
      integer :: iEdge, c1, c2, k
      integer, pointer :: nEdges
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) :: dzVert1, dzVert2, dzEdgeK, dzEdgeKp1, rx1, localMaxRx1
      real (kind=RKIND), pointer :: globalRx1Max
      real (kind=RKIND), dimension(:), pointer :: rx1MaxCell
      real (kind=RKIND), dimension(:,:), pointer :: zMid

      logical, pointer :: config_AM_debugDiagnostics_check_state

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_AM_debugDiagnostics_check_state', &
                                config_AM_debugDiagnostics_check_state)

      localMaxRx1 = 0.0_RKIND

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         !call mpas_pool_get_subpool(block % structs, 'debugDiagnosticsAM', debugDiagnosticsAMPool)

         if ( config_AM_debugDiagnostics_check_state ) then
            call ocn_test_ocean_state(dminfo, meshPool, diagnosticsPool, statePool)
         end if

         ! Here are some example variables which may be needed for your analysis member
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

         call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

         !-----------------------------------------------------------------
         !
         ! Compute Haney number, rx1
         !
         !-----------------------------------------------------------------

         call mpas_pool_get_array(diagnosticsPool, 'rx1MaxCell', rx1MaxCell)
         call mpas_pool_get_array(diagnosticsPool, 'globalRx1Max', globalRx1Max)

        ! These could be included for edge or cell fields with depth:
        ! call mpas_pool_get_array(diagnosticsPool, 'rx1Edge', rx1Edge)
        ! call mpas_pool_get_array(diagnosticsPool, 'rx1Cell', rx1Cell)
        ! call mpas_pool_get_array(diagnosticsPool, 'rx1MaxEdge', rx1MaxEdge)
        ! rx1Edge(:,:) = 0.0_RKIND
        ! rx1Cell(:,:) = 0.0_RKIND
        ! rx1MaxEdge(:) = 0.0_RKIND

        rx1MaxCell(:) = 0.0_RKIND
        do iEdge = 1,nEdges
          c1 = cellsOnEdge(1,iEdge)
          c2 = cellsOnEdge(2,iEdge)
          do k = 1,maxLevelEdgeTop(iEdge)-1
            dzVert1 = zMid(k,c1)-zMid(k+1,c1)
            dzVert2 = zMid(k,c2)-zMid(k+1,c2)
            dzEdgeK = zMid(k,c2)-zMid(k,c1)
            dzEdgeKp1 = zMid(k+1,c2)-zMid(k+1,c1)

            rx1 = abs(dzEdgeK+dzEdgeKp1)/(dzVert1+dzVert2)

            rx1MaxCell(c1) = max(rx1MaxCell(c1),rx1)
            rx1MaxCell(c2) = max(rx1MaxCell(c2),rx1)

            ! These could be included for edge or cell fields with depth:
            ! rx1Edge(k,iEdge) = rx1
            ! rx1Cell(k,c1) = max(rx1Cell(k,c1),rx1)
            ! rx1Cell(k,c2) = max(rx1Cell(k,c2),rx1)
            ! rx1MaxEdge(iEdge) = max(rx1MaxEdge(iEdge),rx1)
          end do
        end do

        localMaxRx1 = max(localMaxRx1,maxval(rx1MaxCell))

        block => block % next
      end do
      call mpas_dmpar_max_real(dminfo, localMaxRx1, globalRx1Max)

   end subroutine ocn_compute_debug_diagnostics!}}}

!***********************************************************************
!
!  routine ocn_restart_debug_diagnostics
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_debug_diagnostics(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_debug_diagnostics!}}}

!***********************************************************************
!
!  routine ocn_finalize_debug_diagnostics
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_debug_diagnostics(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_finalize_debug_diagnostics!}}}

   subroutine ocn_test_ocean_state(dminfo, meshPool, diagnosticsPool, statePool)!{{{

      type (dm_info) :: dminfo
      type (mpas_pool_type), pointer :: statePool, meshPool, diagnosticsPool, tracersPool
      character (len=StrKIND), pointer :: xtime
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      integer, dimension(:), pointer :: maxLevelCell
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
                 workGlobalID(2) = indexToCellID(iCell)
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
                 workGlobalID(2) = indexToCellID(iCell)
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
                 workGlobalID(2) = indexToCellID(iCell)
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
                 workGlobalID(2) = indexToCellID(iCell)
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
                 workGlobalID(2) = indexToCellID(iCell)
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

   end subroutine ocn_test_ocean_state!}}}

end module ocn_debug_diagnostics

! vim: foldmethod=marker
