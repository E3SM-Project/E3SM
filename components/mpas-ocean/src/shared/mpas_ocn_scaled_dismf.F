! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_scaled_dismf
!> \brief MPAS ocean scale data ice-shelf melt fluxes
!> \author Xylar Asay-Davis
!> \date   July 2024
!> \details
!>  This module contains routines for scaling data ice-shelf melt fluxes by
!>  the running mean of the removed ice runoff
!>  Design document located at:
!>    https://acme-climate.atlassian.net/wiki/spaces/PSC/pages/4210098268/Design+Document+Data+iceberg+and+ice-shelf+melt+flux+patterns+for+E3SM+runs
!
!-----------------------------------------------------------------------

module ocn_scaled_dismf

    use mpas_kind_types
    use mpas_derived_types
    use mpas_global_sum_mod
    use mpas_timekeeping
    use mpas_timer

    use ocn_config
    use ocn_mesh

    use shr_const_mod

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

    public :: ocn_init_scaled_dismf, &
              ocn_update_scaled_dismf

    !--------------------------------------------------------------------
    !
    ! Private module variables
    !
    !--------------------------------------------------------------------

    logical :: scaledDISMFOn
    character (len=*), parameter :: alarmID = 'scaledDISMFUpdateAlarm'

!***********************************************************************

contains


!***********************************************************************
!
!  routine ocn_init_scaled_dismf
!
!> \brief   Initialize scaling of data ice-shelf melt fluxes
!> \author Xylar Asay-Davis
!> \date   July 2024
!> \details
!>  Set alarms needed to compute daily and running means of removed ice runoff
!
!-----------------------------------------------------------------------

    subroutine ocn_init_scaled_dismf(domain)!{{{

        type (domain_type), intent(inout) :: domain

        ! Alarm variables
        type (MPAS_Time_Type) :: alarmTime
        type (MPAS_TimeInterval_type) :: alarmTimeStep

        integer :: err_tmp

        if (config_scale_dismf_by_removed_ice_runoff) then
            scaledDISMFOn = .true.
        else
            scaledDISMFOn = .false.
            return
        endif

        if (.not. config_remove_ais_ice_runoff) then
            call mpas_log_write('config_scale_dismf_by_removed_ice_runoff = .true. requires config_remove_ais_ice_runoff = .true.', &
                                MPAS_LOG_CRIT)
        endif

        if (trim(config_land_ice_flux_mode) /= 'data') then
            call mpas_log_write('config_scale_dismf_by_removed_ice_runoff = .true. requires config_land_ice_flux_mode = "data"', &
                                MPAS_LOG_CRIT)
        endif

        ! Setup Alarm for updating of scaled DISMF
        alarmTime = mpas_get_clock_time(domain % clock, &
                                        MPAS_START_TIME, &
                                        ierr=err_tmp)
        call mpas_set_timeInterval(alarmTimeStep, &
                                   timeString='0000-00-01_00:00:00', &
                                   ierr=err_tmp)
        call mpas_add_clock_alarm(domain % clock, alarmID, alarmTime + alarmTimeStep, &
                                  alarmTimeInterval=alarmTimeStep, ierr=err_tmp)

    end subroutine ocn_init_scaled_dismf!}}}


!***********************************************************************
!
!  routine ocn_update_scaled_dismf
!
!> \brief  Update scaled data ice-shelf melt fluxes
!> \author Xylar Asay-Davis
!> \date   August 2024
!> \details
!>  Accumulate daily mean of removed runoff.  If we are at the end of a day,
!>  update the running mean of removed ice runoff, and use it to scale
!>  the ice shelf melt flux based on the pattern from a data file
!
!-----------------------------------------------------------------------

    subroutine ocn_update_scaled_dismf(domain)!{{{

        type (domain_type), intent(inout) :: domain

        integer :: err_tmp

        if (.not. scaledDISMFOn) then
            return
        end if

        ! since the current clock time is for the end of the accumulation
        ! intereval, first, accumulate the daily mean
        call accumulate_mean_removed_ice_runoff(domain)

       ! then, compute update the history and running mean if we are at the
       ! end of a day
       if(mpas_is_alarm_ringing(domain % clock, alarmID, ierr=err_tmp)) then
#ifdef MPAS_DEBUG
            call mpas_log_write('       Computing Scaled DISMF')
#endif
            call update_scaled_dismf(domain)
            call mpas_reset_clock_alarm(domain % clock, alarmID, ierr=err_tmp)
        endif

    end subroutine ocn_update_scaled_dismf!}}}


!***********************************************************************
!
!  routine accumulate_mean_removed_ice_runoff
!
!> \brief  Accumulate mean removed ice runoff
!> \author Xylar Asay-Davis
!> \date   August 2024
!> \details
!>  Accumulate current removed ice runoff into the daily mean value
!
!-----------------------------------------------------------------------

    subroutine accumulate_mean_removed_ice_runoff(domain)!{{{

        type (domain_type), intent(inout) :: domain

        type (block_type), pointer :: block_ptr
        type (mpas_pool_type), pointer :: forcingPool

        real (kind=RKIND), dimension(:), pointer :: removedIceRunoffFlux
        real (kind=RKIND), pointer :: avgRemovedIceRunoff
        integer, pointer :: nCellsSolve, nAccumulated

        real (kind=RKIND) :: totalFlux
        real (kind=RKIND), dimension(:), allocatable :: localArrayForReproSum

        integer :: iCell

        block_ptr => domain % blocklist

        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        ! independent of blocks
        call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoff', avgRemovedIceRunoff)
        call mpas_pool_get_dimension(forcingPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_array(forcingPool, 'nAccumulatedRemovedIceRunoff', nAccumulated)

        do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

            call mpas_pool_get_array(forcingPool, 'removedIceRunoffFlux', removedIceRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoff', avgRemovedIceRunoff)

            call mpas_pool_get_dimension(forcingPool, 'nCellsSolve', nCellsSolve)

            call mpas_pool_get_array(forcingPool, 'nAccumulatedRemovedIceRunoff', nAccumulated)


            allocate(localArrayForReproSum(nCellsSolve))
            localArrayForReproSum(:) = 0.0_RKIND

            !$omp parallel
            !$omp do schedule(runtime)
            do iCell=1,nCellsSolve
                localArrayForReproSum(iCell) = removedIceRunoffFlux(iCell) * areaCell(iCell)
            enddo
            !$omp end do
            !$omp end parallel

            block_ptr => block_ptr % next
        end do ! block_ptr

        totalFlux = mpas_global_sum(localArrayForReproSum, domain % dminfo % comm)

        deallocate (localArrayForReproSum)

        avgRemovedIceRunoff = (avgRemovedIceRunoff * nAccumulated + totalFlux) &
                            / ( nAccumulated + 1 )

        nAccumulated = nAccumulated + 1

    end subroutine accumulate_mean_removed_ice_runoff!}}}



!***********************************************************************
!
!  routine update_scaled_dismf
!
!> \brief  Update scaled data ice-shelf melt fluxes
!> \author Xylar Asay-Davis
!> \date   August 2024
!> \details
!>  Update the running mean of removed ice runoff, and use it to scale
!>  the ice shelf melt flux based on the pattern from a data file
!
!-----------------------------------------------------------------------

    subroutine update_scaled_dismf(domain)!{{{

        type (domain_type), intent(inout) :: domain

        type (block_type), pointer :: block_ptr
        type (mpas_pool_type), pointer :: forcingPool

        real (kind=RKIND), pointer :: avgRemovedIceRunoff, runningMeanRemovedIceRunoff
        real (kind=RKIND), dimension(:), pointer :: totalRemovedIceRunoffHistory
        integer, pointer :: nValidHistory, nAccumulated

        real (kind=RKIND), dimension(:), allocatable :: tmpHistory

        real (kind=RKIND) :: previousTotal, timeInterval
        integer :: nHistory

        block_ptr => domain % blocklist
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        call mpas_pool_get_array(forcingPool, 'nAccumulatedRemovedIceRunoff', nAccumulated)
        call mpas_pool_get_array(forcingPool, 'nValidTotalRemovedIceRunoffHistory', nValidHistory)
        call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoff', avgRemovedIceRunoff)
        call mpas_pool_get_array(forcingPool, 'totalRemovedIceRunoffHistory', totalRemovedIceRunoffHistory)
        call mpas_pool_get_array(forcingPool, 'runningMeanRemovedIceRunoff', runningMeanRemovedIceRunoff)

        nHistory = config_ais_ice_runoff_history_days

        previousTotal = totalRemovedIceRunoffHistory(nValidHistory)

        if (nValidHistory == 0) then
            ! we keep zero as the first entry in totalRemovedIceRunoffHistory
            nValidHistory = 1
        end if

        if (nValidHistory == nHistory) then
            ! we need to shift the history, since it's full
            allocate(tmpHistory(nHistory))
            tmpHistory(:) = totalRemovedIceRunoffHistory(:)
            totalRemovedIceRunoffHistory(1:nHistory - 1) = tmpHistory(2:nHistory)
        else
            ! the history isn't full yet, so we just add to the end
            nValidHistory = nValidHistory + 1
        end if

        ! add the new total, the previous total plus the new daily average
        totalRemovedIceRunoffHistory(nValidHistory) = previousTotal + SHR_CONST_CDAY * avgRemovedIceRunoff

        timeInterval = SHR_CONST_CDAY * (nValidHistory - 1)
        ! the running mean is the difference between the newest and oldest
        ! totals divided by the time between them
        runningMeanRemovedIceRunoff = &
            (totalRemovedIceRunoffHistory(nValidHistory) - totalRemovedIceRunoffHistory(1)) &
            / timeInterval

        ! reset daily averaging of the removed runoff
        nAccumulated = 0
        avgRemovedIceRunoff = 0.0_RKIND

    end subroutine update_scaled_dismf!}}}

end module ocn_scaled_dismf