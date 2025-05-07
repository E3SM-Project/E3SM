! Copyright (c) 2013,  Pacific Northwest National Laboratory, (PNNL)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_scaled_sfwf
!> \brief MPAS ocean scale data freshwater fluxes
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
!>  Computes a precipitation factor to multiply the fresh water flux
!>  due to precipitation uniformly to ensure a balance of fresh water
!>  at the ocean surface.
!>  Designed following precip_adjustment subroutine located at:
!>    https://github.com/ESCOMP/POP2-CESM/blob/master/source/forcing_sfwf.F90
!
!-----------------------------------------------------------------------

module ocn_scaled_sfwf

    use mpas_kind_types
    use mpas_derived_types
    use mpas_global_sum_mod
    use mpas_timekeeping
    use mpas_timer

    use ocn_config
    use ocn_mesh

    use shr_const_mod

    use ocn_constants, only: rho_fw, ocn_ref_salinity

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

    public :: ocn_init_scaled_sfwf, &
              ocn_update_scaling_factor

    !--------------------------------------------------------------------
    !
    ! Private module variables
    !
    !--------------------------------------------------------------------

    logical :: scaledSFWFOn
    character (len=*), parameter :: alarmID = 'scaledSFWFUpdateAlarm'

!***********************************************************************

contains


!***********************************************************************
!
!  routine ocn_init_scaled_sfwf
!
!> \brief   Initialize scaling of data freshwater fluxes
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
!>  Set alarms needed to compute yearly and running means of water fluxes 
!
!-----------------------------------------------------------------------

    subroutine ocn_init_scaled_sfwf(domain)!{{{

        type (domain_type), intent(inout) :: domain

        ! Alarm variables
        type (mpas_Time_Type) :: alarmTime
        type (mpas_TimeInterval_type) :: alarmTimeStep

        integer :: err_tmp

        if (config_sfwf_balance .and. config_sfwf_balance_type > 0) then
            scaledSFWFOn = .true.
        else
            scaledSFWFOn = .false.
            return
        endif

        if (config_sfwf_history_years < 1) then
            call mpas_log_write('scaledSFWFOn = .true. requires config_sfwf_history_years >= 1', &
                                mpas_LOG_CRIT)
        endif

        ! Setup Alarm for updating of scaled SFWF
        alarmTime = mpas_get_clock_time(domain % clock, &
                                        mpas_START_TIME, &
                                        ierr=err_tmp)
        call mpas_set_timeInterval(alarmTimeStep, &
                                   timeString='0001-00-00_00:00:00', &
                                   ierr=err_tmp)
        call mpas_add_clock_alarm(domain % clock, alarmID, alarmTime + alarmTimeStep, &
                                  alarmTimeInterval=alarmTimeStep, ierr=err_tmp)

    end subroutine ocn_init_scaled_sfwf!}}}


!***********************************************************************
!
!  routine ocn_update_scaling_factor
!
!> \brief  Update scaling factors for data freshwater fluxes
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
!>  Accumulate yearly mean of surface freshwater fluxes. 
!>  If we are at the end of a year, used it to derive a factor to scale 
!>  down/up precipitation fluxes and maintain freshwater balance.  
!
!-----------------------------------------------------------------------

    subroutine ocn_update_scaling_factor(domain)!{{{

        type (domain_type), intent(inout) :: domain

        integer :: err_tmp

        if (.not. scaledSFWFOn) then
            return
        end if

        ! since the current clock time is for the end of the accumulation
        ! intereval, first, accumulate the yearly mean
        call accumulate_mean_fw_fluxes(domain)

       ! then, compute update the history and running mean if we are at the
       ! end of a year
       if(mpas_is_alarm_ringing(domain % clock, alarmID, ierr=err_tmp)) then
#ifdef mpas_DEBUG
            call mpas_log_write('       Computing Scaling Factor For SFWF')
#endif
            call update_scaling_factor(domain)
            call mpas_reset_clock_alarm(domain % clock, alarmID, ierr=err_tmp)
        endif

    end subroutine ocn_update_scaling_factor!}}}


!***********************************************************************
!
!  routine accumulate_mean_fw_fluxes
!
!> \brief  Accumulate surface water fluxes 
!> \author Shixuan Zhang
!> \date   May 2025
!> \details
!>  Accumulate current surface freshwater and precipitation fluxes 
!>  into the yearly mean value
!
!-----------------------------------------------------------------------

    subroutine accumulate_mean_fw_fluxes(domain)!{{{

        type (domain_type), intent(inout) :: domain

        type (block_type), pointer :: block_ptr

        type (mpas_pool_type), pointer ::  & 
                forcingPool, & 
                statePool, &
                tracersPool, &
                meshPool

        type (MPAS_timeInterval_type) :: &
                timeStepESMF

        real(kind=RKIND), dimension(:), pointer :: &
                areaCell, &
                rainFlux, &
                snowFlux, &
                sshOld, &
                sshNew

      
        real(kind=RKIND), dimension(:,:), pointer :: &
                layerThickness

        real(kind=RKIND), dimension(:,:,:), pointer :: &
                activeTracersNew, &
                activeTracersOld
  
        integer, dimension(:), pointer :: & 
                minLevelCell, & 
                maxLevelCell

        real (kind=RKIND), pointer :: &
                avgprecFlux,  &
                avgfwFlux, &
                avgsalFlux

        integer, pointer :: & 
                nCellsSolve, & 
                nVertLevels, &
                indexSalinity, &
                nAccumulated

        real (kind=RKIND) :: & 
                precFlux, & 
                fwFlux, &
                salFlux

        real(kind=RKIND) :: dt

        real (kind=RKIND), dimension(:), allocatable :: & 
             localArrayForPrecSum, &
             localArrayForSSHSum, &
             localArrayForSalinitySum

        integer :: &
             iCell, &
             k

        call mpas_get_timeInterval(timeStepESMF, dt=dt)

        block_ptr => domain % blocklist

        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        ! independent of blocks
        call mpas_pool_get_array(forcingPool, 'avgprecFlux', avgprecFlux)
        call mpas_pool_get_array(forcingPool, 'avgfwFlux', avgfwFlux)
        call mpas_pool_get_array(forcingPool, 'avgsalFlux', avgsalFlux)

        call mpas_pool_get_dimension(forcingPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_array(forcingPool, 'nAccumulatedfwFlux', nAccumulated)

        do while (associated(block_ptr))

           call mpas_pool_get_dimension(block_ptr % dimensions, "nCellsSolve", nCellsSolve)
           call mpas_pool_get_subpool(block_ptr % structs, "mesh", meshPool)
           call mpas_pool_get_subpool(block_ptr % structs, "forcing", forcingPool)
           call mpas_pool_get_subpool(block_ptr % structs, "state", statePool)

           call mpas_pool_get_array(meshPool, "areaCell", areaCell)
           call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
           call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

           call mpas_pool_get_array(statePool, "layerThickness", layerThickness, 1)
           call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
           call mpas_pool_get_array(statePool, 'ssh', sshOld, 1)

           call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
           call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
           call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersNew, 2) 
           call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersOld, 1) 

           call mpas_pool_get_array(forcingPool, 'rainFlux',               rainFlux)
           call mpas_pool_get_array(forcingPool, 'snowFlux',               snowFlux)

           call mpas_pool_get_array(forcingPool, 'avgprecFlux', avgprecFlux)
           call mpas_pool_get_array(forcingPool, 'avgfwFlux', avgfwFlux)
           call mpas_pool_get_array(forcingPool, 'avgsalFlux', avgsalFlux)

           call mpas_pool_get_array(forcingPool, 'nAccumulatedfwFlux', nAccumulated)

           allocate(localArrayForPrecSum(nCellsSolve))
           allocate(localArrayForSSHSum(nCellsSolve))
           localArrayForPrecSum(:) = 0.0_RKIND
           localArrayForSSHSum(:) = 0.0_RKIND

           !precipitation flux and ssh 
           !$omp parallel
           !$omp do schedule(runtime)
           do iCell=1,nCellsSolve
               localArrayForPrecSum(iCell) = (rainFlux(iCell)+snowFlux(iCell)) * areaCell(iCell) ! unit: kg/s
               localArrayForSSHSum(iCell) = (sshNew(iCell)-sshOld(iCell))* areaCell(iCell) * rho_fw / dt ! unit: kg/s
           enddo
           !$omp end do
           !$omp end parallel

           allocate(localArrayForSalinitySum(nCellsSolve))
           localArrayForSalinitySum(:) = 0.0_RKIND
           
           !salinity flux
           !convert to 
           !$omp parallel
           !$omp do schedule(runtime) private(k)
           do iCell = 1,nCellsSolve
             do k = minLevelCell(iCell), maxLevelCell(iCell)
                localArrayForSalinitySum(iCell) = localArrayForSalinitySum(iCell) - & 
                     areaCell(iCell) * layerThickness(k,iCell) * & ! cell volume
                     (activeTracersNew(indexSalinity,k,iCell) - activeTracersOld(indexSalinity,k,iCell)) * & ! * ! cell salinity
                     rho_fw * 1.e3_RKIND / ocn_ref_salinity / dt ! F = - rho * h * 1/S * dS/dt 
             end do
           end do
           !$omp end do
           !$omp end parallel

            block_ptr => block_ptr % next
        end do ! block_ptr

        precFlux = mpas_global_sum(localArrayForPrecSum, domain % dminfo % comm)
        fwFlux = mpas_global_sum(localArrayForSSHSum, domain % dminfo % comm)
        salFlux = mpas_global_sum(localArrayForSalinitySum, domain % dminfo % comm)

        deallocate (localArrayForPrecSum)
        deallocate (localArrayForSSHSum)
        deallocate (localArrayForSalinitySum)

        avgprecFlux = (avgprecFlux * nAccumulated + precFlux) / ( nAccumulated + 1 )
        avgfwFlux = (avgfwFlux * nAccumulated + fwFlux) / ( nAccumulated + 1 )
        avgsalFlux = (avgsalFlux * nAccumulated + salFlux) / ( nAccumulated + 1 )

        nAccumulated = nAccumulated + 1

    end subroutine accumulate_mean_fw_fluxes!}}}


!***********************************************************************
!
!  routine update_scaling_factor
!
!> \brief  Update scaling factors for data freshwater fluxes 
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
!>  Update the net mass fluxes, and use it to scale the freshwater 
!>  fluxes based on the water flux balance relationship
!
!-----------------------------------------------------------------------

    subroutine update_scaling_factor(domain)!{{{

        type (domain_type), intent(inout) :: domain

        type (block_type), pointer :: block_ptr
        type (mpas_pool_type), pointer :: forcingPool

        real (kind=RKIND), pointer :: & 
                avgprecFlux, & 
                avgfwFlux, & 
                avgsalFlux, & 
                SFWFScalingFactor

        real (kind=RKIND), dimension(:), pointer :: & 
                totalAccumulatedprecFluxHistory, &
                totalAccumulatedfwFluxHistory, &
                totalAccumulatedsalFluxHistory

        integer, pointer :: nValidHistory, nAccumulated

        real (kind=RKIND), dimension(:), allocatable :: &
                tmpprecHistory, &
                tmpfwHistory, &
                tmpsalHistory 

        real (kind=RKIND) :: previousprecFlux, previousfwFlux, previoussalFlux
        real (kind=RKIND) :: precFluxMean, fwFluxMean, salFluxMean
        integer :: nHistory

        block_ptr => domain % blocklist
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_array(forcingPool, 'nAccumulatedfwFlux', nAccumulated)
        call mpas_pool_get_array(forcingPool, 'nValidTotalfwFluxHistory', nValidHistory)
        call mpas_pool_get_array(forcingPool, 'avgprecFlux', avgprecFlux)
        call mpas_pool_get_array(forcingPool, 'avgfwFlux', avgfwFlux)
        call mpas_pool_get_array(forcingPool, 'avgsalFlux', avgsalFlux)
        call mpas_pool_get_array(forcingPool, 'totalAccumulatedprecFluxHistory', totalAccumulatedprecFluxHistory)
        call mpas_pool_get_array(forcingPool, 'totalAccumulatedfwFluxHistory', totalAccumulatedfwFluxHistory)
        call mpas_pool_get_array(forcingPool, 'totalAccumulatedsalFluxHistory', totalAccumulatedsalFluxHistory)
        call mpas_pool_get_array(forcingPool, 'SFWFScalingFactor', SFWFScalingFactor)

        nHistory = config_sfwf_history_years

        previousprecFlux = totalAccumulatedprecFluxHistory(nValidHistory)
        previousfwFlux = totalAccumulatedfwFluxHistory(nValidHistory)
        previoussalFlux = totalAccumulatedsalFluxHistory(nValidHistory)

        if (nValidHistory == 0) then
            ! here followed the module mpas_ocn_scaled_dismf.F
            ! we keep zero as the first entry in totalAccumulated fluxes 
            nValidHistory = 1
        end if

        if (nHistory > 1) then 
            if (nValidHistory == nHistory) then
                ! we need to shift the history, since it's full
                allocate(tmpprecHistory(nHistory))
                allocate(tmpfwHistory(nHistory))
                allocate(tmpsalHistory(nHistory))
                tmpprecHistory(:) = totalAccumulatedprecFluxHistory(:)
                tmpfwHistory(:) = totalAccumulatedfwFluxHistory(:)
                tmpsalHistory(:) = totalAccumulatedsalFluxHistory(:)
                totalAccumulatedprecFluxHistory(1:nHistory - 1) = tmpprecHistory(2:nHistory)
                totalAccumulatedfwFluxHistory(1:nHistory - 1) = tmpfwHistory(2:nHistory)
                totalAccumulatedsalFluxHistory(1:nHistory - 1) = tmpsalHistory(2:nHistory)
            else
                ! the history isn't full yet, so we just add to the end
                nValidHistory = nValidHistory + 1
            end if
            ! add the new total, the previous total plus the new yearly average
            totalAccumulatedprecFluxHistory(nValidHistory) = previousprecFlux + avgprecFlux
            totalAccumulatedfwFluxHistory(nValidHistory) = previousfwFlux + avgfwFlux
            totalAccumulatedsalFluxHistory(nValidHistory) = previoussalFlux + avgsalFlux
        else 
            totalAccumulatedprecFluxHistory(nValidHistory) = avgprecFlux
            totalAccumulatedfwFluxHistory(nValidHistory) = avgfwFlux
            totalAccumulatedsalFluxHistory(nValidHistory) = avgsalFlux
        end if

        ! Scaling factors are ratio of (freshwater flux + salinity flux) to precipitation flux 
        ! Therefore is independent of the time invervals for the accumulated fluxes 
        precFluxMean = totalAccumulatedprecFluxHistory(nValidHistory) 
        fwFluxMean = totalAccumulatedfwFluxHistory(nValidHistory) 
        salFluxMean = totalAccumulatedsalFluxHistory(nValidHistory) 

        SFWFScalingFactor = SFWFScalingFactor - (fwFluxMean + salFluxMean) / precFluxMean

        ! reset yearly averaging of the precipitation and freshwater fluxes 
        nAccumulated = 0
        avgprecFlux = 0.0_RKIND 
        avgfwFlux = 0.0_RKIND
        avgsalFlux = 0.0_RKIND

    end subroutine update_scaling_factor!}}}

end module ocn_scaled_sfwf
