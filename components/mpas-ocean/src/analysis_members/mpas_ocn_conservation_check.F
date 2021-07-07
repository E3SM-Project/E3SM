! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_conservation_check
!
!> \brief MPAS sea ice analysis mode member: conservation_check
!> \author Mark Petersen, Adrian K. Turner
!> \date   August 2021
!> \details
!>  MPAS sea ice analysis mode member: conservation_check
!>
!-----------------------------------------------------------------------

module ocn_conservation_check

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_log, only: mpas_log_write

   use ocn_constants
   use ocn_config

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

   public :: ocn_bootstrap_conservation_check, &
             ocn_init_conservation_check, &
             ocn_precompute_conservation_check, &
             ocn_compute_conservation_check, &
             ocn_restart_conservation_check, &
             ocn_finalize_conservation_check

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: landIceFluxesOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_bootstrap_conservation_check
!
!> \brief   Bootstrap MPAS-OCN analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all boostraps required for the
!>  MPAS-OCN analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_bootstrap_conservation_check(domain, err)!{{{

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

   end subroutine ocn_bootstrap_conservation_check!}}}

!***********************************************************************
!
!  routine ocn_init_conservation_check
!
!> \brief   Initialize MPAS-OCN analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-OCN analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_conservation_check(domain, err)!{{{

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

   end subroutine ocn_init_conservation_check!}}}

!***********************************************************************
!
!  routine ocn_precompute_conservation_check
!
!> \brief   Precompute MPAS-OCN analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-OCN analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_precompute_conservation_check(domain, timeLevel, err)!{{{

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

      type(MPAS_pool_type), pointer :: &
           conservationCheckAMPool, &
           conservationCheckEnergyAMPool, &
           conservationCheckMassAMPool, &
           conservationCheckSaltAMPool

      integer, pointer :: &
           performConservationPrecompute

      real(kind=RKIND), pointer :: &
           initialEnergy, &
           initialMass, &
           initialSalt

      err = 0

      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
      call MPAS_pool_get_array(conservationCheckAMPool, "performConservationPrecompute", performConservationPrecompute)

      if (performConservationPrecompute == 1) then

         ! zero the accumulated fluxes
         call reset_accumulated_variables(domain)

         ! initial total energy
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckEnergyAM", conservationCheckEnergyAMPool)
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "initialEnergy", initialEnergy)

         call compute_total_energy(domain, initialEnergy)

         ! initial total mass
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckMassAM", conservationCheckMassAMPool)
         call MPAS_pool_get_array(conservationCheckMassAMPool, "initialMass", initialMass)

         call compute_total_mass(domain, initialMass)

         ! initial total salt
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckSaltAM", conservationCheckSaltAMPool)
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "initialSalt", initialSalt)

         call compute_total_salt(domain, initialSalt)

         performConservationPrecompute = 0

      endif

   end subroutine ocn_precompute_conservation_check!}}}

!***********************************************************************
!
!  routine ocn_compute_conservation_check
!
!> \brief   Compute MPAS-OCN analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-OCN analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_conservation_check(domain, timeLevel, err)!{{{

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

      integer :: &
           ierr

      type(MPAS_pool_type), pointer :: &
           conservationCheckAMPool

      integer, pointer :: &
           performConservationPrecompute

      type(MPAS_Time_type) :: &
           currentTime

      character(len=strKIND) :: &
           timeStr

      err = 0

      if (config_AM_conservationCheck_write_to_logfile .and. &
          MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then
         call mpas_log_write('==========================================================')
         currentTime = MPAS_get_clock_time(domain % clock, MPAS_NOW, ierr=ierr)
         call MPAS_get_time(currentTime, dateTimeString=timeStr, ierr=ierr)
         call mpas_log_write(' Conservation checks: '//trim(timeStr))
      endif

      ! energy conservation check
      call energy_conservation(domain, err)

      ! mass conservation check
      call mass_conservation(domain, err)

      ! salt conservation check
      call salt_conservation(domain, err)

      if (config_AM_conservationCheck_write_to_logfile .and. &
         MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then
         call mpas_log_write('==========================================================')

         ! set precompute to happen next timestep
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
         call MPAS_pool_get_array(conservationCheckAMPool, "performConservationPrecompute", performConservationPrecompute)
         performConservationPrecompute = 1

      endif

   end subroutine ocn_compute_conservation_check!}}}

!***********************************************************************
!
!  routine energy_conservation
!
!> \brief   Compute MPAS-OCN analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-OCN analysis member.
!
!-----------------------------------------------------------------------

   subroutine energy_conservation(domain, err)

     use ocn_constants, only: &
          latent_heat_fusion_mks

      type(domain_type), intent(inout) :: &
           domain

      integer, intent(out) :: &
           err !< Output: error flag

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           conservationCheckEnergyAMPool

      real(kind=RKIND), pointer :: &
           initialEnergy, &
           finalEnergy, &
           energyChange, &
           netEnergyFlux, &
           absoluteEnergyError, &
           relativeEnergyError

      real(kind=RKIND), pointer :: &
           accumulatedLatentHeatFlux, &
           accumulatedSensibleHeatFlux, &
           accumulatedLongWaveHeatFluxUp, &
           accumulatedLongWaveHeatFluxDown, &
           accumulatedShortWaveHeatFlux, &
           accumulatedSeaIceHeatFlux, &
           accumulatedMeltingSnowHeatFlux, &
           accumulatedMeltingIceRunoffHeatFlux, &
           accumulatedIcebergHeatFlux, &
           accumulatedLandIceHeatFlux, &
           accumulatedRainTemperatureFlux, &
           accumulatedEvapTemperatureFlux, &
           accumulatedSeaIceTemperatureFlux, &
           accumulatedRiverRunoffTemperatureFlux, &
           accumulatedIcebergTemperatureFlux, &
           accumulatedTotalFreshWaterTemperatureFlux

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           forcingPool, &
           tracersSurfaceFluxPool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           latentHeatFlux, &
           sensibleHeatFlux, &
           longWaveHeatFluxUp, &
           longWaveHeatFluxDown, &
           seaIceHeatFlux, &
           shortWaveHeatFlux, &
           snowFlux, &
           iceRunoffFlux, &
           icebergHeatFlux, &
           landIceHeatFlux, &
           rainTemperatureFlux, &
           evapTemperatureFlux, &
           seaIceTemperatureFlux, &
           icebergTemperatureFlux, &
           totalFreshWaterTemperatureFlux

      real(kind=RKIND), dimension(:,:), pointer :: &
           activeTracersSurfaceFluxRunoff

      type (MPAS_timeInterval_type) :: &
           timeStepESMF

      real(kind=RKIND) :: &
           dt

      integer, pointer :: &
           nCellsSolve, &
           index_temperature_flux

      integer :: &
           iCell, &
           ierr

      integer, parameter :: &
           nSums = 7

      character(len=17) :: &
           formatString

      logical, pointer :: &
           activeTracersBulkRestoringPKG

      call mpas_pool_get_package(ocnPackages, 'activeTracersBulkRestoringPKGActive', activeTracersBulkRestoringPKG)

      call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=dt)

      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckEnergyAM", conservationCheckEnergyAMPool)

      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLatentHeatFlux",                 accumulatedLatentHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSensibleHeatFlux",               accumulatedSensibleHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLongWaveHeatFluxUp",             accumulatedLongWaveHeatFluxUp)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLongWaveHeatFluxDown",           accumulatedLongWaveHeatFluxDown)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSeaIceHeatFlux",                 accumulatedSeaIceHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedShortWaveHeatFlux",              accumulatedShortWaveHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedMeltingSnowHeatFlux",            accumulatedMeltingSnowHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedMeltingIceRunoffHeatFlux",       accumulatedMeltingIceRunoffHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedIcebergHeatFlux",                accumulatedIcebergHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLandIceHeatFlux",                accumulatedLandIceHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedRainTemperatureFlux",            accumulatedRainTemperatureFlux)           
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedEvapTemperatureFlux",            accumulatedEvapTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSeaIceTemperatureFlux",          accumulatedSeaIceTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedRiverRunoffTemperatureFlux",     accumulatedRiverRunoffTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedIcebergTemperatureFlux",         accumulatedIcebergTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedTotalFreshWaterTemperatureFlux", accumulatedTotalFreshWaterTemperatureFlux)

      !-------------------------------------------------------------
      ! Net heat flux to ice
      !-------------------------------------------------------------

      if (activeTracersBulkRestoringPKG) then

         allocate(sumArray(nSums))
         allocate(sumArrayOut(nSums))

         sumArray = 0.0_RKIND

         block => domain % blocklist
         do while (associated(block))

            call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

            call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
            call MPAS_pool_get_subpool(block % structs, "forcing", forcingPool)

            call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

            call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
            call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
            call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)
            call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDown)
            call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)
            call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
            call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
            call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'icebergHeatFlux', icebergHeatFlux)
            call mpas_pool_get_array(forcingPool, 'landIceHeatFlux', landIceHeatFlux)

            call mpas_pool_get_array(forcingPool, 'rainTemperatureFlux',            rainTemperatureFlux)
            call mpas_pool_get_array(forcingPool, 'evapTemperatureFlux',            evapTemperatureFlux)
            call mpas_pool_get_array(forcingPool, 'seaIceTemperatureFlux',          seaIceTemperatureFlux)
            call mpas_pool_get_array(forcingPool, 'icebergTemperatureFlux',         icebergTemperatureFlux)
            call mpas_pool_get_array(forcingPool, 'totalFreshWaterTemperatureFlux', totalFreshWaterTemperatureFlux)
            call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux',tracersSurfaceFluxPool)
            call mpas_pool_get_dimension(tracersSurfaceFluxPool, 'index_temperatureSurfaceFlux', index_temperature_flux)
            call mpas_pool_get_array(tracersSurfaceFluxPool, 'activeTracersSurfaceFluxRunoff', activeTracersSurfaceFluxRunoff)

            do iCell = 1, nCellsSolve

               sumArray(1) = sumArray(1) + areaCell(iCell) * latentHeatFlux(iCell)
               sumArray(2) = sumArray(2) + areaCell(iCell) * sensibleHeatFlux(iCell)
               sumArray(3) = sumArray(3) + areaCell(iCell) * longWaveHeatFluxUp(iCell)
               sumArray(4) = sumArray(4) + areaCell(iCell) * longWaveHeatFluxDown(iCell)
               sumArray(5) = sumArray(5) + areaCell(iCell) * seaIceHeatFlux(iCell)
               sumArray(6) = sumArray(6) + areaCell(iCell) * shortWaveHeatFlux(iCell)
               sumArray(7) = sumArray(7) + areaCell(iCell) * snowFlux(iCell) * latent_heat_fusion_mks
               sumArray(8) = sumArray(8) + areaCell(iCell) * iceRunoffFlux(iCell) * latent_heat_fusion_mks
               sumArray(9) = sumArray(9) + areaCell(iCell) * icebergHeatFlux(iCell)

               sumArray(11) = sumArray(11) + areaCell(iCell) * rainTemperatureFlux(iCell)
               sumArray(12) = sumArray(12) + areaCell(iCell) * evapTemperatureFlux(iCell)
               sumArray(13) = sumArray(13) + areaCell(iCell) * seaIceTemperatureFlux(iCell)
               ! river runoff temperature flux
               sumArray(14) = sumArray(14) + areaCell(iCell) * activeTracersSurfaceFluxRunoff(index_temperature_flux,iCell)
               sumArray(15) = sumArray(15) + areaCell(iCell) * icebergTemperatureFlux(iCell)
               sumArray(16) = sumArray(16) + areaCell(iCell) * totalFreshWaterTemperatureFlux(iCell)

            enddo ! iCell

            if (landIceFluxesOn) then
               do iCell = 1, nCellsSolve
                  sumArray(10) = sumArray(10) + areaCell(iCell) * landIceHeatFlux(iCell)
               enddo ! iCell
            end if

            block => block % next
         enddo

         ! perform the sums over processors
         call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

         ! accumulate fluxes
         accumulatedLatentHeatFlux                 = accumulatedLatentHeatFlux                 + sumArrayOut( 1)
         accumulatedSensibleHeatFlux               = accumulatedSensibleHeatFlux               + sumArrayOut( 2)
         accumulatedLongWaveHeatFluxUp             = accumulatedLongWaveHeatFluxUp             + sumArrayOut( 3)
         accumulatedLongWaveHeatFluxDown           = accumulatedLongWaveHeatFluxDown           + sumArrayOut( 4)
         accumulatedSeaIceHeatFlux                 = accumulatedSeaIceHeatFlux                 + sumArrayOut( 5)
         accumulatedShortWaveHeatFlux              = accumulatedShortWaveHeatFlux              + sumArrayOut( 6)
         accumulatedMeltingSnowHeatFlux            = accumulatedMeltingSnowHeatFlux            + sumArrayOut( 7)
         accumulatedMeltingIceRunoffHeatFlux       = accumulatedMeltingIceRunoffHeatFlux       + sumArrayOut( 8)
         accumulatedIcebergHeatFlux                = accumulatedIcebergHeatFlux                + sumArrayOut( 9)
         accumulatedLandIceHeatFlux                = accumulatedLandIceHeatFlux                + sumArrayOut(10)
         accumulatedRainTemperatureFlux            = accumulatedRainTemperatureFlux            + sumArrayOut(11)
         accumulatedEvapTemperatureFlux            = accumulatedEvapTemperatureFlux            + sumArrayOut(12)
         accumulatedSeaIceTemperatureFlux          = accumulatedSeaIceTemperatureFlux          + sumArrayOut(13)
         accumulatedRiverRunoffTemperatureFlux     = accumulatedRiverRunoffTemperatureFlux     + sumArrayOut(14)
         accumulatedIcebergTemperatureFlux         = accumulatedIcebergTemperatureFlux         + sumArrayOut(15)
         accumulatedTotalFreshWaterTemperatureFlux = accumulatedTotalFreshWaterTemperatureFlux + sumArrayOut(16)

         ! cleanup
         deallocate(sumArray)
         deallocate(sumArrayOut)

      endif ! activeTracersBulkRestoringPKG

      !-------------------------------------------------------------
      ! Energy conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! get initial energy
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "initialEnergy", initialEnergy)

         ! get final energy
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "finalEnergy", finalEnergy)
         call compute_total_energy(domain, finalEnergy)

         ! compute the energy change
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyChange", energyChange)
         energyChange = finalEnergy - initialEnergy

         ! calculate the final net energy flux to the ice
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "netEnergyFlux", netEnergyFlux)

         netEnergyFlux = &
                accumulatedLatentHeatFlux &
              + accumulatedSensibleHeatFlux &
              + accumulatedLongWaveHeatFluxUp &
              + accumulatedLongWaveHeatFluxDown &
              + accumulatedSeaIceHeatFlux &
              + accumulatedShortWaveHeatFlux &
              + accumulatedMeltingSnowHeatFlux &
              + accumulatedMeltingIceRunoffHeatFlux &
              + accumulatedIcebergHeatFlux &
              + accumulatedLandIceHeatFlux &
              + accumulatedTotalFreshWaterTemperatureFlux * latent_heat_fusion_mks

         ! compute the final energy error
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "absoluteEnergyError", absoluteEnergyError)
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "relativeEnergyError", relativeEnergyError)

         absoluteEnergyError = netEnergyFlux * dt - energyChange
         relativeEnergyError = absoluteEnergyError / (finalEnergy - 1.0_RKIND) ! why the minus 1????

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            formatString = "(a32,2x,1pes24.17)"

            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write(' Energy conservation check')
            call mpas_log_write(' ')
            call mpas_log_write(' Initial energy             (J) = $r', realArgs=(/initialEnergy/))
            call mpas_log_write(' Final energy               (J) = $r', realArgs=(/finalEnergy/))
            call mpas_log_write(' Energy change              (J) = $r', realArgs=(/energyChange/))
            call mpas_log_write(' ')
            call mpas_log_write(' Latent heat flux           (W) = $r', realArgs=(/accumulatedLatentHeatFlux/))
            call mpas_log_write(' Sensible heat flux         (W) = $r', realArgs=(/accumulatedSensibleHeatFlux/))
            call mpas_log_write(' Longwave heat flux up      (W) = $r', realArgs=(/accumulatedLongWaveHeatFluxUp/))
            call mpas_log_write(' Longwave heat flux down    (W) = $r', realArgs=(/accumulatedLongWaveHeatFluxDown/))
            call mpas_log_write(' Sea ice heat flux          (W) = $r', realArgs=(/accumulatedSeaIceHeatFlux/))
            call mpas_log_write(' Shortwave heat flux        (W) = $r', realArgs=(/accumulatedShortWaveHeatFlux/))
            call mpas_log_write(' MeltingSnowHeatFlux        (W) = $r', realArgs=(/accumulatedMeltingSnowHeatFlux/))
            call mpas_log_write(' MeltingIceRunoffHeatFlux   (W) = $r', realArgs=(/accumulatedMeltingIceRunoffHeatFlux/))
            call mpas_log_write(' IcebergHeatFlux            (W) = $r', realArgs=(/accumulatedIcebergHeatFlux/))
            call mpas_log_write(' LandIceHeatFlux            (W) = $r', realArgs=(/accumulatedLandIceHeatFlux/))
            call mpas_log_write(' RainTemperatureFlux     (Cm/s) = $r', realArgs=(/accumulatedRainTemperatureFlux/))
            call mpas_log_write(' EvapTemperatureFlux     (Cm/s) = $r', realArgs=(/accumulatedEvapTemperatureFlux/))
            call mpas_log_write(' SeaIceTemperatureFlux   (Cm/s) = $r', realArgs=(/accumulatedSeaIceTemperatureFlux/))
            call mpas_log_write(' RiverRunoffTemperatureFlux(Cm/s)=$r', realArgs=(/accumulatedRiverRunoffTemperatureFlux/))
            call mpas_log_write(' IcebergTemperatureFlux  (Cm/s) = $r', realArgs=(/accumulatedIcebergTemperatureFlux/))
            call mpas_log_write(' TotalFreshWaterTemperatureFlux = $r', realArgs=(/accumulatedTotalFreshWaterTemperatureFlux/))
            call mpas_log_write(' Net energy flux          (W) = $r', realArgs=(/netEnergyFlux/))
            call mpas_log_write(' Net energy flux          (J) = $r', realArgs=(/netEnergyFlux * dt/))
            call mpas_log_write(' ')
            call mpas_log_write(' Absolute energy error    (J) = $r', realArgs=(/absoluteEnergyError/))
            call mpas_log_write(' Absolute energy error    (W) = $r', realArgs=(/absoluteEnergyError / dt/))
            call mpas_log_write(' Relative energy error        = $r', realArgs=(/relativeEnergyError/))

         endif

      endif

    end subroutine energy_conservation

!***********************************************************************
!
!  routine mass_conservation
!
!> \brief   Compute MPAS-OCN analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-OCN analysis member.
!
!-----------------------------------------------------------------------

   subroutine mass_conservation(domain, err)

      type (domain_type), intent(inout) :: &
           domain

      integer, intent(out) :: &
           err !< Output: error flag

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           conservationCheckMassAMPool

      real(kind=RKIND), pointer :: &
           initialMass, &
           finalMass, &
           massChange, &
           netMassFlux, &
           absoluteMassError, &
           relativeMassError

      real(kind=RKIND), pointer :: &
           accumulatedRainFlux, &
           accumulatedSnowFlux, &
           accumulatedEvaporationFlux, &
           accumulatedSeaIceFlux, &
           accumulatedRiverRunoffFlux, &
           accumulatedIceRunoffFlux, &
           accumulatedRemovedRiverRunoffFlux, &
           accumulatedRemovedIceRunoffFlux, &
           accumulatedIcebergFlux, &
           accumulatedLandIceFlux

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           forcingPool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           rainFlux, &
           snowFlux, &
           evaporationFlux, &
           seaIceFreshWaterFlux, &
           riverRunoffFlux, &
           iceRunoffFlux, &
           removedRiverRunoffFlux, &
           removedIceRunoffFlux, &
           icebergFreshWaterFlux, &
           landIceFreshwaterFlux

      type (MPAS_timeInterval_type) :: &
           timeStepESMF

      real(kind=RKIND) :: &
           dt

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell, &
           ierr

      integer, parameter :: &
           nSums = 10

      character(len=17) :: &
           formatString

      logical, pointer :: &
           thicknessBulkPKGActive

      call mpas_pool_get_package(ocnPackages, 'thicknessBulkPKGActive', thicknessBulkPKGActive)

      call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=dt)

      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckMassAM", conservationCheckMassAMPool)

       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedRainFlux",              accumulatedRainFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedSnowFlux",              accumulatedSnowFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedEvaporationFlux",       accumulatedEvaporationFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedSeaIceFlux",            accumulatedSeaIceFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedRiverRunoffFlux",       accumulatedRiverRunoffFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedIceRunoffFlux",         accumulatedIceRunoffFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedRemovedRiverRunoffFlux",accumulatedRemovedRiverRunoffFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedRemovedIceRunoffFlux",  accumulatedRemovedIceRunoffFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedIcebergFlux",           accumulatedIcebergFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedLandIceFlux",           accumulatedLandIceFlux)

      !-------------------------------------------------------------
      ! Net mass flux to ocean
      !-------------------------------------------------------------

      if (thicknessBulkPKGActive) then

         allocate(sumArray(nSums))
         allocate(sumArrayOut(nSums))

         sumArray = 0.0_RKIND

         block => domain % blocklist
         do while (associated(block))

            call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

            call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
            call MPAS_pool_get_subpool(block % structs, "forcing", forcingPool)

            call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

            call mpas_pool_get_array(forcingPool, 'rainFlux',               rainFlux)
            call mpas_pool_get_array(forcingPool, 'snowFlux',               snowFlux)
            call mpas_pool_get_array(forcingPool, 'evaporationFlux',        evaporationFlux)
            call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux',   seaIceFreshWaterFlux)
            call mpas_pool_get_array(forcingPool, 'riverRunoffFlux',        riverRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'iceRunoffFlux',          iceRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'removedRiverRunoffFlux', removedRiverRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'removedIceRunoffFlux',   removedIceRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'icebergFreshWaterFlux',  icebergFreshWaterFlux)
            call mpas_pool_get_array(forcingPool, 'landIceFreshwaterFlux',  landIceFreshwaterFlux)

            do iCell = 1, nCellsSolve

                sumArray( 1) = sumArray( 1) + areaCell(iCell) * rainFlux(iCell)
                sumArray( 2) = sumArray( 2) + areaCell(iCell) * snowFlux(iCell)
                sumArray( 3) = sumArray( 3) + areaCell(iCell) * evaporationFlux(iCell)
                sumArray( 4) = sumArray( 4) + areaCell(iCell) * seaIceFreshWaterFlux(iCell)
                sumArray( 5) = sumArray( 5) + areaCell(iCell) * riverRunoffFlux(iCell)
                sumArray( 6) = sumArray( 6) + areaCell(iCell) * iceRunoffFlux(iCell)
                sumArray( 7) = sumArray( 7) + areaCell(iCell) * removedRiverRunoffFlux(iCell)
                sumArray( 8) = sumArray( 8) + areaCell(iCell) * removedIceRunoffFlux(iCell)
                sumArray( 9) = sumArray( 9) + areaCell(iCell) * icebergFreshWaterFlux(iCell)

            enddo ! iCell

            if (landIceFluxesOn) then
               do iCell = 1, nCellsSolve
                  sumArray(10) = sumArray(10) + areaCell(iCell) * landIceFreshwaterFlux(iCell)
               enddo ! iCell
            end if

            block => block % next
         enddo

         ! perform the sums over processors
         call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

         ! accumulate fluxes
         accumulatedRainFlux               = accumulatedRainFlux               + sumArrayOut( 1)
         accumulatedSnowFlux               = accumulatedSnowFlux               + sumArrayOut( 2)
         accumulatedEvaporationFlux        = accumulatedEvaporationFlux        + sumArrayOut( 3)
         accumulatedSeaIceFlux             = accumulatedSeaIceFlux             + sumArrayOut( 4)
         accumulatedRiverRunoffFlux        = accumulatedRiverRunoffFlux        + sumArrayOut( 5)
         accumulatedIceRunoffFlux          = accumulatedIceRunoffFlux          + sumArrayOut( 6)
         accumulatedRemovedRiverRunoffFlux = accumulatedRemovedRiverRunoffFlux + sumArrayOut( 7)
         accumulatedRemovedIceRunoffFlux   = accumulatedRemovedIceRunoffFlux   + sumArrayOut( 8) 
         accumulatedIcebergFlux            = accumulatedIcebergFlux            + sumArrayOut( 9)
         accumulatedLandIceFlux            = accumulatedLandIceFlux            + sumArrayOut(10)

         ! cleanup
         deallocate(sumArray)
         deallocate(sumArrayOut)

      endif ! thicknessBulkPKGActive

      !-------------------------------------------------------------
      ! Mass conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! get initial mass
         call MPAS_pool_get_array(conservationCheckMassAMPool, "initialMass", initialMass)

         ! get final mass
         call MPAS_pool_get_array(conservationCheckMassAMPool, "finalMass", finalMass)
         call compute_total_mass(domain, finalMass)

         ! compute the energy change
         call MPAS_pool_get_array(conservationCheckMassAMPool, "massChange", massChange)
         massChange = finalMass - initialMass

         ! calculate the final net energy flux to the ice
         call MPAS_pool_get_array(conservationCheckMassAMPool, "netMassFlux", netMassFlux)

         netMassFlux = &
              + accumulatedRainFlux &
              + accumulatedSnowFlux &
              + accumulatedEvaporationFlux &
              + accumulatedSeaIceFlux &
              + accumulatedRiverRunoffFlux &
              + accumulatedIceRunoffFlux &
              + accumulatedRemovedRiverRunoffFlux &
              + accumulatedRemovedIceRunoffFlux &
              + accumulatedIcebergFlux &
              + accumulatedLandIceFlux

         ! compute the final energy error
         call MPAS_pool_get_array(conservationCheckMassAMPool, "absoluteMassError", absoluteMassError)
         call MPAS_pool_get_array(conservationCheckMassAMPool, "relativeMassError", relativeMassError)

         absoluteMassError = netMassFlux * dt - massChange
         relativeMassError = absoluteMassError / (finalmass + 1.0_RKIND) ! why the plus 1????

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            formatString = "(a32,2x,1pe24.17)"

            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write(' Mass conservation check')
            call mpas_log_write(' ')
            call mpas_log_write(' Initial mass            (kg) = $r', realArgs=(/initialMass/))
            call mpas_log_write(' Final mass              (kg) = $r', realArgs=(/finalMass/))
            call mpas_log_write(' Mass change             (kg) = $r', realArgs=(/massChange/))
            call mpas_log_write(' ')
            call mpas_log_write(' Rain flux             (kg/s) = $r', realArgs=(/accumulatedRainFlux/))
            call mpas_log_write(' Snow flux             (kg/s) = $r', realArgs=(/accumulatedSnowFlux/))
            call mpas_log_write(' Evaporative flux      (kg/s) = $r', realArgs=(/accumulatedEvaporationFlux/))
            call mpas_log_write(' Sea ice FW flux       (kg/s) = $r', realArgs=(/accumulatedSeaIceFlux/))
            call mpas_log_write(' River runoff flux     (kg/s) = $r', realArgs=(/accumulatedRiverRunoffFlux/))
            call mpas_log_write(' Ice runoff flux       (kg/s) = $r', realArgs=(/accumulatedIceRunoffFlux/))
            call mpas_log_write(' RemovedRiverRunoffFlux(kg/s) = $r', realArgs=(/accumulatedRemovedRiverRunoffFlux/))
            call mpas_log_write(' RemovedIceRunoffFlux  (kg/s) = $r', realArgs=(/accumulatedRemovedIceRunoffFlux/))
            call mpas_log_write(' IcebergFlux           (kg/s) = $r', realArgs=(/accumulatedIcebergFlux/))
            call mpas_log_write(' LandIceFlux           (kg/s) = $r', realArgs=(/accumulatedLandIceFlux/))
            call mpas_log_write(' Net mass flux         (kg/s) = $r', realArgs=(/netMassFlux/))
            call mpas_log_write(' Net mass flux           (kg) = $r', realArgs=(/netMassFlux * dt/))
            call mpas_log_write(' ')
            call mpas_log_write(' Absolute mass error     (kg) = $r', realArgs=(/absoluteMassError/))
            call mpas_log_write(' Absolute mass error   (kg/s) = $r', realArgs=(/absoluteMassError / dt/))
            call mpas_log_write(' Relative mass error          = $r', realArgs=(/relativeMassError/))

         endif

      endif

    end subroutine mass_conservation

!***********************************************************************
!
!  routine salt_conservation
!
!> \brief   Compute MPAS-OCN analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-OCN analysis member.
!
!-----------------------------------------------------------------------

   subroutine salt_conservation(domain, err)

      type (domain_type), intent(inout) :: &
           domain

      integer, intent(out) :: &
           err !< Output: error flag

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           conservationCheckSaltAMPool, &
           meshPool, &
           forcingPool

      real(kind=RKIND), pointer :: &
           initialSalt, &
           finalSalt, &
           saltChange, &
           netSaltFlux, &
           absoluteSaltError, &
           relativeSaltError

      real(kind=RKIND), pointer :: &
           accumulatedSeaIceSalinityFlux

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           seaIceSalinityFlux

      type (MPAS_timeInterval_type) :: &
           timeStepESMF

      real(kind=RKIND) :: &
           dt

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell, &
           ierr

      integer, parameter :: &
           nSums = 1

      character(len=17) :: &
           formatString

      logical, pointer :: &
           activeTracersBulkRestoringPKG

      call mpas_pool_get_package(ocnPackages, 'activeTracersBulkRestoringPKGActive', activeTracersBulkRestoringPKG)

      call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=dt)

      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckSaltAM", conservationCheckSaltAMPool)

      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedSeaIceSalinityFlux", accumulatedSeaIceSalinityFlux)

      !-------------------------------------------------------------
      ! Net salt flux to ice
      !-------------------------------------------------------------

      if (activeTracersBulkRestoringPKG) then

         allocate(sumArray(nSums))
         allocate(sumArrayOut(nSums))

         sumArray = 0.0_RKIND

         block => domain % blocklist
         do while (associated(block))

            call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

            call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
            call MPAS_pool_get_subpool(block % structs, "forcing", forcingPool)

            call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
            call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)

            do iCell = 1, nCellsSolve

               ! salt flux to ocean
               sumArray(1) = sumArray(1) + areaCell(iCell) * seaIceSalinityFlux(iCell)

            enddo ! iCell

            block => block % next
         enddo

         ! perform the sums over processors
         call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

         ! accumulate fluxes
         accumulatedSeaIceSalinityFlux = accumulatedSeaIceSalinityFlux + sumArrayOut(1)

         ! cleanup
         deallocate(sumArray)
         deallocate(sumArrayOut)

      endif ! activeTracersBulkRestoringPKG

      !-------------------------------------------------------------
      ! Salt conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! get initial salt content
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "initialSalt", initialSalt)

         ! get final salt content
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "finalSalt", finalSalt)
         call compute_total_salt(domain, finalSalt)

         ! compute the salt content change
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "saltChange", saltChange)
         saltChange = finalSalt - initialSalt

         ! calculate the final net salt flux to the ice
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "netSaltFlux", netSaltFlux)

         netSaltFlux = &
              accumulatedSeaIceSalinityFlux

         ! compute the final salt error
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "absoluteSaltError", absoluteSaltError)
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "relativeSaltError", relativeSaltError)

         absoluteSaltError = netSaltFlux * dt - saltChange
         relativeSaltError = absoluteSaltError / (finalSalt - 1.0_RKIND) ! why the minus 1????

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            formatString = "(a32,2x,1pe24.17)"

            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write(' Salt conservation check')
            call mpas_log_write(' ')
            call mpas_log_write(' Initial salt            (kg) = $r', realArgs=(/initialSalt/))
            call mpas_log_write(' Final salt              (kg) = $r', realArgs=(/finalSalt/))
            call mpas_log_write(' Salt change             (kg) = $r', realArgs=(/saltChange/))
            call mpas_log_write(' ')
            call mpas_log_write(' Sea ice salinity flux (kg/s) = $r', realArgs=(/accumulatedSeaIceSalinityFlux/))
            call mpas_log_write(' Net salt flux         (kg/s) = $r', realArgs=(/netSaltFlux/))
            call mpas_log_write(' Net salt flux           (kg) = $r', realArgs=(/netSaltFlux * dt/))
            call mpas_log_write(' ')
            call mpas_log_write(' Absolute mass error     (kg) = $r', realArgs=(/absoluteSaltError/))
            call mpas_log_write(' Absolute mass error   (kg/s) = $r', realArgs=(/absoluteSaltError / dt/))
            call mpas_log_write(' Relative salt error          = $r', realArgs=(/relativeSaltError/))

         endif

      endif

    end subroutine salt_conservation

!***********************************************************************
!
!  routine compute_total_energy
!
!> \brief   Compute total energy of sea-ice system
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  Calculate the total energy of the sea-ice system
!
!-----------------------------------------------------------------------

    subroutine compute_total_energy(domain, totalEnergy)

      use ocn_constants, only: &
           rho_sw, &
           cp_sw

      type (domain_type), intent(inout) :: &
           domain

      real(kind=RKIND), intent(out) :: &
           totalEnergy

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           statePool, &
           tracersPool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell

      real(kind=RKIND), dimension(:,:), pointer :: &
           layerThickness

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           activeTracers

      real(kind=RKIND) :: &
           energy

      integer, pointer :: &
           nCellsSolve, &
           nVertLevels, &
           indexTemperature

      integer :: &
           iCell, &
           k

      energy = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, "nVertLevels", nVertLevels)

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "state", statePool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(statePool, "layerThickness", layerThickness, 1)

         ! temperature
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1) !!!! check time level

         do iCell = 1, nCellsSolve

            do k = 1, nVertLevels

               energy = energy + &
                    areaCell(iCell) * layerThickness(k,iCell) * & ! cell volume
                    activeTracers(indexTemperature,k,iCell) * &  ! cell temperature
                    rho_sw * cp_sw                                ! heat capacity

            enddo ! k

         enddo ! iCell

         block => block % next
      end do

      ! sum across processors
      call MPAS_dmpar_sum_real(domain % dminfo, energy, totalEnergy)

    end subroutine compute_total_energy

!***********************************************************************
!
!  routine compute_total_mass
!
!> \brief   Compute total mass of sea-ice system
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  Calculate the total mass of the sea-ice system
!
!-----------------------------------------------------------------------

    subroutine compute_total_mass(domain, totalMass)

      use ocn_constants, only: &
           rho_sw

      type (domain_type), intent(inout) :: &
           domain

      real(kind=RKIND), intent(out) :: &
           totalMass

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           statePool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell

      real(kind=RKIND), dimension(:,:), pointer :: &
           layerThickness

      real(kind=RKIND) :: &
           mass

      integer, pointer :: &
           nCellsSolve, &
           nVertLevels

      integer :: &
           iCell, &
           k

      mass = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, "nVertLevels", nVertLevels)

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "state", statePool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(statePool, "layerThickness", layerThickness, 1)

         do iCell = 1, nCellsSolve

            do k = 1, nVertLevels

               mass = mass + &
                    areaCell(iCell) * layerThickness(k,iCell) * & ! cell volume
                    rho_sw                                        ! density

            enddo ! k

         enddo ! iCell

         block => block % next
      end do

      ! sum across processors
      call MPAS_dmpar_sum_real(domain % dminfo, mass, totalMass)

    end subroutine compute_total_mass

!***********************************************************************
!
!  routine compute_total_salt
!
!> \brief   Compute total salt of sea-ice system
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  Calculate the total salt of the sea-ice system
!
!-----------------------------------------------------------------------

    subroutine compute_total_salt(domain, totalSalt)

      type (domain_type), intent(inout) :: &
           domain

      real(kind=RKIND), intent(out) :: &
           totalSalt

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           statePool, &
           tracersPool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell

      real(kind=RKIND), dimension(:,:), pointer :: &
           layerThickness

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           activeTracers

      real(kind=RKIND) :: &
           salt

      integer, pointer :: &
           nCellsSolve, &
           nVertLevels, &
           indexSalinity

      integer :: &
           iCell, &
           k

      salt = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, "nVertLevels", nVertLevels)

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "state", statePool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(statePool, "layerThickness", layerThickness, 1)

         ! temperature
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1) !!!! check time level

         do iCell = 1, nCellsSolve

            do k = 1, nVertLevels

               salt = salt + &
                    areaCell(iCell) * layerThickness(k,iCell) * & ! cell volume
                    activeTracers(indexSalinity,k,iCell)! * ! cell temperature
                    !?? ! salt capacity??

            enddo ! k

         enddo ! iCell

         block => block % next
      end do

      ! sum across processors
      call MPAS_dmpar_sum_real(domain % dminfo, salt, totalSalt)

    end subroutine compute_total_salt

!***********************************************************************
!
!  routine reset_accumulated_variables
!
!> \brief   Reset the accumulated fluxes
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details This routine resets accumulated fluxes after the
!> conservation calculation has been performed
!
!-----------------------------------------------------------------------

    subroutine reset_accumulated_variables(domain)

      type(domain_type), intent(inout) :: &
           domain

      type(MPAS_pool_type), pointer :: &
           conservationCheckEnergyAMPool, &
           conservationCheckMassAMPool, &
           conservationCheckSaltAMPool

      real(kind=RKIND), pointer :: &
           accumulatedLatentHeatFlux, &
           accumulatedSensibleHeatFlux, &
           accumulatedLongWaveHeatFluxUp, &
           accumulatedLongWaveHeatFluxDown, &
           accumulatedSeaIceHeatFlux, &
           accumulatedShortWaveHeatFlux, &
           accumulatedMeltingSnowHeatFlux, &
           accumulatedMeltingIceRunoffHeatFlux, &
           accumulatedLandIceHeatFlux, &
           accumulatedEvaporationFlux, &
           accumulatedSnowFlux, &
           accumulatedSeaIceFlux, &
           accumulatedRiverRunoffFlux, &
           accumulatedIceRunoffFlux, &
           accumulatedRainFlux, &
           accumulatedSeaIceSalinityFlux

      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckEnergyAM", conservationCheckEnergyAMPool)

      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLatentHeatFlux",       accumulatedLatentHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSensibleHeatFlux",     accumulatedSensibleHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLongWaveHeatFluxUp",   accumulatedLongWaveHeatFluxUp)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLongWaveHeatFluxDown", accumulatedLongWaveHeatFluxDown)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSeaIceHeatFlux",       accumulatedSeaIceHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedShortWaveHeatFlux",    accumulatedShortWaveHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedMeltingSnowHeatFlux",      accumulatedMeltingSnowHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedMeltingIceRunoffHeatFlux",      accumulatedMeltingIceRunoffHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLandIceHeatFlux",      accumulatedLandIceHeatFlux)

      accumulatedLatentHeatFlux       = 0.0_RKIND
      accumulatedSensibleHeatFlux     = 0.0_RKIND
      accumulatedLongWaveHeatFluxUp   = 0.0_RKIND
      accumulatedLongWaveHeatFluxDown = 0.0_RKIND
      accumulatedSeaIceHeatFlux       = 0.0_RKIND
      accumulatedShortWaveHeatFlux    = 0.0_RKIND
      accumulatedMeltingSnowHeatFlux      = 0.0_RKIND
      accumulatedMeltingIceRunoffHeatFlux      = 0.0_RKIND
      accumulatedLandIceHeatFlux      = 0.0_RKIND

      ! mass
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckMassAM", conservationCheckMassAMPool)

      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedEvaporationFlux",      accumulatedEvaporationFlux)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedSnowFlux",             accumulatedSnowFlux)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedSeaIceFlux", accumulatedSeaIceFlux)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedRiverRunoffFlux",      accumulatedRiverRunoffFlux)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedIceRunoffFlux",        accumulatedIceRunoffFlux)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedRainFlux",             accumulatedRainFlux)

      accumulatedEvaporationFlux      = 0.0_RKIND
      accumulatedSnowFlux             = 0.0_RKIND
      accumulatedSeaIceFlux = 0.0_RKIND
      accumulatedRiverRunoffFlux      = 0.0_RKIND
      accumulatedIceRunoffFlux        = 0.0_RKIND
      accumulatedRainFlux             = 0.0_RKIND

      ! salt
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckSaltAM", conservationCheckSaltAMPool)

      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedSeaIceSalinityFlux", accumulatedSeaIceSalinityFlux)

      accumulatedSeaIceSalinityFlux = 0.0_RKIND

    end subroutine reset_accumulated_variables

!***********************************************************************
!
!  routine ocn_restart_conservation_check
!
!> \brief   Save restart for MPAS-OCN analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-OCN analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_conservation_check(domain, err)!{{{

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

      if ( trim(config_land_ice_flux_mode) == 'standalone' .or. &
           trim(config_land_ice_flux_mode) == 'coupled' ) then
         landIceFluxesOn = .true.
      end if

   end subroutine ocn_restart_conservation_check!}}}

!***********************************************************************
!
!  routine ocn_finalize_conservation_check
!
!> \brief   Finalize MPAS-OCN analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-OCN analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_conservation_check(domain, err)!{{{

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

   end subroutine ocn_finalize_conservation_check!}}}

!-----------------------------------------------------------------------

end module ocn_conservation_check

! vim: foldmethod=marker
