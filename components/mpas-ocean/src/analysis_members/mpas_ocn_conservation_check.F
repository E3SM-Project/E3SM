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

   logical :: landIceFreshwaterFluxesOn
   real (kind=RKIND) :: earthAreaE3SM ! surface area of earth in coupler
   integer :: accumulatedFluxCounter

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_bootstrap_conservation_check
!
!> \brief   Bootstrap MPAS-Ocean analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Ocean analysis member.
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
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
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

      ! taking PI from SHR_CONST_PI in share/util/shr_const_mod.F90 to match coupler
      !real (kind=RKIND), parameter :: piE3SM = 3.14159265358979323846_RKIND  ! pi
      real (kind=RKIND), parameter :: piE3SM = 3.141592653589793_RKIND
      ! taking earth radius from SHR_CONST_REARTH in share/util/shr_const_mod.F90 to match coupler
      real (kind=RKIND), parameter :: earthRadiusE3SM = 6.371229e6_RKIND ! radius of earth, m

      err = 0

      earthAreaE3SM = 4.0_RKIND * piE3SM * earthRadiusE3SM*earthRadiusE3SM

   end subroutine ocn_init_conservation_check!}}}

!***********************************************************************
!
!  routine ocn_precompute_conservation_check
!
!> \brief   Precompute MPAS-Ocean analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-Ocean analysis member.
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
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
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
         call mpas_log_write('CONSERVATION CHECKS')
         call mpas_log_write('date: '//trim(timeStr))
      endif

      accumulatedFluxCounter = accumulatedFluxCounter + 1

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
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
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
           accumulatedRemovedIceRunoffHeatFlux, &
           accumulatedIcebergHeatFlux, &
           accumulatedFrazilHeatFlux, &
           accumulatedLandIceHeatFlux, &
           accumulatedLandIceFrazilHeatFlux, &
           accumulatedRainTemperatureFlux, &
           accumulatedEvapTemperatureFlux, &
           accumulatedSeaIceTemperatureFlux, &
           accumulatedRiverRunoffTemperatureFlux, &
           accumulatedIcebergTemperatureFlux

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           forcingPool, &
           statePool, &
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
           removedIceRunoffFlux, &
           icebergHeatFlux, &
           landIceHeatFlux, &
           rainTemperatureFlux, &
           evapTemperatureFlux, &
           seaIceTemperatureFlux, &
           icebergTemperatureFlux, &
           accumulatedFrazilIceMassNew, &
           accumulatedFrazilIceMassOld, &
           accumulatedLandIceFrazilMassNew, &
           accumulatedLandIceFrazilMassOld


      real(kind=RKIND), dimension(:,:), pointer :: &
           activeTracersSurfaceFluxRunoff

      type (MPAS_timeInterval_type) :: &
           timeStepESMF

      real(kind=RKIND) :: &
           dt, dtAvg, v, A, s, c, explicitHeatFluxSum

      integer, pointer :: &
           nCellsSolve, &
           index_temperature_flux

      integer :: &
           iCell, &
           ierr

      integer, parameter :: &
           nSums = 19

      character(len=160) :: &
           m

      logical, pointer :: &
           activeTracersBulkRestoringPKG

      call mpas_pool_get_package(ocnPackages, 'activeTracersBulkRestoringPKGActive', activeTracersBulkRestoringPKG)

      call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=dt)
      dtAvg = dt * accumulatedFluxCounter

      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckEnergyAM", conservationCheckEnergyAMPool)

      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLatentHeatFlux",                 accumulatedLatentHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSensibleHeatFlux",               accumulatedSensibleHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLongWaveHeatFluxUp",             accumulatedLongWaveHeatFluxUp)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLongWaveHeatFluxDown",           accumulatedLongWaveHeatFluxDown)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSeaIceHeatFlux",                 accumulatedSeaIceHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedShortWaveHeatFlux",              accumulatedShortWaveHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedMeltingSnowHeatFlux",            accumulatedMeltingSnowHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedMeltingIceRunoffHeatFlux",       accumulatedMeltingIceRunoffHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedRemovedIceRunoffHeatFlux",       accumulatedRemovedIceRunoffHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedIcebergHeatFlux",                accumulatedIcebergHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedFrazilHeatFlux",                 accumulatedFrazilHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLandIceHeatFlux",                accumulatedLandIceHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLandIceFrazilHeatFlux",          accumulatedLandIceFrazilHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedRainTemperatureFlux",            accumulatedRainTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedEvapTemperatureFlux",            accumulatedEvapTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSeaIceTemperatureFlux",          accumulatedSeaIceTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedRiverRunoffTemperatureFlux",     accumulatedRiverRunoffTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedIcebergTemperatureFlux",         accumulatedIcebergTemperatureFlux)

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
            call MPAS_pool_get_subpool(block % structs, "state", statePool)

            call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

            call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
            call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
            call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)
            call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDown)
            call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)
            call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
            call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
            call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'removedIceRunoffFlux', removedIceRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'icebergHeatFlux', icebergHeatFlux)
            call mpas_pool_get_array(forcingPool, 'landIceHeatFlux', landIceHeatFlux)

            call mpas_pool_get_array(forcingPool, 'rainTemperatureFlux',            rainTemperatureFlux)
            call mpas_pool_get_array(forcingPool, 'evapTemperatureFlux',            evapTemperatureFlux)
            call mpas_pool_get_array(forcingPool, 'seaIceTemperatureFlux',          seaIceTemperatureFlux)
            call mpas_pool_get_array(forcingPool, 'icebergTemperatureFlux',         icebergTemperatureFlux)
            call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux',tracersSurfaceFluxPool)
            call mpas_pool_get_dimension(tracersSurfaceFluxPool, 'index_temperatureSurfaceFlux', index_temperature_flux)
            call mpas_pool_get_array(tracersSurfaceFluxPool, 'activeTracersSurfaceFluxRunoff', activeTracersSurfaceFluxRunoff)
            call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMassNew, 2)
            call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMassOld, 1)

            do iCell = 1, nCellsSolve

               sumArray(1) = sumArray(1) + areaCell(iCell) * latentHeatFlux(iCell)
               sumArray(2) = sumArray(2) + areaCell(iCell) * sensibleHeatFlux(iCell)
               sumArray(3) = sumArray(3) + areaCell(iCell) * longWaveHeatFluxUp(iCell)
               sumArray(4) = sumArray(4) + areaCell(iCell) * longWaveHeatFluxDown(iCell)
               sumArray(5) = sumArray(5) + areaCell(iCell) * seaIceHeatFlux(iCell)
               sumArray(6) = sumArray(6) + areaCell(iCell) * shortWaveHeatFlux(iCell)

               ! Snow and ice melt are a negative heat flux, because they cool the water
               sumArray(7) = sumArray(7) - areaCell(iCell) * snowFlux(iCell) * latent_heat_fusion_mks
               sumArray(8) = sumArray(8) - areaCell(iCell) * iceRunoffFlux(iCell) * latent_heat_fusion_mks
               sumArray(19) = sumArray(19) - areaCell(iCell) * removedIceRunoffFlux(iCell) * latent_heat_fusion_mks

               sumArray(9) = sumArray(9) + areaCell(iCell) * icebergHeatFlux(iCell)
               sumArray(11) = sumArray(11) + areaCell(iCell) * rainTemperatureFlux(iCell)
               sumArray(12) = sumArray(12) + areaCell(iCell) * evapTemperatureFlux(iCell)
               sumArray(13) = sumArray(13) + areaCell(iCell) * seaIceTemperatureFlux(iCell)
               ! river runoff temperature flux
               sumArray(14) = sumArray(14) + areaCell(iCell) * activeTracersSurfaceFluxRunoff(index_temperature_flux,iCell)
               sumArray(15) = sumArray(15) + areaCell(iCell) * icebergTemperatureFlux(iCell)

            enddo ! iCell

            if (config_use_frazil_ice_formation) then
               do iCell = 1, nCellsSolve
                  ! Frazil ice mass is negative. Negative coefficient makes heat
                  ! flux positive, because freezing ice releases heat.
                  sumArray(10) = sumArray(10) - areaCell(iCell) * config_frazil_heat_of_fusion &
                                 * (accumulatedFrazilIceMassNew(iCell) - accumulatedFrazilIceMassOld(iCell))/dt
               enddo
            end if

            if (landIceFreshwaterFluxesOn) then
               do iCell = 1, nCellsSolve
                  sumArray(17) = sumArray(17) + areaCell(iCell) * landIceHeatFlux(iCell)
               enddo
            end if

            if (landIceFreshwaterFluxesOn &
               .and.config_use_frazil_ice_formation &
               .and.config_frazil_under_land_ice) then
               call mpas_pool_get_array(statePool, 'accumulatedLandIceFrazilMass', accumulatedLandIceFrazilMassOld, 1)
               call mpas_pool_get_array(statePool, 'accumulatedLandIceFrazilMass', accumulatedLandIceFrazilMassNew, 2)
               do iCell = 1, nCellsSolve
                  ! Frazil ice mass is negative. Negative coefficient makes heat
                  ! flux positive, because freezing ice releases heat.
                  sumArray(18) = sumArray(18) - areaCell(iCell) * config_frazil_heat_of_fusion &
                                 * (accumulatedLandIceFrazilMassNew(iCell) - accumulatedLandIceFrazilMassOld(iCell))
               enddo
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
         accumulatedFrazilHeatFlux                 = accumulatedFrazilHeatFlux                 + sumArrayOut(10)
         accumulatedRainTemperatureFlux            = accumulatedRainTemperatureFlux            + sumArrayOut(11)
         accumulatedEvapTemperatureFlux            = accumulatedEvapTemperatureFlux            + sumArrayOut(12)
         accumulatedSeaIceTemperatureFlux          = accumulatedSeaIceTemperatureFlux          + sumArrayOut(13)
         accumulatedRiverRunoffTemperatureFlux     = accumulatedRiverRunoffTemperatureFlux     + sumArrayOut(14)
         accumulatedIcebergTemperatureFlux         = accumulatedIcebergTemperatureFlux         + sumArrayOut(15)
         accumulatedLandIceHeatFlux                = accumulatedLandIceHeatFlux                + sumArrayOut(17)
         accumulatedLandIceFrazilHeatFlux          = accumulatedLandIceFrazilHeatFlux          + sumArrayOut(18)
         accumulatedRemovedIceRunoffHeatFlux       = accumulatedRemovedIceRunoffHeatFlux       + sumArrayOut(19)

         ! cleanup
         deallocate(sumArray)
         deallocate(sumArrayOut)

      endif ! activeTracersBulkRestoringPKG

      !-------------------------------------------------------------
      ! Energy conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! Average the fluxes
         accumulatedLatentHeatFlux                 = accumulatedLatentHeatFlux                /accumulatedFluxCounter
         accumulatedSensibleHeatFlux               = accumulatedSensibleHeatFlux              /accumulatedFluxCounter
         accumulatedLongWaveHeatFluxUp             = accumulatedLongWaveHeatFluxUp            /accumulatedFluxCounter
         accumulatedLongWaveHeatFluxDown           = accumulatedLongWaveHeatFluxDown          /accumulatedFluxCounter
         accumulatedSeaIceHeatFlux                 = accumulatedSeaIceHeatFlux                /accumulatedFluxCounter
         accumulatedShortWaveHeatFlux              = accumulatedShortWaveHeatFlux             /accumulatedFluxCounter
         accumulatedMeltingSnowHeatFlux            = accumulatedMeltingSnowHeatFlux           /accumulatedFluxCounter
         accumulatedMeltingIceRunoffHeatFlux       = accumulatedMeltingIceRunoffHeatFlux      /accumulatedFluxCounter
         accumulatedIcebergHeatFlux                = accumulatedIcebergHeatFlux               /accumulatedFluxCounter
         accumulatedFrazilHeatFlux                 = accumulatedFrazilHeatFlux                /accumulatedFluxCounter
         accumulatedLandIceHeatFlux                = accumulatedLandIceHeatFlux               /accumulatedFluxCounter
         accumulatedRainTemperatureFlux            = accumulatedRainTemperatureFlux           /accumulatedFluxCounter
         accumulatedEvapTemperatureFlux            = accumulatedEvapTemperatureFlux           /accumulatedFluxCounter
         accumulatedSeaIceTemperatureFlux          = accumulatedSeaIceTemperatureFlux         /accumulatedFluxCounter
         accumulatedRiverRunoffTemperatureFlux     = accumulatedRiverRunoffTemperatureFlux    /accumulatedFluxCounter
         accumulatedIcebergTemperatureFlux         = accumulatedIcebergTemperatureFlux        /accumulatedFluxCounter
         accumulatedLandIceFrazilHeatFlux          = accumulatedLandIceFrazilHeatFlux         /accumulatedFluxCounter
         accumulatedRemovedIceRunoffHeatFlux       = accumulatedRemovedIceRunoffHeatFlux      /accumulatedFluxCounter

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
                accumulatedSeaIceHeatFlux           &
              + accumulatedShortWaveHeatFlux        &
              + accumulatedLongWaveHeatFluxUp       &
              + accumulatedLongWaveHeatFluxDown     &
              + accumulatedLatentHeatFlux           &
              + accumulatedSensibleHeatFlux         &
              + accumulatedMeltingSnowHeatFlux      &
              + accumulatedMeltingIceRunoffHeatFlux &
              + accumulatedIcebergHeatFlux          &
              + accumulatedFrazilHeatFlux           &
              + accumulatedLandIceHeatFlux          &
              + accumulatedRainTemperatureFlux   *rho_sw*cp_sw &
              + accumulatedEvapTemperatureFlux   *rho_sw*cp_sw &
              + accumulatedSeaIceTemperatureFlux *rho_sw*cp_sw &
              + accumulatedRiverRunoffTemperatureFlux   *rho_sw*cp_sw &
              + accumulatedIcebergTemperatureFlux*rho_sw*cp_sw
              ! note, accumulatedLandIceFrazilHeatFlux not added because already in accumulatedFrazilHeatFlux

         ! compute the final energy error
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "absoluteEnergyError", absoluteEnergyError)
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "relativeEnergyError", relativeEnergyError)

         absoluteEnergyError = netEnergyFlux * dtAvg - energyChange
         relativeEnergyError = absoluteEnergyError / (finalEnergy - 1.0_RKIND)

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            A = earthAreaE3SM
            call mpas_log_write('Conversion factors:')
            write(m,"('Earth area in E3SM:          A = ',es24.16,' m^2')") A; call mpas_log_write(m)
            write(m,"('Averaging time interval:    dt = ',f12.4,' s, ',f12.4,' days')") dtAvg, dtAvg/86400._RKIND; call mpas_log_write(m)
            write(m,"('Number of time steps         N = ',i12)") accumulatedFluxCounter; call mpas_log_write(m)
            write(m,"('density salt water      rho_sw = ',f12.4,' kg/m^3')") rho_sw; call mpas_log_write(m)
            write(m,"('specific heat salt water cp_sw = ',f12.4,' J/kg/K')") cp_sw; call mpas_log_write(m)
            write(m,"('latent heat of fusion      LHF = ',f12.4,' J/kg/K')") latent_heat_fusion_mks; call mpas_log_write(m)
            call mpas_log_write('')
            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write('ENERGY CONSERVATION CHECK')
            call mpas_log_write('')
            s = 0.0_RKIND
            call mpas_log_write('HEAT FLUXES: explicit')
            call mpas_log_write('MPAS-Ocean name            W (raw sum)    coupler name    short name         W/m^2 (flux/A)')
v=accumulatedFrazilHeatFlux          ; write(m,"('frazilMass * LHF /dt     ',es16.8,' o2x_Fioo_q      hfreeze         ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedSeaIceHeatFlux          ; write(m,"('seaIceHeatFlux           ',es16.8,' x2o_Fioi_melth  hmelt           ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedShortWaveHeatFlux       ; write(m,"('shortWaveHeatFlux        ',es16.8,' x2o_Foxx_swnet  hnetsw          ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedLongWaveHeatFluxDown    ; write(m,"('longWaveHeatFluxDown     ',es16.8,' x2o_Faxa_lwdn   hlwdn           ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedLongWaveHeatFluxUp      ; write(m,"('longWaveHeatFluxUp       ',es16.8,' x2o_Foxx_lwup   hlwup           ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedLatentHeatFlux          ; write(m,"('latentHeatFlux           ',es16.8,' x2o_Foxx_lat    hlatvap         ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedMeltingSnowHeatFlux     ; write(m,"('snowFlux * LHF           ',es16.8,'                 hlatfus         ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedMeltingIceRunoffHeatFlux; write(m,"('iceRunoffFlux * LHF      ',es16.8,'                 hiroff          ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
            if (landIceFreshwaterFluxesOn) then
v=accumulatedRemovedIceRunoffHeatFlux; write(m,"('removedIceRunoffFlux *LHF',es16.8,'                 hiroff                          ',f16.8)") v,v/A; call mpas_log_write(m);
v=accumulatedMeltingIceRunoffHeatFlux+accumulatedRemovedIceRunoffHeatFlux
                                       write(m,"('      SUM: ice runoff*LHF',es16.8,'                 hiroff SUM                      ',f16.8)") v,v/A; call mpas_log_write(m)
            endif
v=accumulatedSensibleHeatFlux        ; write(m,"('sensibleHeatFlux         ',es16.8,' x2o_Foxx_sen    hsen            ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedIcebergHeatFlux         ; write(m,"('icebergHeatFlux          ',es16.8,' x2o_Fioi_bergh  hberg           ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
            if (landIceFreshwaterFluxesOn) then
v=accumulatedLandIceHeatFlux         ; write(m,"('landIceHeatFlux          ',es16.8,'                                 ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
            end if
            if (landIceFreshwaterFluxesOn &
               .and.config_use_frazil_ice_formation &
               .and.config_frazil_under_land_ice) then
v=accumulatedLandIceFrazilHeatFlux   ; write(m,"('   landIceFrazilHeatFlux ',es16.8,' (already in hfreeze, do not sum                )',f16.8)") v,v/A; call mpas_log_write(m); ! no sum: s=s+v
            end if
                                       write(m,"('SUM EXPLICIT HEAT FLUXES ',es16.8,'                                 ',f16.8)") s, s/A; call mpas_log_write(m)
            explicitHeatFluxSum = s

            s = 0.0_RKIND
            call mpas_log_write('')
            call mpas_log_write('HEAT FLUXES: implicit through temperature contributions of mass fluxes')
            call mpas_log_write('MPAS-Ocean name            W (sum*rho*cp) coupler name    short name        W/m^2 (flux/A) ')
v=accumulatedRainTemperatureFlux   *rho_sw*cp_sw; write(m,"('RainTemperatureFlux      ',es16.8,'                                 ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedEvapTemperatureFlux   *rho_sw*cp_sw; write(m,"('EvapTemperatureFlux      ',es16.8,'                                 ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedSeaIceTemperatureFlux *rho_sw*cp_sw; write(m,"('SeaIceTemperatureFlux    ',es16.8,'                                 ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedRiverRunoffTemperatureFlux*rho_sw*cp_sw; write(m,"('RiverRunoffTempFlux      ',es16.8,'                                 ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
v=accumulatedIcebergTemperatureFlux*rho_sw*cp_sw; write(m,"('IcebergTemperatureFlux   ',es16.8,'                                 ',f16.8)") v,v/A; call mpas_log_write(m); s=s+v
            write(m,"('SUM IMPLICIT HEAT FLUXES ',es16.8,'                 hh20temp        ',f16.8,es16.8)") s, s/A; call mpas_log_write(m)
            s = s + explicitHeatFluxSum
            write(m,"('SUM IMP+EXP  HEAT FLUXES ',es16.8,'                                 ',f16.8,es16.8)") s, s/A; call mpas_log_write(m)

            call mpas_log_write(' ')
            call mpas_log_write('CHANGE IN HEAT CONTENT: computed from ocean domain')
            call mpas_log_write('                         J ')
            write(m,"('Initial energy         ',es16.8)") initialEnergy; call mpas_log_write(m)
            write(m,"('Final energy           ',es16.8)") finalEnergy; call mpas_log_write(m)
            write(m,"('Energy change          ',es16.8)") energyChange; call mpas_log_write(m)
            call mpas_log_write(' ')
            call mpas_log_write('ENERGY CONSERVATION SUMMARY')
            call mpas_log_write('                         J               W (J/dt)        W/m^2 (J/dt/A)')
            write(m,"('Energy change          ', 2es16.8,f16.8)")  energyChange, energyChange/dtAvg, energyChange/dtAvg/A ; call mpas_log_write(m)
            write(m,"('Net energy flux        ', 2es16.8,f16.8)")  netEnergyFlux * dtAvg, netEnergyFlux, netEnergyFlux/A ; call mpas_log_write(m)
            write(m,"('Absolute energy error  ', 2es16.8,f16.8)")  absoluteEnergyError, absoluteEnergyError/dtAvg, absoluteEnergyError/dtAvg/A ; call mpas_log_write(m)
            call mpas_log_write(' ')
            write(m,"('RELATIVE ENERGY ERROR =', es16.8)")  relativeEnergyError; call mpas_log_write(m)
            call mpas_log_write(' ')

         endif

      endif

    end subroutine energy_conservation

!***********************************************************************
!
!  routine mass_conservation
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
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
           accumulatedFrazilFlux, &
           accumulatedLandIceFlux, &
           accumulatedLandIceFrazilFlux

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           statePool, &
           forcingPool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           rainFlux, &
           snowFlux, &
           evaporationFlux, &
           seaIceFreshwaterFlux, &
           riverRunoffFlux, &
           iceRunoffFlux, &
           removedRiverRunoffFlux, &
           removedIceRunoffFlux, &
           icebergFreshwaterFlux, &
           accumulatedFrazilIceMassNew, &
           accumulatedFrazilIceMassOld, &
           accumulatedLandIceFrazilMassNew, &
           accumulatedLandIceFrazilMassOld, &
           landIceFreshwaterFlux

      real (kind=RKIND), dimension(:,:), pointer :: frazilLayerThicknessTendency

      type (MPAS_timeInterval_type) :: &
           timeStepESMF

      real(kind=RKIND) :: &
           dt, dtAvg, v, A, s, c

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell, ierr, k

      integer, parameter :: &
           nSums = 12

      integer, dimension(:), pointer :: minLevelCell, maxLevelCell

      character(len=160) :: &
           m

      logical, pointer :: &
           thicknessBulkPKGActive

      call mpas_pool_get_package(ocnPackages, 'thicknessBulkPKGActive', thicknessBulkPKGActive)

      call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=dt)
      dtAvg = dt * accumulatedFluxCounter

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
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedFrazilFlux",            accumulatedFrazilFlux)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedLandIceFlux",           accumulatedLandIceFlux)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedLandIceFrazilFlux",     accumulatedLandIceFrazilFlux)

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
            call MPAS_pool_get_subpool(block % structs, "state", statePool)

            call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
            call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            call mpas_pool_get_array(forcingPool, 'rainFlux',               rainFlux)
            call mpas_pool_get_array(forcingPool, 'snowFlux',               snowFlux)
            call mpas_pool_get_array(forcingPool, 'evaporationFlux',        evaporationFlux)
            call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux',   seaIceFreshwaterFlux)
            call mpas_pool_get_array(forcingPool, 'riverRunoffFlux',        riverRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'iceRunoffFlux',          iceRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'removedRiverRunoffFlux', removedRiverRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'removedIceRunoffFlux',   removedIceRunoffFlux)
            call mpas_pool_get_array(forcingPool, 'icebergFreshWaterFlux',  icebergFreshwaterFlux)
            call mpas_pool_get_array(forcingPool, 'landIceFreshwaterFlux',  landIceFreshwaterFlux)
            call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMassNew, 2)
            call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMassOld, 1)
            call mpas_pool_get_array(forcingPool, 'frazilLayerThicknessTendency', frazilLayerThicknessTendency)

            do iCell = 1, nCellsSolve
                sumArray( 1) = sumArray( 1) + areaCell(iCell) * rainFlux(iCell)
                sumArray( 2) = sumArray( 2) + areaCell(iCell) * snowFlux(iCell)
                sumArray( 3) = sumArray( 3) + areaCell(iCell) * evaporationFlux(iCell)
                sumArray( 4) = sumArray( 4) + areaCell(iCell) * seaIceFreshwaterFlux(iCell)
                sumArray( 5) = sumArray( 5) + areaCell(iCell) * riverRunoffFlux(iCell)
                sumArray( 6) = sumArray( 6) + areaCell(iCell) * iceRunoffFlux(iCell)
                sumArray( 7) = sumArray( 7) + areaCell(iCell) * removedRiverRunoffFlux(iCell)
                sumArray( 8) = sumArray( 8) + areaCell(iCell) * removedIceRunoffFlux(iCell)
                sumArray( 9) = sumArray( 9) + areaCell(iCell) * icebergFreshwaterFlux(iCell)
            enddo

            if (config_use_frazil_ice_formation) then
               do iCell = 1, nCellsSolve
                  do k = minLevelCell(iCell), maxLevelCell(iCell)
                     sumArray(10) = sumArray(10) + areaCell(iCell) * frazilLayerThicknessTendency(k,iCell) * rho_sw
                  enddo
               enddo
            end if

            if (landIceFreshwaterFluxesOn) then
               do iCell = 1, nCellsSolve
                  sumArray(11) = sumArray(11) + areaCell(iCell) * landIceFreshwaterFlux(iCell)
               enddo
            end if

            if (landIceFreshwaterFluxesOn &
               .and.config_use_frazil_ice_formation &
               .and.config_frazil_under_land_ice) then
               call mpas_pool_get_array(statePool, 'accumulatedLandIceFrazilMass', accumulatedLandIceFrazilMassOld, 1)
               call mpas_pool_get_array(statePool, 'accumulatedLandIceFrazilMass', accumulatedLandIceFrazilMassNew, 2)
               do iCell = 1, nCellsSolve
                  ! Frazil ice mass is negative. Negative coefficient makes heat
                  ! flux positive, because freezing ice releases heat.
                  sumArray(12) = sumArray(12) + areaCell(iCell) * &
                                 (accumulatedLandIceFrazilMassNew(iCell) - accumulatedLandIceFrazilMassOld(iCell))/dt
               enddo
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
         accumulatedFrazilFlux             = accumulatedFrazilFlux             + sumArrayOut(10)
         accumulatedLandIceFlux            = accumulatedLandIceFlux            + sumArrayOut(11)
         accumulatedLandIceFrazilFlux      = accumulatedLandIceFrazilFlux      + sumArrayOut(12)

         ! cleanup
         deallocate(sumArray)
         deallocate(sumArrayOut)

      endif ! thicknessBulkPKGActive

      !-------------------------------------------------------------
      ! Mass conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! Average the fluxes
         accumulatedRainFlux               = accumulatedRainFlux               /accumulatedFluxCounter
         accumulatedSnowFlux               = accumulatedSnowFlux               /accumulatedFluxCounter
         accumulatedEvaporationFlux        = accumulatedEvaporationFlux        /accumulatedFluxCounter
         accumulatedSeaIceFlux             = accumulatedSeaIceFlux             /accumulatedFluxCounter
         accumulatedRiverRunoffFlux        = accumulatedRiverRunoffFlux        /accumulatedFluxCounter
         accumulatedIceRunoffFlux          = accumulatedIceRunoffFlux          /accumulatedFluxCounter
         accumulatedRemovedRiverRunoffFlux = accumulatedRemovedRiverRunoffFlux /accumulatedFluxCounter
         accumulatedRemovedIceRunoffFlux   = accumulatedRemovedIceRunoffFlux   /accumulatedFluxCounter
         accumulatedIcebergFlux            = accumulatedIcebergFlux            /accumulatedFluxCounter
         accumulatedFrazilFlux             = accumulatedFrazilFlux             /accumulatedFluxCounter
         accumulatedLandIceFlux            = accumulatedLandIceFlux            /accumulatedFluxCounter
         accumulatedLandIceFrazilFlux      = accumulatedLandIceFrazilFlux      /accumulatedFluxCounter

         ! get initial mass
         call MPAS_pool_get_array(conservationCheckMassAMPool, "initialMass", initialMass)

         ! get final mass
         call MPAS_pool_get_array(conservationCheckMassAMPool, "finalMass", finalMass)
         call compute_total_mass(domain, finalMass)

         ! compute the mass change
         call MPAS_pool_get_array(conservationCheckMassAMPool, "massChange", massChange)
         massChange = finalMass - initialMass

         ! calculate the final net mass flux to the ice
         call MPAS_pool_get_array(conservationCheckMassAMPool, "netMassFlux", netMassFlux)

         netMassFlux = &
              + accumulatedRainFlux &
              + accumulatedSnowFlux &
              + accumulatedEvaporationFlux &
              + accumulatedSeaIceFlux &
              + accumulatedRiverRunoffFlux &
              + accumulatedIceRunoffFlux &
              + accumulatedIcebergFlux &
              + accumulatedFrazilFlux &
              + accumulatedLandIceFlux
              ! note, accumulatedLandIceFrazilFlux not added because already in accumulatedFrazilFlux

         ! compute the final mass error
         call MPAS_pool_get_array(conservationCheckMassAMPool, "absoluteMassError", absoluteMassError)
         call MPAS_pool_get_array(conservationCheckMassAMPool, "relativeMassError", relativeMassError)

         absoluteMassError = netMassFlux * dtAvg - massChange
         relativeMassError = absoluteMassError / (finalmass + 1.0_RKIND)

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write('MASS CONSERVATION CHECK')
            call mpas_log_write(' ')
            s = 0.0_RKIND
            A = earthAreaE3SM
            c = 1.0e6_RKIND/A
            call mpas_log_write('MASS FLUXES')
            call mpas_log_write('MPAS-Ocean name          kg/s (F)       coupler name     short name     kg/m^2/s*1e6 (F/A)')
v=accumulatedFrazilFlux            ; write(m,"('frazilFreshwaterFlux   ',es16.8,' o2x_Fioo_frazil  wfreeze        ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
v=accumulatedSeaIceFlux            ; write(m,"('seaIceFreshwaterFlux   ',es16.8,' x2o_Fioi_meltw   wmelt          ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
v=accumulatedRainFlux              ; write(m,"('rainFlux               ',es16.8,' x2o_Faxa_rain    wrain          ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
v=accumulatedSnowFlux              ; write(m,"('snowFlux               ',es16.8,' x2o_Faxa_snow    wsnow          ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
v=accumulatedIcebergFlux           ; write(m,"('icebergFreshwaterFlux  ',es16.8,' x2o_Fioi_bergw   wberg          ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
v=accumulatedEvaporationFlux       ; write(m,"('evaporationFlux        ',es16.8,' x2o_Foxx_evap    wevap          ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
v=accumulatedRiverRunoffFlux       ; write(m,"('riverRunoffFlux        ',es16.8,' x2o_Foxx_rofl    wrunoff        ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
            if (landIceFreshwaterFluxesOn) then
v=accumulatedRemovedRiverRunoffFlux; write(m,"('removedRiverRunoffFlux ',es16.8,'                  wrunoff                        ',f16.8)") v,v*c; call mpas_log_write(m);
v=accumulatedRiverRunoffFlux+accumulatedRemovedRiverRunoffFlux;
                                     write(m,"('      SUM: river runoff',es16.8,' x2o_Foxx_rofl    wrunoff SUM                    ',f16.8)") v,v*c; call mpas_log_write(m)
            endif
v=accumulatedIceRunoffFlux         ; write(m,"('iceRunoffFlux          ',es16.8,' x2o_Foxx_rofi    wfrzrof        ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
            if (landIceFreshwaterFluxesOn) then
v=accumulatedRemovedIceRunoffFlux  ; write(m,"('removedIceRunoffFlux   ',es16.8,'                  wfrzrof                        ',f16.8)") v,v*c; call mpas_log_write(m);
v=accumulatedIceRunoffFlux+accumulatedRemovedIceRunoffFlux;
                                     write(m,"('      SUM: ice runoff  ',es16.8,' x2o_Foxx_rofi    wfrzrof SUM                    ',f16.8)") v,v*c; call mpas_log_write(m)
v=accumulatedLandIceFlux           ; write(m,"('landIceFreshwaterFlux  ',es16.8,'                                 ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
            endif
            if (landIceFreshwaterFluxesOn &
               .and.config_use_frazil_ice_formation &
               .and.config_frazil_under_land_ice) then
v=accumulatedLandIceFrazilFlux     ; write(m,"('  landIceFrazilFlux    ',es16.8,' (already in wfreeze, do not sum                )',f16.8)") v,v*c; call mpas_log_write(m); ! no sum: s=s+v
            endif
                                     write(m,"('SUM VOLUME FLUXES      ',es16.8,'                                 ',f16.8,es16.8)") s, s*c; call mpas_log_write(m)

            call mpas_log_write(' ')
            call mpas_log_write('CHANGE IN MASS: computed from ocean domain')
            call mpas_log_write('                       kg ')
            write(m,"('Initial mass         ',es16.8)") initialMass; call mpas_log_write(m)
            write(m,"('Final mass           ',es16.8)") finalMass; call mpas_log_write(m)
            write(m,"('Mass change          ',es16.8)") massChange; call mpas_log_write(m)
            call mpas_log_write(' ')
            call mpas_log_write('MASS CONSERVATION SUMMARY')
            call mpas_log_write('                       kg              kg/s            kg/m^2/s*1e6')
            write(m,"('Mass change          ', 2es16.8, f16.8)")  massChange, massChange/dtAvg, massChange/dtAvg*c ; call mpas_log_write(m)
            write(m,"('Net mass flux        ', 2es16.8, f16.8)")  netMassFlux * dtAvg, netMassFlux, netMassFlux*c ; call mpas_log_write(m)
            write(m,"('Absolute mass error  ', 2es16.8, f16.8)")  absoluteMassError, absoluteMassError/dtAvg, absoluteMassError/dtAvg*c ; call mpas_log_write(m)
            call mpas_log_write(' ')
            write(m,"('RELATIVE MASS ERROR =', es16.8)")  relativeMassError; call mpas_log_write(m)
            call mpas_log_write(' ')

         endif

      endif

    end subroutine mass_conservation

!***********************************************************************
!
!  routine salt_conservation
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
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
           statePool, &
           forcingPool

      real(kind=RKIND), pointer :: &
           initialSalt, &
           finalSalt, &
           saltChange, &
           netSaltFlux, &
           absoluteSaltError, &
           relativeSaltError

      real(kind=RKIND), pointer :: &
           accumulatedSeaIceSalinityFlux, &
           accumulatedFrazilSalinityFlux, &
           accumulatedLandIceFrazilSalinityFlux

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           accumulatedFrazilIceSalinityNew, &
           accumulatedFrazilIceSalinityOld, &
           accumulatedLandIceFrazilMassNew, &
           accumulatedLandIceFrazilMassOld, &
           seaIceSalinityFlux

      type (MPAS_timeInterval_type) :: &
           timeStepESMF

      real(kind=RKIND) :: &
           dt, dtAvg, v, A, s, c

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell, &
           ierr

      integer, parameter :: &
           nSums = 3

      logical, pointer :: &
           activeTracersBulkRestoringPKG

      character(len=160) :: &
           m

      call mpas_pool_get_package(ocnPackages, 'activeTracersBulkRestoringPKGActive', activeTracersBulkRestoringPKG)

      call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=dt)
      dtAvg = dt * accumulatedFluxCounter

      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckSaltAM", conservationCheckSaltAMPool)

      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedSeaIceSalinityFlux", accumulatedSeaIceSalinityFlux)
      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedFrazilSalinityFlux", accumulatedFrazilSalinityFlux)
      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedLandIceFrazilSalinityFlux", accumulatedLandIceFrazilSalinityFlux)

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
            call MPAS_pool_get_subpool(block % structs, "state", statePool)

            call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
            call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)
            call mpas_pool_get_array(statePool, 'accumulatedFrazilIceSalinity', accumulatedFrazilIceSalinityOld, 1)
            call mpas_pool_get_array(statePool, 'accumulatedFrazilIceSalinity', accumulatedFrazilIceSalinityNew, 2)

            do iCell = 1, nCellsSolve

               ! salt flux to ocean
               sumArray(1) = sumArray(1) + areaCell(iCell) * seaIceSalinityFlux(iCell)
            enddo ! iCell

            if (config_use_frazil_ice_formation) then
               do iCell = 1, nCellsSolve
                  ! units: accumulatedFrazilIceSalinity is in kg/m^2, so divide by dt
                  sumArray(2) = sumArray(2) + areaCell(iCell) &
                                    * (accumulatedFrazilIceSalinityNew(iCell) - accumulatedFrazilIceSalinityOld(iCell))/dt
               enddo
            end if

            if (landIceFreshwaterFluxesOn &
               .and.config_use_frazil_ice_formation &
               .and.config_frazil_under_land_ice) then
               call mpas_pool_get_array(statePool, 'accumulatedLandIceFrazilMass', accumulatedLandIceFrazilMassOld, 1)
               call mpas_pool_get_array(statePool, 'accumulatedLandIceFrazilMass', accumulatedLandIceFrazilMassNew, 2)
               do iCell = 1, nCellsSolve
                  sumArray(3) = sumArray(3) + areaCell(iCell) &
                                 * (accumulatedLandIceFrazilMassNew(iCell) - accumulatedLandIceFrazilMassOld(iCell))/dt
               enddo
            end if

            block => block % next
         enddo

         ! perform the sums over processors
         call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

         ! accumulate fluxes
         accumulatedSeaIceSalinityFlux = accumulatedSeaIceSalinityFlux + sumArrayOut(1)
         accumulatedFrazilSalinityFlux = accumulatedFrazilSalinityFlux + sumArrayOut(2)
         accumulatedLandIceFrazilSalinityFlux = accumulatedLandIceFrazilSalinityFlux + sumArrayOut(3)

         ! cleanup
         deallocate(sumArray)
         deallocate(sumArrayOut)

      endif ! activeTracersBulkRestoringPKG

      !-------------------------------------------------------------
      ! Salt conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! Average the fluxes
         accumulatedSeaIceSalinityFlux = accumulatedSeaIceSalinityFlux /accumulatedFluxCounter
         accumulatedFrazilSalinityFlux = accumulatedFrazilSalinityFlux /accumulatedFluxCounter
         accumulatedLandIceFrazilSalinityFlux = accumulatedLandIceFrazilSalinityFlux /accumulatedFluxCounter

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

         netSaltFlux = accumulatedSeaIceSalinityFlux &
                     + accumulatedFrazilSalinityFlux &
                     + accumulatedLandIceFrazilSalinityFlux

         ! compute the final salt error
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "absoluteSaltError", absoluteSaltError)
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "relativeSaltError", relativeSaltError)

         absoluteSaltError = netSaltFlux * dtAvg - saltChange
         relativeSaltError = absoluteSaltError / (finalSalt - 1.0_RKIND)

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write('SALT CONSERVATION CHECK')
            call mpas_log_write(' ')
            A = earthAreaE3SM
            c = 1.0e6_RKIND/A
            s = 0.0_RKIND
            call mpas_log_write('SALT FLUXES')
            call mpas_log_write('MPAS-Ocean name          kg/s (F)       coupler name     short name     kg/m^2/s*1e6 (F/A)')
v=accumulatedFrazilSalinityFlux    ; write(m,"('frazilSalinityFlux     ',es16.8,'                                 ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
v=accumulatedSeaIceSalinityFlux    ; write(m,"('seaIceSalinityFlux     ',es16.8,' x2o_Fioi_salt    salt           ',f16.8)") v,v*c; call mpas_log_write(m); s=s+v
            if (landIceFreshwaterFluxesOn &
               .and.config_use_frazil_ice_formation &
               .and.config_frazil_under_land_ice) then
v=accumulatedLandIceFrazilSalinityFlux; write(m,"('LandIceFrazilSalinityFlux',es16.8,' (already in wmelt, do not sum)                  ',f16.8)") v,v*c; call mpas_log_write(m); !no sum: s=s+v
            end if

                                     write(m,"('SUM VOLUME FLUXES      ',es16.8,'                                 ',f16.8,es16.8)") s, s*c; call mpas_log_write(m)
            call mpas_log_write(' ')

            call mpas_log_write('CHANGE IN SALT: computed from ocean domain')
            call mpas_log_write('                       kg ')
            write(m,"('Initial salt         ',es16.8)") initialSalt; call mpas_log_write(m)
            write(m,"('Final salt           ',es16.8)") finalSalt; call mpas_log_write(m)
            write(m,"('Salt change          ',es16.8)") saltChange; call mpas_log_write(m)
            call mpas_log_write(' ')
            call mpas_log_write('SALT CONSERVATION SUMMARY')
            call mpas_log_write('                       kg              kg/s            kg/m^2/s*1e6')
            write(m,"('Salt change          ', 3es16.8)")  saltChange, saltChange/dtAvg, saltChange/dtAvg*c ; call mpas_log_write(m)
            write(m,"('Net salt flux        ', 3es16.8)")  netSaltFlux * dtAvg, netSaltFlux, netSaltFlux*c ; call mpas_log_write(m)
            write(m,"('Absolute salt error  ', 3es16.8)")  absoluteSaltError, absoluteSaltError/dtAvg, absoluteSaltError/dtAvg*c ; call mpas_log_write(m)
            call mpas_log_write(' ')
            write(m,"('RELATIVE SALT ERROR =', es16.8)")  relativeSaltError; call mpas_log_write(m)
            call mpas_log_write(' ')

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

      integer, dimension(:), pointer :: minLevelCell, maxLevelCell

      energy = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, "nVertLevels", nVertLevels)

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "state", statePool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call MPAS_pool_get_array(statePool, "layerThickness", layerThickness, 1)

         ! temperature
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

         do iCell = 1, nCellsSolve

            do k = minLevelCell(iCell), maxLevelCell(iCell)

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

      integer, dimension(:), pointer :: minLevelCell, maxLevelCell

      mass = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, "nVertLevels", nVertLevels)

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "state", statePool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call MPAS_pool_get_array(statePool, "layerThickness", layerThickness, 1)

         do iCell = 1, nCellsSolve

            do k = minLevelCell(iCell), maxLevelCell(iCell)

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

      integer, dimension(:), pointer :: minLevelCell, maxLevelCell

      salt = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, "nVertLevels", nVertLevels)

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "state", statePool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call MPAS_pool_get_array(statePool, "layerThickness", layerThickness, 1)

         ! temperature
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1) !!!! check time level

         do iCell = 1, nCellsSolve

            do k = minLevelCell(iCell), maxLevelCell(iCell)

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
           accumulatedShortWaveHeatFlux, &
           accumulatedSeaIceHeatFlux, &
           accumulatedMeltingSnowHeatFlux, &
           accumulatedMeltingIceRunoffHeatFlux, &
           accumulatedRemovedIceRunoffHeatFlux, &
           accumulatedIcebergHeatFlux, &
           accumulatedFrazilHeatFlux, &
           accumulatedLandIceHeatFlux, &
           accumulatedLandIceFrazilHeatFlux, &
           accumulatedRainTemperatureFlux, &
           accumulatedEvapTemperatureFlux, &
           accumulatedSeaIceTemperatureFlux, &
           accumulatedRiverRunoffTemperatureFlux, &
           accumulatedIcebergTemperatureFlux

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
           accumulatedFrazilFlux, &
           accumulatedLandIceFlux, &
           accumulatedLandIceFrazilFlux

      real(kind=RKIND), pointer :: &
           accumulatedSeaIceSalinityFlux, &
           accumulatedFrazilSalinityFlux, &
           accumulatedLandIceFrazilSalinityFlux

      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckEnergyAM", conservationCheckEnergyAMPool)

      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLatentHeatFlux",                 accumulatedLatentHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSensibleHeatFlux",               accumulatedSensibleHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLongWaveHeatFluxUp",             accumulatedLongWaveHeatFluxUp)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLongWaveHeatFluxDown",           accumulatedLongWaveHeatFluxDown)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSeaIceHeatFlux",                 accumulatedSeaIceHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedShortWaveHeatFlux",              accumulatedShortWaveHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedMeltingSnowHeatFlux",            accumulatedMeltingSnowHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedMeltingIceRunoffHeatFlux",       accumulatedMeltingIceRunoffHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedRemovedIceRunoffHeatFlux",       accumulatedRemovedIceRunoffHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedIcebergHeatFlux",                accumulatedIcebergHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedFrazilHeatFlux",                 accumulatedFrazilHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLandIceHeatFlux",                accumulatedLandIceHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLandIceFrazilHeatFlux",          accumulatedLandIceFrazilHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedRainTemperatureFlux",            accumulatedRainTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedEvapTemperatureFlux",            accumulatedEvapTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSeaIceTemperatureFlux",          accumulatedSeaIceTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedRiverRunoffTemperatureFlux",     accumulatedRiverRunoffTemperatureFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedIcebergTemperatureFlux",         accumulatedIcebergTemperatureFlux)

      accumulatedFluxCounter = 0
      accumulatedLatentHeatFlux                = 0.0_RKIND
      accumulatedSensibleHeatFlux              = 0.0_RKIND
      accumulatedLongWaveHeatFluxUp            = 0.0_RKIND
      accumulatedLongWaveHeatFluxDown          = 0.0_RKIND
      accumulatedSeaIceHeatFlux                = 0.0_RKIND
      accumulatedShortWaveHeatFlux             = 0.0_RKIND
      accumulatedMeltingSnowHeatFlux           = 0.0_RKIND
      accumulatedMeltingIceRunoffHeatFlux      = 0.0_RKIND
      accumulatedIcebergHeatFlux               = 0.0_RKIND
      accumulatedFrazilHeatFlux                = 0.0_RKIND
      accumulatedLandIceHeatFlux               = 0.0_RKIND
      accumulatedRainTemperatureFlux           = 0.0_RKIND
      accumulatedEvapTemperatureFlux           = 0.0_RKIND
      accumulatedSeaIceTemperatureFlux         = 0.0_RKIND
      accumulatedRiverRunoffTemperatureFlux    = 0.0_RKIND
      accumulatedIcebergTemperatureFlux        = 0.0_RKIND
      accumulatedLandIceFrazilHeatFlux         = 0.0_RKIND
      accumulatedRemovedIceRunoffHeatFlux      = 0.0_RKIND

      ! mass
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
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedFrazilFlux",            accumulatedFrazilFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedLandIceFlux",           accumulatedLandIceFlux)
       call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedLandIceFrazilFlux",     accumulatedLandIceFrazilFlux)

       accumulatedRainFlux               = 0.0_RKIND
       accumulatedSnowFlux               = 0.0_RKIND
       accumulatedEvaporationFlux        = 0.0_RKIND
       accumulatedSeaIceFlux             = 0.0_RKIND
       accumulatedRiverRunoffFlux        = 0.0_RKIND
       accumulatedIceRunoffFlux          = 0.0_RKIND
       accumulatedRemovedRiverRunoffFlux = 0.0_RKIND
       accumulatedRemovedIceRunoffFlux   = 0.0_RKIND
       accumulatedIcebergFlux            = 0.0_RKIND
       accumulatedFrazilFlux             = 0.0_RKIND
       accumulatedLandIceFlux            = 0.0_RKIND
       accumulatedLandIceFrazilFlux      = 0.0_RKIND

      ! salt
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckSaltAM", conservationCheckSaltAMPool)

      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedSeaIceSalinityFlux", accumulatedSeaIceSalinityFlux)
      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedFrazilSalinityFlux", accumulatedFrazilSalinityFlux)
      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedLandIceFrazilSalinityFlux", accumulatedLandIceFrazilSalinityFlux)

      accumulatedSeaIceSalinityFlux = 0.0_RKIND
      accumulatedFrazilSalinityFlux = 0.0_RKIND
      accumulatedLandIceFrazilSalinityFlux = 0.0_RKIND

    end subroutine reset_accumulated_variables

!***********************************************************************
!
!  routine ocn_restart_conservation_check
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
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
         landIceFreshwaterFluxesOn = .true.
      end if

   end subroutine ocn_restart_conservation_check!}}}

!***********************************************************************
!
!  routine ocn_finalize_conservation_check
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen, Adrian K. Turner
!> \date    August 2021
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
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
