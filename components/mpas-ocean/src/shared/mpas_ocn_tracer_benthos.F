! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_benthos
!
!> \brief MPAS ocean benthos
!> \author Nicole Jeffery
!> \date   07/15/2020
!> \details
!> This module contains routines for computing the benthic tracer tendencies
!> and the bottom boundary flux on ecosys tracers due to the benthos
!
!-----------------------------------------------------------------------

module ocn_tracer_benthos

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_timekeeping
   use mpas_forcing
   use mpas_stream_manager
   use mpas_log
   use ocn_constants
   use ocn_config
   use ocn_tracer_ecosys
   use ocn_framework_forcing
   use ocn_benthos_mod
   use ocn_benthos_parms

   !BEC
   !use bgc_mod
   !use bgc_parms

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

   public :: benthos_tracer_compute, &
!             ocn_tracer_ecosys_benthos_flux_compute,  &
             ocn_get_benthos_data,  &
             ocn_get_benthos_terra_data,  &
             ocn_benthos_forcing_write_restart,  &
             ocn_benthos_terra_forcing_write_restart,  &
             ocn_benthos_tracer_init, &
             ocn_benthos_tracer_finalize

   integer, public:: &
      numBenthicColumnsMax, &
      nBenthicCells

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!-----------------------------------------------------------------------
!  name the necessary benthos derived types
!  all of these are defined in benthos_mod
!-----------------------------------------------------------------------

  type(benthos_indices_type)    , public :: benthos_indices
  type(benthos_element_indices_type) , public :: benthos_element_indices
  type(benthos_input_type)      , public :: benthos_input
  type(benthos_forcing_type)    , public :: benthos_forcing
  type(benthos_output_type)     , public :: benthos_output
  type(benthos_diagnostics_type), public :: benthos_diagnostic_fields
  type(benthos_flux_diagnostics_type), public :: benthos_flux_diagnostic_fields

! hold indices in tracer pool corresponding to each eco tracer array
  type(benthos_indices_type)        , public :: benthosIndices
  type(benthos_element_indices_type), public :: elementIndices
  type(ocean_bgc_indices_type)      , public :: ecosysBenthosIndices

!***********************************************************************

contains

!***********************************************************************
!
!  routine benthos_tracer_compute
!
!> \brief   computes a seafloor tracer tendency due to benthos biogeochemistry
!> \author  Nicole Jeffery
!> \date    07/29/2020
!> \details
!>  This routine computes a seafloor tracer tendency due to the benthos biogeochemistry,
!>  and vertical transport and updates the benthos tracers
!
!-----------------------------------------------------------------------

   subroutine benthos_tracer_compute(benthosTracers, benthosPhPrev3d,benthosTracersStorage, &
      activeTracers, nBenthicTracers, benthosInterfaceLayerThickness, &
      nEcosysTracers, statePool, tendPool, diagnosticsPool, forcingPool, &
      nCellsSolve, nBenthicVertLevels, bottomDepth, density, maxLevelCell, &
      nVertLevels, layerThickness, indexTemperature, indexSalinity, benthosTracersTend,&
      dt,err, config_use_ecosysTracers, ecosysTracers, ecosysTracersTend)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! 1 dimensional arrays
      real (kind=RKIND), dimension(:), intent(in) :: &
         bottomDepth

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness, &         ! zMid?
         density

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in), optional :: &
         ecosysTracers             ! ocean BGC

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(inout), optional :: &
         ecosysTracersTend         ! ocean BGC tracer tendency

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         activeTracers     ! T and S

      real (kind=RKIND), dimension(:), intent(in) :: &
         benthosInterfaceLayerThickness ! m

      integer, dimension(:), intent(in) :: &
         maxLevelCell      ! last active level

      ! scalars
      integer, intent(in) :: &
         nBenthicTracers,    &
         nEcosysTracers,     &
         nCellsSolve,        &
         nBenthicVertLevels, &
         nVertLevels,        &
         indexTemperature,   &
         indexSalinity

      real (KIND=RKIND), intent(in) :: dt

      logical, intent(in) :: &
         config_use_ecosysTracers

      type (mpas_pool_type), intent(in) :: &
         statePool

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         benthosTracers        ! benthos BGC

      !
      ! two dimensional pointers
      !
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         benthosPhPrev3d, &    ! pH
         benthosTracersStorage  ! Deep storage

      type (mpas_pool_type), intent(inout) :: &
         tendPool,  &
         diagnosticsPool, &
         forcingPool

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
        benthosTracersTend

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer :: benthosAuxiliaryPool  ! additional forcing fields
      type (mpas_pool_type), pointer :: seafloorStoragePool   ! additional forcing fields
      type (mpas_pool_type), pointer :: seafloorTracersPool   ! additional forcing fields
      type (mpas_pool_type), pointer :: benthosDiagFieldsLevel1Pool  ! diagnostics
      type (mpas_pool_type), pointer :: benthosDiagFieldsLevel2Pool  ! diagnostics
      type (mpas_pool_type), pointer :: ecosysBenthosForcingPool
      type (mpas_pool_type), pointer :: oceanBenthosFluxPool
      type (mpas_pool_type), pointer :: oceanBottomBenthosFieldsPool
      type (mpas_pool_type), pointer :: oceanBenthosPrecipitationPool

      integer, pointer :: &
         config_set_benthos_flags

      logical, pointer ::  &
         config_use_carbonate_saturation, &
         config_use_secondary_reactions, &
         config_use_depth_dependent_porosity, &
         config_use_non_zero_diffusivity, &
         config_use_sedimentation, &
         config_use_fast_sedimentation, &
         config_use_positive_sedimentation, &
         config_use_bgc_sinking_flux, &
         config_use_step_initial_profiles, &
         config_use_benthic_reactions, &
         config_use_flux_correction, &
         config_use_ocean_conc, &
         config_use_deep_source, &
         config_use_constant_diffusivity, &
         config_use_benthos_monthly_forcing, &
         config_use_benthos_monthly_terra_forcing, &
         config_benthosTracers_diagnostic_fields_level1, &
         config_benthosTracers_diagnostic_fields_level2

      ! seafloorStorage

      integer, dimension(:), pointer :: &
         indexToBenthosCells

      !
      ! level 1 diagnostics
      !

      real (kind=RKIND), dimension(:,:), pointer :: &
         netBenthosReactionTend,   &
         netBenthosTransportTend,  &
         netBenthosTend,           &
         primarySourceTend,        &
         primarySinkTend,          &
         secondarySourceTend,      &
         secondarySinkTend,        &
         carbonateSourceTend,      &
         carbonateSinkTend,        &
         netBenthosTracersCell,    &
         netElements,              &
         netElementsReactionTend,  &
         netElementsTransportTend, &
         netElementsTend,          &
         netBenthosRelErrorTransport, &
         netBenthosElementsRelError, &
         netBenthosElementsRelErrorReactions

      real (kind=RKIND), dimension(:,:), pointer :: &
         benthos_calcToSed,      &
         benthos_pocToSed,       &
         benthos_ponToSed,       &
         benthos_popToSed,       &
         benthos_pfeToSed,       &
         benthos_bsiToSed,       &
         benthos_dustToSed,      &
         pressure          

      real (kind=RKIND), dimension(:), pointer :: &
         oceanBottomBenthosSio3, &
         oceanBottomBenthosPo4, &
         oceanBottomBenthosNh4, &
         oceanBottomBenthosNo3, &
         oceanBottomBenthosO2, &
         oceanBottomBenthosFe, &
         oceanBottomBenthosCo2, &
         oceanBottomBenthosMn, &
         oceanBottomBenthosHco3, &
         oceanBottomBenthosCo3, &
         oceanBottomBenthosSo4, &
         oceanBottomBenthosDIC, &
         oceanBottomBenthosAlk, &
         oceanBottomBenthosTemperature, &
         oceanBottomBenthosSalinity, &
         oceanBottomBenthosDepth, &
         oceanBottomBenthosPressure, &
         oceanBottomBenthosDensity, &
         oceanBenthosPOCFlux, &
         oceanBenthosPOCaFlux, &
         oceanBenthosPOPaFlux, &
         oceanBenthosPONaFlux, &
         oceanBenthosPOCbFlux, &
         oceanBenthosPOPbFlux, &
         oceanBenthosPONbFlux, &
         oceanBenthosPOCcFlux, &
         oceanBenthosPOPcFlux, &
         oceanBenthosPONcFlux, &
         oceanBenthosPICFlux, &
         oceanBenthosCaco3aFlux, &
         oceanBenthosCaco3bFlux, &
         oceanBenthosMgCaco3Flux, &
         oceanBenthosMno2aFlux, &
         oceanBenthosMno2bFlux, &
         oceanBenthosFeoh3Flux, &
         oceanBenthosFeoh3aFlux, &
         oceanBenthosFeoh3bFlux, &
         oceanBenthosFepaFlux, &
         oceanBenthosFepbFlux, &
         oceanBenthosBiogenicSilicaFlux, &
         oceanBenthosDustFlux, &
         oceanBenthosPrecipMno2a, &
         oceanBenthosPrecipMno2b, &
         oceanBenthosPrecipFeoh3a, &
         oceanBenthosPrecipFeoh3b, &
         oceanBenthosPrecipFepa, &
         oceanBenthosPrecipFepb, &
         oceanTerraLithogenicFlux, &
         oceanTerraBenthosPOCFlux, &
         totalSedimentMassFlux

      real (kind=RKIND) :: convertLengthMKStoCGS = 100.0_RKIND
      real (kind=RKIND) :: convertLengthCGStoMKS = 0.01_RKIND
      real (kind=RKIND) :: zTop, zBot

      integer :: iCell, iLevel, iTracer, numColumns, column, nElements, maxLevel, iCellBenthic

      err = 0

!      call mpas_log_write('benthos_compute start, mpas_ocn_tracer_benthos.F',MPAS_LOG_ERR)

      call mpas_timer_start("benthos source-sink")

      call mpas_pool_get_subpool(forcingPool, 'benthosAuxiliary', benthosAuxiliaryPool)
      call mpas_pool_get_subpool(forcingPool, 'oceanBenthosFlux', oceanBenthosFluxPool)
      call mpas_pool_get_subpool(forcingPool, 'oceanBottomBenthosFields', oceanBottomBenthosFieldsPool)
      call mpas_pool_get_subpool(forcingPool, 'oceanBenthosPrecipitation', oceanBenthosPrecipitationPool)
      call mpas_pool_get_subpool(statePool, 'seafloorStorage', seafloorStoragePool)
      call mpas_pool_get_subpool(statePool, 'seafloorTracers', seafloorTracersPool)
      call mpas_pool_get_subpool(diagnosticsPool, 'benthosDiagFieldsLevel1', benthosDiagFieldsLevel1Pool)
      call mpas_pool_get_subpool(diagnosticsPool, 'benthosDiagFieldsLevel2', benthosDiagFieldsLevel2Pool)
      call mpas_pool_get_array(benthosAuxiliaryPool,'indexToBenthosCells', indexToBenthosCells)

      call mpas_pool_get_config(ocnConfigs, 'config_set_benthos_flags',  config_set_benthos_flags)
      call mpas_pool_get_config(ocnConfigs, 'config_use_carbonate_saturation',  config_use_carbonate_saturation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_secondary_reactions',  config_use_secondary_reactions)
      call mpas_pool_get_config(ocnConfigs, 'config_use_ocean_conc',  config_use_ocean_conc)
      call mpas_pool_get_config(ocnConfigs, 'config_use_depth_dependent_porosity',  config_use_depth_dependent_porosity)
      call mpas_pool_get_config(ocnConfigs, 'config_use_non_zero_diffusivity',  config_use_non_zero_diffusivity)
      call mpas_pool_get_config(ocnConfigs, 'config_use_sedimentation',  config_use_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_fast_sedimentation',  config_use_fast_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_positive_sedimentation',  config_use_positive_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_bgc_sinking_flux',  config_use_bgc_sinking_flux)
      call mpas_pool_get_config(ocnConfigs, 'config_use_step_initial_profiles',  config_use_step_initial_profiles)
      call mpas_pool_get_config(ocnConfigs, 'config_use_benthic_reactions',  config_use_benthic_reactions)
      call mpas_pool_get_config(ocnConfigs, 'config_use_flux_correction',  config_use_flux_correction)
      call mpas_pool_get_config(ocnConfigs, 'config_use_deep_source',  config_use_deep_source)
      call mpas_pool_get_config(ocnConfigs, 'config_use_constant_diffusivity',  config_use_constant_diffusivity)
      call mpas_pool_get_config(ocnConfigs, 'config_use_benthos_monthly_forcing', config_use_benthos_monthly_forcing)
      call mpas_pool_get_config(ocnConfigs, 'config_use_benthos_monthly_terra_forcing', config_use_benthos_monthly_terra_forcing)
      call mpas_pool_get_config(ocnConfigs, 'config_benthosTracers_diagnostic_fields_level1', &
              config_benthosTracers_diagnostic_fields_level1)
      call mpas_pool_get_config(ocnConfigs, 'config_benthosTracers_diagnostic_fields_level2', &
              config_benthosTracers_diagnostic_fields_level2)

      if (config_use_ecosysTracers) then           
         call mpas_pool_get_subpool(forcingPool, 'ecosysBenthosForcing', ecosysBenthosForcingPool)
         call mpas_pool_get_array(ecosysBenthosForcingPool, 'benthos_calcToSed', benthos_calcToSed)
         call mpas_pool_get_array(ecosysBenthosForcingPool, 'benthos_pocToSed', benthos_pocToSed)
         call mpas_pool_get_array(ecosysBenthosForcingPool, 'benthos_pfeToSed', benthos_pfeToSed)
         call mpas_pool_get_array(ecosysBenthosForcingPool, 'benthos_popToSed', benthos_popToSed)
         call mpas_pool_get_array(ecosysBenthosForcingPool, 'benthos_ponToSed', benthos_ponToSed)
         call mpas_pool_get_array(ecosysBenthosForcingPool, 'benthos_bsiToSed', benthos_bsiToSed)
         call mpas_pool_get_array(ecosysBenthosForcingPool, 'benthos_dustToSed', benthos_dustToSed)
      end if

      call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)

      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPOCFlux', oceanBenthosPOCFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPOCaFlux', oceanBenthosPOCaFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPOPaFlux', oceanBenthosPOPaFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPONaFlux', oceanBenthosPONaFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPOCbFlux', oceanBenthosPOCbFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPOPbFlux', oceanBenthosPOPbFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPONbFlux', oceanBenthosPONbFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPOCcFlux', oceanBenthosPOCcFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPOPcFlux', oceanBenthosPOPcFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPONcFlux', oceanBenthosPONcFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosPICFlux', oceanBenthosPICFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosCaco3aFlux', oceanBenthosCaco3aFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosCaco3bFlux', oceanBenthosCaco3bFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosMgCaco3Flux', oceanBenthosMgCaco3Flux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosMno2aFlux', oceanBenthosMno2aFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosMno2bFlux', oceanBenthosMno2bFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosFeoh3Flux', oceanBenthosFeoh3Flux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosFeoh3aFlux', oceanBenthosFeoh3aFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosFeoh3bFlux', oceanBenthosFeoh3bFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosFepaFlux', oceanBenthosFepaFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosFepbFlux', oceanBenthosFepbFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosBiogenicSilicaFlux', oceanBenthosBiogenicSilicaFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanBenthosDustFlux', oceanBenthosDustFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'totalSedimentMassFlux', totalSedimentMassFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanTerraLithogenicFlux', oceanTerraLithogenicFlux)
      call mpas_pool_get_array(oceanBenthosFluxPool, 'oceanTerraBenthosPOCFlux', oceanTerraBenthosPOCFlux)

      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosSio3', oceanBottomBenthosSio3)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosPo4', oceanBottomBenthosPo4)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosNh4', oceanBottomBenthosNh4)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosNo3', oceanBottomBenthosNo3)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosO2', oceanBottomBenthosO2)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosFe', oceanBottomBenthosFe)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosCo2', oceanBottomBenthosCo2)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosMn', oceanBottomBenthosMn)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosHco3', oceanBottomBenthosHco3)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosCo3', oceanBottomBenthosCo3)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosSo4', oceanBottomBenthosSo4)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosDIC', oceanBottomBenthosDIC)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosAlk', oceanBottomBenthosAlk)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosTemperature', oceanBottomBenthosTemperature)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosSalinity', oceanBottomBenthosSalinity)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosPressure', oceanBottomBenthosPressure)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosDensity', oceanBottomBenthosDensity)
      call mpas_pool_get_array(oceanBottomBenthosFieldsPool, 'oceanBottomBenthosDepth', oceanBottomBenthosDepth)

      call mpas_pool_get_array(oceanBenthosPrecipitationPool, 'oceanBenthosPrecipMno2a', oceanBenthosPrecipMno2a)
      call mpas_pool_get_array(oceanBenthosPrecipitationPool, 'oceanBenthosPrecipMno2b', oceanBenthosPrecipMno2b)
      call mpas_pool_get_array(oceanBenthosPrecipitationPool, 'oceanBenthosPrecipFeoh3a', oceanBenthosPrecipFeoh3a)
      call mpas_pool_get_array(oceanBenthosPrecipitationPool, 'oceanBenthosPrecipFeoh3b', oceanBenthosPrecipFeoh3b)
      call mpas_pool_get_array(oceanBenthosPrecipitationPool, 'oceanBenthosPrecipFepa', oceanBenthosPrecipFepa)
      call mpas_pool_get_array(oceanBenthosPrecipitationPool, 'oceanBenthosPrecipFepb', oceanBenthosPrecipFepb)

!      call mpas_log_write('after ocnConfigs',MPAS_LOG_ERR)

      if (config_benthosTracers_diagnostic_fields_level1) then

!         call mpas_log_write('config_benthosTracers_diagnostic_fields_level1 = true')

         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netBenthosReactionTend', netBenthosReactionTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netBenthosTransportTend', netBenthosTransportTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netBenthosTend', netBenthosTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'primarySourceTend', primarySourceTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'primarySinkTend', primarySinkTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'secondarySourceTend', secondarySourceTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'secondarySinkTend', secondarySinkTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'carbonateSourceTend', carbonateSourceTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'carbonateSinkTend', carbonateSinkTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netBenthosTracersCell', netBenthosTracersCell)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netElements', netElements)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netElementsReactionTend', netElementsReactionTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netElementsTransportTend', netElementsTransportTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netElementsTend', netElementsTend)

         if (associated(netElements)) nElements = size(netElements, dim=1)
      else
         nElements = 7  ! default 
      end if

      if (config_benthosTracers_diagnostic_fields_level2) then
         call mpas_pool_get_array(benthosDiagFieldsLevel2Pool, 'netBenthosRelErrorTransport', netBenthosRelErrorTransport)
         call mpas_pool_get_array(benthosDiagFieldsLevel2Pool, 'netBenthosElementsRelErrorReactions', netBenthosElementsRelErrorReactions)
         call mpas_pool_get_array(benthosDiagFieldsLevel2Pool, 'netBenthosElementsRelError', netBenthosElementsRelError)
      end if

      numColumns = 1
      column = 1
!     DWJ 08/05/2016: This loop needs OpenMP added to it.

      ! allocate fields
      allocate ( benthos_input%benthosTracerBulk(nBenthicVertLevels+1, numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate ( benthos_input%deepStorage(numBenthicColumnsMax, benthos_tracer_cnt))
      allocate ( benthos_input%oceanBottomDepth(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomDensity(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomPressure(numBenthicColumnsMax))
      allocate ( benthos_input%PH_PREV_3D(nBenthicVertLevels+1, numBenthicColumnsMax))
      allocate ( benthos_input%bottom_level(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomTemperature(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomSalinity(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomPhosphate(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomAmmonium(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomNitrate(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomOxygen(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomIron(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomDIC(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomManganese(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomAlkalinity(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomSilicate(numBenthicColumnsMax))
      allocate ( benthos_input%oceanPOCFlux(numBenthicColumnsMax))
      allocate ( benthos_input%oceanPOPFlux(numBenthicColumnsMax))
      allocate ( benthos_input%oceanPONFlux(numBenthicColumnsMax))
      allocate ( benthos_input%oceanParticulateIronFlux(numBenthicColumnsMax))
      allocate ( benthos_input%oceanCalciteFlux(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBiogenicSilicaFlux(numBenthicColumnsMax))
      allocate ( benthos_input%oceanDustFlux(numBenthicColumnsMax))
      allocate ( benthos_input%lithogenicFlux(numBenthicColumnsMax))
      allocate ( benthos_input%terraPOCFlux(numBenthicColumnsMax))

      allocate ( benthos_output%benthosTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, benthos_tracer_cnt) )
      if (config_use_couple_ocean_ecosys) then
         allocate ( benthos_output%benthosEcosysTendencies(numBenthicColumnsMax, benthos_tracer_cnt) )
      end if
      allocate ( benthos_output%benthosTransportTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate ( benthos_output%benthosReactionTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate ( benthos_output%elementTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, element_tracer_cnt) )
      allocate ( benthos_output%elementReactionTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, element_tracer_cnt) )
      allocate ( benthos_output%elementTransportTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, element_tracer_cnt) )
      allocate ( benthos_output%elementDeepStorage(numBenthicColumnsMax, element_tracer_cnt) )

      allocate ( benthos_output%PH_PREV_3D(nBenthicVertLevels+1, numBenthicColumnsMax))
      allocate ( benthos_output%deepStorage(numBenthicColumnsMax, benthos_tracer_cnt))

      allocate ( benthos_output%oceanBottomTemperature(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomSalinity(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomSilicate(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomDensity(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomPressure(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomDepth(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomPhosphate(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomAmmonium(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomNitrate(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomOxygen(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomIron(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomDIC(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomAlkalinity(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomCarbonDioxide(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomHydrogenSulfide(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomMethane(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomBicarbonate(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomCarbonate(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomSulfate(numBenthicColumnsMax))
      allocate ( benthos_output%oceanBottomManganese(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPOCaFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPOPaFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPONaFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPOCbFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPOPbFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPONbFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPOCcFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPOPcFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPONcFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanCalciteFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanAragoniteFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanMgCalciteFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanManganeseOxideaFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanManganeseOxidebFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanIronOxyhydroxideaFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanIronOxyhydroxidebFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanIronBoundPhosphorusaFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanIronBoundPhosphorusbFlux(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPrecipManganeseOxidea(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPrecipManganeseOxideb(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPrecipIronOxyhydroxidea(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPrecipIronOxyhydroxideb(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPrecipIronBoundPhosphorusa(numBenthicColumnsMax))
      allocate ( benthos_output%oceanPrecipIronBoundPhosphorusb(numBenthicColumnsMax))
      allocate ( benthos_output%oceanTotalSedimentMassFlux(numBenthicColumnsMax))

      allocate ( benthos_forcing%deepStoragePrescribedFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate ( benthos_forcing%benthosOceanSedimentFlux(numBenthicColumnsMax, benthos_tracer_cnt) )

        !---------------------------------------------------------------------------
        !   allocate flux diagnostic output fields
        !   all fluxes are positive down
        !---------------------------------------------------------------------------

      allocate (benthos_flux_diagnostic_fields%deepStorageDiffVelFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%deepStorageFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%benthosOceanDiffVelFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%benthosOceanFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%elementSedimentation(numBenthicColumnsMax, element_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%elementBurialExchange(numBenthicColumnsMax, element_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%elementOceanExchange(numBenthicColumnsMax, element_tracer_cnt) )

      !---------------------------------------------------------------------------
      !   allocate diagnostic output fields
      !---------------------------------------------------------------------------

      ! 2D column integrated diagnostics

      allocate (benthos_diagnostic_fields%relErrorTransport(numBenthicColumnsMax,benthos_tracer_cnt))
      allocate (benthos_diagnostic_fields%relErrorElementsReactions(numBenthicColumnsMax,element_tracer_cnt))
      allocate (benthos_diagnostic_fields%relErrorElements(numBenthicColumnsMax, element_tracer_cnt))
      allocate (benthos_diagnostic_fields%diag_primarySourceTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_primarySinkTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_secondarySourceTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_secondarySinkTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_carbonateSourceTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_carbonateSinkTend(numBenthicColumnsMax, benthos_tracer_cnt) )

      allocate (benthos_diagnostic_fields%diag_netElements(numBenthicColumnsMax, element_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netBenthosTracersCell(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netBenthosReactionTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netBenthosTransportTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netBenthosTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netElementsReactionTend(numBenthicColumnsMax, element_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netElementsTransportTend(numBenthicColumnsMax, element_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netElementsTend(numBenthicColumnsMax, element_tracer_cnt) )

      do iCellBenthic=1,nBenthicCells
         iCell = indexToBenthosCells(iCellBenthic)
         benthos_input%bottom_level(column) = maxLevelCell(iCell)
         maxLevel = maxLevelCell(iCell)
         benthos_input%oceanBottomDepth(column) =  bottomDepth(iCell)

         if (config_use_benthos_monthly_terra_forcing) then
            benthos_input%lithogenicFlux(column) = oceanTerraLithogenicFlux(iCell)
            benthos_input%terraPOCFlux(column) = oceanTerraBenthosPOCFlux(iCell)
         else
            benthos_input%lithogenicFlux(column) = 0.0_RKIND
            benthos_input%terraPOCFlux(column) = 0.0_RKIND
         end if

         if (config_use_benthos_monthly_forcing) then
            benthos_input%oceanBottomDensity(column) = oceanBottomBenthosDensity(iCell)
            benthos_input%oceanBottomPressure(column) =  oceanBottomBenthosPressure(iCell)
            benthos_input%oceanBottomTemperature(column) = oceanBottomBenthosTemperature(iCell)
            benthos_input%oceanBottomSalinity(column) = oceanBottomBenthosSalinity(iCell)

            benthos_input%oceanBottomSilicate(column) = oceanBottomBenthosSio3(iCell)
            benthos_input%oceanBottomPhosphate(column) = oceanBottomBenthosPo4(iCell)
            benthos_input%oceanBottomAmmonium(column) = oceanBottomBenthosNh4(iCell)
            benthos_input%oceanBottomNitrate(column) = oceanBottomBenthosNo3(iCell)
            benthos_input%oceanBottomOxygen(column) = oceanBottomBenthosO2(iCell)
            benthos_input%oceanBottomIron(column) = oceanBottomBenthosFe(iCell)
            benthos_input%oceanBottomDIC(column) = oceanBottomBenthosDIC(iCell)
            benthos_input%oceanBottomAlkalinity(column) = oceanBottomBenthosAlk(iCell)
            benthos_input%oceanBottomManganese(column) = oceanBottomBenthosMn(iCell)

            benthos_input%oceanPOCFlux(column) = oceanBenthosPOCFlux(iCell)
            benthos_input%oceanParticulateIronFlux(column) = oceanBenthosFeoh3Flux(iCell)
            benthos_input%oceanCalciteFlux(column) = oceanBenthosPICFlux(iCell)
            benthos_input%oceanBiogenicSilicaFlux(column) = oceanBenthosBiogenicSilicaFlux(iCell)
            benthos_input%oceanDustFlux(column) = oceanBenthosDustFlux(iCell)
         else
            benthos_input%oceanBottomDensity(column) = density(maxLevel,iCell)
            benthos_input%oceanBottomPressure(column) =  pressure(maxLevel,iCell)
            benthos_input%oceanBottomTemperature(column) = activeTracers(indexTemperature, maxLevel, iCell)
            benthos_input%oceanBottomSalinity(column) = activeTracers(indexSalinity, maxLevel, iCell)
            benthos_input%oceanBottomManganese(column) = 4.0e-3_RKIND !Middag 2011 (3-5 nM Barents and Laptev)
            !
            ! will only be used if use_ocean_bottom_state = true, set in benthos_mod.F90
            !
            if (present(ecosysTracers)) then
               benthos_input%oceanBottomSilicate(column) = max(0.0_RKIND,ecosysTracers(ecosysBenthosIndices%sio3_ind, maxLevel, iCell)) !mmol Si /m^3
               benthos_input%oceanBottomPhosphate(column) = max(0.0_RKIND,ecosysTracers(ecosysBenthosIndices%po4_ind, maxLevel, iCell))
               benthos_input%oceanBottomAmmonium(column) = max(0.0_RKIND,ecosysTracers(ecosysBenthosIndices%nh4_ind, maxLevel, iCell))
               benthos_input%oceanBottomNitrate(column) = max(0.0_RKIND,ecosysTracers(ecosysBenthosIndices%no3_ind, maxLevel, iCell))

               benthos_input%oceanBottomOxygen(column) = max(0.0_RKIND,ecosysTracers(ecosysBenthosIndices%o2_ind, maxLevel, iCell))
               benthos_input%oceanBottomIron(column) = max(0.0_RKIND,ecosysTracers(ecosysBenthosIndices%fe_ind, maxLevel, iCell))
               benthos_input%oceanBottomDIC(column) = max(0.0_RKIND,ecosysTracers(ecosysBenthosIndices%dic_ind, maxLevel, iCell))
               benthos_input%oceanBottomAlkalinity(column) = max(0.0_RKIND,ecosysTracers(ecosysBenthosIndices%alk_ind, maxLevel, iCell))
               benthos_input%oceanPOCFlux(column) = max(0.0_RKIND,benthos_pocToSed(maxLevel,iCell)) ! units are mmolC/m3 m/s
               benthos_input%oceanParticulateIronFlux(column) = max(0.0_RKIND,benthos_pfeToSed(maxLevel,iCell)) ! units are mmolFe/m3 m/s
               benthos_input%oceanCalciteFlux(column) = max(0.0_RKIND,benthos_calcToSed(maxLevel,iCell)) ! units are mmolC/m3 m/s
               benthos_input%oceanBiogenicSilicaFlux(column) = 0.0_RKIND
               benthos_input%oceanDustFlux(column) = 0.0_RKIND
            else
               benthos_input%oceanBottomSilicate(column) = 0.0_RKIND
               benthos_input%oceanBottomPhosphate(column) = 0.0_RKIND
               benthos_input%oceanBottomAmmonium(column) = 0.0_RKIND
               benthos_input%oceanBottomNitrate(column) = 0.0_RKIND

               benthos_input%oceanBottomOxygen(column) = 0.0_RKIND
               benthos_input%oceanBottomIron(column) = 0.0_RKIND
               benthos_input%oceanBottomDIC(column) = 0.0_RKIND
               benthos_input%oceanBottomAlkalinity(column) = 0.0_RKIND
               benthos_input%oceanPOCFlux(column) = 0.0_RKIND
               benthos_input%oceanParticulateIronFlux(column) = 0.0_RKIND
               benthos_input%oceanCalciteFlux(column) = 0.0_RKIND
               !benthos_input%oceanPONFlux(column) = 0.0_RKIND
               !benthos_input%oceanPOPFlux(column) = 0.0_RKIND
            endif !
         endif  ! benthos_monthly_forcing
         benthos_input%deepStorage(column,benthos_indices%poca_ind) = benthosTracersStorage(benthosIndices%poca_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%pocb_ind) = benthosTracersStorage(benthosIndices%pocb_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%pocc_ind) = benthosTracersStorage(benthosIndices%pocc_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%pona_ind) = benthosTracersStorage(benthosIndices%pona_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%ponb_ind) = benthosTracersStorage(benthosIndices%ponb_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%ponc_ind) = benthosTracersStorage(benthosIndices%ponc_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%popa_ind) = benthosTracersStorage(benthosIndices%popa_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%popb_ind) = benthosTracersStorage(benthosIndices%popb_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%popc_ind) = benthosTracersStorage(benthosIndices%popc_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%o2_ind) = benthosTracersStorage(benthosIndices%o2_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%nh4_ind) = benthosTracersStorage(benthosIndices%nh4_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%h2po4_ind) = benthosTracersStorage(benthosIndices%h2po4_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%co2_ind) = benthosTracersStorage(benthosIndices%co2_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%no3_ind) = benthosTracersStorage(benthosIndices%no3_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%mno2a_ind) = benthosTracersStorage(benthosIndices%mno2a_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%mno2b_ind) = benthosTracersStorage(benthosIndices%mno2b_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%mn_ind) = benthosTracersStorage(benthosIndices%mn_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%feoh3a_ind) = benthosTracersStorage(benthosIndices%feoh3a_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%feoh3b_ind) = benthosTracersStorage(benthosIndices%feoh3b_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%fe_ind) = benthosTracersStorage(benthosIndices%fe_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%fepa_ind) = benthosTracersStorage(benthosIndices%fepa_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%fepb_ind) = benthosTracersStorage(benthosIndices%fepb_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%so4_ind) = benthosTracersStorage(benthosIndices%so4_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%h2s_ind) = benthosTracersStorage(benthosIndices%h2s_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%ch4_ind) = benthosTracersStorage(benthosIndices%ch4_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%dic_ind) = benthosTracersStorage(benthosIndices%dic_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%alk_ind) = benthosTracersStorage(benthosIndices%alk_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%hco3_ind) = benthosTracersStorage(benthosIndices%hco3_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%fes_ind) = benthosTracersStorage(benthosIndices%fes_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%fes2_ind) = benthosTracersStorage(benthosIndices%fes2_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%s_ind) = benthosTracersStorage(benthosIndices%s_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%caco3a_ind) = benthosTracersStorage(benthosIndices%caco3a_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%caco3b_ind) = benthosTracersStorage(benthosIndices%caco3b_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%co3_ind) = benthosTracersStorage(benthosIndices%co3_ind,iCell)
         benthos_input%deepStorage(column,benthos_indices%camgco3_ind) = benthosTracersStorage(benthosIndices%camgco3_ind,iCell)

         ! Check for negative tracers
         if (minval(benthosTracers(:,:,iCell)) .lt. -1.0E-11_RKIND) then
               call mpas_log_write('benthosTracers < 0 before SourceSink',MPAS_LOG_ERR)
               call mpas_log_write('iCell: $i',MPAS_LOG_ERR,intArgs=(/iCell/))
               call mpas_log_write('benthos_input%oceanBottomDensity(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomDensity(column)/))
               call mpas_log_write('benthos_input%oceanBottomTemperature(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomTemperature(column)/))
               call mpas_log_write('benthos_input%oceanBottomSalinity(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomSalinity(column)/))
               call mpas_log_write('benthos_input%oceanBottomPressure(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomPressure(column)/))
               call mpas_log_write('benthos_input%oceanBottomSilicate(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomSilicate(column)/))
               call mpas_log_write('benthos_input%oceanBottomPhosphate(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomPhosphate(column)/))
               call mpas_log_write('benthos_input%oceanBottomAmmonium(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomAmmonium(column)/))
               call mpas_log_write('benthos_input%oceanBottomNitrate(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomNitrate(column)/))
               call mpas_log_write('benthos_input%oceanBottomOxygen(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomOxygen(column)/))
               call mpas_log_write('benthos_input%oceanBottomIron(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomIron(column)/))
               call mpas_log_write('benthos_input%oceanBottomDIC(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomDIC(column)/))
               call mpas_log_write('benthos_input%oceanBottomAlkalinity(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomAlkalinity(column)/))
               call mpas_log_write('benthos_input%oceanBottomManganese(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanBottomManganese(column)/))
               call mpas_log_write('benthos_input%oceanPOCFlux(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanPOCFlux(column)/))
               call mpas_log_write('benthos_input%oceanParticulateIronFlux(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanParticulateIronFlux(column)/))
               call mpas_log_write('benthos_input%oceanCalciteFlux(column): $r',MPAS_LOG_ERR,realArgs=(/benthos_input%oceanCalciteFlux(column)/))
               do iTracer=1,nBenthicTracers
                  if (minval(benthosTracers(iTracer,:,iCell)) .lt. 0.0_RKIND) then
                     call mpas_log_write('iTracer: $i',MPAS_LOG_ERR,intArgs=(/iTracer/))
                     do iLevel=1,nBenthicVertLevels+1
                        call mpas_log_write('iLevel: $i',MPAS_LOG_ERR,intArgs=(/iLevel/))
                        call mpas_log_write('benthosTracers(iTracers,iLevel,iCell)): $r',MPAS_LOG_ERR,realArgs=(/benthosTracers(iTracer,iLevel,iCell)/))
                        call mpas_log_write('benthosTracersTend(iTracers,iLevel,iCell)): $r',MPAS_LOG_ERR,realArgs=(/benthosTracersTend(iTracer,iLevel,iCell)/))
                     end do
                  end if
               end do
         end if

         do iLevel=1,nBenthicVertLevels+1
            benthos_input%PH_PREV_3D(iLevel,column) = max(0.0_RKIND,benthosPhPrev3d(iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%poca_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%poca_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%pocb_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%pocb_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%pocc_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%pocc_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%pona_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%pona_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%ponb_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%ponb_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%ponc_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%ponc_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%popa_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%popa_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%popb_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%popb_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%popc_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%popc_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%o2_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%o2_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%nh4_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%nh4_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%h2po4_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%h2po4_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%co2_ind) = max(0.0_RKIND,max(0.0_RKIND,benthosTracers(benthosIndices%co2_ind,iLevel,iCell)))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%no3_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%no3_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%mno2a_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%mno2a_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%mno2b_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%mno2b_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%mn_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%mn_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%feoh3a_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%feoh3a_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%feoh3b_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%feoh3b_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%fe_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%fe_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%fepa_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%fepa_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%fepb_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%fepb_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%so4_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%so4_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%h2s_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%h2s_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%ch4_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%ch4_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%dic_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%dic_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%alk_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%alk_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%hco3_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%hco3_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%fes_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%fes_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%fes2_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%fes2_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%s_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%s_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%caco3a_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%caco3a_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%caco3b_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%caco3b_ind,iLevel,iCell))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%co3_ind) = max(0.0_RKIND,max(0.0_RKIND,benthosTracers(benthosIndices%co3_ind,iLevel,iCell)))
            benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%camgco3_ind) = max(0.0_RKIND,benthosTracers(benthosIndices%camgco3_ind,iLevel,iCell))
         end do
         call benthos_SourceSink(benthos_input, benthos_forcing, benthos_indices, &
            benthos_output, benthos_diagnostic_fields, nBenthicVertLevels, &
            numBenthicColumnsMax,dt,config_use_benthos_monthly_forcing,err)

         ! Define the tendencies
         benthosTracersStorage(benthosIndices%poca_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%poca_ind)
         benthosTracersStorage(benthosIndices%pocb_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%pocb_ind) 
         benthosTracersStorage(benthosIndices%pocc_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%pocc_ind)
         benthosTracersStorage(benthosIndices%pona_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%pona_ind) 
         benthosTracersStorage(benthosIndices%ponb_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%ponb_ind) 
         benthosTracersStorage(benthosIndices%ponc_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%ponc_ind) 
         benthosTracersStorage(benthosIndices%popa_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%popa_ind) 
         benthosTracersStorage(benthosIndices%popb_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%popb_ind) 
         benthosTracersStorage(benthosIndices%popc_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%popc_ind) 
         benthosTracersStorage(benthosIndices%o2_ind,iCell) =  benthos_output%deepStorage(column,benthos_indices%o2_ind)
         benthosTracersStorage(benthosIndices%nh4_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%nh4_ind)
         benthosTracersStorage(benthosIndices%h2po4_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%h2po4_ind) 
         benthosTracersStorage(benthosIndices%co2_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%co2_ind)
         benthosTracersStorage(benthosIndices%no3_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%no3_ind)
         benthosTracersStorage(benthosIndices%mno2a_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%mno2a_ind)
         benthosTracersStorage(benthosIndices%mno2b_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%mno2b_ind)
         benthosTracersStorage(benthosIndices%mn_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%mn_ind)
         benthosTracersStorage(benthosIndices%feoh3a_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%feoh3a_ind)
         benthosTracersStorage(benthosIndices%feoh3b_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%feoh3b_ind)
         benthosTracersStorage(benthosIndices%fe_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%fe_ind)
         benthosTracersStorage(benthosIndices%fepa_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%fepa_ind) 
         benthosTracersStorage(benthosIndices%fepb_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%fepb_ind)
         benthosTracersStorage(benthosIndices%so4_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%so4_ind)
         benthosTracersStorage(benthosIndices%h2s_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%h2s_ind)
         benthosTracersStorage(benthosIndices%ch4_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%ch4_ind)
         benthosTracersStorage(benthosIndices%dic_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%dic_ind)
         benthosTracersStorage(benthosIndices%alk_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%alk_ind)
         benthosTracersStorage(benthosIndices%hco3_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%hco3_ind)
         benthosTracersStorage(benthosIndices%fes_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%fes_ind)
         benthosTracersStorage(benthosIndices%fes2_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%fes2_ind)
         benthosTracersStorage(benthosIndices%s_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%s_ind)
         benthosTracersStorage(benthosIndices%caco3a_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%caco3a_ind)
         benthosTracersStorage(benthosIndices%caco3b_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%caco3b_ind)
         benthosTracersStorage(benthosIndices%co3_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%co3_ind)
         benthosTracersStorage(benthosIndices%camgco3_ind,iCell) = benthos_output%deepStorage(column,benthos_indices%camgco3_ind)

         oceanBottomBenthosSio3(iCell) = benthos_output%oceanBottomSilicate(column)
         oceanBottomBenthosPo4(iCell) = benthos_output%oceanBottomPhosphate(column)
         oceanBottomBenthosNh4(iCell) = benthos_output%oceanBottomAmmonium(column)
         oceanBottomBenthosNo3(iCell) = benthos_output%oceanBottomNitrate(column)
         oceanBottomBenthosO2(iCell) = benthos_output%oceanBottomOxygen(column)
         oceanBottomBenthosFe(iCell) = benthos_output%oceanBottomIron(column)
         oceanBottomBenthosCo2(iCell) = benthos_output%oceanBottomCarbonDioxide(column)
         oceanBottomBenthosMn(iCell) = benthos_output%oceanBottomManganese(column)
         oceanBottomBenthosHco3(iCell) = benthos_output%oceanBottomBicarbonate(column)
         oceanBottomBenthosCo3(iCell) = benthos_output%oceanBottomCarbonate(column)
         oceanBottomBenthosSo4(iCell) = benthos_output%oceanBottomSulfate(column)
         oceanBottomBenthosDIC(iCell) = benthos_output%oceanBottomDIC(column)
         oceanBottomBenthosAlk(iCell) = benthos_output%oceanBottomAlkalinity(column)

         oceanBottomBenthosTemperature(iCell) = benthos_output%oceanBottomTemperature(column)
         oceanBottomBenthosSalinity(iCell) = benthos_output%oceanBottomSalinity(column)
         oceanBottomBenthosDensity(iCell) = benthos_output%oceanBottomDensity(column)
         oceanBottomBenthosPressure(iCell) = benthos_output%oceanBottomPressure(column)
         oceanBottomBenthosDepth(iCell) = benthos_output%oceanBottomDepth(column)

         oceanBenthosPOCaFlux(iCell) = benthos_output%oceanPOCaFlux(column)
         oceanBenthosPOPaFlux(iCell) = benthos_output%oceanPOPaFlux(column)
         oceanBenthosPONaFlux(iCell) = benthos_output%oceanPONaFlux(column)
         oceanBenthosPOCbFlux(iCell) = benthos_output%oceanPOCbFlux(column)
         oceanBenthosPOPbFlux(iCell) = benthos_output%oceanPOPbFlux(column)
         oceanBenthosPONbFlux(iCell) = benthos_output%oceanPONbFlux(column)
         oceanBenthosPOCcFlux(iCell) = benthos_output%oceanPOCcFlux(column)
         oceanBenthosPOPcFlux(iCell) = benthos_output%oceanPOPcFlux(column)
         oceanBenthosPONcFlux(iCell) = benthos_output%oceanPONcFlux(column)
         oceanBenthosCaco3aFlux(iCell) = benthos_output%oceanCalciteFlux(column)
         oceanBenthosCaco3bFlux(iCell) = benthos_output%oceanAragoniteFlux(column)
         oceanBenthosMgCaco3Flux(iCell) = benthos_output%oceanMgCalciteFlux(column)
         oceanBenthosMno2aFlux(iCell) = benthos_output%oceanManganeseOxideaFlux(column)
         oceanBenthosMno2bFlux(iCell) = benthos_output%oceanManganeseOxidebFlux(column)
         oceanBenthosFeoh3aFlux(iCell) = benthos_output%oceanIronOxyhydroxideaFlux(column)
         oceanBenthosFeoh3bFlux(iCell) = benthos_output%oceanIronOxyhydroxidebFlux(column)
         oceanBenthosFepaFlux(iCell) = benthos_output%oceanIronBoundPhosphorusaFlux(column)
         oceanBenthosFepbFlux(iCell) = benthos_output%oceanIronBoundPhosphorusbFlux(column)

         oceanBenthosPrecipMno2a(iCell) = benthos_output%oceanPrecipManganeseOxidea(column)
         oceanBenthosPrecipMno2b(iCell) = benthos_output%oceanPrecipManganeseOxideb(column)
         oceanBenthosPrecipFeoh3a(iCell) = benthos_output%oceanPrecipIronOxyhydroxidea(column)
         oceanBenthosPrecipFeoh3b(iCell) = benthos_output%oceanPrecipIronOxyhydroxideb(column)
         oceanBenthosPrecipFepa(iCell) = benthos_output%oceanPrecipIronBoundPhosphorusa(column)
         oceanBenthosPrecipFepb(iCell) = benthos_output%oceanPrecipIronBoundPhosphorusb(column)
         totalSedimentMassFlux(iCell) = benthos_output%oceanTotalSedimentMassFlux(column)

         do iLevel=1,nBenthicVertLevels+1

            benthosPhPrev3d(iLevel,iCell) = benthos_output%PH_PREV_3D(iLevel,column)

            benthosTracersTend(benthosIndices%poca_ind,iLevel,iCell) =  &
            benthosTracersTend(benthosIndices%poca_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%poca_ind)*benthosInterfaceLayerThickness(iLevel)
            benthosTracersTend(benthosIndices%pocb_ind,iLevel,iCell) =  &
            benthosTracersTend(benthosIndices%pocb_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%pocb_ind)*benthosInterfaceLayerThickness(iLevel)
            benthosTracersTend(benthosIndices%pocc_ind,iLevel,iCell) =  &
            benthosTracersTend(benthosIndices%pocc_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%pocc_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%pona_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%pona_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%pona_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%ponb_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%ponb_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%ponb_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%ponc_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%ponc_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%ponc_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%popa_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%popa_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%popa_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%popb_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%popb_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%popb_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%popc_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%popc_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%popc_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%o2_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%o2_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%o2_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%nh4_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%nh4_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%nh4_ind)*benthosInterfaceLayerThickness(iLevel)

           benthosTracersTend(benthosIndices%h2po4_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%h2po4_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%h2po4_ind)*benthosInterfaceLayerThickness(iLevel)

           benthosTracersTend(benthosIndices%co2_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%co2_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%co2_ind)*benthosInterfaceLayerThickness(iLevel)

           benthosTracersTend(benthosIndices%no3_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%no3_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%no3_ind)*benthosInterfaceLayerThickness(iLevel)

           benthosTracersTend(benthosIndices%mno2a_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%mno2a_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%mno2a_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%mno2b_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%mno2b_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%mno2b_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%mn_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%mn_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%mn_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%feoh3a_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%feoh3a_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%feoh3a_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%feoh3b_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%feoh3b_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%feoh3b_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%fe_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%fe_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%fe_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%fepa_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%fepa_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%fepa_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%fepb_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%fepb_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%fepb_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%so4_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%so4_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%so4_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%h2s_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%h2s_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%h2s_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%ch4_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%ch4_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%ch4_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%hco3_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%hco3_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%hco3_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%fes_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%fes_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%fes_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%fes2_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%fes2_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%fes2_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%s_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%s_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%s_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%caco3a_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%caco3a_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%caco3a_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%caco3b_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%caco3b_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%caco3b_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%co3_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%co3_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%co3_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%camgco3_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%camgco3_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%camgco3_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%dic_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%dic_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%dic_ind)*benthosInterfaceLayerThickness(iLevel)
           benthosTracersTend(benthosIndices%alk_ind,iLevel,iCell) =  &
           benthosTracersTend(benthosIndices%alk_ind,iLevel,iCell)   &
               + benthos_output%benthosTendencies(iLevel,column,benthos_indices%alk_ind)*benthosInterfaceLayerThickness(iLevel)

         end do
         !
         ! level 1 diagnostics
         !

         if (config_benthosTracers_diagnostic_fields_level1) then
            netBenthosReactionTend(benthosIndices%poca_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%poca_ind)
            netBenthosReactionTend(benthosIndices%pocb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%pocb_ind)
            netBenthosReactionTend(benthosIndices%pocc_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%pocc_ind)
            netBenthosReactionTend(benthosIndices%pona_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%pona_ind)
            netBenthosReactionTend(benthosIndices%ponb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%ponb_ind)
            netBenthosReactionTend(benthosIndices%ponc_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%ponc_ind)
            netBenthosReactionTend(benthosIndices%popa_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%popa_ind)
            netBenthosReactionTend(benthosIndices%popb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%popb_ind)
            netBenthosReactionTend(benthosIndices%popc_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%popc_ind)
            netBenthosReactionTend(benthosIndices%o2_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%o2_ind)
            netBenthosReactionTend(benthosIndices%nh4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%nh4_ind)
            netBenthosReactionTend(benthosIndices%h2po4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%h2po4_ind)
            netBenthosReactionTend(benthosIndices%co2_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%co2_ind)
            netBenthosReactionTend(benthosIndices%no3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%no3_ind)
            netBenthosReactionTend(benthosIndices%mno2a_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%mno2a_ind)
            netBenthosReactionTend(benthosIndices%mno2b_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%mno2b_ind)
            netBenthosReactionTend(benthosIndices%mn_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%mn_ind)
            netBenthosReactionTend(benthosIndices%feoh3a_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%feoh3a_ind)
            netBenthosReactionTend(benthosIndices%feoh3b_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%feoh3b_ind)
            netBenthosReactionTend(benthosIndices%fe_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%fe_ind)
            netBenthosReactionTend(benthosIndices%fepa_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%fepa_ind)
            netBenthosReactionTend(benthosIndices%fepb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%fepb_ind)
            netBenthosReactionTend(benthosIndices%so4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%so4_ind)
            netBenthosReactionTend(benthosIndices%h2s_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%h2s_ind)
            netBenthosReactionTend(benthosIndices%ch4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%ch4_ind)
            netBenthosReactionTend(benthosIndices%hco3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%hco3_ind)
            netBenthosReactionTend(benthosIndices%fes_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%fes_ind)
            netBenthosReactionTend(benthosIndices%fes2_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%fes2_ind)
            netBenthosReactionTend(benthosIndices%s_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%s_ind)
            netBenthosReactionTend(benthosIndices%caco3a_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%caco3a_ind)
            netBenthosReactionTend(benthosIndices%caco3b_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%caco3b_ind)
            netBenthosReactionTend(benthosIndices%co3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%co3_ind)
            netBenthosReactionTend(benthosIndices%camgco3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%camgco3_ind)
            netBenthosReactionTend(benthosIndices%dic_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%dic_ind)
            netBenthosReactionTend(benthosIndices%alk_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosReactionTend(column,benthos_indices%alk_ind)

            netBenthosTransportTend(benthosIndices%poca_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%poca_ind)
            netBenthosTransportTend(benthosIndices%pocb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%pocb_ind)
            netBenthosTransportTend(benthosIndices%pocc_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%pocc_ind)
            netBenthosTransportTend(benthosIndices%pona_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%pona_ind)
            netBenthosTransportTend(benthosIndices%ponb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%ponb_ind)
            netBenthosTransportTend(benthosIndices%ponc_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%ponc_ind)
            netBenthosTransportTend(benthosIndices%popa_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%popa_ind)
            netBenthosTransportTend(benthosIndices%popb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%popb_ind)
            netBenthosTransportTend(benthosIndices%popc_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%popc_ind)
            netBenthosTransportTend(benthosIndices%o2_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%o2_ind)
            netBenthosTransportTend(benthosIndices%nh4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%nh4_ind)
            netBenthosTransportTend(benthosIndices%h2po4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%h2po4_ind)
            netBenthosTransportTend(benthosIndices%co2_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%co2_ind)
            netBenthosTransportTend(benthosIndices%no3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%no3_ind)
            netBenthosTransportTend(benthosIndices%mno2a_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%mno2a_ind)
            netBenthosTransportTend(benthosIndices%mno2b_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%mno2b_ind)
            netBenthosTransportTend(benthosIndices%mn_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%mn_ind)
            netBenthosTransportTend(benthosIndices%feoh3a_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%feoh3a_ind)
            netBenthosTransportTend(benthosIndices%feoh3b_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%feoh3b_ind)
            netBenthosTransportTend(benthosIndices%fe_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%fe_ind)
            netBenthosTransportTend(benthosIndices%fepa_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%fepa_ind)
            netBenthosTransportTend(benthosIndices%fepb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%fepb_ind)
            netBenthosTransportTend(benthosIndices%so4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%so4_ind)
            netBenthosTransportTend(benthosIndices%h2s_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%h2s_ind)
            netBenthosTransportTend(benthosIndices%ch4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%ch4_ind)
            netBenthosTransportTend(benthosIndices%hco3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%hco3_ind)
            netBenthosTransportTend(benthosIndices%fes_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%fes_ind)
            netBenthosTransportTend(benthosIndices%fes2_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%fes2_ind)
            netBenthosTransportTend(benthosIndices%s_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%s_ind)
            netBenthosTransportTend(benthosIndices%caco3a_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%caco3a_ind)
            netBenthosTransportTend(benthosIndices%caco3b_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%caco3b_ind)
            netBenthosTransportTend(benthosIndices%co3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%co3_ind)
            netBenthosTransportTend(benthosIndices%camgco3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%camgco3_ind)
            netBenthosTransportTend(benthosIndices%dic_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%dic_ind)
            netBenthosTransportTend(benthosIndices%alk_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTransportTend(column,benthos_indices%alk_ind)


            netBenthosTend(benthosIndices%poca_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%poca_ind)
            netBenthosTend(benthosIndices%pocb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%pocb_ind)
            netBenthosTend(benthosIndices%pocc_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%pocc_ind)
            netBenthosTend(benthosIndices%pona_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%pona_ind)
            netBenthosTend(benthosIndices%ponb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%ponb_ind)
            netBenthosTend(benthosIndices%ponc_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%ponc_ind)
            netBenthosTend(benthosIndices%popa_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%popa_ind)
            netBenthosTend(benthosIndices%popb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%popb_ind)
            netBenthosTend(benthosIndices%popc_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%popc_ind)
            netBenthosTend(benthosIndices%o2_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%o2_ind)
            netBenthosTend(benthosIndices%nh4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%nh4_ind)
            netBenthosTend(benthosIndices%h2po4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%h2po4_ind)
            netBenthosTend(benthosIndices%co2_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%co2_ind)
            netBenthosTend(benthosIndices%no3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%no3_ind)
            netBenthosTend(benthosIndices%mno2a_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%mno2a_ind)
            netBenthosTend(benthosIndices%mno2b_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%mno2b_ind)
            netBenthosTend(benthosIndices%mn_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%mn_ind)
            netBenthosTend(benthosIndices%feoh3a_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%feoh3a_ind)
            netBenthosTend(benthosIndices%feoh3b_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%feoh3b_ind)
            netBenthosTend(benthosIndices%fe_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%fe_ind)
            netBenthosTend(benthosIndices%fepa_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%fepa_ind)
            netBenthosTend(benthosIndices%fepb_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%fepb_ind)
            netBenthosTend(benthosIndices%so4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%so4_ind)
            netBenthosTend(benthosIndices%h2s_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%h2s_ind)
            netBenthosTend(benthosIndices%ch4_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%ch4_ind)
            netBenthosTend(benthosIndices%hco3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%hco3_ind)
            netBenthosTend(benthosIndices%fes_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%fes_ind)
            netBenthosTend(benthosIndices%fes2_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%fes2_ind)
            netBenthosTend(benthosIndices%s_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%s_ind)
            netBenthosTend(benthosIndices%caco3a_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%caco3a_ind)
            netBenthosTend(benthosIndices%caco3b_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%caco3b_ind)
            netBenthosTend(benthosIndices%co3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%co3_ind)
            netBenthosTend(benthosIndices%camgco3_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%camgco3_ind)
            netBenthosTend(benthosIndices%dic_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%dic_ind)
            netBenthosTend(benthosIndices%alk_ind,iCell) = &
               benthos_diagnostic_fields % diag_netBenthosTend(column,benthos_indices%alk_ind)

            netElements(elementIndices%carbon_ind,iCell) = &
               benthos_diagnostic_fields % diag_netElements(column,benthos_element_indices%carbon_ind)
            netElements(elementIndices%nitrogen_ind,iCell) = &
               benthos_diagnostic_fields % diag_netElements(column,benthos_element_indices%nitrogen_ind)
            netElements(elementIndices%oxygen_ind,iCell) = &
               benthos_diagnostic_fields % diag_netElements(column,benthos_element_indices%oxygen_ind)
            netElements(elementIndices%phosphorus_ind,iCell) = &
               benthos_diagnostic_fields % diag_netElements(column,benthos_element_indices%phosphorus_ind)
            netElements(elementIndices%sulfur_ind,iCell) = &
               benthos_diagnostic_fields % diag_netElements(column,benthos_element_indices%sulfur_ind)
            netElements(elementIndices%manganese_ind,iCell) = &
               benthos_diagnostic_fields % diag_netElements(column,benthos_element_indices%manganese_ind)
            netElements(elementIndices%iron_ind,iCell) = &
               benthos_diagnostic_fields % diag_netElements(column,benthos_element_indices%iron_ind)

         endif  !  config_benthosTracers_diagnostic_fields_level1

         !
         ! level 2 diagnostics
         !

         if (config_benthosTracers_diagnostic_fields_level2) then

            netBenthosRelErrorTransport(benthosIndices%poca_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%poca_ind)
            netBenthosRelErrorTransport(benthosIndices%pocb_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%pocb_ind)
            netBenthosRelErrorTransport(benthosIndices%pocc_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%pocc_ind)
            netBenthosRelErrorTransport(benthosIndices%pona_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%pona_ind)
            netBenthosRelErrorTransport(benthosIndices%ponb_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%ponb_ind)
            netBenthosRelErrorTransport(benthosIndices%ponc_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%ponc_ind)
            netBenthosRelErrorTransport(benthosIndices%popa_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%popa_ind)
            netBenthosRelErrorTransport(benthosIndices%popb_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%popb_ind)
            netBenthosRelErrorTransport(benthosIndices%popc_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%popc_ind)
            netBenthosRelErrorTransport(benthosIndices%o2_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%o2_ind)
            netBenthosRelErrorTransport(benthosIndices%nh4_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%nh4_ind)
            netBenthosRelErrorTransport(benthosIndices%h2po4_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%h2po4_ind)
            netBenthosRelErrorTransport(benthosIndices%co2_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%co2_ind)
            netBenthosRelErrorTransport(benthosIndices%no3_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%no3_ind)
            netBenthosRelErrorTransport(benthosIndices%mno2a_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%mno2a_ind)
            netBenthosRelErrorTransport(benthosIndices%mno2b_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%mno2b_ind)
            netBenthosRelErrorTransport(benthosIndices%mn_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%mn_ind)
            netBenthosRelErrorTransport(benthosIndices%feoh3a_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%feoh3a_ind)
            netBenthosRelErrorTransport(benthosIndices%feoh3b_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%feoh3b_ind)
            netBenthosRelErrorTransport(benthosIndices%fe_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%fe_ind)
            netBenthosRelErrorTransport(benthosIndices%fepa_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%fepa_ind)
            netBenthosRelErrorTransport(benthosIndices%fepb_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%fepb_ind)
            netBenthosRelErrorTransport(benthosIndices%so4_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%so4_ind)
            netBenthosRelErrorTransport(benthosIndices%h2s_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%h2s_ind)
            netBenthosRelErrorTransport(benthosIndices%ch4_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%ch4_ind)
            netBenthosRelErrorTransport(benthosIndices%hco3_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%hco3_ind)
            netBenthosRelErrorTransport(benthosIndices%fes_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%fes_ind)
            netBenthosRelErrorTransport(benthosIndices%fes2_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%fes2_ind)
            netBenthosRelErrorTransport(benthosIndices%s_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%s_ind)
            netBenthosRelErrorTransport(benthosIndices%caco3a_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%caco3a_ind)
            netBenthosRelErrorTransport(benthosIndices%caco3b_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%caco3b_ind)
            netBenthosRelErrorTransport(benthosIndices%co3_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%co3_ind)
            netBenthosRelErrorTransport(benthosIndices%camgco3_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%camgco3_ind)
            netBenthosRelErrorTransport(benthosIndices%dic_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%dic_ind)
            netBenthosRelErrorTransport(benthosIndices%alk_ind,iCell) = &
               benthos_diagnostic_fields % relErrorTransport(column,benthos_indices%alk_ind)

            netBenthosElementsRelError(elementIndices%carbon_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElements(column,benthos_element_indices%carbon_ind)
            netBenthosElementsRelError(elementIndices%nitrogen_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElements(column,benthos_element_indices%nitrogen_ind)
            netBenthosElementsRelError(elementIndices%oxygen_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElements(column,benthos_element_indices%oxygen_ind)
            netBenthosElementsRelError(elementIndices%phosphorus_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElements(column,benthos_element_indices%phosphorus_ind)
            netBenthosElementsRelError(elementIndices%sulfur_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElements(column,benthos_element_indices%sulfur_ind)
            netBenthosElementsRelError(elementIndices%manganese_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElements(column,benthos_element_indices%manganese_ind)
            netBenthosElementsRelError(elementIndices%iron_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElements(column,benthos_element_indices%iron_ind)

            netBenthosElementsRelErrorReactions(elementIndices%carbon_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElementsReactions(column,benthos_element_indices%carbon_ind)
            netBenthosElementsRelErrorReactions(elementIndices%nitrogen_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElementsReactions(column,benthos_element_indices%nitrogen_ind)
            netBenthosElementsRelErrorReactions(elementIndices%oxygen_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElementsReactions(column,benthos_element_indices%oxygen_ind)
            netBenthosElementsRelErrorReactions(elementIndices%phosphorus_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElementsReactions(column,benthos_element_indices%phosphorus_ind)
            netBenthosElementsRelErrorReactions(elementIndices%sulfur_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElementsReactions(column,benthos_element_indices%sulfur_ind)
            netBenthosElementsRelErrorReactions(elementIndices%manganese_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElementsReactions(column,benthos_element_indices%manganese_ind)
            netBenthosElementsRelErrorReactions(elementIndices%iron_ind,iCell) = &
               benthos_diagnostic_fields % relErrorElementsReactions(column,benthos_element_indices%iron_ind)

         endif

         !
         ! level 1 diagnostics with no depth dependence
         !
 !        if (config_ecosysTracers_diagnostic_fields_level1) then
 !           ecosys_diag_photoC_TOT_zint(iCell) =  &
 !                       BGC_diagnostic_fields%diag_photoC_TOT_zint(column)*convertLengthCGStoMKS
 !           ecosys_diag_photoC_NO3_TOT_zint(iCell) =  &
 !                       BGC_diagnostic_fields%diag_photoC_NO3_TOT_zint(column)*convertLengthCGStoMKS
 !           ecosys_diag_O2_ZMIN(iCell) =  &
 !                       BGC_diagnostic_fields%diag_O2_ZMIN(column)
!         endif  !  config_ecosysTracers_diagnostic_fields_level1

         !
         ! level 3 diagnostics with no depth dependence
         !

       end do ! iCell

      ! deallocate fields

      deallocate ( benthos_input%benthosTracerBulk)
      deallocate ( benthos_input%deepStorage)
      deallocate ( benthos_input%oceanBottomDepth)
      deallocate ( benthos_input%oceanBottomDensity)
      deallocate ( benthos_input%oceanBottomPressure)
!      deallocate ( benthos_input%tracerType(benthos_tracer_cnt) )
      deallocate ( benthos_input%PH_PREV_3D)
      deallocate ( benthos_input%bottom_level)
      deallocate ( benthos_input%oceanBottomTemperature)
      deallocate ( benthos_input%oceanBottomSalinity)
      deallocate ( benthos_input%oceanBottomPhosphate)
      deallocate ( benthos_input%oceanBottomAmmonium)
      deallocate ( benthos_input%oceanBottomNitrate)
      deallocate ( benthos_input%oceanBottomOxygen)
      deallocate ( benthos_input%oceanBottomIron)
      deallocate ( benthos_input%oceanBottomDIC)
      deallocate ( benthos_input%oceanBottomManganese)
      deallocate ( benthos_input%oceanBottomAlkalinity)
      deallocate ( benthos_input%oceanBottomSilicate)
      deallocate ( benthos_input%oceanPOCFlux)
      deallocate ( benthos_input%oceanPOPFlux)
      deallocate ( benthos_input%oceanPONFlux)
      deallocate ( benthos_input%oceanParticulateIronFlux)
      deallocate ( benthos_input%oceanCalciteFlux)
      deallocate ( benthos_input%oceanBiogenicSilicaFlux)
      deallocate ( benthos_input%oceanDustFlux)
      deallocate ( benthos_input%lithogenicFlux)
      deallocate ( benthos_input%terraPOCFlux)

      deallocate ( benthos_output%benthosTendencies)
      if (config_use_couple_ocean_ecosys) then
         deallocate ( benthos_output%benthosEcosysTendencies)
      end if
      deallocate ( benthos_output%benthosTransportTendencies )
      deallocate ( benthos_output%benthosReactionTendencies )
      deallocate ( benthos_output%elementTendencies )
      deallocate ( benthos_output%elementReactionTendencies )
      deallocate ( benthos_output%elementTransportTendencies )
      deallocate ( benthos_output%elementDeepStorage)

      deallocate ( benthos_output%PH_PREV_3D)
      deallocate ( benthos_output%deepStorage)

      deallocate ( benthos_output%oceanBottomTemperature)
      deallocate ( benthos_output%oceanBottomSalinity)
      deallocate ( benthos_output%oceanBottomSilicate)
      deallocate ( benthos_output%oceanBottomDensity)
      deallocate ( benthos_output%oceanBottomPressure)
      deallocate ( benthos_output%oceanBottomDepth)
      deallocate ( benthos_output%oceanBottomPhosphate)
      deallocate ( benthos_output%oceanBottomAmmonium)
      deallocate ( benthos_output%oceanBottomNitrate)
      deallocate ( benthos_output%oceanBottomOxygen)
      deallocate ( benthos_output%oceanBottomIron)
      deallocate ( benthos_output%oceanBottomDIC)
      deallocate ( benthos_output%oceanBottomAlkalinity)
      deallocate ( benthos_output%oceanBottomCarbonDioxide)
      deallocate ( benthos_output%oceanBottomHydrogenSulfide)
      deallocate ( benthos_output%oceanBottomMethane)
      deallocate ( benthos_output%oceanBottomBicarbonate)
      deallocate ( benthos_output%oceanBottomCarbonate)
      deallocate ( benthos_output%oceanBottomSulfate)
      deallocate ( benthos_output%oceanBottomManganese)
      deallocate ( benthos_output%oceanPOCaFlux)
      deallocate ( benthos_output%oceanPOPaFlux)
      deallocate ( benthos_output%oceanPONaFlux)
      deallocate ( benthos_output%oceanPOCbFlux)
      deallocate ( benthos_output%oceanPOPbFlux)
      deallocate ( benthos_output%oceanPONbFlux)
      deallocate ( benthos_output%oceanPOCcFlux)
      deallocate ( benthos_output%oceanPOPcFlux)
      deallocate ( benthos_output%oceanPONcFlux)
      deallocate ( benthos_output%oceanCalciteFlux)
      deallocate ( benthos_output%oceanAragoniteFlux)
      deallocate ( benthos_output%oceanMgCalciteFlux)
      deallocate ( benthos_output%oceanManganeseOxideaFlux)
      deallocate ( benthos_output%oceanManganeseOxidebFlux)
      deallocate ( benthos_output%oceanIronOxyhydroxideaFlux)
      deallocate ( benthos_output%oceanIronOxyhydroxidebFlux)
      deallocate ( benthos_output%oceanIronBoundPhosphorusaFlux)
      deallocate ( benthos_output%oceanIronBoundPhosphorusbFlux)
      deallocate ( benthos_output%oceanPrecipManganeseOxidea)
      deallocate ( benthos_output%oceanPrecipManganeseOxideb)
      deallocate ( benthos_output%oceanPrecipIronOxyhydroxidea)
      deallocate ( benthos_output%oceanPrecipIronOxyhydroxideb)
      deallocate ( benthos_output%oceanPrecipIronBoundPhosphorusa)
      deallocate ( benthos_output%oceanPrecipIronBoundPhosphorusb)
      deallocate ( benthos_output%oceanTotalSedimentMassFlux)

      deallocate ( benthos_forcing%deepStoragePrescribedFlux)
      deallocate ( benthos_forcing%benthosOceanSedimentFlux)

        !---------------------------------------------------------------------------
        !   deallocate flux diagnostic output fields
        !   all fluxes are positive down
        !---------------------------------------------------------------------------

      deallocate (benthos_flux_diagnostic_fields%deepStorageDiffVelFlux )
      deallocate (benthos_flux_diagnostic_fields%deepStorageFlux )
      deallocate (benthos_flux_diagnostic_fields%benthosOceanDiffVelFlux )
      deallocate (benthos_flux_diagnostic_fields%benthosOceanFlux )
      deallocate (benthos_flux_diagnostic_fields%elementSedimentation)
      deallocate (benthos_flux_diagnostic_fields%elementBurialExchange)
      deallocate (benthos_flux_diagnostic_fields%elementOceanExchange)

      !---------------------------------------------------------------------------
      !   deallocate diagnostic output fields
      !---------------------------------------------------------------------------

      ! 2D column integrated diagnostics

      deallocate (benthos_diagnostic_fields%relErrorTransport)
      deallocate (benthos_diagnostic_fields%relErrorElementsReactions)
      deallocate (benthos_diagnostic_fields%relErrorElements)
      deallocate (benthos_diagnostic_fields%diag_primarySourceTend )
      deallocate (benthos_diagnostic_fields%diag_primarySinkTend )
      deallocate (benthos_diagnostic_fields%diag_secondarySourceTend )
      deallocate (benthos_diagnostic_fields%diag_secondarySinkTend )
      deallocate (benthos_diagnostic_fields%diag_carbonateSourceTend )
      deallocate (benthos_diagnostic_fields%diag_carbonateSinkTend )

      deallocate (benthos_diagnostic_fields%diag_netElements )
      deallocate (benthos_diagnostic_fields%diag_netBenthosTracersCell )
      deallocate (benthos_diagnostic_fields%diag_netBenthosReactionTend )
      deallocate (benthos_diagnostic_fields%diag_netBenthosTransportTend )
      deallocate (benthos_diagnostic_fields%diag_netBenthosTend )
      deallocate (benthos_diagnostic_fields%diag_netElementsReactionTend )
      deallocate (benthos_diagnostic_fields%diag_netElementsTransportTend )
      deallocate (benthos_diagnostic_fields%diag_netElementsTend )

      call mpas_timer_stop("benthos source-sink")

   !--------------------------------------------------------------------

   end subroutine benthos_tracer_compute!}}}

!***********************************************************************
!
!  routine benthos_tracer_init
!
!> \brief   Initializes
!> \author  Nicole Jeffery
!> \date    07/16/2020
!> \details
!>  This routine initializes indices for benthos tracers
!
!-----------------------------------------------------------------------

   subroutine ocn_benthos_tracer_init(domain,err)!{{{


!NOTE:  called from mpas_ocn_forward_mode.F

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: seafloorTracersPool
      type (mpas_pool_type), pointer :: seafloorStoragePool
      type (mpas_pool_type), pointer :: seafloorTendPool
      type (mpas_pool_type), pointer :: benthosAuxiliaryPool
      type (mpas_pool_type), pointer :: benthosDiagFieldsLevel1Pool
      type (mpas_pool_type), pointer :: benthosMonthlyForcing
      type (mpas_pool_type), pointer :: benthosMonthlyTerraForcing

      real (kind=RKIND), dimension(:,:,:), pointer :: &
          benthosTracers

      ! seafloorStorage

      real (kind=RKIND), dimension(:,:), pointer:: &
          benthosTracersStorage, &
          benthosPhPrev3d

      real (kind=RKIND), dimension(:), pointer :: &
          benthosInterfaceLayerThickness, &
          benthosVertGridPoints, &
          bottomDepth, &
          latCell

      integer, dimension(:), pointer :: &
          indexToBenthosCells

      ! scalar pointers
      integer, pointer :: &
          nBenthicVertLevels,       &
          index_dummy,              &
          config_set_benthos_flags, &
          nCellsSolve

      character(len=strKIND) :: &
         forcingIntervalMonthly,  &
         forcingReferenceTimeMonthly

       logical, pointer :: &
          !config_do_restart, &
          !config_do_restart_benthos, &
          config_use_carbonate_saturation, &
          config_use_secondary_reactions, &
          config_use_depth_dependent_porosity, &
          config_use_non_zero_diffusivity, &
          config_use_sedimentation, &
          config_use_fast_sedimentation, &
          config_use_positive_sedimentation, &
          config_use_bgc_sinking_flux, &
          config_use_step_initial_profiles, &
          config_use_benthic_reactions, &
          config_use_flux_correction, &
          config_use_ocean_conc, &
          config_use_deep_source, &
          config_use_constant_diffusivity, &
          config_use_ecosysTracers, &
          config_use_ocean_bottom_state, &
          config_use_benthos_monthly_forcing, &
          config_use_benthos_monthly_terra_forcing, &
          config_use_couple_ocean_ecosys

        real (kind=RKIND), dimension(:), pointer :: &
           oceanBenthosPOCFluxData, &
           oceanBenthosPICFluxData, &
           oceanBenthosBiogenicSilicaFluxData, &
           oceanBenthosDustFluxData, &
           oceanBenthosMno2FluxData, &
           oceanBenthosFeoh3FluxData, &
           oceanBenthosPrecipMno2Data, &
           oceanBenthosPrecipFeoh3Data, &
           oceanBottomBenthosSio3Data, &
           oceanBottomBenthosPo4Data, &
           oceanBottomBenthosNh4Data, &
           oceanBottomBenthosNo3Data, &
           oceanBottomBenthosO2Data, &
           oceanBottomBenthosFeData, &
           oceanBottomBenthosMnData, &
           oceanBottomBenthosSo4Data, &
           oceanBottomBenthosDICData, &
           oceanBottomBenthosAlkData, &
           oceanBottomBenthosTemperatureData, &
           oceanBottomBenthosSalinityData, &
           oceanBottomBenthosPressureData, &
           oceanBottomBenthosDensityData, &
           oceanTerraBenthosPOCFluxData, &
           oceanTerraLithogenicFluxData

      ! local variables

      integer :: &
         nBenthicTracers,      &
         nElements,            &
         numBenthicColumnsMax, &
         iBenthicTracers,      &
         iVertLevels,          &
         iCell, &
         iCellBenthic

      logical :: &
          test_flag

      !
      ! get tracers pools
      !

      err = 0

      !
      ! Get tracer group so we can get the number of tracers in it
      !

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_subpool(statePool, 'seafloorTracers', seafloorTracersPool)
      call mpas_pool_get_subpool(statePool, 'seafloorStorage', seafloorStoragePool)
      call mpas_pool_get_subpool(tendPool, 'seafloorTend', seafloorTendPool)
      call mpas_pool_get_subpool(diagnosticsPool, 'benthosDiagFieldsLevel1', benthosDiagFieldsLevel1Pool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nBenthicVertLevels', nBenthicVertLevels)

      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(seafloorTracersPool, 'benthosTracers', benthosTracers)
      call mpas_pool_get_array(seafloorTracersPool, 'benthosPhPrev3d', benthosPhPrev3d)
      call mpas_pool_get_array(seafloorStoragePool, 'benthosTracersStorage', benthosTracersStorage)
      call mpas_pool_get_array(seafloorTendPool, 'benthosInterfaceLayerThickness', benthosInterfaceLayerThickness)

!      call mpas_log_write('Error test warning, mpas_ocn_tracer_benthos.F',MPAS_LOG_WARN)
!      call mpas_log_write('Check benthos_tracer_cnt: $i',MPAS_LOG_ERR,intArgs=(/benthos_tracer_cnt/))

      if (associated(benthosTracers)) then   ! associated for benthosGRP

      nBenthicTracers = size(benthosTracers, dim=1)
!      call mpas_log_write('test associated : nBenthicTracers = $i',MPAS_LOG_ERR,intArgs=(/nBenthicTracers/))
      if (benthos_tracer_cnt /= nBenthicTracers) then
         call mpas_log_write('benthos_tracer_cnt /= nBenthicTracers: $i',MPAS_LOG_CRIT,intArgs=(/nBenthicTracers/))
         err = 1
         return
      end if

!-----------------------------------------------------------------------
!  initialize benthos parameters
!-----------------------------------------------------------------------

      allocate( benthos_indices%short_name(benthos_tracer_cnt) )
      allocate( benthos_indices%long_name(benthos_tracer_cnt) )
      allocate( benthos_indices%units(benthos_tracer_cnt) )

! no need to allocate the above fields for ecosysBenthosIndices (only use the index)

!-----------------------------------------------------------------------
!  sets most of benthos parameters
!  sets namelist defaults
!-----------------------------------------------------------------------
 
      call benthos_parms_init   ! initialize namelist defaults
      !
      ! for now only do 1 column at a time
      !
      T0_Kelvin_benthos = T0_Kelvin

      numBenthicColumnsMax = 1

      benthos_indices%poca_ind    =  1
      benthos_indices%pocb_ind    =  2
      benthos_indices%pocc_ind    =  3
      benthos_indices%pona_ind    =  4
      benthos_indices%ponb_ind    =  5
      benthos_indices%ponc_ind    =  6
      benthos_indices%popa_ind    =  7
      benthos_indices%popb_ind    =  8
      benthos_indices%popc_ind    =  9
      benthos_indices%o2_ind      = 10
      benthos_indices%nh4_ind     = 11
      benthos_indices%h2po4_ind   = 12
      benthos_indices%co2_ind     = 13
      benthos_indices%no3_ind     = 14
      benthos_indices%mno2a_ind   = 15
      benthos_indices%mno2b_ind   = 16
      benthos_indices%mn_ind      = 17
      benthos_indices%feoh3a_ind  = 18
      benthos_indices%feoh3b_ind  = 19
      benthos_indices%fe_ind      = 20
      benthos_indices%fepa_ind    = 21
      benthos_indices%fepb_ind    = 22
      benthos_indices%so4_ind     = 23
      benthos_indices%h2s_ind     = 24
      benthos_indices%ch4_ind     = 25
      benthos_indices%dic_ind     = 26
      benthos_indices%alk_ind     = 27
      benthos_indices%hco3_ind    = 28
      benthos_indices%fes_ind     = 29
      benthos_indices%fes2_ind    = 30
      benthos_indices%s_ind       = 31
      benthos_indices%caco3a_ind  = 32
      benthos_indices%caco3b_ind  = 33
      benthos_indices%co3_ind     = 34
      benthos_indices%camgco3_ind = 35

      benthos_element_indices%carbon_ind     = 1
      benthos_element_indices%oxygen_ind     = 2
      benthos_element_indices%nitrogen_ind   = 3
      benthos_element_indices%phosphorus_ind = 4
      benthos_element_indices%sulfur_ind     = 5
      benthos_element_indices%manganese_ind  = 6
      benthos_element_indices%iron_ind       = 7

      nElements = 7

      allocate( benthos_element_indices%short_name(nElements) )
      allocate( benthos_element_indices%long_name(nElements) )
      allocate( benthos_element_indices%units(nElements) )

!      call mpas_log_write('test assign indices: benthos_indices%caco3a_ind = $i',MPAS_LOG_ERR,intArgs=(/benthos_indices%caco3a_ind/))
!      call mpas_log_write('test assign indices: benthos_element_indices%iron_ind = $i',MPAS_LOG_ERR,intArgs=(/benthos_element_indices%iron_ind/))
!      err = 1
!      return

      call mpas_pool_get_dimension(seafloorTracersPool, 'index_pocaBenthos',  index_dummy)
      benthosIndices%poca_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_pocbBenthos',  index_dummy)
      benthosIndices%pocb_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_poccBenthos',  index_dummy)
      benthosIndices%pocc_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_ponaBenthos',  index_dummy)
      benthosIndices%pona_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_ponbBenthos',  index_dummy)
      benthosIndices%ponb_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_poncBenthos',  index_dummy)
      benthosIndices%ponc_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_popaBenthos',  index_dummy)
      benthosIndices%popa_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_popbBenthos',  index_dummy)
      benthosIndices%popb_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_popcBenthos',  index_dummy)
      benthosIndices%popc_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_o2Benthos',  index_dummy)
      benthosIndices%o2_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_nh4Benthos',  index_dummy)
      benthosIndices%nh4_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_h2po4Benthos',  index_dummy)
      benthosIndices%h2po4_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_co2Benthos',  index_dummy)
      benthosIndices%co2_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_no3Benthos',  index_dummy)
      benthosIndices%no3_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_mno2aBenthos',  index_dummy)
      benthosIndices%mno2a_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_mno2bBenthos',  index_dummy)
      benthosIndices%mno2b_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_mnBenthos',  index_dummy)
      benthosIndices%mn_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_feoh3aBenthos',  index_dummy)
      benthosIndices%feoh3a_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_feoh3bBenthos',  index_dummy)
      benthosIndices%feoh3b_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_feBenthos',  index_dummy)
      benthosIndices%fe_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_fepaBenthos',  index_dummy)
      benthosIndices%fepa_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_fepbBenthos',  index_dummy)
      benthosIndices%fepb_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_so4Benthos',  index_dummy)
      benthosIndices%so4_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_h2sBenthos',  index_dummy)
      benthosIndices%h2s_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_ch4Benthos',  index_dummy)
      benthosIndices%ch4_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_dicBenthos',  index_dummy)
      benthosIndices%dic_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_alkBenthos',  index_dummy)
      benthosIndices%alk_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_hco3Benthos',  index_dummy)
      benthosIndices%hco3_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_fesBenthos',  index_dummy)
      benthosIndices%fes_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_fes2Benthos',  index_dummy)
      benthosIndices%fes2_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_sBenthos',  index_dummy)
      benthosIndices%s_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_caco3aBenthos',  index_dummy)
      benthosIndices%caco3a_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_caco3bBenthos',  index_dummy)
      benthosIndices%caco3b_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_co3Benthos',  index_dummy)
      benthosIndices%co3_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_camgco3Benthos',  index_dummy)
      benthosIndices%camgco3_ind  = index_dummy

      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netCarbon',  index_dummy)

!      call mpas_log_write('test get dimension:index_netCarbon = $i',MPAS_LOG_ERR,intArgs=(/index_dummy/))

      elementIndices%carbon_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netOxygen',  index_dummy)
      elementIndices%oxygen_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netNitrogen',  index_dummy)
      elementIndices%nitrogen_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netPhosphorus',  index_dummy)
      elementIndices%phosphorus_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netSulfur',  index_dummy)
      elementIndices%sulfur_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netManganese',  index_dummy)
      elementIndices%manganese_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netIron',  index_dummy)
      elementIndices%iron_ind  = index_dummy

!      call mpas_log_write('test get dimension : elementIndices%iron_ind = $i',MPAS_LOG_ERR,intArgs=(/elementIndices%iron_ind/))

!      call mpas_pool_get_dimension(tracersPool, 'index_DMS',  index_dummy)
!      dmsIndices%dms_ind  = index_dummy
!      call mpas_pool_get_dimension(tracersPool, 'index_DMSP', index_dummy)
!      dmsIndices%dmsp_ind = index_dummy

      call mpas_pool_get_config(ocnConfigs, 'config_use_ecosysTracers',  config_use_ecosysTracers)
      if (config_use_ecosysTracers) then
         call mpas_pool_get_dimension(tracersPool, 'index_PO4',      index_dummy)
         ecosysBenthosIndices%po4_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_NO3',      index_dummy)
         ecosysBenthosIndices%no3_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_SiO3',     index_dummy)
         ecosysBenthosIndices%sio3_ind            = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_NH4',      index_dummy)
         ecosysBenthosIndices%nh4_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_Fe',       index_dummy)
         ecosysBenthosIndices%fe_ind              = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_Lig',       index_dummy)
         ecosysBenthosIndices%lig_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_O2',       index_dummy)
         ecosysBenthosIndices%o2_ind              = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_DIC',      index_dummy)
         ecosysBenthosIndices%dic_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_DIC_ALT_CO2', index_dummy)
         ecosysBenthosIndices%dic_alt_co2_ind     = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_ALK',      index_dummy)
         ecosysBenthosIndices%alk_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_ALK_ALT_CO2',      index_dummy)
         ecosysBenthosIndices%alk_alt_co2_ind     = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_DOC',      index_dummy)
         ecosysBenthosIndices%doc_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_DON',      index_dummy)
         ecosysBenthosIndices%don_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_DOP',      index_dummy)
         ecosysBenthosIndices%dop_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_DOPr',     index_dummy)
         ecosysBenthosIndices%dopr_ind            = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_DONr',     index_dummy)
         ecosysBenthosIndices%donr_ind            = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_DOCr',     index_dummy)
         ecosysBenthosIndices%docr_ind            = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_zooC',     index_dummy)
         ecosysBenthosIndices%zooC_ind            = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_spChl',    index_dummy)
         ecosysBenthosIndices%spChl_ind           = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_spC',      index_dummy)
         ecosysBenthosIndices%spC_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_spFe',     index_dummy)
         ecosysBenthosIndices%spFe_ind            = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_spP',      index_dummy)
         ecosysBenthosIndices%spP_ind             = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_spCaCO3',  index_dummy)
         ecosysBenthosIndices%spCaCO3_ind         = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_diatChl',  index_dummy)
         ecosysBenthosIndices%diatChl_ind         = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_diatC',    index_dummy)
         ecosysBenthosIndices%diatC_ind           = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_diatFe',   index_dummy)
         ecosysBenthosIndices%diatFe_ind          = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_diatP',    index_dummy)
         ecosysBenthosIndices%diatP_ind           = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_diatSi',   index_dummy)
         ecosysBenthosIndices%diatSi_ind          = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_diazChl',  index_dummy)
         ecosysBenthosIndices%diazChl_ind         = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_diazC',    index_dummy)
         ecosysBenthosIndices%diazC_ind           = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_diazFe',   index_dummy)
         ecosysBenthosIndices%diazFe_ind          = index_dummy
         call mpas_pool_get_dimension(tracersPool, 'index_diazP',    index_dummy)
         ecosysBenthosIndices%diazP_ind           = index_dummy
         !call mpas_pool_get_dimension(tracersPool, 'index_phaeoChl', index_dummy)
         !ecosysBenthosIndices%phaeoChl_ind        = index_dummy
         !call mpas_pool_get_dimension(tracersPool, 'index_phaeoC',   index_dummy)
         !ecosysBenthosIndices%phaeoC_ind          = index_dummy
         !call mpas_pool_get_dimension(tracersPool, 'index_phaeoFe',  index_dummy)
         !ecosysBenthosIndices%phaeoFe_ind         = index_dummy
      else
         index_dummy = 1
         ecosysBenthosIndices%po4_ind             = index_dummy
         ecosysBenthosIndices%no3_ind             = index_dummy
         ecosysBenthosIndices%sio3_ind            = index_dummy
         ecosysBenthosIndices%nh4_ind             = index_dummy
         ecosysBenthosIndices%fe_ind              = index_dummy
         ecosysBenthosIndices%lig_ind             = index_dummy
         ecosysBenthosIndices%o2_ind              = index_dummy
         ecosysBenthosIndices%dic_ind             = index_dummy
         ecosysBenthosIndices%dic_alt_co2_ind     = index_dummy
         ecosysBenthosIndices%alk_ind             = index_dummy
         ecosysBenthosIndices%alk_alt_co2_ind     = index_dummy
         ecosysBenthosIndices%doc_ind             = index_dummy
         ecosysBenthosIndices%don_ind             = index_dummy
         ecosysBenthosIndices%dop_ind             = index_dummy
         ecosysBenthosIndices%dopr_ind            = index_dummy
         ecosysBenthosIndices%donr_ind            = index_dummy
         ecosysBenthosIndices%docr_ind            = index_dummy
         ecosysBenthosIndices%zooC_ind            = index_dummy
         ecosysBenthosIndices%spChl_ind           = index_dummy
         ecosysBenthosIndices%spC_ind             = index_dummy
         ecosysBenthosIndices%spFe_ind            = index_dummy
         ecosysBenthosIndices%spP_ind             = index_dummy
         ecosysBenthosIndices%spCaCO3_ind         = index_dummy
         ecosysBenthosIndices%diatChl_ind         = index_dummy
         ecosysBenthosIndices%diatC_ind           = index_dummy
         ecosysBenthosIndices%diatFe_ind          = index_dummy
         ecosysBenthosIndices%diatP_ind           = index_dummy
         ecosysBenthosIndices%diatSi_ind          = index_dummy
         ecosysBenthosIndices%diazChl_ind         = index_dummy
         ecosysBenthosIndices%diazC_ind           = index_dummy
         ecosysBenthosIndices%diazFe_ind          = index_dummy
         ecosysBenthosIndices%diazP_ind           = index_dummy
         !ecosysBenthosIndices%phaeoChl_ind        = index_dummy
         !ecosysBenthosIndices%phaeoC_ind          = index_dummy
         !ecosysBenthosIndices%phaeoFe_ind         = index_dummy
      end if

!      call mpas_log_write('test get dimension : ecosysBenthosIndices%phaeoFe_ind = $i',MPAS_LOG_ERR,intArgs=(/ecosysBenthosIndices%phaeoFe_ind/))

      call mpas_pool_get_config(ocnConfigs, 'config_set_benthos_flags',  config_set_benthos_flags)
      call mpas_pool_get_config(ocnConfigs, 'config_use_carbonate_saturation',  config_use_carbonate_saturation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_secondary_reactions',  config_use_secondary_reactions)
      call mpas_pool_get_config(ocnConfigs, 'config_use_depth_dependent_porosity',  config_use_depth_dependent_porosity)
      call mpas_pool_get_config(ocnConfigs, 'config_use_non_zero_diffusivity',  config_use_non_zero_diffusivity)
      call mpas_pool_get_config(ocnConfigs, 'config_use_sedimentation',  config_use_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_fast_sedimentation',  config_use_fast_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_positive_sedimentation',  config_use_positive_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_bgc_sinking_flux',  config_use_bgc_sinking_flux)
      call mpas_pool_get_config(ocnConfigs, 'config_use_step_initial_profiles',  config_use_step_initial_profiles)
      call mpas_pool_get_config(ocnConfigs, 'config_use_benthic_reactions',  config_use_benthic_reactions)
      call mpas_pool_get_config(ocnConfigs, 'config_use_flux_correction',  config_use_flux_correction)
      call mpas_pool_get_config(ocnConfigs, 'config_use_ocean_conc',  config_use_ocean_conc)
      call mpas_pool_get_config(ocnConfigs, 'config_use_deep_source',  config_use_deep_source)
      call mpas_pool_get_config(ocnConfigs, 'config_use_constant_diffusivity',  config_use_constant_diffusivity)
      call mpas_pool_get_config(ocnConfigs, 'config_use_benthos_monthly_forcing',  config_use_benthos_monthly_forcing)
      call mpas_pool_get_config(ocnConfigs, 'config_use_benthos_monthly_terra_forcing',  config_use_benthos_monthly_terra_forcing)
      call mpas_pool_get_config(ocnConfigs, 'config_use_ocean_bottom_state',  config_use_ocean_bottom_state)
      call mpas_pool_get_config(ocnConfigs, 'config_use_couple_ocean_ecosys',  config_use_couple_ocean_ecosys)

      call mpas_pool_get_subpool(forcingPool, 'benthosAuxiliary', benthosAuxiliaryPool)
      call mpas_pool_get_array(benthosAuxiliaryPool,'benthosVertGridPoints', benthosVertGridPoints)
      call mpas_pool_get_array(benthosAuxiliaryPool,'indexToBenthosCells', indexToBenthosCells)

    ! Set configuration flags

      call benthos_flags_init(config_set_benthos_flags, config_use_carbonate_saturation, &
           config_use_secondary_reactions, config_use_depth_dependent_porosity, &
           config_use_non_zero_diffusivity, config_use_sedimentation, &
           config_use_fast_sedimentation, config_use_positive_sedimentation, &
           config_use_bgc_sinking_flux, config_use_step_initial_profiles, &
           config_use_benthic_reactions, config_use_flux_correction, config_use_ocean_conc, &
           config_use_deep_source, config_use_constant_diffusivity, config_use_ocean_bottom_state,&
           config_use_benthos_monthly_forcing, config_use_couple_ocean_ecosys)

!NJ-TEST: VERIFIED
test_flag = .false.
      if (test_flag) then
         call mpas_log_write('test benthos_flags_init')
         call mpas_log_write('config_set_benthos_flags: $i',MPAS_LOG_ERR,intArgs=(/config_set_benthos_flags/))
         if (config_use_carbonate_saturation) then
            call mpas_log_write('config_use_carbonate_saturation = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_carbonate_saturation = false', MPAS_LOG_ERR)
         end if
         if (config_use_secondary_reactions) then
            call mpas_log_write('config_use_secondary_reactions = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_secondary_reactions = false', MPAS_LOG_ERR)
         end if
         if (config_use_depth_dependent_porosity) then
            call mpas_log_write('config_use_depth_dependent_porosity = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_depth_dependent_porosity = false', MPAS_LOG_ERR)
         end if
         if (config_use_non_zero_diffusivity) then
            call mpas_log_write('config_use_non_zero_diffusivity = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_non_zero_diffusivity = false', MPAS_LOG_ERR)
         end if
         if (config_use_sedimentation) then
            call mpas_log_write('config_use_sedimentation = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_sedimentation = false', MPAS_LOG_ERR)
         end if
         if (config_use_fast_sedimentation) then
            call mpas_log_write('config_use_fast_sedimentation = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_fast_sedimentation = false', MPAS_LOG_ERR)
         end if
         if (config_use_positive_sedimentation) then
            call mpas_log_write('config_use_positive_sedimentation = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_positive_sedimentation = false', MPAS_LOG_ERR)
         end if
         if (config_use_bgc_sinking_flux) then
            call mpas_log_write('config_use_bgc_sinking_flux = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_bgc_sinking_flux = false', MPAS_LOG_ERR)
         end if
         if (config_use_step_initial_profiles) then
            call mpas_log_write('config_use_step_initial_profiles = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_step_initial_profiles = false', MPAS_LOG_ERR)
         end if
         if (config_use_benthic_reactions) then
            call mpas_log_write('config_use_benthic_reactions = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_benthic_reactions = false', MPAS_LOG_ERR)
         end if
         if (config_use_flux_correction) then
            call mpas_log_write('config_use_flux_correction = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_flux_correction = false', MPAS_LOG_ERR)
         end if
         if (config_use_ocean_conc) then
            call mpas_log_write('config_use_ocean_conc = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_ocean_conc = false', MPAS_LOG_ERR)
         end if
         if (config_use_deep_source) then
            call mpas_log_write('config_use_deep_source = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_deep_source = false', MPAS_LOG_ERR)
         end if
         if (config_use_constant_diffusivity) then
            call mpas_log_write('config_use_constant_diffusivity = true', MPAS_LOG_ERR)
         else
            call mpas_log_write('config_use_constant_diffusivity = false', MPAS_LOG_ERR)
         end if
!         err = 1
!         return
     end if ! test_benthos_flags

     ! benthos_init sets short and long names, units in benthos_indices
     ! also sets benthos model parameter fields

      call benthos_init(benthos_indices, benthos_element_indices, nBenthicVertLevels, nBenthicTracers, &
         nElements, benthosInterfaceLayerThickness, benthosVertGridPoints)

!      call mpas_log_write('after benthos_init', MPAS_LOG_ERR)
!     NJ-TEST: VERIFIED
      test_flag = .false.
      if (test_flag) then
         call mpas_log_write('test benthos_init', MPAS_LOG_ERR)
         call mpas_log_write('nBenthicTracers= $i',MPAS_LOG_ERR,intArgs=(/nBenthicTracers/))
         call mpas_log_write('nElements= $i',MPAS_LOG_ERR,intArgs=(/nElements/))
         call mpas_log_write('benthos_indices%units(benthos_indices%poca_ind)='//trim(benthos_indices%units(benthos_indices%poca_ind)))
         call mpas_log_write('benthos_indices%short_name(benthos_indices%poca_ind)='//trim(benthos_indices%short_name(benthos_indices%poca_ind)))
         call mpas_log_write('benthos_indices%long_name(benthos_indices%poca_ind)='//trim(benthos_indices%long_name(benthos_indices%poca_ind)))
         call mpas_log_write('benthos_element_indices%units(benthos_element_indices%iron_ind)='//trim(benthos_element_indices%units(benthos_element_indices%iron_ind)))
         call mpas_log_write('benthos_element_indices%short_name(benthos_element_indices%iron_ind)='//trim(benthos_element_indices%short_name(benthos_element_indices%iron_ind)))
         call mpas_log_write('benthos_element_indices%long_name(benthos_element_indices%iron_ind)='//trim(benthos_element_indices%long_name(benthos_element_indices%iron_ind)))
!         err = 1
!         return
      end if

!NOTES:

!also check short_name with mpas variable name

      allocate ( benthos_input%tracerType(benthos_tracer_cnt) )

      !call MPAS_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
      !call MPAS_pool_get_config(domain % configs, 'config_do_restart_benthos', config_do_restart_benthos)

      if (config_use_benthos_monthly_forcing) then
         ! initialize monthly forcing to be read from file (sediment fluxes)

         forcingIntervalMonthly = "0000-01-00_00:00:00"
         forcingReferenceTimeMonthly = "0000-01-15_00:00:00"

         call MPAS_forcing_init_group( forcingGroupHead,  &
                "benthosMonthlyClimatology", &
                domain, &
                '0000-01-01_00:00:00', &
                '0000-01-01_00:00:00', &
                '0001-00-00_00:00:00', &
                config_do_restart)

         call mpas_pool_get_subpool(domain % blocklist % structs, 'benthosMonthlyForcing', benthosMonthlyForcing)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosPOCFluxData', oceanBenthosPOCFluxData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosPICFluxData', oceanBenthosPICFluxData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosBiogenicSilicaFluxData', oceanBenthosBiogenicSilicaFluxData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosDustFluxData', oceanBenthosDustFluxData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosMno2FluxData', oceanBenthosMno2FluxData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosFeoh3FluxData', oceanBenthosFeoh3FluxData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosPrecipMno2Data', oceanBenthosPrecipMno2Data)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosPrecipFeoh3Data', oceanBenthosPrecipFeoh3Data)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosSio3Data', oceanBottomBenthosSio3Data)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosPo4Data', oceanBottomBenthosPo4Data)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosNh4Data', oceanBottomBenthosNh4Data)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosNo3Data', oceanBottomBenthosNo3Data)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosO2Data', oceanBottomBenthosO2Data)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosFeData', oceanBottomBenthosFeData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosMnData', oceanBottomBenthosMnData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosSo4Data', oceanBottomBenthosSo4Data)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosDICData', oceanBottomBenthosDICData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosAlkData', oceanBottomBenthosAlkData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosTemperatureData', oceanBottomBenthosTemperatureData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosSalinityData', oceanBottomBenthosSalinityData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosDensityData', oceanBottomBenthosDensityData)
         call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosPressureData', oceanBottomBenthosPressureData)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBenthosPOCFluxData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBenthosPOCFluxData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBenthosPICFluxData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBenthosPICFluxData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBenthosBiogenicSilicaFluxData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBenthosBiogenicSilicaFluxData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBenthosDustFluxData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBenthosDustFluxData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBenthosMno2FluxData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBenthosMno2FluxData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBenthosFeoh3FluxData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBenthosFeoh3FluxData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBenthosPrecipMno2Data', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBenthosPrecipMno2Data', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBenthosPrecipFeoh3Data', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBenthosPrecipFeoh3Data', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosSio3Data', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosSio3Data', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosPo4Data', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosPo4Data', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosNh4Data', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosNh4Data', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosNo3Data', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosNo3Data', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosO2Data', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosO2Data', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosFeData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosFeData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosMnData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosMnData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosSo4Data', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosSo4Data', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosDICData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosDICData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosAlkData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosAlkData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosTemperatureData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosTemperatureData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosSalinityData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosSalinityData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosDensityData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosDensityData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyClimatology', &
                'oceanBottomBenthosPressureData', &
                'benthos_monthly_flux', &
                'benthosMonthlyForcing',  &
                'oceanBottomBenthosPressureData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field_data( forcingGroupHead, &
           'benthosMonthlyClimatology',  &
            domain % streamManager,  &
            config_do_restart, &
            .false.)

     end if ! config_use_benthos_monthly_forcing

     if (config_use_benthos_monthly_terra_forcing) then
         ! initialize monthly forcing to be read from file (terrigenous sediment fluxes)

         forcingIntervalMonthly = "0000-01-00_00:00:00"
         forcingReferenceTimeMonthly = "0000-01-15_00:00:00"

         call MPAS_forcing_init_group( forcingGroupHead,  &
                "benthosMonthlyTerraClimatology", &
                domain, &
                '0000-01-01_00:00:00', &
                '0000-01-01_00:00:00', &
                '0001-00-00_00:00:00', &
                config_do_restart)

         call mpas_pool_get_subpool(domain % blocklist % structs, 'benthosMonthlyTerraForcing', benthosMonthlyTerraForcing)
         call mpas_pool_get_array(benthosMonthlyTerraForcing, 'oceanTerraBenthosPOCFluxData', oceanTerraBenthosPOCFluxData)
         call mpas_pool_get_array(benthosMonthlyTerraForcing, 'oceanTerraLithogenicFluxData', oceanTerraLithogenicFluxData)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyTerraClimatology', &
                'oceanTerraBenthosPOCFluxData', &
                'benthos_monthly_terra_flux', &
                'benthosMonthlyTerraForcing',  &
                'oceanTerraBenthosPOCFluxData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'benthosMonthlyTerraClimatology', &
                'oceanTerraLithogenicFluxData', &
                'benthos_monthly_terra_flux', &
                'benthosMonthlyTerraForcing',  &
                'oceanTerraLithogenicFluxData', &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

         call MPAS_forcing_init_field_data( forcingGroupHead, &
        'benthosMonthlyTerraClimatology',  &
         domain % streamManager,  &
         config_do_restart, &
         .false.)

     end if

   ! TracerType: 1 for solutes, 0 for solids

      do iBenthicTracers = 1, nBenthicTracers
         benthos_input%tracerType(iBenthicTracers) = 1
      end do

   ! solids
      benthos_input%tracerType(benthos_indices%poca_ind) = 0
      benthos_input%tracerType(benthos_indices%pocb_ind) = 0
      benthos_input%tracerType(benthos_indices%pocc_ind) = 0
      benthos_input%tracerType(benthos_indices%pona_ind) = 0
      benthos_input%tracerType(benthos_indices%ponb_ind) = 0
      benthos_input%tracerType(benthos_indices%ponc_ind) = 0
      benthos_input%tracerType(benthos_indices%popa_ind) = 0
      benthos_input%tracerType(benthos_indices%popb_ind) = 0
      benthos_input%tracerType(benthos_indices%popc_ind) = 0
      benthos_input%tracerType(benthos_indices%mno2a_ind) = 0
      benthos_input%tracerType(benthos_indices%mno2b_ind) = 0
      benthos_input%tracerType(benthos_indices%feoh3a_ind) = 0
      benthos_input%tracerType(benthos_indices%feoh3b_ind) = 0
      benthos_input%tracerType(benthos_indices%fepa_ind) = 0
      benthos_input%tracerType(benthos_indices%fepb_ind) = 0
      benthos_input%tracerType(benthos_indices%caco3a_ind) = 0
      benthos_input%tracerType(benthos_indices%caco3b_ind) = 0
      benthos_input%tracerType(benthos_indices%camgco3_ind) = 0
      benthos_input%tracerType(benthos_indices%s_ind) = 0
      benthos_input%tracerType(benthos_indices%fes_ind) = 0
      benthos_input%tracerType(benthos_indices%fes2_ind) = 0

!     specify benthos grid points
      nBenthicCells = 0
      indexToBenthosCells(:) = 0
      do iCell = 1, nCellsSolve
         if (bottomDepth(iCell) .le. 400.0_RKIND .and. latCell(iCell) .ge. 50.0_RKIND*deg_to_rad) then
!         if (bottomDepth(iCell) .le. 400.0_RKIND) then
            nBenthicCells = nBenthicCells + 1
            indexToBenthosCells(nBenthicCells) = iCell
            call mpas_log_write('nBenthicCells: $i',MPAS_LOG_ERR,intArgs=(/nBenthicCells/))
            call mpas_log_write('indexToBenthosCells(nBenthicCells): $i',MPAS_LOG_ERR,intArgs=(/indexToBenthosCells(nBenthicCells)/))
         endif
       enddo
   ! prescribed initial concentration array
      call mpas_log_write('After initialize Benthic mask',MPAS_LOG_ERR)
      call mpas_log_write('Final nBenthicCells: $i',MPAS_LOG_ERR,intArgs=(/nBenthicCells/))

      if (.not. config_do_restart_benthos) then
          benthosTracersStorage(:,:) = 0.0_RKIND
          benthosPhPrev3d(:,:) = 0.0_RKIND
          benthosTracers(:,:,:) = 0.0_RKIND

         if (config_use_benthos_monthly_forcing) then

         do iCellBenthic = 1, nBenthicCells !CellsSolve
            iCell = indexToBenthosCells(iCellBenthic)
            do iBenthicTracers = 1, nBenthicTracers
               benthosTracersStorage(iBenthicTracers,iCell) = 0.0_RKIND
            end do
            do iVertLevels = 1, nBenthicVertLevels+1
               benthosPhPrev3d(iVertLevels,iCell) = 8.0_RKIND
               benthosTracers(benthosIndices%poca_ind,iVertLevels,iCell) = 0.158_RKIND
               benthosTracers(benthosIndices%pocb_ind,iVertLevels,iCell) = 1227.2_RKIND
               benthosTracers(benthosIndices%pocc_ind,iVertLevels,iCell) = 3060.1_RKIND
               benthosTracers(benthosIndices%popa_ind,iVertLevels,iCell) = 4.02e-6_RKIND
               benthosTracers(benthosIndices%popb_ind,iVertLevels,iCell) = 4.16_RKIND !*0.15_RKIND*CtoP !
               benthosTracers(benthosIndices%popc_ind,iVertLevels,iCell) = 10.6_RKIND !*0.36_RKIND*CtoP !
               benthosTracers(benthosIndices%pona_ind,iVertLevels,iCell) = 0.024_RKIND !*CtoP/NtoP !
               benthosTracers(benthosIndices%ponb_ind,iVertLevels,iCell) = 118.4_RKIND !*0.15_RKIND*CtoP/NtoP !
               benthosTracers(benthosIndices%ponc_ind,iVertLevels,iCell) = 306.1_RKIND !*0.36_RKIND*CtoP/NtoP !
               benthosTracers(benthosIndices%o2_ind,iVertLevels,iCell) = oceanBottomBenthosO2Data(iCell) !0.0_RKIND   !
               benthosTracers(benthosIndices%nh4_ind,iVertLevels,iCell) = oceanBottomBenthosNh4Data(iCell) !84.3_RKIND  !**
               benthosTracers(benthosIndices%h2po4_ind,iVertLevels,iCell) = oceanBottomBenthosPo4Data(iCell) !36.4_RKIND !** 90.0_RKIND !   umol/L  Reed
               benthosTracers(benthosIndices%co2_ind,iVertLevels,iCell) = 26.5_RKIND  !   mM Krumins
               benthosTracers(benthosIndices%no3_ind,iVertLevels,iCell) = oceanBottomBenthosNo3Data(iCell) !0.0_RKIND  !   mM Krumins
               benthosTracers(benthosIndices%mno2a_ind,iVertLevels,iCell) = 0.0_RKIND!        mmol/kg  Reed et al
               benthosTracers(benthosIndices%mno2b_ind,iVertLevels,iCell) = 0.0_RKIND
               benthosTracers(benthosIndices%mn_ind,iVertLevels,iCell) = 0.0_RKind! 59.0_RKIND !**  umol/L Reed
               benthosTracers(benthosIndices%feoh3a_ind,iVertLevels,iCell) = 0.0_RKIND !    mmol/kg
               benthosTracers(benthosIndices%feoh3b_ind,iVertLevels,iCell) = 0.0_RKIND !    mmol/kg
               benthosTracers(benthosIndices%fe_ind,iVertLevels,iCell) = oceanBottomBenthosFeData(iCell) !0.85_RKIND !**
               benthosTracers(benthosIndices%fepa_ind,iVertLevels,iCell) = 0.0_RKIND !  mmol/kg  Reed
               benthosTracers(benthosIndices%fepb_ind,iVertLevels,iCell) = 0.0_RKIND !  mmol/kg  Reed
               benthosTracers(benthosIndices%so4_ind,iVertLevels,iCell) = 28000.0_RKIND !1080.9_RKIND  !*1000.0_RKIND !  umol/L from mmol/L Reed
               benthosTracers(benthosIndices%h2s_ind,iVertLevels,iCell) = 0.0_RKIND !300.0_RKIND !** !961.3_RKIND !   gets to 600 umol/L at bottom but theres a bottom flux
               benthosTracers(benthosIndices%ch4_ind,iVertLevels,iCell) = 0.0_RKIND !1.3_RKIND ! **
               benthosTracers(benthosIndices%hco3_ind,iVertLevels,iCell) = 2693.2_RKIND !  mM  Krumins et al
               benthosTracers(benthosIndices%fes_ind,iVertLevels,iCell) = 59.7_RKIND !  2000.0 ! 200.0 !  mmol/kg  Reed et al
               benthosTracers(benthosIndices%fes2_ind,iVertLevels,iCell) = 17.8_RKIND !  2000.0 ! 200.0 !  mmol/kg  Reed et al
               benthosTracers(benthosIndices%s_ind,iVertLevels,iCell) = 0.0_RKIND !  2000.0 ! 200.0 !  mmol/kg  Reed et al
               benthosTracers(benthosIndices%co3_ind,iVertLevels,iCell) = 140.4_RKIND !  66.75 ! 65 umol/kg * 1.027 g/m3   Sulpis et al 2018
               benthosTracers(benthosIndices%dic_ind,iVertLevels,iCell) = oceanBottomBenthosDICData(iCell) !2860.1_RKIND 
               benthosTracers(benthosIndices%alk_ind,iVertLevels,iCell) = oceanBottomBenthosAlkData(iCell) !3141.2_RKIND
               benthosTracers(benthosIndices%caco3a_ind,iVertLevels,iCell) = 1000.0_RKIND*.5_RKIND !**  ratio of PIC is 63:13:24...  Roughly  POC:PIC = (:,:,benthosIndices%depth)^(:,:,benthosIndices%-0.56)*64.3   so about 4.8:1
               benthosTracers(benthosIndices%caco3b_ind,iVertLevels,iCell) = 1000.0_RKIND*0.35_RKIND !**
               benthosTracers(benthosIndices%camgco3_ind,iVertLevels,iCell)= 1000.0_RKIND*0.15_RKIND !**
!               call mpas_log_write('initial benthosTracers(benthosIndices%caco3b_ind,iVertLevels,iCell): $r',MPAS_LOG_ERR,realArgs=(/benthosTracers(benthosIndices%caco3b_ind,iVertLevels,iCell)/))

            end do              ! iVertLevels
            if (config_use_step_initial_profiles) then
               do iBenthicTracers = 1,nBenthicTracers
                  do iVertLevels = 1, nBenthicVertLevels+1
                     if (iVertLevels .le. 15 .or. iVertLevels .ge. 26) then
                        benthosTracers(iBenthicTracers,iVertLevels,iCell) = benthosTracers(iBenthicTracers,iVertLevels,iCell)*0.1_RKIND
                     else if (config_set_benthos_flags .eq. 6 .and. iBenthicTracers .eq. benthosIndices%poca_ind) then
                        benthosTracers(iBenthicTracers,iVertLevels,iCell) = 0.0_RKIND
                     end if
                  end do        ! ivertLevels
               end do ! iBenthicTracers
            end if  ! step function
         end do                 ! iCell
         else ! not benthos_monthly_forcing
         do iCellBenthic = 1, nBenthicCells !CellsSolve
            iCell = indexToBenthosCells(iCellBenthic)
            do iBenthicTracers = 1, nBenthicTracers
               benthosTracersStorage(iBenthicTracers,iCell) = 0.0_RKIND
            end do
            do iVertLevels = 1, nBenthicVertLevels+1
               benthosPhPrev3d(iVertLevels,iCell) = 8.0_RKIND
               benthosTracers(benthosIndices%poca_ind,iVertLevels,iCell) = 0.158_RKIND
               benthosTracers(benthosIndices%pocb_ind,iVertLevels,iCell) = 1227.2_RKIND
               benthosTracers(benthosIndices%pocc_ind,iVertLevels,iCell) = 3060.1_RKIND
               benthosTracers(benthosIndices%popa_ind,iVertLevels,iCell) = 4.02e-6_RKIND
               benthosTracers(benthosIndices%popb_ind,iVertLevels,iCell) = 4.16_RKIND !*0.15_RKIND*CtoP !
               benthosTracers(benthosIndices%popc_ind,iVertLevels,iCell) = 10.6_RKIND !*0.36_RKIND*CtoP !
               benthosTracers(benthosIndices%pona_ind,iVertLevels,iCell) = 0.024_RKIND !*CtoP/NtoP !
               benthosTracers(benthosIndices%ponb_ind,iVertLevels,iCell) = 118.4_RKIND !*0.15_RKIND*CtoP/NtoP !
               benthosTracers(benthosIndices%ponc_ind,iVertLevels,iCell) = 306.1_RKIND !*0.36_RKIND*CtoP/NtoP !
               benthosTracers(benthosIndices%o2_ind,iVertLevels,iCell) = 0.0_RKIND   !
               benthosTracers(benthosIndices%nh4_ind,iVertLevels,iCell) = 84.3_RKIND  !**
               benthosTracers(benthosIndices%h2po4_ind,iVertLevels,iCell) =36.4_RKIND !** 90.0_RKIND !   umol/L  Reed
               benthosTracers(benthosIndices%co2_ind,iVertLevels,iCell) = 26.5_RKIND  !   mM Krumins
               benthosTracers(benthosIndices%no3_ind,iVertLevels,iCell) = 0.0_RKIND  !   mM Krumins
               benthosTracers(benthosIndices%mno2a_ind,iVertLevels,iCell) = 0.0_RKIND!        mmol/kg  Reed et al
               benthosTracers(benthosIndices%mno2b_ind,iVertLevels,iCell) = 0.0_RKIND
               benthosTracers(benthosIndices%mn_ind,iVertLevels,iCell) = 59.0_RKIND !**  umol/L Reed
               benthosTracers(benthosIndices%feoh3a_ind,iVertLevels,iCell) = 0.0_RKIND !    mmol/kg
               benthosTracers(benthosIndices%feoh3b_ind,iVertLevels,iCell) = 0.0_RKIND !    mmol/kg
               benthosTracers(benthosIndices%fe_ind,iVertLevels,iCell) = 0.85_RKIND !**  mmol/kg  Reed
               benthosTracers(benthosIndices%fepa_ind,iVertLevels,iCell) = 0.0_RKIND !  mmol/kg  Reed
               benthosTracers(benthosIndices%fepb_ind,iVertLevels,iCell) = 0.0_RKIND !  mmol/kg  Reed
               benthosTracers(benthosIndices%so4_ind,iVertLevels,iCell) = 1080.9_RKIND  !*1000.0_RKIND !  umol/L from mmol/L Reed
               benthosTracers(benthosIndices%h2s_ind,iVertLevels,iCell) = 300.0_RKIND !** !961.3_RKIND !   gets to 600 umol/L at bottom but theres a bottom flux
               benthosTracers(benthosIndices%ch4_ind,iVertLevels,iCell) = 1.3_RKIND ! **
               benthosTracers(benthosIndices%hco3_ind,iVertLevels,iCell) = 2693.2_RKIND !  mM  Krumins et al
               benthosTracers(benthosIndices%fes_ind,iVertLevels,iCell) = 59.7_RKIND !  2000.0 ! 200.0 !  mmol/kg  Reed et al
               benthosTracers(benthosIndices%fes2_ind,iVertLevels,iCell) = 17.8_RKIND !  2000.0 ! 200.0 !  mmol/kg  Reed et al
               benthosTracers(benthosIndices%s_ind,iVertLevels,iCell) = 0.0_RKIND !  2000.0 ! 200.0 !  mmol/kg  Reed et al
               benthosTracers(benthosIndices%co3_ind,iVertLevels,iCell) = 140.4_RKIND !  66.75 ! 65 umol/kg * 1.027 g/m3   Sulpis et al 2018
               benthosTracers(benthosIndices%dic_ind,iVertLevels,iCell) = 2860.1_RKIND !benthosTracers(benthosIndices%co3_ind,iVertLevels,iCell) + benthosTracers(benthosIndices%co2_ind,iVertLevels,iCell) + benthosTracers(benthosIndices%hco3_ind,iVertLevels,iCell)
               benthosTracers(benthosIndices%alk_ind,iVertLevels,iCell) = 3141.2_RKIND !  87 of the ocean conc. (:,:,benthosIndices%mmol/m3) same as ocean conc. Krumins et al
               benthosTracers(benthosIndices%caco3a_ind,iVertLevels,iCell) = 1000.0_RKIND*.5_RKIND !**  ratio of PIC is 63:13:24...  Roughly  POC:PIC = (:,:,benthosIndices%depth)^(:,:,benthosIndices%-0.56)*64.3   so about 4.8:1
               benthosTracers(benthosIndices%caco3b_ind,iVertLevels,iCell) = 1000.0_RKIND*0.35_RKIND !**
               benthosTracers(benthosIndices%camgco3_ind,iVertLevels,iCell)= 1000.0_RKIND*0.15_RKIND !**
!               call mpas_log_write('initial benthosTracers(benthosIndices%caco3b_ind,iVertLevels,iCell): $r',MPAS_LOG_ERR,realArgs=(/benthosTracers(benthosIndices%caco3b_ind,iVertLevels,iCell)/))

            end do              ! iVertLevels
            if (config_use_step_initial_profiles) then
               do iBenthicTracers = 1,nBenthicTracers
                  do iVertLevels = 1, nBenthicVertLevels+1
                     if (iVertLevels .le. 15 .or. iVertLevels .ge. 26) then
                        benthosTracers(iBenthicTracers,iVertLevels,iCell) = benthosTracers(iBenthicTracers,iVertLevels,iCell)*0.1_RKIND
                     else if (config_set_benthos_flags .eq. 6 .and. iBenthicTracers .eq. benthosIndices%poca_ind) then
                        benthosTracers(iBenthicTracers,iVertLevels,iCell) = 0.0_RKIND
                     end if
                  end do        ! ivertLevels
               end do ! iBenthicTracers
            end if  ! step function
         end do                 ! iCell
         endif
      end if                    ! .not. config_do_restart_benthos

      !NJ-TEST: VERIFIED
      test_flag = .false.
      if (test_flag) then
         call mpas_log_write('test benthos_flags_init initialize',MPAS_LOG_ERR)
         do iBenthicTracers = 1, nBenthicTracers
            call mpas_log_write('iBenthicTracers: $i',MPAS_LOG_ERR,intArgs=(/iBenthicTracers/))
            call mpas_log_write('benthos_input%tracerType(iBenthicTracers) : $i',MPAS_LOG_ERR,intArgs=(/benthos_input%tracerType(iBenthicTracers)/))
            call mpas_log_write('maxval(benthosTracers(iBenthicTracers,1,:) : $r',MPAS_LOG_ERR,realArgs=(/maxval(benthosTracers(iBenthicTracers,1,:))/))
            call mpas_log_write('maxval(benthosTracersStorage(iBenthicTracers,:)): $r',MPAS_LOG_ERR,realArgs=(/maxval(benthosTracersStorage(iBenthicTracers,:))/) )
         end do
         call mpas_log_write('maxval(benthosPhPrev3d(nBenthicVertLevels+1,:)): $r',MPAS_LOG_ERR,realArgs=(/maxval(benthosPhPrev3d(nBenthicVertLevels+1,: )) /) )
!         err = 1
!         return
      end if

      else
         call mpas_log_write('not associated benthosTracers', MPAS_LOG_ERR)
      end if  !  associated(benthosTracers)

   !--------------------------------------------------------------------

   end subroutine ocn_benthos_tracer_init!}}}

!***********************************************************************
!
!  routine benthos_tracer_finalize
!
!> \brief   Finalizes Benthos
!> \author  Nicole Jeffery
!> \date    12/04/2023
!> \details
!>  This routine deallocates arrays and indices for benthos tracers
!
!-----------------------------------------------------------------------

   subroutine ocn_benthos_tracer_finalize()!{{{


      deallocate( benthos_indices%short_name )
      deallocate( benthos_indices%long_name )
      deallocate( benthos_indices%units )

      deallocate( benthos_element_indices%short_name )
      deallocate( benthos_element_indices%long_name )
      deallocate( benthos_element_indices%units )

      deallocate ( benthos_input%tracerType)

   end subroutine ocn_benthos_tracer_finalize !}}}

!***********************************************************************

!***********************************************************************
!
!  routine get_benthos_data
!
!> \brief   retrieve data needed to compute benthos ocean climatology
!> \author  Nicole Jeffery
!> \date    04/20/21
!> \details
!>  This routine calls mpas_forcing routines to acquire needed benthos forcing data and interpolates
!>  between time levels.  directly copied from ocn_get_ecosysData
!
!-----------------------------------------------------------------------

    subroutine ocn_get_benthos_data( streamManager, &
        domain, &
        simulationClock, &
        firstTimeStep) !{{{

        type (MPAS_streamManager_type), intent(inout) :: streamManager

        type (domain_type) :: domain
        type (MPAS_timeInterval_type) :: timeStepEco
        type (MPAS_clock_type) :: simulationClock

        logical, intent(in) :: firstTimeStep
        !character(len=strKind), pointer :: config_dt
        real(kind=RKIND) :: dt

        type (mpas_pool_type), pointer :: forcingPool
        type (mpas_pool_type), pointer :: meshPool
        type (mpas_pool_type), pointer :: benthosAuxiliary
        type (mpas_pool_type), pointer :: oceanBenthosPrecipitation, &
           oceanBottomBenthosFields, oceanBenthosFlux
        type (mpas_pool_type), pointer :: benthosMonthlyForcing

! input flux components in oceanBenthosPrecipitation
        real (kind=RKIND), dimension(:), pointer :: &
           oceanBenthosPrecipMno2a, &
           oceanBenthosPrecipMno2b, &
           oceanBenthosPrecipFeoh3a, &
           oceanBenthosPrecipFeoh3b, &
           oceanBenthosPrecipFepa, &
           oceanBenthosPrecipFepb

! input flux components in oceanBenthosFlux
        real (kind=RKIND), dimension(:), pointer :: &
           oceanBenthosPOCFlux, &
           oceanBenthosPICFlux, &
           oceanBenthosMno2Flux, &
           oceanBenthosFeoh3Flux, &
           oceanBenthosBiogenicSilicaFlux, &
           oceanBenthosDustFlux

        ! input flux components in    oceanBottomBenthosFields
        real (kind=RKIND), dimension(:), pointer :: &
           oceanBottomBenthosSio3, &
           oceanBottomBenthosNo3, &
           oceanBottomBenthosNh4, &
           oceanBottomBenthosO2, &
           oceanBottomBenthosDIC, &
           oceanBottomBenthosAlk, &
           oceanBottomBenthosPo4, &
           oceanBottomBenthosFe, &
           oceanBottomBenthosCo2, &
           oceanBottomBenthosHco3, &
           oceanBottomBenthosCo3, &
           oceanBottomBenthosSo4, &
           oceanBottomBenthosMn, &
           oceanBottomBenthosTemperature, &
           oceanBottomBenthosSalinity, &
           oceanBottomBenthosPressure, &
           oceanBottomBenthosDensity

        ! input flux components in oceanBottomBenthosFields
        real (kind=RKIND), dimension(:), pointer :: &
           oceanBenthosPOCFluxData, &
           oceanBenthosPICFluxData, &
           oceanBenthosBiogenicSilicaFluxData, &
           oceanBenthosDustFluxData, &
           oceanBenthosMno2FluxData, &
           oceanBenthosFeoh3FluxData, &
           oceanBenthosPrecipMno2Data, &
           oceanBenthosPrecipFeoh3Data, &
           oceanBottomBenthosSio3Data, &
           oceanBottomBenthosPo4Data, &
           oceanBottomBenthosNh4Data, &
           oceanBottomBenthosNo3Data, &
           oceanBottomBenthosO2Data, &
           oceanBottomBenthosFeData, &
           oceanBottomBenthosMnData, &
           oceanBottomBenthosSo4Data, &
           oceanBottomBenthosDICData, &
           oceanBottomBenthosAlkData, &
           oceanBottomBenthosTemperatureData, &
           oceanBottomBenthosSalinityData, &
           oceanBottomBenthosPressureData, &
           oceanBottomBenthosDensityData

        integer, pointer :: nCells
        integer :: iCell

        !logical, pointer :: config_do_restart, config_do_restart_benthos

        !call MPAS_pool_get_config(domain%configs, 'config_dt', config_dt)
        !call MPAS_pool_get_config(domain%configs, 'config_do_restart', config_do_restart)
        !call MPAS_pool_get_config(domain%configs, 'config_do_restart_benthos', config_do_restart_benthos)

        call mpas_set_timeInterval(timeStepEco,timeString=config_dt)
        call mpas_get_timeInterval(timeStepEco,dt=dt)

        if (firstTimeStep .and. config_do_restart) then
          call MPAS_forcing_get_forcing(forcingGroupHead, &
             'benthosMonthlyClimatology', streamManager, 0.0_RKIND)
        else
          call MPAS_forcing_get_forcing(forcingGroupHead, &
             'benthosMonthlyClimatology', streamManager, dt)
        end if

        call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(forcingPool, 'oceanBenthosFlux', oceanBenthosFlux)
        call mpas_pool_get_subpool(forcingPool, 'oceanBenthosPrecipitation', oceanBenthosPrecipitation)
        call mpas_pool_get_subpool(forcingPool, 'oceanBottomBenthosFields', oceanBottomBenthosFields)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'benthosMonthlyForcing', benthosMonthlyForcing)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosPOCFluxData', oceanBenthosPOCFluxData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosPICFluxData', oceanBenthosPICFluxData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosBiogenicSilicaFluxData', oceanBenthosBiogenicSilicaFluxData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosDustFluxData', oceanBenthosDustFluxData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosFeoh3FluxData', oceanBenthosFeoh3FluxData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosMno2FluxData', oceanBenthosMno2FluxData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosPrecipFeoh3Data', oceanBenthosPrecipFeoh3Data)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBenthosPrecipMno2Data', oceanBenthosPrecipMno2Data)

        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosSio3Data', oceanBottomBenthosSio3Data)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosPo4Data', oceanBottomBenthosPo4Data)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosNh4Data', oceanBottomBenthosNh4Data)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosNo3Data', oceanBottomBenthosNo3Data)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosO2Data', oceanBottomBenthosO2Data)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosFeData', oceanBottomBenthosFeData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosMnData', oceanBottomBenthosMnData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosSo4Data', oceanBottomBenthosSo4Data)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosDICData', oceanBottomBenthosDICData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosAlkData', oceanBottomBenthosAlkData)

        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosDensityData', oceanBottomBenthosDensityData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosTemperatureData', oceanBottomBenthosTemperatureData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosPressureData', oceanBottomBenthosPressureData)
        call mpas_pool_get_array(benthosMonthlyForcing, 'oceanBottomBenthosSalinityData', oceanBottomBenthosSalinityData)

        call mpas_pool_get_array(oceanBenthosFlux, 'oceanBenthosPOCFlux', oceanBenthosPOCFlux )
        call mpas_pool_get_array(oceanBenthosFlux, 'oceanBenthosPICFlux', oceanBenthosPICFlux)
        call mpas_pool_get_array(oceanBenthosFlux, 'oceanBenthosDustFlux', oceanBenthosDustFlux)
        call mpas_pool_get_array(oceanBenthosFlux, 'oceanBenthosBiogenicSilicaFlux', oceanBenthosBiogenicSilicaFlux)
        call mpas_pool_get_array(oceanBenthosFlux, 'oceanBenthosMno2Flux', oceanBenthosMno2Flux)
        call mpas_pool_get_array(oceanBenthosFlux, 'oceanBenthosFeoh3Flux', oceanBenthosFeoh3Flux)

!        call mpas_pool_get_array(oceanBenthosPrecipitation, 'oceanBenthosPrecipMno2a', oceanBenthosPrecipMno2a)
!        call mpas_pool_get_array(oceanBenthosPrecipitation, 'oceanBenthosPrecipMno2b', oceanBenthosPrecipMno2b)
!        call mpas_pool_get_array(oceanBenthosPrecipitation, 'oceanBenthosPrecipFeoh3a', oceanBenthosPrecipFeoh3a)
!        call mpas_pool_get_array(oceanBenthosPrecipitation, 'oceanBenthosPrecipFeoh3b', oceanBenthosPrecipFeoh3b)
!        call mpas_pool_get_array(oceanBenthosPrecipitation, 'oceanBenthosPrecipFepa', oceanBenthosPrecipFepa)
!        call mpas_pool_get_array(oceanBenthosPrecipitation, 'oceanBenthosPrecipFepb', oceanBenthosPrecipFepb)

        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosSio3', oceanBottomBenthosSio3)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosNo3', oceanBottomBenthosNo3)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosNh4', oceanBottomBenthosNh4)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosO2', oceanBottomBenthosO2)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosDIC', oceanBottomBenthosDIC)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosAlk', oceanBottomBenthosAlk)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosPo4', oceanBottomBenthosPo4)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosFe', oceanBottomBenthosFe)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosCo2', oceanBottomBenthosCo2)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosCo3', oceanBottomBenthosCo3)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosHco3', oceanBottomBenthosHco3)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosSo4', oceanBottomBenthosSo4)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosMn', oceanBottomBenthosMn)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosTemperature', oceanBottomBenthosTemperature)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosSalinity', oceanBottomBenthosSalinity)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosPressure', oceanBottomBenthosPressure)
        call mpas_pool_get_array(oceanBottomBenthosFields, 'oceanBottomBenthosDensity', oceanBottomBenthosDensity)

        do iCell = 1, nCells
           ! particulate data is in mmol/m^3*m/s
           oceanBenthosPOCFlux(iCell) = max(0.0_RKIND,oceanBenthosPOCFluxData(iCell))
           oceanBenthosPICFlux(iCell) = max(0.0_RKIND,oceanBenthosPICFluxData(iCell))
           oceanBenthosBiogenicSilicaFlux(iCell) = max(0.0_RKIND,oceanBenthosBiogenicSilicaFluxData(iCell))
           oceanBenthosDustFlux(iCell) = max(0.0_RKIND,oceanBenthosDustFluxData(iCell))
           oceanBenthosMno2Flux(iCell) =  3.6e-13_RKIND ! Arkona Basin
           oceanBenthosFeoh3Flux(iCell) = max(0.0_RKIND,oceanBenthosFeoh3FluxData(iCell))

           oceanBottomBenthosDensity(iCell) = max(1020.0_RKIND,oceanBottomBenthosDensityData(iCell))
           oceanBottomBenthosTemperature(iCell) = oceanBottomBenthosTemperatureData(iCell)
           oceanBottomBenthosSalinity(iCell) = max(30.0_RKIND,oceanBottomBenthosSalinityData(iCell))
           oceanBottomBenthosPressure(iCell) = max(101300.0_RKIND,oceanBottomBenthosPressureData(iCell))
           oceanBottomBenthosMn(iCell) = 4.0e-3_RKIND !Middag 2011 (3-5 nM Barents and Laptev)
           !Colombo 2020 (CAA 3-9 nM; Chukchi Shelf 3-15 nM)
           !(Klunder et al., 2012a; Kondo et al., 2016; Vieira et al., 2019)

           oceanBottomBenthosSio3(iCell) = max(0.0_RKIND,oceanBottomBenthosSio3Data(iCell))
           oceanBottomBenthosNo3(iCell) = max(0.0_RKIND,oceanBottomBenthosNo3Data(iCell))
           oceanBottomBenthosNh4(iCell) = max(0.0_RKIND,oceanBottomBenthosNh4Data(iCell))
           oceanBottomBenthosO2(iCell) = max(0.0_RKIND,oceanBottomBenthosO2Data(iCell))
           oceanBottomBenthosDIC(iCell) = max(0.0_RKIND,oceanBottomBenthosDICData(iCell))
           oceanBottomBenthosAlk(iCell) = max(0.0_RKIND,oceanBottomBenthosAlkData(iCell))
           oceanBottomBenthosPo4(iCell) = max(0.0_RKIND,oceanBottomBenthosPo4Data(iCell))
           oceanBottomBenthosFe(iCell) = max(0.0_RKIND,oceanBottomBenthosFeData(iCell))
           oceanBottomBenthosSo4(iCell) = 28000.0_RKIND

        end do

        !call mpas_log_write('maxval(oceanBottomBenthosTemperature(:)): $r',MPAS_LOG_ERR,realArgs=(/maxval(oceanBottomBenthosTemperature(:))/))
        !call mpas_log_write('minval(oceanBottomBenthosTemperature(:)): $r',MPAS_LOG_ERR,realArgs=(/minval(oceanBottomBenthosTemperature(:))/))

        !call mpas_log_write('maxval(oceanBottomBenthosNO3(:)): $r',MPAS_LOG_ERR,realArgs=(/maxval(oceanBottomBenthosNO3(:))/))
        !call mpas_log_write('minval(oceanBottomBenthosNO3(:)): $r',MPAS_LOG_ERR,realArgs=(/minval(oceanBottomBenthosNO3(:))/))

    end subroutine ocn_get_benthos_data!}}}

!***********************************************************************
!
!  routine get_benthos_terra_data
!
!> \brief   retrieve data from terrigenous sources to the benthos
!> \author  Nicole Jeffery
!> \date    04/20/21
!> \details
!>  This routine calls mpas_forcing routines to acquire for lithogenic and
!> land-derived POC.  directly copied from ocn_get_ecosysData
!
!-----------------------------------------------------------------------

    subroutine ocn_get_benthos_terra_data( streamManager, &
        domain, &
        simulationClock, &
        firstTimeStep) !{{{

        type (MPAS_streamManager_type), intent(inout) :: streamManager

        type (domain_type) :: domain
        type (MPAS_timeInterval_type) :: timeStepEco
        type (MPAS_clock_type) :: simulationClock

        logical, intent(in) :: firstTimeStep
        character(len=strKind), pointer :: config_dt
        real(kind=RKIND) :: dt

        type (mpas_pool_type), pointer :: forcingPool
        type (mpas_pool_type), pointer :: meshPool
        type (mpas_pool_type), pointer :: benthosAuxiliary
        type (mpas_pool_type), pointer :: oceanBenthosFlux
        type (mpas_pool_type), pointer :: benthosMonthlyTerraForcing

! input flux components in oceanBenthosFlux
        real (kind=RKIND), dimension(:), pointer :: &
           oceanTerraBenthosPOCFlux, &
           oceanTerraLithogenicFlux

        ! input flux components in oceanBottomBenthosFields
        real (kind=RKIND), dimension(:), pointer :: &
           oceanTerraBenthosPOCFluxData, &
           oceanTerraLithogenicFluxData

        integer, pointer :: nCells
        integer :: iCell

        logical, pointer :: config_do_restart, config_do_restart_benthos

        call MPAS_pool_get_config(domain%configs, 'config_dt', config_dt)
        call MPAS_pool_get_config(domain%configs, 'config_do_restart', config_do_restart)
        call MPAS_pool_get_config(domain%configs, 'config_do_restart_benthos', config_do_restart_benthos)

        call mpas_set_timeInterval(timeStepEco,timeString=config_dt)
        call mpas_get_timeInterval(timeStepEco,dt=dt)

        if (firstTimeStep .and. config_do_restart) then
          call MPAS_forcing_get_forcing(forcingGroupHead, &
             'benthosMonthlyTerraClimatology', streamManager, 0.0_RKIND)
        else
          call MPAS_forcing_get_forcing(forcingGroupHead, &
             'benthosMonthlyTerraClimatology', streamManager, dt)
        end if

        call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(forcingPool, 'oceanBenthosFlux', oceanBenthosFlux)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'benthosMonthlyTerraForcing', benthosMonthlyTerraForcing)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(benthosMonthlyTerraForcing, 'oceanTerraBenthosPOCFluxData', oceanTerraBenthosPOCFluxData)
        call mpas_pool_get_array(benthosMonthlyTerraForcing, 'oceanTerraLithogenicFluxData', oceanTerraLithogenicFluxData)
        call mpas_pool_get_array(oceanBenthosFlux, 'oceanTerraBenthosPOCFlux', oceanTerraBenthosPOCFlux )
        call mpas_pool_get_array(oceanBenthosFlux, 'oceanTerraLithogenicFlux', oceanTerraLithogenicFlux)

        do iCell = 1, nCells
           ! particulate data is in mmol/m^3*m/s
           oceanTerraBenthosPOCFlux(iCell) = max(0.0_RKIND,oceanTerraBenthosPOCFluxData(iCell))
           oceanTerraLithogenicFlux(iCell) = max(0.0_RKIND,oceanTerraLithogenicFluxData(iCell))
        end do

    end subroutine ocn_get_benthos_terra_data!}}}

!***********************************************************************
!
!  routine ocn_benthos_forcing_write_restart
!
!> \brief   writes restart timestamp for benthos data to be read in on future restart
!> \author  Nicole Jeffery
!> \date    04/30/2021
!
!  Copied after Mathew Maltrud, ocn_ecosys_forcing_write_restart
!-----------------------------------------------------------------------

   subroutine ocn_benthos_forcing_write_restart(domain)!{{{

      type(domain_type) :: domain

      call MPAS_forcing_write_restart_times(forcingGroupHead)

    end subroutine ocn_benthos_forcing_write_restart!}}}

!***********************************************************************
!
!  routine ocn_benthos_terra_forcing_write_restart
!
!> \brief   writes restart timestamp for benthos data to be read in on future restart
!> \author  Nicole Jeffery
!> \date    04/30/2021
!
!  Copied after Mathew Maltrud, ocn_ecosys_forcing_write_restart
!-----------------------------------------------------------------------

   subroutine ocn_benthos_terra_forcing_write_restart(domain)!{{{

      type(domain_type) :: domain

      call MPAS_forcing_write_restart_times(forcingGroupHead)

    end subroutine ocn_benthos_terra_forcing_write_restart!}}}

end module ocn_tracer_benthos

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
