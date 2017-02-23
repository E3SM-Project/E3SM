! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_DMS
!
!> \brief MPAS ocean DMS
!> \author Mathew Maltrud
!> \date   11/01/2015
!> \details
!>  This module contains routines for computing tracer forcing due to DMS
!
!-----------------------------------------------------------------------

module ocn_tracer_DMS

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants

   use DMS_mod
   use DMS_parms
   use BGC_mod
   use BGC_parms

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

   public :: ocn_tracer_DMS_compute, &
             ocn_tracer_DMS_surface_flux_compute,  &
             ocn_tracer_DMS_init

   integer, public:: &
      numColumnsMax

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!-----------------------------------------------------------------------
!  name the necessary DMS derived types
!  all of these are defined in DMS_mod
!-----------------------------------------------------------------------

  type(DMS_indices_type)    , public :: DMS_indices
  type(DMS_input_type)      , public :: DMS_input
  type(DMS_forcing_type)    , public :: DMS_forcing
  type(DMS_output_type)     , public :: DMS_output
  type(DMS_diagnostics_type), public :: DMS_diagnostic_fields
  type(DMS_flux_diagnostics_type), public :: DMS_flux_diagnostic_fields

! hold indices in tracer pool corresponding to each tracer array
  type(DMS_indices_type), public :: dmsIndices
  type(BGC_indices_type) :: ecosysIndices

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_DMS_compute
!
!> \brief   computes a tracer tendency due to DMS
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine computes a tracer tendency due to DMS
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_DMS_compute(activeTracers, DMSTracers, nTracersDMS, ecosysTracers, nTracersEcosys,   &
      forcingPool, nCellsSolve, maxLevelCell,  &
      nVertLevels, layerThickness, indexTemperature, indexSalinity, DMSTracersTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! one dimensional arrays
      integer, dimension(:), intent(in) :: &
         maxLevelCell

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         DMSTracers
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         ecosysTracers
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         activeTracers


      type (mpas_pool_type), intent(in) :: forcingPool

      ! scalars
      integer, intent(in) :: nTracersDMS, nTracersEcosys, nCellsSolve, nVertLevels
      integer, intent(in) :: indexTemperature, indexSalinity


      !
      ! two dimensional pointers
      !
      real (kind=RKIND), dimension(:), pointer :: &
         shortWaveHeatFlux

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         DMSTracersTend

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

      ! source/sink wants cm instead of m
      real (kind=RKIND) :: zTop, zBot, convertLengthScale = 100.0_RKIND

      integer :: iCell, iLevel, iTracer, numColumns, column, iLevelSurface

      call mpas_timer_start("DMS source-sink")

      err = 0

      call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)

      numColumns = 1
      column = 1
      iLevelSurface = 1

      !DWJ 08/05/2016: This loop needs OpenMP added to it.
      do iCell=1,nCellsSolve
         DMS_input%number_of_active_levels(column) = maxLevelCell(iCell)

         DMS_forcing%ShortWaveFlux_surface(column)  = shortWaveHeatFlux(iCell)
         DMS_forcing%SST(column) = activeTracers(indexTemperature,iLevelSurface,iCell)
         DMS_forcing%SSS(column) = activeTracers(indexSalinity,iLevelSurface,iCell)

         do iLevel=1,maxLevelCell(iCell)
            DMS_input%cell_thickness(iLevel,column)    = layerThickness(iLevel,iCell)*convertLengthScale

            DMS_input%DMS_tracers(iLevel,column,DMS_indices%dms_ind)  = DMSTracers(dmsIndices%dms_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%dmsp_ind) = DMSTracers(dmsIndices%dmsp_ind,iLevel,iCell)

            DMS_input%DMS_tracers(iLevel,column,DMS_indices%no3_ind)      = ecosysTracers(ecosysIndices%no3_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%doc_ind)      = ecosysTracers(ecosysIndices%doc_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%zooC_ind)     = ecosysTracers(ecosysIndices%zooC_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%spC_ind)      = ecosysTracers(ecosysIndices%spC_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%spChl_ind)    = ecosysTracers(ecosysIndices%spChl_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%spCaCO3_ind)  = ecosysTracers(ecosysIndices%spCaCO3_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%diatC_ind)    = ecosysTracers(ecosysIndices%diatC_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%diatChl_ind)  = ecosysTracers(ecosysIndices%diatChl_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%phaeoC_ind)   = ecosysTracers(ecosysIndices%phaeoC_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%phaeoChl_ind) = ecosysTracers(ecosysIndices%phaeoChl_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%diazC_ind)    = ecosysTracers(ecosysIndices%diazC_ind,iLevel,iCell)
            DMS_input%DMS_tracers(iLevel,column,DMS_indices%diazChl_ind)  = ecosysTracers(ecosysIndices%diazChl_ind,iLevel,iCell)

         enddo  !  iLevel

         call DMS_SourceSink(DMS_indices, DMS_input, DMS_forcing,   &
                             DMS_output, DMS_diagnostic_fields, nVertLevels,   &
                             numColumnsMax, numColumns)

         do iLevel=1,maxLevelCell(iCell)

            DMSTracersTend(dmsIndices%dms_ind,iLevel,iCell) = DMSTracersTend(dmsIndices%dms_ind,iLevel,iCell)   &
               + DMS_output%DMS_tendencies(iLevel,column,DMS_indices%dms_ind)*layerThickness(iLevel,iCell)
            DMSTracersTend(dmsIndices%dmsp_ind,iLevel,iCell) = DMSTracersTend(dmsIndices%dmsp_ind,iLevel,iCell)   &
               + DMS_output%DMS_tendencies(iLevel,column,DMS_indices%dmsp_ind)*layerThickness(iLevel,iCell)

         enddo

      enddo  !  iCell

      call mpas_timer_stop("DMS source-sink")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_DMS_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_DMS_surface_flux_compute
!
!> \brief   computes a tracer tendency due to DMS
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine computes a tracer tendency due to DMS
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_DMS_surface_flux_compute(activeTracers, DMSTracers, forcingPool,  &
      nTracers, nCellsSolve, zMid, indexTemperature, indexSalinity, DMSSurfaceFlux, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         DMSSurfaceFlux

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         DMSTracers
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         activeTracers

      ! scalars
      integer, intent(in) :: nTracers, nCellsSolve, indexTemperature, indexSalinity

      type (mpas_pool_type), intent(inout) :: forcingPool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

      type (mpas_pool_type), pointer :: ecosysAuxiliary,  &
                                        DMSSeaIceCoupling, &
                                        DMSFluxDiagnostics

      integer :: numColumns, column, iCell, iTracer, iLevelSurface

      real (kind=RKIND), dimension(:), pointer :: &
         atmosphericPressure,   &
         iceFraction,          &
         landIceFraction,      &
         windSpeedSquared10m,  &
         iceFluxDMS,           &
         iceFluxDMSP,          &
         dms_flux_diag_ifrac,  &
         dms_flux_diag_xkw,    &
         dms_flux_diag_atm_press,  &
         dms_flux_diag_pv,     &
         dms_flux_diag_schmidt,&
         dms_flux_diag_sat,    &
         dms_flux_diag_surf,   &
         dms_flux_diag_ws

      real (kind=RKIND) :: &
         renormFluxes = 0.01_RKIND, &
         PascalsToAtmospheres = 1.0_RKIND/101.325e+3_RKIND,  &
         mSquared_to_cmSquared = 1.0e+4_RKIND
!        PascalsToAtmospheres = 1.0_RKIND,  &
!        mSquared_to_cmSquared = 1.0_RKIND
!        PascalsToAtmospheres = 0.0_RKIND,  &
!        mSquared_to_cmSquared = 1.0_RKIND

      err = 0

      call mpas_timer_start("DMS surface flux")

      call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
      call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)
      call mpas_pool_get_array(forcingPool, 'landIceFraction', landIceFraction)

      call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)
      call mpas_pool_get_array(ecosysAuxiliary, 'windSpeedSquared10m', windSpeedSquared10m)

      call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)
      call mpas_pool_get_array(DMSSeaIceCoupling, 'iceFluxDMS', iceFluxDMS)
      call mpas_pool_get_array(DMSSeaIceCoupling, 'iceFluxDMSP', iceFluxDMSP)

      call mpas_pool_get_subpool(forcingPool, 'DMSFluxDiagnostics', DMSFluxDiagnostics)
      call mpas_pool_get_array(DMSFluxDiagnostics, 'dms_flux_diag_xkw', dms_flux_diag_xkw)
      call mpas_pool_get_array(DMSFluxDiagnostics, 'dms_flux_diag_atm_press', dms_flux_diag_atm_press)
      call mpas_pool_get_array(DMSFluxDiagnostics, 'dms_flux_diag_pv', dms_flux_diag_pv)
      call mpas_pool_get_array(DMSFluxDiagnostics, 'dms_flux_diag_schmidt', dms_flux_diag_schmidt)
      call mpas_pool_get_array(DMSFluxDiagnostics, 'dms_flux_diag_sat', dms_flux_diag_sat)
      call mpas_pool_get_array(DMSFluxDiagnostics, 'dms_flux_diag_surf', dms_flux_diag_surf)
      call mpas_pool_get_array(DMSFluxDiagnostics, 'dms_flux_diag_ws', dms_flux_diag_ws)
      call mpas_pool_get_array(DMSFluxDiagnostics, 'dms_flux_diag_ifrac', dms_flux_diag_ifrac)

      DMS_forcing%lcalc_DMS_gas_flux  = .true.

      numColumns = 1
      column = 1
      iLevelSurface = 1

      !DWJ 08/05/2016: This loop needs OpenMP added to it.
      do iCell=1,nCellsSolve

         DMS_forcing%surfacePressure(column) = atmosphericPressure(iCell)*PascalsToAtmospheres
         DMS_forcing%iceFraction(column) = iceFraction(iCell)
!maltrud assume for now that if there is any land ice, it is all land ice
         if (associated(landIceFraction)) then
            if (landIceFraction(iCell) > 0.0_RKIND) DMS_forcing%iceFraction(column) = 1.0_RKIND
         endif
         DMS_forcing%windSpeedSquared10m(column) = windSpeedSquared10m(iCell)*mSquared_to_cmSquared
         DMS_forcing%SST(column) = activeTracers(indexTemperature,iLevelSurface,iCell)
         DMS_forcing%SSS(column) = activeTracers(indexSalinity,iLevelSurface,iCell)

         DMS_input%DMS_tracers(1,column,DMS_indices%dms_ind)  = DMSTracers(dmsIndices%dms_ind,1,iCell)
         DMS_input%DMS_tracers(1,column,DMS_indices%dmsp_ind) = DMSTracers(dmsIndices%dmsp_ind,1,iCell)

         call DMS_SurfaceFluxes(DMS_indices, DMS_input, DMS_forcing,   &
                                DMS_flux_diagnostic_fields,   &
                                numColumnsMax, column)

         DMSSurfaceFlux(dmsIndices%dms_ind,iCell) = DMS_forcing%netFlux(column,DMS_indices%dms_ind)*renormFluxes +  &
            iceFluxDMS(iCell)
         DMSSurfaceFlux(dmsIndices%dmsp_ind,iCell) = DMS_forcing%netFlux(column,DMS_indices%dmsp_ind)*renormFluxes +  &
            iceFluxDMSP(iCell)

         dms_flux_diag_ifrac(iCell)     = DMS_flux_diagnostic_fields%diag_DMS_IFRAC(column)
         dms_flux_diag_xkw(iCell)       = DMS_flux_diagnostic_fields%diag_DMS_XKW(column)
         dms_flux_diag_atm_press(iCell) = DMS_flux_diagnostic_fields%diag_DMS_ATM_PRESS(column)
         dms_flux_diag_pv(iCell)        = DMS_flux_diagnostic_fields%diag_DMS_PV(column)
         dms_flux_diag_schmidt(iCell)   = DMS_flux_diagnostic_fields%diag_DMS_SCHMIDT(column)
         dms_flux_diag_sat(iCell)       = DMS_flux_diagnostic_fields%diag_DMS_SAT(column)
         dms_flux_diag_surf(iCell)      = DMS_flux_diagnostic_fields%diag_DMS_SURF(column)
         dms_flux_diag_ws(iCell)        = DMS_flux_diagnostic_fields%diag_DMS_WS(column)

      enddo  !  iCell

      call mpas_timer_stop("DMS surface flux")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_DMS_surface_flux_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_DMS_init
!
!> \brief   Initializes ocean surface restoring
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine initializes fields required for tracer surface flux restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_DMS_init(domain,err)!{{{

!NOTE:  called from mpas_ocn_forward_mode.F

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool

      ! three dimensional pointers
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        DMSTracers

      ! scalars
      integer :: nTracers, numColumnsMax

      ! scalar pointers
      integer, pointer :: nVertLevels, index_dummy

      !
      ! get tracers pools
      !

      err = 0

      !
      ! Get tracer group so we can get the number of tracers in it
      !

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'DMSTracers', DMSTracers, 1)

      ! make sure DMS is turned on

      if (associated(DMSTracers)) then

      ! cannot use DMS_tracer_cnt since it has dms, dmsp, and 12 ecosys fields

      nTracers = size(DMSTracers, dim=1)
      if (nTracers /= 2) then
         err = 1
         return
      endif

      !
      ! pull nVertLevels out of the mesh structure
      !

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

!-----------------------------------------------------------------------
!  initialize DMS parameters
!-----------------------------------------------------------------------

   allocate( DMS_indices%short_name(DMS_tracer_cnt) )
   allocate( DMS_indices%long_name(DMS_tracer_cnt) )
   allocate( DMS_indices%units(DMS_tracer_cnt) )

! no need to allocate the above fields for dmsIndices (?)

!-----------------------------------------------------------------------
!  sets most of DMS parameters
!  sets namelist defaults
!-----------------------------------------------------------------------

   call DMS_parms_init

      ! modify namelist values here....

      !
      ! for now only do 1 column at a time
      !
      numColumnsMax = 1

      DMS_indices%dms_ind      = 1
      DMS_indices%dmsp_ind     = 2
      DMS_indices%no3_ind      = 3
      DMS_indices%doc_ind      = 4
      DMS_indices%zooC_ind     = 5
      DMS_indices%spC_ind      = 6
      DMS_indices%spCaCO3_ind  = 7
      DMS_indices%diatC_ind    = 8
      DMS_indices%diazC_ind    = 9
      DMS_indices%phaeoC_ind   = 10
      DMS_indices%spChl_ind    = 11
      DMS_indices%diatChl_ind  = 12
      DMS_indices%diazChl_ind  = 13
      DMS_indices%phaeoChl_ind = 14

      call mpas_pool_get_dimension(tracersPool, 'index_DMS',  index_dummy)
      dmsIndices%dms_ind  = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DMSP', index_dummy)
      dmsIndices%dmsp_ind = index_dummy

      call mpas_pool_get_dimension(tracersPool, 'index_NO3',      index_dummy)
      ecosysIndices%no3_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DOC',      index_dummy)
      ecosysIndices%doc_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_zooC',     index_dummy)
      ecosysIndices%zooC_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spChl',    index_dummy)
      ecosysIndices%spChl_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spC',      index_dummy)
      ecosysIndices%spC_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spCaCO3',  index_dummy)
      ecosysIndices%spCaCO3_ind         = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatChl',  index_dummy)
      ecosysIndices%diatChl_ind         = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatC',    index_dummy)
      ecosysIndices%diatC_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diazChl',  index_dummy)
      ecosysIndices%diazChl_ind         = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diazC',    index_dummy)
      ecosysIndices%diazC_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_phaeoChl', index_dummy)
      ecosysIndices%phaeoChl_ind        = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_phaeoC',   index_dummy)
      ecosysIndices%phaeoC_ind          = index_dummy

! DMS_init sets short and long names, units in DMS_indices

      call DMS_init(DMS_indices)

!NOTES:

!also check short_name with mpas variable name

!-----------------------------------------------------------------------
!  allocate input, forcing, diagnostic arrays
!-----------------------------------------------------------------------

      allocate ( DMS_input%DMS_tracers(nVertLevels, numColumnsMax, DMS_tracer_cnt) )
      allocate ( DMS_input%cell_thickness(nVertLevels, numColumnsMax) )
      allocate ( DMS_input%number_of_active_levels(numColumnsMax) )

      allocate ( DMS_forcing%ShortWaveFlux_surface(numColumnsMax) )
      allocate ( DMS_forcing%surfacePressure(numColumnsMax) )
      allocate ( DMS_forcing%iceFraction(numColumnsMax) )
      allocate ( DMS_forcing%windSpeedSquared10m(numColumnsMax) )
      allocate ( DMS_forcing%SST(numColumnsMax) )
      allocate ( DMS_forcing%SSS(numColumnsMax) )

      allocate ( DMS_forcing%netFlux(numColumnsMax, DMS_tracer_cnt) )

      allocate ( DMS_output%DMS_tendencies(nVertLevels, numColumnsMax, DMS_tracer_cnt) )

    !---------------------------------------------------------------------------
    !   allocate flux diagnostic output fields
    !---------------------------------------------------------------------------

    allocate (DMS_flux_diagnostic_fields%diag_DMS_IFRAC(numColumnsMax) )
    allocate (DMS_flux_diagnostic_fields%diag_DMS_XKW(numColumnsMax) )
    allocate (DMS_flux_diagnostic_fields%diag_DMS_ATM_PRESS(numColumnsMax) )
    allocate (DMS_flux_diagnostic_fields%diag_DMS_PV(numColumnsMax) )
    allocate (DMS_flux_diagnostic_fields%diag_DMS_SCHMIDT(numColumnsMax) )
    allocate (DMS_flux_diagnostic_fields%diag_DMS_SAT(numColumnsMax) )
    allocate (DMS_flux_diagnostic_fields%diag_DMS_SURF(numColumnsMax) )
    allocate (DMS_flux_diagnostic_fields%diag_DMS_WS(numColumnsMax) )

    !---------------------------------------------------------------------------
    !   allocate diagnostic output fields
    !---------------------------------------------------------------------------

    allocate (DMS_diagnostic_fields%diag_DMS_S_DMSP(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMS_S_TOTAL(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMS_R_B(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMS_R_PHOT(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMS_R_BKGND(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMS_R_TOTAL(nVertLevels, numColumnsMax) )

    allocate (DMS_diagnostic_fields%diag_DMSP_S_PHAEO(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMSP_S_NONPHAEO(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMSP_S_ZOO(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMSP_S_TOTAL(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMSP_R_B(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMSP_R_BKGND(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_DMSP_R_TOTAL(nVertLevels, numColumnsMax) )

    allocate (DMS_diagnostic_fields%diag_Cyano_frac(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_Cocco_frac(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_Eukar_frac(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_diatS(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_diatN(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_phytoN(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_coccoS(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_cyanoS(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_eukarS(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_diazS(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_phaeoS(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_zooS(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_zooCC(nVertLevels, numColumnsMax) )
    allocate (DMS_diagnostic_fields%diag_RSNzoo(nVertLevels, numColumnsMax) )

    end if  !  associated(DMS_tracers)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_DMS_init!}}}

!***********************************************************************

end module ocn_tracer_DMS

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
