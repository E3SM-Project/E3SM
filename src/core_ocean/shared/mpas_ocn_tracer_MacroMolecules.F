! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_MacroMolecules
!
!> \brief MPAS ocean MacroMolecules
!> \author Mathew Maltrud
!> \date   11/01/2015
!> \details
!>  This module contains routines for computing tracer forcing due to MacroMolecules
!
!-----------------------------------------------------------------------

module ocn_tracer_MacroMolecules

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants

   use MACROS_mod
   use MACROS_parms
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

   public :: ocn_tracer_MacroMolecules_compute, &
             ocn_tracer_MacroMolecules_surface_flux_compute,  &
             ocn_tracer_MacroMolecules_init

   integer, public:: &
      numColumnsMax

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!-----------------------------------------------------------------------
!  name the necessary MacroMolecules derived types
!  all of these are defined in MacroMolecules_mod
!-----------------------------------------------------------------------

  type(MACROS_indices_type)    , public :: MacroMolecules_indices
  type(MACROS_input_type)      , public :: MacroMolecules_input
  type(MACROS_output_type)     , public :: MacroMolecules_output
  type(MACROS_diagnostics_type), public :: MacroMolecules_diagnostic_fields

! hold indices in tracer pool corresponding to each tracer array
  type(MACROS_indices_type), public :: macrosIndices
  type(BGC_indices_type) :: ecosysIndices

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_MacroMolecules_compute
!
!> \brief   computes a tracer tendency due to MacroMolecules
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine computes a tracer tendency due to MacroMolecules
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_MacroMolecules_compute(MacroMoleculesTracers, nTracersMacroMolecules,   &
      ecosysTracers, nTracersEcosys, forcingPool, &
      nCellsSolve, maxLevelCell, nVertLevels, layerThickness, MacroMoleculesTracersTend, err)!{{{

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
         MacroMoleculesTracers
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         ecosysTracers

      type (mpas_pool_type), intent(in) :: forcingPool

      ! scalars
      integer, intent(in) :: nTracersMacroMolecules, nTracersEcosys, nCellsSolve, nVertLevels

      !
      ! two dimensional pointers
      !
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         MacroMoleculesTracersTend

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

      integer :: iCell, iLevel, iTracer, numColumns, column

      call mpas_timer_start("MacroMolecules source-sink")

      err = 0

      numColumns = 1
      column = 1
      !DWJ 08/05/2016: This loop needs OpenMP added to it
      do iCell=1,nCellsSolve
         MacroMolecules_input%number_of_active_levels(column) = maxLevelCell(iCell)
         do iLevel=1,maxLevelCell(iCell)
            MacroMolecules_input%cell_thickness(iLevel,column)    = layerThickness(iLevel,iCell)*convertLengthScale

            MacroMolecules_input%MACROS_tracers(iLevel,column,MacroMolecules_indices%prot_ind)  =   &
               MacroMoleculesTracers(macrosIndices%prot_ind,iLevel,iCell)
            MacroMolecules_input%MACROS_tracers(iLevel,column,MacroMolecules_indices%poly_ind) =   &
               MacroMoleculesTracers(macrosIndices%poly_ind,iLevel,iCell)
            MacroMolecules_input%MACROS_tracers(iLevel,column,MacroMolecules_indices%lip_ind) =   &
               MacroMoleculesTracers(macrosIndices%lip_ind,iLevel,iCell)

            MacroMolecules_input%MACROS_tracers(iLevel,column,MacroMolecules_indices%zooC_ind)     =   &
               ecosysTracers(ecosysIndices%zooC_ind,iLevel,iCell)
            MacroMolecules_input%MACROS_tracers(iLevel,column,MacroMolecules_indices%spC_ind)      =   &
               ecosysTracers(ecosysIndices%spC_ind,iLevel,iCell)
            MacroMolecules_input%MACROS_tracers(iLevel,column,MacroMolecules_indices%diatC_ind)    =   &
               ecosysTracers(ecosysIndices%diatC_ind,iLevel,iCell)
            MacroMolecules_input%MACROS_tracers(iLevel,column,MacroMolecules_indices%phaeoC_ind)   =   &
               ecosysTracers(ecosysIndices%phaeoC_ind,iLevel,iCell)
            MacroMolecules_input%MACROS_tracers(iLevel,column,MacroMolecules_indices%diazC_ind)    =   &
               ecosysTracers(ecosysIndices%diazC_ind,iLevel,iCell)

         enddo  !  iLevel

         call MACROS_SourceSink(MacroMolecules_indices, MacroMolecules_input,  &
                             MacroMolecules_output, MacroMolecules_diagnostic_fields, nVertLevels,   &
                             numColumnsMax, numColumns)

         do iLevel=1,maxLevelCell(iCell)

            MacroMoleculesTracersTend(macrosIndices%prot_ind,iLevel,iCell) = &
               MacroMoleculesTracersTend(macrosIndices%prot_ind,iLevel,iCell)   &
                  + MacroMolecules_output%MACROS_tendencies(iLevel,column,MacroMolecules_indices%prot_ind)
            MacroMoleculesTracersTend(macrosIndices%poly_ind,iLevel,iCell) = &
               MacroMoleculesTracersTend(macrosIndices%poly_ind,iLevel,iCell)   &
                  + MacroMolecules_output%MACROS_tendencies(iLevel,column,MacroMolecules_indices%poly_ind)
            MacroMoleculesTracersTend(macrosIndices%lip_ind,iLevel,iCell) = &
               MacroMoleculesTracersTend(macrosIndices%lip_ind,iLevel,iCell)   &
                  + MacroMolecules_output%MACROS_tendencies(iLevel,column,MacroMolecules_indices%lip_ind)

         enddo

      enddo  !  iCell

      call mpas_timer_stop("MacroMolecules source-sink")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_MacroMolecules_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_MacroMolecules_surface_flux_compute
!
!> \brief   computes a tracer tendency due to MacroMolecules
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine computes a tracer tendency due to MacroMolecules
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_MacroMolecules_surface_flux_compute(activeTracers, MacroMoleculesTracers, forcingPool,  &
      nTracers, nCellsSolve, zMid, indexTemperature, indexSalinity, MacroMoleculesSurfaceFlux, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         MacroMoleculesSurfaceFlux

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         MacroMoleculesTracers
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

      integer :: iCell

      call mpas_timer_start("MacroMolecules surface flux")

      err = 0

      ! fluxes are zero

      !DWJ 08/05/2016: This loop needs OpenMP added to it
      do iCell = 1, nCellsSolve

         MacroMoleculesSurfaceFlux(macrosIndices%prot_ind,iCell) = 0.0_RKIND
         MacroMoleculesSurfaceFlux(macrosIndices%poly_ind,iCell) = 0.0_RKIND
         MacroMoleculesSurfaceFlux(macrosIndices%lip_ind, iCell) = 0.0_RKIND

      enddo  !  iCell

      call mpas_timer_stop("MacroMolecules surface flux")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_MacroMolecules_surface_flux_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_MacroMolecules_init
!
!> \brief   Initializes ocean surface restoring
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine initializes fields required for tracer surface flux restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_MacroMolecules_init(domain,err)!{{{

!NOTE:  called from mpas_ocn_forward_mode.F

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool

      ! three dimensional pointers
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        MacroMoleculesTracers

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
      call mpas_pool_get_array(tracersPool, 'MacroMoleculesTracers', MacroMoleculesTracers, 1)

      ! make sure MacrosMolecules is turned on

      if (associated(MacroMoleculesTracers)) then

      ! cannot use MacroMolecules_tracer_cnt since it has poly, prot, lip and 5 ecosys fields

      nTracers = size(MacroMoleculesTracers, dim=1)
      if (nTracers /= 3) then
         err = 1
         return
      endif

      !
      ! pull nVertLevels out of the mesh structure
      !

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

!-----------------------------------------------------------------------
!  initialize MacroMolecules parameters
!-----------------------------------------------------------------------

   allocate( MacroMolecules_indices%short_name(MACROS_tracer_cnt) )
   allocate( MacroMolecules_indices%long_name(MACROS_tracer_cnt) )
   allocate( MacroMolecules_indices%units(MACROS_tracer_cnt) )

! no need to allocate the above fields for macrosIndices (?)

!-----------------------------------------------------------------------
!  sets most of MacroMolecules parameters
!  sets namelist defaults
!-----------------------------------------------------------------------

   call MACROS_parms_init

! modify namelist values here....

      !
      ! for now only do 1 column at a time
      !
      numColumnsMax = 1

      MacroMolecules_indices%prot_ind     = 1
      MacroMolecules_indices%poly_ind     = 2
      MacroMolecules_indices%lip_ind      = 3
      MacroMolecules_indices%zooC_ind     = 4
      MacroMolecules_indices%spC_ind      = 5
      MacroMolecules_indices%diatC_ind    = 6
      MacroMolecules_indices%diazC_ind    = 7
      MacroMolecules_indices%phaeoC_ind   = 8

      call mpas_pool_get_dimension(tracersPool, 'index_PROT',  index_dummy)
      macrosIndices%prot_ind  = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_POLY', index_dummy)
      macrosIndices%poly_ind = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_LIP', index_dummy)
      macrosIndices%lip_ind = index_dummy

      call mpas_pool_get_dimension(tracersPool, 'index_zooC',     index_dummy)
      ecosysIndices%zooC_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spC',      index_dummy)
      ecosysIndices%spC_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatC',    index_dummy)
      ecosysIndices%diatC_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diazC',    index_dummy)
      ecosysIndices%diazC_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_phaeoC',   index_dummy)
      ecosysIndices%phaeoC_ind          = index_dummy

! MacroMolecules_init sets short and long names, units in MacroMolecules_indices

      call MACROS_init(MacroMolecules_indices)

!NOTES:

!also check short_name with mpas variable name

!-----------------------------------------------------------------------
!  allocate input, forcing, diagnostic arrays
!-----------------------------------------------------------------------

      allocate ( MacroMolecules_input%MACROS_tracers(nVertLevels, numColumnsMax, MACROS_tracer_cnt) )
      allocate ( MacroMolecules_input%cell_thickness(nVertLevels, numColumnsMax) )
      allocate ( MacroMolecules_input%number_of_active_levels(numColumnsMax) )

      allocate ( MacroMolecules_output%MACROS_tendencies(nVertLevels, numColumnsMax, MACROS_tracer_cnt) )

    !---------------------------------------------------------------------------
    !   allocate diagnostic output fields
    !---------------------------------------------------------------------------

    allocate (MacroMolecules_diagnostic_fields%diag_PROT_S_TOTAL(nVertLevels, numColumnsMax) )
    allocate (MacroMolecules_diagnostic_fields%diag_POLY_S_TOTAL(nVertLevels, numColumnsMax) )
    allocate (MacroMolecules_diagnostic_fields%diag_LIP_S_TOTAL(nVertLevels, numColumnsMax) )
    allocate (MacroMolecules_diagnostic_fields%diag_PROT_R_TOTAL(nVertLevels, numColumnsMax) )
    allocate (MacroMolecules_diagnostic_fields%diag_POLY_R_TOTAL(nVertLevels, numColumnsMax) )
    allocate (MacroMolecules_diagnostic_fields%diag_LIP_R_TOTAL(nVertLevels, numColumnsMax) )

    end if  !  associated(MacroMoleculesTracers)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_MacroMolecules_init!}}}

!***********************************************************************

end module ocn_tracer_MacroMolecules

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
