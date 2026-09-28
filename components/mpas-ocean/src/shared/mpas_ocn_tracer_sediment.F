! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_sediment
!
!> \brief MPAS-Ocean sediment tracer module
!> \author Nicole Jeffery, LANL
!> \date   07/2026
!> \details
!>  This module implements suspended sediment transport as an independent
!>  tracer class in MPAS-Ocean.  It supports up to 3 runtime-configurable
!>  sediment size classes, each with distinct sinking rates and critical
!>  shear stresses for erosion.
!>
!>  Physical parameterizations:
!>    - Vertical settling: semi-Lagrangian PPM scheme (CFL-free),
!>      implemented in ocn_sediment_settling (mpas_ocn_sediment_settling.F)
!>    - Bottom erosion: critical shear stress formulation (CROCO/MUSTANG)
!>    - Horizontal transport: existing MPAS monotone tracer advection
!>    - Sources (Phase 3): river runoff, ice-ocean coupling
!>
!>  Backward compatibility:
!>    Setting config_sediment_tracer_count = 1 with
!>    config_sediment_sinking_rate_1  = -0.001
!>    config_sediment_critical_shear_1 = 0.2
!>    config_sediment_erosion_rate     = 4.0e-5
!>    reproduces the toy DOCr test results.
!>
!>  References:
!>    Colella & Woodward 1984 (PPM)
!>    CROCO MUSTANG sediment module
!>    Thomsen & Flach 2000 (critical shear stress)
!
!-----------------------------------------------------------------------

module ocn_tracer_sediment

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

   use ocn_config
   use ocn_mesh
   use ocn_constants

   use ocn_sediment_settling

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_tracer_sediment_init,           &
             ocn_tracer_sediment_tend,            &
             ocn_tracer_sediment_surface_flux_compute

   !--------------------------------------------------------------------
   !
   ! Private module variables: sediment class properties read from namelist
   !
   !--------------------------------------------------------------------

   integer :: nSedimentClasses   !< active number of size classes (1-3)

   !> Sinking rates (m/s, negative = downward) indexed 1:nSedimentClasses
   real (kind=RKIND), dimension(3) :: sinkingRate

   !> Critical bed shear stress for erosion (N/m^2)
   real (kind=RKIND), dimension(3) :: criticalShear

   !> Shared parameters
   real (kind=RKIND) :: particleDensity   !< kg/m^3
   real (kind=RKIND) :: erosionRate       !< ERATE (kg/m^2/s)
   real (kind=RKIND) :: bedPorosity       !< dimensionless
   real (kind=RKIND) :: riverSedimentLoad !< kg/m3

!***********************************************************************

   contains

!***********************************************************************
!
!  routine ocn_tracer_sediment_init
!
!> \brief   Initialises the sediment tracer module
!> \author  Nicole Jeffery
!> \date    07/2026
!> \details
!>  Reads per-class properties from namelist, zeros sediment auxiliary
!>  fields (bed mass, fluxes), and stores sinking rates/critical shear
!>  in module-level arrays for use during each time step.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_sediment_init(domain, err)!{{{

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      integer, intent(out) :: err                 !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: sedimentAuxiliary

      real (kind=RKIND), dimension(:), pointer :: &
         sedimentBedMass1,      &
         sedimentBedMass2,      &
         sedimentBedMass3,      &
         sedimentSettlingFlux1, &
         sedimentSettlingFlux2, &
         sedimentSettlingFlux3, &
         sedimentErosionFlux1,  &
         sedimentErosionFlux2,  &
         sedimentErosionFlux3,  &
         bottomStressCellSediment

      ! namelist config pointers
      integer, pointer   :: config_sediment_tracer_count
      real (kind=RKIND), pointer :: &
         config_sediment_sinking_rate_1, &
         config_sediment_sinking_rate_2, &
         config_sediment_sinking_rate_3, &
         config_sediment_critical_shear_1, &
         config_sediment_critical_shear_2, &
         config_sediment_critical_shear_3, &
         config_sediment_particle_density, &
         config_sediment_erosion_rate,     &
         config_sediment_bed_porosity

      err = 0

      !--- read namelist class count and properties ---

      call mpas_pool_get_config(ocnConfigs, 'config_sediment_tracer_count',    &
                                 config_sediment_tracer_count)
      call mpas_pool_get_config(ocnConfigs, 'config_sediment_sinking_rate_1',  &
                                 config_sediment_sinking_rate_1)
      call mpas_pool_get_config(ocnConfigs, 'config_sediment_sinking_rate_2',  &
                                 config_sediment_sinking_rate_2)
      call mpas_pool_get_config(ocnConfigs, 'config_sediment_sinking_rate_3',  &
                                 config_sediment_sinking_rate_3)
      call mpas_pool_get_config(ocnConfigs, 'config_sediment_critical_shear_1',&
                                 config_sediment_critical_shear_1)
      call mpas_pool_get_config(ocnConfigs, 'config_sediment_critical_shear_2',&
                                 config_sediment_critical_shear_2)
      call mpas_pool_get_config(ocnConfigs, 'config_sediment_critical_shear_3',&
                                 config_sediment_critical_shear_3)
      call mpas_pool_get_config(ocnConfigs, 'config_sediment_particle_density',&
                                 config_sediment_particle_density)
      call mpas_pool_get_config(ocnConfigs, 'config_sediment_erosion_rate',    &
                                 config_sediment_erosion_rate)
      call mpas_pool_get_config(ocnConfigs, 'config_sediment_bed_porosity',    &
                                 config_sediment_bed_porosity)

      nSedimentClasses  = config_sediment_tracer_count

      sinkingRate(1)    = config_sediment_sinking_rate_1
      sinkingRate(2)    = config_sediment_sinking_rate_2
      sinkingRate(3)    = config_sediment_sinking_rate_3

      criticalShear(1)  = config_sediment_critical_shear_1
      criticalShear(2)  = config_sediment_critical_shear_2
      criticalShear(3)  = config_sediment_critical_shear_3

      particleDensity   = config_sediment_particle_density
      erosionRate       = config_sediment_erosion_rate
      bedPorosity       = config_sediment_bed_porosity
      riverSedimentLoad = config_river_sediment_load

      if (nSedimentClasses < 1 .or. nSedimentClasses > 3) then
         call mpas_log_write( &
            'ocn_tracer_sediment_init: config_sediment_tracer_count must be 1, 2, or 3', &
            MPAS_LOG_CRIT)
         err = 1
         return
      endif

      call mpas_log_write('ocn_tracer_sediment_init: initializing $i sediment class(es)', &
                          intArgs=(/nSedimentClasses/))

      !--- zero sediment auxiliary fields ---

      call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(forcingPool, 'sedimentAuxiliary', sedimentAuxiliary)

      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentBedMass1',      sedimentBedMass1)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentBedMass2',      sedimentBedMass2)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentBedMass3',      sedimentBedMass3)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentSettlingFlux1', sedimentSettlingFlux1)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentSettlingFlux2', sedimentSettlingFlux2)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentSettlingFlux3', sedimentSettlingFlux3)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentErosionFlux1',  sedimentErosionFlux1)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentErosionFlux2',  sedimentErosionFlux2)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentErosionFlux3',  sedimentErosionFlux3)
      call mpas_pool_get_array(sedimentAuxiliary, 'bottomStressCellSediment', bottomStressCellSediment)

      sedimentBedMass1(:)       = 0.0_RKIND  ! NOTE: replace with restart read in Phase 2
      sedimentBedMass2(:)       = 0.0_RKIND
      sedimentBedMass3(:)       = 0.0_RKIND
      sedimentSettlingFlux1(:)  = 0.0_RKIND
      sedimentSettlingFlux2(:)  = 0.0_RKIND
      sedimentSettlingFlux3(:)  = 0.0_RKIND
      sedimentErosionFlux1(:)   = 0.0_RKIND
      sedimentErosionFlux2(:)   = 0.0_RKIND
      sedimentErosionFlux3(:)   = 0.0_RKIND
      bottomStressCellSediment(:) = 0.0_RKIND

   end subroutine ocn_tracer_sediment_init!}}}

!***********************************************************************
!
!  routine ocn_tracer_sediment_tend
!
!> \brief   Computes sediment settling and erosion tendencies
!> \author  Nicole Jeffery
!> \date    07/2026
!> \details
!>  For each active sediment size class, this routine:
!>    1. Computes bottom shear stress from near-bottom velocity
!>    2. Calls the semi-Lagrangian PPM settling routine
!>       (ocn_sediment_settling_flx) which returns both the settling
!>       tendency and updates the seafloor bed mass
!>    3. Stores per-class settling and erosion flux diagnostics
!>
!>  Horizontal transport is handled by the standard MPAS tracer advection
!>  framework; this routine provides only the vertical settling tendency.
!>
!>  Phase 2 integration is active through ocn_tend_tracer in
!>  mpas_ocn_tendency.F.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_sediment_tend(statePool, forcingPool, meshPool, &
                                       sedimentTracers, sedimentTend, &
                                       kineticEnergyCell, layerThickEdgeDrag, &
                                       layerThickness, dt, &
                                       timeLevel, nCellsSolve, err) !{{{

      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(inout) :: forcingPool
      type (mpas_pool_type), intent(in) :: meshPool
      integer, intent(in)  :: timeLevel           !< Input: time level
      real (kind=RKIND), intent(in) :: dt         !< Input: time step (s)
      integer, intent(in) :: nCellsSolve
      integer, intent(out) :: err                 !< Output: error flag

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         sedimentTracers   !< (nClass, nVertLevels, nCells)

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         sedimentTend      !< (nClass, nVertLevels, nCells)

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell, &
         layerThickEdgeDrag, &
         layerThickness    !< (nVertLevels, nCells)

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: sedimentAuxiliary

      real (kind=RKIND), dimension(:,:), pointer :: &
         normalVelocity

      real (kind=RKIND), dimension(:), pointer :: &
         bottomStressCellSediment, &
         sedimentBedMass1, sedimentBedMass2, sedimentBedMass3, &
         sedimentSettlingFlux1, sedimentSettlingFlux2, sedimentSettlingFlux3, &
         sedimentErosionFlux1,  sedimentErosionFlux2,  sedimentErosionFlux3

      ! per-cell working arrays for each class
      real (kind=RKIND), dimension(:), allocatable :: &
         workTend, workConc

      integer :: iCell, iClass, kmin, kmax

      err = 0

      call mpas_pool_get_subpool(forcingPool, 'sedimentAuxiliary', sedimentAuxiliary)

      call mpas_pool_get_array(statePool,   'normalVelocity',    normalVelocity,    timeLevel)

      call mpas_pool_get_array(sedimentAuxiliary, 'bottomStressCellSediment', bottomStressCellSediment)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentBedMass1',      sedimentBedMass1)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentBedMass2',      sedimentBedMass2)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentBedMass3',      sedimentBedMass3)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentSettlingFlux1', sedimentSettlingFlux1)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentSettlingFlux2', sedimentSettlingFlux2)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentSettlingFlux3', sedimentSettlingFlux3)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentErosionFlux1',  sedimentErosionFlux1)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentErosionFlux2',  sedimentErosionFlux2)
      call mpas_pool_get_array(sedimentAuxiliary, 'sedimentErosionFlux3',  sedimentErosionFlux3)

      !--- compute bottom shear stress (once, shared across all classes) ---

      call ocn_compute_bottom_shear_stress(normalVelocity, kineticEnergyCell, &
                                            layerThickEdgeDrag, bottomStressCellSediment)

      !--- loop over active sediment classes ---

      allocate(workTend(nVertLevels), workConc(nVertLevels))

      do iCell = 1, nCellsSolve

         kmin = minLevelCell(iCell)
         kmax = maxLevelCell(iCell)

         do iClass = 1, nSedimentClasses

            workConc(:)  = sedimentTracers(iClass, :, iCell)
            workTend(:)  = 0.0_RKIND

            select case (iClass)
               case (1)
                  call ocn_sediment_settling_flx(                        &
                     workConc,  sinkingRate(iClass),  workTend,          &
                     layerThickness(:, iCell),  dt,  kmin,  kmax,        &
                     sedimentSettlingFlux1(iCell),                        &
                     sedimentBedMass1(iCell),                             &
                     bottomStressCellSediment(iCell),                     &
                     sedimentErosionFlux1(iCell))
               case (2)
                  call ocn_sediment_settling_flx(                        &
                     workConc,  sinkingRate(iClass),  workTend,          &
                     layerThickness(:, iCell),  dt,  kmin,  kmax,        &
                     sedimentSettlingFlux2(iCell),                        &
                     sedimentBedMass2(iCell),                             &
                     bottomStressCellSediment(iCell),                     &
                     sedimentErosionFlux2(iCell))
               case (3)
                  call ocn_sediment_settling_flx(                        &
                     workConc,  sinkingRate(iClass),  workTend,          &
                     layerThickness(:, iCell),  dt,  kmin,  kmax,        &
                     sedimentSettlingFlux3(iCell),                        &
                     sedimentBedMass3(iCell),                             &
                     bottomStressCellSediment(iCell),                     &
                     sedimentErosionFlux3(iCell))
            end select

            sedimentTend(iClass, :, iCell) = sedimentTend(iClass, :, iCell) + workTend(:)

         end do ! iClass

      end do ! iCell

      deallocate(workTend, workConc)

   end subroutine ocn_tracer_sediment_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_sediment_surface_flux_compute
!
!> \brief   Computes sediment surface fluxes (river + ice coupling)
!> \author  Nicole Jeffery
!> \date    07/2026
!> \details
!>  Phase 2. Scale river runoff flux with fixed sediment load
!>  Phase 3.  River runoff and ice-ocean sediment fluxes
!>  will be computed here and accumulated into sedimentTracersSurfaceFlux.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_sediment_surface_flux_compute(sedimentTracers, forcingPool,  &
      nCellsSolve, sedimentSurfaceFlux, err)!{{{

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         sedimentTracers           !< (nClass, nVertLevels, nCells)

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         sedimentSurfaceFlux       !< (nClass, nCells)

      ! scalars
      integer, intent(in) :: nCellsSolve

      type (mpas_pool_type), intent(inout) :: forcingPool

      integer, intent(out) :: err                 !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      ! input flux components in forcing pool
      real (kind=RKIND), dimension(:), pointer :: &
         riverRunoffFlux

      integer :: iCell, iClass, column

      err = 0

      ! Zero surface fluxes; Phase 3 will add river and ice-ocean contributions
      sedimentSurfaceFlux(:,:) = 0.0_RKIND

      if (config_use_sediment_river_flux) then
         call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
         do iCell=1,nCellsSolve
            do iClass = 1, nSedimentClasses
               sedimentSurfaceFlux(iClass,iCell) = sedimentSurfaceFlux(iClass,iCell) + &
               max(0.0_RKIND,riverRunoffFlux(iCell))/rho_fw * riverSedimentLoad
            end do
         end do
      end if
      if (config_use_sediment_sea_ice_coupling) then ! placeholder for Phase 3
         do iCell=1,nCellsSolve
            do iClass = 1, nSedimentClasses
               sedimentSurfaceFlux(iClass,iCell) = sedimentSurfaceFlux(iClass,iCell) + 0.0_RKIND
            end do
         end do

      end if

   end subroutine ocn_tracer_sediment_surface_flux_compute!}}}

!***********************************************************************
end module ocn_tracer_sediment

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
