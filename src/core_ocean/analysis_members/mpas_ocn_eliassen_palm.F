! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_eliassen_palm
!
!> \brief MPAS ocean analysis core member: Eliassen-Palm Flux Tensor
!> \author Juan A. Saenz, Todd Ringler
!> \date   May 2015
!> \details
!>  This module contains the routines for computing the Eliassen and Palm Flux Tensor
!>  in buoyancy coordinates, and related quantities.
!
!-----------------------------------------------------------------------

module ocn_eliassen_palm

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use mpas_pool_routines
   use mpas_constants
   use ocn_constants
   use ocn_diagnostics_routines
   use ocn_equation_of_state

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

   public :: ocn_init_eliassen_palm, &
             ocn_compute_eliassen_palm, &
             ocn_restart_eliassen_palm, &
             ocn_finalize_eliassen_palm

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND), parameter :: epsilonEPFT=1.0e-15_RKIND

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_eliassen_palm
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_eliassen_palm(domain, err)!{{{

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
      integer :: err_tmp
      integer :: k
      real (KIND=RKIND) :: deltaBuoyancy, deltaDensity

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: amEPFTPool

      integer :: nBuoyancyLayers
      real (kind=RKIND), dimension(:), pointer :: potentialDensityMidRef
      real (kind=RKIND), dimension(:), pointer :: potentialDensityTopRef
      real (kind=RKIND), dimension(:), pointer :: buoyancyMidRef
      real (kind=RKIND), dimension(:), pointer :: buoyancyInterfaceRef

      logical, pointer :: amEPFTActive
      logical, pointer :: config_AM_eliassenPalm_compute_on_startup
      integer, pointer :: config_AM_eliassenPalm_nBuoyancyLayers
      real (kind=RKIND), pointer :: config_AM_eliassenPalm_rhomax_buoycoor
      real (kind=RKIND), pointer :: config_AM_eliassenPalm_rhomin_buoycoor

      integer, pointer :: nSamplesEA

      real (kind=RKIND), dimension(:,:), pointer :: buoyancyMaskEA
      real (kind=RKIND), dimension(:,:), pointer :: sigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: heightMidBuoyCoorEA
      real (kind=RKIND), dimension(:,:), pointer :: montgPotBuoyCoorEA
      real (kind=RKIND), dimension(:,:), pointer :: montgPotGradZonalEA
      real (kind=RKIND), dimension(:,:), pointer :: montgPotGradMeridEA
      real (kind=RKIND), dimension(:,:), pointer :: heightMidBuoyCoorSqEA
      real (kind=RKIND), dimension(:,:), pointer :: heightMGradZonalEA
      real (kind=RKIND), dimension(:,:), pointer :: heightMGradMeridEA
      real (kind=RKIND), dimension(:,:), pointer :: usigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: vsigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: varpisigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: uusigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: vvsigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: uvsigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: uvarpisigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: vvarpisigmaEA

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_AM_eliassenPalm_compute_on_startup', &
        config_AM_eliassenPalm_compute_on_startup)
      call mpas_pool_get_config(domain % configs, 'config_AM_eliassenPalm_nBuoyancyLayers', &
        config_AM_eliassenPalm_nBuoyancyLayers)
      call mpas_pool_get_config(domain % configs, 'config_AM_eliassenPalm_rhomax_buoycoor', &
        config_AM_eliassenPalm_rhomax_buoycoor)
      call mpas_pool_get_config(domain % configs, 'config_AM_eliassenPalm_rhomin_buoycoor', &
        config_AM_eliassenPalm_rhomin_buoycoor)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'eliassenPalmAM', amEPFTPool)

         !-----------------------------------------------------------------
         ! set up pointers
         !-----------------------------------------------------------------
         call mpas_pool_get_array(amEPFTPool, 'potentialDensityMidRef', potentialDensityMidRef)
         call mpas_pool_get_array(amEPFTPool, 'potentialDensityTopRef', potentialDensityTopRef)
         call mpas_pool_get_array(amEPFTPool, 'buoyancyMidRef', buoyancyMidRef)
         call mpas_pool_get_array(amEPFTPool, 'buoyancyInterfaceRef', buoyancyInterfaceRef)

         !-----------------------------------------------------------------
         ! compute buoyancy and density increment of each layer
         ! at present we use layer interfaces that are evenly-spaced in buoyancy space
         !-----------------------------------------------------------------
         nBuoyancyLayers = config_AM_eliassenPalm_nBuoyancyLayers
         deltaDensity = (config_AM_eliassenPalm_rhomax_buoycoor &
            - config_AM_eliassenPalm_rhomin_buoycoor) / config_AM_eliassenPalm_nBuoyancyLayers
         deltaBuoyancy = -gravity * deltaDensity / rho_sw

         !-----------------------------------------------------------------
         ! compute density/bouyancy at top of each layer
         !-----------------------------------------------------------------
         do k = 1, nBuoyancyLayers
            potentialDensityTopRef(k) = config_AM_eliassenPalm_rhomin_buoycoor + deltaDensity * (k-1)
            buoyancyInterfaceRef(k) = -gravity &
              * (config_AM_eliassenPalm_rhomin_buoycoor - rho_sw) / rho_sw &
              + deltaBuoyancy * (k-1)
         end do
         k=nBuoyancyLayers
         buoyancyInterfaceRef(k+1) = buoyancyInterfaceRef(k) + deltaBuoyancy

         !-----------------------------------------------------------------
         ! compute density/bouyancy for each layer
         !-----------------------------------------------------------------
         do k = 1, nBuoyancyLayers-1
            potentialDensityMidRef(k) = 0.5_RKIND*(potentialDensityTopRef(k) + potentialDensityTopRef(k+1))
            buoyancyMidRef(k) = 0.5_RKIND*(buoyancyInterfaceRef(k) + buoyancyInterfaceRef(k+1))
         end do
         k=nBuoyancyLayers
         potentialDensityMidRef(k) = 0.5*(potentialDensityTopRef(k) + config_AM_eliassenPalm_rhomax_buoycoor)
         buoyancyMidRef(k) = 0.5_RKIND*(buoyancyInterfaceRef(k) + buoyancyInterfaceRef(k+1))

         block => block % next

      end do


   end subroutine ocn_init_eliassen_palm!}}}

!***********************************************************************
!
!  routine ocn_compute_eliassen_palm
!
!> \brief   Compute Eliassen-Palm flux tensor
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details
!>  This routine conducts all computations required for the EPFT analysis member.
!>  Each time this AM is called, the instananeous ocean state is interpolated
!>  onto the target buoyancy values. The state is then accumulated into the
!>  ensemble average variable arrays (varEA). Based on the current
!>  estimate of the ensemble average, thickness-weight velocities are estimated
!>  along with the computation of the Eliassen-Palm flux tensor.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_eliassen_palm(domain, timeLevel, err)!{{{

      use mpas_vector_reconstruction

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

      !-----------------------------------------------------------------
      ! define types that live inside of domain
      !-----------------------------------------------------------------
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: am_epftPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      !-----------------------------------------------------------------
      ! define pointers to namelist config variables local to the EPFT module
      !-----------------------------------------------------------------
      logical, pointer :: config_AM_eliassenPalm_debug
      real (kind=RKIND), pointer :: config_AM_eliassenPalm_rhomin_buoycoor
      real (kind=RKIND), pointer :: config_AM_eliassenPalm_rhomax_buoycoor

      !-----------------------------------------------------------------
      ! define local scalars holding length of dimensions
      !-----------------------------------------------------------------
      integer, pointer :: nVertLevels, nBuoyancyLayers, nBuoyancyLayersP1
      integer, pointer :: nEdges, nCells, nCellsSolve ! nCellsSolve does not include halo

      !-----------------------------------------------------------------
      ! define buoyancy coordinates and fields related to the vertical direction
      !-----------------------------------------------------------------
      integer, dimension(:), pointer :: maxLevelCell
      real(KIND=RKIND), dimension(:), pointer :: potentialDensityMidRef
      real(KIND=RKIND), dimension(:), pointer :: potentialDensityTopRef
      real(KIND=RKIND), dimension(:), pointer :: buoyancyMidRef
      real(KIND=RKIND), dimension(:), pointer :: buoyancyInterfaceRef
      real(KIND=RKIND), dimension(:), pointer :: bottomDepth

      !-----------------------------------------------------------------
      ! define mesh variables
      !-----------------------------------------------------------------
      real(KIND=RKIND), dimension(:), pointer :: fCell
      integer, dimension(:,:), pointer :: cellMask

      !-----------------------------------------------------------------
      ! define fields related to the Ensemble Average (EA)
      !-----------------------------------------------------------------
      integer, pointer :: nSamplesEA
      real(KIND=RKIND), dimension(:,:), pointer :: buoyancyMaskEA
      real(KIND=RKIND), dimension(:,:), pointer :: sigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: heightMidBuoyCoorEA
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradZonalEA
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradMeridEA
      real(KIND=RKIND), dimension(:,:), pointer :: heightMidBuoyCoorSqEA
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotBuoyCoorEA
      real(KIND=RKIND), dimension(:,:), pointer :: heightMGradZonalEA
      real(KIND=RKIND), dimension(:,:), pointer :: heightMGradMeridEA
      real(KIND=RKIND), dimension(:,:), pointer :: usigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: vsigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: varpisigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: uusigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: vvsigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: uvsigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: uvarpisigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: vvarpisigmaEA

      !-----------------------------------------------------------------
      ! define the Thickness-Weighted Average (TWA) velocity
      !-----------------------------------------------------------------
      real(KIND=RKIND), dimension(:,:), pointer :: uTWA
      real(KIND=RKIND), dimension(:,:), pointer :: vTWA
      real(KIND=RKIND), dimension(:,:), pointer :: varpiTWA
      real(KIND=RKIND), dimension(:,:), pointer :: duTWAdz
      real(KIND=RKIND), dimension(:,:), pointer :: dvTWAdz

      !-----------------------------------------------------------------
      ! define Ertel's potential vorticity and related fields
      !-----------------------------------------------------------------
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPV
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPVGradZonal
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPVGradMerid
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPVTendency
      real(KIND=RKIND), dimension(:,:,:), pointer :: ErtelPVFlux
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPVFlux1
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPVFlux2

      !-----------------------------------------------------------------
      ! define the Eliassen-Palm flux tensor and related fields
      !-----------------------------------------------------------------
      real(KIND=RKIND), dimension(:,:,:,:), pointer :: EPFT
      real(KIND=RKIND), dimension(:,:,:), pointer :: divEPFT
      real(KIND=RKIND), dimension(:,:), pointer :: divEPFT1
      real(KIND=RKIND), dimension(:,:), pointer :: divEPFT2
      real(KIND=RKIND), dimension(:,:), pointer :: divEPFTshear1
      real(KIND=RKIND), dimension(:,:), pointer :: divEPFTshear2
      real(KIND=RKIND), dimension(:,:), pointer :: divEPFTdrag1
      real(KIND=RKIND), dimension(:,:), pointer :: divEPFTdrag2
      real(KIND=RKIND), dimension(:,:), pointer :: uuTWACorr
      real(KIND=RKIND), dimension(:,:), pointer :: vvTWACorr
      real(KIND=RKIND), dimension(:,:), pointer :: uvTWACorr
      real(KIND=RKIND), dimension(:,:), pointer :: epeTWA
      real(KIND=RKIND), dimension(:,:), pointer :: eddyFormDragZonal
      real(KIND=RKIND), dimension(:,:), pointer :: eddyFormDragMerid

      !-----------------------------------------------------------------
      ! define scratch fields used as work variables and for testing
      !-----------------------------------------------------------------
      type(field1DInteger), pointer :: firstLayerBuoyCoorField
      type(field1DInteger), pointer :: lastLayerBuoyCoorField
      type(field2DReal), pointer :: heightMidBuoyCoorField
      type(field2DReal), pointer :: heightTopBuoyCoorField
      type(field2DReal), pointer :: heightInterfaceBuoyCoorField
      type(field2DReal), pointer :: sigmaField
      type(field2DReal), pointer :: montgPotBuoyCoorField
      type(field2DReal), pointer :: montgPotNormalGradOnEdgeField
      type(field2DReal), pointer :: uMidBuoyCoorField
      type(field2DReal), pointer :: vMidBuoyCoorField
      type(field2DReal), pointer :: densityMidBuoyCoorField
      type(field2DReal), pointer :: densityTopBuoyCoorField
      type(field2DReal), pointer :: buoyancyMaskField
      type(field2DReal), pointer :: montgPotGradXField
      type(field2DReal), pointer :: montgPotGradYField
      type(field2DReal), pointer :: montgPotGradZField
      type(field2DReal), pointer :: montgPotGradZonalField
      type(field2DReal), pointer :: montgPotGradMeridField
      type(field2DReal), pointer :: wrk3DnVertLevelsP1Field
      type(field2DReal), pointer :: wrk3DnVertLevelsField
      type(field2DReal), pointer :: wrk3DBuoyCoorField
      type(field2DReal), pointer :: ErtelPVNormalGradOnEdgeField
      type(field2DReal), pointer :: ErtelPVGradXField
      type(field2DReal), pointer :: ErtelPVGradYField
      type(field2DReal), pointer :: ErtelPVGradZField
      type(field3DReal), pointer :: wrkVectorField
      type(field4DReal), pointer :: wrkTensorField

      type(field2DReal), pointer :: array1_3DField
      type(field2DReal), pointer :: array2_3DField
      type(field2DReal), pointer :: array3_3DField
      type(field2DReal), pointer :: array1_3DbuoyField
      type(field2DReal), pointer :: array2_3DbuoyField
      type(field2DReal), pointer :: PVMidBuoyCoorField
      type(field2DReal), pointer :: PVMidBuoyCoorEAField
      type(field2DReal), pointer :: uMidBuoyCoorEAField
      type(field2DReal), pointer :: vMidBuoyCoorEAField
      type(field2DReal), pointer :: uPVMidBuoyCoorEAField
      type(field2DReal), pointer :: vPVMidBuoyCoorEAField
      type(field3DReal), pointer :: PVFluxTestField

      !-----------------------------------------------------------------
      ! define pointers to scratch fields
      !-----------------------------------------------------------------
      integer, dimension(:), pointer :: firstLayerBuoyCoor
      integer, dimension(:), pointer :: lastLayerBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: heightMidBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: heightTopBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: heightInterfaceBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: sigma
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotNormalGradOnEdge
      real(KIND=RKIND), dimension(:,:), pointer :: uMidBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: vMidBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: densityMidBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: densityTopBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: buoyancyMask
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradX
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradY
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradZ
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradZonal
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradMerid
      real(KIND=RKIND), dimension(:,:), pointer :: wrk3DnVertLevelsP1
      real(KIND=RKIND), dimension(:,:), pointer :: wrk3DnVertLevels
      real(KIND=RKIND), dimension(:,:), pointer :: wrk3DBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPVNormalGradOnEdge
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPVGradX
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPVGradY
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPVGradZ
      real(KIND=RKIND), dimension(:,:,:), pointer :: wrkVector
      real(KIND=RKIND), dimension(:,:,:,:), pointer :: wrkTensor

      real(KIND=RKIND), dimension(:,:), pointer :: array1_3D
      real(KIND=RKIND), dimension(:,:), pointer :: array2_3D
      real(KIND=RKIND), dimension(:,:), pointer :: array3_3D
      real(KIND=RKIND), dimension(:,:), pointer :: array1_3Dbuoy
      real(KIND=RKIND), dimension(:,:), pointer :: array2_3Dbuoy
      real(KIND=RKIND), dimension(:,:), pointer :: PVMidBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: PVMidBuoyCoorEA
      real(KIND=RKIND), dimension(:,:), pointer :: uMidBuoyCoorEA
      real(KIND=RKIND), dimension(:,:), pointer :: vMidBuoyCoorEA
      real(KIND=RKIND), dimension(:,:), pointer :: uPVMidBuoyCoorEA
      real(KIND=RKIND), dimension(:,:), pointer :: vPVMidBuoyCoorEA
      real(KIND=RKIND), dimension(:,:,:), pointer :: PVFluxTest

      !-----------------------------------------------------------------
      ! define some arrays in z-coordinates, obtained from diagnostics and forcing
      !-----------------------------------------------------------------
      real(KIND=RKIND), dimension(:), pointer :: atmosphericPressure
      real(KIND=RKIND), dimension(:,:), pointer :: zMid
      real(KIND=RKIND), dimension(:,:), pointer :: zTop
      real(KIND=RKIND), dimension(:,:), pointer :: density
      real(KIND=RKIND), dimension(:,:), pointer :: potentialDensity
      real(KIND=RKIND), dimension(:,:), pointer :: pressure
      real(KIND=RKIND), dimension(:,:), pointer :: velocityZonal
      real(KIND=RKIND), dimension(:,:), pointer :: velocityMeridional
      real(KIND=RKIND), dimension(:,:), pointer :: relativeVorticityCell ! jas used for testing
      !real(KIND=RKIND), dimension(:,:), pointer :: wCellCenter

      !-----------------------------------------------------------------
      ! define local test variables
      !-----------------------------------------------------------------
      ! jas to do : move these to scratch in Registry_epft
      integer :: nCellsCum
      real(KIND=RKIND) :: RMSlocal1, RMSglobal1
      real(KIND=RKIND) :: RMSlocal2, RMSglobal2
      real(KIND=RKIND) :: RMSPVFlux1local, RMSPVFlux1global
      real(KIND=RKIND) :: RMSPVFlux2local, RMSPVFlux2global

      !-----------------------------------------------------------------
      ! define local work variables
      !-----------------------------------------------------------------
      integer :: nCellsGlobal, k, i

      err = 0

      nCellsCum = 0
      RMSlocal1  = 0.0_RKIND
      RMSlocal2  = 0.0_RKIND
      RMSglobal1 = 0.0_RKIND
      RMSglobal2 = 0.0_RKIND
      RMSPVFlux1local  = 0.0_RKIND
      RMSPVFlux2local  = 0.0_RKIND
      RMSPVFlux1global = 0.0_RKIND
      RMSPVFlux2global = 0.0_RKIND

      dminfo = domain % dminfo

      !--------------------------------------------------
      ! get config variables
      !--------------------------------------------------
      call mpas_pool_get_config(domain % configs, 'config_AM_eliassenPalm_debug', &
        config_AM_eliassenPalm_debug)
      call mpas_pool_get_config(domain % configs, 'config_AM_eliassenPalm_rhomin_buoycoor', &
        config_AM_eliassenPalm_rhomin_buoycoor)
      call mpas_pool_get_config(domain % configs, 'config_AM_eliassenPalm_rhomax_buoycoor', &
        config_AM_eliassenPalm_rhomax_buoycoor)

      if(config_AM_eliassenPalm_debug) then
         write(stderrUnit, *) ' '
         write(stderrUnit, *) 'starting ocn_compute_epft'
         write(stderrUnit, *) ' '
      end if

      block => domain % blocklist
      do while (associated(block))

         !--------------------------------------------------
         ! assign pointers for each pool
         !--------------------------------------------------
         call mpas_pool_get_subpool(block % structs, 'eliassenPalmAM', am_epftPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'eliassenPalmAMPKGScratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         !--------------------------------------------------
         ! assign pointers for mesh-related variables
         !--------------------------------------------------
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
         call mpas_pool_get_array(meshPool, 'fCell', fCell)
         call mpas_pool_get_array(meshPool, 'cellMask', cellMask) ! used for tests


         !--------------------------------------------------
         ! get scratch field pointers
         !--------------------------------------------------
         call mpas_pool_get_field(scratchPool, 'firstLayerBuoyCoor', firstLayerBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'lastLayerBuoyCoor', lastLayerBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'heightMidBuoyCoor', heightMidBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'heightTopBuoyCoor', heightTopBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'heightInterfaceBuoyCoor', heightInterfaceBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'sigma', sigmaField)
         call mpas_pool_get_field(scratchPool, 'montgPotBuoyCoor', montgPotBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'montgPotNormalGradOnEdge', montgPotNormalGradOnEdgeField)
         call mpas_pool_get_field(scratchPool, 'uMidBuoyCoor', uMidBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'vMidBuoyCoor', vMidBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'densityMidBuoyCoor', densityMidBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'densityTopBuoyCoor', densityTopBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'buoyancyMask', buoyancyMaskField)
         call mpas_pool_get_field(scratchPool, 'montgPotGradX', montgPotGradXField)
         call mpas_pool_get_field(scratchPool, 'montgPotGradY', montgPotGradYField)
         call mpas_pool_get_field(scratchPool, 'montgPotGradZ', montgPotGradZField)
         call mpas_pool_get_field(scratchPool, 'montgPotGradZonal', montgPotGradZonalField)
         call mpas_pool_get_field(scratchPool, 'montgPotGradMerid', montgPotGradMeridField)
         call mpas_pool_get_field(scratchPool, 'wrk3DnVertLevelsP1', wrk3DnVertLevelsP1Field)
         call mpas_pool_get_field(scratchPool, 'wrk3DnVertLevels', wrk3DnVertLevelsField)
         call mpas_pool_get_field(scratchPool, 'wrk3DBuoyCoor', wrk3DBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'ErtelPVNormalGradOnEdge', ErtelPVNormalGradOnEdgeField)
         call mpas_pool_get_field(scratchPool, 'ErtelPVGradX', ErtelPVGradXField)
         call mpas_pool_get_field(scratchPool, 'ErtelPVGradY', ErtelPVGradYField)
         call mpas_pool_get_field(scratchPool, 'ErtelPVGradZ', ErtelPVGradZField)
         call mpas_pool_get_field(scratchPool, 'wrkVector', wrkVectorField)
         call mpas_pool_get_field(scratchPool, 'wrkTensor', wrkTensorField)

         call mpas_pool_get_field(scratchPool, 'array1_3D', array1_3DField)
         call mpas_pool_get_field(scratchPool, 'array2_3D', array2_3DField)
         call mpas_pool_get_field(scratchPool, 'array3_3D', array3_3DField)
         call mpas_pool_get_field(scratchPool, 'array1_3Dbuoy', array1_3DbuoyField)
         call mpas_pool_get_field(scratchPool, 'array2_3Dbuoy', array2_3DbuoyField)
         call mpas_pool_get_field(scratchPool, 'PVMidBuoyCoor', PVMidBuoyCoorField)
         call mpas_pool_get_field(scratchPool, 'PVMidBuoyCoorEA', PVMidBuoyCoorEAField)
         call mpas_pool_get_field(scratchPool, 'uMidBuoyCoorEA', uMidBuoyCoorEAField)
         call mpas_pool_get_field(scratchPool, 'vMidBuoyCoorEA', vMidBuoyCoorEAField)
         call mpas_pool_get_field(scratchPool, 'uPVMidBuoyCoorEA', uPVMidBuoyCoorEAField)
         call mpas_pool_get_field(scratchPool, 'vPVMidBuoyCoorEA', vPVMidBuoyCoorEAField)
         call mpas_pool_get_field(scratchPool, 'PVFluxTest', PVFluxTestField)

         !--------------------------------------------------
         ! allocate scratch field variables
         !--------------------------------------------------
         call mpas_allocate_scratch_field(firstLayerBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(lastLayerBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(heightMidBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(heightTopBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(heightInterfaceBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(sigmaField, .true.)
         call mpas_allocate_scratch_field(montgPotBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(montgPotNormalGradOnEdgeField, .true.)
         call mpas_allocate_scratch_field(uMidBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(vMidBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(densityMidBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(densityTopBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(buoyancyMaskField, .true.)
         call mpas_allocate_scratch_field(montgPotGradXField, .true.)
         call mpas_allocate_scratch_field(montgPotGradYField, .true.)
         call mpas_allocate_scratch_field(montgPotGradZField, .true.)
         call mpas_allocate_scratch_field(montgPotGradZonalField, .true.)
         call mpas_allocate_scratch_field(montgPotGradMeridField, .true.)
         call mpas_allocate_scratch_field(wrk3DnVertLevelsP1Field, .true.)
         call mpas_allocate_scratch_field(wrk3DnVertLevelsField, .true.)
         call mpas_allocate_scratch_field(wrk3DBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(ErtelPVNormalGradOnEdgeField, .true.)
         call mpas_allocate_scratch_field(ErtelPVGradXField, .true.)
         call mpas_allocate_scratch_field(ErtelPVGradYField, .true.)
         call mpas_allocate_scratch_field(ErtelPVGradZField, .true.)
         call mpas_allocate_scratch_field(wrkVectorField, .true.)
         call mpas_allocate_scratch_field(wrkTensorField, .true.)

         call mpas_allocate_scratch_field(array1_3DField, .true.)
         call mpas_allocate_scratch_field(array2_3DField, .true.)
         call mpas_allocate_scratch_field(array3_3DField, .true.)
         call mpas_allocate_scratch_field(array1_3DbuoyField, .true.)
         call mpas_allocate_scratch_field(array2_3DbuoyField, .true.)
         call mpas_allocate_scratch_field(PVMidBuoyCoorField, .true.)
         call mpas_allocate_scratch_field(PVMidBuoyCoorEAField, .true.)
         call mpas_allocate_scratch_field(uMidBuoyCoorEAField, .true.)
         call mpas_allocate_scratch_field(vMidBuoyCoorEAField, .true.)
         call mpas_allocate_scratch_field(uPVMidBuoyCoorEAField, .true.)
         call mpas_allocate_scratch_field(vPVMidBuoyCoorEAField, .true.)
         call mpas_allocate_scratch_field(PVFluxTestField, .true.)

         !--------------------------------------------------
         ! assign pointers for scratch and test variables
         !--------------------------------------------------
         firstLayerBuoyCoor      => firstLayerBuoyCoorField % array
         lastLayerBuoyCoor       => lastLayerBuoyCoorField % array
         heightMidBuoyCoor       => heightMidBuoyCoorField % array
         heightTopBuoyCoor       => heightTopBuoyCoorField % array
         heightInterfaceBuoyCoor => heightInterfaceBuoyCoorField % array
         sigma                   => sigmaField % array
         montgPotBuoyCoor        => montgPotBuoyCoorField % array
         montgPotNormalGradOnEdge=> montgPotNormalGradOnEdgeField % array
         uMidBuoyCoor            => uMidBuoyCoorField % array
         vMidBuoyCoor            => vMidBuoyCoorField % array
         densityMidBuoyCoor      => densityMidBuoyCoorField % array
         densityTopBuoyCoor      => densityTopBuoyCoorField % array
         buoyancyMask            => buoyancyMaskField % array
         montgPotGradX           => montgPotGradXField % array
         montgPotGradY           => montgPotGradYField % array
         montgPotGradZ           => montgPotGradZField % array
         montgPotGradZonal       => montgPotGradZonalField % array
         montgPotGradMerid       => montgPotGradMeridField % array
         wrk3DnVertLevelsP1      => wrk3DnVertLevelsP1Field % array
         wrk3DnVertLevels        => wrk3DnVertLevelsField % array
         wrk3DBuoyCoor           => wrk3DBuoyCoorField % array
         ErtelPVNormalGradOnEdge => ErtelPVNormalGradOnEdgeField % array
         ErtelPVGradX            => ErtelPVGradXField % array
         ErtelPVGradY            => ErtelPVGradYField % array
         ErtelPVGradZ            => ErtelPVGradZField % array
         wrkVector               => wrkVectorField % array
         wrkTensor               => wrkTensorField % array

         array1_3D =>  array1_3DField % array
         array2_3D =>  array2_3DField % array
         array3_3D =>  array3_3DField % array
         array1_3Dbuoy =>  array1_3DbuoyField % array
         array2_3Dbuoy =>  array2_3DbuoyField % array
         PVMidBuoyCoor    =>  PVMidBuoyCoorField % array
         PVMidBuoyCoorEA  =>  PVMidBuoyCoorEAField % array
         uMidBuoyCoorEA   =>  uMidBuoyCoorEAField % array
         vMidBuoyCoorEA   =>  vMidBuoyCoorEAField % array
         uPVMidBuoyCoorEA =>  uPVMidBuoyCoorEAField % array
         vPVMidBuoyCoorEA =>  vPVMidBuoyCoorEAField % array
         PVFluxTest       =>  PVFluxTestField % array

         !--------------------------------------------------
         ! assign pointers used from forcing pool
         !--------------------------------------------------
         call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)

         !--------------------------------------------------
         ! get diagnostic variables
         !--------------------------------------------------
         call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
         call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
         call mpas_pool_get_array(diagnosticsPool, 'density', density)
         call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', potentialDensity)
         call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)

         !--------------------------------------------------
         ! variables that define the vertical coordinate system in density/buoyancy space
         !--------------------------------------------------
         call mpas_pool_get_dimension(block % dimensions, 'nBuoyancyLayers', nBuoyancyLayers)
         call mpas_pool_get_dimension(block % dimensions, 'nBuoyancyLayersP1', nBuoyancyLayersP1)
         call mpas_pool_get_array(am_epftPool, 'potentialDensityMidRef', potentialDensityMidRef)
         call mpas_pool_get_array(am_epftPool, 'potentialDensityTopRef', potentialDensityTopRef)
         call mpas_pool_get_array(am_epftPool, 'buoyancyMidRef', buoyancyMidRef)
         call mpas_pool_get_array(am_epftPool, 'buoyancyInterfaceRef', buoyancyInterfaceRef)

         !--------------------------------------------------
         ! assign pointers for ensemble average (EA) and thickness-weighted averaged state
         !--------------------------------------------------
         call mpas_pool_get_array(am_epftPool, 'nSamplesEA', nSamplesEA)
         call mpas_pool_get_array(am_epftPool, 'buoyancyMaskEA', buoyancyMaskEA)
         call mpas_pool_get_array(am_epftPool, 'sigmaEA', sigmaEA)
         call mpas_pool_get_array(am_epftPool, 'heightMidBuoyCoorEA', heightMidBuoyCoorEA)
         call mpas_pool_get_array(am_epftPool, 'heightMidBuoyCoorSqEA', heightMidBuoyCoorSqEA)
         call mpas_pool_get_array(am_epftPool, 'montgPotBuoyCoorEA', montgPotBuoyCoorEA)
         call mpas_pool_get_array(am_epftPool, 'montgPotGradZonalEA', montgPotGradZonalEA)
         call mpas_pool_get_array(am_epftPool, 'montgPotGradMeridEA', montgPotGradMeridEA)
         call mpas_pool_get_array(am_epftPool, 'heightMGradZonalEA', HeightMGradZonalEA)
         call mpas_pool_get_array(am_epftPool, 'heightMGradMeridEA', HeightMGradMeridEA)
         call mpas_pool_get_array(am_epftPool, 'usigmaEA', usigmaEA)
         call mpas_pool_get_array(am_epftPool, 'vsigmaEA', vsigmaEA)
         call mpas_pool_get_array(am_epftPool, 'varpisigmaEA', varpisigmaEA)
         call mpas_pool_get_array(am_epftPool, 'uusigmaEA', uusigmaEA)
         call mpas_pool_get_array(am_epftPool, 'vvsigmaEA', vvsigmaEA)
         call mpas_pool_get_array(am_epftPool, 'uvsigmaEA', uvsigmaEA)
         call mpas_pool_get_array(am_epftPool, 'uvarpisigmaEA', uvarpisigmaEA)
         call mpas_pool_get_array(am_epftPool, 'vvarpisigmaEA', vvarpisigmaEA)

         !--------------------------------------------------
         ! assign pointers for thickness-weighted averaged state
         !--------------------------------------------------
         call mpas_pool_get_array(am_epftPool, 'uTWA', uTWA)
         call mpas_pool_get_array(am_epftPool, 'vTWA', vTWA)
         call mpas_pool_get_array(am_epftPool, 'varpiTWA', varpiTWA)
         call mpas_pool_get_array(am_epftPool, 'duTWAdz', duTWAdz)
         call mpas_pool_get_array(am_epftPool, 'dvTWAdz', dvTWAdz)

         !--------------------------------------------------
         ! Eliassen-Palm Flux Tensor and related products
         !--------------------------------------------------
         call mpas_pool_get_array(am_epftPool, 'EPFT', EPFT)
         call mpas_pool_get_array(am_epftPool, 'divEPFT', divEPFT)
         call mpas_pool_get_array(am_epftPool, 'divEPFT1', divEPFT1)
         call mpas_pool_get_array(am_epftPool, 'divEPFT2', divEPFT2)
         call mpas_pool_get_array(am_epftPool, 'divEPFTshear1', divEPFTshear1)
         call mpas_pool_get_array(am_epftPool, 'divEPFTshear2', divEPFTshear2)
         call mpas_pool_get_array(am_epftPool, 'divEPFTdrag1', divEPFTdrag1)
         call mpas_pool_get_array(am_epftPool, 'divEPFTdrag2', divEPFTdrag2)
         call mpas_pool_get_array(am_epftPool, 'uuTWACorr', uuTWACorr)
         call mpas_pool_get_array(am_epftPool, 'vvTWACorr', vvTWACorr)
         call mpas_pool_get_array(am_epftPool, 'uvTWACorr', uvTWACorr)
         call mpas_pool_get_array(am_epftPool, 'epeTWA', epeTWA)
         call mpas_pool_get_array(am_epftPool, 'eddyFormDragZonal', eddyFormDragZonal)
         call mpas_pool_get_array(am_epftPool, 'eddyFormDragMerid', eddyFormDragMerid)

         call mpas_pool_get_array(am_epftPool, 'ErtelPVFlux' , ErtelPVFlux)
         call mpas_pool_get_array(am_epftPool, 'ErtelPVFlux1', ErtelPVFlux1)
         call mpas_pool_get_array(am_epftPool, 'ErtelPVFlux2', ErtelPVFlux2)
         call mpas_pool_get_array(am_epftPool, 'ErtelPVTendency', ErtelPVTendency)
         call mpas_pool_get_array(am_epftPool, 'ErtelPV', ErtelPV)
         call mpas_pool_get_array(am_epftPool, 'ErtelPVGradZonal', ErtelPVGradZonal)
         call mpas_pool_get_array(am_epftPool, 'ErtelPVGradMerid', ErtelPVGradMerid)

         ! Compute potentialDensity over the entire block to ensure it's valid for EPFT computation
         call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, scratchPool, &
                                            nCells, 1, 'absolute', potentialDensity, &
                                            err) !, timeLevelIn=1)

         !--------------------------------------------------
         ! Get variables associated to diabatic processes
         !--------------------------------------------------
         !jas issue diabatic terms
         ! diabaticHeating(nVertLevels,nCells)! "vertical velocity" in buoyancy space
         !wCellCenter       = 0.0
         ! Get diabaticTimeTendency of a buoyancy surface, omega with funny hat, if any.
         !call any existing MPAS-O subroutines for this


         !-------------------------------------------------------------
         ! begin computation
         !-------------------------------------------------------------

         !-------------------------------------------------------------
         ! compute firstLayerBuoyCoor and lastLayerBuoyCoor
         !    firstLayerBuoyCoor == top buoyancy coordinate to exist in each column
         !    lastLayerBuoyCoor == bottom buoyancy coordinate to exist in each column
         !    buoyancyMask == 1.0 between layers firstLayerBuoyCoor and lastLayerBuoyCoor
         !-------------------------------------------------------------
         call get_masks_in_buoyancy_coordinates(nVertLevels, nCells, nBuoyancyLayers, &
           maxLevelCell, potentialDensity, potentialDensityMidRef, &
           firstLayerBuoyCoor, lastLayerBuoyCoor, buoyancyMask)

         !-------------------------------------------------------------
         ! TEST for general consistency
         !-------------------------------------------------------------
         if(config_AM_eliassenPalm_debug) then
            print *, ' '
            print *, 'timeLevel:', timeLevel
            print *, 'nBuoyancyLayers:', timeLevel
            print *, ' '
            print *, 'config_AM_eliassenPalm_rhomin_buoycoor, config_AM_eliassenPalm_rhomax_buoycoor'
            print *, config_AM_eliassenPalm_rhomin_buoycoor, config_AM_eliassenPalm_rhomax_buoycoor
            print *, '(config_AM_eliassenPalm_rhomax_buoycoor - config_AM_eliassenPalm_rhomin_buoycoor)/nBuoyancyLayers'
            print *, (config_AM_eliassenPalm_rhomax_buoycoor - config_AM_eliassenPalm_rhomin_buoycoor)/nBuoyancyLayers
            print *, 'potentialDensityTopRef (nBuoyancyLayers)'
            print *, potentialDensityTopRef
            print *, 'potentialDensityTopRef(2:nBuoyancyLayers)-potentialDensityTopRef(:nBuoyancyLayers-1)'
            print *, potentialDensityTopRef(2:nBuoyancyLayers)-potentialDensityTopRef(:nBuoyancyLayers-1)
            print *, 'potentialDensityMidRef (nBuoyancyLayers)'
            print *, potentialDensityMidRef
            print *, 'potentialDensityMidRef(2:nBuoyancyLayers)-potentialDensityMidRef(:nBuoyancyLayers-1)'
            print *, potentialDensityMidRef(2:nBuoyancyLayers)-potentialDensityMidRef(:nBuoyancyLayers-1)
            print *, 'nCells,nBuoyancyLayers', nCells,nBuoyancyLayers
            print *, 'nCells*nBuoyancyLayers', nCells*nBuoyancyLayers
            print *, 'No. valid cells in buoyancy coords sum(buoyancyMask)', sum(buoyancyMask)
            print *, 'nCells*nVertLevels', nCells*nVertLevels
            print *, 'sum(cellMask)', sum(cellMask)
            print *, 'minval(potentialDensity), maxval(potentialDensity)'
            print *, minval(potentialDensity), maxval(potentialDensity)
            print *, 'minval(density), maxval(density)'
            print *, minval(density), maxval(density)
         endif

         !-------------------------------------------------------------
         ! INTERPOLATION TEST 1
         ! stratified, horizontally uniform
         ! Interpolating from z, rho to z, rho
         !-------------------------------------------------------------
         if(config_AM_eliassenPalm_debug) then
           do i = 1, nCells
              array1_3D(:,i) = -zMid(:,nCells/2)
              array2_3D(:,i) = potentialDensity(:,nCells/2)
           end do
           print *, ' '
           print *, 'TEST1: Testing interpolatoin function'
           print *, 'Interpolating from (z, rho) to (z, rho)'
           print *, 'call linear_interp_1d_field_along_column(nVertLevels, nCells, ' &
              // 'nVertLevels, maxLevelCell, array1_3D, array2_3D, array1_3D(:,1), array3_3D)'

           print *, 'sum(array1_3D)/nCells + sum(zMid(:,nCells/2))'
           print *, sum(array1_3D)/nCells + sum(zMid(:,nCells/2))
           print *, 'sum(array1_3D)/nCells - sum(array1_3D(:,1))'
           print *, sum(array1_3D)/nCells - sum(array1_3D(:,1))

           call linear_interp_1d_field_along_column(nVertLevels, nCells, nVertLevels, &
              maxLevelCell, array1_3D, array2_3D, array1_3D(:,1), array3_3D)
           print *, 'array1_3D(:,1)'
           print *, array1_3D(:,1)
           print *, 'zMid(:,nCells/2)'
           print *, zMid(:,nCells/2)
           print *, 'array2_3D(:,1)'
           print *, array2_3D(:,1)
           print *, 'array3_3D(:,1)'
           print *, array3_3D(:,1)
           print *, 'array2_3D(:,1)-array3_3D(:,1)'
           print *, array2_3D(:,1)-array3_3D(:,1)

           do i = 1,nCells
              do k = 1, maxLevelCell(i)
                 RMSlocal1 = RMSlocal1 + &
                    ((array3_3D(k,i) - array2_3D(k,i)))**2
                    !((array3_3D(k,i) - array2_3D(k,i))/array2_3D(k,i))**2
              end do
           end do
         endif

         !-------------------------------------------------------------
         ! INTERPOLATION TEST 2
         ! Define a stratification where potential density varies linearly with depth.
         ! Using reference potential density that varies linearly with index.
         ! Interpolate z from that potential density to reference potential density.
         ! Compare to expected values.
         !-------------------------------------------------------------
         if(config_AM_eliassenPalm_debug) then
            do i = 1,nCells
               do k = 1, nVertLevels
                  array1_3D(k,i) = config_AM_eliassenPalm_rhomin_buoycoor*1.02_RKIND + &
                     (zMid(k,i)-zMid(1,i)) * &
                     (config_AM_eliassenPalm_rhomax_buoycoor*0.98_RKIND &
                     - config_AM_eliassenPalm_rhomin_buoycoor*1.02_RKIND) / &
                     (zMid(nVertLevels,i) - zMid(1,i))
                  array2_3D(k,i) = config_AM_eliassenPalm_rhomin_buoycoor*1.02_RKIND + &
                     (zTop(k,i)-zMid(1,i)) * &
                     (config_AM_eliassenPalm_rhomax_buoycoor*0.98_RKIND &
                     - config_AM_eliassenPalm_rhomin_buoycoor*1.02_RKIND) / &
                     (zMid(nVertLevels,i) - zMid(1,i))
               end do
            end do
            call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
               maxLevelCell, array1_3D, zMid, potentialDensityMidRef, array1_3Dbuoy)

            do i = 1,nCells
               do k = 1, nBuoyancyLayers
                  array2_3Dbuoy(k,i) = zMid(1,i) + &
                     (potentialDensityMidRef(k) - config_AM_eliassenPalm_rhomin_buoycoor*1.02_RKIND) * &
                     (zMid(nVertLevels,i) - zMid(1,i)) / &
                     (config_AM_eliassenPalm_rhomax_buoycoor*0.98_RKIND &
                     - config_AM_eliassenPalm_rhomin_buoycoor*1.02_RKIND)
               end do
            end do
            do i = 1,nCells
               do k = 1, nBuoyancyLayers
                  RMSlocal2 = RMSlocal2 + &
                     ((array1_3Dbuoy(k,i) - array2_3Dbuoy(k,i))/array2_3Dbuoy(k,i))**2
               end do
            end do
            print *, ' '
            print *, 'TEST2: Testing interpolation function'
            print *, 'interpolating a linear function'
            print *, 'array1_3Dbuoy(:,nCells/2)'
            print *, array1_3Dbuoy(:,nCells/2)
            print *, 'array2_3Dbuoy(:,nCells/2)'
            print *, array2_3Dbuoy(:,nCells/2)
            print *, 'array1_3Dbuoy(:,nCells/2) - array2_3Dbuoy(:,nCells/2)'
            print *, array1_3Dbuoy(:,nCells/2) - array2_3Dbuoy(:,nCells/2)
         endif


         !-------------------------------------------------------------
         ! check to see if at any point in the domain:
         !    potentialDensity < potentialDensityTopRef(1)
         !    potentialDensity > potentialDensityTopRef(nBuoyancyLayersP1)
         ! either case means that buoyancy coordinate does not span the fluid domain
         !-------------------------------------------------------------
         call check_potentialDensityRef_range(nVertLevels, nCells, maxLevelCell, potentialDensity)

         !-------------------------------------------------------------
         ! interpolate state variable from z-space into buoyancy-space
         !-------------------------------------------------------------
         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, potentialDensity, zMid, &
            potentialDensityMidRef, heightMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, potentialDensity, zMid, &
            potentialDensityTopRef, heightTopBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, potentialDensity, velocityZonal, &
            potentialDensityMidRef, uMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, potentialDensity, velocityMeridional, &
            potentialDensityMidRef, vMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, potentialDensity, density, &
            potentialDensityMidRef, densityMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, potentialDensity, density, &
            potentialDensityTopRef, densityTopBuoyCoor)

         ! Diabatic terms
         !call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
         !   maxLevelCell, potentialDensity, wCellCenter, potentialDensityTopRef, wMidBuoyCoor)

         !-------------------------------------------------------------
         ! fill in data above firstLayerBuoyCoor and below lastLayerBuoyCoor
         !-------------------------------------------------------------
         do i = 1, nCells
            do k = 1, firstLayerBuoyCoor(i)-1
              heightMidBuoyCoor(k,i) = zTop(1,i)
              heightTopBuoyCoor(k,i) = zTop(1,i)
              uMidBuoyCoor(k,i) = velocityZonal(1,i)
              vMidBuoyCoor(k,i) = velocityMeridional(1,i)
              densityMidBuoyCoor(k,i) = density(1,i)
              densityTopBuoyCoor(k,i) = density(1,i)
              ! diabatic
              !wMidBuoyCoor(k,i) = wCellCenter(1,i)
            end do
            do k = lastLayerBuoyCoor(i) + 1, nBuoyancyLayers
              heightMidBuoyCoor(k,i) = -bottomDepth(i)
              heightTopBuoyCoor(k,i) = -bottomDepth(i)
              uMidBuoyCoor(k,i) = velocityZonal(maxLevelCell(i),i)
              vMidBuoyCoor(k,i) = velocityMeridional(maxLevelCell(i),i)
              densityMidBuoyCoor(k,i) = density(maxLevelCell(i),i)
              densityTopBuoyCoor(k,i) = density(maxLevelCell(i),i)
              ! diabatic
              !wMidBuoyCoor(k,i) = wCellCenter(maxLevelCell(i),i)
            end do
            heightInterfaceBuoyCoor(1:nBuoyancyLayers,i) = heightTopBuoyCoor(1:nBuoyancyLayers,i)
            heightInterfaceBuoyCoor(nBuoyancyLayers+1,i) = -bottomDepth(i)
         end do

         !-------------------------------------------------------------
         ! compute sigma, aka "layer thickness", units of s^2
         !-------------------------------------------------------------
         call computeSigma(nCells, nBuoyancyLayers, &
            heightInterfaceBuoyCoor, buoyancyInterfaceRef, sigma)

         !-------------------------------------------------------------
         ! using data interpolated to buoyancy space, compute Montgomery potential
         !-------------------------------------------------------------
         call computeMontgomeryPotential(nBuoyancyLayers, nCells, atmosphericPressure, &
            densityMidBuoyCoor, potentialDensityMidRef, heightInterfaceBuoyCoor, montgPotBuoyCoor)

         !-------------------------------------------------------------
         ! compute the normal derivative of Montgomery potential at cell edges
         !-------------------------------------------------------------
         call computeNormalGradientOnEdge(nBuoyancyLayers, nCells, nEdges, &
            meshPool, montgPotBuoyCoor, montgPotNormalGradOnEdge)

         !-------------------------------------------------------------
         ! reconstruct full gradient vector at cell centers
         !-------------------------------------------------------------
         call mpas_reconstruct(meshPool, montgPotNormalGradOnEdge, &
            montgPotGradX, montgPotGradY, montgPotGradZ, &
            montgPotGradZonal, montgPotGradMerid, includeHalos=.true.)

         !-------------------------------------------------------------
         ! Increment first-order running ensemble average fields
         !-------------------------------------------------------------
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, buoyancyMask, buoyancyMaskEA)
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, sigma, sigmaEA)
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, heightMidBuoyCoor, heightMidBuoyCoorEA)
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, montgPotBuoyCoor, montgPotBuoyCoorEA)
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, montgPotGradZonal, montgPotGradZonalEA)
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, montgPotGradMerid, montgPotGradMeridEA)

         !-------------------------------------------------------------
         ! Increment second-order running ensemble average fields
         !-------------------------------------------------------------
         wrk3DBuoyCoor = heightMidBuoyCoor * heightMidBuoyCoor
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, heightMidBuoyCoorSqEA)

         wrk3DBuoyCoor = heightMidBuoyCoor * montgPotGradZonal
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, heightMGradZonalEA)

         wrk3DBuoyCoor = heightMidBuoyCoor * montgPotGradMerid
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, heightMGradMeridEA)

         wrk3DBuoyCoor = uMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, usigmaEA)

         wrk3DBuoyCoor = vMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, vsigmaEA)

         ! Diabatic terms
         !wrk3DBuoyCoor = wMidBuoyCoor * sigma
         !call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, varpisigmaEA)
         varpisigmaEA = 0.0_RKIND

         !-------------------------------------------------------------
         ! Increment third-order running ensemble average fields
         !-------------------------------------------------------------
         wrk3DBuoyCoor = uMidBuoyCoor * uMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, uusigmaEA)

         wrk3DBuoyCoor = vMidBuoyCoor * vMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, vvsigmaEA)

         wrk3DBuoyCoor = uMidBuoyCoor * vMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, uvsigmaEA)

         ! Diabatic terms
         !wrk3DBuoyCoor = uMidBuoyCoor * wMidBuoyCoor * sigma
         !call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, uvarpisigmaEA)
         uvarpisigmaEA = 0.0_RKIND

         ! Diabatic terms
         !wrk3DBuoyCoor = vMidBuoyCoor * wMidBuoyCoor* sigma
         !call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, wrk3DBuoyCoor, vvarpisigmaEA)
         vvarpisigmaEA = 0.0_RKIND

         !-------------------------------------------------------------
         ! update number of samples in ensemble average
         !-------------------------------------------------------------
         nSamplesEA = nSamplesEA + 1

         !-------------------------------------------------------------
         ! compute the thickness-weighted average velocity
         ! based on current estimate of ensemble-average state
         !-------------------------------------------------------------
         call calculateTWA(nBuoyancyLayers, nCells, nBuoyancyLayers, sigmaEA, usigmaEA,  uTWA)
         call calculateTWA(nBuoyancyLayers, nCells, nBuoyancyLayers, sigmaEA, vsigmaEA,  vTWA)
         ! Diabatic terms
         !call calculateTWA(nBuoyancyLayers, nCells, nBuoyancyLayers, sigmaEA, varpisigmaEA,  varpiTWA)
         varpiTWA = 0.0_RKIND

         !-------------------------------------------------------------
         ! compute the Eliassen-Palm flux tensor
         ! based on current estimate of ensemble-average state
         !-------------------------------------------------------------
         call calculateEPFTfromTWA(nBuoyancyLayers, nCells, &
            sigmaEA, heightMidBuoyCoorEA, &
            heightMidBuoyCoorSqEA, montgPotGradZonalEA, montgPotGradMeridEA, &
            heightMGradZonalEA, heightMGradMeridEA, uTWA, vTWA, varpiTWA, &
            uusigmaEA, vvsigmaEA, uvsigmaEA, uvarpisigmaEA, vvarpisigmaEA, EPFT)

         !-------------------------------------------------------------
         ! Calculate eddy correlations
         !-------------------------------------------------------------
         call calculateCorrelationfromTWA(nBuoyancyLayers, nCells, &
            sigmaEA, uTWA, uTWA, uuSigmaEA, uuTWACorr)
         call calculateCorrelationfromTWA(nBuoyancyLayers, nCells, &
            sigmaEA, vTWA, vTWA, vvSigmaEA, vvTWACorr)
         call calculateCorrelationfromTWA(nBuoyancyLayers, nCells, &
            sigmaEA, uTWA, vTWA, uvSigmaEA, uvTWACorr)
         call calculateEPEfromTWA(nBuoyancyLayers, nCells, &
            sigmaEA, heightMidBuoyCoorEA, heightMidBuoyCoorSqEA, epeTWA)
         call calculateEddyFormDragFromTWA(nBuoyancyLayers, nCells, sigmaEA, &
            heightMidBuoyCoorEA, montgPotGradZonalEA, heightMGradZonalEA, eddyFormDragZonal)
         call calculateEddyFormDragFromTWA(nBuoyancyLayers, nCells, sigmaEA, &
            heightMidBuoyCoorEA, montgPotGradMeridEA, heightMGradMeridEA, eddyFormDragMerid)

         !-------------------------------------------------------------
         ! compute the total force from the EPFT: div(EPFT)
         !-------------------------------------------------------------
         call calculateDivEPFT(config_AM_eliassenPalm_debug, &
            domain % on_a_sphere, rho_sw, nBuoyancyLayers, nCells, nEdges, &
            meshPool, buoyancyMidRef, sigmaEA, buoyancyMaskEA, EPFT, divEPFT)
         ! decompose the vector into its components for output
         divEPFT1 = divEPFT(1,:,:)
         divEPFT2 = divEPFT(2,:,:)

         !-------------------------------------------------------------
         ! compute the force from horizontal shear component of the EPFT
         !-------------------------------------------------------------
         wrkTensor = 0.0_RKIND
         wrkTensor(1:2,1:2,:,:) = EPFT(1:2,1:2,:,:)
         call calculateDivEPFT(config_AM_eliassenPalm_debug, &
            domain % on_a_sphere, rho_sw, nBuoyancyLayers, nCells, nEdges, &
            meshPool, buoyancyMidRef, sigmaEA, buoyancyMaskEA, wrkTensor, wrkVector)
         divEPFTshear1 = wrkVector(1,:,:)
         divEPFTshear2 = wrkVector(2,:,:)

         !-------------------------------------------------------------
         ! compute the force from vertical form drag component of the EPFT
         !-------------------------------------------------------------
         wrkTensor = 0.0_RKIND
         wrkTensor(3,1:2,:,:) = EPFT(3,1:2,:,:)
         call calculateDivEPFT(config_AM_eliassenPalm_debug, &
            domain % on_a_sphere, rho_sw, nBuoyancyLayers, nCells, nEdges, &
            meshPool, buoyancyMidRef, sigmaEA, buoyancyMaskEA, wrkTensor, wrkVector)
         divEPFTdrag1 = wrkVector(1,:,:)
         divEPFTdrag2 = wrkVector(2,:,:)

         !-------------------------------------------------------------
         ! transform div(EPFT) into a flux of Ertel PV
         !-------------------------------------------------------------
         call calculateErtelPVFlux(nCells, nBuoyancyLayers, sigmaEA, divEPFT, ErtelPVFlux)
         ErtelPVFlux1 = ErtelPVFlux(1,:,:)
         ErtelPVFlux2 = ErtelPVFlux(2,:,:)

         !-------------------------------------------------------------
         ! compute Ertel PV tendency from Ertel PV fluxes, div(ErtelPVFlux)
         !-------------------------------------------------------------
         call calculateErtelPVTendencyFromPVFlux(config_AM_eliassenPalm_debug, &
            domain % on_a_sphere, nBuoyancyLayers, nCells, nEdges, &
            meshPool, sigmaEA, ErtelPVFlux, ErtelPVTendency)

         !-------------------------------------------------------------
         ! compute Ertel PV based on EA and TWA fields
         !-------------------------------------------------------------
         call computeErtelPV(nCells, nBuoyancyLayers, nEdges, meshPool, &
            fCell, uTWA, vTWA, sigmaEA, ErtelPV)

         !-------------------------------------------------------------
         ! compute the normal derivative of EPV at cell edges
         !-------------------------------------------------------------
         call computeNormalGradientOnEdge(nBuoyancyLayers, nCells, nEdges, &
            meshPool, ErtelPV, ErtelPVNormalGradOnEdge)

         !-------------------------------------------------------------
         ! reconstruct full gradient vector at cell centers
         !-------------------------------------------------------------
         call mpas_reconstruct(meshPool, ErtelPVNormalGradOnEdge, &
            ErtelPVGradX, ErtelPVGradY, ErtelPVGradZ, ErtelPVGradZonal, ErtelPVGradMerid, includeHalos=.true.)

         !-------------------------------------------------------------
         ! compute the vertical derivative of uTWA
         !-------------------------------------------------------------
         call computeVerticalDerivative(nCells, nBuoyancyLayers, &
           firstLayerBuoyCoor, lastLayerBuoyCoor, heightMidBuoyCoor, uTWA, duTWAdz)

         !-------------------------------------------------------------
         ! compute the vertical derivative of vTWA
         !-------------------------------------------------------------
         call computeVerticalDerivative(nCells, nBuoyancyLayers, &
           firstLayerBuoyCoor, lastLayerBuoyCoor, heightMidBuoyCoor, vTWA, dvTWAdz)

         !-------------------------------------------------------------
         ! Compute the geometric decomposition in terms of angles and
         ! eccentricities using the entries of EPFT.
         ! (not yet implemented)
         !-------------------------------------------------------------
         !call eddyGeomDecompEPFT(EPFT, ...)



         !-------------------------------------------------------------
         ! TEST:
         ! calculate potential vorticity fluxes using curl of u
         !-------------------------------------------------------------
         if(config_AM_eliassenPalm_debug) then

            call mpas_pool_get_array(diagnosticsPool, 'relativeVorticityCell', relativeVorticityCell)

            ! store relVortMidBuoyCoor in array1_3Dbuoy
            call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
               maxLevelCell, potentialDensity, relativeVorticityCell, &
               potentialDensityMidRef, array1_3Dbuoy)

            do i = 1,nCells
               do k=firstLayerBuoyCoor(i), lastLayerBuoyCoor(i)
                  PVMidBuoyCoor(k,i) = (fCell(i) + array1_3Dbuoy(k,i) ) / sigma(k,i)
               end do
            end do

            call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
               uMidBuoyCoor, uMidBuoyCoorEA)

            call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
               vMidBuoyCoor, vMidBuoyCoorEA)

            call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
               PVMidBuoyCoor, PVMidBuoyCoorEA)

            wrk3DBuoyCoor = uMidBuoyCoor * PVMidBuoyCoor
            call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
               wrk3DBuoyCoor, uPVMidBuoyCoorEA)

            wrk3DBuoyCoor = vMidBuoyCoor * PVMidBuoyCoor
            call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
               wrk3DBuoyCoor, vPVMidBuoyCoorEA)

            PVFluxTest(1,:,:) = uPVMidBuoyCoorEA - uMidBuoyCoorEA * PVMidBuoyCoorEA
            PVFluxTest(2,:,:) = vPVMidBuoyCoorEA - vMidBuoyCoorEA * PVMidBuoyCoorEA

            do i = 1,nCells
               do k = firstLayerBuoyCoor(i), lastLayerBuoyCoor(i)
                  RMSPVFlux1Local = RMSPVFlux1local + &
                     ( ErtelPVFlux(1,k,i) - PVFLuxTest(1,k,i) )**2
                  RMSPVFlux2Local = RMSPVFlux2local + &
                     ( ErtelPVFlux(2,k,i) - PVFLuxTest(2,k,i) )**2
               end do
            end do

         end if


         !-------------------------------------------------------------
         ! deallocate scratch space and test space variables
         !-------------------------------------------------------------
         call mpas_deallocate_scratch_field(firstLayerBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(lastLayerBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(heightMidBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(heightTopBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(heightInterfaceBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(sigmaField, .true.)
         call mpas_deallocate_scratch_field(montgPotBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(montgPotNormalGradOnEdgeField, .true.)
         call mpas_deallocate_scratch_field(uMidBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(vMidBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(densityMidBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(densityTopBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(buoyancyMaskField, .true.)
         call mpas_deallocate_scratch_field(montgPotGradXField, .true.)
         call mpas_deallocate_scratch_field(montgPotGradYField, .true.)
         call mpas_deallocate_scratch_field(montgPotGradZField, .true.)
         call mpas_deallocate_scratch_field(montgPotGradZonalField, .true.)
         call mpas_deallocate_scratch_field(montgPotGradMeridField, .true.)
         call mpas_deallocate_scratch_field(wrk3DnVertLevelsP1Field, .true.)
         call mpas_deallocate_scratch_field(wrk3DnVertLevelsField, .true.)
         call mpas_deallocate_scratch_field(wrk3DBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(ErtelPVNormalGradOnEdgeField, .true.)
         call mpas_deallocate_scratch_field(ErtelPVGradXField, .true.)
         call mpas_deallocate_scratch_field(ErtelPVGradYField, .true.)
         call mpas_deallocate_scratch_field(ErtelPVGradZField, .true.)
         call mpas_deallocate_scratch_field(wrkVectorField, .true.)
         call mpas_deallocate_scratch_field(wrkTensorField, .true.)

         call mpas_deallocate_scratch_field(array1_3DField, .true.)
         call mpas_deallocate_scratch_field(array2_3DField, .true.)
         call mpas_deallocate_scratch_field(array3_3DField, .true.)
         call mpas_deallocate_scratch_field(array1_3DbuoyField, .true.)
         call mpas_deallocate_scratch_field(array2_3DbuoyField, .true.)
         call mpas_deallocate_scratch_field(PVMidBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(PVMidBuoyCoorEAField, .true.)
         call mpas_deallocate_scratch_field(uMidBuoyCoorEAField, .true.)
         call mpas_deallocate_scratch_field(vMidBuoyCoorEAField, .true.)
         call mpas_deallocate_scratch_field(uPVMidBuoyCoorEAField, .true.)
         call mpas_deallocate_scratch_field(vPVMidBuoyCoorEAField, .true.)
         call mpas_deallocate_scratch_field(PVFluxTestField, .true.)

         !-------------------------------------------------------------
         ! update test variables
         !-------------------------------------------------------------
         nCellsCum = nCellsCum + nCells

         !-------------------------------------------------------------
         ! move to the next block
         !-------------------------------------------------------------
         block => block % next

      end do

      !-------------------------------------------------------------
      ! TESTS: tallying up tests across processors.
      !-------------------------------------------------------------
      if(config_AM_eliassenPalm_debug) then
        RMSglobal1 = 1.0e+36_RKIND
        call mpas_dmpar_sum_int(dminfo, nCellsCum, nCellsGlobal)
        call mpas_dmpar_sum_real(dminfo, RMSlocal1, RMSglobal1)
        call mpas_dmpar_sum_real(dminfo, RMSlocal2, RMSglobal2)

        if (dminfo % my_proc_id == IO_NODE) then
           print *, ' '
           print *, 'RKIND=', RKIND
           print *, 'rms relative error interp test1:',sqrt(RMSglobal1/nCellsGlobal)
           print *, 'rms relative error interp test2:',sqrt(RMSglobal2/nCellsGlobal)

           print *, ' '
        endif

        call mpas_dmpar_sum_real(dminfo, sum(abs(ErtelPVFlux(1,:,:))), RMSglobal1)
        call mpas_dmpar_max_real(dminfo, maxval(abs(ErtelPVFlux(1,:,:))), RMSglobal2)
        if (dminfo % my_proc_id == IO_NODE) then
           print *, 'Checking ErtelPVFlux'
           print *, 'Global sum(abs(ErtelPVFlux(1,:,:))) = ', RMSglobal1
           print *, 'Global max(abs(ErtelPVFlux(1,:,:))) = ', RMSglobal2
        endif

        call mpas_dmpar_sum_real(dminfo, sum(abs(ErtelPVFlux(2,:,:))), RMSglobal1)
        call mpas_dmpar_max_real(dminfo, maxval(abs(ErtelPVFlux(2,:,:))), RMSglobal2)
        if (dminfo % my_proc_id == IO_NODE) then
           print *, 'Global sum(abs(ErtelPVFlux(2,:,:))) = ', RMSglobal1
           print *, 'Global max(abs(ErtelPVFlux(2,:,:))) = ', RMSglobal2
        endif

        call mpas_dmpar_sum_real(dminfo, sum(abs(ErtelPVFlux(3,:,:))), RMSglobal1)
        call mpas_dmpar_max_real(dminfo, maxval(abs(ErtelPVFlux(3,:,:))), RMSglobal2)
        if (dminfo % my_proc_id == IO_NODE) then
           print *, 'Global sum(abs(ErtelPVFlux(3,:,:))) = ', RMSglobal1
           print *, 'Global max(abs(ErtelPVFlux(3,:,:))) = ', RMSglobal2
        endif

        call mpas_dmpar_sum_real(dminfo, RMSPVFlux1Local, RMSPVFlux1global)
        call mpas_dmpar_sum_real(dminfo, RMSPVFlux2Local, RMSPVFlux2global)
        if (dminfo % my_proc_id == IO_NODE) then
           print *, 'rms relative error test PVFlux1:',sqrt(RMSPVFlux1global/nCellsGlobal)
           print *, 'rms relative error test PVFlux2:',sqrt(RMSPVFLux2global/nCellsGlobal)

           print *, ' '
        endif

      endif

      if(config_AM_eliassenPalm_debug) then
         write(stderrUnit, *) ' '
         write(stderrUnit, *) 'exiting ocn_compute_epft'
         write(stderrUnit, *) ' '
      end if

    end subroutine ocn_compute_eliassen_palm!}}}

!***********************************************************************
!
!  routine ocn_restart_eliassen_palm
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_eliassen_palm(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

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

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_eliassen_palm!}}}

!***********************************************************************
!
!  routine ocn_finalize_eliassen_palm
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Juan A. Saenz
!> \date    May 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_eliassen_palm(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

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

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_finalize_eliassen_palm!}}}



!***********************************************************************
! Local routines start here
!***********************************************************************


!***********************************************************************
!
!  subroutine get_masks_in_buoyancy_coordinates
!
!> \brief   Get masks in buoyancy coordinates
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details
!>  firstLayerBuoyCoor(iCell): the index of the smallest reference density that
!>    is >= the smallest actual density in a column.
!>  lastLayerBuoyCoor(iCell): the index of the largest reference density that is <= the
!>    largest actual density in a column.
!>  Set masks in buoyancy coordinates:
!>    mask = 1: cell is a valid ocean cell
!>    mask = 0: cell is not a valid ocean cell
!>  Required: potentialDensityMidRef monotonically increases with index value
!
!-----------------------------------------------------------------------
   subroutine get_masks_in_buoyancy_coordinates(nVertLevels, nCells, nBuoyancyLayers, &
      maxLevelCell, potentialDensity, potentialDensityMidRef, &
      firstLayerBuoyCoor, lastLayerBuoyCoor, buoyancyMask)!{{{

      !-----------------------------------------------------------------
      ! intent(in)
      !-----------------------------------------------------------------
      integer, intent(in) :: nVertLevels, nCells, nBuoyancyLayers
      integer, dimension(nCells), intent(in)  :: maxLevelCell
      real (kind=RKIND), dimension(nVertLevels, nCells), intent(in) :: potentialDensity
      real (kind=RKIND), dimension(nBuoyancyLayers), intent(in) :: potentialDensityMidRef

      !-----------------------------------------------------------------
      ! intent(out)
      !-----------------------------------------------------------------
      integer, dimension(nCells), intent(out) :: firstLayerBuoyCoor
      integer, dimension(nCells), intent(out) :: lastLayerBuoyCoor
      real (kind=RKIND), dimension(nBuoyancyLayers, nCells), intent(out) :: buoyancyMask

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      integer :: iCell, maxLevel, kB, kBBottom, kBTop

      !-----------------------------------------------------------------
      ! initialize fields assuming no density layers exist
      !-----------------------------------------------------------------
      firstLayerBuoyCoor = nBuoyancyLayers
      lastLayerBuoyCoor = 1
      buoyancyMask     = 0.0_RKIND

      !-----------------------------------------------------------------
      ! loop over all cells
      !    when searching from the top down
      !      find first target density greater than density in top model layer
      !    when searching from the bottom up
      !      find first target density less than density in bottom model layer
      !-----------------------------------------------------------------
      do iCell = 1, nCells

         ! find the bottom model layer for this cell
         maxLevel = maxLevelCell(iCell)

         ! search top down
         do kB = 1, nBuoyancyLayers
            if (potentialDensityMidRef(kB) >= potentialDensity(1,iCell) ) then
               firstLayerBuoyCoor(iCell) = kB
               exit
            endif
         enddo

         ! search bottom up
         do kB = nBuoyancyLayers, 1, -1
            if (potentialDensityMidRef(kB) <= potentialDensity(maxLevel,iCell) ) then
               lastLayerBuoyCoor(iCell) = kB
               exit
            endif
         enddo

         ! set mask to 1 inside the range
         do kB = firstLayerBuoyCoor(iCell), lastLayerBuoyCoor(iCell)
           buoyancyMask(kB,iCell) = 1.0_RKIND
         enddo

      enddo

   end subroutine get_masks_in_buoyancy_coordinates!}}}


!***********************************************************************
!
!  subroutine check_potentialDensityRef_range
!
!> \brief   Check if the range of values in potentialDensityTopRef contains current state
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details
!>  Check if the range of values in potentialDensityTopRef contains all values in
!>  potentialDensity of the current state.
!>  If not, print a warning.
!
!-----------------------------------------------------------------------
   subroutine check_potentialDensityRef_range(nVertLevels, nCells, maxLevelCell, &
         potentialDensity)!{{{
      !-----------------------------------------------------------------
      ! intent(in)
      !-----------------------------------------------------------------
      integer, intent(in) :: nVertLevels, nCells
      integer, dimension(nCells), intent(in)  :: maxLevelCell
      real (kind=RKIND), dimension(nVertLevels, nCells), intent(in) :: potentialDensity

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      integer :: k, iCell, iCellMinBound, iCellMaxBound
      logical :: printWarning

      real (kind=RKIND), pointer :: config_AM_eliassenPalm_rhomin_buoycoor, &
        config_AM_eliassenPalm_rhomax_buoycoor

      call mpas_pool_get_config(ocnConfigs, 'config_AM_eliassenPalm_rhomin_buoycoor', &
        config_AM_eliassenPalm_rhomin_buoycoor)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_eliassenPalm_rhomax_buoycoor', &
        config_AM_eliassenPalm_rhomax_buoycoor)

      printWarning = .false.
      iCellMinBound = -1
      iCellMaxBound = -1

      do iCell = 1, nCells
         if (potentialDensity(1,iCell) < config_AM_eliassenPalm_rhomin_buoycoor) then
            printWarning = .true.
            iCellMinBound = iCell
            exit
         end if
         if (potentialDensity(maxLevelCell(iCell),iCell) > config_AM_eliassenPalm_rhomax_buoycoor) then
            printWarning = .true.
            iCellMaxBound = iCell
            exit
         end if
      enddo

      if (printWarning) then
         write(stderrUnit,*)
         write(stderrUnit,*) ' *** WARNING: in eliassen_palm analysis member, subroutine check_potentialDensityRef_range.'
         write(stderrUnit,*) '  One or more columns in the ocean domain have densities that are not'
         write(stderrUnit,*) '  contained in the defined buoyancy space of the EPFT module'
         if (iCellMinBound.gt.0) write(stderrUnit,*) '  fluid is lighter than min buoyancy at cell: ',iCellMinBound
         if (iCellMaxBound.gt.0) write(stderrUnit,*) '  fluid is lighter than max buoyancy at cell: ',iCellMaxBound
         write(stderrUnit,*)
      end if

   end subroutine check_potentialDensityRef_range!}}}


!***********************************************************************
!
!  subroutine linear_interp_1d_field_along_column
!
!> \brief   One-dimensional interpolation in buoyancy coordinates
!> \author  Juan A. Saenz, Todd Ringler
!> \date    17 December 2013
!> \details
!>  Interpolate a field yFieldIn residing on xFieldIn onto xColumnOut and store
!>  and return in yFieldOut.
!>  Interpolation is done using one-dimensional interpolation along xColumnOut.
!>  Required: xFieldIn monotonically increases with index value
!
!-----------------------------------------------------------------------

   subroutine linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
      maxLevelCell, xFieldIn, yFieldIn, xColumnOut, yFieldOut)!{{{

      !-----------------------------------------------------------------
      ! intent(in)
      !-----------------------------------------------------------------
      integer, intent(in) :: nVertLevels, nCells, nBuoyancyLayers
      integer, dimension(nCells), intent(in)  :: maxLevelCell
      real (kind=RKIND), dimension(nVertLevels, nCells), intent(in) :: xFieldIn
      real (kind=RKIND), dimension(nVertLevels, nCells), intent(in) :: yFieldIn
      real (kind=RKIND), dimension(nBuoyancyLayers), intent(in) :: xColumnOut

      !-----------------------------------------------------------------
      ! intent(out)
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(nBuoyancyLayers, nCells), intent(out) :: yFieldOut

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      integer :: iCell, maxLevel, kB, kBBottom, kBTop, kDataAbove, kDataBelow, kData
      real (kind=RKIND) :: dx, dy
      real (kind=RKIND), dimension(nVertLevels, nCells) :: xSrc ! source data
      real (kind=RKIND), dimension(nBuoyancyLayers) :: xDst     ! destination data

      ! jas issue: the code below works for monotonically decreasing arrays.
      ! however above it expects arguments that are monotonically increasing arrays.
      xSrc = -xFieldIn
      xDst = -xColumnOut

      !-----------------------------------------------------------------
      ! test for monoticity of xSrc
      !-----------------------------------------------------------------
      ! jas issue : to do

      !-----------------------------------------------------------------
      ! initialize intent(out)
      !-----------------------------------------------------------------
      yFieldOut = 0.0_RKIND

      !-----------------------------------------------------------------
      ! loop over all columns
      !-----------------------------------------------------------------
      do iCell = 1, nCells

         ! find the index of the bottom level of a column
         maxLevel = maxLevelCell(iCell)

         ! Monotonically decreasing xSrc required
         ! Find index of first element in xDst that is inside xSrc(:,iCell)
         kBTop = 1
         do kB = 1, nBuoyancyLayers
            ! the following line ensures that
            ! if all xDst > xSrc(1,iCell) then kBTop = nBuoyancyLayers
            kBTop = kB
            if (xDst(kB) <= xSrc(1,iCell) ) then
               exit
            endif
         enddo

         !find last target buoyancy level inside column
         kBBottom = nBuoyancyLayers
         do kB = nBuoyancyLayers, 1, -1
            ! the following line ensures that
            ! if all xDst < xSrc(1,iCell) then kBBottom = 1
            kBBottom = kB
            if (xDst(kB) >= xSrc(maxLevel,iCell) ) then
               exit
            endif
         enddo

         ! For the target x levels outside the x range in a column:
         ! set data from 1:kBTop-1 to surface values
         do kB = 1, kBTop-1
           yFieldOut(kB,iCell) = yFieldIn(1,iCell)
         enddo
         !set data from kBBottom+1:nBuoyancyLayers to bottom values
         do kB = kBBottom+1, nBuoyancyLayers
           yFieldOut(kB,iCell) = yFieldIn(maxLevel,iCell)
         enddo

         ! The interpolation:
         ! for the target buoyancy levels within the buoyancy range in a column:
         kDataAbove = 1
         kDataBelow = kDataAbove + 1
         do kB = kBTop, kBBottom
            ! for each xDst(kB) value, find the corresponding upper and lower
            ! xSrc value in the field data, then interpolate y between those values.
            if (xDst(kB) < xSrc(kDataBelow,iCell)) then
               do kData = kDataBelow, maxLevel
                  if (xDst(kB) > xSrc(kData,iCell) ) then
                     kDataBelow=kData
                     kDataAbove=kDataBelow-1
                     exit
                  endif
               enddo
            endif

            dx = xSrc(kDataBelow,iCell) - xSrc(kDataAbove,iCell)
            dy = yFieldIn(kDataBelow,iCell) - yFieldIn(kDataAbove,iCell)
            yFieldOut(kB,iCell) = yFieldIn(kDataAbove,iCell) + &
               (xDst(kB)-xSrc(kDataAbove,iCell)) * dy/dx
         enddo

      enddo

   end subroutine linear_interp_1d_field_along_column!}}}


!***********************************************************************
!
!  subroutine computeSigma
!
!> \brief   Calculate the inverse of the derivative of buoy wrt z
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details
!>  This subroutine calculates the inverse of the derivative of buoy wrt z.
!
!-----------------------------------------------------------------------

   subroutine computeSigma(nCells, nLayers, &
         heightInterface, buoyInterface, sigma)!{{{
      !-----------------------------------------------------------------
      ! intent(in)
      !-----------------------------------------------------------------
      integer, intent(in) :: nCells, nLayers
      real (kind=RKIND), dimension(:,:), intent(in) :: heightInterface
      real (kind=RKIND), dimension(:), intent(in) :: buoyInterface

      !-----------------------------------------------------------------
      ! intent(out)
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:,:), intent(out) :: sigma

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: iCell, k

      !-----------------------------------------------------------------
      ! initialize sigma assuming zero thickness layers everywhere
      !-----------------------------------------------------------------
      sigma = 0.0_RKIND

      !-----------------------------------------------------------------
      ! loop over all column, sigam = delta z / delta b
      !  note: positive z points "up", i.e. from k+1 to k
      !  note: positive b points "up", i.e. from k+1 to k
      !-----------------------------------------------------------------
      do iCell = 1, nCells
         do k = 1,nLayers
            sigma(k,iCell) = (heightInterface(k+1,iCell) - heightInterface(k,iCell)) / &
               (buoyInterface(k+1) - buoyInterface(k))
         enddo
      enddo

   end subroutine computeSigma!}}}



!***********************************************************************
!
!  subroutine computeMontgomeryPotential
!
!> \brief   Compute the Montgomery potential
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 20015
!> \details
!>  This subroutine computes the Montgomery potential using eqn 2.10 in
!>  R.L. Higdon and R.A. Szoeke (1997), J. Comp. Phys. 135, 30–53, Article No. CP975733
!
!>  Montgomery Potential (MP) in layer k is MP(k-1) + pInterface(k)*deltaAlpha
!>  where deltaAlpha is (1/potDens(k) - 1/potDens(k-1))
!>  and pInterface(k) is the pressure at interface k, i.e. at top of layer k.
!>
!>  Montgomery potential of a layer is constant across layer
!-----------------------------------------------------------------------

   subroutine computeMontgomeryPotential(nLayers, nCells, pSurface, &
         density, potDens, heightInterface, MontgomeryPotential)!{{{

      !-----------------------------------------------------------------
      ! intent(in)
      !-----------------------------------------------------------------
      integer, intent(in) :: nLayers, nCells
      real (kind=RKIND), dimension(nCells), intent(in) :: pSurface
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: density
      real (kind=RKIND), dimension(nLayers), intent(in) :: potDens
      real (kind=RKIND), dimension(nLayers+1, nCells), intent(in) :: heightInterface

      !-----------------------------------------------------------------
      ! intent(out)
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(nLayers, nCells), intent(out) :: MontgomeryPotential

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: iCell, k
      real (kind=RKIND) :: pInterfacek ! pressure at interface k, i.e. at top of layer k

      !-----------------------------------------------------------------
      ! initialize intent(out)
      !-----------------------------------------------------------------
      MontgomeryPotential = 0.0_RKIND

      !-----------------------------------------------------------------
      ! loop over all columns
      !-----------------------------------------------------------------
      do iCell = 1, nCells

         !-----------------------------------------------------------------
         ! compute Montgomery potential in top buoyancy layer
         ! at present, assume atmosphere surface pressure is zero (or a constant)
         !-----------------------------------------------------------------
         pInterfacek = 0.0_RKIND
         k = 1
         MontgomeryPotential(k,iCell) = pInterfacek/potDens(k) + gravity*heightInterface(k,iCell)

         !-----------------------------------------------------------------
         ! compute Montgomery potential by accumulating jump across each layer interace
         !   Jump == pressure at interface * (alpha (below interface) - alpha (above interface))
         !-----------------------------------------------------------------
         do k = 2, nLayers
            pInterfacek = pInterfacek + &
               gravity * ( heightInterface(k-1,iCell)-heightInterface(k,iCell) ) * density(k-1,iCell)
            MontgomeryPotential(k,iCell) = MontgomeryPotential(k-1,iCell) + &
               pInterfacek * ( 1/potDens(k) - 1/potDens(k-1) )
         enddo

      enddo

   end subroutine computeMontgomeryPotential!}}}



!***********************************************************************
!
!  subroutine computeNormalGradientOnEdge
!
!> \brief   Compute the gradient of a quantity that exists on cell centers
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details
!>  This subroutine computes the normal derivative of a scalar
!>  quantity that exists on cell centers. Routine assumes that
!>  data is valid throughout the entire column, as is the case
!>  when working in buoyancy coordinates
!
!-----------------------------------------------------------------------

   subroutine computeNormalGradientOnEdge(nBLayers, nCells, nEdges, &
         meshPool, field, normalGradOnEdge)!{{{

      !-----------------------------------------------------------------
      ! intent(in)
      !-----------------------------------------------------------------
      integer, intent(in) :: nBLayers, nCells, nEdges
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      real (kind=RKIND), dimension(:,:), intent(in) :: field

      !-----------------------------------------------------------------
      ! intent(out)
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:,:), intent(out) :: normalGradOnEdge

      !-----------------------------------------------------------------
      !local variables
      !-----------------------------------------------------------------
      integer :: iEdge, k, cell1, cell2, kMin, kMax
      integer, pointer :: nBuoyancyLayers
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: boundaryEdge
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND) :: invLength

      !-----------------------------------------------------------------
      ! assign pointers
      !-----------------------------------------------------------------
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'boundaryEdge', boundaryEdge)

      call mpas_pool_get_dimension(meshPool, 'nBuoyancyLayers', nBuoyancyLayers)

      !-----------------------------------------------------------------
      ! initialize intent(out)
      !-----------------------------------------------------------------
      normalGradOnEdge = 0.0_RKIND

      !-----------------------------------------------------------------
      ! loop over edges, compute derivative as (cell2 - cell1) / dc
      !-----------------------------------------------------------------
      do iEdge = 1, nEdges
         ! do not compute the normal derivative at land/sea interface
         if (boundaryEdge(1,iEdge) == 1) then
            normalGradOnEdge(:,iEdge) = 0.0_RKIND
         else
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            invLength = 1.0_RKIND / dcEdge(iEdge)
            do k = 1, nBuoyancyLayers
              normalGradOnEdge(k,iEdge) = ( field(k,cell2) - field(k,cell1) )*invLength
            enddo
         end if
      enddo

   end subroutine computeNormalGradientOnEdge!}}}



!***********************************************************************
!
!  subroutine updateEnsembleAverage
!
!> \brief   Update ensemble average
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details
!>  This subroutine updates the ensemble average
!
!-----------------------------------------------------------------------

      subroutine updateEnsembleAverage(nLayers, nCells, nSamples, A, Abar)!{{{

      !-----------------------------------------------------------------
      ! intent(in)
      !-----------------------------------------------------------------
      integer, intent(in) :: nLayers, nCells, nSamples
      real (kind=RKIND), dimension(nLayers, nCells), intent(in)  :: A

      !-----------------------------------------------------------------
      ! intent(inout)
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(nLayers, nCells), intent(inout)  :: Abar

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: iCell, k

      !-----------------------------------------------------------------
      ! Abar is the current estimate of the ensemble average
      ! on input, Abar was built using nSamples of A
      ! To update Abar, we multiple Abar times nSamples, add in the
      ! current value (A), then normalize by (nSamples + 1)
      !-----------------------------------------------------------------
      do iCell = 1, nCells
         do k = 1, nLayers
            Abar(k,iCell) = (nSamples * Abar(k,iCell) + A(k,iCell)) / (nSamples + 1.0_RKIND)
         enddo
      enddo

   end subroutine updateEnsembleAverage!}}}


!***********************************************************************
!
!  subroutine calculateTWA
!
!> \brief   Calculate the thickness weighted average
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details
!>  This subroutine calculates the thickness weighted average
!
!-----------------------------------------------------------------------
      subroutine calculateTWA(nLayers, nCells, nBuoyancyLayers, sigmaEA, &
            varSigmaEA, varTWA)!{{{

      !-----------------------------------------------------------------
      ! intent(in)
      !-----------------------------------------------------------------
      integer, intent(in) :: nLayers, nCells, nBuoyancyLayers
      real (kind=RKIND), dimension(nLayers, nCells), intent(in)  :: sigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in)  :: varSigmaEA

      !-----------------------------------------------------------------
      ! intent(inout)
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(nLayers, nCells), intent(out) :: varTWA

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: iCell, k

      !-----------------------------------------------------------------
      ! initialize intent(out)
      !-----------------------------------------------------------------
      varTWA = 0.0_RKIND

      do iCell = 1, nCells
         do k = 1,nBuoyancyLayers
            varTWA(k,iCell) = varSigmaEA(k,iCell) / max(epsilonEPFT,sigmaEA(k,iCell))
         enddo
      enddo

   end subroutine calculateTWA!}}}


!***********************************************************************
!
!  subroutine calculateEPFTfromTWA
!
!> \brief   Calculate the Eliassen-Palm flux tensor from TWAs
!> \author  Juan A. Saenz
!> \date    January 2014
!> \details
!>  This subroutine calculates the Eliassen and Palm flux tensor from thickness
!>  weighted averages.
!>  EPTF_pq(x,y,z) is represented as EPFT(p,q,k,i)
!-----------------------------------------------------------------------

   subroutine calculateEPFTfromTWA(nLayers, nCells, &
         sigmaEA, heightEA, heightSqEA, MxEA, MyEA, HMxEA, HMyEA, uTWA, vTWA, varpiTWA, &
         uuSigmaEA, vvSigmaEA, uvSigmaEA, uvarpisigmaEA, vvarpisigmaEA, Etensor)!{{{
      implicit none
      integer, intent(in) :: nLayers, nCells
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: sigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: heightEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: heightSqEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: MxEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: MyEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: HMxEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: HMyEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: uTWA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: vTWA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: varpiTWA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: uuSigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: vvSigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: uvSigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: uvarpisigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: vvarpisigmaEA
      real (kind=RKIND), dimension(3, 3, nLayers, nCells), intent(out) :: Etensor

      ! local variables
      integer :: iCell, kLayer
      real (kind=RKIND) :: sigma
      real (kind=RKIND) :: uppupp, vppvpp, uppvpp, uppwpp, vppwpp
      real (kind=RKIND) :: HpHp, HpMxp, HpMyp
      real (kind=RKIND) :: dummy1, dummy2, dummy3

      Etensor = 0.0_RKIND

      do iCell = 1, nCells
         do kLayer = 1,nLayers

            sigma = max(sigmaEA(kLayer,iCell), epsilonEPFT)

            uppupp = uuSigmaEA(kLayer,iCell) / sigma - uTWA(kLayer,iCell)*uTWA(kLayer,iCell)
            vppvpp = vvSigmaEA(kLayer,iCell) / sigma - vTWA(kLayer,iCell)*vTWA(kLayer,iCell)
            uppvpp = uvSigmaEA(kLayer,iCell) / sigma - uTWA(kLayer,iCell)*vTWA(kLayer,iCell)
            uppwpp = uvarpisigmaEA(kLayer,iCell) / sigma - uTWA(kLayer,iCell)*varpiTWA(kLayer,iCell)
            vppwpp = vvarpisigmaEA(kLayer,iCell) / sigma - vTWA(kLayer,iCell)*varpiTWA(kLayer,iCell)

            HpHp   = heightSqEA(kLayer,iCell) - heightEA(kLayer,iCell)*heightEA(kLayer,iCell)
            HpMxp  = HMxEA(kLayer,iCell) - heightEA(kLayer,iCell) * MxEA(kLayer,iCell)
            HpMyp  = HMyEA(kLayer,iCell) - heightEA(kLayer,iCell) * MyEA(kLayer,iCell)

            !EPTF_pq(x,y,z) is represented as EPFT(p,q,kLayer,iCell)
            !column 1: Eu
            Etensor(1,1,kLayer,iCell) = uppupp + 0.5_RKIND * HpHp / sigma
            Etensor(2,1,kLayer,iCell) = uppvpp
            Etensor(3,1,kLayer,iCell) = uppwpp + HpMxp / sigma

            !column 2: Ev
            Etensor(1,2,kLayer,iCell) = uppvpp
            Etensor(2,2,kLayer,iCell) = vppvpp + 0.5_RKIND * HpHp / sigma
            Etensor(3,2,kLayer,iCell) = vppwpp + HpMyp / sigma

            !column 3: Ew
            Etensor(1,3,kLayer,iCell) = 0.0_RKIND
            Etensor(2,3,kLayer,iCell) = 0.0_RKIND
            Etensor(3,3,kLayer,iCell) = 0.0_RKIND

         enddo
      enddo

   end subroutine calculateEPFTfromTWA!}}}


!***********************************************************************
!
!  subroutine calculateCorrelationfromTWA
!
!> \brief   Calculate the eddy correlation from TWAs
!> \author  Juan A. Saenz
!> \date    July 2015
!> \details
!>  This subroutine calculates the eddy kinetic energy from thickness
!>  weighted averages.
!-----------------------------------------------------------------------

   subroutine calculateCorrelationfromTWA(nLayers, nCells, &
         sigmaEA, uTWA, vTWA, uvSigmaEA, uvCorr)!{{{
      implicit none
      integer, intent(in) :: nLayers, nCells
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: sigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: uTWA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: vTWA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: uvSigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(out) :: uvCorr

      ! local variables
      integer :: iCell, kLayer
      real (kind=RKIND) :: sigma
      real (kind=RKIND) :: uppupp, vppvpp

      uvCorr = 0.0_RKIND ! jas issue: assign a mask value instead

      do iCell = 1, nCells
         do kLayer = 1,nLayers
            sigma = max(sigmaEA(kLayer,iCell), epsilonEPFT)
            uvCorr(kLayer, iCell) = uvSigmaEA(kLayer,iCell) / sigma - uTWA(kLayer,iCell)*vTWA(kLayer,iCell)
         enddo
      enddo

   end subroutine calculateCorrelationfromTWA!}}}


!***********************************************************************
!
!  subroutine calculateEPEfromTWA
!
!> \brief   Calculate the eddy potential energy from TWAs
!> \author  Juan A. Saenz
!> \date    July 2015
!> \details
!>  This subroutine calculates the eddy potential energy from thickness
!>  weighted averages.
!-----------------------------------------------------------------------

   subroutine calculateEPEfromTWA(nLayers, nCells, &
         sigmaEA, heightEA, heightSqEA, epe)!{{{
      implicit none
      integer, intent(in) :: nLayers, nCells
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: sigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: heightEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: heightSqEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(out) :: epe

      ! local variables
      integer :: iCell, kLayer
      real (kind=RKIND) :: sigma
      real (kind=RKIND) :: HpHp

      epe = 0.0_RKIND

      do iCell = 1, nCells
         do kLayer = 1,nLayers

            sigma = max(sigmaEA(kLayer,iCell), epsilonEPFT)
            HpHp   = heightSqEA(kLayer,iCell) - heightEA(kLayer,iCell)*heightEA(kLayer,iCell)

            epe(kLayer,iCell) = 0.5_RKIND * HpHp / sigma

         enddo
      enddo

   end subroutine calculateEPEfromTWA!}}}


!***********************************************************************
!
!  subroutine calculateEddyFormDragfromTWA
!
!> \brief   Calculate the eddy form drag from TWAs
!> \author  Juan A. Saenz
!> \date    July 2015
!> \details
!>  This subroutine calculates the eddy form drag from thickness
!>  weighted averages.
!-----------------------------------------------------------------------

   subroutine calculateEddyFormDragfromTWA(nLayers, nCells, &
         sigmaEA, heightEA, MxEA, HMxEA, formDragX)!{{{
      implicit none
      integer, intent(in) :: nLayers, nCells
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: sigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: heightEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: MxEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: HMxEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(out) :: formDragX

      ! local variables
      integer :: iCell, kLayer
      real (kind=RKIND) :: sigma
      real (kind=RKIND) :: HpMxp

       formDragX = 0.0_RKIND

      do iCell = 1, nCells
         do kLayer = 1,nLayers

            sigma = max(sigmaEA(kLayer,iCell), epsilonEPFT)
            HpMxp  = HMxEA(kLayer,iCell) - heightEA(kLayer,iCell) * MxEA(kLayer,iCell)

            formDragX(kLayer,iCell) = HpMxp / sigma
         enddo
      enddo

   end subroutine calculateEddyFormDragfromTWA!}}}


!***********************************************************************
!
!  subroutine calculateDivEPFT
!
!> \brief   Calculate the divergence of EPFT
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details
!>  This subroutine calculates the divergence of the Eliassen-Palm flux tensor.
!>  This is done by calculating the divergence of the column vectors in the tensor.
!>  The divergence of a vector v in buoyancy coordinates is given by (Young 2012):
!>      div of v = sigma^-1 (sigma * v_i)_xi
!
!-----------------------------------------------------------------------

   subroutine calculateDivEPFT(debugFlag, onASphere, rho0, nLayers, nCells, nEdges, &
        meshPool, buoyancyMidRef, sigmaEA, buoyancyMaskEA, tensorCellIn, vectorCellOut)!{{{

      use mpas_vector_operations

      logical, intent(in) :: debugFlag
      logical, intent(in) :: onASphere
      integer, intent(in) :: nLayers, nCells, nEdges
      type (mpas_pool_type), intent(in) :: meshPool
      real (kind=RKIND), intent(in) :: rho0 ! config_density0
      real (kind=RKIND), dimension(:), intent(in) :: buoyancyMidRef
      real (kind=RKIND), dimension(:,:), intent(in) :: sigmaEA
      real (kind=RKIND), dimension(:,:), intent(in) :: buoyancyMaskEA
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: tensorCellIn
      real (kind=RKIND), dimension(:,:,:), intent(out) :: vectorCellOut

      ! local variables
      logical :: includeHalo
      integer :: q, iCell, kLayer, iComponent
      real (kind=RKIND) :: wrk, wrkAbove, wrkBelow, sigma, db
      real (kind=RKIND), dimension(:), pointer :: latCell
      real (kind=RKIND), dimension(:), pointer :: lonCell
      integer, dimension(:,:), pointer :: edgeSignOnCell
      real (kind=RKIND), dimension(:,:), allocatable :: scalarWrk1
      real (kind=RKIND), dimension(:,:,:), allocatable :: vectorCellWrk1
      real (kind=RKIND), dimension(:,:,:), allocatable :: vectorCellWrk2
      real (kind=RKIND), dimension(:,:,:), allocatable :: vectorEdgeWrk1
      real (kind=RKIND), dimension(:), allocatable :: vertVector

      ! variables used for testing and debugging
      real (kind=RKIND), dimension(:), allocatable :: divExact
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      integer, dimension(:,:), pointer :: boundaryCell

      if (debugFlag) then
        allocate(divExact(nCells+1))
      end if

      allocate(scalarWrk1(nLayers,nCells+1))
      allocate(vectorCellWrk1(3,nLayers,nCells+1))
      allocate(vectorCellWrk2(3,nLayers,nCells+1))
      allocate(vectorEdgeWrk1(3,nLayers,nEdges+1))
      allocate(vertVector(nLayers))

      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)

      includeHalo = .true.

      ! initialize work and intent(out)
      vectorCellOut = 0.0_RKIND

      ! loop over all three column vectors
      do q = 1, 3

         scalarWrk1 = 0.0_RKIND

         ! horizontal derivatives
         vectorCellWrk1 = tensorCellIn(:,q,:,:)

         ! weight the vector with sigmaEA(:,:)
         do iComponent = 1,3
            vectorCellWrk1(iComponent,:,:) = sigmaEA(:,:)*vectorCellWrk1(iComponent,:,:)
         enddo


         ! use q=3 as a test vector
         if (q.eq.3 .and. debugFlag) then
           do iCell = 1,nCells
             vectorCellWrk1(1,:,iCell) = xCell(iCell)
             vectorCellWrk1(2,:,iCell) = yCell(iCell)
             vectorCellWrk1(3,:,iCell) = zCell(iCell)
             ! the analytical divergence:
             divExact(iCell) = -1.0_RKIND*sin(lonCell(iCell)) * &
                (1.0_RKIND + 2.0_RKIND*sin(latCell(iCell))) !+ 3.0*sin(latCell(iCell))
           enddo
         endif

         ! zero the vertical component of vectorCellWrk1
         ! the vertical will be treated seperately below
         vectorCellWrk1(3,:,:) = 0.0_RKIND

         if (onASphere) then

            ! convert from lat/lon to Cartesian (x,y,z)
            do iCell = 1,nCells
               do kLayer = 1,nLayers
                  call mpas_vector_LonLatR_to_R3(vectorCellWrk1(:,kLayer,iCell), &
                     lonCell(iCell), latCell(iCell), vectorCellWrk2(:,kLayer,iCell))
               end do
            end do

            ! copy vector measured in x,y,z back into Wrk1
            vectorCellWrk1 = vectorCellWrk2

         end if

         ! average the vector from cell centers to cell edges
         call mpas_vector_R3Cell_to_Edge(vectorCellWrk1, meshPool, &
            vectorEdgeWrk1)

         ! computed the divergence via the weak, line-integral form
         call mpas_divergence_in_r3_buoyancy(vectorEdgeWrk1, meshPool, &
            edgeSignOnCell, includeHalo, scalarWrk1)

         ! use q=3 as a test vector
         if (q.eq.3 .and. debugFlag) then
           print *, ' '
           print *, 'calculateDivEPFT:'
           print *, 'layer, RMS relative error on layer :'
           do kLayer = 1,nLayers
              wrk = sqrt( &
                sum( &
                ( &
                (divExact(:)-scalarWrk1(kLayer,:))/ max(abs(divExact(:)),1.0e-15_RKIND)  &
                )**2 * &
                (1.0_RKIND - boundaryCell(1,:)) &
                ) / nCells )
              print *, kLayer, wrk
           enddo
         endif


         ! Normalize by sigma, after having taken the derivative of sigma * v_i
         if (q < 3 .or. .not. debugFlag) then
            do iCell = 1,nCells
               do kLayer = 1,nLayers
                  sigma = max(sigmaEA(kLayer,iCell), epsilonEPFT)
                  scalarWrk1(kLayer,iCell) = scalarWrk1(kLayer,iCell) / sigma
               end do
            end do
         end if


         ! vertical derivative
         do iCell = 1,nCells

            ! copy the vertical component of EPFT into a work array
            vertVector(:) = tensorCellIn(3,q,:,iCell)

            ! use q=3 as a test vector
            if (q.eq.3 .and. debugFlag) then
               vertVector(:) = 0.0_RKIND
            endif


            do kLayer = 1,nLayers

               wrk = 0.0_RKIND

               if(kLayer.eq.1) then
                 wrkAbove=sigmaEA(kLayer,iCell)*vertVector(kLayer)
                 wrkBelow=sigmaEA(kLayer+1,iCell)*vertVector(kLayer+1)
                 db = buoyancyMidRef(kLayer)-buoyancyMidRef(kLayer+1)
               else if (kLayer.eq.nLayers) then
                 wrkAbove=sigmaEA(kLayer-1,iCell)*vertVector(kLayer-1)
                 wrkBelow=sigmaEA(kLayer,iCell)*vertVector(kLayer)
                 db = buoyancyMidRef(kLayer-1)-buoyancyMidRef(kLayer)
               else
                 wrkAbove=sigmaEA(kLayer-1,iCell)*vertVector(kLayer-1)
                 wrkBelow=sigmaEA(kLayer+1,iCell)*vertVector(kLayer+1)
                 db = buoyancyMidRef(kLayer-1)-buoyancyMidRef(kLayer+1)
               endif

               sigma = max(sigmaEA(kLayer,iCell), epsilonEPFT)
               wrk = (wrkAbove - wrkBelow) / db / sigma

               scalarWrk1(kLayer,iCell) = scalarWrk1(kLayer,iCell) + wrk

            end do ! kLayer = 1,nLayers

         end do ! iCell = 1,nCells

         vectorCellOut(q,:,:) = scalarWrk1

      end do !do q=1,3

      if (debugFlag) then
        deallocate(divExact)
      end if
      deallocate(scalarWrk1)
      deallocate(vectorCellWrk1)
      deallocate(vectorCellWrk2)
      deallocate(vectorEdgeWrk1)
      deallocate(vertVector)

   end subroutine calculateDivEPFT!}}}



!***********************************************************************
!
!  subroutine calculateErtelPVFlux
!
!> \brief   Calculate the Ertel potential vorticity fluxes
!> \author  Juan A. Saenz
!> \date    January 2014
!> \details
!>  This subroutine calculates the Ertel potential vorticity fluxes
!>  using the divergence of EPFT, as outlined in eqn 129 of Young 2012.
!
!-----------------------------------------------------------------------

   subroutine calculateErtelPVFlux(nCells, nBuoyancyLayers, &
        sigma, divEPFT, ErtelPVFlux)!{{{

      integer, intent(in) :: nCells, nBuoyancyLayers
      real (kind=RKIND), dimension(:,:), intent(in) :: sigma
      real (kind=RKIND), dimension(:,:,:), intent(in) :: divEPFT
      real (kind=RKIND), dimension(:,:,:), intent(out) :: ErtelPVFlux

      ! local variables
      integer :: i, k

      ErtelPVFlux(1,:,:) = divEPFT(2,:,:)
      ErtelPVFlux(2,:,:) = -1.0_RKIND * divEPFT(1,:,:)
      ErtelPVFlux(3,:,:) = 0.0_RKIND

      do i = 1, nCells
         do k = 1,nBuoyancyLayers
            ErtelPVFlux(:,k,i) = ErtelPVFlux(:,k,i) / max(sigma(k,i),epsilonEPFT)
         end do
      end do

    end subroutine calculateErtelPVFlux


!***********************************************************************
!
!  subroutine calculateErtelPVTendencyFromPVFlux
!
!> \brief   Calculate the Ertel PV tendency from Ertel PV flux
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details
!>  This subroutine calculates the Ertel PV tendency as the divergence of
!>  the Ertel PV flux, where the latter only has horizontal components.
!
!-----------------------------------------------------------------------

   subroutine calculateErtelPVTendencyFromPVFlux(debugFlag, onASphere, nLayers, nCells, nEdges, &
         meshPool, sigma, vectorCell, divVectorCell)!{{{

      use mpas_vector_operations

      logical, intent(in) :: debugFlag
      logical, intent(in) :: onASphere
      integer, intent(in) :: nLayers, nCells, nEdges
      type (mpas_pool_type), intent(in) :: meshPool
      real (kind=RKIND), dimension(:,:), intent(in) :: sigma
      real (kind=RKIND), dimension(:,:,:), intent(in) :: vectorCell
      real (kind=RKIND), dimension(:,:), intent(out) :: divVectorCell

      ! local variables
      logical :: includeHalo
      integer :: i, k, iComponent
      real (kind=RKIND), dimension(:), pointer :: latCell
      real (kind=RKIND), dimension(:), pointer :: lonCell
      integer, dimension(:,:), pointer :: edgeSignOnCell, boundaryCell
      real (kind=RKIND), dimension(:,:,:), allocatable :: vectorCellWrk1
      real (kind=RKIND), dimension(:,:,:), allocatable :: vectorCellWrk2
      real (kind=RKIND), dimension(:,:,:), allocatable :: vectorEdgeWrk1

      ! test variables
      real (kind=RKIND) :: wrk
      real (kind=RKIND), dimension(:), allocatable :: divExact
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell

      if (debugFlag) then
        allocate(divExact(nCells+1))
      end if


      allocate(vectorCellWrk1(3,nLayers,nCells+1))
      allocate(vectorCellWrk2(3,nLayers,nCells+1))
      allocate(vectorEdgeWrk1(3,nLayers,nEdges)) !jas issue Todd had set this to nEdges+1

      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)

      includeHalo = .true.

      divVectorCell = 0.0_RKIND

      ! copy vector into work array
      vectorCellWrk1 = vectorCell

      ! weight the vector with sigma(:,:)
      do iComponent = 1,3
         vectorCellWrk1(iComponent,:,:) = sigma(:,:)*vectorCellWrk1(iComponent,:,:)
      enddo


      if (debugFlag) then
        do i= 1,nCells
          vectorCellWrk1(1,:,i) = xCell(i)
          vectorCellWrk1(2,:,i) = yCell(i)
          vectorCellWrk1(3,:,i) = zCell(i)
          ! the analytical divergence:
          divExact(i) = -1.0_RKIND*sin(lonCell(i)) * &
             (1.0_RKIND + 2.0_RKIND*sin(latCell(i))) !+ 3.0*sin(latCell(i))
        enddo
      endif


      if (onASphere) then
      ! copy vector into work array
         do i = 1,nCells
            do k = 1,nLayers
               call mpas_vector_LonLatR_to_R3(vectorCellWrk1(:,k,i), &
                  lonCell(i), latCell(i), vectorCellWrk2(:,k,i))
            end do
         end do

         ! copy transformed vector back to Wrk1
         vectorCellWrk1 = vectorCellWrk2
      end if

      ! average vector from cell centers to cell edges
      call mpas_vector_R3Cell_to_Edge(vectorCellWrk1, meshPool, &
         vectorEdgeWrk1)

      ! computed divergence via weak-form, line integral
      call mpas_divergence_in_r3_buoyancy(vectorEdgeWrk1, meshPool, edgeSignOnCell, &
         includeHalo, divVectorCell)


      if (debugFlag) then
        print *, ' '
        print *, 'calculateErtelPVTendencyFromPVFlux:'
        print *, 'k, RMS relative error on layer:'
        do k= 1,nLayers
           wrk = sqrt( &
             sum( &
             ( &
             (divExact(:)-divVectorCell(k,:))/ max(abs(divExact(:)),1.0e-15_RKIND)  &
             )**2 * &
             (1.0_RKIND - boundaryCell(1,:)) &
             ) / nCells )
           print *, k, wrk
        enddo
      endif


      if (.not. debugFlag) then
         do i = 1,nCells
            do k = 1,nLayers
              divVectorCell(k,i) = divVectorCell(k,i) / max(sigma(k,i), epsilonEPFT)
            end do
         end do
      end if

      if (debugFlag) then
        deallocate(divExact)
      end if
      deallocate(vectorCellWrk1)
      deallocate(vectorCellWrk2)
      deallocate(vectorEdgeWrk1)

   end subroutine calculateErtelPVTendencyFromPVFlux!}}}



!***********************************************************************
!
!  subroutine computeErtelPV
!
!> \brief   Calculate Ertel potential vorticity on buoyancy surfaces
!> \author  Juan A. Saenz
!> \date    January 2014
!> \details
!>  This subroutine calculates Ertel potential voriticity in buoyancy surfaces
!
!-----------------------------------------------------------------------

   subroutine computeErtelPV(nCells, nLayers, nEdges, meshPool, &
      fCell, uCell, vCell, sigma, ErtelPV)

      use mpas_vector_reconstruction

      integer, intent(in) :: nCells, nLayers, nEdges
      type (mpas_pool_type), intent(in) :: meshPool
      real (kind=RKIND), dimension(:), intent(in) :: fCell
      real (kind=RKIND), dimension(:,:), intent(in) :: uCell, vCell
      real (kind=RKIND), dimension(:,:), intent(in) :: sigma
      real (kind=RKIND), dimension(:,:), intent(out) :: ErtelPV

      ! local variables
      integer :: i, k
      real (kind=RKIND), dimension(:,:), allocatable :: velNormalGradOnEdge
      real (kind=RKIND), dimension(:,:), allocatable :: velGradX, velGradY, velGradZ
      real (kind=RKIND), dimension(:,:), allocatable :: velGradZonal, velGradMerid
      real (kind=RKIND), dimension(:,:), allocatable :: vGradZonal, uGradMerid

      allocate(velNormalGradOnEdge(nLayers, nEdges)) ! jas issue this one seems correct
      allocate(velGradX(nLayers, nCells)) ! jas issue gets sent to routine where looping over nCellsSolve
      allocate(velGradY(nLayers, nCells)) ! jas issue gets sent to routine where looping over nCellsSolve
      allocate(velGradZ(nLayers, nCells)) ! jas issue gets sent to routine where looping over nCellsSolve
      allocate(velGradZonal(nLayers, nCells)) ! jas issue gets sent to routine where looping over nCellsSolve
      allocate(velGradMerid(nLayers, nCells)) ! jas issue gets sent to routine where looping over nCellsSolve
      allocate(uGradMerid(nLayers, nCells)) ! jas issue gets assigned array of size is nCells
      allocate(vGradZonal(nLayers, nCells)) ! jas issue gets assigned array of size is nCells

      ! calculate derivative of uTWA with respect to the meridional direction
      call computeNormalGradientOnEdge(nLayers, nCells, nEdges, &
        meshPool, uCell, velNormalGradOnEdge)
      call mpas_reconstruct(meshPool, velNormalGradOnEdge, &
        velGradX, velGradY, velGradZ, &
        velGradZonal, velGradMerid, includeHalos=.true.)
      uGradMerid = velGradMerid

      ! calculate derivative of vTWA with respect to the zonal direction
      call computeNormalGradientOnEdge(nLayers, nCells, nEdges, &
        meshPool, vCell, velNormalGradOnEdge)
      call mpas_reconstruct(meshPool, velNormalGradOnEdge, &
        velGradX, velGradY, velGradZ, &
        velGradZonal, velGradMerid, includeHalos=.true.)
      vGradZonal = velGradZonal

      ErtelPV = 0.0_RKIND

      do i = 1, nCells
        do k = 1,nLayers
          ErtelPV(k,i) = (fCell(i) + vGradZonal(k,i) - uGradMerid(k,i))/max(sigma(k,i),epsilonEPFT)
        end do
      end do

      deallocate(velNormalGradOnEdge)
      deallocate(velGradX)
      deallocate(velGradY)
      deallocate(velGradZ)
      deallocate(velGradZonal)
      deallocate(velGradMerid)
      deallocate(vGradZonal)
      deallocate(uGradMerid)


   end subroutine computeErtelPV


!***********************************************************************
!
!  subroutine computeVerticalDerivative
!
!> \brief   Calculate the the vertical derivative, in depth coordinates, of a scalar
!> \author  Juan A. Saenz
!> \date    July, 2015
!> \details
!>  This subroutine calculates the vertical derivative, in depth coordinates, of a scalar.
!>  The scalar is assumed to exist in the middle of a cell layer.
!>  The vertical derivative in the middle of the cell layer is returned.
!
!-----------------------------------------------------------------------

   subroutine computeVerticalDerivative(nCells, nLayers, &
       firstLayer, lastLayer, heightMid, field, derivativeField)!{{{
      integer, intent(in) :: nCells, nLayers
      integer, dimension(nCells), intent(in)  :: firstLayer, lastLayer
      real (kind=RKIND), dimension(:,:), intent(in) :: heightMid
      real (kind=RKIND), dimension(:,:), intent(in) :: field
      real (kind=RKIND), dimension(:,:), intent(out) :: derivativeField

      ! local variables
      integer :: iCell, kLayer
      real (kind=RKIND) :: wrkAbove, wrkBelow, dz

      derivativeField(nLayers, 1:nCells) = 0.0_RKIND

      do iCell = 1,nCells

         if ( lastLayer(iCell) > firstLayer(iCell) ) then
            wrkAbove = field(firstLayer(iCell),iCell)
            wrkBelow = field(firstLayer(iCell)+1,iCell)
            dz = heightMid(firstLayer(iCell),iCell)-heightMid(firstLayer(iCell)+1,iCell)

            derivativeField(firstLayer(iCell), iCell) = (wrkAbove - wrkBelow) / dz

            do kLayer = firstLayer(iCell)+1, lastLayer(iCell)-1

               wrkAbove = field(kLayer-1,iCell)
               wrkBelow = field(kLayer+1,iCell)
               dz = heightMid(kLayer-1,iCell)-heightMid(kLayer+1,iCell)

               derivativeField(kLayer, iCell) = (wrkAbove - wrkBelow) / dz

            end do ! kLayer = firstLayer(iCell)+1, lastLayer(iCell)-1

            wrkAbove = field(lastLayer(iCell)-1,iCell)
            wrkBelow = field(lastLayer(iCell),iCell)
            dz = heightMid(lastLayer(iCell)-1,iCell)-heightMid(lastLayer(iCell),iCell)

            derivativeField(lastLayer(iCell), iCell) = (wrkAbove - wrkBelow) / dz
         end if
      end do ! iCell = 1,nCells
   end subroutine computeVerticalDerivative!}}}


!***********************************************************************
!
!  subroutine eddyGeomDecompEPFT
!
!> \brief   Calculate the eddy geometric decomposition from EPFT
!> \author  Juan A. Saenz
!> \date    January 2014
!> \details
!>  This subroutine calculates the eddy geometric decomposition from EPFT
!
!-----------------------------------------------------------------------

   subroutine eddyGeomDecompEPFT()!sigmaRef, ErtelPVFlux, ErtelPVTendency)!{{{
         ! Compute the geometric decomposition in terms of angles and eccentricities using
         ! the entries of EPFT.

   end subroutine eddyGeomDecompEPFT!}}}


!***********************************************************************
!
!  routine mpas_divergence_in_r3_buoyancy
!
!> \brief   MPAS 3D divergence routine
!> \author  Todd Ringler
!> \date    02/07/14
!> \details
!> This routine computes the of an input vector.
!-----------------------------------------------------------------------
   subroutine mpas_divergence_in_r3_buoyancy(vectorR3Edge, meshPool, &
      edgeSignOnCell, includeHalo, divCell)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         vectorR3Edge  !< Input: vector at edge, R3, indices (direction,verticalIndex,edgeIndex)

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      integer, dimension(:,:), intent(in) :: &
         edgeSignOnCell        !< Input: Direction of vector connecting cells

      logical, intent(in) :: &
         includeHalo !< Input: If true, halo cells and edges are included in computation

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         divCell          !< Output: scalar divergence, indices (verticalIndex,edgeIndex)

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, iCell, i, k, p
      integer, pointer :: nVertLevels, nCells

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND) :: edgeNormalDotVector
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors

      call mpas_pool_get_dimension(meshPool, 'nBuoyancyLayers', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors)

      divCell(:,:) = 0.0_RKIND
      do iCell = 1, nCells
         invAreaCell = 1.0_RKIND / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, nVertLevels
               edgeNormalDotVector = 0.0_RKIND
               do p=1,3
                 edgeNormalDotVector = edgeNormalDotVector + &
                    edgeNormalVectors(p,iEdge)*vectorR3Edge(p,k,iEdge)
                enddo
               divCell(k,iCell) = divCell(k,iCell) - &
                  edgeSignOnCell(i,iCell) * dvEdge(iEdge) * invAreaCell * &
                  edgeNormalDotVector
            end do
         end do
      end do

  end subroutine mpas_divergence_in_r3_buoyancy!}}}


!***********************************************************************
!
!  routine mpas_vector_R3Cell_to_Edge
!
!> \brief   MPAS 3D divergence routine
!> \author  Todd Ringler
!> \date    02/07/14
!> \details
!> This routine averages a vector field from cells to edges
!-----------------------------------------------------------------------
  subroutine mpas_vector_R3Cell_to_Edge(vectorCell, meshPool, &
        vectorEdge)

      real (kind=RKIND), dimension(:,:,:), intent(in) :: vectorCell
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      real (kind=RKIND), dimension(:,:,:), intent(out) :: vectorEdge

      !local variables
      integer :: iEdge, k, cell1, cell2
      integer, pointer :: nEdges, nBuoyancyLayers
      integer, dimension(:,:), pointer :: cellsOnEdge, boundaryEdge

      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'boundaryEdge', boundaryEdge)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nBuoyancyLayers', nBuoyancyLayers)

      vectorEdge = 0.0_RKIND

      do iEdge = 1, nEdges
         ! Enforce vector value of zero on boundary edges, e.g. no slip for velocities
         if (boundaryEdge(1,iEdge) == 1) then
            vectorEdge(:,:,iEdge) = 0.0_RKIND
         else
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            do k = 1, nBuoyancyLayers
                vectorEdge(:,k,iEdge) = 0.5_RKIND*( vectorCell(:,k,cell2) + vectorCell(:,k,cell1) )
            enddo
         end if
      enddo

   end subroutine mpas_vector_R3Cell_to_Edge!}}}


end module ocn_eliassen_palm

! vim: foldmethod=marker
