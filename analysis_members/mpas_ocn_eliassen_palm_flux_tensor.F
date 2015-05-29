! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_epft
!
!> \brief MPAS ocean analysis core member: epft
!> \author Juan A. Saenz, Todd Ringler
!> \date   May, 2015
!> \details
!>  This module contains the routines for computing the Eliassen and Palm Flux Tensor 
!>  in buoyancy coordinates, and related quantities.
!
!-----------------------------------------------------------------------

module ocn_eliassen_palm_flux_tensor

   use mpas_grid_types
   use mpas_timer
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use mpas_configure
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

   public :: ocn_setup_packages_eliassen_palm_flux_tensor, &
             ocn_init_eliassen_palm_flux_tensor, &
             ocn_compute_eliassen_palm_flux_tensor, &
             ocn_restart_eliassen_palm_flux_tensor, &
             ocn_finalize_eliassen_palm_flux_tensor

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: am_eliassen_palm_flux_tensorTimer
   logical :: amEPFTOn
   real (kind=RKIND), parameter :: epsilonEPFT=1.0e-15

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_setup_packages_eliassen_palm_flux_tensor
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Juan Saenz, Todd Ringler
!> \date    May 2015
!> \details 
!>  This routine is intended to configure the packages for this MPAS
!>   ocean analysis member
!
!-----------------------------------------------------------------------

   subroutine ocn_setup_packages_eliassen_palm_flux_tensor(configPool, packagePool, err)!{{{

      use mpas_packages

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

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
      logical, pointer :: am_eliassen_palm_flux_tensor_Active

      err = 0

      call mpas_pool_get_config(configPool, "config_use_epft", config_use_epft)
      call mpas_pool_get_package(packagePool, &
          'am_eliassen_palm_flux_tensor_Active', am_eliassen_palm_flux_tensor_Active)

      ! turn on package for this analysis member
      am_eliassen_palm_flux_tensor_Active = .true.

   end subroutine ocn_setup_packages_eliassen_palm_flux_tensor!}}}


!***********************************************************************
!
!  routine ocn_init_eliassen_palm_flux_tensor
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_eliassen_palm_flux_tensor(domain, err)!{{{

      use mpas_packages

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

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
      !real(KIND=RKIND) :: global_min, global_max, local_min, local_max

      type (block_type), pointer :: block
      type (amEPFT_type), pointer :: amEPFT

      real(KIND=RKIND), dimension(:), pointer :: potentialDensityMidRef
      real(KIND=RKIND), dimension(:), pointer :: potentialDensityTopRef
      real(KIND=RKIND), dimension(:,:), pointer :: buoyMaskEA



      err = 0

      block => domain % blocklist
      do while (associated(block))

         amEPFT => block % amEPFT

         ! Calculate target values
         potentialDensityMidRef => amEPFT % potentialDensityMidRef % array
         potentialDensityTopRef => amEPFT % potentialDensityTopRef % array

         do k = 1, config_nBuoyancyLayers
            potentialDensityTopRef(k) = config_rhomin_buoycoor + &
               (config_rhomax_buoycoor - config_rhomin_buoycoor) / &
               (config_nBuoyancyLayers) * (k-1)
         end do
         do k = 1, config_nBuoyancyLayers-1
            potentialDensityMidRef(k) = &
               0.5*(potentialDensityTopRef(k) + potentialDensityTopRef(k+1))
         end do
         potentialDensityMidRef(config_nBuoyancyLayers) = &
            0.5*(potentialDensityTopRef(config_nBuoyancyLayers) + config_rhomax_buoycoor)

         if (.not. config_do_restart .or. config_oac_epft_reset) then
            amEPFT % buoyMaskEA % array = 0.0
            amEPFT % sigmaEA % array = 0.0
            amEPFT % nSamplesEA % scalar = 0.0
            amEPFT % heightMidBuoyCoorEA % array = 0.0
            amEPFT % montgPotBuoyCoorEA % array = 0.0
            amEPFT % montgPotGradZonalEA % array = 0.0
            amEPFT % montgPotGradMeridEA % array = 0.0
            amEPFT % heightMidBuoyCoorSqEA % array = 0.0
            amEPFT % HeightMGradZonalEA % array = 0.0
            amEPFT % HeightMGradMeridEA % array = 0.0
            amEPFT % usigmaEA % array = 0.0
            amEPFT % vsigmaEA % array = 0.0
            amEPFT % uusigmaEA % array = 0.0
            amEPFT % vvsigmaEA % array = 0.0
            amEPFT % uvsigmaEA % array = 0.0
            amEPFT % uwsigmaEA % array = 0.0
            amEPFT % vwsigmaEA % array = 0.0
         end if

         block => block % next

      end do


   end subroutine ocn_init_eliassen_palm_flux_tensor!}}}

!***********************************************************************
!
!  routine ocn_compute_eliassen_palm_flux_tensor
!
!> \brief   Compute Eliassen-Palm flux tensor
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details 
!>  This routine conducts all computation required for the EPFT analysis member.
!>  Each time this AM is called, the instananeous ocean state is interpolated
!>  onto the target buoyancy values. The state is then accumulated in the
!>  accumulated into the ensemble average (*EA) arrays. Based on the current
!>  estimate of the ensemble average, thickness-weight velocities are estimates
!>  along with the computation of the Eliassen-Palm flux tensor
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_eliassen_palm_flux_tensor(domain, timeLevel, err)!{{{

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
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: forcingPool ! jas-issue does this exist?
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool 
      type (mpas_pool_type), pointer :: am_epftPool
      type (mpas_pool_type), pointer :: am_eliassen_palm_flux_tensor ! jas issue ?

      logical, pointer :: config_epft_debug

      !-----------------------------------------------------------------
      ! define local scalars holding length of dimensions
      !-----------------------------------------------------------------
      integer, pointer :: nVertLevels, nBuoyLayers, nBuoyLayersP1
      integer, pointer :: nEdges, nCells, nCellsSolve ! nCellsSolve includes halos
      

      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:), pointer :: firstLayerBuoyCoor
      integer, dimension(:), pointer :: lastLayerBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: buoyMask

      integer :: nSamplesEA
      
      real(KIND=RKIND), dimension(:), pointer :: potentialDensityMidRef
      real(KIND=RKIND), dimension(:), pointer :: potentialDensityTopRef
      real(KIND=RKIND), dimension(:), pointer :: buoyancyMidRef
      real(KIND=RKIND), dimension(:), pointer :: buoyancyInterfaceRef
      real(KIND=RKIND), dimension(:,:), pointer :: buoyMaskEA
      real(KIND=RKIND), dimension(:,:), pointer :: sigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: heightMidBuoyCoorEA
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradZonalEA
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradMeridEA
      real(KIND=RKIND), dimension(:,:), pointer :: heightMidBuoyCoorSqEA
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotBuoyCoorEA
      real(KIND=RKIND), dimension(:,:), pointer :: HeightMGradZonalEA
      real(KIND=RKIND), dimension(:,:), pointer :: HeightMGradMeridEA
      real(KIND=RKIND), dimension(:,:), pointer :: usigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: vsigmaEA
      !real(KIND=RKIND), dimension(:,:), pointer :: wsigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: uusigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: vvsigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: uvsigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: uwsigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: vwsigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: uTWA
      real(KIND=RKIND), dimension(:,:), pointer :: vTWA
      real(KIND=RKIND), dimension(:,:), pointer :: wTWA
      real(KIND=RKIND), dimension(:,:,:,:), pointer :: EPFT
      real(KIND=RKIND), dimension(:,:,:), pointer :: divEPFT 
      real(KIND=RKIND), dimension(:,:,:), pointer :: ErtelPVFlux
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPVTendency
      real(KIND=RKIND), dimension(:,:), pointer :: ErtelPV


      real(KIND=RKIND), dimension(:), pointer :: SSH

      real(KIND=RKIND), dimension(:), pointer :: bottomDepth

      real(KIND=RKIND), dimension(:,:), pointer :: heightMidBuoyCoor    
      real(KIND=RKIND), dimension(:,:), pointer :: heightTopBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: heightInterfaceBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: uMidBuoyCoor    
      real(KIND=RKIND), dimension(:,:), pointer :: vMidBuoyCoor    
      real(KIND=RKIND), dimension(:,:), pointer :: densityMidBuoyCoor    
      real(KIND=RKIND), dimension(:,:), pointer :: densityTopBuoyCoor    
      real(KIND=RKIND), dimension(:,:), pointer :: sigma
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotNormalGradOnEdge
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradX         
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradY         
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradZ         
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradZonal     
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradMerid

      real(KIND=RKIND), dimension(:), pointer :: surfacePressure
      real(KIND=RKIND), dimension(:,:), pointer :: zMid
      real(KIND=RKIND), dimension(:,:), pointer :: zTop
      real(KIND=RKIND), dimension(:,:), pointer :: density
      real(KIND=RKIND), dimension(:,:), pointer :: potentialDensity
      real(KIND=RKIND), dimension(:,:), pointer :: pressure
      real(KIND=RKIND), dimension(:,:), pointer :: uCellCenter
      real(KIND=RKIND), dimension(:,:), pointer :: vCellCenter
!      real(KIND=RKIND), dimension(:,:) :: wCellCenter

      ! work variables 
      integer :: k
      real(KIND=RKIND), dimension(:,:), pointer :: wrk3DnVertLevels
      real(KIND=RKIND), dimension(:,:), pointer :: wrk3DBuoyCoor

      ! test variables
      integer :: nCellsGlobal, i
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
      real(KIND=RKIND), dimension(:,:), pointer :: relativeVorticityCell
      real(KIND=RKIND), dimension(:), pointer :: fCell

      err = 0

      dminfo = domain % dminfo

      call mpas_timer_start("compute_eliassen_palm_flux_tensor", .false., &
        am_eliassen_palm_flux_tensorTimer)
      
      call mpas_pool_get_config(domain % configs, 'config_epft_debug', config_epft_debug)

      
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'eliassenPalmFluxTensorScratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'amEliassenPalmFluxTensor', am_epftPool)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nBuoyLayers', nBuoyLayers)
         !call mpas_pool_get_dimension(block % dimensions, 'nBuoyLayersP1', nBuoyLayersP1)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(forcingPool, 'seaSurfacePressure', surfacePressure)

         call mpas_pool_get_field(scratchPool, 'heightMidBuoyCoor', heightMidBuoyCoor)
         call mpas_pool_get_field(scratchPool, 'heightTopBuoyCoor', heightTopBuoyCoor)
         call mpas_pool_get_field(scratchPool, 'heightInterfaceBuoyCoor', heightInterfaceBuoyCoor)
         call mpas_pool_get_field(scratchPool, 'uMidBuoyCoor', uMidBuoyCoor)
         call mpas_pool_get_field(scratchPool, 'vMidBuoyCoor', vMidBuoyCoor)
         call mpas_pool_get_field(scratchPool, 'densityMidBuoyCoor', densityMidBuoyCoor)
         call mpas_pool_get_field(scratchPool, 'densityTopBuoyCoor', densityTopBuoyCoor)
         call mpas_pool_get_field(scratchPool, 'sigma', sigma)
         call mpas_pool_get_field(scratchPool, 'montgPotBuoyCoor', montgPotBuoyCoor)
         call mpas_pool_get_field(scratchPool, 'montgPotNormalGradOnEdge', montgPotNormalGradOnEdge)
         call mpas_pool_get_field(scratchPool, 'firstLayerBuoyCoor', firstLayerBuoyCoor)
         call mpas_pool_get_field(scratchPool, 'lastLayerBuoyCoor', lastLayerBuoyCoor)
         call mpas_pool_get_field(scratchPool, 'buoyMask', buoyMask)
         call mpas_pool_get_field(scratchPool, 'montgPotGradX', montgPotGradX)
         call mpas_pool_get_field(scratchPool, 'montgPotGradY', montgPotGradY)
         call mpas_pool_get_field(scratchPool, 'montgPotGradZ', montgPotGradZ)
         call mpas_pool_get_field(scratchPool, 'montgPotGradZonal', montgPotGradZonal)
         call mpas_pool_get_field(scratchPool, 'montgPotGradMerid', montgPotGradMerid)
         call mpas_pool_get_field(scratchPool, 'wrk3DnVertLevelsP1', wrk3DnVertLevelsP1)
         call mpas_pool_get_field(scratchPool, 'wrk3DnVertLevels', wrk3DnVertLevels)
         call mpas_pool_get_field(scratchPool, 'wrk3DBuoyCoor', wrk3DBuoyCoor)
         
         call mpas_allocate_scratch_field(heightMidBuoyCoor, .true.)
         call mpas_allocate_scratch_field(heightTopBuoyCoor, .true.)
         call mpas_allocate_scratch_field(heightInterfaceBuoyCoor, .true.)
         call mpas_allocate_scratch_field(uMidBuoyCoor, .true.)
         call mpas_allocate_scratch_field(vMidBuoyCoor, .true.)
         call mpas_allocate_scratch_field(densityMidBuoyCoor, .true.)
         call mpas_allocate_scratch_field(densityTopBuoyCoor, .true.)
         call mpas_allocate_scratch_field(sigma, .true.)
         call mpas_allocate_scratch_field(montgPotBuoyCoor, .true.)
         call mpas_allocate_scratch_field(montgPotNormalGradOnEdge, .true.)
         call mpas_allocate_scratch_field(firstLayerBuoyCoor, .true.)
         call mpas_allocate_scratch_field(lastLayerBuoyCoor, .true.)
         call mpas_allocate_scratch_field(buoyMask, .true.)
         call mpas_allocate_scratch_field(montgPotGradX, .true.)
         call mpas_allocate_scratch_field(montgPotGradY, .true.)
         call mpas_allocate_scratch_field(montgPotGradZ, .true.)
         call mpas_allocate_scratch_field(montgPotGradZonal, .true.)
         call mpas_allocate_scratch_field(montgPotGradMerid, .true.)
         call mpas_allocate_scratch_field(wrk3DnVertLevelsP1, .true.)
         call mpas_allocate_scratch_field(wrk3DnVertLevels, .true.)
         call mpas_allocate_scratch_field(wrk3DBuoyCoor, .true.)

         ! test variables
         call mpas_pool_get_field(scratchPool, 'array1_3D', array1_3D)
         call mpas_pool_get_field(scratchPool, 'array2_3D', array2_3D)
         call mpas_pool_get_field(scratchPool, 'array3_3D', array3_3D)
         call mpas_pool_get_field(scratchPool, 'array1_3Dbuoy', array1_3Dbuoy)
         call mpas_pool_get_field(scratchPool, 'array2_3Dbuoy', array2_3Dbuoy)

         call mpas_allocate_scratch_field(array1_3D, .true.)
         call mpas_allocate_scratch_field(array2_3D, .true.)
         call mpas_allocate_scratch_field(array3_3D, .true.)
         call mpas_allocate_scratch_field(array1_3Dbuoy, .true.)
         call mpas_allocate_scratch_field(array2_3Dbuoy, .true.)
        
         call mpas_pool_get_field('PVMidBuoyCoor', PVMidBuoyCoor)
         call mpas_pool_get_field('PVMidBuoyCoorEA', PVMidBuoyCoorEA)
         call mpas_pool_get_field('uMidBuoyCoorEA', uMidBuoyCoorEA)
         call mpas_pool_get_field('vMidBuoyCoorEA', vMidBuoyCoorEA)
         call mpas_pool_get_field('uPVMidBuoyCoorEA', uPVMidBuoyCoorEA)
         call mpas_pool_get_field('vPVMidBuoyCoorEA', vPVMidBuoyCoorEA)
         call mpas_pool_get_field('PVFluxTest', PVFluxTest)

         call mpas_allocate_scratch_field(PVMidBuoyCoor, .true.)
         call mpas_allocate_scratch_field(PVMidBuoyCoorEA, .true.)
         call mpas_allocate_scratch_field(uMidBuoyCoorEA , .true.)
         call mpas_allocate_scratch_field(vMidBuoyCoorEA , .true.)
         call mpas_allocate_scratch_field(uPVMidBuoyCoorEA , .true.)
         call mpas_allocate_scratch_field(vPVMidBuoyCoorEA, .true.)
         call mpas_allocate_scratch_field(PVFluxTest, .true.)

         call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
         call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
         call mpas_pool_get_array(diagnosticsPool, 'density', density)
         call mpas_pool_get_array(diagnosticsPool, 'potentialdensity', potentialDensity)
         call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)
         call mpas_pool_get_array(diagnosticsPool, 'normalVelocityZonal', uCellCenter)
         call mpas_pool_get_array(diagnosticsPool, 'normalVelocityMeridional', vCellCenter)
         
         call mpas_pool_get_array(am_epftPool, 'potentialDensityMidRef', potentialDensityMidRef)
         call mpas_pool_get_array(am_epftPool, 'potentialDensityTopRef', potentialDensityTopRef)
         call mpas_pool_get_array(am_epftPool, 'buoyancyMidRef', buoyancyMidRef)
         call mpas_pool_get_array(am_epftPool, 'buoyancyInterfaceRef', buoyancyInterfaceRef)
         call mpas_pool_get_array(am_epftPool, 'buoyMaskEA', buoyMaskEA)
         call mpas_pool_get_array(am_epftPool, 'sigmaEA', sigmaEA)
         call mpas_pool_get_array(am_epftPool, 'heightMidBuoyCoorEA', heightMidBuoyCoorEA)
         call mpas_pool_get_array(am_epftPool, 'montgPotGradZonalEA', montgPotGradZonalEA)
         call mpas_pool_get_array(am_epftPool, 'montgPotGradMeridEA', montgPotGradMeridEA)
         call mpas_pool_get_array(am_epftPool, 'heightMidBuoyCoorSqEA', heightMidBuoyCoorSqEA)
         call mpas_pool_get_array(am_epftPool, 'montgPotBuoyCoorEA', montgPotBuoyCoorEA)
         call mpas_pool_get_array(am_epftPool, 'HeightMGradZonalEA', HeightMGradZonalEA)
         call mpas_pool_get_array(am_epftPool, 'HeightMGradMeridEA', HeightMGradMeridEA)
         call mpas_pool_get_array(am_epftPool, 'usigmaEA', usigmaEA)
         call mpas_pool_get_array(am_epftPool, 'vsigmaEA', vsigmaEA)
         call mpas_pool_get_array(am_epftPool, 'uusigmaEA', uusigmaEA)
         call mpas_pool_get_array(am_epftPool, 'vvsigmaEA', vvsigmaEA)
         call mpas_pool_get_array(am_epftPool, 'uvsigmaEA', uvsigmaEA)
         call mpas_pool_get_array(am_epftPool, 'uwsigmaEA', uwsigmaEA)
         call mpas_pool_get_array(am_epftPool, 'vwsigmaEA', vwsigmaEA)
         call mpas_pool_get_array(am_epftPool, 'uTWA', uTWA)
         call mpas_pool_get_array(am_epftPool, 'vTWA', vTWA)
         call mpas_pool_get_array(am_epftPool, 'wTWA', wTWA)
         call mpas_pool_get_array(am_epftPool, 'EPFT', EPFT)
         call mpas_pool_get_array(am_epftPool, 'divEPFT', divEPFT)
         call mpas_pool_get_array(am_epftPool, 'ErtelPVFlux', ErtelPVFlux)
         call mpas_pool_get_array(am_epftPool, 'ErtelPVTendency', ErtelPVTendency)
         call mpas_pool_get_array(am_epftPool, 'ErtelPV', ErtelPV)
        
         call mpas_pool_get_array(statePool, 'SSH', SSH)
         
         call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

         nSamplesEA            =  nSamplesEA % scalar
         
         heightMidBuoyCoor       => heightMidBuoyCoor % array
         heightTopBuoyCoor       => heightTopBuoyCoor % array
         heightInterfaceBuoyCoor => heightInterfaceBuoyCoor % array
         uMidBuoyCoor            => uMidBuoyCoor % array
         vMidBuoyCoor            => vMidBuoyCoor % array
         densityMidBuoyCoor      => densityMidBuoyCoor % array
         densityTopBuoyCoor      => densityTopBuoyCoor % array
         sigma                   => sigma % array
         montgPotBuoyCoor        => montgPotBuoyCoor % array
         montgPotNormalGradOnEdge=> montgPotNormalGradOnEdge % array
         firstLayerBuoyCoor      => firstLayerBuoyCoor % array
         lastLayerBuoyCoor       => lastLayerBuoyCoor % array
         buoyMask                => buoyMask % array
         montgPotGradX           => montgPotGradX % array
         montgPotGradY           => montgPotGradY % array
         montgPotGradZ           => montgPotGradZ % array
         montgPotGradZonal       => montgPotGradZonal % array
         montgPotGradMerid       => montgPotGradMerid % array
         wrk3DnVertLevelsP1      => wrk3DnVertLevelsP1 % array
         wrk3DnVertLevels        => wrk3DnVertLevels % array
         wrk3DBuoyCoor           => wrk3DBuoyCoor % array
         
         array1_3D =>  array1_3D % array
         array2_3D =>  array2_3D % array
         array3_3D =>  array3_3D % array
         array1_3Dbuoy =>  array1_3Dbuoy % array
         array2_3Dbuoy =>  array2_3Dbuoy % array

         PVMidBuoyCoor    =>  PVMidBuoyCoor % array
         PVMidBuoyCoorEA  =>  PVMidBuoyCoorEA % array
         uMidBuoyCoorEA   =>  uMidBuoyCoorEA % array
         vMidBuoyCoorEA   =>  vMidBuoyCoorEA % array
         uPVMidBuoyCoorEA =>  uPVMidBuoyCoorEA % array
         vPVMidBuoyCoorEA =>  vPVMidBuoyCoorEA % array
         PVFluxTest       =>  PVFluxTest % array

         
         nBuoyLayersP1 = nBuoyLayers+1

         
         ! jas diabatic terms
         !diabaticHeating(nVertLevels,nCells)! "vertical velocity" in buoyancy space
         !wCellCenter       = 0.0

         !jas issue
         ! Get diabaticTimeTendency of a buoyancy surface, omega with funny hat, if any.
         !call any existing MPAS-O subroutines for this


         !-------------------------------------------------------------
         ! begin computation
         !-------------------------------------------------------------

         call get_masks_in_buoyancy_coordinates(nVertLevels, nCells, nBuoyLayers, &
           maxLevelCell, potentialDensity, potentialDensityMidRef, &
           firstLayerBuoyCoor, lastLayerBuoyCoor, buoyMask)

         if(config_oac_epft_debug) then
            print *, ' '
            print *, 'timeLevel:', timeLevel
            print *, ' '
            print *, 'potentialDensityTopRef'
            print *, potentialDensityTopRef
            print *, 'potentialDensityMidRef'
            print *, potentialDensityMidRef
            print *, 'nCells*nBuoyLayers', nCells*nBuoyLayers
            print *, 'sum(buoyMask)', sum(buoyMask)
            print *, 'nCells*nVertLevels', nCells*nVertLevels
            print *, 'sum(mesh%cellMask%array)', sum(mesh%cellMask%array)
            print *, 'minval(potentialDensity), maxval(potentialDensity)'
            print *, minval(potentialDensity), maxval(potentialDensity)
            print *, 'minval(density), maxval(density)'
            print *, minval(density), maxval(density)
         endif


! INTERPOLATION TEST 1
! stratified, horizontally uniform
! Interpolating from z, rho to z, rho
         if(config_oac_epft_debug) then
           do i = 1, nCells
              array1_3D(:,i) = -zMid(:,nCells/2)
              array2_3D(:,i) = potentialDensity(:,nCells/2)
           end do
           print *, ' '
           print *, 'Testing interpolatoin function'
           print *, 'Interpolating from z, rho to z, rho'
           print *, 'call linear_interp_1d_field_along_column(nVertLevels, nCells, &
              nVertLevels, maxLevelCell, array1_3D, array2_3D, array1_3D(:,1), array3_3D)'
           
           print *, 'sum(array1_3D)/nCells + sum(zMid(:,nCells/2))'
           print *, sum(array1_3D)/nCells + sum(zMid(:,nCells/2))
           print *, 'sum(array1_3D)/nCells - sum(array1_3D(:,1))'
           print *, sum(array1_3D)/nCells - sum(array1_3D(:,1))
  
           call linear_interp_1d_field_along_column(nVertLevels, nCells, nVertLevels, &
              maxLevelCell, array1_3D, array2_3D, array1_3D(:,1), array3_3D)
           print *, 'array1_3D(:,1)'
           print *, array1_3D(:,1)
           print *, '-zMid(:,nCells/2)'
           print *, -zMid(:,nCells/2)
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


! INTERPOLATION TEST 2
! Define a stratification where potential density varies linearly with depth
! Using reference potential density that varies linearly with index
! Interpolate z from that potential density to reference potential density
! compare to expected values
         if(config_oac_epft_debug) then
            do i = 1,nCells
               do k = 1, nVertLevels
                  array1_3D(k,i) = config_rhomin_buoycoor*1.02 + &
                     (zMid(k,i)-zMid(1,i)) * &
                     (config_rhomax_buoycoor*0.98 - config_rhomin_buoycoor*1.02) / &
                     (zMid(nVertLevels,i) - zMid(1,i))
                  array2_3D(k,i) = config_rhomin_buoycoor*1.02 + &
                     (zTop(k,i)-zMid(1,i)) * &
                     (config_rhomax_buoycoor*0.98 - config_rhomin_buoycoor*1.02) / &
                     (zMid(nVertLevels,i) - zMid(1,i))
               end do
            end do
            call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyLayers, &
               maxLevelCell, array1_3D, zMid, potentialDensityMidRef, array1_3Dbuoy)
            
            do i = 1,nCells
               do k = 1, nBuoyLayers
                  array2_3Dbuoy(k,i) = zMid(1,i) + &
                     (potentialDensityMidRef(k) - config_rhomin_buoycoor*1.02) * &
                     (zMid(nVertLevels,i) - zMid(1,i)) / &
                     (config_rhomax_buoycoor*0.98 - config_rhomin_buoycoor*1.02)
               end do
            end do
            do i = 1,nCells
               do k = 1, nBuoyLayers
                  RMSlocal2 = RMSlocal2 + &
                     ((array1_3Dbuoy(k,i) - array2_3Dbuoy(k,i))/array2_3Dbuoy(k,i))**2
               end do
            end do
         endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! end chunk for testing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! start chunk commented during testing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call check_potentialDensityRef_range(nVertLevels, nCells, maxLevelCell, potentialDensity)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyLayers, &
            maxLevelCell, -potentialDensity, zMid, &
            -potentialDensityMidRef, heightMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyLayers, &
            maxLevelCell, -potentialDensity, zMid, &
            -potentialDensityTopRef, heightTopBuoyCoor)
         do i=1,nCells
            !correct the top of heightTopBuoyCoor
            do k=1,firstLayerBuoyCoor(i)
              heightTopBuoyCoor(k,i)=zTop(1,i)
            enddo
            ! correct the bottom of heightTopBuoyCoor
            do k=lastLayerBuoyCoor(i)+1,nBuoyLayers
              heightTopBuoyCoor(k,i)=-bottomDepth(i)
            enddo
            ! copy into interface variable
            heightInterfaceBuoyCoor(1:nBuoyLayers,i)=heightTopBuoyCoor(1:nBuoyLayers,i)
            heightInterfaceBuoyCoor(nBuoyLayers+1,i)=-bottomDepth(i)
         enddo

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyLayers, &
            maxLevelCell, -potentialDensity, uCellCenter, &
            -potentialDensityMidRef, uMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyLayers, &
            maxLevelCell, -potentialDensity, vCellCenter, &
            -potentialDensityMidRef, vMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyLayers, &
            maxLevelCell, -potentialDensity, density, &
            -potentialDensityMidRef, densityMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyLayers, &
            maxLevelCell, -potentialDensity, density, &
            -potentialDensityTopRef, densityTopBuoyCoor)

         !call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyLayers, &
         !   maxLevelCell, -potentialDensity, Q, potentialDensityTopRef, QMidRef)


         call computeBuoyancyColumn(nBuoyLayers, potentialDensityMidRef, buoyancyMidRef)
         call computeBuoyancyColumnP1(nBuoyLayersP1, potentialDensityTopRef, &
            buoyancyInterfaceRef)


         call computeSigma(nCells, nBuoyLayers, firstLayerBuoyCoor, lastLayerBuoyCoor, &
            heightInterfaceBuoyCoor, buoyancyInterfaceRef, sigma)


         call computeMontgomeryPotential(nBuoyLayers, nCells, surfacePressure, &
            firstLayerBuoyCoor, lastLayerBuoyCoor, SSH, densityMidBuoyCoor, &
            potentialDensityMidRef, heightInterfaceBuoyCoor, montgPotBuoyCoor)
         call computeNormalGradientOnEdge(nBuoyLayers, nCells, nEdges, &
            mesh, &
            montgPotBuoyCoor, montgPotNormalGradOnEdge)
         call mpas_reconstruct(mesh, montgPotNormalGradOnEdge, &
            montgPotGradX, montgPotGradY, montgPotGradZ, &
            montgPotGradZonal, montgPotGradMerid)

! jas issue: in some cases it might be cleaner to pass mesh instead of 
! nBuoyLayers, nCells, maxLevelCell...
  
         ! Increment first-order running mean fields:
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            buoyMask, buoyMaskEA)
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            sigma, sigmaEA)
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            heightMidBuoyCoor, heightMidBuoyCoorEA)
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            montgPotBuoyCoor, montgPotBuoyCoorEA)
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            montgPotGradZonal, montgPotGradZonalEA)
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            montgPotGradMerid, montgPotGradMeridEA)
   
   
         ! Increment second-order running mean fields
         wrk3DBuoyCoor = heightMidBuoyCoor * heightMidBuoyCoor
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, heightMidBuoyCoorSqEA)

         wrk3DBuoyCoor = heightMidBuoyCoor * montgPotGradZonal
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, HeightMGradZonalEA)

         wrk3DBuoyCoor = heightMidBuoyCoor * montgPotGradMerid
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, HeightMGradMeridEA)

         wrk3DBuoyCoor = uMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, usigmaEA)

         wrk3DBuoyCoor = vMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, vsigmaEA)

         !wrk3DBuoyCoor = wMidBuoyCoor * sigma
         !call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
         !   wrk3DBuoyCoor, wsigmaEA)
   
   
         ! Increment third-order running mean fields 
         wrk3DBuoyCoor = uMidBuoyCoor * uMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, uusigmaEA)

         wrk3DBuoyCoor = vMidBuoyCoor * vMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, vvsigmaEA)

         wrk3DBuoyCoor = uMidBuoyCoor * vMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, uvsigmaEA)

         !wrk3DBuoyCoor = uMidBuoyCoor * wMidBuoyCoor * sigma
         !call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
         !   wrk3DBuoyCoor, uwsigmaEA)
         uwsigmaEA = 0.0

         !wrk3DBuoyCoor = vMidBuoyCoor * wMidBuoyCoor* sigma
         !call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
         !   wrk3DBuoyCoor, vwsigmaEA)
         vwsigmaEA = 0.0
   
         ! update number of samples in ensemble average
         amEPFT % nSamplesEA % scalar = amEPFT % nSamplesEA % scalar + 1
   
   
   
         ! Calculate the thickness weighted averages
         call calculateTWA(nBuoyLayers, nCells, nBuoyLayers, &
            sigmaEA, usigmaEA,  uTWA)
         call calculateTWA(nBuoyLayers, nCells, nBuoyLayers, &
            sigmaEA, vsigmaEA,  vTWA)
         !call calculateTWA(nBuoyLayers, nCells, nBuoyLayers, &
         !   sigmaEA, wsigmaEA,  wTWA)
         wTWA = 0.0

         
         call calculateEPFTfromTWA(nBuoyLayers, nCells, &
            sigmaEA, heightMidBuoyCoorEA, &
            heightMidBuoyCoorSqEA, montgPotGradZonalEA, montgPotGradMeridEA, &
            HeightMGradZonalEA, HeightMGradMeridEA, uTWA, vTWA, wTWA, &
            uusigmaEA, vvsigmaEA, uvsigmaEA, uwsigmaEA, vwsigmaEA, EPFT)

        call calculateDivEPFT(nBuoyLayers, nCells, nEdges, &
            mesh, buoyancyInterfaceRef, sigmaEA, buoyMaskEA, EPFT, divEPFT)
         
        call calculateErtelPVFlux(nCells, nBuoyLayers, &
            sigmaEA, divEPFT, ErtelPVFlux)

        call calculateErtelPVTendencyFromPVFlux(nBuoyLayers, nCells, nEdges, &
            mesh, sigmaEA, ErtelPVFlux, ErtelPVTendency)


         fCell => mesh % fCell % array
         call computeErtelPV(nCells, nBuoyLayers, nEdges, mesh, &
            fCell, uTWA, vTWA, sigmaEA, ErtelPV)

         ! Compute the geometric decomposition in terms of angles and 
         ! eccentricities using the entries of EPFT.
         !call eddyGeomDecompEPFT(EPFT, ...)




         ! calculate potential vorticity fluxes using curl of u
         if(config_oac_epft_debug) then

            relativeVorticityCell => diagnostics % relativeVorticityCell % array
            
            ! store relVortMidBuoyCoor in array1_3Dbuoy
            call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyLayers, &
               maxLevelCell, -potentialDensity, relativeVorticityCell, &
               -potentialDensityMidRef, array1_3Dbuoy)
            
            do i = 1,nCells
               do k=firstLayerBuoyCoor(i), lastLayerBuoyCoor(i)
                  PVMidBuoyCoor(k,i) = (fCell(i) + array1_3Dbuoy(k,i) ) / sigma(k,i)
               end do
            end do

            call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
               uMidBuoyCoor, uMidBuoyCoorEA)
            
            call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
               vMidBuoyCoor, vMidBuoyCoorEA)
            
            call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
               PVMidBuoyCoor, PVMidBuoyCoorEA)
            
            wrk3DBuoyCoor = uMidBuoyCoor * PVMidBuoyCoor
            call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
               wrk3DBuoyCoor, uPVMidBuoyCoorEA)
            
            wrk3DBuoyCoor = vMidBuoyCoor * PVMidBuoyCoor
            call updateEnsembleAverage(nBuoyLayers, nCells, nSamplesEA, &
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



         ! Clean up
         ! jas issue: make sure I deallocate everything
         call mpas_deallocate_scratch_field(amEPFT % firstLayerBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % lastLayerBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % buoyMask, .true.)
         call mpas_deallocate_scratch_field(amEPFT % heightMidBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % heightTopBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % heightInterfaceBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % uMidBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % vMidBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % densityMidBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % densityTopBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % buoyancyMidRef, .true.)
         call mpas_deallocate_scratch_field(amEPFT % buoyancyInterfaceRef, .true.)
         call mpas_deallocate_scratch_field(amEPFT % sigma, .true.)
         call mpas_deallocate_scratch_field(amEPFT % montgPotBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % montgPotNormalGradOnEdge, .true.)
         call mpas_deallocate_scratch_field(amEPFT % wrk3DnVertLevels, .true.)
         call mpas_deallocate_scratch_field(amEPFT % wrk3DBuoyCoor, .true.)

         call mpas_deallocate_scratch_field(amEPFT % array1_3D, .true.)
         call mpas_deallocate_scratch_field(amEPFT % array2_3D, .true.)
         call mpas_deallocate_scratch_field(amEPFT % array3_3D, .true.)
         call mpas_deallocate_scratch_field(amEPFT % array1_3Dbuoy, .true.)
         call mpas_deallocate_scratch_field(amEPFT % array2_3Dbuoy, .true.)

         call mpas_deallocate_scratch_field(amEPFT % PVMidBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(amEPFT % PVMidBuoyCoorEA, .true.)
         call mpas_deallocate_scratch_field(amEPFT % uPVMidBuoyCoorEA , .true.)
         call mpas_deallocate_scratch_field(amEPFT % vPVMidBuoyCoorEA, .true.)
         call mpas_deallocate_scratch_field(amEPFT % PVFluxTest, .true.)


         nCellsCum = nCellsCum + nCells

         block => block % next
      end do

      ! mpi gather/scatter calls may be placed here.
      if(config_oac_epft_debug) then
        RMSglobal1 = 1.0D36
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

      
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'am_eliassen_palm_flux_tensor', am_epftPool)

         ! assignment of final am_eliassen_palm_flux_tensor variables could occur here.

         block => block % next
      end do

      call mpas_timer_stop("eliassen_palm_flux_tensor", am_eliassen_palm_flux_tensorTimer)

   end subroutine ocn_compute_eliassen_palm_flux_tensor!}}}

!***********************************************************************
!
!  routine ocn_restart_eliassen_palm_flux_tensor
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_eliassen_palm_flux_tensor(domain, err)!{{{

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

   end subroutine ocn_restart_eliassen_palm_flux_tensor!}}}

!***********************************************************************
!
!  routine ocn_finalize_eliassen_palm_flux_tensor
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_eliassen_palm_flux_tensor(domain, err)!{{{

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

   end subroutine ocn_finalize_eliassen_palm_flux_tensor!}}}

end module ocn_eliassen_palm_flux_tensor

! vim: foldmethod=marker
