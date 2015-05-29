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


!***********************************************************************
!
!  subroutine get_masks_in_buoyancy_coordinates
!
!> \brief   Get masks in buoyancy coordinates
!> \author  Juan A. Saenz
!> \date    Jan 2014
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
   subroutine get_masks_in_buoyancy_coordinates(nVertLevels, nCells, nBuoyLayers, &
      maxLevelCell, potentialDensity, potentialDensityMidRef, &
      firstLayerBuoyCoor, lastLayerBuoyCoor, buoyMask)!{{{
      
      integer, intent(in) :: nVertLevels, nCells, nBuoyLayers
      integer, dimension(nCells), intent(in)  :: maxLevelCell
      integer, dimension(nCells), intent(out) :: firstLayerBuoyCoor
      integer, dimension(nCells), intent(out) :: lastLayerBuoyCoor
      real (kind=RKIND), dimension(nBuoyLayers, nCells), intent(out) :: buoyMask
      real (kind=RKIND), dimension(nVertLevels, nCells), intent(in) :: potentialDensity
      real (kind=RKIND), dimension(nBuoyLayers), intent(in) :: potentialDensityMidRef

      ! Local variables
      integer :: iCell, maxLevel, kB, kBBottom, kBTop
      
      firstLayerBuoyCoor = 1
      lastLayerBuoyCoor = nBuoyLayers
      buoyMask     = 0.0


      do iCell = 1, nCells
      
         maxLevel = maxLevelCell(iCell)

         do kB = 1, nBuoyLayers
            if (potentialDensityMidRef(kB) >= potentialDensity(1,iCell) ) then
               firstLayerBuoyCoor(iCell) = kB
               exit
            endif
         enddo
         
         do kB = nBuoyLayers, 1, -1
            if (potentialDensityMidRef(kB) <= potentialDensity(maxLevel,iCell) ) then
               lastLayerBuoyCoor(iCell) = kB
               exit
            endif
         enddo
         
         ! set mask to 1 inside the range
         do kB = firstLayerBuoyCoor(iCell), lastLayerBuoyCoor(iCell)
           buoyMask(kB,iCell) = 1.0
         enddo
      
      enddo
   
   end subroutine get_masks_in_buoyancy_coordinates!}}}



!***********************************************************************
!
!  subroutine check_potentialDensityRef_range
!
!> \brief   Check if the range of values in potentialDensityTopRef contains current state
!> \author  Juan A. Saenz
!> \date    Jan 2014
!> \details 
!>  Check if the range of values in potentialDesnityTopRef contains all values in 
!>  potentialDensity of the current state. 
!>  If not, print a warning.
!
!-----------------------------------------------------------------------
   subroutine check_potentialDensityRef_range(nVertLevels, nCells, maxLevelCell, &
         potentialDensity)!{{{
      integer, intent(in) :: nVertLevels, nCells
      integer, dimension(nCells), intent(in)  :: maxLevelCell
      real (kind=RKIND), dimension(nVertLevels, nCells), intent(in) :: potentialDensity

      ! Local variables
      integer :: k, i
      logical :: printWarning
      
      printWarning = .false.
      
      do i = 1, nCells
         if (potentialDensity(1,i) < config_rhomin_buoycoor) then
            printWarning = .true.
            exit
         end if
         if (potentialDensity(maxLevelCell(i),i) > config_rhomax_buoycoor) then
            printWarning = .true.
            exit
         end if
      enddo
      
      !jas issue: do we want to print a warning once, or at every i,k out of range?
      if (printWarning) then
         write(stderrUnit,*) 'Warning: in EPFT package, reference potential density does &
            not span the values of potentialDensity in the current state'
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
!>  Required: xFieldIn monotonically decreases with index value
!
!-----------------------------------------------------------------------

   subroutine linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyLayers, &
      maxLevelCell, xFieldIn, yFieldIn, xColumnOut, yFieldOut)!{{{
      
      integer, intent(in) :: nVertLevels, nCells, nBuoyLayers
      integer, dimension(nCells), intent(in)  :: maxLevelCell
      real (kind=RKIND), dimension(nVertLevels, nCells), intent(in) :: xFieldIn
      real (kind=RKIND), dimension(nVertLevels, nCells), intent(in) :: yFieldIn
      real (kind=RKIND), dimension(nBuoyLayers), intent(in) :: xColumnOut
      real (kind=RKIND), dimension(nBuoyLayers, nCells), intent(out) :: yFieldOut

      ! Local variables
      integer :: iCell, maxLevel, kB, kBBottom, kBTop, kDataAbove, kDataBelow, kData
      real (kind=RKIND) :: dx, dy

      yFieldOut = 0.0

      do iCell = 1, nCells
      
         ! find the index of the bottom level of a column
         maxLevel = maxLevelCell(iCell)

         ! Monotonically decreasing xFieldIn required
         ! Find index of first element in xColumnOut that is inside xFieldIn(:,iCell)
         kBTop = 1
         do kB = 1, nBuoyLayers
            ! the following line ensures that
            ! if all xColumnOut > xFieldIn(1,iCell) then kBTop = nBuoyLayers
            kBTop = kB
            if (xColumnOut(kB) <= xFieldIn(1,iCell) ) then
               exit
            endif
         enddo
         
         !find last target buoyancy level inside column
         kBBottom = nBuoyLayers
         do kB = nBuoyLayers, 1, -1
            ! the following line ensures that
            ! if all xColumnOut < xFieldIn(1,iCell) then kBBottom = 1
            kBBottom = kB
            if (xColumnOut(kB) >= xFieldIn(maxLevel,iCell) ) then
               exit
            endif
         enddo
         
         
         ! For the target x levels outside the x range in a column:
         ! set data from 1:kBTop-1 to surface values
         do kB = 1, kBTop-1
           yFieldOut(kB,iCell) = yFieldIn(1,iCell)
         enddo
         !set data from kBBottom+1:nBuoyLayers to bottom values
         do kB = kBBottom+1, nBuoyLayers
           yFieldOut(kB,iCell) = yFieldIn(maxLevel,iCell)
         enddo

         ! The interpolation: 
         ! for the target buoyancy levels within the buoyancy range in a column:
         ! jas issue: this can be replaced by a call to 
         !     src/operators/mpas_spline_interpolatoin.F:mpas_interpolate_linear()
         kDataAbove = 1
         kDataBelow = kDataAbove + 1
         do kB = kBTop, kBBottom
            ! for each xColumnOut(kB) value, find the corresponding upper and lower 
            ! xFieldIn value in the field data, then interpolate y between those values.
            if (xColumnOut(kB) < xFieldIn(kDataBelow,iCell)) then
               do kData = kDataBelow, maxLevel
                  if (xColumnOut(kB) > xFieldIn(kData,iCell) ) then
                     kDataBelow=kData
                     kDataAbove=kDataBelow-1
                     exit
                  endif
               enddo
            endif
            
            dx = xFieldIn(kDataBelow,iCell) - xFieldIn(kDataAbove,iCell)
            dy = yFieldIn(kDataBelow,iCell) - yFieldIn(kDataAbove,iCell)
            yFieldOut(kB,iCell) = yFieldIn(kDataAbove,iCell) + &
               (xColumnOut(kB)-xFieldIn(kDataAbove,iCell)) * dy/dx
         enddo
      
      enddo
   
   end subroutine linear_interp_1d_field_along_column!}}}


!***********************************************************************
!
!  subroutine computeBuoyancyColumn
!
!> \brief   Compute buoyancy
!> \author  Juan A. Saenz
!> \date    17 December 2013
!> \details 
!>  This subroutine computes buoyancy
!
!-----------------------------------------------------------------------

   subroutine computeBuoyancyColumn(nLayers, potentialDensity, buoyancy)!{{{
      integer, intent(in) :: nLayers
      real (kind=RKIND), dimension(nLayers), intent(in) :: potentialDensity
      real (kind=RKIND), dimension(nLayers), intent(out) :: buoyancy
   
      !local variables
      integer :: i, k
      real (kind=RKIND) :: rho0
      
      rho0 = config_density0
      
      buoyancy = 0.0
      
      do k = 1, nLayers
         buoyancy(k) = -gravity * (potentialDensity(k)-rho0) / rho0
      enddo
   
   end subroutine computeBuoyancyColumn!}}}


!***********************************************************************
!
!  subroutine computeBuoyancyColumnP1
!
!> \brief   Compute buoyancy
!> \author  Juan A. Saenz
!> \date    Jan 2014
!> \details 
!>  This subroutine computes buoyancy
!
!-----------------------------------------------------------------------

   subroutine computeBuoyancyColumnP1(nLayers, potentialDensity, buoyancy)!{{{
      integer, intent(in) :: nLayers
      real (kind=RKIND), dimension(nLayers-1), intent(in) :: potentialDensity
      real (kind=RKIND), dimension(nLayers), intent(out) :: buoyancy
   
      !local variables
      integer :: i, k
      real (kind=RKIND) :: rho0
      
      rho0 = config_density0
      
      buoyancy = 0.0
      
      do k = 1, nLayers-1
         buoyancy(k) = -gravity * (potentialDensity(k)-rho0) / rho0
      enddo
   
      buoyancy(nLayers) = -gravity * (config_rhomax_buoycoor-rho0) / rho0
      
   
   end subroutine computeBuoyancyColumnP1!}}}



!***********************************************************************
!
!  subroutine computeSigma
!
!> \brief   Calculate the inverse of the derivative of buoy wrt z
!> \author  Juan A. Saenz
!> \date    December 2013
!> \details 
!>  This subroutine calculates the inverse of the derivative of buoy wrt z.
!
!-----------------------------------------------------------------------

   subroutine computeSigma(nCells, nLayers, firstLayerBuoyCoor, &
         lastLayerBuoyCoor, heightInterface, buoyInterface, sigma)!{{{
      integer, intent(in) :: nCells, nLayers
      integer, dimension(:), intent(in) :: firstLayerBuoyCoor
      integer, dimension(:), intent(in) :: lastLayerBuoyCoor
      real (kind=RKIND), dimension(:,:), intent(in) :: heightInterface
      real (kind=RKIND), dimension(:), intent(in) :: buoyInterface
      real (kind=RKIND), dimension(:,:), intent(out) :: sigma
      
      
      ! local variables
      integer :: i, k
      
      sigma = 0.0
      
      do i = 1, nCells
         do k = 1,nLayers
            sigma(k,i) = (heightInterface(k,i) - heightInterface(k+1,i)) / &
               (buoyInterface(k) - buoyInterface(k+1))
         enddo
      enddo
   
   end subroutine computeSigma!}}}



!***********************************************************************
!
!  subroutine computeMontgomeryPotential
!
!> \brief   Compute the Montgomery potential
!> \author  Juan A. Saenz
!> \date    17 December 2013
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

   subroutine computeMontgomeryPotential(nLayers, nCells, pSurface, firstLayer, &
         lastLayer, SSH, density, potDens, heightInterface, MontgomeryPotential)!{{{
      integer, intent(in) :: nLayers, nCells
      integer, dimension(nCells), intent(in) :: firstLayer
      integer, dimension(nCells), intent(in) :: lastLayer
      real (kind=RKIND), dimension(nCells), intent(in) :: pSurface
      real (kind=RKIND), dimension(nCells), intent(in) :: SSH
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: density
      real (kind=RKIND), dimension(nLayers), intent(in) :: potDens
      real (kind=RKIND), dimension(nLayers+1, nCells), intent(in) :: heightInterface
      real (kind=RKIND), dimension(nLayers, nCells), intent(out) :: MontgomeryPotential
      
      ! local variables
      integer :: i, k
      real (kind=RKIND) :: pInterfacek ! pressure at interface k, i.e. at top of layer k
      
      MontgomeryPotential = 0.0
      

      do i = 1, nCells
         
         !pInterfacek = pSurface(i)
         pInterfacek = 0.0
         k = 1
         MontgomeryPotential(k,i) = pInterfacek/potDens(k) + gravity*heightInterface(k,i)
         
         do k = 2, nLayers 
            pInterfacek = pInterfacek + &
               gravity * ( heightInterface(k-1,i)-heightInterface(k,i) ) * density(k-1,i)
            MontgomeryPotential(k,i) = MontgomeryPotential(k-1,i) + &
               pInterfacek * ( 1/potDens(k) - 1/potDens(k-1) )
         enddo
      
      enddo
   
   end subroutine computeMontgomeryPotential!}}}



!***********************************************************************
!
!  subroutine computeNormalGradientOnEdge
!
!> \brief   Compute the gradient of a quantity that exists on cell centers
!> \author  Juan A. Saenz
!> \date    December 2013
!> \details 
!>  This subroutine computes the gradient of a quantity that exists on cell centers
!
!-----------------------------------------------------------------------

   subroutine computeNormalGradientOnEdge(nBLayers, nCells, nEdges, &
         mesh, field, normalGradOnEdge)!{{{
      integer, intent(in) :: nBLayers, nCells, nEdges
      type (mesh_type), intent(in) :: mesh !< Input: mesh information
      real (kind=RKIND), dimension(:,:), intent(in) :: field
      real (kind=RKIND), dimension(:,:), intent(out) :: normalGradOnEdge

      !local variables
      integer :: nEdgesSolve, iEdge, k, cell1, cell2, kMin, kMax
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: boundaryEdge
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND) :: invLength

      cellsOnEdge       => mesh % cellsOnEdge % array
      dcEdge            => mesh % dcEdge % array
      maxLevelEdgeTop   => mesh % maxLevelEdgeTop % array
      boundaryEdge      => mesh % boundaryEdge % array

      normalGradOnEdge = 0.0

      do iEdge = 1, nEdges
         ! enforce enforce zero gradient on boundary edges
         if (boundaryEdge(1,iEdge) == 1) then
            normalGradOnEdge(:,iEdge) = 0.0
         else
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            invLength = 1.0 / dcEdge(iEdge)
            do k = 1, mesh % nBuoyancyLayers
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
!> \author  Juan A. Saenz
!> \date    17 December 2013
!> \details 
!>  This subroutine updates the ensemble average
!
!-----------------------------------------------------------------------

   subroutine updateEnsembleAverage(nLayers, nCells, nSamples, A, Abar)!{{{
      integer, intent(in) :: nLayers, nCells, nSamples
      real (kind=RKIND), dimension(nLayers, nCells), intent(in)  :: A
      real (kind=RKIND), dimension(nLayers, nCells), intent(inout)  :: Abar
   
      !test
      integer :: i, k
   
      do i = 1, nCells
         do k = 1, nLayers
            Abar(k,i) = (nSamples * Abar(k,i) + A(k,i)) / (nSamples + 1.0)
         enddo
      enddo

   end subroutine updateEnsembleAverage!}}}


!***********************************************************************
!
!  subroutine calculateTWA
!
!> \brief   Calculate the thickness weighted average
!> \author  Juan A. Saenz
!> \date    January 2014
!> \details 
!>  This subroutine calculates the thickness weighted average
!
!-----------------------------------------------------------------------
   subroutine calculateTWA(nLayers, nCells, nBuoyancyLayers, sigmaEA, &
         varSigmaEA, varTWA)!{{{
      integer, intent(in) :: nLayers, nCells, nBuoyancyLayers
      real (kind=RKIND), dimension(nLayers, nCells), intent(in)  :: sigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in)  :: varSigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(out) :: varTWA
   
      ! local variables
      integer :: i, k
   
      varTWA = 0.0
      
      do i = 1, nCells
         do k = 1,nBuoyancyLayers
            varTWA(k,i) = varSigmaEA(k,i) / max(1.0e-15,sigmaEA(k,i))
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
         sigmaEA, heightEA, heightSqEA, MxEA, MyEA, HMxEA, HMyEA, uTWA, vTWA, wTWA, &
         uuSigmaEA, vvSigmaEA, uvSigmaEA, uwSigmaEA, vwSigmaEA, Etensor)!{{{
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
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: wTWA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: uuSigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: vvSigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: uvSigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: uwSigmaEA
      real (kind=RKIND), dimension(nLayers, nCells), intent(in) :: vwSigmaEA
      real (kind=RKIND), dimension(3, 3, nLayers, nCells), intent(out) :: Etensor
   
      ! local variables
      integer :: iCell, kLayer
      real (kind=RKIND) :: sigma
      real (kind=RKIND) :: uppupp, vppvpp, uppvpp, uppwpp, vppwpp
      real (kind=RKIND) :: HpHp, HpMxp, HpMyp
      
      Etensor = 0.0

      do iCell = 1, nCells
         do kLayer = 1,nLayers
            
            sigma = max(sigmaEA(kLayer,iCell), 1.0e-15)

            uppupp = uuSigmaEA(kLayer,iCell) / sigma - uTWA(kLayer,iCell)*uTWA(kLayer,iCell)
            vppvpp = vvSigmaEA(kLayer,iCell) / sigma - vTWA(kLayer,iCell)*vTWA(kLayer,iCell)
            uppvpp = uvSigmaEA(kLayer,iCell) / sigma - uTWA(kLayer,iCell)*vTWA(kLayer,iCell)
            uppwpp = uwSigmaEA(kLayer,iCell) / sigma - uTWA(kLayer,iCell)*wTWA(kLayer,iCell)
            vppwpp = vwSigmaEA(kLayer,iCell) / sigma - vTWA(kLayer,iCell)*wTWA(kLayer,iCell)
            HpHp   = heightSqEA(kLayer,iCell) - heightEA(kLayer,iCell)*heightEA(kLayer,iCell)
            HpMxp  = HMxEA(kLayer,iCell) - heightEA(kLayer,iCell)*MxEA(kLayer,iCell)
            HpMyp  = HMyEA(kLayer,iCell) - heightEA(kLayer,iCell)*MyEA(kLayer,iCell)

            !EPTF_pq(x,y,z) is represented as EPFT(p,q,kLayer,iCell)
            !column 1: Eu
            Etensor(1,1,kLayer,iCell) = uppupp + 0.5 * HpHp / sigma
            Etensor(2,1,kLayer,iCell) = uppvpp
            Etensor(3,1,kLayer,iCell) = uppwpp + HpMxp / sigma

            !column 2: Ev
            Etensor(1,2,kLayer,iCell) = uppvpp 
            Etensor(2,2,kLayer,iCell) = vppvpp + 0.5 * HpHp / sigma
            Etensor(3,2,kLayer,iCell) = vppwpp + HpMyp / sigma

            !column 3: Ew
            Etensor(1,3,kLayer,iCell) = 0.0
            Etensor(2,3,kLayer,iCell) = 0.0
            Etensor(3,3,kLayer,iCell) = 0.0
            
         enddo
      enddo
   
   end subroutine calculateEPFTfromTWA!}}}


!***********************************************************************
!
!  subroutine calculateDivEPFT
!
!> \brief   Calculate the divergence of EPFT 
!> \author  Juan A. Saenz
!> \date    January 2014
!> \details 
!>  This subroutine calculates the divergence of the Elliassen-Palm flux tensor
!
!-----------------------------------------------------------------------

   subroutine calculateDivEPFT(nLayers, nCells, nEdges, &
        mesh, buoyancyInterfaceRef, sigmaEA, buoyMaskEA, tensorCellIn, vectorCellOut)!{{{

      use mpas_vector_operations

      integer, intent(in) :: nLayers, nCells, nEdges
      type (mesh_type), intent(in) :: mesh
      real (kind=RKIND), dimension(:), intent(in) :: buoyancyInterfaceRef
      real (kind=RKIND), dimension(:,:), intent(in) :: sigmaEA
      real (kind=RKIND), dimension(:,:), intent(in) :: buoyMaskEA
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: tensorCellIn
      real (kind=RKIND), dimension(:,:,:), intent(out) :: vectorCellOut

      ! local variables
      logical :: includeHalo, on_a_sphere
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
      real (kind=RKIND) :: rho0
      
      ! variables used for testing and debugging
      real (kind=RKIND), dimension(:), allocatable :: divExact
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      xCell => mesh % xCell % array
      yCell => mesh % yCell % array
      zCell => mesh % zCell % array

      if (config_oac_epft_debug) then
        allocate(divExact(nCells+1))
      end if


      rho0 = config_density0

      allocate(scalarWrk1(nLayers,nCells+1))
      allocate(vectorCellWrk1(3,nLayers,nCells+1))
      allocate(vectorCellWrk2(3,nLayers,nCells+1))
      allocate(vectorEdgeWrk1(3,nLayers,nEdges))
      allocate(vertVector(nLayers))

      on_a_sphere = mesh % on_a_sphere
      edgeSignOnCell => mesh % edgeSignOnCell % array
      latCell => mesh % latCell % array
      lonCell => mesh % lonCell % array
      
      includeHalo = .true.
      

      ! initialize work and intent(out)
      vectorCellOut = 0.0

      ! loop over all three column vectors
      do q = 1, 3
         
         scalarWrk1 = 0.0

         ! horizontal derivatives
         vectorCellWrk1 = tensorCellIn(:,q,:,:)

         ! weight the vector with sigmaEA(:,:)
         do iComponent = 1,3
            vectorCellWrk1(iComponent,:,:) = sigmaEA(:,:)*vectorCellWrk1(iComponent,:,:)
         enddo


         ! use q=3 as a test vector
         if (q.eq.3 .and. config_oac_epft_debug) then
           do iCell = 1,nCells
             vectorCellWrk1(1,:,iCell) = xCell(iCell)
             vectorCellWrk1(2,:,iCell) = yCell(iCell)
             vectorCellWrk1(3,:,iCell) = zCell(iCell)
             ! the analytical divergence:
             divExact(iCell) = -1.0*sin(lonCell(iCell)) * &
                (1.0 + 2.0*sin(latCell(iCell))) !+ 3.0*sin(latCell(iCell))
           enddo
         endif

         ! zero the vertical component of vectorCellWrk1
         ! the vertical will be treated seperately below
         vectorCellWrk1(3,:,:) = 0.0

         if (on_a_sphere) then
            
            do iCell = 1,nCells
               do kLayer = 1,nLayers
                  call mpas_vector_LonLatR_to_R3(vectorCellWrk1(:,kLayer,iCell), &
                     lonCell(iCell), latCell(iCell), vectorCellWrk2(:,kLayer,iCell))
               end do
            end do

            vectorCellWrk1 = vectorCellWrk2

         end if

         call mpas_vector_R3Cell_to_Edge(vectorCellWrk1, mesh, &
            vectorEdgeWrk1)

         call mpas_divergence_in_r3_buoyancy(vectorEdgeWrk1, mesh, &
            edgeSignOnCell, includeHalo, scalarWrk1)

         ! use q=3 as a test vector
         if (q.eq.3 .and. config_oac_epft_debug) then
           print *, ' '
           do kLayer = 1,nLayers
              wrk = sqrt( &
                sum( &
                ( & 
                (divExact(:)-scalarWrk1(kLayer,:))/ max(abs(divExact(:)),1.0e-15)  &
                )**2 * &
                (1.0 - mesh % boundaryCell % array(1,:)) &
                ) / nCells )
              print *, 'div RMS relative error on layer:', wrk
           enddo
         endif


         if (q < 3 .or. .not. config_oac_epft_debug) then
            do iCell = 1,nCells
               do kLayer = 1,nLayers
                  sigma = max(sigmaEA(kLayer,iCell), 1.0e-15)
                  scalarWrk1(kLayer,iCell) = scalarWrk1(kLayer,iCell) / sigma
               end do
            end do
         end if


         ! vertical derivative
         do iCell = 1,nCells

            ! copy the vertical component of EPFT into a work array
            vertVector(:) = tensorCellIn(3,q,:,iCell)

            ! use q=3 as a test vector
            if (q.eq.3 .and. config_oac_epft_debug) then
               vertVector(:) = 0.0
            endif


            do kLayer = 1,nLayers

               wrk = 0.0

               ! jas issue: change buoyancyInterfaceRef to buoyancyMidRef and generalize 
               if(kLayer.eq.1) then
                 wrkAbove=sigmaEA(kLayer,iCell)*vertVector(kLayer)
                 wrkBelow=sigmaEA(kLayer+1,iCell)*vertVector(kLayer+1)
                 db = buoyancyInterfaceRef(kLayer)-buoyancyInterfaceRef(kLayer+1) 
               else if (kLayer.eq.nLayers) then
                 wrkAbove=sigmaEA(kLayer-1,iCell)*vertVector(kLayer-1)
                 wrkBelow=sigmaEA(kLayer,iCell)*vertVector(kLayer)
                 db = buoyancyInterfaceRef(kLayer)-buoyancyInterfaceRef(kLayer+1) 
               else
                 wrkAbove=sigmaEA(kLayer-1,iCell)*vertVector(kLayer-1)
                 wrkBelow=sigmaEA(kLayer+1,iCell)*vertVector(kLayer+1)
                 db = 2.0*(buoyancyInterfaceRef(kLayer)-buoyancyInterfaceRef(kLayer+1))
               endif

               ! jas issue: should sigma also be clipped like this for k-1 and k+1?
               sigma = max(sigmaEA(kLayer,iCell), 1.0e-15)
               wrk = (wrkAbove - wrkBelow) / db / sigma

               scalarWrk1(kLayer,iCell) = scalarWrk1(kLayer,iCell) + wrk

               ! temporarily mask divEPFT to add in visualization
               !if(buoyMaskEA(kLayer,iCell).lt.0.5) scalarWrk1(kLayer,iCell) = 0.0

            end do

         end do

         vectorCellOut(q,:,:) = scalarWrk1

      end do


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
      ErtelPVFlux(2,:,:) = -1.0 * divEPFT(1,:,:)
      ErtelPVFlux(3,:,:) = 0.0

      do i = 1, nCells
         do k = 1,nBuoyancyLayers
            ErtelPVFlux(:,k,i) = ErtelPVFlux(:,k,i) / max(sigma(k,i),1.0e-15)
         end do
      end do

    end subroutine calculateErtelPVFlux



!***********************************************************************
!
!  routine mpas_tensor_cell_to_edge_BuoyCoor
!
!> \brief   Interpolate a matrix from cell to edge
!> \author  Mark Petersen, Juan A. Saenz
!> \date    Jan 2014
!> \details 
!>  This routine interpolates a matrix from cell to edge locations,
!>  looping through nBuoyancyLayers.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_cell_to_edge_BuoyCoor(matrixCell, grid, &
        includeHalo, matrixEdge)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: &
         matrixCell   !< Input: matrix located at Cell

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      logical, intent(in) :: & 
         includeHalo !< Input: If true, halo cells and edges are included in computation

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:,:), intent(out) :: &
         matrixEdge   !< Output: matrix located at Edge

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, cell1, cell2, p, q, k
      integer :: nEdgesCompute, nBuoyLayers, nCells
      integer, dimension(:,:), pointer :: cellsOnEdge

      if (includeHalo) then
         nEdgesCompute = grid % nEdges
      else 
         nEdgesCompute = grid % nEdgesSolve
      endif
      nBuoyLayers = grid % nBuoyancyLayers
      nCells = grid % nCells

      cellsOnEdge => grid % cellsOnEdge % array

      ! error check that index 1 of matrixEdge and matrixCell are same length?

      do iEdge=1,nEdgesCompute
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nBuoyLayers
            do q = 1, 3
               do p = 1, 3
                  matrixEdge(p,q,k,iEdge) = &
                     0.5*(matrixCell(p,q,k,cell1) + matrixCell(p,q,k,cell2))
               end do
            end do
         enddo
      enddo

   end subroutine mpas_tensor_cell_to_edge_BuoyCoor!}}}


!***********************************************************************
!
!  subroutine calculateErtelPVTendencyFromPVFlux
!
!> \brief   Calculate the Ertel PV tendency from Ertel PV flux 
!> \author  Juan A. Saenz
!> \date    January 2014
!> \details 
!>  This subroutine calculates the Ertel PV tendency as the divergence of
!>  the Ertel PV flux
!
!-----------------------------------------------------------------------

   subroutine calculateErtelPVTendencyFromPVFlux(nLayers, nCells, nEdges, &
         mesh, sigmaEA, vectorCell, divVectorCell)!{{{

      use mpas_vector_operations

      integer, intent(in) :: nLayers, nCells, nEdges
      type (mesh_type), intent(in) :: mesh
      real (kind=RKIND), dimension(:,:), intent(in) :: sigmaEA
      real (kind=RKIND), dimension(:,:,:), intent(in) :: vectorCell
      real (kind=RKIND), dimension(:,:), intent(out) :: divVectorCell

      ! local variables
      logical :: includeHalo, on_a_sphere
      integer :: i, k, iComponent
      real (kind=RKIND) :: sigma
      real (kind=RKIND), dimension(:), pointer :: latCell
      real (kind=RKIND), dimension(:), pointer :: lonCell
      integer, dimension(:,:), pointer :: edgeSignOnCell
      real (kind=RKIND), dimension(:,:,:), allocatable :: vectorCellWrk1
      real (kind=RKIND), dimension(:,:,:), allocatable :: vectorCellWrk2
      real (kind=RKIND), dimension(:,:,:), allocatable :: vectorEdgeWrk1

      ! test variables
      real (kind=RKIND) :: wrk
      real (kind=RKIND), dimension(:), allocatable :: divExact
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      xCell => mesh % xCell % array
      yCell => mesh % yCell % array
      zCell => mesh % zCell % array

      if (config_oac_epft_debug) then
        allocate(divExact(nCells+1))
      end if


      allocate(vectorCellWrk1(3,nLayers,nCells+1))
      allocate(vectorCellWrk2(3,nLayers,nCells+1))
      allocate(vectorEdgeWrk1(3,nLayers,nEdges))

      on_a_sphere = mesh % on_a_sphere
      edgeSignOnCell => mesh % edgeSignOnCell % array
      latCell => mesh % latCell % array
      lonCell => mesh % lonCell % array
      
      includeHalo = .true.
 
      vectorCellWrk1 = vectorCell

      ! weight the vector with sigmaEA(:,:)
      do iComponent = 1,3
         vectorCellWrk1(iComponent,:,:) = sigmaEA(:,:)*vectorCellWrk1(iComponent,:,:)
      enddo

         
      if (config_oac_epft_debug) then
        do i= 1,nCells
          vectorCellWrk1(1,:,i) = xCell(i)
          vectorCellWrk1(2,:,i) = yCell(i)
          vectorCellWrk1(3,:,i) = zCell(i)
          ! the analytical divergence:
          divExact(i) = -1.0*sin(lonCell(i)) * &
             (1.0 + 2.0*sin(latCell(i))) !+ 3.0*sin(latCell(i))
        enddo
      endif


      if (on_a_sphere) then

         do i = 1,nCells
            do k = 1,nLayers
               call mpas_vector_LonLatR_to_R3(vectorCellWrk1(:,k,i), &
                  lonCell(i), latCell(i), vectorCellWrk2(:,k,i))
            end do
         end do

         vectorCellWrk1 = vectorCellWrk2

      end if

      call mpas_vector_R3Cell_to_Edge(vectorCellWrk1, mesh, &
         vectorEdgeWrk1)

      call mpas_divergence_in_r3_buoyancy(vectorEdgeWrk1, mesh, edgeSignOnCell, &
         includeHalo, divVectorCell)

         
      if (config_oac_epft_debug) then
        print *, ' '
        do k= 1,nLayers
           wrk = sqrt( &
             sum( &
             ( & 
             (divExact(:)-divVectorCell(k,:))/ max(abs(divExact(:)),1.0e-15)  &
             )**2 * &
             (1.0 - mesh % boundaryCell % array(1,:)) &
             ) / nCells )
           print *, 'div RMS relative error on layer:', wrk
        enddo
      endif


      if (.not. config_oac_epft_debug) then
         do i = 1,nCells
            do k = 1,nLayers
              sigma = max(sigmaEA(k,i), 1.0e-15)
              !sigma = 1.0
              divVectorCell(k,i) = divVectorCell(k,i) / sigma
            end do
         end do
      end if


      deallocate(vectorCellWrk1)
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
   
   subroutine computeErtelPV(nCells, nLayers, nEdges, mesh, &
      fCell, uCell, vCell, sigma, ErtelPV)
      
      use mpas_vector_reconstruction
      
      integer, intent(in) :: nCells, nLayers, nEdges
      type (mesh_type), intent(in) :: mesh
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

      allocate(velNormalGradOnEdge(nLayers, nEdges))
      allocate(velGradX(nLayers, nCells))
      allocate(velGradY(nLayers, nCells))
      allocate(velGradZ(nLayers, nCells))
      allocate(velGradZonal(nLayers, nCells))
      allocate(velGradMerid(nLayers, nCells))
      allocate(vGradZonal(nLayers, nCells))
      allocate(uGradMerid(nLayers, nCells))

      ! calculate derivative of uTWA with respect to y
      call computeNormalGradientOnEdge(nLayers, nCells, nEdges, &
        mesh, &
        uCell, velNormalGradOnEdge)
      call mpas_reconstruct(mesh, velNormalGradOnEdge, &
        velGradX, velGradY, velGradZ, &
        velGradZonal, velGradMerid)
      uGradMerid = velGradMerid

      ! calculate derivative of vTWA with respect to x
      call computeNormalGradientOnEdge(nLayers, nCells, nEdges, &
        mesh, &
        vCell, velNormalGradOnEdge)
      call mpas_reconstruct(mesh, velNormalGradOnEdge, &
        velGradX, velGradY, velGradZ, &
        velGradZonal, velGradMerid)
      vGradZonal = velGradZonal

      do i = 1, nCells
        do k = 1,nLayers 
          ErtelPV(k,i) = (fCell(i) + vGradZonal(k,i) - uGradMerid(k,i))/max(sigma(k,i),1.0e-15)
        end do
      end do

   end subroutine computeErtelPV


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
   subroutine mpas_divergence_in_r3_buoyancy(vectorR3Edge, grid, &
      edgeSignOnCell, includeHalo, divCell)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         vectorR3Edge  !< Input: vector at edge, R3, indices (direction,verticalIndex,edgeIndex)

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

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

      integer :: iEdge, iCell, nCellsCompute, i, k, p, nVertLevels

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND) :: edgeNormalDotVector
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors

      if (includeHalo) then
         nCellsCompute = grid % nCells
      else
         nCellsCompute = grid % nCellsSolve
      endif
      nVertLevels = grid % nBuoyancyLayers

      edgesOnCell        => grid % edgesOnCell % array
      nEdgesOnCell       => grid % nEdgesOnCell % array
      dvEdge             => grid % dvEdge % array
      areaCell           => grid % areaCell % array
      edgeNormalVectors  => grid % edgeNormalVectors % array

      divCell(:,:) = 0.0
      do iCell = 1, nCellsCompute
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, nVertLevels
               edgeNormalDotVector = 0.0
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
  subroutine mpas_vector_R3Cell_to_Edge(vectorCell, mesh, &
        vectorEdge)

      real, dimension(:,:,:), intent(in) :: vectorCell
      type (mesh_type), intent(in) :: mesh !< Input: mesh information
      real (kind=RKIND), dimension(:,:,:), intent(out) :: vectorEdge

      !local variables
      integer :: nEdges, iEdge, k, cell1, cell2
      integer, dimension(:,:), pointer :: cellsOnEdge, boundaryEdge

      cellsOnEdge       => mesh % cellsOnEdge % array
      boundaryEdge      => mesh % boundaryEdge % array

      vectorEdge = 0.0

      do iEdge = 1, mesh % nEdges
         ! Enforce vector value of zero on boundary edges, e.g. no slip for velocities
         if (boundaryEdge(1,iEdge) == 1) then
            vectorEdge(:,:,iEdge) = 0.0
         else
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            do k = 1, mesh % nBuoyancyLayers
                vectorEdge(:,k,iEdge) = 0.5*( vectorCell(:,k,cell2) + vectorCell(:,k,cell1) )
            enddo
         end if
      enddo

   end subroutine mpas_vector_R3Cell_to_Edge!}}}


end module ocn_eliassen_palm_flux_tensor

! vim: foldmethod=marker
