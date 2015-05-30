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
!> \date   May 2015
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
      ! input variables
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

      ! turn on package for this analysis member based on configure option
      ! (at present, this routine is only called when true)
      amEPFTACtive = .false.
      if (config_use_epft) amEPFTACtive = .true.

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

      integer :: nBuoyancyLayers
      real (kind=RKIND), dimension(:), pointer :: potentialDensityMidRef
      real (kind=RKIND), dimension(:), pointer :: potentialDensityTopRef
      real (kind=RKIND), dimension(:), pointer :: buoyancyMidRef
      real (kind=RKIND), dimension(:), pointer :: buoyancyInterfaceRef

      logical, pointer :: amEPFTActive, config_do_restart, config_epft_reset
      integer, pointer :: config_epft_nBuoyancyLayers
      real (kind=RKIND), pointer :: config_epft_rhomax_buoycoor 
      real (kind=RKIND), pointer :: config_epft_rhomin_buoycoor
      real (kind=RKIND), pointer :: config_density0

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
      real (kind=RKIND), dimension(:,:), pointer :: uusigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: vvsigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: uvsigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: uwsigmaEA
      real (kind=RKIND), dimension(:,:), pointer :: vwsigmaEA

      err = 0

      call mpas_pool_get_package(domain % packages, 'amEPFTActive', amEPFTActive)

      if(.not.amEPFTActive) return

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(domain % configs, 'config_epft_reset', config_epft_reset)
      call mpas_pool_get_config(domain % configs, 'config_epft_nBuoyancyLayers', config_epft_nBuoyancyLayers)
      call mpas_pool_get_config(domain % configs, 'config_epft_rhomax_buoycoor', config_epft_rhomax_buoycoor)
      call mpas_pool_get_config(domain % configs, 'config_epft_rhomin_buoycoor', config_epft_rhomin_buoycoor)
      call mpas_pool_get_config(domain % configs, 'config_density0', config_density0)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'amEPFT', amEPFTPool)

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
         nBuoyancyLayers = config_epft_nBuoyancyLayers
         deltaDensity = (config_epft_rhomax_buoycoor - config_epft_rhomin_buoycoor) / config_epft_nBuoyancyLayers
         deltaBuoyancy = -gravity * deltaDensity / config_density0

         !-----------------------------------------------------------------
         ! compute density/bouyancy at top of each layer
         !-----------------------------------------------------------------
         do k = 1, nBuoyancyLayers
            potentialDensityTopRef(k) = config_epft_rhomin_buoycoor + deltaDensity * (k-1)
            buoyancyInterfaceRef(k) = -gravity * (config_epft_rhomin_buoycoor - config_density0) / config_density0 + deltaBuoyancy * (k-1)
         end do
         k=nBuoyancyLayers
         buoyancyInterfaceRef(k+1) = buoyancyInterfaceRef(k) + deltaBuoyancy

         !-----------------------------------------------------------------
         ! compute density/bouyancy for each layer
         !-----------------------------------------------------------------
         do k = 1, nBuoyancyLayers-1
            potentialDensityMidRef(k) = 0.5*(potentialDensityTopRef(k) + potentialDensityTopRef(k+1))
            buoyancyMidRef(k) = 0.5*(buoyancyInterfaceRef(k) + buoyancyInterfaceRef(k+1))
         end do
         k=nBuoyancyLayers
         potentialDensityMidRef(k) = 0.5*(potentialDensityTopRef(k-1) + config_epft_rhomax_buoycoor)
         buoyancyMidRef(k) = 0.5*(buoyancyInterfaceRef(k) + buoyancyInterfaceRef(k+1))
         
         !-----------------------------------------------------------------
         ! initialize ensemble averages when it is not a restart or when a reset is specified
         !-----------------------------------------------------------------
         if (.not. config_do_restart .or. config_epft_reset) then
            call mpas_pool_get_array(amEPFTPool, 'buoyancyMaskEA', buoyancyMaskEA)
            call mpas_pool_get_array(amEPFTPool, 'sigmaEA', sigmaEA)
            call mpas_pool_get_array(amEPFTPool, 'nSamplesEA', nSamplesEA)
            call mpas_pool_get_array(amEPFTPool, 'heightMidBuoyCoorEA', heightMidBuoyCoorEA)
            call mpas_pool_get_array(amEPFTPool, 'montgPotBuoyCoorEA', montgPotBuoyCoorEA)
            call mpas_pool_get_array(amEPFTPool, 'montgPotGradZonalEA', montgPotGradZonalEA)
            call mpas_pool_get_array(amEPFTPool, 'montgPotGradMeridEA', montgPotGradMeridEA)
            call mpas_pool_get_array(amEPFTPool, 'heightMidBuoyCoorSqEA', heightMidBuoyCoorSqEA)
            call mpas_pool_get_array(amEPFTPool, 'heightMGradZonalEA', heightMGradZonalEA)
            call mpas_pool_get_array(amEPFTPool, 'heightMGradMeridEA', heightMGradMeridEA)
            call mpas_pool_get_array(amEPFTPool, 'usigmaEA', usigmaEA)
            call mpas_pool_get_array(amEPFTPool, 'vsigmaEA', vsigmaEA)
            call mpas_pool_get_array(amEPFTPool, 'uusigmaEA', uusigmaEA)
            call mpas_pool_get_array(amEPFTPool, 'vvsigmaEA', vvsigmaEA)
            call mpas_pool_get_array(amEPFTPool, 'uvsigmaEA', uvsigmaEA)
            call mpas_pool_get_array(amEPFTPool, 'uwsigmaEA', uwsigmaEA)
            call mpas_pool_get_array(amEPFTPool, 'vwsigmaEA', vwsigmaEA)

            buoyancyMaskEA = 0.0
            sigmaEA = 0.0
            nSamplesEA = 0.0
            heightMidBuoyCoorEA = 0.0
            montgPotBuoyCoorEA = 0.0
            montgPotGradZonalEA = 0.0
            montgPotGradMeridEA = 0.0
            heightMidBuoyCoorSqEA = 0.0
            heightMGradZonalEA = 0.0
            heightMGradMeridEA = 0.0
            usigmaEA = 0.0
            vsigmaEA = 0.0
            uusigmaEA = 0.0
            vvsigmaEA = 0.0
            uvsigmaEA = 0.0
            uwsigmaEA = 0.0
            vwsigmaEA = 0.0
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
      integer, pointer :: nVertLevels, nBuoyancyLayers, nBuoyLayersP1
      integer, pointer :: nEdges, nCells, nCellsSolve ! nCellsSolve includes halos
      

      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:), pointer :: firstLayerBuoyCoor
      integer, dimension(:), pointer :: lastLayerBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: buoyancyMask

      integer :: nSamplesEA
      
      real(KIND=RKIND), dimension(:), pointer :: potentialDensityMidRef
      real(KIND=RKIND), dimension(:), pointer :: potentialDensityTopRef
      real(KIND=RKIND), dimension(:), pointer :: buoyancyMidRef
      real(KIND=RKIND), dimension(:), pointer :: buoyancyInterfaceRef
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
         call mpas_pool_get_dimension(block % dimensions, 'nBuoyancyLayers', nBuoyLayers)
         !call mpas_pool_get_dimension(block % dimensions, 'nBuoyancyLayersP1', nBuoyLayersP1)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'fCell', fCell)


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
         call mpas_pool_get_field(scratchPool, 'buoyancyMask', buoyancyMask)
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
         call mpas_allocate_scratch_field(buoyancyMask, .true.)
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
         call mpas_pool_get_array(am_epftPool, 'buoyancyMaskEA', buoyancyMaskEA)
         call mpas_pool_get_array(am_epftPool, 'sigmaEA', sigmaEA)
         call mpas_pool_get_array(am_epftPool, 'heightMidBuoyCoorEA', heightMidBuoyCoorEA)
         call mpas_pool_get_array(am_epftPool, 'montgPotGradZonalEA', montgPotGradZonalEA)
         call mpas_pool_get_array(am_epftPool, 'montgPotGradMeridEA', montgPotGradMeridEA)
         call mpas_pool_get_array(am_epftPool, 'heightMidBuoyCoorSqEA', heightMidBuoyCoorSqEA)
         call mpas_pool_get_array(am_epftPool, 'montgPotBuoyCoorEA', montgPotBuoyCoorEA)
         call mpas_pool_get_array(am_epftPool, 'heightMGradZonalEA', HeightMGradZonalEA)
         call mpas_pool_get_array(am_epftPool, 'heightMGradMeridEA', HeightMGradMeridEA)
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
         buoyancyMask                => buoyancyMask % array
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

         
         nBuoyancyLayersP1 = nBuoyLayers+1

         
         ! jas diabatic terms
         !diabaticHeating(nVertLevels,nCells)! "vertical velocity" in buoyancy space
         !wCellCenter       = 0.0

         !jas issue
         ! Get diabaticTimeTendency of a buoyancy surface, omega with funny hat, if any.
         !call any existing MPAS-O subroutines for this


         !-------------------------------------------------------------
         ! begin computation
         !-------------------------------------------------------------

         call get_masks_in_buoyancy_coordinates(nVertLevels, nCells, nBuoyancyLayers, &
           maxLevelCell, potentialDensity, potentialDensityMidRef, &
           firstLayerBuoyCoor, lastLayerBuoyCoor, buoyancyMask)

         if(config_oac_epft_debug) then
            print *, ' '
            print *, 'timeLevel:', timeLevel
            print *, ' '
            print *, 'potentialDensityTopRef'
            print *, potentialDensityTopRef
            print *, 'potentialDensityMidRef'
            print *, potentialDensityMidRef
            print *, 'nCells*nBuoyancyLayers', nCells*nBuoyLayers
            print *, 'sum(buoyancyMask)', sum(buoyancyMask)
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
            call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
               maxLevelCell, array1_3D, zMid, potentialDensityMidRef, array1_3Dbuoy)
            
            do i = 1,nCells
               do k = 1, nBuoyancyLayers
                  array2_3Dbuoy(k,i) = zMid(1,i) + &
                     (potentialDensityMidRef(k) - config_rhomin_buoycoor*1.02) * &
                     (zMid(nVertLevels,i) - zMid(1,i)) / &
                     (config_rhomax_buoycoor*0.98 - config_rhomin_buoycoor*1.02)
               end do
            end do
            do i = 1,nCells
               do k = 1, nBuoyancyLayers
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

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, -potentialDensity, zMid, &
            -potentialDensityMidRef, heightMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, -potentialDensity, zMid, &
            -potentialDensityTopRef, heightTopBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, -potentialDensity, uCellCenter, &
            -potentialDensityMidRef, uMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, -potentialDensity, vCellCenter, &
            -potentialDensityMidRef, vMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, -potentialDensity, density, &
            -potentialDensityMidRef, densityMidBuoyCoor)

         call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
            maxLevelCell, -potentialDensity, density, &
            -potentialDensityTopRef, densityTopBuoyCoor)

         !call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
         !   maxLevelCell, -potentialDensity, Q, potentialDensityTopRef, QMidRef)


         ! JAS this can be moved to init ?
         call computeBuoyancyColumn(nBuoyancyLayers, potentialDensityMidRef, buoyancyMidRef)
         call computeBuoyancyColumnP1(nBuoyancyLayersP1, potentialDensityTopRef, &
            buoyancyInterfaceRef)

         !-------------------------------------------------------------
         ! fill in data above firstLayerBuoyCoor and below lastLayerBuoyCoor
         !-------------------------------------------------------------
         do i = 1, nCells
            do k = 1, firstLayerBuoyCoor(i)-1
              heightMidBuoyCoor(k,i) = zTop(1,i)
              heightTopBuoyCoor(k,i) = zTop(1,i)
              uMidBuoyCoor(k,i) = normalVelocityZonal(1,i)
              vMidBuoyCoor(k,i) = normalVelocityMeridional(1,i)
              densityMidBuoyCoor(k,i) = density(1,i)
              densityTopBuoyCoor(k,i) = density(1,i)
              ! TDR: diabatic
              !wMidBuoyCoor(k,i) = wCellCenter(1,i)
            end do
            do k = lastLayerBuoyCoor(i) + 1, nBuoyancyLayers
              heightMidBuoyCoor(k,i) = -bottomDepth(i)
              heightTopBuoyCoor(k,i) = -bottomDepth(i)
              uMidBuoyCoor(k,i) = normalVelocityZonal(maxLevelCell(i),i)
              vMidBuoyCoor(k,i) = normalVelocityMeridional(maxLevelCell(i),i)
              densityMidBuoyCoor(k,i) = density(maxLevelCell(i),i)
              densityTopBuoyCoor(k,i) = density(maxLevelCell(i),i)
              ! TDR: diabatic
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
         call computeMontgomeryPotential(nBuoyancyLayers, nCells, surfacePressure, &
            firstLayerBuoyCoor, lastLayerBuoyCoor, SSH, densityMidBuoyCoor, &
            potentialDensityMidRef, heightInterfaceBuoyCoor, montgPotBuoyCoor)
         
         !-------------------------------------------------------------
         ! compute the normal derivative of Montgomery potential at cell edges
         !-------------------------------------------------------------
        call computeNormalGradientOnEdge(nBuoyancyLayers, nCells, nEdges, &
            mesh, &
            montgPotBuoyCoor, montgPotNormalGradOnEdge)
         
         !-------------------------------------------------------------
         ! reconstruct full gradient vector at cell centers
         !-------------------------------------------------------------
        call mpas_reconstruct(mesh, montgPotNormalGradOnEdge, &
            montgPotGradX, montgPotGradY, montgPotGradZ, &
            montgPotGradZonal, montgPotGradMerid)

         !-------------------------------------------------------------
         ! Increment first-order running mean fields
         !-------------------------------------------------------------
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            buoyancyMask, buoyancyMaskEA)
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            sigma, sigmaEA)
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            heightMidBuoyCoor, heightMidBuoyCoorEA)
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            montgPotBuoyCoor, montgPotBuoyCoorEA)
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            montgPotGradZonal, montgPotGradZonalEA)
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            montgPotGradMerid, montgPotGradMeridEA)
   
         !-------------------------------------------------------------
         ! Increment second-order running mean fields
         !-------------------------------------------------------------
         wrk3DBuoyCoor = heightMidBuoyCoor * heightMidBuoyCoor
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, heightMidBuoyCoorSqEA)

         wrk3DBuoyCoor = heightMidBuoyCoor * montgPotGradZonal
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, heightMGradZonalEA)

         wrk3DBuoyCoor = heightMidBuoyCoor * montgPotGradMerid
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, heightMGradMeridEA)

         wrk3DBuoyCoor = uMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, usigmaEA)

         wrk3DBuoyCoor = vMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, vsigmaEA)

         ! Diabatic terms
         !wrk3DBuoyCoor = wMidBuoyCoor * sigma
         !call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
         !   wrk3DBuoyCoor, wsigmaEA)
   
   
         !-------------------------------------------------------------
         ! Increment third-order running mean fields
         !-------------------------------------------------------------
         wrk3DBuoyCoor = uMidBuoyCoor * uMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, uusigmaEA)

         wrk3DBuoyCoor = vMidBuoyCoor * vMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, vvsigmaEA)

         wrk3DBuoyCoor = uMidBuoyCoor * vMidBuoyCoor * sigma
         call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
            wrk3DBuoyCoor, uvsigmaEA)

         ! Diabatic terms
         !wrk3DBuoyCoor = uMidBuoyCoor * wMidBuoyCoor * sigma
         !call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
         !   wrk3DBuoyCoor, uwsigmaEA)
         uwsigmaEA = 0.0

         ! Diabatic terms
         !wrk3DBuoyCoor = vMidBuoyCoor * wMidBuoyCoor* sigma
         !call updateEnsembleAverage(nBuoyancyLayers, nCells, nSamplesEA, &
         !   wrk3DBuoyCoor, vwsigmaEA)
         vwsigmaEA = 0.0
   
         !-------------------------------------------------------------
         ! update number of samples in ensemble average
         !-------------------------------------------------------------
         nSamplesEA = nSamplesEA + 1
   
         !-------------------------------------------------------------
         ! based on current estimate of ensemble-average state,
         ! compute the thickness-weighted average velocity
         !-------------------------------------------------------------
         call calculateTWA(nBuoyancyLayers, nCells, nBuoyLayers, &
            sigmaEA, usigmaEA,  uTWA)
         call calculateTWA(nBuoyancyLayers, nCells, nBuoyLayers, &
            sigmaEA, vsigmaEA,  vTWA)
         ! Diabatic terms
         !call calculateTWA(nBuoyancyLayers, nCells, nBuoyLayers, &
         !   sigmaEA, wsigmaEA,  wTWA)
         wTWA = 0.0

         !-------------------------------------------------------------
         ! based on current estimate of ensemble-average state,
         ! compute the Eliassen-Palm flux tensor
         !-------------------------------------------------------------
         call calculateEPFTfromTWA(nBuoyancyLayers, nCells, &
            sigmaEA, heightMidBuoyCoorEA, &
            heightMidBuoyCoorSqEA, montgPotGradZonalEA, montgPotGradMeridEA, &
            heightMGradZonalEA, heightMGradMeridEA, uTWA, vTWA, wTWA, &
            uusigmaEA, vvsigmaEA, uvsigmaEA, uwsigmaEA, vwsigmaEA, EPFT)

         !-------------------------------------------------------------
         ! compute the force applied to the momentum equation as div(EPFT)
         !-------------------------------------------------------------
        call calculateDivEPFT(nBuoyancyLayers, nCells, nEdges, &
            mesh, buoyancyInterfaceRef, sigmaEA, buoyancyMaskEA, EPFT, divEPFT)
         
         !-------------------------------------------------------------
         ! transform div(EPFT) into a flux of Ertel's PV
         !-------------------------------------------------------------
        call calculateErtelPVFlux(nCells, nBuoyancyLayers, &
            sigmaEA, divEPFT, ErtelPVFlux)

         !-------------------------------------------------------------
         ! compute div(ErtelPVFlux) to obtain tendency of Ertel's PV
         !-------------------------------------------------------------
        call calculateErtelPVTendencyFromPVFlux(nBuoyancyLayers, nCells, nEdges, &
            mesh, sigmaEA, ErtelPVFlux, ErtelPVTendency)


         !-------------------------------------------------------------
         ! compute Ertel PV based on EA/TWA fields
         !-------------------------------------------------------------
         call computeErtelPV(nCells, nBuoyancyLayers, nEdges, mesh, &
            fCell, uTWA, vTWA, sigmaEA, ErtelPV)

         !-------------------------------------------------------------
         ! Compute the geometric decomposition in terms of angles and 
         ! eccentricities using the entries of EPFT.
         ! (not yet implemented)
         !-------------------------------------------------------------
         !call eddyGeomDecompEPFT(EPFT, ...)




         ! calculate potential vorticity fluxes using curl of u
         if(config_oac_epft_debug) then

            relativeVorticityCell => diagnostics % relativeVorticityCell % array
            
            ! store relVortMidBuoyCoor in array1_3Dbuoy
            call linear_interp_1d_field_along_column(nVertLevels, nCells, nBuoyancyLayers, &
               maxLevelCell, -potentialDensity, relativeVorticityCell, &
               -potentialDensityMidRef, array1_3Dbuoy)
            
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
         ! deallocate scratch space
         !-------------------------------------------------------------
         call mpas_deallocate_scratch_field(firstLayerBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(lastLayerBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(buoyancyMaskField, .true.)
         call mpas_deallocate_scratch_field(sigmaField, .true.)
         call mpas_deallocate_scratch_field(heightMidBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(heightTopBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(heightInterfaceBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(uMidBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(vMidBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(densityMidBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(densityTopBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(montgPotBuoyCoorField, .true.)
         call mpas_deallocate_scratch_field(montgPotNormalGradOnEdgeField, .true.)
         call mpas_deallocate_scratch_field(montgPotGradXField, .true.)
         call mpas_deallocate_scratch_field(montgPotGradYField, .true.)
         call mpas_deallocate_scratch_field(montgPotGradZField, .true.)
         call mpas_deallocate_scratch_field(montgPotGradZonalField, .true.)
         call mpas_deallocate_scratch_field(montgPotGradMeridField, .true.)
         call mpas_deallocate_scratch_field(wrk3DnVertLevelsField, .true.)
         call mpas_deallocate_scratch_field(wrk3DBuoyCoorField, .true.)
         
         call mpas_deallocate_scratch_field(array1_3D, .true.)
         call mpas_deallocate_scratch_field(array2_3D, .true.)
         call mpas_deallocate_scratch_field(array3_3D, .true.)
         call mpas_deallocate_scratch_field(array1_3Dbuoy, .true.)
         call mpas_deallocate_scratch_field(array2_3Dbuoy, .true.)

         call mpas_deallocate_scratch_field(PVMidBuoyCoor, .true.)
         call mpas_deallocate_scratch_field(PVMidBuoyCoorEA, .true.)
         call mpas_deallocate_scratch_field(uPVMidBuoyCoorEA , .true.)
         call mpas_deallocate_scratch_field(vPVMidBuoyCoorEA, .true.)
         call mpas_deallocate_scratch_field(PVFluxTest, .true.)


         nCellsCum = nCellsCum + nCells

         !-------------------------------------------------------------
         ! move to the next block
         !-------------------------------------------------------------
      
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
!> \author  Juan A. Saenz
!> \date    May 2015
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
      buoyancyMask     = 0.0

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
           buoyancyMask(kB,iCell) = 1.0
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

      real (kind=RKIND), pointer :: config_epft_rhomin_buoycoor, config_epft_rhomax_buoycoor

      call mpas_pool_get_config(ocnConfigs, 'config_epft_rhomin_buoycoor', config_epft_rhomin_buoycoor)
      call mpas_pool_get_config(ocnConfigs, 'config_epft_rhomax_buoycoor', config_epft_rhomax_buoycoor)

      printWarning = .false.
      iCellMinBound = -1
      iCellMaxBound = -1

      do iCell = 1, nCells
         if (potentialDensity(1,iCell) < config_epft_rhomin_buoycoor) then
            printWarning = .true.
            iCellMinBound = iCell
            exit
         end if
         if (potentialDensity(maxLevelCell(iCell),iCell) > config_epft_rhomax_buoycoor) then
            printWarning = .true.
            iCellMaxBound = iCell
            exit
         end if
      enddo
      
      if (printWarning) then
         write(stderrUnit,*)
         write(stderrUnit,*) 'Warning: in EPFT package, subroutine check_potentialDensityRef_range'
         write(stderrUnit,*) 'One or more columns in the ocean doman have densities that are not'
         write(stderrUnit,*) 'contained in the defined buoyancy space of the EPFT module'
         if (iCellMinBound.gt.0) write(stderrUnit,*) 'fluid is lighter than min buoyancy at cell: ',iCellMinBound
         if (iCellMaxBound.gt.0) write(stderrUnit,*) 'fluid is lighter than max buoyancy at cell: ',iCellMaxBound
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
!>  Required: xFieldIn monotonically decreases with index value
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

      !-----------------------------------------------------------------
      ! initialize intent(out)
      !-----------------------------------------------------------------
      yFieldOut = 0.0

      !-----------------------------------------------------------------
      ! loop over all columns
      !-----------------------------------------------------------------
      do iCell = 1, nCells
      
         ! find the index of the bottom level of a column
         maxLevel = maxLevelCell(iCell)

         ! Monotonically decreasing xFieldIn required
         ! Find index of first element in xColumnOut that is inside xFieldIn(:,iCell)
         kBTop = 1
         do kB = 1, nBuoyancyLayers
            ! the following line ensures that
            ! if all xColumnOut > xFieldIn(1,iCell) then kBTop = nBuoyancyLayers
            kBTop = kB
            if (xColumnOut(kB) <= xFieldIn(1,iCell) ) then
               exit
            endif
         enddo
         
         !find last target buoyancy level inside column
         kBBottom = nBuoyancyLayers
         do kB = nBuoyancyLayers, 1, -1
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
         !set data from kBBottom+1:nBuoyancyLayers to bottom values
         do kB = kBBottom+1, nBuoyancyLayers
           yFieldOut(kB,iCell) = yFieldIn(maxLevel,iCell)
         enddo

         ! The interpolation: 
         ! for the target buoyancy levels within the buoyancy range in a column:
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
      
      
      ! local variables
      integer :: iCell, k
      
      !-----------------------------------------------------------------
      ! initialize sigma assuming zero thickness layers everywhere
      !-----------------------------------------------------------------
      sigma = 0.0
      
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
      MontgomeryPotential = 0.0
      
      !-----------------------------------------------------------------
      ! loop over all columns
      !-----------------------------------------------------------------
      do iCell = 1, nCells
         
         !-----------------------------------------------------------------
         ! compute Montgomery potential in top buoyancy layer
         ! at present, assume atmosphere surface pressure is zero (or a constant)
         !-----------------------------------------------------------------
         pInterfacek = 0.0
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
      normalGradOnEdge = 0.0

      !-----------------------------------------------------------------
      ! loop over edges, compute derivative as (cell2 - cell1) / dc
      !-----------------------------------------------------------------
      do iEdge = 1, nEdges
         ! do not compute the normal derivative at land/sea interface
         if (boundaryEdge(1,iEdge) == 1) then
            normalGradOnEdge(:,iEdge) = 0.0
         else
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            invLength = 1.0 / dcEdge(iEdge)
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
            Abar(k,iCell) = (nSamples * Abar(k,iCell) + A(k,iCell)) / (nSamples + 1.0)
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
      varTWA = 0.0
      
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
            
            sigma = max(sigmaEA(kLayer,iCell), epsilonEPFT)

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
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details 
!>  This subroutine calculates the divergence of the Eliassen-Palm flux tensor
!
!-----------------------------------------------------------------------

   subroutine calculateDivEPFT(onASphere, nLayers, nCells, nEdges, &
        meshPool, buoyancyMidRef, sigmaEA, buoyancyMaskEA, tensorCellIn, vectorCellOut)!{{{

      use mpas_vector_operations

      logical, intent(in) :: onASphere
      integer, intent(in) :: nLayers, nCells, nEdges
      type (mpas_pool_type), intent(in) :: meshPool
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
      allocate(vectorEdgeWrk1(3,nLayers,nEdges+1)) 
      allocate(vertVector(nLayers))

      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      
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
            if (q.eq.3 .and. config_oac_epft_debug) then
               vertVector(:) = 0.0
            endif


            do kLayer = 1,nLayers

               wrk = 0.0

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

            end do  ! do iCell=1,nCells

         end do  ! do q=1,3

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
            ErtelPVFlux(:,k,i) = ErtelPVFlux(:,k,i) / max(sigma(k,i),epsilonEPFT)
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
      integer :: nEdgesCompute, nBuoyancyLayers, nCells
      integer, dimension(:,:), pointer :: cellsOnEdge

      if (includeHalo) then
         nEdgesCompute = grid % nEdges
      else 
         nEdgesCompute = grid % nEdgesSolve
      endif
      nBuoyancyLayers = grid % nBuoyancyLayers
      nCells = grid % nCells

      cellsOnEdge => grid % cellsOnEdge % array

      ! error check that index 1 of matrixEdge and matrixCell are same length?

      do iEdge=1,nEdgesCompute
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nBuoyancyLayers
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
!> \author  Juan A. Saenz, Todd Ringler
!> \date    May 2015
!> \details 
!>  This subroutine calculates the Ertel PV tendency as the divergence of
!>  the Ertel PV flux
!
!-----------------------------------------------------------------------

   subroutine calculateErtelPVTendencyFromPVFlux(onASphere, nLayers, nCells, nEdges, &
         meshPool, sigma, vectorCell, divVectorCell)!{{{

      use mpas_vector_operations

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
      xCell => mesh % xCell % array
      yCell => mesh % yCell % array
      zCell => mesh % zCell % array

      if (config_oac_epft_debug) then
        allocate(divExact(nCells+1))
      end if


      allocate(vectorCellWrk1(3,nLayers,nCells+1))
      allocate(vectorCellWrk2(3,nLayers,nCells+1))
      allocate(vectorEdgeWrk1(3,nLayers,nEdges+1)) !jas issua check that this is +1

      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      
      includeHalo = .true.
 
      ! copy vector into work array
      vectorCellWrk1 = vectorCell

      ! weight the vector with sigma(:,:)
      do iComponent = 1,3
         vectorCellWrk1(iComponent,:,:) = sigma(:,:)*vectorCellWrk1(iComponent,:,:)
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
              divVectorCell(k,i) = divVectorCell(k,i) / max(sigma(k,i), epsilonEPFT)
            end do
         end do
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

      allocate(velNormalGradOnEdge(nLayers, nEdges+1)) ! jas issue check all these are +1
      allocate(velGradX(nLayers, nCells+1))
      allocate(velGradY(nLayers, nCells+1))
      allocate(velGradZ(nLayers, nCells+1))
      allocate(velGradZonal(nLayers, nCells+1))
      allocate(velGradMerid(nLayers, nCells+1))
      allocate(vGradZonal(nLayers, nCells+1))
      allocate(uGradMerid(nLayers, nCells+1))

      ! calculate derivative of uTWA with respect to the meridional direction
      call computeNormalGradientOnEdge(nLayers, nCells, nEdges, &
        meshPool, uCell, velNormalGradOnEdge)
      call mpas_reconstruct(meshPool, velNormalGradOnEdge, &
        velGradX, velGradY, velGradZ, &
        velGradZonal, velGradMerid)
      uGradMerid = velGradMerid

      ! calculate derivative of vTWA with respect to the zonal direction
      call computeNormalGradientOnEdge(nLayers, nCells, nEdges, &
        meshPool, vCell, velNormalGradOnEdge)
      call mpas_reconstruct(meshPool, velNormalGradOnEdge, &
        velGradX, velGradY, velGradZ, &
        velGradZonal, velGradMerid)
      vGradZonal = velGradZonal

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

      divCell(:,:) = 0.0
      do iCell = 1, nCells
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
  subroutine mpas_vector_R3Cell_to_Edge(vectorCell, meshPool, &
        vectorEdge)

      real, dimension(:,:,:), intent(in) :: vectorCell
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

      vectorEdge = 0.0

      do iEdge = 1, nEdges
         ! Enforce vector value of zero on boundary edges, e.g. no slip for velocities
         if (boundaryEdge(1,iEdge) == 1) then
            vectorEdge(:,:,iEdge) = 0.0
         else
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            do k = 1, nBuoyancyLayers
                vectorEdge(:,k,iEdge) = 0.5*( vectorCell(:,k,cell2) + vectorCell(:,k,cell1) )
            enddo
         end if
      enddo

   end subroutine mpas_vector_R3Cell_to_Edge!}}}


end module ocn_eliassen_palm_flux_tensor

! vim: foldmethod=marker
