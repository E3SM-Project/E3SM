! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  oac_epft
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

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_setup_packages_eliassen_palm_flux_tensor
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
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
!> \author  Juan A. Saenz
!> \date    May 2015
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_eliassen_palm_flux_tensor(domain, err)!{{{

      use mpas_packages

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
!> \author  Juan A. Saenz
!> \date    May 2015
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_eliassen_palm_flux_tensor(domain, timeLevel, err)!{{{

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

      type (mpas_pool_type), pointer :: am_eliassen_palm_flux_tensorPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: am_eliassen_palm_flux_tensor

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nBuoyLayers, nBuoyLayersP1
      integer, pointer :: nEdges, nCells, nCellsSolve, nCellsCum ! nCellsSolve includes halos

      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:), pointer :: firstLayerBuoyCoor
      integer, dimension(:), pointer :: lastLayerBuoyCoor
      real(KIND=RKIND), dimension(:,:), pointer :: buoyMask

      integer :: nSamplesEA
      real(KIND=RKIND), dimension(:,:), pointer :: sigmaEA
      real(KIND=RKIND), dimension(:,:), pointer :: heightMidBuoyCoorEA
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotBuoyCoorEA
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradZonalEA
      real(KIND=RKIND), dimension(:,:), pointer :: montgPotGradMeridEA
      real(KIND=RKIND), dimension(:,:), pointer :: heightMidBuoyCoorSqEA
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


      real(KIND=RKIND), dimension(:), pointer :: potentialDensityMidRef
      real(KIND=RKIND), dimension(:), pointer :: potentialDensityTopRef
      real(KIND=RKIND), dimension(:), pointer :: buoyancyMidRef
      real(KIND=RKIND), dimension(:), pointer :: buoyancyInterfaceRef
      real(KIND=RKIND), dimension(:), pointer :: bottomDepth

      real(KIND=RKIND), dimension(:,:), pointer :: buoyMaskEA
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
      real(KIND=RKIND), dimension(:), pointer :: SSH
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
      real(KIND=RKIND) :: RMSlocal1, RMSglobal1
      real(KIND=RKIND) :: RMSlocal2, RMSglobal2
      real(KIND=RKIND) :: RMSPVFlux1local, RMSPVFlux1global
      real(KIND=RKIND) :: RMSPVFlux2local, RMSPVFlux2global
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

      call mpas_timer_start("compute_eliassen_palm_flux_tensor", .false., am_eliassen_palm_flux_tensorTimer)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'am_eliassen_palm_flux_tensor', am_eliassen_palm_flux_tensorPool)

         ! Here are some example variables which may be needed for your analysis member
         call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)

         ! Computations which are functions of nCells, nEdges, or nVertices 
         ! must be placed within this block loop
         ! Here are some example loops
         do iCell = 1,nCellsSolve
            do k = 1, maxLevelCell(iCell)
               do iTracer = 1, num_tracers
               ! computations on tracers(iTracer,k, iCell)
               end do
            end do
         end do

         block => block % next
      end do

      ! mpi gather/scatter calls may be placed here.
      ! Here are some examples.  See mpas_oac_global_stats.F for further details.
!      call mpas_dmpar_sum_real_array(dminfo, nVariables, sumSquares(1:nVariables), reductions(1:nVariables))
!      call mpas_dmpar_min_real_array(dminfo, nMins, mins(1:nMins), reductions(1:nMins))
!      call mpas_dmpar_max_real_array(dminfo, nMaxes, maxes(1:nMaxes), reductions(1:nMaxes))

      ! Even though some variables do not include an index that is decomposed amongst 
      ! domain partitions, we assign them within a block loop so that all blocks have the
      ! correct values for writing output.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'am_eliassen_palm_flux_tensor', am_eliassen_palm_flux_tensorPool)

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
