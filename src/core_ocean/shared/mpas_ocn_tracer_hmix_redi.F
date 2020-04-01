! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hmix_Redi
!
!> \brief MPAS ocean horizontal tracer mixing driver
!> \details
!>  This module contains the main driver routine for computing
!>  horizontal mixing tendencies.
!>
!>  It provides an init and a tend function. Each are described below.
!
!-----------------------------------------------------------------------

module ocn_tracer_hmix_Redi

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading
   use mpas_constants

   use ocn_config
   use ocn_constants

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

   public :: ocn_tracer_hmix_Redi_tend, &
             ocn_tracer_hmix_Redi_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND) :: term1TaperFactor

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_hmix_Redi_tend
!
!> \brief   Computes Laplacian tendency term for horizontal tracer mixing
!> \details
!>  This routine computes the horizontal mixing tendency for tracers
!>  based on current state using a Laplacian parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_Redi_tend(meshPool, scratchPool, layerThicknessEdge, zMid, tracers, &
                    redi_term1, redi_term2, redi_term3,  &
                    RediKappa, dt, isActiveTracer, slopeTriadUp, slopeTriadDown,  &
                    kappaRediCell, kappaRediEdge, rediLimiterCount, tend, err)!{{{


      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: Scratch information

      real (kind=RKIND), dimension(:), intent(in) :: &
         RediKappa

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge, &!< Input: thickness at edge
         zMid,               &!< Input: Z coordinate at the center of a cell
         kappaRediCell        !< Input: vertical structure of GM/Redi limiter

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         slopeTriadUp, slopeTriadDown, &
         tracers, kappaRediEdge !< Input: tracer quantities

      real (kind=RKIND), intent(in) :: dt

      integer, intent(inout) :: rediLimiterCount

      logical, intent(in) :: isActiveTracer

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend,redi_term1, redi_term2, redi_term3     !< Input/Output: velocity tendency

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

      integer :: iCell, iEdge, cell1, cell2, iCellSelf
      integer :: i, k, iTr, nTracers, nCells, nEdges, km1, kp1, nCellsP1
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      integer, dimension(:,:), allocatable :: boundaryMask

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnCell, maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, edgeSignOnCell, cellsOnCell

      real (kind=RKIND) :: tempTracer, invAreaCell1, invAreaCell2, invAreaCell, areaEdge
      real (kind=RKIND) :: flux, flux_term1, flux_term2, flux_term3, flux_terms12, sign_term1, dTracerDx, coef
      real (kind=RKIND) :: r_tmp, tracer_turb_flux, kappaRediEdgeVal
      real (kind=RKIND), dimension(:), pointer :: areaCell, dvEdge, dcEdge
      integer, dimension(:), allocatable :: limitCount
      real (kind=RKIND), dimension(:), allocatable :: areaSum, minimumVal
      real (kind=RKIND), dimension(:,:), allocatable :: fluxRediZTop
      real (kind=RKIND), dimension(:,:,:), allocatable :: redi_term2_edge, redi_term3_topOfCell
      err = 0

      if (.not.config_use_Redi) return

      call mpas_timer_start("tracer redi")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      nTracers = size(tracers, dim=1)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      nCells = nCellsArray( size(nCellsArray) )
      nEdges = nEdgesArray( size(nEdgesArray) )

      allocate(minimumVal(nTracers))
      allocate(limitCount(nTracers))
      allocate(redi_term2_edge(nTracers,nVertLevels,nEdges))
      allocate(redi_term3_topOfCell(nTracers,nVertLevels+1,nCells))

      redi_term2_edge(:,:,:) = 0.0_RKIND
      redi_term3_topOfCell(:,:,:) = 0.0_RKIND

      ! mrp this is in original del2. Apparently we don't do halo cells for del2.
      nCells = nCellsArray(1)

      redi_term1(:,:,:) = 0.0_RKIND
      ! Term 1: this is the "standard" horizontal del2 term, but with RediKappa coefficient.
      ! \kappa_2 \nabla \phi on edge
      ! mrp performance: merge term 1 into main cell loop later
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         invAreaCell = 1.0_RKIND / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            r_tmp = dvEdge(iEdge) / dcEdge(iEdge)

            do k = 1, maxLevelEdgeTop(iEdge)

               kappaRediEdgeVal = 0.25_RKIND * (kappaRediCell(k,cell1) + kappaRediCell(k,cell2) + &
                        kappaRediCell(k+1,cell1) + kappaRediCell(k+1,cell2)) * &
                        (0.5_RKIND*(kappaRediEdge(k,1,iEdge) + kappaRediEdge(k,2,iEdge)))
!                                kappaRediCell(k+1,cell1) + kappaRediCell(k+1,cell2))
               do iTr = 1, nTracers
                  ! \kappa_2 \nabla \phi on edge
                  tracer_turb_flux = tracers(iTr, k, cell2) - tracers(iTr, k, cell1)

                  ! div(h \kappa_2 \nabla \phi) at cell center
                  flux = layerThicknessEdge(k, iEdge) * tracer_turb_flux * r_tmp * RediKappa(iCell) !* (1.0_RKIND + term1TaperFactor * &
!                          (kappaRediEdge - 1.0_RKIND)) * RediKappa(iCell) * r_tmp

                  redi_term1(iTr,k,iCell) = redi_term1(iTr,k,iCell) - (1.0_RKIND + term1TaperFactor * &
                          (kappaRediEdgeVal - 1.0_RKIND)) *edgeSignOnCell(i, iCell) * flux * invAreaCell
                  tend(iTr,k,iCell) = tend(iTr,k,iCell) - (1.0_RKIND + term1TaperFactor * &
                          (kappaRediEdgeVal - 1.0_RKIND)) * edgeSignOnCell(i, iCell) * flux * invAreaCell

               end do
            end do

         end do
      end do
      !$omp end do

      !mrp performance note: can try to reduce halo width later
      nCells = nCellsArray( 1 )
      nCellsP1 = nCellsArray( size(nCellsArray) ) + 1
      redi_term2(:,:,:) = 0.0_RKIND
      ! Term 2: div( h S dphi/dz)
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         invAreaCell = 1.0_RKIND / areaCell(iCell)

         do i = 1, nEdgesOnCell(iCell)
            ! Check if neighboring cell exists
            if  (cellsOnCell(i,iCell).eq.nCellsP1) cycle
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            coef = edgeSignOnCell(i, iCell) * dvEdge(iEdge) * invAreaCell

           k=1
           kappaRediEdgeVal = 0.25_RKIND * (kappaRediCell(k,cell1) + kappaRediCell(k,cell2) + &
                            kappaRediCell(k+1,cell1) + kappaRediCell(k+1,cell2))
           do iTr = 1, nTracers

               ! For top layer, only use triads pointing down:
               flux_term2 = coef * RediKappa(iCell) * layerThicknessEdge(k, iEdge) * &
                  0.25_RKIND * &
                  (  kappaRediCell(k+1,cell1) * slopeTriadDown(k,1,iEdge) * (tracers(iTr,k,cell1) - tracers(iTr,k+1,cell1)) &
                                  / (zMid(k,cell1) - zMid(k+1,cell1)) &
                   + kappaRediCell(k+1,cell2) * slopeTriadDown(k,2,iEdge)  * (tracers(iTr,k,cell2) - tracers(iTr,k+1,cell2)) &
                                  / (zMid(k,cell2) - zMid(k+1,cell2)) ) 

               redi_term2(iTr,k,iCell) = redi_term2(iTr,k,iCell) - flux_term2
               ! Add tendency for Term 2: div( h S dphi/dz) at top layer
               tend(iTr,k,iCell) = tend(iTr,k,iCell) - flux_term2
            end do

            do k = 2, maxLevelEdgeTop(iEdge) - 1
              kappaRediEdgeVal = 0.25_RKIND * (kappaRediCell(k,cell1) + kappaRediCell(k,cell2) + &
                              kappaRediCell(k+1,cell1) + kappaRediCell(k+1,cell2))
              do iTr = 1, nTracers
                  ! There are four triads associated with each edge:
                  flux_term2 = coef * RediKappa(iCell) * layerThicknessEdge(k, iEdge) * &
                     0.25_RKIND * &
                     (  kappaRediCell(k,cell1) * slopeTriadUp(k,1,iEdge) * (tracers(iTr,k-1,cell1) - tracers(iTr,k,cell1)) &
                                     / (zMid(k-1,cell1) - zMid(k,cell1)) &
                      + kappaRediCell(k,cell2) * slopeTriadUp(k,2,iEdge) * (tracers(iTr,k-1,cell2) - tracers(iTr,k,cell2)) &
                                     / (zMid(k-1,cell2) - zMid(k,cell2)) &
                      + kappaRediCell(k+1,cell1) * slopeTriadDown(k,1,iEdge) * (tracers(iTr,k,cell1) - tracers(iTr,k+1,cell1)) &
                                     / (zMid(k,cell1) - zMid(k+1,cell1)) &
                      + kappaRediCell(k+1,cell2) * slopeTriadDown(k,2,iEdge)  * (tracers(iTr,k,cell2) - tracers(iTr,k+1,cell2)) &
                                     / (zMid(k,cell2) - zMid(k+1,cell2))  ) 

                  redi_term2(iTr,k,iCell) = redi_term2(iTr,k,iCell) - flux_term2
                  ! Add tendency for Term 2: div( h S dphi/dz) at mid layers
                  tend(iTr,k,iCell) = tend(iTr,k,iCell) - flux_term2
               end do
            end do

           k=maxLevelEdgeTop(iEdge)
           kappaRediEdgeVal = 0.25_RKIND * (kappaRediCell(k,cell1) + kappaRediCell(k,cell2) + &
                           kappaRediCell(k+1,cell1) + kappaRediCell(k+1,cell2))
           do iTr = 1, nTracers
               ! For bottom layer, only use triads pointing up:
               flux_term2 = coef * RediKappa(iCell) * layerThicknessEdge(k, iEdge) * &
                  0.25_RKIND * &
                  (  kappaRediCell(k,cell1) * slopeTriadUp(k,1,iEdge) * (tracers(iTr,k-1,cell1) - tracers(iTr,k,cell1)) &
                                  / (zMid(k-1,cell1) - zMid(k,cell1)) &
                   + kappaRediCell(k,cell2) * slopeTriadUp(k,2,iEdge) * (tracers(iTr,k-1,cell2) - tracers(iTr,k,cell2)) &
                                  / (zMid(k-1,cell2) - zMid(k,cell2)) ) 

               if ( maxLevelCell(cell1) > maxLevelEdgeTop(iEdge) ) then
                 flux_term2 = flux_term2 + coef * RediKappa(iCell) *  &
                        layerThicknessEdge(k, iEdge) * 0.25_RKIND * kappaRediCell(k+1,cell1) * slopeTriadDown(k,1,iEdge) * &
                        (tracers(iTr,k,cell1) - tracers(iTr,k+1,cell1)) / (zMid(k,cell1) - zMid(k+1,cell1))
               endif
               if ( maxLevelCell(cell2) > maxLevelEdgeTop(iEdge) ) then
                 flux_term2 = flux_term2 + coef * RediKappa(iCell) *  &
                   layerThicknessEdge(k, iEdge) * 0.25_RKIND * kappaRediCell(k+1,cell2) * slopeTriadDown(k,2,iEdge) * &
                   (tracers(iTr,k,cell2) - tracers(iTr,k+1,cell2)) / (zMid(k,cell2) - zMid(k+1,cell2))
               endif
              redi_term2(iTr,k,iCell) = redi_term2(iTr,k,iCell) - flux_term2

              redi_term2_edge(iTr,k,iEdge) = -flux_term2

              ! Add tendency for Term 2: div( h S dphi/dz) at bottom layer
!               tend(iTr,k,iCell) = tend(iTr,k,iCell) - flux_term2
            end do
         end do
      end do
      !$omp end do

      ! Term 3: d/dz ( h S grad phi) 
      ! note that dz and h cancel, so this is (S grad phi)

      allocate(fluxRediZTop(nTracers,nVertLevels+1))
      allocate(areaSum(nVertLevels+1))

      !mrp performance note: can try to reduce halo width later
      nCells = nCellsArray(1)
      !nCells = nCellsArray( size(nCellsArray) )

      redi_term3(:,:,:) = 0.0_RKIND
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         areaSum(:) = 1.0e-34_RKIND
         fluxRediZTop(:,:) = 0.0_RKIND
         invAreaCell = 1.0_RKIND / areaCell(iCell)

         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            areaEdge = dcEdge(iEdge) * dvEdge(iEdge) * 0.25_RKIND
            ! mrp performance note: Replace with a logical array on init later
            if (cell1==iCell) then
                iCellSelf = 1
            elseif (cell2==iCell) then
                iCellSelf = 2
            else
                print *, 'iCellSelf, something went wrong',cell1,cell2,iCell,nCells
            endif

            do k = 2, maxLevelEdgeTop(iEdge)
!               areaSum(k) = areaSum(k) + areaEdge
!               areaSum(k+1) = areaSum(k+1) + areaEdge

               kappaRediEdgeVal = 0.25_RKIND * (kappaRediCell(k,cell1) + kappaRediCell(k,cell2) + &
                           kappaRediCell(k+1,cell1) + kappaRediCell(k+1,cell2))

               do iTr = 1, nTracers
                  ! Gradient of tracer at edge.
                  ! dvEdge on following line is equal to areaEdge/dcEdge
                  dTracerDx = (tracers(iTr,k,cell2) - tracers(iTr,k,cell1)) * dvEdge(iEdge) * 0.25_RKIND
                  fluxRediZTop(iTr,k) = fluxRediZTop(iTr,k) + &
                     0.5_RKIND*(slopeTriadUp(k,iCellSelf,iEdge) + slopeTriadDown(k-1,iCellSelf,iEdge)) &
                            * dTracerDx * kappaRediEdgeVal * invAreaCell                  
                          !                  fluxRediZTop(iTr,k+1) = fluxRediZTop(iTr,k+1) + &
!                     slopeTriadDown(k,iCellSelf,iEdge) * dTracerDx * kappaRediEdgeVal
               end do
            end do

            do k = 2, maxLevelEdgeTop(iEdge)
              redi_term3_topOfCell(:,k,iCell) = fluxRediZTop(:,k)
            enddo
            redi_term3_topOfCell(:,maxLevelCell(iCell)+1,iCell) = 0.0_RKIND
         end do

         ! impose no-flux boundary conditions at top and bottom of column
         fluxRediZTop(:,1) = 0.0_RKIND
!         fluxRediZTop(:,maxLevelCell(iCell)) = 0.0_RKIND
         fluxRediZTop(:,maxLevelCell(iCell)+1) = 0.0_RKIND
         do k = 1, maxLevelCell(iCell)
            do iTr = 1, nTracers
               ! Add tendency for Term 3: d/dz ( h S grad phi) = ( S grad phi) fluxes
               ! 2.0 in next line is because a dot product on a C-grid
               ! requires a factor of 1/2 to average to the cell center.
               flux_term3 = RediKappa(iCell) * 2.0_RKIND * &
                   (fluxRediZTop(iTr,k) - fluxRediZTop(iTr,k+1))

               redi_term3(iTr,k,iCell) = redi_term3(iTr,k,iCell) + flux_term3
!               tend(iTr,k,iCell) = tend(iTr,k,iCell) + flux_term3
            end do
         end do
      end do
      !$omp end do
      deallocate(fluxRediZTop)
      deallocate(areaSum)

      limitCount(:) = 0
      !loop over cells and check for out of bounds
      if(isActiveTracer) then
        minimumVal(1) = -2.0_RKIND
        minimumVal(2:nTracers) = 0.0_RKIND
        rediLimiterCount = 0
      else
        minimumVal(:) = 0.0_RKIND
      endif

      do iCell = 1, nCells
        do k = 1, maxLevelCell(iCell)
           do iTr=1,ntracers
              tempTracer = tracers(iTr,k,iCell) + dt*(redi_term1(iTr,k,iCell) + redi_term2(iTr,k,iCell) + &
                              redi_term3(iTr,k,iCell))
              if(tempTracer < minimumVal(iTr)) then

                  do i = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(i, iCell)
                     redi_term2_edge(iTr,k,iEdge) = 0.0_RKIND
                  enddo

                  redi_term3_topOfCell(iTr,k,iCell) = 0.0_RKIND
                  redi_term3_topOfCell(iTr,k+1,iCell) = 0.0_RKIND
                  limitCount(iTr) = limitCount(iTr) + 1
              endif
           enddo
        enddo
     enddo

     if(isActiveTracer) then
       rediLimiterCount = limitCount(1)
     endif

     !now go back and reapply all tendencies
     do iCell=1,nCells
        do k = 1, maxLevelCell(iCell)
          do iTr = 1,ntracers
             tend(iTr,k,iCell) = tend(iTr,k,iCell) + redi_term1(iTr,k,iCell)
             tend(iTr,k,iCell) = tend(iTr,k,iCell) + rediKappa(iCell)*2.0_RKIND* &
                      (redi_term3_topOfCell(iTr,k,iCell) - redi_term3_topOfCell(iTr,k+1,iCell))
          enddo
        enddo

        do i = 1, nEdgesOnCell(iCell)
           iEdge = edgesOnCell(i, iCell)
           do k = 1,maxLevelEdgeTop(iEdge)
             do iTr = 1,ntracers
              tend(iTr,k,iCell) = tend(iTr,k,iCell) + redi_term2_edge(iTr,k,iEdge)
            enddo
           enddo
        enddo
     enddo


      deallocate(minimumVal)
      deallocate(redi_term2_edge)
      deallocate(redi_term3_topOfCell)

      call mpas_timer_stop("tracer redi")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_Redi_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hmix_Redi_init
!
!> \brief   Initializes ocean tracer horizontal mixing quantities
!> \details
!>  This routine initializes a variety of quantities related to
!>  Laplacian horizontal velocity mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_Redi_init(domain, err)!{{{

      type (domain_type), intent(in) :: domain

      integer, intent(out) :: err !< Output: error flag

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: forcingPool

      real (kind=RKIND), dimension(:), pointer :: RediKappa, RediKappaData
      real (kind=RKIND), dimension(:), pointer :: areaCell

      integer :: k, iCell
      integer, pointer :: nVertLevels, nCells
      real (kind=RKIND) :: sqrtPiInv2, avgCellDiameter

      err = 0

      if (.not.config_use_Redi) return

      if (config_redi_N2_limit_term1) then
        term1TaperFactor = 1.0_RKIND
      else
        term1TaperFactor = 0.0_RKIND
      endif

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(diagnosticsPool, 'RediKappa', RediKappa)
         call mpas_pool_get_array(forcingPool, 'RediKappaData', RediKappaData)

         ! initialize Redi kappa array
         if (config_Redi_closure == 'constant') then
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               RediKappa(iCell) = config_Redi_kappa
            end do
            !$omp end do
         else if (config_Redi_closure == 'data') then
             ! read RediKappa in from input
            call mpas_log_write( "config_Redi_closure = 'data'. "// &
               "Make sure that the variable RediKappa is read in from an input file.")
         else
            call mpas_log_write( 'Invalid choice of config_Redi_closure.', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain % dminfo)
         end if

         ! Add resolution taper
         if (config_eddying_resolution_taper == 'none') then
            ! Nothing to do, as we just keep the same assignment as above.
         else if (config_eddying_resolution_taper == 'ramp') then
            sqrtPiInv2 = 2.0_RKIND / sqrt(pii)
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               avgCellDiameter = sqrtPiInv2*sqrt(areaCell(iCell))
               if (avgCellDiameter <= config_eddying_resolution_ramp_min) then
                  RediKappa(iCell) = 0.0_RKIND
               else if (avgCellDiameter >= config_eddying_resolution_ramp_max) then
                  ! nothing to do, i.e. RediKappa(iCell) = RediKappa(iCell)
               else
                  RediKappa(iCell) = RediKappa(iCell) &
                     *(avgCellDiameter - config_eddying_resolution_ramp_min) & 
                     /(config_eddying_resolution_ramp_max - config_eddying_resolution_ramp_min)
               end if
            end do
            !$omp end do
         else
            call mpas_log_write( 'Invalid choice of config_eddying_resolution_taper.', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain % dminfo)
         end if
         block => block % next
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_Redi_init!}}}

!***********************************************************************

end module ocn_tracer_hmix_Redi

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
