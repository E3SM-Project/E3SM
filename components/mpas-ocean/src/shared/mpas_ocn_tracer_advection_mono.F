! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection_mono
!
!> \brief MPAS monotonic tracer advection with FCT
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This module contains routines for monotonic advection of tracers
!>  using a Flux Corrected Transport (FCT) algorithm
!
!-----------------------------------------------------------------------

module ocn_tracer_advection_mono

   ! module includes
#ifdef _ADV_TIMERS
   use mpas_timer
#endif
   use mpas_kind_types

   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables
   use ocn_tracer_advection_shared
   use ocn_tracer_advection_vert

   implicit none
   private
   save

   ! module private variables
   real (kind=RKIND) ::  &
      coef3rdOrder        !< high-order horizontal coefficient
   logical ::            &
      monotonicityCheck   !< flag to check monotonicity

   ! public method interfaces
   public :: ocn_tracer_advection_mono_tend, &
             ocn_tracer_advection_mono_init

!**********************************************************************

   contains

!**********************************************************************
!
!  routine ocn_tracer_advection_mono_tend
!
!> \brief MPAS monotonic tracer horizontal advection tendency with FCT
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This routine computes the monotonic tracer horizontal advection
!>  tendency using a flux-corrected transport (FCT) algorithm.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_advection_mono_tend( &
                                        tend, tracers, layerThickness, &
                                        normalThicknessFlux, w, dt,    &
                                        computeBudgets)!{{{

      !-----------------------------------------------------------------
      ! Input/Output parameters
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend    !< [inout] Tracer tendency to which advection added

      !-----------------------------------------------------------------
      ! Input parameters
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers               !< [in] Current tracer values

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness,      &!< [in] Thickness
         normalThicknessFlux, &!< [in] Thichness weighted velocitiy
         w                     !< [in] Vertical velocity

      real (kind=RKIND), intent(in) :: &
         dt                    !< [in] Timestep

      logical, intent(in) :: &
         computeBudgets        !< [in] Flag to compute active tracer budgets

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer ::          &
         i, iCell, iEdge, &! horz indices
         cell1, cell2,    &! neighbor cell indices
         nCells, nEdges,  &! numbers of cells or edges
         k,k2,kmin, kmax, &! vert index variants
         kmin1, kmax1,    &! vert index variants
         iTracer,         &! tracer index
         numTracers        ! total number of tracers

      real (kind=RKIND) ::  &
         signedFactor,      &! temp factor including flux sign
         tracerMinNew,      &! updated tracer minimum
         tracerMaxNew,      &! updated tracer maximum
         tracerUpwindNew,   &! tracer updated with upwind flx
         scaleFactor,       &! factor for normalizing fluxes
         flux,              &! flux temporary
         tracerWeight,      &! tracer weighting temporary
         invAreaCell1,      &! inverse cell area
         invAreaCell2,      &! inverse cell area
         coef1, coef3        ! temporary coefficients

      real (kind=RKIND), dimension(:), allocatable :: &
         wgtTmp,            &! vertical temporaries for
         flxTmp, sgnTmp      !   high-order flux computation

      real (kind=RKIND), dimension(:,:), allocatable :: &
         tracerCur,     &! reordered current tracer
         tracerMax,     &! max tracer in neighbors for limiting
         tracerMin,     &! min tracer in neighbors for limiting
         hNewInv,       &! inverse of new layer thickness
         hProv,         &! provisional layer thickness
         hProvInv,      &! inverse of provisional layer thickness
         flxIn,         &! flux coming into each cell
         flxOut,        &! flux going out of each cell
         workTend,      &! temp for holding some tendency values
         lowOrderFlx,   &! low order flux for FCT
         highOrderFlx    ! high order flux for FCT

      real (kind=RKIND), parameter :: &
         eps = 1.e-10_RKIND  ! small number to avoid numerical difficulties

      ! end of preamble
      !----------------
      ! begin code

#ifdef _ADV_TIMERS
      call mpas_timer_start('allocates')
#endif
      ! Get dimensions
      numTracers  = size(tracers,dim=1)

      ! allocate temporary arrays
      allocate(wgtTmp      (nVertLevels), &
               flxTmp      (nVertLevels), &
               sgnTmp      (nVertLevels), &
               tracerCur   (nVertLevels  ,nCellsAll+1), &
               tracerMin   (nVertLevels  ,nCellsAll), &
               tracerMax   (nVertLevels  ,nCellsAll), &
               hNewInv     (nVertLevels  ,nCellsAll), &
               hProv       (nVertLevels  ,nCellsAll), &
               hProvInv    (nVertLevels  ,nCellsAll), &
               flxIn       (nVertLevels  ,nCellsAll+1), &
               flxOut      (nVertLevels  ,nCellsAll+1), &
               workTend    (nVertLevels  ,nCellsAll+1), &
               lowOrderFlx (nVertLevels+1,max(nCellsAll,nEdgesAll)+1), &
               highOrderFlx(nVertLevels+1,max(nCellsAll,nEdgesAll)+1))

      !$acc enter data create(wgtTmp, flxTmp, sgnTmp, &
      !$acc                   tracerCur, tracerMin, tracerMax, &
      !$acc                   hNewInv, hProv, hProvInv, flxIn, flxOut, &
      !$acc                   workTend, lowOrderFlx, highOrderFlx)

#ifdef _ADV_TIMERS
      call mpas_timer_stop('allocates')
      call mpas_timer_start('prov thickness')
#endif

      ! Compute some provisional layer thicknesses
      ! Note: This assumes we are in the first part of the horizontal/
      ! vertical operator splitting, which is true because currently
      ! we dont flip order and horizontal is always first.
      ! See notes in commit 2cd4a89d.

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(areaCell, dvEdge, minLevelCell, maxLevelCell, &
      !$acc            nEdgesOnCell, edgesOnCell, edgeSignOnCell, &
      !$acc            layerThickness, normalThicknessFlux, w, &
      !$acc            hProv, hProvInv, hNewInv) &
      !$acc    private(k,kmin,kmax,i,iEdge, invAreaCell1, signedFactor) 
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(k,kmin,kmax,i,iEdge, invAreaCell1, signedFactor) 
#endif
      do iCell = 1, nCellsAll
         invAreaCell1 = dt/areaCell(iCell)
         kmin = minLevelCell(iCell)
         kmax = maxLevelCell(iCell)
         do k = kmin, kmax
            hProv(k, iCell) = layerThickness(k, iCell)
         end do
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            signedFactor = invAreaCell1*dvEdge(iEdge)* &
                           edgeSignOnCell(i,iCell)
            ! Provisional layer thickness is after horizontal
            ! thickness flux only
            do k = kmin, kmax
               hProv(k,iCell) = hProv(k,iCell) &
                              + signedFactor*normalThicknessFlux(k,iEdge)
            end do
         end do
         ! New layer thickness is after horizontal and vertical
         ! thickness flux
         do k = kmin, kmax
            hProvInv(k,iCell) = 1.0_RKIND/ hProv(k,iCell)
            hNewInv (k,iCell) = 1.0_RKIND/(hProv(k,iCell) - &
                                dt*w(k,iCell) + dt*w(k+1, iCell))
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef _ADV_TIMERS
      call mpas_timer_stop('prov thickness')
#endif

      ! Loop over tracers. One tracer is advected at a time.
      do iTracer = 1, numTracers

#ifdef _ADV_TIMERS
        call mpas_timer_start('cell init')
#endif

        ! Extract current tracer and change index order to improve locality
#ifdef MPAS_OPENACC
        !$acc parallel loop collapse(2) &
        !$acc    present(tracerCur, tracers)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(k)
#endif
        do iCell = 1, nCellsAll+1
        do k=1, nVertLevels
           tracerCur(k,iCell) = tracers(iTracer,k,iCell)
        end do ! k loop
        end do ! iCell loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

        ! Compute the high and low order horizontal fluxes.
#ifdef _ADV_TIMERS
        call mpas_timer_stop('cell init')
        call mpas_timer_start('tracer bounds')
#endif

        ! set nCells to first halo level
        nCells = nCellsHalo(1)

        ! Determine bounds on tracer (tracerMin and tracerMax) from
        ! surrounding cells for later limiting.

#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(minLevelCell, maxLevelCell, nEdgesOnCell, &
        !$acc            cellsOnCell, tracerCur, tracerMin, tracerMax) &
        !$acc    private(k, kmin, kmax, kmin1, kmax1, i, cell2)
#else
        !$omp parallel
        !$omp do schedule(runtime) &
        !$omp    private(k, kmin, kmax, kmin1, kmax1, i, cell2)
#endif
        do iCell = 1, nCells
           kmin = minLevelCell(iCell)
           kmax = maxLevelCell(iCell)
           do k=kmin,kmax
              tracerMin(k,iCell) = tracerCur(k,iCell)
              tracerMax(k,iCell) = tracerCur(k,iCell)
           end do
           do i = 1, nEdgesOnCell(iCell)
              cell2 = cellsOnCell(i,iCell)
              kmin1 = max(kmin, minLevelCell(cell2))
              kmax1 = min(kmax, maxLevelCell(cell2))
              do k=kmin1,kmax1
                 tracerMax(k,iCell) = max(tracerMax(k,iCell), &
                                          tracerCur(k,cell2))
                 tracerMin(k,iCell) = min(tracerMin(k,iCell), &
                                          tracerCur(k,cell2))
              end do ! k loop
           end do ! i loop over nEdgesOnCell
        end do ! cell loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_stop('tracer bounds')
        call mpas_timer_start('horiz flux')
#endif
        ! Need all the edges around the 1 halo cells and owned cells
        nEdges = nEdgesHalo(2)

        ! Compute the high order horizontal flux

#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(minLevelCell, maxLevelCell, minLevelEdgeBot, &
        !$acc            maxLevelEdgeTop, nAdvCellsForEdge, & 
        !$acc            advCellsForEdge, cellsOnEdge, dvEdge, &
        !$acc            advMaskHighOrder, advCoefs, advCoefs3rd, &
        !$acc            tracerCur, normalThicknessFlux, &
        !$acc            highOrderFlx, lowOrderFlx) &
        !$acc    private(i, k, icell, cell1, cell2, coef1, coef3, &
        !$acc            wgtTmp, sgnTmp, flxTmp, tracerWeight)
#else
        !$omp parallel
        !$omp do schedule(runtime) &
        !$omp    private(i, k, icell, cell1, cell2, coef1, coef3, &
        !$omp            wgtTmp, sgnTmp, flxTmp, tracerWeight)
#endif
        do iEdge = 1, nEdges
           cell1 = cellsOnEdge(1, iEdge)
           cell2 = cellsOnEdge(2, iEdge)

           ! compute some common intermediate factors
           do k = 1, nVertLevels
              wgtTmp(k) = normalThicknessFlux(k,iEdge)* &
                          advMaskHighOrder(k,iEdge)
              sgnTmp(k) = sign(1.0_RKIND, &
                               normalThicknessFlux(k,iEdge))
              flxTmp(k) = 0.0_RKIND
           end do

           ! Compute 3rd or 4th fluxes where requested.
           do i = 1, nAdvCellsForEdge(iEdge)
              iCell = advCellsForEdge(i,iEdge)
              coef1 = advCoefs       (i,iEdge)
              coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder
              do k = minLevelCell(iCell), maxLevelCell(iCell)
                 flxTmp(k) = flxTmp(k) + tracerCur(k,iCell)* &
                             wgtTmp(k)*(coef1 + coef3*sgnTmp(k))
              end do ! k loop
           end do ! i loop over nAdvCellsForEdge

           do k=1,nVertLevels
              highOrderFlx(k,iEdge) = flxTmp(k)
           end do

           ! Compute 2nd order fluxes where needed.
           ! Also compute low order upwind horizontal flux (monotonic)
           ! Remove low order flux from the high order flux
           ! Store left over high order flux in highOrderFlx array
           do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
              tracerWeight = (1.0_RKIND - advMaskHighOrder(k,iEdge)) &
                           * (dvEdge(iEdge) * 0.5_RKIND)             &
                           * normalThicknessFlux(k, iEdge)

              lowOrderFlx(k,iEdge) = dvEdge(iEdge) * &
               (max(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracerCur(k,cell1) &
              + min(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracerCur(k,cell2))

              highOrderFlx(k,iEdge) = highOrderFlx(k,iedge) &
                                    + tracerWeight*(tracerCur(k,cell1) &
                                                  + tracerCur(k,cell2))

              highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) &
                                    -  lowOrderFlx(k,iEdge)
           end do ! k loop
        end do ! iEdge loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_stop('horiz flux')
        call mpas_timer_start('scale factor build')
#endif

        ! Initialize flux arrays for all cells

#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present (workTend, flxIn, flxOut)
#else
        !$omp parallel
        !$omp do schedule(runtime)
#endif
        do iCell = 1, nCellsAll+1
        do k=1, nVertLevels
           workTend(k, iCell) = 0.0_RKIND
           flxIn   (k, iCell) = 0.0_RKIND
           flxOut  (k, iCell) = 0.0_RKIND
        end do ! k loop
        end do ! iCell loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

        ! Need one halo of cells around owned cells
        nCells = nCellsHalo(1)

#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(minLevelEdgeBot, maxLevelEdgeTop, &
        !$acc            minLevelCell, maxLevelCell, areaCell, &
        !$acc            nEdgesOnCell, edgesOnCell, cellsOnEdge, &
        !$acc            edgeSignOnCell, layerThickness, hProvInv, &
        !$acc            tracerCur, tracerMax, tracerMin, workTend, &
        !$acc            flxOut, flxIn, lowOrderFlx, highOrderFlx) &
        !$acc    private(i, k, iEdge, cell1, cell2, &
        !$acc            invAreaCell1, signedFactor, scaleFactor, &
        !$acc            tracerUpwindNew, tracerMinNew, tracerMaxNew)
#else
        !$omp parallel
        !$omp do schedule(runtime) &
        !$omp    private(i, k, iEdge, cell1, cell2, &
        !$omp            invAreaCell1, signedFactor, scaleFactor, &
        !$omp            tracerUpwindNew, tracerMinNew, tracerMaxNew)
#endif
        do iCell = 1, nCells
           invAreaCell1 = 1.0_RKIND / areaCell(iCell)

           ! Finish computing the low order horizontal fluxes
           ! Upwind fluxes are accumulated in workTend
           do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              cell1 = cellsOnEdge(1,iEdge)
              cell2 = cellsOnEdge(2,iEdge)
              signedFactor = edgeSignOnCell(i, iCell) * invAreaCell1

              do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)

                 ! Here workTend is the advection tendency due to the
                 ! upwind (low order) fluxes.
                 workTend(k,iCell) = workTend(k,iCell) &
                                   + signedFactor*lowOrderFlx(k,iEdge)

                 ! Accumulate remaining high order fluxes
                 flxOut(k,iCell) = flxOut(k,iCell) + min(0.0_RKIND,  &
                                   signedFactor*highOrderFlx(k,iEdge))
                 flxIn (k,iCell) = flxIn (k,iCell) + max(0.0_RKIND,  &
                                   signedFactor*highOrderFlx(k,iEdge))

              end do
           end do

           ! Build the factors for the FCT
           ! Computed using the bounds that were computed previously,
           ! and the bounds on the newly updated value
           ! Factors are placed in the flxIn and flxOut arrays
           do k = minLevelCell(iCell), maxLevelCell(iCell)
              ! Here workTend is the upwind tendency
              tracerUpwindNew = (tracerCur(k,iCell)*layerThickness(k,iCell) &
                              + dt*workTend(k,iCell))*hProvInv(k,iCell)
              tracerMinNew = tracerUpwindNew &
                           + dt*flxOut(k,iCell)*hProvInv(k,iCell)
              tracerMaxNew = tracerUpwindNew &
                           + dt*flxIn (k,iCell)*hProvInv(k,iCell)

              scaleFactor = (tracerMax(k,iCell) - tracerUpwindNew)/ &
                            (tracerMaxNew - tracerUpwindNew + eps)
              flxIn (k,iCell) = min(1.0_RKIND, &
                                max(0.0_RKIND, scaleFactor))
              scaleFactor = (tracerUpwindNew - tracerMin(k,iCell))/ &
                            (tracerUpwindNew - tracerMinNew + eps)
              flxOut(k,iCell) = min(1.0_RKIND, &
                                max(0.0_RKIND, scaleFactor))
           end do ! k loop
        end do ! iCell loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_stop('scale factor build')
        call mpas_timer_start('rescale horiz fluxes')
#endif
        ! Need all of the edges around owned cells
        nEdges = nEdgesHalo(1)
        !  rescale the high order horizontal fluxes

#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(minLevelEdgeBot, maxLevelEdgeTop, &
        !$acc            cellsOnEdge, highOrderFlx, flxIn, flxOut) &
        !$acc    private(k, cell1, cell2)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(k, cell1, cell2)
#endif
        do iEdge = 1, nEdges
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
              highOrderFlx(k,iEdge) = max(0.0_RKIND,highOrderFlx(k,iEdge))* &
                                      min(flxOut(k,cell1), flxIn (k,cell2)) &
                                    + min(0.0_RKIND,highOrderFlx(k,iEdge))* &
                                      min(flxIn (k,cell1), flxOut(k,cell2))
           end do ! k loop
        end do ! iEdge loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_stop('rescale horiz fluxes')
        call mpas_timer_start('flux accumulate')
#endif

        ! Accumulate the scaled high order vertical tendencies
        ! and the upwind tendencies
#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(minLevelEdgeBot, maxLevelEdgeTop, &
        !$acc            minLevelCell, maxLevelCell, areaCell, &
        !$acc            nEdgesOnCell, edgesOnCell, edgeSignOnCell, &
        !$acc            workTend, highOrderFlx, layerThickness, &
        !$acc            tracerCur, hProvInv, tend) &
        !$acc    private(i, k, iEdge, invAreaCell1, signedFactor) 
#else
        !$omp parallel
        !$omp do schedule(runtime) &
        !$omp    private(i, k, iEdge, invAreaCell1, signedFactor) 
#endif
        do iCell = 1, nCellsOwned
           invAreaCell1 = 1.0_RKIND / areaCell(iCell)

           ! Accumulate the scaled high order horizontal tendencies
           do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              signedFactor = invAreaCell1*edgeSignOnCell(i,iCell)
              do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                 ! workTend on RHS is upwind tendency
                 ! workTend on LHS is total horiz advect tendency
                 workTend(k,iCell) = workTend(k,iCell) &
                                   + signedFactor*highOrderFlx(k,iEdge)
              end do
           end do

           do k = minLevelCell(iCell), maxLevelCell(iCell)
              ! workTend  on RHS is total horiz advection tendency
              ! tracerCur on LHS is provisional tracer after
              !                     horizontal fluxes only.
              tracerCur(k,iCell) = (tracerCur(k,iCell)* &
                                    layerThickness(k,iCell) &
                                    + dt*workTend(k,iCell)) &
                                 * hProvInv(k,iCell)
              tend(iTracer,k,iCell) = tend(iTracer,k,iCell) &
                                    + workTend(k,iCell)
           end do

        end do ! iCell loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_stop('flux accumulate')
        call mpas_timer_start('advect diags horiz')
#endif
        ! Compute budget and monotonicity diagnostics if needed
        if (computeBudgets) then

           nEdges = nEdgesHalo(2)
#ifdef MPAS_OPENACC
           !$acc parallel loop private(k) &
           !$acc    present(activeTracerHorizontalAdvectionEdgeFlux, &
           !$acc            minLevelEdgeBot, maxLevelEdgeTop, &
           !$acc            lowOrderFlx, highOrderFlx, dvEdge)
#else
           !$omp parallel
           !$omp do schedule(runtime) private(k)
#endif
           do iEdge = 1,nEdges
           do k = minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
              ! Save u*h*T flux on edge for analysis. This variable will be
              ! divided by h at the end of the time step.
              activeTracerHorizontalAdvectionEdgeFlux(iTracer,k,iEdge) = &
                 (lowOrderFlx(k,iEdge) + highOrderFlx(k,iEdge))/dvEdge(iEdge)
           enddo
           enddo
#ifndef MPAS_OPENACC
           !$omp end do
#endif

#ifdef MPAS_OPENACC
           !$acc parallel loop private(k) &
           !$acc    present(activeTracerHorizontalAdvectionTendency, &
           !$acc            minLevelCell, maxLevelCell, workTend)
#else
           !$omp do schedule(runtime) private(k)
#endif
           do iCell = 1, nCellsOwned
           do k = minLevelCell(iCell), maxLevelCell(iCell)
              activeTracerHorizontalAdvectionTendency(iTracer,k,iCell) = &
                                                     workTend(k,iCell)
           end do
           end do ! iCell loop
#ifndef MPAS_OPENACC
           !$omp end do
           !$omp end parallel
#endif

        end if ! computeBudgets

        if (monotonicityCheck) then
           ! Check tracer values against local min,max to detect
           ! non-monotone values and write warning if found

           ! Perform check on host since print involved
           !$acc update host(tracerCur, tracerMin, tracerMax)

           !$omp parallel
           !$omp do schedule(runtime) private(k)
           do iCell = 1, nCellsOwned
           do k = minLevelCell(iCell), maxLevelCell(iCell)
              if(tracerCur(k,iCell) < tracerMin(k, iCell)-eps) then
                 call mpas_log_write( &
                    'Horizontal minimum out of bounds on tracer: $i $r $r ', &
                    MPAS_LOG_WARN, intArgs=(/iTracer/),                      &
                    realArgs=(/ tracerMin(k, iCell), tracerCur(k,iCell) /) )
              end if

              if(tracerCur(k,iCell) > tracerMax(k,iCell)+eps) then
                 call mpas_log_write( &
                    'Horizontal maximum out of bounds on tracer: $i $r $r ', &
                    MPAS_LOG_WARN, intArgs=(/iTracer/),                      &
                    realArgs=(/ tracerMax(k, iCell), tracerCur(k,iCell) /) )
              end if
           end do
           end do
           !$omp end do
           !$omp end parallel
        end if ! monotonicity check
#ifdef _ADV_TIMERS
        call mpas_timer_stop('advect diags horiz')
#endif

        !-----------------------------------------------------------------------
        !
        !  Horizontal advection complete
        !  Begin vertical advection
        !
        !-----------------------------------------------------------------------

#ifdef _ADV_TIMERS
        call mpas_timer_start('cell init vert')
#endif
        ! Initialize variables for use in this iTracer iteration

#ifdef MPAS_OPENACC
        !$acc parallel loop collapse(2) present(workTend)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(k)
#endif
        do iCell = 1, nCellsAll
        do k=1, nVertLevels
           workTend(k, iCell) = 0.0_RKIND
        end do ! k loop
        end do ! iCell loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_stop('cell init vert')
        call mpas_timer_start('vertical bounds')
#endif

        ! Need all owned and 1 halo cells
        nCells = nCellsHalo(1)

        ! Determine bounds on tracerCur from neighbor values for limiting
#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(tracerCur, tracerMin, tracerMax, &
        !$acc            minLevelCell, maxLevelCell) &
        !$acc    private(k, kmin, kmax)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(k, kmin, kmax)
#endif
        do iCell = 1, nCells
           kmin = minLevelCell(iCell)
           kmax = maxLevelCell(iCell)

           ! take care of top cell
           tracerMax(kmin,iCell) = max(tracerCur(1,iCell), &
                                       tracerCur(2,iCell))
           tracerMin(kmin,iCell) = min(tracerCur(1,iCell), &
                                       tracerCur(2,iCell))
           do k=kmin+1,kmax-1
              tracerMax(k,iCell) = max(tracerCur(k-1,iCell), &
                                       tracerCur(k  ,iCell), &
                                       tracerCur(k+1,iCell))
              tracerMin(k,iCell) = min(tracerCur(k-1,iCell), &
                                       tracerCur(k  ,iCell), &
                                       tracerCur(k+1,iCell))
           end do
           ! finish with bottom cell
           tracerMax(kmax,iCell) = max(tracerCur(kmax  ,iCell), &
                                       tracerCur(kmax-1,iCell))
           tracerMin(kmax,iCell) = min(tracerCur(kmax  ,iCell), &
                                       tracerCur(kmax-1,iCell))
        end do ! cell loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_stop('vertical bounds')
        call mpas_timer_start('vertical flux high')
#endif

        ! Compute the high order vertical fluxes based on order selected.

        call ocn_tracer_advection_vert_flx(tracerCur, w, hProv, &
                                           highOrderFlx)

#ifdef _ADV_TIMERS
        call mpas_timer_stop('vertical flux high')
        call mpas_timer_start('vertical flux')
#endif

        ! Compute low order upwind vertical flux (monotonic)
        ! Remove low order flux from the high order flux.
        ! Store left over high order flux in highOrderFlx array.

#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(minLevelCell, maxLevelCell, w, tracerCur, &
        !$acc            lowOrderFlx, highOrderFlx, workTend, &
        !$acc            flxIn, flxOut) &
        !$acc    private(k, kmin, kmax)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(k, kmin, kmax)
#endif
        do iCell = 1, nCells
           kmin = minLevelCell(iCell)
           kmax = maxLevelCell(iCell)

           lowOrderFlx(kmin,iCell) = 0.0_RKIND
           do k = kmin+1, kmax
              lowOrderFlx(k,iCell) = &
                    min(0.0_RKIND,w(k,iCell))*tracerCur(k-1,iCell) + &
                    max(0.0_RKIND,w(k,iCell))*tracerCur(k  ,iCell)
              highOrderFlx(k,iCell) = highOrderFlx(k,iCell) &
                                    -  lowOrderFlx(k,iCell)
           end do ! k loop
           lowOrderFlx(kmax+1,iCell) = 0.0_RKIND

           ! Upwind fluxes are accumulated in workTend
           ! flxIn  contains total remaining high order flux into iCell
           !          it is positive.
           ! flxOut contains total remaining high order flux out of iCell
           !          it is negative
           do k = kmin,kmax
              workTend(k,iCell) = lowOrderFlx(k+1,iCell) &
                                - lowOrderFlx(k  ,iCell)
              flxIn (k,iCell) = max(0.0_RKIND, highOrderFlx(k+1,iCell)) &
                              - min(0.0_RKIND, highOrderFlx(k  ,iCell))
              flxOut(k,iCell) = min(0.0_RKIND, highOrderFlx(k+1,iCell)) &
                              - max(0.0_RKIND, highOrderFlx(k  ,iCell))
           end do ! k Loop
        end do ! iCell Loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_stop('vertical flux')
        call mpas_timer_start('scale factor build')
#endif

        ! Build the scale factors to limit flux for FCT
        ! Computed using the bounds that were computed previously,
        ! and the bounds on the newly updated value
        ! Factors are placed in the flxIn and flxOut arrays

        nCells = nCellsHalo(1) ! Need one halo around owned cells

#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(minLevelCell, maxLevelCell, hProv, hNewInv, &
        !$acc            tracerCur, tracerMax, tracerMin, &
        !$acc            workTend, flxIn, flxOut) &
        !$acc    private(k, tracerMinNew, tracerMaxNew, &
        !$acc            tracerUpwindNew, scaleFactor) 
#else
        !$omp parallel
        !$omp do schedule(runtime) &
        !$omp    private(k, tracerMinNew, tracerMaxNew, &
        !$omp            tracerUpwindNew, scaleFactor)
#endif
        do iCell = 1, nCells

           do k = minLevelCell(iCell), maxLevelCell(iCell)
              ! workTend on the RHS is the upwind tendency
              tracerMinNew = (tracerCur(k,iCell)*hProv(k,iCell) &
                           + dt*(workTend(k,iCell)+flxOut(k,iCell))) &
                           * hNewInv(k,iCell)
              tracerMaxNew = (tracerCur(k,iCell)*hProv(k,iCell) &
                           + dt*(workTend(k,iCell)+flxIn(k,iCell))) &
                           * hNewInv(k,iCell)
              tracerUpwindNew = (tracerCur(k,iCell)*hProv(k,iCell) &
                              + dt*workTend(k,iCell)) * hNewInv(k,iCell)

              scaleFactor = (tracerMax(k,iCell)-tracerUpwindNew)/ &
                            (tracerMaxNew-tracerUpwindNew+eps)
              flxIn (k,iCell) = min(1.0_RKIND, &
                                max(0.0_RKIND, scaleFactor))

              scaleFactor = (tracerUpwindNew-tracerMin(k,iCell))/ &
                            (tracerUpwindNew-tracerMinNew+eps)
              flxOut(k,iCell) = min(1.0_RKIND, &
                                max(0.0_RKIND, scaleFactor))
           end do ! k loop
        end do ! iCell loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_stop('scale factor build')
        call mpas_timer_start('flux accumulate')
#endif

        ! Accumulate the scaled high order vertical tendencies
        ! and the upwind tendencies

#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(minLevelCell, maxLevelCell, highOrderFlx, &
        !$acc            flxIn, flxOut, workTend, tend) &
        !$acc    private(k, kmin, kmax, flux)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(k, kmin, kmax, flux)
#endif
        do iCell = 1, nCellsOwned
           kmin = minLevelCell(iCell)
           kmax = maxLevelCell(iCell)
           ! rescale the high order vertical flux
           do k = kmin+1, kmax
              flux =  highOrderFlx(k,iCell)
              highOrderFlx(k,iCell) = max(0.0_RKIND,flux)*   &
                                      min(flxOut(k  ,iCell), &
                                          flxIn (k-1,iCell)) &
                                    + min(0.0_RKIND,flux)*   &
                                      min(flxOut(k-1,iCell), &
                                          flxIn (k  ,iCell))
           end do ! k loop

           do k = kmin,kmax
              ! workTend on the RHS is upwind tendency
              ! workTend on the LHS is total vertical advection tendency
              workTend(k, iCell) = workTend(k, iCell)       &
                                 + (highOrderFlx(k+1,iCell) &
                                  - highOrderFlx(k  ,iCell))
              tend(iTracer,k,iCell) = tend(iTracer,k,iCell) &
                                    + workTend(k,iCell)
           end do ! k loop

        end do ! iCell loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

        ! Compute advection diagnostics and monotonicity checks if requested
#ifdef _ADV_TIMERS
        call mpas_timer_stop('flux accumulate')
        call mpas_timer_start('advect diags vert')
#endif
        if (computeBudgets) then

#ifdef MPAS_OPENACC
           !$acc parallel loop &
           !$acc    present(lowOrderFlx, highOrderFlx, workTend, &
           !$acc            activeTracerVerticalAdvectionTopFlux, &
           !$acc            activeTracerVerticalAdvectionTendency, &
           !$acc            minLevelCell, maxLevelCell) &
           !$acc    private(k, kmin, kmax)
#else
           !$omp parallel
           !$omp do schedule(runtime) private(k,kmin,kmax)
#endif
           do iCell = 1, nCellsOwned
              kmin = minLevelCell(iCell)
              kmax = maxLevelCell(iCell)
              do k = kmin+1,kmax
                 activeTracerVerticalAdvectionTopFlux(iTracer,k,iCell) = &
                    lowOrderFlx(k,iCell) + highOrderFlx(k,iCell)
              end do
              do k = kmin,kmax
                 activeTracerVerticalAdvectionTendency(iTracer,k,iCell) = &
                    workTend(k,iCell)
              end do
           end do ! iCell loop
#ifndef MPAS_OPENACC
           !$omp end do
           !$omp end parallel
#endif
        end if ! computeBudgets

        if (monotonicityCheck) then

           ! Check for monotonicity of new tracer value
           ! Use flxIn as a temp for new tracer value

#ifdef MPAS_OPENACC
           !$acc parallel loop collapse(2) &
           !$acc    present(flxIn, tracerCur, hProv, hNewInv, &
           !$acc            workTend, cellMask)
#else
           !$omp parallel
           !$omp do schedule(runtime) private(k)
#endif
           do iCell = 1,nCellsOwned
           do k = 1,nVertLevels
              ! workTend on the RHS is total vertical advection tendency
              flxIn(k,iCell) = (tracerCur(k, iCell)*hProv(k, iCell) &
                             + dt*workTend(k, iCell)) &
                             * hNewInv(k,iCell)*cellMask(k,iCell)
           end do ! vert
           end do ! cell loop
#ifndef MPAS_OPENACC
           !$omp end do
           !$omp end parallel
#endif

           ! Transfer new tracer values to host for diagnostic check
           !$acc update host(flxIn)

           ! Print out info for cells that are out of bounds
           ! flxIn holds new tracer values
           do iCell = 1, nCellsOwned
           do k = minLevelCell(iCell), maxLevelCell(iCell)
              if (flxIn(k,iCell) < tracerMin(k, iCell)-eps) then
                 call mpas_log_write( &
                   'Vertical minimum out of bounds on tracer: $i $i $i $r $r ',&
                   MPAS_LOG_WARN, intArgs=(/iTracer, k, iCell/), &
                   realArgs=(/ tracerMin(k,iCell), flxIn(k,iCell) /) )
              end if
              if (flxIn(k,iCell) > tracerMax(k,iCell)+eps) then
                 call mpas_log_write( &
                   'Vertical maximum out of bounds on tracer: $i $i $i $r $r ',&
                    MPAS_LOG_WARN, intArgs=(/iTracer, k, iCell/), &
                    realArgs=(/ tracerMax(k,iCell), flxIn(k,iCell) /) )
              end if
           end do
           end do

        end if ! monotonicity check
#ifdef _ADV_TIMERS
        call mpas_timer_stop('advect diags vert')
#endif
      end do ! iTracer loop

#ifdef _ADV_TIMERS
      call mpas_timer_start('deallocates')
#endif
      !$acc exit data delete(wgtTmp, flxTmp, sgnTmp, &
      !$acc                  tracerCur, tracerMin, tracerMax, &
      !$acc                  hNewInv, hProv, hProvInv, flxIn, flxOut, &
      !$acc                  workTend, lowOrderFlx, highOrderFlx)

      deallocate(wgtTmp,       &
                 flxTmp,       &
                 sgnTmp,       &
                 tracerCur,    &
                 tracerMin,    &
                 tracerMax,    &
                 hNewInv,      &
                 hProv,        &
                 hProvInv,     &
                 flxIn,        &
                 flxOut,       &
                 workTend,     &
                 lowOrderFlx,  &
                 highOrderFlx)

#ifdef _ADV_TIMERS
      call mpas_timer_stop('deallocates')
#endif

   end subroutine ocn_tracer_advection_mono_tend!}}}

!**********************************************************************
!
!  routine ocn_tracer_advection_mono_init
!
!> \brief MPAS initialize monotonic tracer advection tendency with FCT
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This routine initializes monotonic tracer advection quantities for
!>  the flux-corrected transport (FCT) algorithm.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_advection_mono_init(err)!{{{

      !*** output parameters

      integer, intent(out) :: &
         err                   !< [out] error flag

      ! end of preamble
      !----------------
      ! begin code

      err = 0 ! initialize error code to success

      ! Check that the halo is wide enough for FCT
      if (config_num_halos < 3) then
         call mpas_log_write( &
            'Monotonic advection cannot be used with less than 3 halos.', &
            MPAS_LOG_CRIT)
         err = -1
      end if

      ! Set blending coefficient if 3rd order horizontal advection chosen
      select case (config_horiz_tracer_adv_order)
      case (2)
         coef3rdOrder = 0.0_RKIND
      case (3)
         coef3rdOrder = config_coef_3rd_order
      case (4)
         coef3rdOrder = 0.0_RKIND
      case default
         coef3rdOrder = 0.0_RKIND
         call mpas_log_write( &
            'Invalid value for horizontal advection order, defaulting to 2',&
            MPAS_LOG_WARN)
      end select ! horizontal advection order

      ! Set flag for checking monotonicity
      monotonicityCheck = config_check_tracer_monotonicity

   end subroutine ocn_tracer_advection_mono_init!}}}

!***********************************************************************

end module ocn_tracer_advection_mono

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
