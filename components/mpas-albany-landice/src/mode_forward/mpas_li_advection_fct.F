! and the University Corporation for Atmospheric Research (UCAR).
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_tracer_advection_fct
!
!> \brief MPAS monotonic tracer advection with FCT
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones, Trevor Hillebrand
!> \date   October 2017, updated May 2019, ported to MALI Sept 2023
!> \details
!>  This module contains routines for monotonic advection of tracers
!>  using a Flux Corrected Transport (FCT) algorithm
!
!-----------------------------------------------------------------------

module li_tracer_advection_fct

   ! module includes
#ifdef _ADV_TIMERS
   use mpas_timer
#endif
   use mpas_kind_types

   use li_config
   use li_mesh
   use li_tracer_advection_fct_shared

   implicit none
   private
   save

   ! module private variables
   real (kind=RKIND) ::  &
      coef3rdOrder        !< high-order horizontal coefficient
   logical ::            &
      monotonicityCheck   !< flag to check monotonicity

   ! public method interfaces
   public :: li_tracer_advection_fct_tend, &
             li_tracer_advection_fct_init

!**********************************************************************

   contains

!**********************************************************************
!
!  routine li_tracer_advection_fct_tend
!
!> \brief MPAS monotonic tracer horizontal advection tendency with FCT
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones, Trevor Hillebrand
!> \date   October 2017, updated May 2019, ported to MALI Sept 2023
!> \details
!>  This routine computes the monotonic tracer horizontal advection
!>  tendency using a flux-corrected transport (FCT) algorithm.
!
!-----------------------------------------------------------------------

   subroutine li_tracer_advection_fct_tend( &
                                        tend, tracers, layerThickness, &
                                        normalThicknessFlux, w, dt,    &
                                        nTracers, &
                                        layerThicknessEdgeFlux, &
                                        computeBudgets)!{{{
      use li_mesh
      !-----------------------------------------------------------------
      ! Input/Output parameters
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend    !< [inout] Tracer tendency to which advection added
      real (kind=RKIND), dimension(:,:), intent(inout), optional :: &
         layerThicknessEdgeFlux  !< [inout] used to compute higher order normalThicknessFlux
      !-----------------------------------------------------------------
      ! Input parameters
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracers               !< [in] Current tracer values

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness,      &!< [in] Thickness
         normalThicknessFlux, &!< [in] Thichness weighted velocitiy
         w                     !< [in] Vertical velocity

      real (kind=RKIND), intent(in) :: &
         dt                    !< [in] Timestep

      integer, intent(in) :: nTracers

      logical, intent(in) :: &
         computeBudgets    !< [in] Flag to compute active tracer budgets
      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer ::          &
         i, iCell, iEdge, &! horz indices
         cell1, cell2,    &! neighbor cell indices
         k,kmin, kmax, &! vert index variants
         kmin1, kmax1,    &! vert index variants
         iTracer           ! tracer index

      real (kind=RKIND) ::  &
         signedFactor,      &! temp factor including flux sign
         tracerMinNew,      &! updated tracer minimum
         tracerMaxNew,      &! updated tracer maximum
         tracerUpwindNew,   &! tracer updated with upwind flx
         scaleFactor,       &! factor for normalizing fluxes
         flux,              &! flux temporary
         tracerWeight,      &! tracer weighting temporary
         invAreaCell1,      &! inverse cell area
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

      ! allocate temporary arrays
      allocate(wgtTmp      (nVertLevels), &
               flxTmp      (nVertLevels), &
               sgnTmp      (nVertLevels), &
               tracerCur   (nVertLevels  ,nCells+1), &
               tracerMin   (nVertLevels  ,nCells), &
               tracerMax   (nVertLevels  ,nCells), &
               hNewInv     (nVertLevels  ,nCells), &
               hProv       (nVertLevels  ,nCells), &
               hProvInv    (nVertLevels  ,nCells), &
               flxIn       (nVertLevels  ,nCells+1), &
               flxOut      (nVertLevels  ,nCells+1), &
               workTend    (nVertLevels  ,nCells+1), &
               lowOrderFlx (nVertLevels+1,max(nCells,nEdges)+1), &
               highOrderFlx(nVertLevels+1,max(nCells,nEdges)+1))

      ! Initialize variables so you don't get floating point exceptions
      wgtTmp(:) = 0.0_RKIND
      flxTmp(:) = 0.0_RKIND
      sgnTmp(:) = 0.0_RKIND
      tracerCur(:,:) = 0.0_RKIND
      tracerMin(:,:) = 0.0_RKIND
      tracerMax(:,:) = 0.0_RKIND
      hNewInv(:,:) = 0.0_RKIND
      hProv(:,:) = 0.0_RKIND
      hProvInv(:,:) = 0.0_RKIND
      flxIn(:,:) = 0.0_RKIND
      flxOut(:,:) = 0.0_RKIND
      workTend(:,:) = 0.0_RKIND
      lowOrderFlx(:,:) = 0.0_RKIND
      highOrderFlx(:,:) = 0.0_RKIND

#ifdef _ADV_TIMERS
      call mpas_timer_stop('allocates')
      call mpas_timer_start('prov thickness')
#endif

      ! Compute some provisional layer thicknesses
      ! Note: This assumes we are in the first part of the horizontal/
      ! vertical operator splitting, which is true because currently
      ! we dont flip order and horizontal is always first.
      ! See notes in commit 2cd4a89d.

      do iCell = 1, nCells
         invAreaCell1 = dt/areaCell(iCell)
         kmin = 1 !minLevelCell(iCell)
         kmax = nVertLevels !maxLevelCell(iCell)
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
            if (hProv(k, iCell) > 0.0_RKIND) then
               hProvInv(k,iCell) = 1.0_RKIND/ hProv(k,iCell)
               hNewInv (k,iCell) = 1.0_RKIND/(hProv(k,iCell) - &
                                   dt*w(k,iCell) + dt*w(k+1, iCell))
            else
               hProvInv(k, iCell) = 0.0_RKIND
               hNewInv(k, iCell)  = 0.0_RKIND
            endif
         end do
      end do

#ifdef _ADV_TIMERS
      call mpas_timer_stop('prov thickness')
#endif

      ! Loop over tracers. One tracer is advected at a time.
      do iTracer = 1, nTracers

#ifdef _ADV_TIMERS
        call mpas_timer_start('cell init')
#endif

        ! Extract current tracer and change index order to improve locality
        do iCell = 1, nCells+1
        do k=1, nVertLevels
           tracerCur(k,iCell) = tracers(iTracer,k,iCell)
        end do ! k loop
        end do ! iCell loop

        ! Compute the high and low order horizontal fluxes.
#ifdef _ADV_TIMERS
        call mpas_timer_stop('cell init')
        call mpas_timer_start('tracer bounds')
#endif

        ! set nCells to first halo level
        ! nCells = nCellsHalo(1)

        ! Determine bounds on tracer (tracerMin and tracerMax) from
        ! surrounding cells for later limiting.

        do iCell = 1, nCellsHalo(1)
           kmin = 1 !minLevelCell(iCell)
           kmax = nVertLevels !maxLevelCell(iCell)
           do k=kmin,kmax
              tracerMin(k,iCell) = tracerCur(k,iCell)
              tracerMax(k,iCell) = tracerCur(k,iCell)
           end do
           ! TODO: Determine how/if this translates to MALI
           do i = 1, nEdgesOnCell(iCell)
              cell2 = cellsOnCell(i,iCell)
              !kmin1 = max(kmin, minLevelCell(cell2))
              !kmax1 = min(kmax, maxLevelCell(cell2))
              do k=kmin, kmax !kmin1,kmax1
                 tracerMax(k,iCell) = max(tracerMax(k,iCell), &
                                          tracerCur(k,cell2))
                 tracerMin(k,iCell) = min(tracerMin(k,iCell), &
                                          tracerCur(k,cell2))
              end do ! k loop
           end do ! i loop over nEdgesOnCell
        end do ! cell loop

#ifdef _ADV_TIMERS
        call mpas_timer_stop('tracer bounds')
        call mpas_timer_start('horiz flux')
#endif
        ! Need all the edges around the 1 halo cells and owned cells
        ! nEdges = nEdgesHalo(2)

        ! Compute the high order horizontal flux
        do iEdge = 1, nEdgesHalo(2)
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
              do k = 1, nVertLevels !minLevelCell(iCell), maxLevelCell(iCell)
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
           do k = 1, nVertLevels !minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
              tracerWeight = advMask2ndOrder(k,iEdge) &
                           * (1.0_RKIND - advMaskHighOrder(k, iEdge)) &
                           * (dvEdge(iEdge) * 0.5_RKIND)             &
                           * normalThicknessFlux(k, iEdge)

              lowOrderFlx(k,iEdge) = dvEdge(iEdge) * &
               (max(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracerCur(k,cell1) &
              + min(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracerCur(k,cell2))

              highOrderFlx(k,iEdge) = highOrderFlx(k,iedge) &
                                    + tracerWeight*(tracerCur(k,cell1) &
                                                  + tracerCur(k,cell2))

              ! only remove low order flux where high order flux is valid to
              ! avoid introducing nonzero values to highOrderFlx where it is invalid
              highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) &
                                    -  lowOrderFlx(k,iEdge) * &
                                    (advMaskHighOrder(k,iEdge) + advMask2ndOrder(k,iEdge) &
                                     * (1.0_RKIND - advMaskHighOrder(k, iEdge)) )
           end do ! k loop
        end do ! iEdge loop

#ifdef _ADV_TIMERS
        call mpas_timer_stop('horiz flux')
        call mpas_timer_start('scale factor build')
#endif

        ! Initialize flux arrays for all cells
        do iCell = 1, nCells+1
        do k=1, nVertLevels
           workTend(k, iCell) = 0.0_RKIND
           flxIn   (k, iCell) = 0.0_RKIND
           flxOut  (k, iCell) = 0.0_RKIND
        end do ! k loop
        end do ! iCell loop

        ! Need one halo of cells around owned cells
        ! nCells = nCellsHalo(1)

        do iCell = 1, nCellsHalo(1)
           invAreaCell1 = 1.0_RKIND / areaCell(iCell)

           ! Finish computing the low order horizontal fluxes
           ! Upwind fluxes are accumulated in workTend
           do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              cell1 = cellsOnEdge(1,iEdge)
              cell2 = cellsOnEdge(2,iEdge)
              signedFactor = edgeSignOnCell(i, iCell) * invAreaCell1

              do k = 1, nVertLevels !minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)

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
           do k = 1, nVertLevels !minLevelCell(iCell), maxLevelCell(iCell)
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

#ifdef _ADV_TIMERS
        call mpas_timer_stop('scale factor build')
        call mpas_timer_start('rescale horiz fluxes')
#endif
        ! Need all of the edges around owned cells
        ! nEdges = nEdgesHalo(1)
        !  rescale the high order horizontal fluxes
        do iEdge = 1, nEdgesHalo(1)
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           do k = 1, nVertLevels !minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
              highOrderFlx(k,iEdge) = max(0.0_RKIND,highOrderFlx(k,iEdge))* &
                                      min(flxOut(k,cell1), flxIn (k,cell2)) &
                                    + min(0.0_RKIND,highOrderFlx(k,iEdge))* &
                                      min(flxIn (k,cell1), flxOut(k,cell2))
           end do ! k loop
        end do ! iEdge loop

#ifdef _ADV_TIMERS
        call mpas_timer_stop('rescale horiz fluxes')
        call mpas_timer_start('flux accumulate')
#endif

        ! Accumulate the scaled high order vertical tendencies
        ! and the upwind tendencies
        do iCell = 1, nCellsSolve
           invAreaCell1 = 1.0_RKIND / areaCell(iCell)

           ! Accumulate the scaled high order horizontal tendencies
           do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              signedFactor = invAreaCell1*edgeSignOnCell(i,iCell)
              do k = 1, nVertLevels !minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                 ! workTend on RHS is upwind tendency
                 ! workTend on LHS is total horiz advect tendency
                 workTend(k,iCell) = workTend(k,iCell) &
                                   + signedFactor*highOrderFlx(k,iEdge)
              end do
           end do

           do k = 1, nVertLevels !minLevelCell(iCell), maxLevelCell(iCell)
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

#ifdef _ADV_TIMERS
        call mpas_timer_stop('flux accumulate')
        call mpas_timer_start('advect diags horiz')
#endif
        ! Compute budget and monotonicity diagnostics if needed

        ! Use layerThicknessEdgeFlux from the call to fct for thickness
        ! advection as the higher-order normalThicknessFlux for fct tracer advection.
        if (present(layerThicknessEdgeFlux)) then
            do iEdge = 1,nEdges
            do k = 1, nVertLevels
               ! Save u*h*T flux on edge for analysis. This variable will be
               ! divided by h at the end of the time step.
               ! average normal velocities from layer interfaces to layer midpoints
               if (dvEdge(iEdge) > 0.0_RKIND) then
                  layerThicknessEdgeFlux(k,iEdge) = &
                      (lowOrderFlx(k,iEdge) + highOrderFlx(k,iEdge))/dvEdge(iEdge)
               else
                  layerThicknessEdgeFlux(k,iEdge) = 0.0_RKIND
               endif
            enddo
            enddo
        endif

!           ! TODO: Determine if it is necessary to define activeTracerHorizontalAdvectionTendency
!           ! for MALI
!           !do iCell = 1, nCellsSolve
!           !do k = 1, nVertLevels !minLevelCell(iCell), maxLevelCell(iCell)
!           !   activeTracerHorizontalAdvectionTendency(iTracer,k,iCell) = &
!           !                                          workTend(k,iCell)
!           !end do
!           !end do ! iCell loop

        if (monotonicityCheck) then
           ! Check tracer values against local min,max to detect
           ! non-monotone values and write warning if found

           do iCell = 1, nCellsSolve
           do k = 1, nVertLevels !minLevelCell(iCell), maxLevelCell(iCell)
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
        end if ! monotonicity check
#ifdef _ADV_TIMERS
        call mpas_timer_stop('advect diags horiz')
#endif

!TODO: implement vertical advection?
        ! Update tracer array
        tracers(iTracer,:,:) = tracerCur(:,:)
      end do ! iTracer loop

#ifdef _ADV_TIMERS
      call mpas_timer_start('deallocates')
#endif

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

   end subroutine li_tracer_advection_fct_tend!}}}

!**********************************************************************
!
!  routine li_tracer_advection_fct_init
!
!> \brief MPAS initialize monotonic tracer advection tendency with FCT
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This routine initializes monotonic tracer advection quantities for
!>  the flux-corrected transport (FCT) algorithm.
!
!-----------------------------------------------------------------------

   subroutine li_tracer_advection_fct_init(err)!{{{

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
         coef3rdOrder = config_advection_coef_3rd_order
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

   end subroutine li_tracer_advection_fct_init!}}}

!***********************************************************************

end module li_tracer_advection_fct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
