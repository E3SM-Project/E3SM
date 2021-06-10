! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection_std
!
!> \brief MPAS standard tracer advection
!> \author Doug Jacobsen, Phil Jones
!> \date   03/09/12, updated May 2019
!> \details
!>  This module contains routines for advection of tracers using a standard
!>  FV algorithm in MPAS discretization.
!
!-----------------------------------------------------------------------

module ocn_tracer_advection_std

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_log
   use ocn_mesh
   use ocn_tracer_advection_shared

   implicit none
   private
   save

   ! private module variables
   real (kind=RKIND) :: &
      coef3rdOrder       !< coefficient for blending high-order terms

   integer :: vertOrder  !< choice of order for vertical advection
   integer, parameter :: &! enumerator for supported vertical adv order
      vertOrder2=2,      &!< 2nd order
      vertOrder3=3,      &!< 3rd order
      vertOrder4=4        !< 4th order

   logical :: monotonicityCheck !< flag to check monotonicity

   ! public method interfaces
   public :: ocn_tracer_advection_std_tend, &
             ocn_tracer_advection_std_init

!***********************************************************************

   contains

!***********************************************************************
!
!  routine ocn_tracer_advection_std_tend
!
!> \brief MPAS standard tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine computes the standard tracer advection tendencity.
!>  Both horizontal and vertical.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_advection_std_tend(tracers, normalThicknessFlux, &
                                            w, layerThickness, dt, tend)

      !-----------------------------------------------------------------
      ! Input/Output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend            !< [inout] Accumulated tracer tendency

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers               !< [in] current tracer values

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalThicknessFlux, &!< [in] Thickness weighted velocitiy
         w,                   &!< [in] Vertical velocity
         layerThickness        !< [in] Thickness

      real (kind=RKIND), intent(in) :: &
         dt                    !< [in] Timestep

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: &
         i, k,         &! loop indices for neighbors, vertical
         kmin, kmax,   &! min, max active vertical layers
         iCell, iEdge, &! loop indices for cells, edges
         iTracer,      &! tracer index
         cell1, cell2, &! neighbor cell indices across edge
         numTracers

      real (kind=RKIND) :: &
         tracerWgt,        &! local temporary
         invAreaCell1,     &! inverse cell area
         verticalWeightK,  &! weights for vertical advection
         verticalWeightKm1,&! weights for vertical advection
         vertDivFactor      ! vertical divergence factor

      real (kind=RKIND), dimension(:,:), allocatable :: &
         tracerCur,        &! reordered tracer at current time
         highOrderHorzFlx, &! high-order flux in horizontal
         highOrderVertFlx   ! high-order flux in vertical

      real (kind=RKIND), parameter :: &
         eps = 1.e-10_RKIND  ! small value to avoid div by zero

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Get dimensions
      numTracers = size(tracers,dim=1)

      ! Allocate some arrays
      allocate(tracerCur       (nVertLevels, nCellsAll), &
               highOrderHorzFlx(nVertLevels, nEdgesAll), &
               highOrderVertFlx(nVertLevels, nCellsAll))

      ! Loop over tracers. One tracer is advected at a time. It is
      ! copied into a temporary array in order to improve locality
      do iTracer = 1, numTracers

        ! Initialize variables for use in this iTracer iteration
        !$omp parallel
        !$omp do schedule(runtime) private(k)
        do iCell = 1, nCellsAll
        do k=1,nVertLevels
           tracerCur(k,iCell) = tracers(iTracer,k,iCell)
           highOrderVertFlx(k, iCell) = 0.0_RKIND
        end do
        end do
        !$omp end do

        !$omp do schedule(runtime) private(k)
        do iEdge = 1, nEdgesAll
        do k=1,nVertLevels
           highOrderHorzFlx(k, iEdge) = 0.0_RKIND
        end do
        end do
        !$omp end do

        ! Compute the high order vertical flux.

        ! First take care of levels at top, bottom defaulting
        ! to lower order

        !$omp do schedule(runtime) &
        !$omp    private(k, verticalWeightK, verticalWeightKm1)
        do iCell = 1, nCellsAll
           ! top active level
           k = max(minLevelCell(iCell), &
                   min(maxLevelCell(iCell), minLevelCell(iCell)+1))
           verticalWeightK   = layerThickness(k-1,iCell) / &
                              (layerThickness(k  ,iCell) + &
                               layerThickness(k-1,iCell))
           verticalWeightKm1 = layerThickness(k  ,iCell) / &
                              (layerThickness(k  ,iCell) + &
                               layerThickness(k-1,iCell))
           highOrderVertFlx(k,iCell) = w(k,iCell)&
                            *(verticalWeightK  *tracerCur(k  ,iCell) + &
                              verticalWeightKm1*tracerCur(k-1,iCell))

           ! lowest active level
           k = max(minLevelCell(iCell), maxLevelCell(iCell))
           verticalWeightK   = layerThickness(k-1,iCell) / &
                              (layerThickness(k  ,iCell) + &
                               layerThickness(k-1,iCell))
           verticalWeightKm1 = layerThickness(k  ,iCell) / &
                              (layerThickness(k  ,iCell) + &
                               layerThickness(k-1,iCell))
           highOrderVertFlx(k,iCell) = w(k,iCell)* &
                            (verticalWeightK  *tracerCur(k  ,iCell) + &
                             verticalWeightKm1*tracerCur(k-1,iCell))
        end do ! iCell Loop
        !$omp end do

        ! Now computer interior layers at high order based on
        ! user-requested order

        select case (vertOrder)
        case (vertOrder4)

           !$omp do schedule(runtime) &
           !$omp    private(k, kmin, kmax, &
           !$omp            verticalWeightK, verticalWeightKm1)
           do iCell = 1, nCellsAll
              kmin = minLevelCell(iCell)
              kmax = maxLevelCell(iCell)
              do k=kmin+2,kmax-1
                 highOrderVertFlx(k, iCell) = w(k,iCell)* &
                                 ( 7.0_RKIND*(tracerCur(k  ,iCell)+ &
                                              tracerCur(k-1,iCell)) &
                                           - (tracerCur(k+1,iCell)+ &
                                              tracerCur(k-2,iCell)))/&
                                              12.0_RKIND
              end do ! vertical loop
           end do ! iCell Loop
           !$omp end do

        case (vertOrder3)

           !$omp do schedule(runtime) &
           !$omp    private(k, kmin, kmax, &
           !$omp            verticalWeightK, verticalWeightKm1)
           do iCell = 1, nCellsAll
              kmin = minLevelCell(iCell)
              kmax = maxLevelCell(iCell)
              do k=kmin+2,kmax-1
                 highOrderVertFlx(k, iCell) = (w(k,iCell)* &
                                 (7.0_RKIND * (tracerCur(k  ,iCell)+ &
                                               tracerCur(k-1,iCell)) - &
                                              (tracerCur(k+1,iCell)+ &
                                               tracerCur(k-2,iCell))) - &
                              coef3rdOrder*abs(w(k,iCell))* &
                                             ((tracerCur(k+1,iCell)- &
                                               tracerCur(k-2,iCell)) - &
                                    3.0_RKIND*(tracerCur(k  ,iCell)- &
                                               tracerCur(k-1,iCell))))/ &
                                               12.0_RKIND
              end do ! vertical loop
           end do ! iCell Loop
           !$omp end do

        case (vertOrder2)

           !$omp do schedule(runtime) &
           !$omp    private(k, kmin, kmax, &
           !$omp            verticalWeightK, verticalWeightKm1)
           do iCell = 1, nCellsAll
              kmin = minLevelCell(iCell)
              kmax = maxLevelCell(iCell)
              do k=kmin+2,kmax-1
                 verticalWeightK   = layerThickness(k-1,iCell) / &
                                    (layerThickness(k  ,iCell) + &
                                     layerThickness(k-1,iCell))
                 verticalWeightKm1 = layerThickness(k  ,iCell) / &
                                    (layerThickness(k  ,iCell) + &
                                     layerThickness(k-1,iCell))
                 highOrderVertFlx(k,iCell) = w(k, iCell) * &
                           (verticalWeightK  *tracerCur(k  ,iCell) + &
                            verticalWeightKm1*tracerCur(k-1,iCell))
              end do ! vertical loop
           end do ! iCell Loop
           !$omp end do

        end select ! vertOrder

        ! Compute the high order horizontal flux

        !$omp do schedule(runtime) &
        !$omp    private(cell1, cell2, k, kmix, kmax, tracerWgt, i, iCell)
        do iEdge = 1, nEdgesAll
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)
          kmin = minLevelEdgeBot(iEdge)
          kmax = maxLevelEdgeTop(iEdge)

          ! Compute 2nd order fluxes where needed.
          do k = kmin, kmax
            tracerWgt = (1.0_RKIND - advMaskHighOrder(k,iEdge)) &
                      * (dvEdge(iEdge) * 0.5_RKIND) &
                      * normalThicknessFlux(k, iEdge)

            highOrderHorzFlx(k, iEdge) = highOrderHorzFlx(k, iedge) &
                                       + tracerWgt*(tracerCur(k,cell1) &
                                                  + tracerCur(k,cell2))
          end do ! k loop

          ! Compute 3rd or 4th fluxes where requested.
          do i = 1, nAdvCellsForEdge(iEdge)
            iCell = advCellsForEdge(i,iEdge)
            kmin = minLevelCell(iCell)
            kmax = maxLevelCell(iCell)

            do k = kmin,kmax
              tracerWgt = advMaskHighOrder(k,iEdge) &
                        * (advCoefs(i,iEdge) + coef3rdOrder &
                        * sign(1.0_RKIND,normalThicknessFlux(k,iEdge)) &
                        * advCoefs3rd(i,iEdge))

              tracerWgt = normalThicknessFlux(k,iEdge)*tracerWgt
              highOrderHorzFlx(k,iEdge) = highOrderHorzFlx(k,iEdge) &
                                        + tracerWgt*tracerCur(k,iCell)
            end do ! k loop
          end do ! i loop over nAdvCellsForEdge
        end do ! iEdge loop
        !$omp end do

        ! Accumulate the scaled high order horizontal tendencies
        !$omp do schedule(runtime) &
        !$omp    private(invAreaCell1, i, iEdge, k, kmin, kmax)
        do iCell = 1, nCellsAll
          invAreaCell1 = 1.0_RKIND / areaCell(iCell)
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            kmin = minLevelEdgeBot(iEdge)
            kmax = maxLevelEdgeTop(iEdge)

            do k = kmin,kmax
              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) &
                                      + edgeSignOnCell(i, iCell) &
                                      * highOrderHorzFlx(k, iEdge) &
                                      * invAreaCell1
            end do
          end do
        end do
        !$omp end do

        ! Accumulate the scaled high order vertical tendencies.
        vertDivFactor = 1.0_RKIND
        !$omp do schedule(runtime) private(k,kmin,kmax)
        do iCell = 1, nCellsOwned
           kmin = minLevelCell(iCell)
           kmax = maxLevelCell(iCell)

           do k = kmin,kmax
              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + &
                         vertDivFactor*(highOrderVertFlx(k+1, iCell) &
                                      - highOrderVertFlx(k, iCell))
           end do ! k loop
        end do ! iCell loop
        !$omp end do
        !$omp end parallel
      end do ! iTracer loop

      deallocate(tracerCur, highOrderHorzFlx, highOrderVertFlx)

   end subroutine ocn_tracer_advection_std_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_advection_std_init
!
!> \brief MPAS initialize standard tracer advection tendency.
!> \author Doug Jacobsen, Phil Jones
!> \date   03/09/12, updated May 2019
!> \details
!>  This routine initializes constants and choices for the standard
!>  tracer advection tendency.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_advection_std_init(horzAdvOrder, vertAdvOrder, &
                                            inCoef3rdOrder,             &
                                            checkMonotonicity, err) !{{{

      !*** input parameters

      integer, intent(in) :: &
         horzAdvOrder         !< [in] Order for horizontal advection
      integer, intent(in) :: &
         vertAdvOrder         !< [in] Order for vertical advection
      real (kind=RKIND), intent(in) :: &
         inCoef3rdOrder       !< [in] Coefficient for blending advection orders
      logical, intent(in) :: &
         checkMonotonicity    !< [in] Flag to check on monotonicity of tracers

      !*** output parameters

      integer, intent(out) :: err !< [out] Error Flag

      ! end of preamble
      !----------------
      ! begin code

      err = 0 ! set error code to success

      ! set monotonicity flag
      monotonicityCheck = checkMonotonicity

      ! set 3rd order coefficient based on horizontal order choice
      select case (horzAdvOrder)
      case (2)
         coef3rdOrder = 0.0_RKIND
      case (3)
         coef3rdOrder = inCoef3rdOrder
      case (4)
         coef3rdOrder = 0.0_RKIND
      case default
         call mpas_log_write( &
            'Invalid value for horz advection order, defaulting to 2nd order', &
            MPAS_LOG_WARN)
      end select ! horzAdvOrder

      ! set choice of vertical advection order
      select case (vertAdvOrder)
      case (2)
         vertOrder = vertOrder2
      case (3)
         vertOrder = vertOrder3
      case (4)
         vertOrder = vertOrder4
      case default
         vertOrder = vertOrder2
         call mpas_log_write( &
         'Invalid value for vertical advection order, defaulting to 2nd order',&
         MPAS_LOG_WARN)
      end select ! vertAdvOrder

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_advection_std_init!}}}

!***********************************************************************

end module ocn_tracer_advection_std

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
