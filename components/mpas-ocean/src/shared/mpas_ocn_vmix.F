! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix
!
!> \brief MPAS ocean vertical mixing driver
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module is the main driver for
!>  vertical mixing in the ocean.
!>
!
!-----------------------------------------------------------------------

module ocn_vmix

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use ocn_mesh

   use mpas_constants
   use ocn_constants
   use ocn_config
   use ocn_vmix_cvmix
   use ocn_vmix_gotm
   use ocn_vmix_coefs_redi
   use ocn_diagnostics_variables

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

   public :: ocn_vmix_coefs, &
             ocn_vel_vmix_tend_implicit, &
             ocn_tracer_vmix_tend_implicit, &
             ocn_vmix_init, &
             ocn_vmix_implicit

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: velVmixOn, tracerVmixOn
   real (kind=RKIND) :: implicitBottomDragCoef
   real (kind=RKIND) :: rayleighDampingCoef, rayleighBottomDampingCoef, &
                        rayleighDepthVariable

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs
!
!> \brief   Computes coefficients for vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers based user choices of mixing parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs(meshPool, statePool, forcingPool, scratchPool, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      type (mpas_pool_type), intent(in) :: scratchPool !< Input/Output: Scratch structure

      integer, intent(in), optional :: timeLevelIn !< Input: Time level for state pool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool             !< Input/Output: state information

      type (mpas_pool_type), intent(inout) :: &
         forcingPool             !< Input/Output: forcing information

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

      integer :: err1, err2, err3
      integer :: timeLevel

      integer :: iEdge, iCell, k, nEdges, nCells

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing coefficients
      !
      !-----------------------------------------------------------------

      err = 0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if


      nEdges = nEdgesHalo( 1 )

#ifdef MPAS_OPENACC
      !$acc parallel loop collapse(2) present(vertViscTopOfEdge)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iEdge = 1, nEdges
      do k = 1, nVertLevels
         vertViscTopOfEdge(k, iEdge) = 0.0_RKIND
      end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nCells = nCellsHalo( 1 )

#ifdef MPAS_OPENACC
      !$acc parallel loop collapse(2) present(vertDiffTopOfCell)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCells
      do k = 1, nVertLevels
         vertDiffTopOfCell(k, iCell) = 0.0_RKIND
      end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc update host(vertViscTopOfEdge, vertDiffTopOfCell)
#endif

      call ocn_vmix_coefs_cvmix_build(meshPool, statePool, forcingPool, err1, timeLevel)
      call ocn_vmix_coefs_redi_build(meshPool, statePool, err2, timeLevel)
      call ocn_vmix_coefs_gotm_build(statePool, forcingPool, err3, timeLevel)

      err = ior(ior(err1, err2), err3)

#ifdef MPAS_OPENACC
      !$acc update device(vertDiffTopOfCell, vertViscTopOfCell, vertNonLocalFlux)
#endif
   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit_rayleigh
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!> \author  Mark Petersen, Phillip J. Wolfram
!> \date    September 2011, July 2019
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using computed coefficients and includes implicit rayleigh drag.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_implicit_rayleigh(dt, kineticEnergyCell, vertViscTopOfEdge, layerThickness, & !{{{
                                                  layerThickEdgeMean, normalVelocity, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeMean !< Input: mean thickness at edge

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

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

      integer :: iEdge, k, cell1, cell2, Nsurf, N

      real (kind=RKIND), dimension(:), allocatable :: bTemp, C, rTemp
      real (kind=RKIND) :: A, m
      real (kind=RKIND) :: edgeThicknessTotal

      err = 0

      if(.not.velVmixOn) return

      allocate(bTemp(nVertLevels),C(nVertLevels),rTemp(nVertLevels))

#ifdef MPAS_OPENACC
      !$acc enter data create(bTemp, C, rTemp)

      !$acc parallel loop present(maxLevelEdgeTop, minLevelEdgeBot, cellsOnEdge, &
      !$acc    layerThickEdgeMean, layerThickness, vertViscTopOfEdge, normalVelocity, &
      !$acc    kineticEnergyCell) &
      !$acc    private(Nsurf, N, cell1, cell2, edgeThicknessTotal, A, bTemp, C, m, rTemp, k)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(Nsurf, N, cell1, cell2, edgeThicknessTotal, k, A, bTemp, C, m, rTemp)
#endif
      do iEdge = 1, nEdgesOwned
        Nsurf = minLevelEdgeBot(iEdge)
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         edgeThicknessTotal = 0.0_RKIND
         do k = Nsurf, N
            edgeThicknessTotal = edgeThicknessTotal + layerThickEdgeMean(k,iEdge)
         enddo

         ! one active layer
         if (N .eq. Nsurf) then
           normalVelocity(N,iEdge) = normalVelocity(N,iEdge) &
               / (1.0_RKIND + dt*implicitBottomDragCoef &
               * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge) &
               ! added Rayleigh terms
               + dt*(rayleighBottomDampingCoef + rayleighDampingCoef &
               / ((1.0_RKIND - rayleighDepthVariable) + rayleighDepthVariable*edgeThicknessTotal)))
         else

           ! tridiagonal matrix algorithm
           C(Nsurf) = -2.0_RKIND*dt*vertViscTopOfEdge(Nsurf+1,iEdge) &
                  / (layerThickEdgeMean(Nsurf,iEdge) + layerThickEdgeMean(Nsurf+1,iEdge)) &
                  / layerThickEdgeMean(Nsurf,iEdge)
           bTemp(Nsurf) = 1.0_RKIND - C(Nsurf) &
              ! added Rayleigh terms
              + dt*rayleighDampingCoef/((1.0_RKIND - rayleighDepthVariable) + rayleighDepthVariable*edgeThicknessTotal)
           rTemp(Nsurf) = normalVelocity(Nsurf,iEdge)

           ! first pass: set the coefficients
           do k = Nsurf+1, N-1
              A        = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
                     / (layerThickEdgeMean(k-1,iEdge) + layerThickEdgeMean(k,iEdge)) &
                     / layerThickEdgeMean(k,iEdge)
              m        = A/bTemp(k-1)
              C(k)     = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
                     / (layerThickEdgeMean(k,iEdge) + layerThickEdgeMean(k+1,iEdge)) &
                     / layerThickEdgeMean(k,iEdge)
              bTemp(k) = 1.0_RKIND - A - C(k) &
                     ! added Rayleigh terms
                     + dt*(rayleighDampingCoef/((1.0_RKIND - rayleighDepthVariable) + rayleighDepthVariable*edgeThicknessTotal)) &
                     - m*C(k-1)
              rTemp(k) = normalVelocity(k,iEdge) - m*rTemp(k-1)
           enddo

           A = -2.0_RKIND*dt*vertViscTopOfEdge(N,iEdge) &
             / (layerThickEdgeMean(N-1,iEdge) + layerThickEdgeMean(N,iEdge)) &
             / layerThickEdgeMean(N,iEdge)
           m = A/bTemp(N-1)

           ! x(N) = rTemp(N) / bTemp(N)
           ! Apply bottom drag boundary condition on the viscous term
           ! using sqrt(2.0*kineticEnergyEdge(k,iEdge))
           normalVelocity(N,iEdge) = (normalVelocity(N,iEdge) - m*rTemp(N-1)) &
               / (1.0_RKIND - A + dt*implicitBottomDragCoef &
               * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge) &
               ! added Rayleigh terms
               + dt*(rayleighBottomDampingCoef + rayleighDampingCoef &
               / ((1.0_RKIND - rayleighDepthVariable) + rayleighDepthVariable*edgeThicknessTotal)) &
               - m*C(N-1))
           ! second pass: back substitution
           do k = N-1, Nsurf, -1
              normalVelocity(k,iEdge) = (rTemp(k) - C(k)*normalVelocity(k+1,iEdge)) / bTemp(k)
           enddo
           normalVelocity(1:Nsurf-1,iEdge) = 0.0_RKIND
           normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

         end if ! N=Nsurf, i.e. single layer
        end if
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(bTemp, C, rTemp)
#endif

      deallocate(bTemp,C,rTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit_rayleigh!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using computed coefficients.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_implicit(dt, kineticEnergyCell, vertViscTopOfEdge, layerThickness, & !{{{
                                         layerThickEdgeMean, normalVelocity, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeMean !< Input: mean thickness at edge

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

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

      integer :: iEdge, k, cell1, cell2, N, Nsurf

      real (kind=RKIND), dimension(:), allocatable :: bTemp, C, rTemp
      real (kind=RKIND) :: A, m

      err = 0

      if(.not.velVmixOn) return

      allocate(bTemp(nVertLevels),C(nVertLevels),rTemp(nVertLevels))

#ifdef MPAS_OPENACC
      !$acc enter data create(bTemp, C, rTemp)

      !$acc parallel loop present(maxLevelEdgeTop, minLevelEdgeBot, cellsOnEdge, &
      !$acc    layerThickEdgeMean, layerThickness, vertViscTopOfEdge, normalVelocity, &
      !$acc    kineticEnergyCell) &
      !$acc    private(Nsurf, N, cell1, cell2, A, bTemp, C, m, rTemp, k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(Nsurf, N, cell1, cell2, A, bTemp, C, m, rTemp, k)
#endif
      do iEdge = 1, nEdgesOwned
        N = maxLevelEdgeTop(iEdge)
        Nsurf = minLevelEdgeBot(iEdge)
        if (N .gt. 0) then

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! one active layer
         if (N .eq. Nsurf) then
            normalVelocity(N,iEdge) = normalVelocity(N,iEdge)  &
                / (1.0_RKIND + dt*implicitBottomDragCoef &
                   * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge) )
         else

           ! tridiagonal matrix algorithm
           C(Nsurf)     = -2.0_RKIND*dt*vertViscTopOfEdge(Nsurf+1,iEdge) &
                      / (layerThickEdgeMean(Nsurf,iEdge) + layerThickEdgeMean(Nsurf+1,iEdge)) &
                      / layerThickEdgeMean(Nsurf,iEdge)
           bTemp(Nsurf) = 1.0_RKIND - C(Nsurf)
           rTemp(Nsurf) = normalVelocity(Nsurf,iEdge)

           ! first pass: set the coefficients
           do k = Nsurf+1, N-1
              A        = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
                     / (layerThickEdgeMean(k-1,iEdge) + layerThickEdgeMean(k,iEdge)) &
                     / layerThickEdgeMean(k,iEdge)
              m        = A/bTemp(k-1)
              C(k)     = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
                     / (layerThickEdgeMean(k,iEdge) + layerThickEdgeMean(k+1,iEdge)) &
                     / layerThickEdgeMean(k,iEdge)
              bTemp(k) = 1.0_RKIND - A - C(k) - m*C(k-1)
              rTemp(k) = normalVelocity(k,iEdge) - m*rTemp(k-1)
           enddo

           A = -2.0_RKIND*dt*vertViscTopOfEdge(N,iEdge) &
             / (layerThickEdgeMean(N-1,iEdge) + layerThickEdgeMean(N,iEdge)) &
             / layerThickEdgeMean(N,iEdge)
           m = A/bTemp(N-1)

           ! x(N) = rTemp(N) / bTemp(N)
           ! Apply bottom drag boundary condition on the viscous term
           ! using sqrt(2.0*kineticEnergyEdge(k,iEdge))
           normalVelocity(N,iEdge) = (normalVelocity(N,iEdge) - m*rTemp(N-1)) &
               / (1.0_RKIND - A + dt*implicitBottomDragCoef &
               * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge) &
               - m*C(N-1))
           ! second pass: back substitution
           do k = N-1, Nsurf, -1
              normalVelocity(k,iEdge) = (rTemp(k) - C(k)*normalVelocity(k+1,iEdge)) / bTemp(k)
           enddo
           normalVelocity(1:Nsurf-1,iEdge) = 0.0_RKIND
           normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

         end if ! N=Nsurf, i.e. single layer
        end if
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(bTemp, C, rTemp)
#endif

      deallocate(bTemp,C,rTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit_variable
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!>          with spatially-variable bottom drag
!> \author  Mark Petersen, Phillip J. Wolfram
!> \date    September 2011, September 2019
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using computed coefficients from spatially-variable bottom drag.
!>  Except for bottom drag coefficient, routine should be identifcal to
!>  ocn_vel_vmix_tend_implicit above.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_implicit_spatially_variable(bottomDrag, dt, kineticEnergyCell, & !{{{
                                         vertViscTopOfEdge, layerThickness, &
                                         layerThickEdgeMean, normalVelocity, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeMean !< Input: mean thickness at edge

       real (kind=RKIND), dimension(:), intent(in) :: &
         bottomDrag !< Input: bottomDrag at cell centeres

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

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

      integer :: iEdge, k, cell1, cell2, N, Nsurf
      real (kind=RKIND) :: implicitCd

      real (kind=RKIND), dimension(:), allocatable :: bTemp, C, rTemp
      real (kind=RKIND) :: A, m

      err = 0

      if(.not.velVmixOn) return

      allocate(bTemp(nVertLevels),C(nVertLevels),rTemp(nVertLevels))

#ifdef MPAS_OPENACC
      !$acc enter data create(bTemp, C, rTemp)
      !$acc enter data create(bottomDrag)

      !$acc parallel loop present(maxLevelEdgeTop, minLevelEdgeBot, cellsOnEdge, &
      !$acc    layerThickEdgeMean, layerThickness, vertViscTopOfEdge, normalVelocity, &
      !$acc    kineticEnergyCell, bottomDrag) &
      !$acc    private(Nsurf, N, cell1, cell2, implicitCd, A, bTemp, C, m, rTemp, k)
#else
      !$omp do schedule(runtime)
#endif
      do iEdge = 1, nEdgesOwned
        Nsurf = minLevelEdgeBot(iEdge)
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         ! average cell-based implicit bottom drag to edges
         implicitCd = 0.5_RKIND*(bottomDrag(cell1) + bottomDrag(cell2))

         ! one active layer
         if (N .eq. Nsurf) then
            normalVelocity(N,iEdge) = normalVelocity(N,iEdge)  &
                / (1.0_RKIND + dt*implicitCD &
                   * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge) )
         else

           ! tridiagonal matrix algorithm
           C(Nsurf)     = -2.0_RKIND*dt*vertViscTopOfEdge(Nsurf+1,iEdge) &
                      / (layerThickEdgeMean(Nsurf,iEdge) + layerThickEdgeMean(Nsurf+1,iEdge)) &
                      / layerThickEdgeMean(Nsurf,iEdge)
           bTemp(Nsurf) = 1.0_RKIND - C(Nsurf)
           rTemp(Nsurf) = normalVelocity(Nsurf,iEdge)

           ! first pass: set the coefficients
           do k = Nsurf+1, N-1
              A        = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
                     / (layerThickEdgeMean(k-1,iEdge) + layerThickEdgeMean(k,iEdge)) &
                     / layerThickEdgeMean(k,iEdge)
              m        = A/bTemp(k-1)
              C(k)     = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
                     / (layerThickEdgeMean(k,iEdge) + layerThickEdgeMean(k+1,iEdge)) &
                     / layerThickEdgeMean(k,iEdge)
              bTemp(k) = 1.0_RKIND - A - C(k) - m*C(k-1)
              rTemp(k) = normalVelocity(k,iEdge) - m*rTemp(k-1)
           enddo

           A = -2.0_RKIND*dt*vertViscTopOfEdge(N,iEdge) &
             / (layerThickEdgeMean(N-1,iEdge) + layerThickEdgeMean(N,iEdge)) &
             / layerThickEdgeMean(N,iEdge)
           m = A/bTemp(N-1)

           ! x(N) = rTemp(N) / bTemp(N)
           ! Apply bottom drag boundary condition on the viscous term
           ! using sqrt(2.0*kineticEnergyEdge(k,iEdge))
           normalVelocity(N,iEdge) = (normalVelocity(N,iEdge) - m*rTemp(N-1)) &
               / (1.0_RKIND - A + dt*implicitCD &
               * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge) &
               - m*C(N-1))
           ! second pass: back substitution
           do k = N-1, Nsurf, -1
              normalVelocity(k,iEdge) = (rTemp(k) - C(k)*normalVelocity(k+1,iEdge)) / bTemp(k)
           enddo
           normalVelocity(1:Nsurf-1,iEdge) = 0.0_RKIND
           normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND
         end if ! one active layer
        end if
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(bTemp, C, rTemp)
      !$acc exit data delete(bottomDrag)
#endif
      deallocate(bTemp,C,rTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit_spatially_variable!}}}


!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit_variable_mannings
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!>          with spatially-variable bottom drag using Mannings friction
!> \author  Mark Petersen, Phillip J. Wolfram
!> \date    September 2011, September 2019
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using computed coefficients from spatially-variable bottom drag.
!>  Except for bottom drag coefficient, routine should be identifcal to
!>  ocn_vel_vmix_tend_implicit above.  Cd uses Mannings' n values for the
!>  Cd=g*n^2*h^{-1/3}.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_implicit_spatially_variable_mannings(forcingPool, bottomDrag, dt, & !{{{
                                         kineticEnergyCell, vertViscTopOfEdge, layerThickness, &
                                         layerThickEdgeMean, normalVelocity, ssh, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         forcingPool          !< Input: forcing information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeMean !< Input: mean thickness at edge

       real (kind=RKIND), dimension(:), intent(in) :: &
         bottomDrag !< Input: bottomDrag at cell centeres

       real (kind=RKIND), dimension(:), pointer :: ssh

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

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

      integer :: iEdge, k, cell1, cell2, N, Nsurf
      real (kind=RKIND) :: implicitCd

      real (kind=RKIND), dimension(:), allocatable :: bTemp, C, rTemp
      real (kind=RKIND) :: A, m

      ! vegetation_drag
      real (kind=RKIND), dimension(:), pointer :: vegetationHeight
      real (kind=RKIND), dimension(:), pointer :: vegetationDiameter
      real (kind=RKIND), dimension(:), pointer ::vegetationDensity
      real (kind=RKIND), dimension(:), pointer ::vegetationManning
      integer, dimension(:), pointer ::vegetationMask
      real (kind=RKIND) :: old_bottom_Cd, lambda, beta, alpha, total_h
      real (kind=RKIND) :: inundation_depth, von_karman, cff1, cff2, cff3, cff4
      integer :: iCell

      err = 0
      von_karman = 0.4_RKIND

      if(.not.velVmixOn) return

      call mpas_pool_get_array(forcingPool, 'vegetationMask', vegetationMask)
      call mpas_pool_get_array(forcingPool, 'vegetationHeight', vegetationHeight)
      call mpas_pool_get_array(forcingPool, 'vegetationDensity', vegetationDensity)
      call mpas_pool_get_array(forcingPool, 'vegetationDiameter', vegetationDiameter)
      call mpas_pool_get_array(forcingPool, 'vegetationManning', vegetationManning)

      allocate(bTemp(nVertLevels),C(nVertLevels),rTemp(nVertLevels))

#ifdef MPAS_OPENACC
      !$acc enter data create(bTemp, C, rTemp)
      !$acc enter data create(ssh, bottomDrag)
#endif

      ! Compute bottomDrag (Manning roughness) induced by vegetation
      if (config_use_vegetation_drag .AND. config_use_vegetation_manning_equation) then
        do iCell = 1, nCellsOwned
          if (vegetationDensity(iCell) * vegetationHeight(iCell) * vegetationDiameter(iCell) .eq. 0.0_RKIND) then
            vegetationMask(iCell) = 0
          endif
          if (vegetationMask(iCell) .eq. 1) then
            total_h = bottomDepth(iCell) + ssh(iCell)
            old_bottom_Cd = gravity * bottomDrag(iCell)**2 * total_h**(-1.0_RKIND/3.0_RKIND)
            inundation_depth = MIN(vegetationHeight(iCell), total_h)
            inundation_depth = MAX(inundation_depth, 1e-6)
            lambda = vegetationDiameter(iCell) * vegetationDensity(iCell)
            alpha = (config_vegetation_drag_coefficient*lambda/ &
                   (4.0_RKIND*von_karman**2 * inundation_depth**2))**(1.0_RKIND/3.0_RKIND)
            beta = 0.5_RKIND*alpha*inundation_depth*(1.0_RKIND - EXP(-2.0_RKIND*alpha*inundation_depth)) &
                  / (1.0_RKIND - EXP(-alpha*inundation_depth))**2
            cff1 = total_h**(2.0_RKIND/3.0_RKIND) &
                  * SQRT((0.5_RKIND*beta*lambda*config_vegetation_drag_coefficient*inundation_depth &
                  + old_bottom_Cd)/(gravity*total_h))
            cff2 = (alpha*inundation_depth)**2/(1.0_RKIND - EXP(-alpha*inundation_depth))
            cff3 = (1.0_RKIND - EXP(-alpha*inundation_depth))/(alpha**2 * inundation_depth*total_h)
            cff4 = LOG(total_h/inundation_depth) - (1.0_RKIND-inundation_depth/total_h) &
                  * (1.0_RKIND - 1.0_RKIND/(alpha*inundation_depth))
            vegetationManning(iCell) = cff1/(cff2*(cff3+cff4))
            vegetationManning(iCell) = MAX(bottomDrag(iCell), vegetationManning(iCell))
          else
            vegetationManning(iCell) = bottomDrag(iCell)
          endif
        enddo
      endif

#ifdef MPAS_OPENACC
      !$acc enter data create(vegetationManning)

      !$acc parallel loop present(maxLevelEdgeTop, minLevelEdgeBot, cellsOnEdge, &
      !$acc    layerThickEdgeMean, layerThickness, vertViscTopOfEdge, normalVelocity, &
      !$acc    kineticEnergyCell, bottomDrag, vegetationManning, ssh, bottomDepth) &
      !$acc    private(Nsurf, N, cell1, cell2, implicitCd, A, bTemp, C, m, rTemp, k)
#else
      !$omp do schedule(runtime)
#endif
      do iEdge = 1, nEdgesOwned
        Nsurf = minLevelEdgeBot(iEdge)
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! average cell-based implicit bottom drag to edges and convert Mannings n to Cd
         if (config_use_implicit_bottom_roughness) then
           ! NCOM's formula for bottom drag coefficient
           ! The original formula is cb(i,j) = max(cbmin, (vonk/log(0.5*depth(i,j)/z0))**2 ))
           ! where vonk=0.16, z0=0.001, bracketed between [0.0025, 0.1]
           ! The constant below is 250 = 0.5*0.5/z0 with the extra 0.5 from
           ! averaging the two bottomDepth cells.
           implicitCd = max(0.0025_RKIND, &
                        min(0.1_RKIND, &
                          0.16_RKIND/(log(250.0_RKIND*(bottomDepth(cell1)+bottomDepth(cell2))))**2 ))
         elseif (config_use_vegetation_drag .AND. config_use_vegetation_manning_equation) then
           implicitCd = gravity*(0.5_RKIND*(vegetationManning(cell1) + vegetationManning(cell2)))**2 * &
            (0.5_RKIND * (ssh(cell1) + ssh(cell2) + bottomDepth(cell1) + bottomDepth(cell2)))**(-1.0_RKIND/3.0_RKIND)
         else
           implicitCd = gravity*(0.5_RKIND*(bottomDrag(cell1) + bottomDrag(cell2)))**2 * &
            (0.5_RKIND * (ssh(cell1) + ssh(cell2) + bottomDepth(cell1) + bottomDepth(cell2)))**(-1.0_RKIND/3.0_RKIND)
         endif

         ! one active layer
         if (N .eq. Nsurf) then
           normalVelocity(N,iEdge) = normalVelocity(N,iEdge)  &
                / (1.0_RKIND + dt*implicitCD &
                   * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge) )
         else

           ! tridiagonal matrix algorithm
           C(Nsurf)     = -2.0_RKIND*dt*vertViscTopOfEdge(Nsurf+1,iEdge) &
                      / (layerThickEdgeMean(Nsurf,iEdge) + layerThickEdgeMean(Nsurf+1,iEdge)) &
                      / layerThickEdgeMean(Nsurf,iEdge)
           bTemp(Nsurf) = 1.0_RKIND - C(Nsurf)
           rTemp(Nsurf) = normalVelocity(Nsurf,iEdge)

           ! first pass: set the coefficients
           do k = Nsurf+1, N-1
              A        = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
                     / (layerThickEdgeMean(k-1,iEdge) + layerThickEdgeMean(k,iEdge)) &
                     / layerThickEdgeMean(k,iEdge)
              m        = A/bTemp(k-1)
              C(k)     = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
                     / (layerThickEdgeMean(k,iEdge) + layerThickEdgeMean(k+1,iEdge)) &
                     / layerThickEdgeMean(k,iEdge)
              bTemp(k) = 1.0_RKIND - A - C(k) - m*C(k-1)
              rTemp(k) = normalVelocity(k,iEdge) - m*rTemp(k-1)
           enddo

           A = -2.0_RKIND*dt*vertViscTopOfEdge(N,iEdge) &
             / (layerThickEdgeMean(N-1,iEdge) + layerThickEdgeMean(N,iEdge)) &
             / layerThickEdgeMean(N,iEdge)
           m = A/bTemp(N-1)

           ! x(N) = rTemp(N) / bTemp(N)
           ! Apply bottom drag boundary condition on the viscous term
           ! using sqrt(2.0*kineticEnergyEdge(k,iEdge))
           normalVelocity(N,iEdge) = (normalVelocity(N,iEdge) - m*rTemp(N-1)) &
               / (1.0_RKIND - A + dt*implicitCD &
               * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge) &
               - m*C(N-1))
           ! second pass: back substitution
           do k = N-1, Nsurf, -1
              normalVelocity(k,iEdge) = (rTemp(k) - C(k)*normalVelocity(k+1,iEdge)) / bTemp(k)
           enddo
           normalVelocity(1:Nsurf-1,iEdge) = 0.0_RKIND
           normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

         endif ! one active layer
        end if
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(bTemp, C, rTemp)
      !$acc exit data delete(ssh, bottomDrag, vegetationManning)
#endif
      deallocate(bTemp,C,rTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit_spatially_variable_mannings!}}}


!***********************************************************************
!
!  routine ocn_tracer_vmix_tend_implicit
!
!> \brief   Computes tendencies for implicit tracer vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for
!>  tracers using computed coefficients.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vmix_tend_implicit(dt, vertDiffTopOfCell, layerThickness, tracers, &
                  vertNonLocalFlux, tracerGroupSurfaceFlux, config_cvmix_kpp_nonlocal_with_implicit_mix, &
                  err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertDiffTopOfCell !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness, &             !< Input: thickness at cell center
         tracerGroupSurfaceFlux        !< Input: surface flux for tracers nonlocal computation

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         vertNonLocalFlux             !non local flux at interfaces

      logical, intent(in) :: config_cvmix_kpp_nonlocal_with_implicit_mix
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracers        !< Input: tracers

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

      integer :: iCell, k, iTracer, num_tracers, N, Nsurf

      real (kind=RKIND), dimension(:), allocatable :: bTemp, C
      real (kind=RKIND), dimension(:,:), allocatable :: rTemp
      real (kind=RKIND) :: A, m

      err = 0

      if(.not.tracerVmixOn) return

      num_tracers = size(tracers, dim=1)

      allocate(bTemp(nVertLevels),C(nVertLevels))
      allocate(rTemp(num_tracers,nVertLevels))

      call mpas_timer_start('vmix tracers tend imp loop', .false.)

#ifdef MPAS_OPENACC
      !$acc enter data create(bTemp, C, rTemp)

      !$acc parallel loop present(maxLevelCell, minLevelCell, vertDiffTopOfCell, &
      !$acc    layerThickness, tracers, tracerGroupSurfaceFlux, vertNonLocalFlux) &
      !$acc    private(Nsurf, N, A, bTemp, C, m, rTemp, iTracer, k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(Nsurf, N, A, bTemp, C, m, rTemp, iTracer, k)
#endif
      do iCell = 1, nCellsOwned
         Nsurf = minLevelCell(iCell)
         N = maxLevelCell(iCell)

         ! tridiagonal matrix algorithm
         C(Nsurf) = -2.0_RKIND*dt*vertDiffTopOfCell(Nsurf+1,iCell) &
                    / (layerThickness(Nsurf,iCell) + layerThickness(Nsurf+1,iCell)) &
                    / layerThickness(Nsurf,iCell)
         bTemp(Nsurf) = 1.0_RKIND - C(Nsurf)
         if ( config_cvmix_kpp_nonlocal_with_implicit_mix ) then
            do iTracer = 1, num_tracers
               rTemp(iTracer,Nsurf) = tracers(iTracer,Nsurf,iCell) + dt*tracerGroupSurfaceFlux(iTracer,iCell) &
                    * (-vertNonLocalFlux(1, Nsurf+1,iCell) )/ layerThickness(Nsurf,iCell)
            enddo
         else
            do iTracer = 1, num_tracers
               rTemp(iTracer,Nsurf) = tracers(iTracer,Nsurf,iCell)
            enddo
         endif

         ! first pass: set the coefficients
         do k = Nsurf+1, N-1
            A        = -2.0_RKIND*dt*vertDiffTopOfCell(k,iCell) &
                   / (layerThickness(k-1,iCell) + layerThickness(k,iCell)) / layerThickness(k,iCell)
            m        = A/bTemp(k-1)
            C(k)     = -2.0_RKIND*dt*vertDiffTopOfCell(k+1,iCell) &
                   / (layerThickness(k,iCell) + layerThickness(k+1,iCell)) / layerThickness(k,iCell)
            bTemp(k) = 1.0_RKIND - A - C(k) - m*C(k-1)
            if ( config_cvmix_kpp_nonlocal_with_implicit_mix ) then
               do iTracer = 1, num_tracers
                  rTemp(iTracer,k) = tracers(iTracer,k,iCell) + dt*tracerGroupSurfaceFlux(iTracer,iCell) &
                       * (vertNonLocalFlux(1,k,iCell) - vertNonLocalFlux(1,k+1,iCell)) / layerThickness(k,iCell) &
                       - m*rTemp(iTracer,k-1)
               enddo
            else
               do iTracer = 1, num_tracers
                  rTemp(iTracer,k) = tracers(iTracer,k,iCell) - m*rTemp(iTracer,k-1)
               enddo
            endif
         enddo

         A = -2.0_RKIND*dt*vertDiffTopOfCell(N,iCell) &
           / (layerThickness(N-1,iCell) + layerThickness(N,iCell)) / layerThickness(N,iCell)
         m = A/bTemp(N-1)

         ! x(N) = rTemp(N) / bTemp(N)
         if ( config_cvmix_kpp_nonlocal_with_implicit_mix ) then
            do iTracer = 1, num_tracers
               tracers(iTracer,N,iCell) = (tracers(iTracer,N,iCell) + dt*tracerGroupSurfaceFlux(iTracer,iCell) &
                      * vertNonLocalFlux(1,N,iCell) / layerThickness(N,iCell) &
                      - m*rTemp(iTracer,N-1)) / (1.0_RKIND - A - m*C(N-1))
            enddo
         else
            do iTracer = 1, num_tracers
               tracers(iTracer,N,iCell) = (tracers(iTracer,N,iCell) - m*rTemp(iTracer,N-1)) &
                      / (1.0_RKIND - A - m*C(N-1))
            enddo
         endif
         ! second pass: back subsititution
         do k = N-1, Nsurf, -1
            do iTracer = 1, num_tracers
               tracers(iTracer,k,iCell) = (rTemp(iTracer,k) - C(k)*tracers(iTracer,k+1,iCell)) / bTemp(k)
            enddo
         enddo
         tracers(:,1:Nsurf-1,iCell) = -1e34
         tracers(:,N+1:nVertLevels,iCell) = -1e34

      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(bTemp, C, rTemp)
#endif
      call mpas_timer_stop('vmix tracers tend imp loop')

      deallocate(bTemp, C, rTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vmix_tend_implicit!}}}

!***********************************************************************
!
!  routine ocn_vmix_implicit
!
!> \brief   Driver for implicit vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine is a driver for handling implicit vertical mixing
!>  of both momentum and tracers for a block. It's intended to reduce
!>  redundant code.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_implicit(dt, meshPool, statePool, forcingPool, scratchPool, err, timeLevelIn)!{{{
      real (kind=RKIND), intent(in) :: dt
      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool
      type (mpas_pool_type), intent(inout) :: forcingPool
      type (mpas_pool_type), intent(in) :: scratchPool !< Input/Output: Scratch structure
      integer, intent(out) :: err
      integer, intent(in), optional :: timeLevelIn

      type (mpas_pool_type), pointer :: tracersPool, tracersSurfaceFluxPool

      integer :: iCell, timeLevel, k, cell1, cell2, iEdge, iTracer, num_tracers, nCells, nEdges, N, Nsurf
      real (kind=RKIND), dimension(:), pointer :: bottomDrag, ssh
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: tracerGroupSurfaceFlux
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup
      real (kind=RKIND), dimension(:,:,:), allocatable :: nonLocalFluxTend

      type (mpas_pool_iterator_type) :: groupItr

      character (len=StrKIND) :: modifiedGroupName
      err = 0

      call mpas_timer_start('vmix imp')

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

      call mpas_pool_get_array(forcingPool, 'bottomDrag', bottomDrag)

      call mpas_timer_start('vmix coefs', .false.)
      call ocn_vmix_coefs(meshPool, statePool, forcingPool, scratchPool, err, timeLevel)
      call mpas_timer_stop('vmix coefs')

      ! if using CVMix, then viscosity has to be averaged from cell centers to cell edges
      if ( config_use_cvmix ) then

         call mpas_timer_start('CVMix avg', .false.)
#ifdef MPAS_OPENACC
         !$acc parallel loop present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, &
         !$acc    vertViscTopOfCell, vertViscTopOfEdge) &
         !$acc    private(cell1, cell2, k, Nsurf, N)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(cell1, cell2, k, Nsurf, N)
#endif
         do iEdge=1,nEdgesOwned
            Nsurf = minLevelEdgeBot(iEdge)
            N = maxLevelEdgeTop(iEdge)
            vertViscTopOfEdge(1:Nsurf-1, iEdge) = 0.0_RKIND
            vertViscTopOfEdge(N+1:nVertLevels, iEdge) = 0.0_RKIND
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
#ifdef MPAS_OPENACC
            !$acc loop vector
#endif
            do k=Nsurf, N
               vertViscTopOfEdge(k,iEdge) = 0.5_RKIND*(vertViscTopOfCell(k,cell2)+vertViscTopOfCell(k,cell1))
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
         call mpas_timer_stop('CVMix avg')
      endif

      ! if using GOTM, then viscosity has to be averaged from cell centers to cell edges
      if ( config_use_gotm ) then

         call mpas_timer_start('GOTM avg', .false.)
#ifdef MPAS_OPENACC
         !$acc parallel loop present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, &
         !$acc    vertViscTopOfCell, vertViscTopOfEdge) &
         !$acc    private(cell1, cell2, k, Nsurf, N)
#else
         !$omp do schedule(runtime) private(cell1, cell2, k, Nsurf, N)
#endif
         do iEdge=1,nEdgesOwned
            Nsurf = minLevelEdgeBot(iEdge)
            N = maxLevelEdgeTop(iEdge)
            vertViscTopOfEdge(1:Nsurf-1, iEdge) = 0.0_RKIND
            vertViscTopOfEdge(N+1:nVertLevels, iEdge) = 0.0_RKIND
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
#ifdef MPAS_OPENACC
            !$acc loop vector
#endif
            do k=Nsurf, N
               vertViscTopOfEdge(k,iEdge) = 0.5_RKIND*(vertViscTopOfCell(k,cell2)+vertViscTopOfCell(k,cell1))
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
#endif
         call mpas_timer_stop('GOTM avg')
      endif
#ifdef MPAS_OPENACC
      !$acc update host(vertViscTopOfEdge)
#endif

      !
      !  Implicit vertical solve for momentum
      !
      call mpas_timer_start('vmix solve momentum', .false.)
#ifdef MPAS_OPENACC
      !$acc enter data copyin(layerThickness, normalVelocity)
#endif
      if (config_use_implicit_bottom_drag_variable) then
        call ocn_vel_vmix_tend_implicit_spatially_variable(bottomDrag, dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdgeMean, normalVelocity, err)
      else if (config_use_implicit_bottom_drag_variable_mannings.or. &
               config_use_implicit_bottom_roughness) then
        ! update bottomDrag via Cd=g*n^2*h^(-1/3)
        call ocn_vel_vmix_tend_implicit_spatially_variable_mannings(forcingPool, bottomDrag, &
          dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdgeMean, normalVelocity, &
          ssh, err)
      else if (config_Rayleigh_friction.or. &
               config_Rayleigh_bottom_friction.or. &
               config_Rayleigh_damping_depth_variable) then
        call ocn_vel_vmix_tend_implicit_rayleigh(dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdgeMean, normalVelocity, err)
      else
        call ocn_vel_vmix_tend_implicit(dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdgeMean, normalVelocity, err)
      end if
#ifdef MPAS_OPENACC
      !$acc exit data copyout(normalVelocity)
#endif
      call mpas_timer_stop('vmix solve momentum')

      !
      !  Implicit vertical solve for all tracers
      !

      call mpas_timer_start('vmix solve tracers', .false.)
      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )

         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroup, timeLevel)
#ifdef MPAS_OPENACC
            !$acc enter data copyin(tracersGroup)
#endif
            ! store tracers
            if (trim(groupItr % memberName) == 'activeTracers') then
               if (config_compute_active_tracer_budgets) then
                  num_tracers = size(activeTracerVertMixTendency, dim=1)
#ifdef MPAS_OPENACC
                  !$acc parallel loop gang vector collapse(3) present(activeTracerVertMixTendency, tracersGroup)
#else
                  !$omp parallel
                  !$omp do schedule(runtime) private(k, iTracer)
#endif
                  do iCell = 1, nCellsOwned
                  do k = 1, nVertLevels
                  do iTracer = 1, num_tracers
                     activeTracerVertMixTendency(iTracer,k,iCell)=tracersGroup(iTracer,k,iCell)
                  end do
                  end do
                  end do
#ifndef MPAS_OPENACC
                  !$omp end do
                  !$omp end parallel
#endif
               endif
            endif

            if ( associated(tracersGroup) ) then
               if (trim(groupItr % memberName) == 'activeTracers') then
                  call mpas_pool_get_array(tracersSurfaceFluxPool, 'nonLocalSurfaceTracerFlux', &
                           tracerGroupSurfaceFlux)
               else
                  modifiedGroupName = trim(groupItr % memberName) // "SurfaceFlux"
                  call mpas_pool_get_array(tracersSurfaceFluxPool, trim(modifiedGroupName), &
                          tracerGroupSurfaceFlux)
               endif
#ifdef MPAS_OPENACC
               !$acc enter data copyin(tracerGroupSurfaceFlux)
#endif

               call ocn_tracer_vmix_tend_implicit(dt, vertDiffTopOfCell, layerThickness, tracersGroup, &
                        vertNonLocalFlux, tracerGroupSurfaceFlux,  &
                        config_cvmix_kpp_nonlocal_with_implicit_mix, err)
#ifdef MPAS_OPENACC
               !$acc exit data delete(tracerGroupSurfaceFlux)
#endif
            end if

            ! difference tracers to compute influence of vertical mixing and divide by dt
            if (trim(groupItr % memberName) == 'activeTracers') then
               if (config_compute_active_tracer_budgets) then
#ifdef MPAS_OPENACC
                  !$acc parallel loop gang vector collapse(3) present(activeTracerVertMixTendency, tracersGroup)
#else
                  !$omp parallel
                  !$omp do schedule(runtime) private(k, iTracer)
#endif
                  do iCell = 1, nCellsOwned
                  do k = 1, nVertLevels
                  do iTracer = 1, num_tracers
                     activeTracerVertMixTendency(iTracer,k,iCell) = &
                        (tracersGroup(iTracer,k,iCell) - activeTracerVertMixTendency(iTracer,k,iCell)) / dt
                  end do
                  end do
                  end do
#ifndef MPAS_OPENACC
                  !$omp end do
                  !$omp end parallel
#endif

#ifdef MPAS_OPENACC
                  !$acc update host(activeTracerVertMixTendency)
#endif
               endif
            endif

#ifdef MPAS_OPENACC
      !$acc exit data copyout(tracersGroup)
#endif
         end if
      end do
      call mpas_timer_stop('vmix solve tracers')

#ifdef MPAS_OPENACC
      !$acc exit data delete(layerThickness)
#endif

      call mpas_timer_stop('vmix imp')

   end subroutine ocn_vmix_implicit!}}}

!***********************************************************************
!
!  routine ocn_vmix_init
!
!> \brief   Initializes ocean vertical mixing quantities
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_init(domain, err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      integer :: err_tmp

      err = 0

      rayleighDampingCoef = 0.0_RKIND

      if (config_Rayleigh_friction) then
          rayleighDampingCoef = config_Rayleigh_damping_coeff
      endif

      rayleighBottomDampingCoef = 0.0_RKIND

      if (config_Rayleigh_bottom_friction) then
          rayleighBottomDampingCoef = config_Rayleigh_bottom_damping_coeff
      endif

      rayleighDepthVariable = 0.0_RKIND
      if (config_Rayleigh_damping_depth_variable) then
        rayleighDepthVariable = 1.0_RKIND
      end if

      velVmixOn = .true.
      tracerVmixOn = .true.

      if(config_disable_vel_vmix.or.config_disable_vel_all_tend) velVmixOn = .false.
      if(config_disable_tr_vmix.or.config_disable_tr_all_tend) tracerVmixOn = .false.

      implicitBottomDragCoef = 0.0_RKIND

      if (config_use_implicit_bottom_drag) then
          implicitBottomDragCoef = config_implicit_bottom_drag_coeff
      endif

      call ocn_vmix_cvmix_init(domain,err_tmp)
      err = ior(err, err_tmp)
      call ocn_vmix_gotm_init(domain,err_tmp)
      err = ior(err, err_tmp)
      call ocn_vmix_coefs_redi_init(err_tmp)
      err = ior(err, err_tmp)

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_init!}}}
!***********************************************************************

end module ocn_vmix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
