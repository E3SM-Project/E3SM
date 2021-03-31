! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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

   private :: tridiagonal_solve, &
              tridiagonal_solve_mult

   public :: ocn_vmix_coefs, &
             ocn_vel_vmix_tend_implicit, &
             ocn_tracer_vmix_tend_implicit, &
             ocn_vmix_init, &
             ocn_vmix_implicit, &
             ocn_compute_kpp_rhs

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

      integer :: iEdge, iCell, nEdges, nCells
      integer, dimension(:), pointer :: nEdgesArray, nCellsArray

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

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

      nEdges = nEdgesArray( 2 )

      !$omp parallel
      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         vertViscTopOfEdge(:, iEdge) = 0.0_RKIND
      end do
      !$omp end do
      !$omp end parallel

      nCells = nCellsArray( 2 )

      !$omp parallel
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         vertDiffTopOfCell(:, iCell) = 0.0_RKIND
      end do
      !$omp end do
      !$omp end parallel

      call ocn_vmix_coefs_cvmix_build(meshPool, statePool, forcingPool, err1, timeLevel)
      call ocn_vmix_coefs_redi_build(meshPool, statePool, err2, timeLevel)
      call ocn_vmix_coefs_gotm_build(statePool, forcingPool, err3, timeLevel)

      err = ior(ior(err1, err2), err3)

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

   subroutine ocn_vel_vmix_tend_implicit_rayleigh(meshPool, dt, kineticEnergyCell, vertViscTopOfEdge, layerThickness, & !{{{
                                         layerThicknessEdge, normalVelocity, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThicknessEdge        !< Input: thickness at edge

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

      integer :: iEdge, k, cell1, cell2, Nsurf, N, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nEdgesArray

      integer, dimension(:), pointer :: minLevelEdgeBot, maxLevelEdgeTop

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), allocatable :: A, B, C, velTemp
      real (kind=RKIND) :: edgeThicknessTotal

      err = 0

      if(.not.velVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      nEdges = nEdgesArray( 1 )

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(Nsurf, N, cell1, cell2, edgeThicknessTotal, k, A, B, C, velTemp)
      do iEdge = 1, nEdges
        Nsurf = minLevelEdgeBot(iEdge)
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         ! Compute A(k), B(k), C(k)
         ! layerThicknessEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThicknessEdge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         edgeThicknessTotal = 0.0_RKIND
         do k = Nsurf, N
            layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
            edgeThicknessTotal = edgeThicknessTotal + layerThicknessEdge(k,iEdge)
         end do

         ! A is lower diagonal term
         A(1:Nsurf)=0.0_RKIND
         do k = Nsurf+1, N
            A(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
               / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! C is upper diagonal term
         C(1:Nsurf-1)=0.0_RKIND
         do k = Nsurf, N-1
            C(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
               / (layerThicknessEdge(k,iEdge) + layerThicknessEdge(k+1,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! B is diagonal term
         B(1:Nsurf-1)=0.0_RKIND
         B(Nsurf) = 1.0_RKIND - C(Nsurf) &
           + dt*rayleighDampingCoef/((1.0_RKIND - rayleighDepthVariable) + rayleighDepthVariable*edgeThicknessTotal)
         do k = Nsurf+1, N-1
            B(k) = 1.0_RKIND - A(k) - C(k) &
              ! added Rayleigh terms
              + dt*(rayleighDampingCoef/((1.0_RKIND - rayleighDepthVariable) + rayleighDepthVariable*edgeThicknessTotal))
         enddo

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
         B(N) = 1.0_RKIND - A(N) + dt*implicitBottomDragCoef &
              * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThicknessEdge(N,iEdge) &
              ! added Rayleigh terms
              + dt*(rayleighBottomDampingCoef + &
                    rayleighDampingCoef /((1.0_RKIND - rayleighDepthVariable) + rayleighDepthVariable*edgeThicknessTotal))

         call tridiagonal_solve(A(Nsurf+1:N),B(Nsurf:N),C(Nsurf:N-1),normalVelocity(Nsurf:N,iEdge),velTemp(Nsurf:N),N-Nsurf+1)

         normalVelocity(1:Nsurf-1,iEdge) = 0.0_RKIND
         normalVelocity(Nsurf:N,iEdge) = velTemp(Nsurf:N)
         normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

        end if
      end do
      !$omp end do
      !$omp end parallel

      deallocate(A,B,C,velTemp)

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

   subroutine ocn_vel_vmix_tend_implicit(meshPool, dt, kineticEnergyCell, vertViscTopOfEdge, layerThickness, & !{{{
                                         layerThicknessEdge, normalVelocity, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThicknessEdge        !< Input: thickness at edge

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

      integer :: iEdge, k, cell1, cell2, N, Nsurf, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nEdgesArray

      integer, dimension(:), pointer :: maxLevelEdgeTop, minLevelEdgeBot

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), allocatable :: A, B, C, velTemp

      err = 0

      if(.not.velVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      nEdges = nEdgesArray( 1 )

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))

      !$omp parallel
      !$omp do schedule(runtime) private(Nsurf, N, cell1, cell2, A, B, C, velTemp, k)
      do iEdge = 1, nEdges
        N = maxLevelEdgeTop(iEdge)
        Nsurf = minLevelEdgeBot(iEdge)
        if (N .gt. 0) then

         ! Compute A(k), B(k), C(k)
         ! layerThicknessEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThicknessEdge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = Nsurf, N
            layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do

         ! A is lower diagonal term
         A(1:Nsurf)=0.0_RKIND
         do k = Nsurf+1, N
            A(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
               / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! C is upper diagonal term
         C(1:Nsurf-1)=0.0_RKIND
         do k = Nsurf, N-1
            C(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
               / (layerThicknessEdge(k,iEdge) + layerThicknessEdge(k+1,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! B is diagonal term
         B(1:Nsurf-1)=0.0_RKIND
         B(Nsurf) = 1.0_RKIND - C(Nsurf)
         do k = Nsurf+1, N-1
            B(k) = 1.0_RKIND - A(k) - C(k)
         enddo

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
         B(N) = 1.0_RKIND - A(N) + dt*implicitBottomDragCoef &
              * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThicknessEdge(N,iEdge)

         call tridiagonal_solve(A(Nsurf+1:N),B(Nsurf:N),C(Nsurf:N-1),normalVelocity(Nsurf:N,iEdge),velTemp(Nsurf:N),N-Nsurf+1)

         normalVelocity(1:Nsurf-1,iEdge) = 0.0_RKIND
         normalVelocity(Nsurf:N,iEdge) = velTemp(Nsurf:N)
         normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

        end if
      end do
      !$omp end do
      !$omp end parallel

      deallocate(A,B,C,velTemp)

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

   subroutine ocn_vel_vmix_tend_implicit_spatially_variable(meshPool, bottomDrag, dt, kineticEnergyCell, & !{{{
                                         vertViscTopOfEdge, layerThickness, &
                                         layerThicknessEdge, normalVelocity, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

       real (kind=RKIND), dimension(:), intent(in) :: &
         bottomDrag !< Input: bottomDrag at cell centeres

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThicknessEdge        !< Input: thickness at edge

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

      integer :: iEdge, k, cell1, cell2, N, Nsurf, nEdges
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: implicitCd
      integer, dimension(:), pointer :: nEdgesArray

      integer, dimension(:), pointer :: minLevelEdgeBot, maxLevelEdgeTop

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), allocatable :: A, B, C, velTemp

      err = 0

      if(.not.velVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      nEdges = nEdgesArray( 1 )

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
        Nsurf = minLevelEdgeBot(iEdge)
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         ! Compute A(k), B(k), C(k)
         ! layerThicknessEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThicknessEdge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = Nsurf, N
            layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do

         ! average cell-based implicit bottom drag to edges
         implicitCd = 0.5_RKIND*(bottomDrag(cell1) + bottomDrag(cell2))

         ! A is lower diagonal term
         A(1:Nsurf)=0.0_RKIND
         do k = Nsurf+1, N
            A(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
               / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! C is upper diagonal term
         C(1:Nsurf-1)=0.0_RKIND
         do k = Nsurf, N-1
            C(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
               / (layerThicknessEdge(k,iEdge) + layerThicknessEdge(k+1,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! B is diagonal term
         B(1:Nsurf-1)=0.0_RKIND
         B(Nsurf) = 1.0_RKIND - C(Nsurf)
         do k = Nsurf+1, N-1
            B(k) = 1.0_RKIND - A(k) - C(k)
         enddo

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
         ! use implicitCd from spatially variable bottom drag
         B(N) = 1.0_RKIND - A(N) + dt*implicitCd &
              * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThicknessEdge(N,iEdge)

         call tridiagonal_solve(A(Nsurf+1:N),B(Nsurf:N),C(Nsurf:N-1),normalVelocity(Nsurf:N,iEdge),velTemp(Nsurf:N),N-Nsurf+1)

         normalVelocity(1:Nsurf-1,iEdge) = 0.0_RKIND
         normalVelocity(Nsurf:N,iEdge) = velTemp(Nsurf:N)
         normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

        end if
      end do
      !$omp end do

      deallocate(A,B,C,velTemp)

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

   subroutine ocn_vel_vmix_tend_implicit_spatially_variable_mannings(meshPool, forcingPool, bottomDrag, dt, & !{{{
                                         kineticEnergyCell, vertViscTopOfEdge, layerThickness, &
                                         layerThicknessEdge, normalVelocity, ssh, bottomDepth, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

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

       real (kind=RKIND), dimension(:), intent(in) :: &
         bottomDrag !< Input: bottomDrag at cell centeres

       real (kind=RKIND), dimension(:), pointer :: ssh

       real (kind=RKIND), dimension(:), pointer :: bottomDepth

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThicknessEdge        !< Input: thickness at edge

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

      integer :: iEdge, k, cell1, cell2, N, Nsurf, nEdges
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: implicitCd
      integer, dimension(:), pointer :: nEdgesArray

      integer, dimension(:), pointer :: minLevelEdgeBot, maxLevelEdgeTop

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), allocatable :: A, B, C, velTemp

      ! vegetation_drag
      real (kind=RKIND), dimension(:), pointer :: vegetationHeight
      real (kind=RKIND), dimension(:), pointer :: vegetationDiameter
      real (kind=RKIND), dimension(:), pointer ::vegetationDensity
      real (kind=RKIND), dimension(:), pointer ::vegetationManning
      integer, dimension(:), pointer ::vegetationMask
      real (kind=RKIND) :: old_bottom_Cd, lambda, beta, alpha, total_h
      real (kind=RKIND) :: inundation_depth, von_karman, cff1, cff2, cff3, cff4
      integer :: iCell, nCells
      integer, pointer :: nCellsSolve

      err = 0
      von_karman = 0.4_RKIND

      if(.not.velVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(forcingPool, 'vegetationMask', vegetationMask)
      call mpas_pool_get_array(forcingPool, 'vegetationHeight', vegetationHeight)
      call mpas_pool_get_array(forcingPool, 'vegetationDensity', vegetationDensity)
      call mpas_pool_get_array(forcingPool, 'vegetationDiameter', vegetationDiameter)
      call mpas_pool_get_array(forcingPool, 'vegetationManning', vegetationManning)

      nEdges = nEdgesArray( 1 )
      nCells = nCellsSolve

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))

      ! Compute bottomDrag (Manning roughness) induced by vegetation
      if (config_use_vegetation_drag .AND. config_use_vegetation_manning_equation) then
        do iCell = 1, nCells
          if (vegetationDensity(iCell) * vegetationHeight(iCell) * vegetationDiameter(iCell) .eq. 0.0_RKIND) then
            vegetationMask(iCell) = 0
          endif
          if (vegetationMask(iCell) .eq. 1) then
            total_h = bottomDepth(iCell) + ssh(iCell)
            old_bottom_Cd = gravity * bottomDrag(iCell)**2.0_RKIND * total_h**(1.0_RKIND/3.0_RKIND)
            inundation_depth = MIN(vegetationHeight(iCell), total_h)
            inundation_depth = MAX(inundation_depth, 1e-6)
            lambda = vegetationDiameter(iCell) * vegetationDensity(iCell)
            alpha = (config_vegetation_drag_coefficient*lambda/ &
                   (4.0_RKIND*von_karman**2.0_RKIND*inundation_depth**2.0_RKIND))**(1.0_RKIND/3.0_RKIND)
            beta = 0.5_RKIND*alpha*inundation_depth*(1.0_RKIND - EXP(-2.0_RKIND*alpha*inundation_depth)) &
                  / (1.0_RKIND - EXP(-alpha*inundation_depth))**2.0_RKIND
            cff1 = total_h**(2.0_RKIND/3.0_RKIND) &
                  * SQRT((0.5_RKIND*beta*lambda*config_vegetation_drag_coefficient*inundation_depth &
                  + old_bottom_Cd)/(gravity*total_h))
            cff2 = (alpha*inundation_depth)**2.0_RKIND/(1.0_RKIND - EXP(-alpha*inundation_depth))
            cff3 = (1.0_RKIND - EXP(-alpha*inundation_depth))/(alpha**2.0_RKIND*inundation_depth*total_h)
            cff4 = LOG(total_h/inundation_depth) - (1.0_RKIND-inundation_depth/total_h) &
                  * (1.0_RKIND - 1.0_RKIND/(alpha*inundation_depth))
            vegetationManning(iCell) = cff1/(cff2*(cff3+cff4))
            vegetationManning(iCell) = MAX(bottomDrag(iCell), vegetationManning(iCell))
          else
            vegetationManning(iCell) = bottomDrag(iCell)
          endif
        enddo
      endif

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
        Nsurf = minLevelEdgeBot(iEdge)
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         ! Compute A(k), B(k), C(k)
         ! layerThicknessEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThicknessEdge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = Nsurf, N
            layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do

         ! average cell-based implicit bottom drag to edges and convert Mannings n to Cd
         if (config_use_vegetation_drag .AND. config_use_vegetation_manning_equation) then
           implicitCd = gravity*(0.5_RKIND*(vegetationManning(cell1) + vegetationManning(cell2)))**2.0 * &
            (0.5_RKIND * (ssh(cell1) + ssh(cell2) + bottomDepth(cell1) + bottomDepth(cell2)))**(-1.0_RKIND/3.0_RKIND)
         else
           implicitCd = gravity*(0.5_RKIND*(bottomDrag(cell1) + bottomDrag(cell2)))**2.0 * &
            (0.5_RKIND * (ssh(cell1) + ssh(cell2) + bottomDepth(cell1) + bottomDepth(cell2)))**(-1.0_RKIND/3.0_RKIND)
         endif

         ! A is lower diagonal term
         A(1:Nsurf)=0.0_RKIND
         do k = Nsurf+1, N
            A(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
               / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! C is upper diagonal term
         C(1:Nsurf-1)=0.0_RKIND
         do k = Nsurf, N-1
            C(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
               / (layerThicknessEdge(k,iEdge) + layerThicknessEdge(k+1,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! B is diagonal term
         B(1:Nsurf-1)=0.0_RKIND
         B(Nsurf) = 1.0_RKIND - C(Nsurf)
         do k = Nsurf+1, N-1
            B(k) = 1.0_RKIND - A(k) - C(k)
         enddo

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
         ! use implicitCd from spatially variable bottom drag
         B(N) = 1.0_RKIND - A(N) + dt*implicitCd &
              * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThicknessEdge(N,iEdge)

         call tridiagonal_solve(A(Nsurf+1:N),B(Nsurf:N),C(Nsurf:N-1),normalVelocity(Nsurf:N,iEdge),velTemp(Nsurf:N),N-Nsurf+1)

         normalVelocity(1:Nsurf-1,iEdge) = 0.0_RKIND
         normalVelocity(Nsurf:N,iEdge) = velTemp(Nsurf:N)
         normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

        end if
      end do
      !$omp end do

      deallocate(A,B,C,velTemp)

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

   subroutine ocn_tracer_vmix_tend_implicit(meshPool, dt, vertDiffTopOfCell, layerThickness, tracers, &
                  vertNonLocalFlux, tracerGroupSurfaceFlux, config_cvmix_kpp_nonlocal_with_implicit_mix, &
                  err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

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

      integer :: iCell, k, num_tracers, N, Nsurf, nCells
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray

      integer, dimension(:), pointer :: maxLevelCell, minLevelCell

      real (kind=RKIND), dimension(:), allocatable :: A,B,C
      real (kind=RKIND), dimension(:,:), allocatable :: tracersTemp, rhs

      err = 0

      if(.not.tracerVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      num_tracers = size(tracers, dim=1)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),tracersTemp(num_tracers,nVertLevels))
      allocate(rhs(num_tracers,nVertLevels))

      nCells = nCellsArray( 1 )

      call mpas_timer_start('vmix tracers tend imp loop', .false.)
      !$omp parallel
      !$omp do schedule(runtime) private(Nsurf, N, A, B, C, rhs, tracersTemp, k)
      do iCell = 1, nCells
         ! Compute A(k), B(k), C(k) for tracers
         Nsurf = minLevelCell(iCell)
         N = maxLevelCell(iCell)

         ! A is lower diagonal term
         A(1:Nsurf)=0.0_RKIND
         do k = Nsurf+1, N
            A(k) = -2.0_RKIND*dt*vertDiffTopOfCell(k,iCell) &
                 / (layerThickness(k-1,iCell) + layerThickness(k,iCell)) / layerThickness(k,iCell)
         enddo

         ! C is upper diagonal term
         C(1:Nsurf-1)=0.0_RKIND
         do k = Nsurf, N-1
            C(k) = -2.0_RKIND*dt*vertDiffTopOfCell(k+1,iCell) &
                 / (layerThickness(k,iCell) + layerThickness(k+1,iCell)) / layerThickness(k,iCell)
         enddo
         C(N) = 0.0_RKIND

         ! B is diagonal term
         B(1:Nsurf-1)=0.0_RKIND
         do k = Nsurf, N
            B(k) = 1.0_RKIND - A(k) - C(k)
         enddo

         if ( config_cvmix_kpp_nonlocal_with_implicit_mix ) then
            call ocn_compute_kpp_rhs(tracers(:,:,iCell), rhs(:,:), dt, Nsurf, N, num_tracers, &
                             layerThickness(:,iCell), vertNonLocalFlux(:,:,iCell), &
                             tracerGroupSurfaceFlux(:,iCell))
         else
            rhs(:,:) = tracers(:,:,iCell)
         endif

         call tridiagonal_solve_mult(A(Nsurf+1:N), B(Nsurf:N), C(Nsurf:N-1), rhs(:,Nsurf:N), &
              tracersTemp(:,Nsurf:N), N-Nsurf+1, num_tracers)

         tracers(:,Nsurf:N,iCell) = tracersTemp(:,Nsurf:N)
         tracers(:,N+1:nVertLevels,iCell) = -1e34
         tracers(:,1:Nsurf-1,iCell) = -1e34
      end do
      !$omp end do
      !$omp end parallel
      call mpas_timer_stop('vmix tracers tend imp loop')

      deallocate(A, B, C, tracersTemp, rhs)

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

      integer :: iCell, timeLevel, k, cell1, cell2, iEdge, nCells, nEdges
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray
      real (kind=RKIND), dimension(:), pointer :: bottomDrag, ssh, bottomDepth   ! needed for depth-variable computation
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: nonLocalSurfaceTracerFlux, tracerGroupSurfaceFlux
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup
      real (kind=RKIND), dimension(:,:,:), allocatable :: nonLocalFluxTend
      integer, dimension(:), pointer :: maxLevelCell, minLevelCell, maxLevelEdgeTop, minLevelEdgeTop, minLevelEdgeBot
      integer, dimension(:,:), pointer :: cellsOnEdge

      type (mpas_pool_iterator_type) :: groupItr

      character (len=StrKIND) :: modifiedGroupName
      integer, pointer :: indexTempFlux, indexSaltFlux, nVertLevels
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

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'minLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeTop', minLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

      call mpas_pool_get_array(forcingPool, 'bottomDrag', bottomDrag)

      call mpas_timer_start('vmix coefs', .false.)
      call ocn_vmix_coefs(meshPool, statePool, forcingPool, scratchPool, err, timeLevel)
      call mpas_timer_stop('vmix coefs')

      nCells = nCellsArray(1)
      ! if using CVMix, then viscosity has to be averaged from cell centers to cell edges
      if ( config_use_cvmix ) then

         nEdges = nEdgesArray( 1 )
         call mpas_timer_start('CVMix avg', .false.)
         !$omp parallel
         !$omp do schedule(runtime) private(cell1, cell2, k)
         do iEdge=1,nEdges
            vertViscTopOfEdge(:, iEdge) = 0.0_RKIND
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
               vertViscTopOfEdge(k,iEdge) = 0.5_RKIND*(vertViscTopOfCell(k,cell2)+vertViscTopOfCell(k,cell1))
            end do
         end do
         !$omp end do
         !$omp end parallel
         call mpas_timer_stop('CVMix avg')
      endif

      ! if using GOTM, then viscosity has to be averaged from cell centers to cell edges
      if ( config_use_gotm ) then

         nEdges = nEdgesArray( 1 )
         call mpas_timer_start('GOTM avg', .false.)
         !$omp do schedule(runtime) private(cell1, cell2, k)
         do iEdge=1,nEdges
            vertViscTopOfEdge(:, iEdge) = 0.0_RKIND
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
               vertViscTopOfEdge(k,iEdge) = 0.5_RKIND*(vertViscTopOfCell(k,cell2)+vertViscTopOfCell(k,cell1))
            end do
         end do
         !$omp end do
         call mpas_timer_stop('GOTM avg')
      endif

      !
      !  Implicit vertical solve for momentum
      !
      call mpas_timer_start('vmix solve momentum', .false.)
      if (config_use_implicit_bottom_drag_variable) then
        call ocn_vel_vmix_tend_implicit_spatially_variable(meshPool, bottomDrag, dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdge, normalVelocity, err)
      else if (config_use_implicit_bottom_drag_variable_mannings) then
        ! update bottomDrag via Cd=g*n^2*h^(-1/3)
        call ocn_vel_vmix_tend_implicit_spatially_variable_mannings(meshPool, forcingPool, bottomDrag, &
          dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdge, normalVelocity, &
          ssh, bottomDepth, err)
      else if (config_Rayleigh_friction.or. &
               config_Rayleigh_bottom_friction.or. &
               config_Rayleigh_damping_depth_variable) then
        call ocn_vel_vmix_tend_implicit_rayleigh(meshPool, dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdge, normalVelocity, err)
      else
        call ocn_vel_vmix_tend_implicit(meshPool, dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdge, normalVelocity, err)
      end if
      call mpas_timer_stop('vmix solve momentum')

      !
      !  Implicit vertical solve for all tracers
      !

      call mpas_timer_start('vmix solve tracers', .false.)
      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )

         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroup, timeLevel)
            ! store tracers
            if (trim(groupItr % memberName) == 'activeTracers') then
               if (config_compute_active_tracer_budgets) then
                  !$omp parallel
                  !$omp do schedule(runtime)
                  do iCell = 1, nCells
                     activeTracerVertMixTendency(:,:,iCell)=tracersGroup(:,:,iCell)
                  end do
                  !$omp end do
                  !$omp end parallel
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

               call ocn_tracer_vmix_tend_implicit(meshPool, dt, vertDiffTopOfCell, layerThickness, tracersGroup, &
                        vertNonLocalFlux, tracerGroupSurfaceFlux,  &
                        config_cvmix_kpp_nonlocal_with_implicit_mix, err)
            end if

            ! difference tracers to compute influence of vertical mixing and divide by dt
            if (trim(groupItr % memberName) == 'activeTracers') then
               if (config_compute_active_tracer_budgets) then
                  !$omp parallel
                  !$omp do schedule(runtime)
                  do iCell = 1, nCells
                     activeTracerVertMixTendency(:,:,iCell) = &
                        (tracersGroup(:,:,iCell) - activeTracerVertMixTendency(:,:,iCell)) / dt
                  end do
                  !$omp end do
                  !$omp end parallel
               endif
            endif

         end if
      end do
      call mpas_timer_stop('vmix solve tracers')

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
!
!  routine tridiagonal_solve
!
!> \brief   Solve the matrix equation Ax=r for x, where A is tridiagonal.
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  Solve the matrix equation Ax=r for x, where A is tridiagonal.
!>  A is an nxn matrix, with:
!>  a sub-diagonal, filled from 1:n-1 (a(1) appears on row 2)
!>  b diagonal, filled from 1:n
!>  c sup-diagonal, filled from 1:n-1  (c(1) apears on row 1)
!
!-----------------------------------------------------------------------
   subroutine tridiagonal_solve(a,b,c,r,x,n) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer,intent(in) :: n
      real (KIND=RKIND), dimension(n), intent(in) :: a,b,c,r

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (KIND=RKIND), dimension(n), intent(out) :: x

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (KIND=RKIND), dimension(n) :: bTemp,rTemp
      real (KIND=RKIND) :: m
      integer i

      ! Use work variables for b and r
      bTemp(1) = b(1)
      rTemp(1) = r(1)

      ! First pass: set the coefficients
      do i = 2,n
         m = a(i-1)/bTemp(i-1)
         bTemp(i) = b(i) - m*c(i-1)
         rTemp(i) = r(i) - m*rTemp(i-1)
      end do

      x(n) = rTemp(n)/bTemp(n)
       ! Second pass: back-substition
      do i = n-1, 1, -1
         x(i) = (rTemp(i) - c(i)*x(i+1))/bTemp(i)
      end do

   end subroutine tridiagonal_solve !}}}

!***********************************************************************
!
!  routine tridiagonal_solve_mult
!
!> \brief   Solve multiple matrix equations Ax=r for x, where A is tridiagonal.
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  Solve the matrix equation Ax=r for x, where A is tridiagonal.
!>  A is an nxn matrix, with:
!>  a sub-diagonal, filled from 1:n-1 (a(1) appears on row 2)
!>  b diagonal, filled from 1:n
!>  c sup-diagonal, filled from 1:n-1  (c(1) apears on row 1)
!
!-----------------------------------------------------------------------
subroutine tridiagonal_solve_mult(a,b,c,r,x,n,nSystems)!{{{

   integer,intent(in) :: n, nSystems
   real (KIND=RKIND), dimension(n), intent(in) :: a,b,c
   real (KIND=RKIND), dimension(nSystems,n), intent(in) :: r
   real (KIND=RKIND), dimension(nSystems,n), intent(out) :: x
   real (KIND=RKIND), dimension(n) :: bTemp
   real (KIND=RKIND), dimension(nSystems,n) :: rTemp
   real (KIND=RKIND) :: m
   integer i,j

   ! Use work variables for b and r
   bTemp(1) = b(1)
   do j = 1,nSystems
      rTemp(j,1) = r(j,1)
   end do

   ! First pass: set the coefficients
   do i = 2,n
      m = a(i-1)/bTemp(i-1)
      bTemp(i) = b(i) - m*c(i-1)
      do j = 1,nSystems
         rTemp(j,i) = r(j,i) - m*rTemp(j,i-1)
      end do
   end do

   do j = 1,nSystems
      x(j,n) = rTemp(j,n)/bTemp(n)
   end do
   ! Second pass: back-substition
   do i = n-1, 1, -1
      do j = 1,nSystems
         x(j,i) = (rTemp(j,i) - c(i)*x(j,i+1))/bTemp(i)
      end do
   end do

end subroutine tridiagonal_solve_mult!}}}

!***********************************************************************
!
!  subroutine ocn_compute_kpp_rhs
!
!> \brief   Computes the non local flux tendency for KPP
!> \author  Luke Van Roekel
!> \date    October 2017
!> \details
!>   Computes non local flux tendency from KPP when
!>   config_cvmix_kpp_nonlocal_with_implicit_mix = .true.
!>   otherwise this term is computed in ocn_tend_tracer
!
!-----------------------------------------------------------------------

subroutine ocn_compute_kpp_rhs(tracers, rhs, dt, minLevelCell, maxLevelCell, nTracers, &
                         layerThickness, vertNonLocalFlux, tracerGroupSurfaceFlux)!{{{

  real (kind=RKIND), intent(in) :: dt
  real (kind=RKIND), dimension(:,:), intent(in) :: vertNonLocalFlux, tracers
  real (kind=RKIND), dimension(:), intent(in) :: layerThickness, tracerGroupSurfaceFlux
  real (kind=RKIND), dimension(:,:), intent(out) :: rhs
  integer, intent(in) :: minLevelCell, maxLevelCell, nTracers
  integer :: iTracer, k

  rhs(:,1:minLevelCell-1) = 0.0_RKIND
  
  do k=minLevelCell+1,maxLevelCell-1
     do iTracer=1,nTracers
        rhs(iTracer, k) = tracers(iTracer,k) + dt * tracerGroupSurfaceFlux(iTracer) *   &
                     (vertNonLocalFlux(1,k) - vertNonLocalFlux(1,k+1)) / layerThickness(k)
     enddo
  enddo

  k=minLevelCell
  do iTracer=1,nTracers
     rhs(iTracer, k) = tracers(iTracer,k) + dt * tracerGroupSurfaceFlux(iTracer) *  &
                              (-vertNonLocalFlux(1,k+1) )/ layerThickness(k)
  enddo

  k=maxLevelCell
  do iTracer=1,nTracers
     rhs(iTracer,k) = tracers(iTracer,k) + dt * tracerGroupSurfaceFlux(iTracer) * &
                                     vertNonLocalFlux(1,k) / layerThickness(k)
  enddo

end subroutine ocn_compute_kpp_rhs!}}}
!***********************************************************************

end module ocn_vmix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
