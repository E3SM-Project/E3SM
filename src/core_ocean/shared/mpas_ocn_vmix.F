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

   use ocn_constants
   use ocn_vmix_coefs_const
   use ocn_vmix_coefs_tanh
   use ocn_vmix_coefs_rich
   use ocn_vmix_cvmix
   use ocn_vmix_coefs_redi

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

   subroutine ocn_vmix_coefs(meshPool, statePool, forcingPool, diagnosticsPool, scratchPool, err, timeLevelIn)!{{{

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

      type (mpas_pool_type), intent(inout) :: &
         diagnosticsPool             !< Input/Output: diagnostic information

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

      integer :: err1, err2, err3, err4, err5
      integer :: timeLevel

      integer :: iEdge, iCell, nEdges, nCells
      integer, dimension(:), pointer :: nEdgesArray, nCellsArray

      real (kind=RKIND), dimension(:,:), pointer :: vertViscTopOfEdge, vertDiffTopOfCell

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

      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfEdge', vertViscTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)

      nEdges = nEdgesArray( 2 )

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         vertViscTopOfEdge(:, iEdge) = 0.0_RKIND
      end do
      !$omp end do

      nCells = nCellsArray( 2 )

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         vertDiffTopOfCell(:, iCell) = 0.0_RKIND
      end do
      !$omp end do

      call ocn_vmix_coefs_const_build(meshPool, statePool, diagnosticsPool, err1, timeLevel)
      call ocn_vmix_coefs_tanh_build(meshPool, statePool, diagnosticsPool, err2, timeLevel)
      call ocn_vmix_coefs_rich_build(meshPool, statePool, diagnosticsPool, scratchPool, err3, timeLevel)
      call ocn_vmix_coefs_cvmix_build(meshPool, statePool, forcingPool, diagnosticsPool, err4, timeLevel)
      call ocn_vmix_coefs_redi_build(meshPool, statePool, diagnosticsPool, err5, timeLevel)

      err = ior(ior(ior(err1, ior(err2, err3)), err4), err5)

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs!}}}

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

      integer :: iEdge, k, cell1, cell2, N, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nEdgesArray

      integer, dimension(:), pointer :: maxLevelEdgeTop

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), allocatable :: A, B, C, velTemp

      err = 0

      if(.not.velVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      nEdges = nEdgesArray( 1 )

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))
      A(1)=0

      !$omp do schedule(runtime) private(N, cell1, cell2, k)
      do iEdge = 1, nEdges
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         ! Compute A(k), B(k), C(k)
         ! layerThicknessEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThicknessEdge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, N
            layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do

         ! A is lower diagonal term
         do k = 2, N
            A(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
               / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! C is upper diagonal term
         do k = 1, N-1
            C(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
               / (layerThicknessEdge(k,iEdge) + layerThicknessEdge(k+1,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! B is diagonal term
         B(1) = 1 - C(1)
         do k = 2, N-1
            B(k) = 1 - A(k) - C(k)
         enddo

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
         B(N) = 1 - A(N) + dt*implicitBottomDragCoef &
              * sqrt(kineticEnergyCell(k,cell1) + kineticEnergyCell(k,cell2)) / layerThicknessEdge(k,iEdge)

         call tridiagonal_solve(A(2:N),B,C(1:N-1),normalVelocity(:,iEdge),velTemp,N)

         normalVelocity(1:N,iEdge) = velTemp(1:N)
         normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

        end if
      end do
      !$omp end do

      deallocate(A,B,C,velTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit!}}}

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

      integer :: iCell, k, num_tracers, N, nCells
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), dimension(:), allocatable :: A,B,C
      real (kind=RKIND), dimension(:,:), allocatable :: tracersTemp, rhs

      err = 0

      if(.not.tracerVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      num_tracers = size(tracers, dim=1)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),tracersTemp(num_tracers,nVertLevels))
      allocate(rhs(num_tracers,nVertLevels))

      nCells = nCellsArray( 1 )

      call mpas_timer_start('vmix tracers tend imp loop', .false.)
      !$omp do schedule(runtime) private(N, k)
      do iCell = 1, nCells
         ! Compute A(k), B(k), C(k) for tracers
         N = maxLevelCell(iCell)

         ! A is lower diagonal term
         A(1)=0
         do k = 2, N
            A(k) = -2.0_RKIND*dt*vertDiffTopOfCell(k,iCell) &
                 / (layerThickness(k-1,iCell) + layerThickness(k,iCell)) / layerThickness(k,iCell)
         enddo

         ! C is upper diagonal term
         do k = 1, N-1
            C(k) = -2.0_RKIND*dt*vertDiffTopOfCell(k+1,iCell) &
                 / (layerThickness(k,iCell) + layerThickness(k+1,iCell)) / layerThickness(k,iCell)
         enddo
         C(N) = 0.0_RKIND

         ! B is diagonal term
         do k = 1, N
            B(k) = 1 - A(k) - C(k)
         enddo

         if ( config_cvmix_kpp_nonlocal_with_implicit_mix ) then
            call ocn_compute_kpp_rhs(tracers(:,:,iCell), rhs(:,:), dt, N, num_tracers, &
                             layerThickness(:,iCell), vertNonLocalFlux(:,:,iCell), &
                             tracerGroupSurfaceFlux(:,iCell))
         else
            rhs(:,:) = tracers(:,:,iCell)
         endif

         call tridiagonal_solve_mult(A(2:N), B, C(1:N-1), rhs(:,:), &
              tracersTemp, N, nVertLevels, num_tracers)

         tracers(:,1:N,iCell) = tracersTemp(:,1:N)
         tracers(:,N+1:nVertLevels,iCell) = -1e34
      end do
      !$omp end do
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

   subroutine ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, forcingPool, scratchPool, err, timeLevelIn)!{{{
      real (kind=RKIND), intent(in) :: dt
      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(inout) :: diagnosticsPool
      type (mpas_pool_type), intent(inout) :: statePool
      type (mpas_pool_type), intent(inout) :: forcingPool
      type (mpas_pool_type), intent(in) :: scratchPool !< Input/Output: Scratch structure
      integer, intent(out) :: err
      integer, intent(in), optional :: timeLevelIn

      type (mpas_pool_type), pointer :: tracersPool, tracersSurfaceFluxPool

      integer :: iCell, timeLevel, k, cell1, cell2, iEdge, nCells, nEdges
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, layerThickness, layerThicknessEdge, vertViscTopOfEdge, &
                                                    vertDiffTopOfCell, kineticEnergyCell
      real (kind=RKIND), dimension(:,:), pointer :: vertViscTopOfCell, nonLocalSurfaceTracerFlux, tracerGroupSurfaceFlux
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup, vertNonLocalFlux
      real (kind=RKIND), dimension(:,:,:), allocatable :: nonLocalFluxTend
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracerVertMixTendency
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge
      logical, pointer :: config_use_cvmix,config_compute_active_tracer_budgets
      logical, pointer :: config_cvmix_kpp_nonlocal_with_implicit_mix

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
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix', config_use_cvmix)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_nonlocal_with_implicit_mix', &
                                config_cvmix_kpp_nonlocal_with_implicit_mix)
      call mpas_pool_get_config(ocnConfigs, 'config_compute_active_tracer_budgets',  &
                                config_compute_active_tracer_budgets)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfEdge', vertViscTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfCell', vertViscTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, &
         'activeTracerVertMixTendency',activeTracerVertMixTendency)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

      call mpas_timer_start('vmix coefs', .false.)
      call ocn_vmix_coefs(meshPool, statePool, forcingPool, diagnosticsPool, scratchPool, err, timeLevel)
      call mpas_timer_stop('vmix coefs')

      nCells = nCellsArray(1)
      call mpas_pool_get_array(diagnosticsPool, 'vertNonLocalFlux', vertNonLocalFlux)
      ! if using CVMix, then viscosity has to be averaged from cell centers to cell edges
      if ( config_use_cvmix ) then

         nEdges = nEdgesArray( 1 )
         call mpas_timer_start('CVMix avg', .false.)
         !$omp do schedule(runtime) private(cell1, cell2, k)
         do iEdge=1,nEdges
            vertViscTopOfEdge(:, iEdge) = 0.0_RKIND
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,maxLevelEdgeTop(iEdge)
               vertViscTopOfEdge(k,iEdge) = 0.5_RKIND*(vertViscTopOfCell(k,cell2)+vertViscTopOfCell(k,cell1))
            end do
         end do
         !$omp end do
         call mpas_timer_stop('CVMix avg')
      endif

      !
      !  Implicit vertical solve for momentum
      !
      call mpas_timer_start('vmix solve momentum', .false.)
      call ocn_vel_vmix_tend_implicit(meshPool, dt, kineticEnergyCell, vertViscTopOfEdge, layerThickness, layerThicknessEdge, &
                                      normalVelocity, err)
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
                  !$omp do schedule(runtime)
                  do iCell = 1, nCells
                     activeTracerVertMixTendency(:,:,iCell)=tracersGroup(:,:,iCell)
                  end do
                  !$omp end do
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
                  !$omp do schedule(runtime)
                  do iCell = 1, nCells
                     activeTracerVertMixTendency(:,:,iCell) = &
                        (tracersGroup(:,:,iCell) - activeTracerVertMixTendency(:,:,iCell)) / dt
                  end do
                  !$omp end do
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
      logical, pointer :: config_disable_vel_vmix, config_disable_tr_vmix
      logical, pointer :: config_use_implicit_bottom_drag
      real (kind=RKIND), pointer :: config_implicit_bottom_drag_coeff

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_disable_vel_vmix', config_disable_vel_vmix)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_tr_vmix', config_disable_tr_vmix)
      call mpas_pool_get_config(ocnConfigs, 'config_use_implicit_bottom_drag', config_use_implicit_bottom_drag)
      call mpas_pool_get_config(ocnConfigs, 'config_implicit_bottom_drag_coeff', config_implicit_bottom_drag_coeff)

      velVmixOn = .true.
      tracerVmixOn = .true.

      if(config_disable_vel_vmix) velVmixOn = .false.
      if(config_disable_tr_vmix) tracerVmixOn = .false.

      implicitBottomDragCoef = 0.0_RKIND

      if (config_use_implicit_bottom_drag) then
          implicitBottomDragCoef = config_implicit_bottom_drag_coeff
      endif

      call ocn_vmix_coefs_const_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_vmix_coefs_tanh_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_vmix_coefs_rich_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_vmix_cvmix_init(domain,err_tmp)
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
subroutine tridiagonal_solve_mult(a,b,c,r,x,n,nDim,nSystems)!{{{

   integer,intent(in) :: n, nDim, nSystems
   real (KIND=RKIND), dimension(n), intent(in) :: a,b,c
   real (KIND=RKIND), dimension(nSystems,nDim), intent(in) :: r
   real (KIND=RKIND), dimension(nSystems,nDim), intent(out) :: x
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

subroutine ocn_compute_kpp_rhs(tracers, rhs, dt, maxLevelCell, nTracers, &
                         layerThickness, vertNonLocalFlux, tracerGroupSurfaceFlux)!{{{

  real (kind=RKIND), intent(in) :: dt
  real (kind=RKIND), dimension(:,:), intent(in) :: vertNonLocalFlux, tracers
  real (kind=RKIND), dimension(:), intent(in) :: layerThickness, tracerGroupSurfaceFlux
  real (kind=RKIND), dimension(:,:), intent(inout) :: rhs
  integer, intent(in) :: maxLevelCell, nTracers
  integer :: iTracer, k

  do k=2,maxLevelCell-1
     do iTracer=1,nTracers
        rhs(iTracer, k) = tracers(iTracer,k) + dt * tracerGroupSurfaceFlux(iTracer) *   &
                     (vertNonLocalFlux(1,k) - vertNonLocalFlux(1,k+1)) / layerThickness(k)
     enddo
  enddo

  k=1
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
