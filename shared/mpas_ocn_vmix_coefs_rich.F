! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_coefs_rich
!
!> \brief MPAS ocean vertical mixing coefficients
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module contains the routines for computing
!>  richardson vertical mixing coefficients.
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_coefs_rich

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timer
   use mpas_threading

   use ocn_constants
   use ocn_equation_of_state

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

   public :: ocn_vmix_coefs_rich_build, &
             ocn_vmix_coefs_rich_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: richViscOn, richDiffOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs_rich_build
!
!> \brief   Computes coefficients for vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine computes the vertical mixing coefficients for momentum
!>  and activeTracers based user choices of mixing parameterization.
!
!-----------------------------------------------------------------------
   subroutine ocn_vmix_coefs_rich_build(meshPool, statePool, diagnosticsPool, scratchPool, err, timeLevelIn)!{{{

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

      integer :: err1, err2, err3, timeLevel, nCells
      integer, pointer :: indexTemperature, indexSalinity, nCellsSolve

      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfEdge, vertDiffTopOfCell, normalVelocity, layerThickness, layerThicknessEdge, density, displacedDensity

      real (kind=RKIND), dimension(:,:), pointer :: RiTopOfEdge, RiTopOfCell

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if ( .not. richViscOn .and. .not. richDiffOn ) then
         return
      end if

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfEdge', vertViscTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'RiTopOfEdge', RiTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'RiTopOfCell', RiTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'displacedDensity', displacedDensity)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)

      nCells = nCellsSolve

      ! compute in-place density
      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, scratchPool, nCells, 0, 'relative', density, err, &
                                         timeLevelIn=timeLevel)

      ! compute displacedDensity, density displaced adiabatically to the mid-depth one layer deeper.
      ! That is, layer k has been displaced to the depth of layer k+1.
      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, scratchPool, nCells, 1, 'relative', &
         displacedDensity, err, timeLevelIn=timeLevel)

      call ocn_vmix_get_rich_numbers(meshPool, scratchPool, indexTemperature, indexSalinity, normalVelocity, layerThickness, &
                                     layerThicknessEdge, density, displacedDensity, activeTracers, RiTopOfEdge, RiTopOfCell, err1)

      call ocn_vel_vmix_coefs_rich(meshPool, RiTopOfEdge, layerThicknessEdge, vertViscTopOfEdge, err2)
      call ocn_tracer_vmix_coefs_rich(meshPool, RiTopOfCell, layerThickness, vertDiffTopOfCell, err3)

      err = ior(err1, ior(err2, err3))

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_rich_build!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_coefs_rich
!
!> \brief   Computes coefficients for vertical momentum mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine computes the richardson vertical mixing coefficients for momentum
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_coefs_rich(meshPool, RiTopOfEdge, layerThicknessEdge, vertViscTopOfEdge, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge        !< Input: thickness at edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         RiTopOfEdge   !< Richardson number at top of edge

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

      real (kind=RKIND), dimension(:,:), intent(out) :: vertViscTopOfEdge !< Output: vertical viscosity

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, k, nEdges
      integer, pointer :: nEdgesSolve

      integer, dimension(:), pointer :: maxLevelEdgeTop

      real (kind=RKIND), pointer :: config_rich_mix, config_bkrd_vert_visc, config_convective_visc

      err = 0

      if(.not.richViscOn) return

      call mpas_timer_start('vel rich coef')

      call mpas_pool_get_config(ocnConfigs, 'config_rich_mix', config_rich_mix)
      call mpas_pool_get_config(ocnConfigs, 'config_bkrd_vert_visc', config_bkrd_vert_visc)
      call mpas_pool_get_config(ocnConfigs, 'config_convective_visc', config_convective_visc)

      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      nEdges = nEdgesSolve

      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 2, maxLevelEdgeTop(iEdge)
            ! efficiency note: these if statements are inside iEdge and k loops.
            ! Perhaps there is a more efficient way to do this.
            if (RiTopOfEdge(k,iEdge)>0.0_RKIND) then
               vertViscTopOfEdge(k,iEdge) = vertViscTopOfEdge(k, iEdge) + config_bkrd_vert_visc &
                  + config_rich_mix / (1.0_RKIND + 5.0_RKIND*RiTopOfEdge(k,iEdge))**2
               if (vertViscTopOfEdge(k,iEdge) > config_convective_visc) then
                  vertViscTopOfEdge(k,iEdge) = config_convective_visc
               end if
            else
               ! for Ri<0 use the convective value for the viscosity
               vertViscTopOfEdge(k,iEdge) = config_convective_visc
            end if
         end do
      end do
      !$omp end do

      call mpas_timer_stop('vel rich coef')


   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_coefs_rich!}}}

!***********************************************************************
!
!  routine ocn_tracer_vmix_coefs_rich
!
!> \brief   Computes coefficients for vertical tracer mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine computes the richardson vertical mixing coefficients for activeTracers
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vmix_coefs_rich(meshPool, RiTopOfCell, layerThickness, vertDiffTopOfCell, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness             !< Input: thickness at cell center

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         RiTopOfCell   !< Input: Richardson number at top of cell

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

      real (kind=RKIND), dimension(:,:), intent(out) :: vertDiffTopOfCell !< Output: vertical diffusions

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, k, nCells
      integer, pointer :: nCellsSolve

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: coef
      real (kind=RKIND), pointer :: config_bkrd_vert_diff, config_bkrd_vert_visc, config_rich_mix, config_convective_diff

      err = 0

      if(.not.richDiffOn) return

      call mpas_timer_start('tracer rich coef')

      call mpas_pool_get_config(ocnConfigs, 'config_bkrd_vert_diff', config_bkrd_vert_diff)
      call mpas_pool_get_config(ocnConfigs, 'config_bkrd_vert_visc', config_bkrd_vert_visc)
      call mpas_pool_get_config(ocnConfigs, 'config_rich_mix', config_rich_mix)
      call mpas_pool_get_config(ocnConfigs, 'config_convective_diff', config_convective_diff)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      nCells = nCellsSolve

      coef = -gravity / rho_sw / 2.0_RKIND
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         do k = 2, maxLevelCell(iCell)
            ! efficiency note: these if statements are inside iEdge and k loops.
            ! Perhaps there is a more efficient way to do this.
            if (RiTopOfCell(k,iCell)>0.0_RKIND) then
               vertDiffTopOfCell(k,iCell) = vertDiffTopOfCell(k, iCell) + config_bkrd_vert_diff &
                  + (config_bkrd_vert_visc &
                     + config_rich_mix / (1.0_RKIND + 5.0_RKIND*RiTopOfCell(k,iCell))**2) &
                  / (1.0_RKIND + 5.0_RKIND*RiTopOfCell(k,iCell))
               if (vertDiffTopOfCell(k,iCell) > config_convective_diff) then
                  vertDiffTopOfCell(k,iCell) = config_convective_diff
               end if
             else
               ! for Ri<0 use the convective value for the diffusion
               vertDiffTopOfCell(k,iCell) = config_convective_diff
            end if
         end do
      end do
      !$omp end do

      call mpas_timer_stop('tracer rich coef')

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vmix_coefs_rich!}}}

!***********************************************************************
!
!  routine ocn_vmix_get_rich_numbers
!
!> \brief   Build richardson numbers for vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine builds the arrays needed for richardson number vertical
!>  mixing coefficients.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_get_rich_numbers(meshPool, scratchPool, indexTemperature, indexSalinity, normalVelocity, & !{{{
                                        layerThickness, layerThicknessEdge, density, displacedDensity, activeTracers, &
                                        RiTopOfEdge, RiTopOfCell, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      type (mpas_pool_type), intent(in) :: scratchPool !< Input: scratch variables

      integer, intent(in) :: indexTemperature !< Input: index for temperature
      integer, intent(in) :: indexSalinity !< Input: index for salinity

      real (kind=RKIND), dimension(:,:), intent(in) :: normalVelocity       !< Input: horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness       !< Input: thickness
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThicknessEdge  !< Input: thickness at edge

      real (kind=RKIND), dimension(:,:,:), intent(in) :: activeTracers !< Input: activeTracers

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

      real (kind=RKIND), dimension(:,:), intent(inout) :: density    !< Input/output: density
      real (kind=RKIND), dimension(:,:), intent(inout) :: displacedDensity    !< Input/output: displaced density
      real (kind=RKIND), dimension(:,:), intent(inout) :: RiTopOfEdge     !< Input/output: Richardson number top of cell
      real (kind=RKIND), dimension(:,:), intent(inout) :: RiTopOfCell     !< Input/output: Richardson number top of cell

      integer, intent(inout) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, iEdge, k, i
      integer :: cell1, cell2, nCells, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOncell, edgeSignOnCell

      real (kind=RKIND) :: coef, invAreaCell
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: ddensityTopOfCell, du2TopOfCell, &
                                                    ddensityTopOfEdge, du2TopOfEdge
      type (field2DReal), pointer :: ddensityTopOfCellField, du2TopOfCellField, &
                                      ddensityTopOfEdgeField, du2TopOfEdgeField

      err = 0

      if ( ( .not. richViscOn ) .and. ( .not. richDiffOn ) ) return

      call mpas_timer_start('get rich nums')

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      call mpas_pool_get_field(scratchPool, 'ddensityTopOfCell', ddensityTopOfCellField)
      call mpas_pool_get_field(scratchPool, 'ddensityTopOfEdge', ddensityTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'du2TopOfCell', du2TopOfCellField)
      call mpas_pool_get_field(scratchPool, 'du2TopOfEdge', du2TopOfEdgeField)
      call mpas_allocate_scratch_field(ddensityTopOfCellField, .true., .false.)
      call mpas_allocate_scratch_field(ddensityTopOfEdgeField, .true., .false.)
      call mpas_allocate_scratch_field(du2TopOfCellField, .true., .false.)
      call mpas_allocate_scratch_field(du2TopOfEdgeField, .true., .false.)
      call mpas_threading_barrier()

      ddensityTopOfCell => ddensityTopOfCellField % array
      ddensityTopOfEdge => ddensityTopOfEdgeField % array
      du2TopOfCell => du2TopOfCellField % array
      du2TopOfEdge => du2TopOfEdgeField % array

      nCells = nCellsArray( size(nCellsArray) )

      ! ddensityTopOfCell(k) = $\rho^*_{k-1}-\rho_k$, where $\rho^*$ has been adiabatically displaced to level k.
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         ddensityTopOfCell(:, iCell) = 0.0_RKIND
         du2TopOfCell(:, iCell) = 0.0_RKIND
         RiTopOfCell(:, iCell) = 0.0_RKIND
      end do
      !$omp end do

      nEdges = nEdgesArray( size(nEdgesArray) )

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         ddensityTopOfEdge(:, iEdge) = 0.0_RKIND
         du2TopOfEdge(:, iEdge) = 0.0_RKIND
         RiTopOfEdge(:, iEdge) = 0.0_RKIND
      end do
      !$omp end do

      nCells = nCellsArray( 2 )

      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         do k = 2, maxLevelCell(iCell)
            ddensityTopOfCell(k,iCell) = displacedDensity(k-1,iCell) - density(k,iCell)
          end do
      end do
      !$omp end do

      nEdges = nEdgesArray( 2 )

      !$omp do schedule(runtime) private(cell1, cell2, k)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 2, maxLevelEdgeTop(iEdge)
            ! interpolate ddensityTopOfCell to ddensityTopOfEdge
            ddensityTopOfEdge(k,iEdge) = &
               (ddensityTopOfCell(k,cell1) + &
                ddensityTopOfCell(k,cell2))/2

            ! du2TopOfEdge(k) = $u_{k-1}-u_k$
            du2TopOfEdge(k,iEdge) = (normalVelocity(k-1,iEdge) - normalVelocity(k,iEdge))**2
         end do
       end do
       !$omp end do

       nCells = nCellsArray( 1 )

      ! interpolate du2TopOfEdge to du2TopOfCell
      !$omp do schedule(runtime) private(invAreaCell, i, iEdge, k)
      do iCell = 1, nCells
        invAreaCell = 1.0_RKIND / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)

          do k = 2, maxLevelEdgeBot(iEdge)
            du2TopOfCell(k, iCell) = du2TopOfCell(k, iCell) + 0.5_RKIND * dcEdge(iEdge) * dvEdge(iEdge) &
                                   * du2TopOfEdge(k, iEdge) * invAreaCell
          end do
        end do
      end do
      !$omp end do

      ! compute RiTopOfEdge using ddensityTopOfEdge and du2TopOfEdge
      ! coef = -g/density_0/2
      coef = -gravity / rho_sw / 2.0_RKIND

      nEdges = nEdgesArray( 2 )

      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 2, maxLevelEdgeTop(iEdge)
            RiTopOfEdge(k,iEdge) = coef * ddensityTopOfEdge(k,iEdge) &
               * ( layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge) ) &
               / ( du2TopOfEdge(k,iEdge) + 1e-20_RKIND )
         end do
      end do
      !$omp end do

      nCells = nCellsArray( 1 )

      ! compute RiTopOfCell using ddensityTopOfCell and du2TopOfCell
      ! coef = -g/density_0/2
      !$omp do schedule(runtime) private(k)
      do iCell = 1,nCells
         do k = 2,maxLevelCell(iCell)
            RiTopOfCell(k,iCell) = coef * ddensityTopOfCell(k,iCell) &
               * (layerThickness(k-1,iCell) + layerThickness(k,iCell)) &
               / (du2TopOfCell(k,iCell) + 1e-20_RKIND)
         end do
      end do
      !$omp end do

      call mpas_threading_barrier()
      call mpas_deallocate_scratch_field(ddensityTopOfCellField, .true.)
      call mpas_deallocate_scratch_field(ddensityTopOfEdgeField, .true.)
      call mpas_deallocate_scratch_field(du2TopOfCellField, .true.)
      call mpas_deallocate_scratch_field(du2TopOfEdgeField, .true.)

      call mpas_timer_stop('get rich nums')

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_get_rich_numbers!}}}

!***********************************************************************
!
!  routine ocn_vmix_coefs_rich_init
!
!> \brief   Initializes ocean momentum vertical mixing quantities
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical velocity mixing in the ocean. Since a variety of
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization.
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_coefs_rich_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_use_rich_visc, config_use_rich_diff

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_rich_visc', config_use_rich_visc)
      call mpas_pool_get_config(ocnConfigs, 'config_use_rich_diff', config_use_rich_diff)

      richViscOn = config_use_rich_visc
      richDiffOn = config_use_rich_diff

!     if (config_vert_visc_type.eq.'rich') then
!         richViscOn = .true.
!     endif

!     if (config_vert_diff_type.eq.'rich') then
!         richDiffOn = .true.
!     endif


   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_rich_init!}}}

!***********************************************************************

end module ocn_vmix_coefs_rich

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
