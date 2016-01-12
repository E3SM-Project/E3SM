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
!>  and tracers based user choices of mixing parameterization.
!
!-----------------------------------------------------------------------
   subroutine ocn_vmix_coefs_rich_build(meshPool, statePool, diagnosticsPool, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

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

      integer :: err1, err2, err3, timeLevel
      integer, pointer :: indexT, indexS

      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfEdge, vertDiffTopOfCell, normalVelocity, layerThickness, layerThicknessEdge, density, displacedDensity

      real (kind=RKIND), dimension(:,:), pointer :: RiTopOfEdge, RiTopOfCell

      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_dimension(statePool, 'index_temperature', indexT)
      call mpas_pool_get_dimension(statePool, 'index_salinity', indexS)

      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfEdge', vertViscTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'RiTopOfEdge', RiTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'RiTopOfCell', RiTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'displacedDensity', displacedDensity)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)

      call mpas_timer_start("eos rich")

      ! compute in-place density
      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, 0, 'relative', density, err, timeLevelIn=timeLevel)

      ! compute displacedDensity, density displaced adiabatically to the mid-depth one layer deeper.  
      ! That is, layer k has been displaced to the depth of layer k+1.
      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, 1, 'relative', displacedDensity, err, timeLevelIn=timeLevel)

      call mpas_timer_stop("eos rich")

      call ocn_vmix_get_rich_numbers(meshPool, indexT, indexS, normalVelocity, layerThickness, layerThicknessEdge, & 
                                  density, displacedDensity, tracers, RiTopOfEdge, RiTopOfCell, err1)

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

      integer :: iEdge, k
      integer, pointer :: nEdges

      integer, dimension(:), pointer :: maxLevelEdgeTop

      real (kind=RKIND), pointer :: config_rich_mix, config_bkrd_vert_visc, config_convective_visc

      err = 0

      if(.not.richViscOn) return

      call mpas_pool_get_config(ocnConfigs, 'config_rich_mix', config_rich_mix)
      call mpas_pool_get_config(ocnConfigs, 'config_bkrd_vert_visc', config_bkrd_vert_visc)
      call mpas_pool_get_config(ocnConfigs, 'config_convective_visc', config_convective_visc)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      do iEdge = 1, nEdges
         do k = 2, maxLevelEdgeTop(iEdge)
            ! efficiency note: these if statements are inside iEdge and k loops.
            ! Perhaps there is a more efficient way to do this.
            if (RiTopOfEdge(k,iEdge)>0.0) then
               vertViscTopOfEdge(k,iEdge) = vertViscTopOfEdge(k, iEdge) + config_bkrd_vert_visc &
                  + config_rich_mix / (1.0 + 5.0*RiTopOfEdge(k,iEdge))**2
               if (vertViscTopOfEdge(k,iEdge) > config_convective_visc) then
                  vertViscTopOfEdge(k,iEdge) = config_convective_visc
               end if
            else
               ! for Ri<0 use the convective value for the viscosity
               vertViscTopOfEdge(k,iEdge) = config_convective_visc
            end if
         end do
      end do


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
!>  This routine computes the richardson vertical mixing coefficients for tracers
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

      integer :: iCell, k
      integer, pointer :: nCells

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: coef
      real (kind=RKIND), pointer :: config_density0, config_bkrd_vert_diff, config_bkrd_vert_visc, config_rich_mix, config_convective_diff

      err = 0

      if(.not.richDiffOn) return

      call mpas_pool_get_config(ocnConfigs, 'config_density0', config_density0)
      call mpas_pool_get_config(ocnConfigs, 'config_bkrd_vert_diff', config_bkrd_vert_diff)
      call mpas_pool_get_config(ocnConfigs, 'config_bkrd_vert_visc', config_bkrd_vert_visc)
      call mpas_pool_get_config(ocnConfigs, 'config_rich_mix', config_rich_mix)
      call mpas_pool_get_config(ocnConfigs, 'config_convective_diff', config_convective_diff)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      coef = -gravity / config_density0 / 2.0
      do iCell = 1, nCells
         do k = 2, maxLevelCell(iCell)
            ! efficiency note: these if statements are inside iEdge and k loops.
            ! Perhaps there is a more efficient way to do this.
            if (RiTopOfCell(k,iCell)>0.0) then
               vertDiffTopOfCell(k,iCell) = vertDiffTopOfCell(k, iCell) + config_bkrd_vert_diff &
                  + (config_bkrd_vert_visc & 
                     + config_rich_mix / (1.0 + 5.0*RiTopOfCell(k,iCell))**2) &
                  / (1.0 + 5.0*RiTopOfCell(k,iCell))
               if (vertDiffTopOfCell(k,iCell) > config_convective_diff) then
                  vertDiffTopOfCell(k,iCell) = config_convective_diff
               end if
             else
               ! for Ri<0 use the convective value for the diffusion
               vertDiffTopOfCell(k,iCell) = config_convective_diff
            end if
         end do
      end do


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

   subroutine ocn_vmix_get_rich_numbers(meshPool, indexT, indexS, normalVelocity, layerThickness, layerThicknessEdge, & !{{{
                                 density, displacedDensity, tracers, RiTopOfEdge, RiTopOfCell, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      integer, intent(in) :: indexT !< Input: index for temperature
      integer, intent(in) :: indexS !< Input: index for salinity

      real (kind=RKIND), dimension(:,:), intent(in) :: normalVelocity       !< Input: horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness       !< Input: thickness
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThicknessEdge  !< Input: thickness at edge 

      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: tracers

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
      integer :: cell1, cell2
      integer, pointer :: nVertLevels, nCells, nEdges

      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOncell, edgeSignOnCell

      real (kind=RKIND) :: coef, invAreaCell
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), allocatable :: ddensityTopOfCell, du2TopOfCell, &
                                                        ddensityTopOfEdge, du2TopOfEdge

      real (kind=RKIND), pointer :: config_density0

      err = 0

      if ( ( .not. richViscOn ) .and. ( .not. richDiffOn ) ) return

      call mpas_pool_get_config(ocnConfigs, 'config_density0', config_density0)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

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

      allocate( &
         ddensityTopOfCell(nVertLevels+1,nCells+1), ddensityTopOfEdge(nVertLevels+1,nEdges), &
         du2TopOfCell(nVertLevels+1,nCells+1), du2TopOfEdge(nVertLevels+1,nEdges))

      ! ddensityTopOfCell(k) = $\rho^*_{k-1}-\rho_k$, where $\rho^*$ has been adiabatically displaced to level k.
      ddensityTopOfCell = 0.0
      do iCell = 1, nCells
         do k = 2, maxLevelCell(iCell)
            ddensityTopOfCell(k,iCell) = displacedDensity(k-1,iCell) - density(k,iCell)
          end do
      end do

      ! interpolate ddensityTopOfCell to ddensityTopOfEdge
      ddensityTopOfEdge = 0.0
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 2, maxLevelEdgeTop(iEdge)
            ddensityTopOfEdge(k,iEdge) = &
               (ddensityTopOfCell(k,cell1) + &
                ddensityTopOfCell(k,cell2))/2  
         end do
       end do

      ! du2TopOfEdge(k) = $u_{k-1}-u_k$
      du2TopOfEdge=0.0
      do iEdge = 1, nEdges
         do k = 2, maxLevelEdgeTop(iEdge)
            du2TopOfEdge(k,iEdge) = (normalVelocity(k-1,iEdge) - normalVelocity(k,iEdge))**2
         end do
      end do

      ! interpolate du2TopOfEdge to du2TopOfCell
      du2TopOfCell = 0.0
      do iCell = 1, nCells
        invAreaCell = 1.0 / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)

          do k = 2, maxLevelEdgeBot(iEdge)
            du2TopOfCell(k, iCell) = du2TopOfCell(k, iCell) + 0.5 * dcEdge(iEdge) * dvEdge(iEdge) * du2TopOfEdge(k, iEdge) * invAreaCell
          end do
        end do
      end do

      ! compute RiTopOfEdge using ddensityTopOfEdge and du2TopOfEdge
      ! coef = -g/density_0/2
      RiTopOfEdge = 0.0
      coef = -gravity / config_density0 / 2.0
      do iEdge = 1, nEdges
         do k = 2, maxLevelEdgeTop(iEdge)
            RiTopOfEdge(k,iEdge) = coef * ddensityTopOfEdge(k,iEdge) &
               * ( layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge) ) &
               / ( du2TopOfEdge(k,iEdge) + 1e-20 )
         end do
      end do

      ! compute RiTopOfCell using ddensityTopOfCell and du2TopOfCell
      ! coef = -g/density_0/2
      RiTopOfCell = 0.0
      do iCell = 1,nCells
         do k = 2,maxLevelCell(iCell)
            RiTopOfCell(k,iCell) = coef * ddensityTopOfCell(k,iCell) &
               * (layerThickness(k-1,iCell) + layerThickness(k,iCell)) &
               / (du2TopOfCell(k,iCell) + 1e-20)
         end do
      end do

      deallocate(ddensityTopOfCell, ddensityTopOfEdge, &
        du2TopOfCell, du2TopOfEdge)

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
