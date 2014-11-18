! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection
!
!> \brief MPAS ocean tracer advection driver
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains driver routine for tracer advection tendencys
!>  as well as the routines for setting up advection coefficients and 
!>  initialization of the advection routines.
!
!-----------------------------------------------------------------------

module ocn_tracer_advection

   use mpas_kind_types
   use mpas_grid_types
   use mpas_sort
   use mpas_hash

   use mpas_tracer_advection_std
   use mpas_tracer_advection_mono

   use ocn_constants
     
   implicit none
   private
   save

   public :: ocn_tracer_advection_init,         &
             ocn_tracer_advection_tend

   logical :: tracerAdvOn
   logical :: monotonicOn

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_tracer_advection_tend
!
!> \brief MPAS ocean tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine is the driver routine for computing the tendency for
!>  advection of tracers.
!
!-----------------------------------------------------------------------
   subroutine ocn_tracer_advection_tend(tracers, normalThicknessFlux, w, layerThickness, verticalCellSize, dt, meshPool, tend_layerThickness, tend)!{{{

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: tracer tendency
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input/Output: tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFlux !< Input: Thickness weighted horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: w  !< Input: Vertical velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness  !< Input: Thickness field
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of a cell
      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      real (kind=RKIND), dimension(:,:), intent(in) :: tend_layerThickness !< Input: Thickness tendency information

      real (kind=RKIND), dimension(:,:), pointer :: advCoefs, advCoefs3rd

      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, nAdvCellsForEdge
      integer, dimension(:,:), pointer :: highOrderAdvectionMask, edgeSignOnCell, advCellsForEdge

      if(.not. tracerAdvOn) return

      call mpas_pool_get_array(meshPool, 'advCoefs', advCoefs)
      call mpas_pool_get_array(meshPool, 'advCoefs3rd', advCoefs3rd)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask', highOrderAdvectionMask)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(meshPool, 'advCellsForEdge', advCellsForEdge)

      if(monotonicOn) then
         call mpas_tracer_advection_mono_tend(tracers, advCoefs, advCoefs3rd, &
            nAdvCellsForEdge, advCellsForEdge, normalThicknessFlux, w, layerThickness, &
            verticalCellSize, dt, meshPool, tend_layerThickness, tend, maxLevelCell, maxLevelEdgeTop, &
            highOrderAdvectionMask, edgeSignOnCell_in = edgeSignOnCell)
      else
         call mpas_tracer_advection_std_tend(tracers, advCoefs, advCoefs3rd, &
            nAdvCellsForEdge, advCellsForEdge, normalThicknessFlux, w, layerThickness, &
            verticalCellSize, dt, meshPool, tend_layerThickness, tend, maxLevelCell, maxLevelEdgeTop, &
            highOrderAdvectionMask, edgeSignOnCell_in = edgeSignOnCell)
      endif
   end subroutine ocn_tracer_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_tracer_advection_init
!
!> \brief MPAS ocean tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine is the driver routine for initialization of 
!>  the tracer advection routines.
!
!-----------------------------------------------------------------------
   subroutine ocn_tracer_advection_init(err)!{{{

      integer, intent(inout) :: err !< Input/Output: Error flag

      integer :: err_tmp

      logical, pointer :: config_disable_tr_adv, config_dzdk_positive, config_check_tracer_monotonicity, config_monotonic
      integer, pointer :: config_horiz_tracer_adv_order, config_vert_tracer_adv_order, config_num_halos
      real (kind=RKIND), pointer :: config_coef_3rd_order

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_num_halos', config_num_halos)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_tr_adv', config_disable_tr_adv)
      call mpas_pool_get_config(ocnConfigs, 'config_dzdk_positive', config_dzdk_positive)
      call mpas_pool_get_config(ocnConfigs, 'config_check_tracer_monotonicity', config_check_tracer_monotonicity)
      call mpas_pool_get_config(ocnConfigs, 'config_horiz_tracer_adv_order', config_horiz_tracer_adv_order)
      call mpas_pool_get_config(ocnConfigs, 'config_vert_tracer_adv_order', config_vert_tracer_adv_order)
      call mpas_pool_get_config(ocnConfigs, 'config_coef_3rd_order', config_coef_3rd_order)
      call mpas_pool_get_config(ocnConfigs, 'config_monotonic', config_monotonic)

      tracerAdvOn = .true.

      if(config_disable_tr_adv) tracerAdvOn = .false.

      call mpas_tracer_advection_std_init(config_horiz_tracer_adv_order, config_vert_tracer_adv_order, config_coef_3rd_order, config_dzdk_positive, config_check_tracer_monotonicity, err_tmp)
      call mpas_tracer_advection_mono_init(config_num_halos, config_horiz_tracer_adv_order, config_vert_tracer_adv_order, config_coef_3rd_order, config_dzdk_positive, config_check_tracer_monotonicity, err_tmp)

      err = ior(err, err_tmp)

      monotonicOn = .false.

      if(config_monotonic) then
         monotonicOn = .true.
      endif

   end subroutine ocn_tracer_advection_init!}}}

end module ocn_tracer_advection
