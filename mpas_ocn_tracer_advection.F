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
   use mpas_configure
   use mpas_sort
   use mpas_hash

   use mpas_tracer_advection_std
   use mpas_tracer_advection_mono
     
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
   subroutine ocn_tracer_advection_tend(tracers, normalThicknessFlux, w, layerThickness, verticalCellSize, dt, mesh, tend_layerThickness, tend)!{{{

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: tracer tendency
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input/Output: tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFlux !< Input: Thickness weighted horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: w  !< Input: Vertical velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness  !< Input: Thickness field
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of a cell
      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (mesh_type), intent(in) :: mesh !< Input: mesh information
      real (kind=RKIND), dimension(:,:), intent(in) :: tend_layerThickness !< Input: Thickness tendency information

      if(.not. tracerAdvOn) return

      if(monotonicOn) then
         call mpas_tracer_advection_mono_tend(tracers, mesh % advCoefs % array, mesh % advCoefs3rd % array, &
            mesh % nAdvCellsForEdge % array, mesh % advCellsForEdge % array, normalThicknessFlux, w, layerThickness, &
            verticalCellSize, dt, mesh, tend_layerThickness, tend, mesh % maxLevelCell % array, mesh % maxLevelEdgeTop % array, &
            mesh % highOrderAdvectionMask % array, edgeSignOnCell_in = mesh % edgeSignOnCell % array)
      else
         call mpas_tracer_advection_std_tend(tracers, mesh % advCoefs % array, mesh % advCoefs3rd % array, &
            mesh % nAdvCellsForEdge % array, mesh % advCellsForEdge % array, normalThicknessFlux, w, layerThickness, &
            verticalCellSize, dt, mesh, tend_layerThickness, tend, mesh % maxLevelCell % array, mesh % maxLevelEdgeTop % array, &
            mesh % highOrderAdvectionMask % array, edgeSignOnCell_in = mesh % edgeSignOnCell % array)
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

      err = 0

      tracerAdvOn = .true.

      if(config_disable_tr_adv) tracerAdvOn = .false.

      call mpas_tracer_advection_std_init(config_horiz_tracer_adv_order, config_vert_tracer_adv_order, config_coef_3rd_order, config_dzdk_positive, config_check_tracer_monotonicity, err_tmp)
      call mpas_tracer_advection_mono_init(config_horiz_tracer_adv_order, config_vert_tracer_adv_order, config_coef_3rd_order, config_dzdk_positive, config_check_tracer_monotonicity, err_tmp)

      err = ior(err, err_tmp)

      monotonicOn = .false.

      if(config_monotonic) then
         monotonicOn = .true.
      endif

   end subroutine ocn_tracer_advection_init!}}}

end module ocn_tracer_advection
