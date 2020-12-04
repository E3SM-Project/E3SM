! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vvel_advection
!
!> \brief MPAS ocean vertical velocity advection driver
!> \author Luke Van Roekel
!> \date   January 2021
!> \details
!>  This module contains initialization and driver routines for computing vertical
!>  velocity advection tendencies. It is based, in part, on the tracer advection code
!>
!
!-------------------------------------------------------------------------------

module ocn_vvel_advection

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_sort
   use mpas_hash
   use mpas_timer

   use ocn_vvel_horiz_advection_vector_invariant
   use ocn_vvel_horiz_advection
   use ocn_vvel_vert_advection

   use ocn_constants
   use ocn_config

   implicit none
   private
   save

   ! public module method interfaces
   public :: ocn_vvel_advection_init,         &
             ocn_vvel_advection_tend

   ! privat module variables
   logical :: vvelAdvOn !< flag to turn on tracer advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_tracer_advection_tend
!
!> \brief MPAS ocean tracer advection tendency
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This routine is the driver routine for computing tracer advection
!>  tendencies. It simply calls submodule tendency routines based on choice of
!>  algorithm.
!
!-------------------------------------------------------------------------------

   subroutine ocn_vvel_advection_tend(verticalVelocity, normalThicknessFlux,       &
                                        layerThickness, layerThicknessEdge, dt, meshPool,          &
                                        scratchPool, diagnosticsPool, tend)!{{{

      !*** Input/Output parameters

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend            !< [in,out] vertical velocity tendency to which advection added

      !*** Input parameters

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity             !< [in] Current tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalThicknessFlux !< [in] Thickness weighted horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness,     &      !< [in] Thickness field at cell centers
         layerThicknessEdge
      real (kind=RKIND), intent(in) :: &
         dt                  !< [in] Time step
      type (mpas_pool_type), intent(in) :: &
         meshPool            !< [in] Mesh information
      type (mpas_pool_type), intent(in) :: &
         scratchPool         !< [in] Scratch fields
      type (mpas_pool_type), intent(in) :: &
         diagnosticsPool     !< [in] Diagnostic fields

      !*** Local variables

      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         advCoefs, advCoefs3rd ! advection coefficients

      integer, dimension(:), pointer, contiguous :: &
         maxLevelCell,    &! index of max level at each cell
         maxLevelEdgeTop, &! max level at edge with both cells active
         nAdvCellsForEdge  ! number of advective cells for each edge
      integer, dimension(:,:), pointer, contiguous :: &!
         highOrderAdvectionMask, &! mask for higher order contributions
         edgeSignOnCell,  &! sign at cell edge for fluxes
         advCellsForEdge   ! index of advective cells for each edge

      ! end of preamble
      !----------------
      ! begin code

      ! immediate return if tracer advection not selected
      if(.not. vvelAdvOn) return

      call mpas_timer_start("Vertical Velocity adv")

      ! extract pool variables
      call mpas_pool_get_array(meshPool, 'advCoefs', advCoefs)
      call mpas_pool_get_array(meshPool, 'advCoefs3rd', advCoefs3rd)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask', &
                                          highOrderAdvectionMask)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(meshPool, 'advCellsForEdge', advCellsForEdge)


      if(config_vvel_advection_use_vector_invariant_form) then
         call ocn_vvel_horiz_vector_invariant_advection_tend(verticalVelocity, advCoefs, advCoefs3rd, &
               nAdvCellsForEdge, advCellsForEdge, normalThicknessFlux, layerThicknessEdge, &
               dt, meshPool, tend, maxLevelCell, maxLevelEdgeTop, &
               highOrderAdvectionMask, edgeSignOnCell)
      else
         call ocn_vvel_horiz_advection_tend(verticalVelocity, advCoefs, advCoefs3rd, &
               nAdvCellsForEdge, advCellsForEdge, normalThicknessFlux, layerThickness, &
               dt, meshPool, tend, maxLevelCell, maxLevelEdgeTop, &
               highOrderAdvectionMask, edgeSignOnCell)
      end if

      call ocn_vvel_vert_advection_tend( verticalVelocity, layerThickness, dt, meshPool, &
                                                tend, maxLevelCell, maxLevelEdgeTop)

      call mpas_timer_stop("Vertical Velocity adv")

   end subroutine ocn_vvel_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_advection_init
!
!> \brief MPAS ocean vertical velocity advection tendency
!> \author Luke Van Roekel
!> \date   January 2021
!> \details
!>  This routine is the driver routine for initializing various vertical velocity
!>  advection choices and variables.
!
!-------------------------------------------------------------------------------

   subroutine ocn_vvel_advection_init(err)!{{{

      !*** output parameters

      integer, intent(out) :: err !< [out] Error flag

      !*** local variables

      integer :: errTmp ! temporary error flag

      ! end preamble
      !-------------
      ! begin code

      err = 0 ! initialize error code to success

      ! set some basic flags for options
      vvelAdvOn = .not. config_disable_vvel_adv

      ! set all other options from submodule initialization routines
      if(config_vvel_advection_use_vector_invariant_form) then
         call ocn_vvel_horiz_vector_invariant_advection_init(config_vert_vvel_adv_order,   &
                                         config_coef_3rd_order,  errTmp)
      else
         call ocn_vvel_horiz_advection_init(config_vert_vvel_adv_order, &
                                         config_coef_3rd_order, errTmp)
      end if

      call ocn_vvel_horiz_advection_init(config_horiz_vvel_adv_order,    &
                                         config_coef_3rd_order, err)

      ! if an error is returned from init routines, write an error
      ! message and return a non-zero error code
      if (err /= 0 .or. errTmp /= 0) then
         err = 1
         call mpas_log_write(                                 &
            'Error encountered during vertical velocity advection init', &
            MPAS_LOG_ERR, masterOnly=.true.)
      endif

   end subroutine ocn_vvel_advection_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_vvel_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
