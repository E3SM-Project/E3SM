! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hmix_del2
!
!> \brief Ocean horizontal mixing - Laplacian parameterization
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains routines for computing horizontal mixing
!>  tendencies using a Laplacian formulation.
!
!-----------------------------------------------------------------------

module ocn_hmix_del2_constant

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading
   use mpas_vector_operations
   use mpas_matrix_operations
   use ocn_constants
   use ocn_config

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

   public :: ocn_hmix_del2_constant_coeffs, &
             ocn_hmix_del2_constant_coeffs_init

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical ::  hmixDel2On  !< integer flag to determine whether del2 chosen
         
!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hmix_del2_tend
!
!> \brief   Computes tendency term for Laplacian horizontal momentum mixing
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    22 August 2011
!> \details
!>  This routine computes the horizontal mixing tendency for momentum
!>  based on a Laplacian form for the mixing, \f$\nu_2 \nabla^2 u\f$
!>  This tendency takes the
!>  form \f$\nu( \nabla divergence + k \times \nabla relativeVorticity )\f$,
!>  where \f$\nu\f$ is a viscosity and \f$k\f$ is the vertical unit vector.
!>  This form is strictly only valid for constant \f$\nu\f$ .
!
!-----------------------------------------------------------------------

   subroutine ocn_hmix_del2_constant(meshPool, del2Coefficient, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool            !< Input: mesh information

      !------ -----------------------------------------------------------
      !
      ! input /output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         del2Coefficient       !< Input: viscosity

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

      integer :: iEdge, cell1, cell2, vertex1, vertex2, k, nEdges
      integer, dimension(:), pointer :: nEdgesArray
      integer, dimension(:), pointer :: maxLevelEdgeTop

      real (kind=RKIND) :: visc2
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2

      !-----------------------------------------------------------------
      !
      ! exit if this mixing is not selected
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.hmixDel2ContantOn) return

      call mpas_timer_start("vel del2 constant")

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)

      nEdges = nEdgesArray( 1 )

      !$omp parallel
      !$omp do schedule(runtime) private(visc2)
      do iEdge = 1, nEdges
         do k = 1, maxLevelEdgeTop(iEdge)

            visc2 =  config_mom_del2 * meshScalingDel2(iEdge)

            del2Coefficient(k,iEdge) = visc2

         end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("vel del2 constant")

   !--------------------------------------------------------------------

   end subroutine ocn_hmix_del2_constant_coeffs!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_del2_init
!
!> \brief   Initializes ocean momentum Laplacian horizontal mixing
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  Laplacian horizontal momentum mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_hmix_del2_constant_coeffs_init(err)!{{{


   integer, intent(out) :: err !< Output: error flag

   !--------------------------------------------------------------------
   !
   ! set some local module variables based on input config choices
   !
   !--------------------------------------------------------------------

   err = 0

   hmixDel2ConstantOn = .false.

   if ( config_use_del2_constant ) hmixDel2ContantOn = .true.


   !--------------------------------------------------------------------

   end subroutine ocn_hmix_del2_constant_coeffs_init!}}}

!***********************************************************************

end module ocn_hmix_del2_constant

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
