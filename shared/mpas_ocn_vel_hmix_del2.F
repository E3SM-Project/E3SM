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

module ocn_vel_hmix_del2

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading
   use mpas_vector_operations
   use mpas_matrix_operations
   use ocn_constants

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

   public :: ocn_vel_hmix_del2_tend, &
             ocn_vel_hmix_del2_init

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

   subroutine ocn_vel_hmix_del2_tend(meshPool, divergence, relativeVorticity, viscosity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         divergence      !< Input: velocity divergence

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         relativeVorticity       !< Input: relative vorticity

      type (mpas_pool_type), intent(in) :: &
         meshPool            !< Input: mesh information

      !------ -----------------------------------------------------------
      !
      ! input /output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         viscosity       !< Input: viscosity

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend             !< Input/Output: velocity tendency

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
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgeMask

      real (kind=RKIND) :: u_diffusion, invLength1, invLength2, visc2
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2, &
              dcEdge, dvEdge

      real (kind=RKIND), pointer :: config_mom_del2

      !-----------------------------------------------------------------
      !
      ! exit if this mixing is not selected
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.hmixDel2On) return

      call mpas_timer_start("vel del2")

      call mpas_pool_get_config(ocnConfigs, 'config_mom_del2', config_mom_del2)

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)

      nEdges = nEdgesArray( 1 )

      !$omp do schedule(runtime) private(cell1, cell2, vertex1, vertex2, invLength1, invLength2, k, u_diffusion, visc2)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         invLength1 = 1.0_RKIND / dcEdge(iEdge)
         invLength2 = 1.0_RKIND / dvEdge(iEdge)

         do k = 1, maxLevelEdgeTop(iEdge)

            ! Here -( relativeVorticity(k,vertex2) - relativeVorticity(k,vertex1) ) / dvEdge(iEdge)
            ! is - \nabla relativeVorticity pointing from vertex 2 to vertex 1, or equivalently
            !    + k \times \nabla relativeVorticity pointing from cell1 to cell2.

            u_diffusion = ( divergence(k,cell2)  - divergence(k,cell1) ) * invLength1 &
                         -( relativeVorticity(k,vertex2) - relativeVorticity(k,vertex1) ) * invLength2

            visc2 =  config_mom_del2 * meshScalingDel2(iEdge)

            tend(k,iEdge) = tend(k,iEdge) + edgeMask(k, iEdge) * visc2 * u_diffusion

            viscosity(k,iEdge) = viscosity(k,iEdge) + visc2

         end do
      end do
      !$omp end do

      call mpas_timer_stop("vel del2")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del2_tend!}}}

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

   subroutine ocn_vel_hmix_del2_init(err)!{{{


   integer, intent(out) :: err !< Output: error flag

   real (kind=RKIND), pointer :: config_mom_del2
   logical, pointer :: config_use_mom_del2

   !--------------------------------------------------------------------
   !
   ! set some local module variables based on input config choices
   !
   !--------------------------------------------------------------------

   err = 0

   call mpas_pool_get_config(ocnConfigs, 'config_mom_del2', config_mom_del2)
   call mpas_pool_get_config(ocnConfigs, 'config_use_mom_del2', config_use_mom_del2)

   hmixDel2On = .false.

   if ( config_mom_del2 > 0.0_RKIND ) then
      hmixDel2On = .true.
   endif

   if ( .not. config_use_mom_del2 ) hmixDel2On = .false.


   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del2_init!}}}

!***********************************************************************

end module ocn_vel_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
