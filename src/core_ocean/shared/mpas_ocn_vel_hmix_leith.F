! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hmix_leith
!
!> \brief Ocean horizontal mixing - Leith parameterization 
!> \author Mark Petersen
!> \date   22 October 2012
!> \details
!>  This module contains routines for computing horizontal mixing 
!>  tendencies using the Leith parameterization.
!
!-----------------------------------------------------------------------

module ocn_vel_hmix_leith

   use mpas_grid_types
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

   public :: ocn_vel_hmix_leith_tend, &
             ocn_vel_hmix_leith_init

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical ::  hmixLeithOn  !< integer flag to determine whether leith chosen

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hmix_leith_tend
!
!> \brief  Computes tendency term for horizontal momentum mixing with Leith parameterization
!> \author Mark Petersen, Todd Ringler
!> \date   22 October 2012
!> \details 
!> This routine computes the horizontal mixing tendency for momentum
!> based on the Leith closure.  The Leith closure is the
!> enstrophy-cascade analogy to the Smagorinsky (1963) energy-cascade
!> closure, i.e. Leith (1996) assumes an inertial range of enstrophy flux
!> moving toward the mesh scale. The assumption of an enstrophy cascade
!> and dimensional analysis produces right-hand-side dissipation,
!> $\bf{D}$, of velocity of the form
!> $ {\bf D} = \nabla \cdot \left( \nu_\ast \nabla {\bf u} \right) 
!>    = \nabla \cdot \left( \gamma \left| \nabla \omega  \right| 
!>      \left( \Delta x \right)^3 \nabla \bf{u} \right)
!> where $\omega$ is the relative vorticity and $\gamma$ is a non-dimensional, 
!> $O(1)$ parameter. We set $\gamma=1$.

!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_leith_tend(meshPool, divergence, relativeVorticity, viscosity, tend, err)!{{{

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

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend             !< Input/Output: velocity tendency

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         viscosity       !< Input: viscosity

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

      integer :: iEdge, cell1, cell2, vertex1, vertex2, k
      integer, pointer :: nEdgesSolve
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgeMask

      real (kind=RKIND) :: u_diffusion, invLength1, invLength2, visc2
      real (kind=RKIND), dimension(:), pointer :: meshScaling, &
              dcEdge, dvEdge

      real (kind=RKIND), pointer :: config_leith_parameter, config_leith_dx, config_leith_visc2_max

      !-----------------------------------------------------------------
      !
      ! exit if this mixing is not selected
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.hmixLeithOn) return

      call mpas_pool_get_config(ocnConfigs, 'config_Leith_parameter', config_leith_parameter)
      call mpas_pool_get_config(ocnConfigs, 'config_Leith_dx', config_leith_dx)
      call mpas_pool_get_config(ocnConfigs, 'config_Leith_visc2_max', config_leith_visc2_max)

      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'meshScaling', meshScaling)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)

      do iEdge = 1, nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         invLength1 = 1.0 / dcEdge(iEdge)
         invLength2 = 1.0 / dvEdge(iEdge)

         do k = 1, maxLevelEdgeTop(iEdge)

            ! Here -( relativeVorticity(k,vertex2) - relativeVorticity(k,vertex1) ) / dvEdge(iEdge)
            ! is - \nabla relativeVorticity pointing from vertex 2 to vertex 1, or equivalently 
            !    + k \times \nabla relativeVorticity pointing from cell1 to cell2.

            u_diffusion = ( divergence(k,cell2)  - divergence(k,cell1) ) * invLength1 &
                         -( relativeVorticity(k,vertex2) - relativeVorticity(k,vertex1) ) * invLength2

            ! Here the first line is (\delta x)^3
            ! the second line is |\nabla \omega|
            ! and u_diffusion is \nabla^2 u (see formula for $\bf{D}$ above).
            visc2 = ( config_leith_parameter * config_leith_dx * meshScaling(iEdge) / 3.14)**3 &
                     * abs( relativeVorticity(k,vertex2) - relativeVorticity(k,vertex1) ) * invLength1 * sqrt(3.0)
            visc2 = min(visc2, config_leith_visc2_max)

            tend(k,iEdge) = tend(k,iEdge) + edgeMask(k, iEdge) * visc2 * u_diffusion

            viscosity(k,iEdge) = viscosity(k,iEdge) + visc2

         end do
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_leith_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_leith_init
!
!> \brief   Initializes ocean momentum horizontal mixing with Leith parameterization
!> \author Mark Petersen
!> \date   22 October 2012
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  Leith parameterization for horizontal momentum mixing in the ocean.  
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_leith_init(err)!{{{


   integer, intent(out) :: err !< Output: error flag

   !--------------------------------------------------------------------
   !
   ! set some local module variables based on input config choices
   !
   !--------------------------------------------------------------------
   logical, pointer :: config_use_Leith_del2

   err = 0

   call mpas_pool_get_config(ocnConfigs, 'config_use_Leith_del2', config_use_Leith_del2)

   hmixLeithOn = .false.

   if (config_use_leith_del2) then
      hmixLeithOn = .true.
   endif

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_leith_init!}}}

!***********************************************************************

end module ocn_vel_hmix_leith

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
