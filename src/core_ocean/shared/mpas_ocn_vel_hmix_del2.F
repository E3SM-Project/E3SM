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

   use mpas_grid_types
   use mpas_vector_operations
   use mpas_matrix_operations
   use mpas_tensor_operations
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
             ocn_vel_hmix_del2_tensor_tend, &
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

      integer :: iEdge, cell1, cell2, vertex1, vertex2, k
      integer, pointer :: nEdgesSolve
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

      call mpas_pool_get_config(ocnConfigs, 'config_mom_del2', config_mom_del2)

      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)
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

            visc2 =  config_mom_del2 * meshScalingDel2(iEdge)

            tend(k,iEdge) = tend(k,iEdge) + edgeMask(k, iEdge) * visc2 * u_diffusion

            viscosity(k,iEdge) = viscosity(k,iEdge) + visc2

         end do
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del2_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_del2_tensor_tend
!
!> \brief   Computes tendency term for Laplacian horizontal momentum mixing
!> \author  Mark Petersen
!> \date    July 2013
!> \details 
!>  This routine computes the horizontal mixing tendency for momentum
!>  using tensor operations, 
!>  based on a Laplacian form for the mixing, \f$\nabla\cdot( \nu_2 \nabla(u))\f$
!>  where \f$\nu_2\f$ is a viscosity.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del2_tensor_tend(meshPool, normalVelocity, tangentialVelocity, viscosity, scratchPool, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: velocity normal to an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tangentialVelocity     !< Input: velocity, tangent to an edge

      type (mpas_pool_type), intent(in) :: &
         meshPool            !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         viscosity       !< Input/Output: viscosity

      type (mpas_pool_type), intent(inout) :: &
         scratchPool !< Input/Output: Scratch structure

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

      integer :: iEdge, k
      integer, pointer :: nEdgesSolve, nEdges, nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: edgeMask, edgeSignOnCell

      real (kind=RKIND) :: visc2
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2
      real (kind=RKIND), dimension(:,:), pointer :: normalVectorEdge, edgeTangentVectors
      real (kind=RKIND), dimension(:,:,:), pointer :: &
         strainRateR3Cell, strainRateR3Edge, divTensorR3Cell, outerProductEdge

      type (field2DReal), pointer :: normalVectorEdgeField
      type (field3DReal), pointer :: strainRateR3CellField, strainRateR3EdgeField, divTensorR3CellField, outerProductEdgeField

      logical, pointer :: config_use_mom_del2_tensor
      real (kind=RKIND), pointer :: config_mom_del2_tensor 

      !-----------------------------------------------------------------
      !
      ! exit if this mixing is not selected
      !
      !-----------------------------------------------------------------

      err = 0
      call mpas_pool_get_config(ocnConfigs, 'config_use_mom_del2_tensor', config_use_mom_del2_tensor)

      if ( .not. config_use_mom_del2_tensor ) return

      call mpas_pool_get_config(ocnConfigs, 'config_mom_del2_tensor ', config_mom_del2_tensor )

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors)

      call mpas_pool_get_field(scratchPool, 'strainRateR3Cell',strainRateR3CellField)
      call mpas_pool_get_field(scratchPool, 'strainRateR3Edge',strainRateR3EdgeField)
      call mpas_pool_get_field(scratchPool, 'divTensorR3Cell', divTensorR3CellField)
      call mpas_pool_get_field(scratchPool, 'outerProductEdge',outerProductEdgeField)
      call mpas_pool_get_field(scratchPool, 'normalVectorEdge',normalVectorEdgeField)

      call mpas_allocate_scratch_field(strainRateR3CellField, .true.)
      call mpas_allocate_scratch_field(strainRateR3EdgeField, .true.)
      call mpas_allocate_scratch_field(divTensorR3CellField, .true.)
      call mpas_allocate_scratch_field(outerProductEdgeField, .true.)
      call mpas_allocate_scratch_field(normalVectorEdgeField, .true.)

      strainRateR3Cell => strainRateR3CellField % array
      strainRateR3Edge => strainRateR3EdgeField % array
      divTensorR3Cell  => divTensorR3CellField % array
      outerProductEdge => outerProductEdgeField % array
      normalVectorEdge => normalVectorEdgeField % array

      call mpas_strain_rate_R3Cell(normalVelocity, tangentialVelocity, &
         meshPool, edgeSignOnCell, edgeTangentVectors, .true., &
         outerProductEdge, strainRateR3Cell)

      call mpas_matrix_cell_to_edge(strainRateR3Cell, meshPool, .true., strainRateR3Edge)

      ! The following loop could possibly be reduced to nEdgesSolve
      do iEdge = 1, nEdges  
         visc2 = config_mom_del2_tensor * meshScalingDel2(iEdge)
         do k = 1, maxLevelEdgeTop(iEdge)
            strainRateR3Edge(:,k,iEdge) = visc2 * strainRateR3Edge(:,k,iEdge) 
            viscosity(k,iEdge) = viscosity(k,iEdge) + visc2
         end do
         ! Impose zero strain rate at land boundaries
         do k = maxLevelEdgeTop(iEdge)+1, nVertLevels
            strainRateR3Edge(:,k,iEdge) = 0.0
         end do
      end do

      ! may change boundaries to false later
      call mpas_divergence_of_tensor_R3Cell(strainRateR3Edge, meshPool, edgeSignOnCell, .true., divTensorR3Cell)

      call mpas_vector_R3Cell_to_normalVectorEdge(divTensorR3Cell, meshPool, .true., normalVectorEdge)

      ! The following loop could possibly be reduced to nEdgesSolve
      do iEdge = 1, nEdges 
         do k = 1, maxLevelEdgeTop(iEdge)
            tend(k,iEdge) = tend(k,iEdge) + edgeMask(k, iEdge) * normalVectorEdge(k,iEdge)
         end do
      end do

      call mpas_deallocate_scratch_field(strainRateR3CellField, .true.)
      call mpas_deallocate_scratch_field(strainRateR3EdgeField, .true.)
      call mpas_deallocate_scratch_field(divTensorR3CellField, .true.)
      call mpas_deallocate_scratch_field(outerProductEdgeField, .true.)
      call mpas_deallocate_scratch_field(normalVectorEdgeField, .true.)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del2_tensor_tend!}}}

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

   if ( config_mom_del2 > 0.0 ) then
      hmixDel2On = .true.
   endif

   if ( .not. config_use_mom_del2 ) hmixDel2On = .false.


   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del2_init!}}}

!***********************************************************************

end module ocn_vel_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
