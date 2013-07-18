! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
   use mpas_configure
   use mpas_vector_operations
   use mpas_matrix_operations
   use mpas_tensor_operations

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

   subroutine ocn_vel_hmix_del2_tend(mesh, divergence, relativeVorticity, viscosity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         divergence      !< Input: velocity divergence

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         relativeVorticity       !< Input: relative vorticity

      type (mesh_type), intent(in) :: &
         mesh            !< Input: mesh information

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

      integer :: iEdge, nEdgesSolve, cell1, cell2, vertex1, vertex2, k
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgeMask

      real (kind=RKIND) :: u_diffusion, invLength1, invLength2, visc2
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2, &
              dcEdge, dvEdge

      !-----------------------------------------------------------------
      !
      ! exit if this mixing is not selected
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.hmixDel2On) return

      nEdgesSolve = mesh % nEdgesSolve
      maxLevelEdgeTop => mesh % maxLevelEdgeTop % array
      cellsOnEdge => mesh % cellsOnEdge % array
      verticesOnEdge => mesh % verticesOnEdge % array
      meshScalingDel2 => mesh % meshScalingDel2 % array
      edgeMask => mesh % edgeMask % array
      dcEdge => mesh % dcEdge % array
      dvEdge => mesh % dvEdge % array

      do iEdge=1,nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         invLength1 = 1.0 / dcEdge(iEdge)
         invLength2 = 1.0 / dvEdge(iEdge)

         do k=1,maxLevelEdgeTop(iEdge)

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

   subroutine ocn_vel_hmix_del2_tensor_tend(grid, normalVelocity, tangentialVelocity, viscosity, tend, scratch, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: velocity normal to an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tangentialVelocity     !< Input: velocity, tangent to an edge

      type (mesh_type), intent(in) :: &
         grid            !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend             !< Input/Output: velocity tendency

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         viscosity       !< Input/Output: viscosity

      type (scratch_type), intent(inout) :: &
         scratch !< Input/Output: Scratch structure

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

      integer :: iEdge, nEdgesSolve, nEdges, k
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: edgeMask, edgeSignOnCell

      real (kind=RKIND) :: visc2
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2
      real (kind=RKIND), dimension(:,:), pointer :: normalVectorEdge, edgeTangentVectors
      real (kind=RKIND), dimension(:,:,:), pointer :: &
         strainRateR3Cell, strainRateR3Edge, divTensorR3Cell, outerProductEdge

      !-----------------------------------------------------------------
      !
      ! exit if this mixing is not selected
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.config_use_mom_del2_tensor) return

      nEdges = grid % nEdges
      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      meshScalingDel2 => grid % meshScalingDel2 % array
      edgeMask => grid % edgeMask % array
      edgeSignOnCell => grid % edgeSignOnCell % array
      edgeTangentVectors => grid % edgeTangentVectors % array

      call mpas_allocate_scratch_field(scratch % strainRateR3Cell, .true.)
      call mpas_allocate_scratch_field(scratch % strainRateR3Edge, .true.)
      call mpas_allocate_scratch_field(scratch % divTensorR3Cell, .true.)
      call mpas_allocate_scratch_field(scratch % outerProductEdge, .true.)
      call mpas_allocate_scratch_field(scratch % normalVectorEdge, .true.)

      strainRateR3Cell => scratch % strainRateR3Cell % array
      strainRateR3Edge => scratch % strainRateR3Edge % array
      divTensorR3Cell  => scratch % divTensorR3Cell % array
      outerProductEdge => scratch % outerProductEdge % array
      normalVectorEdge => scratch % normalVectorEdge % array

      call mpas_strain_rate_R3Cell(normalVelocity, tangentialVelocity, &
         grid, edgeSignOnCell, edgeTangentVectors, .true., &
         outerProductEdge, strainRateR3Cell)

      call mpas_matrix_cell_to_edge(grid, strainRateR3Cell, .true., strainRateR3Edge)

      do iEdge=1,nEdges  ! nEdgesSolve?
         visc2 = config_mom_del2_tensor * meshScalingDel2(iEdge)

         do k=1,maxLevelEdgeTop(iEdge)

            strainRateR3Edge(:,k,iEdge) = visc2 * strainRateR3Edge(:,k,iEdge) 

            viscosity(k,iEdge) = viscosity(k,iEdge) + visc2

         end do
      end do

     ! may change boundaries to false later
     call mpas_divergence_of_tensor_R3Cell(strainRateR3Edge, grid, edgeSignOnCell, .true., divTensorR3Cell)

     call mpas_vector_R3Cell_to_normalVectorEdge(divTensorR3Cell, grid, .true., normalVectorEdge)

      do iEdge=1,nEdges  !Solve?
         do k=1,maxLevelEdgeTop(iEdge)

            tend(k,iEdge) = tend(k,iEdge) + edgeMask(k, iEdge) * normalVectorEdge(k,iEdge)

         end do
      end do

      call mpas_deallocate_scratch_field(scratch % strainRateR3Cell, .true.)
      call mpas_deallocate_scratch_field(scratch % strainRateR3Edge, .true.)
      call mpas_deallocate_scratch_field(scratch % divTensorR3Cell, .true.)
      call mpas_deallocate_scratch_field(scratch % outerProductEdge, .true.)
      call mpas_deallocate_scratch_field(scratch % normalVectorEdge, .true.)

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

   !--------------------------------------------------------------------
   !
   ! set some local module variables based on input config choices
   !
   !--------------------------------------------------------------------

   err = 0

   hmixDel2On = .false.

   if ( config_mom_del2 > 0.0 ) then
      hmixDel2On = .true.
   endif

   if(.not.config_use_mom_del2) hmixDel2On = .false.


   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del2_init!}}}

!***********************************************************************

end module ocn_vel_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
