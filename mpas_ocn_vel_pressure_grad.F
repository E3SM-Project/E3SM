! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_pressure_grad
!
!> \brief MPAS ocean pressure gradient module
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module contains the routine for computing 
!>  tendencie from the horizontal pressure gradient.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_pressure_grad

   use mpas_grid_types
   use mpas_configure
   use mpas_constants

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

   public :: ocn_vel_pressure_grad_tend, &
             ocn_vel_pressure_grad_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: pgradOn
   real (kind=RKIND) :: density0Inv, gdensity0Inv


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_pressure_grad_tend
!
!> \brief   Computes tendency term for horizontal pressure gradient
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine computes the pressure gradient tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_pressure_grad_tend(grid, pressure, zMid, density, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         pressure, & !< Input: Pressure field or Mongomery potential
         zMid, &     !< Input: z-coordinate at mid-depth of layer
         density         !< Input: density

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

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

      integer :: nEdgesSolve, iEdge, k, cell1, cell2
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask

      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND) :: invdcEdge

      err = 0

      if(.not.pgradOn) return

      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array
      dcEdge => grid % dcEdge % array
      edgeMask => grid % edgeMask % array

      ! pressure for generalized coordinates
      ! -1/density_0 (grad p_k + density g grad z_k^{mid})

      ! For pure isopycnal coordinates, we are still using 
      ! grad(M), the gradient of Montgomery Potential, because
      ! we have set density0Inv=1 and gdensity0Inv=0 in the init routine,
      ! and pressure is passed in as montgomeryPotential.

      do iEdge=1,nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invdcEdge = 1.0 / dcEdge(iEdge)

         do k=1,maxLevelEdgeTop(iEdge)
            tend(k,iEdge) = tend(k,iEdge)     &
              - edgeMask(k,iEdge) * density0Inv*(  pressure(k,cell2) &
                         - pressure(k,cell1) )* invdcEdge &
              - edgeMask(k,iEdge) * gdensity0Inv*  0.5*(density(k,cell1)+density(k,cell2)) &
                        *(  zMid(k,cell2) &
                          - zMid(k,cell1) )* invdcEdge
                      
         end do

      end do



   !--------------------------------------------------------------------

   end subroutine ocn_vel_pressure_grad_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_pressure_grad_init
!
!> \brief   Initializes ocean momentum horizontal pressure gradient
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine initializes parameters required for the computation of the
!>  horizontal pressure gradient.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_pressure_grad_init(err)!{{{

   !--------------------------------------------------------------------


      !-----------------------------------------------------------------
      !
      ! Output Variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag


      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      err = 0

      pgradOn = .true.

      if (config_pressure_gradient_type.eq.'MontgomeryPotential') then
        density0Inv = 1.0
        gdensity0Inv = 0.0
      else 
        density0Inv = 1.0/config_density0
        gdensity0Inv = gravity/config_density0
      end if

      if(config_disable_vel_pgrad) pgradOn = .false.



   !--------------------------------------------------------------------

   end subroutine ocn_vel_pressure_grad_init!}}}

!***********************************************************************

end module ocn_vel_pressure_grad

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
