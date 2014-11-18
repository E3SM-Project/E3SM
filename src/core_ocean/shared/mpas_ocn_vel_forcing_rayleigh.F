! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_rayleigh
!
!> \brief MPAS ocean Rayleigh Friction (to be used to smooth "shocks" from cold starts)
!> \author Todd Ringler
!> \date   5 January 2012
!> \details
!>  This module contains the routine for computing 
!>  tendencies based on linear Rayleigh friction.
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_rayleigh

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

   public :: ocn_vel_forcing_rayleigh_tend, &
             ocn_vel_forcing_rayleigh_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: rayleighFrictionOn
   real (kind=RKIND) :: rayleighDampingCoef


!***********************************************************************

contains

!***********************************************************************
!
!  ocn_vel_forcing_rayleigh_tend
!
!> \brief   Computes tendency term from Rayleigh friction
!> \author  Todd Ringler
!> \date    5 January 2012
!> \details 
!>  This routine computes the Rayleigh friction tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_rayleigh_tend(meshPool, normalVelocity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity    !< Input: velocity 

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

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

      integer :: iEdge, k
      integer, pointer :: nEdgesSolve
      integer, dimension(:), pointer :: maxLevelEdgeTop

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if ( .not. rayleighFrictionOn ) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      do iEdge = 1, nEdgesSolve
        do k = 1, maxLevelEdgeTop(iEdge)

           tend(k,iEdge) = tend(k,iEdge) - rayleighDampingCoef * normalVelocity(k,iEdge)

        enddo
      enddo


   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_rayleigh_tend!}}}

!***********************************************************************
!
!  ocn_vel_forcing_rayleigh_init
!
!> \brief   Initializes ocean Rayleigh friction
!> \author  Todd Ringler
!> \date    5 January 2012
!> \details 
!>  This routine initializes quantities related to 
!>  in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_rayleigh_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_Rayleigh_friction
      real (kind=RKIND), pointer :: config_Rayleigh_damping_coeff

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_Rayleigh_friction', config_Rayleigh_friction)
      call mpas_pool_get_config(ocnConfigs, 'config_Rayleigh_damping_coeff', config_Rayleigh_damping_coeff)

      rayleighDampingCoef = 0.0

      if (config_Rayleigh_friction) then
          rayleighFrictionOn = .true.
          rayleighDampingCoef = config_Rayleigh_damping_coeff
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_rayleigh_init!}}}

!***********************************************************************

end module ocn_vel_forcing_rayleigh

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
