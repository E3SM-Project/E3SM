! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hmix
!
!> \brief MPAS ocean horizontal momentum mixing driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for computing 
!>  horizontal mixing tendencies.  
!>
!>  It provides an init and a tend function. Each are described below.
!
!-----------------------------------------------------------------------

module ocn_vel_hmix

   use mpas_grid_types
   use mpas_configure
   use mpas_timer
   use ocn_vel_hmix_del2
   use ocn_vel_hmix_leith
   use ocn_vel_hmix_del4

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

   public :: ocn_vel_hmix_tend, &
             ocn_vel_hmix_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: hmixOn
   type (timer_node), pointer :: del2Timer, leithTimer, del4Timer


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hmix_tend
!
!> \brief   Computes tendency term for horizontal momentum mixing
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine computes the horizontal mixing tendency for momentum
!>  based on current state and user choices of mixing parameterization.
!>  Multiple parameterizations may be chosen and added together.  These
!>  tendencies are generally computed by calling the specific routine
!>  for the chosen parameterization, so this routine is primarily a
!>  driver for managing these choices.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_tend(grid, divergence, relativeVorticity, viscosity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         divergence    !< Input: velocity divergence

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         relativeVorticity     !< Input: relative vorticity

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         viscosity     !< Input: viscosity

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

      integer :: err1, err2, err3

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      if(.not.hmixOn) return

      viscosity = 0.0

      call mpas_timer_start("del2", .false., del2Timer)
      call ocn_vel_hmix_del2_tend(grid, divergence, relativeVorticity, viscosity, tend, err1)
      call mpas_timer_stop("del2", del2Timer)

      call mpas_timer_start("leith", .false., leithTimer)
      call ocn_vel_hmix_leith_tend(grid, divergence, relativeVorticity, viscosity, tend, err2)
      call mpas_timer_stop("leith", leithTimer)

      call mpas_timer_start("del4", .false., del4Timer)
      call ocn_vel_hmix_del4_tend(grid, divergence, relativeVorticity, tend, err3)
      call mpas_timer_stop("del4", del4Timer)

      err = ior(ior(err1, err2),err3)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_init
!
!> \brief   Initializes ocean momentum horizontal mixing quantities
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  horizontal velocity mixing in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      integer :: err1, err2, err3

      hmixOn = .true.

      call ocn_vel_hmix_del2_init(err1)
      call ocn_vel_hmix_leith_init(err2)
      call ocn_vel_hmix_del4_init(err3)

      err = ior(ior(err1, err2),err3)

      if(config_disable_vel_hmix) hmixOn = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_init!}}}

!***********************************************************************

end module ocn_vel_hmix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
