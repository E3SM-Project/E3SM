! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing
!
!> \brief MPAS ocean forcing driver
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for computing
!>  tendencies from forcings.
!
!-----------------------------------------------------------------------

module ocn_vel_forcing

   use mpas_kind_types
   use mpas_derived_types
   use mpas_log
   use ocn_constants
   use ocn_vel_forcing_surface_stress
   use ocn_vel_forcing_explicit_bottom_drag
   use ocn_vel_forcing_topographic_wave_drag

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

   public :: ocn_vel_forcing_tend, &
             ocn_vel_forcing_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_tend
!
!> \brief   Computes tendency term from forcings
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    15 September 2011
!> \details
!>  This routine computes the forcing tendency for momentum
!>  based on current state and user choices of forcings.
!>  Multiple forcings may be chosen and added together.  These
!>  tendencies are generally computed by calling the specific routine
!>  for the chosen forcing, so this routine is primarily a
!>  driver for managing these choices.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_tend(normalVelocity, sfcFlxAttCoeff, &
                                   surfaceStress, kineticEnergyCell, &
                                   layerThickEdgeMean, tend, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity,    &!< [in] Normal velocity at edges
         kineticEnergyCell, &!< [in] kinetic energy at cell
         layerThickEdgeMean  !< [in] mean thickness at edge

      real (kind=RKIND), dimension(:), intent(in) :: &
         sfcFlxAttCoeff, &!< [in] attenuation coefficient for sfc fluxes
         surfaceStress    !< [in] surface stress at edges

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend            !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: err1 ! local error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** initialize error flag, then call individual forcing routines

      err = 0

      call ocn_vel_forcing_surface_stress_tend(sfcFlxAttCoeff, &
                   surfaceStress, layerThickEdgeMean, tend, err1)
      err = ior(err, err1)

      call ocn_vel_forcing_explicit_bottom_drag_tend(normalVelocity, &
                   kineticEnergyCell, layerThickEdgeMean, tend, err1)
      err = ior(err, err1)

      call ocn_vel_forcing_topographic_wave_drag_tend(normalVelocity, &
                   tend, err1)

      err = ior(err, err1)
   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_init
!
!> \brief   Initializes ocean forcings
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes quantities related to forcings
!>  in the ocean. Since a multiple forcings are available,
!>  this routine primarily calls the
!>  individual init routines for each forcing.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_init(err)!{{{

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: err1, err2 ! local error flags

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      err = 0

      call ocn_vel_forcing_surface_stress_init(err1)
      if (err1 /= 0) then
         call mpas_log_write( &
            'vel_forcing_init encountered error in sfc stress init', &
             MPAS_LOG_ERR)
         err = err1
      endif

      call ocn_vel_forcing_explicit_bottom_drag_init(err2)
      if (err2 /= 0) then
         call mpas_log_write( &
            'vel_forcing_init encountered error in bot drag init', &
             MPAS_LOG_ERR)
         err = err2
      endif

      call ocn_vel_forcing_topographic_wave_drag_init(err2)
      if (err2 /= 0) then
         call mpas_log_write( &
            'vel_forcing_init encountered error in topographic wave drag', &
             MPAS_LOG_ERR)
         err = err2
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_init!}}}

!***********************************************************************

end module ocn_vel_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
