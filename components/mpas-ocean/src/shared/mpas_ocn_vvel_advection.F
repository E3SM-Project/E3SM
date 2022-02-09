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
!> \authors Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date   November 2021
!> \details
!>  This module contains initialization and driver routines for computing (non-hydrostatic) 
!>  vertical velocity advection tendencies. It is based, in part, on the tracer advection 
!>  code.
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
!  routine ocn_vvel_advection_tend
!
!> \brief MPAS ocean vertical velocity advection tendency
!> \author Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date   November 2021
!> \details
!>  This routine is the driver routine for computing vvel advection
!>  tendencies. It simply calls submodule tendency routines based on choice of
!>  algorithm.
!
!-------------------------------------------------------------------------------

   subroutine ocn_vvel_advection_tend(meshPool, verticalVelocity, normalThicknessFlux,       &
                                        layerThickness, vertAleTransportTop, dt,          &
                                        tend)!{{{

      implicit none

      !*** Input/Output parameters

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend            !< [in,out] vertical velocity tendency to which advection added

      !*** Input parameters

      type (mpas_pool_type), intent(in) :: &
         meshPool                     !< [in] mesh information
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity             !< [in] non-hydrostatic vertical velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalThicknessFlux          !< [in] thickness weighted horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness               !< [in] thickness field at cell centers
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertAleTransportTop          !< [in] velocity through layer top levels
      real (kind=RKIND), intent(in) :: &
         dt                           !< [in] time-step

      ! end of preamble
      !----------------
      ! begin code

      ! immediate return if tracer advection not selected
      if(.not. vvelAdvOn) return

      call mpas_timer_start("Vertical Velocity adv")

      call ocn_vvel_horiz_advection_tend(verticalVelocity, normalThicknessFlux, &
         layerThickness, dt, meshPool, tend)

      call ocn_vvel_vert_advection_tend(verticalVelocity, layerThickness, dt, &
         vertAleTransportTop, tend)

      call mpas_timer_stop("Vertical Velocity adv")

   end subroutine ocn_vvel_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_advection_init
!
!> \brief MPAS ocean vertical velocity advection tendency
!> \author Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date   November 2021
!> \details
!>  This routine is the driver for initializing various (non-hydrostatic) 
!>  vertical velocity advection choices and variables.
!
!-------------------------------------------------------------------------------

   subroutine ocn_vvel_advection_init(err)!{{{

      implicit none

      !*** output parameters

      integer, intent(out) :: err !< [out] Error flag

      ! end preamble
      !-------------
      ! begin code

      err = 0 ! initialize error code to success

      ! set some basic flags for options
      vvelAdvOn = .not. config_disable_vvel_adv

      ! set all other options from submodule initialization routines
      if (err.eq.0) then     
         call ocn_vvel_horiz_advection_init(err)
      end if 

      if (err.eq.0) then
         call ocn_vvel_vert_advection_init (err)
      end if 

      ! if an error is returned from init routines, write an error
      ! message and return a non-zero error code
      if (err /= 0) then 
         err = 1
         call mpas_log_write(                                 &
            'Error encountered during vertical velocity advection init', &
            MPAS_LOG_ERR, masterOnly=.true.)
      endif

   end subroutine ocn_vvel_advection_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_vvel_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
