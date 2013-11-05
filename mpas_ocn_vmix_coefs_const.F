! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_coefs_const
!
!> \brief MPAS ocean vertical mixing coefficients
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module contains the routines for computing 
!>  constant vertical mixing coefficients.  
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_coefs_const

   use mpas_grid_types
   use mpas_configure
   use mpas_timer

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

   private :: ocn_vel_vmix_coefs_const, &
              ocn_tracer_vmix_coefs_const

   public :: ocn_vmix_coefs_const_build, &
             ocn_vmix_coefs_const_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: constViscOn, constDiffOn

   real (kind=RKIND) :: constVisc, constDiff


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs_const_build
!
!> \brief   Computes coefficients for vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers based user choices of mixing parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs_const_build(mesh, s, d, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(inout) :: &
         s             !< Input/Output: state information

      type (diagnostics_type), intent(inout) :: &
         d             !< Input/Output: diagnostic information

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

      integer :: err1, err2

      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfEdge, vertDiffTopOfCell

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      vertViscTopOfEdge => d % vertViscTopOfEdge % array
      vertDiffTopOfCell => d % vertDiffTopOfCell % array

      call ocn_vel_vmix_coefs_const(mesh, vertViscTopOfEdge, err1)
      call ocn_tracer_vmix_coefs_const(mesh, vertDiffTopOfCell, err2)

      err = ior(err1, err2)

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_const_build!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_coefs_const
!
!> \brief   Computes coefficients for vertical momentum mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine computes the constant vertical mixing coefficients for momentum
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_coefs_const(mesh, vertViscTopOfEdge, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: vertViscTopOfEdge !< Output: vertical viscosity

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.constViscOn) return

      vertViscTopOfEdge = vertViscTopOfEdge + constVisc

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_coefs_const!}}}

!***********************************************************************
!
!  routine ocn_tracer_vmix_coefs_const
!
!> \brief   Computes coefficients for vertical tracer mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine computes the constant vertical mixing coefficients for tracers
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vmix_coefs_const(mesh, vertDiffTopOfCell, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: vertDiffTopOfCell !< Output: Vertical diffusion

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.constDiffOn) return

      vertDiffTopOfCell = vertDiffTopOfCell + constDiff

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vmix_coefs_const!}}}


!***********************************************************************
!
!  routine ocn_vmix_coefs_const_init
!
!> \brief   Initializes ocean momentum vertical mixing quantities
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  vertical velocity mixing in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_coefs_const_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      constViscOn = config_use_const_visc
      constDiffOn = config_use_const_diff
      constVisc = config_vert_visc
      constDiff = config_vert_diff

!     if (config_vert_visc_type.eq.'const') then
!         constViscOn = .true.
!         constVisc = config_vert_visc
!     endif

!     if (config_vert_diff_type.eq.'const') then
!         constDiffOn = .true.
!         constDiff = config_vert_diff
!     endif


   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_const_init!}}}

!***********************************************************************

end module ocn_vmix_coefs_const

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
