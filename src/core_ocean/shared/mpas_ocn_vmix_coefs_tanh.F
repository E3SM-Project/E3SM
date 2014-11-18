! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_coefs_tanh
!
!> \brief MPAS ocean vertical mixing coefficients
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module contains the routines for computing 
!>  tanhant vertical mixing coefficients.  
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_coefs_tanh

   use mpas_grid_types
   use mpas_timer
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

   public :: ocn_vmix_coefs_tanh_build, &
             ocn_vmix_coefs_tanh_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: tanhViscOn, tanhDiffOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs_tanh_build
!
!> \brief   Computes coefficients for vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers based user choices of mixing parameterization.
!
!-----------------------------------------------------------------------

   ! DWJ-POOL -- Remove Statepool?
   subroutine ocn_vmix_coefs_tanh_build(meshPool, statePool, diagnosticsPool, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information
      integer, intent(in), optional :: timeLevelIn

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool             !< Input/Output: state information

      type (mpas_pool_type), intent(inout) :: &
         diagnosticsPool             !< Input/Output: diagnostic information

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

      integer :: err1, err2, timeLevel

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

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfEdge', vertViscTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)

      call ocn_vel_vmix_coefs_tanh(meshPool, vertViscTopOfEdge, err1)
      call ocn_tracer_vmix_coefs_tanh(meshPool, vertDiffTopOfCell, err2)

      err = ior(err1, err2)

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_tanh_build!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_coefs_tanh
!
!> \brief   Computes coefficients for vertical momentum mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine computes the tanh vertical mixing coefficients for momentum
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_coefs_tanh(meshPool, vertViscTopOfEdge, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

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

      real (kind=RKIND), dimension(:,:), intent(out) :: vertViscTopOfEdge !< Output: Vertical viscosity

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: k
      integer, pointer :: nVertLevels

      real (kind=RKIND), dimension(:), pointer :: refBottomDepth
      real (kind=RKIND), pointer :: config_max_visc_tanh, config_min_visc_tanh, config_ZMid_tanh
      real (kind=RKIND), pointer :: config_zWidth_tanh

      err = 0

      if(.not.tanhViscOn) return

      call mpas_pool_get_config(ocnConfigs, 'config_max_visc_tanh', config_max_visc_tanh)
      call mpas_pool_get_config(ocnConfigs, 'config_min_visc_tanh', config_min_visc_tanh)
      call mpas_pool_get_config(ocnConfigs, 'config_ZMid_tanh', config_ZMid_tanh)
      call mpas_pool_get_config(ocnConfigs, 'config_zWidth_tanh', config_zWidth_tanh)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

      ! refBottomDepth is used here for simplicity.  Using zMid and h, which 
      ! vary in time, would give the exact location of the top, but it
      ! would only change the diffusion value very slightly.
      do k = 2, nVertLevels
         vertViscTopOfEdge(k,:) = vertViscTopOfEdge(k,:) - (config_max_visc_tanh - config_min_visc_tanh) / 2.0 &
            * tanh((refBottomDepth(k-1) + config_ZMid_tanh) &
                  / config_zWidth_tanh) &
            + (config_max_visc_tanh + config_min_visc_tanh) / 2
      end do


   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_coefs_tanh!}}}

!***********************************************************************
!
!  routine ocn_tracer_vmix_coefs_tanh
!
!> \brief   Computes coefficients for vertical tracer mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine computes the tanh vertical mixing coefficients for tracers
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vmix_coefs_tanh(meshPool, vertDiffTopOfCell, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

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

      integer :: k
      integer, pointer :: nVertLevels

      real (kind=RKIND), dimension(:), pointer :: refBottomDepth

      real (kind=RKIND), pointer :: config_max_diff_tanh, config_min_diff_tanh, config_ZMid_tanh
      real (kind=RKIND), pointer :: config_zWidth_tanh

      err = 0

      if ( .not. tanhDiffOn ) return

      call mpas_pool_get_config(ocnConfigs, 'config_max_diff_tanh', config_max_diff_tanh)
      call mpas_pool_get_config(ocnConfigs, 'config_min_diff_tanh', config_min_diff_tanh)
      call mpas_pool_get_config(ocnConfigs, 'config_ZMid_tanh', config_ZMid_tanh)
      call mpas_pool_get_config(ocnConfigs, 'config_zWidth_tanh', config_zWidth_tanh)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

      ! refBottomDepth is used here for simplicity.  Using zMid and h, which 
      ! vary in time, would give the exact location of the top, but it
      ! would only change the diffusion value very slightly.
      do k=2,nVertLevels
         vertDiffTopOfCell(k,:) = vertDiffTopOfCell(k,:) - (config_max_diff_tanh - config_min_diff_tanh) / 2.0 &
            * tanh((refBottomDepth(k-1) + config_ZMid_tanh) &
                  / config_zWidth_tanh) &
            + (config_max_diff_tanh + config_min_diff_tanh) / 2
      end do


   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vmix_coefs_tanh!}}}


!***********************************************************************
!
!  routine ocn_vmix_coefs_tanh_init
!
!> \brief   Initializes ocean vertical mixing quantities
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  tanh vertical mixing in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs_tanh_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_use_tanh_visc, config_use_tanh_diff

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_tanh_visc', config_use_tanh_visc)
      call mpas_pool_get_config(ocnConfigs, 'config_use_tanh_diff', config_use_tanh_diff)

      tanhViscOn = config_use_tanh_visc
      tanhDiffOn = config_use_tanh_diff

!     if (config_vert_visc_type.eq.'tanh') then
!         tanhViscOn = .true.
!     endif

!     if (config_vert_diff_type.eq.'tanh') then
!         tanhDiffOn = .true.
!     endif

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_tanh_init!}}}

!***********************************************************************

end module ocn_vmix_coefs_tanh

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
