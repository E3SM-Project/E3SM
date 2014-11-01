!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_coefs_const
!
!> \brief MPAS ocean vertical mixing coefficients
!> \author Doug Jacobsen
!> \date   19 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for compounding
!>  the Redi vertical mixing coefficients.  
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_coefs_redi

   use mpas_grid_types
   use mpas_configure
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

   private :: ocn_tracer_vmix_coefs_redi

   public :: ocn_vmix_coefs_redi_build, &
             ocn_vmix_coefs_redi_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: rediDiffOn
   logical, pointer :: config_use_standardGM
   real (kind=RKIND), pointer :: config_Redi_kappa

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs_redi_build
!
!> \brief   Computes coefficients for vertical mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers based user choices of mixing parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs_redi_build(meshPool, statePool, diagnosticsPool, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      integer, intent(in), optional :: timeLevelIn !< Input: Time level for state pool

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

      real (kind=RKIND), dimension(:,:), pointer :: &
        vertDiffTopOfCell, k33

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'k33',k33)

      if (config_use_standardGM) then
          call ocn_tracer_vmix_coefs_redi(meshPool, vertDiffTopOfCell, k33, err)
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_redi_build!}}}

!***********************************************************************
!
!  routine ocn_tracer_vmix_coefs_redi
!
!> \brief   Computes coefficients for vertical tracer mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the rediant vertical mixing coefficients for tracers
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vmix_coefs_redi(meshPool, vertDiffTopOfCell, vertRediDiff, err)!{{{

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

      real (kind=RKIND), dimension(:,:), intent(inout) :: vertDiffTopOfCell !< Output: Vertical diffusion

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in)    :: vertRediDiff !< Output: Vertical diffusion

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.rediDiffOn) return

      vertDiffTopOfCell = vertDiffTopOfCell + vertRediDiff

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vmix_coefs_redi!}}}


!***********************************************************************
!
!  routine ocn_vmix_coefs_redi_init
!
!> \brief   Initializes ocean momentum vertical mixing quantities
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  vertical velocity mixing in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_coefs_redi_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_standardGM',config_use_standardGM)
      call mpas_pool_get_config(ocnConfigs, 'config_Redi_kappa',config_Redi_kappa)

      if (config_use_standardGM) then
         rediDiffOn = .True.
      else
         rediDiffOn = .False.
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_redi_init!}}}

!***********************************************************************

end module ocn_vmix_coefs_redi

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
