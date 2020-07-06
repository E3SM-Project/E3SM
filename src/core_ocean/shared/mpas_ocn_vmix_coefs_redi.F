!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_coefs_redi
!
!> \brief MPAS ocean vertical mixing coefficients
!> \author Doug Jacobsen
!> \date   19 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for compounding
!>  the Redi vertical mixing coefficients.
!>
!-----------------------------------------------------------------------

module ocn_vmix_coefs_redi

   use mpas_derived_types
   use mpas_pool_routines
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
   logical, pointer :: config_use_Redi

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
      real (kind=RKIND), dimension(:), pointer :: &
        RediKappa

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'k33', k33)
      call mpas_pool_get_array(diagnosticsPool, 'RediKappa', RediKappa)

      if (config_use_Redi) then
          call ocn_tracer_vmix_coefs_redi(meshPool, vertDiffTopOfCell, k33, RediKappa, err)
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

   subroutine ocn_tracer_vmix_coefs_redi(meshPool, vertDiffTopOfCell, k33, RediKappa, err)!{{{

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

      real (kind=RKIND), dimension(:), intent(in) :: &
         RediKappa

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         k33

      real (kind=RKIND) :: rediKappaCell
      real (kind=RKIND),dimension(:),pointer :: dvEdge, dcEdge, areaCell

      integer :: i, iCell, nCells, iEdge

      integer, dimension(:), pointer :: nEdgesOnCell, nCellsArray
      integer, dimension(:,:), pointer :: edgesOnCell

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.rediDiffOn) return

      call mpas_timer_start('tracer redi coef')

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

      nCells = nCellsArray(1)

      !$omp do schedule(runtime) private(i, rediKappaCell, iEdge)
      do iCell = 1, nCells
         rediKappaCell = 0.0_RKIND
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            rediKappaCell = rediKappaCell + 0.25_RKIND*RediKappa(iEdge) * dvEdge(iEdge) * dcEdge(iEdge)
         enddo

         vertDiffTopOfCell(:, iCell) = vertDiffTopOfCell(:, iCell) + rediKappaCell * k33(:, iCell) &
                        / areaCell(iCell)
      end do
      !$omp end do

      call mpas_timer_stop('tracer redi coef')

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

      call mpas_pool_get_config(ocnConfigs, 'config_use_Redi',config_use_Redi)

      if (config_use_Redi) then
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
