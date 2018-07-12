! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! This module is for packaging output quantities from RRTMGP based on spectral flux profiles
!    This implementation reports the g-point fluxes
!
module mo_fluxes_bygpoint
  use mo_rte_kind,      only: wp
  use mo_fluxes,        only: ty_fluxes
  use mo_optical_props, only: ty_optical_props
  implicit none

  ! Output from radiation calculations
  !   Data components are pointers so results can be written directly into memory
  !   reduce() function accepts spectral flux profiles
  type, extends(ty_fluxes) :: ty_fluxes_bygpoint
    real(wp), dimension(:,:,:), pointer :: gpt_flux_up => NULL(), & ! Band-by-band fluxes
                                           gpt_flux_dn => NULL()    ! (ncol, nlev, nband)
    real(wp), dimension(:,:,:), pointer :: gpt_flux_net => NULL()   ! Net (down - up)
    real(wp), dimension(:,:,:), pointer :: gpt_flux_dn_dir => NULL() ! Direct flux down
  contains
    procedure :: reduce => reduce_bygpoint
    procedure :: are_desired => are_desired_bygpoint
  end type ty_fluxes_bygpoint
contains
  ! --------------------------------------------------------------------------------------
  function reduce_bygpoint(this, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir) result(error_msg)
    class(ty_fluxes_bygpoint),           intent(inout) :: this
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_up ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_dn ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    class(ty_optical_props),           intent(in   ) :: spectral_disc  !< derived type with spectral information
    logical,                           intent(in   ) :: top_at_1
    real(kind=wp), dimension(:,:,:), optional, &
                                       intent(in   ) :: gpt_flux_dn_dir! Direct flux down
    character(len=128)                               :: error_msg
    ! ------
    integer :: ncol, nlev, ngpt, nbnd
    ! ------
    error_msg = ""
    ncol = size(gpt_flux_up, DIM=1)
    nlev = size(gpt_flux_up, DIM=2)
    ngpt = size(gpt_flux_up, DIM=3)

    if(associated(this%gpt_flux_up)) then
      if(any([size(this%gpt_flux_up, 1) /= ncol,  &
              size(this%gpt_flux_up, 2) /= nlev,  &
              size(this%gpt_flux_up, 3) /= ngpt])) then
        error_msg = "reduce: gpt_flux_up array incorrectly sized (can't compute net flux either)"
      else
        this%gpt_flux_up(:,:,:) = gpt_flux_up(:,:,:)
      end if
    end if
    if(associated(this%gpt_flux_dn)) then
      if(any([size(this%gpt_flux_dn, 1) /= ncol,  &
              size(this%gpt_flux_dn, 2) /= nlev,  &
              size(this%gpt_flux_dn, 3) /= ngpt])) then
        error_msg = "reduce: gpt_flux_dn array incorrectly sized (can't compute net flux either)"
      else
        this%gpt_flux_dn(:,:,:) = gpt_flux_dn(:,:,:)
      end if
    end if
    if(associated(this%gpt_flux_net)) then
      if(any([size(this%gpt_flux_net, 1) /= ncol,  &
              size(this%gpt_flux_net, 2) /= nlev,  &
              size(this%gpt_flux_net, 3) /= ngpt])) then
        error_msg = "reduce: gpt_flux_net array incorrectly sized (can't compute net flux either)"
      else
        this%gpt_flux_net(:,:,:) = gpt_flux_dn(:,:,:) - gpt_flux_up(:,:,:)
      end if
    end if
    if(associated(this%gpt_flux_dn_dir)) then
      if(any([size(this%gpt_flux_dn_dir, 1) /= ncol,  &
              size(this%gpt_flux_dn_dir, 2) /= nlev,  &
              size(this%gpt_flux_dn_dir, 3) /= ngpt])) then
        error_msg = "reduce: gpt_flux_dn_dir array incorrectly sized (can't compute net flux either)"
      else if(present(gpt_flux_dn_dir)) then
        this%gpt_flux_dn_dir(:,:,:) = gpt_flux_dn_dir(:,:,:)
      end if
    end if
  end function reduce_bygpoint

    ! --------------------------------------------------------------------------------------
    ! Are any fluxes desired from this set of g-point fluxes? We can tell because memory will
    !   be allocated for output
    !
    function are_desired_bygpoint(this)
      class(ty_fluxes_bygpoint), intent(in   ) :: this
      logical                                  :: are_desired_bygpoint

      are_desired_bygpoint = any([associated(this%gpt_flux_up),     &
                                  associated(this%gpt_flux_dn),     &
                                  associated(this%gpt_flux_dn_dir), &
                                  associated(this%gpt_flux_net)])
    end function are_desired_bygpoint

end module mo_fluxes_bygpoint
