module eos_interface

  use iso_c_binding,  only: c_int, c_bool, c_ptr, c_f_pointer
  use dimensions_mod, only: nlev, nlevp, np
  use kinds,          only: real_kind 
  use hybvcoord_mod,  only: hvcoord_t 
  use parallel_mod,   only: abortmp

  implicit none

  type(hvcoord_t) :: hvcoord

  public :: init_f90
  public :: pnh_and_exner_from_eos_f90
contains

  subroutine init_f90 (hyai, ps0, hydrostatic) bind(c)
    use control_mod, only: theta_hydrostatic_mode
    !
    ! Inputs
    !
    real (kind=real_kind), intent(in) :: hyai(nlevp)
    real (kind=real_kind), intent(in) :: ps0
    logical (kind=C_BOOL), intent(in) :: hydrostatic

    hvcoord%hyai = hyai
    hvcoord%ps0 = ps0

    theta_hydrostatic_mode = hydrostatic

  end subroutine init_f90

  subroutine pnh_and_exner_from_eos_f90(num_elems, vtheta_dp_ptr, dp3d_ptr, phi_i_ptr, pnh_ptr, exner_ptr, dpnh_dp_i_ptr) bind(c)
    use eos, only: pnh_and_exner_from_eos
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: num_elems
    type (c_ptr), intent(in) :: vtheta_dp_ptr, dp3d_ptr, phi_i_ptr, pnh_ptr, exner_ptr, dpnh_dp_i_ptr

    !
    ! Locals
    !
    real (kind=real_kind), dimension(:,:,:,:),  pointer :: vtheta_dp, dp3d, phi_i, pnh, exner, dpnh_dp_i
    integer :: ie

    call c_f_pointer(vtheta_dp_ptr, vtheta_dp, [np,np,nlev, num_elems])
    call c_f_pointer(dp3d_ptr,      dp3d,      [np,np,nlev, num_elems])
    call c_f_pointer(phi_i_ptr,     phi_i,     [np,np,nlevp,num_elems])
    call c_f_pointer(pnh_ptr,       pnh,       [np,np,nlev, num_elems])
    call c_f_pointer(exner_ptr,     exner,     [np,np,nlev,num_elems])
    call c_f_pointer(dpnh_dp_i_ptr, dpnh_dp_i, [np,np,nlev,num_elems])

    do ie=1,num_elems
      call pnh_and_exner_from_eos(hvcoord,               &
                                  vtheta_dp(:,:,:,ie),   &
                                  dp3d(:,:,:,ie),        &
                                  phi_i(:,:,:,ie),       &
                                  pnh(:,:,:,ie),         &
                                  exner(:,:,:,ie),       &
                                  dpnh_dp_i(:,:,:,ie))
    enddo

  end subroutine pnh_and_exner_from_eos_f90

end module eos_interface
