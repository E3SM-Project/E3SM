module elem_ops_interface

  use iso_c_binding,  only: c_int, c_bool, c_ptr, c_f_pointer
  use dimensions_mod, only: nlev, nlevp, np
  use kinds,          only: real_kind 
  use hybvcoord_mod,  only: hvcoord_t 
  use parallel_mod,   only: abortmp

  implicit none

  type(hvcoord_t) :: hvcoord

  public :: init_f90
  public :: set_theta_ref_f90
contains

  subroutine init_f90 (hyai, ps0) bind(c)
    !
    ! Inputs
    !
    real (kind=real_kind), intent(in) :: hyai(nlevp)
    real (kind=real_kind), intent(in) :: ps0

    hvcoord%hyai = hyai
    hvcoord%ps0 = ps0

  end subroutine init_f90

  subroutine set_theta_ref_f90(num_elems, dp_ptr, theta_ref_ptr) bind(c)
    use element_ops, only: set_theta_ref
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: num_elems
    type (c_ptr), intent(in) :: dp_ptr, theta_ref_ptr

    !
    ! Locals
    !
    real (kind=real_kind), dimension(:,:,:,:),  pointer :: dp, theta_ref
    integer :: ie

    call c_f_pointer(dp_ptr,        dp,        [np,np,nlev,num_elems])
    call c_f_pointer(theta_ref_ptr, theta_ref, [np,np,nlev,num_elems])

    do ie=1,num_elems
      call set_theta_ref(hvcoord,             &
                         dp(:,:,:,ie),        &
                         theta_ref(:,:,:,ie))
    enddo
  end subroutine set_theta_ref_f90

end module elem_ops_interface
