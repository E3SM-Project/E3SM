module elem_ops_interface

  use iso_c_binding,  only: c_int, c_bool, c_ptr, c_f_pointer
  use dimensions_mod, only: nlev, nlevp, np
  use kinds,          only: real_kind
  use hybvcoord_mod,  only: hvcoord_t
  use parallel_mod,   only: abortmp

  implicit none

  type(hvcoord_t) :: hvcoord
  public :: init_f90
  public :: compute_r_star_f90
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

  subroutine compute_r_star_f90(num_elems, moist, Q_ptr, R_ptr) bind(c)
    use element_ops, only: get_R_star
    use control_mod, only: use_moisture
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: num_elems
    logical (kind=c_bool), intent(in) :: moist
    type (c_ptr), intent(in) :: Q_ptr, R_ptr

    !
    ! Locals
    !
    real (kind=real_kind), dimension(:,:,:,:),  pointer :: Q, R
    integer :: ie

    use_moisture = moist

    call c_f_pointer(Q_ptr, Q, [np,np,nlev,num_elems])
    call c_f_pointer(R_ptr, R, [np,np,nlev,num_elems])

    do ie=1,num_elems
      call get_R_star(R(:,:,:,ie),Q(:,:,:,ie))
    enddo
  end subroutine compute_R_star_f90

end module elem_ops_interface
