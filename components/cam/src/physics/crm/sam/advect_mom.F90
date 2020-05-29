module advect_mom_mod
  use advect2_mom_xy_mod
  use advect2_mom_z_mod
  implicit none

contains

  subroutine advect_mom(ncrms)
    use vars
    implicit none
    integer, intent(in) :: ncrms
    integer i,j,k

    call advect2_mom_xy(ncrms)
    call advect2_mom_z(ncrms)

  end subroutine advect_mom

end module advect_mom_mod
