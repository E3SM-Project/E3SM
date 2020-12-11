module diffuse_mom_mod
  use diffuse_mom2D_mod
  use diffuse_mom3D_mod
  implicit none

contains

  subroutine diffuse_mom(ncrms,grdf_x, grdf_y, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)

    !  Interface to the diffusion routines

    use vars
    implicit none
    integer, intent(in) :: ncrms
    integer i,j,k
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tk(ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(ncrms,nzm)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z

    !call t_startf ('diffuse_mom')

    if(RUN3D) then
      call diffuse_mom3D(ncrms,grdf_x, grdf_y, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)
    else
      call diffuse_mom2D(ncrms,grdf_x,         grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)
    endif

    !call t_stopf ('diffuse_mom')

  end subroutine diffuse_mom

end module diffuse_mom_mod
