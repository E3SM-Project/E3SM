module varspecs_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  include 'netcdf.inc'

  type varspecs
    character*(nf_max_name) :: name
    character*8 :: vshape

    integer :: varid
    integer :: totsiz
    integer :: xtype
    integer :: nx, ny, nz
    integer :: tpos
    integer :: count(nf_max_var_dims)

    integer, pointer :: numx(:)
    real(r8), pointer :: xpos(:,:)
    real(r8), pointer :: ypos(:)
    real(r8), pointer :: zpos(:)

  end type varspecs
end module varspecs_mod
