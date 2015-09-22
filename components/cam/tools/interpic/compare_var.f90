subroutine compare_var (vari, varo)
  use varspecs_mod, only: varspecs

  implicit none

  include 'netcdf.inc'
!
! Input arguments
!
  type(varspecs) :: vari
  type(varspecs) :: varo
!
! Local workspace
!
  logical isfloatingi
  logical isfloatingo

  if (vari%name /= varo%name) then
    write(6,*)'compare_var: names do not match: ', trim(vari%name), trim(varo%name)
    stop 999
  end if

  isfloatingi = vari%xtype == nf_float .or. vari%xtype == nf_double
  isfloatingo = varo%xtype == nf_float .or. varo%xtype == nf_double

  if (vari%xtype /= varo%xtype .and. .not. (isfloatingi .and. isfloatingo)) then
    write(6,*)'compare_var: types are incompatible for: ', &
              trim(vari%name), ' ', trim(varo%name)
  end if

  return
end subroutine compare_var
