subroutine cpvar (ncidi, ncido, vi, vo, name, &
                  totsiz, xtype, start, count)

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  include 'netcdf.inc'
!
! Input arguments
!
  integer ncidi, ncido
  integer vi, vo
  integer xtype
  integer start(nf_max_var_dims)
  integer count(nf_max_var_dims)
  integer totsiz

  character*(nf_max_name) :: name
!
! Local workspace
!
  character, allocatable :: cbuf(:)
  real(r8), allocatable :: buf(:)
  integer, allocatable :: ibuf(:)

  if (xtype == NF_FLOAT .or. xtype == NF_DOUBLE) then

    allocate (buf(totsiz))
    call wrap_get_vara_double (ncidi, vi, start, count, buf)
    call wrap_put_vara_double (ncido, vo, start, count, buf)
    deallocate (buf)

  else if (xtype == NF_INT) then

    allocate (ibuf(totsiz))
    call wrap_get_vara_int (ncidi, vi, start, count, ibuf)
    call wrap_put_vara_int (ncido, vo, start, count, ibuf)
    deallocate (ibuf)

  else if (xtype == NF_CHAR) then
    
    allocate (cbuf(totsiz))
    call wrap_get_vara_text (ncidi, vi, start, count, cbuf)
    call wrap_put_vara_text (ncido, vo, start, count, cbuf)
    deallocate (cbuf)
    
  else

    write(6,*)'Unknown type for variable ',name
    stop 999

  end if

  return
end subroutine cpvar
