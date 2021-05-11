subroutine handle_special_cases (ncidt, ncido, nvars, ntime, unlimdimid)
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  include 'netcdf.inc'
!
! Input arguments
!
  integer ncidt, ncido
  integer nvars, ntime
  integer unlimdimid
!
! Local workspace
!
  character*(nf_max_name) :: name

  integer i, n    ! indices
  integer vt, vo
  integer natts
  integer xtype
  integer nvdims
  integer dimlen
  integer totsiz
  integer vardids(nf_max_var_dims) ! variable dimension ids
  integer start(nf_max_var_dims)
  integer count(nf_max_var_dims)
  integer ret

  character*(nf_max_name) :: attname
!
! Externals
!
  logical, external :: is_special_case
!
! Copy all the special case variables from template to output
! Probably should have separate "define" and "put" loops to optimize redef/enddef
!
  start(:) = 1

  do vt=1,nvars
    call wrap_inq_var (ncidt, vt, name, xtype, nvdims, vardids, natts)

    if (is_special_case (name, ncidt)) then
      totsiz = 1                       ! Init to size that can be multiplied
      do n=1,nvdims
        if (vardids(n) == unlimdimid) then
          count(n) = ntime
          totsiz = totsiz * ntime
        else
          call wrap_inq_dimlen (ncidt, vardids(n), dimlen)
          count(n) = dimlen
          totsiz = totsiz * dimlen
        end if
      end do

      call wrap_def_var (ncido, name, xtype, nvdims, vardids, vo)
!
! Copy attributes from input to output, then the variable itself
!
      do n=1,natts
        call wrap_inq_attname (ncidt, vt, n, attname)
        call wrap_copy_att (ncidt, vt, attname, ncido, vo)
      end do

      if (nf_enddef (ncido) /= NF_NOERR) stop 999
      
      call cpvar (ncidt, ncido, vt, vo, name, &
                  totsiz, xtype, start, count )

      ret = nf_redef (ncido)

    end if
  end do
end subroutine handle_special_cases
