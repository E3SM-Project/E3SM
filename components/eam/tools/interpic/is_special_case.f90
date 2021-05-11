logical function is_special_case (name, ncid)

  implicit none

  include 'netcdf.inc'
!
! Input arguments
!
  integer ncid
  character*(*) name
!
! Local workspace
!
  character*(nf_max_name) dimname
  integer n
  integer ndims
!
! Hardwire the names known to be functions of coordinate variables
!
  if (name == 'rlon' .or. name == 'nlon' .or. name == 'wnummax' .or. &
      name == 'hyai' .or. name == 'hybi' .or. name == 'hyam' .or. &
      name == 'hybm' .or. name == 'gw') then
    is_special_case = .true.
    return
  end if
!
! Loop through the dimensions and see if "name" matches a dimension
!
  if (nf_inq_ndims (ncid, ndims) /= nf_noerr) then
    call err_exit ('is_special_case: nf_inq_ndims failure')
  end if

  do n=1,ndims
    call wrap_inq_dimname (ncid, n, dimname)
    if (dimname == name) then
      is_special_case = .true.
      return
    end if
  end do

  is_special_case = .false.
  return
end function is_special_case
