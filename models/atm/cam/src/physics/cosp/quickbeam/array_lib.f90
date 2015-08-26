! ARRAY_LIB: Array procedures for F90
! Compiled/Modified:
!   07/01/06  John Haynes (haynes@atmos.colostate.edu)
!
! infind (function)
! lin_interpolate (function)
  
  module array_lib
  implicit none

  contains

! ----------------------------------------------------------------------------
! function INFIND
! ----------------------------------------------------------------------------
  function infind(list,val,sort,dist)
  use m_mrgrnk
  implicit none
!
! Purpose:
!   Finds the index of an array that is closest to a value, plus the
!   difference between the value found and the value specified
!
! Inputs:
!   [list]   an array of sequential values
!   [val]    a value to locate
! Optional input:
!   [sort]   set to 1 if [list] is in unknown/non-sequential order
!
! Returns:
!   index of [list] that is closest to [val]
!
! Optional output:
!   [dist]   set to variable containing [list([result])] - [val]
!
! Requires:
!   mrgrnk library
!
! Created:
!   10/16/03  John Haynes (haynes@atmos.colostate.edu)
! Modified:
!   01/31/06  IDL to Fortran 90
 
! ----- INPUTS -----
  real*8, dimension(:), intent(in) :: list
  real*8, intent(in) :: val  
  integer, intent(in), optional :: sort
  
! ----- OUTPUTS -----
  integer*4 :: infind
  real*8, intent(out), optional :: dist

! ----- INTERNAL -----
  real*8, dimension(size(list)) :: lists
  integer*4 :: nlist, result, tmp(1), sort_list
  integer*4, dimension(size(list)) :: mask, idx

  if (present(sort)) then
    sort_list = sort
  else
    sort_list = 0
  endif  

  nlist = size(list)
  if (sort_list == 1) then
    call mrgrnk(list,idx)
    lists = list(idx)
  else
    lists = list
  endif

  if (val >= lists(nlist)) then
    result = nlist
  else if (val <= lists(1)) then
    result = 1
  else
    mask(:) = 0
    where (lists < val) mask = 1
      tmp = minloc(mask,1)
      if (abs(lists(tmp(1)-1)-val) < abs(lists(tmp(1))-val)) then
        result = tmp(1) - 1
      else
        result = tmp(1)
      endif
  endif
  if (present(dist)) dist = lists(result)-val
  if (sort_list == 1) then
    infind = idx(result)
  else
    infind = result
  endif

  end function infind

! ----------------------------------------------------------------------------
! function LIN_INTERPOLATE
! ----------------------------------------------------------------------------  
  subroutine lin_interpolate(yarr,xarr,yyarr,xxarr,tol)
  use m_mrgrnk
  implicit none
!
! Purpose:
!   linearly interpolate a set of y2 values given a set of y1,x1,x2
!
! Inputs:
!   [yarr]    an array of y1 values
!   [xarr]    an array of x1 values
!   [xxarr]   an array of x2 values
!   [tol]     maximum distance for a match
!
! Output:
!   [yyarr]   interpolated array of y2 values
!
! Requires:
!   mrgrnk library
!
! Created:
!   06/07/06  John Haynes (haynes@atmos.colostate.edu)

! ----- INPUTS -----
  real*8, dimension(:), intent(in) :: yarr, xarr, xxarr
  real*8, intent(in) :: tol

! ----- OUTPUTS -----
  real*8, dimension(size(xxarr)), intent(out) :: yyarr

! ----- INTERNAL -----
  real*8, dimension(size(xarr)) :: ysort, xsort
  integer*4, dimension(size(xarr)) :: ist
  integer*4 :: nx, nxx, i, iloc
  real*8 :: d, m

  nx = size(xarr)
  nxx = size(xxarr)

! // xsort, ysort are sorted versions of xarr, yarr  
  call mrgrnk(xarr,ist)
  ysort = yarr(ist)
  xsort = xarr(ist)
  
  do i=1,nxx
    iloc = infind(xsort,xxarr(i),dist=d)
    if (d > tol) then
      print *, 'interpolation error'
      stop
    endif
    if (iloc == nx) then
!     :: set to the last value
      yyarr(i) = ysort(nx)
    else
!     :: is there another closeby value?
      if (abs(xxarr(i)-xsort(iloc+1)) < 2*tol) then
!       :: yes, do a linear interpolation      
        m = (ysort(iloc+1)-ysort(iloc))/(xsort(iloc+1)-xsort(iloc))
        yyarr(i) = ysort(iloc) + m*(xxarr(i)-xsort(iloc))
      else
!       :: no, set to the only nearby value
        yyarr(i) = ysort(iloc)
      endif
    endif
  enddo
  
  end subroutine lin_interpolate

  end module array_lib
