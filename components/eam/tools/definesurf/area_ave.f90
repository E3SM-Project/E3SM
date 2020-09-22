subroutine area_ave (nlat_i , nlon_i , numlon_i, fld_i , &
                     nlat_o , nlon_o , numlon_o, fld_o , &
                     i_ovr  , j_ovr  , w_ovr   , nmax  )

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
! ------------------------ code history ---------------------------
! source file:       area_ave.F
! purpose:           area averaging of field from input to output grids
! date last revised: November 1996
! author:            Gordon Bonan
! standardized:
! reviewed:
! -----------------------------------------------------------------

! ------------------- input variables -----------------------------
  integer nlat_i   ! number of latitude points for input grid
  integer nlat_o   ! number of latitude points for output grid
  integer nlon_i   ! maximum number of longitude points for input grid
  integer nlon_o   ! maximum number of longitude points for output grid
  integer nmax     ! maximum number of overlapping cells
  integer numlon_i(nlat_i)  ! input grid number of lon points at each lat
  integer numlon_o(nlat_o)  ! input grid number of lon points at each lat
  integer i_ovr(nlon_o,nlat_o,nmax)  ! lon index, overlapping input cell
  integer j_ovr(nlon_o,nlat_o,nmax)  ! lat index, overlapping input cell

  real(r8) fld_i(nlon_i,nlat_i) !field for input grid
  real(r8) w_ovr(nlon_o,nlat_o,nmax)  ! overlap weights for input cells
! -----------------------------------------------------------------

! ------------------- output variables ----------------------------
  real(r8) fld_o(nlon_o,nlat_o) !field for output grid
! -----------------------------------------------------------------

! ------------------- local variables -----------------------------
  integer jo,ji             !latitude index for output,input grids
  integer io,ii             !longitude index for output,input grids
  integer n                 !overlapping cell index
! -----------------------------------------------------------------

  do jo = 1, nlat_o
    do io =1, numlon_o(jo)
      fld_o(io,jo) = 0.
    end do
  end do

  do n = 1, nmax
    do jo = 1, nlat_o
      do io =1, numlon_o(jo)
        ii = i_ovr(io,jo,n)
        ji = j_ovr(io,jo,n)
        fld_o(io,jo) = fld_o(io,jo) + w_ovr(io,jo,n)*fld_i(ii,ji)
      end do
    end do
  end do

  return
end subroutine area_ave
