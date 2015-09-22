subroutine ao_i(nlon_i , nlat_i , numlon_i, lon_i , lat_i , &
                nlon_o , nlat_o , numlon_o, lon_o , lat_o , &
                mx_ovr , i_ovr  , j_ovr   , w_ovr , re    , &
                area_o , relerr )

  use shr_kind_mod, only: r8 => shr_kind_r8

! -----------------------------------------------------------------
  implicit none
! ------------------------ code history ---------------------------
! source file:       ao_i.F
! purpose:           area averaging initialization: indices and weights
! date last revised: November 1996 
! author:            Gordon Bonan
! standardized: 
! reviewed:    
! -----------------------------------------------------------------

! ------------------------ notes ----------------------------------
! get indices and weights for area-averaging between input and output grids

! o input grid does not have to be finer resolution than output grid

! o both grids must be oriented south to north, i.e., cell(lat+1)
!   must be north of cell(lat). the southern edge of the first row 
!   must be -90 (south pole) and the northern edge of the last row
!   must be +90 (north pole)

! o both grids must be oriented eastwards, i.e., cell(lon+1) must be
!   east of cell(lon). but the two grids do not have to start at the 
!   same longitude, i.e., one grid can start at dateline and go east;
!   the other grid can start at greenwich and go east. longitudes for
!   the western edge of the cells must increase continuously and span
!   360 degrees. examples
!   dateline            :  -180 to 180         (- longitudes west of greenwich)
!   greenwich           :     0 to 360
!   greenwich (centered): -dx/2 to -dx/2 + 360 (- longitudes west of greenwich)

! for each output grid cell
! o number of input grid cells that overlap with output grid cell (n_ovr)
! o longitude index (1 <= i_ovr <= nlon_i) of the overlapping input grid cell
! o latitude index  (1 <= j_ovr <= nlat_i) of the overlapping input grid cell

! for field values fld_i on an  input grid with dimensions nlon_i and nlat_i
!     field values fld_o on an output grid with dimensions nlon_o and nlat_o are
! fld_o(io,jo) = 
! fld_i(i_ovr(io,jo,     1),j_ovr(io,jo,     1)) * w_ovr(io,jo,     1) +
!                             ... + ... +
! fld_i(i_ovr(io,jo,mx_ovr),j_ovr(io,jo,mx_ovr)) * w_ovr(io,jo,mx_ovr) 

! error check: overlap weights of input cells sum to 1 for each output cell
! -----------------------------------------------------------------

! ------------------- input variables -----------------------------
  integer nlon_i              !input grid max number of input longitude points
  integer nlat_i              !input grid number of input  latitude points
  integer numlon_i(nlat_i)    !input grid number of lon points for each lat
  integer nlon_o              !output grid max number of output lon points
  integer nlat_o              !output grid number of output latitude points
  integer numlon_o(nlat_o)    !output grid number of lon points for each lat
  integer mx_ovr              !max num of input cells that overlap output cell

  real(r8) lon_i(nlon_i+1,nlat_i) !input  grid cell lon, western edge  (degrees)
  real(r8) lon_o(nlon_o+1,nlat_o) !output grid cell lon, western edge  (degrees)
  real(r8) lat_i(nlat_i+1)        !input  grid cell lat,  southern edge (degrees)
  real(r8) lat_o(nlat_o+1)        !output grid cell lat,  southern edge (degrees)
  real(r8) area_o(nlon_o,nlat_o)  !cell area on output grid 
  real(r8) re                     !radius of earth
  real(r8) relerr                 !max error: sum overlap weights ne 1
! -----------------------------------------------------------------

! ------------------- output variables ----------------------------
  integer i_ovr(nlon_o,nlat_o,mx_ovr)  !lon index, overlapping input cell
  integer j_ovr(nlon_o,nlat_o,mx_ovr)  !lat index, overlapping input cell
  real(r8)  w_ovr(nlon_o,nlat_o,mx_ovr)  !overlap weights for input cells
! -----------------------------------------------------------------

! ------------------- local variables -----------------------------
  integer io,ii        !input and output grids longitude loop index
  integer jo,ji        !input and output grids latitude  loop index
  integer n            !overlapping cell index

  real(r8) offset          !used to shift x-grid 360 degrees
  real(r8) f_ovr           !sum of overlap weights for cells on output grid
!
! Dynamic
!
  integer n_ovr(nlon_o,nlat_o) !number of overlapping input cells

! -----------------------------------------------------------------
! initialize overlap weights on output grid to zero for maximum 
! number of overlapping points. set lat and lon indices of overlapping
! input cells to dummy values. set number of overlapping cells to zero
! -----------------------------------------------------------------

  do n = 1, mx_ovr 
    do jo = 1, nlat_o      
      do io = 1, numlon_o(jo)  
        i_ovr(io,jo,n) = 1
        j_ovr(io,jo,n) = 1
        w_ovr(io,jo,n) = 0. 
      end do
    end do
  end do

  do jo = 1, nlat_o      
    do io = 1, numlon_o(jo)  
      n_ovr(io,jo) = 0
    end do
  end do

! -----------------------------------------------------------------
! first pass to find cells that overlap, area of overlap, and weights
! -----------------------------------------------------------------

  call ao (nlon_i , nlat_i , numlon_i, lon_i  , lat_i  , &
           nlon_o , nlat_o , numlon_o, lon_o  , lat_o  , &
           area_o , re     , mx_ovr  , n_ovr  , i_ovr  , &
           j_ovr  , w_ovr  )

! -----------------------------------------------------------------
! second pass to find cells that overlap, area of overlap, and weights
! -----------------------------------------------------------------

! shift x-grid to locate periodic grid intersections
! the following assumes that all lon_i(1,:) have the same value
! independent of latitude and that the same holds for lon_o(1,:)

  if (lon_i(1,1) .lt. lon_o(1,1)) then
    offset = 360.0
  else
    offset = -360.0
  end if

  do ji = 1,nlat_i
    do ii = 1, numlon_i(ji) + 1
      lon_i(ii,ji) = lon_i(ii,ji) + offset
    end do
  end do

! find overlap

  call ao (nlon_i , nlat_i , numlon_i , lon_i , lat_i  , &
           nlon_o , nlat_o , numlon_o , lon_o , lat_o  , &
           area_o , re     , mx_ovr   , n_ovr , i_ovr  , &
           j_ovr  , w_ovr  )

! restore x-grid (un-shift x-grid)

  do ji = 1,nlat_i
    do ii = 1, numlon_i(ji) + 1
      lon_i(ii,ji) = lon_i(ii,ji) - offset
    end do
  end do

! -----------------------------------------------------------------
! error check: overlap weights for input grid cells must sum to 1
! -----------------------------------------------------------------

  do jo = 1, nlat_o
    do io = 1, numlon_o(jo)
      f_ovr = 0.
      
      do n = 1, mx_ovr
        f_ovr = f_ovr + w_ovr(io,jo,n)
      end do

      if (abs(f_ovr-1.) .gt. relerr) then
        write (6,*) 'AO_I error: area not conserved for',' lon,lat = ', io,jo
        write (6,'(a30,e20.10)') '   sum of overlap weights = ', f_ovr
        call endrun
      end if

    end do
  end do

  return
end subroutine ao_i
