subroutine map_i (nlon_i   , nlat_i  , numlon_i, lon_i   , lat_i, &
                  nlon_o   , nlat_o  , numlon_o, lon_o   , lat_o, &
                  mxovr_i2o, iovr_i2o, jovr_i2o, wovr_i2o)

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

! ------------------------ code history ---------------------------
! source file:       map_i.F
! purpose:           driver for area averaging initialization
! date last revised: July 2000
! author:            Mariana Vertenstein
! -----------------------------------------------------------------

! ------------------------ notes ----------------------------------
! o get indices and weights for area-averaging:
!
!   from input surface grid to output model grid
!
! o input surface and output model grids can be any resolution BUT:
!
!   both grids must be oriented south to north, i.e., cell(lat+1)
!   must be north of cell(lat). the southern edge of the first row
!   must be -90 (south pole) and the northern edge of the last row
!   must be +90 (north pole)
!
!   both grids must be oriented eastwards, i.e., cell(lon+1) must be
!   east of cell(lon). but the two grids do not have to start at the
!   same longitude, i.e., one grid can start at dateline and go east;
!   the other grid can start at greenwich and go east. longitudes for
!   the western edge of the cells must increase continuously and span
!   360 degrees. examples
!   dateline            :  -180 to 180         (- longitudes west of greenwich)
!   greenwich           :     0 to 360
!   greenwich (centered): -dx/2 to -dx/2 + 360 (- longitudes west of greenwich)
!
! o field values fld_i on an  input grid with dimensions nlon_i and nlat_i =>
!   field values fld_o on an output grid with dimensions nlon_o and nlat_o as
!
!   fld_o(io,jo) =
!   fld_i(i_ovr(io,jo,     1 ),j_ovr(io,jo,     1 )) * w_ovr(io,jo,     1 ) +
!   fld_i(i_ovr(io,jo,mxovr_i),j_ovr(io,jo,mxovr_i)) * w_ovr(io,jo,mxovr_i)
!
! o error checks:
!   overlap weights of input cells sum to 1 for each output cell
!   global sums of dummy fields are conserved for input => model area-averaging
! -----------------------------------------------------------------

! ------------------- arguments -----------------------------------
  integer , intent(in) :: nlon_i                  !input grid max number of longitude points
  integer , intent(in) :: nlat_i                  !input grid number of latitude  points
  integer , intent(in) :: numlon_i(nlat_i)        !input grid number of longitude points at each lat
  real(r8), intent(in) :: lon_i(nlon_i+1,nlat_i)  !input grid cell longitude, west edge (degrees)
  real(r8), intent(in) :: lat_i(nlat_i+1)         !input grid cell latitude, south edge (degrees)
  integer , intent(in) :: nlon_o                  !model grid max number of longitude points
  integer , intent(in) :: nlat_o                  !model grid number of latitude  points
  integer , intent(in) :: numlon_o(nlat_o)        !model grid number of longitude points at each lat
  real(r8), intent(in) :: lon_o(nlon_o+1,nlat_o)  !model grid cell longitude, west edge  (degrees)
  real(r8), intent(in) :: lat_o(nlat_o+1)         !model grid cell latitude, south edge (degrees)
  integer , intent(in) :: mxovr_i2o               !max number of input cells that overlap model cell
  integer , intent(out):: iovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lon index of overlap input cell
  integer , intent(out):: jovr_i2o(nlon_o,nlat_o,mxovr_i2o) !lat index of overlap input cell
  real(r8), intent(out):: wovr_i2o(nlon_o,nlat_o,mxovr_i2o) !weight    of overlap input cell
! -----------------------------------------------------------------
!
! ------------------- local variables -----------------------------
!
  real(r8) fld_i(nlon_i,nlat_i)   !dummy input grid field
  real(r8) fld_o(nlon_o,nlat_o)   !dummy model grid field
  real(r8) area_i(nlon_i,nlat_i)  !input grid cell area
  real(r8) area_o(nlon_o,nlat_o)  !model grid cell area
  real(r8) re                     !radius of earth
  real(r8) sum_fldo               !global sum of dummy model field
  real(r8) sum_fldi               !global sum of dummy input field
  integer io,ii                   !model and input longitude loop indices
  integer jo,ji                   !model and input latitude  loop indices
  real(r8), parameter :: relerr = 0.000001 !relative error for error checks
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! get cell areas
! -----------------------------------------------------------------

  call cell_area (nlat_i, nlon_i, numlon_i, lon_i, lat_i, re, area_i)
                                                          
  call cell_area (nlat_o, nlon_o, numlon_o, lon_o, lat_o, re, area_o)

! -----------------------------------------------------------------
! get indices and weights for mapping from input grid to model grid
! -----------------------------------------------------------------

  call ao_i (nlon_i   , nlat_i   , numlon_i, lon_i    , lat_i , &
             nlon_o   , nlat_o   , numlon_o, lon_o    , lat_o , &
             mxovr_i2o, iovr_i2o , jovr_i2o, wovr_i2o , re    , &
             area_o   , relerr   )

! -----------------------------------------------------------------
! error check: global sum fld_o = global sum fld_i
! -----------------------------------------------------------------
!
! make dummy input field and sum globally
!
  sum_fldi = 0.
  do ji = 1, nlat_i      
    do ii = 1, numlon_i(ji)
      fld_i(ii,ji) = (ji-1)*nlon_i + ii
      sum_fldi = sum_fldi + area_i(ii,ji)*fld_i(ii,ji)
    end do
  end do
!
! area-average model field from input field
!
  call area_ave (nlat_i   , nlon_i   , numlon_i ,fld_i    , &
                 nlat_o   , nlon_o   , numlon_o ,fld_o    , &
                 iovr_i2o , jovr_i2o , wovr_i2o , mxovr_i2o)
!
! global sum of model field
!
  sum_fldo = 0.
  do jo = 1, nlat_o
    do io = 1, numlon_o(jo)
      sum_fldo = sum_fldo + area_o(io,jo)*fld_o(io,jo)
    end do
  end do
!
! check for conservation
!
  if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
    write (6,*) 'map_i error srf => model: srf field not conserved'
    write (6,'(a23,e20.10)') 'global sum model field = ',sum_fldo
    write (6,'(a23,e20.10)') 'global sum srf field   = ',sum_fldi
    call endrun
  end if

  return
end subroutine map_i
