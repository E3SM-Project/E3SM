subroutine ao (nlon_i , nlat_i , numlon_i, lon_i  , lat_i  , &
               nlon_o , nlat_o , numlon_o, lon_o  , lat_o  , &
               area_o , re     , mx_ovr  , n_ovr  , i_ovr  , &
               j_ovr  , w_ovr  )

  use precision

! -----------------------------------------------------------------
  implicit none
! ------------------------ code history ---------------------------
! source file:       ao.F
! purpose:           weights and indices for area of overlap between 
!                    input and output grids
! date last revised: March 1996 
! author:            Gordon Bonan
! standardized: 
! reviewed:    
! -----------------------------------------------------------------

! ------------------- input variables -----------------------------
  integer nlon_i              !maximum number of input  longitude points
  integer nlat_i              !number of input  latitude points
  integer numlon_i(nlat_i)    !number of input lon pts for each latitude 
  integer nlon_o              !maximum number of output longitude points
  integer nlat_o              !number of output latitude points
  integer numlon_o(nlat_o)    !number of output lon pts for each latitude 
  integer mx_ovr              !maximum number of overlapping input cells

  real(r8) lon_i(nlon_i+1,nlat_i) !input  grid cell longitude, w. edge (deg)
  real(r8) lon_o(nlon_o+1,nlat_o) !output grid cell longitude, w. edge (deg)
  real(r8) lat_i(nlat_i+1)        !input  grid cell latitude, s. edge (deg)
  real(r8) lat_o(nlat_o+1)        !output grid cell latitude, s. edge (deg)
  real(r8) area_o(nlon_o,nlat_o)  !area of output grid cell
  real(r8) re                     !radius of earth
! -----------------------------------------------------------------

! ------------------- input/output variables ----------------------
  integer n_ovr(nlon_o,nlat_o       ) !number of overlapping input cells
  integer i_ovr(nlon_o,nlat_o,mx_ovr) !lon index, overlapping input cell
  integer j_ovr(nlon_o,nlat_o,mx_ovr) !lat index, overlapping input cell

  real(r8)  w_ovr(nlon_o,nlat_o,mx_ovr) !overlap weights for input cells
! -----------------------------------------------------------------

! ------------------- local variables -----------------------------
  integer io,ii        !output and input grids longitude loop index
  integer jo,ji        !output and input grids latitude  loop index

  real(r8) lonw,lone,dx    !west, east longitudes of overlap and difference
  real(r8) lats,latn,dy    !south, north latitudes of overlap and difference
  real(r8) deg2rad         !pi/180
  real(r8) a_ovr           !area of overlap
  real(r8) zero,one
  parameter (zero=0.0)   ! Needed as arg to "max"
  parameter (one=1.)     ! Needed as arg to "atan"
! -----------------------------------------------------------------

  deg2rad = (4.*atan(one)) / 180.

! -----------------------------------------------------------------
! for each output grid cell: find overlapping input grid cell and area of
! input grid cell that overlaps with output grid cell. cells overlap if:
!
! southern edge of input grid < northern edge of output grid AND
! northern edge of input grid > southern edge of output grid
!
! western edge of input grid < eastern edge of output grid AND
! eastern edge of input grid > western edge of output grid
!
!           lon_o(io,jo)      lon_o(io+1,jo)
!
!              |                   |
!              --------------------- lat_o(jo+1)
!              |                   |
!              |                   |
!    xxxxxxxxxxxxxxx lat_i(ji+1)   |
!    x         |   x               |
!    x  input  |   x   output      |
!    x  cell   |   x    cell       |
!    x  ii,ji  |   x   io,jo       |
!    x         |   x               |
!    x         ----x---------------- lat_o(jo  )
!    x             x
!    xxxxxxxxxxxxxxx lat_i(ji  )
!    x             x
! lon_i(ii,ji) lon_i(ii+1,ji)
! -----------------------------------------------------------------

! note that code does not vectorize but is only called during 
! initialization.

  do jo = 1, nlat_o
    do io = 1, numlon_o(jo)

! loop through all input grid cells to find overlap with output grid.

      do ji = 1, nlat_i                            
        if ( lat_i(ji  ).lt.lat_o(jo+1) .and. &
             lat_i(ji+1).gt.lat_o(jo  ) ) then                !lat ok

          do ii = 1, numlon_i(ji)
            if ( lon_i(ii  ,ji).lt.lon_o(io+1,jo) .and. &
                 lon_i(ii+1,ji).gt.lon_o(io  ,jo) ) then    !lon okay

! increment number of overlapping cells. make sure 0 < n_ovr < mx_ovr

              n_ovr(io,jo) = n_ovr(io,jo) + 1
!             if (n_ovr(io,jo) .gt. mx_ovr) then
!               write (6,*) 'AO error: n_ovr= ',n_ovr(io,jo), &
!                           ' exceeded mx_ovr = ',mx_ovr, &
!                           ' for output lon,lat = ',io,jo
!               call endrun
!             end if

! determine area of overlap

              lone = min(lon_o(io+1,jo),lon_i(ii+1,ji))*deg2rad !e edge
              lonw = max(lon_o(io  ,jo),lon_i(ii  ,ji))*deg2rad !w edge
              dx = max(zero,(lone-lonw))
              latn = min(lat_o(jo+1),lat_i(ji+1))*deg2rad    !n edge
              lats = max(lat_o(jo  ),lat_i(ji  ))*deg2rad    !s edge
              dy = max(zero,(sin(latn)-sin(lats)))
              a_ovr = dx*dy*re*re

! determine indices and weights. re cancels in the division by area

              i_ovr(io,jo,n_ovr(io,jo)) = ii
              j_ovr(io,jo,n_ovr(io,jo)) = ji
              w_ovr(io,jo,n_ovr(io,jo)) = a_ovr/area_o(io,jo)

            end if
          end do

        end if
      end do

    end do
  end do
      
  return
end subroutine ao
