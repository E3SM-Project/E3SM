subroutine interplandm (plono, nlato, nlono, lato, rlono, &
                        landmfile, landmo) 
!
! Read LANDM_COSLAT from input file and interpolate to output grid.
! The input grid is assumed rectangular, but the output grid may
! be reduced.
!
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  include 'netcdf.inc'
!
! Input arguments
!
  integer , intent(in) :: plono                ! output longitude dimension
  integer , intent(in) :: nlato                ! number of latitudes
  integer , intent(in) :: nlono(nlato)         ! number of reduced latitudes
  real(r8), intent(in) :: lato(nlato)          ! latitude at center of grid cell
  real(r8), intent(in) :: rlono(plono,nlato)   ! longitude on (potentially reduced) output grid
  character(len=*), intent(in) :: landmfile    ! file containing input LANDM_COSLAT
!
! Output arguments
!
  real(r8), intent(out) :: landmo(plono,nlato) ! landm on reduced grid

! Local variables

  integer :: nloni
  integer :: nlati
  integer :: i,j                    ! spatial indices
  integer :: ret                    ! return code

  integer :: landmfileid                     ! netcdf file id for landm file
  integer :: londimid, latdimid              ! lon, lat dimension ids
  integer :: lonid, latid                    ! lon, lat var ids
  integer :: landmid                         ! landm variable id

  real(r8), allocatable :: landmi(:,:)       ! landm on full grid
  real(r8), allocatable :: lati(:)
  real(r8), allocatable :: loni(:)
  real(r8), allocatable :: xtemp(:,:)        ! temporary for interpolation

  ret = nf_open (landmfile, nf_nowrite, landmfileid)
  if (ret /= nf_noerr) then
     write(6,*)nf_strerror(ret)
     write(6,*)'Unable to open input file ', trim (landmfile)
     stop 999
  end if
!
! Retrieve grid info and LANDM_COSLAT field from from offline file.
!
  call wrap_inq_dimid (landmfileid, 'lat', latdimid)
  call wrap_inq_dimlen (landmfileid, latdimid, nlati)

  call wrap_inq_dimid (landmfileid, 'lon', londimid)
  call wrap_inq_dimlen (landmfileid, londimid, nloni)

  allocate (lati(nlati))
  allocate (loni(nloni))
  allocate (landmi(nloni,nlati))

  call wrap_inq_varid (landmfileid, 'lat', latid)
  call wrap_get_var8 (landmfileid, latid, lati)

  call wrap_inq_varid (landmfileid, 'lon', lonid)
  call wrap_get_var8 (landmfileid, lonid, loni)

  call wrap_inq_varid (landmfileid, 'LANDM_COSLAT', landmid)
  call wrap_get_var8 (landmfileid, landmid, landmi)

  allocate (xtemp(nloni,nlato))
!
! For rectangular -> reduced, interpolate first in latitude, then longitude
!
  do i=1,nloni
     call lininterp (landmi(i,1), nlati, nloni, lati, &
                     xtemp(i,1), nlato, nloni, lato, .false.)
  end do

  do j=1,nlato
    call lininterp (xtemp(1,j), nloni, 1, loni, &
                    landmo(1,j), nlono(j), 1, rlono(1,j), .true.)
  end do

  deallocate (xtemp)
  deallocate (lati)
  deallocate (loni)
  deallocate (landmi)

  return
end subroutine interplandm
