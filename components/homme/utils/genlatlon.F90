! ============================================
! Construct a gaussian lat lon grid
! ============================================

program genlatlon
  ! =========================
  use kinds
  ! =========================
  use math_constants
  ! =========================
  use quadrature_mod
  ! =========================
implicit none


  integer nlon,nlat,ngrid
  character(len=80)instring
  character(len=80)gridfname
  character(len=80)gridtype
  type (quadrature_t) :: gs

  integer i,j
  real*8, dimension(:), allocatable :: lat
  real*8, dimension(:), allocatable :: lon

  gridtype = "latlong"
  print *,"enter  number of equally spaced longitude pts:"
  read(5,*)nlon

  print *,"enter  number of Gaussian latitude pts:"
  read(5,*)nlat

  allocate(lon(nlon))
  allocate(lat(nlat))

  ngrid=nlon*nlat

  print *,"enter grid file name:"
  read(5,*)instring

  gridfname = TRIM(ADJUSTL(instring))
  open(unit=7,file=gridfname,form="FORMATTED",status="REPLACE")

  gs=gauss(nlat)

  do i=1,nlon
     lon(i)=(i-1)*(2.0D0*DD_PI/REAL(nlon,kind=real_kind))
  end do        

  do j=1,nlat
     lat(j) = ASIN(gs%points(j))
  end do

#if 0
  write(7,*)ngrid
  do j=1,nlat
     do i=1,nlon
        write(7,10)lon(i),lat(j)
 10     format(e22.16,1x,e22.16)
     end do
  end do
#else
  write(7,*)gridtype
  write(7,*)nlon
  write(7,*)nlat

  do i=1,nlon
     write(7,10)lon(i)
 10  format(e22.16)
  end do

  do j=1,nlat
     write(7,10)lat(j)
  end do
#endif

  close(7)

end program genlatlon
