
program diff
  use netcdf, only: nf90_open, nf90_inq_varid, nf90_get_var, NF90_NOWRITE, NF90_NOERR, nf90_strerror, nf90_close
  implicit none
  integer, parameter :: nlev = 26, ntime = 2
  integer            :: nlon = 128, nlat = 64
  integer :: c_ncid, g_ncid, c_varid, g_varid
  double precision, dimension(:,:,:), allocatable :: c_dat, g_dat, res, nrm
  double precision :: l1, l2, linf
  character(len=256) :: c_fname = 'asp_baroclinic2.cpuref.nc'
  character(len=256) :: g_fname = '../../../../asp-baroclinic/HOMME-2-0-0-high-L26-nu1e14-sub4/movies/asp_baroclinic2.nc'
  character(len=32)  :: grid    = 'ne8'
  character(len=16) :: vnames(5) = (/ 'Q' , 'Q2' , 'Q3' , 'Q4' , 'Q5' /)
  integer, dimension(4) :: dbeg, dend
  double precision, parameter :: pi = 3.1415926535897932385D0
  integer :: i

  if (command_argument_count().ge.1) call get_command_argument(1,grid)
  if (command_argument_count().ge.2) call get_command_argument(2,c_fname)
  if (command_argument_count().ge.3) call get_command_argument(3,g_fname)

  if     (trim(grid) == 'ne8'  ) then
    nlon = 128
    nlat = 64
  elseif (trim(grid) == 'ne15' ) then
    nlon = 256
    nlat = 128
  elseif (trim(grid) == 'ne30' ) then
    nlon = 512
    nlat = 256
  elseif (trim(grid) == 'ne60' ) then
    nlon = 1024
    nlat = 512
  elseif (trim(grid) == 'ne96' ) then
    nlon = 1536
    nlat = 768
  elseif (trim(grid) == 'ne120') then
    nlon = 2048
    nlat = 1024
  endif

  dbeg = (/ 1 , 1 , 1 , ntime /)
  dend = (/ nlon , nlat , nlev , 1 /)

  allocate(c_dat(nlon,nlat,nlev))
  allocate(g_dat(nlon,nlat,nlev))
  allocate(res  (nlon,nlat,nlev))
  allocate(nrm  (nlon,nlat,nlev))

  call nc_check_err( __LINE__ , nf90_open(trim(c_fname) , NF90_NOWRITE , c_ncid ) )
  call nc_check_err( __LINE__ , nf90_open(trim(g_fname) , NF90_NOWRITE , g_ncid ) )
  do i = 1 , 5
    call nc_check_err( __LINE__ , nf90_inq_varid( c_ncid , trim(vnames(i)) , c_varid ) )
    call nc_check_err( __LINE__ , nf90_inq_varid( g_ncid , trim(vnames(i)) , g_varid ) )
    call nc_check_err( __LINE__ , nf90_get_var( c_ncid , c_varid  , c_dat , dbeg , dend ) )
    call nc_check_err( __LINE__ , nf90_get_var( g_ncid , g_varid  , g_dat , dbeg , dend ) )
    res = abs( c_dat - g_dat )
    nrm = abs( c_dat )
    l1 = sum(res) / sum(nrm)
    l2 = sqrt( sum(res**2) / sum(nrm**2) )
    linf = maxval(res) / maxval(nrm)
    write(*,fmt='(A4,3ES20.10)') trim(vnames(i))//': ' , l1 , l2 , linf
  enddo
  call nc_check_err( __LINE__ , nf90_close( c_ncid ) )
  call nc_check_err( __LINE__ , nf90_close( g_ncid ) )

contains

  subroutine nc_check_err(id,status)
    implicit none
    integer, intent(in   ) :: id, status
    if (status /= NF90_NOERR) then
      write(*,*) 'ID: ',id,' Error: ',nf90_strerror(status)
      stop
    end if
  end subroutine nc_check_err 

end program diff

