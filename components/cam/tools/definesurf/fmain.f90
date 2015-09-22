program fmain

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  include 'netcdf.inc'
!
! Local workspace
!
  real(r8), parameter :: fillvalue = 1.d36
  real(r8), parameter :: filter_coefficient = 0.25D0

  character(len=128) :: topofile  = ' ' ! input high resolution (10 min) file name
  character(len=128) :: landmfile = ' ' ! input land mask file name
  character(len=128) :: gridfile  = ' ' ! input initial condition file with grid definition
  character(len=128) :: outbcfile = ' ' ! output boundary condition file with PHIS, SGH, etc.
  character(len= 80) :: arg             ! used for parsing command line arguments
  character(len=256) :: cmdline         ! input command line
  character(len=256) :: history         ! history attribute text
  character(len=  8) :: datestring
  character(len= 32) :: z_filter_type   ! type of filter applied to height
  character(len= 32) :: s_filter_type   ! type of filter applied to standard deviations

  logical verbose                       ! Add print statements
  logical make_ross                     ! Make Ross ice shelf south of -79
  logical filter_del2                   ! Execute SJ Lin's del2 terrain filter
  logical filter_remap                  ! Execute SJ Lin's newer remapping terrain filter
  logical filter_sgh                    ! Filter SGH and SGH30 in addition to height
  logical reduced_grid                  ! reduced grid defined
  logical have_sgh30                    ! input topofile has sgh30, output will also

  integer cmdlen                        ! character array lengths
  integer gridid
  integer foutid                        ! output file id
  integer lonid, londimid, rlonid       ! longitude dimension variable ids
  integer latid, latdimid               ! latitude  dimension variable ids
  integer sghid, phisid, landfid, nlonid, landmid, sgh30id ! output variable netcdf ids
  integer start(4), count(4)
  integer plon, nlat
  integer i, j
  integer ret
  integer nargs                        ! input arg
  integer n                            ! index loops thru input args

  integer dim(2)                        ! dimension list for output variables

  integer , allocatable :: nlon(:)
  real(r8), allocatable :: mlatcnts(:)    ! model cell center latitudes
  real(r8), allocatable :: mloncnts(:,:)  ! model cell center longitudes
  real(r8), allocatable :: sgh(:,:)
  real(r8), allocatable :: sgh30(:,:)
  real(r8), allocatable :: phis(:,:)
  real(r8), allocatable :: fland(:,:)
  real(r8), allocatable :: landm(:,:)

  integer iargc
  external iargc
!
! Default settings before parsing argument list
!
  verbose      = .false.
  make_ross    = .true.
  filter_del2  = .false.
  filter_remap = .false.
  filter_sgh   = .false.
  reduced_grid = .false.

! parse input arguments

  nargs = iargc()
  n = 1
  cmdline = char(10) // 'definesurf '
  do while (n .le. nargs)
    arg = ' '
    call getarg (n, arg)
    n = n + 1

    select case (arg)
! topography file name (10')
    case ('-t')
      call getarg (n, arg)
      n = n + 1
      topofile = arg
      cmdline = trim(cmdline) // ' -t ' // trim(topofile)
! grid file name
    case ('-g')
      call getarg (n, arg)
      n = n + 1
      gridfile = arg
      cmdline = trim(cmdline) // ' -g ' // trim(gridfile)
! verbose mode
    case ('-v')
      verbose = .true.
      cmdline = trim(cmdline) // ' -v'
! landmask file name
    case ('-l')
      call getarg (n, arg)
      n = n + 1
      landmfile = arg
      cmdline = trim(cmdline) // ' -l ' // trim(landmfile)
! extend Ross Sea
    case ('-r')
      make_ross = .false.
      cmdline = trim(cmdline) // ' -r'
! use del2 filter on heights
    case ('-del2')
      filter_del2 = .true.
      cmdline = trim(cmdline) // ' -del2'
! use remap filter on heights
    case ('-remap')
      filter_remap = .true.
      cmdline = trim(cmdline) // ' -remap'
! apply filter to sgh (and sgh30) in addition to height
    case ('-sgh')
      filter_sgh = .true.
      cmdline = trim(cmdline) // ' -sgh'
! not one of the above, must be output file name
    case default
      if (outbcfile .eq. ' ') then
        outbcfile = arg
      else
        write (6,*) 'Argument ', arg,' is not known'
        call usage_exit (' ')
      end if
      cmdline = trim(cmdline) // ' ' // trim(arg)
    end select
  end do
  
  if (outbcfile == ' ') then
    call usage_exit ('Must enter an output file name')
  end if
  
  if (gridfile == ' ') then
    call usage_exit ('Must enter gridfile name via -g arg (can use a model history file)')
  end if

  if (topofile == ' ') then
    call usage_exit ('Must enter topofile name via -t arg')
  end if

  if (filter_remap .and. filter_del2) then
     write(6,*)'Both filter_remap and filter_del2 set: using filter_remap'
  end if

  if (.not. filter_remap .and. .not. filter_del2) then
     write(6,*)'No filter being applied to height field'
     if (filter_sgh) call usage_exit ('Must filter height to filter sgh')
  end if

  if (landmfile == ' ') then
    call usage_exit ('Must enter landmfile name via -l arg')
  end if

! Open the grid file 
  ret = nf_open (trim(gridfile), nf_nowrite, gridid) 
  if (ret /= nf_noerr) then
    write(6,*)nf_strerror(ret)
    write(6,*)'Unable to open input file ', trim(gridfile), ' for writing'
    stop 999
  end if

! Get the grid dimensions from the grid file
  call wrap_inq_dimid  (gridid, 'lon',    londimid)
  call wrap_inq_dimlen (gridid, londimid, plon    )
  call wrap_inq_dimid  (gridid, 'lat',    latdimid)
  call wrap_inq_dimlen (gridid, latdimid, nlat    )
!
! Get longitude and latitude arrays for model grid.  
! If reduced grid, 2-d variable containing lon values for each lat is called "rlon".
! First allocate space for dynamic arrays now that sizes are known
!
  allocate (nlon(nlat))
  allocate (mlatcnts(nlat))
  allocate (mloncnts(plon,nlat))

  if (nf_inq_varid (gridid, 'nlon', nlonid) == nf_noerr) then
    if (nf_get_var_int (gridid, nlonid, nlon) /= nf_noerr) then
      write(6,*)'nf_get_var_int() failed for nlon'
      call endrun
    end if
    reduced_grid = .true.
  else
    nlon(:) = plon
  end if
    
  do j=1,nlat 
    if (nlon(j)<1 .or. nlon(j)>plon) then
      write(6,*)'nlon(',j,')=',nlon(j),' is invalid.'
      write(6,*)'Must be between 1 and ',plon
      call endrun
    end if
  end do

  call wrap_inq_varid (gridid, 'lat', latid)
  call wrap_get_var8 (gridid, latid, mlatcnts)

  if (nf_inq_varid (gridid, 'lon', lonid) == nf_noerr) then
     call wrap_get_var8 (gridid, lonid, mloncnts(1,1))
     do j=2,nlat
        mloncnts(:,j) = mloncnts(:,1)
     end do
  else
     call wrap_inq_varid (gridid, 'rlon', rlonid)
     call wrap_get_var8 (gridid, rlonid, mloncnts)
  end if

! Close the grid file
  if (nf_close (gridid) == nf_noerr) then
    write(6,*) 'close grid file ', trim(gridfile)
  else
    write(6,*) 'ERROR CLOSING NETCDF FILE ',trim(gridfile)
  end if
!
! Allocate space for variables
!
  allocate (sgh(plon,nlat))
  allocate (sgh30(plon,nlat))
  allocate (phis(plon,nlat))
  allocate (fland(plon,nlat))
  allocate (landm(plon,nlat))
!
! Determine model topographic height and 2 standard deviations
!
  call sghphis (plon, nlat, nlon, mlatcnts, mloncnts, topofile, &
                verbose, sgh, sgh30, have_sgh30, phis, fland)

! Do the terrain filter.
! Note: not valid if a reduced grid is used.
  if (filter_remap) then
     z_filter_type = 'remap'
     write(6,*)'Remapping terrain filtering'
! 7 and 3 are the recommended mapping accuracy settings 
     call map2f (plon, nlat, phis, 7, 3, .true.)
     if (filter_sgh) then
        s_filter_type = 'remap'
        write(6,*)'Filtering standard deviation'
        call map2f (plon, nlat, sgh, 7, 3, .true.)
        if(have_sgh30) call map2f(plon, nlat, sgh30, 7, 3, .true.)
     else
        s_filter_type = 'none (2x[1-2-1])'
        write(6,*)'Not filtering standard deviation'
     end if
  else if (filter_del2) then
     z_filter_type = 'del2'
     write(6,*) 'Del2 Terrain filtering'
     call sm2(plon, nlat, phis, plon/12, filter_coefficient)
     if (filter_sgh) then
        s_filter_type = 'del2'
        write(6,*)'Filtering standard deviation'
        call sm2(plon, nlat, sgh,  plon/12, filter_coefficient)
        if(have_sgh30) call sm2(plon, nlat, sgh30,  plon/12, filter_coefficient)
     else
        s_filter_type = 'none (2x[1-2-1])'
        write(6,*)'Not filtering standard deviation'
     end if
  else
     z_filter_type = 'none'
     s_filter_type = 'none (2x[1-2-1])'
  endif
!
! Adjustments to land fraction:
! 1. Extend land fraction for Ross Ice shelf
! 2. Set land fractions < .001 to 0.0
! 3. flag regions outside reduced grid
!
  do j=1,nlat
    do i=1,nlon(j)
!
! Overwrite FLAND flag as land for Ross ice shelf
      if (make_ross .and. mlatcnts(j) < -79.) then
        fland(i,j) = 1.
      end if

      if (fland(i,j) < .001_r8) fland(i,j) = 0.0

    end do
!
! Fill region outside reduced grid with flag values
    do i=nlon(j)+1,plon
      sgh(i,j)   = fillvalue
      if(have_sgh30) sgh30(i,j)   = fillvalue
      phis(i,j)  = fillvalue
      fland(i,j)   = fillvalue
      landm(i,j) = fillvalue
    end do
  end do
!
! Calculate LANDM field required by cloud water.  
!
!JR Replace original resolution-dependent calculation with interpolation.
!JR
!JR  call inimland (plon, nlat, nlon, mlatcnts, mloncnts, topofile, &
!JR                 verbose, make_ross, landm)
!
  call interplandm (plon, nlat, nlon, mlatcnts, mloncnts, &
                    landmfile, landm)

! Create NetCDF file for output
  ret = nf_create (outbcfile, NF_CLOBBER, foutid)
  if (ret .ne. NF_NOERR) call handle_error(ret)

! Create dimensions for output
  call wrap_def_dim (foutid, 'lon', plon, lonid)
  call wrap_def_dim (foutid, 'lat', nlat, latid)
  dim(1)=lonid
  dim(2)=latid

! Create latitude dimension variable for output
  ret = nf_def_var (foutid,'lat', NF_DOUBLE, 1, latid, latdimid)
  if (ret .ne. NF_NOERR) call handle_error(ret)
  call wrap_put_att_text (foutid,latdimid,'long_name', 'latitude')
  call wrap_put_att_text (foutid,latdimid,'units'    , 'degrees_north')

! Create longitude dimension variable for output
  if (.not.reduced_grid) then
     ret = nf_def_var (foutid,'lon', NF_DOUBLE, 1, lonid, londimid)
     if (ret .ne. NF_NOERR) call handle_error(ret)
     call wrap_put_att_text (foutid,londimid,'long_name', 'longitude')
     call wrap_put_att_text (foutid,londimid,'units'    , 'degrees_east')

! For reduced grid, add longitude limits (nlon) and lons (rlon)
  else
     ret = nf_def_var (foutid,'nlon', NF_INT, 1, lonid, londimid)
     if (ret .ne. NF_NOERR) call handle_error(ret)
     ret = nf_def_var (foutid,'rlon', NF_DOUBLE, 2, dim, rlonid)
     if (ret .ne. NF_NOERR) call handle_error(ret)
     call wrap_put_att_text (foutid,rlonid,'long_name', 'longitude')
     call wrap_put_att_text (foutid,rlonid,'units'    , 'degrees_east')
  end if

! Create variables for output
  ret = nf_def_var (foutid,'PHIS'        , NF_DOUBLE, 2, dim, phisid)
  if (ret .ne. NF_NOERR) call handle_error(ret)
  call wrap_put_att_double (foutid, phisid, '_FillValue', nf_double, 1, fillvalue)
  call wrap_put_att_double (foutid, phisid, 'missing_value', nf_double, 1, fillvalue)
  call wrap_put_att_text   (foutid, phisid, 'long_name' , 'surface geopotential')
  call wrap_put_att_text   (foutid, phisid, 'units'     , 'm2/s2')
  call wrap_put_att_text   (foutid, phisid, 'from_hires', 'true')
  call wrap_put_att_text   (foutid, phisid, 'filter'    , z_filter_type)

  ret = nf_def_var (foutid,'SGH'         , NF_DOUBLE, 2, dim, sghid)
  if (ret .ne. NF_NOERR) call handle_error(ret)
  call wrap_put_att_double (foutid, sghid, '_FillValue', nf_double, 1, fillvalue)
  call wrap_put_att_double (foutid, sghid, 'missing_value', nf_double, 1, fillvalue)
  call wrap_put_att_text   (foutid, sghid, 'long_name' , 'standard deviation of 10-min elevations')
  call wrap_put_att_text   (foutid, sghid, 'units'     , 'm')
  call wrap_put_att_text   (foutid, sghid, 'from_hires', 'true')
  call wrap_put_att_text   (foutid, sghid, 'filter'    , s_filter_type)

  if (have_sgh30) then
     ret = nf_def_var (foutid,'SGH30'       , NF_DOUBLE, 2, dim, sgh30id)
     if (ret .ne. NF_NOERR) call handle_error(ret)
     call wrap_put_att_double (foutid, sgh30id, '_FillValue', nf_double, 1, fillvalue)
     call wrap_put_att_double (foutid, sgh30id, 'missing_value', nf_double, 1, fillvalue)
     call wrap_put_att_text   (foutid, sgh30id, 'long_name' , 'standard deviation of elevation from 30s to 10m')
     call wrap_put_att_text   (foutid, sgh30id, 'units'     , 'm')
     call wrap_put_att_text   (foutid, sgh30id, 'from_hires', 'true')
     call wrap_put_att_text   (foutid, sgh30id, 'filter'    , s_filter_type)
  endif

  ret = nf_def_var (foutid,'LANDFRAC'    , NF_DOUBLE, 2, dim, landfid)
  if (ret .ne. NF_NOERR) call handle_error(ret)
  call wrap_put_att_double (foutid, landfid, '_FillValue', nf_double, 1, fillvalue)
  call wrap_put_att_double (foutid, landfid, 'missing_value', nf_double, 1, fillvalue)
  call wrap_put_att_text   (foutid, landfid, 'long_name' , 'gridbox land fraction')
  call wrap_put_att_text   (foutid, landfid, 'from_hires', 'true')

  ret = nf_def_var (foutid,'LANDM_COSLAT', NF_DOUBLE, 2, dim, landmid)
  if (ret .ne. NF_NOERR) call handle_error(ret)
  call wrap_put_att_double (foutid, landmid, '_FillValue', nf_double, 1, fillvalue)
  call wrap_put_att_double (foutid, landmid, 'missing_value', nf_double, 1, fillvalue)
  call wrap_put_att_text   (foutid, landmid, 'long_name' , &
                          'land ocean transition mask: ocean (0), continent (1), transition (0-1)')
  call wrap_put_att_text   (foutid, landmid, 'from_hires', 'true')

! Define history attribute.
  call DATE_AND_TIME(DATE=datestring)
  history = 'Written on date: ' // datestring // cmdline
  call wrap_put_att_text (foutid, nf_global, 'history', history)

! Define Ross Sea attribute
  if (make_ross) then
    write (6,*) 'Extending Ross ice shelf south of -79 degrees'
    call wrap_put_att_text (foutid, nf_global, 'make_ross', 'true')
  else
    write (6,*) 'Not doing anything special for Ross ice shelf'
    call wrap_put_att_text (foutid, nf_global, 'make_ross', 'false')
  end if

! Define source file attributes
  call wrap_put_att_text (foutid, nf_global, 'topofile', topofile)
  cmdlen = len_trim (gridfile)
  call wrap_put_att_text (foutid, nf_global, 'gridfile', gridfile)
  cmdlen = len_trim (landmfile)
  call wrap_put_att_text (foutid, nf_global, 'landmask', landmfile)


! End definition of netCDF file
  ret = nf_enddef (foutid)
  if (ret/=NF_NOERR) call handle_error (ret)


! Write data to file
  write(6,*) 'Writing surface quantities'

! Write dimension variables
  call wrap_put_var8 (foutid, latdimid, mlatcnts)
  if (.not.reduced_grid) then
     call wrap_put_var8 (foutid, londimid, mloncnts(:,1))
  else
     ret = nf_put_var_int (foutid, nlonid, nlon)
     if (ret/=NF_NOERR) call handle_error (ret)
     call wrap_put_vara8 (foutid, rlonid, start, count, mloncnts)
  end if

  start(:) = 1
  count(1) = plon
  count(2) = nlat
  count(3:) = 1

  call wrap_put_vara8 (foutid, sghid,   start, count, sgh)
  if(have_sgh30) call wrap_put_vara8 (foutid, sgh30id,   start, count, sgh30)
  call wrap_put_vara8 (foutid, phisid , start, count, phis)
  call wrap_put_vara8 (foutid, landfid, start, count, fland)
  call wrap_put_vara8 (foutid, landmid, start, count, landm)

  if (nf_close (foutid) == nf_noerr) then
    write(6,*) 'Successfully defined surface quantities on ', trim(outbcfile)
  else
    write(6,*) 'ERROR CLOSING NETCDF FILE ',trim(outbcfile)
  end if

  deallocate (nlon)
  deallocate (mlatcnts)
  deallocate (mloncnts)
  deallocate (sgh) 
  deallocate (sgh30)
  deallocate (phis)
  deallocate (fland)
  deallocate (landm)

  stop 0
end program fmain

subroutine usage_exit (arg)
  implicit none
  character*(*) arg
  
  if (arg /= ' ') write (6,*) arg
  write (6,*) 'Usage: definesurf -t topofile -g gridfile -l landmfile [-v] [-r] [-del2] [-remap] outfile'
  write (6,*) '       -v     verbose mode'
  write (6,*) '       -r     Do *not* extend Ross Ice Shelf as land ice'
  write (6,*) '       -del2  use del2 terrain filter (not a valid option for reduced grid)'
  write (6,*) '       -remap use remapping filter (not a valid option for reduced grid)'
  write (6,*) '       -sgh   filter sgh and sgh30 using same terrain filter'
  stop 999
end subroutine usage_exit
