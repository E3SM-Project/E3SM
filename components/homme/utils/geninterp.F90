program geninterp
implicit none

  integer, parameter :: END_OF_FILE = -1

  character(len=80)  :: instring
  character(len=80)  :: fname
  character(len=80)  :: gridfname
  character(len=80)  :: gridtype
  character(len=80)  :: batchname

  logical ex

  integer hobi     ! A switch for high-order OR bilinear interpolation
 
  integer zavg
  integer ios
  integer bunit,gunit

  bunit = 17
  gunit = 19
  batchname = "geninterp.batch"

  inquire(file=batchname,EXIST=ex)

  if (ex) then
     open(unit=bunit,file=batchname,form="FORMATTED",status="OLD")
     do while(.true.)

        ! ===============================
        ! Read in the grid file of lats and
        ! lons to interpolate
        ! values to...
        ! ===============================


        read(bunit,*,iostat=ios)instring
        if (ios == END_OF_FILE) then
            close(bunit)
            stop
        else
            gridfname = TRIM(ADJUSTL(instring))
        endif

        ! ================================
        ! Peak at grid file to determine
        ! the type...
        ! ================================
        if(gridfname.eq."NATIVE") then
           gridtype = "NATIVE"
        else
           open(unit=gunit,file=gridfname,form="FORMATTED")
           read(gunit,*)instring
           gridtype = TRIM(ADJUSTL(instring))
           close(gunit)
        
        end if
        ! ===============================
        ! Read in cube grid filename to 
        ! interpolate from...
        ! ===============================

        read(bunit,*,iostat=ios)instring
        if (ios == END_OF_FILE) then
           close(bunit)
           stop
        else
           fname = TRIM(ADJUSTL(instring))
        endif

        read(bunit,*,iostat=ios)zavg
        if (ios == END_OF_FILE) then
           close(bunit)
           stop
        endif

        ! Nair
        read(bunit,*,iostat=ios)hobi
        if (ios == END_OF_FILE) then
           close(bunit)
           stop
        endif
        if (gridtype == "NATIVE") then

        else if (gridtype == "latlong") then
           call cube_interp_latlong(gridfname,fname,hobi)
           if (zavg==1) call zonal_avg(gridfname,fname)
        else
           call cube_interp_generic(gridfname,fname,hobi)
        end if

     end do
     
  else

    ! =====================================
    ! Read in the grid to interpolate to from
    ! std in...
    ! =====================================

    print *,"Welcome to Geninterp!"
    print *,"---------------------"
    print *,"file name with new grid coordinates?"
    read(5,*)instring
    gridfname = TRIM(ADJUSTL(instring))

    ! ================================
    ! Peak at grid file to determine
    ! the type...
    ! ================================

    open(unit=gunit,file=gridfname,form="FORMATTED")
    read(gunit,*)instring
    gridtype = TRIM(ADJUSTL(instring))
    close(gunit)

    ! =============================================
    ! Read in the field file to interpolate from
    ! std in...
    ! =============================================


    print *,"enter file to interpolate"
    read(5,*)instring
    fname = TRIM(ADJUSTL(instring))

    print *,"zonally average (1=yes)"
    read(5,*)zavg

    print *,"Interpolation Type (0=high-order OR  1=bilinear)"
    read(5,*)hobi

    ! ==============================
    ! Interpolations default to 
    ! generic...
    ! ==============================

    if (gridtype == "latlong") then
       call cube_interp_latlong(gridfname,fname,hobi)
       if (zavg==1) call zonal_avg(gridfname,fname)
    else
       call cube_interp_generic(gridfname,fname,hobi)
    end if

  end if
  
  stop

end program geninterp

! ==============================================
! Generic interpolation routine from cube grid to
! list of lat-lon points...
! ==============================================

subroutine cube_interp_generic(gridfname,fname,itype)

! =========================================
  use kinds
! =========================================
  use math_constants
! =========================================
  use quadrature_mod
! =========================================
  use interpolate_mod
! =========================================
  use coordinate_systems_mod
! =========================================
implicit none

  character(len=*)gridfname
  character(len=*)fname

  ! ============================
  ! Local variable..
  ! ============================

  integer, intent(in)                                  ::  itype

  real (kind=real_kind), dimension(:,:,:), allocatable :: fld
  type (spherical_polar_t), dimension(:),  allocatable :: grid_new
  real (kind=real_kind), dimension(:),     allocatable :: fld_new
  real (kind=real_kind), dimension(:),     allocatable :: deglat
  real (kind=real_kind), dimension(:),     allocatable :: deglon
 

  type (quadrature_t)   :: gp
  type (interpolate_t)  :: interp

  character(len=80)  :: quadtype
  character(len=80)  :: newfname
  
  integer iface,i,j,k
  integer npts,ne,nface,nlev
  integer n1,lev
  
  integer ngrid
  integer fldunit,gunit,ounit

  fldunit = 8
  gunit   = 9
  ounit   = 10

  ! ===============================================
  ! Open grid file read in grid ...
  ! ===============================================

  open(unit=gunit,file=gridfname,form="formatted")

  ! =============================
  ! Open cube grid field...
  ! =============================

  open(unit=fldunit,file=fname,form="UNFORMATTED",access="SEQUENTIAL",status="old")

  read(fldunit)quadtype(1:2)
  read(fldunit)ne,npts,nface,nlev

  if (quadtype == "gs") then
     gp = gauss(npts)
  else
     gp = gausslobatto(npts)
  end if

   !! Interpolation type selection
    ! itype = 0  -> high-order
    ! itype = 1  -> bilinear 
  
  call interpolate_create(gp,interp)

     interp%itype = itype

  n1=npts*ne  
  allocate(fld(n1,n1,nface))

  ! =============================
  ! Open the output file ...
  ! =============================

  newfname = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(gridfname))
  open(unit=ounit,file=newfname,form="FORMATTED")   

  read(gunit,*)ngrid
  allocate(grid_new(ngrid))
  allocate(fld_new(ngrid))
  allocate(deglat(ngrid))
  allocate(deglon(ngrid))

  do i=1,ngrid
     read(gunit,*)grid_new(i)%lon,grid_new(i)%lat
  end do

  ! ===========================================
  ! Lat long data is assumed to be in radians
  ! ===========================================

  do i=1,ngrid
     deglon(i) = grid_new(i)%lon/degrad - 180.0D0
     deglat(i) = grid_new(i)%lat/degrad
  end do

  write(ounit)deglon
  write(ounit)deglat

  ! =================================
  ! Interpolate the old field to the 
  ! new grid...
  ! =================================

  do k=1,nlev
     read(fldunit)fld
     call gen_grid_assemble(grid_new,fld,fld_new,interp)     
     write(ounit)fld_new
  end do

  close(gunit)
  close(fldunit)
  close(ounit)

end subroutine cube_interp_generic

! ==============================================
! Special interpolation routine from cube grid to
! a regular lat-lon tensor product grid...
! ==============================================

subroutine cube_interp_latlong(gridfname,fnamein,itype)
! =========================================
  use kinds
! =========================================
  use math_constants
! =========================================
  use quadrature_mod
! =========================================
  use interpolate_mod
! =========================================
  use coordinate_systems_mod
! =========================================
  use typeSizes
! =========================================
  use netcdf
! =========================================
implicit none

  character(len=*)gridfname
  character(len=*)fnamein
  character(len=80)fname

  ! ============================
  ! Local variable..
  ! ============================

  integer, intent(in)                                  :: itype

  real (kind=real_kind), dimension(:,:,:), allocatable :: fld
  type (spherical_polar_t), dimension(:),  allocatable :: grid_new
  real (kind=real_kind), dimension(:),     allocatable :: fld_new
  real (kind=real_kind), dimension(:),     allocatable :: lat
  real (kind=real_kind), dimension(:),     allocatable :: lon
  real (kind=real_kind), dimension(:),     allocatable :: deglat
  real (kind=real_kind), dimension(:),     allocatable :: deglon

  type (quadrature_t)   :: gp
  type (interpolate_t)  :: interp

  character(len=80)  :: quadtype
  character(len=80)  :: newfname
  character(len=80)  :: tail
  character(len=80)  :: gridtype
  character(len=80)  :: fldname ! gotten by pruning numbers off end of field name
  
  integer iface,i,j,k
  integer npts,ne,nface,nlev
  integer n1,lev
  integer nlon,nlat
  integer fldunit,gunit,ounit, ltime, ierr

  ! netCDF related variables
  integer :: ncFileID
  integer :: latDimID, lonDimID, levDimID, timeDimID
  integer :: intVarID, latVarID, lonVarID, levVarID
  real (kind = real_kind) :: fillVal
  logical series 
  fldunit = 8
  gunit   = 9
!!!  ounit   = 10

  ! ===============================================
  ! Open grid file read in grid ...
  ! ===============================================

  open(unit=gunit,file=gridfname,form="FORMATTED")

  ! =============================
  ! Open cube grid field...
  ! =============================

  
  i = len(trim(fnamein))
  if(fnamein(i:i).eq.".") then
     print *, "interpolate series"
     fname = trim(fnamein)//"0"
     series = .true.
     ltime = 1
   else
      fname = fnamein
      series = .false.
  endif
  
  open(unit=fldunit,file=fname,form="UNFORMATTED",access="SEQUENTIAL",status="old")

  read(fldunit)quadtype(1:2)
  read(fldunit)ne,npts,nface,nlev

  print *,"quadtype:",quadtype(1:2)
  print *,"ne,npts,nface,nlev=",ne,npts,nface,nlev

  if (quadtype == "gs") then
     gp = gauss(npts)
  else
     gp = gausslobatto(npts)
  end if

  call interpolate_create(gp,interp)

       interp%itype = itype

  n1=npts*ne  
  allocate(fld(n1,n1,nface))

  ! =============================
  ! Open the output file ...
  ! =============================

!!!  newfname = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(gridfname))
!!!  open(unit=ounit,file=newfname,form="UNFORMATTED",access="SEQUENTIAL")   
  
  newfname = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(gridfname))//".nc"

  ! The filename fname is assumed to be composed of
  ! test.fldname.n, where test is the experiment and n is a frame number in a movie.
  ! fldname is gotten by stripping the leading and trailing names pieces off

  call tailstrip(fname,tail)
  print *,"tail=",tail
  call prune_trailing_numbers(tail,fldname)

  print *,"newfname=",newfname
  print *,"tail=",tail
  print *,"fld=",fldname

  ! Create the netCDF file
  call check(nf90_create(path = newfname, cmode = nf90_clobber, ncid = ncFileID))

  read(gunit,*)gridtype
  read(gunit,*)nlon
  read(gunit,*)nlat

  allocate(grid_new(nlon))
  allocate(fld_new(nlon))
  allocate(lat(nlat))
  allocate(deglat(nlat))
  allocate(lon(nlon))
  allocate(deglon(nlon))

  ! Define the dimensions
  call check(nf90_def_dim(ncid = ncFileID, name = "lat", len = nlat, dimid = latDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "lon", len = nlon, dimid = lonDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "lev", len = nlev, dimid = levDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "time", len =nf90_unlimited , dimid = timeDimID))

  ! Create variables and attributes
  call check(nf90_def_var(ncid = ncFileID, name = TRIM(ADJUSTL(fldname)), xtype = nf90_double,     &
                     dimids = (/ lonDimID, latDimID, levDimID, timeDimID /), varID = intVarID) )

  ! Use a 8-byte float constant, to match variable type
  fillVal = -9999.D0
  call check(nf90_put_att(ncFileID, intVarID,  "_FillValue", fillVal ) )
                      
  call check(nf90_def_var(ncFileID, "lat", nf90_double, dimids = latDimID, varID = latVarID) )
  call check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"))
  call check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"))

  call check(nf90_def_var(ncFileID, "lon", nf90_double, lonDimID, lonVarID) )
  call check(nf90_put_att(ncFileID, lonVarID, "long_name", "longitude"))
  call check(nf90_put_att(ncFileID, lonVarID, "units",     "degrees_east"))

  ! Global attributes
  call check(nf90_put_att(ncFileID, nf90_global, "history", "created by geninterp"))
  
  ! Leave define mode
  call check(nf90_enddef(ncfileID))

  ! ==========================================
  ! Lat and long are assumed to be in Radians
  ! ==========================================

  do i=1,nlon
     read(gunit,*)lon(i)
     deglon(i) = lon(i)/degrad - 180.0D0
  end do

  do j=1,nlat
     read(gunit,*)lat(j)
     deglat(j) = lat(j)/degrad
  end do

!!!  write(ounit)nlon
!!!  write(ounit)nlat
!!!  write(ounit)nlev
!!!  write(ounit)deglon
!!!  write(ounit)deglat
  
  ! Write the dimension variables
  call check(nf90_put_var(ncFileID, latVarId, deglat) )
  call check(nf90_put_var(ncFileID, lonVarId, deglon) )

  ! =================================
  ! Interpolate the old field to the 
  ! new grid...
  ! =================================

101 do k=1,nlev
     read(fldunit)fld
     do j=1,nlat
        grid_new(:)%lon = lon(:)
        grid_new(:)%lat = lat(j)
        call gen_grid_assemble(grid_new,fld,fld_new,interp)     
!!!        write(ounit)fld_new
        if(series) then
           call check(nf90_put_var(ncFileID, intVarID, values = fld_new, &
                start = (/ 1, j, k, ltime /), count = (/ nlon, 1, 1, 1 /)) )
        else
           call check(nf90_put_var(ncFileID, intVarID, values = fld_new, &
                start = (/ 1, j, k /), count = (/ nlon, 1, 1/)) )
        end if
     end do
  end do

  if(series) then
     close(fldunit)
     if(ltime<10) then
        write(fname,"(a,i1)") trim(fnamein),ltime
     else
        write(fname,"(a,i2)") trim(fnamein),ltime
     end if
     ltime=ltime+1
     print *, fname
     open(unit=fldunit,file=fname,form="UNFORMATTED",access="SEQUENTIAL",status="old",&
          iostat=ierr)
     if(ierr.eq.0) then
        read(fldunit)quadtype(1:2)
        read(fldunit)ne,npts,nface,nlev
        goto 101
     end if
  end if
  close(gunit)
  close(fldunit)
!!!  close(ounit)

  ! Close netCDF file
  call check(nf90_close(ncFileID))

end subroutine cube_interp_latlong

subroutine zonal_avg(gridfname,fname)
! =========================================
  use kinds
! =========================================
implicit none

  character(len=*)gridfname
  character(len=*)fname

  ! ============================
  ! Local variable..
  ! ============================

  real (kind=real_kind), dimension(:),     allocatable :: fld_new
  real (kind=real_kind), dimension(:),     allocatable :: deglat
  real (kind=real_kind), dimension(:),     allocatable :: deglon
  real (kind=real_kind), dimension(:),     allocatable :: sigmam
  real (kind=real_kind), dimension(:),     allocatable :: sigmah
  real (kind=real_kind) :: fld_zonal

  character(len=80)  :: infname
  character(len=80)  :: zfname

  integer i,j,k
  integer nlon,nlat,nlev
  integer zunit,iunit

  iunit = 21
  zunit = 22

  ! =============================
  ! Open the output file ...
  ! =============================

  infname = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(gridfname))
  open(unit=iunit,file=infname,form="UNFORMATTED",access="SEQUENTIAL")   

  ! ====================================
  ! Zonally average contents of infname 
  ! new grid...
  ! ====================================

  read(iunit)nlon
  read(iunit)nlat
  read(iunit)nlev

  allocate(fld_new(nlon))
  allocate(deglat(nlat))
  allocate(deglon(nlon))
  allocate(sigmam(nlev))
  allocate(sigmah(nlev+1))

  read(iunit)deglon
  read(iunit)deglat

  do k=1,nlev+1
     sigmah(k) = REAL(k-1,kind(0.0D0)) / REAL(nlev,kind(0.0D0))
  end do

  do k=1,nlev
     sigmam(k) = 0.50D0*(sigmah(k+1) + sigmah(k))
  end do

  zfname = TRIM(ADJUSTL(infname))//".zonal"
  open(unit=zunit,file=zfname,form="FORMATTED")

  do k=1,nlev
     do j=1,nlat
        read(iunit)fld_new
        fld_zonal=SUM(fld_new(:))/nlon
        write(zunit,10)deglat(j),sigmam(k),fld_zonal
 10 format(e22.15,1x,e22.15,1x,e22.15)
     end do
  end do

  close(iunit)
  close(zunit)

end subroutine zonal_avg

subroutine check(status)
  use typeSizes
  use netcdf
  integer, intent ( in) :: status
    
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
  end if
end subroutine check
