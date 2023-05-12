module mkprocdata_map

contains

subroutine mkmap(filei, fileo, fmap, ftemplate)

  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use fileutils
  use gridmapMod

  implicit none
  character(len=256), intent(in) :: filei   ! input  dataset
  character(len=256), intent(in) :: fileo   ! output mapped dataset
  character(len=256), intent(in) :: fmap    ! mapping file
  character(len=256), intent(in) :: ftemplate ! template file, containing lat & lon arrays desired in output file


  integer :: nDimensions    ! number of dimensions defined for this netCDF dataset.
  integer :: nVariables     ! number of variables defined for this netCDF dataset.
  integer :: output_n       ! number of variables in the output file that are obtained by
                            ! regridding variables from the input file
  integer :: nAttributes    ! number of global attributes defined for this netCDF dataset.
  integer :: unlimitedDimID ! ID of the unlimited dimension, if there is one
                            ! If no unlimited length dimension has been defined, 
                            ! -1 is returned.
  integer :: ier            ! returned error code
  integer :: nlen
  integer :: nlon, nlat
  integer :: dimid
  integer :: dimid_lndgrid
  integer :: dimid_lon, dimid_lat
  integer :: ncidi, ncido, ncidt
  integer :: varidi, varido
  integer :: varid_area
  integer :: dimlen
  integer :: xtype
  integer :: ndimsi,ndimso
  integer :: dimidsi(4)     ! dimension id array
  integer :: dimidso(4)     ! dimension id array
  integer :: nAtts          ! number of variable attributes
  integer :: n,nv,nt,nd,na  ! indices
  integer :: attlen

  ! input_ids & output_ids: arrays defining mapping between input & output variables:
  integer, dimension(:), allocatable :: input_ids
  integer, dimension(:), allocatable :: output_ids

  character(len=256):: locfn
  character(len=128):: dimname
  character(len=128):: varname
  character(len=128):: attname 
  character(len=128):: cattvalue
  real(r8)          :: dattvalue
  real(r4)          :: fattvalue
  integer           :: iattvalue
  type(gridmap_type):: tgridmap
  logical :: mapvar 

  !--------------------------------------------------------
  ! Read in mapping file - will have frac in it - and input and output domain
  !--------------------------------------------------------

  call getfil (fmap, locfn, 0)

  call handle_ncerr(nf90_open(locfn, NF90_NOWRITE, ncidi))

  call get_dimlen(ncidi, 'n_b', nlen)

  call handle_ncerr(nf90_close(ncidi))
  call gridmap_mapread(tgridmap, locfn)

  !--------------------------------------------------------
  ! Read in template file to get nlon & nlat
  !--------------------------------------------------------

  call getfil (ftemplate, locfn, 0)

  call handle_ncerr(nf90_open(locfn, NF90_NOWRITE, ncidt))

  call get_dimlen(ncidt, 'lon', nlon)
  call get_dimlen(ncidt, 'lat', nlat)

  call handle_ncerr(nf90_close(ncidt))

  write(6,*) 'nlon = ', nlon
  write(6,*) 'nlat = ', nlat

  if (nlon*nlat /= nlen) then
     write(6,*) 'must have nlon*nlat == nlen'
     write(6,*) 'nlon = ', nlon
     write(6,*) 'nlat = ', nlat
     write(6,*) 'nlen = ', nlen
     stop 1
  end if

  call getfil (filei, locfn, 0)

  !--------------------------------------------------------
  ! Create output file (put it in define mode)
  !--------------------------------------------------------

  call handle_ncerr(nf90_create(fileo, NF90_64BIT_OFFSET, ncido))

  !--------------------------------------------------------
  ! Define output dimensions - creating file puts it in define mode 
  !--------------------------------------------------------

  call handle_ncerr(nf90_open(locfn, NF90_NOWRITE, ncidi))
  call handle_ncerr(nf90_inquire(ncidi, nDimensions, nVariables, &
       nAttributes, unlimitedDimId))

  do nd = 1,nDimensions
     ! Determine input dimensions
     call handle_ncerr(nf90_inquire_dimension(ncidi, dimid=nd, name=dimname, len=dimlen))
 
     ! Define output variables
     ! Assume that input dimensions are time, lndgrid
     ! 2d lon,lat       <=> 1d lndgrid
     if (dimname == 'time') then
        call handle_ncerr(nf90_def_dim(ncido, name=dimname, len=nf90_unlimited, dimid=dimid))
     else
        if (trim(dimname) == 'lndgrid') then
           dimid_lndgrid= nd
           call handle_ncerr(nf90_def_dim(ncido, name='lon', len=nlon, dimid=dimid_lon)) 
           call handle_ncerr(nf90_def_dim(ncido, name='lat', len=nlat, dimid=dimid_lat)) 
        else
           call handle_ncerr(nf90_def_dim(ncido, name=dimname, len=dimlen, dimid=dimid))
        end if
     end if
     write(6,*)'n = ',nd,' dimname= ',trim(dimname)
  end do

  !--------------------------------------------------------
  ! Define output variables
  !--------------------------------------------------------

  allocate(input_ids(nVariables), output_ids(nVariables))

  ! Loop over input variables
  output_n = 0
  do nv = 1,nVariables

     ! Determine input variable
     call handle_ncerr(nf90_Inquire_Variable(ncid=ncidi, varid=nv, natts=natts, &
          name=varname, ndims=ndimsi, dimids=dimidsi, xtype=xtype)) 

     if (ignore_var(varname)) then
        write(6,*)'skipping writing out variable ',trim(varname)         

     else
        output_n = output_n + 1
        
        ! Determine output dimension ids
        if (dimidsi(1) == dimid_lndgrid) then
           ndimso = ndimsi + 1
           dimidso(1) = dimid_lon
           dimidso(2) = dimid_lat
           do n = 2,ndimsi
              call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(n), dimname))
              call handle_ncerr(nf90_inq_dimid(ncido, dimname, dimidso(n+1)))
           end do
        else
           ndimso = ndimsi
           do n = 1,ndimsi
              call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(n), dimname))
              call handle_ncerr(nf90_inq_dimid(ncido, dimname, dimidso(n))) 
           end do
        end if

        ! Define output variable and attributes 
        call handle_ncerr(nf90_def_var(ncido, name=varname, xtype=xtype, &
             dimids=dimidso(1:ndimso), varid=varido))

        input_ids(output_n) = nv
        output_ids(output_n) = varido

        do na = 1,natts
           call handle_ncerr(nf90_inq_attname(ncidi, varid=nv, attnum=na, name=attname))
           call handle_ncerr(nf90_inquire_attribute(ncidi, varid=nv, name=attname, &
                len=attlen, xtype=xtype))
           if (xtype == nf90_char) then
              call handle_ncerr(nf90_get_att(ncidi, varid=nv, name=attname, values=cattvalue))
              call handle_ncerr(nf90_put_att(ncido, varid=varido, name=attname, &
                   values=cattvalue(1:attlen)))
           else if (xtype == nf90_double) then
              call handle_ncerr(nf90_get_att(ncidi, varid=nv, name=attname, values=dattvalue))
              call handle_ncerr(nf90_put_att(ncido, varid=varido, name=attname, values=dattvalue))
           else if (xtype == nf90_float) then
              call handle_ncerr(nf90_get_att(ncidi, varid=nv, name=attname, values=fattvalue))
              call handle_ncerr(nf90_put_att(ncido, varid=varido, name=attname, values=fattvalue))
           else if (xtype == nf90_int) then
              call handle_ncerr(nf90_get_att(ncidi, varid=nv, name=attname, values=iattvalue))
              call handle_ncerr(nf90_put_att(ncido, varid=varido, name=attname, values=iattvalue))
           end if
        end do

     end if  ! .not. ignore_var
  end do
  ! now output_n is the number of output variables that are regridded from input variables

  call def_area(ncido, dimid_lon, dimid_lat, varid_area)

  ! End define mode
  call handle_ncerr(nf90_enddef(ncido))  

  !--------------------------------------------------------
  ! Read in variables and write them out
  !--------------------------------------------------------

  do nv = 1,output_n
     call outputvar(input_ids(nv), output_ids(nv), nlon, nlat, ncidi, ncido, tgridmap)
  end do

  call write_area(ncido, varid_area, tgridmap, nlon, nlat)

  call handle_ncerr(nf90_close(ncidi))
  call handle_ncerr(nf90_close(ncido))

  deallocate(input_ids, output_ids)

contains

   ! return true if we should ignore this variable, false otherwise
   ! (this allows us to ignore variables that aren't handled properly in the regridding)
   logical function ignore_var(varname)
      character(len=*), intent(in) :: varname

      ! We do NOT exclude 'area' because it may be useful to be able to refer to this in
      ! the regridded file. See also the 'gw' variable created through def_area /
      ! write_area for the actual grid cell areas of the output grid.
      ignore_var = (trim(varname) == 'date_written' .or. &
                    trim(varname) == 'time_written' .or. &
                    trim(varname) == 'lon' .or. &
                    trim(varname) == 'lat' .or. &
                    trim(varname) == 'landmask' .or. &
                    trim(varname) == 'pftmask')
   end function ignore_var

   ! return true if we should weight this variable by landfrac, false otherwise
   ! note that most variables are weighted by landfrac, so this function just lists those
   ! variables that are NOT weighted by landfrac
   logical function weight_by_landfrac(varname)
      character(len=*), intent(in) :: varname

      weight_by_landfrac = .not. (trim(varname) == 'area' .or. &
                                  trim(varname) == 'landfrac')
   end function weight_by_landfrac

   ! Get the length of a dimension from an open netcdf file
   subroutine get_dimlen(ncid, dimname, len)
      use netcdf

      implicit none

      ! Subroutine arguments
      integer, intent(in) :: ncid ! ID of an open netcdf file
      character(len=*), intent(in) :: dimname ! name of dimension of interest
      integer, intent(out) :: len ! length of dimension

      ! Local variables
      integer :: dimid
      integer :: ier

      call handle_ncerr(nf90_inq_dimid(ncid, dimname, dimid))
      call handle_ncerr(nf90_inquire_dimension(ncid, dimid, len=len))

   end subroutine get_dimlen
 

  subroutine outputvar(varidi, varido, nlon, nlat, ncidi, ncido, tgridmap)

    integer, intent(in) :: varidi
    integer, intent(in) :: varido
    integer, intent(in) :: nlon
    integer, intent(in) :: nlat
    integer, intent(in) :: ncidi
    integer, intent(in) :: ncido
    type(gridmap_type), intent(inout) :: tgridmap

    integer :: len1,len2,len3,len4
    integer :: n1,n2,n3,n4
    integer :: dimidsi(4)     
    integer :: dimidso(4)     
    real(r8), allocatable :: rarrayi(:),rarrayo(:)
    character(len=128):: dimname
    character(len=128):: varname
    character(len=128):: attname 
    integer :: varid_landfrac
    logical :: first_time = .true.
    real(r8), allocatable, save :: landfraci(:) 
    real(r8), allocatable, save :: src_mask(:)
    integer :: vid
    real(r8):: spval = 1.e36 

    ! Allocate landfraci and src_mask; determine landfrac on input file
    if (first_time) then
       dimname = 'lndgrid'
       call handle_ncerr(nf90_inq_dimid(ncidi, dimname, dimidsi(1)))
       call handle_ncerr(nf90_inquire_dimension(ncidi, dimid=dimidsi(1), len=len1))
       allocate(landfraci(len1))
       allocate(src_mask(len1))
       call handle_ncerr(nf90_inq_varid(ncid=ncidi, name='landfrac', varid=varid_landfrac))
       call handle_ncerr(nf90_get_var(ncidi, varid_landfrac, landfraci))
       first_time = .false.
    end if

    call handle_ncerr(nf90_Inquire_Variable(ncid=ncidi, varid=varidi, &
         name=varname, ndims=ndimsi, dimids=dimidsi, xtype=xtype))
    write(6,*)'varidi = ',varidi,' varido = ', varido, ' varname= ',trim(varname),' ndimsi= ',ndimsi 
    
    call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(1), &
         len=len1, name=dimname))
    allocate(rarrayi(len1))
    if (trim(dimname)=='lndgrid') then
       mapvar = .true.
    else
       mapvar = .false.
    end if

    call handle_ncerr(nf90_Inquire_Variable(ncid=ncido, varid=varido, &
         name=varname, ndims=ndimso, dimids=dimidso, xtype=xtype))
    if (mapvar) then
       call handle_ncerr(nf90_inquire_dimension(ncido, dimidso(1), &
            len=len1, name=dimname))
       call handle_ncerr(nf90_inquire_dimension(ncido, dimidso(2), &
            len=len2, name=dimname))
       allocate(rarrayo(len1*len2))
    else
       len1 = size(rarrayi)
       allocate(rarrayo(len1))
    end if

    ! src_mask will give the relative weight of each grid cell. i.e., if two source grid
    ! cells have the same area of overlap with a destination grid cell, then their
    ! relative weightings are given by their src_mask values.
    ! For most grid cells, this relative weighting is landfrac, but there are a few
    ! fields for which we do not want to weight by landfrac.
    ! Note that, for some fields (fields that are NA over certain landunits) we should
    ! really also be weighting by an additional factor saying the fraction of land area
    ! over which that field applies; but more metadata need to be added to the clm
    ! history files before that will be possible.
    if (mapvar .and. weight_by_landfrac(varname)) then
       src_mask(:) = landfraci(:)
    else
       ! note that, if we get here because mapvar is false, we currently don't use
       ! src_mask; but we set it to 1 anyway to be safe
       src_mask(:) = 1
    end if
    
    if (ndimsi == 1) then
       call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(1), len=len1))
       call handle_ncerr(nf90_get_var(ncidi, varidi, rarrayi))
       if (mapvar) then
          if (xtype == nf90_int) then
             where (rarrayi == -9999) src_mask = 0
             call gridmap_areaave(tgridmap, rarrayi, rarrayo, src_mask, &
                  spval=-9999._r8)
          else
             where (rarrayi == spval) src_mask = 0
             call gridmap_areaave(tgridmap, rarrayi, rarrayo, src_mask)
          end if
          call handle_ncerr(nf90_put_var(ncido, varido, &
               reshape(rarrayo,(/nlon,nlat/))))
       else
          rarrayo(:)= rarrayi(:)
          call handle_ncerr(nf90_put_var(ncido, varido, rarrayo))
       end if
    end if
    
    if (ndimsi == 2) then
       call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(1), len=len1))
       call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(2), len=len2))
       do n2 = 1,len2
          call handle_ncerr(nf90_get_var(ncidi, varidi, rarrayi, &
               start=(/1,n2/), count=(/len1,1/)))
          if (mapvar) then
             if (xtype == nf90_int) then
                where (rarrayi == -9999) src_mask = 0
                call gridmap_areaave(tgridmap, rarrayi, rarrayo, src_mask, &
                     spval=-9999._r8)
             else
                where (rarrayi == spval) src_mask = 0
                call gridmap_areaave(tgridmap, rarrayi, rarrayo, src_mask)
             end if

             call handle_ncerr(nf90_put_var(ncido, varido, &
                  reshape(rarrayo,(/nlon,nlat/)), &
                  start=(/1,1,n2/), count=(/nlon,nlat/)))
          else
             rarrayo(:)= rarrayi(:)
             call handle_ncerr(nf90_put_var(ncido, varido, rarrayo, &
                  start=(/1,n2/), count=(/len1,1/)))
           end if
        end do
     end if
     
     if (ndimsi == 3) then
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(1), len=len1))
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(2), len=len2))
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(3), len=len3))
        do n2 = 1,len2
        do n3 = 1,len3           
           call handle_ncerr(nf90_get_var(ncidi, varidi, rarrayi, &
                start=(/1,n2,n3/), count=(/len1,1,1/)))
           if (mapvar) then
              where (rarrayi == spval) src_mask = 0
              call gridmap_areaave(tgridmap, rarrayi, rarrayo, src_mask)
              call handle_ncerr(nf90_put_var(ncido, varido, &
                   reshape(rarrayo,(/nlon,nlat/)), &
                   start=(/1,1,n2,n3/), count=(/nlon,nlat/)))
           else
              rarrayo(:)= rarrayi(:)
              call handle_ncerr(nf90_put_var(ncido, varido, rarrayo, &
                   start=(/1,n2,n3/), count=(/len1,1,1/)))
           end if
        end do
        end do
     end if

     if (ndimsi == 4) then
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(1), len=len1))
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(2), len=len2))
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(3), len=len3))
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(4), len=len4))
        do n2 = 1,len2
        do n3 = 1,len3           
        do n4 = 1,len4           
           call handle_ncerr(nf90_get_var(ncidi, varidi, rarrayi, &
                start=(/1,n2,n3,n4/), count=(/len1,1,1,1/)))
           if (mapvar) then
              src_mask(:) = 0.
              where (rarrayi == spval) src_mask = spval
              call gridmap_areaave(tgridmap, rarrayi, rarrayo)
              call handle_ncerr(nf90_put_var(ncido, varido, &
                   reshape(rarrayo,(/nlon,nlat/)), &
                   start=(/ 1,1,n2,n3,n4 /), count=(/nlon,nlat,1,1,1/)))
           else
              rarrayo(:)= rarrayi(:)
              call handle_ncerr(nf90_put_var(ncido, varido, rarrayo, &
                   start=(/1,n2,n3,n4/), count=(/len1,1,1,1/)))
           end if
        end do
        end do
        end do
     end if

     deallocate(rarrayi)
     deallocate(rarrayo)

  end subroutine outputvar


  ! Define the area variable 'gw' in the output file, return the variable id of this
  ! variable. This variable will give the actual grid cell areas of the output grid. The
  ! name ('gw') is consistent with the similar variable in the cam homme regridding
  ! routine.
  subroutine def_area(ncido, dimid_lon, dimid_lat, varid)
     implicit none

     ! Subroutine arguments
     integer, intent(in) :: ncido      ! netcdf id of output netcdf file, in define mode
     integer, intent(in) :: dimid_lon  ! id of lon dimension in the file given by ncido
     integer, intent(in) :: dimid_lat  ! id of lat dimension in the file given by ncido
     integer, intent(out) :: varid     ! variable id of the newly-created area variable
     
     call handle_ncerr(nf90_def_var(ncido, name='gw', xtype=NF90_FLOAT, &
          dimids=(/dimid_lon,dimid_lat/), varid=varid))

     call handle_ncerr(nf90_put_att(ncido, varid=varid, name='long_name', &
          values='grid cell areas'))

     call handle_ncerr(nf90_put_att(ncido, varid=varid, name='units', &
          values='km^2'))

  end subroutine def_area

  ! Write the area variable to the output file. This assumes that def_area has already
  ! been called, and varid is the variable id returned by that subroutine.
  subroutine write_area(ncido, varid, tgridmap, nlon, nlat)
     use constMod, only : re_km  ! radius of earth (km)

     implicit none

     ! Subroutine arguments
     integer, intent(in) :: ncido   ! netcdf id of output netcdf file, open for writing
     integer, intent(in) :: varid   ! variable id of the area variable in which we'll put data
     type(gridmap_type), intent(in) :: tgridmap ! grid information; we're interested in
                                                ! the area_dst component, assumed to be in
                                                ! square radians
     integer, intent(in) :: nlon    ! number of longitudes in the output map
     integer, intent(in) :: nlat    ! number of latitudes in the output map

     ! Local variables
     real(r8), dimension(:,:), allocatable :: area_2d  ! grid cell areas, on a 2-d array,
                                                       ! in km^2


     allocate(area_2d(nlon,nlat))
     area_2d = reshape(tgridmap%area_dst,(/nlon,nlat/)) * re_km**2

     call handle_ncerr(nf90_put_var(ncido, varid, area_2d))

     deallocate(area_2d)
  end subroutine write_area

  subroutine handle_err(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
       write(6,*)trim(nf90_strerror(status))
       stop 1
    end if
  end subroutine handle_err

  subroutine handle_ncerr(ret)
     integer         , intent(in) :: ret  ! return code from netCDF library routine
     if ( ret .ne. NF90_NOERR ) then
        write(6,*) nf90_strerror( ret )
        stop 1
     endif
   end subroutine handle_ncerr

end subroutine mkmap

end module mkprocdata_map
