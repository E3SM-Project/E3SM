!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_example_clim.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glint_example_clim

  ! Subroutines used to initialize and compute forcing for the glint_example climate driver

  use glimmer_global, only: dp, fname_length
  use glint_global_grid

  implicit none

  type glex_climate

     ! Mass-balance coupling timing parameters --------------------------
     integer                 :: total_years  = 10  ! Length of run in years
     integer                 :: climate_tstep = 6  ! Climate time-step in hours

     ! Filenames --------------------------------------------------------
     character(fname_length) :: precip_file = '' !*FD Name of precip file
     character(fname_length) :: stemp_file  = '' !*FD Name of surface temp file
     character(fname_length) :: orog_file   = '' !*FD Name of orography file

     ! Variable names ---------------------------------------------------
     character(fname_length) :: precip_varname = '' !*FD precip variable name
     character(fname_length) :: stemp_varname  = '' !*FD temperature variable name
     character(fname_length) :: orog_varname   = '' !*FD orography variable name

     ! Arrays for holding climatology -----------------------------------
     real(dp),dimension(:,:),  pointer :: orog_clim     => null()  !*FD Orography
     real(dp),dimension(:,:,:),pointer :: precip_clim   => null()  !*FD Precip
     real(dp),dimension(:,:,:),pointer :: surftemp_clim => null()  !*FD Surface temperature

     ! Grid variables ---------------------------------------------------
     type(global_grid) :: clim_grid

     ! Time variables ---------------------------------------------------
     real(dp),dimension(:),pointer :: pr_time => null() !*FD Time in precip climatology
     real(dp),dimension(:),pointer :: st_time => null() !*FD Time in surftemp climatology

     ! Other parameters -------------------------------------------------
     integer  :: days_in_year  = 365
     integer  :: hours_in_year = 365*24
     real(dp) :: precip_scale  = 1.d0  ! Factor for scaling precip
     logical  :: temp_in_kelvin=.true. ! Set if temperature field is in Kelvin

     !NOTE: The glint_example driver assumes we will read in precip, surface air temp,
     !       and surface orography as required for a PDD scheme.
     !      But this module includes a subroutine (compute_smb_gcm) that computes a crude SMB 
     !       based on these inputs, so we can run Glint in SMB mode instead of PDD mode. 
     logical :: gcm_smb = .false.   ! if true, pass SMB to glint instead of PDD info

  end type glex_climate

  interface read_ncdf
     module procedure read_ncdf_1d,read_ncdf_2d,read_ncdf_3d
  end interface

contains

  subroutine glex_clim_init(params,filename)

    use glimmer_config
    use glimmer_log

    type(glex_climate) :: params   !*FD Climate parameters
    character(*)       :: filename !*FD config filename

    type(ConfigSection),pointer :: config !*FD structure holding sections of configuration file   
    type(global_grid) :: pgrid,sgrid,ogrid 
    character(20) :: sttu,prtu ! Units
    integer :: ierr,i

    call ConfigRead(filename,config)
    call glex_clim_readconfig(params,config)
    call glex_clim_printconfig(params)
    call CheckSections(config)

    ! Read in global grids

    call read_ncdf_ggrid(params%precip_file,pgrid)
    call read_ncdf_ggrid(params%stemp_file, sgrid)
    call read_ncdf_ggrid(params%orog_file,  ogrid)

    ! Check all grids are the same, and copy
    
    call check_ggrids(pgrid,sgrid,ogrid)
    params%clim_grid=pgrid

    ! Read in time axes

    call read_ncdf(params%precip_file,'time',params%pr_time,units=prtu)
    call read_ncdf(params%stemp_file, 'time',params%st_time,units=sttu)

    ! Scale as fractions of a year if necessary
    
    call scale_time(params,params%pr_time,prtu)
    call scale_time(params,params%st_time,sttu)

    ! Read in data

    call read_ncdf(params%precip_file,params%precip_varname,params%precip_clim)
    call read_ncdf(params%stemp_file, params%stemp_varname, params%surftemp_clim)
    call read_ncdf(params%orog_file,  params%orog_varname,  params%orog_clim)

    ! Scale precip 
    params%precip_clim = params%precip_clim * params%precip_scale

    ! Convert temps to degrees C if necessary
    if (params%temp_in_kelvin) then
       params%surftemp_clim=params%surftemp_clim-273.15
    end if

  end subroutine glex_clim_init

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glex_clim_readconfig(params,config)

    use glimmer_config
    use glimmer_log

    type(glex_climate)           :: params !*FD Climate parameters
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file   
    type(ConfigSection), pointer :: section

    call GetSection(config,section,'GLEX climate')
    if (associated(section)) then
       call GetValue(section,'days_in_year',params%days_in_year)
       call GetValue(section,'total_years',params%total_years)
       call GetValue(section,'climate_tstep',params%climate_tstep)
       call GetValue(section,'gcm_smb',params%gcm_smb)

       params%hours_in_year = params%days_in_year*24
    end if

    call GetSection(config,section,'GLEX precip')
    if (associated(section)) then
       call GetValue(section,'filename',params%precip_file)
       call GetValue(section,'variable',params%precip_varname)
       call GetValue(section,'scaling',params%precip_scale)
    end if

    call GetSection(config,section,'GLEX temps')
    if (associated(section)) then
       call GetValue(section,'filename',params%stemp_file)
       call GetValue(section,'variable',params%stemp_varname)
       call GetValue(section,'kelvin',  params%temp_in_kelvin)
    end if

    call GetSection(config,section,'GLEX orog')
    if (associated(section)) then
       call GetValue(section,'filename',params%orog_file)
       call GetValue(section,'variable',params%orog_varname)
    end if

    if (params%precip_file=='') &
       call write_log('GLINT Example: precip filename must be supplied',GM_FATAL)
    if (params%stemp_file=='') &
       call write_log('GLINT Example: temperature filename must be supplied',GM_FATAL)
    if (params%orog_file=='') &
       call write_log('GLINT Example: orography filename must be supplied',GM_FATAL)
    if (params%precip_varname=='') &
       call write_log('GLINT Example: precip variable must be specified',GM_FATAL)
    if (params%stemp_varname=='') &
       call write_log('GLINT Example: temperature variable must be specified',GM_FATAL)
    if (params%orog_varname=='') &
       call write_log('GLINT Example: orography variable must be specified',GM_FATAL)

  end subroutine glex_clim_readconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glex_clim_printconfig(params)

    use glimmer_log

    type(glex_climate) :: params
    character(100) :: message
    
    call write_log('GLINT Example configuration')
    call write_log('---------------------------')

    call write_log('Precip: '//trim(params%precip_varname)//' in file '// &
         trim(params%precip_file))
    call write_log('Surface temperature: '//trim(params%stemp_varname)//' in file '// &
         trim(params%stemp_file))
    call write_log('Orography: '//trim(params%orog_varname)//' in file '// &
         trim(params%orog_file))

    if (params%temp_in_kelvin) then
       call write_log('Temperatures in Kelvin')
    else
       call write_log('Temperatures in degC')
    end if
    
    if (params%precip_scale /= 1.d0) then
       write(message,*)'Precipitation scaled by ', params%precip_scale
       call write_log(message)
    end if

    if (params%gcm_smb) then
       call write_log ('Will pass surface mass balance (not PDD info) to Glint')
    endif
       
    call write_log('')

  end subroutine glex_clim_printconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_ggrid(fname,ggrid)

    use netcdf
    use glimmer_log

   !*FD Constructs a global grid data type from file

    character(*),intent(in)       :: fname
    type(global_grid),intent(out) :: ggrid

    integer :: ncerr     ! NetCDF error 
    integer :: ncid      ! NetCDF file id
    integer :: varid     ! NetCDF variable id
    integer :: ndims     ! Number of dimensions

    integer :: lon_id,lon_nd,lat_id,lat_nd
    character(30) :: lon_varn,lat_varn

    integer :: nx,ny
    integer,dimension(1) :: lldimids
    real(dp),dimension(:),allocatable :: lons,lats

    ! Open file
    ncerr = nf90_open(fname,0,ncid)
    call handle_err(ncerr,__LINE__)

    ! Look for desired dimension names - the standard name attribute is the
    ! place to look.
    call ncdf_find_var(ncid,(/'longitude'/),lon_varn,lon_id,lon_nd)
    call ncdf_find_var(ncid,(/'latitude' /),lat_varn,lat_id,lat_nd)

    ! Check they're only 1D arrays
    if (lon_nd/=1 .or. lat_nd/=1) &
         call write_log('Latitude and longitude variables must be 1D',GM_FATAL)

    ! Find out the sizes
    ncerr = nf90_inquire_variable(ncid,lon_id,dimids=lldimids)
    call handle_err(ncerr,__LINE__)
    ncerr = nf90_inquire_dimension(ncid,lldimids(1),len=nx)
    call handle_err(ncerr,__LINE__)
    ncerr = nf90_inquire_variable(ncid,lat_id,dimids=lldimids)
    call handle_err(ncerr,__LINE__)
    ncerr = nf90_inquire_dimension(ncid,lldimids(1),len=ny)
    call handle_err(ncerr,__LINE__)

    ! Allocate temporary arrays
    allocate(lons(nx),lats(ny))

    ! Read in lats and lons
    ncerr = nf90_get_var(ncid,lon_id,lons)
    call handle_err(ncerr,__LINE__)
    ncerr = nf90_get_var(ncid,lat_id,lats)
    call handle_err(ncerr,__LINE__)

    ! NB we are ignoring cell boundaries here.
    ! Construct global grid type
    call new_global_grid(ggrid,real(lons,dp),real(lats,dp),correct=.false.)

    ! Close file
    ncerr = nf90_close(ncid)
    call handle_err(ncerr,__LINE__)

  end subroutine read_ncdf_ggrid

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  subroutine check_ggrids(g1,g2,g3)

    use glimmer_log

    !*FD Compares three grids to make sure they are all the same

    type(global_grid),intent(in) :: g1,g2,g3
    logical :: fail

    fail=.false.

    if (g1%nx/=g2%nx .or. g1%nx/=g3%nx) fail = .true.
    if (g1%ny/=g2%ny .or. g1%ny/=g3%ny) fail = .true.

    if (any(g1%lons/=g2%lons) .or. any(g1%lons/=g3%lons)) fail = .true.
    if (any(g1%lats/=g2%lats) .or. any(g1%lats/=g3%lats)) fail = .true.

    if (fail) &
         call write_log('GLINT Example: All three grids must be the same',GM_FATAL)

  end subroutine check_ggrids

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  subroutine ncdf_find_var(ncid,stdnames,varname,varid,ndims,fatal)

    use netcdf
    use glimmer_log

    !*FD Returns the name and id of the first found variable which has
    !*FD the requested standard name attribute

    integer,intent(in) :: ncid !*FD ID of an open netcdf file
    character(*),dimension(:),intent(in)  :: stdnames !*FD standard names sought
    character(*),intent(out) :: varname !*FD variable name
    integer,intent(out) :: varid !*FD variable ID
    integer,intent(out) :: ndims !*FD Number of dimensions of variable
    logical,intent(in),optional :: fatal !*FD set true to halt if name not found

    integer :: nvars,iv,natts,ia,nd,nsn,in
    integer :: ncerr
    character(50) :: an,sn
    character(100) :: message
    logical :: ft

    if (present(fatal)) then
       ft = fatal
    else
       ft = .true.
    end if

    nsn=size(stdnames)

    ncerr = nf90_inquire(ncid,nVariables=nvars)
    call handle_err(ncerr,__LINE__)

    ! Loop over variables
    do iv = 1,nvars
       ncerr = nf90_inquire_variable(ncid,iv,varname,ndims=nd,nAtts=natts)
       call handle_err(ncerr,__LINE__)

       ! Loop over attributes
       do ia = 1,natts
          an = ''
          ncerr = nf90_inq_attname(ncid, iv, ia, an)
          call handle_err(ncerr,__LINE__)

          ! If standard name, get value and check 
          ! against targets in turn
          if (trim(an)=='standard_name' .or. trim(an)=='long_name') then
             sn = ''
             ncerr = nf90_get_att(ncid, iv, an, sn)
             call handle_err(ncerr,__LINE__)
             do in = 1,nsn
                if (trim(sn)==trim(stdnames(in))) then
                   varid = iv
                   ndims = nd
                   return
                end if
             end do
          end if

       end do
    end do

    ! If we get to here, we've failed to find the right std name anywhere.
    if (ft) then
       message = 'Failed to find standard names: '
       do in = 1,nsn
          message = trim(message)//' '//trim(stdnames(in))
          if (in/=nsn) message = trim(message)//','
       end do
       call write_log(message,GM_FATAL)
    end if

  end subroutine ncdf_find_var

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! READ_NCDF routines
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_1d(filename,varname,array,units)

    use netcdf
    use glimmer_log

    character(*), intent(in)      :: filename,varname
    real(dp),dimension(:),pointer :: array
    character(*),optional,intent(out) :: units

    real(dp),dimension(:),allocatable :: dim1

    integer  :: ncerr     ! NetCDF error 
    integer  :: ncid      ! NetCDF file id
    integer  :: varid     ! NetCDF variable id
    integer  :: ndims     ! Number of dimensions
    real(dp) :: offset,scale
    integer,      dimension(1) :: dimids,dimlens
    character(20),dimension(1) :: dimnames
    character(100) :: message
    integer :: u_attlen

    if (associated(array)) deallocate(array)
    offset = 0.d0
    scale  = 1.d0

    call read_ncdf_findvar(filename,ncid,varid,ndims,varname)

    ! If not a 1d variable, flag and error and exit ----

    if (ndims/=1) then
       write(message,*)'NetCDF: Requested variable has ',ndims,' dimensions, 1 required'
       call write_log(message,GM_FATAL)
    end if

    call read_ncdf_dimnames(ncid,varid,ndims,dimids,dimlens,dimnames)

    ! Allocate output and dimension arrays -------------

    allocate(array(dimlens(1)))

    ! Retrieve variable contents -----------------------

    ncerr=nf90_get_var(ncid, varid, array)
    call handle_err(ncerr,__LINE__)

    call read_ncdf_scaling(ncid,varid,offset,scale)

    array=offset+(array*scale)

    ! Find units if necessary
    if (present(units)) then
       ncerr=nf90_inquire_attribute(ncid, varid, 'units', len=u_attlen)
       call handle_err(ncerr,__LINE__)
       ncerr=nf90_get_att(ncid, varid, 'units', units)
       call handle_err(ncerr,__LINE__)
       units=units(1:u_attlen)
    end if

  end subroutine read_ncdf_1d

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_2d(filename,varname,array)

    use netcdf

    character(*), intent(in)        :: filename,varname
    real(dp),dimension(:,:),pointer :: array

    integer  :: ncerr     ! NetCDF error 
    integer  :: ncid      ! NetCDF file id
    integer  :: varid     ! NetCDF variable id
    integer  :: ndims     ! Number of dimensions
    real(dp) :: offset,scale
    integer,      dimension(2) :: dimids,dimlens
    character(20),dimension(2) :: dimnames

    if (associated(array)) deallocate(array)
    offset = 0.d0
    scale  = 1.d0

    call read_ncdf_findvar(filename,ncid,varid,ndims,varname)

    ! If not a 2d variable, flag and error and exit ----

    if (ndims /= 2) then
       print*,'NetCDF: Requested variable only has ',ndims,' dimensions'
       stop
    end if

    call read_ncdf_dimnames(ncid,varid,ndims,dimids,dimlens,dimnames)

    ! Allocate output -------------

    allocate(array(dimlens(1),dimlens(2)))

    ! Retrieve variable contents -----------------------

    ncerr=nf90_get_var(ncid, varid, array)
    call handle_err(ncerr,__LINE__)

    call read_ncdf_scaling(ncid,varid,offset,scale)

    array = offset + (array*scale)

  end subroutine read_ncdf_2d

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_3d(filename,varname,array)

    use netcdf

    character(*), intent(in)          :: filename,varname
    real(dp),dimension(:,:,:),pointer :: array

    integer  :: ncerr     ! NetCDF error 
    integer  :: ncid      ! NetCDF file id
    integer  :: varid     ! NetCDF variable id
    integer  :: ndims     ! Number of dimensions
    real(dp) :: offset,scale
    integer,      dimension(3) :: dimids,dimlens
    character(20),dimension(3) :: dimnames

    if (associated(array)) deallocate(array)
    offset = 0.d0
    scale  = 1.d0

    call read_ncdf_findvar(filename,ncid,varid,ndims,varname)

    ! If not a 3d variable, flag and error and exit ----

    if (ndims /= 3) then
       print*,'NetCDF: Requested variable only has ',ndims,' dimensions'
       stop
    end if

    call read_ncdf_dimnames(ncid,varid,ndims,dimids,dimlens,dimnames)

    ! Allocate output -------------

    allocate(array(dimlens(1),dimlens(2),dimlens(3)))

    ! Retrieve variable contents -----------------------

    ncerr=nf90_get_var(ncid, varid, array)
    call handle_err(ncerr,__LINE__)

    call read_ncdf_scaling(ncid,varid,offset,scale)

    array = offset + (array*scale)

  end subroutine read_ncdf_3d

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_findvar(filename,ncid,varid,ndims,varname)

    use netcdf

    character(*) :: filename
    integer :: ncid,varid,ndims
    character(*) :: varname
    integer :: ncerr

    ! Open file

    ncerr=nf90_open(filename,0,ncid)
    call handle_err(ncerr,__LINE__)

    ! Find out the id of variable and its dimensions

    ncerr=nf90_inq_varid(ncid,varname,varid)
    call handle_err(ncerr,__LINE__)
    ncerr=nf90_inquire_variable(ncid, varid, ndims=ndims)
    call handle_err(ncerr,__LINE__)

  end subroutine read_ncdf_findvar

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_dimnames(ncid,varid,ndims,dimids,dimlens,dimnames)

    use netcdf

    integer :: ncid,varid,ndims
    integer,dimension(:) :: dimids,dimlens
    character(*),dimension(:) :: dimnames

    integer :: ncerr,i

    ! Get dimensions ids 

    ncerr=nf90_inquire_variable(ncid, varid, dimids=dimids)
    call handle_err(ncerr,__LINE__)

    ! Retrieve dimension names

    do i=1,ndims
       ncerr=nf90_inquire_dimension(ncid, dimids(i),name=dimnames(i),len=dimlens(i))
       call handle_err(ncerr,__LINE__)
    end do

  end subroutine read_ncdf_dimnames

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_scaling(ncid,varid,offset,scale)

    use netcdf

    integer :: ncid,varid
    real(dp) :: offset,scale
    integer :: ncerr

    ! Get scaling and offset, if present, and apply ----

    ncerr=nf90_get_att(ncid, varid, 'add_offset', offset)
    if (ncerr /= NF90_NOERR) then
       offset = 0.d0
       ncerr = NF90_NOERR
    end if

    ncerr=nf90_get_att(ncid, varid, 'scale_factor', scale)
    if (ncerr/=NF90_NOERR) then
       scale = 1.d0
       ncerr = NF90_NOERR
    end if

  end subroutine read_ncdf_scaling

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine scale_time(params,time,units)
    
    use glimmer_log
    
    type(glex_climate),   intent(in)    :: params
    real(dp),dimension(:),intent(inout) :: time
    character(*),         intent(in)    :: units

    select case(trim(units))

    case('years')
       ! Do nothing

    case('hours')
       time = time/real(params%hours_in_year,dp)

    !WHL - Added 'months' and 'days'
    !TODO - Test these options
    case('months')
       time = time/12.d0      

    case('days')
       time = time/real(params%days_in_year,dp)

    case default
       call write_log('Time units '//trim(units)//' unrecognised',GM_FATAL)

    end select

  end subroutine scale_time

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine handle_err(status,line)

    use netcdf

    integer, intent (in) :: status
    integer, intent (in) :: line

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       print *, 'Line:',line
       stop "Stopped"
    end if
  end subroutine handle_err

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine example_climate(params,precip,temp,time)

    use glimmer_log

    type(glex_climate) :: params
    real(dp),dimension(:,:),intent(out)  :: precip,temp
    real(dp),intent(in) :: time ! Time (hours)

    integer :: ntemp,nprecip
    real(dp) :: tsp,tst
    real(dp) :: pos
    integer :: lower,upper

    real(dp) :: fyear

    ! Calculate fraction of year
    fyear = real(mod(time,real(params%hours_in_year,dp))) / real(params%hours_in_year,dp)

    ! Do temperature interpolation
    call bracket_point(fyear, params%st_time, lower, upper, pos)
    temp = linear_interp(params%surftemp_clim(:,:,lower), params%surftemp_clim(:,:,upper), pos)

    ! precip
    call bracket_point(fyear, params%pr_time, lower, upper, pos)
    precip = linear_interp(params%precip_clim(:,:,lower), params%precip_clim(:,:,upper), pos)

  end subroutine example_climate

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function linear_interp(a,b,pos)

    real(dp),dimension(:,:),intent(in) :: a,b
    real(dp),dimension(size(a,1),size(a,2)) :: linear_interp
    real(dp),               intent(in) :: pos

    linear_interp = a*(1.d0-pos) + b*pos

  end function linear_interp

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine bracket_point(n,a,lower,upper,frac)

    real(dp),             intent(in)  :: n      ! current fraction of year
    real(dp),dimension(:),intent(in)  :: a      ! array of fractional year values
    integer,              intent(out) :: lower
    integer,              intent(out) :: upper
    real(dp),             intent(out) :: frac

    real(dp),dimension(0:size(a)+1) :: aa
    integer :: na

    ! Array bounds
    na = size(a)
    aa(1:na) = a
    aa(0)    = -1 + a(na)
    aa(na+1) =  1 + aa(1)

    lower = 0
    upper = 1
    do
       if (n >= aa(lower) .and. n < aa(upper)) then
          exit
       end if
       lower = lower + 1
       upper = upper + 1
    end do
    frac = (n-aa(lower)) / (aa(upper)-aa(lower))

    call fixbounds(lower,1,na)
    call fixbounds(upper,1,na)

  end subroutine bracket_point

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine fixbounds(in,bottom,top)

    integer :: in,top,bottom

    do
       if (in<=top) exit
       in = in - (top-bottom+1)
    end do

    do
       if (in>=bottom) exit
       in = in + (top-bottom+1)
    end do

  end subroutine fixbounds

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_gcm_smb(temp,        precip,   &
                             orog,                  &
                             qsmb,        tsfc,     &
                             topo,                  &
                             glc_nec,     glc_topomax)

     use glimmer_physcon, only: scyr

     ! This is a crude parameterization for estimating qsmb and tsfc in each elevation class,
     !  given temp and precip at a given mean surface elevation.
     ! With these estimates we can run the standalone model (glint_example) as if we were
     !  getting qsmb and tsfc from a GCM.
     ! By tuning ablt_const, we can get an SMB that is not so different from what the PDD scheme computes.

     ! input fields on global grid
     real(dp), dimension(:,:), intent(in) :: temp     ! 2 m air temp (deg C)
     real(dp), dimension(:,:), intent(in) :: precip   ! precip rate  (mm/s = kg/m2/s)
     real(dp), dimension(:,:), intent(in) :: orog     ! global orography (m)

     ! output fields on global grid, with elevation class index
     real(dp), dimension(:,:,:), intent(inout) :: qsmb    ! ice sfc mass balance (kg/m2/s)
     real(dp), dimension(:,:,:), intent(inout) :: tsfc    ! ice sfc temp (deg C)
     real(dp), dimension(:,:,:), intent(inout) :: topo    ! ice sfc elevation (m)

     integer, intent(in) :: glc_nec                       ! number of elevation classes

     real(dp), dimension(0:glc_nec), intent(in) :: glc_topomax  ! upper elevation of each class (m)

     integer :: nx, ny, nec
     integer :: i, j, k

     real(dp), dimension(glc_nec) :: glc_topomid   ! midrange elevation of each class (m)

     real(dp) :: ablt                              ! ablation rate (kg/m2/s)

     real(dp), parameter :: lapse_rate = 0.006d0   ! temp lapse rate (deg/m)

     real(dp), parameter :: ablt_const = 5000.d0/scyr  ! ablation rate per degree above 0 C (converted from kg/m2/yr to kg/m2/s)
                                                       ! can be tuned to agree (more or less) with acab from PDD scheme

     ! get global grid size
     nx = size(temp,1)
     ny = size(temp,2)

     ! compute mid-range elevation of each class
     do k = 1, glc_nec-1
        glc_topomid(k) = 0.5d0 * (glc_topomax(k-1) + glc_topomax(k))
     enddo
     k = glc_nec
     glc_topomid(k) = 2.d0*glc_topomax(k-1) - glc_topomax(k-2)

     do k = 1, glc_nec

        do j = 1, ny
        do i = 1, nx

           ! set topo to midrange value for this elevation class
           topo(i,j,k) = glc_topomid(k)

           ! set tsfc assuming a fixed lapse rate
           tsfc(i,j,k) = temp(i,j) - lapse_rate * (topo(i,j,k) - orog(i,j))

           ! simple parameterization for ablation as function of temperature
           if (tsfc(i,j,k) > 0.d0) then
              ablt = ablt_const * tsfc(i,j,k)
           else
              ablt = 0.d0
           endif

           ! set smb as function of precip and ablation
           qsmb(i,j,k) = precip(i,j) - ablt

        enddo

     enddo   ! i
     enddo   ! j

  end subroutine compute_gcm_smb

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glint_example_clim
 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
