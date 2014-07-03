
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glex_ebm_clim.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module glint_ebm_ex_clim

  !*FD Module containing the glint example climate driver, for use with
  !*FD the energy balance mass-balance model. Unlike the main glint_example code,
  !*FD this one assumes that the data are all on the same grid.

  use glimmer_global
  use glint_global_grid

  implicit none

  type glex_ebm_climate
     ! Mass-balance coupling timing parameters --------------------------
     integer                 :: total_years=10    ! Length of run in years
     integer                 :: climate_tstep=6   ! Climate time-step in hours
     real(rk)                :: diurnal_cycle=0.0 ! Imposed Diurnal cycle (degC)
     ! Filenames --------------------------------------------------------
     character(fname_length) :: precip_file = 'monthly_precip_mean_1974-2003.nc' !*FD Name of precip file
     character(fname_length) :: stemp_file  = 'surf_temp_6h_1974-2003.nc'        !*FD Name of surface temp file
     character(fname_length) :: orog_file   = 'global_orog.nc'                   !*FD Name of orography file
     character(fname_length) :: zonwind_file  = ''      !*FD Name of zonal wind file
     character(fname_length) :: merwind_file  = ''      !*FD Name of meridional wind file
     character(fname_length) :: humid_file    = ''      !*FD Name of humidity file
     character(fname_length) :: lwdown_file   = ''      !*FD Name of lwdown file
     character(fname_length) :: swdown_file   = ''      !*FD Name of swdown file
     character(fname_length) :: airpress_file = ''      !*FD Name of air pressure file
     ! Variable names ---------------------------------------------------
     character(fname_length) :: precip_varname   = 'prate'
     character(fname_length) :: stemp_varname    = 'air_temperature'
     character(fname_length) :: orog_varname     = 'hgt'
     character(fname_length) :: zonwind_varname  = ''
     character(fname_length) :: merwind_varname  = ''
     character(fname_length) :: humid_varname    = ''
     character(fname_length) :: lwdown_varname   = ''
     character(fname_length) :: swdown_varname   = ''
     character(fname_length) :: airpress_varname = ''
     ! Arrays for holding climatology -----------------------------------
     real(rk),dimension(:,:),  pointer :: orog_clim     => null()  !*FD Orography
     real(rk),dimension(:,:,:),pointer :: precip_clim   => null()  !*FD Precip
     real(rk),dimension(:,:,:),pointer :: surftemp_clim => null()  !*FD Surface temperature
     real(rk),dimension(:,:,:),pointer :: zonwind_clim  => null()  !*FD Zonal wind
     real(rk),dimension(:,:,:),pointer :: merwind_clim  => null()  !*FD Meridional Wind
     real(rk),dimension(:,:,:),pointer :: humid_clim    => null()  !*FD Humidity
     real(rk),dimension(:,:,:),pointer :: lwdown_clim   => null()  !*FD Longwave down
     real(rk),dimension(:,:,:),pointer :: swdown_clim   => null()  !*FD Shortwave down
     real(rk),dimension(:,:,:),pointer :: airpress_clim => null()  !*FD Air pressure
     ! Grid variables ---------------------------------------------------
     type(global_grid),pointer :: grid   => null()
     type(global_grid),pointer :: grid_check => null()
     ! Other parameters -------------------------------------------------
     integer :: days_in_year=365
     integer :: hours_in_year=365*24
  end type glex_ebm_climate

  interface read_ncdf
     module procedure read_ncdf_2d,read_ncdf_3d
  end interface

contains

  subroutine glex_ebm_clim_init(params,filename)
    
    use glimmer_config
    use glint_global_interp
    use glimmer_log

    type(glex_ebm_climate) :: params   !*FD Climate parameters
    character(*)           :: filename !*FD config filename

    type(ConfigSection),pointer :: config !*FD structure holding sections of configuration file   
    integer :: ierr,i

    call ConfigRead(filename,config)
    call glex_ebm_clim_readconfig(params,config)
    call CheckSections(config)

    ! Read in climate data

    call read_ncdf(params%precip_file,   params%precip_varname,   params%precip_clim,    params%grid)
    call read_ncdf(params%stemp_file,    params%stemp_varname,    params%surftemp_clim,  params%grid_check)
    if (params%grid/=params%grid_check) call write_log('global data grid mis-match',GM_FATAL,&
         __FILE__,__LINE__)
    call read_ncdf(params%orog_file,     params%orog_varname,     params%orog_clim,      params%grid_check)
    if (params%grid/=params%grid_check) call write_log('global data grid mis-match',GM_FATAL,&
         __FILE__,__LINE__)
    call read_ncdf(params%zonwind_file,  params%zonwind_varname,  params%zonwind_clim,   params%grid_check)
    if (params%grid/=params%grid_check) call write_log('global data grid mis-match',GM_FATAL,&
         __FILE__,__LINE__)
    call read_ncdf(params%merwind_file,  params%merwind_varname,  params%merwind_clim,   params%grid_check)
    if (params%grid/=params%grid_check) call write_log('global data grid mis-match',GM_FATAL,&
         __FILE__,__LINE__)
    call read_ncdf(params%humid_file,    params%humid_varname,    params%humid_clim,     params%grid_check)
    if (params%grid/=params%grid_check) call write_log('global data grid mis-match',GM_FATAL,&
         __FILE__,__LINE__)
    call read_ncdf(params%lwdown_file,   params%lwdown_varname,   params%lwdown_clim,    params%grid_check)
    if (params%grid/=params%grid_check) call write_log('global data grid mis-match',GM_FATAL,&
         __FILE__,__LINE__)
    call read_ncdf(params%swdown_file,   params%swdown_varname,   params%swdown_clim,    params%grid_check)
    if (params%grid/=params%grid_check) call write_log('global data grid mis-match',GM_FATAL,&
         __FILE__,__LINE__)
    call read_ncdf(params%airpress_file, params%airpress_varname, params%airpress_clim,  params%grid_check)
    if (params%grid/=params%grid_check) call write_log('global data grid mis-match',GM_FATAL,&
         __FILE__,__LINE__)

    ! Fix up a few things

    params%surftemp_clim=params%surftemp_clim-273.15       ! Convert temps to degreesC
    params%precip_clim=params%precip_clim*1000.0/21600.0   ! Convert precip to mm/s
    params%lwdown_clim=params%lwdown_clim/21600.0          ! Convert fluxes to W/m^2
    params%swdown_clim=params%swdown_clim/21600.0          ! Convert fluxes to W/m^2
    params%orog_clim=params%orog_clim/9.81                 ! Convert geopotential to m

  end subroutine glex_ebm_clim_init

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glex_ebm_clim_readconfig(params,config)

    use glimmer_config

    type(glex_ebm_climate)       :: params !*FD Climate parameters
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file   
    type(ConfigSection), pointer :: section

    call GetSection(config,section,'GLEX precip')
    if (associated(section)) then
       call GetValue(section,'filename',params%precip_file)
       call GetValue(section,'variable',params%precip_varname)
    end if
    call GetSection(config,section,'GLEX temps')
    if (associated(section)) then
       call GetValue(section,'filename',params%stemp_file)
       call GetValue(section,'variable',params%stemp_varname)
    end if
    call GetSection(config,section,'GLEX orog')
    if (associated(section)) then
       call GetValue(section,'filename',params%orog_file)
       call GetValue(section,'variable',params%orog_varname)
    end if
    call GetSection(config,section,'GLEX humid')
    if (associated(section)) then
       call GetValue(section,'filename',params%humid_file)
       call GetValue(section,'variable',params%humid_varname)
    end if
    call GetSection(config,section,'GLEX zonwind')
    if (associated(section)) then
       call GetValue(section,'filename',params%zonwind_file)
       call GetValue(section,'variable',params%zonwind_varname)
    end if
    call GetSection(config,section,'GLEX merwind')
    if (associated(section)) then
       call GetValue(section,'filename',params%merwind_file)
       call GetValue(section,'variable',params%merwind_varname)
    end if
    call GetSection(config,section,'GLEX lwdown')
    if (associated(section)) then
       call GetValue(section,'filename',params%lwdown_file)
       call GetValue(section,'variable',params%lwdown_varname)
    end if
    call GetSection(config,section,'GLEX swdown')
    if (associated(section)) then
       call GetValue(section,'filename',params%swdown_file)
       call GetValue(section,'variable',params%swdown_varname)
    end if
    call GetSection(config,section,'GLEX airpress')
    if (associated(section)) then
       call GetValue(section,'filename',params%airpress_file)
       call GetValue(section,'variable',params%airpress_varname)
    end if
    call GetSection(config,section,'GLEX climate')
    if (associated(section)) then
       call GetValue(section,'days_in_year',params%days_in_year)
       call GetValue(section,'total_years',params%total_years)
       call GetValue(section,'climate_tstep',params%climate_tstep)
       call GetValue(section,'diurnal_cycle',params%diurnal_cycle)
       params%hours_in_year=params%days_in_year*24
    end if

  end subroutine glex_ebm_clim_readconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_2d(filename,varname,array,grid)

    use netcdf

    character(*)                    :: filename,varname
    real(rk),dimension(:,:),pointer :: array
    type(global_grid),      pointer :: grid

    real(rk),dimension(:),allocatable :: dim1,dim2
    real(rk),dimension(:),pointer :: lonbound => NULL()
    real(rk),dimension(:),pointer :: latbound => NULL()

    integer  :: ncerr     ! NetCDF error 
    integer  :: ncid      ! NetCDF file id
    integer  :: varid     ! NetCDF variable id
    integer  :: ndims     ! Number of dimensions
    integer  :: i,args
    real(rk) :: offset=0.0,scale=1.0
    integer,      dimension(2) :: dimids,dimlens
    character(20),dimension(2) :: dimnames
    logical :: lonb_present,latb_present

    if (associated(array)) deallocate(array)
    if (associated(grid))  deallocate(grid)

    call read_ncdf_common1(filename,ncid,varid,ndims,varname)
    
    ! If not a 3d variable, flag and error and exit ----

    if (ndims/=2) then
       print*,'NetCDF: Requested variable only has ',ndims,' dimensions'
       stop
    end if

    call read_ncdf_common2(ncid,varid,ndims,dimids,dimlens,dimnames)

    ! Allocate output and dimension arrays -------------

    allocate(array(dimlens(1),dimlens(2)))
    allocate(dim1(dimlens(1)))
    allocate(dim2(dimlens(2)))
    
    ! Retrieve variable contents -----------------------

    ncerr=nf90_get_var(ncid, varid, array)
    call handle_err(ncerr,__LINE__)

    call read_ncdf_common3(ncid,varid,offset,scale)

    array=offset+(array*scale)

    ! Get dimension variables --------------------------

    call read_ncdf_getdim(ncid,varid,dimnames(1),dim1)
    call read_ncdf_getdim(ncid,varid,dimnames(2),dim2)

    ! Get boundary arrays, if present ------------------

    call read_ncdf_getbound(ncid,varid,'bounds_lon',lonb_present,lonbound)
    call read_ncdf_getbound(ncid,varid,'bounds_lat',latb_present,latbound)

    ! Construct global grid ----------------------------

    call read_ncdf_common4(grid,lonb_present,latb_present,dim1,dim2,lonbound,latbound)

    ! Tidy up ------------------------------------------

    deallocate(dim1,dim2)
    if (associated(latbound)) deallocate(latbound)
    if (associated(lonbound)) deallocate(lonbound)

  end subroutine read_ncdf_2d

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_3d(filename,varname,array,grid)

    use netcdf

    character(*)                      :: filename,varname
    real(rk),dimension(:,:,:),pointer :: array
    type(global_grid),        pointer :: grid

    real(rk),dimension(:),allocatable :: dim1,dim2,dim3
    real(rk),dimension(:),pointer :: lonbound => NULL()
    real(rk),dimension(:),pointer :: latbound => NULL()

    integer  :: ncerr     ! NetCDF error 
    integer  :: ncid      ! NetCDF file id
    integer  :: varid     ! NetCDF variable id
    integer  :: ndims     ! Number of dimensions
    integer  :: i,args
    real(rk) :: offset=0.0,scale=1.0
    integer,      dimension(3) :: dimids,dimlens
    character(20),dimension(3) :: dimnames
    real(dp),     dimension(:,:),allocatable :: lnb,ltb
    logical :: lonb_present,latb_present

    if (associated(array)) deallocate(array)
    if (associated(grid))  deallocate(grid)

    call read_ncdf_common1(filename,ncid,varid,ndims,varname)
    
    ! If not a 3d variable, flag and error and exit ----

    if (ndims/=3) then
       print*,'NetCDF: Requested variable only has ',ndims,' dimensions'
       stop
    end if

    call read_ncdf_common2(ncid,varid,ndims,dimids,dimlens,dimnames)

    ! Allocate output and dimension arrays -------------

    allocate(array(dimlens(1),dimlens(2),dimlens(3)))
    allocate(dim1(dimlens(1)))
    allocate(dim2(dimlens(2)))
    allocate(dim3(dimlens(3)))
    
    ! Retrieve variable contents -----------------------

    ncerr=nf90_get_var(ncid, varid, array)
    call handle_err(ncerr,__LINE__)

    call read_ncdf_common3(ncid,varid,offset,scale)

    array=offset+(array*scale)

    ! Get dimension variables --------------------------

    call read_ncdf_getdim(ncid,varid,dimnames(1),dim1)
    call read_ncdf_getdim(ncid,varid,dimnames(2),dim2)
    call read_ncdf_getdim(ncid,varid,dimnames(3),dim3)

    ! Get boundary arrays, if present ------------------

    call read_ncdf_getbound(ncid,varid,'bounds_lon',lonb_present,lonbound)
    call read_ncdf_getbound(ncid,varid,'bounds_lat',latb_present,latbound)

    ! Construct global grid ----------------------------

    call read_ncdf_common4(grid,lonb_present,latb_present,dim1,dim2,lonbound,latbound)

    ! Tidy up ------------------------------------------

    deallocate(dim1,dim2,dim3)
    if (associated(latbound)) deallocate(latbound)
    if (associated(lonbound)) deallocate(lonbound)

  end subroutine read_ncdf_3d

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_common1(filename,ncid,varid,ndims,varname)

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

  end subroutine read_ncdf_common1

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_common2(ncid,varid,ndims,dimids,dimlens,dimnames)

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

  end subroutine read_ncdf_common2

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_common3(ncid,varid,offset,scale)

    use netcdf

    integer :: ncid,varid
    real(rk) :: offset,scale
    integer :: ncerr

    ! Get scaling and offset, if present, and apply ----

    ncerr=nf90_get_att(ncid, varid, 'add_offset', offset)
    if (ncerr/=NF90_NOERR) then
       offset=0.0
       ncerr=NF90_NOERR
    end if

    ncerr=nf90_get_att(ncid, varid, 'scale_factor', scale)
    if (ncerr/=NF90_NOERR) then
       scale=1.0
       ncerr=NF90_NOERR
    end if

  end subroutine read_ncdf_common3

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_getdim(ncid,varid,dimname,dim)

    use netcdf

    integer :: ncid,varid
    character(*) :: dimname
    real(rk),dimension(:) :: dim
    integer :: ncerr

    ! Get dimension variables 

    ncerr=nf90_inq_varid(ncid,dimname,varid)
    call handle_err(ncerr,__LINE__)
    ncerr=nf90_get_var(ncid, varid, dim)
    call handle_err(ncerr,__LINE__)

  end subroutine read_ncdf_getdim

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_getbound(ncid,varid,boundname,present,bound)

    use netcdf

    integer :: ncid,varid
    character(*) :: boundname
    logical :: present
    real(rk),dimension(:),pointer :: bound

    integer :: ncerr
    integer,dimension(2) :: dimids,dimlens
    character(20),dimension(2) :: dimnames
    real(rk),dimension(:,:),allocatable :: b
    integer :: i

    ncerr=nf90_inq_varid(ncid,boundname,varid)
    if (ncerr/=NF90_NOERR) then
       present=.false.
       ncerr=NF90_NOERR
    else
       ncerr=nf90_inquire_variable(ncid, varid, dimids=dimids)
       call handle_err(ncerr,__LINE__)
       do i=1,2
          ncerr=nf90_inquire_dimension(ncid, dimids(i), name=dimnames(i),len=dimlens(i))
          call handle_err(ncerr,__LINE__)
       end do
       if (associated(bound)) deallocate(bound)
       allocate(b(dimlens(1),dimlens(2)),bound(dimlens(2)+1))
       ncerr=nf90_get_var(ncid, varid,b)
       call handle_err(ncerr,__LINE__)
       do i=1,dimlens(2)
          bound(i)=b(1,i)
          bound(i+1)=b(2,i)
       end do
       deallocate(b)
       present=.true.
    end if

  end subroutine read_ncdf_getbound

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_ncdf_common4(grid,lonb_present,latb_present,dim1,dim2,lonbound,latbound)

    type(global_grid),        pointer :: grid
    logical :: lonb_present,latb_present
    real(rk),dimension(:) :: dim1,dim2
    real(rk), dimension(:), pointer :: lonbound,latbound
    integer :: args

    ! Construct grid type

    allocate(grid)
    args=0
    if (lonb_present) args=args+1
    if (latb_present) args=args+2

    select case(args)
    case(0)
       call new_global_grid(grid,dim1,dim2,correct=.false.)
    case(1)
       call new_global_grid(grid,dim1,dim2,lonb=lonbound)
    case(2)
       call new_global_grid(grid,dim1,dim2,latb=latbound,correct=.false.)
    case(3)
       call new_global_grid(grid,dim1,dim2,lonb=lonbound,latb=latbound)
    end select

  end subroutine read_ncdf_common4

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

  subroutine ebm_ex_climate(params,precip,temp,zonwind,merwind,humid,lwdown,swdown,airpress,time)

    use glimmer_log

    type(glex_ebm_climate) :: params
    real(rk),dimension(:,:),intent(out)  :: precip,temp,zonwind,merwind,humid,lwdown,swdown,airpress
    real(rk),intent(in) :: time

    integer :: ntemp,nprecip,nzonwind,nmerwind,nhumid,nlwdown,nswdown,nairpress
    real(rk) :: tsp,tst,tsz,tsm,tsh,tsl,tss,tsa
    real(rk) :: pos
    integer :: lower,upper
    character(150) :: msg

    ntemp     = size(params%surftemp_clim,3)
    nprecip   = size(params%precip_clim,3)
    nzonwind  = size(params%zonwind_clim,3)
    nmerwind  = size(params%merwind_clim,3)
    nhumid    = size(params%humid_clim,3)
    nlwdown   = size(params%lwdown_clim,3)
    nswdown   = size(params%swdown_clim,3)
    nairpress = size(params%airpress_clim,3)

    tst=params%hours_in_year/ntemp
    tsp=params%hours_in_year/nprecip
    tsz=params%hours_in_year/nzonwind
    tsm=params%hours_in_year/nmerwind
    tsh=params%hours_in_year/nhumid
    tsl=params%hours_in_year/nlwdown
    tss=params%hours_in_year/nswdown
    tsa=params%hours_in_year/nairpress
    
    temp    =interp_field(time,tst,ntemp,    params%surftemp_clim)
    precip  =interp_field(time,tsp,nprecip,  params%precip_clim)
    zonwind =interp_field(time,tsz,nzonwind, params%zonwind_clim)
    merwind =interp_field(time,tsm,nmerwind, params%merwind_clim)
    humid   =interp_field(time,tsh,nhumid,   params%humid_clim)
    lwdown  =interp_field(time,tsl,nlwdown,  params%lwdown_clim)
    swdown  =interp_field(time,tss,nswdown,  params%swdown_clim)
    airpress=interp_field(time,tsa,nairpress,params%airpress_clim)

    ! Add diurnal cycle to temperature. We assume that
    ! the lowest temperature is at midnight, and the highest at midday.
    ! Obviously, this probably isn't true...

    if (mod(int(time),24)==0)  temp=temp-params%diurnal_cycle
    if (mod(int(time),24)==12) temp=temp+params%diurnal_cycle

  end subroutine ebm_ex_climate

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function interp_field(time,ts,n,clim)
 
    real(rk),intent(in) :: time
    real(rk),intent(in) :: ts
    integer, intent(in) :: n
    real(rk),dimension(:,:,:),intent(in) :: clim
    real(rk),dimension(size(clim,1),size(clim,2)) :: interp_field

    integer :: lower,upper
    real(rk) :: pos

    lower=int(time/ts)
    upper=lower+1
    pos=mod(time,ts)/ts
    call fixbounds(lower,1,n)
    call fixbounds(upper,1,n)
    interp_field=linear_interp(clim(:,:,lower),clim(:,:,upper),pos)

  end function interp_field

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function linear_interp(a,b,pos)

    real(rk),dimension(:,:),intent(in) :: a,b
    real(rk),dimension(size(a,1),size(a,2)) :: linear_interp
    real(rk),               intent(in) :: pos

    linear_interp=a*(1.0-pos)+b*pos

  end function linear_interp

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine fixbounds(in,bottom,top)

    integer :: in,top,bottom

    do
       if (in<=top) exit
       in=in-(top-bottom+1)
    end do

    do
       if (in>=bottom) exit
       in=in+(top-bottom+1)
    end do

  end subroutine fixbounds

end module glint_ebm_ex_clim
