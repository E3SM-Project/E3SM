!
! Interpolates data in both horizontal and vertical directions
! and writes it to a netcdf file using PIO (which uses pNetCDF)
!
! For questions, please contact: david.hall@colorado.edu
!_______________________________________________________________________

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module netcdf_interp_mod

  !$ use omp_lib                                                                ! load openmp library code
  use control_mod,      only: test_case!, vinterp_nlev
  use dimensions_mod,   only: nelemd, np, ne, nlev
  use element_mod,      only: element_t
  use hybrid_mod,       only: hybrid_t,hybrid_create
  use hybvcoord_mod,    only: hvcoord_t
  use kinds,            only: real_kind
  use element_ops,      only: get_field, get_vector_field
  use parallel_mod,     only: parallel_t, haltmp, syncmp
  use shr_const_mod,    only: shr_const_cday
  use time_mod,         only: timelevel_t,time_at
  use vertical_se,      only: v_interpolate, L_interp

  use common_io_mod,    only: &
    max_output_streams, output_frequency, &
    output_start_time, output_end_time, varname_len, piofs, num_io_procs

  use interpolate_mod,  only: &
    interpolate_t,interpdata_t, setup_latlon_interp, interpolate_scalar, &
    get_interp_parameter, get_interp_lat, get_interp_lon, interpolate_vector

  use pio_io_mod,       only: &
    nf_output_init_begin, nf_handle, &
    nf_output_register_variables, nf_output_register_dims, &
    nf_variable_attributes, nf_output_init_complete, pio_double, &
    get_current_varnames, nf_selectedvar, nf_get_frame, nf_put_var => nf_put_var_pio, &
    nfsizekind, get_varindex, pio_put_var, nf_close_all, pio_setDebugLevel, &
    pio_put_att, pio_put_var, pio_global, nf_init_decomp, nf_advance_frame,&
    pio_syncfile

  implicit none

  integer, parameter :: n_dims = 4                                      ! number of dimensions
  character*(*),parameter :: dim_names(n_dims)=(/'lon ','lat ','lev ','time'/)

  !_____________________________________________________________________
  ! scalar field descriptor

  type nc_sfield_t
    integer         :: var_type           ! pio type of variable
    integer         :: dim_ids(n_dims)    ! dimension ids
    integer         :: tensor_dim         ! tensor dim:0=scalar,1=vector,2=tensor
    logical         :: is_required        ! variable required flag
    character*(5)   :: short_name         ! short name of variable
    character*(256) :: units              ! physical units
    character*(256) :: long_name          ! full name of variable
  endtype

  !_____________________________________________________________________
  ! vector field descriptor

  type nc_vfield_t
    character*(4)   :: vec_name           ! vector name
    integer         :: n_comp             ! number of components
    character*(4)   :: short_name(2)      ! name of each component
  endtype

  !_____________________________________________________________________
  !  possible output variables

  integer, parameter :: n_vars= 38 ! number of possible output fields

  type(nc_sfield_t), parameter :: nc_vars(n_vars) =&
  (/&
    nc_sfield_t(pio_double, (/1,0,0,0/),0, .true. , 'lon  ', 'degrees_east',    'column longitude'           ),&
    nc_sfield_t(pio_double, (/2,0,0,0/),0, .true. , 'lat  ', 'degrees_north',   'column latitude'            ),&
    nc_sfield_t(pio_double, (/3,0,0,0/),0, .true. , 'lev  ', '',                'vertical coordinate'        ),&
    nc_sfield_t(pio_double, (/4,0,0,0/),0, .true. , 'time ', 'days',            'model elapsed time'         ),&
    nc_sfield_t(pio_double, (/3,0,0,0/),0, .true. , 'hyam ', '',                'hybrid A at midpoints'      ),&
    nc_sfield_t(pio_double, (/3,0,0,0/),0, .true. , 'hybm ', '',                'hybrid B at midpoints'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'ps   ', 'Pa',              'hydrostatic surf pressure'  ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'z    ', 'meters',          'altitiude'                  ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'geo  ', 'meters^2/sec^2',  'geopotential height'        ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'phi  ', 'meters^2/sec^2',  'geopotential height'        ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'phis ', 'meters^2/sec^2',  'surface geopotential'       ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'p    ', 'Pa',              'hydrostatic atmo pressure'  ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'pt   ', 'Pa',              'total atmospheric pressure' ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'T    ', 'degrees Kelvin',  'temperature'                ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'theta', 'degrees Kelvin',  'potential temperature'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),1, .false., 'u    ', 'meters/second',   'longitudinal wind component'),&
    nc_sfield_t(pio_double, (/1,2,3,4/),1, .false., 'v    ', 'meters/second',   'latitudinal wind component' ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'w    ', 'meters/second',   'vertical wind component'    ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'pp   ', 'Pa',              'pressure deviation'         ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'Q    ', 'kg/kg',           'tracer 1 mixing ratio'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'Q2   ', 'kg/kg',           'tracer 2 mixing ratio'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'Q3   ', 'kg/kg',           'tracer 3 mixing ratio'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'Q4   ', 'kg/kg',           'tracer 4 mixing ratio'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'omega', '',                ''                           ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'D3   ', '',                '3d divergence of velocity'  ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'speed', 'meters/second',   'wind speed'                 ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'KE   ', 'Joules',          'kinetic energy'             ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'alpha', '',                'specific volume'            ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'rho  ', '',                'density'                    ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'dpdn ', '',                'pseudo-density'             ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'm    ', '',                'pseudo-density'             ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'cflux', '',                'column int horiz mass flux' ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'ndot ', '',                'vertical eta velocity'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'nflux', '',                'vertical eta flux'          ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'O1   ', '',                'debugging output 1'         ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'O2   ', '',                'debugging output 2'         ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'O3   ', '',                'debugging output 3'         ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'O4   ', '',                'debugging output 4'         ) &
    !nc_sfield_t(type,       dim_ids,    n_spatial_dims, required?, short_name, units, long_name)
  /)

  type(nc_vfield_t),parameter::nc_vfields(1)=&                          ! list vector fields
  (/&
    !nc_vfield_t(vec_name, n_comps, (/component_names/))
    nc_vfield_t ('v'     , 2      , (/'u   ','v   '/)) &
  /)

  integer :: ni                                                         ! number of vertical levels
  integer :: dim_sizes(n_dims)                                          ! size of each dimension

  type(interpdata_t), allocatable :: interpdata(:)                      ! interpolated data-locations
  type(nf_handle),    save        :: files(max_output_streams)          ! netcdf file-handles

  integer ::  &
    n_interp, &                                                         ! num horiztonal interp points on this proc
    n_lat,    &                                                         ! num points in lat direction
    n_lon                                                               ! num points in lon direction

  integer(kind=nfsizekind)::  &
    start2d(3), &                                                       ! netcdf offset for 2d data
    start3d(4), &                                                       ! netcdf offset for 3d data
    count2d(3), &                                                       ! netcdf data size for 2d data
    count3d(4)                                                          ! netcdf data size for 3d data

  character(len=varname_len), pointer :: output_varnames(:)             ! list of variables to write


  ! convert nc_vars to 1d-arrays

  integer,        parameter :: n_vfields          = size(nc_vfields)            ! count num vector fields
  ! using *(*) gives segfault with gnu fortran:
  character*(5),  parameter :: var_names  (n_vars)= nc_vars(:)%short_name       ! get array of variable names
  integer,        parameter :: var_type   (n_vars)= nc_vars(:)%var_type         ! get array of variable types
  logical,        parameter :: var_req    (n_vars)= nc_vars(:)%is_required      ! get array of requirement flags
  integer                   :: var_dim_ids(n_dims,n_vars)                       ! get array of dimension ids

  contains

  !_____________________________________________________________________________
  function index_of(var_name) result(index)
    character*(*), intent(in) :: var_name

    integer :: i, index

    index=-1
    do i=1,n_vars
      if(var_name == var_names(i)) index=i
    enddo
  end function

  !_____________________________________________________________________________
  integer function size_of_dim(index) result(r)
    integer, intent(in) :: index
    r = 1
    if(index>0) r = dim_sizes(index)
    if (r==0) r=1
  end function

  !_____________________________________________________________________________
  subroutine netcdf_interp_init(elem, hybrid, hvcoord)

    ! Define netcdf variables, dimensions, and attributes

    type(element_t),        intent(in):: elem(:)
    type(hybrid_t), target, intent(in):: hybrid
    type(hvcoord_t),        intent(in):: hvcoord

    type(parallel_t), pointer:: par

    integer, parameter :: run_type = 0                                  ! set run_type = 'new run'

    integer i                                                           ! loop variable
    integer sz(3)                                                       ! buffers for dimension sizes
    integer iorank
    integer :: ios                                                      ! io stream number
    integer result

    call pio_setDebugLevel(0)                                           ! disable pio debugging
    par=> hybrid%par

    if (par%masterproc) then
      print *,'initializing netcdf output =',var_names
     ! print *,"vinterp_type = ",vinterp_type
     ! print *,"vinterp_nlev = ",vinterp_nlev
    endif

    ni =nlev
    !if(vinterp_type==1) ni = vinterp_nlev

    ! convert nc_vars dimensions into 2d-array

    do i=1,n_vars
      var_dim_ids(:,i)= nc_vars(i)%dim_ids(:)
    enddo

    ! allocate memory for interpolated points

    allocate(interpdata(nelemd))                                        ! allocate one interpdata struct per element

    ! initialize lat-lon interpolation

    call setup_latlon_interp(elem, interpdata, par)

    n_lat    = get_interp_parameter('nlat')                             ! get num lat points in interpolation
    n_lon    = get_interp_parameter('nlon')                             ! get num lon points in interpolation
    n_interp = sum(interpdata(1:nelemd)%n_interp)                       ! count horizontal interp points on this proc

    print *,par%rank,': n_lat=',n_lat,' n_lon=',n_lon,' n_interp=',n_interp

    ! set dimension sizes

    dim_sizes=(/n_lon, n_lat, ni, 0/)                                   ! set size of each dimension

    ! initialize netcdf file: add dimensions and variables

    if (par%masterproc) print *,'known output variables =',var_names

    call nf_output_init_begin(files, par%masterproc, par%nprocs, par%rank, par%comm, test_case, run_type)
    call nf_output_register_dims(files, n_dims, dim_names, dim_sizes)         
    call nf_output_register_variables(files, n_vars, var_names, var_dim_ids, var_type, var_req)

    ! set variable attributes: long name and dimensional units

    do i=1,n_vars
      call nf_variable_attributes(files, nc_vars(i)%short_name, nc_vars(i)%long_name, nc_vars(i)%units)
    enddo

    ! save horizontal NP and NE as global attributes

    do ios = 1, max_output_streams
      if((output_frequency(ios) .gt. 0) ) then
        result = pio_put_att(files(ios)%fileid, pio_global, 'np', np)
        result = pio_put_att(files(ios)%fileid, pio_global, 'ne', ne)
      endif
    enddo

    ! exit netcdf definition mode

    call nf_output_init_complete(files)

    ! write lat and lon coordinates to file

    call write_coordinates_to_file(hybrid, hvcoord)

    ! assign sets of lat-lon points to the different output streams

    call decompose_output()

  end subroutine

  !_____________________________________________________________________
  subroutine netcdf_interp_finish()

    ! Flush data and close netcdf output streams

    call nf_close_all(files)
  end subroutine

  !_____________________________________________________________________
  subroutine write_coordinates_to_file(hybrid,hvcoord)

    ! Write latitude and longitude interpolation coordinates to file

    type(hybrid_t), intent(in):: hybrid
    type(hvcoord_t),intent(in):: hvcoord

    integer :: ios                                                      ! io-stream index
    integer :: iorank                                                   ! rank of iostream
    integer :: varid                                                    ! variable id
    integer :: vindex                                                   ! variable index
    integer :: ierr                                                     ! error number

    real(real_kind), allocatable :: lat(:),lon(:),lev(:),hyam(:),hybm(:)              ! array of lat and lon values

    iorank=piofs%io_rank                                                ! get iostream rank

    allocate(lon(n_lon)); lon  = get_interp_lon()                       ! allocate and fill longitude array
    allocate(lat(n_lat)); lat  = get_interp_lat()                       ! allocate and fill latitiude array
    allocate(lev(ni));    lev  = L_interp                               ! allocate and fill level array
    allocate(hyam(ni));   hyam = v_interpolate(hvcoord%hyam,ni)         ! allocate and fill hybrid A array
    allocate(hybm(ni));   hybm = v_interpolate(hvcoord%hybm,ni)         ! allocate and fill hybrid B array

    do ios=1,max_output_streams
      if((output_frequency(ios) .gt. 0) ) then
        if(iorank==0) print *,"writing coordinates to ios=",ios

        ! write lon values to file

        vindex  = get_varindex('lon', files(ios)%varlist)               ! look up index of longitude var
        varid   = files(ios)%varlist(vindex)%vardesc%varid              ! get longitude variable id
        ierr    = pio_put_var(files(ios)%fileID, varid, lon)            ! write longitude array to file

        ! write lat values to file
      
        vindex  = get_varindex('lat',files(ios)%varlist)                ! look up index of latitide var
        varid   = files(ios)%varlist(vindex)%vardesc%varid              ! get latitide variable id
        ierr    = pio_put_var(files(ios)%fileID, varid, lat)            ! write latitude array to file

        ! write lev values to file

        vindex  = get_varindex('lev',files(ios)%varlist)                ! look up index of lev var
        varid   = files(ios)%varlist(vindex)%vardesc%varid              ! get lev variable id
        ierr    = pio_put_var(files(ios)%fileID, varid, lev)            ! write level array to file

        ! write hyam values to file

        vindex  = get_varindex('hyam',files(ios)%varlist)               ! look up index of lev var
        varid   = files(ios)%varlist(vindex)%vardesc%varid              ! get lev variable id
        ierr    = pio_put_var(files(ios)%fileID, varid, hyam)           ! write level array to file

        ! write hyam values to file

        vindex  = get_varindex('hybm',files(ios)%varlist)               ! look up index of lev var
        varid   = files(ios)%varlist(vindex)%vardesc%varid              ! get lev variable id
        ierr    = pio_put_var(files(ios)%fileID, varid, hybm)           ! write level array to file

        call pio_syncfile(files(ios)%fileid)

      endif
    end do

    call syncmp(hybrid%par)                                             ! synchronize openmp threads

  end subroutine

  !_____________________________________________________________________________
  subroutine decompose_output()

    ! convert lon,lat,lev positions to 1d indices and distrubute amongsts procs

    integer, pointer ::&
      ldof_2d(:)  , &
      ldof_3d(:)  , &
      io_dof_2d(:), &
      io_dof_3d(:)
  
    integer ie,i,i1d,k, io_dof , iorank                                 ! loop indices

    !map 2d(lon,lat) indices to 1d offsets

    allocate(ldof_2d(n_interp))                                         ! allocate space for 1d index array
    i1d=0                                                               ! initialize 1d index
    do ie=1,nelemd                                                      ! loop over each element
       do i=1,interpdata(ie)%n_interp                                   ! loop over each interp-point in element
          i1d = i1d + 1                                                 ! increment 1d index
          ldof_2d(i1d)= interpdata(ie)%ilon(i)+ (interpdata(ie)%ilat(i)-1)*n_lon! covert 2d lat-lon index to 1d index
       end do
    end do

    ! map 3d(lon,lat,level) indices to 1d-offsets

    allocate(ldof_3d(n_interp * ni))                                    ! allocate space for 1d index array
    i1d=0                                                               ! initialize 1d index
    do k=1,ni                                                           ! loop over each vertical level
       do ie=1,nelemd                                                   ! loop over each element
          do i=1,interpdata(ie)%n_interp                                ! loop over each interp point in element
              i1d=i1d+1
              ldof_3d(i1d)=interpdata(ie)%ilon(i)+(interpdata(ie)%ilat(i)-1)*n_lon+(k-1)*n_lat*n_lon
          end do
       end do
    end do

    ! divide output data into slices

    iorank = piofs%io_rank
    call getiodof(2,(/n_lon,n_lat/),   iorank, io_dof_2d, start2d(1:2), count2d(1:2))
    call getiodof(3,(/n_lon,n_lat,ni/),iorank, io_dof_3d, start3d(1:3), count3d(1:3))

    ! initialize output decomposition

    call nf_init_decomp(files, (/1,2/),   ldof_2d, io_dof_2d, start2d(1:2), count2d(1:2))
    call nf_init_decomp(files, (/1,2,3/), ldof_3d, io_dof_3d, start3d(1:3), count3d(1:3))

  end subroutine

  !_____________________________________________________________________________
  subroutine write_scalar_field_1d(elem, var, ios, n0, hybrid, hvcoord)

    ! Get 3d scalar field data by name, and write to file

    type(element_t),  intent(in) :: elem(:)
    type(nc_sfield_t),intent(in) :: var                                 ! variable descriptor
    integer,          intent(in) :: ios                                 ! io stream number
    integer,          intent(in) :: n0                                  ! time-level index
    type(hybrid_t),   intent(in) :: hybrid
    type(hvcoord_t),  intent(in) :: hvcoord

    real(real_kind) :: element_data(np,np,ni)                           ! native grid data
    real(real_kind), allocatable :: result(:,:)                         ! interpolated data
    integer :: st, en                                                   ! start and end index for interpolation
    integer :: ie, var_index, vertical_dim                              ! loop indicies
    integer :: n_vert                                                   ! number of vertical levels

    output_varnames=>get_current_varnames(ios)                          ! get list of output variables

    ! check to see if field is selected and if it is a scalar field

    if(.not. nf_selectedvar(var%short_name, output_varnames)) return    ! skip if not selected
    if(var%tensor_dim .ne. 0) return                                    ! skip if part of a vector or tensor
    if (hybrid%par%masterproc) print *,'writing interpolated scalar field "',var%short_name,'"'

    var_index     = index_of(var%short_name)
    vertical_dim  = var_dim_ids(3,var_index)
    n_vert        = size_of_dim(vertical_dim)
    !if (hybrid%par%masterproc) print *,"index_of(",var%short_name,")=",var_index," vertical_dim=",vertical_dim," n_vert=",n_vert

    ! interpolate scalar field data to result field

    allocate(result(n_interp,n_vert))                                     ! allocate buffer for interpolated data
    result = 0.0d0

    st = 1                                                              ! set start index
    do ie=1,nelemd                                                      ! loop over elem
      en = st + interpdata(ie)%n_interp-1                               ! set end index
      call get_field(elem(ie),var%short_name,element_data,hvcoord,n0,n0)    ! get scalar field by name
      call interpolate_scalar(interpdata(ie),element_data, np, n_vert, result(st:en,:))
      st=st+interpdata(ie)%n_interp                                     ! shift start index
    enddo

    if(n_vert>1) then
      ! write 3d interpolated scalar field to file

      start3d(4)   = nf_get_frame(files(ios))
      count3d(4)   = 1
      count3d(1:3) = dim_sizes( (/1,2,3/) )                             ! set size of output
      call nf_put_var(files(ios), result, start3d, count3d, name=var%short_name)
    endif

    if(n_vert==1) then
    ! write 2d interpolated scalar field to file

      start2d(3)   = nf_get_frame(files(ios))
      count2d(3)   = 1
      count2d(1:2) = dim_sizes( (/1,2/) )                               ! set size of output
      call nf_put_var(files(ios), result, start2d, count2d, name=var%short_name)
    endif

    deallocate(result)

  end subroutine

  !_____________________________________________________________________________
  ! write_vector_field: write vector field components to file, if selected

  subroutine write_vector_field(elem, var, ios, n0, hybrid, hvcoord)

    type (element_t), intent(in) :: elem(:)
    type(nc_vfield_t),intent(in) :: var                                 ! variable descriptor
    integer,          intent(in) :: ios                                 ! io stream number
    integer,          intent(in) :: n0                                  ! time-level index
    type(hybrid_t),   intent(in) :: hybrid
    type(hvcoord_t),  intent(in) :: hvcoord

    real(real_kind) :: vdata(np,np,2,ni)                                ! nodal vector data

    real(real_kind), allocatable :: result(:,:,:)                       ! interpolated lat lon data

    character*(10) :: short_name                                        ! buffer for current var short name

    integer :: st, en                                                   ! start and end index for interpolation
    integer :: ie,i                                                     ! loop indicies
    logical :: selected                                                 ! indicates if vector is selected
    logical :: comp_selected

    output_varnames=>get_current_varnames(ios)                          ! get list of output variables

    ! check to see if any vector component is selected

    selected = .false.
    do i=1,var%n_comp
      selected = selected .or. nf_selectedvar(var%short_name(i), output_varnames)
    end do
    if(.not. selected) return                                           ! skip if no components selected

    ! Interpolate vector data

    allocate(result(n_interp,ni,2))                                     ! allocate buffer for 3d interpolated data
    st = 1                                                              ! set start index
    do ie=1,nelemd                                                      ! loop over elem
      en = st + interpdata(ie)%n_interp-1                               ! set end index
      vdata =  get_vector_field(elem(ie), var%vec_name, n0, ni)         ! get vector field data by name

      call interpolate_vector(interpdata(ie), elem(ie), vdata, ni, result(st:en,:,:), 0)
      st=st+interpdata(ie)%n_interp                                             ! shift start index
    enddo

    ! Write selected components to file
    start3d(4)=nf_get_frame(files(ios))
    count3d(4)=1
    count3d(1:3) = dim_sizes( (/1,2,3/) )                                       ! set size of output

    do i=1,var%n_comp                                                           ! loop over each vec component
      if( .not. nf_selectedvar(var%short_name(i), output_varnames) ) cycle      ! skip output if not selected

      if (hybrid%par%masterproc) print *,'writing interpolated vector component ',var%short_name(i)
      call nf_put_var(files(ios), result(:,:,i), start3d, count3d, name=var%short_name(i))
    enddo

    deallocate(result)

  end subroutine

  !_____________________________________________________________________________
  subroutine netcdf_interp_write(elem, tl, hybrid, hvcoord)

    ! Check each variable against output list and write it if present

    type (element_t),   intent(in) :: elem(:)
    type (TimeLevel_t), intent(in) :: tl
    type(hybrid_t),     intent(in) :: hybrid
    type(hvcoord_t),    intent(in) :: hvcoord

    integer :: n0                                                       ! time-level to write
    integer :: ios, i, nf_frame                                         ! loop indicies

    n0 = tl%n0

    do ios=1,max_output_streams                                         ! loop over all output streams
      if(output_frequency(ios)==0) cycle

      if(( output_start_time(ios) <= tl%nstep)  .and. &
         ( output_end_time(ios)   >= tl%nstep)  .and. &
         ( modulo(tl%nstep, output_frequency(ios)) .eq. 0)) then

        nf_frame   = nf_get_frame(files(ios))                           ! set netcdf record index
        start2d(3) = nf_frame; count2d(3)=1
        start3d(4) = nf_frame; count3d(4)=1

        if (hybrid%par%masterproc) print *,'current frame = ',start3d(4)

        ! write interpolated scalar fields

        do i=1,n_vars                                                   ! loop over all possible fields
          call write_scalar_field_1d(elem, nc_vars(i), ios, n0, hybrid, hvcoord)    ! write if scalar field and is-selected
        enddo

        ! write interpolated vector fields

        do i=1,n_vfields                                                ! loop over all vector fields
          call write_vector_field(elem,nc_vfields(i), ios, n0, hybrid, hvcoord)  ! write if at least one component selected
        enddo

        ! write simulation time

        call nf_put_var(files(ios),real(time_at(tl%nstep)/shr_const_cday,kind=real_kind),start2d(3:3),count2d(3:3),name='time')

        ! advance to next data-frame

        call nf_advance_frame(files(ios))                               ! advance output-frame number

        ! flush data to file
        call pio_syncfile(files(ios)%fileid)

      endif
    enddo

  end subroutine

  !_____________________________________________________________________________
  ! getiodof: copied unmodified from interp_movie_mod 

  subroutine GetIODOF(ndims, gdims, iorank, iodof, start, count)
    integer, intent(in) :: gdims(ndims)
    integer, intent(in) :: iorank
    integer(kind=nfsizekind), intent(out) :: start(ndims), count(ndims)
    integer, pointer :: iodof(:) ! gcc4.2 didn't like intent(out)

    integer :: nzrank, nxrank, nx, k, i, j, ndims, icnt

    if(iorank>=0) then
       nx=num_io_procs
       count(ndims)=max(1,gdims(ndims)/num_io_procs)
       nx=max(1,num_io_procs/gdims(ndims))
       nzrank=iorank/nx

       k=num_io_procs-gdims(ndims)*nx

       if(iorank>num_io_procs-k-1.and.k>0) then
          nzrank=gdims(ndims)-1
       end if

       start(ndims)=nzrank*count(ndims)+1
       if(gdims(ndims)>num_io_procs) then
          k=gdims(ndims)-num_io_procs*count(ndims)
          if(k>iorank) then
             count(ndims)=count(ndims)+1
          end if
          if(k>=iorank) then
             start(ndims)=start(ndims)+iorank
          else
             start(ndims)=start(ndims)+k
          end if

       end if

       if(k>0 .and.nzrank==gdims(ndims)-1 ) then
          nx=nx+k
       end if
       nxrank=iorank

       do i=ndims-1,1,-1
          !           print *, nxrank, nx

          nxrank=mod(nxrank,nx)
          count(i)=gdims(i)/nx
          k=gdims(i)-count(i)*nx
          if(nxrank<k) then
             count(i)=count(i)+1
             start(i)=count(i)*nxrank+1
          else
             start(i)=count(i)*nxrank+k+1
          end if
          nx=max(1,num_io_procs/gdims(i))

       end do

       icnt=0
       if(ndims.eq.1) then
          allocate(iodof(count(1)))
          do i=1,count(1)
             iodof(i)=start(1)+i-1
          end do
       else if(ndims.eq.2) then
          allocate(iodof(count(1)*count(2)))

          do j=1,count(2)
             do i=1,count(1)
                icnt=icnt+1
                iodof(icnt)=start(1)+gdims(1)*(start(2)-1)+icnt-1
             end do
          end do

       else
          allocate(iodof(count(1)*count(2)*count(3)))

          do k=1,count(3)
             do j=1,count(2)
                do i=1,count(1)
                   icnt=icnt+1
                   iodof(icnt)=start(1)+gdims(1)*(start(2)-1)+ &
                        gdims(1)*gdims(2)*(start(3)-1)+icnt-1
                end do
             end do
          end do
       end if

    else	
       allocate(iodof(1))
       iodof=-1
    end if

  end subroutine GetIODOF

end module netcdf_interp_mod
