module mo_load_coefficients
  use mo_rte_kind,              only: wp
  use mo_gas_concentrations,       only: ty_gas_concs
  use mo_gas_optics, only: ty_gas_optics
  use netcdf
  implicit none
  private
  public :: load_and_init

  ! Overload read routines
  interface read_field
    module procedure &
          read_char_vec, read_logical_vec, &
          read_real_scalar, &
          read_real_1d_field, read_real_2d_field, &
          read_real_3d_field, read_real_4d_field, &
          read_integer_scalar, &
          read_integer_1d_field, read_integer_2d_field, &
          read_integer_3d_field, read_integer_4d_field
  end interface

  ! Define a constant to use for maximum character length for consistency
  integer, parameter :: max_char_length = 256

contains

  subroutine stop_on_err(msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg

    if(msg /= "") then
      write(error_unit, *) msg
      stop
    end if
  end subroutine


  !--------------------------------------------------------------------------------------------------------------------
  ! read optical coefficients from NetCDF file
  subroutine load_and_init(kdist, filename, available_gases)
    class(ty_gas_optics), intent(inout) :: kdist
    character(len=*),                   intent(in   ) :: filename
    class(ty_gas_concs),                intent(in   ) :: available_gases ! Which gases does the host model have available?
    character(len=max_char_length), dimension(:), allocatable :: gas_names
    integer,  dimension(:,:,:), allocatable :: key_species
    integer,  dimension(:,:  ), allocatable :: band2gpt
    real(wp), dimension(:,:  ), allocatable :: band_lims
    real(wp)                                :: press_ref_trop, temp_ref_p, temp_ref_t
    real(wp), dimension(:      ), allocatable :: press_ref
    real(wp), dimension(:      ), allocatable :: temp_ref
    real(wp), dimension(:,:,:  ), allocatable :: vmr_ref
    real(wp), dimension(:,:,:,:), allocatable :: kmajor
    real(wp), dimension(:,:,:), allocatable :: kminor_lower, kminor_upper
    character(len=max_char_length), dimension(:), allocatable :: gas_minor, identifier_minor
    character(len=max_char_length), dimension(:), allocatable :: minor_gases_lower, minor_gases_upper
    integer, dimension(:,:), allocatable :: minor_limits_gpt_lower, minor_limits_gpt_upper
    logical, dimension(:), allocatable :: minor_scales_with_density_lower, minor_scales_with_density_upper
    character(len=max_char_length), dimension(:), allocatable :: scaling_gas_lower, scaling_gas_upper
    logical, dimension(:), allocatable :: scale_by_complement_lower, scale_by_complement_upper
    integer, dimension(:), allocatable :: kminor_start_lower, kminor_start_upper
    real(wp), dimension(:,:,:  ), allocatable :: rayl_lower, rayl_upper
    real(wp), dimension(:      ), allocatable :: solar_src
    real(wp), dimension(:,:    ), allocatable :: totplnk
    real(wp), dimension(:,:,:,:), allocatable :: planck_frac
    ! -----------------
    integer :: ncid, i

    ! -----------------
    ! open coefficient file
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) then 
      call stop_on_err("load_and_init(): can't open file " // trim(fileName))
    end if

    ! Read data from file. Note that the read routines will allocate the output
    ! arrays if they are not already allocated, so we should not have to
    ! allocate ahead fo time.
    call read_field(ncid, 'gas_names', gas_names)
    call read_field(ncid, 'key_species', key_species)
    call read_field(ncid, 'bnd_limits_wavenumber', band_lims)
    call read_field(ncid, 'bnd_limits_gpt', band2gpt)
    call read_field(ncid, 'press_ref', press_ref)
    call read_field(ncid, 'temp_ref',  temp_ref)
    call read_field(ncid, 'absorption_coefficient_ref_P', temp_ref_p)
    call read_field(ncid, 'absorption_coefficient_ref_T', temp_ref_t)
    call read_field(ncid, 'press_ref_trop', press_ref_trop)
    call read_field(ncid, 'kminor_lower', kminor_lower) 
    call read_field(ncid, 'kminor_upper', kminor_upper)
    call read_field(ncid, 'gas_minor', gas_minor)
    call read_field(ncid, 'identifier_minor', identifier_minor)
    call read_field(ncid, 'minor_gases_lower', minor_gases_lower)
    call read_field(ncid, 'minor_gases_upper', minor_gases_upper)
    call read_field(ncid, 'minor_limits_gpt_lower', minor_limits_gpt_lower)
    call read_field(ncid, 'minor_limits_gpt_upper', minor_limits_gpt_upper)
    call read_field(ncid, 'minor_scales_with_density_lower', minor_scales_with_density_lower)
    call read_field(ncid, 'minor_scales_with_density_upper', minor_scales_with_density_upper)
    call read_field(ncid, 'scale_by_complement_lower', scale_by_complement_lower)
    call read_field(ncid, 'scale_by_complement_upper', scale_by_complement_upper)
    call read_field(ncid, 'scaling_gas_lower', scaling_gas_lower)
    call read_field(ncid, 'scaling_gas_upper', scaling_gas_upper)
    call read_field(ncid, 'kminor_start_lower', kminor_start_lower)
    call read_field(ncid, 'kminor_start_upper', kminor_start_upper)
    call read_field(ncid, 'vmr_ref', vmr_ref)
    call read_field(ncid, 'kmajor', kmajor)

    if(var_exists(ncid, 'rayl_lower') .and. var_exists(ncid, 'rayl_upper')) then
      call read_field(ncid, 'rayl_lower', rayl_lower)
      call read_field(ncid, 'rayl_upper', rayl_upper)
    end if

    if(var_exists(ncid, 'totplnk') .and. var_exists(ncid, 'plank_fraction')) then
      !
      ! If there's a totplnk variable in the file it's a longwave (internal sources) type
      !
      call read_field(ncid, 'totplnk', totplnk)
      call read_field(ncid, 'plank_fraction', planck_frac)
      call stop_on_err(kdist%load(available_gases, &
                                  gas_names,   &
                                  key_species, &
                                  band2gpt,    &
                                  band_lims,   &
                                  press_ref,   &
                                  press_ref_trop, &
                                  temp_ref,    &
                                  temp_ref_p, temp_ref_t,     &
                                  vmr_ref, kmajor,            &
                                  kminor_lower, kminor_upper, &
                                  gas_minor,identifier_minor, &
                                  minor_gases_lower, minor_gases_upper, &
                                  minor_limits_gpt_lower, &
                                  minor_limits_gpt_upper, &
                                  minor_scales_with_density_lower, &
                                  minor_scales_with_density_upper, &
                                  scaling_gas_lower, scaling_gas_upper, &
                                  scale_by_complement_lower, &
                                  scale_by_complement_upper, &
                                  kminor_start_lower, &
                                  kminor_start_upper, &
                                  totplnk, planck_frac,       &
                                  rayl_lower, rayl_upper))
      deallocate(totplnk, planck_frac)
    else if (var_exists(ncid, 'solar_source')) then
      !
      ! Solar source doesn't have an dependencies yet
      !
      call read_field(ncid, 'solar_source', solar_src)
      call stop_on_err(kdist%load(available_gases, &
                                  gas_names,   &
                                  key_species, &
                                  band2gpt,    &
                                  band_lims,   &
                                  press_ref,   &
                                  press_ref_trop, &
                                  temp_ref,    &
                                  temp_ref_p, temp_ref_t,     &
                                  vmr_ref, kmajor,            &
                                  kminor_lower, kminor_upper, &
                                  gas_minor,identifier_minor,&
                                  minor_gases_lower, minor_gases_upper, &
                                  minor_limits_gpt_lower, &
                                  minor_limits_gpt_upper, &
                                  minor_scales_with_density_lower, &
                                  minor_scales_with_density_upper, &
                                  scaling_gas_lower, scaling_gas_upper, &
                                  scale_by_complement_lower, &
                                  scale_by_complement_upper, &
                                  kminor_start_lower, &
                                  kminor_start_upper, &
                                  solar_src, &
                                  rayl_lower, rayl_upper))
       deallocate(solar_src)
    else
       call stop_on_err( &
         'load_and_init: no totplnk, no planck_fraction, no solar_src.' &
       )
    end if

    ! Free up allocated memory
    deallocate(gas_names, key_species, band2gpt, band_lims, &
               press_ref, temp_ref, vmr_ref, &
               kmajor, kminor_lower, kminor_upper, gas_minor, &
               identifier_minor, minor_gases_lower, minor_gases_upper, &
               minor_limits_gpt_lower, minor_limits_gpt_upper, &
               minor_scales_with_density_lower, &
               minor_scales_with_density_upper, &
               scaling_gas_lower, scaling_gas_upper, &
               scale_by_complement_lower, scale_by_complement_upper, &
               kminor_start_lower, kminor_start_upper)

    if (allocated(rayl_lower)) deallocate(rayl_lower)
    if (allocated(rayl_upper)) deallocate(rayl_upper)

    ! Close file
    ncid = nf90_close(ncid)
  end subroutine load_and_init
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Ancillary functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_integer_scalar(ncid, varName, vardata)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer, intent(inout)       :: vardata 

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end subroutine read_integer_scalar
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_integer_1d_field(ncid, varName, vardata)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer, allocatable, intent(inout)  :: vardata(:)

    integer :: varsizes(1)
    integer :: varid

    ! Check sizes and allocate if needed
    varsizes = get_data_size(ncid, varName, 1)
    if (allocated(vardata)) then
       if (any(shape(vardata) /= varsizes)) then
          call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent." )
       end if
    else
       allocate(vardata(varsizes(1)))
    end if

    ! Get variable ID from name
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))

    ! Read data
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end subroutine read_integer_1d_field
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_integer_2d_field(ncid, varName, vardata)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer, allocatable, intent(inout) :: vardata(:,:)

    integer :: varsizes(2)
    integer :: varid

    ! Check sizes and allocate if needed
    varsizes = get_data_size(ncid, varName, 2)
    if (allocated(vardata)) then
       if (any(shape(vardata) /= varsizes)) then
          call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent." )
       end if
    else
       allocate(vardata(varsizes(1), varsizes(2)))
    end if

    ! Get variable ID from name
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))

    ! Read data
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end subroutine read_integer_2d_field
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_integer_3d_field(ncid, varName, vardata)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer, allocatable, intent(inout) :: vardata(:,:,:)

    integer :: varsizes(3)
    integer :: varid

    ! Check sizes and allocate if needed
    varsizes = get_data_size(ncid, varName, 3)
    if (allocated(vardata)) then
       if (any(shape(vardata) /= varsizes)) then
          print *, 'shape(vardata) = ', shape(vardata), 'varsizes = ', varsizes
          call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent." )
       end if
    else
       allocate(vardata(varsizes(1), varsizes(2), varsizes(3)))
    end if

    ! Get variable ID from name
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))

    ! Read data
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end subroutine read_integer_3d_field
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_integer_4d_field(ncid, varName, vardata)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer, allocatable, intent(inout) :: vardata(:,:,:,:)

    integer :: varsizes(4)
    integer :: varid

    ! Check sizes and allocate if needed
    varsizes = get_data_size(ncid, varName, 4)
    if (allocated(vardata)) then
       if (any(shape(vardata) /= varsizes)) then
          call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent." )
       end if
    else
       allocate(vardata(varsizes(1), varsizes(2), varsizes(3), varsizes(4)))
    end if

    ! Get variable ID from name
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))

    ! Read data
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end subroutine read_integer_4d_field
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_real_scalar(ncid, varName, vardata)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    real(wp), intent(inout) :: vardata 

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end subroutine read_real_scalar
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_real_1d_field(ncid, varName, vardata)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    real(wp), allocatable, intent(inout) :: vardata(:)

    integer :: varsizes(1)
    integer :: varid

    ! Check sizes and allocate if needed
    varsizes = get_data_size(ncid, varName, 1)
    if (allocated(vardata)) then
       if (any(shape(vardata) /= varsizes)) then
          call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent." )
       end if
    else
       allocate(vardata(varsizes(1)))
    end if

    ! Get variable ID from name
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))

    ! Read data
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end subroutine read_real_1d_field
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_real_2d_field(ncid, varName, vardata)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    real(wp), allocatable, intent(inout) :: vardata(:,:)

    integer :: varsizes(2)
    integer :: varid

    ! Check sizes and allocate if needed
    varsizes = get_data_size(ncid, varName, 2)
    if (allocated(vardata)) then
       if (any(shape(vardata) /= varsizes)) then
          call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent." )
       end if
    else
       allocate(vardata(varsizes(1), varsizes(2)))
    end if

    ! Get variable ID from name
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))

    ! Read data
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end subroutine read_real_2d_field
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_real_3d_field(ncid, varName, vardata)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    real(wp), allocatable, intent(inout) :: vardata(:,:,:)

    integer :: varsizes(3)
    integer :: varid

    ! Check sizes and allocate if needed
    varsizes = get_data_size(ncid, varName, 3)
    if (allocated(vardata)) then
       if (any(shape(vardata) /= varsizes)) then
          call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent." )
       end if
    else
       allocate(vardata(varsizes(1), varsizes(2), varsizes(3)))
    end if

    ! Get variable ID from name
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))

    ! Read data
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end subroutine read_real_3d_field
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_real_4d_field(ncid, varName, vardata)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    real(wp), allocatable, intent(inout) :: vardata(:,:,:,:)

    integer :: varsizes(4)
    integer :: varid

    ! Check sizes and allocate if needed
    varsizes = get_data_size(ncid, varName, 4)
    if (allocated(vardata)) then
       if (any(shape(vardata) /= varsizes)) then
          call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent." )
       end if
    else
       allocate(vardata(varsizes(1), varsizes(2), varsizes(3), varsizes(4)))
    end if

    ! Get variable ID from name
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))

    ! Read data
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end subroutine read_real_4d_field
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_logical_vec(ncid, varName, vardata)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    logical, allocatable, intent(inout) :: vardata(:)
    integer, allocatable :: vardata_tmp(:)

    integer :: varid
    integer :: nx, ix
    integer :: var_sizes(1)

    ! Get variable ID from name
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_logical_vec: can't find variable " // trim(varName))

    ! Get variable dimension sizes
    var_sizes = get_data_size(ncid, varName, 1)

    ! Allocate space for temporary variable
    allocate(vardata_tmp(var_sizes(1)))

    ! Read temporary variable
    if(nf90_get_var(ncid, varid, vardata_tmp)  /= NF90_NOERR) &
      call stop_on_err("read_logical_vec: can't read variable " // trim(varName))

    ! Check if vardata is already allocated; if it is, check sizes. If not,
    ! allocate now.
    if (allocated(vardata)) then
       if (any(shape(vardata) /= var_sizes)) then
          call stop_on_err("read_logical_vec: inconsistent sizes for " // trim(varName))
       end if
    else
       allocate(vardata(var_sizes(1)))
    end if

    ! Convert temporary variable to logical
    do ix = 1, var_sizes(1)
      if (vardata_tmp(ix) .eq. 0) then
        vardata(ix) = .false.
      else
        vardata(ix) = .true.
      endif
    enddo

    deallocate(vardata_tmp)
  end subroutine read_logical_vec
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_char_vec(ncid, varName, vardata)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    character(len=*), allocatable, intent(inout) :: vardata(:)

    ! var_sizes needs to be length 2, because one dimension is reserved for the
    ! character length.
    integer :: var_sizes(2)
    integer :: varid

    ! Get variable ID from name
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) then
      call stop_on_err("read_char_vec: can't find variable " // trim(varName))
    end if

    ! Get dimension sizes; 1st dimension is string length
    var_sizes = get_data_size(ncid, varName, 2)

    ! Check if variable is allocated; if not, allocate. If allocated, check that
    ! sizes match variable sizes on disk. Note that var_sizes(1) is the string
    ! length, so we want to compare against var_sizes(2).
    if (allocated(vardata)) then
       if (size(vardata, 1) /= var_sizes(2)) then
          call stop_on_err("read_char_vec: inconsistent sizes for " // trim(varName))
       end if
    else
       allocate(vardata(var_sizes(2)))
    end if

    ! Finally, read data
    if(nf90_get_var(ncid, varid, vardata)  /= NF90_NOERR) &
      call stop_on_err("read_char_vec: can't read variable " // trim(varName))

  end subroutine read_char_vec
  !--------------------------------------------------------------------------------------------------------------------
  function var_exists(ncid, varName)
    !
    ! Does this variable exist (have a valid var_id) in the open netCDF file?
    !
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    logical :: var_exists

    integer :: varId
    var_exists = nf90_inq_varid(ncid, trim(varName), varid) == NF90_NOERR
  end function var_exists
  !--------------------------------------------------------------------------------------------------------------------
  function get_dim_length(ncid, dimname)
    !
    ! Get the length of a dimension from an open netCDF file
    !  This is unfortunately a two-step process
    !
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer :: get_dim_length

    integer :: dimid

    if(nf90_inq_dimid(ncid, trim(dimname), dimid) == NF90_NOERR) then
      if(nf90_inquire_dimension(ncid, dimid, len=get_dim_length) /= NF90_NOERR) get_dim_length = 0
    else
      get_dim_length = 0
    end if

  end function get_dim_length
  !--------------------------------------------------------------------------------------------------------------------
  function get_data_size(ncid, varName, n)
    !
    ! Returns the extents of a netcdf variable on disk
    !
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer, intent(in) :: n
    integer :: get_data_size(n)

    integer :: i
    integer :: varid, ndims, dimids(n)

    get_data_size(n) = -1
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't find variable " // trim(varName))
    if(nf90_inquire_variable(ncid, varid, ndims = ndims) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't get information for variable " // trim(varName))
    if(ndims /= n) &
      call stop_on_err("get_data_size:  variable " // trim(varName) // " has the wrong number of dimensions" )
    if(nf90_inquire_variable(ncid, varid, dimids = dimids) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't read dimension ids for variable " // trim(varName))
    do i = 1, n
      if(nf90_inquire_dimension(ncid, dimids(i), len = get_data_size(i)) /= NF90_NOERR) &
        call stop_on_err("get_data_size: can't get dimension lengths for variable " // trim(varName))
    end do

  end function get_data_size
  !--------------------------------------------------------------------------------------------------------------------
end module
