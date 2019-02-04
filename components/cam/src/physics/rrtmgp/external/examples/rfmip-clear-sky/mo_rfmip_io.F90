! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! This module reads an example file containing atomspheric conditions (temperature, pressure, gas concentrations)
!   and surface properties (emissivity, temperature), defined on nlay layers across a set of ncol columns subject to
!   nexp perturbations, and returns them in data structures suitable for use in rte and rrtmpg. The input data
!   are partitioned into a user-specified number of blocks.
! For the moment only quantities relevant to longwave calculations are provided.
!
! The example files comes from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
!   The protocol for this experiment allows for different specifications of which gases to consider:
! all gases, (CO2, CH4, N2O) + {CFC11eq; CFC12eq + HFC-134eq}. Ozone is always included
! The protocol does not specify the treatmet of gases like CO
!
! -------------------------------------------------------------------------------------------------
module mo_rfmip_io
  use mo_rte_kind,   only: wp
  use mo_gas_concentrations, &
                        only: ty_gas_concs
  use mo_util_string,   only: lower_case, string_in_array, string_loc_in_array
  use mo_simple_netcdf, only: read_field, write_field, get_dim_size
  use netcdf
  implicit none
  private
  public :: read_kdist_gas_names, read_size, read_and_block_pt, &
            read_and_block_sw_bc, read_and_block_lw_bc, read_and_block_gases_ty, &
            unblock_and_write

  integer :: ncol_l = 0, nlay_l = 0, nexp_l = 0 ! Local copies
contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read the names of the gases known to the k-distribution
  !
  !
  subroutine read_kdist_gas_names(fileName, kdist_gas_names)
    character(len=*),          intent(in   ) :: fileName
    character(len=32), dimension(:), allocatable, &
                               intent(  out) :: kdist_gas_names
    ! ---------------------------
    integer :: ncid, varid
    character(len=9), parameter :: varName = "gas_names"
    ! ---------------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_kdist_gas_names: can't open file " // trim(fileName))

    allocate(kdist_gas_names(get_dim_size(ncid, 'absorber')))

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_kdist_gas_names: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, kdist_gas_names)  /= NF90_NOERR) &
      call stop_on_err("read_kdist_gas_names: can't read variable " // trim(varName))

    ncid = nf90_close(ncid)
  end subroutine read_kdist_gas_names
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Find the size of the problem: columns, layers, perturbations (experiments)
  !
  subroutine read_size(fileName, ncol, nlay, nexp)
    character(len=*),          intent(in   ) :: fileName
    integer,         optional, intent(  out) :: ncol, nlay, nexp
    ! ---------------------------
    integer :: ncid
    ! ---------------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_size: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'site')
    nlay = get_dim_size(ncid, 'layer')
    nexp = get_dim_size(ncid, 'expt')
    if(get_dim_size(ncid, 'level') /= nlay+1) call stop_on_err("read_size: number of levels should be nlay+1")
    ncid = nf90_close(ncid)

    ncol_l = ncol
    nlay_l = nlay
    nexp_l = nexp
  end subroutine read_size
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Return layer and level pressures and temperatures as arrays dimensioned (ncol, nlay/+1, nblocks)
  !   Input arrays are dimensioned (nlay/+1, ncol, nexp)
  !   Output arrays are allocated within this routine
  !
  subroutine read_and_block_pt(fileName, blocksize, &
                               p_lay, p_lev, t_lay, t_lev)
    character(len=*),           intent(in   ) :: fileName
    integer,                    intent(in   ) :: blocksize
    real(wp), dimension(:,:,:), allocatable, & ! [blocksize, nlay/+1, nblocks]
                                intent(  out) :: p_lay, p_lev, t_lay, t_lev
    ! ---------------------------
    integer :: ncid
    integer :: b, nblocks
    real(wp), dimension(:,:  ), allocatable :: temp2d
    real(wp), dimension(:,:,:), allocatable :: temp3d
    ! ---------------------------
    if(any([ncol_l, nlay_l, nexp_l]  == 0)) call stop_on_err("read_and_block_pt: Haven't read problem size yet.")
    if(mod(ncol_l*nexp_l, blocksize) /= 0 ) call stop_on_err("read_and_block_pt: number of columns doesn't fit evenly into blocks.")
    nblocks = (ncol_l*nexp_l)/blocksize
    allocate(p_lay(blocksize, nlay_l,   nblocks), t_lay(blocksize, nlay_l,   nblocks), &
             p_lev(blocksize, nlay_l+1, nblocks), t_lev(blocksize, nlay_l+1, nblocks))

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_and_block_pt: can't find file " // trim(fileName))
    !
    ! Read p, T data; reshape to suit RRTMGP dimensions
    !
    temp3d = reshape(spread(read_field(ncid, "pres_layer", nlay_l,   ncol_l), dim = 3, ncopies = nexp_l), &
                     shape = [nlay_l, blocksize, nblocks])
    do b = 1, nblocks
      p_lay(:,:,b) = transpose(temp3d(:,:,b))
    end do
    temp3d = reshape(       read_field(ncid, "temp_layer", nlay_l,   ncol_l, nexp_l), &
                     shape = [nlay_l, blocksize, nblocks])
    do b = 1, nblocks
      t_lay(:,:,b) = transpose(temp3d(:,:,b))
    end do

    deallocate(temp3d)
    temp3d = reshape(spread(read_field(ncid, "pres_level", nlay_l+1, ncol_l),  dim = 3, ncopies = nexp_l), &
                    shape = [nlay_l+1, blocksize, nblocks])
    do b = 1, nblocks
      p_lev(:,:,b) = transpose(temp3d(:,:,b))
    end do
    temp3d = reshape(       read_field(ncid, "temp_level", nlay_l+1, ncol_l, nexp_l), &
                    shape = [nlay_l+1, blocksize, nblocks])
    do b = 1, nblocks
      t_lev(:,:,b) = transpose(temp3d(:,:,b))
    end do

    ncid = nf90_close(ncid)
  end subroutine read_and_block_pt
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read and reshape shortwave boundary conditions
  !
  subroutine read_and_block_sw_bc(fileName, blocksize, &
                               surface_albedo, total_solar_irradiance, solar_zenith_angle)
    character(len=*),           intent(in   ) :: fileName
    integer,                    intent(in   ) :: blocksize
    real(wp), dimension(:,:), allocatable, &
                                intent(  out) :: surface_albedo, total_solar_irradiance, solar_zenith_angle
    ! ---------------------------
    integer :: ncid
    integer :: nblocks
    ! ---------------------------
    if(any([ncol_l, nlay_l, nexp_l]  == 0)) call stop_on_err("read_and_block_sw_bc: Haven't read problem size yet.")
    if(mod(ncol_l*nexp_l, blocksize) /= 0 ) call stop_on_err("read_and_block_sw_bc: number of columns doesn't fit evenly into blocks.")
    nblocks = (ncol_l*nexp_l)/blocksize
    !
    ! Check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
    !

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_and_block_sw_bc: can't find file " // trim(fileName))

    surface_albedo         = reshape(spread(read_field(ncid, "surface_albedo",         ncol_l), dim=2, ncopies=nexp_l), &
                                     shape = [blocksize, nblocks])
    total_solar_irradiance = reshape(spread(read_field(ncid, "total_solar_irradiance", ncol_l), dim=2, ncopies=nexp_l), &
                                    shape = [blocksize, nblocks])
    solar_zenith_angle     = reshape(spread(read_field(ncid, "solar_zenith_angle",     ncol_l), dim=2, ncopies=nexp_l), &
                                     shape = [blocksize, nblocks])

    ncid = nf90_close(ncid)
  end subroutine read_and_block_sw_bc
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read and reshape longwave boundary conditions
  !
  subroutine read_and_block_lw_bc(fileName, blocksize, &
                                  surface_emissivity, surface_temperature)
    character(len=*),           intent(in   ) :: fileName
    integer,                    intent(in   ) :: blocksize
    real(wp), dimension(:,:), allocatable, &
                                intent(  out) :: surface_emissivity, surface_temperature
    ! ---------------------------
    integer :: ncid
    integer :: nblocks
    ! ---------------------------
    if(any([ncol_l, nlay_l, nexp_l]  == 0)) &
      call stop_on_err("read_and_block_lw_bc: Haven't read problem size yet.")
    if(mod(ncol_l*nexp_l, blocksize) /= 0 ) &
      call stop_on_err("read_and_block_lw_bc: number of columns doesn't fit evenly into blocks.")
    nblocks = (ncol_l*nexp_l)/blocksize
    !
    ! Check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
    !

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_and_block_lw_bc: can't find file " // trim(fileName))

    surface_emissivity  = reshape(spread(read_field(ncid, "surface_emissivity",  ncol_l), dim=2, ncopies=nexp_l), &
                                  shape = [blocksize, nblocks])
    surface_temperature = reshape(spread(read_field(ncid, "surface_temperature", ncol_l), dim=2, ncopies=nexp_l), &
                                  shape = [blocksize, nblocks])

    ncid = nf90_close(ncid)
  end subroutine read_and_block_lw_bc
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read and reshape gas concentrations. RRTMGP requires gas concentrations to be supplied via a class
  !   (ty_gas_concs). Gas concentrations are set via a call to gas_concs%set_vmr(name, values)
  !   where `name` is nominally the chemical formula for the gas in question and `values` may be
  !   a scalar, a 1-d profile assumed to apply to all columns, or an array of dimension (ncol, nlay).
  ! This routine outputs a vector nblocks long of these types so each element of the array can be passed to
  !   the rrtmgp gas optics calculation in turn.
  !
  ! This routine exploits RFMIP conventions: only water vapor and ozone vary by column within
  !   each experiment.
  ! Fields in the RFMIP file have a trailing _GM (global mean); some fields use a chemical formula and other
  !   a descriptive name, so a map is provided between these.
  !
  subroutine read_and_block_gases_ty(fileName, blocksize, gas_names, gas_conc_array)
    character(len=*),           intent(in   ) :: fileName
    integer,                    intent(in   ) :: blocksize
    character(len=*),  dimension(:), &
                                intent(in   ) :: gas_names ! Provided by gas_optics -- which gases do we want to read
    type(ty_gas_concs), dimension(:), allocatable, &
                                intent(  out) :: gas_conc_array

    ! ---------------------------
    integer :: ncid
    integer :: nblocks
    integer :: b, g
    integer,  dimension(:,:),   allocatable :: exp_num
    real(wp), dimension(:),     allocatable :: gas_conc_temp_1d
    real(wp), dimension(:,:,:), allocatable :: gas_conc_temp_3d
    character(len=32)                       :: gas_name_in_file
    character(len=32), dimension(10) :: &
      chem_name = ['co   ', &
                   'ch4  ', &
        				   'o2   ', &
        				   'n2o  ', &
        				   'n2   ', &
        				   'co2  ', &
        				   'CCl4 ', &
        				   'ch4  ', &
        				   'CH3Br', &
   			           'CH3Cl'], &
      desc_name = ['carbon_monoxide     ', &
                   'methane             ', &
                   'oxygen              ', &
          			   'nitrous_oxide       ', &
          			   'nitrogen            ', &
        				   'carbon_dioxide      ', &
        				   'carbon_tetrachloride', &
        				   'methane             ', &
        				   'methyl_bromide      ', &
        				   'methyl_chloride     ']
    ! ---------------------------
    if(any([ncol_l, nlay_l, nexp_l]  == 0)) &
      call stop_on_err("read_and_block_lw_bc: Haven't read problem size yet.")
    if(mod(ncol_l*nexp_l, blocksize) /= 0 ) &
      call stop_on_err("read_and_block_lw_bc: number of columns doesn't fit evenly into blocks.")
    nblocks = (ncol_l*nexp_l)/blocksize
    allocate(gas_conc_array(nblocks))
    ! Experiment index for each colum
    exp_num = reshape(spread([(b, b = 1, nexp_l)], 1, ncopies = ncol_l), shape = [blocksize, nblocks], order=[1,2])

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_and_block_gases_ty: can't find file " // trim(fileName))
    !
    ! Water vapor and ozone depend on col, lay, exp: look just like other fields
    !
    gas_conc_temp_3d = reshape(read_field(ncid, "water_vapor", nlay_l, ncol_l, nexp_l), &
                               shape = [nlay_l, blocksize, nblocks]) * read_scaling(ncid, "water_vapor")
    do b = 1, nblocks
      call stop_on_err(gas_conc_array(b)%set_vmr('h2o', transpose(gas_conc_temp_3d(:,:,b))))
    end do

    gas_conc_temp_3d = reshape(read_field(ncid, "ozone", nlay_l, ncol_l, nexp_l), &
                               shape = [nlay_l, blocksize, nblocks]) * read_scaling(ncid, "ozone")
    do b = 1, nblocks
      call stop_on_err(gas_conc_array(b)%set_vmr('o3', transpose(gas_conc_temp_3d(:,:,b))))
    end do

    !
    ! All other gases are a function of experiment only
    !
    do g = 1, size(gas_names)
      gas_name_in_file = trim(lower_case(gas_names(g)))
      if(gas_name_in_file == 'h2o' .or. gas_name_in_file == 'o3') cycle
      !
      ! Use a mapping between chemical formula and name if it exists
      !
      if(string_in_array(gas_name_in_file, chem_name)) &
        gas_name_in_file = desc_name(string_loc_in_array(gas_name_in_file, chem_name))
      gas_name_in_file = trim(gas_name_in_file) // "_GM"

      ! Read the values as a function of experiment
      gas_conc_temp_1d = read_field(ncid, gas_name_in_file, nexp_l) * read_scaling(ncid, gas_name_in_file)

	  do b = 1, nblocks
        ! Does every value in this block belong to the same experiment?
	    if(all(exp_num(1,b) == exp_num(2:,b))) then
	      ! Provide a scalar value
		    call stop_on_err(gas_conc_array(b)%set_vmr(gas_names(g), gas_conc_temp_1d(exp_num(1,b))))
		  else
		  ! Create 2D field, blocksize x nlay, with scalar values from each experiment
		  call stop_on_err(gas_conc_array(b)%set_vmr(gas_names(g), &
		                                             spread(gas_conc_temp_1d(exp_num(:,b)), 2, ncopies = nlay_l)))
		  end if
	  end do

    end do
    ncid = nf90_close(ncid)
  end subroutine read_and_block_gases_ty
  !--------------------------------------------------------------------------------------------------------------------
  function read_scaling(ncid, varName)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    real(wp)                     :: read_scaling

    integer           :: varid
    character(len=16) :: charUnits

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_scaling: can't find variable " // trim(varName))
    if(nf90_get_att(ncid, varid, "units", charUnits)  /= NF90_NOERR) &
      call stop_on_err("read_scaling: can't read attribute 'units' from variable " // trim(varName))
    read(charUnits, *) read_scaling
    return

  end function read_scaling
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Reshape and reorder values (nominally fluxes) from RTE order (ncol, nlev, nblocks)
  !   to RFMIP order (nlev, ncol, nexp), then write them to a user-specified variable
  !   in a netCDF file.
  !
  subroutine unblock_and_write(fileName, varName, values)
    character(len=*),           intent(in   ) :: fileName, varName
    real(wp), dimension(:,:,:),  & ! [blocksize, nlay/+1, nblocks]
                                intent(in   ) :: values
    ! ---------------------------
    integer :: ncid
    integer :: b, blocksize, nlev, nblocks
    real(wp), dimension(:,:), allocatable :: temp2d
    ! ---------------------------
    if(any([ncol_l, nlay_l, nexp_l]  == 0)) call stop_on_err("unblock_and_write: Haven't read problem size yet.")
    blocksize = size(values,1)
    nlev      = size(values,2)
    nblocks   = size(values,3)
    if(nlev /= nlay_l+1)                   call stop_on_err('unblock_and_write: array values has the wrong number of levels')
    if(blocksize*nblocks /= ncol_l*nexp_l) call stop_on_err('unblock_and_write: array values has the wrong number of blocks/size')

    allocate(temp2D(nlev, ncol_l*nexp_l))
    do b = 1, nblocks
      temp2D(1:nlev, ((b-1)*blocksize+1):(b*blocksize)) = transpose(values(1:blocksize,1:nlev,b))
    end do
    !
    ! Check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
    !

    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("unblock_and_write: can't find file " // trim(fileName))
    call stop_on_err(write_field(ncid, varName,  &
                                 reshape(temp2d, shape = [nlev, ncol_l, nexp_l])))

    ncid = nf90_close(ncid)
  end subroutine unblock_and_write
  !--------------------------------------------------------------------------------------------------------------------
  subroutine stop_on_err(msg)
    !
    ! Print error message and stop
    !
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg
    if(len_trim(msg) > 0) then
      write(error_unit,*) trim(msg)
      stop
    end if
  end subroutine
end module mo_rfmip_io
