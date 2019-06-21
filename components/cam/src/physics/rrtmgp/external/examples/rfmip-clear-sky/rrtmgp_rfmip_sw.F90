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
! Example program to demonstrate the calculation of shortwave radiative fluxes in clear, aerosol-free skies.
!   The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
!   The large problem (1800 profiles) is divided into blocks
!
! Program is invoked as rrtmgp_rfmip_sw [block_size input_file  coefficient_file upflux_file downflux_file]
!   All arguments are optional but need to be specified in order.
!
! -------------------------------------------------------------------------------------------------
!
! Error checking: Procedures in rte+rrtmgp return strings which are empty if no errors occured
!   Check the incoming string, print it out and stop execution if non-empty
!
subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "rrtmgp_rfmip_sw stopping"
    stop
  end if
end subroutine stop_on_err
! -------------------------------------------------------------------------------------------------
!
! Main program
!
! -------------------------------------------------------------------------------------------------
program rrtmgp_rfmip_sw
  ! --------------------------------------------------
  !
  ! Modules for working with rte and rrtmgp
  !
  ! Working precision for real variables
  !
  use mo_rte_kind,           only: wp
  !
  ! Optical properties of the atmosphere as array of values
  !   In the longwave we include only absorption optical depth (_1scl)
  !   Shortwave calculations use optical depth, single-scattering albedo, asymmetry parameter (_2str)
  !
  use mo_optical_props,      only: ty_optical_props_2str
  !
  ! Gas optics: maps physical state of the atmosphere to optical properties
  !
  use mo_gas_optics,         only: ty_gas_optics
  !
  ! Gas optics uses a derived type to represent gas concentrations compactly
  !
  use mo_gas_concentrations, only: ty_gas_concs
  !
  ! RTE shortwave driver
  !
  use mo_rte_sw,             only: rte_sw
  !
  ! RTE driver uses a derived type to reduce spectral fluxes to whatever the user wants
  !   Here we're just reporting broadband fluxes
  !
  use mo_fluxes,             only: ty_fluxes_broadband
  ! --------------------------------------------------
  !
  ! modules for reading and writing files
  !
  ! RRTMGP's gas optics class needs to be initialized with data read from a netCDF files
  !
  use mo_load_coefficients,  only: load_and_init
  use mo_rfmip_io,           only: read_size, read_and_block_pt, read_and_block_gases_ty, unblock_and_write, &
                                   read_and_block_sw_bc, read_kdist_gas_names
#IFDEF USE_TIMING
  !
  ! Timing library
  !
  use gptl,                  only: gptlstart, gptlstop, gptlinitialize, gptlpr, gptlfinalize, gptlsetoption, &
                                   gptlpercent, gptloverhead
#ENDIF
  implicit none
  ! --------------------------------------------------
  !
  ! Local variables
  !
  character(len=132)         :: rfmip_file = 'multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-0-4_none.nc', &
                                kdist_file = 'coefficients_sw.nc', &
                                flxdn_file = 'rsd_template.nc', flxup_file = 'rsu_template.nc'
  integer                    :: nargs, ncol, nlay, nexp, nblocks, block_size
  logical                    :: top_at_1
  integer                    :: b, icol, igpt
  character(len=6)           :: block_size_char

  character(len=32 ), &
            dimension(:),             allocatable :: kdist_gas_names, gases_to_use
  real(wp), dimension(:,:,:),         allocatable :: p_lay, p_lev, t_lay, t_lev ! block_size, nlay, nblocks
  real(wp), dimension(:,:,:), target, allocatable :: flux_up, flux_dn
  real(wp), dimension(:,:  ),         allocatable :: surface_albedo, total_solar_irradiance, solar_zenith_angle
                                                     ! block_size, nblocks
  !
  ! Classes used by rte+rrtmgp
  !
  type(ty_gas_optics)                            :: k_dist
  type(ty_optical_props_2str)                    :: optical_props
  type(ty_fluxes_broadband)                      :: fluxes
  real(wp), dimension(:,:), allocatable          :: toa_flux ! block_size, ngpt
  real(wp), dimension(:  ), allocatable          :: def_tsi, mu0 ! block_size
  logical , dimension(:  ), allocatable          :: usecol ! block_size, ngpt
  !
  ! ty_gas_concentration holds multiple columns; we make an array of these objects to
  !   leverage what we know about the input file
  !
  type(ty_gas_concs), dimension(:), allocatable  :: gas_conc_array

  integer :: ret
  ! -------------------------------------------------------------------------------------------------
  !
  ! Code starts
  ! Argument list:
  !   block size, input file, coefficient file, upflux file, downflux file
  !   all arguments are optional
  !
  nargs = command_argument_count()
  if(nargs >= 2) call get_command_argument(2, rfmip_file)
  if(nargs >= 3) call get_command_argument(3, kdist_file)
  if(nargs >= 4) call get_command_argument(4, flxup_file)
  if(nargs >= 5) call get_command_argument(5, flxdn_file)

  ! How big is the problem? Does it fit into blocks of the size we've specified?
  !
  call read_size(rfmip_file, ncol, nlay, nexp)
  if(nargs >= 1) then
    call get_command_argument(1, block_size_char)
    read(block_size_char, '(i6)') block_size
  else
    block_size = ncol
  end if
  if(mod(ncol*nexp, block_size) /= 0 ) call stop_on_err("rrtmgp_rfmip_sw: number of columns doesn't fit evenly into blocks.")
  nblocks = (ncol*nexp)/block_size
  print *, "Doing ",  nblocks, "blocks of size ", block_size

  !
  ! Names of gases known to the k-distribution.
  !
  call read_kdist_gas_names(kdist_file, kdist_gas_names)
  !
  ! Which gases will be included in the calculation?
  !    By default we'll use all the gases the k-distribution can handle, but
  !    we could provide variants i.e. using equivalent concentrations per RFMIP
  !
  gases_to_use = kdist_gas_names
  print *, "Radiation calculation uses gases "
  print *, "  ", [(trim(gases_to_use(b)) // " ", b = 1, size(gases_to_use))]

  ! --------------------------------------------------
  !
  ! Prepare data for use in rte+rrtmgp
  !
  !
  ! Allocation on assignment within reading routines
  !
  call read_and_block_pt(rfmip_file, block_size, p_lay, p_lev, t_lay, t_lev)
  !
  ! Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
  !
  top_at_1 = p_lay(1, 1, 1) < p_lay(1, nlay, 1)

  !
  ! Read the gas concentrations and surface properties
  !
  call read_and_block_gases_ty(rfmip_file, block_size, gases_to_use, gas_conc_array)
  call read_and_block_sw_bc(rfmip_file, block_size, surface_albedo, total_solar_irradiance, solar_zenith_angle)
  !
  ! Read k-distribution information. load_and_init() reads data from netCDF and calls
  !   k_dist%init(); users might want to use their own reading methods
  !
  call load_and_init(k_dist, trim(kdist_file), gas_conc_array(1))
  if(.not. k_dist%source_is_external()) &
    stop "rrtmgp_rfmip_sw: k-distribution file isn't SW"

  allocate(toa_flux(block_size, k_dist%get_ngpt()), &
           def_tsi(block_size), mu0(block_size), usecol(block_size))
  !
  ! RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
  !   is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
  !   This introduces an error but shows input sanitizing.
  !
  if(top_at_1) then
    p_lev(:,1,:) = k_dist%get_press_ref_min() + epsilon(k_dist%get_press_ref_min())
  else
    p_lev(:,nlay+1,:) &
                 = k_dist%get_press_ref_min() + epsilon(k_dist%get_press_ref_min())
  end if

  !
  ! Allocate space for output fluxes (accessed via pointers in ty_fluxes_broadband),
  !   gas optical properties, and source functions. The %alloc() routines carry along
  !   the spectral discretization from the k-distribution.
  !
  allocate(flux_up(block_size, nlay+1, nblocks), &
           flux_dn(block_size, nlay+1, nblocks))
  call stop_on_err(optical_props%alloc_2str(block_size, nlay, k_dist))
  ! --------------------------------------------------
#IFDEF USE_TIMING
  !
  ! Initialize timers
  !
  ret = gptlsetoption (gptlpercent, 1)        ! Turn on "% of" print
  ret = gptlsetoption (gptloverhead, 0)       ! Turn off overhead estimate
  ret =  gptlinitialize()
#ENDIF
  !
  ! Loop over blocks
  !
  do b = 1, nblocks
    fluxes%flux_up => flux_up(:,:,b)
    fluxes%flux_dn => flux_dn(:,:,b)
    !
    ! Compute the optical properties of the atmosphere and the Planck source functions
    !    from pressures, temperatures, and gas concentrations...
    !
#IFDEF USE_TIMING
    ret =  gptlstart('gas_optics (SW)')
#ENDIF
    call stop_on_err(k_dist%gas_optics(p_lay(:,:,b), &
                                       p_lev(:,:,b),       &
                                       t_lay(:,:,b),       &
                                       gas_conc_array(b),  &
                                       optical_props,      &
                                       toa_flux))
#IFDEF USE_TIMING
    ret =  gptlstop('gas_optics (SW)')
#ENDIF
    !
    ! Normalize incoming solar flux to match RFMIP specification
    !
    def_tsi(1:block_size) = sum(toa_flux, dim=2)
    do igpt = 1, k_dist%get_ngpt()
      do icol = 1, block_size
        toa_flux(icol,igpt) = toa_flux(icol,igpt) * total_solar_irradiance(icol,b)/def_tsi(icol)
      end do
    end do
    !
    ! RTE will fail if passed solar zenith angles greater than 90 degree. We replace any with
    !   nighttime columns with a default solar zenith angle. We'll mask these out later, of
    !   course, but this gives us more work and so a better measure of timing.
    !
    usecol(1:block_size)  = solar_zenith_angle(1:block_size,b) < 90._wp - 2._wp * spacing(90._wp)
    mu0(1:block_size) = merge(cos(solar_zenith_angle(:,b) * acos(-1._wp)/180._wp), 1._wp, usecol)
    !
    ! ... and compute the spectrally-resolved fluxes, providing reduced values
    !    via ty_fluxes_broadband
    !
#IFDEF USE_TIMING
    ret =  gptlstart('rte_sw')
#ENDIF
    call stop_on_err(rte_sw(optical_props,   &
                            top_at_1,        &
                            mu0,             &
                            toa_flux,        &
                            spread(surface_albedo(:,b), 1, ncopies = k_dist%get_nband()), &
                            spread(surface_albedo(:,b), 1, ncopies = k_dist%get_nband()), &
                            fluxes))
#IFDEF USE_TIMING
    ret =  gptlstop('rte_sw')
#ENDIF
    !
    ! Zero out fluxes for which the original solar zenith angle is > 90 degrees.
    !
    do icol = 1, block_size
      if(.not. usecol(icol)) then
        flux_up(icol,:,b)  = 0._wp
        flux_dn(icol,:,b)  = 0._wp
      end if
    end do
  end do
  !
  ! End timers
  !
#IFDEF USE_TIMING
  ret = gptlpr(block_size)
  ret = gptlfinalize()
#ENDIF
  ! --------------------------------------------------
  call unblock_and_write(trim(flxup_file), 'rsu', flux_up)
  call unblock_and_write(trim(flxdn_file), 'rsd', flux_dn)
end program rrtmgp_rfmip_sw
