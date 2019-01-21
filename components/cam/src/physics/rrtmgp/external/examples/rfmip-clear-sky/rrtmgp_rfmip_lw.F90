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
! Example program to demonstrate the calculation of longwave radiative fluxes in clear, aerosol-free skies.
!   The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
!   The large problem (1800 profiles) is divided into blocks
!
! Program is invoked as rrtmgp_rfmip_lw [block_size input_file  coefficient_file upflux_file downflux_file]
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
    write (error_unit,*) "rrtmgp_rfmip_lw stopping"
    stop
  end if
end subroutine stop_on_err
! -------------------------------------------------------------------------------------------------
!
! Main program
!
! -------------------------------------------------------------------------------------------------
program rrtmgp_rfmip_lw
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
  !   Shortwave calculations would use optical depth, single-scattering albedo, asymmetry parameter (_2str)
  !
  use mo_optical_props,      only: ty_optical_props_1scl
  !
  ! Gas optics: maps physical state of the atmosphere to optical properties
  !
  use mo_gas_optics,         only: ty_gas_optics
  !
  ! Gas optics uses a derived type to represent gas concentrations compactly...
  !
  use mo_gas_concentrations, only: ty_gas_concs
  !
  ! ... and another type to encapsulate the longwave source functions.
  !
  use mo_source_functions,   only: ty_source_func_lw
  !
  ! RTE longwave driver
  !
  use mo_rte_lw,             only: rte_lw
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
                                   read_and_block_lw_bc, read_kdist_gas_names
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
                                kdist_file = 'coefficients_lw.nc', &
                                flxdn_file = 'rld_template.nc', flxup_file = 'rlu_template.nc'
  integer                    :: nargs, ncol, nlay, nexp, nblocks, block_size
  logical                    :: top_at_1
  integer                    :: b
  character(len=6)           :: block_size_char

  character(len=32 ), &
            dimension(:),             allocatable :: kdist_gas_names, gases_to_use
  real(wp), dimension(:,:,:),         allocatable :: p_lay, p_lev, t_lay, t_lev ! block_size, nlay, nblocks
  real(wp), dimension(:,:,:), target, allocatable :: flux_up, flux_dn
  real(wp), dimension(:,:  ),         allocatable :: sfc_emis, sfc_t            ! block_size, nblocks

  !
  ! Classes used by rte+rrtmgp
  !
  type(ty_gas_optics)                            :: k_dist
  type(ty_source_func_lw)                        :: source
  type(ty_optical_props_1scl)                    :: optical_props
  type(ty_fluxes_broadband)                      :: fluxes
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
  if(mod(ncol*nexp, block_size) /= 0 ) call stop_on_err("rrtmgp_rfmip_lw: number of columns doesn't fit evenly into blocks.")
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
  call read_and_block_lw_bc(rfmip_file, block_size, sfc_emis, sfc_t)

  !
  ! Read k-distribution information. load_and_init() reads data from netCDF and calls
  !   k_dist%init(); users might want to use their own reading methods
  !
  call load_and_init(k_dist, trim(kdist_file), gas_conc_array(1))
  if(.not. k_dist%source_is_internal()) &
    stop "rrtmgp_rfmip_lw: k-distribution file isn't LW"

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
  allocate(flux_up(    block_size, nlay+1, nblocks), &
           flux_dn(    block_size, nlay+1, nblocks))
  call stop_on_err(source%alloc            (block_size, nlay, k_dist))
  call stop_on_err(optical_props%alloc_1scl(block_size, nlay, k_dist))
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
    ret =  gptlstart('gas_optics (LW)')
#ENDIF
    call stop_on_err(k_dist%gas_optics(p_lay(:,:,b), &
                                       p_lev(:,:,b),       &
                                       t_lay(:,:,b),       &
                                       sfc_t(:  ,b),       &
                                       gas_conc_array(b),  &
                                       optical_props,      &
                                       source,             &
                                       tlev = t_lev(:,:,b)))
#IFDEF USE_TIMING
    ret =  gptlstop('gas_optics (LW)')
#ENDIF
    !
    ! ... and compute the spectrally-resolved fluxes, providing reduced values
    !    via ty_fluxes_broadband
    !
#IFDEF USE_TIMING
    ret =  gptlstart('rte_lw')
#ENDIF
    call stop_on_err(rte_lw(optical_props,   &
                            top_at_1,        &
                            source,          &
                            spread(sfc_emis(:,b), 1, ncopies = k_dist%get_nband()), &
                            fluxes))
#IFDEF USE_TIMING
    ret =  gptlstop('rte_lw')
#ENDIF
  end do
#IFDEF USE_TIMING
  !
  ! End timers
  !
  ret = gptlpr(block_size)
  ret = gptlfinalize()
#ENDIF
  ! --------------------------------------------------
  call unblock_and_write(trim(flxup_file), 'rlu', flux_up)
  call unblock_and_write(trim(flxdn_file), 'rld', flux_dn)
end program rrtmgp_rfmip_lw
