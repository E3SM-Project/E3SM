module radiation_utils

   use shr_kind_mod, only: r8=>shr_kind_r8, cl=>shr_kind_cl
   use assertions, only: assert

   implicit none
   private

   public :: compress_day_columns, expand_day_columns, &
             calculate_heating_rate, clip_values, &
             handle_error, &
             fluxes_t, initialize_fluxes, reset_fluxes, free_fluxes, &
             expand_day_fluxes, get_gas_vmr

   ! Interface blocks for overloaded procedures
   interface compress_day_columns
      module procedure compress_day_columns_1d, compress_day_columns_2d
   end interface compress_day_columns

   interface expand_day_columns
      module procedure expand_day_columns_1d, expand_day_columns_2d
   end interface expand_day_columns

   interface clip_values
      module procedure clip_values_1d, clip_values_2d, clip_values_3d
   end interface clip_values

   ! Type to hold fluxes
   type fluxes_t
      real(r8), allocatable :: flux_up(:,:)
      real(r8), allocatable :: flux_dn(:,:)
      real(r8), allocatable :: flux_net(:,:)
      real(r8), allocatable :: flux_dn_dir(:,:)
      real(r8), allocatable :: bnd_flux_up(:,:,:)
      real(r8), allocatable :: bnd_flux_dn(:,:,:)
      real(r8), allocatable :: bnd_flux_net(:,:,:)
      real(r8), allocatable :: bnd_flux_dn_dir(:,:,:)
   end type


   ! Max length for character strings
   integer, parameter :: max_char_len = 512

   ! Name of this module for error messages
   character(len=*), parameter :: module_name = 'radiation_utils'

contains
   !-------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------
   subroutine initialize_fluxes(ncol, nlevels, nbands, fluxes, do_direct)

      integer, intent(in) :: ncol, nlevels, nbands
      type(fluxes_t), intent(inout) :: fluxes
      logical, intent(in), optional :: do_direct

      logical :: do_direct_local

      if (present(do_direct)) then
         do_direct_local = .true.
      else
         do_direct_local = .false.
      end if

      ! Allocate flux arrays
      ! NOTE: fluxes defined at interfaces, so need to either pass nlevels as
      ! number of model levels plus one, or allocate as nlevels+1 if nlevels
      ! represents number of model levels rather than number of interface levels.

      ! Broadband fluxes
      allocate(fluxes%flux_up(ncol, nlevels))
      allocate(fluxes%flux_dn(ncol, nlevels))
      allocate(fluxes%flux_net(ncol, nlevels))
      if (do_direct_local) allocate(fluxes%flux_dn_dir(ncol, nlevels))

      ! Fluxes by band
      allocate(fluxes%bnd_flux_up(ncol, nlevels, nbands))
      allocate(fluxes%bnd_flux_dn(ncol, nlevels, nbands))
      allocate(fluxes%bnd_flux_net(ncol, nlevels, nbands))
      if (do_direct_local) allocate(fluxes%bnd_flux_dn_dir(ncol, nlevels, nbands))

      ! Initialize
      call reset_fluxes(fluxes)

   end subroutine initialize_fluxes
   !-------------------------------------------------------------------------------
   subroutine reset_fluxes(fluxes)

      type(fluxes_t), intent(inout) :: fluxes

      ! Reset broadband fluxes
      fluxes%flux_up(:,:) = 0._r8
      fluxes%flux_dn(:,:) = 0._r8
      fluxes%flux_net(:,:) = 0._r8
      if (allocated(fluxes%flux_dn_dir)) fluxes%flux_dn_dir(:,:) = 0._r8

      ! Reset band-by-band fluxes
      fluxes%bnd_flux_up(:,:,:) = 0._r8
      fluxes%bnd_flux_dn(:,:,:) = 0._r8
      fluxes%bnd_flux_net(:,:,:) = 0._r8
      if (allocated(fluxes%bnd_flux_dn_dir)) fluxes%bnd_flux_dn_dir(:,:,:) = 0._r8

   end subroutine reset_fluxes
   !-------------------------------------------------------------------------------
   subroutine free_fluxes(fluxes)
      type(fluxes_t), intent(inout) :: fluxes
      if (allocated(fluxes%flux_up)) deallocate(fluxes%flux_up)
      if (allocated(fluxes%flux_dn)) deallocate(fluxes%flux_dn)
      if (allocated(fluxes%flux_net)) deallocate(fluxes%flux_net)
      if (allocated(fluxes%flux_dn_dir)) deallocate(fluxes%flux_dn_dir)
      if (allocated(fluxes%bnd_flux_up)) deallocate(fluxes%bnd_flux_up)
      if (allocated(fluxes%bnd_flux_dn)) deallocate(fluxes%bnd_flux_dn)
      if (allocated(fluxes%bnd_flux_net)) deallocate(fluxes%bnd_flux_net)
      if (allocated(fluxes%bnd_flux_dn_dir)) deallocate(fluxes%bnd_flux_dn_dir)
   end subroutine free_fluxes
   !-------------------------------------------------------------------------------
   subroutine expand_day_fluxes(daytime_fluxes, expanded_fluxes, day_indices)
      type(fluxes_t), intent(in) :: daytime_fluxes
      type(fluxes_t), intent(inout) :: expanded_fluxes
      integer, intent(in) :: day_indices(:)
      integer :: nday, iday, icol

      ! Reset fluxes in expanded_fluxes object to zero
      call reset_fluxes(expanded_fluxes)

      ! Number of daytime columns is number of indices greater than zero
      nday = count(day_indices > 0)

      ! Loop over daytime indices and map daytime fluxes into expanded arrays
      do iday = 1,nday

         ! Map daytime index to proper column index
         icol = day_indices(iday)

         ! Expand broadband fluxes
         expanded_fluxes%flux_up(icol,:) = daytime_fluxes%flux_up(iday,:)
         expanded_fluxes%flux_dn(icol,:) = daytime_fluxes%flux_dn(iday,:)
         expanded_fluxes%flux_net(icol,:) = daytime_fluxes%flux_net(iday,:)
         if (allocated(daytime_fluxes%flux_dn_dir)) then
            expanded_fluxes%flux_dn_dir(icol,:) = daytime_fluxes%flux_dn_dir(iday,:)
         end if

         ! Expand band-by-band fluxes
         expanded_fluxes%bnd_flux_up(icol,:,:) = daytime_fluxes%bnd_flux_up(iday,:,:)
         expanded_fluxes%bnd_flux_dn(icol,:,:) = daytime_fluxes%bnd_flux_dn(iday,:,:)
         expanded_fluxes%bnd_flux_net(icol,:,:) = daytime_fluxes%bnd_flux_net(iday,:,:)
         if (allocated(daytime_fluxes%bnd_flux_dn_dir)) then
            expanded_fluxes%bnd_flux_dn_dir(icol,:,:) = daytime_fluxes%bnd_flux_dn_dir(iday,:,:)
         end if

      end do

   end subroutine expand_day_fluxes
   !-------------------------------------------------------------------------------
   subroutine compress_day_columns_1d(xcol, xday, day_indices)

      real(r8), intent(in) :: xcol(:)
      real(r8), intent(inout) :: xday(:)
      integer, intent(in) :: day_indices(:)
      integer :: icol, iday
      character(len=32) :: subname = 'compress_day_columns'

      ! Check dimensions
      call assert(size(xday, 1) == size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xcol, 1) >= size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      
      ! Do mapping
      do iday = 1,count(day_indices>0)
         icol = day_indices(iday)
         xday(iday) = xcol(icol)
      end do

   end subroutine
   !-------------------------------------------------------------------------------
   subroutine compress_day_columns_2d(xcol, xday, day_indices)

      real(r8), intent(in) :: xcol(:,:)
      real(r8), intent(inout) :: xday(:,:)
      integer, intent(in) :: day_indices(:)
      integer :: icol, iday
      character(len=32) :: subname = 'compress_day_columns'

      ! Check dimensions
      call assert(size(xday, 1) == size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xcol, 1) >= size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xday, 2) == size(xcol, 2), trim(subname) // ': inconsistent sizes')

      ! Do mapping
      do iday = 1,count(day_indices>0)
         icol = day_indices(iday)
         xday(iday,:) = xcol(icol,:)
      end do

   end subroutine
   !-------------------------------------------------------------------------------
   subroutine expand_day_columns_1d(xday, xcol, day_indices)

      real(r8), intent(in) :: xday(:)
      real(r8), intent(inout) :: xcol(:)
      integer, intent(in) :: day_indices(:)
      integer :: icol, iday
      character(len=32) :: subname = 'expand_day_columns_1d'

      ! Check dimensions
      call assert(size(xday, 1) == size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xcol, 1) >= size(day_indices, 1), trim(subname) // ': inconsistent sizes')

      ! We need to reset to zero because we only populate the daytime columns
      xcol(:) = 0._r8

      ! Do mapping
      do iday = 1,count(day_indices>0)
         icol = day_indices(iday)
         xcol(icol) = xday(iday)
      end do
   end subroutine expand_day_columns_1d
   !-------------------------------------------------------------------------------
   subroutine expand_day_columns_2d(xday, xcol, day_indices)

      real(r8), intent(in) :: xday(:,:)
      real(r8), intent(inout) :: xcol(:,:)
      integer, intent(in) :: day_indices(:)
      integer :: icol, iday
      character(len=32) :: subname = 'expand_day_columns_2d'

      ! Check dimensions
      call assert(size(xday, 1) == size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xcol, 1) >= size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xday, 2) == size(xcol, 2), trim(subname) // ': inconsistent sizes')

      ! We need to reset to zero because we only populate the daytime columns
      xcol(:,:) = 0._r8

      ! Do mapping
      do iday = 1,count(day_indices>0)
         icol = day_indices(iday)
         xcol(icol,:) = xday(iday,:)
      end do

   end subroutine expand_day_columns_2d

   !-------------------------------------------------------------------------------

   subroutine get_gas_vmr(icall, state, pbuf, gas_names, gas_vmr) 

      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc
      use rad_constituents, only: rad_cnst_get_gas

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      character(len=*), intent(in), dimension(:) :: gas_names
      real(r8), intent(out), dimension(:,:,:) :: gas_vmr

      ! Mass mixing ratio
      real(r8), pointer :: mmr(:,:)

      ! Gases and molecular weights. Note that we do NOT have CFCs yet (I think
      ! this is coming soon in RRTMGP). RRTMGP also allows for absorption due to
      ! CO and N2, which RRTMG did not have.
      character(len=3), dimension(8) :: gas_species = (/ &
         'H2O', 'CO2', 'O3 ', 'N2O', &
         'CO ', 'CH4', 'O2 ', 'N2 ' &
      /)
      real(r8), dimension(8) :: mol_weight_gas = (/ &
         18.01528, 44.0095, 47.9982, 44.0128, &
         28.0101, 16.04246, 31.998, 28.0134 &
      /)  ! g/mol

      ! Molar weight of air
      real(r8), parameter :: mol_weight_air = 28.97  ! g/mol
                                       
      ! Defaults for gases that are not available (TODO: is this still accurate?)
      real(r8), parameter :: co_vol_mix_ratio = 1.0e-7_r8
      real(r8), parameter :: n2_vol_mix_ratio = 0.7906_r8

      ! Loop indices
      integer :: igas

      ! Number of columns
      integer :: ncol
      integer :: nlev

      ! Name of subroutine for error messages
      character(len=32) :: subname = 'get_gas_vmr'

      ! Number of columns in chunk
      ncol = state%ncol
      nlev = size(gas_vmr,3)

      ! initialize
      gas_vmr(:,:,:) = 0._r8

      ! For each gas species needed for RRTMGP, read the mass mixing ratio from the
      ! CAM rad_constituents interface, convert to volume mixing ratios, and
      ! subset for daytime-only indices if needed.
      do igas = 1,size(gas_names)

         select case(trim(gas_names(igas)))

            case('CO')

               ! CO not available, use default
               gas_vmr(igas,1:ncol,1:nlev) = co_vol_mix_ratio

            case('N2')

               ! N2 not available, use default
               gas_vmr(igas,1:ncol,1:nlev) = n2_vol_mix_ratio

            case('H2O')

               ! Water vapor is represented as specific humidity in CAM, so we
               ! need to handle water a little differently
               call rad_cnst_get_gas(icall, trim(gas_species(igas)), state, pbuf, mmr)

               ! Convert to volume mixing ratio by multiplying by the ratio of
               ! molecular weight of dry air to molecular weight of gas. Note that
               ! first specific humidity (held in the mass_mix_ratio array read
               ! from rad_constituents) is converted to an actual mass mixing
               ! ratio.
               gas_vmr(igas,1:ncol,1:nlev) = mmr(1:ncol,1:nlev) / ( &
                  1._r8 - mmr(1:ncol,1:nlev) &
               )  * mol_weight_air / mol_weight_gas(igas)

            case DEFAULT

               ! Get mass mixing ratio from the rad_constituents interface
               call rad_cnst_get_gas(icall, trim(gas_species(igas)), state, pbuf, mmr)

               ! Convert to volume mixing ratio by multiplying by the ratio of
               ! molecular weight of dry air to molecular weight of gas
               gas_vmr(igas,1:ncol,1:nlev) = mmr(1:ncol,1:nlev) &
                                            * mol_weight_air / mol_weight_gas(igas)

         end select

      end do  ! igas

   end subroutine get_gas_vmr

   !-------------------------------------------------------------------------------

   subroutine calculate_heating_rate(flux_up, flux_dn, pint, heating_rate)

      use physconst, only: gravit

      ! Inputs
      real(r8), intent(in), dimension(:,:) :: flux_up, flux_dn, pint

      ! Output heating rate; same size as pint with one fewer vertical level
      real(r8), intent(out) :: heating_rate(:,:)

      ! Loop indices
      integer :: icol, ilev

      ! Everyone needs a name
      character(len=32) :: subname = 'calculate_heating_rate'

      ! Get dimension sizes and make sure arrays conform
      call assert(size(pint,1) == size(flux_up,1), subname // ': sizes do not conform.')
      call assert(size(pint,2) == size(flux_up,2), subname // ': sizes do not conform.')
      call assert(size(heating_rate,1) == size(flux_up,1), subname // ': sizes do not conform.')
      call assert(size(heating_rate,2) == size(flux_up,2)-1, subname // ': sizes do not conform.')

      ! Loop over levels and calculate heating rates; note that the fluxes *should*
      ! be defined at interfaces, so the loop ktop,kbot and grabbing the current
      ! and next value of k should be safe. ktop should be the top interface, and
      ! kbot + 1 should be the bottom interface.
      !
      ! NOTE: to get heating rate in K/day, normally we would use:
      !
      !     H = dF / dp * g * (sec/day) * (1e-5) / (cpair)
      !
      ! Here we just use
      !
      !     H = dF / dp * g
      !
      ! Why? Something to do with convenience with applying the fluxes to the
      ! heating tendency?
      do ilev = 1,size(pint,2)-1
         do icol = 1,size(pint,1)
            heating_rate(icol,ilev) = ( &
               flux_up(icol,ilev+1) - flux_up(icol,ilev) - &
               flux_dn(icol,ilev+1) + flux_dn(icol,ilev) &
            ) * gravit / (pint(icol,ilev+1) - pint(icol,ilev))
         end do
      end do

   end subroutine calculate_heating_rate

   !----------------------------------------------------------------------------

   ! Routines to clip array values if they are outside of an expected range,
   ! defined by min_x and max_x. Allow passing a varname argument, which will be
   ! appended to warning message and tolerance argument that sets a maximum
   ! tolerance above which a warning or error will be thrown. That is, these 
   ! routines will *always* clip values outside expected range and allow the 
   ! simulation to continue, but a warning will be issued if the values fall 
   ! outside the expected range plus the tolerance. In other words, no warning 
   ! is issued when clipping values outside the valid range plus/minus the 
   ! tolerance.
   function clip_values_1d(x, min_x, max_x, varname, tolerance) result(error_message)
      real(r8), intent(inout) :: x(:)
      real(r8), intent(in) :: min_x
      real(r8), intent(in) :: max_x
      character(len=*), intent(in) :: varname
      real(r8), intent(in), optional :: tolerance
      real(r8) :: tolerance_local
      character(max_char_len) :: error_message

      error_message = ''
      tolerance_local = 0._r8
      if (present(tolerance)) then
         tolerance_local = tolerance
      end if

      ! look for values less than threshold
      if (any(x < min_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x < (min_x - tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x < (min_x - tolerance_local)), ' values below threshold ', &
                     '; min = ', minval(x)
         end if
         ! Clip values
         where (x < min_x)
            x = min_x
         endwhere
      end if

      ! Look for values greater than threshold
      if (any(x > max_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x > (max_x + tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x > (max_x + tolerance_local)), ' values above threshold ', &
                     '; max = ', maxval(x)
         end if
         ! Clip values
         where (x > max_x)
            x = max_x
         end where
      end if
   end function clip_values_1d
   !-------------------------------------------------------------------------------
   function clip_values_2d(x, min_x, max_x, varname, tolerance) result(error_message)
      real(r8), intent(inout) :: x(:,:)
      real(r8), intent(in) :: min_x
      real(r8), intent(in) :: max_x
      character(len=*), intent(in) :: varname
      real(r8), intent(in), optional :: tolerance
      real(r8) :: tolerance_local
      character(max_char_len) :: error_message

      error_message = ''
      tolerance_local = 0._r8
      if (present(tolerance)) then
         tolerance_local = tolerance
      end if

      ! look for values less than threshold
      if (any(x < min_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x < (min_x - tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x < (min_x - tolerance_local)), ' values below threshold ', &
                     '; min = ', minval(x)
         end if
         ! Clip values
         where (x < min_x)
            x = min_x
         endwhere
      end if

      ! Look for values greater than threshold
      if (any(x > max_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x > (max_x + tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x > (max_x + tolerance_local)), ' values above threshold ', &
                     '; max = ', maxval(x)
         end if
         ! Clip values
         where (x > max_x)
            x = max_x
         end where
      end if
   end function clip_values_2d
   !-------------------------------------------------------------------------------
   function clip_values_3d(x, min_x, max_x, varname, tolerance) result(error_message)
      real(r8), intent(inout) :: x(:,:,:)
      real(r8), intent(in) :: min_x
      real(r8), intent(in) :: max_x
      character(len=*), intent(in) :: varname
      real(r8), intent(in), optional :: tolerance
      real(r8) :: tolerance_local
      character(max_char_len) :: error_message

      error_message = ''
      tolerance_local = 0._r8
      if (present(tolerance)) then
         tolerance_local = tolerance
      end if

      ! look for values less than threshold
      if (any(x < min_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x < (min_x - tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x < (min_x - tolerance_local)), ' values below threshold ', &
                     '; min = ', minval(x)
         end if
         ! Clip values
         where (x < min_x)
            x = min_x
         endwhere
      end if

      ! Look for values greater than threshold
      if (any(x > max_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x > (max_x + tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x > (max_x + tolerance_local)), ' values above threshold ', &
                     '; max = ', maxval(x)
         end if
         ! Clip values
         where (x > max_x)
            x = max_x
         end where
      end if
   end function clip_values_3d
   !-------------------------------------------------------------------------------
   subroutine handle_error(error_message, fatal, warn)
      use cam_abortutils, only: endrun
      character(len=*), intent(in) :: error_message
      logical, intent(in), optional :: fatal
      logical, intent(in), optional :: warn
      logical :: fatal_local
      logical :: warn_local

      ! Allow passing of an optional flag to not stop the run if an error is
      ! encountered. This allows this subroutine to be used when inquiring if a
      ! variable exists without failing.
      if (present(fatal)) then
         fatal_local = fatal
      else
         fatal_local = .true.
      end if

      ! Allow optional flag to disable warning messages.
      ! Useful for avoiding low temperature messages in aquaplanet cases.
      if (present(warn)) then
         warn_local = warn
      else
         warn_local = .true.
      end if

      ! If we encounter an error, fail if we require success. Otherwise do
      ! nothing and return silently.
      if (len(trim(error_message)) > 0) then
         if (fatal_local) then
            call endrun(trim(error_message))
         else if (warn_local) then
            print *, trim(error_message)
         end if
      end if
   end subroutine handle_error
   !----------------------------------------------------------------------------
end module radiation_utils
