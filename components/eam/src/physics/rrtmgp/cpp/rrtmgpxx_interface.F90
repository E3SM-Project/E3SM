! Module to bridge the gap between the Fortran and C++ implemenations of
! RRTMGP. Remove class references from function calls, and handle all of that
! here. This is necessary because radiation_tend will remain in F90 (to deal
! with E3SM data types), but we will switch to C++ for the underlying RRTMGP
! code.
module rrtmgpxx_interface

   use perf_mod, only: t_startf, t_stopf
   use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk
   use radiation_utils, only: compress_day_columns, expand_day_columns
   use radiation_state, only: ktop, kbot, nlev_rad

   use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
   use mo_gas_concentrations, only: ty_gas_concs
   use mo_load_coefficients, only: load_and_init
   use mo_rte_kind, only: wp
   use mo_optical_props, only: ty_optical_props_2str, ty_optical_props_1scl
   use mo_fluxes_byband, only: ty_fluxes_byband
   use mo_rrtmgp_clr_all_sky, only: rte_sw, rte_lw
   use assertions, only: assert
   use iso_c_binding

   implicit none

   private

   ! Gas optics objects that hold k-distribution information. These are made
   ! module variables because we only want to initialize them once at init time.
   type(ty_gas_optics_rrtmgp) :: k_dist_sw, k_dist_lw

   ! Make these module variables so that we do not have to provide access to
   ! k_dist objects; this just makes it easier to switch between F90 and C++
   ! interfaces.
   integer, public :: nswbands, nlwbands, nswgpts, nlwgpts

   public :: &
      rrtmgpxx_initialize, rrtmgpxx_finalize, &
      rrtmgpxx_run_sw, rrtmgpxx_run_lw, &
      get_nband_sw, get_nband_lw, &
      get_ngpt_sw, get_ngpt_lw, &
      get_gpoint_bands_sw, get_gpoint_bands_lw, &
      get_min_temperature, get_max_temperature
   
   interface 

      function get_nband_sw() bind(C,name="get_nband_sw")
         use iso_c_binding
         implicit none
         integer(c_int) :: get_nband_sw
      end function

      function get_nband_lw() bind(C,name="get_nband_lw")
         use iso_c_binding
         implicit none
         integer(c_int) :: get_nband_lw
      end function

      function get_ngpt_sw() bind(C, name="get_ngpt_sw")
         use iso_c_binding
         implicit none
         integer(c_int) :: get_ngpt_sw
      end function

      function get_ngpt_lw() bind(C, name="get_ngpt_lw")
         use iso_c_binding
         implicit none
         integer(c_int) :: get_ngpt_lw
      end function

      function get_min_temperature() bind(C, name="get_min_temperature")
         use iso_c_binding
         implicit none
         real(c_double) :: get_min_temperature
      end function

      function get_max_temperature() bind(C, name="get_max_temperature")
         use iso_c_binding
         implicit none
         real(c_double) :: get_max_temperature
      end function

      subroutine get_gpoint_bands_sw(gpoint_bands) bind(C, name="get_gpoint_bands_sw")
         use iso_c_binding
         implicit none
         integer(c_int), dimension(*) :: gpoint_bands
      end subroutine

      subroutine get_gpoint_bands_lw(gpoint_bands) bind(C, name="get_gpoint_bands_lw")
         use iso_c_binding
         implicit none
         integer(c_int), dimension(*) :: gpoint_bands
      end subroutine

      subroutine rrtmgpxx_initialize_cpp(coefficients_file_sw, coefficients_file_lw) bind(C, name="rrtmgpxx_initialize_cpp")
         use iso_c_binding, only: C_CHAR, C_NULL_CHAR
         implicit none
         character(kind=c_char) :: coefficients_file_sw(*)
         character(kind=c_char) :: coefficients_file_lw(*)
      end subroutine rrtmgpxx_initialize_cpp

      subroutine rrtmgpxx_finalize() bind(C, name="rrtmgpxx_finalize")
      end subroutine rrtmgpxx_finalize

      subroutine rrtmgpxx_run_sw( &
         ngas, ncol, nlev, &
         gas_vmr, &
         pmid, tmid, pint, coszrs, &
         albedo_dir, albedo_dif, &
         cld_tau_gpt, cld_ssa_gpt, cld_asm_gpt, &
         aer_tau_bnd, aer_ssa_bnd, aer_asm_bnd, &
         allsky_flux_up, allsky_flux_dn, allsky_flux_net, allsky_flux_dn_dir, &
         allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, allsky_bnd_flux_dn_dir, &
         clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net, clrsky_flux_dn_dir, &
         clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, clrsky_bnd_flux_dn_dir, &
         tsi_scaling &
         ) bind(C, name="rrtmgpxx_run_sw")
         use iso_c_binding
         implicit none
         integer(kind=c_int), value :: ngas, ncol, nlev
         real(kind=c_double), dimension(*) :: &
            gas_vmr, pmid, tmid, pint, coszrs, albedo_dir, albedo_dif, &
            cld_tau_gpt, cld_ssa_gpt, cld_asm_gpt, &
            aer_tau_bnd, aer_ssa_bnd, aer_asm_bnd, &
            allsky_flux_up, allsky_flux_dn, allsky_flux_net, allsky_flux_dn_dir, &
            allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, allsky_bnd_flux_dn_dir, &
            clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net, clrsky_flux_dn_dir, &
            clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, clrsky_bnd_flux_dn_dir
         real(kind=c_double), value :: tsi_scaling
      end subroutine rrtmgpxx_run_sw

      subroutine rrtmgpxx_run_lw ( &
         ngas, ncol, nlev, &
         gas_vmr, &
         pmid, tmid, pint, tint, &
         surface_emissivity, &
         cld_tau, aer_tau, &
         allsky_flux_up_cxx    , allsky_flux_dn_cxx    , allsky_flux_net_cxx, &
         allsky_bnd_flux_up_cxx, allsky_bnd_flux_dn_cxx, allsky_bnd_flux_net_cxx, &
         clrsky_flux_up_cxx    , clrsky_flux_dn_cxx    , clrsky_flux_net_cxx, &
         clrsky_bnd_flux_up_cxx, clrsky_bnd_flux_dn_cxx, clrsky_bnd_flux_net_cxx &
         ) bind(C, name="rrtmgpxx_run_lw")
         use iso_c_binding
         implicit none
         integer(kind=c_int), value :: ngas, ncol, nlev
         real(kind=c_double), dimension(*) :: &
            gas_vmr, &
            pmid, tmid, pint, tint, surface_emissivity, &
            cld_tau, aer_tau, &
            allsky_flux_up_cxx, allsky_flux_dn_cxx, allsky_flux_net_cxx, &
            allsky_bnd_flux_up_cxx, allsky_bnd_flux_dn_cxx, allsky_bnd_flux_net_cxx, &
            clrsky_flux_up_cxx, clrsky_flux_dn_cxx, clrsky_flux_net_cxx, &
            clrsky_bnd_flux_up_cxx, clrsky_bnd_flux_dn_cxx, clrsky_bnd_flux_net_cxx
      end subroutine rrtmgpxx_run_lw

      subroutine add_gas_name(gas_name) bind(C, name="add_gas_name")
         use iso_c_binding, only: C_CHAR
         character(kind=c_char) :: gas_name
      end subroutine add_gas_name

   end interface

contains

   subroutine rrtmgpxx_initialize(active_gases, coefficients_file_sw, coefficients_file_lw)
      use iso_c_binding, only: C_CHAR, C_NULL_CHAR
      character(len=*), intent(in) :: active_gases(:)
      character(len=*), intent(in) :: coefficients_file_sw, coefficients_file_lw
      type(ty_gas_concs) :: available_gases
      ! Read gas optics coefficients from file
      ! Need to initialize available_gases here! The only field of the
      ! available_gases type that is used int he kdist initialize is
      ! available_gases%gas_name, which gives the name of each gas that would be
      ! present in the ty_gas_concs object. So, we can just set this here, rather
      ! than trying to fully populate the ty_gas_concs object here, which would be
      ! impossible from this initialization routine because I do not thing the
      ! rad_cnst objects are setup yet.
      ! the other tasks!
      ! TODO: This needs to be fixed to ONLY read in the data if masterproc, and then broadcast
      call set_available_gases(active_gases, available_gases)
      call load_and_init(k_dist_sw, coefficients_file_sw, available_gases)
      call load_and_init(k_dist_lw, coefficients_file_lw, available_gases)
      ! Add active gases
      call add_gases(active_gases)
      ! Initialize RRTMGP
      call rrtmgpxx_initialize_cpp( &
         C_CHAR_""//trim(coefficients_file_sw)//C_NULL_CHAR, &
         C_CHAR_""//trim(coefficients_file_lw)//C_NULL_CHAR &
      )
      ! Set number of bands based on what we read in from input data
      nswbands = get_nband_sw()
      nlwbands = get_nband_lw()
      ! Number of gpoints depend on inputdata, so initialize here
      nswgpts = get_ngpt_sw()
      nlwgpts = get_ngpt_lw()
   end subroutine rrtmgpxx_initialize

   ! --------------------------------------------------------------------------
   ! Private routines
   ! --------------------------------------------------------------------------

   subroutine set_available_gases(gases, gas_concentrations)

      use mo_gas_concentrations, only: ty_gas_concs
      use mo_rrtmgp_util_string, only: lower_case

      type(ty_gas_concs), intent(inout) :: gas_concentrations
      character(len=*), intent(in) :: gases(:)
      character(len=32), dimension(size(gases)) :: gases_lowercase
      integer :: igas

      ! Initialize with lowercase gas names; we should work in lowercase
      ! whenever possible because we cannot trust string comparisons in RRTMGP
      ! to be case insensitive
      do igas = 1,size(gases)
         gases_lowercase(igas) = trim(lower_case(gases(igas)))
      end do
      call handle_error(gas_concentrations%init(gases_lowercase))

   end subroutine set_available_gases

   subroutine add_gases(gases)
      use mo_rrtmgp_util_string, only: lower_case
      use iso_c_binding, only: C_CHAR, C_NULL_CHAR
      character(len=*), intent(in) :: gases(:)
      integer :: igas
      do igas = 1,size(gases)
         call add_gas_name(trim(lower_case(gases(igas)))//C_NULL_CHAR)
      end do
   end subroutine add_gases

   !----------------------------------------------------------------------------

   subroutine free_optics_sw(optics)
      use mo_optical_props, only: ty_optical_props_2str
      type(ty_optical_props_2str), intent(inout) :: optics
      if (allocated(optics%tau)) deallocate(optics%tau)
      if (allocated(optics%ssa)) deallocate(optics%ssa)
      if (allocated(optics%g)) deallocate(optics%g)
      call optics%finalize()
   end subroutine free_optics_sw

   !----------------------------------------------------------------------------

   subroutine free_optics_lw(optics)
      use mo_optical_props, only: ty_optical_props_1scl
      type(ty_optical_props_1scl), intent(inout) :: optics
      if (allocated(optics%tau)) deallocate(optics%tau)
      call optics%finalize()
   end subroutine free_optics_lw

   !----------------------------------------------------------------------------

   ! Compress optics arrays to smaller arrays containing only daytime columns.
   ! This is to work with the RRTMGP shortwave routines that will fail if they
   ! encounter non-sunlit columns, and also allows us to perform less
   ! computations. This routine is primarily a convenience routine to do all of
   ! the shortwave optics arrays at once, as we do this for individual arrays
   ! elsewhere in the code.
   subroutine compress_optics_sw(day_indices, tau, ssa, asm, tau_day, ssa_day, asm_day)
      integer, intent(in), dimension(:) :: day_indices
      real(wp), intent(in), dimension(:,:,:) :: tau, ssa, asm
      real(wp), intent(out), dimension(:,:,:) :: tau_day, ssa_day, asm_day
      integer :: nday, iday, ilev, ibnd
      nday = count(day_indices > 0)
      do ibnd = 1,size(tau,3)
         do ilev = 1,size(tau,2)
            do iday = 1,nday
               tau_day(iday,ilev,ibnd) = tau(day_indices(iday),ilev,ibnd)
               ssa_day(iday,ilev,ibnd) = ssa(day_indices(iday),ilev,ibnd)
               asm_day(iday,ilev,ibnd) = asm(day_indices(iday),ilev,ibnd)
            end do
         end do
      end do
   end subroutine compress_optics_sw


   subroutine set_gas_concentrations(ncol, gas_names, gas_vmr, gas_concentrations)
      use mo_gas_concentrations, only: ty_gas_concs
      use mo_rrtmgp_util_string, only: lower_case

      integer, intent(in) :: ncol
      character(len=*), intent(in), dimension(:) :: gas_names
      real(wp), intent(in), dimension(:,:,:) :: gas_vmr
      type(ty_gas_concs), intent(out) :: gas_concentrations

      ! Local variables
      real(wp), dimension(ncol,nlev_rad) :: vol_mix_ratio_out

      ! Loop indices
      integer :: igas

      ! Character array to hold lowercase gas names
      character(len=32), allocatable :: gas_names_lower(:)

      ! Name of subroutine for error messages
      character(len=32) :: subname = 'set_gas_concentrations'

      ! Initialize gas concentrations with lower case names
      allocate(gas_names_lower(size(gas_names)))
      do igas = 1,size(gas_names)
         gas_names_lower(igas) = trim(lower_case(gas_names(igas)))
      end do
      call handle_error(gas_concentrations%init(gas_names_lower))

      ! For each gas, add level above model top and set values in RRTMGP object
      do igas = 1,size(gas_names)
         vol_mix_ratio_out = 0
         ! Map to radiation grid
         vol_mix_ratio_out(1:ncol,ktop:kbot) = gas_vmr(igas,1:ncol,1:pver)
         ! Copy top-most model level to top-most rad level (which could be above
         ! the top of the model)
         vol_mix_ratio_out(1:ncol,1) = gas_vmr(igas,1:ncol,1)
         ! Set volumn mixing ratio in gas concentration object for just columns
         ! in this chunk
         call handle_error(gas_concentrations%set_vmr( &
            trim(lower_case(gas_names(igas))), vol_mix_ratio_out(1:ncol,1:nlev_rad)) &
         )
      end do

   end subroutine set_gas_concentrations

   !----------------------------------------------------------------------------

   ! Stop run ungracefully since we don't want dependencies on E3SM abortutils
   ! here
   subroutine handle_error(msg)
      character(len=*), intent(in) :: msg
      if (trim(msg) .ne. '') then
         print *, trim(msg)
         stop
      end if
   end subroutine handle_error

end module rrtmgpxx_interface
