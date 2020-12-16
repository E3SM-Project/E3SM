! Module to bridge the gap between the Fortran and C++ implemenations of
! RRTMGP. Remove class references from function calls, and handle all of that
! here. This is necessary because radiation_tend will remain in F90 (to deal
! with E3SM data types), but we will switch to C++ for the underlying RRTMGP
! code.
module rrtmgpxx_interface

   use iso_c_binding

   implicit none

   private

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
