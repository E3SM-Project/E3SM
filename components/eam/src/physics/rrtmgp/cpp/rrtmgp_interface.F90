! Module to bridge the gap between the Fortran and C++ implemenations of
! RRTMGP. Remove class references from function calls, and handle all of that
! here. This is necessary because radiation_tend will remain in F90 (to deal
! with E3SM data types), but we will switch to C++ for the underlying RRTMGP
! code.
module rrtmgp_interface

   use iso_c_binding

   implicit none

   private

   ! Make these module variables so that we do not have to provide access to
   ! k_dist objects; this just makes it easier to switch between F90 and C++
   ! interfaces.
   integer, public :: nswbands, nlwbands, nswgpts, nlwgpts

   public :: &
      rrtmgp_initialize, rrtmgp_finalize, &
      rrtmgp_run_sw, rrtmgp_run_lw, &
      get_nband_sw, get_nband_lw, &
      get_ngpt_sw, get_ngpt_lw, &
      get_gpoint_bands_sw, get_gpoint_bands_lw, &
      get_min_temperature, get_max_temperature, &
      c_strarr
   
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

      subroutine rrtmgp_initialize_cxx(ngas, gas_names, coefficients_file_sw, coefficients_file_lw) bind(C, name="rrtmgp_initialize_cxx")
         use iso_c_binding
         implicit none
         integer(kind=c_int), value :: ngas
         type(c_ptr), dimension(*) :: gas_names
         character(kind=c_char) :: coefficients_file_sw(*)
         character(kind=c_char) :: coefficients_file_lw(*)
      end subroutine rrtmgp_initialize_cxx

      subroutine rrtmgp_finalize() bind(C, name="rrtmgp_finalize")
      end subroutine rrtmgp_finalize

      subroutine rrtmgp_run_sw( &
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
         ) bind(C, name="rrtmgp_run_sw")
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
      end subroutine rrtmgp_run_sw

      subroutine rrtmgp_run_lw ( &
         ngas, ncol, nlev, &
         gas_vmr, &
         pmid, tmid, pint, tint, &
         surface_emissivity, &
         cld_tau, aer_tau, &
         allsky_flux_up    , allsky_flux_dn    , allsky_flux_net, &
         allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, &
         clrsky_flux_up    , clrsky_flux_dn    , clrsky_flux_net, &
         clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net &
         ) bind(C, name="rrtmgp_run_lw")
         use iso_c_binding
         implicit none
         integer(kind=c_int), value :: ngas, ncol, nlev
         real(kind=c_double), dimension(*) :: &
            gas_vmr, &
            pmid, tmid, pint, tint, surface_emissivity, &
            cld_tau, aer_tau, &
            allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
            allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, &
            clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net, &
            clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net
      end subroutine rrtmgp_run_lw

   end interface

contains

   ! Need to wrap rrtmgp_initialize_cxx to make sure we pass c-style strings
   subroutine rrtmgp_initialize(ngas, gas_names, coefficients_file_sw, coefficients_file_lw)
      implicit none
      integer :: ngas
      character(len=*), dimension(:) :: gas_names
      character(len=*) :: coefficients_file_sw
      character(len=*) :: coefficients_file_lw
      ! Null-terminated C-compatible version gas names
      character(len=len(gas_names)+1), dimension(size(gas_names)), target :: gas_names_c
      call rrtmgp_initialize_cxx( &
         ngas, c_strarr(gas_names, gas_names_c), &
         trim(coefficients_file_sw)//C_NULL_CHAR, &
         trim(coefficients_file_lw)//C_NULL_CHAR &
      )
      ! Set number of gpoints
      nswgpts = get_ngpt_sw()
      nlwgpts = get_ngpt_lw()
   end subroutine rrtmgp_initialize

   ! Utility function to convert F90 string arrays to C-compatible string
   ! pointers; NOTE: str_c seems to need to be intent(out), or else the first
   ! element in the pointer array is messed up for some reason.
   function c_strarr(str, str_c) result(str_p)
      use iso_c_binding
      implicit none
      character(len=*), dimension(:), intent(in) :: str
      character(len=*), dimension(:), target, intent(out) :: str_c
      type(c_ptr), dimension(size(str)) :: str_p
      integer :: istr
      do istr = 1,size(str)
         str_c(istr) = trim(str(istr))//C_NULL_CHAR
         str_p(istr) = c_loc(str_c(istr))
      end do
   end function c_strarr

end module rrtmgp_interface
