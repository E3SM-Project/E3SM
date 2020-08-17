! Module to bridge the gap between the Fortran and C++ implemenations of
! RRTMGP. Remove class references from function calls, and handle all of that
! here. This is necessary because radiation_tend will remain in F90 (to deal
! with E3SM data types), but we will switch to C++ for the underlying RRTMGP
! code.
module rrtmgp_interface

   use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
   use mo_gas_concentrations, only: ty_gas_concs
   use mo_load_coefficients, only: load_and_init
   use mo_rte_kind, only: wp

   implicit none

   private

   ! Gas optics objects that hold k-distribution information. These are made
   ! module variables because we only want to initialize them once at init time.
   type(ty_gas_optics_rrtmgp), public :: k_dist_sw, k_dist_lw

   ! Make these module variables so that we do not have to provide access to
   ! k_dist objects; this just makes it easier to switch between F90 and C++
   ! interfaces.
   integer, public :: nswbands, nlwbands, nswgpts, nlwgpts

   public :: rrtmgp_initialize, &
      get_nbnds_sw, get_nbnds_lw, &
      get_ngpts_sw, get_ngpts_lw, &
      get_min_temperature, get_max_temperature

contains

   integer function get_nbnds_sw()
      get_nbnds_sw = k_dist_sw%get_nband()
   end function get_nbnds_sw

   integer function get_nbnds_lw()
      get_nbnds_lw = k_dist_lw%get_nband()
   end function get_nbnds_lw

   integer function get_ngpts_sw()
      get_ngpts_sw = k_dist_sw%get_ngpt()
   end function get_ngpts_sw

   integer function get_ngpts_lw()
      get_ngpts_lw = k_dist_lw%get_ngpt()
   end function get_ngpts_lw

   subroutine rrtmgp_initialize(active_gases, coefficients_file_sw, coefficients_file_lw)
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
      ! Set number of bands based on what we read in from input data
      nswbands = k_dist_sw%get_nband()
      nlwbands = k_dist_lw%get_nband()
      ! Number of gpoints depend on inputdata, so initialize here
      nswgpts = k_dist_sw%get_ngpt()
      nlwgpts = k_dist_lw%get_ngpt()
   end subroutine rrtmgp_initialize

   real(wp) function get_min_temperature()
      get_min_temperature = min(k_dist_sw%get_temp_min(), k_dist_lw%get_temp_min())
   end function get_min_temperature
     
   real(wp) function get_max_temperature()
      get_max_temperature = max(k_dist_sw%get_temp_max(), k_dist_lw%get_temp_max())
   end function get_max_temperature

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

   ! Stop run ungracefully since we don't want dependencies on E3SM abortutils
   ! here
   subroutine handle_error(msg)
      character(len=*), intent(in) :: msg
      if (trim(msg) .ne. '') then
         print *, trim(msg)
         stop
      end if
   end subroutine handle_error

end module rrtmgp_interface
