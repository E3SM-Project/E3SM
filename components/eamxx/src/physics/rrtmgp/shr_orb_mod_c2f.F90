module shr_orb_mod_c2f

   use iso_c_binding
   use shr_orb_mod, only: shr_orb_params, shr_orb_decl, shr_orb_decl_gregorian, shr_orb_cosz, SHR_ORB_UNDEF_INT
   use shr_orb_mod, only: shr_orb_decl_fixed_year
   implicit none
   public :: shr_orb_params_c2f, shr_orb_decl_c2f, shr_orb_cosz_c2f
   public :: shr_orb_params_continuous_c2f
   public :: shr_orb_decl_gregorian_c2f
   integer(c_int), bind(C) :: shr_orb_undef_int_c2f = SHR_ORB_UNDEF_INT

   abstract interface
      logical(c_bool) function solar_decl_callback(ymd, tod, delta, eccf) bind(C)
         import :: c_bool, c_int, c_double
         integer(c_int), value :: ymd
         real(c_double), value :: tod
         real(c_double), intent(out) :: delta, eccf
      end function solar_decl_callback
      logical(c_bool) function solar_geometry_callback(ymd, tod, delta, eccf, eot_seconds) bind(C)
         import :: c_bool, c_int, c_double
         integer(c_int), value :: ymd
         real(c_double), value :: tod
         real(c_double), intent(out) :: delta, eccf, eot_seconds
      end function solar_geometry_callback
   end interface
   type(c_funptr), private, save :: coupled_solar_callback = c_null_funptr
   type(c_funptr), private, save :: coupled_geometry_callback = c_null_funptr
   logical(c_bool), private, save :: coupled_use_equation_of_time = .false._c_bool

contains

   subroutine shr_orb_set_geometry_callback_c2f(callback, enabled) bind(C)
      type(c_funptr), value :: callback
      logical(c_bool), value :: enabled
      coupled_geometry_callback = callback
      coupled_use_equation_of_time = enabled .and. c_associated(callback)
   end subroutine shr_orb_set_geometry_callback_c2f

   logical(c_bool) function shr_orb_has_geometry_callback_c2f() bind(C)
      shr_orb_has_geometry_callback_c2f = c_associated(coupled_geometry_callback)
   end function shr_orb_has_geometry_callback_c2f

   logical(c_bool) function shr_orb_use_equation_of_time_c2f() bind(C)
      shr_orb_use_equation_of_time_c2f = coupled_use_equation_of_time
   end function shr_orb_use_equation_of_time_c2f

   logical(c_bool) function shr_orb_coupled_geometry_c2f(ymd, tod, delta, eccf, eot_seconds) bind(C)
      integer(c_int), value :: ymd
      real(c_double), value :: tod
      real(c_double), intent(out) :: delta, eccf, eot_seconds
      procedure(solar_geometry_callback), pointer :: callback

      shr_orb_coupled_geometry_c2f = .false._c_bool
      if (.not. c_associated(coupled_geometry_callback)) return
      call c_f_procpointer(coupled_geometry_callback, callback)
      shr_orb_coupled_geometry_c2f = callback(ymd, tod, delta, eccf, eot_seconds)
   end function shr_orb_coupled_geometry_c2f

   subroutine shr_orb_set_decl_callback_c2f(callback) bind(C)
      type(c_funptr), value :: callback
      coupled_solar_callback = callback
   end subroutine shr_orb_set_decl_callback_c2f

   logical(c_bool) function shr_orb_coupled_decl_c2f(ymd, tod, delta, eccf) bind(C)
      integer(c_int), value :: ymd
      real(c_double), value :: tod
      real(c_double), intent(out) :: delta, eccf
      procedure(solar_decl_callback), pointer :: callback

      shr_orb_coupled_decl_c2f = .false._c_bool
      if (.not. c_associated(coupled_solar_callback)) return
      call c_f_procpointer(coupled_solar_callback, callback)
      shr_orb_coupled_decl_c2f = callback(ymd, tod, delta, eccf)
   end function shr_orb_coupled_decl_c2f

   subroutine shr_orb_params_c2f(           &
           orb_year, eccen, mvelpp, lambm0, &
           obliqr, delta, eccf            &
           ) bind(C, name='shr_orb_params_c2f')
      integer(c_int) :: orb_year
      real(c_double) :: eccen
      real(c_double) :: mvelpp
      real(c_double) :: lambm0
      real(c_double) :: obliqr
      real(c_double) :: delta
      real(c_double) :: eccf
      logical :: print_flag = .false.
      call shr_orb_params(orb_year, eccen, mvelpp, lambm0, obliqr, delta, eccf, print_flag)
   end subroutine shr_orb_params_c2f

   subroutine shr_orb_params_continuous_c2f(orb_year, year_fraction, eccen, obliq, mvelp, &
         obliqr, lambm0, mvelpp) bind(C, name='shr_orb_params_continuous_c2f')
      integer(c_int), value :: orb_year
      real(c_double), value :: year_fraction
      real(c_double), intent(inout) :: eccen, obliq, mvelp
      real(c_double), intent(out) :: obliqr, lambm0, mvelpp

      call shr_orb_params(orb_year, eccen, obliq, mvelp, obliqr, lambm0, mvelpp, &
           .false., year_fraction=year_fraction)
   end subroutine shr_orb_params_continuous_c2f

   subroutine shr_orb_decl_c2f(     &
         calday, eccen, mvelpp, lambm0, &
         obliqr, delta, eccf            &
         ) bind(C, name='shr_orb_decl_c2f')
      real(c_double), VALUE :: calday
      real(c_double), VALUE :: eccen
      real(c_double), VALUE :: mvelpp
      real(c_double), VALUE :: lambm0
      real(c_double), VALUE :: obliqr
      real(c_double) :: delta
      real(c_double) :: eccf
      call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
   end subroutine shr_orb_decl_c2f

   subroutine shr_orb_decl_fixed_year_c2f(calday, days_per_year, eccen, mvelpp, lambm0, obliqr, &
         delta, eccf, eot_seconds, use_eot) bind(C)
      real(c_double), value :: calday, days_per_year, eccen, mvelpp, lambm0, obliqr
      real(c_double), intent(out) :: delta, eccf, eot_seconds
      logical(c_bool), value :: use_eot

      eot_seconds = 0.0_c_double
      if (use_eot) then
         call shr_orb_decl_fixed_year(calday, days_per_year, eccen, mvelpp, lambm0, obliqr, &
              delta, eccf, eot_seconds)
      else
         call shr_orb_decl_fixed_year(calday, days_per_year, eccen, mvelpp, lambm0, obliqr, delta, eccf)
      end if
   end subroutine shr_orb_decl_fixed_year_c2f

   subroutine shr_orb_decl_gregorian_c2f(days_since_20000101, eccen, mvelpp, lambm0, &
         obliqr, delta, eccf) bind(C, name='shr_orb_decl_gregorian_c2f')
      real(c_double), value :: days_since_20000101, eccen, mvelpp, lambm0, obliqr
      real(c_double), intent(out) :: delta, eccf

      call shr_orb_decl_gregorian(days_since_20000101, eccen, mvelpp, lambm0, obliqr, delta, eccf)
   end subroutine shr_orb_decl_gregorian_c2f

   subroutine shr_orb_decl_eot_c2f(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf, eot_seconds) bind(C)
      real(c_double), value :: calday, eccen, mvelpp, lambm0, obliqr
      real(c_double), intent(out) :: delta, eccf, eot_seconds
      call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf, eot_seconds=eot_seconds)
   end subroutine shr_orb_decl_eot_c2f

   subroutine shr_orb_decl_gregorian_eot_c2f(elapsed_days, eccen, mvelpp, lambm0, &
         obliqr, delta, eccf, eot_seconds) bind(C)
      real(c_double), value :: elapsed_days, eccen, mvelpp, lambm0, obliqr
      real(c_double), intent(out) :: delta, eccf, eot_seconds
      call shr_orb_decl_gregorian(elapsed_days, eccen, mvelpp, lambm0, obliqr, delta, eccf, eot_seconds)
   end subroutine shr_orb_decl_gregorian_eot_c2f

   real(c_double) function shr_orb_cosz_eot_c2f(jday, lat, lon, declin, dt_avg, eot_seconds) bind(C)
      real(c_double), value :: jday, lat, lon, declin, dt_avg, eot_seconds
      shr_orb_cosz_eot_c2f = shr_orb_cosz(jday, lat, lon, declin, dt_avg, eot_seconds=eot_seconds)
   end function shr_orb_cosz_eot_c2f

   real(c_double) function shr_orb_cosz_c2f( &
         jday, lat, lon, declin, dt_avg   &
         ) bind(C, name='shr_orb_cosz_c2f')
      real(c_double), VALUE :: jday
      real(c_double), VALUE :: lat
      real(c_double), VALUE :: lon
      real(c_double), VALUE :: declin
      real(c_double), VALUE :: dt_avg
      shr_orb_cosz_c2f = shr_orb_cosz(jday, lat, lon, declin, dt_avg)
      return
   end function shr_orb_cosz_c2f

end module shr_orb_mod_c2f
