module shr_orb_mod_c2f

   use iso_c_binding
   use shr_orb_mod, only: shr_orb_params, shr_orb_decl, shr_orb_cosz, SHR_ORB_UNDEF_INT
   implicit none
   public :: shr_orb_params_c2f, shr_orb_decl_c2f, shr_orb_cosz_c2f
   integer(c_int), bind(C) :: shr_orb_undef_int_c2f = SHR_ORB_UNDEF_INT

contains

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
