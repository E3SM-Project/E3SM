module shr_orb_mod_c

   use iso_c_binding
   use shr_orb_mod, only: shr_orb_params, shr_orb_decl, shr_orb_cosz
   implicit none
   public :: shr_orb_params_c, shr_orb_decl_c, shr_orb_cosz_c

contains

   subroutine shr_orb_params_c(           &
           orb_year, eccen, mvelpp, lambm0, &
           obliqr, delta, eccf            &
           ) bind(C, name='shr_orb_params_c')
      integer(c_int) :: orb_year
      real(c_double) :: eccen
      real(c_double) :: mvelpp
      real(c_double) :: lambm0
      real(c_double) :: obliqr
      real(c_double) :: delta
      real(c_double) :: eccf
      logical :: print_flag
      call shr_orb_params(orb_year, eccen, mvelpp, lambm0, obliqr, delta, eccf, print_flag)
   end subroutine shr_orb_params_c
        
   subroutine shr_orb_decl_c(     &
         calday, eccen, mvelpp, lambm0, &
         obliqr, delta, eccf            &
         ) bind(C, name='shr_orb_decl_c')
      real(c_double) :: calday
      real(c_double) :: eccen
      real(c_double) :: mvelpp
      real(c_double) :: lambm0
      real(c_double) :: obliqr
      real(c_double) :: delta
      real(c_double) :: eccf
      call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
   end subroutine shr_orb_decl_c

   real(c_double) function shr_orb_cosz_c( &
         jday, lat, lon, declin, dt_avg   &
         ) bind(C, name='shr_orb_cosz_c')
      real(c_double), intent(in) :: jday
      real(c_double), intent(in) :: lat
      real(c_double), intent(in) :: lon
      real(c_double), intent(in) :: declin
      real(c_double), intent(in) :: dt_avg
      shr_orb_cosz_c = shr_orb_cosz(jday, lat, lon, declin, dt_avg)
      return
   end function shr_orb_cosz_c

end module shr_orb_mod_c
