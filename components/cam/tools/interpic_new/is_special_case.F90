logical function is_special_case (name)

   implicit none

   include 'netcdf.inc'

   ! Input arguments

   character(len=*), intent(in) :: name
   !---------------------------------------------------------------------

   is_special_case = .false.

   ! Hardwire variable names known to be functions of spatial dimensions

   if (name == 'rlon' .or. name == 'nlon'   .or. name == 'wnummax' .or. &
       name == 'hyai' .or. name == 'hybi'   .or. &
       name == 'hyam' .or. name == 'hybm'   .or. &
       name == 'gw'   .or. name == 'w_stag' .or. &
       name == 'lat'  .or. name == 'lon'    .or. &
       name == 'slat' .or. name == 'slon'   .or. &
       name == 'lev'  .or. name == 'ilev') then

      is_special_case = .true.

   end if

end function is_special_case
