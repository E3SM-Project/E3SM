module types
   use prec
!
! 1 for sst, 1 for ice
!
type icesstparms
   character(len=8) :: inname               ! variable name on input file
   integer :: invarid                       ! netcdf varid on input file
!
! Climatological
!
   character(len=32) :: climname            ! output variable name
   character(len=32) :: climname_pre        ! output variable name prior to time diddling
   character(len=128) :: climlongname       ! output variable long name
   character(len=128) :: climlongname_pre   ! output variable long name prior to time diddling
!
! AMIP
!  
   character(len=32) :: amipname            ! output variable name
   character(len=32) :: amipname_pre        ! output variable name prior to time diddling     
   character(len=128) :: amiplongname       ! output variable long name                       
   character(len=128) :: amiplongname_pre   ! output variable long name prior to time diddling

   character(len=8) :: units                ! "fraction" for sea ice, "deg_C" for sst

   integer :: varidclim                     ! variable id goes with climname
   integer :: varidclim_pre                 ! variable id goes with climname_pre
   integer :: varidamip                     ! variable id goes with amipname
   integer :: varidamip_pre                 ! variable id goes with amipname_pre

   real(r8) :: conv   ! convergence criteria
   real(r8) :: tmin   ! min limit to which clipping will be applied
   real(r8) :: tmax   ! max limit to which clipping will be applied
   real(r8) :: varmin ! min variance for grid cell to be included in area means
   real(r8) :: dt     ! smoothing when values go from mininum to maximum
                      ! maximum jump in monthly means allowed is tmax-tmin-dt  (in output units)

!JR correl: How to decay toward climatology in buffer zones
!JR what if buffer zone is longer than 1 year?  Karl will check that the
!           following are from AMIP II data on original 1x1 degree grid.
!JR ice items 6-12 are NOT taken from observations
!JR sst items 9-12 are NOT taken from observations

   real(r8) :: correl(12)
end type icesstparms
end module types
