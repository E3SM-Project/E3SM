!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module SCRIP_NetCDFMod

!BOP
! !MODULE: SCRIP_NetCDFMod
!
! !DESCRIPTION:
! This module contains netCDF error handling functions and the
! use of the netCDF module for all netCDF functions.
!
! !REVISION HISTORY:
! SVN:$Id: $
!
! !USES:

   use SCRIP_KindsMod
   use SCRIP_ErrorMod
   use netcdf

   implicit none
   private
   save

! !DEFINED PARAMETERS:

! !PUBLIC MEMBER FUNCTIONS:

   public :: SCRIP_NetcdfErrorCheck

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: SCRIP_NetcdfErrorCheck
! !INTERFACE:

   function SCRIP_NetcdfErrorCheck(netcdfStat, errorCode, rtnName, &
                                   errorMsg)

! !DESCRIPTION:
! This routine checks netCDF status flags from netCDF routines.
! On error, it adds an error message to the error log and returns
! both a true value to the calling routines as well as setting the
! error code to fail. It is used in the same manner is the
! SCRIP\_ErrorCheck function.
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

   logical (SCRIP_logical) :: &
      SCRIP_NetcdfErrorCheck

   integer (SCRIP_i4), intent(out) :: &
      errorCode ! returned SCRIP error flag

! !INPUT PARAMETERS:

   integer (SCRIP_i4), intent(in) :: &
      netcdfStat ! status flag from netCDF call

   character (*), intent(in) :: &
      rtnName, &! name of calling routine
      errorMsg ! error message for logging

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   character (SCRIP_charLength) :: &
      ncErrMsg ! netCDF error message

!-----------------------------------------------------------------------
!
! if no error, return false and a successful errorCode
!
!-----------------------------------------------------------------------

   if (netcdfStat == NF90_NOERR) then
      errorCode = SCRIP_Success
      SCRIP_NetcdfErrorCheck = .false.

!-----------------------------------------------------------------------
!
! if an error is detected, return a true value and call the SCRIP
! error handlers to log the error. Log both the netCDF error msg
! as well as the error passed by the calling routine.
!
!-----------------------------------------------------------------------

   else

      SCRIP_NetcdfErrorCheck = .true.
      ncErrMsg = nf90_strerror(netcdfStat)
      call SCRIP_ErrorSet(errorCode, rtnName, ncErrMsg)
      call SCRIP_ErrorSet(errorCode, rtnName, errorMsg)

   endif

!-----------------------------------------------------------------------
!EOC

 end function SCRIP_NetcdfErrorCheck

!***********************************************************************

 end module SCRIP_NetCDFMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
