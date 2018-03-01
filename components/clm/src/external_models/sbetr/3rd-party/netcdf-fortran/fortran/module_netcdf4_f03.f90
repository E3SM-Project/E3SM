Module netcdf4_f03

! Module to provide inheritance of Netcdf4 data and interfaces from nf_data
! and nf_interfaces modules. Not really needed by netCDF but provided for
! folks writing new code or updating old code who would prefer using
! modules instead of the old netcdf.inc include file.

! Written by: Richard Weed, Ph.D.
!             Center for Advanced Vehicular Systems 
!             Mississippi State University
!             rweed@cavs.msstate.edu
          

! License (and other Lawyer Language)
 
! This software is released under the Apache 2.0 Open Source License. The
! full text of the License can be viewed at :
!
!   http:www.apache.org/licenses/LICENSE-2.0.html
!
! The author grants to the University Corporation for Atmospheric Research
! (UCAR), Boulder, CO, USA the right to revise and extend the software
! without restriction. However, the author retains all copyrights and
! intellectual property rights explicitly stated in or implied by the
! Apache license

! Version 1.0 - April 2009 - Separate module for netcdf 4 for now to
!                            to simplify build process

! This module can be used as a replacement for an include to netcdf.inc
! The module is named netcdf4_f03 to differentiate it from the
! netcdf_f03 module


 USE netcdf_nf_data        ! Brings in the nf interface error flags etc. 
 USE netcdf_nf_interfaces  ! Brings in explicit interfaces to nf_ routines
 USE netcdf4_nf_interfaces ! Bring in netcdf4 interfaces

 Implicit NONE

!------------------------------- End of Module netcdf4_f03 ---------------------
 End Module netcdf4_f03
