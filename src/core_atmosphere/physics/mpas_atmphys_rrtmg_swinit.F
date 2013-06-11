! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
#define DM_BCAST_INTEGER(A) call mpas_dmpar_bcast_int(dminfo,A)
#define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_reals(dminfo,size(A),A)
#define DM_BCAST_REAL(A) call mpas_dmpar_bcast_real(dminfo,A)

!=============================================================================================
 module mpas_atmphys_rrtmg_swinit
 use mpas_dmpar
 use mpas_kind_types
 use mpas_grid_types
 use mpas_atmphys_constants
 use mpas_atmphys_utilities

!wrf physics
 use module_ra_rrtmg_sw

 implicit none
 private
 public:: rrtmg_initsw_forMPAS

 contains

!=============================================================================================
 subroutine rrtmg_initsw_forMPAS(dminfo)
!=============================================================================================
!input arguments:
 type(dm_info):: dminfo

!---------------------------------------------------------------------------------------------

!read in absorption coefficients and other data:
 call rrtmg_swlookuptable(dminfo)

!Perform g-point reduction and other initializations: specific heat of dry air (cp) used in
!flux to heating rate conversion factor.
 call rrtmg_sw_ini(cp)

 end subroutine rrtmg_initsw_forMPAS

!=============================================================================================
 subroutine rrtmg_swlookuptable(dminfo)
!=============================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo

!local variables:
 integer:: i,istat,rrtmg_unit
 logical:: opened
 character(len=StrKIND):: errmess
!---------------------------------------------------------------------------------------------

!get a unit to open init file:
 if(dminfo % my_proc_id == IO_NODE) then
    do i = 10,99
       inquire(i,opened = opened,iostat=istat)
       if(.not. opened) then
          rrtmg_unit = i
          exit
       endif
    enddo
    if(istat /= 0) &
       call physics_error_fatal('module_ra_rrtmg_sw: rrtm_swlookuptable: Can not '// &
                                'find unused fortran unit to read in lookup table.' )
 endif

!distribute unit to other processors:
 call mpas_dmpar_bcast_int(dminfo,rrtmg_unit)

!open init file:
 if(dminfo % my_proc_id == IO_NODE) then
    open(rrtmg_unit,file='RRTMG_SW_DATA.DBL',form='UNFORMATTED',iostat=istat)

    if(istat /= 0) then
         write(0,*) 'rrtmg_unit=',rrtmg_unit
         write(0,*) 'istat=',istat
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
    endif
 endif

 call sw_kgb16(rrtmg_unit,dminfo)
 call sw_kgb17(rrtmg_unit,dminfo)
 call sw_kgb18(rrtmg_unit,dminfo)
 call sw_kgb19(rrtmg_unit,dminfo)
 call sw_kgb20(rrtmg_unit,dminfo)
 call sw_kgb21(rrtmg_unit,dminfo)
 call sw_kgb22(rrtmg_unit,dminfo)
 call sw_kgb23(rrtmg_unit,dminfo)
 call sw_kgb24(rrtmg_unit,dminfo)
 call sw_kgb25(rrtmg_unit,dminfo)
 call sw_kgb26(rrtmg_unit,dminfo)
 call sw_kgb27(rrtmg_unit,dminfo)
 call sw_kgb28(rrtmg_unit,dminfo)
 call sw_kgb29(rrtmg_unit,dminfo)

 if(dminfo % my_proc_id == IO_NODE) close(rrtmg_unit)

 end subroutine rrtmg_swlookuptable

!=============================================================================================

! **************************************************************************
!  RRTMG Shortwave Radiative Transfer Model
!  Atmospheric and Environmental Research, Inc., Cambridge, MA
!
!  Original by J.Delamere, Atmospheric & Environmental Research.
!  Reformatted for F90: JJMorcrette, ECMWF
!  Revision for GCMs:  Michael J. Iacono, AER, July 2002
!  Further F90 reformatting:  Michael J. Iacono, AER, June 2006
!
!  This file contains 14 READ statements that include the 
!  absorption coefficients and other data for each of the 14 shortwave
!  spectral bands used in RRTMG_SW.  Here, the data are defined for 16
!  g-points, or sub-intervals, per band.  These data are combined and
!  weighted using a mapping procedure in module RRTMG_SW_INIT to reduce
!  the total number of g-points from 224 to 112 for use in the GCM.
! **************************************************************************

! **************************************************************************
      subroutine sw_kgb16(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg16, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            rayl, strrat1, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array rayl contains the Rayleigh extinction coefficient at v = 2925 cm-1.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         rayl, strrat1, layreffr, kao, kbo, selfrefo, forrefo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_REAL(rayl)
      DM_BCAST_REAL(strrat1)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb16

! **************************************************************************
      subroutine sw_kgb17(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg17, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            rayl, strrat, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array rayl contains the Rayleigh extinction coefficient at v = 3625 cm-1.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         rayl, strrat, layreffr, kao, kbo, selfrefo, forrefo, sfluxrefo


      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_REAL(rayl)
      DM_BCAST_REAL(strrat)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb17

! **************************************************************************
      subroutine sw_kgb18(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg18, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            rayl, strrat, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array rayl contains the Rayleigh extinction coefficient at v = 4325 cm-1.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         rayl, strrat, layreffr, kao, kbo, selfrefo, forrefo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_REAL(rayl)
      DM_BCAST_REAL(strrat)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb18 

! **************************************************************************
      subroutine sw_kgb19(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg19, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            rayl, strrat, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array rayl contains the Rayleigh extinction coefficient at v = 4900 cm-1.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         rayl, strrat, layreffr, kao, kbo, selfrefo, forrefo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_REAL(rayl)
      DM_BCAST_REAL(strrat)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb19

! **************************************************************************
      subroutine sw_kgb20(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg20, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            absch4o, rayl, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array rayl contains the Rayleigh extinction coefficient at v = 5670 cm-1.

!     Array absch4o contains the absorption coefficients for methane.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         rayl, layreffr, absch4o, kao, kbo, selfrefo, forrefo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_REAL(rayl)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(absch4o)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb20

! **************************************************************************
      subroutine sw_kgb21(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg21, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            rayl, strrat, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array rayl contains the Rayleigh extinction coefficient at v = 6925 cm-1.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         rayl, strrat, layreffr, kao, kbo, selfrefo, forrefo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_REAL(rayl)
      DM_BCAST_REAL(strrat)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb21

! **************************************************************************
      subroutine sw_kgb22(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg22, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            rayl, strrat, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array rayl contains the Rayleigh extinction coefficient at v = 8000 cm-1.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296_rb,260_rb,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         rayl, strrat, layreffr, kao, kbo, selfrefo, forrefo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_REAL(rayl)
      DM_BCAST_REAL(strrat)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb22

! **************************************************************************
      subroutine sw_kgb23(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg23, only : kao, selfrefo, forrefo, sfluxrefo, &
                            raylo, givfac, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array raylo contains the Rayleigh extinction coefficient at all v for this band

!     Array givfac is the average Giver et al. correction factor for this band. 

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         raylo, givfac, layreffr, kao, selfrefo, forrefo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(raylo)
      DM_BCAST_REAL(givfac)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb23

! **************************************************************************
      subroutine sw_kgb24(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg24, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            raylao, raylbo, abso3ao, abso3bo, strrat, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Arrays raylao and raylbo contain the Rayleigh extinction coefficient at 
!     all v for this band for the upper and lower atmosphere.

!     Arrays abso3ao and abso3bo contain the ozone absorption coefficient at 
!     all v for this band for the upper and lower atmosphere.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         raylao, raylbo, strrat, layreffr, abso3ao, abso3bo, kao, kbo, selfrefo, &
         forrefo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(raylao)
      DM_BCAST_MACRO(raylbo)
      DM_BCAST_REAL(strrat)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(abso3ao)
      DM_BCAST_MACRO(abso3bo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb24

! **************************************************************************
      subroutine sw_kgb25(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg25, only : kao, sfluxrefo, &
                            raylo, abso3ao, abso3bo, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array raylo contains the Rayleigh extinction coefficient at all v = 2925 cm-1.

!     Arrays abso3ao and abso3bo contain the ozone absorption coefficient at 
!     all v for this band for the upper and lower atmosphere.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         raylo, layreffr, abso3ao, abso3bo, kao, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(raylo)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(abso3ao)
      DM_BCAST_MACRO(abso3bo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb25

! **************************************************************************
      subroutine sw_kgb26(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg26, only : sfluxrefo, raylo
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array raylo contains the Rayleigh extinction coefficient at all v for this band.

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         raylo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(raylo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb26

! **************************************************************************
      subroutine sw_kgb27(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg27, only : kao, kbo, sfluxrefo, raylo, &
                            scalekur, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 
!     The values in array sfluxrefo were obtained using the "low resolution"
!     version of the Kurucz solar source function.  For unknown reasons,
!     the total irradiance in this band differs from the corresponding
!     total in the "high-resolution" version of the Kurucz function.
!     Therefore, these values are scaled by the factor SCALEKUR.

!     Array raylo contains the Rayleigh extinction coefficient at all v = 2925 cm-1.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         raylo, scalekur, layreffr, kao, kbo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(raylo)
      DM_BCAST_REAL(scalekur)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb27

! **************************************************************************
      subroutine sw_kgb28(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg28, only : kao, kbo, sfluxrefo, &
                            rayl, strrat, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array raylo contains the Rayleigh extinction coefficient at all v = ???? cm-1.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         rayl, strrat, layreffr, kao, kbo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_REAL(rayl)
      DM_BCAST_REAL(strrat)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(sfluxrefo)

      end subroutine sw_kgb28

! **************************************************************************
      subroutine sw_kgb29(rrtmg_unit,dminfo)
! **************************************************************************

      use rrsw_kg29, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            absh2oo, absco2o, rayl, layreffr
      save

!input arguments:
      type(dm_info),intent(in):: dminfo
      integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Array sfluxrefo contains the Kurucz solar source function for this band. 

!     Array rayl contains the Rayleigh extinction coefficient at all v = 2200 cm-1.

!     Array absh2oo contains the water vapor absorption coefficient for this band.

!     Array absco2o contains the carbon dioxide absorption coefficient for this band.

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

!     write(0,*)
!     write(0,*) '--- enter subroutine sw_kgb29:'
!     write(0,*) '--- rrtmg_unit= ', rrtmg_unit

      istat = 0
      if(dminfo % my_proc_id == IO_NODE) read (unit=rrtmg_unit,iostat=istat) &
         rayl, layreffr, absh2oo, absco2o, kao, kbo, selfrefo, forrefo, sfluxrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_sw: error reading RRTMG_SW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_REAL(rayl)
      DM_BCAST_INTEGER(layreffr)
      DM_BCAST_MACRO(absh2oo)
      DM_BCAST_MACRO(absco2o)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)
      DM_BCAST_MACRO(sfluxrefo)
      
      end subroutine sw_kgb29

!=============================================================================================
 end module mpas_atmphys_rrtmg_swinit
!=============================================================================================
