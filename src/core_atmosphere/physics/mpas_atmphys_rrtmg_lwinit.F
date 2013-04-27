! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
#define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_reals(dminfo,size(A),A)

!==================================================================================================
 module mpas_atmphys_rrtmg_lwinit
 use mpas_dmpar
 use mpas_kind_types
 use mpas_grid_types
 use mpas_atmphys_constants
 use mpas_atmphys_utilities

!wrf physics
 use module_ra_rrtmg_lw

 implicit none
 private
 public:: rrtmg_initlw_forMPAS


!>\brief main initialization module for the RRTMG long wave radiation code, based on the MPI
!> decomposition used in MPAS. It replaces the initialization in module_ra_rrtgm_lw.F.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_rrtmg_lwinit:
!> -----------------------------------------
!> rrtmg_lwinit_forMPAS: call rrtmg_lwlookuptable from subroutine init_radiation_lw.
!> rrtmg_lwlookuptable : read and broadcast all input data on MPAS nodes using MPI decomposition.
!>
 

 contains

!==================================================================================================
 subroutine rrtmg_initlw_forMPAS(dminfo)
!==================================================================================================

!input arguments:
 type(dm_info):: dminfo

!--------------------------------------------------------------------------------------------------

!read in absorption coefficients and other data:
 call rrtmg_lwlookuptable(dminfo)

!Perform g-point reduction and other initializations: specific heat of dry air (cp) used in
!flux to heating rate conversion factor.
 call rrtmg_lw_ini(cp)

 end subroutine rrtmg_initlw_forMPAS

!==================================================================================================
 subroutine rrtmg_lwlookuptable(dminfo)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo

!local variables:
 integer:: i,istat,rrtmg_unit
 logical:: opened
 character(len=StrKIND):: errmess
!--------------------------------------------------------------------------------------------------

!get a unit to open init file:
 istat = -999
 if(dminfo % my_proc_id == IO_NODE) then
    do i = 10,99
       inquire(i,opened = opened,iostat=istat)
       if(.not. opened ) then
          rrtmg_unit = i
          exit
       endif
    enddo
    if(istat /= 0) &
       call physics_error_fatal('module_ra_rrtmg_lw: rrtm_lwlookuptable: Can not '// &
                                'find unused fortran unit to read in lookup table.' )

 endif

!distribute unit to other processors:
 call mpas_dmpar_bcast_int(dminfo,rrtmg_unit)

!open init file:
 if(dminfo % my_proc_id == IO_NODE) then
#ifdef SINGLE_PRECISION
    open(rrtmg_unit,file='RRTMG_LW_DATA.TBL',form='UNFORMATTED',status='OLD',iostat=istat)
#else
    open(rrtmg_unit,file='RRTMG_LW_DATA.DBL',form='UNFORMATTED',status='OLD',iostat=istat)
#endif

    if(istat /= 0) then
       write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit', &
             rrtmg_unit
       call physics_error_fatal(errmess)
    endif
 endif

 call lw_kgb01(rrtmg_unit,dminfo)
 call lw_kgb02(rrtmg_unit,dminfo)
 call lw_kgb03(rrtmg_unit,dminfo)
 call lw_kgb04(rrtmg_unit,dminfo)
 call lw_kgb05(rrtmg_unit,dminfo)
 call lw_kgb06(rrtmg_unit,dminfo)
 call lw_kgb07(rrtmg_unit,dminfo)
 call lw_kgb08(rrtmg_unit,dminfo)
 call lw_kgb09(rrtmg_unit,dminfo)
 call lw_kgb10(rrtmg_unit,dminfo)
 call lw_kgb11(rrtmg_unit,dminfo)
 call lw_kgb12(rrtmg_unit,dminfo)
 call lw_kgb13(rrtmg_unit,dminfo)
 call lw_kgb14(rrtmg_unit,dminfo)
 call lw_kgb15(rrtmg_unit,dminfo)
 call lw_kgb16(rrtmg_unit,dminfo)

 if(dminfo % my_proc_id == IO_NODE) close(rrtmg_unit)

 end subroutine rrtmg_lwlookuptable

!==================================================================================================

! **************************************************************************     
!  RRTMG Longwave Radiative Transfer Model
!  Atmospheric and Environmental Research, Inc., Cambridge, MA
!
!  Original version:   E. J. Mlawer, et al.
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90 formatting:  Michael J. Iacono; June 2006
!
!  This file contains 16 READ statements that include the 
!  absorption coefficients and other data for each of the 16 longwave
!  spectral bands used in RRTMG_LW.  Here, the data are defined for 16
!  g-points, or sub-intervals, per band.  These data are combined and
!  weighted using a mapping procedure in module RRTMG_LW_INIT to reduce
!  the total number of g-points from 256 to 140 for use in the GCM.
! **************************************************************************     

! **************************************************************************
      subroutine lw_kgb01(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg01, only : fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
                           absa, absb, &
                      selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: P = 212.7250 mbar, T = 223.06 K

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  The first
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the corresponding TREF for this  pressure level, 
!     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30, 
!     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second 
!     index, JP, runs from 1 to 13 and refers to the corresponding 
!     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).  
!     The third index, IG, goes from 1 to 16, and tells us which 
!     g-interval the absorption coefficients are for.

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

!     The arrays kao_mn2 and kbo_mn2 contain the coefficients of the 
!     nitrogen continuum for the upper and lower atmosphere.
!     Minor gas mapping levels: 
!     Lower - n2: P = 142.5490 mbar, T = 215.70 K
!     Upper - n2: P = 142.5490 mbar, T = 215.70 K

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
         fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, selfrefo, forrefo

      if(istat /= 0) then
         write(0,*)
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(kao_mn2)
      DM_BCAST_MACRO(kbo_mn2)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb01

! **************************************************************************
      subroutine lw_kgb02(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg02, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P = 1053.630 mbar, T = 294.2 K
!     Upper: P = 3.206e-2 mb, T = 197.92 K

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  The first
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the corresponding TREF for this  pressure level, 
!     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30, 
!     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second 
!     index, JP, runs from 1 to 13 and refers to the corresponding 
!     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).  
!     The third index, IG, goes from 1 to 16, and tells us which 
!     g-interval the absorption coefficients are for.

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
         fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb02

! **************************************************************************
      subroutine lw_kgb03(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg03, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
                            kbo_mn2o, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P = 212.7250 mbar, T = 223.06 K
!     Upper: P = 95.8 mbar, T = 215.7 k

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

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
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array KAO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.  The second index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The third index 
!     runs over the g-channel (1 to 16).

!     The array KBO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amounts ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 to 
!     that of gas2.  The second index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The third index 
!     runs over the g-channel (1 to 16).

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
         fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(kao_mn2o)
      DM_BCAST_MACRO(kbo_mn2o)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb03 

! **************************************************************************
      subroutine lw_kgb04(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg04, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower : P = 142.5940 mbar, T = 215.70 K
!     Upper : P = 95.58350 mb, T = 215.70 K

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels  < ~100mb, temperatures, and ratios
!     of H2O to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.  The second index, JT, which
!     runs from 1 to 5, corresponds to different temperatures.  More 
!     specifically, JT = 3 means that the data are for the corresponding 
!     reference temperature TREF for this  pressure level, JT = 2 refers 
!     to the TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and
!     JT = 5 is for TREF+30.  The third index, JP, runs from 13 to 59 and
!     refers to the corresponding pressure level in PREF (e.g. JP = 13 is
!     for a pressure of 95.5835 mb).  The fourth index, IG, goes from 1 to
!     16, and tells us which g-interval the absorption coefficients are for.

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
         fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb04

! **************************************************************************
      subroutine lw_kgb05(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg05, only : fracrefao, fracrefbo, kao, kbo, kao_mo3, &
                            selfrefo, forrefo, ccl4o
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P = 473.42 mb, T = 259.83
!     Upper: P = 0.2369280 mbar, T = 253.60 K

!     The arrays kao_mo3 and ccl4o contain the coefficients for
!     ozone and ccl4 in the lower atmosphere.
!     Minor gas mapping level:
!     Lower - o3: P = 317.34 mbar, T = 240.77 k
!     Lower - ccl4:

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array KBO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels  < ~100mb, temperatures, and ratios
!     of H2O to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.  The second index, JT, which
!     runs from 1 to 5, corresponds to different temperatures.  More 
!     specifically, JT = 3 means that the data are for the corresponding 
!     reference temperature TREF for this  pressure level, JT = 2 refers 
!     to the TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and
!     JT = 5 is for TREF+30.  The third index, JP, runs from 13 to 59 and
!     refers to the corresponding pressure level in PREF (e.g. JP = 13 is
!     for a pressure of 95.5835 mb).  The fourth index, IG, goes from 1 to
!     16, and tells us which g-interval the absorption coefficients are for.

!     The array KAO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.  The second index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The third index 
!     runs over the g-channel (1 to 16).

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
         fracrefao, fracrefbo, kao, kbo, kao_mo3, ccl4o, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(kao_mo3)
      DM_BCAST_MACRO(ccl4o)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb05

! **************************************************************************
      subroutine lw_kgb06(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg06, only : fracrefao, kao, kao_mco2, selfrefo, forrefo, &
                            cfc11adjo, cfc12o
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: : P = 473.4280 mb, T = 259.83 K

!     The arrays kao_mco2, cfc11adjo and cfc12o contain the coefficients for
!     carbon dioxide in the lower atmosphere and cfc11 and cfc12 in the upper
!     atmosphere.
!     Original cfc11 is multiplied by 1.385 to account for the 1060-1107 cm-1 band.
!     Minor gas mapping level:
!     Lower - co2: P = 706.2720 mb, T = 294.2 k
!     Upper - cfc11, cfc12

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  The first
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the corresponding TREF for this  pressure level, 
!     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30, 
!     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second 
!     index, JP, runs from 1 to 13 and refers to the corresponding 
!     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).  
!     The third index, IG, goes from 1 to 16, and tells us which 
!     g-interval the absorption coefficients are for.

!     The array KAO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

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
         fracrefao, kao, kao_mco2, cfc11adjo, cfc12o, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kao_mco2)
      DM_BCAST_MACRO(cfc11adjo)
      DM_BCAST_MACRO(cfc12o)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb06

! **************************************************************************
      subroutine lw_kgb07(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg07, only : fracrefao, fracrefbo, kao, kbo, kao_mco2, &
                            kbo_mco2, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower : P = 706.27 mb, T = 278.94 K
!     Upper : P = 95.58 mbar, T= 215.70 K

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

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

!     The array KAO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.  The second index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The third index 
!     runs over the g-channel (1 to 16).

!     The array KBO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

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
         fracrefao, fracrefbo, kao, kbo, kao_mco2, kbo_mco2, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(kao_mco2)
      DM_BCAST_MACRO(kbo_mco2)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb07

! **************************************************************************
      subroutine lw_kgb08(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg08, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
                            kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
                            cfc12o, cfc22adjo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P=473.4280 mb, T = 259.83 K
!     Upper: P=95.5835 mb, T= 215.7 K

!     The arrays kao_mco2, kbo_mco2, kao_mn2o, kbo_mn2o contain the coefficients for
!     carbon dioxide and n2o in the lower and upper atmosphere.
!     The array kao_mo3 contains the coefficients for ozone in the lower atmosphere,
!     and arrays cfc12o and cfc12adjo contain the coefficients for cfc12 and cfc22.
!     Original cfc22 is multiplied by 1.485 to account for the 780-850 cm-1 
!     and 1290-1335 cm-1 bands.
!     Minor gas mapping level:
!     Lower - co2: P = 1053.63 mb, T = 294.2 k
!     Lower - o3: P = 317.348 mb, T = 240.77 k
!     Lower - n2o: P = 706.2720 mb, T= 278.94 k
!     Lower - cfc12, cfc22
!     Upper - co2: P = 35.1632 mb, T = 223.28 k
!     Upper - n2o: P = 8.716e-2 mb, T = 226.03 k

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  The first
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the corresponding TREF for this  pressure level, 
!     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30, 
!     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second 
!     index, JP, runs from 1 to 13 and refers to the corresponding 
!     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).  
!     The third index, IG, goes from 1 to 16, and tells us which 
!     g-interval the absorption coefficients are for.

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

!     The array KAO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

!     The array KBO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

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
         fracrefao, fracrefbo, kao, kbo, kao_mco2, kbo_mco2, kao_mn2o, &
         kbo_mn2o, kao_mo3, cfc12o, cfc22adjo, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(kao_mco2)
      DM_BCAST_MACRO(kbo_mco2)
      DM_BCAST_MACRO(kao_mn2o)
      DM_BCAST_MACRO(kbo_mn2o)
      DM_BCAST_MACRO(kao_mo3)
      DM_BCAST_MACRO(cfc12o)
      DM_BCAST_MACRO(cfc22adjo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb08

! **************************************************************************
      subroutine lw_kgb09(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg09, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
                            kbo_mn2o, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P=212.7250 mb, T = 223.06 K
!     Upper: P=3.20e-2 mb, T = 197.92 k

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

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

!     The array KAO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.  The second index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The third index 
!     runs over the g-channel (1 to 16).

!     The array KBO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

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
         fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(kao_mn2o)
      DM_BCAST_MACRO(kbo_mn2o)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb09

! **************************************************************************
      subroutine lw_kgb10(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg10, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P = 212.7250 mb, T = 223.06 K
!     Upper: P = 95.58350 mb, T = 215.70 K

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  The first
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the corresponding TREF for this  pressure level, 
!     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30, 
!     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second 
!     index, JP, runs from 1 to 13 and refers to the corresponding 
!     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).  
!     The third index, IG, goes from 1 to 16, and tells us which 
!     g-interval the absorption coefficients are for.

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
         fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb10

! **************************************************************************
      subroutine lw_kgb11(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg11, only : fracrefao, fracrefbo, kao, kbo, kao_mo2, &
                            kbo_mo2, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P=1053.63 mb, T= 294.2 K
!     Upper: P=0.353 mb, T = 262.11 K

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  The first
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the corresponding TREF for this  pressure level, 
!     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30, 
!     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second 
!     index, JP, runs from 1 to 13 and refers to the corresponding 
!     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).  
!     The third index, IG, goes from 1 to 16, and tells us which 
!     g-interval the absorption coefficients are for.

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

!     The array KAO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

!     The array KBO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

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
         fracrefao, fracrefbo, kao, kbo, kao_mo2, kbo_mo2, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(kao_mo2)
      DM_BCAST_MACRO(kbo_mo2)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb11

! **************************************************************************
      subroutine lw_kgb12(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg12, only : fracrefao, kao, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P = 174.1640 mbar, T= 215.78 K

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
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
         fracrefao, kao, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb12

! **************************************************************************
      subroutine lw_kgb13(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg13, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
                            kbo_mo3, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in):: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P=473.4280 mb, T = 259.83 K      
!     Upper: P=4.758820 mb, T = 250.85 K

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array KAO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.  The second index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The third index 
!     runs over the g-channel (1 to 16).

!     The array KBO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

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
         fracrefao, fracrefbo, kao, kao_mco2, kao_mco, kbo_mo3, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kao_mco2)
      DM_BCAST_MACRO(kao_mco)
      DM_BCAST_MACRO(kbo_mo3)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb13

! **************************************************************************
      subroutine lw_kgb14(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg14, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in) :: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P = 142.5940 mb, T = 215.70 K
!     Upper: P = 4.758820 mb, T = 250.85 K

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

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
         fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb14

! **************************************************************************
      subroutine lw_kgb15(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg15, only : fracrefao, kao, kao_mn2, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer, intent(in) :: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P = 1053. mb, T = 294.2 K

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

!     The array KA_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.  The second index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The third index 
!     runs over the g-channel (1 to 16).

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
         fracrefao, kao, kao_mn2, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kao_mn2)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb15

! **************************************************************************
      subroutine lw_kgb16(rrtmg_unit,dminfo)
! **************************************************************************

      use rrlw_kg16, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      save

!input arguments:
 type(dm_info),intent(in):: dminfo
 integer,intent(in) :: rrtmg_unit

!local variables:                                    
      character(len=StrKIND):: errmess
      integer:: istat

!     Arrays fracrefao and fracrefbo are the Planck fractions for the lower
!     and upper atmosphere.
!     Planck fraction mapping levels: 
!     Lower: P = 387.6100 mbar, T = 250.17 K
!     Upper: P=95.58350 mb, T = 215.70 K

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference MLS column amount value of gas 1 
!     to that of gas2.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this  pressure 
!     level, JT = 2 refers to the temperature
!     TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.

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
         fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo

      if(istat /= 0) then
         write(errmess,'(A,I4)') 'module_ra_rrtmg_lw: error reading RRTMG_LW_DATA on unit ', &
               rrtmg_unit
         call physics_error_fatal(errmess)
      endif

      DM_BCAST_MACRO(fracrefao)
      DM_BCAST_MACRO(fracrefbo)
      DM_BCAST_MACRO(kao)
      DM_BCAST_MACRO(kbo)
      DM_BCAST_MACRO(selfrefo)
      DM_BCAST_MACRO(forrefo)

      end subroutine lw_kgb16

!==================================================================================================
 end module mpas_atmphys_rrtmg_lwinit
!==================================================================================================
