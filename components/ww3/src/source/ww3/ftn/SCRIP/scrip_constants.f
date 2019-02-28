!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This module defines common constants used in many routines.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: constants.f,v 1.2 2000/04/19 21:56:25 pwjones Exp $
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
!
!     This software and ancillary information (herein called software) 
!     called SCRIP is made available under the terms described here.  
!     The software has been approved for release with associated 
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with 
!     the version available from Los Alamos National Laboratory.
!
!     This code has been modified from the version available from 
!     Los Alamos National Laboratory, for the purpose of running it
!     within WW3.
!
!***********************************************************************

      module SCRIP_constants

!-----------------------------------------------------------------------

      use SCRIP_KindsMod  ! defines common data types

      implicit none

      save

!-----------------------------------------------------------------------

      real (kind = SCRIP_r8), parameter :: 
     &                        zero   = 0.0_SCRIP_r8,
     &                        one    = 1.0_SCRIP_r8,
     &                        two    = 2.0_SCRIP_r8,
     &                        three  = 3.0_SCRIP_r8,
     &                        four   = 4.0_SCRIP_r8,
     &                        five   = 5.0_SCRIP_r8,
     &                        half   = 0.5_SCRIP_r8,
     &                        quart  = 0.25_SCRIP_r8,
     &                        bignum = 1.e+20_SCRIP_r8,
     &                        tiny   = 1.e-14_SCRIP_r8,
     &                        pi     = 3.14159265359_SCRIP_r8,
     &                        pi2    = two*pi,
     &                        pih    = half*pi

!-----------------------------------------------------------------------

      end module SCRIP_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
