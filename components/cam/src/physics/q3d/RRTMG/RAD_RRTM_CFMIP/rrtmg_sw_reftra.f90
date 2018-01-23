!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_reftra.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.5 $
!     created:   $Date: 2009/05/22 22:22:22 $

      module rrtmg_sw_reftra

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use parkind, only : im => kind_im, rb => kind_rb
      use rrsw_tbl, only : tblint, bpade, od_lo, exp_tbl
      use rrsw_vsn, only : hvrrft, hnamrft

      implicit none

      contains

! --------------------------------------------------------------------
      subroutine reftra_sw(nlayers, lrtchk, pgg, prmuz, ptau, pw, &
                           pref, prefd, ptra, ptrad)
! --------------------------------------------------------------------
  
! Purpose: computes the reflectivity and transmissivity of a clear or 
!   cloudy layer using a choice of various approximations.
!
! Interface:  *rrtmg_sw_reftra* is called by *rrtmg_sw_spcvrt*
!
! Description:
! explicit arguments :
! --------------------
! inputs
! ------ 
!      lrtchk  = .t. for all layers in clear profile
!      lrtchk  = .t. for cloudy layers in cloud profile 
!              = .f. for clear layers in cloud profile
!      pgg     = assymetry factor
!      prmuz   = cosine solar zenith angle
!      ptau    = optical thickness
!      pw      = single scattering albedo
!
! outputs
! -------
!      pref    : collimated beam reflectivity
!      prefd   : diffuse beam reflectivity 
!      ptra    : collimated beam transmissivity
!      ptrad   : diffuse beam transmissivity
!
!
! Method:
! -------
!      standard delta-eddington, p.i.f.m., or d.o.m. layer calculations.
!      kmodts  = 1 eddington (joseph et al., 1976)
!              = 2 pifm (zdunkowski et al., 1980)
!              = 3 discrete ordinates (liou, 1973)
!
!
! Modifications:
! --------------
! Original: J-JMorcrette, ECMWF, Feb 2003
! Revised for F90 reformatting: MJIacono, AER, Jul 2006
! Revised to add exponential lookup table: MJIacono, AER, Aug 2007
!
! ------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------

      integer(kind=im), intent(in) :: nlayers

      logical, intent(in) :: lrtchk(:)                         ! Logical flag for reflectivity and
                                                               ! and transmissivity calculation; 
                                                               !   Dimensions: (nlayers)

      real(kind=rb), intent(in) :: pgg(:)                      ! asymmetry parameter
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: ptau(:)                     ! optical depth
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: pw(:)                       ! single scattering albedo 
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: prmuz                       ! cosine of solar zenith angle

! ------- Output -------

      real(kind=rb), intent(inout) :: pref(:)                  ! direct beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: prefd(:)                 ! diffuse beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: ptra(:)                  ! direct beam transmissivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: ptrad(:)                 ! diffuse beam transmissivity
                                                               !   Dimensions: (nlayers+1)

! ------- Local -------

      integer(kind=im) :: jk, jl, kmodts
      integer(kind=im) :: itind

      real(kind=rb) :: tblind
      real(kind=rb) :: za, za1, za2
      real(kind=rb) :: zbeta, zdend, zdenr, zdent
      real(kind=rb) :: ze1, ze2, zem1, zem2, zemm, zep1, zep2
      real(kind=rb) :: zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4, zgt
      real(kind=rb) :: zr1, zr2, zr3, zr4, zr5
      real(kind=rb) :: zrk, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
      real(kind=rb) :: zsr3, zt1, zt2, zt3, zt4, zt5, zto1
      real(kind=rb) :: zw, zwcrit, zwo

      real(kind=rb), parameter :: eps = 1.e-08_rb

!     ------------------------------------------------------------------

! Initialize

      hvrrft = '$Revision: 1.5 $'

      zsr3=sqrt(3._rb)
      zwcrit=0.9999995_rb
      kmodts=2

      do jk=1, nlayers
         if (.not.lrtchk(jk)) then
            pref(jk) =0._rb
            ptra(jk) =1._rb
            prefd(jk)=0._rb
            ptrad(jk)=1._rb
         else
            zto1=ptau(jk)
            zw  =pw(jk)
            zg  =pgg(jk)  

! General two-stream expressions

            zg3= 3._rb * zg
            if (kmodts == 1) then
               zgamma1= (7._rb - zw * (4._rb + zg3)) * 0.25_rb
               zgamma2=-(1._rb - zw * (4._rb - zg3)) * 0.25_rb
               zgamma3= (2._rb - zg3 * prmuz ) * 0.25_rb
            else if (kmodts == 2) then  
               zgamma1= (8._rb - zw * (5._rb + zg3)) * 0.25_rb
               zgamma2=  3._rb *(zw * (1._rb - zg )) * 0.25_rb
               zgamma3= (2._rb - zg3 * prmuz ) * 0.25_rb
            else if (kmodts == 3) then  
               zgamma1= zsr3 * (2._rb - zw * (1._rb + zg)) * 0.5_rb
               zgamma2= zsr3 * zw * (1._rb - zg ) * 0.5_rb
               zgamma3= (1._rb - zsr3 * zg * prmuz ) * 0.5_rb
            end if
            zgamma4= 1._rb - zgamma3
    
! Recompute original s.s.a. to test for conservative solution

            zwo= zw / (1._rb - (1._rb - zw) * (zg / (1._rb - zg))**2)
    
            if (zwo >= zwcrit) then
! Conservative scattering

               za  = zgamma1 * prmuz 
               za1 = za - zgamma3
               zgt = zgamma1 * zto1
        
! Homogeneous reflectance and transmittance,
! collimated beam

               ze1 = min ( zto1 / prmuz , 500._rb)
!               ze2 = exp( -ze1 )

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               if (ze1 .le. od_lo) then 
                  ze2 = 1._rb - ze1 + 0.5_rb * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_rb
                  ze2 = exp_tbl(itind)
               endif
!

               pref(jk) = (zgt - za1 * (1._rb - ze2)) / (1._rb + zgt)
               ptra(jk) = 1._rb - pref(jk)

! isotropic incidence

               prefd(jk) = zgt / (1._rb + zgt)
               ptrad(jk) = 1._rb - prefd(jk)        

! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
! table returns a transmittance of 1.0.
               if (ze2 .eq. 1.0_rb) then 
                  pref(jk) = 0.0_rb
                  ptra(jk) = 1.0_rb
                  prefd(jk) = 0.0_rb
                  ptrad(jk) = 1.0_rb
               endif

            else
! Non-conservative scattering

               za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
               za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
               zrk = sqrt ( zgamma1**2 - zgamma2**2)
               zrp = zrk * prmuz               
               zrp1 = 1._rb + zrp
               zrm1 = 1._rb - zrp
               zrk2 = 2._rb * zrk
               zrpp = 1._rb - zrp*zrp
               zrkg = zrk + zgamma1
               zr1  = zrm1 * (za2 + zrk * zgamma3)
               zr2  = zrp1 * (za2 - zrk * zgamma3)
               zr3  = zrk2 * (zgamma3 - za2 * prmuz )
               zr4  = zrpp * zrkg
               zr5  = zrpp * (zrk - zgamma1)
               zt1  = zrp1 * (za1 + zrk * zgamma4)
               zt2  = zrm1 * (za1 - zrk * zgamma4)
               zt3  = zrk2 * (zgamma4 + za1 * prmuz )
               zt4  = zr4
               zt5  = zr5

! mji - reformulated code to avoid potential floating point exceptions
!               zbeta = - zr5 / zr4
               zbeta = (zgamma1 - zrk) / zrkg
!!
        
! Homogeneous reflectance and transmittance

               ze1 = min ( zrk * zto1, 500._rb)
               ze2 = min ( zto1 / prmuz , 500._rb)
!
! Original
!              zep1 = exp( ze1 )
!              zem1 = exp(-ze1 )
!              zep2 = exp( ze2 )
!              zem2 = exp(-ze2 )
!
! Revised original, to reduce exponentials
!              zep1 = exp( ze1 )
!              zem1 = 1._rb / zep1
!              zep2 = exp( ze2 )
!              zem2 = 1._rb / zep2
!
! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               if (ze1 .le. od_lo) then 
                  zem1 = 1._rb - ze1 + 0.5_rb * ze1 * ze1
                  zep1 = 1._rb / zem1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_rb
                  zem1 = exp_tbl(itind)
                  zep1 = 1._rb / zem1
               endif

               if (ze2 .le. od_lo) then 
                  zem2 = 1._rb - ze2 + 0.5_rb * ze2 * ze2
                  zep2 = 1._rb / zem2
               else
                  tblind = ze2 / (bpade + ze2)
                  itind = tblint * tblind + 0.5_rb
                  zem2 = exp_tbl(itind)
                  zep2 = 1._rb / zem2
               endif

! collimated beam

! mji - reformulated code to avoid potential floating point exceptions
!               zdenr = zr4*zep1 + zr5*zem1
!               pref(jk) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
!               zdent = zt4*zep1 + zt5*zem1
!               ptra(jk) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent

               zdenr = zr4*zep1 + zr5*zem1
               zdent = zt4*zep1 + zt5*zem1
               if (zdenr .ge. -eps .and. zdenr .le. eps) then
                  pref(jk) = eps
                  ptra(jk) = zem2
               else 
                  pref(jk) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
                  ptra(jk) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent
               endif
!!

! diffuse beam

               zemm = zem1*zem1
               zdend = 1._rb / ( (1._rb - zbeta*zemm ) * zrkg)
               prefd(jk) =  zgamma2 * (1._rb - zemm) * zdend
               ptrad(jk) =  zrk2*zem1*zdend

            endif

         endif         

      enddo    

      end subroutine reftra_sw

      end module rrtmg_sw_reftra

