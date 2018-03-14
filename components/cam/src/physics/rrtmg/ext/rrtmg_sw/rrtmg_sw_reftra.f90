!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_reftra.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:14 $

      module rrtmg_sw_reftra

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb
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

      integer, intent(in) :: nlayers

      logical, intent(in) :: lrtchk(:)                           ! Logical flag for reflectivity and
                                                                 ! and transmissivity calculation; 
                                                                 !   Dimensions: (nlayers)

      real(kind=r8), intent(in) :: pgg(:)                      ! asymmetry parameter
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: ptau(:)                     ! optical depth
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: pw(:)                       ! single scattering albedo 
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: prmuz                       ! cosine of solar zenith angle

! ------- Output -------

      real(kind=r8), intent(inout) :: pref(:)                    ! direct beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prefd(:)                   ! diffuse beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: ptra(:)                    ! direct beam transmissivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: ptrad(:)                   ! diffuse beam transmissivity
                                                                 !   Dimensions: (nlayers+1)

! ------- Local -------

      integer :: jk, jl, kmodts
      integer :: itind

      real(kind=r8) :: tblind
      real(kind=r8) :: za, za1, za2
      real(kind=r8) :: zbeta, zdend, zdenr, zdent
      real(kind=r8) :: ze1, ze2, zem1, zem2, zemm, zep1, zep2
      real(kind=r8) :: zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4, zgt
      real(kind=r8) :: zr1, zr2, zr3, zr4, zr5
      real(kind=r8) :: zrk, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
      real(kind=r8) :: zsr3, zt1, zt2, zt3, zt4, zt5, zto1
      real(kind=r8) :: zw, zwcrit, zwo

      real(kind=r8), parameter :: eps = 1.e-08_r8

!     ------------------------------------------------------------------

! Initialize

      hvrrft = '$Revision: 1.2 $'

      zsr3=sqrt(3._r8)
      zwcrit=0.9999995_r8
      kmodts=2

      do jk=1, nlayers
         if (.not.lrtchk(jk)) then
            pref(jk) =0._r8
            ptra(jk) =1._r8
            prefd(jk)=0._r8
            ptrad(jk)=1._r8
         else
            zto1=ptau(jk)
            zw  =pw(jk)
            zg  =pgg(jk)  

! General two-stream expressions

            zg3= 3._r8 * zg
            if (kmodts == 1) then
               zgamma1= (7._r8 - zw * (4._r8 + zg3)) * 0.25_r8
               zgamma2=-(1._r8 - zw * (4._r8 - zg3)) * 0.25_r8
               zgamma3= (2._r8 - zg3 * prmuz ) * 0.25_r8
            else if (kmodts == 2) then  
               zgamma1= (8._r8 - zw * (5._r8 + zg3)) * 0.25_r8
               zgamma2=  3._r8 *(zw * (1._r8 - zg )) * 0.25_r8
               zgamma3= (2._r8 - zg3 * prmuz ) * 0.25_r8
            else if (kmodts == 3) then  
               zgamma1= zsr3 * (2._r8 - zw * (1._r8 + zg)) * 0.5_r8
               zgamma2= zsr3 * zw * (1._r8 - zg ) * 0.5_r8
               zgamma3= (1._r8 - zsr3 * zg * prmuz ) * 0.5_r8
            end if
            zgamma4= 1._r8 - zgamma3
    
! Recompute original s.s.a. to test for conservative solution

            zwo= zw / (1._r8 - (1._r8 - zw) * (zg / (1._r8 - zg))**2)
    
            if (zwo >= zwcrit) then
! Conservative scattering

               za  = zgamma1 * prmuz 
               za1 = za - zgamma3
               zgt = zgamma1 * zto1
        
! Homogeneous reflectance and transmittance,
! collimated beam

               ze1 = min ( zto1 / prmuz , 500._r8)
               ze2 = exp( -ze1 )

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
!               if (ze1 .le. od_lo) then 
!                  ze2 = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
!               else
!                  tblind = ze1 / (bpade + ze1)
!                  itind = tblint * tblind + 0.5_r8
!                  ze2 = exp_tbl(itind)
!               endif
!

               pref(jk) = (zgt - za1 * (1._r8 - ze2)) / (1._r8 + zgt)
               ptra(jk) = 1._r8 - pref(jk)

! isotropic incidence

               prefd(jk) = zgt / (1._r8 + zgt)
               ptrad(jk) = 1._r8 - prefd(jk)        

! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
! table returns a transmittance of 1.0.
               if (ze2 .eq. 1.0_r8) then 
                  pref(jk) = 0.0_r8
                  ptra(jk) = 1.0_r8
                  prefd(jk) = 0.0_r8
                  ptrad(jk) = 1.0_r8
               endif

            else
! Non-conservative scattering

               za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
               za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
               zrk = sqrt ( zgamma1**2 - zgamma2**2)
               zrp = zrk * prmuz               
               zrp1 = 1._r8 + zrp
               zrm1 = 1._r8 - zrp
               zrk2 = 2._r8 * zrk
               zrpp = 1._r8 - zrp*zrp
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
               zbeta = (zgamma1 - zrk) / zrkg !- zr5 / zr4
        
! Homogeneous reflectance and transmittance

               ze1 = min ( zrk * zto1, 500._r8)
               ze2 = min ( zto1 / prmuz , 500._r8)
!
! Original
!              zep1 = exp( ze1 )
!              zem1 = exp(-ze1 )
!              zep2 = exp( ze2 )
!              zem2 = exp(-ze2 )
!
! Revised original, to reduce exponentials
               zep1 = exp( ze1 )
               zem1 = 1._r8 / zep1
               zep2 = exp( ze2 )
               zem2 = 1._r8 / zep2
!
! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
!               if (ze1 .le. od_lo) then 
!                  zem1 = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
!                  zep1 = 1._r8 / zem1
!               else
!                  tblind = ze1 / (bpade + ze1)
!                  itind = tblint * tblind + 0.5_r8
!                  zem1 = exp_tbl(itind)
!                  zep1 = 1._r8 / zem1
!               endif
!
!               if (ze2 .le. od_lo) then 
!                  zem2 = 1._r8 - ze2 + 0.5_r8 * ze2 * ze2
!                  zep2 = 1._r8 / zem2
!               else
!                  tblind = ze2 / (bpade + ze2)
!                  itind = tblint * tblind + 0.5_r8
!                  zem2 = exp_tbl(itind)
!                  zep2 = 1._r8 / zem2
!               endif

! collimated beam

               zdenr = zr4*zep1 + zr5*zem1
               zdent = zt4*zep1 + zt5*zem1
               if (zdenr .ge. -eps .and. zdenr .le. eps) then
                  pref(jk) = eps
                  ptra(jk) = zem2
               else
                  pref(jk) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
                  ptra(jk) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent
               endif 

! diffuse beam

               zemm = zem1*zem1
               zdend = 1._r8 / ( (1._r8 - zbeta*zemm ) * zrkg)
               prefd(jk) =  zgamma2 * (1._r8 - zemm) * zdend
               ptrad(jk) =  zrk2*zem1*zdend

            endif

         endif         

      enddo    

      end subroutine reftra_sw

      end module rrtmg_sw_reftra

