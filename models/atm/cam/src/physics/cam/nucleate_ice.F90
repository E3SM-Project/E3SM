module nucleate_ice

!---------------------------------------------------------------------------------
! Purpose:
!   Ice nucleation code.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use wv_saturation,  only: svp_water, svp_ice
use cam_logfile,    only: iulog

implicit none
private
save

public :: nucleati

!===============================================================================
contains
!===============================================================================

subroutine nucleati(  &
   wbar, tair, relhum, cldn, qc,              &
   nfice, rhoair, so4_num, dst_num, soot_num, &
   nuci, onihf, oniimm, onidep, onimey) 

   !---------------------------------------------------------------
   ! Purpose:
   !  The parameterization of ice nucleation.
   !
   ! Method: The current method is based on Liu & Penner (2005)
   !  It related the ice nucleation with the aerosol number, temperature and the
   !  updraft velocity. It includes homogeneous freezing of sulfate, immersion
   !  freezing of soot, and Meyers et al. (1992) deposition nucleation
   !
   ! Authors: Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
   !----------------------------------------------------------------

   ! Input Arguments
   real(r8), intent(in) :: wbar        ! grid cell mean vertical velocity (m/s)
   real(r8), intent(in) :: tair        ! temperature (K)
   real(r8), intent(in) :: relhum      ! relative humidity with respective to liquid
   real(r8), intent(in) :: cldn        ! new value of cloud fraction    (fraction)
   real(r8), intent(in) :: qc          ! liquid water mixing ratio (kg/kg)
   real(r8), intent(in) :: nfice       ! ice mass fraction
   real(r8), intent(in) :: rhoair      ! air density (kg/m3)
   real(r8), intent(in) :: so4_num     ! so4 aerosol number (#/cm^3)
   real(r8), intent(in) :: dst_num     ! total dust aerosol number (#/cm^3)
   real(r8), intent(in) :: soot_num    ! soot (hydrophilic) aerosol number (#/cm^3)

   ! Output Arguments
   real(r8), intent(out) :: nuci       ! ice number nucleated (#/kg)
   real(r8), intent(out) :: onihf      ! nucleated number from homogeneous freezing of so4
   real(r8), intent(out) :: oniimm     ! nucleated number from immersion freezing
   real(r8), intent(out) :: onidep     ! nucleated number from deposition nucleation
   real(r8), intent(out) :: onimey     ! nucleated number from deposition nucleation  (meyers: mixed phase)

   ! Local workspace
   real(r8) :: nihf                      ! nucleated number from homogeneous freezing of so4
   real(r8) :: niimm                     ! nucleated number from immersion freezing
   real(r8) :: nidep                     ! nucleated number from deposition nucleation
   real(r8) :: nimey                     ! nucleated number from deposition nucleation (meyers)
   real(r8) :: n1, ni                    ! nucleated number
   real(r8) :: tc, A, B, C, regm, RHw    ! work variable
   real(r8) :: esl, esi, deles           ! work variable
   real(r8) :: subgrid
   !-------------------------------------------------------------------------------

   ni = 0._r8
   tc = tair - 273.15_r8

   ! initialize
   niimm = 0._r8
   nidep = 0._r8
   nihf  = 0._r8

   if(so4_num.ge.1.0e-10_r8 .and. (soot_num+dst_num).ge.1.0e-10_r8 .and. cldn.gt.0._r8) then

      !-----------------------------
      ! RHw parameterization for heterogeneous immersion nucleation
      A = 0.0073_r8
      B = 1.477_r8
      C = 131.74_r8
      RHw=(A*tc*tc+B*tc+C)*0.01_r8   ! RHi ~ 120-130%

      subgrid = 1.2_r8

      if((tc.le.-35.0_r8) .and. ((relhum*svp_water(tair)/svp_ice(tair)*subgrid).ge.1.2_r8)) then ! use higher RHi threshold

         A = -1.4938_r8 * log(soot_num+dst_num) + 12.884_r8
         B = -10.41_r8  * log(soot_num+dst_num) - 67.69_r8
         regm = A * log(wbar) + B

         if(tc.gt.regm) then    ! heterogeneous nucleation only
            if(tc.lt.-40._r8 .and. wbar.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation
               call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
               niimm=0._r8
               nidep=0._r8
               n1=nihf
            else
               call hetero(tc,wbar,soot_num+dst_num,niimm,nidep)
               nihf=0._r8
               n1=niimm+nidep
            endif
         elseif (tc.lt.regm-5._r8) then ! homogeneous nucleation only
            call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
            niimm=0._r8
            nidep=0._r8
            n1=nihf
         else        ! transition between homogeneous and heterogeneous: interpolate in-between
            if(tc.lt.-40._r8 .and. wbar.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation
               call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
               niimm=0._r8
               nidep=0._r8
               n1=nihf
            else

               call hf(regm-5._r8,wbar,relhum,subgrid,so4_num,nihf)
               call hetero(regm,wbar,soot_num+dst_num,niimm,nidep)

               if(nihf.le.(niimm+nidep)) then
                  n1=nihf
               else
                  n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5._r8)
               endif
            endif
         endif

         ni=n1

      endif
   endif

   ! deposition/condensation nucleation in mixed clouds (-40<T<0C) (Meyers, 1992)
   if(tc.lt.0._r8 .and. tc.gt.-37._r8 .and. qc.gt.1.e-12_r8) then
      esl = svp_water(tair)     ! over water in mixed clouds
      esi = svp_ice(tair)     ! over ice
      deles = (esl - esi)
      nimey=1.e-3_r8*exp(12.96_r8*deles/esi - 0.639_r8) 
   else
      nimey=0._r8
   endif

   nuci=ni+nimey
   if(nuci.gt.9999._r8.or.nuci.lt.0._r8) then
      write(iulog, *) 'Warning: incorrect ice nucleation number (nuci reset =0)'
      write(iulog, *) ni, tair, relhum, wbar, nihf, niimm, nidep,deles,esi,dst_num,so4_num
      nuci=0._r8
   endif

   nuci   = nuci*1.e+6_r8/rhoair    ! change unit from #/cm3 to #/kg
   onimey = nimey*1.e+6_r8/rhoair
   onidep = nidep*1.e+6_r8/rhoair
   oniimm = niimm*1.e+6_r8/rhoair
   onihf  = nihf*1.e+6_r8/rhoair

end subroutine nucleati

!===============================================================================

subroutine hetero(T,ww,Ns,Nis,Nid)

    real(r8), intent(in)  :: T, ww, Ns
    real(r8), intent(out) :: Nis, Nid

    real(r8) A11,A12,A21,A22,B11,B12,B21,B22
    real(r8) A,B,C

!---------------------------------------------------------------------
! parameters

      A11 = 0.0263_r8
      A12 = -0.0185_r8
      A21 = 2.758_r8
      A22 = 1.3221_r8
      B11 = -0.008_r8
      B12 = -0.0468_r8
      B21 = -0.2667_r8
      B22 = -1.4588_r8

!     ice from immersion nucleation (cm^-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)

      Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
      Nis = min(Nis,Ns)

      Nid = 0.0_r8    ! don't include deposition nucleation for cirrus clouds when T<-37C

end subroutine hetero

!===============================================================================

subroutine hf(T,ww,RH,subgrid,Na,Ni)

      real(r8), intent(in)  :: T, ww, RH, subgrid, Na
      real(r8), intent(out) :: Ni

      real(r8)    A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
      real(r8)    A2_fast,B2_fast
      real(r8)    C1_fast,C2_fast,k1_fast,k2_fast
      real(r8)    A1_slow,A2_slow,B1_slow,B2_slow,B3_slow
      real(r8)    C1_slow,C2_slow,k1_slow,k2_slow
      real(r8)    regm
      real(r8)    A,B,C
      real(r8)    RHw

!---------------------------------------------------------------------
! parameters

      A1_fast  =0.0231_r8
      A21_fast =-1.6387_r8  !(T>-64 deg)
      A22_fast =-6.045_r8   !(T<=-64 deg)
      B1_fast  =-0.008_r8
      B21_fast =-0.042_r8   !(T>-64 deg)
      B22_fast =-0.112_r8   !(T<=-64 deg)
      C1_fast  =0.0739_r8
      C2_fast  =1.2372_r8

      A1_slow  =-0.3949_r8
      A2_slow  =1.282_r8
      B1_slow  =-0.0156_r8
      B2_slow  =0.0111_r8
      B3_slow  =0.0217_r8
      C1_slow  =0.120_r8
      C2_slow  =2.312_r8

      Ni = 0.0_r8

!----------------------------
!RHw parameters
      A = 6.0e-4_r8*log(ww)+6.6e-3_r8
      B = 6.0e-2_r8*log(ww)+1.052_r8
      C = 1.68_r8  *log(ww)+129.35_r8
      RHw=(A*T*T+B*T+C)*0.01_r8

      if((T.le.-37.0_r8) .and. ((RH*subgrid).ge.RHw)) then

        regm = 6.07_r8*log(ww)-55.0_r8

        if(T.ge.regm) then    ! fast-growth regime

          if(T.gt.-64.0_r8) then
            A2_fast=A21_fast
            B2_fast=B21_fast
          else
            A2_fast=A22_fast
            B2_fast=B22_fast
          endif

          k1_fast = exp(A2_fast + B2_fast*T + C2_fast*log(ww))
          k2_fast = A1_fast+B1_fast*T+C1_fast*log(ww)

          Ni = k1_fast*Na**(k2_fast)
          Ni = min(Ni,Na)

        else       ! slow-growth regime

          k1_slow = exp(A2_slow + (B2_slow+B3_slow*log(ww))*T + C2_slow*log(ww))
          k2_slow = A1_slow+B1_slow*T+C1_slow*log(ww)

          Ni = k1_slow*Na**(k2_slow)
          Ni = min(Ni,Na)

        endif

      end if

end subroutine hf

!===============================================================================

end module nucleate_ice

