!------------------------------------------------------------------------------
!     'fjx_sub_mod.f90'  for Cloud-J v7.7     (2/2020, mjp)
!------------------------------------------------------------------------------
!
! DESCRIPTION:
!     v7.7  (02/2020) Final synch with Solar-J v7.7
!        Corrects problem with MAX-RAN that was caused by MAX-COR fixes
!        New calling sequence of FPs, added OD18
!     v7.6c (06/2019) Adds geometric option:
!        Corrects (optical) mass to the value for a true spherical atmosphere.
!        NOTE that inferred atmospheric mass for absorption & scattering will be
!           larger than is used in the std planar geopotential hydrostatic atmosphere.
!     (1) Z geometric replaces Z geopotential:  Z-geom = Z-geop / ( 1 - Z-geop/RAD)
!     (2) Layer (optical) mass increases by AMG(L) = (1 + Z-geom-midpt(L)/RAD)**2  > 1
!           to account for rho * dz going from geopotential to geometric
!     (3) Layer (optical) mass increases by a second factor of AMG(L) = (1 + Z/RAD)**2
!           but this is spread over a larger area (also AMG(L)) and so this factor
!           does NOT apply to optical mass in expanding area that used for ray-tracing.
!           It is recouped in (4) below.
!     (4) Flux depositon (absorbed or scattered) in each expanded layer is
!           increased by factor AMG(L) when used in multiple scattering to reflect
!           the smaller grid area at the surface assumed for that calculation.
!     (5) The optical mass grid used in spherical ray-tracing original must be
!           increased by the expansion factor AMG(L) for the scattering code so that
!           the total scattered/absorbed flux in an optically thin layer is included.
!
!     v7.6  (12/2018) CORRECTS the calc of deposition of direct beam (FLXD)
!         >>>> for conservative atmos, now predicts no spurious atmos absorption.
!         >>>> error in incident for clear skies:  <0.00% up to 80 sza, +0.04% at 88 sza
!         >>>> error for cloudy atmos depends on OD & extra layers (ATAU, ATAU0) in top of clouds
!              ATAU/ATAU0 = 1.10/.010  cld OD=38: +0.15% to +0.49% (sza = 0-80, 88) = 2x cost
!              ATAU/ATAU0 = 1.05/.005  cld OD=38: +0.05% to +0.16% (sza = 0-80, 88) = 3x cost
!              typical use before was 1.20/.020 which ads 50% more layers with OD=40 (1.5x cost)
!
!     v7.6  (07/2018) adds refraction to the SPHERE1R calculation.
!            A new version of SPHERE1N is also available (cleaner algorithm)
!            also a flat-Earth version SPHERE1F is available
!
!     v7.6  (06/2018) makes major changes in core scattering routine
!           to allow for angle-dependent albedos, specifically
!     The ALBEDO is now specified for each wavelength at the 4 quad angles
!          AND for the incident SZA (stored as the 5th albedo here).
!     New simpler way of interpolating TAU and F for inserted cloud layers
!     OD600 uses only clouds and aerosols, not O3 and Rayleigh as before
!     Dropped mid-layer odd-points (J's) to cut cost, interpolated cloud layers remain
!
!     v7.4  (08/2015) consistent with 7.1 data and results
!          variables in call to PHOTO_JX are same as in 7.1,
!          but a logical(out) LDARK is added to to count the number of J calcs
!          Works new ver 7.3 tfor cloud-J
!              + data sets for spectra, clouds and SS aerosols, new aerosol format
!          Extended to v7.4 to allow for Solar-J
!          v7.4d fixed the deposition of sunlight for SAZ>90, J's are unchanged.

      MODULE FJX_SUB_MOD

      USE FJX_CMN_MOD

!SJ! note that Solar-J uses these
!   USE CMN_FJX_MOD
!   USE CMN_H_MOD, only: iyear
!   USE RRSW_FASTJ_CMN, only: MXLAY, ngptsw



      IMPLICIT NONE

      PUBLIC :: PHOTO_JX, SOLAR_JX, JP_ATM0, ACLIM_FJX, ACLIM_GEO, ACLIM_RH, EXITC

      CONTAINS

!<<<<<<<<<<<<<<<<<<<<<<<<begin fastJX core code<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<all outside calls go through PHOTO_JX<<<<<<<<<<<<<<<<<<<<

!-----------------------------------------------------------------------
      subroutine PHOTO_JX (U0,SZA,RFL,SOLF, LPRTJ, PPP,ZZZ,TTT,HHH,    &
                       DDD,RRR,OOO, CCC, LWP,IWP,REFFL,REFFI,AERSP,    &
                       NDXAER, L1U,ANU,NJXU, VALJXX,SKPERD,SWMSQ,OD18, LDARK)
!  PHOTO_JX is the gateway to fast-JX calculations:
!    calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
!    needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
!    needs day-of-year for sun distance, SZA (not lat or long)
!-----------------------------------------------------------------------
      implicit none

!---calling sequence variables
      integer, intent(in)                    :: L1U,ANU,NJXU
      real*8,  intent(in)                    :: U0,SZA,SOLF
      real*8,  intent(in), dimension(5,W_+W_r) :: RFL
      logical, intent(in)                    :: LPRTJ
      real*8,  intent(in), dimension(L1U+1)  :: PPP,ZZZ
      real*8,  intent(in), dimension(L1U  )  :: TTT,HHH,DDD,RRR,OOO,CCC
      real*8,  intent(in), dimension(L1U  )  :: LWP,IWP,REFFL,REFFI
      real*8,  intent(in), dimension(L1U,ANU):: AERSP
      integer, intent(in), dimension(L1U,ANU):: NDXAER
! reports out the JX J-values, upper level program converts to CTM chemistry J's
      real*8,  intent(out), dimension(L1U-1,NJXU) ::  VALJXX
      real*8,  intent(out), dimension(S_+2,L1U)   ::  SKPERD
      real*8,  intent(out), dimension(L1U)        ::  OD18
      real*8,  intent(out), dimension(6)     :: SWMSQ
      logical, intent(out)                   :: LDARK

!-----------------------------------------------------------------------
!---key LOCAL atmospheric data needed to solve plane-parallel J & Heating
!-----these are dimensioned L1_
      real*8, dimension(L1_+1) :: PPJ,ZZJ
      real*8, dimension(L1_)   :: TTJ,HHJ,DDJ,RRJ,OOJ,CCJ
      integer,dimension(L1_)   :: JXTRA
!
      real*8, dimension(W_+W_r)       :: FJTOP,FJBOT,FSBOT,FLXD0
      real*8, dimension(5,W_+W_r)     :: FIBOT
      real*8, dimension(L1_,W_+W_r)   :: AVGF, FJFLX
      real*8, dimension(L1_,W_+W_r)   :: DTAUX, FLXD
      real*8, dimension(8,L1_,W_+W_r) :: POMEGAX
!
      real*8, dimension(L1_)        ::  DTAU600
      real*8, dimension(8,L1_)      ::  POMG600
      real*8, dimension(S_,L1_)     ::  FFX
      real*8, dimension(S_,8)       ::  FFXNET
      real*8, dimension(S_,4)       ::  FFXTAU
      real*8, dimension(8)          ::  FFXNETS
!--special Solar-J heating arrays for RRTMG, CLIRAD or LLNL
      real*8, dimension(mxlay, 0:ngptsw-1)  :: TAUG_RRTMG
      real*8, dimension(L1U, 0:30)          :: TAUG_CLIRAD
      real*8, dimension(L1U, 0:21)          :: TAUG_LLNL
!---flux/heating arrays (along with FJFLX,FLXD,FLXD0)
      real*8  FLXJ(L1_),FFX0,FXBOT,FABOT
      real*8  ODABS,ODRAY,ODI,ODL
      real*8  FREFS,FREFL,FREFI, FREF1,FREF2,PREF1,PREF2
      real*8  AMF(L1_+1,L1_+1), AMG(L1_)
!------------key SCATTERING arrays for clouds+aerosols------------------
      real*8  QQEXT(S_),SSALB(S_),SSLEG(8,S_),OPTX(S_),DDENS
      real*8  OD(S_,L1_),SSA(S_,L1_),SLEG(8,S_,L1_)
      real*8  OD600(L1_)
      real*8  PATH,RH,XTINCT,RE_LIQ,RE_ICE,TE_ICE
!------------key arrays AFTER solving for J's---------------------------
      real*8  FFF(W_,L1_),VALJ(X_)
      real*8  FLXUP(S_),FLXDN(S_),DIRUP(S_),DIRDN(S_)
      real*8  VALJL(L1_,X_) !2-D array of J_s returned by JRATET
!
      integer  LU,I,J,K,KR,KR0,KG,JG,L,M,NAER,RATIO(S_), NSS2
      real*8   XQO3,XQO2,TTTX, ODKL,DPKL
      real*8   ODRRTM,FRRTM
      real*8   ZMID
!-----------------------------------------------------------------------
      LU = L1U - 1
      VALJXX(:,:) = 0.d0
      FFXTAU(:,:) = 0.d0
      SKPERD(:,:)=0.d0
      SWMSQ(:)=0.d0
      OD18(:)=0.d0

!---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
!                        or         99.0                 80 km
      if (SZA .gt. 98.d0) then
         LDARK = .true.
         return
      else
         LDARK = .false.
      endif

!---load the atmospheric column data
      do L = 1,L1U
         PPJ(L) = PPP(L) !pressure in mb   bottom edge of layer L
         ZZJ(L) = ZZZ(L) !alt (cm)         ditto
         TTJ(L) = TTT(L) !Temp(K)
         DDJ(L) = DDD(L) !dry air
         OOJ(L) = OOO(L) !O3
         HHJ(L) = HHH(L) !h2o in molecules/cm2
         CCJ(L) = CCC(L) !methane
      enddo
      PPJ(L1U+1) = 0.d0
      ZZJ(L1U+1) = ZZZ(L1U) + ZZHT

      if (ATM0 .ge. 3) then
! v76c >>> correction for spherical geometric atmosphere
! Change altitude from geopotential to geometric /(1 - Z/R)
!     this will affect the ray-tracing geometry and air mass factors (AMF)
        do L = 2,L1U+1
           ZZJ(L) = ZZJ(L)/(1.d0 - ZZJ(L)/RAD)
        enddo
! Scale factor (1+Z/RAD)**2 = AMG(L) for area expansion & mass increase (geop dz)
        do L = 1,L1U
           ZMID = 0.5d0*(ZZJ(L) + ZZJ(L+1))
           AMG(L) = (1.d0 + ZMID/RAD)**2
        enddo
      else
           AMG(:) = 1.d0
      endif

!---calculate spherical weighting functions (AMF: Air Mass Factor)
!   indep of wavelength (refracted path assumes visible index of refraction)
!-----------------------------------------------------------------------
      if (ATM0 .eq. 0) then
         call SPHERE1F (U0,RAD,ZZJ,ZZHT,AMF, L1U)  ! flat Earth, AMF=1/u0
      elseif (ATM0 .eq. 1) then
         call SPHERE1N (U0,RAD,ZZJ,ZZHT,AMF, L1U)  ! spherical straight-line paths
      else     ! 2 or 3
         call SPHERE1R (U0,RAD,ZZJ,ZZHT,AMF, L1U)  ! spherical w/refraction
      endif
!-----------------------------------------------------------------------

      TAUG_RRTMG(:,:)= 0.d0
      TAUG_CLIRAD(:,:)=0.d0
      TAUG_LLNL(:,:) = 0.d0
!SJ! needed in Solar-J version
!      if (W_r .ne. 0)then
!         if (W_r .eq. W_rrtmg)  call RRTMG_SW_INP(IYEAR,L1U,PPJ,ZZJ,DDJ,TTJ,HHJ,OOJ,CCJ,TAUG_RRTMG)
!         if (W_r .eq. W_Clirad) call FJX_CLIRAD_H2O(L1U,PPJ,TTJ,HHJ,TAUG_CLIRAD)
!         if (W_r .eq. W_LLNL)   call FJX_GGLLNL_H2O(L1U,PPJ,TTJ,HHJ,TAUG_LLNL)
!         write(6,'(a,I5,a,I5,a,I5,a,I5)')'W_r=    ', W_r, 'W_rrtmg=', W_rrtmg, 'W_clirad=', W_clirad, 'W_LLNL=', W_LLNL

  !SJ    endif

      if (LPRTJ) then
         write(6,*)'Fast-J v7.6c ---PHOTO_JX internal print: fjx_sub_mod.f90'
         write(6,'(a,3i8)') ' NWBIN / NSBIN:',NWBIN,NSBIN
         write(6,'(a,3i8)') ' g-bin super-bin  L-flux: W_  S_  W_+W_r ',W_,S_,W_+W_r
         write(6,'(3i8)')   (KR, KDOKR(KR), LDOKR(KR), KR = 1,W_+W_r)
         write(6,*) '     L   P1      P2        T  //WPath, Reff, OD'
      endif

! >>>> major loop over standard levels:
      OD(:,:) = 0.d0
      SSA(:,:)= 0.d0
      SLEG(:,:,:)=0.d0
      do L = 1,L1U
         OD600(L) = 0.d0
! initialize scattering/absoprtion data with Rayleigh scattering (always non-zero)
! NB. SLEG(8,Kwavel,Llayer) includes the single-scattering albedo
         do K = 1,S_
            do I = 1,8
               SLEG(I,K,L) = 0.d0
            enddo
            ODRAY  = DDJ(L)*QRAYL(K)
            OD(K,L)  = ODRAY
            SSA(K,L) = ODRAY
            SLEG(1,K,L) = 1.0d0*ODRAY
            SLEG(3,K,L) = 0.5d0*ODRAY
         enddo
!>>>diagnostic print of Rayleigh data:
         if (LPRTJ) then
            write(6,'(a,i3,2f8.2,f8.1,1p,e12.4)') &
                 'Rayl',L,PPP(L),PPP(L+1),TTT(L),OD(18,L)
         endif
!---Liquid Water Cloud
         if (LWP(L) .gt. 1.d-5 .and. REFFL(L) .gt. 0.1d0) then
            RE_LIQ = REFFL(L)
            TE_ICE = TTT(L)
            call OPTICL (RE_LIQ,TE_ICE, DDENS, QQEXT,SSALB,SSLEG)
!---extinction K(m2/g) = 3/4 * Q / [Reff(micron) * density(g/cm3)]
            do K = 1,S_
               ODL = LWP(L) * 0.75d0 * QQEXT(K) / (RE_LIQ * DDENS)
               OD(K,L)  = OD(K,L)  + ODL
               SSA(K,L) = SSA(K,L) + SSALB(K)*ODL
               FFXTAU(K,3) = FFXTAU(K,3) + ODL*(1.d0-SSALB(K))
               FFXTAU(K,4) = FFXTAU(K,4) + ODL
               do I = 1,8
                  SLEG(I,K,L) = SLEG(I,K,L) + SSLEG(I,K)*SSALB(K)*ODL
               enddo
               if (K.eq.18) then
                  OD600(L) = OD600(L) + ODL
               endif
            enddo
!>>>diagnostic print of cloud data:
            if (LPRTJ) then
               write(6,'(a,i3,2f8.2,f8.1,f8.4,f8.2,f8.4)') &
                  'Rayl+Liq Cld',L,PPP(L),PPP(L+1),TTT(L),LWP(L),REFFL(L),OD(18,L)
            endif
         endif
!---Ice Water Cloud
         if (IWP(L) .gt. 1.d-5 .and. REFFI(L) .gt. 0.1d0) then
            RE_ICE = REFFI(L)
            TE_ICE = TTT(L)
            call OPTICI (RE_ICE,TE_ICE, DDENS, QQEXT,SSALB,SSLEG)
!---extinction K(m2/g) = 3/4 * Q / [Reff(micron) * density(g/cm3)]
            do K = 1,S_
               ODL = IWP(L) * 0.75d0 * QQEXT(K) / (RE_ICE * DDENS)
               OD(K,L)  = OD(K,L)  + ODL
               SSA(K,L) = SSA(K,L) + SSALB(K)*ODL
               FFXTAU(K,3) = FFXTAU(K,3) + ODL*(1.d0-SSALB(K))
               FFXTAU(K,4) = FFXTAU(K,4) + ODL
               do I = 1,8
                  SLEG(I,K,L) = SLEG(I,K,L) + SSLEG(I,K)*SSALB(K)*ODL
               enddo
               if (K.eq.18) then
                  OD600(L) = OD600(L) + ODL
               endif
            enddo
!>>>diagnostic print of cloud data:
            if (LPRTJ) then
               write(6,'(a,i3,2f8.2,f8.1,f8.4,f8.2,f8.4)') &
                  'Rayl+Liq+Ice Cld',L,PPP(L),PPP(L+1),TTT(L),IWP(L),REFFI(L),OD(18,L)
            endif
         endif
!---Strat Sulfate Aerosol Cloud: first aerosol index = 1 (bkgrd) or 2 (volcanic)
         do M = 1,ANU
            NAER = NDXAER(L,M)
            if ((NAER.eq.1) .or. (NAER.eq.2)) then
               PATH = AERSP(L,M)
               if (PATH .gt. 0.d0) then
                  call OPTICS (OPTX,SSALB,SSLEG, PATH,NAER)
                  do K = 1,S_
                     OD(K,L)  = OD(K,L)  + OPTX(K)
                     SSA(K,L) = SSA(K,L) + SSALB(K)*OPTX(K)
                     do I = 1,8
                        SLEG(I,K,L) = SLEG(I,K,L) + SSLEG(I,K)*SSALB(K)*OPTX(K)
                     enddo
                  enddo
                  OD600(L) = OD600(L) + OPTX(18)
!>>>diagnostic print of SSA data:
                  if (LPRTJ) then
                     write(6,'(a,i3,2f8.2,f8.1,2f8.5,1x,a12)') &
                       'StratSA',L,PPP(L),PPP(L+1),TTT(L),PATH,OPTX(18),TITLSS(NAER)
                  endif
               endif
            endif
         enddo
!---GEOMIP enhanced Strat Sulfate Aerosols:  index = 1001 to 1000+NGG
         do M = 1,ANU
            NAER = NDXAER(L,M)
            if (NAER .gt. 1000) then
               PATH = AERSP(L,M)
               if (PATH .gt. 0.d0) then
                  call OPTICG (OPTX,SSALB,SSLEG, PATH,NAER)
                  do K = 1,S_
                     OD(K,L)  = OD(K,L)  + OPTX(K)
                     SSA(K,L) = SSA(K,L) + SSALB(K)*OPTX(K)
                     do I = 1,8
                        SLEG(I,K,L) = SLEG(I,K,L) + SSLEG(I,K)*SSALB(K)*OPTX(K)
                     enddo
                  enddo
                  OD600(L)= OD600(L)+ OPTX(18)
!>>>diagnostic print of GEOMIP data:
                  if (LPRTJ) then
                     write(6,'(a,i3,2f8.2,f8.1,3f8.5,i10)') &
                        'SSA-GEO',L,PPP(L),PPP(L+1),TTT(L),PATH,OPTX(18),OD(18,L),NAER
                  endif
               endif
            endif
         enddo
!---OTHER aerosols in layer: check aerosol index
!---this uses data from climatology OR from current CTM (STT of aerosols)
!---subroutines OPTICA & OPTICM return the same information:
!---  PATH is the g/m2 in the layer, NAER in the cloud/aerosol index
!---  UMich aerosols use relative humidity (RH)
         RH = RRR(L)
         do M = 1,ANU
            NAER = NDXAER(L,M)
            PATH = AERSP(L,M)
            if (PATH .gt. 0.d0) then
               if (NAER.gt.2 .and. NAER.lt.1000) then
                  call OPTICA (OPTX,SSALB,SSLEG, PATH,RH, NAER)
                  do K = 1,S_
                     OD(K,L)  = OD(K,L)  + OPTX(K)
                     SSA(K,L) = SSA(K,L) + SSALB(K)*OPTX(K)
                     do I = 1,8
                        SLEG(I,K,L)=SLEG(I,K,L) + SSLEG(I,K)*SSALB(K)*OPTX(K)
                     enddo
                  enddo
                  OD600(L) = OD600(L) + OPTX(18)
!>>>diagnostic print of OPTICA data:
                  if (LPRTJ) then
                     write(6,'(a,i3,2f8.2,f8.1,2f8.5,i5)') &
                        'aerosol',L,PPP(L),PPP(L+1),TTT(L),PATH,OPTX(18),NAER
                  endif
               endif
            endif
         enddo
         do M = 1,ANU
            NAER = NDXAER(L,M)
            PATH = AERSP(L,M)
            if (PATH .gt. 0.d0) then
               if (NAER .lt. 0) then
                  call OPTICM (OPTX,SSALB,SSLEG, PATH,RH, -NAER)
                  do K = 1,S_
                     OD(K,L)  = OD(K,L)  + OPTX(K)
                     SSA(K,L) = SSA(K,L) + SSALB(K)*OPTX(K)
                     do I = 1,8
                        SLEG(I,K,L)=SLEG(I,K,L) + SSLEG(I,K)*SSALB(K)*OPTX(K)
                     enddo
                  enddo
                  OD600(L) = OD600(L) + OPTX(18)
!>>>diagnostic print of OPTICM data:
                  if (LPRTJ) then
                     write(6,'(a,i3,2f8.2,f8.1,2f8.5,i5)') &
                        'aerosol',L,PPP(L),PPP(L+1),TTT(L),PATH,OPTX(18),NAER
                  endif
               endif
            endif
         enddo

!---Add O2 & O3 absorbers to get final optical properties (Fast-J bins only 1:18)
         do K = 1,W_
            TTTX = TTJ(L)
            call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1),QO2(K,2),TQQ(3,1),QO2(K,3), LQQ(1))
            call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2),QO3(K,2),TQQ(3,2),QO3(K,3), LQQ(2))
            ODABS = XQO3*OOJ(L) + XQO2*DDJ(L)*0.20948d0
            OD(K,L)  = OD(K,L)  + ODABS
!SJ!       if (LPRTJ) then
!SJ!          write(6,'(a,i3,2f8.2,f8.1,f12.4)') &
!SJ!               'L/K/PPP/TTT/OD_O2/3',L,K,PPP(L),TTT(L),ODABS
!SJ!       endif
         enddo
!---renormalize the SLEG array by OD - note that SSA is included in SLEG and not used further
         do K = 1,S_
            do I = 1,8
               SLEG(I,K,L) = SLEG(I,K,L)/OD(K,L)
            enddo
            FFXTAU(K,1) = FFXTAU(K,1) + OD(K,L)*(1.d0 - SLEG(1,K,L))
            FFXTAU(K,2) = FFXTAU(K,2) + OD(K,L)
         enddo
      enddo   ! end of 'do L = 1,L1U'

! >>> now transform matrix OD(K,L) & SLEG (I,K,L) ==> DTAUX(L,KR) & POMEGAX(I,L,KR)
!     needed for good caching in the solver
! >>> also expand the K=1:S_ wavelengths of OD & SLEG to the KR=1:W_+W_r
!     for the full sub-bins of RRTMG.
      if (LRRTMG .or. LCLIRAD .or. LGGLLNL) then
         KR0 = W_
         KR = 0
         do K = 1,S_
            do J = 1,NGC(K)  ! for 1:17 just one gbin/bin, for 18:27 there are NGC gbins/bin
               KR = KR+1
               if (KR .lt. KR0) then
                  do L = 1,L1U
                     DTAUX(L,KR) = OD(K,L)
                     do I = 1,8
                        POMEGAX(I,L,KR) = SLEG(I,K,L)
                     enddo
                  enddo
               else
                  do L = 1,L1U
                     ODRRTM=0.D0
                     if (LCLIRAD) then
                        ODRRTM = TAUG_CLIRAD(L,KR-KR0)
                     elseif(LRRTMG)then
                        ODRRTM = TAUG_RRTMG(L,KR-KR0)
                     else if(LGGLLNL)then
                        ODRRTM = TAUG_LLNL(L,KR-KR0)
                     endif
                     DTAUX(L,KR) = OD(K,L) + ODRRTM
                     FRRTM = OD(K,L)/DTAUX(L,KR)
                     do I = 1,8
                        POMEGAX(I,L,KR) = SLEG(I,K,L)*FRRTM
                     enddo
                  enddo
               endif
            enddo
         enddo
      else   ! below = no gas abs bins 18:S_, but cloud/aers abs. can be Cloud-J
         do K=1,S_
            do L=1,L1U
               DTAUX(L,K)= OD(K,L)
               do I=1,8
                  POMEGAX(I,L,K)= SLEG(I,K,L)
               enddo
            enddo
         enddo
      endif
      do L=1, L1U
         OD18(L)= OD600(L)
      enddo


!---Using aerosol+cloud OD/layer in visible (600 nm) calculate how to add layers
!-----------------------------------------------------------------------
      call EXTRAL1(OD600,L1U,N_,ATAU,ATAU0, JXTRA)
!-----------------------------------------------------------------------
!---complete calculation of actinic and net fluxes for all L & wavelengths (incl W_+W_r)
!-----------------------------------------------------------------------
      call OPMIE (DTAUX,POMEGAX,U0,RFL,AMF,AMG,JXTRA, &
              AVGF,FJTOP,FJBOT,FIBOT,FSBOT,FJFLX,FLXD,FLXD0, LDOKR,LU)
!-----------------------------------------------------------------------
      FFF(:,:) = 0.d0
      FREFI = 0.d0
      FREFL = 0.d0
      FREFS = 0.d0
      FLXUP(:) = 0.d0
      DIRUP(:) = 0.d0
      FLXDN(:) = 0.d0
      DIRDN(:) = 0.d0
      FLXJ(:) = 0.d0
      FFX(:,:) = 0.d0
      FFXNET(:,:) = 0.d0
      FFXNETS(:) = 0.d0
      FREF1 = 0.d0
      FREF2 = 0.d0
      PREF1 = 0.d0
      PREF2 = 0.d0
! accumulate data on solar fluxes:  actinic = J-values (1:W_), users FL(K) [photons]
! note FFF(18,L) uses only the first sub-bin of K=18 (97% of flux) and applies 100% of FL(K)
      do K = 1,W_
         if (LDOKR(K) .gt. 0) then
            do L = 1,LU
               FFF(K,L) = SOLF*FL(K)*AVGF(L,K)
            enddo
            PREF1 = PREF1 + FSBOT(K)*SOLF*FL(K)*FP(K)  ! PAR direct
            PREF2 = PREF2 + FJBOT(K)*SOLF*FL(K)*FP(K)  ! PAR diffuse
         endif
      enddo
!---use the FFF() values  in photons/cm2/sec to calculate J's
!---mapping J-values from fast-JX species onto CTM chemistry reactins is done in main code
!-----------------------------------------------------------------------
      call JRATET(PPJ,TTJ,FFF, VALJXX, LU,NJXU)
!-----------------------------------------------------------------------
! accumulate data on solar fluxes:  energy and solar heating (!:S_), uses FW(K) [Watts]
      KG = 0
      do K = 1,S_
         do JG = 1,NSJSUB(K) ! NSJSUB could be NGC or ones set determined at INIT
            KG = KG+1
            if (LDOKR(KG) .gt. 0) then
!  direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
!  also at bottom (DN), does not include diffuse reflected flux.
               FLXUP(K) = FLXUP(K) + FJTOP(KG)*SJSUB(K,JG)
               DIRUP(K) = DIRUP(K) - FLXD0(KG)*SJSUB(K,JG)
               FLXDN(K) = FLXDN(K) - FJBOT(KG)*SJSUB(K,JG)
               DIRDN(K) = DIRDN(K) - FSBOT(KG)*SJSUB(K,JG)
               FREFI = FREFI + FLXD0(KG)*SOLF*FW(K)*SJSUB(K,JG)
               FREFL = FREFL + FJTOP(KG)*SOLF*FW(K)*SJSUB(K,JG)
               FREFS = FREFS + SOLF*FW(K)*SJSUB(K,JG)
               FABOT = (FJBOT(KG) + FSBOT(KG)) - FIBOT(5,KG)
               FXBOT = FSBOT(KG) - FABOT
               FLXJ(1) = FJFLX(1,KG) - FXBOT
               do L = 2,LU
                  FLXJ(L) = FJFLX(L,KG) - FJFLX(L-1,KG)
               enddo
               FLXJ(LU+1) = FJTOP(KG) - FJFLX(LU,KG)
               FFX0 = 0.d0
               do L = 1,L1U
                  FFX0 =     FFX0     + (FLXD(L,KG) - FLXJ(L))*SJSUB(K,JG)
                  FFX(K,L) = FFX(K,L) + (FLXD(L,KG) - FLXJ(L))*SJSUB(K,JG)
               enddo
               FFXNET(K,1) = FFXNET(K,1) + (FLXD0(KG)             )*SJSUB(K,JG)  ! direct(solar) flux dep into atmos (spherical)
               FFXNET(K,2) = FFXNET(K,2) + (FSBOT(KG)             )*SJSUB(K,JG)  ! direct(solar) flux dep onto LB (surface)
               FFXNET(K,3) = FFXNET(K,3) + (FLXD0(KG)  + FSBOT(KG))*SJSUB(K,JG)  ! total solar into atmopshere+surface
               FFXNET(K,4) = FFXNET(K,4) + (FJTOP(KG)             )*SJSUB(K,JG)  ! diffuse flux leaving top-of-atmos
               FFXNET(K,5) = FFXNET(K,5) + (FFX0                  )              ! diffuse flux absorbed in atmos
               FFXNET(K,6) = FFXNET(K,6) + (FABOT                 )*SJSUB(K,JG)  ! total (dir+dif) absorbed at LB (surface)
               FFXNET(K,7) = FFXNET(K,7) + (FSBOT(KG)             )*SJSUB(K,JG)  ! direct flux dep onto LB (surface diags)
               FFXNET(K,8) = FFXNET(K,8) + (FJBOT(KG)             )*SJSUB(K,JG)  ! diffuse flux dep onto LB (surface)
            endif
         enddo ! end JG/KG loop over g-bins embedded in the S_=27 super bins
      enddo  ! end loop over wavelength super bins K
!-----------------------------------------------------------
      FREFL = FREFL/FREFS     !calculate fraction reflected flux (energy weighted)
      FREFI = FREFI/FREFS
! calc K/day & other fluxes
      do L = 1,L1U
         SKPERD(S_+1,L) = 0.d0
         SKPERD(S_+2,L) = 0.d0
         DPKL = HeatFac_/(PPP(L)-PPP(L+1))
         do K = 1,S_
            SKPERD(K,L) = FFX(K,L)*FW(K)*SOLF*DPKL
         enddo
         do K = 1,W_
            SKPERD(S_+1,L) = SKPERD(S_+1,L) + SKPERD(K,L)
         enddo
         do K = W_+1,S_
            SKPERD(S_+2,L) = SKPERD(S_+2,L) + SKPERD(K,L)
         enddo
      enddo
      do J = 1,8
         do K = 1,S_
            FFXNETS(J) = FFXNETS(J) + FFXNET(K,J)*SOLF*FW(K)
         enddo
      enddo
      SWMSQ(:) = 0.d0
      do K=1,S_
         SWMSQ(1) = SWMSQ(1) + FFXNET(K,3)*SOLF*FW(K)
         SWMSQ(2) = SWMSQ(2) + FFXNET(K,4)*SOLF*FW(K)
         SWMSQ(3) = SWMSQ(3) + FFXNET(K,5)*SOLF*FW(K)
         SWMSQ(4) = SWMSQ(4) + FFXNET(K,6)*SOLF*FW(K)
      enddo
      do K=1, W_
         SWMSQ(5) = SWMSQ(5) + FFXNET(K,7)*FL(K)*FP(K)*SOLF
         SWMSQ(6) = SWMSQ(6) + FFXNET(K,8)*FL(K)*FP(K)*SOLF
      enddo

!---diagnostics/variables below are JUST for PRINT and NOT returned to the CTM code
      if (LPRTJ) then
         do L=1,L1U
            DTAU600(L) = DTAUX(L,18)
            do I=1,8
               POMG600(I,L) = POMEGAX(I,L,18)
            enddo
         enddo
         write(6,'(a)') 'Fast-J  v7.6 ---PHOTO_JX internal print: Atmosphere--'
         call JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU600,POMG600,JXTRA, LU)
!SJ!         if (LRRTMG .or. LCLIRAD .or. LGGLLNL) then
!SJ!            write(ParaSummary(26:200),'(a, 10f10.4)') &
!SJ!               ' RFL(,18)/SZA/u0/maxOD600/F-incd/F-refl/: ', &
!SJ!               (RFL(i,18),i=1,5),SZA, U0, MAXVAL(OD600), FREFI, FREFL
!SJ!            write(6,'(a)') ParaSummary(1:200)
!SJ!         endif
         write(6,'(a)') 'Fast-J  v7.6 ---PHOTO_JX internal print: Solar fluxes (W/m2)--'
         write(6,'(a11,f12.4)')    ' inc TOTAL ',SWMSQ(1)
         write(6,'(a11,f12.4)')    ' rfl outtop',SWMSQ(2)
         write(6,'(a11,f12.4)')    ' abs in atm',SWMSQ(3)
         write(6,'(a11,f12.4)')    ' abs at srf',SWMSQ(4)
         write(6,'(a11,1p,e12.4)') ' PAR direct',SWMSQ(5)
         write(6,'(a11,1p,e12.4)') ' PAR diffus',SWMSQ(6)
         NW1=1
         NS1=1
         NW2=W_
         NS2=S_
         write(6,'(a)') 'Spectral budget scaled to 1.0, wavelengths 27:-1:1'
         write(6,'(a11/(9f8.4))') ' inc TOTAL ',(FFXNET(K,3), K=NS2,NS1,-1)
         write(6,'(a11/(9f8.4))') ' rfl outtop',(FFXNET(K,4), K=NS2,NS1,-1)
         write(6,'(a11/(9f8.4))') ' abs in atm',(FFXNET(K,5), K=NS2,NS1,-1)
         write(6,'(a11/(9f8.4))') ' abs at srf',(FFXNET(K,6), K=NS2,NS1,-1)
         do J = 1,8
            do K = 1,S_
               FFXNET(K,J) = FFXNET(K,J)*SOLF*FW(K)
            enddo
         enddo
         do K = 1,W_
            FREF1 = FREF1 + FFXNET(K,3)
         enddo
         do K = W_+1,S_
            FREF2 = FREF2 + FFXNET(K,3)
         enddo
         if (.not.(LRRTMG .or. LCLIRAD)) then
            write(6,'(a,f10.2,f10.6,3f10.4)') 'FJX: SZA/u0/F-incd/F-refl/', &
               SZA,U0,FREFI,FREFL,FREFS
            write(6,'(a,5f10.4)') 'FJX80: albedos@600nm u(1:4) & u0(5)', &
               RFL(1:5,18)
            write(6,'(a5,20i8)')   ' bin:',(K, K=NW1,NW2)
            write(6,'(a5,20f8.1)') ' wvl:',(WL(K), K=NW1,NW2)
            write(6,'(a)') ' ---- 100000=Fsolar   MEAN INTENSITY per wvl bin'
                 RATIO(:) = 0.d0
            do L = LU,1,-1
               do K=NW1,NW2
                  if (LDOKR(K) .gt. 0) then
                     RATIO(K) = (1.d5*FFF(K,L)/(SOLF*FL(K)))
                  endif
               enddo
               write(6,'(i3,2x,18i8)') L,(RATIO(K),K=NW1,NW2)
            enddo
            write(6,'(a)') 'specific intensity onto surface @ 4 quad angles'
            do I = 1,4
               write(6,'(f5.2,20f8.5)') EMU(I),(FIBOT(I,K), K=NW1,NW2)
            enddo
            write(6,'(a)') 'specific intensity up from lower bndry'
            write(6,'(f5.2,20f8.5)') U0,(FIBOT(5,K), K=NW1,NW2)
            write(6,'(a)') 'direct flux at lower bndry'
            write(6,'(f5.2,20f8.5)') U0,(FSBOT(K), K=NW1,NW2)
            write(6,'(a)') 'diffus flux at lower bndry'
            write(6,'(f5.2,20f8.5)') U0,(FJBOT(K), K=NW1,NW2)
            NSS2 = min(S_,NW2+32)  !  limt for large # of super bins: high-res clouds
            write(6,'(a5,32i8)')   ' bin:',(K, K=NW2+1,NSS2)
            write(6,'(a5,32f8.1)') ' wvl:',(WL(K), K=NW2+1,NSS2)
            write(6,'(a)') 'specific intensity to surface at the 4 quad angles'
            do I = 1,4
               write(6,'(f5.2,32f8.5)') EMU(I),(FIBOT(I,K), K=NW2+1,NSS2)
            enddo
            write(6,'(a)') 'specific intensity up from lower bndry'
            write(6,'(f5.2,32f8.5)') U0,(FIBOT(5,K), K=NW2+1,NSS2)
            write(6,'(a)') 'direct flux at lower bndry'
            write(6,'(f5.2,32f8.5)') U0,(FSBOT(K), K=NW2+1,NSS2)
            write(6,'(a)') 'diffus flux at lower bndry'
            write(6,'(f5.2,32f8.5)') U0,(FJBOT(K), K=NW2+1,NSS2)
            write(6,*)
            write(6,*)'Fast-J v7.6 ---PHOTO_JX Net Fluxes include SZA & solar dist'
            write(6,'(a,2f8.2)') ' ---NET FLUXES--- solar < 700 or 778 nm ',FREF1,FREF1+FREF2
            write(6,'(a11,18i8)')   'bins:',(K, K=NW1,NW2)
            write(6,'(a11,18f8.1)') 'wavl:',(WL(K), K=NW1,NW2)
            write(6,'(a11,18f8.2)') 'watt:',(FW(K), K=NW1,NW2)
            write(6,'(a11,18f8.2)') ' sol atm+sf',(FFXNET(K,3), K=NW1,NW2)
            write(6,'(a11,18f8.2)') ' sol in atm',(FFXNET(K,1), K=NW1,NW2)
            write(6,'(a11,18f8.2)') ' sol at srf',(FFXNET(K,2), K=NW1,NW2)
            write(6,'(a11,18f8.2)') ' dif outtop',(FFXNET(K,4), K=NW1,NW2)
            write(6,'(a11,18f8.2)') ' abs in atm',(FFXNET(K,5), K=NW1,NW2)
            write(6,'(a11,18f8.2)') ' abs at srf',(FFXNET(K,6), K=NW1,NW2)
            write(6,'(a11,18f8.2)') ' srf direct',(FFXNET(K,7), K=NW1,NW2)
            write(6,'(a11,18f8.2)') ' srf diffus',(FFXNET(K,8), K=NW1,NW2)
            write(6,'(a11,11f8.2,7f8.3)') ' tau absorb',(FFXTAU(K,1), K=NW1,NW2)
            write(6,'(a11,11f8.2,7f8.3)') ' tau total ',(FFXTAU(K,2), K=NW1,NW2)
            write(6,'(a11,11f8.2,7f8.3)') ' cld absorb',(FFXTAU(K,3), K=NW1,NW2)
            write(6,'(a11,11f8.2,7f8.3)') ' cld total ',(FFXTAU(K,4), K=NW1,NW2)
            write(6,'(a11,1p,e10.3)') ' PAR direct',PREF1
            write(6,'(a11,1p,e10.3)') ' PAR diffus',PREF2
            write(6,'(a,3f8.2)') ' ---NET FLUXES--- solar <700, >700, total: ',&
                             FREF1,FREF2,FREF1+FREF2
            write(6,'(a11,32i8)')   'bins:',(K, K=NW2+1,NSS2)
            write(6,'(a11,32f8.1)') 'wavl:',(WL(K), K=NW2+1,NSS2)
            write(6,'(a11,32f8.2)') 'watt:',(FW(K), K=NW2+1,NSS2)
            write(6,'(a11,33f8.2)') ' sol atm+sf',(FFXNET(K,3), K=NW2+1,NSS2),FFXNETS(3)
            write(6,'(a11,33f8.2)') ' sol in atm',(FFXNET(K,1), K=NW2+1,NSS2),FFXNETS(1)
            write(6,'(a11,33f8.2)') ' sol at srf',(FFXNET(K,2), K=NW2+1,NSS2),FFXNETS(2)
            write(6,'(a11,33f8.2)') ' dif outtop',(FFXNET(K,4), K=NW2+1,NSS2),FFXNETS(4)
            write(6,'(a11,33f8.2)') ' abs in atm',(FFXNET(K,5), K=NW2+1,NSS2),FFXNETS(5)
            write(6,'(a11,33f8.2)') ' abs at srf',(FFXNET(K,6), K=NW2+1,NSS2),FFXNETS(6)
            write(6,'(a11,33f8.2)') ' srf direct',(FFXNET(K,7), K=NW2+1,NSS2),FFXNETS(7)
            write(6,'(a11,33f8.2)') ' srf diffus',(FFXNET(K,8), K=NW2+1,NSS2),FFXNETS(8)
            write(6,'(a11,32f8.3)') ' tau absorb',(FFXTAU(K,1), K=NW2+1,NSS2)
            write(6,'(a11,32f8.3)') ' tau total ',(FFXTAU(K,2), K=NW2+1,NSS2)
            write(6,'(a11,32f8.3)') ' cld absorb',(FFXTAU(K,3), K=NW2+1,NSS2)
            write(6,'(a11,32f8.3)') ' cld total ',(FFXTAU(K,4), K=NW2+1,NSS2)
            write(6,'(a)') 'heating rate profiles in K/day v7.6  180-778nm '
            write(6, '(a4, 32f7.1)')'wvl ',(WL(I),I=NW1,NW2)
            do L = LU,1,-1
               write(6,'(i4,32f7.2)') L,(SKPERD(I,L), I=NW1,NW2)
            enddo
            write(6,'(a)') 'heating rate profiles in K/day v7.6 778-...nm plus 1:18 19:27 1:27'
            write(6, '(a4,32f7.1)')'wvl ',(WL(I),I=NW2+1,NSS2)
            do L = LU,1,-1
               write(6,'(i4,35f7.2)') L,(SKPERD(I,L), I=NW2+1,NSS2+2),SKPERD(S_+1,L)+SKPERD(S_+2,L)
            enddo
            write(6,'(a)') ' Fast-J  v7.6 ----J-values----'
            write(6,'(1x,a,72(a6,3x))') 'L=  ',(TITLEJX(K), K=1,NJX)
            do L = LU,1,-1
               write(6,'(i3,1p, 72e9.2)') L,(VALJXX(L,K),K=1,NJX)
            enddo
       endif  ! end of CLOUDJ print
      endif   ! end of LPRTJ if

      END SUBROUTINE PHOTO_JX


!<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------
      subroutine OPMIE (DTAUX,POMEGAX,U0,RFL,AMF,AMG,JXTRA, &
              FJACT,FJTOP,FJBOT,FIBOT,FSBOT,FJFLX,FLXD,FLXD0, LDOKR,LU)
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in)  ::  DTAUX(L1_,W_+W_r),POMEGAX(M2_,L1_,W_+W_r)
      real*8, intent(in)  ::  AMF(L1_+1,L1_+1),AMG(L1_)
      real*8, intent(in)  ::  U0, RFL(5,W_+W_r)
      integer, intent(in) ::  JXTRA(L1_), LDOKR(W_+W_r),LU
      real*8, intent(out) ::  FJACT(L1_,W_+W_r),FIBOT(5,W_+W_r)
      real*8, intent(out) ::  FJTOP(W_+W_r),FJBOT(W_+W_r),FSBOT(W_+W_r)
      real*8, intent(out) ::  FJFLX(L1_,W_+W_r),FLXD(L1_,W_+W_r),FLXD0(W_+W_r)

      integer JADDTO,L2LEV(L1_+1)
      integer I,II,J,K,L,LL,LL0,IX,JK,   L2,L22,LZ,LZZ,ND
      integer L1U,  LZ0,LZ1,LZMID
      real*8   SUMT,SUMJ,DIVT, FBTMLOG

      real*8  TTAU(L1_+1)
      real*8  ATAUA,ATAUZ,XLTAU,XLTAU1,FJFLX0
      real*8  TAUBTM,TAUTOP,FBTM,FTOP,POMTOP(M2_),POMBTM(M2_)
      real*8  DTAU1(L1_+1), FTAU(L1_+1)
      real*8  DTAUEXT,DTAUSCA,DTAUABS,DECAY
      real*8  POMEGA1(M2_,L1_+1)
!--- variables used in mie code-----------------------------------------
      real*8, dimension(W_+W_r)           :: ZFLUX
      real*8, dimension(N_,W_+W_r)        :: FJ,FZ,ZTAU
      real*8, dimension(M2_,N_,W_+W_r)    :: POMEGA
!
! fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts

! in:
!     DTAUX(1:L1_,1:W_+W_r) = optical depth of each layer
!     POMEGAX(1:8,1:L1_,1:W_+W_r) = scattering phase fn (multiplied by s-s abledo)
!     U0  = cos (SZA)
!  ** RFL(5,1:W_) = Lambertian albedo of surface for angles 1:4 & U0 (#5)
!     AMF(1:L1_+1,1:L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
!        AMF now back to NOT inserting layers
!     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
! out:
!     FJACT(1:L_,1:W_) = mean actinic flux(diff+direct) at std CTM levels(mid-lyr)
!  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
!     FJTOP(1:W_) = diffuse flux out top-of-atmosphere (TAU=0 above top model lyr)
!     FJBOT(1:W_) = diffuse flux onto surface (<0 by definition)
!     FSBOT(1:W_) = direct/solar flux onto surface  (<0 by definition)
!     FIBOT(5,1:W_) = mean intensities onto surface, #5 = I-plus refelected from surface
!     FJFLX(1:L_,1:W_) = diffuse flux across top of model layer L
!        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
!     FLXD(1:L_+1,1:W_) = solar flux deposited in layer L (includes lyr above CTM)
!        this should take into account sphericity, and is not just = mu0
!     FLXD0(1:W_) = sum of solar flux deposited in atmos
!        does NOT include flux on lower surface, does NOT mean absorbed!
!-----------------------------------------------------------------------
!
!     DTAU     Local optical depth of each CTM level
!     TTAU     Optical depth of air vertically above each point (to 0 at top of atm)
!     FTAU     Attenuation of solar beam
!     POMEGAJ  Scattering phase function
!
!---------------------SET UP FOR MIE CODE-------------------------------
!
!-----------------wavelength independent--------------------------------
!
!  Transpose the ascending TTAU grid to a descending ZTAU grid.
!  Double the resolution - TTAU points become the odd points on the
!  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
!  Odd point added at top of grid for unattenuated beam   (Z='inf')
!
!  The following mapping holds for JXTRA=0
!        Surface:   TTAU(1)    ==> ZTAU(2*L1_+1)
!        Top:       TTAU(L1_)  ==> ZTAU(3)
!        Infinity:     0.0     ==> ZTAU(1)
!        index: 2*L1_ - 2*L + 3 ==> LZ
!
!------------------------------------------------------------------------
!  Insert new levels, working downwards from the top of the atmosphere
!  to the surface (increasing in LZ, decreasing in 'L')
!------------------------------------------------------------------------
!  To include the odd/even finite difference eqns for 'j' and 'h' alternating,
!  there are twice the number of layers in the LZ arrays (2*L1_ + 2*JADDTO)
!    and an extra point in the RT solution (at the edges) = 2*L1_+2*JADDTO+1
!
!  must tansfer L=1:L1_ (TTAU,FTAU,POMEGAJ) values onto the reverse
!    order, expanded, doubled-level scatter grid.
!  Note that we need to include the expansion by JXTRA levels (L2L).
!
!----------------------re-grid data---------------------------------------------
!  Calculate cumulative total and define levels we want J-values at.
!  Sum upwards for levels, and then downwards for Mie code readjustments.
!
!     JXTRA(L)  Number of new levels to add between edge (L) and edge (L+1)
!           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
!     JADDTO is the cumulative number of JXTRA levels to be added
!---these should be fixed for all wavelengths to lock-in the array sizes

      L1U = LU + 1   ! uppermost edge is LU+2
        JADDTO = 0  ! no added layers in the topmost added layer (goes to TAU=0)
      do L = 1,L1U
        JADDTO = JADDTO + JXTRA(L)
      enddo
      ND = 2*L1U + 2*JADDTO + 1
      if(ND .gt. N_) then
        call EXITC (' overflow of scatter arrays: ND > N_')
      endif
!---L2LEV(L) = L-index for old layer-edge L in the expanded JXTRA-grid
!     in absence of JXTRA,  L2LEV(L) = L
        L2LEV(1)  = 1
      do L = 2,L1U+1
        L2LEV(L) = L2LEV(L-1) + 1 + JXTRA(L-1)
      enddo

!----------------begin wavelength dependent set up------------------------------
!---Reinitialize arrays
        ZTAU(:,:)     = 0.d0
        POMEGA(:,:,:) = 0.d0
        FJACT(:,:) = 0.d0
        FJTOP(:) = 0.d0
        FJBOT(:) = 0.d0
        FSBOT(:) = 0.d0
        FJFLX(:,:) = 0.d0
        FLXD(:,:) = 0.d0
        FLXD0(:) = 0.d0
        FJ(:,:) = 0.d0
        FZ(:,:) = 0.d0
!---PRIMARY loop over wavelengths
      do K = 1,W_+W_r
      if (LDOKR(K) .gt. 0) then
!---FTAU(L) = solar beam (attenuated) at layer edge L, thus need L1U+1 for uppermost layer.
!   FTAU(1:L1U+1) = 1.0 when u0>0, but for u0<0 must be decayed from FTAU(L1U) to FTAU(L1U+1)
!   To calculate FTAU(L1U+1) we need to have an artificial layer above (DTAU1(L1U+1) = 0)
!        and a finite air mass AMF(L1U+1,L1U+1) = 1.0.  This works fine as long as
!        AMF(L<L1U+1,L1U+1) =0 then U0>0 and FTAU(L1U+1+1) = 1.0, but if
!        AMF(L<L1U+1,L1U+1) >0 then U0<0 and sum XLTAU for L<L1U+1, get >0, & FTAU(L1U+1+1) <0
!---  FTAU and DTAU1 are local, non-K

! new v7.6c >>> correction for spherical geometric atmosphere
!     DTAU1 is increased by AMG = (1+z/R)**2 for solar ray because of geopotential dz
!     if not 'geom' model then AMG = 1
! >>>> DTAU1 is used ONLY to calculate solar ray attenuation
       do L = 1,L1U
          DTAU1(L) = DTAUX(L,K) * AMG(L)
       enddo
          DTAU1(L1U+1) = 0.d0
          LL0 = 0  ! LL0 is the lowest AMF(L,L) = 0 = shadow ht
       do LL = 1,L1U+1
          if (AMF(LL,LL) .le. 0.0d0) then
             LL0 = LL
          endif
       enddo
          FTAU(:) = 0.d0
       do LL = LL0+1,L1U+1
! there is sunlight thru layer LL to layer edge/radius LL (=1:L1U)
! AMF(I,L) includes air mass effective (AMF ~ 1/U0) for all layers I at edge L
! when U0<0 layers for I<L are included (w/double weights).and AMF(L1U,L1U+1) > 0
             XLTAU = 0.0d0
          do II = 1,L1U
             XLTAU = XLTAU + DTAU1(II)*AMF(II,LL)
          enddo
!---FTAU(LL) = solar beam at layer-edge LL which lies below layer with DTAU1(LL)
          if (XLTAU .lt. 82.d0) then
             FTAU(LL) = exp(-XLTAU)
          endif
       enddo
          FSBOT(K) = 0.d0
       if (LL0 .eq. 0) then      ! direct beam on the ground, u0 ~ 1/AMF
          FSBOT(K) = FTAU(1)/AMF(1,1)
       endif

!---calculate column optical depths above each edge L, TTAU(1:L1_+1)
!---v76c optical depth for 1D scattering needs to increase DTAUX(L,K)
!      by AMG**2: 1 for the geop dz (as above for spherical shells)
!                 2 for the squeezing expanded area into the surface column.
!>>>>> Note TTAU is column optical depth used in the multiple scattering code
!>>>>>     DTAUX(L,K) is the layer optical depth from the std planar geopot model
          TTAU(L1U+1) = 0.0d0
       do L = L1U,1,-1
          TTAU(L) = TTAU(L+1) + DTAUX(L,K) * AMG(L)**2
       enddo
!---calculate scattering phase fn (POMEGA(1:8,:,:) at the edge points
!     POMEGAX is a layer average and must be interpolated to POMEGA1
!     Lower boundary value POMEGAX(I,L) is unchanged
       do I = 1,M2_
          POMEGA1(I,1)     = POMEGAX(I,1,K)
          POMEGA1(I,L1U+1) = POMEGAX(I,L1U,K)
       enddo
       do L = 2,L1U
          do I = 1,M2_
             POMEGA1(I,L) = (POMEGAX(I,L,K)*DTAUX(L,K) + &
                  POMEGAX(I,L-1,K)*DTAUX(L-1,K)) / (DTAUX(L,K)+DTAUX(L-1,K))
          enddo
       enddo

!---prepare for interpolating semi-log TAU-grid at top of clouds, L2LEV(L) = added levels
!---prepare master solution grid with both even (j) and odd (h) levels
!---master solution grid (N_) is reversed with TAU=0 at top of atmosphere
!---  LZ is reverse of L with added mid-point (even LZ)

!---Move everything onto the LZ arrays, interpolate the added j-levels (odd pts)
       do L = 1,L1U+1          ! L = index of CTM edge- and mid-layers
          L2 = L2LEV(L)        ! L2 = i# of interp layers
          LZ  = ND + 2 - 2*L2  ! LZ = index for master solution arrays (LZ=1 = top)
          ZTAU(LZ,K) = TTAU(L)
          FZ(LZ,K)   = FTAU(L)
          do I=1,M2_
             POMEGA(I,LZ,K) = POMEGA1(I,L)
          enddo
       enddo

!---Now go thru the pairs of L levels to see if we need JADD levels
       do L = 1,L1U              ! L = index of CTM edge- and mid-layers
          L2 = L2LEV(L)          ! L2 = index for L2 in expanded scale(JXTRA)
          LZ  = ND +   2 - 2*L2  ! LZ = index for scatt arrays
          L22 = L2LEV(L+1) - L2LEV(L) - 1   ! L22 = 0 if no added levels
          if (L22 .gt. 0) then
             TAUBTM = TTAU(L)
             TAUTOP = TTAU(L+1)
             FBTM   = FTAU(L)
             FTOP   = FTAU(L+1)
             do I=1,M2_
                POMBTM(I) = POMEGA1(I,L)
                POMTOP(I) = POMEGA1(I,L+1)
             enddo
! ATAU = 1.05 / 0.005 (most accurate version) each successive delta_TAU increase 5%
                DIVT = 1.d0/(ATAU**(L22+1) - 1.d0)
             do LL = 1,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
                LZZ = LZ - 2*LL      ! LZZ = index(odd) of added level in scatt arrays
!--- for L=1 this is the odd (j) layer just above TAUBTM @ L2 on other grid
                SUMT = (ATAU**(L22+1-LL) - 1.d0)*DIVT
                ZTAU(LZZ,K) = TAUTOP + SUMT*(TAUBTM-TAUTOP)
                if (AMF(1,1) .gt. 0.d0) then
!--- sun is up, top-lit, flux decays downward  - reduce TAU for solar interpolation
! v7.6c ZTAU & TAUTOP/TAUBTM have been increased by AMG**2 for 1D scattering,
!     reduce this for direct beam in a spherical atmosphere
                   DTAUEXT = (ZTAU(LZZ,K)-TAUTOP) / AMG(L)
                   FZ(LZZ,K) =   FTOP*exp(-AMF(L,L)*DTAUEXT)
                else
!--- surface dark, bottom-lit, flux decays upward
                   DTAUEXT = (TAUBTM - ZTAU(LZZ,K)) / AMG(L)
                   FZ(LZZ,K) =   FBTM*exp(-AMF(L,L)*DTAUEXT)
                endif
                do I = 1,M2_
                   POMEGA(I,LZZ,K) = POMTOP(I) + SUMT*(POMBTM(I)-POMTOP(I))
                enddo
             enddo
          endif
       enddo

!---Diagnose direct solar-beam flux deposited in each standard L=1:L1U layer
!---    this is critical as it determines the atmospheric heating rates.

!---Need go thru all the odd (j) LZ layers and sum into the L layers where above shaddow ht
       do L = LL0+1,L1U              ! L = index of CTM edge- and mid-layers
          L2 = L2LEV(L)          ! L2 = index for L2 in expanded scale(JXTRA)
          LZ  = ND +   2 - 2*L2  ! LZ = index for scatt arrays
          L22 = L2LEV(L+1) - L2LEV(L) - 1   ! L22 = 0 if no added levels
        if (L22 .eq. 0) then     !--- Standard L layers, no interpoalted ones!
!--- flux deposited in each layer L(bottom) to L+1(top) must be calculated for both
!    for conservative scattering extinction and separately for absorption extinction
!    For scattering, the correct flux deposited is Dtau*(Ftop + Fbot)/2
!    For absorption, the correct flux needs to be the amount attenuated (down and up)
!      Note use of POMEGAX(1:8,L,K) = original mean layer value L to L+1   (not edge)
!--- calculate & add scattering flux in layer:
            DTAUSCA = DTAU1(L)*POMEGAX(1,L,K)
            FLXD(L,K) = FLXD(L,K) + 0.5d0*(FTAU(L)+FTAU(L+1))*DTAUSCA
!--- calculate absorption decay along solar ray path thru half of layer
            DTAUABS = DTAU1(L) - DTAUSCA
            DECAY = exp(-0.5d0*DTAUABS*AMF(L,L))
          if (DECAY .gt. 0.01d0) then           ! decay no more than 99%
            FLXD(L,K) = FLXD(L,K) + (FTAU(L+1)*(1.d0-DECAY) &
                             + FTAU(L)*(1.d0/DECAY-1.d0))/AMF(L,L)
          else
           if (LL0 .eq. 0) then      ! sun coming from upper layer
            FLXD(L,K) = FLXD(L,K) + FTAU(L+1)/AMF(L,L)
           else                      ! twilight, sun from layer below
            FLXD(L,K) = FLXD(L,K) + FTAU(L)/AMF(L,L)
           endif
          endif
        else

!--- Interpolated layers, L22 > 0 for this L
         do LL = 0,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
            LZZ = LZ - 2*LL      ! LZZ = index(odd) of added level in scatt arrays
!--- calculate & add scattering flux in layer:  assume ZTAU is for 1D RT, not geom-scaled
            DTAUEXT = (ZTAU(LZZ,K)-ZTAU(LZZ-2,K)) / AMG(L)
            DTAUSCA = DTAUEXT*POMEGAX(1,L,K)
            FLXD(L,K) = FLXD(L,K) + 0.5d0*(FZ(LZZ,K)+FZ(LZZ-2,K))*DTAUSCA
!--- calculate absorption decay along solar ray path thru half of layer
            DTAUABS = DTAUEXT - DTAUSCA
            DECAY = exp(-0.5d0*DTAUABS*AMF(L,L))
          if (DECAY .gt. 0.01d0) then           ! decay no more than 99%
            FLXD(L,K) = FLXD(L,K) + (FZ(LZZ-2,K)*(1.d0-DECAY) &
                              + FZ(LZZ,K)*(1.d0/DECAY-1.d0))/AMF(L,L)
          else
           if (LL0 .eq. 0) then      ! sun coming from upper layer
            FLXD(L,K) = FLXD(L,K) + FZ(LZZ-2,K)/AMF(L,L)
           else                      ! twilight, sun from layer below
            FLXD(L,K) = FLXD(L,K) + FZ(LZZ,K)/AMF(L,L)
           endif
          endif
         enddo

        endif
! v76c final correction:
!    increase flux deposited in expanded spherical shells for 1D planar scattering code
        FLXD(L,K) = FLXD(L,K)*AMG(L)
       enddo

!---sum the solar flux depositied in CTM layers
       do L = 1,L1U
          FLXD0(K) = FLXD0(K) + FLXD(L,K)
       enddo

!   Now fill in the even 'h' points with simple interpolation in the scatter arrays:
       do LZ = 2,ND-1,2
          ZTAU(LZ,K) = 0.5d0*(ZTAU(LZ-1,K)+ZTAU(LZ+1,K))
          FZ(LZ,K)   = sqrt(FZ(LZ-1,K)*FZ(LZ+1,K))
          do I=1,M2_
             POMEGA(I,LZ,K) = 0.5d0*(POMEGA(I,LZ-1,K)+POMEGA(I,LZ+1,K))
          enddo
       enddo

      endif  ! k wavelength if end
      enddo  ! k wavelength loop end

!-----------------------------------------------------------------------
       call MIESCT(FJ,FJTOP,FJBOT,FIBOT, POMEGA,FZ,ZTAU,FSBOT,RFL,U0,LDOKR,ND)
!-----------------------------------------------------------------------

!---Integrate average std layer-L intensity from scatter array FJ(LZ=1:ND)
!  The following mapping holds for JXTRA=0
!        Surface:   TTAU(1)    ==> ZTAU(2*L1_+1)
!        Top:       TTAU(L1_)  ==> ZTAU(3)
!        Infinity:     0.0     ==> ZTAU(1)
!        index: 2*L1_ - 2*L + 3 ==> LZ
!
      do K = 1,W_+W_r
      if (LDOKR(K) .gt. 0) then

!---mean direct + diffuse intensity averaged throughout std model layer:
!           L to L+1 edge = L2LEV(L) to L2LEV(L+1)
!           corresponds to LZ1 (bottom edge) to LZ0 (upper edge)
       do L = 1,L1U
         LZ0 = ND+2 - 2*L2LEV(L+1)
         LZ1 = ND+2 - 2*L2LEV(L)
         SUMJ = (4.d0*FJ(LZ0,K)+FZ(LZ0,K))*(ZTAU(LZ0+1,K)-ZTAU(LZ0,K)) &
               +(4.d0*FJ(LZ1,K)+FZ(LZ1,K))*(ZTAU(LZ1,K)-ZTAU(LZ1-1,K))
         SUMT = ZTAU(LZ0+1,K)-ZTAU(LZ0,K) + ZTAU(LZ1,K)-ZTAU(LZ1-1,K)
        do LZ = LZ0+2,LZ1-2,2
         SUMJ =SUMJ+(4.d0*FJ(LZ,K)+FZ(LZ,K))*(ZTAU(LZ+1,K)-ZTAU(LZ-1,K))
         SUMT =SUMT + ZTAU(LZ+1,K)-ZTAU(LZ-1,K)
        enddo
         FJACT(L,K) = SUMJ/SUMT
       enddo

!---mean diffuse-only flux at the edge of the layer-L: 4<I*mu> (not solar)
!---derive from h's above/below LZ edge, inverse interpolated with tau
!--- does not do top or bottom
       do L = 2,L1U
        LZ  = ND+2 - 2*L2LEV(L)
        FJFLX0 = (ZTAU(LZ+1,K)-ZTAU(LZ,K))/(ZTAU(LZ+1,K)-ZTAU(LZ-1,K))
        FJFLX(L-1,K)=4.d0*(FJ(LZ-1,K)*FJFLX0 + FJ(LZ+1,K)*(1.d0-FJFLX0))
       enddo

!---diffuse fluxes reflected at top, incident at bottom:
!        FJTOP(K),FJBOT(K),FIBOT(1:5,K)

      endif
      enddo  ! wavelength loop!

      END SUBROUTINE OPMIE


!-----------------------------------------------------------------------
      subroutine MIESCT(FJ,FJT,FJB,FIB, POMEGA,FZ,ZTAU,FSBOT,RFL,U0,LDOKR,ND)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in)  ::  LDOKR(W_+W_r),ND
      real*8,  intent(in)  ::  POMEGA(M2_,N_,W_+W_r),FZ(N_,W_+W_r), &
                               ZTAU(N_,W_+W_r),RFL(5,W_+W_r),U0,FSBOT(W_+W_r)
      real*8,  intent(out) ::  FJ(N_,W_+W_r),FJT(W_+W_r),FJB(W_+W_r),FIB(5,W_+W_r)
      real*8  PM(M_,M2_),PM0(M2_)
      integer I, IM  ,K
!-----------------------------------------------------------------------
!   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Solution of inhomogeneous Rayleigh scattering atmosphere.
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!
!    takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field! (M=1)
!-----------------------------------------------------------------------
      do I = 1,M_
         call LEGND0 (EMU(I),PM0,M2_)
         do IM = 1,M2_
            PM(I,IM) = PM0(IM)
         enddo
      enddo


!---Note that U0 scattering does not change with altitude
      call LEGND0 (-U0,PM0,M2_)
      do IM=1,M2_
         PM0(IM) = 0.25d0*PM0(IM)
      enddo

!---BLKSLV now called with all the wavelength arrays (K=1:W_)

      call BLKSLV(FJ,POMEGA,FZ,ZTAU,FSBOT,RFL,PM,PM0,FJT,FJB,FIB,LDOKR,ND)

      END SUBROUTINE MIESCT


!-----------------------------------------------------------------------
      subroutine LEGND0 (X,PL,N)
!-----------------------------------------------------------------------
!---Calculates ORDINARY Legendre fns of X (real)
!---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
      implicit none
      integer, intent(in) :: N
      real*8, intent(in)  :: X
      real*8, intent(out) :: PL(N)
      integer I
      real*8  DEN
!---Always does PL(2) = P[1]
      PL(1) = 1.d0
      PL(2) = X
      do I = 3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.0/DEN) - PL(I-2)*(1.d0-1.d0/DEN)
      enddo

      END SUBROUTINE LEGND0


!-----------------------------------------------------------------------
      subroutine BLKSLV &
         (FJ,POMEGA,FZ,ZTAU,FSBOT,RFL,PM,PM0,FJTOP,FJBOT,FIBOT,LDOKR,ND)
!-----------------------------------------------------------------------
!  Sets up and solves the block tri-diagonal system:
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!  This goes back to the old, dumb, fast version 5.3
!-----------------------------------------------------------------------
!SJ!      USE IEEE_ARITHMETIC
      implicit none

      integer, intent(in) ::  LDOKR(W_+W_r),ND
      real*8, intent(in)  ::  POMEGA(M2_,N_,W_+W_r),FZ(N_,W_+W_r), &
                              ZTAU(N_,W_+W_r) ,PM(M_,M2_),PM0(M2_), &
                              RFL(5,W_+W_r),FSBOT(W_+W_r)
      real*8, intent(out) ::  FJ(N_,W_+W_r),FJTOP(W_+W_r),FJBOT(W_+W_r), &
                              FIBOT(5,W_+W_r)

      real*8, dimension(M_,N_,W_+W_r)    ::  A,C,H,   RR

      real*8, dimension(M_,M_,N_,W_+W_r) ::  B,AA,CC,  DD
      real*8, dimension(M_,M_) ::  E
      real*8  SUMB,SUMBX,SUMT,SUMBR,SUMRF
      integer I, J, K, L

      do K = 1,W_+W_r
      if (LDOKR(K) .gt. 0) then
       call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),FSBOT(K),RFL(1,K), &
             PM,PM0, B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K), &
                     A(1,1,K),H(1,1,K),C(1,1,K), ND)
      endif
      enddo

      do K = 1,W_+W_r
      if (LDOKR(K) .gt. 0) then
!-----------UPPER BOUNDARY L=1
       L = 1
        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,1,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,1,K) = -E(I,1)*CC(1,J,1,K)-E(I,2)*CC(2,J,1,K) &
                        -E(I,3)*CC(3,J,1,K)-E(I,4)*CC(4,J,1,K)
         enddo
          RR(J,1,K) = E(J,1)*H(1,1,K)+E(J,2)*H(2,1,K) &
                    +E(J,3)*H(3,1,K)+E(J,4)*H(4,1,K)
        enddo

!----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
       do L = 2,ND-1

        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) + A(I,L,K)*DD(I,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K) - A(J,L,K)*RR(J,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,L,K) = - E(I,J)*C(J,L,K)
         enddo
          RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

       enddo

!---------FINAL DEPTH POINT: L=ND
       L = ND
        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) &
           + AA(I,1,L,K)*DD(1,J,L-1,K) + AA(I,2,L,K)*DD(2,J,L-1,K) &
           + AA(I,3,L,K)*DD(3,J,L-1,K) + AA(I,4,L,K)*DD(4,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K) &
           - AA(J,1,L,K)*RR(1,L-1,K) - AA(J,2,L,K)*RR(2,L-1,K) &
           - AA(J,3,L,K)*RR(3,L-1,K) - AA(J,4,L,K)*RR(4,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    +E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

!-----------BACK SOLUTION
       do L = ND-1,1,-1
        do J = 1,M_
         RR(J,L,K) = RR(J,L,K) &
          + DD(J,1,L,K)*RR(1,L+1,K) + DD(J,2,L,K)*RR(2,L+1,K) &
          + DD(J,3,L,K)*RR(3,L+1,K) + DD(J,4,L,K)*RR(4,L+1,K)
        enddo
       enddo

!----------mean J & H
       do L = 1,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1) + RR(2,L,K)*WT(2) &
                + RR(3,L,K)*WT(3) + RR(4,L,K)*WT(4)
       enddo
       do L = 2,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1)*EMU(1) + RR(2,L,K)*WT(2)*EMU(2) &
                + RR(3,L,K)*WT(3)*EMU(3) + RR(4,L,K)*WT(4)*EMU(4)
       enddo

!---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
       SUMT = RR(1, 1,K)*WT(1)*EMU(1) + RR(2, 1,K)*WT(2)*EMU(2) &
            + RR(3, 1,K)*WT(3)*EMU(3) + RR(4, 1,K)*WT(4)*EMU(4)
       FJTOP(K) = 4.d0*SUMT

!---FJBOT = scaled diffuse flux onto surface
!---FSBOT = direct flux onto surface
!---FABOT = total flux absorbed at surface = FJBOT+FSBOT - reflected flux:
       SUMB = RR(1,ND,K)*WT(1)*EMU(1) + RR(2,ND,K)*WT(2)*EMU(2) &
            + RR(3,ND,K)*WT(3)*EMU(3) + RR(4,ND,K)*WT(4)*EMU(4)
       SUMBR= RR(1,ND,K)*WT(1)*EMU(1)*RFL(1,K) + &
              RR(2,ND,K)*WT(2)*EMU(2)*RFL(2,K) + &
              RR(3,ND,K)*WT(3)*EMU(3)*RFL(3,K) + &
              RR(4,ND,K)*WT(4)*EMU(4)*RFL(4,K)
       SUMRF= WT(1)*EMU(1)*RFL(1,K) + WT(2)*EMU(2)*RFL(2,K)+ &
              WT(3)*EMU(3)*RFL(3,K) + WT(4)*EMU(4)*RFL(4,K)
!---SUMBX = flux from Lambert reflected I+
       SUMBX = (4.d0*SUMBR + FSBOT(K)*RFL(5,K))/(1.d0 + 2.d0*SUMRF)

       FJBOT(K) = 4.d0*SUMB - SUMBX
! now keep the I+ (up from boundary) and the I- (incident diffuse)
       FIBOT(5,K) = SUMBX
       do J = 1,4
          FIBOT(J,K) = 2.d0*RR(J,ND,K) - SUMBX
       enddo

      endif
      enddo

      END SUBROUTINE BLKSLV


!-----------------------------------------------------------------------
      subroutine GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0 &
                    ,B,CC,AA,A,H,C,  ND)
!-----------------------------------------------------------------------
!  Generates coefficient matrices for the block tri-diagonal system:
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  ND
      real*8, intent(in)  ::  POMEGA(M2_,N_),PM(M_,M2_),PM0(M2_)
      real*8, intent(in)  ::  ZFLUX,RFL(5)
      real*8, intent(in),dimension(N_) :: FZ,ZTAU

      real*8, intent(out),dimension(M_,M_,N_) ::  B,AA,CC
      real*8, intent(out),dimension(M_,N_) ::  A,C,H

      integer I, J, K, L1,L2,LL
      real*8  SUM0, SUM1, SUM2, SUM3
      real*8  DELTAU, D1, D2, SURFAC,SUMRFL
!
      real*8, dimension(M_,M_) :: S,T,U,V,W
!---------------------------------------------

!---------upper boundary:  2nd-order terms
       L1 = 1
       L2 = 2
       do I = 1,M_
        SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
       + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
       + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
       + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
       + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
       + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
       + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
       + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
       + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
                + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
                + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
                + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo
!-------------upper boundary, 2nd-order, C-matrix is full (CC)
         DELTAU = ZTAU(L2) - ZTAU(L1)
         D2 = 0.25d0*DELTAU
       do I = 1,M_
        do J = 1,M_
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J)
         CC(I,J,L1) = D2*U(I,J)
        enddo
         H(I,L1) = H(I,L1) + 2.0d0*D2*C(I,L1)
         A(I,L1) = 0.0d0
       enddo
       do I = 1,M_
        D1 = EMU(I)/DELTAU
        B(I,I,L1)  = B(I,I,L1) + D1
        CC(I,I,L1) = CC(I,I,L1) - D1
       enddo

!------------intermediate points:  can be even or odd, A & C diagonal
!---mid-layer h-points, Legendre terms 2,4,6,8
       do LL=2,ND-1,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
           POMEGA(2,LL)*PM(I,2)*PM0(2) + POMEGA(4,LL)*PM(I,4)*PM0(4) &
         + POMEGA(6,LL)*PM(I,6)*PM0(6) + POMEGA(8,LL)*PM(I,8)*PM0(8))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 = &
           POMEGA(2,LL)*PM(I,2)*PM(J,2) + POMEGA(4,LL)*PM(I,4)*PM(J,4) &
          +POMEGA(6,LL)*PM(I,6)*PM(J,6) + POMEGA(8,LL)*PM(I,8)*PM(J,8)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        enddo
       enddo

!---odd-layer j-points, Legendre terms 1,3,5,7
       do LL=3,ND-2,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
           POMEGA(1,LL)*PM(I,1)*PM0(1) + POMEGA(3,LL)*PM(I,3)*PM0(3) &
         + POMEGA(5,LL)*PM(I,5)*PM0(5) + POMEGA(7,LL)*PM(I,7)*PM0(7))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 = &
           POMEGA(1,LL)*PM(I,1)*PM(J,1) + POMEGA(3,LL)*PM(I,3)*PM(J,3) &
          +POMEGA(5,LL)*PM(I,5)*PM(J,5) + POMEGA(7,LL)*PM(I,7)*PM(J,7)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        enddo
       enddo

!---------lower boundary:  2nd-order terms
       L1 = ND
       L2 = ND-1
       do I = 1,M_
        SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
       + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
       + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
       + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
       + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 = &
          POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
        + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 = &
          POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
        + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 = &
          POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
        + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 = &
          POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
        + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
                + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
                + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
                + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo

!------------lower boundary, 2nd-order, A-matrix is full (AA)
! v 7.6 now has separate albedos for each incident angle (4) & direct beam (u0)
!     downward flux at surface (u0*F) = ZFLUX
! new LBC with albedo depending on angle: flux up is weighted alebdo of incident
!
!      I-up * sum[EMU(i) * WT(i)] = 1/2 I-up
!                =sum[ I-down(i) * REFL(i) * EMU(i) * WT(i)]
!
!      J(i) = 1/2 * (I-down(i) + I-up)  ==>  I-down(i) = 2*J(i) - I-up
!  so
!      1/2 * I-up = sum[ (2*J(i) - I-up) * REFL(i)*EMU(i)*WT(i) ]
!
!     I-up * (1 + 2*sum[REFL*EMU*WT]) = 4 * sum[J(i)*REFL*EMU*WT]
! if REFL = constant then we get the 4*REFL/(1+REFL) scaling from Fast-J.

         DELTAU = ZTAU(L1) - ZTAU(L2)
         D2 = 0.25d0*DELTAU
         SUMRFL = 0.d0
       do J = 1,M_
         SUMRFL = SUMRFL + RFL(J)*EMU(J)*WT(J)  ! = 1/2 avg(RFL)
       ENDDO
         SURFAC = 4.d0/(1.d0 + 2.d0*SUMRFL)

       do I = 1,M_
          D1 = EMU(I)/DELTAU
          SUM0 = D1 + D2*(W(I,1)+W(I,2)+W(I,3)+W(I,4))
        do J = 1,M_
         AA(I,J,L1) = - D2*U(I,J)
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J) - SUM0*SURFAC*RFL(J)*EMU(J)*WT(J)
        enddo
         H(I,L1) = H(I,L1) -2.0d0*D2*C(I,L1) +SUM0*SURFAC*0.25*RFL(5)*ZFLUX
       enddo

       do I = 1,M_
          D1 = EMU(I)/DELTAU
        AA(I,I,L1) = AA(I,I,L1) + D1
        B(I,I,L1)  = B(I,I,L1) + D1
        C(I,L1) = 0.0d0
       enddo

      END SUBROUTINE GEN_ID


!<<<<<<<<<<<<<<<<<<<<<<<<<<end fastJX core code<<<<<<<<<<<<<<<<<<<<<<<<<



!<<<<<begin fastJX subroutines called from PHOTO_JX or OPMIE<<<<<<<<<<<<

!------------------------------------------------------------------------------
      subroutine OPTICL (REFF,TEFF, DDENS,QQEXT,SSALB,SSLEG)
!------------------------------------------------------------------------------
! new for FJ v7.5  for LIQUID water clouds only  interpolate properties to R_eff
! every S-bin has its own optical properties
! water clouds for (C1/Deir) GAMMA FN:alf=6  Reff from 1.5 to 48 microns

      implicit none

      real*8, intent(in) ::    REFF         ! effective radius of liq water cloud
      real*8, intent(in) ::    TEFF         ! effective temperature of ice water cloud
      real*8, intent(out) ::   DDENS        ! density of cloud particle (g/cm^3)
      real*8, intent(out)::    QQEXT(S_)    ! optical depth of layer
      real*8, intent(out)::    SSALB(S_)    ! single-scattering albedo
      real*8, intent(out)::    SSLEG(8,S_)  ! scatt phase fn (Leg coeffs)

      integer I,J,K,L, NR
      real*8  FNR

      K = 1   ! liquid water Mie clouds
      DDENS = DCC(K)
          I = 1      !must have at least 2 Reff bins, interpolate in Reff
      do NR = 2,MCC-1
        if (REFF .gt. RCC(NR,K)) then
          I = NR
        endif
      enddo
        FNR = (REFF - RCC(I,K)) / (RCC(I+1,K) - RCC(I,K))
        FNR = min(1.d0, max(0.d0, FNR))
! new - each wavelength S-bins J has its own indexed optical properties
      do J=1,S_
        QQEXT(J) = QCC(J,I,K) + FNR*(QCC(J,I+1,K)-QCC(J,I,K))
        SSALB(J) = SCC(J,I,K) + FNR*(SCC(J,I+1,K)-SCC(J,I,K))
       do L=1,8
        SSLEG(L,J) = PCC(L,J,I,K) + FNR*(PCC(L,J,I+1,K)-PCC(L,J,I,K))
       enddo
      enddo

      END SUBROUTINE OPTICL


!------------------------------------------------------------------------------
      subroutine OPTICI (REFF,TEFF, DDENS,QQEXT,SSALB,SSLEG)
!------------------------------------------------------------------------------
! new for FJ v7.5, parallel with liquid water, but two types of ice-water
! phase functions from a single calculation of Mishchenko, other opticals
!     based on Mie calc. for liquid water to get best possible abosorption in IR

      implicit none

      real*8, intent(in) ::    REFF         ! effective radius of liq water cloud
      real*8, intent(in) ::    TEFF         ! effective temperature of ice water cloud
      real*8, intent(out) ::   DDENS        ! density of cloud particle (g/cm^3)
      real*8, intent(out)::    QQEXT(S_)    ! optical depth of layer
      real*8, intent(out)::    SSALB(S_)    ! single-scattering albedo
      real*8, intent(out)::    SSLEG(8,S_)  ! scatt phase fn (Leg coeffs)

      integer I,J,K,L, NR
      real*8  FNR

        if (TEFF .ge. 233.15d0) then
      K = 2  ! ice irreg (warm)
        else
      K = 3  ! ice hexag (cold)
        endif
      DDENS = DCC(K)
          I = 1      !must have at least 2 Reff bins, interpolate in Reff
      do NR = 2,MCC-1
        if (REFF .gt. RCC(NR,K)) then
          I = NR
        endif
      enddo
        FNR = (REFF - RCC(I,K)) / (RCC(I+1,K) - RCC(I,K))
        FNR = min(1.d0, max(0.d0, FNR))

! new - each wavelength S-bins J has its own indexed optical properties
      do J=1,S_
        QQEXT(J) = QCC(J,I,K) + FNR*(QCC(J,I+1,K)-QCC(J,I,K))
        SSALB(J) = SCC(J,I,K) + FNR*(SCC(J,I+1,K)-SCC(J,I,K))
       do L=1,8
        SSLEG(L,J) = PCC(L,J,I,K) + FNR*(PCC(L,J,I+1,K)-PCC(L,J,I,K))
       enddo
      enddo

      END SUBROUTINE OPTICI


!------------------------------------------------------------------------------
      subroutine OPTICS (OPTD,SSALB,SLEG, PATH,K)
!------------------------------------------------------------------------------
!---for the UCI SSA (stratospheric sulfate aerosol) data sets
!---UCI aersols optical data  v-7.4+
! >>>special for Strat Sulfate Aerosols (SSA)!
! K = 01 S-Bkg   just use 220K & 70 wt%  values for now:  KK = 3 + 1 = 4
! K = 02 S-Vol   LOGN:r=.080 s=.800 n=1.514/.../1.435  refKK = 9 + 4 = 13
!>>> but the output OPTD, SSALB,SLEG now has a full SX-=27 wavelengths, not 5 (200-300-..-999mm)

      implicit none

      real*8, intent(in)::     PATH         ! path (g/m2) of aerosol/cloud
      integer,intent(inout)::     K            ! index of cloud/aerosols
      real*8, intent(out)::    OPTD(S_)    ! optical depth of layer
      real*8, intent(out)::    SSALB(S_)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,S_) ! scatt phase fn (Leg coeffs)

      integer I,J, KK
      real*8  XTINCT, REFF,RHO

      if (K .eq. 1) then
        KK = 4    ! background, 220K, 70 wt%
      elseif (K .eq. 2) then
        KK = 13   ! volcanic,   220K, 70 wt%
      else
        call EXITC ('OPTICS: SSA index out-of-range')
      endif

         REFF = RSS(KK)
         RHO =  DSS(KK)
      do J=1,S_
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
          XTINCT = 0.75d0*QSS(J,KK)/(REFF*RHO)
         OPTD(J) = PATH*XTINCT
         SSALB(J) = SSS(J,KK)
       do I=1,8
         SLEG(I,J) =  PSS(I,J,KK)
       enddo
      enddo
         K = KK

      END SUBROUTINE OPTICS


!------------------------------------------------------------------------------
      subroutine OPTICG (OPTD,SSALB,SLEG, PATH,K)
!------------------------------------------------------------------------------
!---for the GEOMIP SSA (stratospheric sulfate aerosol) data sets
! K = 1001:1015 corresponds to R-eff = 0.02 0.04 0.08 0.10 ...  1.4 2.0 3.0 5.0 microns
!     output OPTD, SSALB,SLEG now has a full SX-=27 wavelengths
      implicit none

      real*8, intent(in)::     PATH        ! path (g/m2) of aerosol/cloud
      integer,intent(in)::     K           ! index of cloud/aerosols
      real*8, intent(out)::    OPTD(S_)    ! optical depth of layer
      real*8, intent(out)::    SSALB(S_)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,S_)  ! scatt phase fn (Leg coeffs)

      integer I,J, KK
      real*8  XTINCT, REFF,RHO

      KK = max(1, min(NGG, K-1000))
      REFF = RGG(KK)
      RHO =  DGG(KK)
      do J=1,S_
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
         XTINCT = 0.75d0*QGG(J,KK)/(REFF*RHO)
         OPTD(J) = PATH*XTINCT
         SSALB(J) = SGG(J,KK)
         do I=1,8
            SLEG(I,J) =  PGG(I,J,KK)
         enddo
      enddo

      END SUBROUTINE OPTICG


!------------------------------------------------------------------------------
      subroutine OPTICA (OPTD,SSALB,SLEG, PATH,RELH,K)
!------------------------------------------------------------------------------
!---v-7.6+ no StratSulfAers (use OPTICS)  also interp/extrap a 1/wavelength
!         std 5 wavelengths:200-300-400-600-999nm
        implicit none
        real*8, intent(out)::    OPTD(S_)    ! optical depth of layer
        real*8, intent(out)::    SSALB(S_)   ! single-scattering albedo
        real*8, intent(out)::    SLEG(8,S_)  ! scatt phase fn (Leg coeffs)
        real*8, intent(in)::     PATH        ! path (g/m2) of aerosol/cloud
        real*8, intent(in)::     RELH        ! relative humidity (0.00->1.00+)
        integer,intent(inout)::     K        ! index of cloud/aerosols
        integer I,J,JMIE
        real*8  XTINCT, REFF,RHO,WAVE, QAAX,SAAX,WAAX
! K=1&2 are the SSA values, not used here any more, make sure they are not asked for.
        if (K.gt.NAA .or. K.lt.3) call EXITC ('OPTICA: aerosol index out-of-range')
        REFF = RAA(K)
        RHO = DAA(K)
        do J = 1,S_
           WAVE =  WL(J)      ! WL(1:S_=1:27) is in common = mean wavelength (nm)
!---Pick pair of Mie wavelength to get scattering properites--sorry for the hardwire here
           JMIE = 1
           WAAX = (WAVE - 200.d0)*0.010d0
           if( WAVE .gt. 300.d0 ) then
              JMIE = 2
              WAAX = (WAVE - 300.d0)*0.010d0
           endif
           if( WAVE .gt. 400.d0 ) then
              JMIE=3
              WAAX = (WAVE - 400.d0)*0.005d0
           endif
           if( WAVE .gt. 600.d0 ) then
              JMIE=4
              WAAX = (WAVE - 600.d0)*0.0025d0
           endif
           if( WAVE .gt. 999.d0 ) then
              QAAX = QAA(5,K) *999.d0/WAVE  ! Q gets smaller
              SSALB(J) = SAA(5,K)        ! single scat albedo & P1-P7 unchanged
              do I=1,8
                 SLEG(I,J) =  PAA(I,5,K)
              enddo
           else
              WAAX = min(1.d0,max(0.d0, WAAX))
              QAAX = QAA(JMIE,K)*(1.d0-WAAX) + QAA(JMIE+1,K)*WAAX
              SSALB(J) = SAA(JMIE,K)*(1.d0-WAAX) + SAA(JMIE+1,K)*WAAX
              do I=1,8
                 SLEG(I,J)= PAA(I,JMIE,K)*(1.d0-WAAX) + PAA(I,JMIE+1,K)*WAAX
              enddo
           endif
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
           XTINCT = 0.75d0*QAAX/(REFF*RHO)
           OPTD(J) = PATH*XTINCT
        enddo

      END SUBROUTINE OPTICA


!------------------------------------------------------------------------------
      subroutine OPTICM (OPTD,SSALB,SLEG, PATH,RELH,LL)
!------------------------------------------------------------------------------
!---U Michigan aerosol data sets, this generate fast-JX data formats.
!---Approximates the Legendre expansion(L) of the scattering phase fn as (2*L+1)*g**L
!---UMAER(I,J,K,L):
!   I=1:3 = [SSAbldeo, g, k-ext(m2/g)]
!   J=1:5 = [200, 300, 400, (550,) 600 , 1000 nm]
!   K=1:21= [0, 5, 10, 15, ..., 90, 95, 99 %RelHum]
!   L=1:33= UM aerosol types [SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4, FF00(0%BC),
!                      FF02, ...FF14(14%BC), BB00, BB02, ...BB30(30%BC)]
      implicit none

      real*8, intent(out)::    OPTD(S_)    ! optical depth of layer
      real*8, intent(out)::    SSALB(S_)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,S_)  ! scatt phase fn (Leg coeffs)
      real*8, intent(in)::     PATH       ! path (g/m2) of aerosol/cloud
      real*8, intent(in)::     RELH       ! relative humidity (0.00->1.00)
      integer,intent(in)::     LL         ! index of cloud/aerosols

      integer KR,J,L, JMIE
      real*8  R,FRH, GCOS, XTINCT, WAVE

!---calculate fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
!---extrapolate phase fn from first term (g)
      L = LL
      if (L.lt.1 .or. L.gt.33)  &
          call EXITC ('OPTICM: aerosol index out-of-range')
!---pick nearest Relative Humidity
      KR =  20.d0*RELH  + 1.5d0
      KR = max(1, min(21, KR))

      do J = 1,S_
         WAVE =  WL(J)
!---Pick nearest Mie wavelength to get scattering properites------------
                               JMIE=1  ! use 200 nm prop for <255 nm
        if( WAVE .gt. 255.d0 ) JMIE=2  ! use 300 nm prop for 255-355 nm
        if( WAVE .gt. 355.d0 ) JMIE=3  ! use 400 nm prop for 355-500 nm
        if( WAVE .gt. 500.d0 ) JMIE=4
        if( WAVE .gt. 800.d0 ) JMIE=5
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
          XTINCT = UMAER(3,JMIE,KR,L)
! rescale/reduce optical depth as 1/WL for > 1000 nm
        if( WAVE .gt. 1000.d0) XTINCT = XTINCT*1000.d0/WAVE
       OPTD(J) = PATH*XTINCT
       SSALB(J) = UMAER(1,JMIE,KR,L)
         GCOS   = UMAER(2,JMIE,KR,L)
       SLEG(1,J) =  1.d0
       SLEG(2,J) =  3.d0*GCOS
       SLEG(3,J) =  5.d0*GCOS**2
       SLEG(4,J) =  7.d0*GCOS**3
       SLEG(5,J) =  9.d0*GCOS**4
       SLEG(6,J) = 11.d0*GCOS**5
       SLEG(7,J) = 13.d0*GCOS**6
       SLEG(8,J) = 15.d0*GCOS**7
      enddo

      END SUBROUTINE OPTICM


!-----------------------------------------------------------------------
      subroutine JRATET(PPJ,TTJ,FFF, VALJL,LU,NJXU)
!-----------------------------------------------------------------------
! in:
!        PPJ(L_+1) = pressure profile at edges
!        TTJ(L_+1) = = temperatures at mid-level
!        FFF(K=1:NW, L=1:L_) = mean actinic flux
! out:
!        VALJL(L_,JX_)  JX_ = no of dimensioned J-values in CTM code
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  :: LU,NJXU
      real*8, intent(in)  ::  PPJ(LU+1),TTJ(LU+1)
      real*8, intent(inout)  ::  FFF(W_,LU)
      real*8, intent(out), dimension(LU,NJXU) ::  VALJL

      real*8  VALJ(X_)
      real*8  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
      real*8  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV

      if (NJXU .lt. NJX) then
        write(6,'(A,2I5)')  'NJXU<NJX',NJXU,NJX
        call EXITC(' JRATET:  CTM has not enough J-values dimensioned')
      endif

      do L = 1,LU

!---need temperature, pressure, and density at mid-layer (for some quantum yields):
        TT   = TTJ(L)
        if (L .eq. 1) then
          PP = PPJ(1)
        else
          PP  = (PPJ(L)+PPJ(L+1))*0.5d0
        endif
         DD = 7.24e18*PP/TT
!---must zero bin-11 (216-222 & 287-291 nm) below 100 hPa since O2 e-fold is too weak
        if (PP .gt. 100.d0) then
          FFF(11,L) = 0.d0
        endif
        do J = 1,NJX
          VALJ(J) = 0.d0
        enddo

!     for J=1:3  O2, O3(total), & O3(O1D)
        do K = 1,W_
          call X_interp (TT,QO2TOT, TQQ(1,1),QO2(K,1),TQQ(2,1),QO2(K,2), TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TT,QO3TOT, TQQ(1,2),QO3(K,1),TQQ(2,2),QO3(K,2), TQQ(3,2),QO3(K,3), LQQ(2))
          call X_interp (TT,QO31DY, TQQ(1,3),Q1D(K,1),TQQ(2,3),Q1D(K,2), TQQ(3,3),Q1D(K,3), LQQ(3))
          QO31D  = QO31DY*QO3TOT
          VALJ(1) = VALJ(1) + QO2TOT*FFF(K,L)
          VALJ(2) = VALJ(2) + QO3TOT*FFF(K,L)
          VALJ(3) = VALJ(3) + QO31D *FFF(K,L)
        enddo

        do J = 4,NJX
          do K = 1,W_
!---also need to allow for Pressure interpolation if SQQ(J) = 'p'
            if (SQQ(J) .eq.'p') then
              call X_interp (PP,QQQT, TQQ(1,J),QQQ(K,1,J), &
                   TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
            else
              call X_interp (TT,QQQT, TQQ(1,J),QQQ(K,1,J), &
                   TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
            endif
              VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
           enddo
        enddo

        do J=1,NJX
          VALJL(L,J) = VALJ(J)
        enddo

      enddo

      END SUBROUTINE JRATET


!-----------------------------------------------------------------------
      subroutine X_interp (TINT,XINT, T1,X1, T2,X2, T3,X3, L123)
!-----------------------------------------------------------------------
!  up-to-three-point linear interpolation function for X-sections
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in)::  TINT,T1,T2,T3, X1,X2,X3
      integer,intent(in)::  L123
      real*8, intent(out)::  XINT

      real*8  TFACT

      if (L123 .le. 1) then
           XINT = X1
      elseif (L123 .eq. 2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
      else
        if (TINT.le. T2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
        else
             TFACT = max(0.d0,min(1.d0,(TINT-T2)/(T3-T2) ))
           XINT = X2 + TFACT*(X3 - X2)
        endif
      endif

      END SUBROUTINE X_interp


!-----------------------------------------------------------------------
      subroutine JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU6,POMEG6,JXTRA,LU)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!---the CTM has L_ = LU layers and fast-JX adds layer LU+1
!---the pressure and altitude(Z) are on layer edge (LU+2)
      integer,intent(in)                  :: LU
      real*8, intent(in), dimension(LU+2) :: PPJ,ZZJ
      real*8, intent(in), dimension(LU+1) :: TTJ,DDJ,OOJ,DTAU6
      real*8, intent(in), dimension(8,LU+1) :: POMEG6
      integer,intent(in), dimension(LU+1) :: JXTRA
!-----------------------------------------------------------------------
      integer  I,J,K,L
      real*8   XCOLO2,XCOLO3,ZKM,DELZ,ZTOP,DAIR,DOZO

      write(6,'(4a)') '   L z(km)     p      T   ', &
       '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb', &
       '  g(cos) CTM lyr=>'
      L = LU+2
      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZZJ(L)*1.d-5,PPJ(L)
      XCOLO2 = 0.d0
      XCOLO3 = 0.d0
      ZTOP = ZZJ(LU+2)
      do L = LU+1,1,-1
        XCOLO2 = XCOLO2 + DDJ(L)*0.20948d0
        XCOLO3 = XCOLO3 + OOJ(L)
        DELZ = ZTOP-ZZJ(L)
        ZTOP = ZZJ(L)
        ZKM = ZZJ(L)*1.d-5
        DAIR = DDJ(L)/DELZ
        DOZO = OOJ(L)/DELZ
        write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZKM,PPJ(L),TTJ(L),DAIR,DOZO,XCOLO2,XCOLO3,DTAU6(L), &
            POMEG6(1,L),POMEG6(2,L)/3.d0, JXTRA(L)
      enddo

      END SUBROUTINE JP_ATM


!-----------------------------------------------------------------------
      subroutine JP_ATM0(PPJ,TTJ,DDJ,OOJ,ZZJ, LU)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!---Atmosphere print, called from outside fjx_sub_mod.f90
!---CTM layers are 1:LU, + top layer (to P=0) added
!---pressure and altitude are on layer edge (1:LU+2)
      integer,intent(in)                  :: LU
      real*8, intent(in), dimension(LU+2) :: PPJ,ZZJ
      real*8, intent(in), dimension(LU+1) :: TTJ,DDJ,OOJ
!-----------------------------------------------------------------------
      integer  I,J,K,L
      real*8   XCOLO2,XCOLO3,ZKM,DELZ,ZTOP
      write(6,'(4a)') '   L z(km)     p      T   ', &
       '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb', &
       '  g(cos) CTM lyr=>'
      L = LU+2
      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZZJ(L)*1.d-5,PPJ(L)
          XCOLO2 = 0.d0
          XCOLO3 = 0.d0
          ZTOP = ZZJ(LU+2)
        do L = LU+1,1,-1
          XCOLO2 = XCOLO2 + DDJ(L)*0.20948d0
          XCOLO3 = XCOLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5
      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,OOJ(L)/DELZ, &
            XCOLO2,XCOLO3
        enddo

      END SUBROUTINE JP_ATM0


!-----------------------------------------------------------------------
      subroutine SPHERE1R (U0,RAD,ZHL,ZZHT,AMF, L1U)
!-----------------------------------------------------------------------
!  version 7.6  - SPHERE1N = drops the mid-layer (v6.2) for comp cost
!     also 7.6  - SPHERE1R = adds refraction (complex ray tracing)
!  for computation cost, not called for SZA > 98 deg (full night)
!-----------------------------------------------------------------------
! in:
!     U0      cos(solar zenith angle)
!     RAD     radius of Earth mean sea level (cm)
!     ZHL(L)  height (cm) of the bottome edge of CTM level L
!     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1))
!     L1U     dimension of CTM = levels +1 (L+1 = above-CTM level)
! out:
!     AMF(K,L) = air mass factor for solar path for CTM layer K at radius L
!          notation is that layer K is bound by radius K and radius K+1
!     AMF(L,L) >0 is used to that there is a path to the sun at radius L
!          even on the dark side, the path will start thru layers K < L
!          but must eventually go up through layer L to reach the sun.
!     AMF(LTOP+1,LTOP+1) is meaningless since there is no layer LTOP+1, and is set =0
!          but when we need to calculate solar intensity at radius LTOP+1
!          on the darkside, we set this = 1. to trigger the AMF(L,L)>0 criterion.
! local:
!     RZ(L) = radius vector length to bottom edge of layer L  = RAD+ZHL(L)
!     ZHL(L1U) = top radius of CTM
!     ZHL(L1U+1) = top radius of atmosphere = RZ(L1U) + ZZHT (in cm)
!     LTOP = L1U = top radius of CTM layers, LTOP+1 = top of atmosphere
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in) ::   L1U
      real*8, intent(in)  ::   U0,RAD,ZHL(L1_+1),ZZHT
      real*8, intent(out) ::   AMF(L1_+1,L1_+1)

      integer  L,L0, K,K0, LTOP
      real*8   ZA0,SZA0,CZA0
      real*8   ZA1,SZA1,SRN0,SA0,A0,CA0,SRN1,SA1,A1,CA1,SA2,A2,CA2
      real*8   DDHT,REF0,F0, C90
      real*8, dimension(L1_+1,L1_+1) :: ZANG,ZAMF
      real*8, dimension(L1_+1) :: RZ,DIVZ,RATZ,RD,RN, PATH1,PATH2,ZANG1
!-----------------------------------------------------------------------
!  this versions sets a density scale ht of DDHT=8km, and a
!        refractive index of 1.000300 at radius = RAD,
!        the 0.000300 scales with atmospheric density
      DDHT = 8.0d5
      REF0 = 300.d-6
      C90 = 1.570796326794897d0
!-----------------------------------------------------------------------
      LTOP = L1U
      do L = 1,LTOP
        RZ(L) = RAD + ZHL(L)
      enddo
        RZ(LTOP+1) = RZ(LTOP) + ZZHT
      do L = 1,LTOP
        DIVZ(L) = 1.d0/(RZ(L+1)-RZ(L))
        RATZ(L) = RZ(L)/RZ(L+1)
        RD(L) = exp(-(RZ(L)-RAD)/DDHT)    ! assume layer density is at the bottom edge
        RN(L) = 1.d0 + REF0*RD(L)
      enddo
        RD(LTOP+1) = 0.d0
        RN(LTOP+1) = 1.d0
        CZA0 = U0
        ZA0 = acos(CZA0)
        SZA0 = sin(ZA0)
      AMF(:,:) = 0.d0
      AMF(LTOP+1,LTOP+1) = 1.d0
!-----------------------------------------------------------------------
        ZAMF(:,:) = 0.d0
        ZANG(:,:) = 0.d0

      if (U0 .lt. 0.d0) goto 1111

! ZA0 .le. 90 deg. = Solar path directly on surface *** even without refraction ***
!  do first downward integration with refraction,
!  calculate elevation angle at each layer edge,
!  redo each path for each layer edge
!  with the elevation angle found in first integration.
! Loop over each lowermost point L
      do L=1,LTOP
! First time thru, just get elevation angle correction (no paths)
         SRN0 = SZA0*RZ(L)*RN(L)    ! invariant path that hits edge L at angle ZA0
         SA0 = SRN0/(RZ(LTOP+1)*RN(LTOP+1)) ! = sin of zenith angle at top of atmosphere
         ZANG1(:) = 0.d0
         ZANG1(LTOP+1) = asin(SA0)
         PATH1(:) = 0.d0
       do K=LTOP,L,-1
!  A1 = zenith angle at bottom layer K (from invariant)
         SA1 = SRN0/(RZ(K)*RN(K))
!  A2 = zenith angle at top of layer K (from triangle geometry)
         SA2 = SA1*RATZ(K)
          A1 = asin(SA1)
          A2 = asin(SA2)
         ZANG1(K) = ZANG1(K+1) + A1 - A2
       enddo
! correct zenith angle at lower edge L by subtracting elevation angle = ZANG1(L)-ZA0 > 0
         ZA1 = ZA0 - (ZANG1(L) - ZA0)
! Second time with elev angle correction. Calculate path. but cumulative angle not needed
         SZA1 = sin(ZA1)
         SRN1 = SZA1*RZ(L)*RN(L)    ! invariant path that hits surface at angle ZA0
         SA0 = SRN0/(RZ(LTOP+1)*RN(LTOP+1)) ! A0 = zenith angle at top of atmosphere
         PATH1(:) = 0.d0
       do K=LTOP,L,-1
!  A1 = zenith angle at bottom layer K (from invariant)
         SA1 = SRN1/(RZ(K)*RN(K))
!  A2 = zenith angle at top of layer K (from triangle geometry)
         SA2 = SA1*RATZ(K)
          A1 = asin(SA1)
          A2 = asin(SA2)
          CA1 = cos(A1)
          CA2 = cos(A2)
         PATH1(K) = RZ(K+1)*CA2 - RZ(K)*CA1  ! path length thru layer K (edges K:K+1)
       enddo
       do K=1,LTOP
        AMF(K,L) = PATH1(K)*DIVZ(K)
       enddo
      enddo
      goto 2222

 1111 continue

! integrate the refracted rath paths that are tangent at each radius RZ
! calculate the zenith angle at the tangent point (Atan(L) > 90)
      do L=1,LTOP+1
         SRN0 = RZ(L)*RN(L)    ! invariant path that hits tangent at RZ(L)
         SA0 = SRN0/(RZ(LTOP+1)*RN(LTOP+1)) ! = sin of zenith angle at top of atmosphere
        ZANG1(:) = 0.d0
        ZANG1(LTOP+1) = asin(SA0)
        PATH1(:) = 0.d0
! begin going downward to tangent layer L (RZ(L) = lower layer edge of layer)
       do K=LTOP,L,-1
!  A1 = zenith angle at bottom layer K (from invariant)
         SA1 = SRN0/(RZ(K)*RN(K))
!  A2 = zenith angle at top of layer K (not refracted)
         SA2 = SA1*RATZ(K)          !  RZ(K)/RZ(K+1)
          A1 = asin(SA1)
          A2 = asin(SA2)
          CA1 = cos(A1)
          CA2 = cos(A2)
         PATH1(K) = RZ(K+1)*CA2 - RZ(K)*CA1  ! path length thru layer K (edges K:K+1)
         ZANG1(K) = ZANG1(K+1) + A1 - A2
       enddo
! this back-up integration is not necessary, symmetric in PATH and in angle about ZANG1(L)
       do K=1,LTOP
        ZAMF(L,K) = PATH1(K)*DIVZ(K)
        ZANG(L,K) = ZANG1(L)+ZANG1(L)-ZANG1(K)
       enddo
        K = LTOP+1
        ZANG(L,K) = ZANG1(L)+ZANG1(L)-ZANG1(K)
      enddo

! there are only LTOP layers with PATHs (air mass factors)
! each layer has its own a path for each edge point
! and layer LTOP is layer (w/thickness) above the model layers to reach top-of-atmos
! define 2 tangent paths to interpolate AMF (PATH2) between K0 and K0+1
! the PATH lengths on the tangent arcs are symmetric (but not about 90 deg)
      do L0 = 1,LTOP+1

       if (ZA0 .lt. ZANG(L0,L0)) then
! do NOT calculate AMF's for top radius/edge LTOP+1 unless lit from below, else AMF=0
        if (L0 .le. LTOP) then
! the column atmosphere zenith angle is less than the tangent point angle
! correct zenith angle at lower edge L by subtracting elevation angle = ZANG(L0,L0) - C90
          ZA1 = ZA0 - (ZANG(L0,L0) - C90)
! Second time with elev angle correction. Calculate path. but cumulative angle not needed
          SZA1 = sin(ZA1)
          SRN1 = SZA1*RZ(L0)*RN(L0)    ! invariant path that hits surface at angle ZA0
          SA0 = SRN0/(RZ(LTOP+1)*RN(LTOP+1)) ! A0 = zenith angle at top of atmosphere
          PATH1(:) = 0.d0
         do K=LTOP,L0,-1
!  A1 = zenith angle at bottom layer K (from refracted invariant)
          SA1 = SRN1/(RZ(K)*RN(K))
!  A2 = zenith angle at top of layer K (from triangle geometry)
          SA2 = SA1*RATZ(K)
          A1 = asin(SA1)
          A2 = asin(SA2)
          CA1 = cos(A1)
          CA2 = cos(A2)
          PATH1(K) = RZ(K+1)*CA2 - RZ(K)*CA1  ! path length thru layer K (edges K:K+1)
         enddo
         do K=1,LTOP
          AMF(K,L0) = PATH1(K)*DIVZ(K)
         enddo
        endif
       else
! post-terminator path: ZA0 > ZANG(1,1) and hence for all ZANG(L0,L0) to top
           K0 = 0
        do K=1,L0-1
         if (ZA0 .le. ZANG(K,L0)) then   ! check if still in Earth shadow
           K0 = K
         endif
        enddo
        if (K0 .gt. 0) then
          F0 = (ZA0 - ZANG(K0+1,L0))/(ZANG(K0,L0)-ZANG(K0+1,L0))
         do L=1,LTOP
          PATH2(L) = F0*ZAMF(K0,L) + (1-F0)*ZAMF(K0+1,L)
         enddo
! load all the paths for each level (pre-terminator)
         do L=1,LTOP
          AMF(L,L0) = PATH2(L)
         enddo
! add on the post-terminator paths *below* level L0=1:LTOP+1 (1:L0-1)
         do L=1,L0-1
          AMF(L,L0) = AMF(L,L0) + PATH2(L)
         enddo
        endif
       endif
      enddo

 2222 continue

      END SUBROUTINE SPHERE1R


!-----------------------------------------------------------------------
      subroutine SPHERE1N (U0,RAD,ZHL,ZZHT,AMF, L1U)
!-----------------------------------------------------------------------
!  version 7.6a  - SPHERE1N = drops the mid-layer (v6.2) for comp cost
!     also 7.6b  - SPHERE1R = adds refraction (complex ray tracing)
!  for computation cost, not called for SZA > 98 deg (full night)
!-----------------------------------------------------------------------
! in:
!     U0      cos(solar zenith angle)
!     RAD     radius of Earth mean sea level (cm)
!     ZHL(L)  height (cm) of the bottome edge of CTM level L
!     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1))
!     L1U     dimension of CTM = levels +1 (L+1 = above-CTM level)
! out:
!     AMF(K,L) = air mass factor for solar path for CTM layer K at radius L
!          notation is that layer K is bound by radius K and radius K+1
! local:
!     RZ(L) = radius vector length to bottom edge of layer L  = RAD+ZHL(L)
!     ZHL(L1U) = top radius of CTM
!     ZHL(L1U+1) = top radius of atmosphere = RZ(L1U) + ZZHT (in cm)
!     LTOP = L1U = top radius of CTM layers, LTOP+1 = top of atmosphere
!-----------------------------------------------------------------------

      implicit none
      integer, intent(in) ::   L1U
      real*8, intent(in)  ::   U0,RAD,ZHL(L1_+1),ZZHT
      real*8, intent(out) ::   AMF(L1_+1,L1_+1)

      integer  L, J, JUP, LTOP, K
      real*8   A0,A1,A2,SA0,SA1,SA2,CA0,CA1,CA2,R0, PATH,ZA0,CZA0,SZA0

      real*8, dimension(L1_+1) :: RZ,DIVZ,RATZ

      LTOP = L1U
      do L = 1,LTOP
        RZ(L) = RAD + ZHL(L)
      enddo
        RZ(LTOP+1) = RZ(LTOP) + ZZHT
      do L = 1,LTOP
        DIVZ(L) = 1.d0/(RZ(L+1)-RZ(L))
        RATZ(L) = RZ(L)/RZ(L+1)
      enddo
        CA0 = U0
        A0 = acos(CA0)
        SA0 = sin(A0)
        R0 = RZ(1)
      AMF(:,:) = 0.d0
      AMF(LTOP+1,LTOP+1) = 1.d0
!-----------------------------------------------------------------------

      if (CA0 .ge. 0.d0) then
! surface in direct sunlight, follow radius vector through layer LTOP
! for each point along the radius L, calculate AMF(J,L)
! AMF(J,L) = optical layer J contribution to solar path at radius level L
        do L=1,LTOP
           SA1 = SA0    ! starting angle at radius point L
           A1 = asin(SA1)
           CA1 = cos(A1)
         do J=L,LTOP
           SA2 = SA1*RATZ(J)        !  RZ(J)/RZ(J+1)
           A2 = asin(SA2)
           CA2 = cos(A2)
           PATH = RZ(J+1)*CA2 - RZ(J)*CA1
          AMF(J,L) = PATH*DIVZ(J)
           SA1 = SA2
           CA1 = CA2
           A1 = A2
         enddo
        enddo
      else
! surface dark, search upward in radius to find a point in sunlight
       do L=2,LTOP+1
        if (SA0*RZ(L) .gt. R0) then
           SA1 = SA0
           CA1 = CA0
           A1 = A0
         do J=L-1,1,-1
           if (SA1*RZ(J+1) .lt. RZ(J)) then  ! path from R(J+1) down to R(J)
             SA2 = SA1/RATZ(J)
             A2 = asin(SA2)
             CA2 = -cos(A2)
            PATH = RZ(J+1)*CA2 - RZ(J)*CA1
            AMF(J,L) = AMF(J,L) + PATH*DIVZ(J)  ! AMF for layer J (to radius L)
             SA1 = SA2
             CA1 = CA2
             A1 = A2
           else                   ! path across terminator (CA=0) in layer J
            PATH = -2.d0*CA1*RZ(J+1)
            AMF(J,L) = AMF(J,L) + PATH*DIVZ(J)   ! AMF for layer J (to radius L)
             CA1 = -CA1
             JUP = J+1
            goto 2        ! end this downward ray, it appears on the other side
           endif
         enddo
     2   continue
         do J=JUP,LTOP             ! start back up SA1 is now lower radius
           SA2 = SA1*RATZ(J)
           A2 = asin(SA2)
           CA2 = cos(A2)
           PATH = RZ(J+1)*CA2 - RZ(J)*CA1
          AMF(J,L) = AMF(J,L) + PATH*DIVZ(J)
           SA1 = SA2
           CA1 = CA2
           A1 = A2
         enddo
        else
        endif
       enddo  ! end of loop over L
      endif    !  end of if CA0 < 0 test, finished

      END SUBROUTINE SPHERE1N


!-----------------------------------------------------------------------
      subroutine SPHERE1F (U0,RAD,ZHL,ZZHT,AMF, L1U)
!-----------------------------------------------------------------------
!     needed for testing flat-disk errors
!  version 7.6a  - SPHERE1N = drops the mid-layer (v6.2) for comp cost
!     also 7.6b  - SPHERE1R = adds refraction (complex ray tracing)
!  for computation cost, not called for SZA > 98 deg (full night)
!-----------------------------------------------------------------------
! in:
!     U0      cos(solar zenith angle)
!     RAD     radius of Earth mean sea level (cm)
!     ZHL(L)  height (cm) of the bottome edge of CTM level L
!     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1))
!     L1U     dimension of CTM = levels +1 (L+1 = above-CTM level)
! out:
!     AMF(K,L) = air mass factor for solar path for CTM layer K at radius L
!          notation is that layer K is bound by radius K and radius K+1
!        = 1/U0
!-----------------------------------------------------------------------

      implicit none
      integer, intent(in) ::   L1U
      real*8, intent(in)  ::   U0,RAD,ZHL(L1_+1),ZZHT
      real*8, intent(out) ::   AMF(L1_+1,L1_+1)

      integer  L, J, LTOP
      real*8   PATH0

      LTOP = L1U
      AMF(:,:) = 0.d0
      AMF(LTOP+1,LTOP+1) = 1.d0

      if (U0 .gt. 0.d0) then
        PATH0 = 1.d0/U0

        do L=1,LTOP
         do J=L,LTOP
          AMF(J,L) = PATH0
         enddo
        enddo

      endif

      END SUBROUTINE SPHERE1F


!-----------------------------------------------------------------------
      subroutine EXTRAL1(DTAU600,L1X,NX,ATAU,ATAU0, JXTRA)
!-----------------------------------------------------------------------
!     version 7.6 replaces v 6.2 and drops back to no mid-layer J(odd) points.
!   Purpose:  reduce spurious negative heating at top of thick clouds.
!---Divide thick layers to achieve better accuracy in the scattering code
!---The key parameters are:
!---        ATAU = factor increase from one layer to the next
!---        ATAU0 = delta-TAU cut-off for cloud OD to insert a layer
!--- Best values are 0.02 & 1.20, doubles # total layers for thick clouds,
!       can use 1.40 with 30% reduction in cost, but 3x larger negative heating,
!       likewise, and use 1.12 but 30% increase in cost, factor of 3 reduction in neg.
!  If cloud top heating (K=1:27) is 24 K/day, (high sun, thick cloud)
!    the visible (k=1:18) error is -2.3 (1.4 -30% layers),
!                                  -0.6 (1.2 best case),
!                                  -0.2 (1.12, +30% layers)
!
!     DTAU600(L=1:L1X) = Optical Depth in layer L, generally 600 nm OD
!        v7.6 = aerosols + cloud ONLY
!     JXTRA(L=1:L1x)  the number in levels to insert in each layer L
!-----------------------------------------------------------------------

      implicit none
      integer, intent(in) ::  L1X            !# layers
      integer, intent(in) ::  NX             !Mie scattering array size (max)
      real*8,  intent(in) ::  DTAU600(L1X)     !cloud+3aerosol OD in each layer
      real*8,  intent(in) ::  ATAU,ATAU0
      integer, intent(out)::  JXTRA(L1X)    !number of sub-layers to be added
      integer JTOTL,JX,L,LL
      real*8  ATAULN,ATAU0X,AJX,DTAU0X

!  need to divide DTAU600 into JX layers such that DTAU600/ATAU0 = ratio =
!       1 + ATAU + ATAU^2 + ATAU^3 + ATAU^(JX-1)  = [ATAU^JX - 1]/[ATAU - 1]
!  then JX = ln(1 + ratio*(ATAU-1)) / ln(ATAU) round off suitably
!       and the last layer added has optical depth = ATAU0 * ATAU^JX
!  note that there are JX+1 sub layers of DTAU inserted.  JX=0, means DTAU600 stays same
         ATAULN = log(ATAU)
         ATAU0X = ATAU0
      do L = L1X,1,-1
         JXTRA(L) = 0
         if (DTAU600(L) .gt. ATAU0X) then    ! ratio DTAU600/ATAU0 > 1
            AJX = log(1.d0 + (ATAU-1.d0)*DTAU600(L)/ATAU0X) / ATAULN
            JX = min(100, max(0, int(AJX + 0.5d0)))
            JXTRA(L) = JX
            ATAU0X = ATAU0X*(ATAU**(JX))
         endif
      enddo
!---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL = L1X + 2
      do L = L1X,1,-1
         JTOTL  = JTOTL + JXTRA(L)
         if (JTOTL*2 .gt. NX)  then
            write(6,'(A,3I5)') 'cutoff EXTRAL1 @L: NX/L1X/L: ',NX,L1X,L
            do LL = L,1,-1
               JXTRA(LL) = 0
            enddo
!           call exitc('STOP at EXTRAL') !not necessary, a warning is OK
            go to 10
         endif
      enddo
  10  continue

      END SUBROUTINE EXTRAL1


!-----------------------------------------------------------------------
      subroutine SOLAR_JX(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
!-----------------------------------------------------------------------
! >>>>>>>> warning tnot specific for SOLAR-J, is it old FAST_J call
!     GMTIME = UT for when J-values are wanted
!           (for implicit solver this is at the end of the time step)
!     NDAY   = integer day of the year (used for solar lat and declin)
!     YGRDJ  = laitude (radians) for grid (I,J)
!     XGDRI  = longitude (radians) for grid (I,J)
!
!     SZA = solar zenith angle in degrees
!     COSSZA = U0 = cos(SZA)
!-----------------------------------------------------------------------
      implicit none

      real*8,  intent(in)  ::  GMTIME,YGRDJ,XGRDI
      integer, intent(in)  ::  NDAY
      real*8,  intent(out) ::  SZA,COSSZA,SOLFX
!
      real*8  LOCT
      real*8  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ
!
      SINDEC = 0.3978d0*sin(0.9863d0*(dble(NDAY)-80.d0)*CPI180)
      SOLDEK = asin(SINDEC)
      COSDEC = cos(SOLDEK)
      SINLAT = sin(YGRDJ)
      SOLLAT = asin(SINLAT)
      COSLAT = cos(SOLLAT)
!
      LOCT   = (((GMTIME)*15.d0)-180.d0)*CPI180 + XGRDI
      COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
      SZA    = acos(COSSZA)/CPI180
!
      SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*C2PI/365.d0))
!
      END SUBROUTINE SOLAR_JX


!SJ!    !!!!!!!!!!!!!!!!! SOLAR-J specific subroutines

!---------------------------------------------------------------------
      subroutine  FJX_CLIRAD_H2O(nlayers, PPP, TTT, HHH, TAUG_CLIRAD)
!---------------------------------------------------------------------
      implicit none

      integer,  intent(in):: nlayers
      real*8 ,  intent(in) :: PPP(nlayers+1), TTT(nlayers), HHH(nlayers)
      real*8 ,  intent(out):: TAUG_CLIRAD(nlayers, 0:30)

      integer G, K, INDKG, L
      real*8, dimension(nlayers):: WT
! heating rate (K/day) = W/m2 deposited in layer * HeatFac_ / delta-P of layer (hPa)
      real*8,  parameter:: HeatFac_ = 86400.d0*9.80616d0/1.00464d5, S0=1360.8
      real*8   HHX, pavg
!-----------------------------------------------------------------------
!    Based on  Chu and Lee (JAS, 1996 Parameterizations for the Absorption of Solar Radiation by Water Vapor and Ozone)
!    Three NIR band wavelength range = 0.70-1.22  1.22-2.27  2.27-10.0 microns
!    WGF = fraction in sub-bin, WGK = H2O Xsection (cm2/g)
!    WGTOT (is redundant here); note that unlike the code in Grant and Grossman, the sum of WGF (g) for each band doesn't add up to 1 and is a fracton of total S0 input
!       real*8, dimension(3), parameter  :: WGTOT =[1, 1, 1]
! WGF is recalculated to add to 1 for each solar bin and tabulated in FJX_spec_Clirad.dat
! Delta-g/WGF (fractions of S0 for 10 sub-bins)
!      real*8, dimension(10,3), parameter  :: WGF = [ &
!       0.20673, 0.03497, 0.03011, 0.02260, 0.01336, 0.00696, 0.00441, 0.00115, 0.00026, 0.00000, &
!       0.08236, 0.01157, 0.01133, 0.01143, 0.01240, 0.01258, 0.01381, 0.00650, 0.00244, 0.00094, &
!       0.01074, 0.00360, 0.00411, 0.00421, 0.00389, 0.00326, 0.00499, 0.00465, 0.00245, 0.00145]
!2012 HITRAN H2O
      real*8, dimension(10,3), parameter  :: WGF =[ &
          .14983, .04105, .04049, .03372, .02792, .01904, .00862, .00235, .00051, .00002, &
          .06266, .01465, .01286, .01165, .01311, .01728, .01949, .01032, .00358, .00122, &
          .00692, .00283, .00350, .00375, .00471, .00534, .00499, .00580, .00339, .00176]
! k-interval
!      real*8, dimension(10,1), parameter :: WGK  =[&
!       0.0010, 0.0133, 0.0422, 0.1334, 0.4217, 1.3340, 5.6230, 31.620, 177.80, 1000.0] !cm2/g
! revised according to 2012 HITRAN H2O molecular line data; unit cm2/gram
      real*8, dimension(10), parameter :: WGK  =[&
       0.0010, 0.0032, 0.0102, 0.0328, 0.1049, 0.4194, 2.5166, 17.616, 123.31, 839.19]
      real*8, parameter :: CMF = 2.98897027277D-23 ! 18.d0 divided by Avagado number

! 0:0 will assign to bin 18 which is 0 for CLIRAD
      TAUG_CLIRAD(:,0:0)= 0.d0
      do L=1, nlayers
         pavg= 0.5d0* (PPP(L)+PPP(L+1))
! Weighting function Pr=300 hPa, Tr=240K;
         WT(L)= (pavg/300.d0)**(0.8) * (1+ 0.00135*(TTT(L)- 240.0d0)) + 1.D-9 !exp(x) ~= 1+x
! In literatue, the original eq. looks like the one below. So the above expression is an approximation.
!        WT(L) =  (pavg/300.d0)**(0.8) * exp(0.00135*(tavg- 240.0d0)) !from Chou (1986, Journal of  Clim. and Appl. Meteo. )
      enddo
! we want g/cm2 of h2o for each box
      INDKG = 0
      do K=1,3
         do G = 1,10
            INDKG = INDKG + 1
!input solar flux
            do L = nlayers, 1,-1
!after L box  sol2 left
               HHX = HHH(L)* CMF !convert molecules/cm2 to g/cm2
               TAUG_CLIRAD(L, INDKG)= WGK(G)*HHX *WT(L) ! optical depth for H2O WGK cm2/gram; HHH (g/cm2); WT (no unit)
            enddo
         enddo
      enddo

      return

      END SUBROUTINE  FJX_CLIRAD_H2O


!!!!!!!!!!!!!!!!!!! SOLAR-J specific subroutines
!---------------------------------------------------------------------
      subroutine  FJX_GGLLNL_H2O(nlayers, PPP, TTT, HHH, TAUG_LLNL)
!---------------------------------------------------------------------
      implicit none

      integer,  intent(in):: nlayers
      real*8 ,  intent(in) :: PPP(nlayers+1), TTT(nlayers), HHH(nlayers)
      real*8 ,  intent(out):: TAUG_LLNL(nlayers, 0:21)

      integer G, K, INDKG, L
      real*8, dimension(nlayers):: WT
! heating rate (K/day) = W/m2 deposited in layer * HeatFac_ / delta-P of layer (hPa)
      real*8,  parameter:: HeatFac_ = 86400.d0*9.80616d0/1.00464d5, S0=1360.8
      real*8   HHX, pavg
!-----------------------------------------------------------------------
!   Grant & Grossman 1998 3-band Solar IR, based on Chu 1992
!     wavelengths = 0.69¡V0.86  0.86¡V2.27  2.27¡V3.85 microns
!     WGTOT = W/m2 in band, WGF = fraction in sub-bin, WGK = H2O Xsection (cm2/g)
!  WGF and WGTOT are tabulated FJX_spec_GGLLNL.dat
      real*8, dimension(3), parameter  :: WGTOT =[209.77, 472.71, 46.788]
      real*8, dimension(7,3), parameter  :: WGF = [ &
       0.948551, 0.051449, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, &
       0.703232, 0.085079, 0.098956, 0.046725, 0.049153, 0.016855, 0.000000, &
       0.389153, 0.106278, 0.142613, 0.118942, 0.151376, 0.068448, 0.023190]
!WGK cross sections unit cm2 g-1
      real*8, dimension(7,3), parameter  :: WGK = [ &
           4.3980E-03, 2.4676E-01,        0.0,        0.0,        0.0,        0.0, 0.0, &
           7.6655E-03, 1.3370E-01, 5.3350E-01, 2.3126E+00, 1.0536E+01, 1.3122E+02, 0.0, &
           1.4989E-02, 1.3525E-01, 5.3707E-01, 3.1426E+00, 2.1238E+01, 1.8492E+02, 1.6292E+03]

      real*8, parameter :: CMF = 2.98897027277D-23 ! 18.d0 divided by Avagado number

! 0:0 will assign to bin 18 which is 0 for CLIRAD
      TAUG_LLNL(:,0)= 0.d0
      do L=1, nlayers
         pavg= 0.5d0* (PPP(L)+PPP(L+1))
! Weighting function Pr=300 hPa, Tr=240K;
         WT(L)= (pavg/300.d0)**(0.8) * (1+ 0.00135*(TTT(L)- 240.0d0)) + 1.D-9 !exp(x) ~= 1+x
! In literatue, the original eq. looks like the one below. So the above expression is an approximation.
!        WT(L) =  (pavg/300.d0)**(0.8) * exp(0.00135*(tavg- 240.0d0)) !from Chou (1986, Journal of  Clim. and Appl. Meteo. )
      enddo
! we want g/cm2 of h2o for each box
      INDKG=0
      do K=1,3
         do G = 1,7
            INDKG = INDKG + 1
!input solar flux
            do L = nlayers, 1,-1
!after L box  sol2 left
               HHX = HHH(L)* CMF !convert molecules/cm2 to g/cm2
               TAUG_LLNL(L, INDKG)= WGK(G,K)*HHX *WT(L) ! optical depth for H2O WGK cm2/gram; HHH (g/cm2); WT (no unit)
!               write(6,'(A10, 2I5, f12.8)')'LLNL tau', L, INDKG, TAUG_LLNL(L,indkg)
            enddo
         enddo
      enddo

      return

      END SUBROUTINE  FJX_GGLLNL_H2O




!-----------------------------------------------------------------------
      subroutine ACLIM_FJX (YLATD,MONTH,PPP, TTT,O3,CH4, L1U)
!-----------------------------------------------------------------------
!  Load fast-JX climatology - T & O3 - for latitude & month & pressure grid
!-----------------------------------------------------------------------
      implicit none
      real*8,  intent(in)  :: YLATD
      integer, intent(in)  :: MONTH, L1U
      real*8,  intent(in),  dimension(L1U+1) :: PPP
      real*8,  intent(out), dimension(L1U)   :: TTT,O3,CH4
      real*8, dimension(LREF)   :: OREF2,TREF2,HREF2,CREF2
      real*8, dimension(LREF+1) :: PSTD
      integer  K, L, M, N
      real*8   DDDL,DLOGP,F0,T0,H0,C0,PB,PC,XC
!  Select appropriate month
      M = max(1,min(12,MONTH))
!  Select appropriate latitudinal profiles
      N = max(1, min(18, (int(YLATD+99)/10 )))
      do K = 1,LREF
         OREF2(K) = O_REF(K,N,M)
         TREF2(K) = T_REF(K,N,M)
         HREF2(K)=  H2O_REF(K,N,M)   ! H2O and CH4 not returned for _FJX
         CREF2(K)=  CH4_REF(K,N,M)
      enddo
!  Apportion O3 and T on supplied climatology z levels onto CTM levels +1
!  with mass (pressure) weighting, assuming constant mixing ratio and
!  temperature half a layer on either side of the point supplied.
!   PPP(L=1:L1_)=edge-pressure of CTM layer, PPP(L1_+1)=0 (top-of-atmos)
!  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
!  Set up pressure levels for O3/T climatology - assume that value
!  given for each 2 km z* level applies from 1 km below to 1 km above,
!  so select pressures at these boundaries. Surface level values at
!  1000 mb are assumed to extend down to the actual PSURF (if > 1000)
      PSTD(1) = max(PPP(1),1000.d0)
      PSTD(2) = 1000.d0 * 10.d0**(-1.d0/16.d0)
      DLOGP   = 10.d0**(-2.d0/16.d0)
      do K = 3,LREF
         PSTD(K) = PSTD(K-1)*DLOGP
      enddo
      PSTD(LREF+1)  = 0.d0
      do L = 1,L1U
         F0 = 0.d0
         T0 = 0.d0
         H0=  0.d0
         C0=  0.d0
         do K = 1,LREF
            PC   = min(PPP(L),PSTD(K))
            PB   = max(PPP(L+1),PSTD(K+1))
            if (PC .gt. PB) then
               XC = (PC-PB)/(PPP(L)-PPP(L+1))
               F0 = F0 + OREF2(K)*XC
               T0 = T0 + TREF2(K)*XC
               H0 = H0 + HREF2(K)*XC
               C0 = C0 + CREF2(K)*XC
            endif
         enddo
         TTT(L)= T0  ! K
         O3(L) = F0  ! ppm
         CH4(L)= C0  ! ppb
      enddo

      END SUBROUTINE ACLIM_FJX


!-----------------------------------------------------------------------
      subroutine ACLIM_RH (PL, TL, QL, RH, L1U)
!-----------------------------------------------------------------------
!  Calculates RH profile given PL(mid-pressure), TL(K), QL (spec hum)
!  May nee RH @ L1U (top layer, not CTM) so aerosol calls are stable
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in):: L1U
      real*8,  intent(in),  dimension(L1U) :: PL,TL,QL
      real*8,  intent(out), dimension(L1U) :: RH
! local variables
      real*8  T, eps, es, qs
      integer  L

      eps=287.04d0/461.50d0
      do L = 1,L1U-1
         if (TL(L) .gt. 273.15) then
            T = TL(L)- 273.15
            es = 6.112*exp(17.67*T/(T+243.50))
         else
            T= TL(L) !in Kelvin
            es= 23.33086- 6111.72784/T + 0.15215*log(T)
            es =exp(es)
         endif
         qs = (eps *es)/(PL(L)-es*(1-eps))
         RH(L)= min(max(QL(L)/qs, 0.d0), 1.d0)
      enddo
      RH(L1U)= RH(L1U-1)

      END SUBROUTINE ACLIM_RH


!-----------------------------------------------------------------------
      subroutine ACLIM_GEO (YLATD,MONTH,PPP, AERS,NAER, L1U)
!-----------------------------------------------------------------------
!  Load GEOMIP SSA climatology (vs P) for latitude & month given pressure grid
!-----------------------------------------------------------------------
!---note that trigger for using GEOMIP aerosol properties is NAER(L) = 1001:1015
!   numbers 1:nnn are reserved for standard ssa and trop aerosols.
      implicit none
      real*8,  intent(in)  :: YLATD
      integer, intent(in)  :: MONTH, L1U
      real*8,  intent(in),  dimension(L1U+1) :: PPP
      real*8,  intent(out), dimension(L1U)   :: AERS
      integer, intent(out), dimension(L1U)   :: NAER

      real*8, dimension(LGREF+2) :: RREF2,XREF2,PREF2    ! param LGREF=19
      real*8, dimension(LGREF+3) :: PSTD2
      integer  I, IGG, K, L, M, N
      real*8   R0,X0,RX0,PB,PC,XC,YN,REDGE(GGA_)

!  Select appropriate month
      M = max(1, min(12, MONTH))
!  Select appropriate latitudinal profiles J=1:64, delta = 2.7906 close enough
      N = 1
      YN = -86.5806d0
      do while (YLATD .gt. YN)
         N = N+1
         YN = YN + 2.7906d0
      enddo
      N = max(1, min(64, N))
!---P_GREF = 2.7,.... 339 hPa (reverse order) ensure ZERO ssa above and below
      RREF2(:) = 0.d0
      XREF2(:) = 0.d0
      do K = 1,LGREF
         PREF2(K+1) = P_GREF(LGREF+1-K)
         RREF2(K+1) = R_GREF(N,LGREF+1-K,M)
         XREF2(K+1) = X_GREF(N,LGREF+1-K,M)
      enddo
      PREF2(1) = PREF2(2) * 1.001d0
      PREF2(LGREF+2) = PREF2(LGREF+1) * 0.999d0
!---re-set PSTD2 to the boundaries between the PREF2 points
      PSTD2(1) = max(PPP(1),PREF2(1))
      do L = 2,LGREF+2
         PSTD2(L) = 0.5d0*(PREF2(L-1)+PREF2(L))
      enddo
      PSTD2(LGREF+3)  = 0.d0
      do I = 2,NGG-1
         REDGE(I) = 0.5d0*(RGG(I)+RGG(I-1))
      enddo
!---integrate for pressure-wtd averages in each layer L = (PPP(L) to PPP(L+1)
      do L = 1,L1U
         X0 = 0.d0
         RX0 = 0.d0
         do K = 1,LGREF
            PC   = min(PPP(L),PSTD2(K))
            PB   = max(PPP(L+1),PSTD2(K+1))
            if (PC .gt. PB) then
               XC = (PC-PB)/(PPP(L)-PPP(L+1))
               X0 = X0 + XREF2(K)*XC
               RX0 = RX0 + RREF2(K)*XREF2(K)*XC
            endif
         enddo
!---for each model layer L, calculate
!   aerosol path in each layer AERS(L) (g/m2)
!   index IGG of closest effective radius to R0 in optical tables
!   AERSL path from X=microg-H2SO4/kg-air (ppb) with 75 wt% of H2SO4 (mass*1/.75)
!     rescale it so that when RGG is used for OD, we get same as if R0 used.
         AERS(L) =  G100*(PPP(L)-PPP(L+1)) * X0 * 1.3333d-6
         NAER(L) = 1001
!---for each model layer L, pick the nearest R-eff and rescale the mass
         if (X0 .gt. 0.d0) then
            R0 = RX0/X0
            IGG = 1
            do I = 2,NGG-1
               if (R0 .gt. REDGE(I)) then
                  IGG = I
               endif
            enddo
            NAER(L) = min(1000+IGG, 1000+GGA_)
            AERS(L) = AERS(L) * RGG(IGG)/R0
         endif
      enddo

      END SUBROUTINE ACLIM_GEO


!-----------------------------------------------------------------------
      subroutine EXITC(T_EXIT)
!-----------------------------------------------------------------------
      character(len=*), intent(in) ::  T_EXIT
      write(6,'(a)') T_EXIT
      stop

      END SUBROUTINE EXITC


      END MODULE FJX_SUB_MOD
