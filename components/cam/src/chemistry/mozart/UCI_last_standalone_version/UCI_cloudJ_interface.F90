module UCI_cloudJ_interface
  !>>>>>>>>  Cloud-J version 7.7
  !
  ! It just uses read in SZA, and LatxMONTH for a climatology of T & O3
  !       CLDFLAG = 1  :  Clear sky J's
  !       CLDFLAG = 2  :  Averaged cloud cover
  !       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
  !       CLDFLAG = 4  :  ****not used
  !       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
  !       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
  !       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin)
  !       CLDFLAG = 8  :  Calculate J's for ALL ICAs (up to 20,000 per cell!)

  use ppgrid,       only : pcols, pver, pverp, begchunk, endchunk
  use chem_mods,    only : ncol_abs => nabscol, phtcnt, pcnstm1 => gas_pcnst
  use dust_model,   only : ndst => dust_nbin
  use mo_constants, only : pi,r2d,boltz,d2r
  use shr_kind_mod, only : r8 => shr_kind_r8

  use cam_history,      only : fieldname_len
  use cam_logfile,      only : iulog
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : iam, masterproc

!!! I/O utilities.  Used by Xactive_Photo.   Use for FastJ?
! use pio
! use cam_pio_utils,only : cam_pio_openfile

  implicit none
  private

!---------------------------------------------------------------------------------
! Public interfaces
!---------------------------------------------------------------------------------
  public :: cloudJ_interface
 
!================================================================================================
contains
!================================================================================================
!---------------------------------------------------------------------------------
! cloudJ_interface routine
!---------------------------------------------------------------------------------

  subroutine cloudJ_interface(photos, vmr, temper, cwat, cldfr, &
                            pmid, zmid, col_dens, zen_angle, srf_alb, &
                            tdens, ps, ts, esfact, relhum, dust_vmr, &
 !                           dt_diag, fracday, &
                            ncol, lchnk)
     
     USE FJX_CMN_MOD
     USE FJX_SUB_MOD
     USE FJX_INIT_MOD
     USE CLD_SUB_MOD, ONLY : CLOUD_JX
     USE OSA_SUB_MOD

     !----------- CloudJ_interface arguments -------------------

     integer,  intent(in)    :: ncol, lchnk
     real(r8), intent(in)    :: esfact                       ! earth sun distance factor
     real(r8), intent(in)    :: ps(pcols)                    ! surface pressure (Pa)
     real(r8), intent(in)    :: ts(ncol)                     ! surface temperature (K)
     real(r8), intent(in)    :: col_dens(ncol,pver,ncol_abs) ! column densities (molecules/cm^2)
     real(r8), intent(in)    :: zen_angle(ncol)              ! solar zenith angle (radians)
     real(r8), intent(in)    :: srf_alb(pcols)               ! surface albedo
     real(r8), intent(in)    :: tdens(ncol,pver)             ! total atms density (molecules/cm^3)
     real(r8), intent(in)    :: vmr(ncol,pver,pcnstm1)       ! species concentration (mol/mol)
     real(r8), intent(in)    :: pmid(pcols,pver)             ! midpoint pressure (Pa)
     real(r8), intent(in)    :: zmid(ncol,pver)              ! midpoint height (m)
     real(r8), intent(in)    :: temper(pcols,pver)           ! midpoint temperature (K)
     real(r8), intent(in)    :: relhum(ncol,pver)            ! relative humidity
     real(r8), intent(in)    :: cwat(ncol,pver)              ! cloud water (kg/kg)
     real(r8), intent(in)    :: cldfr(ncol,pver)             ! cloud fraction
     real(r8), intent(in)    :: dust_vmr(ncol,pver,ndst)     ! dust concentration (mol/mol)
     real(r8), intent(inout) :: photos(ncol,pver,phtcnt)     ! photodissociation rates (1/s)
!     real(r8), intent(out)   :: dt_diag(pcols,8)             ! od diagnostics
!     real(r8), intent(out)   :: fracday(pcols)               ! fraction of day
     

!---------------key params in/out of CLOUD_J-------------------------
     logical                    :: LPRTJ, LDARK
     integer                    :: IRAN
     integer                    :: JVNU,ANU,L1U
     integer                    :: NICA,JCOUNT
     real*8                     :: U0,SZA,SOLF
     real*8,  dimension(L2_  )  :: PPP,ZZZ
     real*8,  dimension(L1_  )  :: TTT,HHH,DDD,RRR,OOO,CCC
     real*8,  dimension(L1_  )  :: O3,CH4,H2O, OD18
     real*8,  dimension(S_+2,L1_):: SKPERD
     real*8,  dimension(6)       :: SWMSQ
     real*8,  dimension(L2_)    :: CLF,LWP,IWP,REFFL,REFFI
     integer, dimension(L2_)    :: CLDIW
     real*8,  dimension(L2_,AN_):: AERSP
     integer, dimension(L2_,AN_):: NDXAER
     real*8,  dimension(L_,JVN_) :: VALJXX
     real*8,  dimension(5,S_)    :: RFL
     real*8,  dimension(NQD_)    :: WTQCA
     
!-------------local use-----------------------
     integer :: NSZA,MONTH,ILAT
! beware the OSA code uses single R*4 variables
     real*4  :: OWAVEL,OWIND,OCHLR,OSA_dir(5)
     real*8  :: YLAT,PSURF,ALBEDO(5),WIND,CHLR
     real*8, dimension(L2_) :: ETAA,ETAB, ZOFL,RI,TI,AER1,AER2,PPPX
     integer,dimension(L2_) :: NAA1,NAA2
     real*8, dimension(L_)  :: WLC,WIC
     real*8, dimension(LWEPAR) :: CLDFRW,CLDIWCW,CLDLWCW
     real*8  SCALEH,CF,PMID_UCI,PDEL,ZDEL,ICWC,F1
     integer I,J,K,L,N
     integer LTOP, NJXX,JP04,JP09
     character*6,  dimension(JVN_)  ::  TITLJXX
     character*11, dimension(4)     ::  TITJX
     real*8 VJOSA(L2_,2),VJSTD(L2_,2)
     
     integer, dimension(10) ::  &
          SZAscan = [ 0, 30, 60, 80, 86, 88, 90, 92, 94, 96]

     integer, parameter :: num_E3SM_atm_levels = 57   ! PJC: this will need to be changed to get from E3SM.
 
     if ( .not.masterproc ) return     !pjc

!$OMP MASTER
     
     write(iulog,*) 'cloudJ_interface, iam = ',iam   !pjc

     ANU = AN_     ! max # FJX aerosols in layer
     JVNU = JVN_   ! max # of J-values.  Currently set in UCI_fjx_cmn_mod.f90.  Can it be set by E3SM?
     L1U = L1_     ! model levels +1

     
!!! Fast-J sets the number of atmospheric levels via the LPAR
!!! parameter in the FJX_CMN_MOD module. This must be consistent with
!!! E3SM.  The easiest solution is to just abort if the # of levels
!!! are inconsistent.   [PJC]
     
     IF ( num_E3SM_atm_levels .NE. LPAR ) THEN
        write (iulog,*) 'cloudj_interface: num_E3SM_atm_levels =',num_E3SM_atm_levels,',  LPAR =',LPAR
        write (iulog,*) 'cloudj_interface: LPAR may need to be modified to match model levels.'
        call endrun('cloudj_interface: ERROR. # of model layers in E3SM and Cloud-J are different.')
     ENDIF
     
!---read in & store all fast-JX data:   single call at set up
!-----------------------------------------------------------------------

!!!  JVNU passed in to INIT_FJX
!!!  NJXX passed back from INIT_FJX. Number of species to calculate J-values for.
!!!  TITLJXX passed back from INIT_FJX.  

     call INIT_FJX (TITLJXX,JVNU,NJXX)
!-----------------------------------------------------------------------
     WRITE (iulog,*) '*** PJC: After INIT_FJX ***'
!--P, T, Cld & Aersl profiles, simple test input case
     open (77,file='tables/atmos_PTClds.dat',status='old',SHARED) !UNIX

     MONTH = 8       ! month
     YLAT  = 20.     ! Latitude
     PSURF = 1013.   ! Surface Pressure


     ALBEDO(1) = 0.05   ! Default albedos for quad angles = 86 deg zenith angle
     ALBEDO(2) = 0.05   ! Default albedos for quad angles = 71 deg zenith angle
     ALBEDO(3) = 0.05   ! Default albedos for quad angles = 48 deg zenith angle
     ALBEDO(4) = 0.05   ! Default albedos for quad angles = 22 deg zenith angle
     ALBEDO(5) = 0.05   ! Default albedo for SZA incident ray

     WIND = 6.0     ! Surface wind in m/s (at what altitude?)
     CHLR = 0.08    ! Chlorophyl mg/m3
     
     write(iulog,'(a,2i5,5x,a,i5)') 'Atmosphere:',LPAR,LWEPAR, 'LPAR / LWEPAR', L1_
     write(iulog,'(a,f10.4)') 'P surface', PSURF
     write(iulog,'(a,3i4)') 'MONTH/ LAT',MONTH,NINT(YLAT)
     write(iulog,'(a,5f8.4)') 'Albedos 1:4 & 5=SZA', ALBEDO
     write(iulog,'(a,2f8.3)') 'OSA: wind & chlor-a',WIND,CHLR
     read (77,*)     ! Read 7 header lines that aren't used anymore.
     read (77,*)
     read (77,*)
     read (77,*)
     read (77,*)
     read (77,*)
     read (77,*)
     do L = 1,LPAR+1      ! Level interfaces
        read (77,'(i3,1x,2f11.7,2x,f5.1,f5.2,f11.2,2(f7.3,i4))') &
             J,ETAA(L),ETAB(L),TI(L),RI(L),ZOFL(L) &     
             ,AER1(L),NAA1(L),AER2(L),NAA2(L)
        ! eta-A   eta-B   Temp(K)   RH(fraction)  equivalent_altitude(m)   AER-P NDA  AER-P NDA
     enddo
     read (77,*)
     do L = LWEPAR,1,-1   ! Clouds in levels   
        read (77,'(i3,1p,e14.5,28x,2e14.5)') &
             J,CLDFRW(L),CLDLWCW(L),CLDIWCW(L)
        ! Cloud_fraction?    WLC(g/g)    WIC(g/g)
     enddo
     close(77)
     
     do L = 1,L1_
        PPP(L) = ETAA(L) + ETAB(L)*PSURF    ! Pressure on level interfaces (mbar)
     enddo
!!!! just for print out and levels          ! PJC: Not used in this code
!!!     do L = 1,L1_
!!!        PPPX(L) = 0.5d0*(PPP(L)+PPP(L+1))   ! Mid-layer pressure (mbar)
!!!     enddo
!---sets climatologies for O3, T, D & Z
!-----------------------------------------------------------------------
     call ACLIM_FJX (YLAT,MONTH,PPP, TTT,O3,CH4, L1_)
        ! Returns Temperature(K), O3(ppm), CH4(ppb) (in layers plus above model top?)
!-----------------------------------------------------------------------
     do L = 1,L_
!!!       TTT(L) = TI(L)  keep climatology T's and O3's
        RRR(L) = RI(L)
     enddo
     ZZZ(1)  = 16.d5*log10(1013.25d0/PPP(1))        ! zzz in cm
     do L = 1,L_
        DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
!!! geopotential since assumes g = constant
        SCALEH      = 1.3806d-19*MASFAC*TTT(L)
        ZZZ(L+1) = ZZZ(L) -( log(PPP(L+1)/PPP(L)) * SCALEH )
        OOO(L) = DDD(L)*O3(L)*1.d-6
        CCC(L) = DDD(L)*CH4(L)*1.d-9
     enddo
     L = L_+1
     ZZZ(L+1)= ZZZ(L) + ZZHT
     DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
     OOO(L) = DDD(L)*O3(L)*1.d-6
     CCC(L) = DDD(L)*CH4(L)*1.d-9
!-----------------------------------------------------------------------
!       call ACLIM_RH (PL, TL, QL, HHH, L1U)
!-----------------------------------------------------------------------
! quick fix Rel Humidity
     HHH(:) = 0.50d0
!!! set up clouds and aerosols
     AERSP(:,:)  = 0.d0
     NDXAER(:,:) = 0
     do L = 1,L_
        NDXAER(L,1) = NAA1(L)
        AERSP(L,1)  = AER1(L)
        NDXAER(L,2) = NAA2(L)
        AERSP(L,2)  = AER2(L)
     enddo
     LTOP  = LWEPAR
     if (maxval(CLDFRW) .le. 0.005d0) then
        IWP(:) = 0.d0
        REFFI(:) = 0.d0
        LWP(:) = 0.d0
        REFFL(:) = 0.d0
     endif
     do L = 1,LTOP
        CLDIW(L) = 0
        CF  = CLDFRW(L)
        if (CF .gt. 0.005d0) then
           CLF(L) = CF
           WLC(L) = CLDLWCW(L) / CF
           WIC(L) = CLDIWCW(L) / CF
!  CLDIW is an integer flag: 1 = water cld, 2 = ice cloud, 3 = both
           if (WLC(L) .gt. 1.d-11) CLDIW(L) = 1
           if (WIC(L) .gt. 1.d-11) CLDIW(L) = CLDIW(L) + 2
        else
           CLF(L) = 0.d0
           WLC(L) = 0.d0
           WIC(L) = 0.d0
        endif
     enddo
!---derive R-effective for clouds:  the current UCI algorithm - use your own
     do L = 1,LTOP
!---ice clouds
        if (WIC(L) .gt. 1.d-12) then
           PDEL = PPP(L) - PPP(L+1)
           ZDEL = (ZZZ(L+1) - ZZZ(L))*0.01d0  ! m
           IWP(L) = 1000.d0*WIC(L)*PDEL*G100 /CLF(L)   ! g/m2
           ICWC =        IWP(L) / ZDEL          ! g/m3
           REFFI(L) = 164.d0 * (ICWC**0.23d0)
        else
           IWP(L) = 0.d0
           REFFI(L) = 0.d0
        endif
!---water clouds
        if (WLC(L) .gt. 1.d-12) then
           PMID_UCI = 0.5d0*(PPP(L)+PPP(L+1))
           PDEL = PPP(L) - PPP(L+1)
           F1   = 0.005d0 * (PMID_UCI - 610.d0)
           F1   = min(1.d0, max(0.d0, F1))
           LWP(L) = 1000.d0*WLC(L)*PDEL*G100     ! g/m2
           REFFL(L) = 9.6d0*F1 + 12.68d0*(1.d0-F1)
        else
           LWP(L) = 0.d0
           REFFL(L) = 0.d0
        endif
     enddo
     do L = 1,LTOP
        CLDFRW(L) = CLF(L)
     enddo
!!!  end of atmosphere setup


!!! begin call to Cloud_J

     do I=1,10    !!!! begin of SZA scan (if wanted)
        NSZA = SZAscan(I)
        SZA = NSZA
        
        write(iulog,'(i5,a)') CLDFLAG, ' CLDFLAG'
        do L = 1,LTOP
           CLF(L) = CLDFRW(L)
        enddo
        IRAN = 1
        SOLF = 1.d0
        U0 = cos(SZA*CPI180)       !  U0 = cos(sza)

! beware the OSA code uses single R*4 variables
        ANGLES(1) = sngl(EMU(1))
        ANGLES(2) = sngl(EMU(2))
        ANGLES(3) = sngl(EMU(3))
        ANGLES(4) = sngl(EMU(4))
        ANGLES(5) = sngl(U0)
        OWIND = sngl(WIND)
        OCHLR = sngl(CHLR)
        do K = 1,NS2
           OWAVEL = sngl(WL(K))
           call OSA(OWAVEL,OWIND,OCHLR, ANGLES,OSA_dir)
           do J = 1,5
              RFL(J,K) = dble(OSA_dir(J))
! this overwrite the OSA with the readin values above
              RFL(J,K) = ALBEDO(J)
           enddo
        enddo

        LPRTJ = .true.
        if (LPRTJ) then
           write(iulog,'(a,f8.3,3f8.5)')'SZA SOLF U0 albedo' &
                ,SZA,SOLF,U0,RFL(5,18)
           call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)
           write(iulog,*) ' wvl  albedo u1:u4 & u0'
           do K=1,NS2
              write(iulog,'(i5,f8.1,5f8.4)') K,WL(K), (RFL(J,K), J=1,5)
           enddo
        endif

        SKPERD(:,:)=0.d0
        SWMSQ(:)= 0.d0
        OD18(:) =0.d0
        WTQCA(:)= 0.d0
        
!=======================================================================
        call CLOUD_JX (U0,SZA,RFL,SOLF,LPRTJ,PPP,ZZZ,TTT,HHH,DDD,       &
               RRR,OOO,CCC,  LWP,IWP,REFFL,REFFI, CLF,CLDCOR,CLDIW,    &
               AERSP,NDXAER,L1U,ANU,JVNU, VALJXX,SKPERD,SWMSQ,OD18,    &
               CLDFLAG,NRANDO,IRAN,LNRG,NICA, JCOUNT,LDARK,WTQCA)
!=======================================================================



! summary print, can change unit to N=7 for file without print from Fast-J
        N=7
        write(N,'(a,2i5)') ' v7.7 CLDFLAG/NSZA=', CLDFLAG,NSZA
        write(N,*) ' LDARK WTQCA',LDARK,WTQCA
        write(N,'(a)') ' AVERAGE Fast-J  v7.7 ----J-values----'
        write(N,'(1x,a,72(a6,3x))') 'L=  ',(TITLEJX(K), K=1,NJX)
        do L = L_,1,-1
           write(N,'(i3,1p, 72e9.2)') L,(VALJXX(L,K),K=1,NJX)
        enddo

        SKPERD(:,:)=NINT( SKPERD(:,:)*100 ) / 100. !Rounding to 2 decimal places for consistent output.
        
        write(N,'(a)') 'heating rate profiles in K/day v7.6  180-778nm '
        write(N, '(a4, 32f7.1)')'wvl ',(WL(K),K=NW1,NW2)
        do L = L_,1,-1
           write(N,'(i4,32f7.2)') L,(SKPERD(K,L), K=NW1,NW2)
        enddo
        write(N,'(a)') 'heating rate profiles in K/day v7.6 778-...nm plus 1:18 19:27 1:27'
        write(N, '(a4,32f7.1)')'wvl ',(WL(K),K=NW2+1,NS2)
        do L = L_,1,-1
           write(N,'(i4,35f7.2)') L,(SKPERD(K,L), K=NW2+1,NS2+2),SKPERD(S_+1,L)+SKPERD(S_+2,L)
        enddo

        write(N,'(a)') 'Fast-J  v7.7 ---PHOTO_JX internal print: Solar fluxes (W/m2)--'
        write(N,'(a11,f12.4)')    ' inc TOTAL ',SWMSQ(1)
        write(N,'(a11,f12.4)')    ' rfl outtop',SWMSQ(2)
        write(N,'(a11,f12.4)')    ' abs in atm',SWMSQ(3)
        write(N,'(a11,f12.4)')    ' abs at srf',SWMSQ(4)
        write(N,'(a11,1p,e12.4)') ' PAR direct',SWMSQ(5)
        write(N,'(a11,1p,e12.4)') ' PAR diffus',SWMSQ(6)
        
!---map the J-values from fast-JX onto CTM (ZPJQUAD) using JIND & JFACTA
!--- from the 'FJX_j2j.dat' tables
!      do J = 1,NRATJ
!         JP = JIND(J)
!         if (JP .gt. 0) then
!            do L = 1,L_
!               ZPJQUAD(L,J) = ZPJQUAD(L,J) + VALJXX(L,JP)*JFACTA(J)
!          enddo
!        endif
!      enddo


     enddo !!!! end of SZA scan (if wanted)


     
!pjc      goto 92
!pjc   91 stop 'error in opening .dat file'
!pjc   92 stop

!$OMP END MASTER

   end subroutine cloudJ_interface

 end module UCI_cloudJ_interface

    
