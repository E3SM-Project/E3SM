!<<<<<<<<<<<<<<<<<<fastJX initialization codes:  need to be called only once
! Solar/Cloud/Fast-J v 7.6+

! INTERFACE:

      MODULE FJX_INIT_MOD

      USE FJX_CMN_MOD
      USE FJX_SUB_MOD, ONLY : EXITC
!SJ!  USE CMN_FJX_MOD
      use cam_history,      only : fieldname_len
      use cam_logfile,      only : iulog
      use cam_abortutils,   only : endrun
      use spmd_utils,       only : iam, masterproc
      use shr_kind_mod, only : r8 => shr_kind_r8,r4 => shr_kind_r4, shr_kind_cl, SHR_KIND_CS
      use shr_file_mod, only : shr_file_getunit, shr_file_freeunit

      IMPLICIT NONE

      logical, parameter :: LPRTJ2 = .FALSE. !Turn On/off diagnostic output to atm.log
      
      PUBLIC  :: INIT_FJX
      
      CONTAINS

!-----------------------------------------------------------------------
      subroutine INIT_FJX (TITLEJXX,NJXU,NJXX)
!-----------------------------------------------------------------------


      implicit none

      integer, intent(in)  ::NJXU
      integer, intent(out) ::NJXX
      character*6, intent(out), dimension(NJXU) :: TITLEJXX
      character*120  TIT_SPEC
      integer  JXUNIT,I, J, K, KR, RANSEED, NUN

      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string

      if (LPRTJ2) write(iulog,*) ' Solar/Cloud-J  ver-7.7 initialization'

! Use channel 8 to read fastJX data files:
!      JXUNIT  = 8
!      NUN = JXUNIT
      NUN = shr_file_getUnit()

!pjc! open (NUN,FILE='tables\CJ77_inp.dat',status='old',form='formatted')  !DOS
!pjc! open (NUN,FILE='tables/CJ77_inp.dat',status='old',form='formatted')  !Unix
      filepath = 'tables/CJ77_inp.dat' 
      open( unit=NUN, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J driver file (iostat='//ios_string//'): '//trim(filepath))
      endif

      read (NUN,'(a120)',iostat=ios) TIT_SPEC
      if (LPRTJ2) write(iulog,'(a)') TIT_SPEC
      call Check_Read_OK(ios)
      read (NUN,'(e10.3)',iostat=ios) RAD
      call Check_Read_OK(ios)
      read (NUN,'(e10.3)',iostat=ios) ZZHT
      call Check_Read_OK(ios)
      read (NUN,'(e10.3)',iostat=ios) ATAU
      call Check_Read_OK(ios)
      read (NUN,'(e10.3)',iostat=ios) ATAU0
      call Check_Read_OK(ios)
      read (NUN,'(e10.3)',iostat=ios) CLDCOR
      call Check_Read_OK(ios)
      read (NUN,'(i10  )',iostat=ios) NWBIN
      call Check_Read_OK(ios)
      read (NUN,'(i10  )',iostat=ios) NSBIN
      call Check_Read_OK(ios)
      read (NUN,'(i10  )',iostat=ios) LNRG
      call Check_Read_OK(ios)
      read (NUN,'(i10  )',iostat=ios) NRANDO
      call Check_Read_OK(ios)
      read (NUN,'(i10  )',iostat=ios) ATM0
      call Check_Read_OK(ios)
      read (NUN,'(i10  )',iostat=ios) CLDFLAG
      call Check_Read_OK(ios)

      if (LPRTJ2) write(iulog,'(a,3i5)') ' finish params LNRG ATM0 CLDFLAG',LNRG,ATM0,CLDFLAG
      close (NUN)
      if (LPRTJ2) write(iulog,'(a,i4)')'PJC 1: NUN =',NUN
      call shr_file_freeunit(NUN)
      NSJSUB(:) = 0
      SJSUB(:,:) = 0.d0   ! default set up for wavelengths when no sub-bins
      SJSUB(:,1) = 1.d0

! Cloud-J default with added near IR bins (if S_ > W)) but no sub bins
      LRRTMG = .false.
      LCLIRAD =.false.
      LGGLLNL =.false.
      NSJSUB(1:S_)= NGC(1:S_)

! while CLIRAD could be configured to run in Cloud-J, it cannot without custom fixes.
!sJ!      if (W_rrtmg .gt. 0) then
!sJ!! use RRTMG gas absorption/NGC is set at cmn_fjx_
!sJ!         NSJSUB(1:SX_)= NGC(1:SX_)
!sJ!         W_r=   W_rrtmg
!sJ!         LRRTMG  =.true.
!sJ!         LCLIRAD =.false.
!sJ!         LGGLLNL =.false.
!sJ!      elseif (W_clirad .gt. 0) then
!sJ!! use CLIRAD gas absorption
!sJ!         NSJSUB(1:W_)=1
!sJ!         NSJSUB(W_+1:SX_)=10
!sJ!         W_r = W_clirad
!sJ!         LRRTMG  =.false.
!sJ!         LCLIRAD =.true.
!sJ!         LGGLLNL =.false.
!sJ!      elseif (W_LLNL .gt. 0) then
!sJ!! use GG-LLNL gas absorption
!sJ!         NSJSUB(1:W_)=1
!sJ!         NSJSUB(W_+1:SX_)= 7
!sJ!         W_r = W_LLNL
!sJ!         LRRTMG  =.false.
!sJ!         LCLIRAD =.false.
!sJ!         LGGLLNL =.true.
!sJ!      endif
      if (LPRTJ2) write(iulog,'(a,3l2)')'LRRTMG/LCLIRAD/LGGLLNL=', LRRTMG, LCLIRAD, LGGLLNL

!  inital RRTMG setup is done in subrotine CHEM_IN of p-input.f
!  note that (if(LRRTMG) call RRTMG_SW_INI(cpdair)) is in CHEM_in of p-input.f
!  lock indexing of RRTMg superbins (1:W_+W_r) onto std bins fluxes (1:S_)
      if (LPRTJ2) write(iulog,'(a,i3,a,i3,a,i3,a,i3)') 'W_rrtmg= ',W_rrtmg,'  S_=',S_,'  W_r=',W_r,'  W_+ W_r= ',W_+W_r
      if (LPRTJ2) write(iulog,'(a,f8.4,a,f8.4,a,i2)') 'ATAU0=',ATAU0,'  ATAU=',ATAU,'   option(ATM0)= ', ATM0

! with Cloud-J v7.6, NO wavelength truncation for trop only, internal fixes remain
      if (W_ .ne. 18) then
        call EXITC(' INIT_JX: invalid no. wavelengths')
      endif

! set up angles of diffuse radiance at ocean surface
      ANGLES(1) = sngl(EMU(1))
      ANGLES(2) = sngl(EMU(2))
      ANGLES(3) = sngl(EMU(3))
      ANGLES(4) = sngl(EMU(4))
      ANGLES(5) = 0.e0 ! assgin U0 in photol_mod.f90

! Read in Fast/Solar-J X-sections (spectral data)
      call RD_XXX('tables/FJX_spec.dat')

      if (.not.(LRRTMG .or. LCLIRAD .or. LGGLLNL)) then
         do I = W_,  S_
            SJSUB(I,1)   = 1.d0
            SJSUB(I,2:16)= 0.d0
         enddo
      endif
      KDOKR(:)=0
      KR = 0
      do K = 1,S_
         do J = 1,NSJSUB(K)
            KR = KR+1
            KDOKR(KR) = K
            if (LPRTJ2) write(iulog,'(A,2I5)')'KR/KDOKR(KR)',KR, KDOKR(KR)
         enddo
      enddo
      if (KR .ne. W_+W_r) CALL EXITC('>>>error w/ RRTM sub bins: KDOKR')
!      IF (LPRTJ2) WRITE (iulog,*) 'PJC: KDOKR =',KDOKR
!      IF (LPRTJ2) WRITE (iulog,*) 'PJC: W_=',W_,'   W_r=',W_r 
      do KR = 1, W_+W_r
         K = KDOKR(KR)
         if (FL(K) .gt. 0.d0) then ! FL is read in call RD_XXX
            LDOKR(KR) = 1
         else
            LDOKR(KR) = 0
         endif
      enddo

! Read in cloud scattering data
      call RD_CLD('tables/FJX_scat-cld.dat')

! Read in strat sulf aerosols scattering data
      call RD_SSA('tables/FJX_scat-ssa.dat')

! Read in aerosols scattering data
      call RD_MIE('tables/FJX_scat-aer.dat')

! Read in UMich aerosol scattering data
      call RD_UM ('tables/FJX_scat-UMa.dat')

! Read in GEOMIP aerosol scattering data
      call RD_GEO ('tables/FJX_scat-geo.dat')

! Read in T & O3 climatology used to fill e.g. upper layers or if O3 not calc.
      call RD_PROF('tables/atmos_std.dat')

! Read in H2O and CH4 profiles for Solar-J
      call RD_TRPROF('tables/atmos_h2och4.dat')

! Read in zonal mean Strat-Sulf-Aerosol monthly data
      call RD_SSAPROF('tables/atmos_geomip.dat')

      NJXX = NJX
      do J = 1,NJXX
        TITLEJXX(J) = TITLEJX(J)
      enddo

! Read in photolysis rates used in chemistry code and mapping onto FJX J's
!---CTM call:  read in J-values names and link to fast-JX names
!     call RD_JS_JX('tables/FJX_j2j.dat', TITLEJXX,NJXX)
!     call RD_JS_JX('tables/pjc_modified/FJX_j2j_superfast.dat', TITLEJXX,NJXX)
      call RD_JS_JX('tables/pjc_modified/FJX_j2j_chemUCI.dat', TITLEJXX,NJXX)

!---for full ASAD:
!     call RD_JS('tables/ratj.d', TITLEJXX,NJXX,TSPECI,JPSPEC  &
!                ,MJVAL,TJVAL,MJX)

!---setup the random number sequence RAN4
      RANSEED = 66
      call RANSET (NRAN_,RAN4,RANSEED)

!      goto 1
!    4 continue
!        call EXITC(' FJX_INIT: error in read')
!    1 continue

      END SUBROUTINE INIT_FJX


!-----------------------------------------------------------------------
      subroutine RD_XXX(NAMFIL)
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections.
!
!>>>>NEW v-7.6+ added Solar-J bins for some to expand to S_
!     NOTE:  W_=18, use NWBIN 8,12,18 to zero flux for wavelengths only in strat
!>>>>NEW v-7.3  expanded input, full names & notes
!>>>>NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!-----------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (JX_spec.dat) >> j2 for fast-J2
!     NUN      Channel number for reading data file
!
!     NJX    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) ::  NAMFIL
      integer  I, J, JJ, K, IW, NQRD, LQ, NWWW, NSSS
      character*120  TIT_SPEC, TIT_J1N
      character*16 TIT_J1L
      character*6  TIT_J1S,TIT_J2S
      real*8  FWSUM
      integer :: NUN
      
      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string

      TQQ(:,:) = 0.d0

!----------spectral data----set for new format data------------------
!   note that X_ = max # Xsects read in
!           NJX = # fast-JX J-values derived from this (.le. X_)
      if (W_ .ne. 18) then
       call EXITC(' no. wavelengths wrong: W_ .ne. 18')
      endif


      NUN = shr_file_getUnit()

      filepath = trim(NAMFIL) 
      open( unit=NUN, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J file (iostat='//ios_string//'): '//trim(filepath))
      endif

!      open (NUN,FILE=trim(NAMFIL),status='old',form='formatted')

      read (NUN,'(a120)',iostat=ios) TIT_SPEC
      call Check_Read_OK(ios)
      read (NUN,*,iostat=ios)
      call Check_Read_OK(ios)
      read (NUN,'(i5,5x,i5)',iostat=ios) NWWW, NSSS
      call Check_Read_OK(ios)

      if (LPRTJ2) write(iulog,'(a120)') TIT_SPEC
      if (LPRTJ2) write(iulog,'(i5,A20,i5,A20)')  NWWW, ' photo-chem wl bins ', NSSS, ' solar heating bins '

      if (NWWW.gt.WX_ .or. NSSS.gt.SX_) then
       call EXITC(' WX_ or SX_ not large enough')
      endif

      NW1 = 1
      NW2 = NWWW
      NS1 = 1
      NS2 = NSSS

!----w-params:  1=w-eff  2=w-bins, 3=solar(photons), 4=solar(W/m2), 5=Y-PAR,  6=Rayleigh, 7=SJ sub-bins
      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
       call Check_Read_OK(ios)
        if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
      read (NUN,'(5x,6e10.3)',iostat=ios)    (WL(IW),IW=1,NSSS)
      call Check_Read_OK(ios)

      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
         if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
      read (NUN,'(5x,6e10.3)',iostat=ios)    (WBIN(IW),IW=1,NSSS)
      call Check_Read_OK(ios)

      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
         if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
      read (NUN,'(5x,6e10.3)',iostat=ios)    (FL(IW),IW=1,NSSS)
      call Check_Read_OK(ios)

      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
         if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
      read (NUN,'(5x,6e10.3)',iostat=ios)    (FW(IW),IW=1,NSSS)
      call Check_Read_OK(ios)
         FWSUM=0.d0
         do IW=1, NSSS
            FWSUM= FWSUM + FW(IW)
         enddo
         if (LPRTJ2) write(iulog,*) 'total Solar flux=', FWSUM

      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
         if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
      read (NUN,'(5x,6e10.3)',iostat=ios)    (FP(IW),IW=1,NSSS)
      call Check_Read_OK(ios)

      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
         if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
      read (NUN,'(5x,6e10.3)',iostat=ios)    (QRAYL(IW),IW=1,NSSS)
      call Check_Read_OK(ios)

!7 SJ-sub-bins
      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
         if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
      do I = 1,NSSS
        SJSUB(I,1)    = 1.0d0
        SJSUB(I,2:15) = 0.0d0
      enddo
      if ((NSSS .eq. 27) .or. (NSSS .eq. 21)) then
!SJ! this is different in SJ, allows for SJSBU(:,1:16)  ? RRTM
         do I= NWWW, NSSS  ! fraction of solar radiation for each sub-bin
            read  (NUN,'(5x,5f10.5)',iostat=ios) (SJSUB(I,IW),IW=1,15)
            call Check_Read_OK(ios)
               if (LPRTJ2) write(iulog,'(5x,5f10.6)') (SJSUB(I,IW),IW=1,15)
         enddo
      endif

!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
!---NB the O3 and q-O3-O1D are at different temperatures and cannot be combined
      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
          SQQ(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      call Check_Read_OK(ios)
      read (NUN,'(a6)',iostat=ios) TIT_J2S
      call Check_Read_OK(ios)
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
                 TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      call Check_Read_OK(ios)
      read (NUN,'(a6)',iostat=ios) TIT_J2S
      call Check_Read_OK(ios)
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
                 TQQ(3,1), (QO2(IW,3),IW=1,NWWW)
      call Check_Read_OK(ios)
         TITLEJX(1) = TIT_J1S
         TITLEJL(1) = TIT_J1L
         LQQ(1) = 3
         if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print

      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
          SQQ(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      call Check_Read_OK(ios)
      read (NUN,'(a6)',iostat=ios) TIT_J2S
      call Check_Read_OK(ios)
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
                 TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      call Check_Read_OK(ios)
      read (NUN,'(a6)',iostat=ios) TIT_J2S
      call Check_Read_OK(ios)
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
                 TQQ(3,2), (QO3(IW,3),IW=1,NWWW)
      call Check_Read_OK(ios)
        TITLEJX(2) = TIT_J1S
        TITLEJL(2) = TIT_J1L
        LQQ(2) = 3
        if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print

      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
          SQQ(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      call Check_Read_OK(ios)
      read (NUN,'(a6)',iostat=ios) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
                 TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      call Check_Read_OK(ios)
      read (NUN,'(a6)',iostat=ios) TIT_J2S
      call Check_Read_OK(ios)
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
                 TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)
      call Check_Read_OK(ios)
        TITLEJX(3) = TIT_J1S
        TITLEJL(3) = TIT_J1L
        LQQ(3) = 3
        if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print

!---Read remaining species:  X-sections at 1-2-3 T_s
!---read in 1 to 3 X-sects per J-value (JJ)
        JJ = 3
!-- read new Xsection block
    3 continue
      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
        if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
        if (TIT_J1S .eq. 'endofJ') goto 1
!---try to add a new Xsect
    2 continue
       JJ = JJ+1
       LQ = 1
         if (JJ .gt. X_) call EXITC(' RD_XXX: X_ not large enough')
       TITLEJX(JJ) = TIT_J1S
       TITLEJL(JJ) = TIT_J1L
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
          SQQ(JJ),TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
      call Check_Read_OK(ios)
        LQQ(JJ) = LQ
!try to read a 2nd Temperature or Pressure
      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
        if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
        if (TIT_J1S .eq. 'endofJ') goto 1
      if (TIT_J1S .eq. TITLEJX(JJ)) then
        LQ = 2
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
        TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
      call Check_Read_OK(ios)
        LQQ(JJ) = LQ
!try to read a 3rd Temperature or Pressure
      read (NUN,'(a6,1x,a16,1x,a120)',iostat=ios) TIT_J1S,TIT_J1L,TIT_J1N
      call Check_Read_OK(ios)
         if (LPRTJ2) write(iulog,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
         if (TIT_J1S .eq. 'endofJ') goto 1
       if (TIT_J1S .eq. TITLEJX(JJ)) then
        LQ = 3
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',iostat=ios)    &
        TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
      call Check_Read_OK(ios)
        LQQ(JJ) = LQ
       else
        goto 2
       endif
      else
        goto 2
      endif
      goto 3
    4 continue
        call EXITC(' RD_XXX: error in read')
    1 continue
      NJX = JJ

!---read in complete, process Xsects for reduced wavelengths (Trop-Only)
!---    possibly also for WACCM >200nm-only version.
!---TROP-ONLY (W_ = 12 or 8) then drop the strat Xsects (labeled 'x')
      if (NWBIN .eq. 12 .or. NWBIN .eq. 8) then
         if (LPRTJ2) write(iulog,'(a)')  &
              ' >>>TROP-ONLY reduced wavelengths, drop strat X-sects'
         JJ = 3
         do J = 4,NJX
            if (SQQ(J) .ne. 'x') then
!---collapse Xsects
               JJ = JJ+1
               if (JJ .lt. J) then
                  TITLEJX(JJ) = TITLEJX(J)
                  LQQ(JJ) = LQQ(J)
                  SQQ(JJ) = SQQ(J)
                  do LQ = 1,LQQ(J)
                     TQQ(LQ,JJ) = TQQ(LQ,J)
                     do IW = 1,NWWW
                        QQQ(IW,LQ,JJ) = QQQ(IW,LQ,J)
                     enddo
                  enddo
               endif
            endif
         enddo
         NJX = JJ
      endif

!print-----
      do J = 1,NJX
         if (LPRTJ2) write(iulog,'(a8,i5,2x,a6,2x,a16,2x,a1,i3,2x,3f6.1)') &
           'X-sects', J,TITLEJX(J),TITLEJL(J), SQQ(J),LQQ(J),(TQQ(I,J),I=1,LQQ(J))
      enddo
!---need to check that TQQ (= T(K) or p(hPa)) is monotonically increasing:
      do J = 1,NJX
         if ((LQQ(J) .eq. 3) .and. (TQQ(2,J) .ge. TQQ(3,J))) then
            call EXITC ('TQQ out of order')
         endif
         if ((LQQ(J) .eq. 2) .and. (TQQ(1,J) .ge. TQQ(2,J))) then
            call EXITC ('TQQ out of order')
         endif
      enddo

!---if FL(K) =0, then scattering skipped, method for dropping to 8 or 12 trop-only bins
      if (NWBIN .eq. 12) then
         do IW = 1,4
            FL(IW) = 0.d0
         enddo
         do IW = 9,10
            FL(IW) = 0.d0
         enddo
      endif
      if (NWBIN .eq. 8) then
         do IW = 1,4
            FL(IW) = 0.d0
         enddo
         FL(5) = FL(5) * 2.d0
         do IW = 6,11
            FL(IW) = 0.d0
         enddo
      endif

      close(NUN)
      if (LPRTJ2) write(iulog,'(a,i4)')'PJC 2: NUN =',NUN
      call shr_file_freeunit(NUN)

      END SUBROUTINE RD_XXX


!-----------------------------------------------------------------------
      subroutine RD_CLD(NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX ver 7.5
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat-cld.dat)
!     NUN      Channel number for reading data file
!     NCC      Number of categories for cloud scattering phase functions
!     QCC      Cloud scattering phase functions
!     WCC      5 Wavelengths for supplied phase functions
!     PCC      Phase function: first 8 terms of expansion
!     RCC      Effective radius associated with cloud type
!     SCC      Single scattering albedo
!     DCC      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) ::  NAMFIL
      integer :: NUN

      integer  I,J,K,L, JCC
      character*120 TITLE0
      real*8     GCCJ, XNDR,XNDI

      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string

      NUN = shr_file_getUnit()
      filepath = NAMFIL 
      open( unit=NUN, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J file (iostat='//ios_string//'): '//trim(filepath))
      endif

!      open (NUN,FILE=NAMFIL,status='old',form='formatted',iostat=ios)

        read (NUN,'(a80)',iostat=ios) TITLE0
        call Check_Read_OK(ios)
          if (LPRTJ2) write(iulog,'(a80)') TITLE0                                  !print----
        read (NUN,'(i4)',iostat=ios)  NCC
        call Check_Read_OK(ios)
        read (NUN,'(i4)',iostat=ios)  MCC
        call Check_Read_OK(ios)
          if (LPRTJ2) write(iulog,'(3i6,a)') NCC,MCC,SX_,' types of clouds & #Reff'   !print----
        read (NUN,*,iostat=ios)
        call Check_Read_OK(ios)
        read (NUN,*,iostat=ios)
        call Check_Read_OK(ios)
        read (NUN,*,iostat=ios)
        call Check_Read_OK(ios)
        read (NUN,*,iostat=ios)
        call Check_Read_OK(ios)
        read (NUN,*,iostat=ios)
        call Check_Read_OK(ios)

        do K = 1, NCC
           read (NUN,'(a12,f8.5)',iostat=ios) TITLCC(K),DCC(K)
           call Check_Read_OK(ios)
           if (LPRTJ2) write(iulog,'(a,i4,1x,a12,f8.5)') 'Cloud#',K,TITLCC(K),DCC(K)  !print----
           do J = 12, SX_
              do I = 1,MCC
                 read (NUN,'(i2, 1x, f5.2, f5.1, f7.1, f5.3, e8.1,f6.3,f8.5,7f6.3)',iostat=ios) &
                      JCC,WCC(J,K),RCC(I,K),GCC(I,K),XNDR,XNDI,                &
                      QCC(J,I,K),SCC(J,I,K),(PCC(L,J,I,K),L=2,8)
                 call Check_Read_OK(ios)
                 if (JCC .ne. J) goto 4
                 PCC(1,J,I,K) = 1.d0
              enddo
              read (NUN,*,iostat=ios)
              call Check_Read_OK(ios)
           enddo
        enddo

! replicate all cloud data for w < 295 nm from J=12 (= 295 nm), OK since mostly trop clouds.
        do K = 1,NCC
           do J = 1,11
              WCC(J,K) = WCC(12,K)
              do I = 1,MCC
                 QCC(J,I,K) = QCC(12,I,K)
                 SCC(J,I,K) = SCC(12,I,K)
                 do L = 1,8
                    PCC(L,J,I,K) = PCC(L,12,I,K)
                 enddo
              enddo
           enddo
        enddo

        goto 2

    4 continue
        call EXITC(' RD_CLD: error in read')

    2 continue
        close(NUN)
        if (LPRTJ2) write(iulog,'(a,i4)')'PJC 3: NUN =',NUN
        call shr_file_freeunit(NUN)

         if (LPRTJ2) write(iulog,'(a,2f9.5,i5)') ' ATAU/ATAU0',ATAU,ATAU0   !print----

      END SUBROUTINE RD_CLD


      subroutine RD_SSA(NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX ver 7.4
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat-ssa.dat)
!     NUN      Channel number for reading data file
!     NSS      Number of categories for cloud scattering phase functions
!     QSS      Cloud scattering phase functions
!     WSS      5 Wavelengths for supplied phase functions
!     PSS      Phase function: first 8 terms of expansion
!     RSS      Effective radius associated with cloud type:  Integ(r^3 dr)/Integ(r^2 dr)
!     GSS      Effective geometric cross section:  Integ(pi r^2 dr)
!     SSS      Single scattering albedo
!     DSS      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) ::  NAMFIL
      integer :: NUN

      integer  I, J, JSS, K, JCC, NSX_
      character*120 TITLE0
      real*8     WJSS,XNDR,XNDI

      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string

      NUN = shr_file_getUnit()

      filepath = NAMFIL
      open( unit=NUN, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J file (iostat='//ios_string//'): '//trim(filepath))
      endif

!      open (NUN,FILE=NAMFIL,status='old',form='formatted',iostat=ios)

      read (NUN,'(a120)',iostat=ios) TITLE0
      call Check_Read_OK(ios)
      if (LPRTJ2) write(iulog,'(a120)') TITLE0                              !print----
      read (NUN,*,iostat=ios)
      call Check_Read_OK(ios)
      read (NUN,'(i4,i4)',iostat=ios)  NSS, NSX_
      call Check_Read_OK(ios)
      read (NUN,*,iostat=ios)
      call Check_Read_OK(ios)
      if (LPRTJ2) write(iulog,'(i6,a)') NSS, ' types of strat sulf aerosols'
      do K = 1,NSS
!SJ! *** the SSA file for LCLIRAD in SJ has different format  '(a12, 3f8,4,...
         read (NUN,'(a10, 3f8.4, 2f8.1)',iostat=ios)    &
              TITLSS(K),RSS(K),GSS(K),DSS(K),TSS(K),WSS(K)
         call Check_Read_OK(ios)
         if (LPRTJ2) write(iulog,'(i4,1x,a12,2f10.4,2f8.1)') K,TITLSS(K),RSS(K),DSS(K),TSS(K),WSS(K)
         do J = 5, NSX_
            read(NUN,'(i2,2f8.4,e8.1,2f8.5,7f6.3)',iostat=ios)      &
                 JSS,WJSS,XNDR,XNDI,QSS(J,K),SSS(J,K),(PSS(I,J,K), I=2,8)
            call Check_Read_OK(ios)
            PSS(1,J,K) = 1.d0
         enddo
      enddo
! reproduce all SSA data for J=1:4 with J=5
      do K = 1,NSS
         do J = 1,4
            QSS(J,K) = QSS(5,K)
            SSS(J,K) = SSS(5,K)
            do I = 1,8
               PSS(I,J,K) = PSS(I,5,K)
            enddo
         enddo
      enddo
      goto 2

    4 continue
        call EXITC(' RD_SSA: error in read')

    2 continue
        close(NUN)
        if (LPRTJ2) write(iulog,'(a,i4)')'PJC 4: NUN =',NUN
        call shr_file_freeunit(NUN)

      END SUBROUTINE RD_SSA


!-----------------------------------------------------------------------
      subroutine RD_MIE(NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols scattering data set for fast-JX ver 7.3+
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     WAA      5 Wavelengths for the supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!     DAA      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) ::  NAMFIL
      integer :: NUN

      integer  I, J, K , JAA
      character*120 TITLE0
!      character*12 TITLAA(A_)   ! TITLAA: Title for scattering data    NEEDS to be in COMMON
      Character*12 TITLAAJ
      real*8   RAAJ, DAAJ

      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string

      NUN = shr_file_getUnit()

      filepath = NAMFIL
      open( unit=NUN, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J file (iostat='//ios_string//'): '//trim(filepath))
      endif
      
!      open (NUN,FILE=NAMFIL,status='old',form='formatted',iostat=ios)

      read (NUN,'(a120)',iostat=ios) TITLE0
      call Check_Read_OK(ios)
!print----
          if (LPRTJ2) write(iulog,'(a120)') TITLE0
        read (NUN,*)
        read (NUN,*)
      do J = 1, A_
        read (NUN,'(i4,1x,a12,1x,2f6.3,1x,a120)',iostat=ios) &
         JAA,TITLAAJ,RAAJ,DAAJ,TITLE0
       if (JAA.gt.0) then
         TITLAA(J) = TITLAAJ
         RAA(J) = RAAJ
         DAA(J) = DAAJ
        do K = 1, 5
         read (NUN,'(f4.0,f7.4,f7.4,7f6.3)',iostat=ios) &
          WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
          PAA(1,K,J) = 1.d0
        enddo
         NAA = J
!print----
          if (LPRTJ2) write(iulog,'(i5,1x,a12,1x,7f7.3,1x,a80)')   &
           J,TITLAAJ,RAAJ,DAAJ,(QAA(K,J),K=1,5),TITLE0
       else
          goto 2
       endif
      enddo
          goto 2

    4 continue
        call EXITC(' RD_MIE: error in read')
    2 continue
        close(NUN)
        if (LPRTJ2) write(iulog,'(a,i4)')'PJC 5: NUN =',NUN
        call shr_file_freeunit(NUN)

      END SUBROUTINE RD_MIE


!-----------------------------------------------------------------------
      subroutine RD_UM(NAMFIL)
!-----------------------------------------------------------------------
!-------UMich aerosol optical data for fast-JX (ver 6.1+)
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) ::  NAMFIL
      integer :: NUN

      integer  I, J, K, L
      character*78 TITLE0
      character*20 TITLUM(33)   ! TITLUM: Title for U Michigan aerosol data set

      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string

      NUN = shr_file_getUnit()

      filepath = NAMFIL
      open( unit=NUN, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J file (iostat='//ios_string//'): '//trim(filepath))
      endif
      
!      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(a78)') TITLE0
      if (LPRTJ2) write(iulog,*) 'UMichigan Aerosols', TITLE0
      read(NUN,'(5x,10f5.0)') WMM
      if (LPRTJ2) write(iulog,'(a,10f7.1)') ' UMIchigan aerosol wavelengths:',WMM

!---33 Different UM Aerosol Types:  SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4,
!---      FF00(0%BC), FF02, ...FF14(14%BC),  BB00, BB02, ...BB30(30%BC)
      do L=1,33
         read(NUN,'(a4)') TITLUM(L)
!---21 Rel Hum:    K=1=0%, =2=5%, ... =20=95%, =21=99%
         do K=1,21
!---6 wavelengths: J=1=200nm, 2=300nm, 3=400nm, (4'=550nm) 5=600nm, 6=1000nm
!---3 optic vars:  I=1=SSAlbedo,  =2=g,  =3=k-ext
            read(NUN,'(18f9.5)')  ((UMAER(I,J,K,L),I=1,3),J=1,6)
         enddo
      enddo
      close(NUN)
      if (LPRTJ2) write(iulog,'(a,i4)')'PJC 6: NUN =',NUN
      call shr_file_freeunit(NUN)

      if (LPRTJ2) write(iulog,'(a)') 'collapse UM wavelengths, drop 550 nm'
      WMM(4) = WMM(5)
      WMM(5) = WMM(6)
      do L=1,33
         do K=1,21
            do I=1,3
               UMAER(I,4,K,L) = UMAER(I,5,K,L)
               UMAER(I,5,K,L) = UMAER(I,6,K,L)
            enddo
         enddo
      enddo

      if (LPRTJ2) write(iulog,'(7(i5,1x,a4))') (L,TITLUM(L), L=1,33)

      END SUBROUTINE RD_UM


!-----------------------------------------------------------------------
      subroutine RD_PROF(NAMFIL)
!-----------------------------------------------------------------------
!  Routine to input T and O3 reference profiles 'atmos_std.dat'
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) ::  NAMFIL
      integer ::  NJ2
!
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
      real*8  OFAC, OFAK

      character*78 TITLE0

      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string
!
      NJ2 = shr_file_getUnit()
      filepath = NAMFIL
      open( unit=NJ2, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J file (iostat='//ios_string//'): '//trim(filepath))
      endif

!      open (NJ2,file=NAMFIL,status='old',form='formatted')

      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
!      if (LPRTJ2) write(iulog,'(1X,A)') TITLE0
      if (LPRTJ2) write(iulog,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11F7.1)') (T_REF(I,L,M), I=1,41)
! volume mixing ratio from 0km to 60 km for every 2km resolution in pressure altitude z*
        read (NJ2,'(3X,11F7.4)') (O_REF(I,L,M), I=1,31)
      enddo
      close (NJ2)
      if (LPRTJ2) write(iulog,'(a,i4)')'PJC 7: NJ2 =',NJ2
      call shr_file_freeunit(NJ2)

!  Extend climatology to 100 km
!  LREF =51 in cmn_fjx_mod.f90
      OFAC = exp(-2.d5/5.d5)
      do I = 32,LREF
        OFAK = OFAC**(I-31)
        do M = 1,NTMONS
        do L = 1,NTLATS
          O_REF(I,L,M) = O_REF(31,L,M)*OFAK
        enddo
        enddo
      enddo
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 42,LREF
        T_REF(I,L,M) = T_REF(41,L,M)
      enddo
      enddo
      enddo

!      close(NJ2)
!      if (LPRTJ2) write(iulog,'(a,i4)')'PJC 8: NJ2 =',NJ2
!      call shr_file_freeunit(NJ2)

 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')

      END SUBROUTINE RD_PROF


!-----------------------------------------------------------------------
      subroutine RD_TRPROF(NAMFIL)
!-----------------------------------------------------------------------
!  Routine to input H2O and CH4 reference profiles 'atmos_h2och4.dat'
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) ::  NAMFIL
      integer ::  NJ2
!
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216

      character*78 TITLE0

      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string
!
      NJ2 = shr_file_getUnit()
      filepath = NAMFIL
      open( unit=NJ2, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J file (iostat='//ios_string//'): '//trim(filepath))
      endif

!      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
!      if (LPRTJ2) write(iulog,'(1X,A)') TITLE0
      if (LPRTJ2) write(iulog,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
!        if (LPRTJ2) write(iulog,'(1X,I3,3X,I2)')   LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11E9.2)') (H2O_REF(I,L,M), I=1,31)
        read (NJ2,'(3X,11F9.2)') (CH4_REF(I,L,M), I=1,31)
!        if (LPRTJ2) write (iulog,'(3X,11E9.2)') (H2O_REF(I,L,M), I=1,31)
!        if (LPRTJ2) write(iulog,'(3X,11F9.2)') (CH4REF(I,L,M), I=1,31)

      enddo
      close (NJ2)
      if (LPRTJ2) write(iulog,'(a,i4)')'PJC 9: NJ2 =',NJ2
      call shr_file_freeunit(NJ2)

!  Extend climatology to 100 km
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 32,LREF
        H2O_REF(I,L,M) = H2O_REF(31,L,M)
        CH4_REF(I,L,M) = CH4_REF(31,L,M)
      enddo
      enddo
      enddo

!      close(NJ2)    ! Unnecessary?
 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')

      END SUBROUTINE RD_TRPROF


!-----------------------------------------------------------------------
      subroutine RD_JS_JX(NAMFIL,TITLEJX,NJX)
!-----------------------------------------------------------------------
!  Read 'FJX_j2j.dat' that defines mapping of fast-JX J's (TITLEJX(1:NJX))
!    onto the CTM reactions:  react# JJ, named T_REACT, uses fast-JX's JVMAP
!    including scaling factor JFACTA
!-----------------------------------------------------------------------
!---mapping variables stored in  block /jvchem/JFACTA,JIND,NRATJ,JLABEL,JVMAP
!           real*8  JFACTA(JVN_)          integer JIND(JVN_), NRATJ
!           character*50 JLABEL(JVN_)     character*6  JVMAP(JVN_)
!     JFACTA    multiplication factor for fast-JX calculated J
!     JLABEL    label(*50) of J-value used in the main chem model
!     JVMAP     label(*6) of J-value used to match with fast-JX J's
!     NRATJ     number of Photolysis reactions in CTM chemistry, derived here
!                   NRATJ must be .le. JVN_
!-----------------------------------------------------------------------
      implicit none
!
      character(*), intent(in)               ::  NAMFIL
      character*6, intent(in),dimension(NJX) ::  TITLEJX
      integer, intent(in)                    ::  NJX
      integer                                ::  NUNIT
      integer   J,JJ,K
      character*120 CLINE
      character*50 T_REACT
      character*6  T_FJX
      real*8 F_FJX

      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string

! Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
! The chemistry code title describes fully the reaction (a50)
! Blank (unfilled) chemistry J's are unmapped
! The number NRATJ is the last JJ readin that is .le. JVN
!   include fractional quantum yield for the fast-JX J's

      JLABEL(:) = '------'
      JVMAP(:) = '------'
      JFACTA(:) = 0.d0

      NUNIT = shr_file_getUnit()
      filepath = NAMFIL
      open( unit=NUNIT, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         if (LPRTJ2) write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J file (iostat='//ios_string//'): '//trim(filepath))
      endif
      
!      open (NUNIT,file=NAMFIL,status='old',form='formatted')

       read (NUNIT,'(a)') CLINE
         if (LPRTJ2) write(iulog,'(a)') CLINE
      do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       if (JJ .gt. JVN_) exit
        JLABEL(JJ) = T_REACT
        JFACTA(JJ) = F_FJX
        JVMAP(JJ) = T_FJX
        NRATJ = JJ
      enddo

      close(NUNIT)
      if (LPRTJ2) write(iulog,'(a,i4)')'PJC 10: NUNIT =',NUNIT
      call shr_file_freeunit(NUNIT)

!---Zero / Set index arrays that map Jvalue(j) onto rates
      do K = 1,NRATJ
         JIND(K) = 0
       do J = 1,NJX
        if (JVMAP(K) .eq. TITLEJX(J)) then
         JIND(K) = J
        endif
       enddo
      enddo

      if (LPRTJ2) write(iulog,'(a,i4,a)')' Photochemistry Scheme with',NRATJ,' J-values'
      do K=1,NRATJ
       if (JVMAP(K) .ne. '------' ) then
        J = JIND(K)
        if (J.eq.0) then
         if (LPRTJ2) write(iulog,'(i5,a50,f6.3,a,1x,a6)') K,JLABEL(K),JFACTA(K), &
               ' no mapping onto fast-JX',JVMAP(K)
        else
         if (LPRTJ2) write(iulog,'(i5,a50,f6.3,a,i4,1x,a6)') K,JLABEL(K),JFACTA(K), &
               ' mapped to FJX:',J,TITLEJX(J)
        endif
       endif
      enddo

!      close(NUNIT)
!      if (LPRTJ2) write(iulog,'(a,i4)')'PJC 11: NUNIT =',NUNIT
!      call shr_file_freeunit(NUNIT)

    END SUBROUTINE RD_JS_JX


!-----------------------------------------------------------------------
      subroutine RD_GEO(NAMFIL)
!-----------------------------------------------------------------------
!-------GEOMIP SSA scattering data set for fast-JX ver 7.5 ONLY RRTMG 27 bins
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat-geo.dat)
!     NUN      Channel number for reading data file
!     NGG      Number of sequentially increasing R-eff's for GEOMIP aerosols
!     RGG      Effective radius associated with cloud type
!     DGG      density (g/cm^3)
!     QGG      ratio optical to geometric X-section
!     SGG      Single scattering albedo
!     PGG      Phase function: first 8 terms of expansion
!-----------------------------------------------------------------------
      implicit none
      character(*), intent(in) ::  NAMFIL
      integer :: NUN

      integer  I, J, K
      character*120 TITLE0
      real*8     WGGJ,XNDR,XNDI,G1,G2,G3

      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string

      NUN = shr_file_getUnit()
      filepath = NAMFIL 
      open( unit=NUN, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J file (iostat='//ios_string//'): '//trim(filepath))
      endif
      
!      open (NUN,FILE=NAMFIL,status='old',form='formatted',iostat=ios)

      read (NUN,'(a120)',iostat=ios) TITLE0
      if (LPRTJ2) write(iulog,'(a120)') TITLE0                              !print----
      read (NUN,'(i4)')  NGG
      if (LPRTJ2) write(iulog,'(i6,a)')  NGG, ' Reff-s for GEO SSA'          !print----
      read (NUN,*)
      read (NUN,*)
      do K = 1,NGG
         read(NUN,'(10x,5f8.4)') RGG(K),G1,DGG(K),G2,G3
         if (LPRTJ2) write(iulog,'(i4,1x,3f8.4,2f8.1)') K,RGG(K),DGG(K), G1,G2,G3
         do J = 5, 27
            read (NUN,'(2x,2f8.4,e8.1,2f8.5,7f6.3)',iostat=ios) &
                 WGGJ,XNDR,XNDI,QGG(J,K),SGG(J,K),(PGG(I,J,K),I=2,8)
            PGG(1,J,K) = 1.d0
         enddo
      enddo
! reproduce all GEO SSA data for w < 202 nm (J=1:4)
      do K = 1,NGG
         do J = 1,4
            QGG(J,K) = QGG(5,K)
            SGG(J,K) = SGG(5,K)
            do I = 1,8
               PGG(I,J,K) = PGG(I,5,K)
            enddo
         enddo
      enddo
      goto 2

    4 continue
      call EXITC(' RD_GEO: error in read')

    2 continue
      close(NUN)
      if (LPRTJ2) write(iulog,'(a,i4)')'PJC 12: NUN =',NUN
      call shr_file_freeunit(NUN)

      END SUBROUTINE RD_GEO


!-----------------------------------------------------------------------
      subroutine RD_SSAPROF(NAMFIL)
!-----------------------------------------------------------------------
!  Routine to input SSA-GEO reference profiles for 'atmos_geomip.dat'
!      R_GEO = effective radius (microns)
!      X_GEO = mass fraction (1e-9 kg-H2SO4/kg-air)
!-----------------------------------------------------------------------
      implicit none
      character(*), intent(in) ::  NAMFIL
      integer ::  NJ2
!
      integer J,L,M
      character*78 TITLE0

      character(len=shr_kind_cl) :: filepath
      integer :: ios
      character(len=4) :: ios_string
!
      NJ2 = shr_file_getUnit()

      filepath = NAMFIL
      open( unit=NJ2, file=filepath, iostat=ios, status="OLD",form='formatted')
      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error opening cloud-J file (iostat='//ios_string//'): '//trim(filepath))
      endif
      
 !     open (NJ2,file=NAMFIL,status='old',form='formatted')

      read (NJ2,'(a)') TITLE0
         if (LPRTJ2) write(iulog,'(1x,a)') TITLE0
      read (NJ2,*)
      read (NJ2,*)
! only specify 19 pressure levels from 2.7 hPa to 340 hPa
      read (NJ2,'(19f7.2)') (P_GREF(L),L=1,19)
      if (LPRTJ2) write(iulog,'(19f7.2)')(P_GREF(L),L=1,19)
      read (NJ2,*)
      read (NJ2,*)
      read (NJ2,*)
! latitude bins 1:36 are Gauss, but approx 1.3953 + (J-33)*2.7906 deg
      read (NJ2,'(32f5.1)') (Y_GREF(L),L=1,32)
      if (LPRTJ2) write (iulog,'(32f5.1)') (Y_GREF(L),L=1,32)

      read (NJ2,*)
      read (NJ2,'(32f5.1)') (Y_GREF(L),L=64,33,-1)
      if (LPRTJ2) write (iulog,'(32f5.1)') (Y_GREF(L),L=33,64)
      read (NJ2,'(a)') TITLE0
         if (LPRTJ2) write(iulog,'(1x,a)') TITLE0
      do M = 1,12
          read (NJ2,*)
        do J = 1,64
          read (NJ2,'(11x,18f6.3)') (R_GREF(J,L,M), L=1,18)
                   R_GREF(J,19,M) = 0.d0
        enddo
      enddo
      read (NJ2,'(a)') TITLE0
         if (LPRTJ2) write(iulog,'(1x,a)') TITLE0
      do M = 1,12
          read (NJ2,*)
        do J = 1,64
          read (NJ2,'(11x,18f6.3)') (X_GREF(J,L,M), L=1,18)
                   X_GREF(J,19,M) = 0.d0
        enddo
      enddo
      read (NJ2,'(a)') TITLE0
         if (LPRTJ2) write(iulog,'(1x,a)') TITLE0
      do M = 1,12
          read (NJ2,*)
        do J = 1,64
          read (NJ2,'(11x,18f6.3)') (A_GREF(J,L,M), L=1,18)
                   A_GREF(J,19,M) = 0.d0
        enddo
      enddo

      close(NJ2)
      if (LPRTJ2) write(iulog,'(a,i4)')'PJC 13: NJ2 =',NJ2
      call shr_file_freeunit(NJ2)

      END SUBROUTINE RD_SSAPROF


!-----------------------------------------------------------------------
      SUBROUTINE RANSET (ND,RAN4L,ISTART)
!-----------------------------------------------------------------------
!  generates a sequence of real*4 pseudo-random numbers RAN4L(1:ND)
!     program RAN3 from Press, based on Knuth
      implicit none
      integer, parameter ::  MBIG=1000000000
      integer, parameter ::  MSEED=161803398
      integer, parameter ::  MZ=0
      real*4 , parameter ::  FAC=1.e-9
      integer,intent(in)    :: ND
      real*4, intent(out)   :: RAN4L(ND)
      integer,intent(inout) :: ISTART
      integer :: MA(55),MJ,MK,I,II,J,K,INEXT,INEXTP
      integer :: ios

!---initialization and/or fix of ISEED < 0
        MJ = MSEED - abs(ISTART)
        MJ = mod(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        do I=1,54
          II = mod(21*I,55)
          MA(II) = MK
          MK = MJ-MK
          if (MK.lt.MZ) then
            MK=MK+MBIG
          endif
          MJ = MA(II)
        enddo
        do K=1,4
         do I=1,55
           MA(I)=MA(I)-MA(1+MOD(I+30,55))
           if (MA(I) .lt. MZ) then
             MA(I) = MA(I)+MBIG
           endif
         enddo
        enddo
        INEXT = 0
        INEXTP = 31
        ISTART = 1
!---generate next ND pseudo-random numbers
      do J=1,ND
         INEXT = mod(INEXT,55) +1
         INEXTP = mod(INEXTP,55) +1
         MJ = MA(INEXT) - MA(INEXTP)
        if (MJ .lt. MZ) then
          MJ=MJ+MBIG
        endif
         MA(INEXT) = MJ
         RAN4L(J) = MJ*FAC
      enddo

      END SUBROUTINE RANSET


!-----------------------------------------------------------------------
      SUBROUTINE Check_Read_OK (ios)
!-----------------------------------------------------------------------
      integer,intent(in) :: ios
      character(len=4) :: ios_string

      if (ios /= 0) then
         write(ios_string, "(I4)") ios
         call endrun('Error reading from Cloud-J file (iostat='//ios_string//')')
      endif
      END SUBROUTINE Check_Read_OK


      END MODULE FJX_INIT_MOD
