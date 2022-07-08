!------------------------------------------------------------------------------
!   'cld_sub_mod.f90' for Cloud-J v7.7 (02/2020) - last change = fixes for MAXRAN
!------------------------------------------------------------------------------
! Determines how to treat clouds, fractional overlap, etc.

      MODULE CLD_SUB_MOD

      USE FJX_CMN_MOD
      USE FJX_SUB_MOD,  ONLY: PHOTO_JX, EXITC

      IMPLICIT NONE

      PUBLIC  :: CLOUD_JX

      CONTAINS

      SUBROUTINE CLOUD_JX (U0,SZA,RFL,SOLF,LPRTJ,PPP,ZZZ,TTT,HHH,DDD,  &
             RRR,OOO,CCC, LWP,IWP,REFFL,REFFI, CLDF,CLDCOR,CLDIW,      &
             AERSP,NDXAER,L1U,ANU,NJXU, VALJXX,SKPERD,SWMSQ,OD18,      &
             CLDFLAG,NRANDO,IRAN,LNRG,NICA, JCOUNT,LDARK,WTQCA)

!---Current recommendation for best average J's is
!     1) cloud decorellation w/ max-overlap blocks:  LNRG = 6 and CLDCOR = 0.33
!     2) cloud quadrature w/ average over eqch QCA:   CLDFALG = 7
!
!   Problem with correlation is that each new cloud layer generates 2x combinations
!      Thus the possibilities are for 2**Lcloudtop ICAs - that is too many (2**35)
!   Using  correlation lengths (3 km for tropical and high clouds, 1.5 km for stratus)
!      Choose 6 bins (of thickness = correl length) as Max-Overlap, then have
!      these bins be randomly or correlated with bins above
!   For now just assume these are random as with the other 2 max-ran groups above.
!
! GRP1 = 0 - 1.5km,  GRP2 = 1.5 - 3.5km,  GRP3 = 3.5 - 6km
! GRP4 =  6 - 9km,   GRP5 = 9 - 13km,     GRP6 = 13km+   (GRP7 = separate cirrus shields)
!
!Key Refs
! Kato, S., et al (2010), Relationships among cloud occurrence frequency, overlap,
!        and effective thickness derived from CALIPSO and CloudSat merged cloud vertical
!        profiles, J. Geophys. Res., 115, D00H28, doi:10.1029/2009JD012277.
! Pincus, R., et al. (2005), Overlap assumptions for assumed probability distribution
!        function cloud schemes in large-scale models, J. Geophys. Res., 110, D15S09,
!        doi:10.1029/2004JD005100.
! Oreopoulos, L., et al (2012) Radiative impacts of cloud heterogeneity and overlap
!        in an atmospheric General Circulation Model,  Atmos. Chem. Phys., 12, 9097-9111,
!        doi:10.5194/acp-12-9097-2012
! Naud, C., & A. D.  DelGenio (2006) Cloud Overlap Dependence on Atmospheric Dynamics,
!        16th ARM Science Team Meeting Proceedings, Albuquerque, NM, March 27 - 31, 2006.
!
!  CLOUD_JX is fractional cloud cover driver for subroutine PHOTO_JX
!    calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
!    needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
!    needs day-of-year for sun distance, SZA (does not directly use lat or long)
!
!--CLOUD_JX:   different cloud schemes (4:8 require max-ran overlap algorithm)
!       CLDFLAG = 1  :  Clear sky J's
!       CLDFLAG = 2  :  Averaged cloud cover
!       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
!       CLDFLAG = 4  :  NO LONGER VALID
!       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
!       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
!       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin) ***recommended
!       CLDFLAG = 8  :  Calcluate J's for ALL ICAs (up to 20,000 per cell!)

!--CLDIW = index for each lcoud layer:
!     = 0 = no cloud
!     = 1 = water cloud only
!     = 2 = ice cloud only
!     = 3 = liquid+ice cloud mix
!-----------------------------------------------------------------------
      implicit none
!---calling sequence variables
      integer, intent(in)                    :: L1U,ANU,NJXU, CLDFLAG,IRAN,NRANDO,LNRG
      real*8,  intent(in)                    :: U0,SZA,SOLF  !v7.7
      real*8,  intent(inout)                 :: CLDCOR       !v7.7
      real*8,  intent(in), dimension(5,W_+W_r) :: RFL
      logical, intent(in)                    :: LPRTJ
      real*8,  intent(in), dimension(L1U+1)  :: PPP,ZZZ
      real*8,  intent(in), dimension(L1U  )  :: TTT,HHH,DDD,RRR,OOO,CCC
      real*8,  intent(in), dimension(L1U  )  :: LWP,IWP,REFFL,REFFI
      real*8,  intent(in), dimension(L1U,ANU):: AERSP
      integer, intent(in), dimension(L1U,ANU):: NDXAER
      real*8,  intent(in), dimension(L1U  )  :: CLDF
      integer, intent(in), dimension(L1U  )  :: CLDIW
! reports out the JX J-values, upper level program converts to CTM chemistry J's
      real*8,  intent(out), dimension(L1U-1,NJXU):: VALJXX
      real*8,  intent(out), dimension(S_+2,L1U)  :: SKPERD
      real*8,  intent(out), dimension(6)         :: SWMSQ
      real*8,  intent(out), dimension(L1U)       :: OD18
      real*8,  intent(out), dimension(NQD_)      :: WTQCA
      integer, intent(out)                       :: NICA,JCOUNT
      logical, intent(out)                       :: LDARK
!-----------------------------------------------------------------------
      logical  LPRTJ0
      integer  I,II,J,K,L,M,N, LTOP, NRG,IRANX
      real*8   CLDFR, XRAN, FSCALE, QCAOD, WTRAN
      real*8,  dimension(L1U)     :: LWPX,IWPX,REFFLX,REFFIX
      real*8,  dimension(LWEPAR)  :: CLTL,CLTI, CLT,CLDX
      integer, dimension(LWEPAR)  :: NCLDF
      integer, dimension(9)       :: GBOT,GTOP,GLVL,GNR,GCMX  ! # max-verlap groups set at 9
      integer, dimension(9,CBIN_+1) :: GFNR
      real*8,  dimension(CBIN_)   ::  CFBIN
      real*8,  dimension(ICA_)    ::    WCOL,OCOL, OCDFS
      integer, dimension(ICA_)    :: ISORT
      real*8,  dimension(LWEPAR+1)   :: TCLD,TTCOL,SOLT,TAUG
      integer, dimension(NQD_)       :: NQ1,NQ2,NDXQS
      real*8,  dimension(L1U-1,NJXU) ::  VALJXXX
      real*8,  dimension(S_+2,L1U)  :: SKPERDD
      real*8,  dimension(6)         :: SWMSQQ
      real*8,  dimension(L1U)       :: OD18Q

!-----------------------------------------------------------------------
      LPRTJ0 = LPRTJ
      JCOUNT = 0
      NICA = 0
      do L = LWEPAR+1, L1U
         LWPX(L) = 0.d0
         IWPX(L) = 0.d0
         REFFLX(L) = 0.d0
         REFFIX(L) = 0.d0
      enddo
      VALJXX(:,:) = 0.d0    ! zero J's Heating R's in case LDARK is returned
      SKPERD(:,:) = 0.d0
      SWMSQ(:)    = 0.d0
      VALJXXX(:,:) = 0.d0   ! zero the PHOTOJ equivalents for wtd averaging
      SKPERDD(:,:) = 0.d0
      SWMSQQ(:)    = 0.d0
      OD18(:)      = 0.d0
      OD18Q(:)     = 0.d0

!---CLOUD_JX:   different cloud schemes
!-----------------------------------------------------------------------
      if (CLDFLAG.lt.1 .or. CLDFLAG.gt.8)then
         call EXITC ('>>>stop, incorrect cloud index')
      endif

!--------------------CLDFLAG =  1, 2, 3---------------------------------
      if (CLDFLAG.le.3) then  !  no ICAs needed
         if (CLDFLAG.eq.2) then
! 2 = average cloud cover
           do L = 1, LWEPAR
              CLDFR = CLDF(L)
              LWPX(L) = LWP(L) * CLDFR
              IWPX(L) = IWP(L) * CLDFR
              REFFLX(L) = REFFL(L)
              REFFIX(L) = REFFI(L)
           enddo

!-----------------------------------------------------------------------
         elseif (CLDFLAG.eq.3) then
! 3 = average cloud cover, adjust cloud fraction **3/2
            do L = 1, LWEPAR
               CLDFR = CLDF(L) * sqrt(CLDF(L))
               LWPX(L) = LWP(L) * CLDFR
               IWPX(L) = IWP(L) * CLDFR
               REFFLX(L) = REFFL(L)
               REFFIX(L) = REFFI(L)
            enddo

!-----------------------------------------------------------------------
         elseif (CLDFLAG.eq.1) then
! 1 = clear sky - no clouds
            do L = 1, LWEPAR
               LWPX(L) = 0.d0
               IWPX(L) = 0.d0
            enddo
         endif
!-----------------------------------------------------------------------

!----all above have only a single, simple call for fast_JX------------
         if(LPRTJ0) then
            write(6,'(2a)') ' cloud_J v7.7 Internal print: clouds = ',&
                         TITCLD(CLDFLAG)
         endif
!-----------------------------------------------------------------------
         call PHOTO_JX (U0,SZA,RFL,SOLF, LPRTJ0, PPP,ZZZ,TTT,HHH,       &
                  DDD,RRR,OOO,CCC, LWPX,IWPX,REFFLX,REFFIX,AERSP,       &
                  NDXAER, L1U,ANU,NJXU, VALJXX,SKPERD,SWMSQ,OD18, LDARK)
         if (.not.LDARK) then
            JCOUNT = JCOUNT + 1
         endif
!-----------------------------------------------------------------------

      else
!--------------------CLDFLAG =  5, 6, 7, 8------------------------------
!-----------------------------------------------------------------------
!  All below = CLDFLAG = 5:8 need to set up cloud overlap & ICAs
!-----------------------------------------------------------------------
         do L = 1, LWEPAR
           CLDX(L) = CLDF(L)
           CLT(L) = 0.d0
           CLTI(L) = 0.d0
           CLTL(L) = 0.d0
           LWPX(L) = LWP(L)
           IWPX(L) = IWP(L)
           REFFLX(L) = REFFL(L)
           REFFIX(L) = REFFI(L)
         enddo
! do cloud fraction binning here and rescale IWP/LWP to conserve layer WP


! generate approx cloud visible optical depths for quadrature and sorting
!   true wavelength dependence will be recalculated in PHOTO_JX
         do L = 1,LWEPAR
            if (REFFIX(L) .gt. 0.d0) then
              CLTI(L) = IWPX(L)*0.75d0*2.d0/(REFFIX(L)*0.917d0)
              CLT(L) = CLT(L) + CLTI(L)
            endif
            if (REFFLX(L) .gt. 0.d0) then
              CLTL(L) = LWPX(L)*0.75d0*2.1d0/REFFLX(L)
              CLT(L) = CLT(L) + CLTL(L)
            endif
         enddo
         LTOP  = LWEPAR
!-------------------------------------------------------------------------
!---Generate max-ran cloud overlap groups used for CLDFLAG = 4:8
!---CLT(cloud ice+liq OD) & IWPX & LWPX adjusted to quantized cld fr
!-------------------------------------------------------------------------
         call ICA_NR(CLDX,CLT,IWPX,LWPX,ZZZ, CLDIW,LTOP,LNRG,CBIN_,ICA_, &
             CFBIN,CLDCOR,NCLDF, GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA)

!---call ICA_ALL to generate the weight and cloud total OD of each ICA
!-------------------------------------------------------------------------
         call ICA_ALL(CLDX,CLT,LTOP,CBIN_,ICA_, CFBIN,     &
            CLDCOR,NCLDF,GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA,  WCOL,OCOL)

         if(LPRTJ0) then
            write(6,*) ' cloud-J v7.7  internal print:  #ICAs = ',NICA
         endif

!-----------------------------------------------------------------------
! 4 = average direct beam over all ICAs  DISCONTINUED
         if (CLDFLAG .eq. 4) then
            call EXITC(' CLD FLAG = 4 not allowed')
         endif

!-----------------------------------------------------------------------
! 5 = random pick of NRANDO(#) ICAs (selected based on fractional area)
         if (CLDFLAG .eq. 5) then

            if(LPRTJ0) then
               write(6,*) ' Average over random selection of ICAs:',NRANDO
            endif
            WTRAN = 1.d0/float(NRANDO)
            OCDFS(1) = WCOL(1)
            do I = 2,NICA
               OCDFS(I) = OCDFS(I-1) + WCOL(I)
            enddo
            do N=1,NRANDO
               IRANX = mod (IRAN+N-1, NRAN_) + 1
               XRAN = RAN4(IRANX)
               I = 1
               do while (XRAN .gt. OCDFS(I) .and. I .lt. NICA)
                  I = I+1
               enddo

               call ICA_III(CLDX,CLT,LTOP,CBIN_,ICA_, I, &
                  CLDCOR,NCLDF,GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, TTCOL)

!---zero out cloud water paths which are not in the selected random ICA
               do L = 1, LTOP
                  LWPX(L) = LWP(L)
                  IWPX(L) = IWP(L)
               enddo
               do L = 1,LTOP
                  if(TTCOL(L) .lt. 1.d-8) then
                     IWPX(L) = 0.d0
                     LWPX(L) = 0.d0
                  endif
               enddo

               if(LPRTJ0) then
                  write(6,'(a,2i6,f8.3)') ' pick random ICA:',N,I,OCOL(I)
                  do L = 1,LTOP
                     if (TTCOL(L) .ge. 1.d-8) then
                        write(6,'(i5,f9.4,2f8.3)') L,TTCOL(L),LWPX(L),IWPX(L)
                     endif
                  enddo
               endif
!-----------------------------------------------------------------------
               call PHOTO_JX (U0,SZA,RFL,SOLF, LPRTJ0, PPP,ZZZ,TTT,HHH,     &
                     DDD,RRR,OOO,CCC, LWPX,IWPX,REFFLX,REFFIX,AERSP,        &
                     NDXAER, L1U,ANU,NJXU, VALJXXX,SKPERDD,SWMSQQ,OD18Q, LDARK)
               if (.not.LDARK) then
                  JCOUNT = JCOUNT + 1
               endif
!-----------------------------------------------------------------------
               LPRTJ0 = .false.
               do J = 1,NJXU
                  do L = 1,L1U-1
                     VALJXX(L,J) = VALJXX(L,J) + VALJXXX(L,J)*WTRAN
                  enddo
               enddo
               do K = 1,S_+2
                  do L = 1,L1U
                     SKPERD(K,L) = SKPERD(K,L) + SKPERDD(K,L)*WTRAN
                  enddo
               enddo
               do M = 1,6
                  SWMSQ(M) = SWMSQ(M) + SWMSQQ(M)*WTRAN
               enddo
               do L = 1,L1U
                  OD18(L) = OD18(L) + OD18Q(L)*WTRAN
               enddo
            enddo
         endif

!-----------------------------------------------------------------------
! 6 = calculate quadrature QCAs, use up to 4 mid points
         if (CLDFLAG .eq. 6) then

            call ICA_QUD(WCOL,OCOL,LTOP,ICA_,NQD_,NICA, &
                         WTQCA, ISORT,NQ1,NQ2,NDXQS)

            if (LPRTJ0) then
               write(6,'(a)') ' quadrature QCAs(mid-pt): wt/range/index/OD'
               do N=1,NQD_
                  if (WTQCA(N).gt.0.d0) then
                     write(6,'(i5,f8.4,3i8,2f10.3)')   &
                     N,WTQCA(N),NQ1(N),NQ2(N),NDXQS(N),OCOL(ISORT(NDXQS(N)))
                  endif
               enddo
            endif
            do N = 1, NQD_
               if (WTQCA(N) .gt. 0.d0) then
                  I = ISORT(NDXQS(N))

                  call ICA_III(CLDX,CLT,LTOP,CBIN_,ICA_, I, CLDCOR,NCLDF, &
                     GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, TTCOL)

!---zero out cloud water paths which are not in the selected QCA
                  do L = 1, LTOP
                     LWPX(L) = LWP(L)
                     IWPX(L) = IWP(L)
                  enddo
                  do L = 1,LTOP
                     if (TTCOL(L) .lt. 1.d-8) then
                        IWPX(L) = 0.d0
                        LWPX(L) = 0.d0
                     endif
                  enddo
!-----------------------------------------------------------------------
              call PHOTO_JX (U0,SZA,RFL,SOLF, LPRTJ0, PPP,ZZZ,TTT,HHH, &
                   DDD,RRR,OOO,CCC, LWPX,IWPX,REFFLX,REFFIX,AERSP,     &
                   NDXAER, L1U,ANU,NJXU, VALJXXX,SKPERDD,SWMSQQ,OD18Q, LDARK)
                  if (.not.LDARK) then
                     JCOUNT = JCOUNT + 1
                  endif
!-----------------------------------------------------------------------
                  LPRTJ0 = .false.
                  do J = 1,NJXU
                     do L = 1,L1U-1
                        VALJXX(L,J) = VALJXX(L,J) + VALJXXX(L,J)*WTQCA(N)
                     enddo
                  enddo
                  do K = 1,S_+2
                     do L = 1,L1U
                        SKPERD(K,L) = SKPERD(K,L) + SKPERDD(K,L)*WTQCA(N)
                     enddo
                  enddo
                  do M = 1,6
                     SWMSQ(M) = SWMSQ(M) + SWMSQQ(M)*WTQCA(N)
                  enddo
                  do L = 1,L1U
                     OD18(L) = OD18(L) + OD18Q(L)*WTQCA(N)
                  enddo
               endif
            enddo
         endif

!-----------------------------------------------------------------------
! 7 = calculate quadrature atmosphere - average cloud within each QCA bin.
         if (CLDFLAG .eq. 7) then

            call ICA_QUD(WCOL,OCOL,LTOP,ICA_,NQD_,NICA, &
                         WTQCA, ISORT,NQ1,NQ2,NDXQS)

            if (LPRTJ0) then
               write(6,'(a)') ' quadrature QCAs(avg-cld): wt/range/index/OD'
               write(6,'(4f8.4)') WTQCA
               write(6,'(4i5)') NQ1
               write(6,'(4i5)') NQ2
            endif
            do N = 1, NQD_
               if (WTQCA(N) .gt. 0.d0) then
                  if (NQ2(N) .ge. NQ1(N)) then
                     IWPX(:) = 0.d0
                     LWPX(:) = 0.d0
                     QCAOD = 0.d0
                     do II = NQ1(N),NQ2(N)
                        I = ISORT(II)

                        call ICA_III(CLDX,CLT,LTOP,CBIN_,ICA_, I, CLDCOR,NCLDF, &
                           GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, TTCOL)

                        if (LPRTJ0) then
                           write(6,'(a,3i5,2f8.4,f9.3)') ' N(QCA)/II/I WCOL,OCOL',  &
                           N,II,I,WCOL(I),WTQCA(N),OCOL(I)
                        endif
                        do L = 1,LTOP
                           if (TTCOL(L) .gt. 1.d-8) then
                              IWPX(L) = IWPX(L) + IWP(L)*WCOL(I)
                              LWPX(L) = LWPX(L) + LWP(L)*WCOL(I)
                              QCAOD = QCAOD + TTCOL(L)*WCOL(I)
                           endif
                        enddo
                     enddo
                     do L = 1,LTOP
                        IWPX(L) = IWPX(L)/WTQCA(N)
                        LWPX(L) = LWPX(L)/WTQCA(N)
                     enddo
                     QCAOD = QCAOD/WTQCA(N)
                     if (LPRTJ0) then
                        write(6,'(a,i3,a,f10.5,f10.3)') &
                        'Quad Atmos Avg #',N,' wt, tot-OD:',WTQCA(N),QCAOD
                        write(6,'(a)') 'L / LWP / IWP'
                        do L=1,LTOP
                           if ((LWPX(L)+IWPX(L)) .gt. 1.d-8)  then
                              write(6,'(i4,2f10.3)') L,LWPX(L),IWPX(L)
                           endif
                        enddo
                     endif

!-----------------------------------------------------------------------
               call PHOTO_JX (U0,SZA,RFL,SOLF, LPRTJ0, PPP,ZZZ,TTT,HHH,     &
                    DDD,RRR,OOO,CCC, LWPX,IWPX,REFFLX,REFFIX,AERSP,         &
                    NDXAER, L1U,ANU,NJXU, VALJXXX,SKPERDD,SWMSQQ,OD18Q, LDARK)
                     if (.not.LDARK) then
                        JCOUNT = JCOUNT + 1
                     endif
!-----------------------------------------------------------------------
                     LPRTJ0 = .false.
                     do J = 1,NJXU
                        do L = 1,L1U-1
                           VALJXX(L,J) = VALJXX(L,J) + VALJXXX(L,J)*WTQCA(N)
                        enddo
                     enddo
                     do K = 1,S_+2
                        do L = 1,L1U
                           SKPERD(K,L) = SKPERD(K,L) + SKPERDD(K,L)*WTQCA(N)
                        enddo
                     enddo
                     do M = 1,6
                        SWMSQ(M) = SWMSQ(M) + SWMSQQ(M)*WTQCA(N)
                     enddo
                     do L = 1,L1U
                        OD18(L) = OD18(L) + OD18Q(L)*WTQCA(N)
                     enddo
                  endif
               endif
            enddo
         endif

!-----------------------------------------------------------------------
! 8 = average J's over all ICAs
         if (CLDFLAG .eq. 8) then

            if(LPRTJ0) then
               write(6,*) ' Average Js over all ICAs: I/ODcol/WTcol'
               write(6,'(i5,2f9.4)') (L,OCOL(L),WCOL(L), L=1,min(12,NICA-1))
               if (NICA.gt.12) write(6,'(a)') '. . .'
                  write(6,'(i5,2f9.4)') NICA,OCOL(NICA),WCOL(NICA)
            endif
            do I = 1, NICA
               call ICA_III(CLDX,CLT,LTOP,CBIN_,ICA_, I, CLDCOR,NCLDF, &
                  GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, TTCOL)
!---zero out cloud water paths which are not in the selected random ICA
               do L = 1, LTOP
                  LWPX(L) = LWP(L)
                  IWPX(L) = IWP(L)
               enddo
               do L = 1,LTOP
                  if(TTCOL(L) .lt. 1.d-8) then
                     IWPX(L) = 0.d0
                     LWPX(L) = 0.d0
                  endif
               enddo
!-----------------------------------------------------------------------
               call PHOTO_JX (U0,SZA,RFL,SOLF, LPRTJ0, PPP,ZZZ,TTT,HHH,     &
                    DDD,RRR,OOO,CCC, LWPX,IWPX,REFFLX,REFFIX,AERSP,         &
                    NDXAER, L1U,ANU,NJXU, VALJXXX,SKPERDD,SWMSQQ,OD18Q, LDARK)
               if (.not.LDARK) then
                  JCOUNT = JCOUNT + 1
               endif
!-----------------------------------------------------------------------
               LPRTJ0 = .false.
               do J = 1,NJXU
                  do L = 1,L1U-1
                     VALJXX(L,J) = VALJXX(L,J) + VALJXXX(L,J)*WCOL(I)
                  enddo
               enddo
               do K = 1,S_+2
                  do L = 1,L1U
                     SKPERD(K,L) = SKPERD(K,L) + SKPERDD(K,L)*WCOL(I)
                  enddo
               enddo
               do M = 1,6
                  SWMSQ(M) = SWMSQ(M) + SWMSQQ(M)*WCOL(I)
               enddo
               do L = 1,L1U
                  OD18(L) = OD18(L) + OD18Q(L)*WCOL(I)
               enddo
            enddo
         endif
!-----------------------------------------------------------------------

      endif

      END SUBROUTINE CLOUD_JX


!-----------------------------------------------------------------------
      SUBROUTINE ICA_NR(CLDF,CLTAU,IWPX,LWPX,ZZZ,CLDIW,LTOP,LNRG,CBIN_, &
            ICA_,CFBIN,CLDCOR,NCLDF, GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA)
!-----------------------------------------------------------------------
!---revised in v7.7 (02/2020) fixed MAX-RAN (#0 & #3) set CLDCOR=0 if need be
!---Read in the cloud fraction (CLDF), cloud OD (CLTAU), cloud index (CLDIW)
!---Derive max-corr cloud overlaps.  Set up the ICAs (Independent Column Atmos)
!   NCLDF(L) has value 1:CBIN_+1 = quantized cloud fraction (0:CBIN_).
!   NCLDF = 0 means NO clouds
!   CFBIN(J) = cloud fraction assumed for bin L=0:CBIN_
!    (e.g., 1(=0-5%) = 0.025)
!   CLTAU(J) = is readjusted for quantum bins to preserve CLDF*CLTAU

!---Definition of ICAs is carefully laid out with key parameters below:
!-----------------------------------------------------------------------
!   NICA = no. of ICAs
!   NRG = no. of sub-groups that are randomly overlapped with each other,
!     but are maximally overlapped among the contiguous layers in the sub-group.
!         Technically, NRG can equal LTOP generating up to 2**LTOP ICAs.
!   GBOT(G=1:NRG) = lower CTM layer of max-overlap group G
!   GTOP(G=1:NRG) = upper CTM layer of max-overlap group G
!          All layers, cloudy or clear are placed in one NRG group.
!   GNR(G=1:NRG) = no. of unique quantized cloud fractions in group G
!                  (.le.CBIN_)
!          Defines the number of uniques fractions in a maximally overlapped
!          group.
!   GFNR(G=1:NRG,1:GNR(G)) = cloud fraction quantum no (value = 0 to NCBIN)
!          Stores the specific cloud fractions counted in GNR.
!-----------------------------------------------------------------------
      implicit none

!---Cloud Cover parameters (in fjx_cmn_mod.f90)
!      integer, parameter ::  NQD_  = 4
!      integer, parameter ::  NRAN_ = 10007  ! dimension for random number
!      integer, parameter ::  CBIN_ = 10     ! # of quantized cld fraction bins
!  may need to reduce quantum number for LNRG6 to be CBIN_ = 10
!      integer, parameter ::  ICA_  = 20000  ! # of Indep Colm Atmospheres

!---Local Cloud Cover parameters
!---define break between randomly or decorrelated overlap groups
      integer, parameter  :: NG_BRK = 0
!---set up for correlated max-overlap groups based on observations
      integer, parameter ::  NRG6_ = 6
      real*8, dimension(NRG6_), parameter:: Zbin =                 &
          [0.d5, 1.5d5, 3.5d5, 6.0d5, 9.0d5, 13.d5]

      integer,intent(in) :: LTOP, LNRG, CBIN_, ICA_
      integer,intent(in),dimension(LTOP) :: CLDIW
      real*8, intent(in),dimension(LTOP) :: CLDF, ZZZ
      real*8, intent(inout)                 :: CLDCOR     !v7.7
      real*8, intent(inout),dimension(LTOP) :: CLTAU,IWPX,LWPX

      integer, intent(out) ::  NRG, NICA
      integer, intent(out), dimension(LTOP) :: NCLDF
      integer, intent(out), dimension(9) :: GBOT,GTOP,GLVL,GNR,GCMX
      integer, intent(out), dimension(9,CBIN_+1) :: GFNR
      real*8,  intent(out), dimension(CBIN_) :: CFBIN

      real*8   FBIN, FSCALE, CLF_MIN, CLF_MAX, FSCALE2
      integer                   ::  NRGX, NICAX
      integer, dimension(9)  :: GMIN,GMAX
      integer, dimension(CBIN_) :: NSAME
      integer  I,K,L,LL,N,NC, L1,L2,L3,  LCLTOP,LCIRRUS
      logical  L1GRP,L2GRP,L3GRP, L6GRP
!-----------------------------------------------------------------------

!---quantize cloud fractions into bins to avoid excessive calculations for
!---  nearly identical maximally overlapping cloud fractions.
!---  CBIN_=20 => N=0=[0%], N=1=[0.001-5%],N=2=[5-10%], .
!                 N=19=[90-95%],N=20=[95-100%]
!---assume the upper end of the range and renormalize in-cloud TAU to preserve
!     CLDF*TAU
        FBIN = CBIN_
      do K = 1,CBIN_
        CFBIN(K) = float(K) / FBIN
      enddo

!---quantize cloud fractions into bins & adjust to conserve TAU*CF
! round up/down <2%=>0% and >98%=>100%  works for FBIN = 10 or 25 (not 40)
      CLF_MIN = 0.02d0
      CLF_MAX = 1.0d0 - CLF_MIN

! clear out any small fraction clouds and quantize the cloud fraction as NCLDF
      do L = 1,LTOP
       if (CLDF(L).lt.CLF_MIN .or. CLDIW(L).eq.0) then
         NCLDF(L) = 0
         CLTAU(L) = 0.d0
         IWPX(L) = 0.d0
         LWPX(L) = 0.d0
       elseif (CLDF(L) .gt. CLF_MAX) then
         NCLDF(L) = FBIN
       else
          FSCALE2 = CLDF(L)*FBIN + 0.4999d0
          NCLDF(L) = max(FSCALE2,1.00001d0)
       endif
      enddo

!---find the cloud-top layer (water or ice) or identify clear sky
         LCLTOP = 0
      do L = 1,LTOP
        if (NCLDF(L).gt.0) then
          LCLTOP = L
        endif
      enddo
      if (LCLTOP .eq. 0) then
          NRG = 0
          NICA = 1
        goto 1
      endif

! rescale LWPX, IWPX, CLTAU
      do L = 1,LCLTOP
        if (NCLDF(L) .gt. 0) then
          FSCALE = CLDF(L) / CFBIN(NCLDF(L))
          CLTAU(L) = CLTAU(L) * FSCALE
          IWPX(L) = IWPX(L) * FSCALE
          LWPX(L) = LWPX(L) * FSCALE
        endif
      enddo

!---define maximally overlapping sub-groups by set levels (LNRG) or min
!   cloud fraction
      if (LNRG .eq. 0) then
!-----------------------------------------------------------------------------
!---Identify the maximally overlapped groups by breaking at a minimun NCLDF(L)
!-----------------------------------------------------------------------------
!---search from bottom to top, finding 1st level in group with cloud fraction
!    .ge. threshold, and then first level above that at which the cld fraction
!    is .lt. threshold. NRG = number of such groups.
        CLDCOR = 0.d0         ! v7.7
        L = 1
        NRG = 0
        do while (L.lt.LTOP)  ! v7.7
          if (NCLDF(L) .gt. NG_BRK) then
            NRG = NRG+1
           if (NRG.gt.9) exit
            GMIN(NRG) = L
            GMAX(NRG) = LCLTOP
            do LL = L+1,LCLTOP
!  look for first layer to drop below CLDF threshold = NGRBRK
              if (NCLDF(LL) .le. NG_BRK) then
                GMAX(NRG) = LL
                exit
              endif
            enddo
            L = GMAX(NRG)+1
          else
            L = L+1
          endif
        enddo
!!! v7.7  avoid gaps in the max overlap groups, set GLVL
        NRG = max(NRG,1)
        GMIN(1)   = 1
        GMAX(NRG) = LTOP
        do N = 2,NRG
           GMIN(N) = GMAX(N-1) + 1
        enddo
        do N=1,NRG
          GLVL(N) = N
        enddo
!!!

      elseif (LNRG .eq. 3) then
!-----------------------------------------------------------------------------
!---Alternative approach to fix a maximum of 3 random-overlap groups
!    (for L60/L57 CTM)
!---GRP=1 (if at all) is L=1:8 (surf to +1 km)
!---GRP=2 (if at all) is L=9 to last LWCloud
!---GRP=3 (if at all) is L=last-LWCld+1 to LTOP
!-----------------------------------------------------------------------------
        CLDCOR = 0.d0    ! v7.7
        L1 = 1
        L2 = 9
!----- L3-1 = uppermost liquid water cloud,  L3 = first of only ice-clouds
        L3 = L2
        do L = LCLTOP,L2,-1
          if (CLDIW(L).eq.1 .or. CLDIW(L).eq.3) then
            L3 = L+1
            exit
          endif
        enddo
           L1GRP = .false.
           L2GRP = .false.
           L3GRP = .false.
        do L = L1,L2-1
           L1GRP = L1GRP .or. (NCLDF(L).gt.0)
        enddo
        do L = L2,L3-1
           L2GRP = L2GRP .or. (NCLDF(L).gt.0)
        enddo
        do L = L3,LCLTOP
           L3GRP = L3GRP .or. (NCLDF(L).gt.0)
        enddo
          NRG = 0
        if (L1GRP) then
          NRG = NRG+1
          GMIN(NRG) = L1
          GMAX(NRG) = L2-1
        endif
        if (L2GRP) then
          NRG = NRG+1
          GMIN(NRG) = L2
          GMAX(NRG) = L3-1
        endif
        if (L3GRP) then
          NRG = NRG+1
          GMIN(NRG) = L3
          GMAX(NRG) = LCLTOP
        endif
           NRG = max(NRG,1)
          GMIN(1)   = 1
!!!v7.7  see above
          GMAX(NRG) = LTOP
        do N=1,NRG
          GLVL(N) = N
        enddo
!!!

      else
!-----------------------------------------------------------------------------
!---Newest recommended approach (v7.3) to use cloud correlation lengths
!   Problem with correlation is that each new cloud layer generates 2x combinations
!      Thus the possibilities are for 2**Lcloudtop ICAs - that is too many (2**35)
!   Using  correlation lengths (3 km for tropical and high clouds, 1.5 km for stratus)
!      Choose 6 bins (of thickness = correl length) as Max-Overlap, then have
!      these bins be randomly or correlated with bins above
!   For now just assume these are random as with the other 2 max-ran groups above.
!
! GRP1 = 0 - 1.5km,  GRP2 = 1.5 - 3.5km,  GRP3 = 3.5 - 6km
! GRP4 =  6 - 9km,   GRP5 = 9 - 13km,     GRP6 = 13km+
!
!Key Refs
!  Kato, S., et al (2010), Relationships among cloud occurrence frequency, overlap,
!        and effective thickness derived from CALIPSO and CloudSat merged cloud vertical
!        profiles, J. Geophys. Res., 115, D00H28, doi:10.1029/2009JD012277.
! Pincus, R., et al. (2005), Overlap assumptions for assumed probability distribution
!        function cloud schemes in large-scale models, J. Geophys. Res., 110, D15S09,
!        doi:10.1029/2004JD005100.
! Oreopoulos, L., et al (2012) Radiative impacts of cloud heterogeneity and overlap
!        in an atmospheric General Circulation Model,  Atmos. Chem. Phys., 12, 9097-9111,
!        doi:10.5194/acp-12-9097-2012
! Naud, C., & A. D.  DelGenio (2006) Cloud Overlap Dependence on Atmospheric Dynamics,
!        16th ARM Science Team Meeting Proceedings, Albuquerque, NM, March 27 - 31, 2006.
!
!---Find the levels in each of the NRG6_  altitude-defined groups
       do L = 1,LCLTOP
        do N = 2,NRG6_
         if (ZZZ(L)-ZZZ(1) .lt. Zbin(N)) then
           GMAX(N-1) = L
         endif
        enddo
       enddo
         GMIN(1) = 1
       do N = 2,NRG6_
         GMIN(N) = GMAX(N-1) + 1
       enddo
       GMAX(NRG6_) = LCLTOP

!---find out if there are any clouds in each of the 6 correlated groups
          NRG = 0
       do N = 1,NRG6_
             L6GRP = .false.
        do L = GMIN(N),GMAX(N)
          if (NCLDF(L) .gt. 0) then
             L6GRP = .true.
          endif
        enddo
        if (L6GRP) then
           NRG = NRG + 1
           GMIN(NRG) = GMIN(N)
           GMAX(NRG) = GMAX(N)
           GLVL(NRG) = N
        endif
       enddo

!---pull off cirrus shields from top MAX-GRP as separate MAX-GRP to
!        allow better correlation between cumulus towers below them
       if (NRG .gt. 0) then
          LCIRRUS = 0
        do L = GMAX(NRG),GMIN(NRG),-1
         if (NCLDF(L).gt.FBIN/2 .and. CLDIW(L).eq.2) then
          LCIRRUS = L
         endif
        enddo
        if (LCIRRUS .gt.GMIN(NRG)) then
!---split the uppermost MAX-GRP
          GMIN(NRG+1) = LCIRRUS
          GMAX(NRG+1) = GMAX(NRG)
          GLVL(NRG+1) = 7
          GMAX(NRG) = LCIRRUS-1
          NRG = NRG+1
        endif
       endif
!!! v7.7  avoid gaps in levels of the max groups
        NRG = max(NRG,1)
        GMIN(1)   = 1
        GMAX(NRG) = LTOP
        do N = 2,NRG
           GMIN(N) = GMAX(N-1) + 1
        enddo
!!!
      endif
!---finished selection of max-overlap groups

!---simplify groups if no clouds with NCLDF > NG_BRK
      GBOT(:) = 0
      GTOP(:) = 0
      if (NRG .eq. 0) then
        NRG = 1
        GBOT(1) = 1
        GTOP(1) = LTOP         ! v7.7
      else
!---assign levels between maximum overlap groups to group above.
        GBOT(1) = 1
        GTOP(1) = GMAX(1)
        do N=2,NRG
          GBOT(N) = max(GTOP(N-1)+1, GMIN(N))
!         GBOT(N) = min(GTOP(N-1)+1, GMIN(N))
          GTOP(N) = GMAX(N)
        enddo
        GTOP(NRG) = LTOP         ! v7.7
      endif
!---for each max-overlap group calculate number of unique cloud fractions
      do N = 1,NRG
        NSAME(:) = 0
        GCMX(N) = 0
        do L = GBOT(N),GTOP(N)
          if (NCLDF(L) .gt. 0) then
            NSAME(NCLDF(L)) = 1
            GCMX(N) = max(GCMX(N),NCLDF(L))
          endif
        enddo
!---sort cloud fractions in deceasing order for each max-overlap group
!---  note that largest bin N=CBIN_ (eg, 95-100%) will be treated as 100%
        GFNR(N,1) = CBIN_
        NC = 1
        do I = CBIN_-1,1,-1
          if(NSAME(I) .gt. 0) then
            NC = NC+1
            GFNR(N,NC) = I
          endif
        enddo
        GNR(N) = NC
        GFNR(N,NC+1) = 0
      enddo
!---number of unique columns in group:  if too many ICAs, drop upper groups!
      NICA = 1
      do N = 1,NRG
        NICA = NICA*GNR(N)
        if (NICA .le. ICA_) then
          NICAX = NICA
          NRGX = N
        endif
      enddo
      if (NICA .gt. ICA_) then
        write(6,*) 'NICA greater than ICA_',NICA,ICA_,NICAX,NRG,NRGX
        NICA = NICAX
        NRG = NRGX
      endif

    1 continue

      END SUBROUTINE ICA_NR


!-----------------------------------------------------------------------
      SUBROUTINE ICA_ALL(CLF,CLT,LTOP,CBINU,ICAU, CFBIN,CLDCOR,NCLDF,  &
           GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA,  WCOL,OCOL)
!-----------------------------------------------------------------------
!    OCOL() = cloud optical depth (total) in each ICA
!    WCOL() = weight(fract area) of ICA,
!    TTCOL() profile of cloud OD is not calculated here.
!    ISORT() = index no. of ICA sorted by column OD from smallest to largest
!---Using the information on max-ran cloud overlap generated by ICA_NR,
!---  this usbroutine generates all the ICAs
!   GBOT(I=1:NRG) = lower CTM layer of max-overlap group I
!   GTOP(I=1:NRG) = upper CTM layer of max-overlap group I
!   GNR(I=1:NRG) = no. of unique quantized cloud fractions in group I
!    (.le.NCBINS)
!   GFNR(I=1:NRG,1:GNR(I)) = cloud fraction quantum no (value = 1 to NCBIN)
!---See JL Neu, MJ Prather, JE Penner (2007), Global atmospheric chemistry:
!      Integrating over fractional cloud cover,J. Geophys. Res., 112, D11306,
!       doi:10.1029/2006JD008007
      implicit none
      integer, intent(in) :: LTOP, CBINU, ICAU, NRG, NICA
      integer, intent(in), dimension(LTOP)  :: NCLDF
      integer, intent(in), dimension(9)     :: GBOT,GTOP,GLVL,GNR,GCMX
      integer, intent(in), dimension(9,CBINU+1) :: GFNR
      real*8,  intent(in), dimension(LTOP)  :: CLF,CLT
      real*8,  intent(in), dimension(CBINU) :: CFBIN
      real*8,  intent(in)                   :: CLDCOR
      real*8,  intent(out),dimension(ICAU)  :: WCOL,OCOL

      real*8  ODCOL,WTCOL,CF0(51),  FWT(10,51),FWTC(10,51),FWTCC(10,51)
      real*8  FIG1,FIG2,GCORR,GCOWT,CORRFAC, FCMX(10) ,CLTOT(100)
      integer I, II, IG1,IG2, G, L,  IGNR(10),GCLDY(10),GRP1,GRP2
      logical L_CLR1,L_CLR2  ,LSKIP   ,LGR_CLR(10)
!-----------------------------------------------------------------------
        CLTOT(:) = 0.d0

        CF0(1) = 0.d0
      do L = 1,CBINU
        CF0(L+1) = CFBIN(L)
      enddo

      do G = 1,NRG
        FCMX(G) = CF0(GCMX(G)+1)           ! max cloud-fraction in MAX-GRP
       if (FCMX(G) .lt. 0.99d0) then
          LGR_CLR(G) = .true.          ! 1st member of MAX-GRP G = clear sky
          GCLDY(G) = 2
       else
          LGR_CLR(G) = .false.
          GCLDY(G) = 1
       endif
       do I = 1,GNR(G)        ! std weighting for each member of each MAX-Group
          FWT(G,I) = CF0(GFNR(G,I)+1) - CF0(GFNR(G,I+1)+1)
          FWTC(G,I) = FWT(G,I)
          FWTCC(G,I) = FWT(G,I)
       enddo
      enddo
        FCMX(NRG+1) = 0.d0

!  pre-calculate correl factors here:  no change if G = 100% cloud or G+1 = 100% cloud or clear
!   also no correlation fix for top MAX-GRP
      do G = 1,NRG-1
       LSKIP =   GCMX(G+1).eq.0 .or. GCMX(G+1).eq.CBINU   &
           .or.  GCMX(G).eq.0 .or. GCMX(G).eq.CBINU     !must have cloudy&clear in both MAX-GRPs
       if (.not.LSKIP) then
          FIG2 = FCMX(G+1)       ! cloudy fract of MAX-GRP just above (G+1)
          GRP2 = GLVL(G+1)       ! upper G6 group for NRG # G+1
          FIG1 = FCMX(G)         ! cloudy fract of current MAX-GRP (sum of cloudy fracts)
          GRP1 = GLVL(G)         ! current G6 group for NRG # G
          CORRFAC = CLDCOR**(GRP2-GRP1)  ! Cloud Correl Factor decreases with gap in G6 groups
! correlation factor: increase fract-area of cloudy member under cloudy section of upper layer (FIG2)
! Note that limits to increase depend on fract of cloud area above and the layer being increased
          GCORR = min(1.d0 + CORRFAC*(1.d0/FIG2 - 1.d0), 1.d0/FIG2, 1.d0/FIG1)   !v7.7 taken out of loop below
        do I = 2,GNR(G)
! enhance weighting for cloudy members below a cloud, reduce weighting below clear sky
         FWTC(G,I) = GCORR * FWT(G,I)
         FWTCC(G,I) = FWT(G,I)*(1.d0-GCORR*FIG2)/(1.d0-FIG2)
        enddo
         FWTC(G,1) = 1.d0 - FIG1*GCORR
         FWTCC(G,1) = 1.d0 - FIG1*(1.d0-GCORR*FIG2)/(1.d0-FIG2)
       endif
      enddo

      do I = 1,NICA
          WTCOL = 1.d0
          ODCOL = 0.d0
! for each ICA locate the members of each GROUP that contributes to it
          II = I
        do G = 1,NRG
          IGNR(G) = mod(II-1, GNR(G)) + 1
          II = (II-1)/GNR(G) + 1
        enddo
        do G = 1,NRG
           IG1 = IGNR(G)          ! working on MAX-GRP = G, member IG1
           L_CLR1 = GFNR(G,  IG1) .gt. GCMX(G)      ! member IG1 is clear
         if (G .eq. NRG) then     ! fix of indexing error in Cloud-J 7.3 did not affect results
           L_CLR2 = .true.
         else
           IG2 = IGNR(G+1)        ! member of MAX-GRP = G+1 for this ICA
           L_CLR2 = GFNR(G+1,IG2) .gt. GCMX(G+1)    ! member above (IG2) is clear
         endif
! all of these combinations should preserve the total weighting for layer member IG1
          if (.not.L_CLR2) then
! upper layer GRP member IG2 is a cloud layer (maybe one of several)
            if (.not.L_CLR1) then
! immediate GRP layer member IG1 is a cloudy one
               GCOWT = FWTC(G,IG1)
            else
! immediate layer GRP member IG1 is the clear one (if it exists)
               GCOWT = FWTC(G,1)
            endif
          else
! upper layer GRP member IG2 is a clear layer
            if (.not.L_CLR1) then
! immediate GRP layer member IG1 is a cloudy one
               GCOWT = FWTCC(G,IG1)
            else
! immediate layer GRP member IG1 is the clear one (if it exists)
               GCOWT = FWTCC(G,1)
             endif
          endif
            WTCOL = WTCOL*GCOWT
          do L = GBOT(G),GTOP(G)
            if (NCLDF(L) .ge. GFNR(G,IG1)) then
              ODCOL = ODCOL + CLT(L)
! could store the full 2-D array of atmospheres if needed:  TCOL(L,I) = CLT(L)
            endif
          enddo
        enddo
          WCOL(I) = WTCOL
          OCOL(I) = ODCOL
      enddo

      END SUBROUTINE ICA_ALL


!-----------------------------------------------------------------------
      SUBROUTINE ICA_III(CLF,CLT,LTOP,CBINU,ICAU, III, &
               CLDCOR,NCLDF, GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, TTCOL)
!-----------------------------------------------------------------------
!    see ICA_ALL, this subroutine picks out the ICA atmosphere #III
!      and loads the REFF/WPs for a FAST_JX calculation.
      implicit none
      integer, intent(in) :: LTOP, CBINU, ICAU, NRG, NICA, III
      integer, intent(in), dimension(LTOP)  :: NCLDF
      integer, intent(in), dimension(9)     :: GBOT,GTOP,GLVL,GNR,GCMX
      integer, intent(in), dimension(9,CBINU+1) :: GFNR
      real*8,  intent(in), dimension(LTOP)  :: CLF,CLT
      real*8,  intent(in)                   :: CLDCOR
      real*8,  intent(out),dimension(LTOP)  :: TTCOL

      integer II, IG, G, L
!-----------------------------------------------------------------------
         TTCOL(:) = 0.d0
      II = max(1, min(NICA,III))
      do G = 1,NRG
          IG = mod(II-1, GNR(G)) + 1
          II = (II-1)/GNR(G) + 1
        do L = GBOT(G),GTOP(G)
          if (NCLDF(L) .ge. GFNR(G,IG)) then
            TTCOL(L) = CLT(L)
          endif
        enddo
      enddo

      END SUBROUTINE ICA_III


!-----------------------------------------------------------------------
      SUBROUTINE ICA_QUD(WCOL,OCOL, LTOP,ICAU,NQDU,NICA, &
                         WTQCA, ISORT,NQ1,NQ2,NDXQS)
!-----------------------------------------------------------------------
!---Take the full set of ICAs and group into the NQD_ ranges of total OD
!---Create the Cumulative Prob Fn and select the mid-point ICA for each group
!---The Quad atmospheres have weights WTQCA
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in)        :: LTOP,ICAU,NQDU,NICA
      real*8,  intent(in), dimension(ICAU)      :: WCOL,OCOL

      real*8, intent(out), dimension(NQDU)      :: WTQCA
      integer, intent(out), dimension(ICAU)     :: ISORT
      integer, intent(out), dimension(NQDU)     :: NQ1,NQ2,NDXQS

      real*8,  dimension(ICA_) :: OCDFS, OCOLS
      integer I, II, J, L, N, N1, N2

      real*8, parameter:: OD_QUAD(4) =[0.5d0, 4.0d0, 30.d0, 1.d9]
!-----------------------------------------------------------------------
      ISORT(:) = 0
      WTQCA(:)  = 0.d0
      NDXQS(:) = 0

!---sort all the Indep Column Atmos (ICAs) in order of increasing column OD
!--- ISORT is the key, giving the ICA number from smallest to largest column OD
!--- OCOLS is the column OD sorted = OCOL(ISORT(I))
!--- OCDFS is the Cum.Prob.Fn. of the successive, sorted ICA
      if (NICA .eq. 1)  then
        ISORT(1) = 1
        OCOLS(1) = OCOL(1)
      else
        call HEAPSORT_A (NICA,OCOL,OCOLS,ISORT,ICA_)
      endif
        OCDFS(1) = WCOL(ISORT(1))
      do I = 2,NICA
        OCDFS(I) = OCDFS(I-1) + WCOL(ISORT(I))
      enddo
!---find beginning/end of quad range, note NQ2 < NQ1 means nothing in that range
          I = 1
      do N = 1,NQDU
       do while (OCOLS(I).lt.OD_QUAD(N) .and. I.le.NICA)
          I = I+1
       enddo
        NQ2(N) = I-1
      enddo
        NQ1(1) = 1
      do N = 2,NQDU
        NQ1(N) = NQ2(N-1) + 1
      enddo
!---define QCA wts from cum prob, pick middle ICA as representative
      do N = 1,NQDU
          N1 = NQ1(N)
          N2 = NQ2(N)
       if (N2 .ge. N1) then
          NDXQS(N) = (N1+N2)/2
         if (N1 .gt. 1) then
            WTQCA(N) = OCDFS(N2)-OCDFS(N1-1)
         else
            WTQCA(N) = OCDFS(N2)
         endif
       endif
      enddo

      END SUBROUTINE ICA_QUD


!-----------------------------------------------------------------------
!      SUBROUTINE ICA_DIRECT -- not used, finally removed in v7.7



!-----------------------------------------------------------------------
      SUBROUTINE HEAPSORT_A (N,A,AX,IX,ND)
!-----------------------------------------------------------------------
!  classic heapsort, sorts real*8 array A(N) into ASCENDING order,
!     places sorted array AX(N):   AX(1) .le. AX(N)
!     returns indexing IX(N) that records the location of A in sequence:
!           A(IX(J)) ==> AX(J), s.t. IX(1) = orig location of smallest A
!                           and IX(N) = original loc. of largest A
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in)  :: N, ND
      real*8, dimension(ND),intent(in)  :: A
      real*8, dimension(ND),intent(out) :: AX
      integer,dimension(ND),intent(out) :: IX
      integer :: I,J,L,IR,IA
      real*8 :: RA
!-----------------------------------------------------------------------
      do I = 1,N
        IX(I) = I
        AX(I) = A(I)
      enddo
      L  = N/2+1
      IR = N
   10 continue
      if (L .gt. 1) then
        L = L-1
        RA = AX(L)
        IA = IX(L)
      else
        RA = AX(IR)
        IA = IX(IR)
        AX(IR) = AX(1)
        IX(IR) = IX(1)
        IR = IR-1
        if (IR .eq. 1) then
          AX(1) = RA
          IX(1) = IA
          return
        endif
      endif
      I = L
      J = L+L
   20 continue
      if (J .le. IR) then
        if (J .lt. IR) then
          if (AX(J) .lt. AX(J+1)) then
            J = J+1
          endif
        endif
        if (RA .lt. AX(J)) then
          AX(I) = AX(J)
          IX(I) = IX(J)
          I = J
          J = J+J
        else
          J = IR+1
        endif
        goto 20
      endif
        AX(I) = RA
        IX(I) = IA
      goto 10

      END SUBROUTINE HEAPSORT_A


      END MODULE CLD_SUB_MOD
