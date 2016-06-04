!=================================================================================================================
 module module_mp_thompson_cldfra3

!module_mp_thompson_cldfra3 contains the subroutine cal_cldfra3 which calculates the cloud fraction as a function
!of relative humidity. The subroutine cal_cldfra3 has been tested in WRF 3.8 with the Thompson cloud microphysics
!and should not be used with other cloud microphysics schemes.

!subroutine cal_cldfra3 was originally copied from ./phys/module_radiation_driver.F from WRF version 3.8.
!Laura D. Fowler (laura@ucar.edu) / 2016-06-04.

 implicit
 private
 public:: cal_cldfra3


 contains


!=================================================================================================================

!+---+-----------------------------------------------------------------+
!..Cloud fraction scheme by G. Thompson (NCAR-RAL), not intended for
!.. combining with any cumulus or shallow cumulus parameterization
!.. scheme cloud fractions.  This is intended as a stand-alone for
!.. cloud fraction and is relatively good at getting widespread stratus
!.. and stratoCu without caring whether any deep/shallow Cu param schemes
!.. is making sub-grid-spacing clouds/precip.  Under the hood, this
!.. scheme follows Mocko and Cotton (1995) in applicaiton of the
!.. Sundqvist et al (1989) scheme but using a grid-scale dependent
!.. RH threshold, one each for land v. ocean points based on
!.. experiences with HWRF testing.
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+

      SUBROUTINE cal_cldfra3(CLDFRA, qv, qc, qi, qs,                    &
     &                 p,t,rho, XLAND, gridkm,                          &
     &                 ids,ide, jds,jde, kds,kde,                       &
     &                 ims,ime, jms,jme, kms,kme,                       &
     &                 its,ite, jts,jte, kts,kte)
!
      USE module_mp_thompson   , ONLY : rsif, rslf
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: ids,ide, jds,jde, kds,kde,                  &
     &                      ims,ime, jms,jme, kms,kme,                  &
     &                      its,ite, jts,jte, kts,kte

      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN):: qv,p,t,rho
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT):: qc,qi,qs
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN):: XLAND

      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT):: cldfra
#if defined(mpas)
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN):: gridkm
#else
      REAL, INTENT(IN):: gridkm
#endif

!..Local vars.
#if defined(mpas)
      REAL:: RH_00, RHI_max
      REAL, DIMENSION(its:ite,jts:jte):: RH_00L, RH_00O
#else
      REAL:: RH_00L, RH_00O, RH_00, RHI_max
#endif
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme):: qvsat
      INTEGER:: i,j,k
      REAL:: TK, TC, qvsi, qvsw, RHUM
      REAL, DIMENSION(kms:kme):: qvs1d, cfr1d, T1d,                     &
     &                           P1d, R1d, qc1d, qi1d, qs1d

      character*512 dbg_msg
      LOGICAL:: debug_flag

!+---+

!..First cut scale-aware. Higher resolution should require closer to
!.. saturated grid box for higher cloud fraction.  Simple functions
!.. chosen based on Mocko and Cotton (1995) starting point and desire
!.. to get near 100% RH as grid spacing moves toward 1.0km, but higher
!.. RH over ocean required as compared to over land.

#if defined(mpas)
      do j = jts,jte
      do i = its,ite
         RH_00L(i,j) = 0.839 + SQRT(1./(50.0+gridkm(i,j)*gridkm(i,j)*gridkm(i,j)*0.5))
         RH_00O(i,j) = 0.879 + SQRT(1./(100.0+gridkm(i,j)*gridkm(i,j)))
      enddo
      enddo
#else
      RH_00L = 0.839 + SQRT(1./(50.0+gridkm*gridkm*gridkm*0.5))
      RH_00O = 0.879 + SQRT(1./(100.0+gridkm*gridkm))
#endif

      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite
         RHI_max = 0.0
         CLDFRA(I,K,J) = 0.0

         if (qc(i,k,j)+qi(i,k,j) .gt. 1.E-8) then
            CLDFRA(I,K,J) = 1.0
            qvsat(i,k,j) = qv(i,k,j)
         else

            TK   = t(i,k,j)
            TC   = TK - 273.16

            qvsw = rslf(P(i,k,j), TK)
            qvsi = rsif(P(i,k,j), TK)

            if (tc .ge. -12.0) then
               qvsat(i,k,j) = qvsw
            elseif (tc .lt. -30.0) then
               qvsat(i,k,j) = qvsi
            else
               qvsat(i,k,j) = qvsw - (qvsw-qvsi)*(-12.0-tc)/(-12.0+30.)
            endif
            RHUM = qv(i,k,j)/qvsat(i,k,j)

#if defined(mpas)
            IF ((XLAND(I,J)-1.5).GT.0.) THEN                             !--- Ocean
               RH_00 = RH_00O(i,j)
            ELSE                                                         !--- Land
               RH_00 = RH_00L(i,j)
            ENDIF
#else
            IF ((XLAND(I,J)-1.5).GT.0.) THEN                             !--- Ocean
               RH_00 = RH_00O
            ELSE                                                         !--- Land
               RH_00 = RH_00L
            ENDIF
#endif

            if (tc .ge. -12.0) then
               RHUM = MIN(0.999, RHUM)
               CLDFRA(I,K,J) = MAX(0.0, 1.0-SQRT((1.0-RHUM)/(1.-RH_00)))
#if defined(mpas)
            elseif (tc.lt.-12..and.tc.gt.-70. .and. RHUM.gt.RH_00O(i,j)) then
               RHI_max = MAX(RHUM+1.E-6, qvsw/qvsi)
               CLDFRA(I,K,J) = MAX(0., ((RH_00O(i,j)-RHUM)/(RH_00O(i,j)-RHI_max)) &
     &                                *((RH_00O(i,j)-RHUM)/(RH_00O(i,j)-RHI_max)))
#else
            elseif (tc.lt.-12..and.tc.gt.-70. .and. RHUM.gt.RH_00O) then
               RHI_max = MAX(RHUM+1.E-6, qvsw/qvsi)
               CLDFRA(I,K,J) = MAX(0., ((RH_00O-RHUM)/(RH_00O-RHI_max)) &
     &                                *((RH_00O-RHUM)/(RH_00O-RHI_max)))
#endif
            endif

            CLDFRA(I,K,J) = MAX(0.0, MIN(CLDFRA(I,K,J), 0.95))

         endif
      ENDDO
      ENDDO
      ENDDO

!..Prepare for a 1-D column to find various cloud layers.

      DO j = jts,jte
      DO i = its,ite
!        if (i.gt.10.and.i.le.20 .and. j.gt.10.and.j.le.20) then
!          debug_flag = .true.
!        else
!           debug_flag = .false.
!        endif
         DO k = kts,kte
            qvs1d(k) = qvsat(i,k,j)
            cfr1d(k) = cldfra(i,k,j)
            T1d(k) = t(i,k,j)
            P1d(k) = p(i,k,j)
            R1d(k) = rho(i,k,j)
            qc1d(k) = qc(i,k,j)
            qi1d(k) = qi(i,k,j)
            qs1d(k) = qs(i,k,j)
         ENDDO

!     if (debug_flag) then
!       WRITE (dbg_msg,*) 'DEBUG-GT: finding cloud layers at point  (', i, ', ', j, ')'
!       CALL wrf_debug (150, dbg_msg)
!     endif
         call find_cloudLayers(qvs1d, cfr1d, T1d, P1d, R1d,             &
     &                         debug_flag, qc1d, qi1d, qs1d, kts,kte)

         DO k = kts,kte
            cldfra(i,k,j) = cfr1d(k)
            qc(i,k,j) = qc1d(k)
            qi(i,k,j) = qi1d(k)
         ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE cal_cldfra3

!+---+-----------------------------------------------------------------+
!..From cloud fraction array, find clouds of multi-level depth and compute
!.. a reasonable value of LWP or IWP that might be contained in that depth,
!.. unless existing LWC/IWC is already there.

      SUBROUTINE find_cloudLayers(qvs1d, cfr1d, T1d, P1d, R1d,          &
     &                            debugfl, qc1d, qi1d, qs1d, kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: kts, kte
      LOGICAL, INTENT(IN):: debugfl
      REAL, DIMENSION(kts:kte), INTENT(IN):: qvs1d,T1d,P1d,R1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: cfr1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: qc1d, qi1d, qs1d

!..Local vars.
      REAL, DIMENSION(kts:kte):: theta, dz
      REAL:: Z1, Z2, theta1, theta2, ht1, ht2
      INTEGER:: k, k2, k_tropo, k_m12C, k_m40C, k_cldb, k_cldt, kbot
      LOGICAL:: in_cloud
      character*512 dbg_msg

!+---+

      k_m12C = 0
      k_m40C = 0
      DO k = kte, kts, -1
         theta(k) = T1d(k)*((100000.0/P1d(k))**(287.05/1004.))
         if (T1d(k)-273.16 .gt. -40.0) k_m40C = MAX(k_m40C, k)
         if (T1d(k)-273.16 .gt. -12.0) k_m12C = MAX(k_m12C, k)
      ENDDO
      if (k_m40C .le. kts) k_m40C = kts
      if (k_m12C .le. kts) k_m12C = kts

      Z2 = 44307.692 * (1.0 - (P1d(kte)/101325.)**0.190)
      DO k = kte-1, kts, -1
         Z1 = 44307.692 * (1.0 - (P1d(k)/101325.)**0.190)
         dz(k+1) = Z2 - Z1
         Z2 = Z1
      ENDDO
      dz(kts) = dz(kts+1)

!..Find tropopause height, best surrogate, because we would not really
!.. wish to put fake clouds into the stratosphere.  The 10/1500 ratio
!.. d(Theta)/d(Z) approximates a vertical line on typical SkewT chart
!.. near typical (mid-latitude) tropopause height.  Since messy data
!.. could give us a false signal of such a transition, do the check over 
!.. three K-level change, not just a level-to-level check.  This method
!.. has potential failure in arctic-like conditions with extremely low
!.. tropopause height, as would any other diagnostic, so ensure resulting
!.. k_tropo level is above 4km.

      DO k = kte-3, kts, -1
         theta1 = theta(k)
         theta2 = theta(k+2)
         ht1 = 44307.692 * (1.0 - (P1d(k)/101325.)**0.190)
         ht2 = 44307.692 * (1.0 - (P1d(k+2)/101325.)**0.190)
         if ( (((theta2-theta1)/(ht2-ht1)) .lt. 10./1500. ) .AND.       &
     &                       (ht1.lt.19000.) .and. (ht1.gt.4000.) ) then 
            goto 86
         endif
      ENDDO
 86   continue
      k_tropo = MAX(kts+2, k+2)

!     if (debugfl) then
!     print*, ' FOUND TROPOPAUSE ', k_tropo, ' near ', ht2, ' m'
!       WRITE (dbg_msg,*) 'DEBUG-GT: FOUND TROPOPAUSE ', k_tropo, ' near ', ht2, ' m'
!       CALL wrf_debug (150, dbg_msg)
!     endif

!..Eliminate possible fractional clouds above supposed tropopause.
      DO k = k_tropo+1, kte
         if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.0.999) then
            cfr1d(k) = 0.
         endif
      ENDDO

!..We would like to prevent fractional clouds below LCL in idealized
!.. situation with deep well-mixed convective PBL, that otherwise is
!.. likely to get clouds in more realistic capping inversion layer.

      kbot = kts+2
      DO k = kbot, k_m12C
         if ( (theta(k)-theta(k-1)) .gt. 0.05E-3*dz(k)) EXIT
      ENDDO
      kbot = MAX(kts+1, k-1)
      DO k = kts, kbot
         if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.0.999) cfr1d(k) = 0.
      ENDDO


!..Starting below tropo height, if cloud fraction greater than 1 percent,
!.. compute an approximate total layer depth of cloud, determine a total 
!.. liquid water/ice path (LWP/IWP), then reduce that amount with tuning 
!.. parameter to represent entrainment factor, then divide up LWP/IWP
!.. into delta-Z weighted amounts for individual levels per cloud layer. 

      k_cldb = k_tropo
      in_cloud = .false.
      k = k_tropo
      DO WHILE (.not. in_cloud .AND. k.gt.k_m12C)
         k_cldt = 0
         if (cfr1d(k).ge.0.01) then
            in_cloud = .true.
            k_cldt = MAX(k_cldt, k)
         endif
         if (in_cloud) then
            DO k2 = k_cldt-1, k_m12C, -1
               if (cfr1d(k2).lt.0.01 .or. k2.eq.k_m12C) then
                  k_cldb = k2+1
                  goto 87
               endif
            ENDDO
 87         continue
            in_cloud = .false.
         endif
         if ((k_cldt - k_cldb + 1) .ge. 2) then
!     if (debugfl) then
!           print*, 'An ice cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       WRITE (dbg_msg,*) 'DEBUG-GT: An ice cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       CALL wrf_debug (150, dbg_msg)
!     endif
            call adjust_cloudIce(cfr1d, qi1d, qs1d, qvs1d, T1d,R1d,dz,  &
     &                           k_cldb,k_cldt,kts,kte)
            k = k_cldb
         endif
         k = k - 1
      ENDDO


      k_cldb = k_tropo
      in_cloud = .false.
      k = k_m12C
      DO WHILE (.not. in_cloud .AND. k.gt.kbot)
         k_cldt = 0
         if (cfr1d(k).ge.0.01) then
            in_cloud = .true.
            k_cldt = MAX(k_cldt, k)
         endif
         if (in_cloud) then
            DO k2 = k_cldt-1, kbot, -1
               if (cfr1d(k2).lt.0.01 .or. k2.eq.kbot) then
                  k_cldb = k2+1
                  goto 88
               endif
            ENDDO
 88         continue
            in_cloud = .false.
         endif
         if ((k_cldt - k_cldb + 1) .ge. 2) then
!     if (debugfl) then
!           print*, 'A water cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       WRITE (dbg_msg,*) 'DEBUG-GT: A water cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       CALL wrf_debug (150, dbg_msg)
!     endif
            call adjust_cloudH2O(cfr1d, qc1d, qvs1d, T1d,R1d,dz,        &
     &                           k_cldb,k_cldt,kts,kte)
            k = k_cldb
         endif
         k = k - 1
      ENDDO

!..Do a final total column adjustment since we may have added more than 1mm
!.. LWP/IWP for multiple cloud decks.

      call adjust_cloudFinal(cfr1d, qc1d, qi1d, R1d,dz, kts,kte,k_tropo)

!     if (debugfl) then
!     print*, ' Made-up fake profile of clouds'
!     do k = kte, kts, -1
!        write(*,'(i3, 2x, f8.2, 2x, f9.2, 2x, f6.2, 2x,  f15.7, 2x, f15.7)') &
!    &        K, T1d(k)-273.15, P1d(k)*0.01, cfr1d(k)*100., qc1d(k)*1000.,qi1d(k)*1000.
!     enddo
!       WRITE (dbg_msg,*) 'DEBUG-GT:  Made-up fake profile of clouds'
!       CALL wrf_debug (150, dbg_msg)
!       do k = kte, kts, -1
!          write(dbg_msg,'(f8.2, 2x, f9.2, 2x, f6.2, 2x,  f15.7, 2x, f15.7)') &
!    &          T1d(k)-273.15, P1d(k)*0.01, cfr1d(k)*100., qc1d(k)*1000.,qi1d(k)*1000.
!          CALL wrf_debug (150, dbg_msg)
!       enddo
!     endif

      END SUBROUTINE find_cloudLayers

!+---+-----------------------------------------------------------------+

      SUBROUTINE adjust_cloudIce(cfr,qi,qs,qvs, T,Rho,dz, k1,k2,kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: k1,k2, kts,kte
      REAL, DIMENSION(kts:kte), INTENT(IN):: cfr, qvs, T, Rho, dz
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: qi, qs
      REAL:: iwp, xiwp, max_iwp, tdz, this_iwp, iwp_exists
      INTEGER:: k
      REAL, PARAMETER:: entr=0.35

      max_iwp = ABS(qvs(k2-1)-qvs(k2))*Rho(k2-1)*dz(k2-1)

      tdz = 0.
      iwp = 0.
      iwp_exists = 0.
      do k = k1, k2
         tdz = tdz + dz(k)
         iwp = iwp + MAX(0., (qvs(k-1)-qvs(k))*Rho(k)*dz(k))
         iwp_exists = iwp_exists + (qi(k)+qs(k))*Rho(k)*dz(k)
      enddo
      if (iwp_exists .gt. 1.0) RETURN
      max_iwp = MAX(max_iwp*(1.-entr), MIN(1.0, iwp*(1.0-entr)))

      do k = k1, k2
         this_iwp = max_iwp*dz(k)/tdz
         if (cfr(k).gt.0.01.and.cfr(k).lt.0.99.and.T(k).ge.203.16) then
            qi(k) = qi(k) + cfr(k)*cfr(k)*this_iwp/Rho(k)/dz(k)
         endif
      enddo

      END SUBROUTINE adjust_cloudIce

!+---+-----------------------------------------------------------------+

      SUBROUTINE adjust_cloudH2O(cfr, qc, qvs, T,Rho,dz, k1,k2,kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: k1,k2, kts,kte
      REAL, DIMENSION(kts:kte):: cfr, qc, qvs, T, Rho, dz
      REAL:: lwp, xlwp, max_lwp, tdz, this_lwp, lwp_exists
      INTEGER:: k
      REAL, PARAMETER:: entr=0.35

      max_lwp = ABS(qvs(k2-1)-qvs(k2))*Rho(k2-1)*dz(k2-1)

      tdz = 0.
      lwp = 0.
      lwp_exists = 0.
      do k = k1, k2
         tdz = tdz + dz(k)
         lwp = lwp + MAX(0., (qvs(k-1)-qvs(k))*Rho(k)*dz(k))
         lwp_exists = lwp_exists + qc(k)*Rho(k)*dz(k)
      enddo
      if (lwp_exists .gt. 1.0) RETURN
      max_lwp = MAX(max_lwp*(1.-entr), MIN(1.0, lwp*(1.0-entr)))

      do k = k1, k2
         this_lwp = max_lwp*dz(k)/tdz
         if (cfr(k).gt.0.95.and.qc(k).lt.1.E-7.and.T(k).lt.253.16) then
            qc(k) = qc(k) + 0.05*this_lwp/Rho(k)/dz(k)
         elseif (cfr(k).gt.0.01.and.cfr(k).lt.0.99.and.T(k).lt.273.16.and.T(k).ge.253.16) then
            qc(k) = qc(k) + cfr(k)*this_lwp/Rho(k)/dz(k)
         elseif (cfr(k).gt.0.01.and.cfr(k).lt.0.99.and.T(k).le.298.16.and.T(k).ge.273.16) then
            qc(k) = qc(k) + cfr(k)*cfr(k)*this_lwp/Rho(k)/dz(k)
         endif
      enddo

      END SUBROUTINE adjust_cloudH2O

!+---+-----------------------------------------------------------------+
!..Do not alter any grid-explicitly resolved hydrometeors, rather only
!.. the supposed amounts due to the cloud fraction scheme.

      SUBROUTINE adjust_cloudFinal(cfr, qc, qi, Rho,dz, kts,kte,k_tropo)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: kts,kte,k_tropo
      REAL, DIMENSION(kts:kte), INTENT(IN):: cfr, Rho, dz
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: qc, qi
      REAL:: lwp, iwp, xfac
      INTEGER:: k

      lwp = 0.
      do k = kts, k_tropo
         if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
            lwp = lwp + qc(k)*Rho(k)*dz(k)
         endif
      enddo

      iwp = 0.
      do k = kts, k_tropo
         if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
            iwp = iwp + qi(k)*Rho(k)*dz(k)
         endif
      enddo

      if (lwp .gt. 1.0) then
         xfac = 1./lwp
         do k = kts, k_tropo
            if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
               qc(k) = qc(k)*xfac
            endif
         enddo
      endif

      if (iwp .gt. 1.0) then
         xfac = 1./iwp
         do k = kts, k_tropo
            if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
               qi(k) = qi(k)*xfac
            endif
         enddo
      endif

      END SUBROUTINE adjust_cloudFinal

!=================================================================================================================
  end module module_mp_thompson_cldfra3
!=================================================================================================================
