MODULE update_thermo_module
! Update thermodynamic variables & Store the diabatic effect to be used for coupling

USE shr_kind_mod, only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld,      only: ntracer,nk1,nk2,nk3,nhalo,nhalo_adv,netsz,nVGCM_seg,nk2_ext
USE constld,      only: a,b,dt,dxsq,dysq,dxdy,crad,rxtau_q, &
                        nomap,notopo,physics,radcode,turbulence, &
                        rad2deg   ! debug

! Subroutines being called
USE bound_channel_module, only: bound_channel
USE bound_extra,    only: bound_normal,bound_vert
USE halo_q,         only: halo_correc_q
USE damping,        only: damping_therm

IMPLICIT NONE
PRIVATE

PUBLIC :: update_thermodynamics,revise_thermodynamics

CONTAINS

! Local Subroutines:
!-------------------------------------------------------------------
! SUBROUTINE update_thermodynamics : update thermodynamic variables
! SUBROUTINE revise_thermodynamics : revise thermodynamic variables
!
! SUBROUTINE RELAXATION  : relaxation toward BG fields (Q3D algorithm)
!                          SUBROUTINE BOUND_G1 : used in RELAXATION
!                          SUBROUTINE BOUND_G2 : used in RELAXATION
!                          SUBROUTINE BOUND_G3 : used in RELAXATION
! SUBROUTINE LDIFFUSION2 : numerical diffusion
!                          FUNCTION TURB_H : used in LDIFFUSION2
!-------------------------------------------------------------------

!=======================================================================
      SUBROUTINE update_thermodynamics (N1, N2, channel)
!=======================================================================
!     Update the thermodynamic variables
      
      INTEGER, INTENT(IN) :: N1, &  ! AB forcing time index for previous timestep
                             N2     ! AB forcing time index for current timestep
      type(channel_t), INTENT(INOUT) :: channel   ! Channel data

      INTEGER :: I, J, K, nt, klow
      INTEGER :: num_seg,mi1,mim,mip,mj1,mjm,mjp
      
      LOGICAL :: DEBUG_PO = .TRUE.    ! JUNG_DEBUG 

      INTEGER :: CHN, NPO, IPO, JPO           ! JUNG_DEBUG 
      INTEGER :: KPO1, KPO2                   ! JUNG_DEBUG 
      REAL (KIND=r8) :: TEMP_TH               ! JUNG_DEBUG 
      
!**********************************************************************
      DO num_seg = 1, 4
!**********************************************************************
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim    = channel%seg(num_seg)%mim
      mip    = channel%seg(num_seg)%mip

      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment
      mjm    = channel%seg(num_seg)%mjm
      mjp    = channel%seg(num_seg)%mjp

!=====================================
      DO J = 1, mj1
      DO I = 1, mi1
      klow = channel%seg(num_seg)%KLOWQ_IJ(I,J)
!=====================================
!     Update theta, qv, and qt from advection tendency
      DO K = klow, nk2
        IF (DEBUG_PO) TEMP_TH = channel%seg(num_seg)%TH3D(I,J,K)
         
        channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,K)       &
                                         + A*channel%seg(num_seg)%FTH3D(I,J,K,N2) &
                                         + B*channel%seg(num_seg)%FTH3D(I,J,K,N1)

        channel%seg(num_seg)%QV3D(I,J,K) = channel%seg(num_seg)%QV3D(I,J,K)       &
                                         + A*channel%seg(num_seg)%FQV3D(I,J,K,N2) &
                                         + B*channel%seg(num_seg)%FQV3D(I,J,K,N1)

        DO nt = 1,ntracer
         channel%seg(num_seg)%QT3D(I,J,K,nt) = channel%seg(num_seg)%QT3D(I,J,K,nt)       &
                                             + A*channel%seg(num_seg)%FQT3D(I,J,K,N2,nt) &
                                             + B*channel%seg(num_seg)%FQT3D(I,J,K,N1,nt)
        ENDDO
       
       IF (DEBUG_PO) THEN
        TEMP_TH = channel%seg(num_seg)%TH3D(I,J,K) - TEMP_TH
        if (TEMP_TH .gt. 10.0_r8) then
          CALL CHK_title ('ADV: TH',channel%num_chn)
          CALL CHK_VALUE (channel%num_chn,I,J,K,num_seg,channel%seg(num_seg)%RLON_T(i,j), &
                          channel%seg(num_seg)%RLAT_T(i,j),    &
                          channel%seg(num_seg)%TH3D(i,j,k), &
                          var1=temp_th)
        endif        
       ENDIF                         
       
      ENDDO  ! k-loop     

!------------------------
      IF (PHYSICS) THEN
!------------------------
      
      DO K = klow, nk2
!      Update theta and qv from microphysics tendency

       IF (DEBUG_PO) TEMP_TH = channel%seg(num_seg)%TH3D(I,J,K)
       
       channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,K) &
                                        + DT*channel%seg(num_seg)%THAD_MICRO(I,J,K)
       channel%seg(num_seg)%QV3D(I,J,K) = channel%seg(num_seg)%QV3D(I,J,K) &
                                        + DT*channel%seg(num_seg)%QVAD_MICRO(I,J,K)

!      Update water species from dynamics and microphysics tendencies
       channel%seg(num_seg)%QC3D(I,J,K) = channel%seg(num_seg)%QC3D(I,J,K)          &
                                        + DT*channel%seg(num_seg)%QCAD_MICRO(I,J,K) &
                                        + A*channel%seg(num_seg)%FQC3D(I,J,K,N2)    &
                                        + B*channel%seg(num_seg)%FQC3D(I,J,K,N1)

       channel%seg(num_seg)%QI3D(I,J,K) = channel%seg(num_seg)%QI3D(I,J,K)          &
                                        + DT*channel%seg(num_seg)%QIAD_MICRO(I,J,K) &
                                        + A*channel%seg(num_seg)%FQI3D(I,J,K,N2)    &
                                        + B*channel%seg(num_seg)%FQI3D(I,J,K,N1)

       channel%seg(num_seg)%QR3D(I,J,K) = channel%seg(num_seg)%QR3D(I,J,K)          &
                                        + DT*channel%seg(num_seg)%QRAD_MICRO(I,J,K) &
                                        + A*channel%seg(num_seg)%FQR3D(I,J,K,N2)    &
                                        + B*channel%seg(num_seg)%FQR3D(I,J,K,N1)

       channel%seg(num_seg)%QS3D(I,J,K) = channel%seg(num_seg)%QS3D(I,J,K)          &
                                        + DT*channel%seg(num_seg)%QSAD_MICRO(I,J,K) &
                                        + A*channel%seg(num_seg)%FQS3D(I,J,K,N2)    &
                                        + B*channel%seg(num_seg)%FQS3D(I,J,K,N1)

       channel%seg(num_seg)%QG3D(I,J,K) = channel%seg(num_seg)%QG3D(I,J,K)          &
                                        + DT*channel%seg(num_seg)%QGAD_MICRO(I,J,K) &
                                        + A*channel%seg(num_seg)%FQG3D(I,J,K,N2)    &
                                        + B*channel%seg(num_seg)%FQG3D(I,J,K,N1)
      
       IF (DEBUG_PO) THEN
        TEMP_TH = channel%seg(num_seg)%TH3D(I,J,K) - TEMP_TH
        if (TEMP_TH .gt. 10.0_r8) then
          CALL CHK_title ('MICRO: TH',channel%num_chn)
          CALL CHK_VALUE (channel%num_chn,I,J,K,num_seg,channel%seg(num_seg)%RLON_T(i,j), &
                          channel%seg(num_seg)%RLAT_T(i,j),    &
                          channel%seg(num_seg)%TH3D(i,j,k), &
                          var1=temp_th)
        endif                                        
       ENDIF                         
      
      ENDDO  ! k-loop              

      IF (RADCODE) THEN
!     Update theta from radiation tendency
      DO K = klow, nk2
        channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,K)        &
                                         + DT*channel%seg(num_seg)%FTHRAD(I,J,K)
      ENDDO
      ENDIF
!------------------------
      ENDIF  ! PHYSICS
!------------------------

      IF (.NOT.NOTOPO) THEN
!--------------------------------
!     Below the lowest air layer
!--------------------------------
      DO K = 2, klow-1
       channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,KLOW)
       channel%seg(num_seg)%QV3D(I,J,K) = 0.0_r8
       DO nt = 1, ntracer
        channel%seg(num_seg)%QT3D(I,J,K,nt) = 0.0_r8
       ENDDO
      ENDDO

      IF (PHYSICS) THEN
        DO K = 2, klow-1
         channel%seg(num_seg)%QC3D(I,J,K) = 0.0_r8
         channel%seg(num_seg)%QI3D(I,J,K) = 0.0_r8
         channel%seg(num_seg)%QR3D(I,J,K) = 0.0_r8
         channel%seg(num_seg)%QS3D(I,J,K) = 0.0_r8
         channel%seg(num_seg)%QG3D(I,J,K) = 0.0_r8
        ENDDO
      ENDIF
!-------------------------
      ENDIF

!=====================================
      ENDDO  ! I-loop
      ENDDO  ! J-loop
!=====================================

!     Diabatic Tendency: radiation + microphysics
      IF (PHYSICS) THEN
       
       DO K = 2,nk2
        DO J = 1,mj1
         DO I = 1,mi1
          channel%seg(num_seg)%FTH3D_DIA(I,J,K) = &
                            channel%seg(num_seg)%FTH3D_DIA(I,J,K)  &
                          + channel%seg(num_seg)%THAD_MICRO(I,J,K) 

          channel%seg(num_seg)%FQV3D_DIA(I,J,K) = &
                            channel%seg(num_seg)%FQV3D_DIA(I,J,K)  &
                          + channel%seg(num_seg)%QVAD_MICRO(I,J,K)

          channel%seg(num_seg)%FQC3D_DIA(I,J,K) = &
                            channel%seg(num_seg)%FQC3D_DIA(I,J,K)  &
                          + channel%seg(num_seg)%QCAD_MICRO(I,J,K)

          channel%seg(num_seg)%FQI3D_DIA(I,J,K) = &
                            channel%seg(num_seg)%FQI3D_DIA(I,J,K)  &
                          + channel%seg(num_seg)%QIAD_MICRO(I,J,K)

          channel%seg(num_seg)%FQR3D_DIA(I,J,K) = &
                            channel%seg(num_seg)%FQR3D_DIA(I,J,K)  &
                          + channel%seg(num_seg)%QRAD_MICRO(I,J,K)

          channel%seg(num_seg)%FQS3D_DIA(I,J,K) = &
                            channel%seg(num_seg)%FQS3D_DIA(I,J,K)  &
                          + channel%seg(num_seg)%QSAD_MICRO(I,J,K)

          channel%seg(num_seg)%FQG3D_DIA(I,J,K) = &
                            channel%seg(num_seg)%FQG3D_DIA(I,J,K)  &
                          + channel%seg(num_seg)%QGAD_MICRO(I,J,K)
         ENDDO
        ENDDO
       ENDDO
       
       IF (RADCODE) THEN
       ! Add the radiation effects at the extended layers       
       DO K = 2,nk2_ext
        DO J = 1,mj1
         DO I = 1,mi1
          channel%seg(num_seg)%FTH3D_DIA(I,J,K) = &
                            channel%seg(num_seg)%FTH3D_DIA(I,J,K)  &
                          + channel%seg(num_seg)%FTHRAD(I,J,K)
         ENDDO
        ENDDO
       ENDDO
       ENDIF       
       
      ENDIF
!**********************************************************************
      ENDDO   ! num_seg
!**********************************************************************

!------------------------
      IF (PHYSICS) THEN
!------------------------
      CALL BOUND_NORMAL  (nhalo,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)

      CALL BOUND_CHANNEL (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
      IF (.not.nomap) &
      CALL HALO_CORREC_Q (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)

      CALL BOUND_VERT (channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                       QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
!------------------------
      ELSE
!------------------------
      CALL BOUND_NORMAL  (nhalo,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)
      CALL BOUND_CHANNEL (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)

      IF (.not.nomap) &
      CALL HALO_CORREC_Q (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)

      CALL BOUND_VERT (channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)
!------------------------
      ENDIF  ! PHYSICS
!------------------------

   END SUBROUTINE update_thermodynamics

!=======================================================================
   SUBROUTINE REVISE_THERMODYNAMICS (channel)
!=======================================================================
!  Revise the thermodynamic variables from turbulence, upper layer damping,
!  relaxation (coupling), and numerical diffusion (if necessary)

      type(channel_t), INTENT(INOUT) :: channel   ! channel data

!     Local variables
      LOGICAL :: LDIFFU = .TRUE.     ! ldiffusion is set only for th3d for now.
      INTEGER :: I, J, K, NT, KLOW
      INTEGER :: mi1, mj1, num_seg

!     TURBULENCE EFFECT (th, qv, qt, qc, qi)
!------------------------------
      IF (TURBULENCE) THEN
!------------------------------
!**********************************************************************
      DO num_seg = 1, 4
!**********************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment

      DO J = 1, mj1
      DO I = 1, mi1
      KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)

      DO K = klow, nk2
        channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,K) &
                                         + DT*channel%seg(num_seg)%THAD3(I,J,K)
        channel%seg(num_seg)%QV3D(I,J,K) = channel%seg(num_seg)%QV3D(I,J,K) &
                                         + DT*channel%seg(num_seg)%QVAD3(I,J,K)

        DO nt = 1,ntracer
         channel%seg(num_seg)%QT3D(I,J,K,nt) = channel%seg(num_seg)%QT3D(I,J,K,nt) &
                                             + DT*channel%seg(num_seg)%QTAD3(I,J,K,nt)
        ENDDO
      ENDDO

      IF (PHYSICS) THEN
        DO K = klow, nk2
          channel%seg(num_seg)%QC3D(I,J,K) = channel%seg(num_seg)%QC3D(I,J,K) &
                                           + DT*channel%seg(num_seg)%QCAD3(I,J,K)
          channel%seg(num_seg)%QI3D(I,J,K) = channel%seg(num_seg)%QI3D(I,J,K) &
                                           + DT*channel%seg(num_seg)%QIAD3(I,J,K)
        ENDDO
      ENDIF

!     DIABATIC TENDENCY:  TURBULENCE EFFECT (QR, QS, QG are not included)
      DO K = klow, nk2
        channel%seg(num_seg)%FTH3D_DIA(I,J,K) = channel%seg(num_seg)%FTH3D_DIA(I,J,K) &
                                              + channel%seg(num_seg)%THAD3(I,J,K)
        channel%seg(num_seg)%FQV3D_DIA(I,J,K) = channel%seg(num_seg)%FQV3D_DIA(I,J,K) &
                                              + channel%seg(num_seg)%QVAD3(I,J,K)

        DO nt = 1,ntracer
         channel%seg(num_seg)%FQT3D_DIA(I,J,K,nt) = &
            channel%seg(num_seg)%FQT3D_DIA(I,J,K,nt)+channel%seg(num_seg)%QTAD3(I,J,K,nt)
        ENDDO

      IF (PHYSICS) THEN
        channel%seg(num_seg)%FQC3D_DIA(I,J,K) = channel%seg(num_seg)%FQC3D_DIA(I,J,K) &
                                              + channel%seg(num_seg)%QCAD3(I,J,K)
        channel%seg(num_seg)%FQI3D_DIA(I,J,K) = channel%seg(num_seg)%FQI3D_DIA(I,J,K) &
                                              + channel%seg(num_seg)%QIAD3(I,J,K)
      ENDIF

      ENDDO

      IF (.NOT.NOTOPO) THEN
        ! Below the lowest air layer (This can change later)
        DO K = 2, KLOW-1
         channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,KLOW)
        ENDDO
      ENDIF      

      ENDDO  ! I-loop
      ENDDO  ! J-loop

!*****************************
      ENDDO  ! num_seg
!*****************************
!------------------------------
      ENDIF  ! TURBULENCE
!------------------------------

!     DAMPING EFFECT (th, qv, qc, and qi)
!     DIABATIC TENDENCY is calculated inside DAMPING_THERM
      IF (CRAD.NE.0.0_r8) CALL DAMPING_THERM (channel)
      
!---------------------------------------------------
!            RELAXATION EFFECT (th, qv, qt)
!---------------------------------------------------

      CALL RELAXATION (channel)

!     Synchronize boundary points after updating the values
!     with turbulence, damping & relaxation
!------------------------
      IF (PHYSICS) THEN
!------------------------
      CALL BOUND_NORMAL  (nhalo,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)

      CALL BOUND_CHANNEL (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
      IF (.not.nomap) &
      CALL HALO_CORREC_Q (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)

      CALL BOUND_VERT (channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                       QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
!------------------------
      ELSE
!------------------------
      CALL BOUND_NORMAL  (nhalo,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)

      CALL BOUND_CHANNEL (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)
      IF (.not.nomap) &
      CALL HALO_CORREC_Q (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)

      CALL BOUND_VERT (channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)
!------------------------
      ENDIF
!------------------------

!     Numerical diffusion only for th
!---------------------------------------------------
      IF (LDIFFU) THEN
!---------------------------------------------------
       CALL LDIFFUSION2 ( channel )     

!     Synchronize boundary points after updating the values with LDIFFUSION2
      CALL BOUND_NORMAL  (nhalo,channel,TH3D=.TRUE.)

      CALL BOUND_CHANNEL (nhalo_adv,channel,TH3D=.TRUE.)
      IF (.not.nomap) &
      CALL HALO_CORREC_Q (nhalo_adv,channel,TH3D=.TRUE.)

      CALL BOUND_VERT (channel,TH3D=.TRUE.)
!---------------------------------------------------
      ENDIF
!---------------------------------------------------

      END SUBROUTINE revise_thermodynamics

!=======================================================================      
      SUBROUTINE CHK_title (title,chn)
!=======================================================================
      CHARACTER (LEN=*), INTENT(IN)  :: title
      INTEGER, INTENT(IN) :: chn
              
      ! Local        
      INTEGER :: ifortw 
      
      ifortw = chn
      if (ifortw == 5) ifortw=55
      if (ifortw == 6) ifortw=65
      
      write(ifortw,*) TRIM(title),',I,J,K,num_seg,lon_deg,lat_deg,var'
        
      END SUBROUTINE chk_title         
      
!=======================================================================      
      SUBROUTINE CHK_VALUE (chn,ival,jval,kval,nval,lon,lat,var,var1,var2)
!=======================================================================
      INTEGER, INTENT(IN) :: chn,ival,jval,kval,nval
      REAL (KIND=r8), INTENT(IN) :: lon,lat,var
      REAL (KIND=r8), OPTIONAL, INTENT(IN) :: var1,var2
              
      ! Local        
      INTEGER :: ifortw 
      REAL (KIND=r8) :: lon_deg,lat_deg
      
      ifortw = chn
      if (ifortw == 5) ifortw=55
      if (ifortw == 6) ifortw=65
      
      lon_deg = lon*rad2deg
      lat_deg = lat*rad2deg

      IF (present(var1)) THEN
        if (present(var2)) then
        write(ifortw,'(4i5,2f8.1,2x,3E15.5)') ival,jval,kval,nval,lon_deg,lat_deg,var,var1,var2
        else
        write(ifortw,'(4i5,2f8.1,2x,2E15.5)') ival,jval,kval,nval,lon_deg,lat_deg,var,var1
        endif
      ELSE 
        write(ifortw,'(4i5,2f8.1,2x,E15.5)') ival,jval,kval,nval,lon_deg,lat_deg,var
      ENDIF
        
      END SUBROUTINE chk_value                
   
!=======================================================================
      SUBROUTINE RELAXATION (channel)
!=======================================================================
!     Relax the net-size (vGCM-cell-size) mean of deviation toward zero.
!
!     (prognostic) channel width is assumed to be 1 (CHN=1)
!
!     PLAN: In the averaging, TOPOGRAPHY effect should be included.
!------------------------------------------------------------------------
      type(channel_t), INTENT(INOUT) :: channel   ! channel data

!     Local variables
      LOGICAL :: INTERPOL = .TRUE.

      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4) :: MN_TH
      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4) :: MN_QV
      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4) :: MN_QC
      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4) :: MN_QI
      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4) :: MN_QR
      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4) :: MN_QS
      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4) :: MN_QG
      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4,ntracer) :: MN_QT

      REAL (KIND=r8) :: XMN(netsz,7+ntracer)

      REAL (KIND=r8), DIMENSION(:,:,:), ALLOCATABLE :: TEMP_TH
      REAL (KIND=r8), DIMENSION(:,:,:), ALLOCATABLE :: TEMP_QV
      REAL (KIND=r8), DIMENSION(:,:,:), ALLOCATABLE :: TEMP_QC
      REAL (KIND=r8), DIMENSION(:,:,:), ALLOCATABLE :: TEMP_QI
      REAL (KIND=r8), DIMENSION(:,:,:), ALLOCATABLE :: TEMP_QR
      REAL (KIND=r8), DIMENSION(:,:,:), ALLOCATABLE :: TEMP_QS
      REAL (KIND=r8), DIMENSION(:,:,:), ALLOCATABLE :: TEMP_QG
      REAL (KIND=r8), DIMENSION(:,:,:,:), ALLOCATABLE :: TEMP_QT

      INTEGER :: mi1,mj1,i,j,k,nt,klow
      INTEGER :: num_seg,ng,lsta,lend,chl,chn,npo

      REAL (KIND=r8) :: lcen,dist,dist1,dist2
      
      chn = 1    !  Assume (prognostic) channel width is 1.

! Calculate the vGCM-cell average of deviation.
!**********************************************************************
      DO num_seg = 1, 4
!**********************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment

!----------------------------
      IF (mi1.GT.mj1) THEN
!---------------------------- x-array
      DO ng = 1, nVGCM_seg
        DO K = 2, nk2

        DO npo = 1, netsz
         chl = channel%seg(num_seg)%lsta(ng) - 1 + npo
        
         XMN(npo,1) = channel%seg(num_seg)%TH3D(chl,chn,K) &
                    - channel%seg(num_seg)%TH3D_BG(chl,chn,K)
                    
         XMN(npo,2) = channel%seg(num_seg)%QV3D(chl,chn,K) &
                    - channel%seg(num_seg)%QV3D_BG(chl,chn,K)
         DO nt = 1,ntracer
          XMN(npo,2+nt) = channel%seg(num_seg)%QT3D(chl,chn,K,nt) &
                        - channel%seg(num_seg)%QT3D_BG(chl,chn,K,nt)
         ENDDO               
        ENDDO
        
        MN_TH(K-1,ng,num_seg) = SUM(XMN(:,1))/FLOAT(netsz)
        MN_QV(K-1,ng,num_seg) = SUM(XMN(:,2))/FLOAT(netsz)

        DO nt = 1,ntracer
         MN_QT(K-1,ng,num_seg,nt) = SUM(XMN(:,2+nt))/FLOAT(netsz)
        ENDDO

        IF (PHYSICS) THEN
          DO npo = 1, netsz
           chl = channel%seg(num_seg)%lsta(ng) - 1 + npo
           
           XMN(npo,3+ntracer) = channel%seg(num_seg)%QC3D(chl,chn,K) &
                              - channel%seg(num_seg)%QC3D_BG(chl,chn,K)
           XMN(npo,4+ntracer) = channel%seg(num_seg)%QI3D(chl,chn,K) &
                              - channel%seg(num_seg)%QI3D_BG(chl,chn,K)
           XMN(npo,5+ntracer) = channel%seg(num_seg)%QR3D(chl,chn,K) &
                              - channel%seg(num_seg)%QR3D_BG(chl,chn,K) 
           XMN(npo,6+ntracer) = channel%seg(num_seg)%QS3D(chl,chn,K) &
                              - channel%seg(num_seg)%QS3D_BG(chl,chn,K) 
           XMN(npo,7+ntracer) = channel%seg(num_seg)%QG3D(chl,chn,K) &
                              - channel%seg(num_seg)%QG3D_BG(chl,chn,K)         
          ENDDO
          
          MN_QC(K-1,ng,num_seg) = SUM(XMN(:,3+ntracer))/FLOAT(netsz)
          MN_QI(K-1,ng,num_seg) = SUM(XMN(:,4+ntracer))/FLOAT(netsz)
          MN_QR(K-1,ng,num_seg) = SUM(XMN(:,5+ntracer))/FLOAT(netsz)
          MN_QS(K-1,ng,num_seg) = SUM(XMN(:,6+ntracer))/FLOAT(netsz)
          MN_QG(K-1,ng,num_seg) = SUM(XMN(:,7+ntracer))/FLOAT(netsz)
        ENDIF

        ENDDO  ! k-loop
      ENDDO    ! ng-loop
!----------------------------
      ELSE
!----------------------------  y-array
      DO ng = 1, nVGCM_seg
        DO K = 2, nk2

        DO npo = 1, netsz
         chl = channel%seg(num_seg)%lsta(ng) - 1 + npo
        
         XMN(npo,1) = channel%seg(num_seg)%TH3D(chn,chl,K) &
                    - channel%seg(num_seg)%TH3D_BG(chn,chl,K)
                    
         XMN(npo,2) = channel%seg(num_seg)%QV3D(chn,chl,K) &
                    - channel%seg(num_seg)%QV3D_BG(chn,chl,K)
         DO nt = 1,ntracer
          XMN(npo,2+nt) = channel%seg(num_seg)%QT3D(chn,chl,K,nt) &
                        - channel%seg(num_seg)%QT3D_BG(chn,chl,K,nt)
         ENDDO               
        ENDDO
        
        MN_TH(K-1,ng,num_seg) = SUM(XMN(:,1))/FLOAT(netsz)
        MN_QV(K-1,ng,num_seg) = SUM(XMN(:,2))/FLOAT(netsz)

        DO nt = 1,ntracer
         MN_QT(K-1,ng,num_seg,nt) = SUM(XMN(:,2+nt))/FLOAT(netsz)
        ENDDO

        IF (PHYSICS) THEN
          DO npo = 1, netsz
           chl = channel%seg(num_seg)%lsta(ng) - 1 + npo
           
           XMN(npo,3+ntracer) = channel%seg(num_seg)%QC3D(chn,chl,K) &
                              - channel%seg(num_seg)%QC3D_BG(chn,chl,K)
           XMN(npo,4+ntracer) = channel%seg(num_seg)%QI3D(chn,chl,K) &
                              - channel%seg(num_seg)%QI3D_BG(chn,chl,K)
           XMN(npo,5+ntracer) = channel%seg(num_seg)%QR3D(chn,chl,K) &
                              - channel%seg(num_seg)%QR3D_BG(chn,chl,K) 
           XMN(npo,6+ntracer) = channel%seg(num_seg)%QS3D(chn,chl,K) &
                              - channel%seg(num_seg)%QS3D_BG(chn,chl,K) 
           XMN(npo,7+ntracer) = channel%seg(num_seg)%QG3D(chn,chl,K) &
                              - channel%seg(num_seg)%QG3D_BG(chn,chl,K)         
          ENDDO
          
          MN_QC(K-1,ng,num_seg) = SUM(XMN(:,3+ntracer))/FLOAT(netsz)
          MN_QI(K-1,ng,num_seg) = SUM(XMN(:,4+ntracer))/FLOAT(netsz)
          MN_QR(K-1,ng,num_seg) = SUM(XMN(:,5+ntracer))/FLOAT(netsz)
          MN_QS(K-1,ng,num_seg) = SUM(XMN(:,6+ntracer))/FLOAT(netsz)
          MN_QG(K-1,ng,num_seg) = SUM(XMN(:,7+ntracer))/FLOAT(netsz)
        ENDIF

        ENDDO  ! k-loop
      ENDDO    ! ng-loop
!----------------------------
      ENDIF
!----------------------------
!****************************
      ENDDO  ! num_seg
!****************************

      IF (INTERPOL) THEN
      ! Filling the halos along the channel: no mapping

      SELECT CASE (channel%num_chg)
      CASE(1)
        CALL BOUND_G1 (MN_TH)
        CALL BOUND_G1 (MN_QV)
        DO nt = 1, ntracer
         CALL BOUND_G1 (MN_QT(:,:,:,nt))
        ENDDO
        IF (PHYSICS) THEN
         CALL BOUND_G1 (MN_QC)
         CALL BOUND_G1 (MN_QI)
         CALL BOUND_G1 (MN_QR)
         CALL BOUND_G1 (MN_QS)
         CALL BOUND_G1 (MN_QG)
        ENDIF

      CASE(2)
        CALL BOUND_G2 (MN_TH)
        CALL BOUND_G2 (MN_QV)
        DO nt = 1, ntracer
         CALL BOUND_G2 (MN_QT(:,:,:,nt))
        ENDDO
        IF (PHYSICS) THEN
         CALL BOUND_G2 (MN_QC)
         CALL BOUND_G2 (MN_QI)
         CALL BOUND_G2 (MN_QR)
         CALL BOUND_G2 (MN_QS)
         CALL BOUND_G2 (MN_QG)
        ENDIF

      CASE(3)
        CALL BOUND_G3 (MN_TH)
        CALL BOUND_G3 (MN_QV)
        DO nt = 1, ntracer
         CALL BOUND_G3 (MN_QT(:,:,:,nt))
        ENDDO
        IF (PHYSICS) THEN
         CALL BOUND_G3 (MN_QC)
         CALL BOUND_G3 (MN_QI)
         CALL BOUND_G3 (MN_QR)
         CALL BOUND_G3 (MN_QS)
         CALL BOUND_G3 (MN_QG)
        ENDIF
      END SELECT

      ENDIF

!**********************************************************************
      DO num_seg = 1, 4
!**********************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment

      ALLOCATE(TEMP_TH(mi1,mj1,nk1))
      ALLOCATE(TEMP_QV(mi1,mj1,nk1))
      ALLOCATE(TEMP_QT(mi1,mj1,nk1,ntracer))

      IF (PHYSICS) THEN
       ALLOCATE(TEMP_QC(mi1,mj1,nk1))
       ALLOCATE(TEMP_QI(mi1,mj1,nk1))
       ALLOCATE(TEMP_QR(mi1,mj1,nk1))
       ALLOCATE(TEMP_QS(mi1,mj1,nk1))
       ALLOCATE(TEMP_QG(mi1,mj1,nk1))
      ENDIF

! Interpolation of vGCM values to CRM grid points
!----------------------------
      IF (mi1.GT.mj1) THEN
!---------------------------- x-array
      DO ng = 1, nVGCM_seg

        lsta = channel%seg(num_seg)%lsta(ng)
        lend = channel%seg(num_seg)%lend(ng)
        lcen = channel%seg(num_seg)%lcen(ng)

        DO chl = lsta, lend

        IF (INTERPOL) THEN

         DIST1 = lcen - float(chl)
         IF (DIST1.GE.0.0_r8) THEN
           DIST2 = FLOAT(netsz) - DIST1
           DIST  = DIST1
           NPO   = ng - 1
         ELSE
           DIST2 = FLOAT(netsz) + DIST1
           DIST  = ABS(DIST1)
           NPO   = ng + 1
         ENDIF

         DO K = 1, nk1
          TEMP_TH(chl,chn,K) = &
            (DIST*MN_TH(K,NPO,num_seg)+DIST2*MN_TH(K,ng,num_seg))/FLOAT(netsz)
          TEMP_QV(chl,chn,K) = &
            (DIST*MN_QV(K,NPO,num_seg)+DIST2*MN_QV(K,ng,num_seg))/FLOAT(netsz)
          DO nt = 1, ntracer
            TEMP_QT(chl,chn,K,nt) = &
              (DIST*MN_QT(K,NPO,num_seg,nt)+DIST2*MN_QT(K,ng,num_seg,nt))/FLOAT(netsz)
          ENDDO
         ENDDO

         IF (PHYSICS) THEN
           DO K = 1, nk1
            TEMP_QC(chl,chn,K) = &
              (DIST*MN_QC(K,NPO,num_seg)+DIST2*MN_QC(K,ng,num_seg))/FLOAT(netsz)
            TEMP_QI(chl,chn,K) = &
              (DIST*MN_QI(K,NPO,num_seg)+DIST2*MN_QI(K,ng,num_seg))/FLOAT(netsz)
            TEMP_QR(chl,chn,K) = &
              (DIST*MN_QR(K,NPO,num_seg)+DIST2*MN_QR(K,ng,num_seg))/FLOAT(netsz)
            TEMP_QS(chl,chn,K) = &
              (DIST*MN_QS(K,NPO,num_seg)+DIST2*MN_QS(K,ng,num_seg))/FLOAT(netsz)
            TEMP_QG(chl,chn,K) = &
              (DIST*MN_QG(K,NPO,num_seg)+DIST2*MN_QG(K,ng,num_seg))/FLOAT(netsz)
           ENDDO
         ENDIF

        ELSE  ! INTERPOL

         DO K = 1, nk1
          TEMP_TH(chl,chn,K) = MN_TH(K,ng,num_seg)
          TEMP_QV(chl,chn,K) = MN_QV(K,ng,num_seg)
          DO nt = 1, ntracer
            TEMP_QT(chl,chn,K,nt) = MN_QT(K,ng,num_seg,nt)
          ENDDO
         ENDDO

         IF (PHYSICS) THEN
           DO K = 1, nk1
            TEMP_QC(chl,chn,K) = MN_QC(K,ng,num_seg)
            TEMP_QI(chl,chn,K) = MN_QI(K,ng,num_seg)
            TEMP_QR(chl,chn,K) = MN_QR(K,ng,num_seg)
            TEMP_QS(chl,chn,K) = MN_QS(K,ng,num_seg)
            TEMP_QG(chl,chn,K) = MN_QG(K,ng,num_seg)
           ENDDO
         ENDIF

        ENDIF ! INTERPOL

        ENDDO  ! chl-loop

      ENDDO  ! ng-loop
!----------------------------
      ELSE
!---------------------------- y-array
      DO ng = 1, nVGCM_seg

        lsta = channel%seg(num_seg)%lsta(ng)
        lend = channel%seg(num_seg)%lend(ng)
        lcen = channel%seg(num_seg)%lcen(ng)

        DO chl = lsta, lend

        IF (INTERPOL) THEN

         DIST1 = lcen - float(chl)
         IF (DIST1.GE.0.0_r8) THEN
           DIST2 = FLOAT(netsz) - DIST1
           DIST  = DIST1
           NPO   = ng - 1
         ELSE
           DIST2 = FLOAT(netsz) + DIST1
           DIST  = ABS(DIST1)
           NPO   = ng + 1
         ENDIF

         DO K = 1, nk1
          TEMP_TH(chn,chl,K) = &
            (DIST*MN_TH(K,NPO,num_seg)+DIST2*MN_TH(K,ng,num_seg))/FLOAT(netsz)
          TEMP_QV(chn,chl,K) = &
            (DIST*MN_QV(K,NPO,num_seg)+DIST2*MN_QV(K,ng,num_seg))/FLOAT(netsz)
          DO nt = 1, ntracer
            TEMP_QT(chn,chl,K,nt) = &
              (DIST*MN_QT(K,NPO,num_seg,nt)+DIST2*MN_QT(K,ng,num_seg,nt))/FLOAT(netsz)
          ENDDO
         ENDDO

         IF (PHYSICS) THEN
           DO K = 1, nk1
            TEMP_QC(chn,chl,K) = &
              (DIST*MN_QC(K,NPO,num_seg)+DIST2*MN_QC(K,ng,num_seg))/FLOAT(netsz)
            TEMP_QI(chn,chl,K) = &
              (DIST*MN_QI(K,NPO,num_seg)+DIST2*MN_QI(K,ng,num_seg))/FLOAT(netsz)
            TEMP_QR(chn,chl,K) = &
              (DIST*MN_QR(K,NPO,num_seg)+DIST2*MN_QR(K,ng,num_seg))/FLOAT(netsz)
            TEMP_QS(chn,chl,K) = &
              (DIST*MN_QS(K,NPO,num_seg)+DIST2*MN_QS(K,ng,num_seg))/FLOAT(netsz)
            TEMP_QG(chn,chl,K) = &
              (DIST*MN_QG(K,NPO,num_seg)+DIST2*MN_QG(K,ng,num_seg))/FLOAT(netsz)
           ENDDO
         ENDIF

        ELSE  ! INTERPOL

         DO K = 1, nk1
          TEMP_TH(chn,chl,K) = MN_TH(K,ng,num_seg)
          TEMP_QV(chn,chl,K) = MN_QV(K,ng,num_seg)
          DO nt = 1, ntracer
            TEMP_QT(chn,chl,K,nt) = MN_QT(K,ng,num_seg,nt)
          ENDDO
         ENDDO

         IF (PHYSICS) THEN
           DO K = 1, nk1
            TEMP_QC(chn,chl,K) = MN_QC(K,ng,num_seg)
            TEMP_QI(chn,chl,K) = MN_QI(K,ng,num_seg)
            TEMP_QR(chn,chl,K) = MN_QR(K,ng,num_seg)
            TEMP_QS(chn,chl,K) = MN_QS(K,ng,num_seg)
            TEMP_QG(chn,chl,K) = MN_QG(K,ng,num_seg)
           ENDDO
         ENDIF

        ENDIF ! INTERPOL

        ENDDO  ! chl-loop

      ENDDO  ! ng-loop
!----------------------------
      ENDIF
!----------------------------

! Apply the relaxation (topography is not considered)

      !--------------------------------
      IF (RXTAU_q .EQ. 0.0_r8) THEN 
      !--------------------------------
      DO K = 2, nk2
       DO J = 1, mj1
        DO I = 1, mi1

          IF (channel%seg(num_seg)%TAU_RX_TQ(I,J,K).NE.0.0_r8) THEN

           channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,K) &
                    - DT*TEMP_TH(I,J,K-1)/channel%seg(num_seg)%TAU_RX_TQ(I,J,K)
           channel%seg(num_seg)%QV3D(I,J,K) = channel%seg(num_seg)%QV3D(I,J,K) &
                    - DT*TEMP_QV(I,J,K-1)/channel%seg(num_seg)%TAU_RX_TQ(I,J,K)

           DO nt = 1,ntracer
            channel%seg(num_seg)%QT3D(I,J,K,nt) = channel%seg(num_seg)%QT3D(I,J,K,nt) &
                    - DT*TEMP_QT(I,J,K-1,nt)/channel%seg(num_seg)%TAU_RX_TQ(I,J,K)
           ENDDO

           IF (PHYSICS) THEN
             channel%seg(num_seg)%QC3D(I,J,K) = channel%seg(num_seg)%QC3D(I,J,K) &
                    - DT*TEMP_QC(I,J,K-1)/channel%seg(num_seg)%TAU_RX_TQ(I,J,K)
             channel%seg(num_seg)%QI3D(I,J,K) = channel%seg(num_seg)%QI3D(I,J,K) &
                    - DT*TEMP_QI(I,J,K-1)/channel%seg(num_seg)%TAU_RX_TQ(I,J,K)
             channel%seg(num_seg)%QR3D(I,J,K) = channel%seg(num_seg)%QR3D(I,J,K) &
                    - DT*TEMP_QR(I,J,K-1)/channel%seg(num_seg)%TAU_RX_TQ(I,J,K)
             channel%seg(num_seg)%QS3D(I,J,K) = channel%seg(num_seg)%QS3D(I,J,K) &
                    - DT*TEMP_QS(I,J,K-1)/channel%seg(num_seg)%TAU_RX_TQ(I,J,K)
             channel%seg(num_seg)%QG3D(I,J,K) = channel%seg(num_seg)%QG3D(I,J,K) &
                    - DT*TEMP_QG(I,J,K-1)/channel%seg(num_seg)%TAU_RX_TQ(I,J,K)
           ENDIF

          ENDIF  ! (TAU_RX_TQ(I,J,K).NE.0.0_r8)

        ENDDO
       ENDDO
      ENDDO
      !------------
      ELSE
      !------------
      DO K = 2, nk2
       DO J = 1, mj1
        DO I = 1, mi1

           channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,K) &
                    - DT*TEMP_TH(I,J,K-1)/RXTAU_q
           channel%seg(num_seg)%QV3D(I,J,K) = channel%seg(num_seg)%QV3D(I,J,K) &
                    - DT*TEMP_QV(I,J,K-1)/RXTAU_q

           DO nt = 1,ntracer
            channel%seg(num_seg)%QT3D(I,J,K,nt) = channel%seg(num_seg)%QT3D(I,J,K,nt) &
                    - DT*TEMP_QT(I,J,K-1,nt)/RXTAU_q
           ENDDO

           IF (PHYSICS) THEN
             channel%seg(num_seg)%QC3D(I,J,K) = channel%seg(num_seg)%QC3D(I,J,K) &
                    - DT*TEMP_QC(I,J,K-1)/RXTAU_q
             channel%seg(num_seg)%QI3D(I,J,K) = channel%seg(num_seg)%QI3D(I,J,K) &
                    - DT*TEMP_QI(I,J,K-1)/RXTAU_q
             channel%seg(num_seg)%QR3D(I,J,K) = channel%seg(num_seg)%QR3D(I,J,K) &
                    - DT*TEMP_QR(I,J,K-1)/RXTAU_q
             channel%seg(num_seg)%QS3D(I,J,K) = channel%seg(num_seg)%QS3D(I,J,K) &
                    - DT*TEMP_QS(I,J,K-1)/RXTAU_q
             channel%seg(num_seg)%QG3D(I,J,K) = channel%seg(num_seg)%QG3D(I,J,K) &
                    - DT*TEMP_QG(I,J,K-1)/RXTAU_q
           ENDIF

        ENDDO
       ENDDO
      ENDDO      
      !------------
      ENDIF
      !------------

      IF (.NOT.NOTOPO) THEN
!     Below the lowest air layer
      DO J = 1, mj1
       DO I = 1, mi1
        KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)
        DO K = 2, KLOW-1
         channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,KLOW)
        ENDDO
       ENDDO
      ENDDO
      ENDIF

      DEALLOCATE(TEMP_TH)
      DEALLOCATE(TEMP_QV)
      DEALLOCATE(TEMP_QT)

      IF (PHYSICS) THEN
       DEALLOCATE(TEMP_QC)
       DEALLOCATE(TEMP_QI)
       DEALLOCATE(TEMP_QR)
       DEALLOCATE(TEMP_QS)
       DEALLOCATE(TEMP_QG)
      ENDIF

!**********************************************************************
      ENDDO  ! num_seg
!**********************************************************************

      CONTAINS

      SUBROUTINE BOUND_G1 (A)
!     Filling the halos along the channel (Group1)
      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4),INTENT(INOUT) :: A

      INTEGER :: chl_p,chl_m,K

        chl_p = nVGCM_seg + 1
        chl_m = 0

        DO K = 1, nk1
          A(K,chl_p,1) = A(K,1,2)
          A(K,chl_m,1) = A(K,nVGCM_seg,4)

          A(K,chl_p,2) = A(K,1,3)
          A(K,chl_m,2) = A(K,nVGCM_seg,1)

          A(K,chl_p,3) = A(K,1,4)
          A(K,chl_m,3) = A(K,nVGCM_seg,2)

          A(K,chl_p,4) = A(K,1,1)
          A(K,chl_m,4) = A(K,nVGCM_seg,3)
        ENDDO

      END SUBROUTINE bound_g1

      SUBROUTINE BOUND_G2 (A)
!     Filling the halos along the channel (Group2)
      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4),INTENT(INOUT) :: A

      INTEGER :: chl_p,chl_m,K

        chl_p = nVGCM_seg + 1
        chl_m = 0

        DO K = 1, nk1
          A(K,chl_p,1) = A(K,1,2)
          A(K,chl_m,1) = A(K,1,4)

          A(K,chl_p,2) = A(K,1,3)
          A(K,chl_m,2) = A(K,nVGCM_seg,1)

          A(K,chl_p,3) = A(K,nVGCM_seg,4)
          A(K,chl_m,3) = A(K,nVGCM_seg,2)

          A(K,chl_p,4) = A(K,nVGCM_seg,3)
          A(K,chl_m,4) = A(K,1,1)
        ENDDO

      END SUBROUTINE bound_g2

      SUBROUTINE BOUND_G3 (A)
!     Filling the halos along the channel (Group3)
      REAL (KIND=r8),DIMENSION(nk1,0:nVGCM_seg+1,4),INTENT(INOUT) :: A

      INTEGER :: chl_p,chl_m,K

        chl_p = nVGCM_seg + 1
        chl_m = 0

        DO K = 1, nk1
          A(K,chl_p,1) = A(K,1,2)
          A(K,chl_m,1) = A(K,1,4)

          A(K,chl_p,2) = A(K,nVGCM_seg,3)
          A(K,chl_m,2) = A(K,nVGCM_seg,1)

          A(K,chl_p,3) = A(K,nVGCM_seg,2)
          A(K,chl_m,3) = A(K,nVGCM_seg,4)

          A(K,chl_p,4) = A(K,1,3)
          A(K,chl_m,4) = A(K,1,1)
        ENDDO

      END SUBROUTINE bound_g3

      END SUBROUTINE relaxation

!=======================================================================
      SUBROUTINE LDIFFUSION2 ( channel )
!=======================================================================
!     Linear 2nd-order Horizontal Diffusion

      type(channel_t), INTENT(INOUT) :: channel   ! Channel data

      ! Local variables

      ! # of variables, to which diffusion is applied. Currently, only th3d is used.
      INTEGER, PARAMETER   :: nvar = 1
      REAL (KIND=r8) :: COEF_FAC = 10.0_r8

      REAL (KIND=r8), DIMENSION(:,:,:,:), ALLOCATABLE :: AVAL
      REAL (KIND=r8), DIMENSION(:,:,:,:), ALLOCATABLE :: TEND

      REAL (KIND=r8) :: COEFX,COEFY,COEFA

      REAL (KIND=r8) :: KVAL,KVAL0
      REAL (KIND=r8) :: FXP,FXM,FYP,FYM
      REAL (KIND=r8) :: FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0
      REAL (KIND=r8) :: FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0

      INTEGER :: I,J,K,KLOW,KLOWQ_MAX
      INTEGER :: mi1,mj1,mim,mip,mjm,mjp,num_seg,nt

      KVAL0 = DXSQ/(4.0_r8*DT)
      KVAL  = KVAL0/COEF_FAC

!**********************************************************************
      DO num_seg = 1, 4
!**********************************************************************
      KLOWQ_MAX = channel%seg(num_seg)%KLOWQ_MAX_GLOB

      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim = channel%seg(num_seg)%mim
      mip = channel%seg(num_seg)%mip

      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment
      mjm = channel%seg(num_seg)%mjm
      mjp = channel%seg(num_seg)%mjp

!-------------------------
!     Allocate Local Data
!-------------------------
      ALLOCATE(AVAL(mim:mip,mjm:mjp,nk3,nvar))
      ALLOCATE(TEND(mi1,mj1,nk2,nvar))

!-------------------------
!     Arrange Input Data
!-------------------------
      AVAL(:,:,:,1) = channel%seg(num_seg)%TH3D(:,:,:)
      ! Add another variables here if necessary (if nvar > 1).

!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      DO J = 1, mj1
      DO I = 1, mi1

      COEFX = channel%seg(num_seg)%COEFX(I,J)
      COEFY = channel%seg(num_seg)%COEFY(I,J)
      COEFA = channel%seg(num_seg)%COEFA(I,J)

      KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)

      DO K=KLOW,KLOWQ_MAX-1
!     Influenced by topography

        IF(channel%seg(num_seg)%DM_Q_TE(I,J,K).EQ.1) THEN
         FXP    = 0.0_r8
         FYP_XP = 0.0_r8
         FYM_XP = 0.0_r8
        ELSE
         FXP    = channel%seg(num_seg)%RGG_U(1,I,J)
         FYP_XP = channel%seg(num_seg)%RGG_V(2,I+1,J)
         FYM_XP = channel%seg(num_seg)%RGG_V(2,I+1,J-1)
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TW(I,J,K).EQ.1) THEN
         FXM    = 0.0_r8
         FYP_XM = 0.0_r8
         FYM_XM = 0.0_r8
        ELSE
         FXM    = channel%seg(num_seg)%RGG_U(1,I-1,J)
         FYP_XM = channel%seg(num_seg)%RGG_V(2,I-1,J)
         FYM_XM = channel%seg(num_seg)%RGG_V(2,I-1,J-1)
        ENDIF

         FYP_X0 = 0.0_r8
         FYM_X0 = 0.0_r8

        IF(channel%seg(num_seg)%DM_Q_TN(I,J,K).EQ.1) THEN
         FYP    = 0.0_r8
         FXP_YP = 0.0_r8
         FXM_YP = 0.0_r8
        ELSE
         FYP    = channel%seg(num_seg)%RGG_V(3,I,J)
         FXP_YP = channel%seg(num_seg)%RGG_U(2,I,J+1)
         FXM_YP = channel%seg(num_seg)%RGG_U(2,I-1,J+1)
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TS(I,J,K).EQ.1) THEN
         FYM    = 0.0_r8
         FXP_YM = 0.0_r8
         FXM_YM = 0.0_r8
        ELSE
         FYM    = channel%seg(num_seg)%RGG_V(3,I,J-1)
         FXP_YM = channel%seg(num_seg)%RGG_U(2,I,J-1)
         FXM_YM = channel%seg(num_seg)%RGG_U(2,I-1,J-1)
        ENDIF

         FXP_Y0 = 0.0_r8
         FXM_Y0 = 0.0_r8

      do nt = 1, nvar
        TEND(I,J,K,nt) = &
          TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,                       &
                 FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,               &
                 FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,               &
                 AVAL(I-1,J-1,K,nt),AVAL(I,J-1,K,nt),AVAL(I+1,J-1,K,nt),  &
                 AVAL(I-1,J,K,nt),AVAL(I,J,K,nt),AVAL(I+1,J,K,nt),        &
                 AVAL(I-1,J+1,K,nt),AVAL(I,J+1,K,nt),AVAL(I+1,J+1,K,nt))
      enddo

      ENDDO   ! K-LOOP (Low level: mountainous region)

      DO K=KLOWQ_MAX,NK2
!     Not influenced by topography (upper atmosphere)

         FXP    = channel%seg(num_seg)%RGG_U(1,I,J)
         FYP_XP = channel%seg(num_seg)%RGG_V(2,I+1,J)
         FYM_XP = channel%seg(num_seg)%RGG_V(2,I+1,J-1)

         FXM    = channel%seg(num_seg)%RGG_U(1,I-1,J)
         FYP_XM = channel%seg(num_seg)%RGG_V(2,I-1,J)
         FYM_XM = channel%seg(num_seg)%RGG_V(2,I-1,J-1)

         FYP_X0 = 0.0_r8
         FYM_X0 = 0.0_r8

         FYP    = channel%seg(num_seg)%RGG_V(3,I,J)
         FXP_YP = channel%seg(num_seg)%RGG_U(2,I,J+1)
         FXM_YP = channel%seg(num_seg)%RGG_U(2,I-1,J+1)

         FYM    = channel%seg(num_seg)%RGG_V(3,I,J-1)
         FXP_YM = channel%seg(num_seg)%RGG_U(2,I,J-1)
         FXM_YM = channel%seg(num_seg)%RGG_U(2,I-1,J-1)

         FXP_Y0 = 0.0_r8
         FXM_Y0 = 0.0_r8

      do nt = 1, nvar
        TEND(I,J,K,nt) = &
           TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,                       &
                  FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,               &
                  FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,               &
                  AVAL(I-1,J-1,K,nt),AVAL(I,J-1,K,nt),AVAL(I+1,J-1,K,nt),  &
                  AVAL(I-1,J,K,nt),AVAL(I,J,K,nt),AVAL(I+1,J,K,nt),        &
                  AVAL(I-1,J+1,K,nt),AVAL(I,J+1,K,nt),AVAL(I+1,J+1,K,nt))
      enddo

      ENDDO  ! K-LOOP (high level: mountain-free region)

      ENDDO  ! J-loop
      ENDDO  ! I-loop

      DO nt =1, nvar
       DO K = 2, NK2
        DO J = 1, mj1
         DO I = 1, mi1
          AVAL(I,J,K,nt) = AVAL(I,J,K,nt) + DT*KVAL*TEND(I,J,K,nt)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

!-------------------------
!     Arrange Output Data
!-------------------------
      channel%seg(num_seg)%TH3D(:,:,:) = AVAL(:,:,:,1)
      ! Add another variables here if necessary (if nvar > 1).

!-------------------------
!     Deallocate Local Data
!-------------------------
      DEALLOCATE(AVAL)
      DEALLOCATE(TEND)

!*************************
      ENDDO  ! num_seg
!*************************

   END SUBROUTINE ldiffusion2

!=======================================================================
      FUNCTION TURB_H(CX,CY,CA,KXP,KXM,KYP,KYM, &
                      KYP_XP,KYM_XP,KYP_XM,KYM_XM,KYP_X0,KYM_X0, &
                      KXP_YP,KXM_YP,KXP_YM,KXM_YM,KXP_Y0,KXM_Y0, &
                      A1,A2,A3,A4,A5,A6,A7,A8,A9) &
               RESULT(TURB_H_result)
!=======================================================================
!     Calculate the turbulent effects (horizontal)
!     CX, CY, CA: coefficients
!     KXP and others: C (see the document)
!     A: target variable
!     A7(I-1,J+1,K),A8(I,J+1,K),A9(I+1,J+1,K)
!     A4(I-1,J  ,K),A5(I,J  ,K),A6(I+1,J  ,K)
!     A1(I-1,J-1,K),A2(I,J-1,K),A3(I+1,J-1,K)
!------------------------------------------------------------------
      IMPLICIT NONE

      REAL (KIND=r8), INTENT(IN) :: CX,CY,CA,KXP,KXM,KYP,KYM
      REAL (KIND=r8), INTENT(IN) :: KYP_XP,KYM_XP,KYP_XM,KYM_XM,KYP_X0,KYM_X0
      REAL (KIND=r8), INTENT(IN) :: KXP_YP,KXM_YP,KXP_YM,KXM_YM,KXP_Y0,KXM_Y0
      REAL (KIND=r8), INTENT(IN) :: A1,A2,A3,A4,A5,A6,A7,A8,A9
      REAL (KIND=r8) :: TURB_H_result

      TURB_H_result = CX*(KXP*(A6-A5)-KXM*(A5-A4))      &
                    + CY*(KYP*(A8-A5)-KYM*(A5-A2))      &
                    + CA*(KYP_XP*(A9-A6)+KYM_XP*(A6-A3) &
                         +KYP_X0*(A8-A5)+KYM_X0*(A5-A2) &
                         -KYP_XM*(A7-A4)-KYM_XM*(A4-A1) &
                         +KXP_YP*(A9-A8)+KXM_YP*(A8-A7) &
                         +KXP_Y0*(A6-A5)+KXM_Y0*(A5-A4) &
                         -KXP_YM*(A3-A2)-KXM_YM*(A2-A1))

      END FUNCTION TURB_H

END MODULE update_thermo_module
