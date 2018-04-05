!=============================================================================
      SUBROUTINE physics_interface (N2, ITT, channel)
!=============================================================================                                    
!     This is the interface subroutine between the Vector Vorticity Model (VVM)
!     and the Lorenz grid Microphysics (also, call RADIATION_RRTMG).
!  
!     channel_width = 1 is assumed: chn = 1
!     If channel_width is wider than 1, then additional do-loop should be added.
!
!     For column physics (microphysics & radiation), 
!     the whole channel array is internally formed and used for calculation. 
!
!     OUTPUT: THAD_MICRO,QVAD_MICRO,
!             QCAD_MICRO,QIAD_MICRO, 
!             QRAD_MICRO,QSAD_MICRO,QGAD_MICRO (tendencies due to microphysics)
!             SPREC (sfc precipitation)
!
!    (IN/OUT) FQR3D,FQS3D,FQG3D (advection tendencies of qr, qs, qg)
!             (FTHRAD can be prescribed here if RADCODE=.F.)
!-----------------------------------------------------------------------------
      USE shr_kind_mod,   only: r8 => shr_kind_r8
      USE vvm_data_types, only: channel_t   
      
      USE parmsld, only: channel_seg_l,nk1,nk2
      USE constld, only: radcode,microcode
      
      USE mphysics_variables, only: theta,qv,qc,qi,qr,qs,qg,     &
                                    tendency_microphysics_theta, &
                                    tendency_microphysics_qv,    &
                                    tendency_microphysics_qc,    &
                                    tendency_microphysics_qi,    &
                                    tendency_microphysics_qr,    &
                                    tendency_microphysics_qs,    &
                                    tendency_microphysics_qg,    &
                                    tendency_rain,tendency_snow,tendency_graupel, &
                                    surface_rain, &
                                    latent_heating_rate,hx_sub
      
      IMPLICIT NONE

      ! Input arguments 
      INTEGER, INTENT(IN) :: N2,   &  ! AB forcing time index for current timestep    
                             ITT      ! Current model time step
                             
      type(channel_t), INTENT(INOUT) :: channel   ! Channel data
      
     ! Local variables 
      INTEGER :: L,mi1,mj1,num_seg,chp,chn,chl,K
!-----------------------------------------------------------------------------      
      L   = N2
      CHN = 1      ! Assume: channel_width = 1 

! Radiation
!*****************************************************************************      
      IF (RADCODE) THEN
!***************************************************************************** 
      
        CALL RADIATION_RRTMG(ITT, channel)    
        
!**************************************         
      ENDIF 
!**************************************      

! Microphysics
!*****************************************************************************
      IF (MICROCODE) THEN    
!*****************************************************************************      

       ! Initialize tendency terms
       TENDENCY_MICROPHYSICS_THETA = 0.0_r8
       TENDENCY_MICROPHYSICS_QV    = 0.0_r8
       TENDENCY_MICROPHYSICS_QC    = 0.0_r8
       TENDENCY_MICROPHYSICS_QI    = 0.0_r8
       TENDENCY_MICROPHYSICS_QR    = 0.0_r8
       TENDENCY_MICROPHYSICS_QS    = 0.0_r8
       TENDENCY_MICROPHYSICS_QG    = 0.0_r8
       
       TENDENCY_RAIN    = 0.0_r8
       TENDENCY_SNOW    = 0.0_r8
       TENDENCY_GRAUPEL = 0.0_r8
       
       LATENT_HEATING_RATE = 0.0_r8 
       
       SURFACE_RAIN = 0.0_r8
       
!      Assign input arrays from model fields 
!      (combine 4 channel segments to set up one channel array)        
!===============================   
       DO num_seg = 1, 4
!=============================== 
       mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
       mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment    
      
!      ---------------------      
       IF (mi1.GT.mj1) THEN
!      ---------------------  ! x-array
        DO chp = 1,mi1
         chl = channel_seg_l*(num_seg-1) + chp
       
         HX_SUB(chl) = channel%seg(num_seg)%HX(chp,chn)
         DO K = 1, nk1 
          THETA(chl,K) = channel%seg(num_seg)%TH3D(chp,chn,K+1)
          QV(chl,K)    = channel%seg(num_seg)%QV3D(chp,chn,K+1)
          QC(chl,K)    = channel%seg(num_seg)%QC3D(chp,chn,K+1)
          QI(chl,K)    = channel%seg(num_seg)%QI3D(chp,chn,K+1)
          QR(chl,K)    = channel%seg(num_seg)%QR3D(chp,chn,K+1)
          QS(chl,K)    = channel%seg(num_seg)%QS3D(chp,chn,K+1)
          QG(chl,K)    = channel%seg(num_seg)%QG3D(chp,chn,K+1)
         ENDDO
        ENDDO   
!      ---------------------                     
       ELSE
!      ---------------------  ! y-array   
        DO chp = 1,mj1
         chl = channel_seg_l*(num_seg-1) + chp
        
         HX_SUB(chl) = channel%seg(num_seg)%HX(chn,chp)
         DO K = 1, nk1
          THETA(chl,K) = channel%seg(num_seg)%TH3D(chn,chp,K+1)
          QV(chl,K)    = channel%seg(num_seg)%QV3D(chn,chp,K+1)
          QC(chl,K)    = channel%seg(num_seg)%QC3D(chn,chp,K+1)
          QI(chl,K)    = channel%seg(num_seg)%QI3D(chn,chp,K+1)
          QR(chl,K)    = channel%seg(num_seg)%QR3D(chn,chp,K+1)
          QS(chl,K)    = channel%seg(num_seg)%QS3D(chn,chp,K+1)
          QG(chl,K)    = channel%seg(num_seg)%QG3D(chn,chp,K+1)
         ENDDO
        ENDDO 
       
!      ---------------------       
       ENDIF  
!      ---------------------              

!===============================   
       ENDDO  ! num_seg 
!=============================== 

       CALL MICROPHYSICS
       
!      Assign output arrays to model fields  
!      (decompose one channel array to 4 channel segments)       
!===============================   
       DO num_seg = 1, 4
!=============================== 
       mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
       mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment    
      
!      ---------------------      
       IF (mi1.GT.mj1) THEN
!      ---------------------  ! x-array

       DO chp = 1,mi1
        chl = channel_seg_l*(num_seg-1) + chp
        
        ! Surface precipitation
        channel%seg(num_seg)%SPREC(chp,chn) = SURFACE_RAIN(chl)

        ! Microphysics tendency terms & latent heating rate
        DO K = 2, nk2
         channel%seg(num_seg)%THAD_MICRO(chp,chn,K) = &
                TENDENCY_MICROPHYSICS_THETA(chl,k-1)
         channel%seg(num_seg)%QVAD_MICRO(chp,chn,K) = &
                   TENDENCY_MICROPHYSICS_QV(chl,k-1)
         channel%seg(num_seg)%QCAD_MICRO(chp,chn,K) = &
                   TENDENCY_MICROPHYSICS_QC(chl,k-1)
         channel%seg(num_seg)%QIAD_MICRO(chp,chn,K) = &
                   TENDENCY_MICROPHYSICS_QI(chl,k-1)
         channel%seg(num_seg)%QSAD_MICRO(chp,chn,K) = &
                   TENDENCY_MICROPHYSICS_QS(chl,k-1)
         channel%seg(num_seg)%QGAD_MICRO(chp,chn,K) = &
                   TENDENCY_MICROPHYSICS_QG(chl,k-1)
         channel%seg(num_seg)%QRAD_MICRO(chp,chn,K) = &
                   TENDENCY_MICROPHYSICS_QR(chl,k-1)
        ENDDO 

        ! Rain, snow and graupel tendencies returned from Microphysics 
        ! contain the vertical advection by the terminal velocity
        DO K = 2, nk2
         channel%seg(num_seg)%FQR3D(chp,chn,K,L) = &
            channel%seg(num_seg)%FQR3D(chp,chn,K,L) + TENDENCY_RAIN(chl,k-1)
         channel%seg(num_seg)%FQS3D(chp,chn,K,L) = &
            channel%seg(num_seg)%FQS3D(chp,chn,K,L) + TENDENCY_SNOW(chl,k-1)
         channel%seg(num_seg)%FQG3D(chp,chn,K,L) = &
            channel%seg(num_seg)%FQG3D(chp,chn,K,L) + TENDENCY_GRAUPEL(chl,k-1)
        ENDDO 

       ENDDO  ! chp_loop
              
!      ---------------------                     
       ELSE
!      ---------------------  ! y-array   

       DO chp = 1,mj1
        chl = channel_seg_l*(num_seg-1) + chp
        
        ! Surface precipitation
        channel%seg(num_seg)%SPREC(chn,chp) = SURFACE_RAIN(chl)

        ! Microphysics tendency terms & latent heating rate
        DO K = 2, nk2
         channel%seg(num_seg)%THAD_MICRO(chn,chp,K) = &
                TENDENCY_MICROPHYSICS_THETA(chl,k-1)
         channel%seg(num_seg)%QVAD_MICRO(chn,chp,K) = &
                   TENDENCY_MICROPHYSICS_QV(chl,k-1)
         channel%seg(num_seg)%QCAD_MICRO(chn,chp,K) = &
                   TENDENCY_MICROPHYSICS_QC(chl,k-1)
         channel%seg(num_seg)%QIAD_MICRO(chn,chp,K) = &
                   TENDENCY_MICROPHYSICS_QI(chl,k-1)
         channel%seg(num_seg)%QSAD_MICRO(chn,chp,K) = &
                   TENDENCY_MICROPHYSICS_QS(chl,k-1)
         channel%seg(num_seg)%QGAD_MICRO(chn,chp,K) = &
                   TENDENCY_MICROPHYSICS_QG(chl,k-1)
         channel%seg(num_seg)%QRAD_MICRO(chn,chp,K) = &
                   TENDENCY_MICROPHYSICS_QR(chl,k-1)
        ENDDO 

        ! Rain, snow and graupel tendencies returned from Microphysics 
        ! contain the vertical advection by the terminal velocity
        DO K = 2, nk2
         channel%seg(num_seg)%FQR3D(chn,chp,K,L) = &
            channel%seg(num_seg)%FQR3D(chn,chp,K,L) + TENDENCY_RAIN(chl,k-1)
         channel%seg(num_seg)%FQS3D(chn,chp,K,L) = &
            channel%seg(num_seg)%FQS3D(chn,chp,K,L) + TENDENCY_SNOW(chl,k-1)
         channel%seg(num_seg)%FQG3D(chn,chp,K,L) = &
            channel%seg(num_seg)%FQG3D(chn,chp,K,L) + TENDENCY_GRAUPEL(chl,k-1)
        ENDDO 

       ENDDO  ! chp_loop
              
!      ---------------------       
       ENDIF  
!      ---------------------              

!===============================   
       ENDDO  ! num_seg 
!===============================        

!**************************************
      ENDIF  ! MICROCODE
!**************************************    

      RETURN
      END SUBROUTINE physics_interface
