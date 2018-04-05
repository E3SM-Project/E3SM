MODULE q3d_mapping_module
! Contains the programs related to preparing the mapping information

   USE shr_kind_mod,    only: r8 => shr_kind_r8
   USE vGCM_data_types, only: vGCM_state_t
   USE vvm_data_types,  only: channel_t
      
   USE parmsld, only: channel_seg_l,netsz,nhalo
   USE constld, only: dx,dy,dxsq,dysq,dxdy,dx_gcm,dy_gcm,wrxmu, &
                      pi,rearth,omega,rad2deg,nomap

   ! Subroutines being called
   USE utils, only: indexr,find_face,latval,lonval,alphav,betav,jaco, &
                    gcont11,gcont12,gcont22,amtx11,amtx12,amtx21,amtx22,imatrix 

IMPLICIT NONE
PRIVATE

PUBLIC :: find_channel_info, mapping_channel

CONTAINS

!===================================================================================
   SUBROUTINE FIND_CHANNEL_INFO (vGCM_st,channel)
!===================================================================================
!  Find nface (face number) of each channel segment based on lon/lat of vGCM grids,
!       channel-segment data size, and channel group ID.
!
!  INPUT : vGCM_st(:)%lon & lat
!
!  OUTPUT: vGCM_st(:)%nface, vGCM_st(:)%alpha1, vGCM_st(:)%beta1,
!          channel%num_chg, 
!          channel%seg(:)%nface, channel%seg(:)%nx_size, channel%seg(:)%ny_size 
!-----------------------------------------------------------------------------------
    type(vGCM_state_t), intent(inout) :: vGCM_st(4)  ! vGCM state data for one channel
    type(channel_t),    intent(inout) :: channel     ! CRM data for one channel
    
    ! Local
    INTEGER :: num_seg,MAX0,MIN0
    REAL (KIND=r8) :: LON_VAL,LAT_VAL
    
    DO num_seg = 1, 4
      LON_VAL = vGCM_st(num_seg)%lon(0,1)  ! Longitude [rad] of the 1st vGCM grid in each segment
      LAT_VAL = vGCM_st(num_seg)%lat(0,1)  ! Latitude [rad] of the 1st vGCM grid in each segment
    
      vGCM_st(num_seg)%nface = FIND_FACE ( LON_VAL, LAT_VAL )
      
      vGCM_st(num_seg)%alpha1 = ALPHAV ( vGCM_st(num_seg)%nface, LON_VAL, LAT_VAL, PI )
      vGCM_st(num_seg)%beta1  = BETAV ( vGCM_st(num_seg)%nface, LON_VAL, LAT_VAL, PI )
    ENDDO    
    
    MAX0 = MAXVAL(vGCM_st(:)%nface)
    MIN0 = MINVAL(vGCM_st(:)%nface)
    
    IF (MAX0.EQ.4) THEN
      channel%num_chg = 1        ! channel group id = 1 (F4, F1, F2, F3)
    ELSE IF (MIN0.EQ.1) THEN
      channel%num_chg = 2        ! channel group id = 2 (F6, F1, F5, F3)
    ELSE 
      channel%num_chg = 3        ! channel group id = 3 (F4, F5, F2, F6)
    ENDIF 
    
    DO num_seg = 1, 4
      channel%seg(num_seg)%nface = vGCM_st(num_seg)%nface
    ENDDO 
    
    ! # of prognostic grids
    SELECT CASE (channel%num_chg)
    CASE(1)
      DO num_seg = 1, 4
        channel%seg(num_seg)%nx_size = channel_seg_l
        channel%seg(num_seg)%ny_size = 1
      ENDDO  
    
    CASE(2)
      DO num_seg = 1, 4
        channel%seg(num_seg)%nx_size = 1
        channel%seg(num_seg)%ny_size = channel_seg_l     
      ENDDO     
    CASE(3)
      DO num_seg = 1, 4
       if (channel%seg(num_seg)%nface .GE. 5) then
         channel%seg(num_seg)%nx_size = channel_seg_l
         channel%seg(num_seg)%ny_size = 1
       else
         channel%seg(num_seg)%nx_size = 1
         channel%seg(num_seg)%ny_size = channel_seg_l  
       endif
      ENDDO  
    END SELECT   
    
   END SUBROUTINE find_channel_info
   
!===================================================================================
   SUBROUTINE MAPPING_CHANNEL(num_ch,vGCM_st,channel)
!===================================================================================
!  Initialize some parameters in constld (dx,dy,dxsq,dysq,dxdy,dx_gcm,dy_gcm,wrxmu)  
!  Obtain mapping information for a channel
!-----------------------------------------------------------------------------------
    INTEGER, intent(in) :: num_ch                 ! channel index
    type(vGCM_state_t), intent(in) :: vGCM_st(4)  ! vGCM state data for one channel
    type(channel_t), intent(inout) :: channel     ! Channel data

    ! Local variables 
    REAL (kind=r8) :: dalpha,dbeta,alpha_lowb,beta_lowb
    INTEGER :: num_seg,mi1,mim,mip,mj1,mjm,mjp
    INTEGER :: I,J,NPO
    
    REAL (KIND=r8), PARAMETER :: LAT_FIX = 15.0_r8   ! [deg] 

    dalpha = PI/(2.0_r8*channel_seg_l)
    dbeta  = dalpha
    
    IF (num_ch .EQ. 1) THEN 
!     Initialize the horizontal grid size of CRM (constld parameters)    
      DX = REARTH*DALPHA
      DY = REARTH*DBETA
    
      DXSQ = DX*DX
      DYSQ = DY*DY
      DXDY = DX*DY 
    
      DX_GCM = netsz*DX
      DY_GCM = netsz*DY
    
      WRXMU  = 2.0_r8/(DX*DX) 
    ENDIF
    
!**************************************************************************** 
      DO num_seg = 1, 4
!**************************************************************************** 
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim = channel%seg(num_seg)%mim    
      mip = channel%seg(num_seg)%mip    
    
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      mjm = channel%seg(num_seg)%mjm     
      mjp = channel%seg(num_seg)%mjp    

! Set up central angle coordinates for each channel segment
!=============================      
      IF (mi1.GT.mj1) THEN
!============================= x-channel-segment      

      channel%seg(num_seg)%CALPHA_T(1) = - 0.25_r8*PI + 0.5_r8*dalpha
      DO I = 2, mi1
       channel%seg(num_seg)%CALPHA_T(I) = channel%seg(num_seg)%CALPHA_T(I-1) + dalpha 
      ENDDO
      DO I = 1, mi1
       channel%seg(num_seg)%CALPHA_U(I) = channel%seg(num_seg)%CALPHA_T(I) + 0.5_r8*dalpha
      ENDDO

!     The x-channel-segment is centered at vGCM_st(num_seg)%beta1
      BETA_lowb = vGCM_st(num_seg)%beta1 - (0.5_r8*mj1*dbeta)
          
      DO J = 1, mj1
       channel%seg(num_seg)%CBETA_T(J) = BETA_lowb + (J-0.5_r8)*dbeta
      ENDDO
      DO J = 1, mj1
       channel%seg(num_seg)%CBETA_V(J) = channel%seg(num_seg)%CBETA_T(J) + 0.5_r8*dbeta
      ENDDO 
!=============================      
      ELSE
!============================= y-channel-segment  

      channel%seg(num_seg)%CBETA_T(1) = - 0.25_r8*PI + 0.5_r8*dbeta       
      DO J = 2, mj1
       channel%seg(num_seg)%CBETA_T(J) = channel%seg(num_seg)%CBETA_T(J-1) + dbeta 
      ENDDO
      DO J = 1, mj1
       channel%seg(num_seg)%CBETA_V(J) = channel%seg(num_seg)%CBETA_T(J) + 0.5_r8*dbeta 
      ENDDO
            
!     The y-channel-segment is centered at vGCM_st(num_seg)%alpha1
      alpha_lowb = vGCM_st(num_seg)%alpha1 - (0.5_r8*mi1*dalpha)
             
      DO I = 1, mi1
       channel%seg(num_seg)%CALPHA_T(I) = alpha_lowb + (I-0.5_r8)*dalpha
      ENDDO
      DO I = 1, mi1
       channel%seg(num_seg)%CALPHA_U(I) = channel%seg(num_seg)%CALPHA_T(I) + 0.5_r8*dalpha
      ENDDO 
!=============================      
      ENDIF
!=============================

!     Determine the coordinates of halo points    
      DO NPO = 1, nhalo
       channel%seg(num_seg)%CALPHA_T(  1-NPO) = channel%seg(num_seg)%CALPHA_T(  1)-dalpha*NPO
       channel%seg(num_seg)%CALPHA_T(mi1+NPO) = channel%seg(num_seg)%CALPHA_T(mi1)+dalpha*NPO 
       
       channel%seg(num_seg)%CALPHA_U(  1-NPO) = channel%seg(num_seg)%CALPHA_U(  1)-dalpha*NPO
       channel%seg(num_seg)%CALPHA_U(mi1+NPO) = channel%seg(num_seg)%CALPHA_U(mi1)+dalpha*NPO 
      ENDDO     
            
      DO NPO = 1, nhalo
       channel%seg(num_seg)%CBETA_T(  1-NPO) = channel%seg(num_seg)%CBETA_T(  1)-dbeta*NPO
       channel%seg(num_seg)%CBETA_T(mj1+NPO) = channel%seg(num_seg)%CBETA_T(mj1)+dbeta*NPO 
       
       channel%seg(num_seg)%CBETA_V(  1-NPO) = channel%seg(num_seg)%CBETA_V(  1)-dbeta*NPO
       channel%seg(num_seg)%CBETA_V(mj1+NPO) = channel%seg(num_seg)%CBETA_V(mj1)+dbeta*NPO 
      ENDDO 
      
      DO I = mim, mip
       channel%seg(num_seg)%CALPHA_V(I) = channel%seg(num_seg)%CALPHA_T(I) 
       channel%seg(num_seg)%CALPHA_Z(I) = channel%seg(num_seg)%CALPHA_U(I)     
      ENDDO

      DO J = mjm, mjp
       channel%seg(num_seg)%CBETA_U(J) = channel%seg(num_seg)%CBETA_T(J) 
       channel%seg(num_seg)%CBETA_Z(J) = channel%seg(num_seg)%CBETA_V(J)     
      ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (NOMAP) THEN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%FVAL(I,J) = 2.0_r8*OMEGA*DSIN(LAT_FIX/RAD2DEG) 
       
         channel%seg(num_seg)%RG_T(I,J) = 1.0_r8     
         channel%seg(num_seg)%RG_U(I,J) = 1.0_r8  
         channel%seg(num_seg)%RG_V(I,J) = 1.0_r8     
         channel%seg(num_seg)%RG_Z(I,J) = 1.0_r8 
       
         channel%seg(num_seg)%GCONT_T(1,I,J) = 1.0_r8
         channel%seg(num_seg)%GCONT_U(1,I,J) = 1.0_r8
         channel%seg(num_seg)%GCONT_V(1,I,J) = 1.0_r8
         channel%seg(num_seg)%GCONT_Z(1,I,J) = 1.0_r8
        
         channel%seg(num_seg)%GCONT_T(2,I,J) = 0.0_r8
         channel%seg(num_seg)%GCONT_U(2,I,J) = 0.0_r8
         channel%seg(num_seg)%GCONT_V(2,I,J) = 0.0_r8
         channel%seg(num_seg)%GCONT_Z(2,I,J) = 0.0_r8
        
         channel%seg(num_seg)%GCONT_T(3,I,J) = 1.0_r8
         channel%seg(num_seg)%GCONT_U(3,I,J) = 1.0_r8
         channel%seg(num_seg)%GCONT_V(3,I,J) = 1.0_r8
         channel%seg(num_seg)%GCONT_Z(3,I,J) = 1.0_r8
       ENDDO
      ENDDO       
      
      DO J=mjm,mjp
       DO I=mim,mip 
        channel%seg(num_seg)%GG_T(1,I,J) = 1.0_r8
        channel%seg(num_seg)%GG_U(1,I,J) = 1.0_r8
        channel%seg(num_seg)%GG_V(1,I,J) = 1.0_r8
        channel%seg(num_seg)%GG_Z(1,I,J) = 1.0_r8
        
        channel%seg(num_seg)%GG_T(2,I,J) = 0.0_r8
        channel%seg(num_seg)%GG_U(2,I,J) = 0.0_r8
        channel%seg(num_seg)%GG_V(2,I,J) = 0.0_r8
        channel%seg(num_seg)%GG_Z(2,I,J) = 0.0_r8
        
        channel%seg(num_seg)%GG_T(3,I,J) = 1.0_r8
        channel%seg(num_seg)%GG_U(3,I,J) = 1.0_r8
        channel%seg(num_seg)%GG_V(3,I,J) = 1.0_r8
        channel%seg(num_seg)%GG_Z(3,I,J) = 1.0_r8

        channel%seg(num_seg)%RGG_T(1,I,J) = 1.0_r8
        channel%seg(num_seg)%RGG_U(1,I,J) = 1.0_r8
        channel%seg(num_seg)%RGG_V(1,I,J) = 1.0_r8
        channel%seg(num_seg)%RGG_Z(1,I,J) = 1.0_r8
        
        channel%seg(num_seg)%RGG_T(2,I,J) = 0.0_r8
        channel%seg(num_seg)%RGG_U(2,I,J) = 0.0_r8
        channel%seg(num_seg)%RGG_V(2,I,J) = 0.0_r8
        channel%seg(num_seg)%RGG_Z(2,I,J) = 0.0_r8
        
        channel%seg(num_seg)%RGG_T(3,I,J) = 1.0_r8
        channel%seg(num_seg)%RGG_U(3,I,J) = 1.0_r8
        channel%seg(num_seg)%RGG_V(3,I,J) = 1.0_r8
        channel%seg(num_seg)%RGG_Z(3,I,J) = 1.0_r8        
       ENDDO
      ENDDO      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ELSE  ! NOMAP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Calculate mapping information

!---------------------------------------------------------------------------
!        CALCULATION OF Latitude and Longitude in [rad]
!        (Each face has different rules. Each process has different values.)
!---------------------------------------------------------------------------
      DO J = mjm, mjp
        DO I = mim, mip
          channel%seg(num_seg)%RLAT_T(I,J) = LATVAL(channel%seg(num_seg)%NFACE,       &
                                                    channel%seg(num_seg)%CALPHA_T(I), &
                                                    channel%seg(num_seg)%CBETA_T(J))
                                                    
          channel%seg(num_seg)%RLAT_U(I,J) = LATVAL(channel%seg(num_seg)%NFACE,       &
                                                    channel%seg(num_seg)%CALPHA_U(I), &
                                                    channel%seg(num_seg)%CBETA_U(J))
                                                    
          channel%seg(num_seg)%RLAT_V(I,J) = LATVAL(channel%seg(num_seg)%NFACE,       &
                                                    channel%seg(num_seg)%CALPHA_V(I), &
                                                    channel%seg(num_seg)%CBETA_V(J))
                                                    
          channel%seg(num_seg)%RLAT_Z(I,J) = LATVAL(channel%seg(num_seg)%NFACE,       &
                                                    channel%seg(num_seg)%CALPHA_Z(I), &
                                                    channel%seg(num_seg)%CBETA_Z(J)) 
        ENDDO
      ENDDO 
              
      DO J = mjm, mjp
        DO I = mim, mip
          channel%seg(num_seg)%RLON_T(I,J) = LONVAL(channel%seg(num_seg)%NFACE,       &
                                                    channel%seg(num_seg)%CALPHA_T(I), &
                                                    channel%seg(num_seg)%CBETA_T(J),PI) 
                                                    
          channel%seg(num_seg)%RLON_U(I,J) = LONVAL(channel%seg(num_seg)%NFACE,       &
                                                    channel%seg(num_seg)%CALPHA_U(I), &
                                                    channel%seg(num_seg)%CBETA_U(J),PI) 
                                                    
          channel%seg(num_seg)%RLON_V(I,J) = LONVAL(channel%seg(num_seg)%NFACE,       &
                                                    channel%seg(num_seg)%CALPHA_V(I), &
                                                    channel%seg(num_seg)%CBETA_V(J),PI) 
                                                    
          channel%seg(num_seg)%RLON_Z(I,J) = LONVAL(channel%seg(num_seg)%NFACE,       &
                                                    channel%seg(num_seg)%CALPHA_Z(I), &
                                                    channel%seg(num_seg)%CBETA_Z(J),PI) 
        ENDDO
      ENDDO

!---------------------------------------------------------------------------
!     Determine the coefficient of Coriolis Force defined at zeta-point 
!---------------------------------------------------------------------------  
      DO J = mjm, mjp
       DO I = mim, mip     
        channel%seg(num_seg)%FVAL(I,J) = &
           2.0_r8*OMEGA*DSIN(channel%seg(num_seg)%RLAT_Z(I,J)+PI/4.0_r8) 
       ENDDO
      ENDDO

!---------------------------------------------------------------------------     
!     Determine the Jacobian of the transformation: SQRT(G)  
!---------------------------------------------------------------------------
!     RG=1 at alpha=beta=0, RG~0.71 at alpha=pi/4 & beta=0, RG~0.77 at alpha=beta=pi/4

      DO J = mjm, mjp
       DO I = mim, mip 
        channel%seg(num_seg)%RG_T(I,J) = JACO(channel%seg(num_seg)%CALPHA_T(I), &
                                              channel%seg(num_seg)%CBETA_T(J))   
        channel%seg(num_seg)%RG_U(I,J) = JACO(channel%seg(num_seg)%CALPHA_U(I), &
                                              channel%seg(num_seg)%CBETA_U(J))       
        channel%seg(num_seg)%RG_V(I,J) = JACO(channel%seg(num_seg)%CALPHA_V(I), &
                                              channel%seg(num_seg)%CBETA_V(J))      
        channel%seg(num_seg)%RG_Z(I,J) = JACO(channel%seg(num_seg)%CALPHA_Z(I), &
                                              channel%seg(num_seg)%CBETA_Z(J))  
       ENDDO
      ENDDO
      
!---------------------------------------------------------------------------
!     Determine the coefficients of the contravariant metric tensor
!     1: g11   2: g12=g21  3: g22
!---------------------------------------------------------------------------
      DO J = mjm, mjp
       DO I = mim, mip 
        channel%seg(num_seg)%GCONT_T(1,I,J) = GCONT11(channel%seg(num_seg)%CALPHA_T(I), &
                                                      channel%seg(num_seg)%CBETA_T(J))
        channel%seg(num_seg)%GCONT_U(1,I,J) = GCONT11(channel%seg(num_seg)%CALPHA_U(I), &
                                                      channel%seg(num_seg)%CBETA_U(J))
        channel%seg(num_seg)%GCONT_V(1,I,J) = GCONT11(channel%seg(num_seg)%CALPHA_V(I), &
                                                      channel%seg(num_seg)%CBETA_V(J))
        channel%seg(num_seg)%GCONT_Z(1,I,J) = GCONT11(channel%seg(num_seg)%CALPHA_Z(I), &
                                                      channel%seg(num_seg)%CBETA_Z(J))
        
        channel%seg(num_seg)%GCONT_T(2,I,J) = GCONT12(channel%seg(num_seg)%CALPHA_T(I), &
                                                      channel%seg(num_seg)%CBETA_T(J))
        channel%seg(num_seg)%GCONT_U(2,I,J) = GCONT12(channel%seg(num_seg)%CALPHA_U(I), &
                                                      channel%seg(num_seg)%CBETA_U(J))
        channel%seg(num_seg)%GCONT_V(2,I,J) = GCONT12(channel%seg(num_seg)%CALPHA_V(I), &
                                                      channel%seg(num_seg)%CBETA_V(J))
        channel%seg(num_seg)%GCONT_Z(2,I,J) = GCONT12(channel%seg(num_seg)%CALPHA_Z(I), &
                                                      channel%seg(num_seg)%CBETA_Z(J))
        
        channel%seg(num_seg)%GCONT_T(3,I,J) = GCONT22(channel%seg(num_seg)%CALPHA_T(I), &
                                                      channel%seg(num_seg)%CBETA_T(J))
        channel%seg(num_seg)%GCONT_U(3,I,J) = GCONT22(channel%seg(num_seg)%CALPHA_U(I), &
                                                      channel%seg(num_seg)%CBETA_U(J))
        channel%seg(num_seg)%GCONT_V(3,I,J) = GCONT22(channel%seg(num_seg)%CALPHA_V(I), &
                                                      channel%seg(num_seg)%CBETA_V(J))
        channel%seg(num_seg)%GCONT_Z(3,I,J) = GCONT22(channel%seg(num_seg)%CALPHA_Z(I), &
                                                      channel%seg(num_seg)%CBETA_Z(J))
       ENDDO
      ENDDO
      
      DO J = mjm, mjp
       DO I = mim, mip 
        channel%seg(num_seg)%RGG_T(1,I,J) = channel%seg(num_seg)%RG_T(I,J)      &
                                           *channel%seg(num_seg)%GCONT_T(1,I,J)
        channel%seg(num_seg)%RGG_U(1,I,J) = channel%seg(num_seg)%RG_U(I,J)      &
                                           *channel%seg(num_seg)%GCONT_U(1,I,J)
        channel%seg(num_seg)%RGG_V(1,I,J) = channel%seg(num_seg)%RG_V(I,J)      &
                                           *channel%seg(num_seg)%GCONT_V(1,I,J)
        channel%seg(num_seg)%RGG_Z(1,I,J) = channel%seg(num_seg)%RG_Z(I,J)      &
                                           *channel%seg(num_seg)%GCONT_Z(1,I,J)
        
        channel%seg(num_seg)%RGG_T(2,I,J) = channel%seg(num_seg)%RG_T(I,J)      &
                                           *channel%seg(num_seg)%GCONT_T(2,I,J)
        channel%seg(num_seg)%RGG_U(2,I,J) = channel%seg(num_seg)%RG_U(I,J)      &
                                           *channel%seg(num_seg)%GCONT_U(2,I,J)
        channel%seg(num_seg)%RGG_V(2,I,J) = channel%seg(num_seg)%RG_V(I,J)      &
                                           *channel%seg(num_seg)%GCONT_V(2,I,J)
        channel%seg(num_seg)%RGG_Z(2,I,J) = channel%seg(num_seg)%RG_Z(I,J)      &
                                           *channel%seg(num_seg)%GCONT_Z(2,I,J)
        
        channel%seg(num_seg)%RGG_T(3,I,J) = channel%seg(num_seg)%RG_T(I,J)      & 
                                           *channel%seg(num_seg)%GCONT_T(3,I,J)
        channel%seg(num_seg)%RGG_U(3,I,J) = channel%seg(num_seg)%RG_U(I,J)      &
                                           *channel%seg(num_seg)%GCONT_U(3,I,J)
        channel%seg(num_seg)%RGG_V(3,I,J) = channel%seg(num_seg)%RG_V(I,J)      &
                                           *channel%seg(num_seg)%GCONT_V(3,I,J)
        channel%seg(num_seg)%RGG_Z(3,I,J) = channel%seg(num_seg)%RG_Z(I,J)      &
                                           *channel%seg(num_seg)%GCONT_Z(3,I,J)
       ENDDO
      ENDDO      

      DO J = mjm, mjp
       DO I = mim, mip 
        channel%seg(num_seg)%GG_T(1,I,J) = channel%seg(num_seg)%RG_T(I,J)       & 
                                          *channel%seg(num_seg)%RGG_T(1,I,J)
        channel%seg(num_seg)%GG_U(1,I,J) = channel%seg(num_seg)%RG_U(I,J)       &  
                                          *channel%seg(num_seg)%RGG_U(1,I,J)
        channel%seg(num_seg)%GG_V(1,I,J) = channel%seg(num_seg)%RG_V(I,J)       &
                                          *channel%seg(num_seg)%RGG_V(1,I,J)
        channel%seg(num_seg)%GG_Z(1,I,J) = channel%seg(num_seg)%RG_Z(I,J)       &
                                          *channel%seg(num_seg)%RGG_Z(1,I,J)
        
        channel%seg(num_seg)%GG_T(2,I,J) = channel%seg(num_seg)%RG_T(I,J)       &
                                          *channel%seg(num_seg)%RGG_T(2,I,J)
        channel%seg(num_seg)%GG_U(2,I,J) = channel%seg(num_seg)%RG_U(I,J)       &
                                          *channel%seg(num_seg)%RGG_U(2,I,J)
        channel%seg(num_seg)%GG_V(2,I,J) = channel%seg(num_seg)%RG_V(I,J)       &
                                          *channel%seg(num_seg)%RGG_V(2,I,J)
        channel%seg(num_seg)%GG_Z(2,I,J) = channel%seg(num_seg)%RG_Z(I,J)       &
                                          *channel%seg(num_seg)%RGG_Z(2,I,J)
        
        channel%seg(num_seg)%GG_T(3,I,J) = channel%seg(num_seg)%RG_T(I,J)       &
                                          *channel%seg(num_seg)%RGG_T(3,I,J)
        channel%seg(num_seg)%GG_U(3,I,J) = channel%seg(num_seg)%RG_U(I,J)       &
                                          *channel%seg(num_seg)%RGG_U(3,I,J)
        channel%seg(num_seg)%GG_V(3,I,J) = channel%seg(num_seg)%RG_V(I,J)       &
                                          *channel%seg(num_seg)%RGG_V(3,I,J)
        channel%seg(num_seg)%GG_Z(3,I,J) = channel%seg(num_seg)%RG_Z(I,J)       &
                                          *channel%seg(num_seg)%RGG_Z(3,I,J)
       ENDDO
      ENDDO            
       
!---------------------------------------------------------------------------
!        CALCULATION OF Transformation matrix: A 
!        (Each face has different rules. Each process has different values.)
!---------------------------------------------------------------------------

      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%AM_T(1,I,J) = AMTX11(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_T(I), &
                                                  channel%seg(num_seg)%CBETA_T(J))
                                                  
        channel%seg(num_seg)%AM_U(1,I,J) = AMTX11(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_U(I), &
                                                  channel%seg(num_seg)%CBETA_U(J))
                                                  
        channel%seg(num_seg)%AM_V(1,I,J) = AMTX11(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_V(I), &
                                                  channel%seg(num_seg)%CBETA_V(J))
                                                  
        channel%seg(num_seg)%AM_Z(1,I,J) = AMTX11(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_Z(I), &
                                                  channel%seg(num_seg)%CBETA_Z(J))
        
        channel%seg(num_seg)%AM_T(2,I,J) = AMTX12(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_T(I), &
                                                  channel%seg(num_seg)%CBETA_T(J))
                                                  
        channel%seg(num_seg)%AM_U(2,I,J) = AMTX12(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_U(I), &
                                                  channel%seg(num_seg)%CBETA_U(J))
                                                  
        channel%seg(num_seg)%AM_V(2,I,J) = AMTX12(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_V(I), &
                                                  channel%seg(num_seg)%CBETA_V(J))
                                                  
        channel%seg(num_seg)%AM_Z(2,I,J) = AMTX12(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_Z(I), &
                                                  channel%seg(num_seg)%CBETA_Z(J))
        
        channel%seg(num_seg)%AM_T(3,I,J) = AMTX21(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_T(I), &
                                                  channel%seg(num_seg)%CBETA_T(J))
                                                  
        channel%seg(num_seg)%AM_U(3,I,J) = AMTX21(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_U(I), &
                                                  channel%seg(num_seg)%CBETA_U(J))
                                                  
        channel%seg(num_seg)%AM_V(3,I,J) = AMTX21(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_V(I), &
                                                  channel%seg(num_seg)%CBETA_V(J))
                                                  
        channel%seg(num_seg)%AM_Z(3,I,J) = AMTX21(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_Z(I), &
                                                  channel%seg(num_seg)%CBETA_Z(J))
        
        channel%seg(num_seg)%AM_T(4,I,J) = AMTX22(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_T(I), &
                                                  channel%seg(num_seg)%CBETA_T(J))
                                                  
        channel%seg(num_seg)%AM_U(4,I,J) = AMTX22(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_U(I), &
                                                  channel%seg(num_seg)%CBETA_U(J))
                                                  
        channel%seg(num_seg)%AM_V(4,I,J) = AMTX22(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_V(I), &
                                                  channel%seg(num_seg)%CBETA_V(J))
                                                  
        channel%seg(num_seg)%AM_Z(4,I,J) = AMTX22(channel%seg(num_seg)%NFACE,       &
                                                  channel%seg(num_seg)%CALPHA_Z(I), &
                                                  channel%seg(num_seg)%CBETA_Z(J))
       ENDDO
      ENDDO            

!---------------------------------------------------------------------------
!     CALCULATION OF INVERSE Transformation matrix: A^-1 
!     In each face, AA^-1 must be Identity matrix. Thus, function cal_inverse 
!     is used, instead of the analytical expression.    
!---------------------------------------------------------------------------   

      CALL CAL_INVERSE (mim,mip,mjm,mjp, &
                        channel%seg(num_seg)%AM_T,channel%seg(num_seg)%AMI_T)
      CALL CAL_INVERSE (mim,mip,mjm,mjp, &
                        channel%seg(num_seg)%AM_U,channel%seg(num_seg)%AMI_U)
      CALL CAL_INVERSE (mim,mip,mjm,mjp, &
                        channel%seg(num_seg)%AM_V,channel%seg(num_seg)%AMI_V)
      CALL CAL_INVERSE (mim,mip,mjm,mjp, &
                        channel%seg(num_seg)%AM_Z,channel%seg(num_seg)%AMI_Z)
      
!---------------------------------------------------------------------------
!     Calculate alpha & beta values of halo points in the neighboring face.
!---------------------------------------------------------------------------

      CALL CAL_HALO_LOC (mim,mip,mjm,mjp,channel%seg(num_seg)%NFACE,                  &
                         DALPHA,DBETA,                                                &
                         channel%seg(num_seg)%CALPHA_T,channel%seg(num_seg)%CALPHA_Z, &
                         channel%seg(num_seg)%CBETA_T, channel%seg(num_seg)%CBETA_Z,  &
                         channel%seg(num_seg)%HALPHA_T,channel%seg(num_seg)%HALPHA_Z, &
                         channel%seg(num_seg)%HBETA_T, channel%seg(num_seg)%HBETA_Z)

!---------------------------------------------------------------------------
!     Find i- and j- index for 1-D halo interpolation
!---------------------------------------------------------------------------

      CALL FIND_HALO_LOC (mi1,mj1,mim,mip,mjm,mjp,                                     &
                          channel%seg(num_seg)%CALPHA_T,channel%seg(num_seg)%CALPHA_Z, &
                          channel%seg(num_seg)%CBETA_T, channel%seg(num_seg)%CBETA_Z,  &
                          channel%seg(num_seg)%HALPHA_T,channel%seg(num_seg)%HALPHA_Z, &
                          channel%seg(num_seg)%HBETA_T, channel%seg(num_seg)%HBETA_Z,  &
                          channel%seg(num_seg)%IHALO_LOC_T, &
                          channel%seg(num_seg)%IHALO_LOC_Z, &
                          channel%seg(num_seg)%JHALO_LOC_T, &
                          channel%seg(num_seg)%JHALO_LOC_Z)       

!---------------------------------------------------------------------------
!     CALCULATION of transformation matrix for vector halos at q-points 
!--------------------------------------------------------------------------- 

      CALL CAL_HALO_AMATRIX (mim,mip,mjm,mjp,channel%seg(num_seg)%NFACE,DALPHA,DBETA,    &
                             channel%seg(num_seg)%HALPHA_T,channel%seg(num_seg)%HBETA_T, &
                             channel%seg(num_seg)%AM_VORT_ALPHA, &
                             channel%seg(num_seg)%AM_VORT_BETA)  
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ENDIF  ! NOMAP  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
     
!********************************   
       ENDDO  ! num_seg 
!********************************       

       CONTAINS      

!===============================================================================
      SUBROUTINE CAL_INVERSE (mim,mip,mjm,mjp,DumA,DumB)
!===============================================================================    
!     Calculate the inverse matrix  
      
      INTEGER, INTENT(IN) :: mim,mip,mjm,mjp
      REAL (kind=r8), DIMENSION(4,mim:mip,mjm:mjp), INTENT(IN)  :: DumA
      REAL (kind=r8), DIMENSION(4,mim:mip,mjm:mjp), INTENT(OUT) :: DumB
      
      ! Local
      REAL (kind=r8) :: TEMP(2,2),TEMPI(2,2) 
      INTEGER :: I,J   
      
      DO J = mjm, mjp
       DO I = mim, mip
        TEMP(1,1) = DumA(1,I,J)
        TEMP(1,2) = DumA(2,I,J)
        TEMP(2,1) = DumA(3,I,J)
        TEMP(2,2) = DumA(4,I,J)
        
        TEMPI = IMATRIX(TEMP)
        
        DumB(1,I,J) = TEMPI(1,1)
        DumB(2,I,J) = TEMPI(1,2)
        DumB(3,I,J) = TEMPI(2,1)
        DumB(4,I,J) = TEMPI(2,2)
       ENDDO
      ENDDO
            
      END SUBROUTINE cal_inverse

!===============================================================================
      SUBROUTINE CAL_HALO_LOC (mim,mip,mjm,mjp,NFACE,             &
                               DALPHA,DBETA,                      &
                               CALPHA_T,CALPHA_Z,CBETA_T,CBETA_Z, &
                               HALPHA_T,HALPHA_Z,HBETA_T,HBETA_Z)
!===============================================================================                               
!     Determine the locations of halo points near the face edges
!     HALPHA_T, HALPHA_Z, HBETA_T, HBETA_Z
      
      INTEGER, INTENT(IN) :: mim,mip,mjm,mjp,NFACE
      REAL (kind=r8), INTENT(IN) :: DALPHA,DBETA
      REAL (kind=r8), DIMENSION(mim:mip), INTENT(IN) :: CALPHA_T,CALPHA_Z
      REAL (kind=r8), DIMENSION(mjm:mjp), INTENT(IN) :: CBETA_T,CBETA_Z
      
      REAL (kind=r8), DIMENSION(mim:mip,nhalo), INTENT(OUT) :: HALPHA_T,HALPHA_Z
      REAL (kind=r8), DIMENSION(nhalo,mjm:mjp), INTENT(OUT) :: HBETA_T,HBETA_Z
      
      ! Local variables
      INTEGER :: I,J,IFACE
      REAL (kind=r8) :: RLON,RLAT,ALPHA,BETA

!===========================================
      SELECT CASE (NFACE)
!===========================================   
!---------------     
        CASE(1)
!---------------   
! OUT(FACE5):       LAT < 35   vs.   OUT(FACE5): -45 < LON < 45   
!  IN(FACE1): -35 < LAT < 35   vs.    IN(FACE1): -45 < LON < 45    
!  no correction of the range of lat & lon
 
        IFACE = 5
        DO I = mim, mip
         ALPHA = CALPHA_T(I)
         DO J = 1, nhalo
          BETA = -0.25_r8*PI + 0.5_r8*DBETA - DBETA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_T(I,J) = ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO I=mim,mip
         ALPHA = CALPHA_Z(I)
         DO J = 1, nhalo
          BETA = -0.25_r8*PI - DBETA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_Z(I,J)= ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO            

! OUT(FACE2): -35 < LAT < 35     vs. OUT(FACE2):       LON < 45   
!  IN(FACE1): -35 < LAT < 35     vs.  IN(FACE1): -45 < LON < 45   
!  no correction of the range of lat & lon
        
        IFACE = 2
        DO J = mjm, mjp
         BETA = CBETA_T(J)
         DO I = 1, nhalo
          ALPHA = -0.25_r8*PI + 0.5_r8*DALPHA - DALPHA*I            
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_T(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO J = mjm, mjp
         BETA = CBETA_Z(J)
         DO I = 1, nhalo
          ALPHA = -0.25_r8*PI - DALPHA*I         
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_Z(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO       
!---------------            
        CASE(2)
!--------------- 
! OUT(FACE5):       LAT < 35   vs.   OUT(FACE5): 45 < LON < 135
!  IN(FACE2): -35 < LAT < 35   vs.    IN(FACE2): 45 < LON < 135 
!  no correction of the range of lat & lon
 
        IFACE = 5
        DO I = mim, mip
         BETA = CALPHA_T(I)
         DO J = 1, nhalo
          ALPHA = 0.25_r8*PI - 0.5_r8*DALPHA + DALPHA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_T(I,J) = ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO I = mim, mip
         BETA = CALPHA_Z(I)
         DO J = 1, nhalo
          ALPHA = 0.25_r8*PI + DALPHA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_Z(I,J)= ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO            

! OUT(FACE3): -35 < LAT < 35     vs. OUT(FACE3):      LON < 135   
!  IN(FACE1): -35 < LAT < 35     vs.  IN(FACE1): 45 < LON < 135   
!  no correction of the range of lat & lon
        
        IFACE = 3
        DO J = mjm, mjp
         BETA = CBETA_T(J)
         DO I = 1, nhalo
          ALPHA = -0.25_r8*PI + 0.5_r8*DALPHA - DALPHA*I            
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_T(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO J = mjm, mjp
         BETA = CBETA_Z(J)
         DO I = 1, nhalo
          ALPHA = -0.25_r8*PI - DALPHA*I         
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_Z(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO    
!---------------            
        CASE(3)
!---------------  
! OUT(FACE5):       LAT < 35   vs.   OUT(FACE5): 135 < LON < 180; -135 > LON > -180
!  IN(FACE3): -35 < LAT < 35   vs.    IN(FACE3): 135 < LON < 225
!  correction of the range of lon

        IFACE = 5
        DO I = mim, mip
         ALPHA = - CALPHA_T(I)
         DO J = 1, nhalo
          BETA = 0.25_r8*PI - 0.5_r8*DBETA + DBETA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          IF (RLON.LT.0.0_r8) RLON = RLON+2.*PI
          
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_T(I,J) = ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO I = mim, mip
         ALPHA = - CALPHA_Z(I)
         DO J = 1, nhalo
          BETA = 0.25_r8*PI + DBETA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI)
          IF (RLON.LT.0.0_r8) RLON = RLON+2.0_r8*PI 
          
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_Z(I,J)= ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO            

! OUT(FACE4): -35 < LAT < 35     vs. OUT(FACE4):       LON < 225  
!  IN(FACE3): -35 < LAT < 35     vs.  IN(FACE3): 135 < LON < 225   
!  no correction of the range of lat & lon
        
        IFACE = 4
        DO J = mjm, mjp
         BETA = CBETA_T(J)
         DO I = 1, nhalo
          ALPHA = -0.25_r8*PI + 0.5_r8*DALPHA - DALPHA*I            
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_T(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO J = mjm, mjp
         BETA = CBETA_Z(J)
         DO I = 1, nhalo
          ALPHA = -0.25_r8*PI - DALPHA*I         
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_Z(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO    
!---------------            
        CASE(4)
!--------------- 
! OUT(FACE5):       LAT < 35   vs.   OUT(FACE5): -135 < LON < -45
!  IN(FACE4): -35 < LAT < 35   vs.    IN(FACE4):  225 < LON < 315
!  correction of the range of lon

        IFACE = 5
        DO I = mim, mip
         BETA = - CALPHA_T(I)
         DO J = 1, nhalo
          ALPHA = -0.25_r8*PI + 0.5_r8*DALPHA - DALPHA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLON = RLON+2.*PI
          
          IF (RLON.LT.0.0_r8) RLON = RLON+2.*PI
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_T(I,J) = ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO I = mim, mip
         BETA = - CALPHA_Z(I)
         DO J = 1, nhalo
          ALPHA = -0.25_r8*PI - DALPHA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLON = RLON+2.0_r8*PI
          
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_Z(I,J)= ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO            

! OUT(FACE1): -35 < LAT <  35   vs. OUT(FACE1):       LON < -45 
!  IN(FACE4): -35 < LAT <  35   vs.  IN(FACE4): 225 < LON < 315  
!  correction of the range of lon
        
        IFACE = 1
        DO J = mjm, mjp
         BETA = CBETA_T(J)
         DO I = 1, nhalo
          ALPHA = -0.25_r8*PI + 0.5_r8*DALPHA - DALPHA*I            
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLON  = RLON+2.0_r8*PI
          
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_T(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO J = mjm, mjp
         BETA = CBETA_Z(J)
         DO I = 1, nhalo
          ALPHA = -0.25_r8*PI - DALPHA*I         
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLON  = RLON+2.0_r8*PI
          
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_Z(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO    
!---------------            
        CASE(5)
!--------------- 
! OUT(FACE3):  35 < LAT        vs.   OUT(FACE3): 135 < LON < 225
!  IN(FACE5):  35 < LAT < 90   vs.    IN(FACE5): 135 < LON < 180; -135 > LON > -180
!  correction of the range of lon

        IFACE = 3 
        DO I = mim, mip
         ALPHA = - CALPHA_T(I)
         DO J = 1, nhalo
          BETA = 0.25_r8*PI - 0.5_r8*DBETA + DBETA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          IF (RLON.GT.PI) RLON = RLON-2.*PI
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_T(I,J) = ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO I = mim, mip
         ALPHA = - CALPHA_Z(I) 
         DO J = 1, nhalo
          BETA = 0.25_r8*PI + DBETA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          IF (RLON.GT.PI) RLON = RLON-2.0_r8*PI
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_Z(I,J)= ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO            
        
! OUT(FACE2):  35 < LAT        vs.   OUT(FACE2): 45 < LON < 135
!  IN(FACE5):  35 < LAT < 90   vs.    IN(FACE5): 45 < LON < 135
!  no correction of the range of lat & lon
        
        IFACE = 2 
        DO J = mjm, mjp
         ALPHA = CBETA_T(J)
         DO I = 1, nhalo
          BETA = 0.25_r8*PI - 0.5_r8*DBETA + DBETA*I            
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_T(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO J = mjm, mjp
         ALPHA = CBETA_Z(J)
         DO I = 1, nhalo
          BETA = 0.25_r8*PI + DBETA*I         
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_Z(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO    
!---------------            
        CASE(6)
!---------------     
! OUT(FACE1):        LAT < -35   vs.  OUT(FACE1): -45 < LON < 45
!  IN(FACE6):  -90 < LAT < -35   vs.   IN(FACE6): -45 < LON < 45
!  no correction of the range of lat & lon
   
        IFACE = 1 
        DO I = mim, mip
         ALPHA = CALPHA_T(I)
         DO J = 1, nhalo
          BETA = -0.25_r8*PI + 0.5_r8*DBETA - DBETA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_T(I,J) = ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO I = mim, mip
         ALPHA = CALPHA_Z(I)
         DO J = 1, nhalo
          BETA = -0.25_r8*PI - DBETA*J         
          RLON = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT = LATVAL(IFACE,ALPHA,BETA)
          HALPHA_Z(I,J)= ALPHAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO            
        
! OUT(FACE2):        LAT < -35   vs.  OUT(FACE2): 45 < LON < 135
!  IN(FACE6):  -90 < LAT < -35   vs.   IN(FACE6): 45 < LON < 135
!  no correction of the range of lat & lon
        
        IFACE = 2 
        DO J = mjm, mjp
         ALPHA = - CBETA_T(J)
         DO I=1,nhalo
          BETA = -0.25_r8*PI + 0.5_r8*DBETA - DBETA*I            
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_T(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO  
        DO J = mjm, mjp
         ALPHA = - CBETA_Z(J)
         DO I = 1, nhalo
          BETA = -0.25_r8*PI - DBETA*I         
          RLON  = LONVAL(IFACE,ALPHA,BETA,PI) 
          RLAT  = LATVAL(IFACE,ALPHA,BETA)
          HBETA_Z(I,J)= BETAV(NFACE,RLON,RLAT,PI)
         ENDDO
        ENDDO              
      
!===========================================          
      END SELECT   
!=========================================== 
      
      END SUBROUTINE cal_halo_loc 

!===============================================================================      
      SUBROUTINE FIND_HALO_LOC (mi1,mj1,mim,mip,mjm,mjp,           &
                                CALPHA_T,CALPHA_Z,CBETA_T,CBETA_Z, &
                                HALPHA_T,HALPHA_Z,HBETA_T,HBETA_Z, &
                                IHALO_LOC_T,IHALO_LOC_Z,           &
                                JHALO_LOC_T,JHALO_LOC_Z) 
!===============================================================================      
!     Calculate mapping rules   

      INTEGER, INTENT(IN) :: mi1,mj1,mim,mip,mjm,mjp

      REAL (kind=r8), DIMENSION(mim:mip), INTENT(IN) :: CALPHA_T,CALPHA_Z
      REAL (kind=r8), DIMENSION(mjm:mjp), INTENT(IN) :: CBETA_T,CBETA_Z
      
      REAL (kind=r8), DIMENSION(mim:mip,nhalo), INTENT(IN) :: HALPHA_T,HALPHA_Z
      REAL (kind=r8), DIMENSION(nhalo,mjm:mjp), INTENT(IN) :: HBETA_T,HBETA_Z
      
      INTEGER, DIMENSION(mim:mip,nhalo), INTENT(OUT) :: IHALO_LOC_T,IHALO_LOC_Z
      INTEGER, DIMENSION(nhalo,mjm:mjp), INTENT(OUT) :: JHALO_LOC_T,JHALO_LOC_Z
            
      ! Logical variables      
      LOGICAL :: LF
      INTEGER :: I,J,NE_JV,NE_IV  
      
      NE_JV = mj1 + nhalo*2
      NE_IV = mi1 + nhalo*2 
      
      DO J = mjm, mjp
       DO I = 1, nhalo
        JHALO_LOC_T(I,J) = INDEXR(HBETA_T(I,J), NE_JV, CBETA_T, LF)
        JHALO_LOC_Z(I,J) = INDEXR(HBETA_Z(I,J), NE_JV, CBETA_Z, LF)
       ENDDO
       DO I = 1, nhalo
        JHALO_LOC_T(I,J) = JHALO_LOC_T(I,J) - nhalo 
        JHALO_LOC_Z(I,J) = JHALO_LOC_Z(I,J) - nhalo 
       ENDDO
      ENDDO       

      DO J = 1, nhalo
       DO I = mim, mip
        IHALO_LOC_T(I,J) = INDEXR(HALPHA_T(I,J), NE_IV, CALPHA_T, LF)
        IHALO_LOC_Z(I,J) = INDEXR(HALPHA_Z(I,J), NE_IV, CALPHA_Z, LF)
       ENDDO
       DO I = mim, mip
        IHALO_LOC_T(I,J) = IHALO_LOC_T(I,J) - nhalo
        IHALO_LOC_Z(I,J) = IHALO_LOC_Z(I,J) - nhalo
       ENDDO
      ENDDO       
            
      END SUBROUTINE find_halo_loc

      SUBROUTINE CAL_HALO_AMATRIX (mim,mip,mjm,mjp,NFACE,         &
                                   DALPHA,DBETA,HALPHA_T,HBETA_T, &
                                   AM_VORT_ALPHA,AM_VORT_BETA)
!***************************************************************************
!     CALCULATION of transformation matrix for vector halos
!     
!     Transformation matrix for internal halo points (used in VORT_COMM_PRE)
!     Internal halo points, which are not the coordinate grids (coordinates of other face):    
!                  | 1 -> nhalo      nhalo+1 -> 2*nhalo |
!*************************************************************************** 
      INTEGER, INTENT(IN) :: mim,mip,mjm,mjp,NFACE
      
      REAL (kind=r8), INTENT(IN) :: DALPHA,DBETA
      REAL (kind=r8), DIMENSION(mim:mip,nhalo), INTENT(IN) :: HALPHA_T
      REAL (kind=r8), DIMENSION(nhalo,mjm:mjp), INTENT(IN) :: HBETA_T
      
      REAL (kind=r8), DIMENSION(4,mim:mip,nhalo*2), INTENT(OUT) :: AM_VORT_ALPHA
      REAL (kind=r8), DIMENSION(4,nhalo*2,mjm:mjp), INTENT(OUT) :: AM_VORT_BETA
      
      ! Local variables
      INTEGER :: I,J,NPO
      REAL (kind=r8) :: ALPHA,BETA

!-------------------------------------------------------       
      DO J = mjm, mjp
!-------------------------------------------------------       
      
       DO I = 1, nhalo
       BETA = HBETA_T(I,J)
       
!       internal halo points near alpha = -pi/4        
        ALPHA = -0.25_r8*PI - 0.5_r8*DALPHA + I*DALPHA
        AM_VORT_BETA(1,I,J) = AMTX11(NFACE,ALPHA,BETA)
        AM_VORT_BETA(2,I,J) = AMTX12(NFACE,ALPHA,BETA)
        AM_VORT_BETA(3,I,J) = AMTX21(NFACE,ALPHA,BETA)
        AM_VORT_BETA(4,I,J) = AMTX22(NFACE,ALPHA,BETA)

        NPO = (2*nhalo+1) - I
!       internal halo points near alpha = pi/4       
        ALPHA = 0.25_r8*PI + 0.5_r8*DALPHA - I*DALPHA  
        AM_VORT_BETA(1,NPO,J) = AMTX11(NFACE,ALPHA,BETA)
        AM_VORT_BETA(2,NPO,J) = AMTX12(NFACE,ALPHA,BETA)
        AM_VORT_BETA(3,NPO,J) = AMTX21(NFACE,ALPHA,BETA)
        AM_VORT_BETA(4,NPO,J) = AMTX22(NFACE,ALPHA,BETA)   
       ENDDO

!---------------------------        
      ENDDO  ! j-loop      
!---------------------------       

!------------------------------------------------------- 
      DO I = mim, mip
!-------------------------------------------------------       

       DO J = 1, nhalo
       ALPHA = HALPHA_T(I,J)
        
!       internal halo points near beta = -pi/4 
        BETA = -0.25_r8*PI - 0.5_r8*DBETA + J*DBETA
        AM_VORT_ALPHA(1,I,J) = AMTX11(NFACE,ALPHA,BETA)
        AM_VORT_ALPHA(2,I,J) = AMTX12(NFACE,ALPHA,BETA)
        AM_VORT_ALPHA(3,I,J) = AMTX21(NFACE,ALPHA,BETA)
        AM_VORT_ALPHA(4,I,J) = AMTX22(NFACE,ALPHA,BETA)

        NPO = (2*nhalo+1) - J
!       internal halo points near beta = pi/4       
        BETA =  0.25_r8*PI + 0.5_r8*DBETA - J*DBETA  
        AM_VORT_ALPHA(1,I,NPO) = AMTX11(NFACE,ALPHA,BETA)
        AM_VORT_ALPHA(2,I,NPO) = AMTX12(NFACE,ALPHA,BETA)
        AM_VORT_ALPHA(3,I,NPO) = AMTX21(NFACE,ALPHA,BETA)
        AM_VORT_ALPHA(4,I,NPO) = AMTX22(NFACE,ALPHA,BETA)   
       ENDDO
!---------------------------        
      ENDDO  ! i-loop        
!---------------------------             
            
      END SUBROUTINE cal_halo_amatrix 

   END SUBROUTINE mapping_channel

END MODULE q3d_mapping_module
