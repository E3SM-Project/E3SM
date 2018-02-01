MODULE q3d_interface_module
! Subprograms control the VVM channels

USE shr_kind_mod,    only: dbl_kind => shr_kind_r8
USE vGCM_data_types, only: vGCM_state_t,vGCM_tend_t,vGCM_out_t,channel_vGCM_t
USE vvm_data_types,  only: channel_t,allocate_channel_data

USE parmsld, only: nVGCM_seg,netsz,nhalo_cal,nhalo_adv,nhalo,ntracer
USE constld, only: d0_0,d1_0,d2_0,d3_0,d0_5,itt_max,a,b,dt,physics      

#ifdef FIXTHIS
USE time_manager_module, only: time_manager
#endif

USE q3d_init_module,    only: init_common,init_channel
USE q3d_mapping_module, only: find_channel_info,mapping_channel
USE q3d_bg_module,      only: cal_bg
USE q3d_effect_module,  only: ini_tendency,ini_crm_eft,cal_eddy_trans, &
                              cal_tendency,ini_sfc_eft,cal_sfc_eft, &
                              eddy_prepare
USE q3d_rtime_module,   only: cal_rtime

IMPLICIT NONE
PRIVATE

public :: crm_init
public :: crm_run
public :: crm_final

CONTAINS
!==============================================================================
  subroutine crm_init (isRestart, nch, channel, vGCM)
!==============================================================================
    LOGICAL, INTENT(IN) :: isRestart  ! .true. iff restart run
    INTEGER, INTENT(IN) :: nch        ! number of channels per node
     
    type(channel_t),      intent(inout) :: channel(nch) ! CRM data for multiple channels
    type(channel_vGCM_t), intent(inout) :: vGCM(nch)    ! vGCM data for multiple channels

!   Local 
    INTEGER num_seg,num_gcm,isize(4,nch),jsize(4,nch),k 

!   Specify model start time in fractional Julian day: determine rjday0 in timeinfo.F90
!   (JUNG_PLAN) Add it here.    

!   Initialize common parameters and vertical profiles in CONSTLD.f90
    CALL INIT_common (isRestart) 

!   Find the channel group-ID, each segment size & cube-face number
    do k = 1, nch
    
      CALL FIND_CHANNEL_INFO (vGCM(k)%vGCM_state,channel(k))  
        
      do num_seg = 1, 4
       isize(num_seg,k) = channel(k)%seg(num_seg)%nx_size
       jsize(num_seg,k) = channel(k)%seg(num_seg)%ny_size
      enddo     
    enddo

!   Allocate CRM channel segments                         
    CALL allocate_channel_data(channel, nch, isize, jsize)
            
!   Set up segment parameters: channel(:)%seg(:)
    DO k = 1, nch
    DO num_seg = 1, 4
    
    channel(k)%seg(num_seg)%mi1 = isize(num_seg,k)
    channel(k)%seg(num_seg)%mj1 = jsize(num_seg,k)
    
    channel(k)%seg(num_seg)%mim    = 1 - nhalo
    channel(k)%seg(num_seg)%mim_c  = 1 - nhalo_cal
    channel(k)%seg(num_seg)%mim_ce = channel(k)%seg(num_seg)%mim_c - 1
    channel(k)%seg(num_seg)%mim_a  = 1 - nhalo_adv
    
    channel(k)%seg(num_seg)%mip    = channel(k)%seg(num_seg)%mi1 + nhalo
    channel(k)%seg(num_seg)%mip_c  = channel(k)%seg(num_seg)%mi1 + nhalo_cal
    channel(k)%seg(num_seg)%mip_ce = channel(k)%seg(num_seg)%mip_c + 1
    channel(k)%seg(num_seg)%mip_a  = channel(k)%seg(num_seg)%mi1 + nhalo_adv
    
    channel(k)%seg(num_seg)%mjm    = 1 - nhalo
    channel(k)%seg(num_seg)%mjm_c  = 1 - nhalo_cal
    channel(k)%seg(num_seg)%mjm_ce = channel(k)%seg(num_seg)%mjm_c - 1
    channel(k)%seg(num_seg)%mjm_a  = 1 - nhalo_adv
    
    channel(k)%seg(num_seg)%mjp    = channel(k)%seg(num_seg)%mj1 + nhalo
    channel(k)%seg(num_seg)%mjp_c  = channel(k)%seg(num_seg)%mj1 + nhalo_cal
    channel(k)%seg(num_seg)%mjp_ce = channel(k)%seg(num_seg)%mjp_c + 1
    channel(k)%seg(num_seg)%mjp_a  = channel(k)%seg(num_seg)%mj1 + nhalo_adv 
    
    DO num_gcm = 1, nVGCM_seg 
     channel(k)%seg(num_seg)%lsta(num_gcm) = 1 + (num_gcm-1)*netsz
     channel(k)%seg(num_seg)%lend(num_gcm) = netsz + (num_gcm-1)*netsz
    ENDDO
    DO num_gcm = 1, nVGCM_seg 
     channel(k)%seg(num_seg)%lcen(num_gcm) = &
                D0_5*(channel(k)%seg(num_seg)%lsta(num_gcm) &
               +channel(k)%seg(num_seg)%lend(num_gcm))
    ENDDO
    
    ENDDO
    ENDDO

    do k = 1, nch
    
      ! Prepare mapping information.
      CALL MAPPING_CHANNEL (k, vGCM(k)%vGCM_state, channel(k))   
      
      ! Calculate background fields for initialization 
      CALL CAL_BG (vGCM(k)%vGCM_state, channel(k))
      
      ! Initialize the channel data.
      CALL INIT_channel (channel(k)) 
      
    enddo  
    
  end subroutine crm_init

!==============================================================================
  subroutine crm_run (itt, nch, channel, vGCM)
!==============================================================================
!   main subprogram handling the VVM channels (called by the main program)
!------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: itt     ! vGCM time step count
    INTEGER, INTENT(IN) :: nch     ! number of channels per node
    
    type(channel_t),      intent(inout) :: channel(nch)  ! CRM data for multiple channels
    type(channel_vGCM_t), intent(inout) :: vGCM(nch)     ! vGCM data for multiple channels

    ! Local
    LOGICAL :: FIRST = .TRUE.
    INTEGER K
    
#ifdef FIXTHIS
    INTEGER :: ITT_CRM_COUNT       ! vGCM time step in CRM's CRM's time-scale 
        
    ITT_CRM_COUNT = (ITT-1)*itt_max

    do k = 1, nch
    
      ! Calculate background fields 
      IF (.NOT. FIRST) CALL CAL_BG (vGCM(k)%vGCM_state, channel(k))
      ! Since BG-fields were already calculated in the initialization.

      ! Control the CRM integration
      CALL CRM_PREDICTION (ITT_CRM_COUNT, channel(k))

      ! Prepare the feedback (from channel data to vGCM_tend data)
      CALL CRM_OUT_TEND (channel(k), vGCM(k)%vGCM_tend)
      
      ! Prepare the vGCM-scale output (from channel data to vGCM_out data)
      CALL CRM_OUT (channel(k), vGCM(k)%vGCM_out)
         
    enddo
    
    FIRST = .FALSE.
    
#endif

  end subroutine crm_run

!==============================================================================
  subroutine crm_final (itt, channel, vGCM)
!==============================================================================
    integer,              intent(in)    :: itt           ! vGCM time step count
    type(channel_t),      intent(inout) :: channel(nch)  ! CRM data for multiple channels
    type(channel_vGCM_t), intent(inout) :: vGCM(nch)     ! vGCM data for multiple channels
    
!   What is this program for?       

  end subroutine crm_final

#ifdef FIXTHIS
!===================================================================================
    SUBROUTINE CRM_PREDICTION (ITT_CRM_COUNT, channel)
!===================================================================================
    integer,              intent(in)    :: ITT_CRM_COUNT
    type(channel_t),      intent(inout) :: channel       ! CRM data for one channel
    
    ! Local
    LOGICAL :: MAKE_AVG
    INTEGER :: ITT_CRM,ITT_NEW,N1,N2

    MAKE_AVG = .FALSE.

!   Initialize the variables related to the calculation of mean eddy effects (for time average)
    CALL INI_CRM_EFT (channel)

    IF (physics) THEN
!    Initialize the variables related to the calculation of mean sfc fluxes (for time average)
     CALL INI_SFC_EFT (channel)
    ENDIF

!   Calculate the relaxation timescale (CRM toward GCM)
    CALL CAL_RTIME (channel)

!-----------------------------------------------------------------------
      DO ITT_CRM = 1, itt_max         ! CRM subcycle time marching
!-----------------------------------------------------------------------
      ITT_NEW = ITT_CRM_COUNT + ITT_CRM

      IF(ITT_CRM.eq.1) THEN
        a = dt                 ! forward scheme to initialize the CRM
        b = D0_0               ! forward scheme to initialize the CRM
      ELSE
        a =   D3_0 / D2_0 * dt     ! AB second-order scheme
        b = - D1_0 / D2_0 * dt     ! AB second-order scheme
      ENDIF

      N1 = MOD ( ITT_CRM    , 2 ) + 1
      N2 = MOD ( ITT_CRM - 1, 2 ) + 1

#ifdef FIXTHIS
!     Defines the year and fractional Julian day
!     rjday, iyr in in timeinfo.F90
      CALL TIME_MANAGER(ITT_NEW)
#endif

      CALL INI_TENDENCY (channel)
      
!     Individual CRM integration
      CALL AB_3D (N1, N2, ITT_NEW, channel)

      IF (ITT_CRM.eq.itt_max) MAKE_AVG=.TRUE.

!     Determine eddy components
      CALL EDDY_PREPARE (channel)

!     Calculate the mean eddy effcts
      CALL CAL_EDDY_TRANS (itt_max,MAKE_AVG,channel)
      
      CALL CAL_TENDENCY (itt_max,MAKE_AVG,channel)

      IF (physics) THEN
!      Calculate the mean sfc fluxes
       CALL CAL_SFC_EFT (itt_max,MAKE_AVG,channel)
      ENDIF

!-----------------------------------
      ENDDO   ! CRM time marching
!-----------------------------------  

    END SUBROUTINE crm_prediction
    
!===================================================================================
    SUBROUTINE CRM_OUT_TEND (channel,vGCM_tend)
!===================================================================================
    type(channel_t),   intent(in)  :: channel       ! CRM data for one channel
    type(vGCM_tend_t), intent(out) :: vGCM_tend(4)  ! vGCM tendency
    
    ! Local
    INTEGER :: num_seg,mi1,mj1,nt

!*******************************  
    DO num_seg = 1, 4
!*******************************
    mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
    mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 

!--------------------------------------------------
   IF (mi1.GT.mj1) THEN
!-------------------------------------------------- x-channel segment  
     CALL TEND_X (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%TH_E1(2:nk2,:), &
                  channel%seg(num_seg)%TH_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dTH)

     CALL TEND_X (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QV_E1(2:nk2,:), &
                  channel%seg(num_seg)%QV_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQV)            

     if (physics) then
     CALL TEND_X (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QC_E1(2:nk2,:), &
                  channel%seg(num_seg)%QC_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQC)

     CALL TEND_X (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QI_E1(2:nk2,:), &
                  channel%seg(num_seg)%QI_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQI)            
                  
     CALL TEND_X (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QR_E1(2:nk2,:), &
                  channel%seg(num_seg)%QR_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQR)

     CALL TEND_X (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QS_E1(2:nk2,:), &
                  channel%seg(num_seg)%QS_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQS)            
                  
     CALL TEND_X (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QG_E1(2:nk2,:), &
                  channel%seg(num_seg)%QG_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQG)
     endif 
     
     DO NT = 1,ntracer 
     CALL TEND_X (channel%seg(num_seg)%NFACE,nk1,         & 
                  channel%seg(num_seg)%QT_E1(2:nk2,:,nt), &
                  channel%seg(num_seg)%QT_E2(2:nk2,:,nt), & 
                  vGCM_tend(num_seg)%dQT(:,:,nt))            
     ENDDO
                  
     CALL TEND_X (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%U_E1(2:nk2,:),  &
                  channel%seg(num_seg)%U_E2(2:nk2,:),  & 
                  vGCM_tend(num_seg)%dU)

     CALL TEND_X (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%V_E1(2:nk2,:),  &
                  channel%seg(num_seg)%V_E2(2:nk2,:),  & 
                  vGCM_tend(num_seg)%dV)                                                                                          
!--------------------------------------------------
   ELSE
!-------------------------------------------------- y-channel segment  
     CALL TEND_Y (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%TH_E1(2:nk2,:), &
                  channel%seg(num_seg)%TH_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dTH)
                  
     CALL TEND_Y (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QV_E1(2:nk2,:), &
                  channel%seg(num_seg)%QV_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQV)
     
     if (physics) then             
     CALL TEND_Y (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QC_E1(2:nk2,:), &
                  channel%seg(num_seg)%QC_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQC)                                    

     CALL TEND_Y (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QI_E1(2:nk2,:), &
                  channel%seg(num_seg)%QI_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQI)
                  
     CALL TEND_Y (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QR_E1(2:nk2,:), &
                  channel%seg(num_seg)%QR_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQR)
                  
     CALL TEND_Y (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QS_E1(2:nk2,:), &
                  channel%seg(num_seg)%QS_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQS)                                    

     CALL TEND_Y (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%QG_E1(2:nk2,:), &
                  channel%seg(num_seg)%QG_E2(2:nk2,:), & 
                  vGCM_tend(num_seg)%dQG)
     endif 
     
     DO NT = 1,ntracer              
     CALL TEND_Y (channel%seg(num_seg)%NFACE,nk1,         & 
                  channel%seg(num_seg)%QT_E1(2:nk2,:,nt), &
                  channel%seg(num_seg)%QT_E2(2:nk2,:,nt), & 
                  vGCM_tend(num_seg)%dQT(:,:,nt))
     ENDDO
                  
     CALL TEND_Y (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%U_E1(2:nk2,:),  &
                  channel%seg(num_seg)%U_E2(2:nk2,:),  & 
                  vGCM_tend(num_seg)%dU)                                    

     CALL TEND_Y (channel%seg(num_seg)%NFACE,nk1,      & 
                  channel%seg(num_seg)%V_E1(2:nk2,:),  &
                  channel%seg(num_seg)%V_E2(2:nk2,:),  & 
                  vGCM_tend(num_seg)%dV)         
!--------------------------------------------------
   ENDIF
!-------------------------------------------------- 
       
!*******************************  
   ENDDO  ! num_seg 
!*******************************

    CONTAINS
    
    SUBROUTINE TEND_X (NFACE,kdim,E1,E2,EOUT)
!   Combine the tendencies due to transport and diabatic effects
!   arrange the data for the vGCM_tend grid structure     
    
    INTEGER, INTENT(IN) :: NFACE     ! # of cube-face, to which the data belongs
    INTEGER, INTENT(IN) :: kdim      ! Vertical size of input data  
    
    REAL(kind=dbl_kind),DIMENSION(kdim,nVGCM_seg),INTENT(IN)  :: E1,E2
    REAL(kind=dbl_kind),DIMENSION(nVGCM_seg,kdim),INTENT(OUT) :: EOUT
    
    ! LOCAL
    INTEGER I,K,NPO 
    
    IF (NFACE .NE. 6) THEN
      DO I = 1,nVGCM_seg
       DO K = 1, kdim
        EOUT(I,K) = E1(K,I) + E2(K,I)
       ENDDO 
      ENDDO  
    ELSE
      DO I = 1,nVGCM_seg
       npo = nVGCM_seg - I + 1
       DO K = 1, kdim  
        EOUT(I,K) = E1(K,npo) + E2(K,npo)
       ENDDO 
      ENDDO  
    ENDIF
     
    END SUBROUTINE tend_x

    SUBROUTINE TEND_Y (NFACE,kdim,E1,E2,EOUT)
!   Combine the tendencies due to transport and diabatic effects
!   arrange the data for the vGCM_tend grid structure  
    
    INTEGER, INTENT(IN) :: NFACE     ! # of cube-face, to which the data belongs
    INTEGER, INTENT(IN) :: kdim      ! Vertical size of input data  
    
    REAL(kind=dbl_kind),DIMENSION(kdim,nVGCM_seg),INTENT(IN)  :: E1,E2
    REAL(kind=dbl_kind),DIMENSION(nVGCM_seg,kdim),INTENT(OUT) :: EOUT
    
    ! LOCAL
    INTEGER I,K,NPO 
    
    IF (NFACE.EQ.2 .OR. NFACE.EQ.3) THEN
      DO I = 1,nVGCM_seg
       npo = nVGCM_seg - I + 1
       DO K = 1, kdim
        EOUT(I,K) = E1(K,npo) + E2(K,npo)
       ENDDO 
      ENDDO  
    ELSE
      DO I = 1,nVGCM_seg
       DO K = 1, kdim  
        EOUT(I,K) = E1(K,I) + E2(K,I)
       ENDDO 
      ENDDO  
    ENDIF
     
    END SUBROUTINE tend_y   
          
    END SUBROUTINE crm_out_tend

!===================================================================================
    SUBROUTINE CRM_OUT (channel,vGCM_out)
    ! Currently, it has only 2D fields. 
!===================================================================================
    type(channel_t),   intent(in) :: channel      ! CRM data for one channel
    type(vGCM_out_t), intent(out) :: vGCM_out(4)  ! vGCM tendency
    
    ! Local
    INTEGER :: num_seg,mi1,mj1

!*******************************  
    DO num_seg = 1, 4
!*******************************
    mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
    mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 

!--------------------------------------------------
   IF (mi1.GT.mj1) THEN
!-------------------------------------------------- x-channel segment  

     if (physics) then
     CALL FIELD_X (channel%seg(num_seg)%NFACE,     & 
                   channel%seg(num_seg)%SPREC_E0,  &
                   vGCM_out(num_seg)%SPREC)

     CALL FIELD_X (channel%seg(num_seg)%NFACE,     & 
                   channel%seg(num_seg)%WTH_E0,    &
                   vGCM_out(num_seg)%WTH)            
                  
     CALL FIELD_X (channel%seg(num_seg)%NFACE,     & 
                  channel%seg(num_seg)%WQV_E0,     &
                  vGCM_out(num_seg)%WQV)

     CALL FIELD_X (channel%seg(num_seg)%NFACE,     & 
                  channel%seg(num_seg)%UW_E0,      &
                  vGCM_out(num_seg)%UW)            
                  
     CALL FIELD_X (channel%seg(num_seg)%NFACE,     & 
                  channel%seg(num_seg)%WV_E0,      &
                  vGCM_out(num_seg)%WV)
     endif 
                                                                      
!--------------------------------------------------
   ELSE
!-------------------------------------------------- y-channel segment  
     
     if (physics) then             
     CALL FIELD_Y (channel%seg(num_seg)%NFACE,     & 
                   channel%seg(num_seg)%SPREC_E0,  &
                   vGCM_out(num_seg)%SPREC)                                    

     CALL FIELD_Y (channel%seg(num_seg)%NFACE,     & 
                   channel%seg(num_seg)%WTH_E0,    &
                   vGCM_out(num_seg)%WTH)
                  
     CALL FIELD_Y (channel%seg(num_seg)%NFACE,     & 
                   channel%seg(num_seg)%WQV_E0,    &
                   vGCM_out(num_seg)%WQV)
                  
     CALL FIELD_Y (channel%seg(num_seg)%NFACE,     & 
                   channel%seg(num_seg)%UW_E0,     &
                   vGCM_out(num_seg)%UW)                                    

     CALL FIELD_Y (channel%seg(num_seg)%NFACE,     & 
                   channel%seg(num_seg)%WV_E0,     &
                   vGCM_out(num_seg)%WV)
     endif 
     
!--------------------------------------------------
   ENDIF
!-------------------------------------------------- 
       
!*******************************  
   ENDDO  ! num_seg 
!*******************************

    CONTAINS
    
    SUBROUTINE FIELD_X (NFACE,E1,EOUT)
!   arrange the data for the vGCM_out grid structure     
    
    INTEGER, INTENT(IN) :: NFACE     ! # of cube-face, to which the data belongs
    
    REAL(kind=dbl_kind),DIMENSION(nVGCM_seg),INTENT(IN)  :: E1
    REAL(kind=dbl_kind),DIMENSION(nVGCM_seg),INTENT(OUT) :: EOUT
    
    ! LOCAL
    INTEGER I,NPO 
    
    IF (NFACE .NE. 6) THEN
      DO I = 1,nVGCM_seg
        EOUT(I) = E1(I) 
      ENDDO  
    ELSE
      DO I = 1,nVGCM_seg
        npo = nVGCM_seg - I + 1
        EOUT(I) = E1(npo)
      ENDDO  
    ENDIF
     
    END SUBROUTINE field_x

    SUBROUTINE FIELD_Y (NFACE,E1,EOUT)
!   arrange the data for the vGCM_tend grid structure  
    
    INTEGER, INTENT(IN) :: NFACE     ! # of cube-face, to which the data belongs
    
    REAL(kind=dbl_kind),DIMENSION(nVGCM_seg),INTENT(IN)  :: E1
    REAL(kind=dbl_kind),DIMENSION(nVGCM_seg),INTENT(OUT) :: EOUT
    
    ! LOCAL
    INTEGER I,NPO 
    
    IF (NFACE.EQ.2 .OR. NFACE.EQ.3) THEN
      DO I = 1,nVGCM_seg
        npo = nVGCM_seg - I + 1
        EOUT(I) = E1(npo)
      ENDDO  
    ELSE
      DO I = 1,nVGCM_seg
        EOUT(I) = E1(I) 
      ENDDO  
    ENDIF
     
    END SUBROUTINE field_y   
          
    END SUBROUTINE crm_out
        
#endif

END MODULE q3d_interface_module
