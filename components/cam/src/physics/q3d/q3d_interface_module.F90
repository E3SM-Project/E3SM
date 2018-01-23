MODULE q3d_interface_module
! Subprograms control the VVM channels

USE shr_kind_mod,    only: dbl_kind => shr_kind_r8
USE vGCM_data_types, only: vGCM_state_t, vGCM_tend_t, channel_vGCM_t
USE vvm_data_types,  only: channel_t,allocate_channel_data

USE parmsld, only: nVGCM_seg,netsz,nhalo_cal,nhalo_adv,nhalo
USE constld, only: d0_0,d1_0,d2_0,d3_0,d0_5,itt_max,a,b,dt,physics      

! subprograms being called
USE q3d_eddy_module, only: eddy_prepare

#ifdef FIXTHIS
USE time_manager_module, only: time_manager
#endif

USE q3d_init_module, only: init_common,init_channel
USE q3d_mapping_module, only: find_channel_info,mapping_channel
USE q3d_bg_module, only: cal_bg
USE q3d_effect_module, only: ini_tendency,ini_crm_eft,cal_eddy_trans, &
                             cal_tendency,ini_sfc_eft,cal_sfc_eft
USE q3d_rtime_module, only: cal_rtime

IMPLICIT NONE
PRIVATE

public :: crm_init
public :: crm_run
public :: crm_final

CONTAINS
!==============================================================================
  subroutine crm_init (isRestart, nch, channel, vGCM)
!==============================================================================
    logical, intent(in) :: isRestart  ! .true. iff restart run
    integer, intent(in) :: nch        ! number of channels per node
     
    type(channel_t),      intent(inout) :: channel(nch) ! CRM data for multiple channels
    type(channel_vGCM_t), intent(inout) :: vGCM(nch)    ! vGCM data for multiple channels

! JH: Assume that channel(:) & vGCM(:) are already allocated before calling crm_init.  

!   Local 
    integer num_seg,num_gcm,isize(4,nch),jsize(4,nch),k 

!   Specify model start time in fractional Julian day: determine rjday0 in timeinfo.F90
!   (Q: Get it from DYCORE ?)

!   Initialize common parameters and vertical profiles
    CALL INIT_common 

!   Find channel group-ID, segment size, cube face number
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
      
      ! Initialize the channel data.
      CALL INIT_channel (channel(k)) 
      
    enddo  
    
  end subroutine crm_init

!==============================================================================
  subroutine crm_run (itt, nch, channel, vGCM)
!==============================================================================
!   main subprogram handling the VVM channels (called by the main program)
!------------------------------------------------------------------------

    integer, intent(in) :: itt     ! vGCM time step count
    integer, intent(in) :: nch     ! number of channels per node
    
    type(channel_t),      intent(inout) :: channel(nch)  ! CRM data for multiple channels
    type(channel_vGCM_t), intent(inout) :: vGCM(nch)     ! vGCM data for multiple channels

    integer k
    
#ifdef FIXTHIS
    INTEGER :: ITT_CRM_COUNT       ! vGCM time step in CRM's view 
        
!!XXgoldyXX: Change this to be Julian day and seconds
    ITT_CRM_COUNT = INT((ITT-1)*itt_max)

    do k = 1, nch
    
      ! Calculate background fields 
      CALL CAL_BG (vGCM(k)%vGCM_state, channel(k))

      ! Control the CRM integration
      CALL CRM_PREDICTION (ITT_CRM_COUNT, channel(k))

!     CALL FEEDBACK (channel(k), vGCM(k))
!     Prepare the feedback (from channel to vGCM)
         
    enddo
    
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

    MAKE_AVG=.FALSE.

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

#endif

END MODULE q3d_interface_module
