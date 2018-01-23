MODULE rcalc_module 

USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld,        only: ntracer,nk2,nk3
USE constld,        only: d1_0,physics

! Subroutines being called
USE advec_3d_module,      only: advec_3d
USE update_thermo_module, only: update_thermodynamics
! also call PHYSICS_INTERFACE (not a module)

IMPLICIT NONE
PRIVATE

PUBLIC :: rcalc_3d

CONTAINS

! Local Subroutines:
!------------------------
! SUBROUTINE rcalc_3d
!------------------------
!    D1_0  = 1.0_dbl_kind

!=======================================================================
   SUBROUTINE RCALC_3D (N1, N2, ITT, channel)
!=======================================================================
!  Call Advection and Physics (Microphysics & Radiation)

      INTEGER, INTENT(IN) :: &
         ITT,       & ! time step count (used for radiation and surface flux calculation)
         N1,        & ! AB forcing time index for previous timestep
         N2           ! AB forcing time index for current timestep
         
      type(channel_t), intent(inout) :: channel   ! channel data
      
      ! Local
      REAL (KIND=dbl_kind), PARAMETER :: TP_HEIGHT = 12000._dbl_kind
      
      integer mi1,mim,mim_a,mim_c,mip,mip_c,mj1,mjm,mjm_a,mjm_c,mjp,mjp_c    
      integer L,num_seg,nt 
                              
!     temporary use of channel%seg(num_seg)%TERM1 and channel%seg(num_seg)%Z3DX0
      
      L = N2

!===============================   
      DO num_seg = 1, 4
!=============================== 
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment
      
      mim    = channel%seg(num_seg)%mim    
      mim_c  = channel%seg(num_seg)%mim_c  
      mim_a  = channel%seg(num_seg)%mim_a  
    
      mip    = channel%seg(num_seg)%mip    
      mip_c  = channel%seg(num_seg)%mip_c  
    
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      
      mjm    = channel%seg(num_seg)%mjm    
      mjm_c  = channel%seg(num_seg)%mjm_c  
      mjm_a  = channel%seg(num_seg)%mjm_a  
    
      mjp    = channel%seg(num_seg)%mjp    
      mjp_c  = channel%seg(num_seg)%mjp_c         
      
      channel%seg(num_seg)%TERM1(:,:,:) = D1_0
      
      CALL ADVEC_3D (mi1,mim,mim_a,mim_c,mip,mip_c,       &
                     mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,       &
                     channel%seg(num_seg)%TERM1,          &
                     channel%seg(num_seg)%Z3DX0(1:mi1,1:mj1,:), &
                     channel%seg(num_seg)%U3DX,           &
                     channel%seg(num_seg)%U3DY,           &
                     channel%seg(num_seg)%W3D,            &
                     channel%seg(num_seg)%RG_T,           &
                     channel%seg(num_seg)%RG_U,           &
                     channel%seg(num_seg)%RG_V,           &
                     channel%seg(num_seg)%KLOWQ_IJ,       &
                     .FALSE.,DIV_CAL=.TRUE.)

      ! TH3D
      CALL ADVEC_3D (mi1,mim,mim_a,mim_c,mip,mip_c,       &
                     mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,       &
                     channel%seg(num_seg)%TH3D,           &
                     channel%seg(num_seg)%FTH3D(1,1,1,L), &
                     channel%seg(num_seg)%U3DX,           &
                     channel%seg(num_seg)%U3DY,           &
                     channel%seg(num_seg)%W3D,            &
                     channel%seg(num_seg)%RG_T,           &
                     channel%seg(num_seg)%RG_U,           &
                     channel%seg(num_seg)%RG_V,           &
                     channel%seg(num_seg)%KLOWQ_IJ,       &
                     .FALSE.,TPHGT=TP_HEIGHT,             &
                     TERM_ADD=channel%seg(num_seg)%Z3DX0(1:mi1,1:mj1,:))
      ! QV3D                                    
      CALL ADVEC_3D (mi1,mim,mim_a,mim_c,mip,mip_c,       &
                     mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,       &
                     channel%seg(num_seg)%QV3D,           &
                     channel%seg(num_seg)%FQV3D(1,1,1,L), &
                     channel%seg(num_seg)%U3DX,           &
                     channel%seg(num_seg)%U3DY,           &
                     channel%seg(num_seg)%W3D,            &
                     channel%seg(num_seg)%RG_T,           &
                     channel%seg(num_seg)%RG_U,           &
                     channel%seg(num_seg)%RG_V,           &
                     channel%seg(num_seg)%KLOWQ_IJ,       &
                     .FALSE.)

      IF (PHYSICS) THEN

!     Positive Definite = .TRUE. (Must be True since a vertical filling is used)            
      
      ! QC3D
      CALL ADVEC_3D (mi1,mim,mim_a,mim_c,mip,mip_c,       &
                     mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,       &
                     channel%seg(num_seg)%QC3D,           &
                     channel%seg(num_seg)%FQC3D(1,1,1,L), &
                     channel%seg(num_seg)%U3DX,           &
                     channel%seg(num_seg)%U3DY,           &
                     channel%seg(num_seg)%W3D,            &
                     channel%seg(num_seg)%RG_T,           &
                     channel%seg(num_seg)%RG_U,           &
                     channel%seg(num_seg)%RG_V,           &
                     channel%seg(num_seg)%KLOWQ_IJ,       &
                     .TRUE.)
      
      ! QI3D
      CALL ADVEC_3D (mi1,mim,mim_a,mim_c,mip,mip_c,       &
                     mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,       &
                     channel%seg(num_seg)%QI3D,           &
                     channel%seg(num_seg)%FQI3D(1,1,1,L), &
                     channel%seg(num_seg)%U3DX,           &
                     channel%seg(num_seg)%U3DY,           &
                     channel%seg(num_seg)%W3D,            &
                     channel%seg(num_seg)%RG_T,           &
                     channel%seg(num_seg)%RG_U,           &
                     channel%seg(num_seg)%RG_V,           &
                     channel%seg(num_seg)%KLOWQ_IJ,       &
                     .TRUE.)
      
      ! QR3D                                          
      CALL ADVEC_3D (mi1,mim,mim_a,mim_c,mip,mip_c,       &
                     mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,       &
                     channel%seg(num_seg)%QR3D,           &
                     channel%seg(num_seg)%FQR3D(1,1,1,L), &
                     channel%seg(num_seg)%U3DX,           &
                     channel%seg(num_seg)%U3DY,           &
                     channel%seg(num_seg)%W3D,            &
                     channel%seg(num_seg)%RG_T,           &
                     channel%seg(num_seg)%RG_U,           &
                     channel%seg(num_seg)%RG_V,           &
                     channel%seg(num_seg)%KLOWQ_IJ,       &
                     .TRUE.)
      
      ! QS3D               
      CALL ADVEC_3D (mi1,mim,mim_a,mim_c,mip,mip_c,       &
                     mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,       &
                     channel%seg(num_seg)%QS3D,           &
                     channel%seg(num_seg)%FQS3D(1,1,1,L), &
                     channel%seg(num_seg)%U3DX,           &
                     channel%seg(num_seg)%U3DY,           &
                     channel%seg(num_seg)%W3D,            &
                     channel%seg(num_seg)%RG_T,           &
                     channel%seg(num_seg)%RG_U,           &
                     channel%seg(num_seg)%RG_V,           &
                     channel%seg(num_seg)%KLOWQ_IJ,       &
                     .TRUE.)
      
      ! QG3D               
      CALL ADVEC_3D (mi1,mim,mim_a,mim_c,mip,mip_c,       &
                     mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,       &
                     channel%seg(num_seg)%QG3D,           &
                     channel%seg(num_seg)%FQG3D(1,1,1,L), &
                     channel%seg(num_seg)%U3DX,           &
                     channel%seg(num_seg)%U3DY,           &
                     channel%seg(num_seg)%W3D,            &
                     channel%seg(num_seg)%RG_T,           &
                     channel%seg(num_seg)%RG_U,           &
                     channel%seg(num_seg)%RG_V,           &
                     channel%seg(num_seg)%KLOWQ_IJ,       &
                     .TRUE.)  
                                                             
      ENDIF  ! PHYSICS

!     Positive Definite can be set to True or False (no filling is used)
      DO nt = 1,ntracer
       CALL ADVEC_3D (mi1,mim,mim_a,mim_c,mip,mip_c,          &
                      mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,          &
                      channel%seg(num_seg)%QT3D(:,:,:,nt),    &
                      channel%seg(num_seg)%FQT3D(1,1,1,L,nt), &
                      channel%seg(num_seg)%U3DX,              &
                      channel%seg(num_seg)%U3DY,              &
                      channel%seg(num_seg)%W3D,               &
                      channel%seg(num_seg)%RG_T,              &
                      channel%seg(num_seg)%RG_U,              &
                      channel%seg(num_seg)%RG_V,              &
                      channel%seg(num_seg)%KLOWQ_IJ,          &
                      .FALSE.)  
      ENDDO

!===============================
      ENDDO   ! num_seg 
!=============================== 

!     Physics calculations
      IF (PHYSICS)  CALL PHYSICS_INTERFACE(N2, ITT, channel)

!     Update thermodynamic variables
      CALL UPDATE_THERMODYNAMICS (N1, N2, channel)

   END SUBROUTINE rcalc_3d

END MODULE rcalc_module
