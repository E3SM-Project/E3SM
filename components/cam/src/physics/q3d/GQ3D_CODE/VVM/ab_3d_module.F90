MODULE ab_3d_module
! Controls the integration of individual CRM (one timestep calculation)

USE shr_kind_mod,   only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t

USE CONSTLD,        only: notherm,buoy,physics,turbulence,iver_wind, &
                          rad2deg   ! debug
      
! Subroutines being called
USE vort_3d_module,       only: vort_3d,vort_3d_corec
USE wind_module,          only: wind_3d_v2,wind_3d_v1,wind_3d_v0
USE turb_3d_module,       only: turb_3d,turb_3d_therm,turb_3d_vort
USE buoyf_module,         only: buoyf_3d
USE q_chk_module,         only: q_chk_3d
USE update_thermo_module, only: update_thermodynamics,revise_thermodynamics
USE rcalc_module,         only: rcalc_3d

IMPLICIT NONE
PRIVATE

PUBLIC :: ab_3d

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE ab_3d
!-----------------------------------------------------------------------
!=======================================================================
   SUBROUTINE AB_3D (N1, N2, ITT, channel)
!=======================================================================
      type(channel_t), intent(inout) :: channel   ! channel data
      
      INTEGER, INTENT(IN) :: &
         ITT,       & ! time step count
         N1,        & ! AB forcing time index for previous timestep
         N2           ! AB forcing time index for current timestep
         
      LOGICAL :: DEBUG = .FALSE.      ! JUNG_DEBUG 
      CHARACTER (LEN=50) :: title    ! JUNG_DEBUG 

        IF (DEBUG) THEN
          write(unit=title,fmt='(a9,i3)') 'CRM_ITT =',ITT  
          CALL CHK_WRITE (title,channel) 
        ENDIF
!     Predict thermodynamic fields (advection, microphysics, radiation)
!     ITT is related to how often radiation is updated.  
      IF (.NOT.NOTHERM) CALL RCALC_3D (N1, N2, ITT, channel)
      
        IF (DEBUG) CALL CHK_WRITE ('A4_RCALC',channel,TH3D=.TRUE.,QV3D=.TRUE.)    
      
!     Calculate the nonlinear turbulence coefficients
!     ITT is related to how often surface flux is calculated.
      IF (TURBULENCE) CALL TURB_3D (ITT, channel)
      
        IF (DEBUG) CALL CHK_WRITE ('A4_TURB',channel)
        
      IF (.NOT.NOTHERM) THEN

!     Calculate the thermodynamic tendencies due to turbulence
      IF (TURBULENCE) CALL TURB_3D_THERM(channel)

        IF (DEBUG) CALL CHK_WRITE ('A4_TURB_THERM',channel,TH3D=.TRUE.,QV3D=.TRUE.)   
         
!     Re-update the thermodynamic variables
      CALL REVISE_THERMODYNAMICS(channel)

        IF (DEBUG) CALL CHK_WRITE ('A4_RE_THERM',channel,TH3D=.TRUE.,QV3D=.TRUE.)      

      IF (PHYSICS) THEN
!       Fill negative values & Saturation adjustment
        CALL Q_CHK_3D (channel)
        
        IF (DEBUG) CALL CHK_WRITE ('A4_Q_CHK',channel,TH3D=.TRUE.,QV3D=.TRUE.)
      ENDIF

      ENDIF

!     Calculate the buoyancy based on the predicted thermodynamic fields 
      IF (BUOY) CALL BUOYF_3D (channel)

!     Predict the vorticity components
      CALL VORT_3D  (N1, N2, channel) 
      
        IF (DEBUG) CALL CHK_WRITE ('A4_VORT',channel,Z3DX=.TRUE.,Z3DY=.TRUE.,Z3DZ=.TRUE.)

!     Calculate the vorticity tendencies due to turbulence
      IF (TURBULENCE) CALL TURB_3D_VORT(channel)
      
        IF (DEBUG) CALL CHK_WRITE ('A4_TURB_VORT',channel)

!     Re-update the vorticity components
      CALL VORT_3D_COREC(channel)   

        IF (DEBUG) CALL CHK_WRITE ('A4_VORT_COREC',channel,Z3DX=.TRUE.,Z3DY=.TRUE.,Z3DZ=.TRUE.)      

!     Diagnose the wind components
!      CALL WIND_3D  (.TRUE., .TRUE., channel)  ! UCHANGE=.T., WCHANGE=.T. 

      SELECT CASE (iver_wind)
      CASE(0)
        CALL WIND_3D_V0  (.TRUE., .TRUE., channel)  ! UCHANGE=.T., WCHANGE=.T. 
      CASE(1)
        CALL WIND_3D_V1  (.TRUE., .TRUE., channel)  ! UCHANGE=.T., WCHANGE=.T. 
      CASE(2)
        CALL WIND_3D_V2  (.TRUE., .TRUE., channel)  ! UCHANGE=.T., WCHANGE=.T. 
      END SELECT      

        IF (DEBUG) CALL CHK_WRITE ('A4_WIND',channel,W3D=.TRUE.,U3DX=.TRUE.,U3DY=.TRUE.)   

   END SUBROUTINE ab_3d

!=======================================================================      
      SUBROUTINE CHK_WRITE (title,channel, &
                            TH3D,QV3D,W3D,U3DX,U3DY,Z3DX,Z3DY,Z3DZ,PSI,CHI)
!=======================================================================
!     NUM_HALO can be different among the variables,
!     but, the horizontal array size of the variables is fixed, (mim:mip,mjm:mjp)
     
      CHARACTER (LEN=*), INTENT(IN)  :: title
      type(channel_t), INTENT(INOUT) :: channel   ! Channel data
      
      LOGICAL, INTENT(IN), OPTIONAL :: TH3D,QV3D,W3D,U3DX,U3DY, &
                                       Z3DX,Z3DY,Z3DZ,PSI,CHI 
              
      ! Local        
      INTEGER :: num_seg,mi1,mj1,ifortw 
      INTEGER :: ilo1,ilo2,ilo3,klo
      
      ifortw = channel%num_chn
      if (ifortw == 5) ifortw=55
      if (ifortw == 6) ifortw=65

      write(ifortw,*) TRIM(title)
       
!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      
      klo  = 2
      ilo1 = 1 
      ilo2 = 187
      ilo3 = 375

!---------------------------------
! TH3D
!---------------------------------
      IF (PRESENT(TH3D)) THEN
        if (TH3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'TH', &
                             channel%seg(num_seg)%TH3D(ilo1,1,klo), &
                             channel%seg(num_seg)%TH3D(ilo2,1,klo), &
                             channel%seg(num_seg)%TH3D(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'TH', &
                             channel%seg(num_seg)%TH3D(1,ilo1,klo), &
                             channel%seg(num_seg)%TH3D(1,ilo2,klo), &
                             channel%seg(num_seg)%TH3D(1,ilo3,klo)
          
        ENDIF                      
        endif                        
      ENDIF
!---------------------------------
! QV3D
!---------------------------------
      IF (PRESENT(QV3D)) THEN
        if (QV3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QV', &
                             channel%seg(num_seg)%QV3D(ilo1,1,klo), &
                             channel%seg(num_seg)%QV3D(ilo2,1,klo), &
                             channel%seg(num_seg)%QV3D(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QV', &
                             channel%seg(num_seg)%QV3D(1,ilo1,klo), &
                             channel%seg(num_seg)%QV3D(1,ilo2,klo), &
                             channel%seg(num_seg)%QV3D(1,ilo3,klo)
          
        ENDIF                      
        endif                        
      ENDIF
      
!---------------------------------
! W3D
!---------------------------------
     IF (PRESENT(W3D)) THEN
        if (W3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'W ', &
                             channel%seg(num_seg)%W3D(ilo1,1,klo), &
                             channel%seg(num_seg)%W3D(ilo2,1,klo), &
                             channel%seg(num_seg)%W3D(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'W ', &
                             channel%seg(num_seg)%W3D(1,ilo1,klo), &
                             channel%seg(num_seg)%W3D(1,ilo2,klo), &
                             channel%seg(num_seg)%W3D(1,ilo3,klo)
          
        ENDIF                      
        endif                        
      ENDIF
      
!---------------------------------
! U3DX
!---------------------------------
     IF (PRESENT(U3DX)) THEN
        if (U3DX) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'U ', &
                             channel%seg(num_seg)%U3DX(ilo1,1,klo), &
                             channel%seg(num_seg)%U3DX(ilo2,1,klo), &
                             channel%seg(num_seg)%U3DX(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'U ', &
                             channel%seg(num_seg)%U3DX(1,ilo1,klo), &
                             channel%seg(num_seg)%U3DX(1,ilo2,klo), &
                             channel%seg(num_seg)%U3DX(1,ilo3,klo)
          
        ENDIF                      
        endif                        
      ENDIF      

!---------------------------------
! U3DY
!---------------------------------
     IF (PRESENT(U3DY)) THEN
        if (U3DY) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'V ', &
                             channel%seg(num_seg)%U3DY(ilo1,1,klo), &
                             channel%seg(num_seg)%U3DY(ilo2,1,klo), &
                             channel%seg(num_seg)%U3DY(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'V ', &
                             channel%seg(num_seg)%U3DY(1,ilo1,klo), &
                             channel%seg(num_seg)%U3DY(1,ilo2,klo), &
                             channel%seg(num_seg)%U3DY(1,ilo3,klo)
          
        ENDIF                      
        endif                        
      ENDIF

!---------------------------------
! Z3DX
!---------------------------------
     IF (PRESENT(Z3DX)) THEN
        if (Z3DX) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'ZX', &
                             channel%seg(num_seg)%Z3DX(ilo1,1,klo), &
                             channel%seg(num_seg)%Z3DX(ilo2,1,klo), &
                             channel%seg(num_seg)%Z3DX(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'ZX', &
                             channel%seg(num_seg)%Z3DX(1,ilo1,klo), &
                             channel%seg(num_seg)%Z3DX(1,ilo2,klo), &
                             channel%seg(num_seg)%Z3DX(1,ilo3,klo)
          
        ENDIF                      
        endif                        
      ENDIF          
      
!---------------------------------
! Z3DY
!---------------------------------
     IF (PRESENT(Z3DY)) THEN
        if (Z3DY) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'ZY', &
                             channel%seg(num_seg)%Z3DY(ilo1,1,klo), &
                             channel%seg(num_seg)%Z3DY(ilo2,1,klo), &
                             channel%seg(num_seg)%Z3DY(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'ZY', &
                             channel%seg(num_seg)%Z3DY(1,ilo1,klo), &
                             channel%seg(num_seg)%Z3DY(1,ilo2,klo), &
                             channel%seg(num_seg)%Z3DY(1,ilo3,klo)
          
        ENDIF                      
        endif                        
      ENDIF      
      
!---------------------------------
! Z3DZ
!---------------------------------
     IF (PRESENT(Z3DZ)) THEN
        if (Z3DZ) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'ZZ', &
                             channel%seg(num_seg)%Z3DZ(ilo1,1,klo), &
                             channel%seg(num_seg)%Z3DZ(ilo2,1,klo), &
                             channel%seg(num_seg)%Z3DZ(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'ZZ', &
                             channel%seg(num_seg)%Z3DZ(1,ilo1,klo), &
                             channel%seg(num_seg)%Z3DZ(1,ilo2,klo), &
                             channel%seg(num_seg)%Z3DZ(1,ilo3,klo)
          
        ENDIF                      
        endif                        
      ENDIF      
      
!---------------------------------
! PSI
!---------------------------------
     IF (PRESENT(PSI)) THEN
        if (PSI) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'PS',   &
                             channel%seg(num_seg)%PSI(ilo1,1), &
                             channel%seg(num_seg)%PSI(ilo2,1), &
                             channel%seg(num_seg)%PSI(ilo3,1)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'PS',   &
                             channel%seg(num_seg)%PSI(1,ilo1), &
                             channel%seg(num_seg)%PSI(1,ilo2), &
                             channel%seg(num_seg)%PSI(1,ilo3)
          
        ENDIF                      
        endif                        
      ENDIF      
      
!---------------------------------
! CHI
!---------------------------------
     IF (PRESENT(CHI)) THEN
        if (CHI) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'CH',   &
                             channel%seg(num_seg)%CHI(ilo1,1), &
                             channel%seg(num_seg)%CHI(ilo2,1), &
                             channel%seg(num_seg)%CHI(ilo3,1)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'CH',   &
                             channel%seg(num_seg)%CHI(1,ilo1), &
                             channel%seg(num_seg)%CHI(1,ilo2), &
                             channel%seg(num_seg)%CHI(1,ilo3)
          
        ENDIF                      
        endif                        
      ENDIF                                      

!**********************************************************************   
      ENDDO   
!**********************************************************************

      END SUBROUTINE chk_write    

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
   
END MODULE ab_3d_module
