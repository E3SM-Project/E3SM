    SUBROUTINE RADIATION_RRTMG(ITT, channel)                               
!========================================================================
!  This is the interface code between the VVM and the RRTMG radiation calculation.
!
!  channel_width = 1 is assumed (e.g., j=1 is the prediction point for x-array): chn = 1
!  If channel_width is wider than 1, then additional do-loop should be added.
!
!  For column physics (microphysics & radiation), 
!  the whole channel array is internally formed and used for calculation. 
!  
!  INPUT: TH3D,QV3D,QC3D,QI3D and arguments
!         rjday0,rjday,iyr (module: timeinfo)
!
!  OUTPUT:
!      FTHRAD,   tendency of potential temperature due to radiation (K/s)
!      OLR,      outgoing long wave radiation (W/m**2)
!      FULWO,    upward longwave flux (W/m**2)
!      FDLWO,    downward longwave flux (W/m**2)
!      FUSWO,    upward shortwave flux (W/m**2)
!      FDSWO,    downward shortwave flux (W/m**2)
!      DTRADLW,  longwave heating rate (K/s)
!      DTRADSW,  shortwave heating rate (K/s)
!      WPLIQ,    liquid water path (g/m**2)
!      WPICE,    ice water path (g/m**2)
!      RELIQ,    effective radius of water (microns)
!      REICE,    effective radius of ice (microns)    
!      FULWTOA,  upward longwave flux at the top of the atmosphere (W/m**2)
!      FUSWTOA,  upward shortwave flux at the top of the atmosphere (W/m**2)
!      FDSWTOA   downward shortwave flux at the top of the atmosphere (W/m**2)
!------------------------------------------------------------------
      USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
      USE parmsld,        only: channel_seg_l
      USE constld,        only: pi,pibar,rad2deg
      USE timeinfo,       only: rjday0,rjday,iyr
      
      USE vvm_data_types, only: channel_t   
      
!     Radiation input parameters and internal variables
      USE rrtm_params, only: latitude,longitude
      USE rrtm_grid,   only: nz,nzm,day,day0,nstep,icycle,nrad,iyear,masterproc    
      USE rrtm_vars,   only: tabs,qv,qcl,qci,sstxy,pres,presi
      USE rad,         only: qrad, &
                             lwUp_3d,lwDown_3d,swUp_3d,swDown_3d, &
                             lwHeatingRate_3d,swHeatingRate_3d,   &
                             lwp_3d,iwp_3d,reliq_3d,reice_3d,     &
                             NetlwUpToa,NetswDownToa,NetswUpToa
                          
      IMPLICIT NONE
      
      ! Input arguments
      INTEGER, INTENT(IN) :: ITT    ! Current dynamic time step
      type(channel_t), intent(inout) :: channel   ! Channel data
 
      ! Local variables
      INTEGER :: mi1,mj1,num_seg,chp,chn,chl,K
      
      REAL (KIND=dbl_kind), PARAMETER :: &
         QVMIN = 1.E-06, &   ! Minimum value of qv used in radiation calculation
         QCMIN = 1.E-07, &   ! Minimum value of qc used in radiation calculation
         QIMIN = 1.E-08      ! Minimum value of qi used in radiation calculation
            
! Timing variables
      NSTEP  = ITT
      ICYCLE = 1
     
      DAY0 = REAL(RJDAY0) 
      DAY  = REAL(RJDAY)
      IYEAR = IYR
   
      CHN = 1      ! Assume: channel_width = 1 
      
!**********************************************************************************
      IF ( (itt .EQ. 1) .OR. mod(itt,nrad)==1 ) THEN
!**********************************************************************************
!!!        if(masterproc) PRINT*, 'RRTMG radiation is called'

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
       
         SSTXY(chl,chn) = channel%seg(num_seg)%TG(chp,chn)
         
         ! Latitude and longitude of column grid points [deg]
         longitude(chl,chn) = channel%seg(num_seg)%RLON_T(chp,chn) * rad2deg
         latitude (chl,chn) = channel%seg(num_seg)%RLAT_T(chp,chn) * rad2deg

         ! Define radiation input variables.  Omit subsurface and top layer 
         ! values (indices 1 and NZ+1, respectively).
         DO K = 2, NZ
          TABS(chl,chn,K-1) = channel%seg(num_seg)%TH3D(chp,chn,K)*PIBAR(K)
          QV  (chl,chn,K-1) = MAX(channel%seg(num_seg)%QV3D(chp,chn,K), QVMIN)
          QCL (chl,chn,K-1) = MAX(channel%seg(num_seg)%QC3D(chp,chn,K), QCMIN)
          QCI (chl,chn,K-1) = MAX(channel%seg(num_seg)%QI3D(chp,chn,K), QIMIN)
         ENDDO
         
        ENDDO  ! chp_loop    
!      ---------------------                     
       ELSE
!      ---------------------  ! y-array   
        DO chp = 1,mj1
         chl = channel_seg_l*(num_seg-1) + chp
       
         SSTXY(chl,chn) = channel%seg(num_seg)%TG(chn,chp)
         
         ! Latitude and longitude of column grid points [deg]
         longitude(chl,chn) = channel%seg(num_seg)%RLON_T(chn,chp) * rad2deg
         latitude (chl,chn) = channel%seg(num_seg)%RLAT_T(chn,chp) * rad2deg

         ! Define radiation input variables.  Omit subsurface and top layer 
         ! values (indices 1 and NZ+1, respectively).
         DO K = 2, NZ
          TABS(chl,chn,K-1) = channel%seg(num_seg)%TH3D(chn,chp,K)*PIBAR(K)
          QV  (chl,chn,K-1) = MAX(channel%seg(num_seg)%QV3D(chn,chp,K), QVMIN)
          QCL (chl,chn,K-1) = MAX(channel%seg(num_seg)%QC3D(chn,chp,K), QCMIN)
          QCI (chl,chn,K-1) = MAX(channel%seg(num_seg)%QI3D(chn,chp,K), QIMIN)
         ENDDO

       ENDDO  ! chp_loop        
!      ---------------------       
       ENDIF  
!      ---------------------              
!===============================   
       ENDDO  ! num_seg 
!===============================    
  
!----------------------------------
      ! Call the RRTMG radiation driver
      CALL rad_full()
!!!      if(masterproc) print *,'calling rad_full'      
!----------------------------------

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
        
        ! Outgoing longwave radiation and other TOA fluxes
        channel%seg(num_seg)%OLR(chp,chn)     = FULWTOA(chl,chn)
        channel%seg(num_seg)%FUSWTOA(chp,chn) = NetswUpToa(chl,chn)
        channel%seg(num_seg)%FDSWTOA(chp,chn) = NetswDownToa(chl,chn)
        channel%seg(num_seg)%FULWTOA(chp,chn) = NetlwUpToa(chl,chn)
     

        ! Calculate potential temperature tendency term
        DO K = 2, NZ
         channel%seg(num_seg)%FTHRAD(chp,chn,K)  = qrad(chl,chn,K-1) / PIBAR(K)
         channel%seg(num_seg)%FULWO(chp,chn,K)   = lwUp_3d(chl,chn,K-1)
         channel%seg(num_seg)%FDLWO(chp,chn,K)   = lwDown_3d(chl,chn,K-1)
         channel%seg(num_seg)%FUSWO(chp,chn,K)   = swUp_3d(chl,chn,K-1)
         channel%seg(num_seg)%FDSWO(chp,chn,K)   = swDown_3d(chl,chn,K-1)
         channel%seg(num_seg)%DTRADLW(chp,chn,K) = lwHeatingRate_3d(chl,chn,K-1)
         channel%seg(num_seg)%DTRADSW(chp,chn,K) = swHeatingRate_3d(chl,chn,K-1)
        ENDDO
     
        DO K = 1, NZM
         channel%seg(num_seg)%WPLIQ(chp,chn,K) = lwp_3d(chl,chn,K)
         channel%seg(num_seg)%WPICE(chp,chn,K) = iwp_3d(chl,chn,K)
         channel%seg(num_seg)%RELIQ(chp,chn,K) = reliq_3d(chl,chn,K)
         channel%seg(num_seg)%REICE(chp,chn,K) = reice_3d(chl,chn,K)
        ENDDO

       ENDDO  ! chp_loop
              
!      ---------------------                     
       ELSE
!      ---------------------  ! y-array   

       DO chp = 1,mj1
        chl = channel_seg_l*(num_seg-1) + chp
        
        ! Outgoing longwave radiation and other TOA fluxes
        channel%seg(num_seg)%OLR(chn,chp)     = FULWTOA(chl,chn)
        channel%seg(num_seg)%FUSWTOA(chn,chp) = NetswUpToa(chl,chn)
        channel%seg(num_seg)%FDSWTOA(chn,chp) = NetswDownToa(chl,chn)
        channel%seg(num_seg)%FULWTOA(chn,chp) = NetlwUpToa(chl,chn)
     

        ! Calculate potential temperature tendency term
        DO K = 2, NZ
         channel%seg(num_seg)%FTHRAD(chn,chp,K)  = qrad(chl,chn,K-1) / PIBAR(K)
         channel%seg(num_seg)%FULWO(chn,chp,K)   = lwUp_3d(chl,chn,K-1)
         channel%seg(num_seg)%FDLWO(chn,chp,K)   = lwDown_3d(chl,chn,K-1)
         channel%seg(num_seg)%FUSWO(chn,chp,K)   = swUp_3d(chl,chn,K-1)
         channel%seg(num_seg)%FDSWO(chn,chp,K)   = swDown_3d(chl,chn,K-1)
         channel%seg(num_seg)%DTRADLW(chn,chp,K) = lwHeatingRate_3d(chl,chn,K-1)
         channel%seg(num_seg)%DTRADSW(chn,chp,K) = swHeatingRate_3d(chl,chn,K-1)
        ENDDO
     
        DO K = 1, NZM
         channel%seg(num_seg)%WPLIQ(chn,chp,K) = lwp_3d(chl,chn,K)
         channel%seg(num_seg)%WPICE(chn,chp,K) = iwp_3d(chl,chn,K)
         channel%seg(num_seg)%RELIQ(chn,chp,K) = reliq_3d(chl,chn,K)
         channel%seg(num_seg)%REICE(chn,chp,K) = reice_3d(chl,chn,K)
        ENDDO

       ENDDO  ! chp_loop
       
!      ---------------------       
       ENDIF  
!      ---------------------              

!===============================   
       ENDDO  ! num_seg 
!===============================    


!**********************************************************************************
      ENDIF  ! End (nstep .EQ. 1) .OR. (nradsteps .GE. nrad)
!**********************************************************************************

      RETURN
      END SUBROUTINE RADIATION_RRTMG
