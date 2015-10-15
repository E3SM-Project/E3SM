!ccc                                                                    
!ccc computes difference bwtween two crcp history files                 
!ccc using unit 16 and unit 17, reports l2, linf of a-b                 
!ccc for each time sample                                               
!ccc                                                                    
      PROGRAM crcp_norm 
      PARAMETER (ntime0 = 1 * 24 * 60 * 4) 
!      parameter(ntime0=1*24*10*4)                                      
      PARAMETER (ntime1 = ntime0 / (24 * 60) ) 
      PARAMETER (nx = 101, nz = 61) 
                                                                        
      DIMENSION theta0 (nx, nz) 
      DIMENSION qv0 (nx, nz) 
      DIMENSION qc0 (nx, nz), qr0 (nx, nz) 
      DIMENSION ux0 (nx, nz), uz0 (nx, nz) 
      DIMENSION heat_av0 (nz) 
      DIMENSION slat_av0 (nx) 
      DIMENSION ssen_av0 (nx) 
      DIMENSION sprec_av0 (nx) 
      DIMENSION radc_av0 (nz) 
                                                                        
      DIMENSION theta (nx, nz) 
      DIMENSION qv (nx, nz) 
      DIMENSION qc (nx, nz), qr (nx, nz) 
      DIMENSION ux (nx, nz), uz (nx, nz) 
      DIMENSION heat_av (nz) 
      DIMENSION slat_av (nx) 
      DIMENSION ssen_av (nx) 
      DIMENSION sprec_av (nx) 
      DIMENSION radc_av (nz) 
                                                                        
      REAL theta_l2, theta_linf 
      REAL qv_l2, qv_linf 
      REAL qc_l2, qc_linf 
      REAL qr_l2, qr_linf 
      REAL ux_l2, ux_linf 
      REAL uz_l2, uz_linf 
      REAL heat_av_l2, heat_av_linf 
      REAL slat_av_l2, slat_av_linf 
      REAL ssen_av_l2, ssen_av_linf 
      REAL sprec_av_l2, sprec_av_linf 
      REAL radc_av_l2, radc_av_linf 
                                                                        
      DO i = 0, ntime1 
                                                                        
      READ (16) time0, nx0, nz0 
      READ (16) ux0, uz0, theta0, qv0, qc0, qr0 
      READ (16) heat_av0, radc_av0, slat_av0, ssen_av0, sprec_av0 
                                                                        
      READ (17) time1, nx1, nz1 
      READ (17) ux, uz, theta, qv, qc, qr 
      READ (17) heat_av, radc_av, slat_av, ssen_av, sprec_av 
                                                                        
      IF (time0.ne.time1) then 
         WRITE ( * , * ) '**** time mismatch' 
         STOP 
      ENDIF 
                                                                        
      IF ( (nx0.ne.nx1) .or. (nz0.ne.nz1) ) then 
         WRITE ( * , * ) '**** grid mismatch' 
         STOP 
      ENDIF 

      WRITE ( * , * ) 'time=', time0 
                                                                        
      CALL error_2d (nx, nz, theta0, theta, theta_l2, theta_linf,'theta') 
      CALL error_2d (nx, nz, qv0, qv, qv_l2, qv_linf, 'qv') 
      CALL error_2d (nx, nz, qc0, qc, qc_l2, qc_linf, 'qc') 
      CALL error_2d (nx, nz, qr0, qr, qr_l2, qr_linf, 'qr') 
      CALL error_2d (nx, nz, ux0, ux, ux_l2, ux_linf, 'ux') 
      CALL error_2d (nx, nz, uz0, uz, uz_l2, uz_linf, 'uz') 

      CALL error_1d (nx, slat_av0, slat_av, slat_av_l2, slat_av_linf) 
      CALL error_1d (nx, ssen_av0, ssen_av, ssen_av_l2, ssen_av_linf) 
      CALL error_1d (nx, sprec_av0, sprec_av, sprec_av_l2,sprec_av_linf)
      CALL error_1d (nz, heat_av0, heat_av, heat_av_l2, heat_av_linf) 
      CALL error_1d (nz, radc_av0, radc_av, radc_av_l2, radc_av_linf) 
                                                                        
!      WRITE ( * , * ) 'theta error: ', theta_l2, theta_linf 
!      WRITE ( * ,  * ) 'qv    error: ', qv_l2, qv_linf 
!      WRITE ( * ,  * ) 'qc    error: ', qc_l2, qc_linf 
!      WRITE ( * ,  * ) 'qr    error: ', qr_l2, qr_linf 
!      WRITE ( * ,  * ) 'ux    error: ', ux_l2, ux_linf 
!      WRITE ( * ,  * ) 'uz    error: ', uz_l2, uz_linf 
      WRITE ( * ,  * ) 'slat  error: ', slat_av_l2, slat_av_linf 
      WRITE ( * ,  * ) 'ssen  error: ', ssen_av_l2, ssen_av_linf 
      WRITE ( * , * ) 'sprec error: ', sprec_av_l2, sprec_av_linf 
      WRITE ( * ,  * ) 'heat  error: ', heat_av_l2, heat_av_linf 
      WRITE ( * ,  * ) 'radc  error: ', radc_av_l2, radc_av_linf 
      WRITE ( *, * ) 
                                                                        
      enddo 
                                                                        
      END PROGRAM crcp_norm                         
                                                                        
      SUBROUTINE error_2d (nx, nz, f0, f1, err_l2, err_linf, name) 
      INTEGER nx, nz, i, k , imax, kmax
      character *(*) name
      DIMENSION f0 (nx, nz), f1 (nx, nz) 
      REAL err_l2, err_linf, d , dmax

      dmax=0.0
      err_linf = 0.0 
      err_l2 = 0.0 
      d0 = abs (f0 (1, 1) ) 
      DO k = 1, nz 
      DO i = 1, nx 
      d = abs (f1 (i, k) - f0 (i, k) ) 

      if(d>dmax) then
         dmax=d
         imax=i
         kmax=k
      endif

      err_linf = max (err_linf, d) 
      err_l2 = err_l2 + d * d 
      enddo 
      enddo 
      err_l2 = sqrt (err_l2) / float (nx * nz) 

      write(*,'(a6,a, 2e18.10,2i4)') name,' error:',err_l2,err_linf

      if(dmax>0.0) then
         print *, 'max at i=',imax,' k=',kmax,f0(imax,kmax),f1(imax,kmax)
      endif

      RETURN 
      END SUBROUTINE error_2d                       
                                                                        
      SUBROUTINE error_1d (n, f0, f1, err_l2, err_linf) 
      INTEGER n, i 
      DIMENSION f0 (n), f1 (n) 
      REAL err_l2, err_linf, d 
      err_linf = 0.0 
      err_l2 = 0.0 
      DO i = 1, n 
      d = abs (f1 (i) - f0 (i) ) 
      err_linf = max (err_linf, d) 
      err_l2 = err_l2 + d * d 
      enddo 
      err_l2 = sqrt (err_l2) / float (n) 
      RETURN 
      END SUBROUTINE error_1d                       
