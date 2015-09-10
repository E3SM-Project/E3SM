module prof_init_mod
  use grid_init_mod, only : dz, height

  use moist_init_mod, only : gg,cp,rg,rv, pi, tt0
  real, allocatable :: rho0(:), th0(:), th_e(:), ux_e(:),      &
       uy_e(:)  , qv_e(:), tm_e(:), qr_e(:), qc_e(:)

    !c convert the sounding to have winds along low-level winds:            
!    eps_wind = 1.e-5 
    ! surface to 1600 m                                      
!    kaver = 5 
  real :: eps_wind
  integer :: kaver
contains                                                                        
  SUBROUTINE prof_init(press,templs,zinls,vap,uuls,vvls,st,nz,npin) 
    implicit none
    integer, intent(in) :: nz, npin
    real, intent(in) :: press (npin), zinls(npin), &
         vap(npin), uuls(npin), vvls(npin),  templs(npin), st

    real :: temp(npin), uu(npin), vv(npin), zin(npin)

    real :: relhum (npin), zinkm (npin)
    real ::  coe2, deltz, th00, presnl, tavi, tempk, delt
    real :: tempkm, tvirt, tt00, wind_cos, wind_sin, deg, cap, capi
    real :: ymean, hmean, xmean, sum, sum1, sum2, sum3, rh00, cs, exs, pr00
    real :: thetme

    integer :: iisn, kk, l, km, k


    !c statement functions:                                                 
    !      alim01(fi)=amax1(0.,amin1(1.,fi))                                
    !      comb(tm,td,tu,ad,au)=                                            
    !     1  alim01((tm-td)/(tu-td))*au + alim01((tu-tm)/(tu-td))*ad        


    temp = templs
    zin = zinls
    uu= uuls
    vv=vvls

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
    DO k = 1, npin 
       !onvert from relative humidity to mixing ratio:                         
       thetme = (1.e3 / press (k) ) ** (rg / cp) 
       !onvert from temperature (deg C or K) into potential temperature        
       temp (k) = (temp (k) + tt0) * thetme 
    enddo
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      allocate(rho0(nz),th0(nz),th_e(nz),ux_e (nz),      &
           uy_e(nz),qv_e(nz),tm_e(nz),qr_e(nz),qc_e(nz))



    sum1 = 0. 
    sum2 = 0. 
    sum3 = 0. 


    DO k = 1, kaver 
       sum1 = sum1 + (press (k) - press (k + 1) ) 
       sum2 = sum2 + uu (k) * (press (k) - press (k + 1) ) 
       sum3 = sum3 + vv (k) * (press (k) - press (k + 1) ) 
    enddo
    xmean = sum2 / sum1 
    ymean = sum3 / sum1 
    hmean = sqrt (xmean**2 + ymean**2) + eps_wind 
    !c consider rotating here...                                            
    wind_sin = (ymean + eps_wind) / hmean 
    wind_cos = (xmean) / hmean 
    !c just E-W                                                             
    !       wind_sin=0.                                                     
    !       wind_cos=1.                                                     
 
    deg = atan (wind_sin / wind_cos) * 180. / pi
    
    PRINT * , 'low trop wind: mean,deg: ', hmean, deg , wind_sin, wind_cos,eps_wind

    DO k = 1, npin 
       uu (k) = uu (k) * wind_cos + vv (k) * wind_sin 
       vv (k) = - uu (k) * wind_sin + vv (k) * wind_cos 
    enddo
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

    !ompute approximated height of pressure levels:                         
    zin (1) = 0. 
    DO k = 2, npin 
       km = k - 1 
       tempk = temp (k) * (1.e3 / press (k) ) ** ( - rg / cp) * &
            (1. + .6e-3 * vap (k) )                                                 
       tempkm = temp (km) * (1.e3 / press (km) ) ** ( - rg / cp) *       &
            (1. + .6e-3 * vap (km) )                                          
       delt = tempk - tempkm 
       IF (delt.gt.1.e-4) then 
          tavi = alog (tempk / tempkm) / delt 
       ELSE 
          tavi = 1. / tempk 
       ENDIF
       deltz = - rg / (tavi * gg) * alog (press (k) / press (km) ) 
       zin (k) = zin (km) + deltz 
    enddo
    !      print*,'INPUT SOUNDING'                                          
    !      do k=1,npin                                                      
    !       print 921,k,zin(k),press(k),temp(k),vap(k)                      
    !921    format(1x,'k,z,p,theta,qv: ',i4,4e17.5)                         
    !      enddo                                                            

    !ompute environmental profiles from sounding assuming no topography:    
    !c surface data:                                                        
    iisn = 1 
    th_e (1) = temp (iisn) 
    tm_e (1) = th_e (1) * (1000. / press (iisn) ) ** ( - rg / cp) 
    qv_e (1) = vap (iisn) * 1.e-3 
    ux_e (1) = uu (1) 
    uy_e (1) = vv (1) 
    !      print*,'DETERMINED SURFACE DATA'                                 
    !c higher levels - interpolate:                                         
    !     print*,'INTERPOLATION TO HIGHER LEVELS'                           
    l = nz 
    DO k = 2, l 
!       zz (k) = height (k) 
       !       print*,'k,z= ',k,height(k)                                          
       DO kk = 2, npin 
          iisn = kk - 1 
          IF (zin (kk) .ge.height (k) ) goto 665 
       enddo
       !       print*,'INPUT SOUNDING DOES NOT GO HIGH ENOUGH. STOP.'          
       STOP 'SOUNDING' 
665    CONTINUE 
       !       print*,'iisn=',iisn                                             
       coe2 = (height (k) - zin (iisn) ) / (zin (iisn + 1) - zin (iisn) ) 
       th_e (k) = coe2 * temp (iisn + 1) + (1. - coe2) * temp (iisn) 
       qv_e (k) = (coe2 * vap (iisn + 1) + (1. - coe2) * vap (iisn) )    &
            * 1.e-3                                                           
       presnl = coe2 * press (iisn + 1) + (1. - coe2) * press (iisn) 
       tm_e (k) = th_e (k) * (1000. / presnl) ** ( - rg / cp) 
       ux_e (k) = coe2 * uu (iisn + 1) + (1. - coe2) * uu (iisn) 
       uy_e (k) = coe2 * vv (iisn + 1) + (1. - coe2) * vv (iisn) 
    enddo
    qr_e=0.
    qc_e=0.

    !ompute th00,tt00,pr00,rh00 and average stability for base state profile
    th00 = th_e (1) 
    tt00 = tm_e (1) 
    tvirt = tm_e (1) * (1. + .6 * qv_e (1) ) 
    rh00 = press (1) * 100. / (rg * tvirt) 
    pr00 = press (1) * 100. 
    sum = 0. 
    DO k = 2, nz - 1 
       sum = sum + (th_e (k + 1) - th_e (k - 1) ) / th_e (k) 
    enddo
!    st = sum / (float (nz - 2) * 2. * dz * gac (k) ) 

    !ccc overwrite with st calculated for deep atmosphere:                  
    ! taken from global model...                       

    !ccc overwrite with st calculated for deep atmosphere:                  
!    print*,'th00,tt00,pr00,rh00,st: ',th00,tt00,pr00,rh00,st, height

    !ompute reference state vertical profiles                               
    cap = rg / cp 
    capi = 1. / cap 
    cs = gg / (cp * tt00 * st) 
    DO k = 1, nz 
!       print *, k, height(k)
       exs = exp ( - st * height (k) ) 
       th0 (k) = th00 / exs 
       rho0 (k) = rh00 * exs * (1. - cs * (1. - exs) ) ** (capi - 1.) 
    enddo

     print*,'PROFILES'                                                 
      do k=1,nz                                                          
         print 200,height(k)/1.e3,th0(k),rho0(k),th_e(k),                    & 
         tm_e(k),qv_e(k)*1.e3,ux_e(k),uy_e(k)                         
      enddo                                                             

    RETURN  
 200     format(1x,'z,th0,rho0,the,tme,qve,ue:',f7.1,f7.2,2f7.1,e10.3,f6.1,2f6.1)    

  END SUBROUTINE prof_init
end module prof_init_mod
