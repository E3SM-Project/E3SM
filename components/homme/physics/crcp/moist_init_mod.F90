module moist_init_mod
! JPE _PRIM identifies the homme dynamics
#ifdef _PRIM 
  use physical_constants, only : gg=>g, cp, rg=>rgas, rv=>rwater_vapor, pi=>DD_PI    
#else
  real, parameter :: gg=9.806 , cp=1005., rg=287., rv=461.
  real, parameter :: pi = 3.141592653589793238462643383279 
#endif
  real, parameter :: hlatv=2.53e6, hlats=2.84e6 
  real, parameter :: ar=5.2e2, br=3., cr=130., dr=0.5, er=0.8
  real, parameter :: alphr=1.0, betr=2.0, gamb1r=6.0
  real, parameter :: gambd1r=11.7, anor=1.e7, anos=1.e7, dconc=200. 
  real, parameter :: tup=268., tdn=253., tt0=273.16, ee0=611.
  real, parameter :: as=2.5e-2, bs=2.0, cs=4.0, ds=0.25
  real, parameter :: es=0.2, alphs=0.3, bets=3.0
  real, parameter :: gamb1s=2.0, gambd1s=2.56
  real :: ddisp
contains
  SUBROUTINE moist_init 
    implicit none
! Create parametrs for the model:                                         
                                                                        
!      COMMON / rain_p0 / ar, br, cr, dr, er, alphr, betr, gamb1r,       &
!      gambd1r, anor                                                     
!      COMMON / rain_p1 / dconc, ddisp 
!      COMMON / snow_p0 / as, bs, cs, ds, es, alphs, bets, gamb1s,       &
!      gambd1s, anos                                                     
!      COMMON / temp_p / tup, tdn 
!      COMMON / latent / hlatv, hlats 
!      COMMON / reference / tt0, ee0 
!                                                                        
!      COMMON / const / gg, cp, rg, rv 
                                                                        
!c mass, terminal velocity, diameter                                    
!      DATA ar, br, cr, dr / 5.2e2, 3., 130., 0.5 / 
!      DATA as, bs, cs, ds / 2.5e-2, 2., 4., 0.25 / 
!c collection ef., alpha, beta                                          
!      DATA er, alphr, betr / 0.8, 1., 2. / 
!      DATA es, alphs, bets / 0.2, .3, 3. / 
!c No                                                                   
!      DATA anor, anos / 2 * 1.e7 / 
!c latent heats:                                                        
!      DATA hlatv, hlats / 2.53e6, 2.84e6 / 
!c cloud droplet concentration (per cc)                                 
                        ! must be between 50 and 2000                   
!      DATA dconc / 200. / 
!!       data dconc /100./ ! must be between 50 and 2000                 
!!      data dconc /1000./ ! must be between 50 and 2000                 
!c limiting temperatures                                                
!      DATA tup, tdn / 268., 253. / 
!c gammas:                                                              
!      DATA gamb1r, gambd1r / 6.0, 11.7 / 
!      DATA gamb1s, gambd1s / 2.0, 2.56 / 
!c reference temperature and saturated vapor pressure:                  
!      DATA tt0, ee0 / 273.16, 611. / 
!      DATA gg, cp, rg, rv / 9.72, 1005., 287., 461. / 
                                                                        
!c check consistency                                                    
      IF (dconc.lt.50..or.dconc.gt.2000.) then 
         PRINT * , ' *** inconsistent droplet concentration. stop.' 
         STOP 'dconc' 
      ENDIF 
!c calculate relative dispersion for Berry's autoconversion:            

      ddisp = 0.146 - 5.964e-2 * alog (dconc / 2000.) 

      PRINT 2075, anor, anos, dconc, ddisp 
 2075 FORMAT(1x,' N0 in raindrop distr.: ',e15.3/                       &
     & 1x,' N0 in snowflake distr.: ',e15.3/                            &
     & 1x,' Berry parameters of cloud droplet spectrum:'/               &
     & 1x,' droplet conc., relative disp.: ',2f12.4)                    
                                                                        
      RETURN 
      END SUBROUTINE moist_init                     
    end module moist_init_mod
