module icecld_fraction_module

! --------------------------------------------------------------
! Purpose:  This module is used to calculate ice cloud fraction
! Minghuai Wang (minghuai.wang@pnnl.gov), 2013-01
!--------------------------------------------------------------

implicit none

public  aist_single

private 

   ! -------------- !
   ! Set Parameters !
   ! -------------- !

   ! ----------------------------------------!
   ! Parameters for Ice Stratus              !
   ! The default value are taken from CAM5.2 !
   ! -------------------------- -------------!
   real*4, private, parameter :: rhmini       = 0.90       ! Minimum rh for ice cloud fraction > 0.
   real*4, private, parameter :: rhmaxi       = 1.00        ! rhi at which ice cloud fraction = 1.
   real*4, private, parameter :: qist_min     = 1.e-7      ! Minimum in-stratus ice IWC constraint [ kg/kg ]
   real*4, private, parameter :: qist_max     = 5.e-3      ! Maximum in-stratus ice IWC constraint [ kg/kg ]
   real*4, private            :: rhminl = 0.8875                      ! Critical RH for low-level  liquid stratus clouds
   real*4, private            :: rhminl_adj_land = 0.100             ! rhminl adjustment for snowfree land
   real*4, private            :: rhminh = 0.800               ! Critical RH for high-level liquid stratus clouds
   real*4, private            :: premit = 40000.0                      ! Top    height for mid-level liquid stratus fraction
   real*4, private            :: premib = 700.0                      ! Bottom height for mid-level liquid stratus fraction
   integer,  private            :: iceopt = 5                    ! option for ice cloud closure 
                                                                ! 1=wang & sassen 2=schiller (iciwc)  
                                                                ! 3=wood & field, 4=Wilson (based on smith)
                                                                ! 5=modified slingo (ssat & empyt cloud)        
   real*4, private            :: icecrit=0.93                   ! Critical RH for ice clouds in Wilson & Ballard closure ( smaller = more ice clouds )


contains

!===================================================================
subroutine icecld_diag(   )

end subroutine icecld_diag 
!------------------------------------------------------------------

!====================================================================
subroutine aist_single( qv, T, p, qi, landfrac, snowh, aist )

   ! --------------------------------------------------------- !
   ! Compute non-physical ice stratus fraction                 ! 
   ! 
   ! This subroutine is adopted from cldwat2m_macro.F90 in CAM5.2,
   ! and was originally used in CAM5.2 to diagnose ice cloud fraction   
   ! --------------------------------------------------------- !

   use params,     only: rair => rgas, rv 
!   use wv_saturation, only: vqsatd2_water_single
   use module_mp_GRAUPEL, only: POLYSVP

   implicit none
  
   real*4, intent(in)  :: qv              ! Grid-mean water vapor[kg/kg]
   real*4, intent(in)  :: T               ! Temperature
   real*4, intent(in)  :: p               ! Pressure [Pa]
   real*4, intent(in)  :: qi              ! Grid-mean ice water content [kg/kg]
   real*4, intent(in)  :: landfrac        ! Land fraction
   real*4, intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

   real*4, intent(out) :: aist            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   ! Local variables
   real*4 rhmin                           ! Critical RH
   real*4 rhwght

   real*4 a,b,c,as,bs,cs                  ! Fit parameters
   real*4 Kc                              ! Constant for ice cloud calc (wood & field)
   real*4 ttmp                            ! Limited temperature
   real*4 icicval                         ! Empirical IWC value [ kg/kg ]
   real*4 rho                             ! Local air density
   real*4 esl                             ! Liq sat vapor pressure
   real*4 esi                             ! Ice sat vapor pressure
   real*4 ncf,phi                         ! Wilson and Ballard parameters
   real*4 esat, qsat, dqsdT

   real*4 rhi                             ! grid box averaged relative humidity over ice
   real*4 minice                          ! minimum grid box avg ice for having a 'cloud'
   real*4 mincld                          ! minimum ice cloud fraction threshold
   real*4 icimr                           ! in cloud ice mixing ratio
 ! real*4 qist_min                        ! minimum in cloud ice mixing ratio
 ! real*4 qist_max                        ! maximum in cloud ice mixing ratio                
   real*4 rhdif                           ! working variable for slingo scheme
   
   real*4 EP_2 

   ! Statement functions
   logical land
   land = nint(landfrac) == 1

   ! --------- !
   ! Constants !
   ! --------- !

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87
     b = 0.569
     c = 0.002892
   ! Schiller parameters ( Option.2 )
     as = -68.4202
     bs = 0.983917
     cs = 2.81795
   ! Wood and Field parameters ( Option.3 )
     Kc = 75.
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.e-12 
     mincld = 1.e-4 
   ! qist_min = 1.e-7  
   ! qist_max = 5.e-3 

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

!     call vqsatd2_water_single(T,p,esat,qsat,dqsdT)
     EP_2 = rair/rv
     esl = min(0.99*p, polysvp(T,0))
     esi = min(0.99*p, polysvp(T,1))
     esi = min(esi, esl)
     qsat= EP_2*esl/(p-esl)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195. ,min(T,253. )) - 273.16 
             icicval = a + b * ttmp + c * ttmp**2. 
             rho = p/(rair*T)
             icicval = icicval * 1.e-6  / rho 
         else
             ttmp = max(190. ,min(T,273.16 ))
             icicval = 10.  **(as * bs**ttmp + cs)
             icicval = icicval * 1.e-6  * 18.  / 28.97 
         endif
         aist =  max(0. ,min(qi/icicval,1. )) 
     elseif( iceopt.eq.3 ) then
         aist = 1.  - exp(-Kc*qi/(qsat*(esi/esl)))
         aist = max(0. ,min(aist,1. ))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land .and. (snowh.le.0.000001 ) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land .and. (snowh.le.0.000001 ) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0 -rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0 -rhwght)
           ! endif
         endif
         ncf = qi/((1.  - icecrit)*qsat)
         if( ncf.le.0.  ) then 
             aist = 0. 
         elseif( ncf.gt.0.  .and. ncf.le.1. /6.  ) then 
             aist = 0.5 *(6.  * ncf)**(2. /3. )
         elseif( ncf.gt.1. /6.  .and. ncf.lt.1.  ) then
             phi = (acos(3. *(1. -ncf)/2. **(3. /2. ))+4. *3.1415927 )/3. 
             aist = (1.  - 4.  * cos(phi) * cos(phi))
         else
             aist = 1. 
         endif
             aist = max(0. ,min(aist,1. ))
     elseif (iceopt.eq.5) then 
! set rh ice cloud fraction
             rhi= (qv+qi)/qsat * (esl/esi)
             rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
             aist = min(1.0 , max(rhdif,0. )**2)

! limiter to remove empty cloud and ice with no cloud
! and set icecld fraction to mincld if ice exists

!             if (qi.lt.minice) then
!                aist=0. 
!             else
!                aist=max(mincld,aist)
!             endif

! enforce limits on icimr
!             if (qi.ge.minice) then
!                icimr=qi/aist

!minimum
!                if (icimr.lt.qist_min) then
!                   aist = max(0. ,min(1. ,qi/qist_min))
!                endif
!maximum
!                if (icimr.gt.qist_max) then
!                   aist = max(0. ,min(1. ,qi/qist_max))
!                endif

!             endif
     endif 

   ! 0.999  is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0. ,min(aist,0.999 ))

   return
   end subroutine aist_single

end module icecld_fraction_module
