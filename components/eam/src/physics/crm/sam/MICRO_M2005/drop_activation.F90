module drop_activation
#ifdef MODAL_AERO
!----------------------------------------------------------------------------------------------------
!
! Purposes: calcualte dropelt number concentration activated from aerosol particle, used 
!           in Morrison's two-moment microphysics in SAM. It treats multimode aerosol population, 
!           and aerosol fields are taken from the modal aerosol treatment in CAM. 
!
! Method: This module is adopted from the module of ndrop used in CAM, originally writted by 
!         Steven Ghan.
! 
! Revision history: 
! July, 2009:  adopted from the module of ndrop used in CAM. 
!
!---------------------------------------------------------------------------------------------------- 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use modal_aero_data, only: ntot_amode
   use params, only: crm_rknd

   implicit none
   private
   save

   public :: drop_activation_init,  drop_activation_Ghan

   real(r8) :: npv(ntot_amode) ! number per volume concentration
   real(r8) :: alogsig(ntot_amode) ! natl log of geometric standard dev of aerosol
   real(r8) :: exp45logsig(ntot_amode)
   real(r8) :: argfactor(ntot_amode)
   real(r8) :: f1(ntot_amode),f2(ntot_amode)  ! abdul-razzak functions of width
   real(r8) :: t0 ! reference temperature
   real(r8) :: aten
   real(r8) :: surften       ! surface tension of water w/respect to air (N/m)
   real(r8) :: alogten,alog2,alog3,alogaten
   real(r8) :: third, twothird, sixth, zero
   real(r8) :: sq2, sqpi, pi

contains
!----------------------------------------------------------------------------------

!==================================================================================
subroutine drop_activation_init
!------------------------------------------------------------------------
! Initialize constants, and prescribed parameters. 
!-----------------------------------------------------------------------
  use modal_aero_data
  use physconst, only: rhoh2o, mwh2o, r_universal
  implicit none

  integer l,m
  real(r8) arg

!  mathematical constants

  zero=0._r8
  third=1./3._r8
  twothird=2.*third
  sixth=1./6._r8
  sq2=sqrt(2._r8)
  pi=4._r8*atan(1.0_r8)
  sqpi=sqrt(pi)

  t0=273.
  surften=0.076_r8
  aten=2.*mwh2o*surften/(r_universal*t0*rhoh2o)
  alogaten=log(aten)
  alog2=log(2._r8)
  alog3=log(3._r8)

  do m=1,ntot_amode
!    use only if width of size distribution is prescribed
     alogsig(m)=log(sigmag_amode(m))
     exp45logsig(m)=exp(4.5*alogsig(m)*alogsig(m))
     argfactor(m)=2./(3.*sqrt(2.)*alogsig(m))
     f1(m)=0.5*exp(2.5*alogsig(m)*alogsig(m))
     f2(m)=1.+0.25*alogsig(m)
  end do

  return
end subroutine drop_activation_init
!-------------------------------------------------------------------------------------------------------

!=======================================================================================================
subroutine drop_activation_Ghan(ncrms,icrm,wnuc4, tair4, rhoair4,  &
                          ndrop4, ines, smaxinout4, k)
!-------------------------------------------------------------------------------------------------------
!
!  Purpose and method:  calculates number, surface, and mass fraction of aerosols activated as CCN
!      calculates flux of cloud droplets, surface area, and aerosol mass into cloud
!      assumes an internal mixture within each of up to pmode multiple aerosol modes
!      a gaussiam spectrum of updrafts can be treated.

!      mks units

!      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
!      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.
!
! Revision history: 
! 2009-07-17:  Originally written by Gteven Ghan, and adopted by Minghuai Wang. 
!
!------------------------------------------------------------------------------------------------------------

      use physconst, only: rair, epsilo, cpair, rh2o, latvap, gravit,   &
                                 rhoh2o, mwh2o, r_universal
      use wv_saturation, only: estblf
      use physconst, only: epsqs => epsilo
      use error_function, only: erf
      use modal_aero_data
      use vars,   only: naer, vaer, hgaer

      implicit none


!   Input
      integer, intent(in) :: ncrms,icrm
      real(crm_rknd), intent (in)  ::    wnuc4          ! updraft velocity (m/s)
      real(crm_rknd), intent (in)  ::    tair4          ! air temperature (K)
      real(crm_rknd), intent (in)  ::    rhoair4        ! air density (kg/m3)
      integer, intent(in)  ::  ines           !  whether non-equillium saturation is used (ines=1: used).  
      real(crm_rknd), intent (inout)  ::    smaxinout4      ! For ines=1,  it is non-equlibrium saturation ratio (input)
                                                 ! for ines=0,  it is smax calculted from the activation parameterizaiton (output).
      integer, intent(in)  ::  k              ! the index of vertical levels. 

!   Output 
      real(crm_rknd), intent (out) ::   ndrop4          ! activated droplet number concentration  


!  Local
      real(r8)  ::    wnuc          ! updraft velocity (m/s)
      real(r8)  ::    tair          ! air temperature (K)
      real(r8)  ::    rhoair        ! air density (kg/m3)
      real(r8) na(ntot_amode)           ! aerosol number concentration (/m3)
      integer nmode      ! number of aerosol modes
      real(r8) volume(ntot_amode)     ! aerosol volume concentration (m3/m3)
      real(r8) hygro(ntot_amode)  ! hygroscopicity of aerosol mode

      real(r8) fn(ntot_amode)      ! number fraction of aerosols activated
      real(r8) fm(ntot_amode)      ! mass fraction of aerosols activated
      real(r8) fluxn(ntot_amode)   ! flux of activated aerosol number fraction into cloud (cm/s)
      real(r8) fluxm(ntot_amode)   ! flux of activated aerosol mass fraction into cloud (cm/s)
      real(r8) flux_fullact   ! flux of activated aerosol fraction assuming 100% activation (cm/s)
                              !    rce-comment
                              !    used for consistency check -- this should match (ekd(k)*zs(k))
                              !    also, fluxm/flux_fullact gives fraction of aerosol mass flux
                              !       that is activated
!      local

      real(r8), parameter :: p0 = 1013.25e2_r8    ! reference pressure (Pa)
      real(r8) sign(ntot_amode)    ! geometric standard deviation of size distribution
      real(r8) pres ! pressure (Pa)
      real(r8) diff0 ! diffusivity (m2/s)
      real(r8) conduct0 ! thermal conductivity (Joule/m/sec/deg)
      real(r8) es ! saturation vapor pressure
      real(r8) qs ! water vapor saturation mixing ratio
      real(r8) dqsdt ! change in qs with temperature
      real(r8) dqsdp ! change in qs with pressure
      real(r8) g ! thermodynamic function (m2/s)
      real(r8) zeta(ntot_amode), eta(ntot_amode)
      real(r8) lnsmax ! ln(smax)
      real(r8) alpha
      real(r8) gamma
      real(r8) beta
      real(r8) sqrtg(ntot_amode)
      real(r8) :: amcube(ntot_amode) ! cube of dry mode radius (m)
      real(r8) :: lnsm(ntot_amode) ! ln(smcrit)
      real(r8) smc(ntot_amode) ! critical supersaturation for number mode radius
      real(r8) alw,sqrtalw
      real(r8) smax
      real(r8) x,arg
      real(r8) xmincoeff
      real(r8) z
      real(r8) etafactor1,etafactor2(ntot_amode),etafactor2max
      real(r8) wmaxf           ! maximum update velocity   [m/s] 
      real(crm_rknd)    ndrop_act       
      integer m,n
!      numerical integration parameters
      real(r8), parameter :: eps=0.3_r8,fmax=0.99_r8,sds=3._r8

      real(r8), parameter :: namin=1.e6_r8   ! minimum aerosol number concentration (/m3)

      wnuc = wnuc4
      tair = tair4
      rhoair = rhoair4

! Set aerosol fields
      na = naer(icrm,k, :) 
      volume = vaer(icrm,k, :)
      hygro = hgaer(icrm,k, :) 

      nmode = ntot_amode
      wmaxf = 10.0

      fn(:)=0._r8
      fm(:)=0._r8
      fluxn(:)=0._r8
      fluxm(:)=0._r8
      flux_fullact=0._r8
      ndrop4 = 0.
      ndrop_act = 0.

      if(nmode.eq.1.and.na(1).lt.1.e-20_r8)return

      pres=rair*rhoair*tair
      diff0=0.211e-4_r8*(p0/pres)*(tair/t0)**1.94
      conduct0=(5.69_r8+0.017_r8*(tair-t0))*4.186e2_r8*1.e-5_r8 ! convert to J/m/s/deg
      es = estblf(tair)
      qs = epsilo*es/(pres-(1.0_r8 - epsqs)*es)
      dqsdt=latvap/(rh2o*tair*tair)*qs
      alpha=gravit*(latvap/(cpair*rh2o*tair*tair)-1./(rair*tair))
      gamma=(1+latvap/cpair*dqsdt)/(rhoair*qs)
      etafactor2max=1.e10/(alpha*wmaxf)**1.5 ! this should make eta big if na is very small.

      do m=1,nmode
         if(volume(m).gt.1.e-39_r8.and.na(m).gt.1.e-39_r8)then
!            number mode radius (m)
!           write(6,*)'alogsig,volc,na=',alogsig(m),volc(m),na(m)
            amcube(m)=(3.*volume(m)/(4.*pi*exp45logsig(m)*na(m)))  ! only if variable size dist
!           growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
!           should depend on mean radius of mode to account for gas kinetic effects
!           see Fountoukis and Nenes, JGR2005 and Meskhidze et al., JGR2006
!           for approriate size to use for effective diffusivity.
            g=1._r8/(rhoh2o/(diff0*rhoair*qs)                                    &
              +latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair)-1._r8))
            sqrtg(m)=sqrt(g)
            beta=2._r8*pi*rhoh2o*g*gamma
            etafactor2(m)=1._r8/(na(m)*beta*sqrtg(m))
            if(hygro(m).gt.1.e-10)then
               smc(m)=2.*aten*sqrt(aten/(27.*hygro(m)*amcube(m))) ! only if variable size dist
            else
               smc(m)=100.
            endif
!           write(6,*)'sm,hygro,amcube=',smcrit(m),hygro(m),amcube(m)
         else
            g=1._r8/(rhoh2o/(diff0*rhoair*qs)                                    &
              +latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair)-1._r8))
            sqrtg(m)=sqrt(g)
            smc(m)=1._r8
            etafactor2(m)=etafactor2max ! this should make eta big if na is very small.
         endif
         lnsm(m)=log(smc(m)) ! only if variable size dist
!        write(6,'(a,i4,4g12.2)')'m,na,amcube,hygro,sm,lnsm=', &
!                   m,na(m),amcube(m),hygro(m),sm(m),lnsm(m)
      enddo

!        single updraft

         if(wnuc.gt.0._r8)then

            alw=alpha*wnuc
            sqrtalw=sqrt(alw)
            etafactor1=alw*sqrtalw

            do m=1,nmode
               eta(m)=etafactor1*etafactor2(m)
               zeta(m)=twothird*sqrtalw*aten/sqrtg(m)
            enddo

            call maxsat(zeta,eta,nmode,smc,smax)

            lnsmax=log(smax)
            xmincoeff=alogaten-twothird*(lnsmax-alog2)-alog3

            do m=1,nmode
!               modal
               x=twothird*(lnsm(m)-lnsmax)/(sq2*alogsig(m))
               fn(m)=0.5_r8*(1._r8-erf(x))
               arg=x-1.5_r8*sq2*alogsig(m)
               fm(m)=0.5_r8*(1._r8-erf(arg))
               if(wnuc.gt.0._r8)then
                  fluxn(m)=fn(m)*wnuc
                  fluxm(m)=fm(m)*wnuc
               endif
               ndrop_act = ndrop_act + fn(m) * na (m)
            enddo
            flux_fullact = wnuc

            if(ines.eq.0) then
              ndrop4 = ndrop_act
              smaxinout4 = smax
            else if(ines.eq.1) then
! for non-equlibrium ss
              smax = smaxinout4
              lnsmax=log(smax)
              xmincoeff=alogaten-twothird*(lnsmax-alog2)-alog3

              do m=1,nmode
!               modal
                x=twothird*(lnsm(m)-lnsmax)/(sq2*alogsig(m))
                fn(m)=0.5_r8*(1._r8-erf(x))
                arg=x-1.5_r8*sq2*alogsig(m)
                fm(m)=0.5_r8*(1._r8-erf(arg))
                if(wnuc.gt.0._r8)then
                   fluxn(m)=fn(m)*wnuc
                   fluxm(m)=fm(m)*wnuc
                endif
                ndrop4 = ndrop4 + fn(m) * na (m)
              enddo
              flux_fullact = wnuc
              ndrop4 = min(ndrop4, ndrop_act)
            end if
            
          endif

! sensitivity tests:
!          ndrop4 = max(ndrop4, 100.*1.0e6)  ! the minimum activated droplet number is 100 /cm3

      return
end subroutine drop_activation_Ghan
!----------------------------------------------------------------------------------------

!=======================================================================================
      subroutine maxsat(zeta,eta,nmode,smc,smax)

!      calculates maximum supersaturation for multiple
!      competing aerosol modes.

!      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
!      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

      implicit none

      integer nmode ! number of modes
      real(r8) smc(ntot_amode) ! critical supersaturation for number mode radius
      real(r8) zeta(ntot_amode), eta(ntot_amode)
      real(r8) smax ! maximum supersaturation
      integer m  ! mode index
      real(r8) sum, g1, g2, g1sqrt, g2sqrt

      do m=1,nmode
         if(zeta(m).gt.1.e5_r8*eta(m).or.smc(m)*smc(m).gt.1.e5_r8*eta(m))then
!            weak forcing. essentially none activated
            smax=1.e-20_r8
         else
!            significant activation of this mode. calc activation all modes.
            go to 1
         endif
      enddo

      return

  1   continue

      sum=0
      do m=1,nmode
         if(eta(m).gt.1.e-20_r8)then
            g1=zeta(m)/eta(m)
            g1sqrt=sqrt(g1)
            g1=g1sqrt*g1
            g2=smc(m)/sqrt(eta(m)+3._r8*zeta(m))
            g2sqrt=sqrt(g2)
            g2=g2sqrt*g2
            sum=sum+(f1(m)*g1+f2(m)*g2)/(smc(m)*smc(m))
         else
            sum=1.e20_r8
         endif
      enddo

      smax=1._r8/sqrt(sum)

      return

end subroutine maxsat
!--------------------------------------------------------------------------------------

#endif
end module drop_activation

