module activate_drop_mam

!---------------------------------------------------------------------------------
!
! Routines for droplet activation by modal aerosols
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
#ifndef HAVE_ERF_INTRINSICS
use shr_spfn_mod,   only: erf => shr_spfn_erf
#endif
use wv_saturation,  only: qsat

use cam_abortutils,       only: endrun

implicit none
private       
save

public &
   actdrop_mam_init,  &
   actdrop_mam_calc

integer :: iulog           ! fortran unit to use for log messages

real(r8) :: pi
real(r8) :: rair
real(r8) :: rgas
real(r8) :: rh2o
real(r8) :: rhoh2o
real(r8) :: mwh2o
real(r8) :: epsilo      ! Ratio of h2o to dry air molecular weights 
real(r8) :: latvap
real(r8) :: cpair
real(r8) :: gravit

real(r8), parameter :: t0       = 273._r8
real(r8), parameter :: p0       = 1013.25e2_r8    ! reference pressure (Pa)
real(r8), parameter :: surften  = 0.076_r8
real(r8), parameter :: zero     = 0._r8
real(r8), parameter :: third    = 1._r8/3._r8
real(r8), parameter :: twothird = 2._r8*third
real(r8), parameter :: sixth    = 1._r8/6._r8

real(r8) :: sq2, alog2, alog3, sqpi
real(r8) :: aten, alogaten

integer :: nmode
real(r8), allocatable :: alogsig(:)     ! natl log of geometric standard dev of aerosol
real(r8), allocatable :: exp45logsig(:)
real(r8), allocatable :: f1(:)          ! abdul-razzak functions of width
real(r8), allocatable :: f2(:)          ! abdul-razzak functions of width

!===============================================================================
contains
!===============================================================================

subroutine actdrop_mam_init( &
   iulog_in, pi_in, rair_in, rgas_in, rh2o_in,          &
   rhoh2o_in, mwh2o_in, epsilo_in, latvap_in, cpair_in, &
   gravit_in, ntot_amode, sigmag_amode)

   integer,              intent(in)  :: iulog_in         ! fortran unit to use for log messages
   real(r8),             intent(in)  :: pi_in
   real(r8),             intent(in)  :: rair_in
   real(r8),             intent(in)  :: rgas_in
   real(r8),             intent(in)  :: rh2o_in
   real(r8),             intent(in)  :: rhoh2o_in
   real(r8),             intent(in)  :: mwh2o_in
   real(r8),             intent(in)  :: epsilo_in       ! Ratio of h2o to dry air molecular weights
   real(r8),             intent(in)  :: latvap_in
   real(r8),             intent(in)  :: cpair_in
   real(r8),             intent(in)  :: gravit_in
   integer,              intent(in)  :: ntot_amode
   real(r8),             intent(in)  :: sigmag_amode(ntot_amode)

   integer :: m

   character(len=*), parameter :: subname='actdrop_mam_init'
   !---------------------------------------------------------------------------------
    
   iulog  = iulog_in

   pi     = pi_in
   rair   = rair_in
   rgas   = rgas_in
   rh2o   = rh2o_in
   rhoh2o = rhoh2o_in
   mwh2o  = mwh2o_in
   epsilo = epsilo_in
   latvap = latvap_in
   cpair  = cpair_in
   gravit = gravit_in

   sq2    = sqrt(2._r8)
   alog2  = log(2._r8)
   alog3  = log(3._r8)
   sqpi   = sqrt(pi)
   aten   = 2._r8*mwh2o*surften/(rgas*t0*rhoh2o)
   alogaten = log(aten)

   allocate( &
      alogsig(ntot_amode),      &
      exp45logsig(ntot_amode),  &
      f1(ntot_amode),           &
      f2(ntot_amode) )

   do m = 1, ntot_amode
      ! use only if width of size distribution is prescribed
      alogsig(m)     = log(sigmag_amode(m))
      exp45logsig(m) = exp(4.5_r8*alogsig(m)*alogsig(m))
      f1(m)          = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m))
      f2(m)          = 1._r8 + 0.25_r8*alogsig(m)
   end do

end subroutine actdrop_mam_init

!---------------------------------------------------------------------------------

subroutine actdrop_mam_calc(        &
   wbar, sigw, wdiab, wminf, wmaxf, &
   tair, rhoair, na, nmode, volume, &
   hygro, in_cloud, smax_f, fn, fm, &
   fluxn, fluxm, flux_fullact)

   ! calculates number, surface, and mass fraction of aerosols activated as CCN
   ! calculates flux of cloud droplets, surface area, and aerosol mass into cloud
   ! assumes an internal mixture within each of up to nmode multiple aerosol modes
   ! a gaussiam spectrum of updrafts can be treated.
   !
   ! **** NOTE **** The modification to the in-cloud calculation of smax has only
   !                implemented in the single updraft branch of the sigw conditional
   !
   ! mks units
   !
   ! Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   ! 2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

   ! Input Arguments
   real(r8), intent(in) :: wbar       ! grid cell mean vertical velocity (m/s)
   real(r8), intent(in) :: sigw       ! subgrid standard deviation of vertical vel (m/s)
   real(r8), intent(in) :: wdiab      ! diabatic vertical velocity (0 if adiabatic)
   real(r8), intent(in) :: wminf      ! minimum updraft velocity for integration (m/s)
   real(r8), intent(in) :: wmaxf      ! maximum updraft velocity for integration (m/s)
   real(r8), intent(in) :: tair       ! air temperature (K)
   real(r8), intent(in) :: rhoair     ! air density (kg/m3)
   real(r8), intent(in) :: na(:)      ! (nmode) aerosol number concentration (/m3)
   integer,  intent(in) :: nmode      ! number of aerosol modes
   real(r8), intent(in) :: volume(:)  ! (nmode) aerosol volume concentration (m3/m3)
   real(r8), intent(in) :: hygro(:)   ! (nmode) hygroscopicity of aerosol mode
   logical,  intent(in) :: in_cloud   ! switch to modify calculations when above cloud base
   real(r8), intent(in) :: smax_f     ! droplet and rain size distr factor in the smax calculation
                                      ! used when in_cloud=.true.

   ! Output Arguments
   real(r8), intent(out) :: fn(:)      ! (nmode) number fraction of aerosols activated
   real(r8), intent(out) :: fm(:)      ! (nmode) mass fraction of aerosols activated
   real(r8), intent(out) :: fluxn(:)   ! (nmode) flux of activated aerosol number fraction into cloud (cm/s)
   real(r8), intent(out) :: fluxm(:)   ! (nmode) flux of activated aerosol mass fraction into cloud (cm/s)
   real(r8), intent(out) :: flux_fullact   ! flux of activated aerosol fraction assuming 100% activation (cm/s)

   !    rce-comment
   !    used for consistency check -- this should match (ekd(k)*zs(k))
   !    also, fluxm/flux_fullact gives fraction of aerosol mass flux
   !       that is activated

   !      local

   integer, parameter:: nx=200
   integer :: m, n

   real(r8), parameter :: eps=0.3_r8, fmax=0.99_r8, sds=3._r8

   real(r8) :: pres ! pressure (Pa)
   real(r8) :: diff0, conduct0

   real(r8) :: es ! saturation vapor pressure
   real(r8) :: qs ! water vapor saturation mixing ratio
   real(r8) :: dqsdt ! change in qs with temperature

   real(r8) :: alpha
   real(r8) :: gamma
   real(r8) :: etafactor2(nmode), etafactor2max

   real(r8) :: amcube(nmode) ! cube of dry mode radius (m)
   real(r8) :: grow, beta
   real(r8) :: sqrtg(nmode)
   real(r8) :: smc(nmode)    ! critical supersaturation for number mode radius
   real(r8) :: lnsm(nmode)   ! ln(smc)

   real(r8) :: wmin, wmax, w, dw, dwmax, dwmin, wnuc, dwnew
   real(r8) :: dfmin, dfmax
   real(r8) :: sumflxn(nmode)
   real(r8) :: sumflxm(nmode)
   real(r8) :: sumfn(nmode)
   real(r8) :: sumfm(nmode)
   real(r8) :: fnold(nmode)   ! number fraction activated
   real(r8) :: fmold(nmode)   ! mass fraction activated
   real(r8) :: sumflx_fullact
   real(r8) :: fold, wold, gold
   real(r8) :: alw, sqrtalw, etafactor1
   real(r8) :: zeta(nmode), eta(nmode)
   real(r8) :: smax
   real(r8) :: lnsmax ! ln(smax)
   real(r8) :: x, fnew
   real(r8) :: g, z, fnmin

   real(r8) :: arg
   real(r8) :: fnbar, fmbar, wb

   real(r8) :: z1, z2
   real(r8) :: integ,integf
   real(r8) :: wf1, wf2, zf1, zf2, gf1, gf2, gf

   character(len=*), parameter :: subname='actdrop_mam_calc'
   !---------------------------------------------------------------------------------

   fn(:)        = 0._r8
   fm(:)        = 0._r8
   fluxn(:)     = 0._r8
   fluxm(:)     = 0._r8
   flux_fullact = 0._r8

   if (nmode == 1 .and. na(1) < 1.e-20_r8) return

   if (sigw <= 1.e-5_r8 .and. wbar <= 0._r8) return

   pres     = rair*rhoair*tair
   diff0    = 0.211e-4_r8*(p0/pres)*(tair/t0)**1.94_r8
   conduct0 = (5.69_r8 + 0.017_r8*(tair-t0))*4.186e2_r8*1.e-5_r8 ! convert to J/m/s/deg
   call qsat(tair, pres, es, qs)
   dqsdt = latvap/(rh2o*tair*tair)*qs
   alpha = gravit*(latvap/(cpair*rh2o*tair*tair) - 1._r8/(rair*tair))
   gamma = (1.0_r8 + latvap/cpair*dqsdt)/(rhoair*qs)

   etafactor2max=1.e10_r8/(alpha*wmaxf)**1.5_r8 ! this should make eta big if na is very small.

   ! growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
   ! should depend on mean radius of mode to account for gas kinetic effects
   ! see Fountoukis and Nenes, JGR2005 and Meskhidze et al., JGR2006
   ! for approriate size to use for effective diffusivity.
   grow = 1._r8/(rhoh2o/(diff0*rhoair*qs)  &
          + latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair) - 1._r8))

   do m = 1, nmode

      sqrtg(m) = sqrt(grow)

      if (volume(m) > 1.e-39_r8 .and. na(m) > 1.e-39_r8) then

         ! number mode radius (m)
         ! write(iulog,*)'alogsig,volc,na=',alogsig(m),volc(m),na(m)
         amcube(m) = (3._r8*volume(m)/(4._r8*pi*exp45logsig(m)*na(m)))  ! only if variable size dist

         beta = 2._r8*pi*rhoh2o*grow*gamma
         etafactor2(m) = 1._r8/(na(m)*beta*sqrtg(m))

         if (hygro(m) > 1.e-10_r8) then
            smc(m) = 2._r8*aten*sqrt(aten/(27._r8*hygro(m)*amcube(m))) ! only if variable size dist
         else
            smc(m) = 100._r8
         endif
         ! write(iulog,*)'hygro,amcube=',hygro(m),amcube(m)
      else
         etafactor2(m) = etafactor2max ! this should make eta big if na is very small.
         smc(m) = 1._r8
      endif

      lnsm(m) = log(smc(m))

      ! write(iulog,'(a,i4,4g12.2)')'m,na,amcube,hygro,sm,lnsm=', &
      !                   m,na(m),amcube(m),hygro(m),sm(m),lnsm(m)
   end do

   if (sigw > 1.e-5_r8) then ! spectrum of updrafts

      wmax  = min(wmaxf, wbar+sds*sigw)
      wmin  = max(wminf, -wdiab)
      wmin  = max(wmin, wbar-sds*sigw)
      w     = wmin
      dwmax = eps*sigw
      dw    = dwmax
      dfmax = 0.2_r8
      dfmin = 0.1_r8

      ! *** return ***
      if (wmax <= w) return

      sumflxn(:)     = 0._r8
      sumfn(:)       = 0._r8
      fnold(:)       = 0._r8
      sumflxm(:)     = 0._r8
      sumfm(:)       = 0._r8
      fmold(:)       = 0._r8
      sumflx_fullact = 0._r8

      fold = 0._r8
      wold = 0._r8
      gold = 0._r8

      dwmin = min(dwmax, 0.01_r8)

      do n = 1, nx

100      wnuc       = w + wdiab
         ! write(iulog,*)'wnuc=',wnuc
         alw        = alpha*wnuc
         sqrtalw    = sqrt(alw)
         etafactor1 = alw*sqrtalw

         do m = 1, nmode
            eta(m)  = etafactor1*etafactor2(m)
            zeta(m) = twothird*sqrtalw*aten/sqrtg(m)
         enddo

         call maxsat(zeta, eta, nmode, smc, smax)
         ! write(iulog,*)'w,smax=',w,smax

         lnsmax = log(smax)
         x      = twothird*(lnsm(nmode) - lnsmax)/(sq2*alogsig(nmode))
         fnew   = 0.5_r8*(1._r8 - erf(x))

         dwnew = dw
         if (fnew-fold > dfmax .and. n > 1) then
            ! reduce updraft increment for greater accuracy in integration
            if (dw .gt. 1.01_r8*dwmin) then
               dw = 0.7_r8*dw
               dw = max(dw, dwmin)
               w  = wold + dw
               go to 100
            else
               dwnew = dwmin
            end if
         end if

         if (fnew-fold < dfmin) then
            ! increase updraft increment to accelerate integration
            dwnew = min(1.5_r8*dw, dwmax)
         end if
         fold = fnew

         z = (w - wbar)/(sigw*sq2)
         g = exp(-z*z)
         fnmin = 1._r8

         do m = 1, nmode
            ! modal
            x     = twothird*(lnsm(m) - lnsmax)/(sq2*alogsig(m))
            fn(m) = 0.5_r8*(1._r8 - erf(x))
            fnmin = min(fn(m), fnmin)
            ! integration is second order accurate
            ! assumes linear variation of f*g with w
            fnbar = (fn(m)*g + fnold(m)*gold)
            arg   = x - 1.5_r8*sq2*alogsig(m)
            fm(m) = 0.5_r8*(1._r8 - erf(arg))
            fmbar = (fm(m)*g + fmold(m)*gold)
            wb    = (w + wold)
            if (w > 0._r8) then
               sumflxn(m) = sumflxn(m) + sixth*(wb*fnbar           &
                            + (fn(m)*g*w + fnold(m)*gold*wold))*dw
               sumflxm(m) = sumflxm(m) + sixth*(wb*fmbar           &
                            + (fm(m)*g*w + fmold(m)*gold*wold))*dw
            endif
            sumfn(m) = sumfn(m) + 0.5_r8*fnbar*dw
            ! write(iulog,'(a,9g10.2)')'lnsmax,lnsm(m),x,fn(m),fnold(m),g,gold,fnbar,dw=',lnsmax,lnsm(m),x,fn(m),fnold(m),g,gold,fnbar,dw
            fnold(m) = fn(m)
            sumfm(m) = sumfm(m) + 0.5_r8*fmbar*dw
            fmold(m) = fm(m)
         end do

         ! same form as sumflxm but replace the fm with 1.0
         sumflx_fullact = sumflx_fullact + sixth*(wb*(g + gold) + (g*w + gold*wold))*dw
         ! sumg=sumg+0.5_r8*(g+gold)*dw

         gold = g
         wold = w
         dw   = dwnew

         if (n > 1 .and. (w > wmax .or. fnmin > fmax)) exit

         w = w + dw

         if (n == nx) then
            write(iulog,*)'do loop is too short in activate'
            write(iulog,*)'wmin=',wmin,' w=',w,' wmax=',wmax,' dw=',dw
            write(iulog,*)'wbar=',wbar,' sigw=',sigw,' wdiab=',wdiab
            write(iulog,*)'wnuc=',wnuc
            write(iulog,*)'na=',(na(m),m=1,nmode)
            write(iulog,*)'fn=',(fn(m),m=1,nmode)
            !   dump all subr parameters to allow testing with standalone code
            !   (build a driver that will read input and call activate)
            write(iulog,*)'wbar,sigw,wdiab,tair,rhoair,nmode='
            write(iulog,*) wbar,sigw,wdiab,tair,rhoair,nmode
            write(iulog,*)'na=',na
            write(iulog,*)'volume=', (volume(m),m=1,nmode)
            write(iulog,*)'hydro='
            write(iulog,*) hygro
            call endrun(subname)
         end if

      end do  ! n = 1, nx


      if (w < wmaxf) then

         ! contribution from all updrafts stronger than wmax
         ! assuming constant f (close to fmax)
         wnuc = w + wdiab

         z1 = (w - wbar)/(sigw*sq2)
         z2 = (wmaxf - wbar)/(sigw*sq2)
         g  = exp(-z1*z1)
         integ = sigw*0.5_r8*sq2*sqpi*(erf(z2) - erf(z1))
         ! consider only upward flow into cloud base when estimating flux
         wf1 = max(w, zero)
         zf1 = (wf1 - wbar)/(sigw*sq2)
         gf1 = exp(-zf1*zf1)
         wf2 = max(wmaxf, zero)
         zf2 = (wf2 - wbar)/(sigw*sq2)
         gf2 = exp(-zf2*zf2)
         gf  = (gf1 - gf2)
         integf = wbar*sigw*0.5_r8*sq2*sqpi*(erf(zf2) - erf(zf1)) + sigw*sigw*gf

         do m = 1, nmode
            sumflxn(m) = sumflxn(m) + integf*fn(m)
            sumfn(m)   = sumfn(m) + fn(m)*integ
            sumflxm(m) = sumflxm(m) + integf*fm(m)
            sumfm(m)   = sumfm(m) + fm(m)*integ
         end do
         ! same form as sumflxm but replace the fm with 1.0
         sumflx_fullact = sumflx_fullact + integf
         ! sumg=sumg+integ
      end if

      do m = 1, nmode

         fn(m) = sumfn(m)/(sq2*sqpi*sigw)
         ! fn(m)=sumfn(m)/(sumg)

         if (fn(m) > 1.01_r8) then
            write(iulog,*)'fn=',fn(m),' > 1 in activate'
            write(iulog,*)'w,m,na,amcube=',w,m,na(m),amcube(m)
            write(iulog,*)'integ,sumfn,sigw=',integ,sumfn(m),sigw
            call endrun(subname)
         end if

         fluxn(m) = sumflxn(m)/(sq2*sqpi*sigw)
         fm(m)    = sumfm(m)/(sq2*sqpi*sigw)
         ! fm(m)=sumfm(m)/(sumg)
         if (fm(m) > 1.01_r8) then
            write(iulog,*)'fm=',fm(m),' > 1 in activate'
         end if
         fluxm(m) = sumflxm(m)/(sq2*sqpi*sigw)
      end do
      ! same form as fluxm
      flux_fullact = sumflx_fullact/(sq2*sqpi*sigw)

   else ! single updraft

      wnuc = wbar + wdiab

      if (wnuc > 0._r8) then

         w = wbar

         if (in_cloud) then

            if (smax_f > 0._r8) then
               smax = alpha*w/(2.0_r8*pi*rhoh2o*grow*gamma*smax_f)
            else
               smax = 1.e-20_r8
            end if

         else ! at cloud base

            alw        = alpha*wnuc
            sqrtalw    = sqrt(alw)
            etafactor1 = alw*sqrtalw

            do m = 1, nmode
               eta(m)  = etafactor1*etafactor2(m)
               zeta(m) = twothird*sqrtalw*aten/sqrtg(m)
            end do

            call maxsat(zeta, eta, nmode, smc, smax)

         end if

         lnsmax    = log(smax)

         do m = 1, nmode

            x     = twothird*(lnsm(m) - lnsmax)/(sq2*alogsig(m))
            fn(m) = 0.5_r8*(1._r8 - erf(x))
            arg   = x - 1.5_r8*sq2*alogsig(m)
            fm(m) = 0.5_r8*(1._r8 - erf(arg))

            if (wbar > 0._r8) then
               fluxn(m) = fn(m)*w
               fluxm(m) = fm(m)*w
            end if

         end do

         flux_fullact = w

      end if  ! wnuc>0

   endif  ! single updraft

end subroutine actdrop_mam_calc

!===============================================================================

subroutine maxsat(zeta, eta, nmode, smc, smax)

   ! calculates maximum supersaturation for multiple
   ! competing aerosol modes.
   !
   ! Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   ! 2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

   integer,  intent(in)  :: nmode       ! number of modes
   real(r8), intent(in)  :: smc(nmode)  ! critical supersaturation for number mode radius
   real(r8), intent(in)  :: zeta(nmode)
   real(r8), intent(in)  :: eta(nmode)
   real(r8), intent(out) :: smax        ! maximum supersaturation

   integer  :: m
   real(r8) :: sum, g1, g2, g1sqrt, g2sqrt
   !---------------------------------------------------------------------------------

   do m = 1, nmode
      if (zeta(m) > 1.e5_r8*eta(m) .or. smc(m)*smc(m) > 1.e5_r8*eta(m)) then
         ! weak forcing. essentially none activated
         smax = 1.e-20_r8
      else
         ! significant activation of this mode. calc activation all modes.
         exit
      end if

      ! No significant activation in any mode.  Do nothing.
      if (m == nmode) return

   end do

   sum = 0.0_r8
   do m = 1, nmode
      if (eta(m) > 1.e-20_r8) then
         g1     = zeta(m)/eta(m)
         g1sqrt = sqrt(g1)
         g1     = g1sqrt*g1
         g2     = smc(m)/sqrt(eta(m) + 3._r8*zeta(m))
         g2sqrt = sqrt(g2)
         g2     = g2sqrt*g2
         sum    = sum + (f1(m)*g1 + f2(m)*g2)/(smc(m)*smc(m))
      else
         sum = 1.e20_r8
      end if
   end do

   smax = 1._r8/sqrt(sum)

end subroutine maxsat

!===============================================================================

end module activate_drop_mam
