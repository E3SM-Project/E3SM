
module grmult
!-----------------------------------------------------------------------
!
! Purpose:
! Compute non-linear dynamics terms in grid point space (in preparation
! for SLD interpolation). This is broken up into two interfaces an
! initialization phase for constants that will be used and the run phase
! where the non-linear terms are actually computed.
!
!---------------------------Code history--------------------------------
!
! Author:  J. Olson
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use pmgrid,          only: plon, plev, plevp
   use physconst,       only: rair
   use hycoef,          only : hyam, hybm, ps0, hyai, hybi, hypd, hybd, hypm, hypi, nprlev
   implicit none

   private
!
! Public interfaces
!
   public grmult_init   ! Initialization
   public grmult_run    ! Run method
!
! Private module data
!
   real(r8) :: hortalc(plev)           ! analytic Hortal temperature correction term
   real(r8) :: hdel   (0:plev-1)       ! del-"hortalc" in vertical
   real(r8) :: onemeps                 ! 1 - epssld (SLD decentering coefficient)
   real(r8) :: onepeps                 ! 1 + epssld (SLD decentering coefficient)
   real(r8) :: detai(plevp)            ! interval between interfaces
   real(r8) :: facm1                   ! interpolation factor for time n-1
   real(r8) :: facm2                   ! interpolation factor for time n-2

contains

subroutine grmult_init()
!-----------------------------------------------------------------------
!
! Purpose:
! Initialization of constants used for grmult.
!
!-----------------------------------------------------------------------
   use physconst, only: rga
   use scanslt,   only: epssld

!-----------------------------------------------------------------------
!
!  Local variables
!
!
   real(r8) tmp            ! tmp variable
   real(r8) tmp1           ! tmp variable
   real(r8) psref          ! reference surf pres in Hortal calc.
   real(r8) pref           ! reference pres at level k in Hortal calc
   real(r8) tnot           ! reference surf temp in Hortal calc.
   real(r8) tref           ! reference temp at level k in Hortal calc.
   real(r8) lapse          ! standard atmosphere lapse rate
   real(r8) tcut           ! cutoff temperature below which "Hortal"
   integer :: k            ! Index
!
! Initialize array of constants for Hortal temperature correction
!
   psref = 101325._r8
   tnot  = 288._r8
   lapse = 0.0065_r8
   tcut  = 216.5_r8
   tmp   = rair*lapse*rga
!
   do k = 1,plev
      pref = hyam(k)*ps0 + hybm(k)*psref
      tmp1 = pref/psref
      tref = tnot*tmp1**tmp
      if (tref > tcut) then
         hortalc(k) = hybm(k)*tmp*tref/tmp1
      else
         hortalc(k) = 0._r8
      endif
   end do
!
! Compute vertical difference
!
   hdel(0) = hortalc(1)
   do k = 1,plev-1
      hdel(k) = hortalc(k+1) - hortalc(k)
   end do
   onemeps = 1._r8 - epssld
   onepeps = 1._r8 + epssld
   facm1   =  3._r8/2._r8
   facm2   = -1._r8/2._r8
!
   do k = 1,plev
      detai  (k) = (hyai(k+1) + hybi(k+1)) - (hyai(k) + hybi(k))
   end do

end subroutine grmult_init

!
!-----------------------------------------------------------------------
!

subroutine grmult_run (ztodt   ,rcoslat ,coslat  ,coriol  ,div     , &
                   q3      ,t3      ,u3      ,v3      , &
                   ed1     ,ql      ,qm      ,tl      , &
                   tm      ,phis    ,phisl   ,phism   ,dpsl    , &
                   dpsm    ,omga    ,pmid    ,pdel    ,ps      , &
                   logps   ,rpmid   ,etamid  ,fu      ,fv      , &
                   t2      ,t3m1    ,u3m1    ,v3m1    ,etadot  , &
                   dpslon  ,dpslat  ,urhs    ,vrhs    ,trhs    , &
                   prhs    ,lnpssld ,prhssld ,tarrsld , &
                   parrsld ,ktoop   ,nlon    )
!-----------------------------------------------------------------------
!
! Purpose:
! Compute non-linear dynamics terms in grid point space (in preparation
! for SLD interpolation)
!
!
!-----------------------------------------------------------------------
  use physconst,    only: zvir, cappa, cpvir
  use time_manager, only: is_first_step
  use commap,       only: bps, t0, href


!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: ztodt                 ! delta-t
  real(r8), intent(in)   :: rcoslat               ! 1/(cos lat)
  real(r8), intent(in)   :: coslat                ! cos(lat)
  real(r8), intent(in)   :: coriol                ! Coriolis parameter
  real(r8), intent(in)   :: div     (plon,plev)   ! divergence
  real(r8), intent(in)   :: q3      (plon,plev)   ! specific humidity
  real(r8), intent(in)   :: t3      (plon,plev)   ! temperature
  real(r8), intent(in)   :: u3      (plon,plev)   ! zonal wind
  real(r8), intent(in)   :: v3      (plon,plev)   ! meridional wind
  real(r8), intent(in)   :: ed1     (plon,plev)   ! (1/ps)etadot(dp/deta) (time n-1)
  real(r8), intent(in)   :: ql      (plon,plev)   ! longitudinal derivative of q
  real(r8), intent(in)   :: qm      (plon,plev)   ! latitudinal  derivative of q
  real(r8), intent(in)   :: tl      (plon,plev)   ! zonal derivative of T
  real(r8), intent(in)   :: tm      (plon,plev)   ! meridional derivative of T
  real(r8), intent(in)   :: phis    (plon)        ! Phi at surface
  real(r8), intent(in)   :: phisl   (plon)        ! longitudinal derivative of phis
  real(r8), intent(in)   :: phism   (plon)        ! latitudinal  devivative of phis
  real(r8), intent(in)   :: dpsl    (plon)        ! longitudinal component of grad ln(ps)
  real(r8), intent(in)   :: dpsm    (plon)        ! latitudinal  component of grad ln(ps)
  real(r8), intent(in)   :: omga    (plon,plev)   ! vertical pressure velocity
  real(r8), intent(in)   :: pmid    (plon,plev)   ! pressure at full levels
  real(r8), intent(in)   :: pdel    (plon,plev)   ! layer thicknesses (pressure)
  real(r8), intent(in)   :: ps      (plon)        ! Surface pressure (n)
  real(r8), intent(in)   :: logps   (plon)        ! log(ps)
  real(r8), intent(in)   :: rpmid   (plon,plev)   ! 1./pmid
  real(r8), intent(in)   :: etamid  (plev)        ! midpoint values of eta (a+b)
  real(r8), intent(inout):: fu      (plon,plev)   ! nonlinear term - u momentum eqn
  real(r8), intent(inout):: fv      (plon,plev)   ! nonlinear term - v momentum eqn
  real(r8), intent(inout):: t2      (plon,plev)   ! nonlinear term - temperature
  real(r8), intent(out)  :: t3m1    (plon,plev)   ! T at time n-1
  real(r8), intent(inout):: u3m1    (plon,plev)   ! U at previous time step
  real(r8), intent(inout):: v3m1    (plon,plev)   ! V at previous time step
  real(r8), intent(out)  :: etadot  (plon,plevp)  ! vertical velocity in eta coordinates
  real(r8), intent(out)  :: dpslon  (plon,plev)   ! Pressure gradient term
  real(r8), intent(out)  :: dpslat  (plon,plev)   ! Pressure gradient term
  real(r8), intent(inout):: urhs    (plon,plev)   ! RHS of U  eqn valid for mid-point
  real(r8), intent(inout):: vrhs    (plon,plev)   ! RHS of V  eqn valid for mid-point
  real(r8), intent(inout):: trhs    (plon,plev)   ! RHS of T  eqn valid for mid-point
  real(r8), intent(inout):: prhs    (plon,plev)   ! RHS of Ps eqn valid for mid-point
  real(r8), intent(out)  :: lnpssld (plon,plev)   ! RHS Ps term for SLD
  real(r8), intent(out)  :: prhssld (plon,plev)   ! RHS Ps term for SLD
  real(r8), intent(out)  :: tarrsld (plon,plev)   ! T  at arr. pt. (SLD)
  real(r8), intent(out)  :: parrsld (plon,plev)   ! Ps at arr. pt. (SLD)
  real(r8), intent(out)  :: ktoop   (plon,plev)   ! (Kappa*T)*(omega/P)
  integer , intent(in)   :: nlon                  ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  real(r8) tmp1                  ! temporary workspace
  real(r8) tmp2                  ! temporary workspace
  real(r8) tmp                   ! temporary workspace
  real(r8) tmpk                  ! workspace
  real(r8) tmpkp1                ! workspace
  real(r8) u3l     (plon,plev)   ! u-wind used locally only
  real(r8) v3l     (plon,plev)   ! v-wind used locally only
  real(r8) tv      (plon,plev)   ! virtual temperature
  real(r8) ddpk    (plon)        ! partial sum of div*delta p
  real(r8) ddpn    (plon)        ! complete sum of div*delta p
  real(r8) vkdp    (plon,plev)   ! V dot grad(ln(ps))
  real(r8) vpdsk   (plon)        ! partial sum  V dot grad(ln(ps)) delta b
  real(r8) vpdsn   (plon)        ! complete sum V dot grad(ln(ps)) delta b
  real(r8) lpsstar (plon)        ! Reference ln(Ps) (used to define a new 
!                                ! perturbation Ps)
  real(r8) lpsstarl(plon)        ! long. grad of reference ln(Ps)
  real(r8) lpsstarm(plon)        ! lat.  grad of reference ln(Ps)
  real(r8) rtv     (plon,plev)   ! rair*(tv+t0)
  real(r8) dt                    ! time step
  real(r8) hsl     (plon,plev)   ! zonal      deriv of hydrostatic term
  real(r8) hsm     (plon,plev)   ! meridional deriv of hydrostatic term
  real(r8) pspsl   (plon)        ! Ps*d(lnPs)/d(long.)
  real(r8) pspsm   (plon)        ! Ps*d(lnPs)/d(lat. )
  real(r8) ed1p    (plon,plevp)  ! (1/ps)etadot(dp/deta) (time n-1)
  real(r8) rtvl                  ! zonal      derivative of R*Tv
  real(r8) rtvm                  ! meridional derivative of R*Tv
  real(r8) abp0                  ! constant for grad(H(n)) matrix
!
! Arrays which hold results from the RHS of the prognostic equations.
! Results will be interpolated to trajectory departure points in the routine
! SCANSLT.
!
  real(r8) tsld0a  (plon,plev)   ! RHS of T  eqn valid for mid-point
  real(r8) usldm   (plon,plev)   ! RHS of U  eqn valid for departure pt (n)
  real(r8) vsldm   (plon,plev)   ! RHS of V  eqn valid for departure pt (n)
  real(r8) tsldm   (plon,plev)   ! RHS of T  eqn valid for departure pt (n)
  real(r8) psldm   (plon,plev)   ! RHS of Ps eqn valid for departure pt (n)
  real(r8) urhsl   (plon,plev)   ! RHS of U  eqn valid for mid-point (n+1/2)
  real(r8) vrhsl   (plon,plev)   ! RHS of V  eqn valid for mid-point (n+1/2)
  real(r8) trhsl   (plon,plev)   ! RHS of T  eqn valid for mid-point (n+1/2)
  real(r8) prhsl   (plon,plev)   ! RHS of Ps eqn valid for mid-point (n+1/2)


  integer i,l,k                  ! longitude, level indices
  integer npr                    ! index
!
!-----------------------------------------------------------------------
!
!
! Compute U/V for time n + 1/2
! This will be used locally in the Hortal Temperature correction
!
  do k = 1,plev
     do i = 1,nlon
        u3l  (i,k) = facm1*u3(i,k) + facm2*u3m1(i,k)
        v3l  (i,k) = facm1*v3(i,k) + facm2*v3m1(i,k)
     end do
  end do
!
! Zero auxiliary fields
!
  tmp = 1._r8/(rair*t0(plev))
  do i=1,nlon
     ddpk    (i)     = 0._r8
     ddpn    (i)     = 0._r8
     vpdsk   (i)     = 0._r8
     vpdsn   (i)     = 0._r8
     lpsstar (i)     = -phis (i)*tmp
     lpsstarl(i)     = -phisl(i)*tmp
     lpsstarm(i)     = -phism(i)*tmp
     pspsl   (i)     = ps(i)*dpsl(i)
     pspsm   (i)     = ps(i)*dpsm(i)
     etadot(i,1)     = 0._r8
     ed1p  (i,1)     = 0._r8
     etadot(i,plevp) = 0._r8
  end do
!
! Virtual temperature
!
tv(:nlon,:) = t3(:nlon,:) * (1.0_r8 + zvir*q3(:nlon,:))
!
! calculate some auxiliary quantities
!
  do k=1,plev
     do i=1,nlon
        ed1p(i,k+1) = ed1(i,k)
        rtv(i,k) = rair*tv(i,k)
!
! sum(plev)(div(k)*dp(k))
!
        ddpn(i) = ddpn(i) + div(i,k)*pdel(i,k)
     end do
  end do
!
! sum(plev)(v(k)*grad(lnps)*db(k))
!
  do k=nprlev,plev
     do i=1,nlon
        vkdp(i,k)= rcoslat*(u3(i,k)*pspsl(i) + v3(i,k)*pspsm(i))
        vpdsn(i) = vpdsn(i) + vkdp(i,k)*hybd(k)
     end do
  end do
!
! Compute etadot (top and bottom = 0.)
!
  do k = 1,plev-1
!
! Compute etadot(dp/deta)(k+1/2) and sum(k)(div(j)*dp(j))
!
     do i=1,nlon
        ddpk(i) = ddpk(i) + div(i,k)*pdel(i,k)
        etadot(i,k+1) = -ddpk(i)
     end do
!
! sum(k)(v(j)*grad(ps)*db(j))
!
     if (k.ge.nprlev) then
        do i=1,nlon
           vpdsk(i) = vpdsk(i) + vkdp(i,k)*hybd(k)
           etadot(i,k+1) = etadot(i,k+1) - vpdsk(i) + hybi(k+1)*(ddpn(i)+vpdsn(i))
        end do
     end if
!
! Convert eta-dot(dp/deta) to eta-dot
!
     tmp = etamid(k+1) - etamid(k)
     do i = 1,nlon
        etadot(i,k+1) = etadot(i,k+1)*tmp/(pmid(i,k+1) - pmid(i,k))
     end do
  end do
!
! Zonal and meridional derivatives of the hydrostatic term in the
! momentum equations.
!
  do k = 1,plev
     do i = 1,nlon
        tmp1     = (1._r8 + zvir*q3(i,k))
        tmp2     = t3(i,k)*zvir
        rtvl     = rair*( tmp1*tl(i,k) + tmp2*ql(i,k) )
        rtvm     = rair*( tmp1*tm(i,k) + tmp2*qm(i,k) )
!
        tmp      = rpmid(i,k)*pdel(i,k)
        hsl(i,k) = -rtvl*tmp
        hsm(i,k) = -rtvm*tmp
     end do
     if(k .ge. nprlev) then
        abp0 = ps0*(hyam(k)*(hybi(k+1) - hybi(k)) - hybm(k)*(hyai(k+1) - hyai(k)))
        do i = 1,nlon
           tmp      = rtv(i,k)*abp0/(pmid(i,k)*pmid(i,k))
           hsl(i,k) = hsl(i,k) - pspsl(i)*tmp
           hsm(i,k) = hsm(i,k) - pspsm(i)*tmp
        end do
     endif
  end do
!
! Calculate RHS of all prognostic eqns
!
  do k = 1,plev
     tmp  = cappa*t0(k)*hypi(plevp)/hypm(k)
     tmp1 = hypd(k)/hypi(plevp)
     do i = 1,nlon
!
! Surface Pressure eqn
!
        prhsl(i,k) =  div(i,k)*(tmp1 - pdel(i,k)/ps(i))
        psldm(i,k) = -div(i,k)*tmp1 - ed1p(i,k+1) + ed1p(i,k)
!
! Temperature eqn
!
        ktoop (i,k) = cappa*tv(i,k)/(1._r8 + cpvir*q3(i,k))*omga(i,k)*rpmid(i,k)
        trhsl (i,k) = ktoop(i,k)
        trhsl (i,k) = trhsl(i,k) - omga(i,k)/ps(i)*tmp
        tsldm (i,k) = 0.5_r8*(ed1p(i,k+1) + ed1p(i,k))*tmp
!
! ... horizontal advection portion of Hortal Temperature correction
!
        trhsl (i,k) = trhsl(i,k) - &
                      rcoslat*( u3(i,k)*lpsstarl(i) + v3(i,k)*lpsstarm(i) )*hortalc(k)
!
! ... Ritchie damping term for Temperature eqn.
!
        tsld0a(i,k) = rcoslat*( u3l(i,k)*lpsstarl(i) + v3l(i,k)*lpsstarm(i) )
!
! U/V eqns (includes only the diagonal portion of the hydrostatic term)
!
        urhsl(i,k) = 0.5_r8*hsl(i,k) + href(k,k)*tl(i,k) + bps(k)*dpsl(i)
        vrhsl(i,k) = 0.5_r8*hsm(i,k) + href(k,k)*tm(i,k) + bps(k)*dpsm(i)
        usldm(i,k) = -phisl(i) + v3(i,k)*coriol*coslat -href(k,k)*tl(i,k) - bps(k)*dpsl(i)
        vsldm(i,k) = -phism(i) - u3(i,k)*coriol*coslat -href(k,k)*tm(i,k) - bps(k)*dpsm(i)
     end do
!
! Add pressure gradient terms to momentum tendencies
!
     if (k.ge.nprlev) then
        do i=1,nlon
           tmp        = rtv(i,k)*hybm(k)*rpmid(i,k)
           dpslon(i,k) = rcoslat*tmp*pspsl(i)
           dpslat(i,k) = rcoslat*tmp*pspsm(i)
           urhsl(i,k) = urhsl(i,k) - dpslon(i,k)*coslat
           vrhsl(i,k) = vrhsl(i,k) - dpslat(i,k)*coslat
        end do
     else
        do i = 1,nlon
           dpslon(i,k) = 0._r8
           dpslat(i,k) = 0._r8
        end do
     end if
  end do
!
! Interior levels of the hydrostatic term
!
  do k=1,plev-1
     do l=k+1,plev
        do i=1,nlon
           urhsl(i,k) = urhsl(i,k) + href(l,k)*tl(i,l) + hsl(i,l)
           vrhsl(i,k) = vrhsl(i,k) + href(l,k)*tm(i,l) + hsm(i,l)
!
           usldm(i,k) = usldm(i,k) - href(l,k)*tl(i,l)
           vsldm(i,k) = vsldm(i,k) - href(l,k)*tm(i,l)
        end do
     end do
  end do
!
! Compute vertical advection portion of Hortal Temperature correction
!
  npr = nprlev - 1
  if(npr .lt. 1) npr = 1
  do k = npr,plev-1
     tmpk   = 0.5_r8*hdel(k-1)/detai(k)
     tmpkp1 = 0.5_r8*hdel(k  )/detai(k)
     do i = 1,nlon
        trhsl(i,k) = trhsl(i,k) - (etadot(i,k  )*tmpk + etadot(i,k+1)*tmpkp1)*lpsstar(i)
     end do
  end do
!
! ... bottom level
!
  tmpk = 0.5_r8*hdel(plev-1)/detai(plev)
  do i = 1,nlon
     trhsl(i,plev) = trhsl(i,plev) - etadot(i,plev)*tmpk*lpsstar(i)
  end do
!
! Compute (by extrapolation) RHS terms for time n + 1/2
!
  if(is_first_step()) then
     do k = 1,plev
        do i = 1,nlon
           urhs  (i,k) = urhsl (i,k)
           vrhs  (i,k) = vrhsl (i,k)
           trhs  (i,k) = trhsl (i,k)
           prhs  (i,k) = prhsl (i,k)
        end do
     end do
  endif
!
  do k = 1,plev
     do i = 1,nlon
        tmp         =       urhsl (i,k)
        urhsl (i,k) = facm1*urhsl (i,k) + facm2*urhs  (i,k)
        urhs  (i,k) = tmp
        tmp         =       vrhsl (i,k)
        vrhsl (i,k) = facm1*vrhsl (i,k) + facm2*vrhs  (i,k)
        vrhs  (i,k) = tmp
        tmp         =       trhsl (i,k)
        trhsl (i,k) = facm1*trhsl (i,k) + facm2*trhs  (i,k)
        trhs  (i,k) = tmp
        tmp         =       prhsl (i,k)
        prhsl (i,k) = facm1*prhsl (i,k) + facm2*prhs  (i,k)
        prhs  (i,k) = tmp
     end do
  end do
!
! Prepare SLD arrays for interpolation by SCANSLT
!
! Append appropriate coefficients (including decentering epsilon values -- See
! Notes 13,14,15)
!
  dt   = 0.5_r8*ztodt
  tmp1 = 0.5_r8*ztodt/coslat
  do k = 1,plev
     tmp = cappa*t0(k)*hypi(plevp)/hypm(k)
     do i = 1,nlon
        urhsl (i,k) = tmp1*urhsl (i,k)
        vrhsl (i,k) = tmp1*vrhsl (i,k)
        trhsl (i,k) = dt  *trhsl (i,k)
        tsld0a(i,k) = dt  *tsld0a(i,k)
        prhsl (i,k) = dt  *prhsl (i,k)
        t2    (i,k) = ztodt*t2   (i,k)
        fu    (i,k) = ztodt*fu   (i,k)
        fv    (i,k) = ztodt*fv   (i,k)
!
! Combine terms.
! (Time split the midpoint results between arrival and departure
!  points)
!
        tmp2          = lpsstar(i)*hortalc(k)
        u3m1    (i,k) = (u3  (i,k)*coslat + onemeps*usldm(i,k)*dt)/coslat + &
                                            onemeps*urhsl(i,k) + &
                                            onemeps*fu   (i,k)*0.5_r8
        v3m1    (i,k) = (v3  (i,k)*coslat + onemeps*vsldm(i,k)*dt)/coslat + &
                                            onemeps*vrhsl(i,k) + &
                                            onemeps*fv   (i,k)*0.5_r8
        t3m1    (i,k) = t3  (i,k) +         onemeps*tsldm(i,k)*dt + &
                                            onemeps*trhsl(i,k) - tmp2 + &
                                            onemeps*t2(i,k)*0.5_r8
        prhssld (i,k) = (psldm(i,k)*dt + prhsl(i,k))*onemeps
        tarrsld (i,k) = (trhsl(i,k) + &
                         hybm(k)*tmp*tsld0a(i,k))*onepeps + tmp2 + onepeps*t2(i,k)*0.5_r8
        parrsld (i,k) = (prhsl(i,k) - hybd(k)*tsld0a(i,k))*onepeps
        lnpssld (i,k) = logps(i) - lpsstar(i) - tsld0a(i,k)*onemeps
        fu      (i,k) = urhsl(i,k)*coslat*onepeps + fu(i,k)*coslat*onepeps*0.5_r8
        fv      (i,k) = vrhsl(i,k)*coslat*onepeps + fv(i,k)*coslat*onepeps*0.5_r8
     end do
  end do
!
  return
end subroutine grmult_run

!
!-----------------------------------------------------------------------
!

end module grmult
