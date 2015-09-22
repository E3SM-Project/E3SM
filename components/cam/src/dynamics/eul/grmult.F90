
subroutine grmult(rcoslat ,d       ,qm1     ,tm1     ,um1     ,&
                  vm1     ,z       ,tm2     ,phis    ,dpsl    ,&
                  dpsm    ,omga    ,pdel    ,pbot    ,logpsm2 ,&
                  logpsm1 ,rpmid   ,rpdel   ,fu      ,fv      ,&
                  t2      ,ut      ,vt      ,drhs    ,pmid    ,&
                  etadot  ,etamid  ,engy    ,ddpn    ,vpdsn   ,&
                  dpslon  ,dpslat  ,vat     ,ktoop   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Non-linear dynamics calculations in grid point space
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, J. Hack, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev, plevp, plon
   use pspect
   use commap
   use physconst, only: rair, cappa, cpvir, zvir
   use hycoef, only : hybi, hybm, hybd, nprlev

   implicit none

!
! Input arguments
!
   real(r8), intent(in) :: rcoslat              ! 1./cosine(latitude)
   real(r8), intent(in) :: d(plon,plev)         ! divergence
   real(r8), intent(in) :: qm1(plon,plev)       ! specific humidity
   real(r8), intent(in) :: tm1(plon,plev)       ! temperature
   real(r8), intent(in) :: um1(plon,plev)       ! zonal wind * cos(lat)
   real(r8), intent(in) :: vm1(plon,plev)       ! meridional wind * cos(lat)
   real(r8), intent(in) :: z(plon,plev)         ! vorticity
   real(r8), intent(in) :: phis(plon)           ! surface geopotential
   real(r8), intent(in) :: dpsl(plon)           ! longitudinal component of grad ln(ps)
   real(r8), intent(in) :: dpsm(plon)           ! latitudinal component of grad ln(ps)
   real(r8), intent(in) :: omga(plon,plev)      ! vertical pressure velocity
   real(r8), intent(in) :: pdel(plon,plev)      ! layer thicknesses (pressure)
   real(r8), intent(in) :: pbot(plon)           ! bottom interface pressure
   real(r8), intent(in) :: logpsm2(plon)        ! log(psm2)
   real(r8), intent(in) :: logpsm1(plon)        ! log(ps)
   real(r8), intent(in) :: rpmid(plon,plev)     ! 1./pmid
   real(r8), intent(in) :: rpdel(plon,plev)     ! 1./pdel
   real(r8), intent(in) :: tm2(plon,plev)       ! temperature at previous time step
   integer, intent(in) :: nlon
!
! Input/Output arguments
!
   real(r8), intent(inout) :: fu(plon,plev)        ! nonlinear term - u momentum eqn
   real(r8), intent(inout) :: fv(plon,plev)        ! nonlinear term - v momentum eqn
   real(r8), intent(inout) :: t2(plon,plev)        ! nonlinear term - temperature
   real(r8), intent(inout) :: ut(plon,plev)        ! (u*TM1) - heat flux - zonal 
   real(r8), intent(inout) :: vt(plon,plev)        ! (u*TM1) - heat flux - meridional
   real(r8), intent(inout) :: drhs(plon,plev)      ! RHS of divergence eqn (del^2 term)
   real(r8), intent(inout) :: pmid(plon,plev)      ! pressure at full levels
   real(r8), intent(inout) :: etadot(plon,plevp)   ! vertical velocity in eta coordinates
   real(r8), intent(in)    :: etamid(plev)         ! midpoint values of eta (a+b)
   real(r8), intent(inout) :: engy(plon,plev)      ! kinetic energy
!
! Output arguments
!
   real(r8), intent(out) :: ddpn(plon)           ! complete sum of d*delta p
   real(r8), intent(out) :: vpdsn(plon)          ! complete sum V dot grad(ln(ps)) delta b
   real(r8), intent(out) :: dpslat(plon,plev)    ! ln(ps) component of lon press gradient
   real(r8), intent(out) :: dpslon(plon,plev)    ! ln(ps) component of lat press gradient
   real(r8), intent(out) :: vat   (plon,plev)    ! Vertical advection of temperature
   real(r8), intent(out) :: ktoop (plon,plev)    ! (Kappa*T)*(omega/P)

!
!---------------------------Local workspace-----------------------------
!
   real(r8) tv(plon,plev)        ! virtual temperature
   real(r8) ddpk(plon)           ! partial sum of d*delta p
   real(r8) vkdp                 ! V dot grad(ln(ps))
   real(r8) vpdsk(plon)          ! partial sum  V dot grad(ln(ps)) delta b
   real(r8) tk0(plon)            ! tm1 at phony level 0
   real(r8) uk0(plon)            ! u at phony level 0
   real(r8) vk0(plon)            ! v at phone level 0
   real(r8) rtv(plon,plev)       ! rair*tv
   real(r8) pterm(plon,plev)     ! intermediate term for hydrostatic eqn
   real(r8) tterm(plon,plev)     ! intermediate term for hydrostatic eqn
   real(r8) tmp                  ! temporary workspace
   real(r8) tmpk                 ! temporary workspace
   real(r8) tmpkp1               ! temporary workspace
   real(r8) edotdpde(plon,plevp) ! etadot*dp/deta
   real(r8) udel(plon,0:plev-1)  ! vertical u difference
   real(r8) vdel(plon,0:plev-1)  ! vertical v difference
   real(r8) tdel(plon,0:plev-1)  ! vertical TM1 difference

   integer i,k,kk             ! longitude, level indices
!
! Initialize arrays which represent vertical sums (ddpk, ddpn, vpdsk,
! vpdsn).  Set upper boundary condition arrays (k=0: tk0, uk0, vk0).
!
   ddpk  = 0.0_r8
   ddpn  = 0.0_r8
   vpdsk = 0.0_r8
   vpdsn = 0.0_r8
   tk0   = 0.0_r8
   uk0   = 0.0_r8
   vk0   = 0.0_r8
!
! Virtual temperature
!
tv(:nlon,:) = tm1(:nlon,:) * (1.0_r8 + zvir * qm1(:nlon,:))

!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         rtv(i,k) = rair*tv(i,k)
      end do
   end do
!
!$OMP PARALLEL DO PRIVATE (I, K, VKDP)
   do i=1,nlon
!
! sum(plev)(d(k)*dp(k))
!
      do k=1,plev
         ddpn(i) = ddpn(i) + d(i,k)*pdel(i,k)
      end do
!
! sum(plev)(v(k)*grad(lnps)*db(k))
!
      do k=nprlev,plev
         vkdp = rcoslat*(um1(i,k)*dpsl(i) + vm1(i,k)*dpsm(i))*pbot(i)
         vpdsn(i) = vpdsn(i) + vkdp*hybd(k)
      end do
!
! Compute etadot (dp/deta) (k+1/2).  Note: sum(k)(d(j)*dp(j)) required in
! pressure region. sum(k)(d(j)*dp(j)) and sum(k)(v(j)*grad(ps)*db(j)) 
! required in hybrid region
!
      edotdpde(i,1) = 0._r8
      do k=1,nprlev-1
         ddpk(i) = ddpk(i) + d(i,k)*pdel(i,k)
         edotdpde(i,k+1) = -ddpk(i)
      end do
!
      do k=nprlev,plev-1
         ddpk(i) = ddpk(i) + d(i,k)*pdel(i,k)
         vkdp = rcoslat*(um1(i,k)*dpsl(i) + vm1(i,k)*dpsm(i))*pbot(i)
         vpdsk(i) = vpdsk(i) + vkdp*hybd(k)
         edotdpde(i,k+1) = -ddpk(i) - vpdsk(i) + hybi(k+1)*(ddpn(i)+vpdsn(i))
      end do
      edotdpde(i,plevp) = 0._r8
!
!
   end do

!
! Nonlinear advection terms.  u*tm1, v*tm1, kinetic energy first
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         ut(i,k) = um1(i,k)*tm1(i,k)
         vt(i,k) = vm1(i,k)*tm1(i,k)
         engy(i,k) = 0.5_r8*(um1(i,k)**2 + vm1(i,k)**2)
      end do
   end do
!
! Compute workspace arrays for delta-u, delta-v, delta-tm1 (k)
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=0,plev-1
      if (k == 0) then
         do i=1,nlon
            udel(i,0) = um1(i,1) - uk0(i)
            vdel(i,0) = vm1(i,1) - vk0(i)
            tdel(i,0) = tm1(i,1) - tk0(i)
         end do
      else
         do i=1,nlon
            udel(i,k) = um1(i,k+1) - um1(i,k)
            vdel(i,k) = vm1(i,k+1) - vm1(i,k)
            tdel(i,k) = tm1(i,k+1) - tm1(i,k)
         end do
      endif
   end do
!
!$OMP PARALLEL DO PRIVATE (K, I, TMPK, TMPKP1, TMP)
   do k=1,plev
!
    if (k < nprlev) then
!
! Horizontal advection: u*z, v*z, energy conversion term (omega/p),
! vertical advection for interface above.  Pure pressure region first.
!
      do i=1,nlon
         dpslat(i,k) = 0._r8
         dpslon(i,k) = 0._r8
         tmpk   = 0.5_r8*rpdel(i,k)*edotdpde(i,k  )
         tmpkp1 = 0.5_r8*rpdel(i,k)*edotdpde(i,k+1)
         fu(i,k) = fu(i,k) + vm1(i,k)*z(i,k) - udel(i,k-1)*tmpk - udel(i,k  )*tmpkp1
         fv(i,k) = fv(i,k) - um1(i,k)*z(i,k) - vdel(i,k-1)*tmpk - vdel(i,k  )*tmpkp1
         vat  (i,k) = - (tdel(i,k-1)*tmpk + tdel(i,k)*tmpkp1)
         ktoop(i,k) = cappa*tv(i,k)/(1._r8 + cpvir*qm1(i,k))* &
                      omga(i,k)*rpmid(i,k)
         t2   (i,k) = t2(i,k) + d(i,k)*tm1(i,k) - tdel(i,k-1)*tmpk + &
                      ktoop(i,k) - tdel(i,k)*tmpkp1
      end do
!
    else if (k < plev) then
!
! Hybrid region above bottom level: Computations are the same as in pure
! pressure region, except that pressure gradient terms are added to 
! momentum tendencies.
!
      do i=1,nlon
         tmpk   = 0.5_r8*rpdel(i,k)*edotdpde(i,k  )
         tmpkp1 = 0.5_r8*rpdel(i,k)*edotdpde(i,k+1)
         tmp = rtv(i,k)*hybm(k)*rpmid(i,k)*pbot(i)
         dpslon(i,k) = rcoslat*tmp*dpsl(i)
         dpslat(i,k) = rcoslat*tmp*dpsm(i)
         fu(i,k) = fu(i,k) + vm1(i,k)*z(i,k) - udel(i,k-1)*tmpk - &
            udel(i,k  )*tmpkp1 - dpslon(i,k)
         fv(i,k) = fv(i,k) - um1(i,k)*z(i,k) - vdel(i,k-1)*tmpk - &
            vdel(i,k  )*tmpkp1 - dpslat(i,k)
         vat  (i,k) = - (tdel(i,k-1)*tmpk + tdel(i,k)*tmpkp1)
         ktoop(i,k) = cappa*tv(i,k)/(1._r8 + cpvir*qm1(i,k))* &
                      omga(i,k)*rpmid(i,k)
         t2   (i,k) = t2(i,k) + d(i,k)*tm1(i,k) - tdel(i,k-1)*tmpk + &
                      ktoop(i,k) - tdel(i,k)*tmpkp1
      end do
!
    else
!
! Bottom level
!
      do i=1,nlon
         tmpk = 0.5_r8*rpdel(i,plev)*edotdpde(i,plev  )
         tmp  = rtv(i,plev)*hybm(plev)*rpmid(i,plev)*pbot(i)
         dpslon(i,plev) = rcoslat*tmp*dpsl(i)
         dpslat(i,plev) = rcoslat*tmp*dpsm(i)
         fu(i,plev) = fu(i,plev) + vm1(i,plev)*z(i,plev) - &
            udel(i,plev-1)*tmpk - dpslon(i,plev)
         fv(i,plev) = fv(i,plev) - um1(i,plev)*z(i,plev) - &
            vdel(i,plev-1)*tmpk - dpslat(i,plev)
         vat  (i,plev) = -(tdel(i,plev-1)*tmpk)
         ktoop(i,plev) = cappa*tv(i,plev)/(1._r8 + cpvir*qm1(i,plev))* &
                         omga(i,plev)*rpmid(i,plev)
         t2   (i,plev) = t2(i,plev) + d(i,plev)*tm1(i,plev) - &
                         tdel(i,plev-1)*tmpk + ktoop(i,plev)
      end do
!
    end if
!
   enddo
!
! Convert eta-dot(dp/deta) to eta-dot (top and bottom = 0.)
!
   etadot(:,1) = 0._r8
   etadot(:,plevp) = 0._r8
!$OMP PARALLEL DO PRIVATE (K, TMP, I)
   do k=2,plev
      tmp = etamid(k) - etamid(k-1)
      do i=1,nlon
         etadot(i,k) = edotdpde(i,k)*tmp/(pmid(i,k) - pmid(i,k-1))
      end do
   end do
!
!-----------------------------------------------------------------------
!
! Divergence and hydrostatic equations
!
! Del squared part of RHS of divergence equation.
! Kinetic energy and diagonal term of hydrostatic equation.
! Total temperature as opposed to  perturbation temperature is acceptable
! since del-square operator will operate on this term.
! (Also store some temporary terms.)
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         tterm(i,k) = 0.5_r8*tm2(i,k) - tm1(i,k)
         pterm(i,k) = rtv(i,k)*rpmid(i,k)*pdel(i,k)
         drhs(i,k) = phis(i) + engy(i,k) + rtv(i,k)*0.5_r8* &
            rpmid(i,k)*pdel(i,k) + href(k,k)*tterm(i,k) + &
            bps(k)*(0.5_r8*logpsm2(i) - logpsm1(i))
      end do
   end do

!
! Bottom level term of hydrostatic equation
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev-1
      do i=1,nlon
         drhs(i,k) = drhs(i,k) + rtv(i,plev)* &
            rpmid(i,plev)*pdel(i,plev) + &
            href(plev,k)*tterm(i,plev)
      end do
   end do
!
! Interior terms of hydrostatic equation
!
!$OMP PARALLEL DO PRIVATE (K, KK, I)
   do k=1,plev-2
      do kk=k+1,plev-1
         do i=1,nlon
            drhs(i,k) = drhs(i,k) + pterm(i,kk) + href(kk,k)*tterm(i,kk)
         end do
      end do
   end do
!    
   return
end subroutine grmult
