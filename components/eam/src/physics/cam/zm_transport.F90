module zm_transport
   !----------------------------------------------------------------------------
   ! 
   ! Transport routines for the Zhang-McFarlane deep convection scheme
   !
   !----------------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use cam_abortutils,  only: endrun
   use constituents,    only: cnst_get_type_byind
   use zm_conv,         only: zm_microp
   use cam_logfile,     only: iulog

   implicit none

   ! public methods
   public :: zm_transport_tracer    ! convective tracer transport
   public :: zm_transport_momentum  ! convective momentum transport

   private

   real(r8), parameter :: mbsth = 1.e-15_r8  ! threshold below which we treat the mass fluxes as zero (in mb/s)
   
contains

!===================================================================================================

subroutine zm_transport_tracer( doconvtran, q, ncnst, &
                                mu, md, du, eu, ed, dp, &
                                jt, mx, ideep, il1g, il2g, &
                                fracis, dqdt, dpdry, dt ) 
   !---------------------------------------------------------------------------- 
   ! Purpose: Convective transport of tracer species
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                               intent(in)  :: ncnst       ! number of tracers to transport
   logical,  dimension(ncnst),            intent(in)  :: doconvtran  ! flag for doing convective transport
   real(r8), dimension(pcols,pver,ncnst), intent(in)  :: q           ! tracer array (including water vapor)
   real(r8), dimension(pcols,pver),       intent(in)  :: mu          ! mass flux up
   real(r8), dimension(pcols,pver),       intent(in)  :: md          ! mass flux down
   real(r8), dimension(pcols,pver),       intent(in)  :: du          ! mass detraining from updraft
   real(r8), dimension(pcols,pver),       intent(in)  :: eu          ! mass entraining from updraft
   real(r8), dimension(pcols,pver),       intent(in)  :: ed          ! mass entraining from downdraft
   real(r8), dimension(pcols,pver),       intent(in)  :: dp          ! delta pressure between interfaces
   real(r8), dimension(pcols,pver,ncnst), intent(in)  :: fracis      ! fraction of tracer that is insoluble
   integer,  dimension(pcols),            intent(in)  :: jt          ! index of cloud top for each column
   integer,  dimension(pcols),            intent(in)  :: mx          ! index of cloud top for each column
   integer,  dimension(pcols),            intent(in)  :: ideep       ! gathering array
   integer,                               intent(in)  :: il1g        ! gathered min ncol index
   integer,                               intent(in)  :: il2g        ! gathered max ncol index
   real(r8), dimension(pcols,pver),       intent(in)  :: dpdry       ! delta pressure between interfaces
   real(r8),                              intent(in)  :: dt          ! model time increment
   real(r8), dimension(pcols,pver,ncnst), intent(out) :: dqdt        ! output tendency array
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i,k                   ! loop indeces
   integer  :: kbm                   ! Highest altitude index of cloud base
   integer  :: kk                    ! Work index
   integer  :: kkp1                  ! Work index
   integer  :: km1                   ! Work index
   integer  :: kp1                   ! Work index
   integer  :: ktm                   ! Highest altitude index of cloud top
   integer  :: m                     ! Work index
   real(r8) :: cabv                 ! Mix ratio of constituent above
   real(r8) :: cbel                 ! Mix ratio of constituent below
   real(r8) :: cdifr                ! Normalized diff between cabv and cbel
   real(r8) :: chat(pcols,pver)     ! Mix ratio in env at interfaces
   real(r8) :: cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
   real(r8) :: const(pcols,pver)    ! Gathered tracer array
   real(r8) :: fisg(pcols,pver)     ! gathered insoluble fraction of tracer
   real(r8) :: conu(pcols,pver)     ! Mix ratio in updraft at interfaces
   real(r8) :: dcondt(pcols,pver)   ! Gathered convective tendency array
   real(r8) :: mupdudp              ! A work variable
   real(r8) :: minc                 ! A work variable
   real(r8) :: maxc                 ! A work variable
   real(r8) :: fluxin               ! A work variable
   real(r8) :: fluxout              ! A work variable
   real(r8) :: netflux              ! A work variable
   real(r8) :: dutmp(pcols,pver)    ! Mass detraining from updraft
   real(r8) :: eutmp(pcols,pver)    ! Mass entraining into updraft
   real(r8) :: edtmp(pcols,pver)    ! Mass entraining into downdraft
   real(r8) :: dptmp(pcols,pver)    ! Delta pressure between interfaces
   real(r8) :: negadt               ! for Conservation check
   real(r8) :: qtmp                 ! for Conservation check
   ! constants
   real(r8), parameter :: small        = 1.e-36_r8 ! a small number to avoid division by zero
   real(r8), parameter :: cdifr_min    = 1.e-6_r8  ! minimum layer difference for geometric averaging
   real(r8), parameter :: maxc_factor  = 1.e-12_r8
   real(r8), parameter :: flux_factor  = 1.e-12_r8
   !----------------------------------------------------------------------------
   
   ! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

   ! Loop ever each constituent (skip water vapor at m=1)
   do m = 2, ncnst

      if (doconvtran(m)) then

         if (cnst_get_type_byind(m).eq.'dry') then
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dpdry(i,k)
                  dutmp(i,k) = du(i,k)*dp(i,k)/dpdry(i,k)
                  eutmp(i,k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                  edtmp(i,k) = ed(i,k)*dp(i,k)/dpdry(i,k)
               end do
            end do
         else
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dp(i,k)
                  dutmp(i,k) = du(i,k)
                  eutmp(i,k) = eu(i,k)
                  edtmp(i,k) = ed(i,k)
               end do
            end do
         endif

         ! Gather up the constituent and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
               fisg(i,k) = fracis(ideep(i),k,m)
            end do
         end do

         ! From now on work only with gathered data

         ! Interpolate environment tracer values to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc < 0) then
                  cdifr = 0._r8
               else
                  cdifr = abs( const(i,k) - const(i,km1) )/max(maxc,small)
               endif
               ! If the two layers differ significantly use a geometric averaging
               if (cdifr > cdifr_min) then
                  cabv = max(const(i,km1),maxc*maxc_factor)
                  cbel = max(const(i,k  ),maxc*maxc_factor)
                  chat(i,k) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel
               else ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5_r8*( const(i,k) + const(i,km1) )
               end if
               ! Provisional updraft and downdraft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)
               ! provisional tendencies
               dcondt(i,k) = 0._r8
            end do
         end do

         ! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         do i = il1g,il2g
            mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
            if (mupdudp > mbsth) then
               conu(i,kk) = ( +eutmp(i,kk)*fisg(i,kk)*const(i,kk)*dptmp(i,kk) )/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) = ( -edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1) )/md(i,k)
            endif
         end do

         ! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
               if (mupdudp > mbsth) then
                  conu(i,kk) = ( mu(i,kkp1)*conu(i,kkp1) + eutmp(i,kk)*fisg(i,kk)*const(i,kk)*dptmp(i,kk) )/mupdudp
               endif
            end do
         end do

         ! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                  cond(i,k) =  ( md(i,km1)*cond(i,km1) - edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1) )/md(i,k)
               endif
            end do
         end do

         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g

               ! limit fluxes outside convection to mass in appropriate layer
               ! these limiters are probably only safe for positive definite quantitities
               ! it assumes that mu and md already satify a courant number limit of 1
               fluxin  =   mu(i,kp1)*conu(i,kp1) + mu(i,k  )*min(chat(i,k  ),const(i,km1)) &
                         -(md(i,k  )*cond(i,k  ) + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
               fluxout =   mu(i,k  )*conu(i,k  ) + mu(i,kp1)*min(chat(i,kp1),const(i,k  )) &
                         -(md(i,kp1)*cond(i,kp1) + md(i,k  )*min(chat(i,k  ),const(i,k  )))

               netflux = fluxin - fluxout
               if (abs(netflux) < max(fluxin,fluxout)*flux_factor) then
                  netflux = 0._r8
               endif
               dcondt(i,k) = netflux/dptmp(i,k)
            end do
         end do

#ifdef CPRCRAY
!DIR$ NOINTERCHANGE
#endif
         do k = kbm,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k == mx(i)) then
                  fluxin  = mu(i,k)*min(chat(i,k),const(i,km1)) - md(i,k)*cond(i,k)
                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*min(chat(i,k),const(i,k))
                  netflux = fluxin - fluxout
                  if (abs(netflux) < max(fluxin,fluxout)*flux_factor) then
                     netflux = 0._r8
                  endif
                  dcondt(i,k) = netflux/dptmp(i,k)
               else if (k > mx(i)) then
                  dcondt(i,k) = 0._r8
               end if
            end do
         end do

         ! Conservation check for ZM microphysics
         if (zm_microp) then
            do i = il1g,il2g
               do k = jt(i),mx(i)
                  if (dcondt(i,k)*dt+const(i,k)<0._r8) then
                     negadt = dcondt(i,k)+const(i,k)/dt
                     dcondt(i,k) = -const(i,k)/dt
                     do kk= k+1, mx(i)
                        if (negadt<0._r8 .and. dcondt(i,kk)*dt+const(i,kk)>0._r8 ) then
                           qtmp = dcondt(i,kk)+negadt*dptmp(i,k)/dptmp(i,kk)
                           if (qtmp*dt+const(i,kk)>0._r8) then
                              dcondt(i,kk)= qtmp
                              negadt=0._r8
                           else
                              negadt= negadt+(const(i,kk)/dt+dcondt(i,kk))*dptmp(i,kk)/dptmp(i,k)
                              dcondt(i,kk)= -const(i,kk)/dt
                           end if
                        end if
                     end do
                     do kk= k-1, jt(i), -1
                        if (negadt<0._r8 .and. dcondt(i,kk)*dt+const(i,kk)>0._r8 ) then
                           qtmp = dcondt(i,kk)+negadt*dptmp(i,k)/dptmp(i,kk)
                           if (qtmp*dt+const(i,kk)>0._r8) then
                              dcondt(i,kk)= qtmp
                              negadt=0._r8
                           else
                              negadt= negadt+(const(i,kk)/dt+dcondt(i,kk))*dptmp(i,kk)/dptmp(i,k)
                              dcondt(i,kk)= -const(i,kk)/dt
                           end if
                        end if
                     end do
                     if (negadt<0._r8) then
                        dcondt(i,k) = dcondt(i,k) - negadt
                     end if
                  end if
               end do
            end do
         end if

         ! Initialize output tendency to zero, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8
         do k = 1,pver
            kp1 = min(pver,k+1)
#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif
            do i = il1g,il2g
               dqdt(ideep(i),k,m) = dcondt(i,k)
            end do
         end do

      end if ! for doconvtran

   end do ! m = 2, ncnst

   return

end subroutine zm_transport_tracer

!===================================================================================================

subroutine zm_transport_momentum( ncol, wind_in, nwind, &
                                  mu, md, du, eu, ed, dp, &
                                  jt, mx, ideep, il1g, il2g, &
                                  wind_tend, pguall, pgdall, icwu, icwd, dt, seten )
   !---------------------------------------------------------------------------- 
   ! Purpose: Convective transport of momentum
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                               intent(in)  :: ncol        ! number of atmospheric columns
   integer,                               intent(in)  :: nwind       ! number of tracers to transport
   real(r8), dimension(pcols,pver,nwind), intent(in)  :: wind_in     ! input Momentum array
   real(r8), dimension(pcols,pver),       intent(in)  :: mu          ! mass flux up
   real(r8), dimension(pcols,pver),       intent(in)  :: md          ! mass flux down
   real(r8), dimension(pcols,pver),       intent(in)  :: du          ! mass detraining from updraft
   real(r8), dimension(pcols,pver),       intent(in)  :: eu          ! mass entraining from updraft
   real(r8), dimension(pcols,pver),       intent(in)  :: ed          ! mass entraining from downdraft
   real(r8), dimension(pcols,pver),       intent(in)  :: dp          ! gathered pressure delta between interfaces
   real(r8),                              intent(in)  :: dt          ! time step in seconds : 2*delta_t
   integer,  dimension(pcols),            intent(in)  :: jt          ! index of cloud top for each column
   integer,  dimension(pcols),            intent(in)  :: mx          ! index of cloud top for each column
   integer,  dimension(pcols),            intent(in)  :: ideep       ! gathering array
   integer,                               intent(in)  :: il1g        ! gathered min ncol index
   integer,                               intent(in)  :: il2g        ! gathered max ncol index
   real(r8), dimension(pcols,pver,nwind), intent(out) :: wind_tend   ! output momentum tendency
   real(r8), dimension(pcols,pver,nwind), intent(out) :: pguall      ! apparent force from  updraft PG
   real(r8), dimension(pcols,pver,nwind), intent(out) :: pgdall      ! apparent force from  downdraft PG
   real(r8), dimension(pcols,pver,nwind), intent(out) :: icwu        ! in-cloud winds in updraft
   real(r8), dimension(pcols,pver,nwind), intent(out) :: icwd        ! in-cloud winds in downdraft
   real(r8), dimension(pcols,pver),       intent(out) :: seten       ! dry static energy tendency
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i,k,m                      ! loop indices
   integer  :: kbm                        ! Highest altitude index of cloud base
   integer  :: kk                         ! Work index
   integer  :: kkp1                       ! Work index
   integer  :: kkm1                       ! Work index
   integer  :: km1                        ! Work index
   integer  :: kp1                        ! Work index
   integer  :: ktm                        ! Highest altitude index of cloud top
   real(r8) :: wind0(pcols,pver,nwind)    ! gathered initial wind
   real(r8) :: windf(pcols,pver,nwind)    ! gathered final wind
   real(r8) :: wind_mid(pcols,pver)       ! gathered momentum in environment at mid-points
   real(r8) :: wind_int(pcols,pver)       ! gathered momentum in environment at interfaces
   real(r8) :: wind_int_d(pcols,pver)     ! gathered momentum in downdrafts at interfaces
   real(r8) :: wind_int_u(pcols,pver)     ! gathered momentum in updrafts at interfaces
   real(r8) :: wind_tend_tmp(pcols,pver)  ! gathered provisional momentum tendency
   real(r8) :: mupdudp                    ! work variable
   real(r8) :: mududp(pcols,pver)         ! work variable
   real(r8) :: mddudp(pcols,pver)         ! work variable
   real(r8) :: pgu(pcols,pver)            ! pressure gradient term for updraft
   real(r8) :: pgd(pcols,pver)            ! pressure gradient term for downdraft
   real(r8) :: gseten(pcols,pver)         ! gathered dry static energy tendency
   real(r8) :: mflux(pcols,pverp,nwind)   ! gathered momentum flux
   real(r8) :: ubot,utop                  ! U momentum at top/bot of layer
   real(r8) :: vbot,vtop                  ! V momentum at top/bot of layer
   real(r8) :: fkeb,fket                  ! flux of KE at top/bot of layer
   real(r8) :: ketend_cons                ! KE tendency for energy fixer
   real(r8) :: ketend                     ! KE tendency for energy fixer
   ! constants
   real(r8), parameter :: momcu = 0.4_r8  ! pressure gradient term constant for updrafts
   real(r8), parameter :: momcd = 0.4_r8  ! pressure gradient term constant for downdrafts
   !----------------------------------------------------------------------------
   
   ! Initialize variables
   pguall(1:ncol,1:pver, 1:nwind) = 0.0_r8
   pgdall(1:ncol,1:pver, 1:nwind) = 0.0_r8
   mflux (1:ncol,1:pverp,1:nwind) = 0.0_r8
   wind0 (1:ncol,1:pver, 1:nwind) = 0.0_r8
   windf (1:ncol,1:pver, 1:nwind) = 0.0_r8
   seten (1:ncol,1:pver)          = 0.0_r8
   gseten(1:ncol,1:pver)          = 0.0_r8

   ! Initialize in-cloud winds to environmental wind
   icwu(1:ncol,1:pver,1:nwind) = wind_in(1:ncol,1:pver,1:nwind)
   icwd(1:ncol,1:pver,1:nwind) = wind_in(1:ncol,1:pver,1:nwind)

   ! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

   ! Loop ever each wind component
   do m = 1, nwind

      ! Gather up the winds and set tend to zero
      do k = 1,pver
         do i =il1g,il2g
            wind_mid(i,k) = wind_in(ideep(i),k,m)
            wind0(i,k,m) = wind_mid(i,k)
         end do
      end do

      ! From now on work only with gathered data

      ! Interpolate winds to interfaces
      do k = 1,pver
         km1 = max(1,k-1)
         do i = il1g, il2g
            ! use arithmetic mean
            wind_int(i,k) = 0.5_r8*( wind_mid(i,k) + wind_mid(i,km1) )
            ! Provisional up and down draft values
            wind_int_u(i,k) = wind_int(i,k)
            wind_int_d(i,k) = wind_int(i,k)
            ! provisional tendency
            wind_tend_tmp(i,k) = 0._r8
         end do
      end do

      !-------------------------------------------------------------------------
      ! Calculate pressure perturbation terms

      ! upper boundary - assume mu is zero
      k=1
      pgu(:il2g,k) = 0.0_r8
      pgd(:il2g,k) = 0.0_r8

      ! interior points
      do k=2,pver-1
         km1 = max(1,k-1)
         kp1 = min(pver,k+1)
         do i = il1g,il2g
            mududp(i,k) = ( mu(i,k)  *(wind_mid(i,k)  -wind_mid(i,km1))/dp(i,km1) &
                           +mu(i,kp1)*(wind_mid(i,kp1)-wind_mid(i,k)  )/dp(i,k))
            mddudp(i,k) = ( md(i,k)  *(wind_mid(i,k)  -wind_mid(i,km1))/dp(i,km1) &
                           +md(i,kp1)*(wind_mid(i,kp1)-wind_mid(i,k)  )/dp(i,k))
            pgu(i,k) = -momcu * 0.5_r8 * mududp(i,k)
            pgd(i,k) = -momcd * 0.5_r8 * mddudp(i,k)
         end do
      end do

      ! bottom boundary 
      k = pver
      km1 = max(1,k-1)
      do i=il1g,il2g
         mududp(i,k) = mu(i,k) * (wind_mid(i,k)-wind_mid(i,km1))/dp(i,km1)
         mddudp(i,k) = md(i,k) * (wind_mid(i,k)-wind_mid(i,km1))/dp(i,km1) 
         pgu(i,k) = -momcu * mududp(i,k)
         pgd(i,k) = -momcd * mddudp(i,k)
      end do
       
      !-------------------------------------------------------------------------
      ! Calculate in-cloud velocity

      ! levels adjacent to top and bottom
      k = 2
      km1 = 1
      kk = pver
      kkm1 = max(1,kk-1)
      do i = il1g,il2g
         mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
         if (mupdudp > mbsth) then
            wind_int_u(i,kk) = ( +eu(i,kk)*wind_mid(i,kk)*dp(i,kk) + pgu(i,kk)*dp(i,kk) )/mupdudp
         endif
         if (md(i,k) < -mbsth) then
            wind_int_d(i,k) = ( -ed(i,km1)*wind_mid(i,km1)*dp(i,km1) ) - pgd(i,km1)*dp(i,km1)/md(i,k)
         endif                     
      end do

      ! Updraft from bottom to top
      do kk = pver-1,1,-1
         kkm1 = max(1,kk-1)
         kkp1 = min(pver,kk+1)
         do i = il1g,il2g
            mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
            if (mupdudp > mbsth) then
               wind_int_u(i,kk) = (  mu(i,kkp1)*wind_int_u(i,kkp1) + eu(i,kk)*wind_mid(i,kk)*dp(i,kk) + pgu(i,kk)*dp(i,kk) )/mupdudp
            endif
         end do
      end do

      ! Downdraft from top to bottom
      do k = 3,pver
         km1 = max(1,k-1)
         do i = il1g,il2g
            if (md(i,k) < -mbsth) then
               wind_int_d(i,k) = ( md(i,km1)*wind_int_d(i,km1) - ed(i,km1)*wind_mid(i,km1)*dp(i,km1) - pgd(i,km1)*dp(i,km1) )/md(i,k)
            endif
         end do
      end do

      !-------------------------------------------------------------------------
      ! Calculate momentum tendency

      do k = ktm,pver
         km1 = max(1,k-1)
         kp1 = min(pver,k+1)
         do i = il1g,il2g
            wind_tend_tmp(i,k) = ( mu(i,kp1)* (wind_int_u(i,kp1)-wind_int(i,kp1)) &
                                  -mu(i,k)  * (wind_int_u(i,k)  -wind_int(i,k)  ) &
                                  +md(i,kp1)* (wind_int_d(i,kp1)-wind_int(i,kp1)) &
                                  -md(i,k)  * (wind_int_d(i,k)  -wind_int(i,k)  ) )/dp(i,k)
         end do
      end do

#ifdef CPRCRAY
!DIR$ NOINTERCHANGE
#endif

      ! dcont for bottom layer
      do k = kbm,pver
         km1 = max(1,k-1)
         do i = il1g,il2g
            if (k==mx(i)) then
               wind_tend_tmp(i,k) = (-mu(i,k)*(wind_int_u(i,k)-wind_int(i,k)) &
                                     -md(i,k)*(wind_int_d(i,k)-wind_int(i,k)) )*(1._r8/dp(i,k))
            end if
         end do
      end do

      ! Initialize to zero everywhere, then scatter tendency back to full array
      wind_tend(1:ncol,1:pver,m) = 0._r8

      do k = 1,pver
         do i = il1g,il2g
            wind_tend(ideep(i),k,m) = wind_tend_tmp(i,k)
            ! Output apparent force on the mean flow from pressure gradient
            pguall(ideep(i),k,m) = -pgu(i,k)
            pgdall(ideep(i),k,m) = -pgd(i,k)
            icwu  (ideep(i),k,m) =  wind_int_u(i,k)
            icwd  (ideep(i),k,m) =  wind_int_d(i,k)
         end do
      end do

      ! Calculate momentum flux in units of mb*m/s2 
      do k = ktm,pver
         do i = il1g,il2g
            mflux(i,k,m) = -mu(i,k)*( wind_int_u(i,k) - wind_int(i,k) ) &
                           -md(i,k)*( wind_int_d(i,k) - wind_int(i,k) )
         end do
      end do

      ! Calculate winds at the end of the time step 
      do k = ktm,pver
         do i = il1g,il2g
            km1 = max(1,k-1)
            kp1 = k+1
            windf(i,k,m) = wind_mid(i,k) - ( mflux(i,kp1,m) - mflux(i,k,m) ) * dt/dp(i,k)
         end do
      end do

   end do ! m = 1, nwind

   !----------------------------------------------------------------------------
   ! Need to add an energy fix to account for the dissipation of kinetic energy
   ! Formulation follows from Boville and Bretherton (2003) - modified by Phil Rasch
   do k = ktm,pver
      km1 = max(1,k-1)
      kp1 = min(pver,k+1)
      do i = il1g,il2g
         ! calculate the KE fluxes at top and bot of layer 
         ! based on a discrete approximation to b&b eq(35) F_KE = u*F_u + v*F_v at interface
         utop = ( wind0(i,k  ,1) + wind0(i,km1,1) )/2._r8
         vtop = ( wind0(i,k  ,2) + wind0(i,km1,2) )/2._r8
         ubot = ( wind0(i,kp1,1) + wind0(i,k,  1) )/2._r8
         vbot = ( wind0(i,kp1,2) + wind0(i,k,  2) )/2._r8
         fket = utop*mflux(i,k  ,1) + vtop*mflux(i,k  ,2) ! top of layer
         fkeb = ubot*mflux(i,k+1,1) + vbot*mflux(i,k+1,2) ! bot of layer
         ! divergence of these fluxes should give a conservative redistribution of KE
         ketend_cons = (fket-fkeb)/dp(i,k)
         ! tendency in kinetic energy resulting from the momentum transport
         ketend = ((windf(i,k,1)**2 + windf(i,k,2)**2) - (wind0(i,k,1)**2 + wind0(i,k,2)**2))*0.5_r8/dt
         ! the difference should be the dissipation
         gseten(i,k) = ketend_cons - ketend
      end do
   end do

   ! Scatter dry static energy to full array
   do k = 1,pver
      do i = il1g,il2g
         seten(ideep(i),k) = gseten(i,k)
      end do
   end do

   return

end subroutine zm_transport_momentum

!===================================================================================================

end module zm_transport
