#ifndef ZM_ZM_TRANSPORT_TRACER_IMPL_HPP
#define ZM_ZM_TRANSPORT_TRACER_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_transport_tracer. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_transport_tracer(
  // Inputs
  const MemberType& team,
  const Int& pver,                        // number of mid-point levels
  const uview_1d<const bool>& doconvtran, // flag for doing convective transport
  const uview_2d<const Real>& q,          // tracer array (including water vapor)
  const Int& ncnst,                       // number of tracers to transport
  const uview_1d<const Real>& mu,         // mass flux up
  const uview_1d<const Real>& md,         // mass flux down
  const uview_1d<const Real>& du,         // mass detraining from updraft
  const uview_1d<const Real>& eu,         // mass entraining from updraft
  const uview_1d<const Real>& ed,         // mass entraining from downdraft
  const uview_1d<const Real>& dp,         // delta pressure between interfaces
  const Int& jt,                          // index of cloud top for each column
  const Int& mx,                          // index of cloud top for each column
  const Int& ideep,                       // gathering array
  const Int& il1g,                        // gathered min ncol index
  const Int& il2g,                        // gathered max ncol index
  const uview_2d<const Real>& fracis,     // fraction of tracer that is insoluble
  const uview_1d<const Real>& dpdry,      // delta pressure between interfaces
  const Real& dt,                         // model time increment)
  // Outputs
  const uview_2d<Real>& dqdt)             // output tendency array
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
  /*
       !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i,k                   ! loop indices
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

#ifndef SCREAM_CONFIG_IS_CMAKE
         if (cnst_get_type_byind(m).eq.'dry') then
#else
         if (.false.) then
#endif
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
         if (zm_param%zm_microp) then
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
  */
}

} // namespace zm
} // namespace scream

#endif
