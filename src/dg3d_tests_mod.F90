#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dg3d_tests_mod
!=======================================================================================================!
! R. D. Nair 04/04/06 
!=======================================================================================================!
  ! ------------------------
  use kinds
  ! ------------------------
  use physical_constants
  ! ------------------------
  use dimensions_mod
  ! ------------------------
  use derivative_mod
  ! ------------------------
 !use edge_mod
  ! ------------------------
 !use bndry_mod
  ! ------------------------
  use coordinate_systems_mod
  ! ------------------------
  use quadrature_mod
  ! ------------------------
 !use ref_state_mod
  ! ------------------------
 !use global_norms_mod
  ! ------------------------
  use time_mod
  ! ------------------------
  !use state_mod
  ! ------------------------
 !use hybrid_mod
  ! ------------------------
 !use control_mod     ! Added to support output_prefix 
  ! ------------------------
  use dg3d_vertical_mod, only: eta_levels, eta_all_levels, pt2temp, temp2pt
!
!=======================================================================================================!
!=======================================================================================================!
implicit none
private
!=======================================================================================================!
!	A Baroclinic test by Jablonowski and Williamson (MWR,2005)
!=======================================================================================================!
   real (kind=real_kind), private, parameter  :: u0 = 35.0D0        ! Zonal Mean wind 
   real (kind=real_kind), private, parameter  :: t0 = 288.0D0       ! Mean temp
   real (kind=real_kind), private, parameter  :: p_0 = 100000.0D0   ! Initial Surface pressure 
   real (kind=real_kind), private, parameter  :: gama = 0.005D0     ! Lapse rate
   real (kind=real_kind), private, parameter  :: r_d = 287.04D0     ! Gas const (dry)
   real (kind=real_kind), private, parameter  :: c_p = 1004.64D0    ! Cp        
   real (kind=real_kind), private, parameter  :: grv = 9.80616D0    ! Gravity   
   real (kind=real_kind), private, parameter  :: omg = 7.29212D-05  ! Omega     
   real (kind=real_kind), private, parameter  :: erad = 6.371229D06 ! Earth radius
   real (kind=real_kind), private, parameter  :: ddt = 4.8D05       ! Temp  gradient 
   real (kind=real_kind), private, parameter  :: eta_t  = 0.2D0     ! eta level 
   real (kind=real_kind), private, parameter  :: eta_s  = 1.0D0     ! eta level 
   real (kind=real_kind), private, parameter  :: eta_0  = 0.252D0   ! eta level 
!------------------------------------------------------------------------------------------------
!       Held-Suarez Specific 
!------------------------------------------------------------------------------------------------
   real (kind=real_kind), private, parameter :: secpd   =   24.0D0 * 3600.0D0 
   real (kind=real_kind), private, parameter :: sigma_b = 0.70D0
   real (kind=real_kind), private, parameter :: k_a     = 1.0D0/(40.0D0*secpd)
   real (kind=real_kind), private, parameter :: k_f     = 1.0D0/(1.0D0*secpd)
   real (kind=real_kind), private, parameter :: k_s     = 1.0D0/(4.0D0*secpd)
   real (kind=real_kind), private, parameter :: dT_y    = 60.0D0
   real (kind=real_kind), private, parameter :: dtheta_z= 10.0D0
!=======================================================================================================!  
  public  :: jw_baroclinic
  public  :: heldsuarez_initial
  public  :: heldsuarez_uv_forcing
  public  :: heldsuarez_th_forcing
  public  :: heldsuarez_th_correction
  public  :: heldsuarez_pt_forcing
!=======================================================================================================!
 contains
!=======================================================================================================!
!=======================================================================================================!
subroutine jw_baroclinic(ie,sphere,D,cori,sgp,ptop,tbar,v3d,pt3d,dp3d, qt3d)
!=======================================================================================================!
    integer, intent(in) :: ie
    type (spherical_polar_t), intent(in):: sphere(np,np)
    real (kind=real_kind), intent(in)   :: D(2,2,np,np)
    real (kind=real_kind), intent(out)  :: sgp(np,np),ptop(np,np),cori(np,np)
    real (kind=real_kind), intent(out)  :: tbar(nlev)     
    real (kind=real_kind), intent(out)  :: v3d(np,np,2,nlev)
    real (kind=real_kind), intent(out)  :: pt3d(np,np,nlev),dp3d(np,np,nlev)
    real (kind=real_kind), intent(out)  :: qt3d(np,np,nlev)    !new var for moist

!=======================================================================================================!
    real (kind=real_kind), dimension(np,np,nlev+1):: pr3d, gp3d
    real (kind=real_kind), dimension(np,np,nlev)  :: t3d, dgp3d, dgp 
    real (kind=real_kind) ::  eta(nlev), etv(nlev), etk(nlev+1)
    real (kind=real_kind) ::  ptbar(nlev), ak(nlev+1),bk(nlev+1)
    real (kind=real_kind) ::  v(np,np,2), vv3(np,np,nlev+1)
    real (kind=real_kind) ::  latc, lonc, rc,  aa, lat,lon, snlat, cslat , v1,v2 
    real (kind=real_kind) ::  trm1,trm2,trm3,trm4, term, rdcp, cprd
    integer :: i,j,k
!=======================================================================================================!
! Variables used computing initialization of tracer variable "q3d"

    real (kind=real_kind)  :: p1, p0, z0, h0,  zz, rr, z1 
!=======================================================================================================!


      Call eta_levels(ak,bk,eta,etk)

        rdcp = r_d / c_p
        cprd = c_p / r_d


       do k = 1, nlev
          etv(k) = (eta(k) - eta_0) * DD_PI * 0.5D0 
       end do
       
       latc = DD_PI * (2.0D0/9.0D0) 
       lonc = DD_PI * (1.0D0/9.0D0) 
       

    do k=1, nlev 
       do j=1,np
          do i=1,np

          snlat=SIN(sphere(i,j)%lat)
          cslat=COS(sphere(i,j)%lat)
          lon  =sphere(i,j)%lon
          lat  =sphere(i,j)%lat

          aa = SIN(latc)*snlat + COS(latc)*cslat*COS(lon - lonc) 

          rc =  10.0D0  * ACOS(aa)

           v1 = u0 * (cos(etv(k)))**1.5D0 * (sin(2.0D0 * lat))**2  +  exp(-rc*rc) 

           v2 =  0.0D0

             ! =====================================================
             ! map sphere velocities onto the contravariant cube velocities
             ! =====================================================
       
          v3d(i,j,1,k)= v1*D(1,1,i,j) + v2*D(1,2,i,j)
          v3d(i,j,2,k)= v1*D(2,1,i,j) + v2*D(2,2,i,j)

         !v3d(i,j,1,k) = v1       
         !v3d(i,j,2,k) = v2      

          cori(i,j) = 2.0D0 * omg * snlat

          end do
       end do
    end do

  ! Initial 3D tracer field (ASP: q1)

    lonc = DD_PI / 9.0D0
    latc = DD_PI * (11.0D0/18.0D0) 
    z0 = 0.6D0
    h0 = 0.1D0
    
    do k=1,nlev
       do j=1,np
          do i=1,np
             snlat=SIN(sphere(i,j)%lat)
             cslat=COS(sphere(i,j)%lat)
              lon = sphere(i,j)%lon
              lat = sphere(i,j)%lat

                aa = SIN(latc)*snlat + COS(latc)*cslat*COS(lon - lonc) 
                rc = 10.0D0  * ACOS(aa)

                zz = exp(-(rc*rc + ((eta(k) - z0)/h0)**2))  
        !! ASP-q1     
          !   if (abs(zz) >=  1.0D-08 ) then 
          !       qt3d(i,j,k) = zz  
          !    else 
          !       qt3d(i,j,k) = 0.0D0 
          !   endif 
             
         !! ASP-q4    
                 qt3d(i,j,k) = (tanh(3.0D0*abs(lat) -DD_PI) + 1.0D0)*0.5D0 
          end do
       end do
    end do

 ! Initial Pressure levels (nlev+1)

    do k=1,nlev+1
       do j=1,np
          do i=1,np
            pr3d(i,j,k) = p_0*etk(k)
          end do
       end do
    end do

 ! Pressure thickness  (for dP thickness case)

   do k=1,nlev
      do j=1,np
         do i=1,np
           dp3d(i,j,k) =  pr3d(i,j,k+1) - pr3d(i,j,k)
         end do
      end do
   end do

 ! Top level pressure (constant)

      ptop(:,:) = pr3d(:,:,1)

 ! Layer-mean Temperature fields 

   do k = 1, nlev
    if (eta(k) <= eta_t) then 
       tbar(k) = t0 * eta(k)**(r_d*gama/grv) + ddt * (eta_t - eta(k))**5 
    else 
       tbar(k) = t0 * eta(k)**(r_d*gama/grv)
    endif
       ptbar(k) = tbar(k) * (1.0D0 / eta(k))**(r_d/c_p)
   end do


    do k=1,nlev
       do j=1,np
          do i=1,np
             snlat=SIN(sphere(i,j)%lat)
             cslat=COS(sphere(i,j)%lat)

             trm1 = 0.75D0 * (eta(k) * DD_PI*u0 /r_d) * sin(etv(k)) *sqrt(cos(etv(k)))
             trm2 = -2.0D0 *snlat**6 *(cslat**2 + 1.0D0/3.0D0) + 10.0D0/63.0D0 
             trm3 =  2.0D0 * u0* (cos(etv(k)))**1.5D0
             trm4 = (1.60D0 *cslat**3 *(snlat**2 + 2.0D0/3.0D0) - DD_PI *0.25D0)* erad*omg

             t3d(i,j,k) = tbar(k) + trm1 *(trm2 * trm3 + trm4 )
 
             pt3d(i,j,k) =  t3d(i,j,k) * (1.0D0/eta(k))**(r_d/c_p)
          end do
       end do
    end do
!=======================================================================================================!
! if  GP height thickness is used
!     do k=1,nlev
!        do j=1,np
!           do i=1,np
!             dgp(i,j,k) =   c_p * pt3d(i,j,k)* ((pr3d(i,j,k+1)/p_0)**(rdcp) - &
!                                                    (pr3d(i,j,k)/p_0)**(rdcp))
!             dp3d(i,j,k) = dgp(i,j,k)
!           end do
!        end do
!     end do
!=======================================================================================================!
    !Surface geopotential

    do j=1,np
       do i=1,np
          snlat=SIN(sphere(i,j)%lat)
          cslat=COS(sphere(i,j)%lat)

          trm1 =  u0* ( cos((eta_s - eta_0)*DD_PI*0.5D0) )**1.5D0

          trm2 = -2.0D0 *snlat**6 *(cslat**2 + 1.0D0/3.0D0) + 10.0D0/63.0D0 
          trm3 = (1.60D0 *cslat**3 *(snlat**2 + 2.0D0/3.0D0) - DD_PI *0.25D0)* erad*omg

          sgp(i,j) =   trm1 *(trm2 * trm1 + trm3) 
       end do
    end do

    !=======================================================================================================!
    
 end subroutine jw_baroclinic 





!=======================================================================================================!
!++++ HELD-SUAREZ  Zone! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
!=======================================================================================================!
subroutine heldsuarez_initial(ie,sphere,D,cori,sgp,ptop,tbar,v3d,pt3d,dp3d)
!=======================================================================================================!
    integer, intent(in) :: ie
    type (spherical_polar_t), intent(in):: sphere(np,np)
    real (kind=real_kind), intent(in)   :: D(2,2,np,np)
    real (kind=real_kind), intent(out)  :: sgp(np,np),ptop(np,np),cori(np,np)
    real (kind=real_kind), intent(out)  :: tbar(nlev)     
    real (kind=real_kind), intent(out)  :: v3d(np,np,2,nlev)
    real (kind=real_kind), intent(out)  :: pt3d(np,np,nlev),dp3d(np,np,nlev)
!=======================================================================================================!
    real (kind=real_kind), dimension(np,np,nlev+1):: pr3d, gp3d
    real (kind=real_kind), dimension(np,np,nlev)  :: t3d, dgp3d, dgp 
    real (kind=real_kind) ::  eta(nlev), etv(nlev), etk(nlev+1)
    real (kind=real_kind) ::  ptbar(nlev), ak(nlev+1),bk(nlev+1)
    real (kind=real_kind) ::  v(np,np,2), vv3(np,np,nlev+1), sigm(nlev+1)
    real (kind=real_kind) ::  latc, lonc, rr,rc,  aa, lat,lon, snlat, cslat , v1,v2 
    real (kind=real_kind) ::  trm1,trm2,trm3,trm4, term, rdcp, dsig , sigk, pk
    integer :: i,j,k
!=======================================================================================================!
    Call eta_levels(ak,bk,eta,etk)

    rdcp = r_d/c_p

    sigm(:) = 0.0D0 

    sigm(1) = 5.0D0/p_0
    sigm(nlev+1) = 1.0D0       
    dsig = (sigm(nlev+1) - sigm(1))/dble(nlev)

  !Equi-spaced sigma 
   do k= 1, nlev
     sigm(k+1) =  sigm(k) +  dsig
    ! print*, 'My sigma-------> ', sigm(k) 
   enddo
      
  !Initial state
    do k=1, nlev 
       do j=1,np
          do i=1,np

          v1 = 0.0D0
          v2 = 0.0D0
       
          v3d(i,j,1,k)= v1*D(1,1,i,j) + v2*D(1,2,i,j)
          v3d(i,j,2,k)= v1*D(2,1,i,j) + v2*D(2,2,i,j)

          snlat=SIN(sphere(i,j)%lat)
          cori(i,j)= 2.0D0*omg*snlat

          end do
       end do
    end do

 ! Initial Pressure levels (nlev+1)
    do k=1,nlev+1
       do j=1,np
          do i=1,np
            pr3d(i,j,k)= p_0*etk(k)
           !pr3d(i,j,k)= p_0*sigm(k)
          end do
       end do
    end do

 ! Pressure thickness 

   do k=1,nlev
      do j=1,np
         do i=1,np
           dp3d(i,j,k)=  pr3d(i,j,k+1)-pr3d(i,j,k)
         end do
      end do
   end do

 ! Top level pressure (constant)
   ptop(:,:) = pr3d(:,:,1)

!=======================================================================================================!
! Layer-mean Temperature fields  (Initial Isothermal atmosphere with 300K)
   do k = 1, nlev
       tbar(k) = 300.0D0 
       ptbar(k)= tbar(k) * (1.0D0 / eta(k))**rdcp
   end do

   do k=1,nlev
       sigk = (sigm(k) + sigm(k+1))*0.5D0 
       do j=1,np
          do i=1,np
             pk = (pr3d(i,j,k+1) + pr3d(i,j,k))*0.5D0 
             t3d(i,j,k) = tbar(k) 
             pt3d(i,j,k)= t3d(i,j,k)*(1.0D0/eta(k))**rdcp
            !pt3d(i,j,k)= t3d(i,j,k)*(1.0D0/sigk)**rdcp
            !pt3d(i,j,k)= t3d(i,j,k)*(pk/p_0)**rdcp
          end do
       end do
    end do
!=======================================================================================================!
! if  GP height thickness is used
!    do k=1,nlev
!       do j=1,np
!          do i=1,np
!            dgp(i,j,k) = c_p * pt3d(i,j,k)* ((pr3d(i,j,k+1)/p_0)**(rdcp) - (pr3d(i,j,k)/p_0)**(rdcp))
!            dp3d(i,j,k)= dgp(i,j,k)
!          end do
!       end do
!    end do
!=======================================================================================================!
 !Surface geopotential
    do j=1,np
       do i=1,np
          sgp(i,j)= 0.0D0
       end do
    end do
!=======================================================================================================!
 end subroutine heldsuarez_initial 
!=======================================================================================================!
!    HS Momentum forcings
!=======================================================================================================!
  function heldsuarez_uv_forcing(pr3d,v) result(hs_vfrc)

    real (kind=real_kind), intent(in) :: pr3d(np,np,nlev+1)
    real (kind=real_kind), intent(in) :: v(np,np,2,nlev)
   !real (kind=real_kind), intent(in) :: ps(np,np)
    real (kind=real_kind)             :: hs_vfrc(np,np,2,nlev)
    real (kind=real_kind) ::  etk(nlev+1),ak(nlev+1),bk(nlev+1)
    real (kind=real_kind) ::  eta(nlev),am(nlev),bm(nlev)
    
    ! Local variables

    integer i,j,k
    real (kind=real_kind) :: k_v
    real (kind=real_kind) :: rps(np,np)
    real (kind=real_kind) :: p,peta,sigma

    call eta_all_levels(ak,bk,etk,am,bm,eta)

    do j=1,np
       do i=1,np
         !rps(i,j)= p_0 /ps(i,j)
         rps(i,j) = 1.0D0/pr3d(i,j,nlev+1)
       end do
    end do

    do k=1,nlev
       do j=1,np
          do i=1,np
        !    p    = vert%Amid(k) + vert%Bmid(k)*ps(i,j)
        !    etam = p*rps(i,j)
        !    sigma=  ak(k)*rps(i,j) + bk(k)
        !     sigma= am(k) + bm(k)
        !     peta =  ak(k)*p_0 + bk(k)* pr3d(i,j,nlev+1)
              peta = (pr3d(i,j,k) + pr3d(i,j,k+1)) * 0.5D0 
              sigma=  peta *rps(i,j) 
	    ! sigma= eta(k)
	     k_v = k_f*MAX(0.0D0,(sigma-sigma_b)/(1.0D0-sigma_b))
             hs_vfrc(i,j,1,k) = -k_v*v(i,j,1,k)
             hs_vfrc(i,j,2,k) = -k_v*v(i,j,2,k)
          end do
       end do
    end do

  end function heldsuarez_uv_forcing
!=======================================================================================================!
!    HS Temparature forcings
!=======================================================================================================!
  function heldsuarez_th_forcing(dt,sphere,pr3d,t3d) result(th_force)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pr3d(np,np,nlev+1)
    real (kind=real_kind), intent(in) :: t3d(np,np,nlev)

  ! Local variables
  ! real (kind=real_kind), intent(out) :: k_th(np,np,nlev), t_eq(np,np,nlev)
    real (kind=real_kind) :: k_th(np,np,nlev), t_eq(np,np,nlev)
    real (kind=real_kind) :: etk(nlev+1),ak(nlev+1),bk(nlev+1)
    real (kind=real_kind) :: eta(nlev),am(nlev),bm(nlev)
    real (kind=real_kind) :: th_force(np,np,nlev)
    real (kind=real_kind) :: k_t(np,np)
    real (kind=real_kind) :: snlatsq(np,np)
    real (kind=real_kind) :: cslatsq(np,np)
    real (kind=real_kind) :: rps(np,np),ps(np,np)

    real (kind=real_kind) :: p,logprat,pratk,Teq
    real (kind=real_kind) :: logps0,etam,temp,tnew
    real (kind=real_kind) :: lat,snlat,kapa,peta,sigma 
    real (kind=real_kind) :: rec_one_minus_sigma_b,den,num
    integer i,j,k
!==============================================================================================================
    call eta_all_levels(ak,bk,etk,am,bm,eta)

    logps0= LOG(p_0)
    kapa  = 2.0D0/7.0D0
    ps(:,:)=  pr3d(:,:,nlev+1)

    do j=1,np
       do i=1,np
         snlat       = SIN(sphere(i,j)%lat)
         snlatsq(i,j)= snlat*snlat
         cslatsq(i,j)= 1.0D0 - snlatsq(i,j)
         rps(i,j)    = 1.0D0/ps(i,j)
       end do
    end do

    rec_one_minus_sigma_b = 1.0D0/(1.0D0 - sigma_b)

    do k=1,nlev
       do j=1,np
          do i=1,np

          temp= t3d(i,j,k)             
          p = am(k)*p_0 + bm(k)*ps(i,j)
             logprat= LOG(p)-logps0
             pratk  = EXP(kappa*(logprat))
          sigma  = am(k)+ bm(k)
             k_th(i,j,k)= k_a + (k_s-k_a)*cslatsq(i,j)*cslatsq(i,j)* &
                           MAX(0.0D0,(sigma - sigma_b)/(1.0D0 - sigma_b))
            t_eq(i,j,k)= MAX(200.0D0,(315.0D0 - dT_y*snlatsq(i,j) - dtheta_z*logprat*cslatsq(i,j))*pratk)
    !======================================================================================
            th_force(i,j,k)= -k_th(i,j,k)*(temp - t_eq(i,j,k))
    !======================================================================================
         ! num= temp + dt*k_th(i,j,k)*t_eq(i,j,k)
         ! den= 1 + dt*k_th(i,j,k)
         ! th_force(i,j,k)= num/den
  
          end do
       end do
    end do

  end function heldsuarez_th_forcing 
!==============================================================================================================
!==============================================================================================================
function heldsuarez_th_correction(dt,oldt3d,sphere,pr3d,t3d) result(pt3d)

    type (spherical_polar_t), intent(in):: sphere(np,np)
    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pr3d(np,np,nlev+1)
    real (kind=real_kind), intent(in) :: t3d(np,np,nlev), oldt3d(np,np,nlev)
    real (kind=real_kind) :: pt3d(np,np,nlev)

  ! Local variables
    real (kind=real_kind) :: t3dnew(np,np,nlev)
    real (kind=real_kind) :: k_th(np,np,nlev),t_eq(np,np,nlev)
    real (kind=real_kind) :: etk(nlev+1),ak(nlev+1),bk(nlev+1)
    real (kind=real_kind) :: eta(nlev),am(nlev),bm(nlev)
    real (kind=real_kind) :: k_t(np,np)
    real (kind=real_kind) :: snlatsq(np,np)
    real (kind=real_kind) :: cslatsq(np,np)
    real (kind=real_kind) :: rps(np,np),ps(np,np)

    real (kind=real_kind) :: p,logprat,pratk,Teq
    real (kind=real_kind) :: logps0,etam,temp,tnew
    real (kind=real_kind) :: lat,snlat,kapa,peta,sigma 
    real (kind=real_kind) :: rec_one_minus_sigma_b,den,num
    integer i,j,k
!==============================================================================================================
    call eta_all_levels(ak,bk,etk,am,bm,eta)

    logps0= LOG(p_0)
    kapa  = 2.0D0/7.0D0
    ps(:,:)= pr3d(:,:,nlev+1)

    do j=1,np
       do i=1,np
         snlat       = SIN(sphere(i,j)%lat)
         snlatsq(i,j)= snlat*snlat
         cslatsq(i,j)= 1.0D0 - snlatsq(i,j)
         rps(i,j)    = 1.0D0/ps(i,j)
       end do
    end do

    do k=1,nlev 
       ! pott(:,:) = oldt3d(:,:,k)
       do j=1,np
          do i=1,np
	     !p = am(k)*p_0 + bm(k)*ps(i,j)
	     !logprat= LOG(p)-logps0
             !pratk  = EXP(kappa*(logprat))
	     !sigma  = am(k) + bm(k)
	     !sigma  = eta(k)
	     !sigma  = p/p_0

	     !temp= t3d(i,j,k)
	     temp= oldt3d(i,j,k)
             p= (pr3d(i,j,k) + pr3d(i,j,k+1))*0.5D0
	     !sigma= eta(k)
	     sigma= p / pr3d(i,j,nlev+1)
	     pratk= (p/p_0)**kapa	     
	     logprat= LOG(p/p_0)
	     k_th(i,j,k)= k_a + (k_s-k_a)*cslatsq(i,j)*cslatsq(i,j)* &
                          MAX(0.0D0,(sigma - sigma_b)/(1.0D0 - sigma_b))
             !t_eq(i,j,k)= MAX(200.0D0,(315.0D0 - dT_y*snlatsq(i,j) - dtheta_z*logprat*cslatsq(i,j))*pratk)
             t_eq(i,j,k)= MAX(200.0D0/pratk,(315.0D0 - dT_y*snlatsq(i,j) - dtheta_z*logprat*cslatsq(i,j)))
	    !======================================================================================
              num= temp  + dt*k_th(i,j,k)*t_eq(i,j,k) *0.5
	      den= 1.0D0 + dt*k_th(i,j,k)  *0.5
	      pt3d(i,j,k)= num/den     
	    ! t3dnew(i,j,k)= num/den
	    ! t3dnew(i,j,k)=  temp*(1.0D0 - k_th(i,j,k) * dt)  + dt*k_th(i,j,k)*t_eq(i,j,k) 
  
          end do
       end do
    end do
   
  ! pt3d(:,:,:)=  temp2pt(pr3d,t3dnew)
    
end function heldsuarez_th_correction
!==============================================================================================================
!==============================================================================================================
function heldsuarez_pt_forcing(dt,sphere,pr3d,t3d) result(pt3d_force)

    type (spherical_polar_t), intent(in):: sphere(np,np)
    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pr3d(np,np,nlev+1)
    real (kind=real_kind), intent(in) :: t3d(np,np,nlev)
    real (kind=real_kind) :: pt3d_force(np,np,nlev)

  ! Local variables
    real (kind=real_kind) :: t3dnew(np,np,nlev)
    real (kind=real_kind) :: k_th(np,np,nlev),t_eq(np,np,nlev)
    real (kind=real_kind) :: etk(nlev+1),ak(nlev+1),bk(nlev+1)
    real (kind=real_kind) :: eta(nlev),am(nlev),bm(nlev)
    real (kind=real_kind) :: k_t(np,np)
    real (kind=real_kind) :: snlatsq(np,np)
    real (kind=real_kind) :: cslatsq(np,np)
    real (kind=real_kind) :: rps(np,np),ps(np,np)

    real (kind=real_kind) :: p,logprat,pratk,Teq
    real (kind=real_kind) :: logps0,etam,temp,tnew
    real (kind=real_kind) :: lat,snlat,kapa,peta,sigma 
    real (kind=real_kind) :: rec_one_minus_sigma_b,den,num
    integer i,j,k
!==============================================================================================================
    call eta_all_levels(ak,bk,etk,am,bm,eta)

    logps0= LOG(p_0)
    kapa  = 2.0D0/7.0D0
    ps(:,:)= pr3d(:,:,nlev+1)

    do j=1,np
       do i=1,np
         snlat       = SIN(sphere(i,j)%lat)
         snlatsq(i,j)= snlat*snlat
         cslatsq(i,j)= 1.0D0 - snlatsq(i,j)
         rps(i,j)    = 1.0D0/ps(i,j)
       end do
    end do

    rec_one_minus_sigma_b = 1.0D0/(1.0D0 - sigma_b)

  !======================================================================================
    do k=1,nlev
       do j=1,np
          do i=1,np
           temp= t3d(i,j,k)
             p= (pr3d(i,j,k) + pr3d(i,j,k+1))*0.5D0
            pratk= (p/p_0)**kapa
            logprat= LOG(p/p_0)
            sigma  = eta(k)
            k_th(i,j,k)= k_a + (k_s-k_a)*cslatsq(i,j)*cslatsq(i,j)* &
                          MAX(0.0D0,(sigma - sigma_b)/(1.0D0 - sigma_b))
            t_eq(i,j,k)= MAX(200.0D0,(315.0D0 - dT_y*snlatsq(i,j) - dtheta_z*logprat*cslatsq(i,j))*pratk)
            !t3dnew(i,j,k)= k_th(i,j,k)*(t3d(i,j,k)-t_eq(i,j,k))
            t3dnew(i,j,k)= t3d(i,j,k) - t_eq(i,j,k)
  
          end do
       end do
    end do
   
    pt3d_force(:,:,:)=  temp2pt(pr3d,t3dnew)
    
    do k=1,nlev
       do j=1,np
          do i=1,np
             pt3d_force(i,j,k) = -k_th(i,j,k) * pt3d_force(i,j,k)  
          end do
       end do
    end do

    return
end function heldsuarez_pt_forcing
!==============================================================================================================
end module dg3d_tests_mod
