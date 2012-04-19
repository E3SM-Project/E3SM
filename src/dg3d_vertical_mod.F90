#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dg3d_vertical_mod
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
   real (kind=real_kind), private, parameter  :: gama= 0.005D0     ! Lapse rate
   real (kind=real_kind), private, parameter  :: r_d = 287.04D0      ! Gas const (dry)
   real (kind=real_kind), private, parameter  :: c_p = 1004.64D0    ! Cp        
   real (kind=real_kind), private, parameter  :: grv = 9.80616D0    ! Gravity   
   real (kind=real_kind), private, parameter  :: omg = 7.29212D-05  ! Omega     
   real (kind=real_kind), private, parameter  :: erad = 6.371229D06 ! Earth radius
   real (kind=real_kind), private, parameter  :: ddt = 4.8D05       ! Temp  gradient 
   real (kind=real_kind), private, parameter  :: eta_t  = 0.2D0     ! eta level 
   real (kind=real_kind), private, parameter  :: eta_s  = 1.0D0     ! eta level 
   real (kind=real_kind), private, parameter  :: eta_0  = 0.252D0   ! eta level 
   Integer, private, parameter :: mono=2
!=======================================================================================================! 
  public  :: pt2temp
  public  :: temp2pt 
  public  :: lagrangian_surfvars
  public  :: etalevel_prs_gp
  public  :: etalevel_temp
  public  :: eta_levels  
  public  :: eta_all_levels
  public  :: intp_cubiclag
  public  :: linear_int
  public  :: index_search

  private  :: fluxjacobian_bcl
!=======================================================================================================!
 contains
!=======================================================================================================!
  subroutine lagrangian_surfvars(ie,sgp,ptop,pt3d,dp3d,pr3d,gp3d,ht3d,peta,t3d)
!=======================================================================================================!
    integer, intent(in)                  :: ie
    real (kind=real_kind), intent(in)    :: sgp(np,np), ptop(np,np) 
    real (kind=real_kind), intent(in)    :: pt3d(np,np,nlev),dp3d(np,np,nlev)
    real (kind=real_kind), intent(out)   :: pr3d(np,np,nlev+1),gp3d(np,np,nlev+1)
    real (kind=real_kind), intent(out)   :: ht3d(np,np,nlev),peta(np,np,nlev), t3d(np,np,nlev)
!=======================================================================================================!
    real (kind=real_kind), dimension(np,np,nlev)   :: dgp3d, dlnp
    real (kind=real_kind), dimension(np,np,nlev+1) :: lnpr
    real (kind=real_kind), dimension(np,np,nlev+1)   :: pprs
    real (kind=real_kind) ::  eta(nlev), ak(nlev+1), bk(nlev+1), etk(nlev+1)
    real (kind=real_kind) ::  xx(4), yy(4)          
    real (kind=real_kind) ::  rdcp,p1,p2,pp,cprd,pp_0,dpr,alphak,pk
    integer :: i,j,k 
!=======================================================================================================!
    rdcp = r_d / c_p
    cprd = 1.0D0 / rdcp 
    pp_0 = 1.0D0 / (p_0)**rdcp

 !Note: dp3d --> GP thickness (for GP-thickness case otherwise pressure thickness)

! Top and bottom boundary 
    pr3d(:,:,1)= ptop(:,:)
    gp3d(:,:,nlev+1)= sgp(:,:)

!Pressure at etk(k) levels (interfaces)   for dP thickness case
    do k=1,nlev
       do j=1,np
          do i=1,np
            pr3d(i,j,k+1) = pr3d(i,j,k) + dp3d(i,j,k)
          end do
       end do
    end do

   do k=1, nlev+1
      do j=1,np
         do i=1,np
           pprs(i,j,k) = (pr3d(i,j,k))**rdcp 
         end do
      end do
   end do

! !Computed geopotential thickness [d(phi) = -Cp*Theta_0*d(P^kappa)  (dP-Thickness case)
   do k=nlev, 1, -1
      do j=1,np
         do i=1,np
             dgp3d(i,j,k) =   c_p * pp_0 * pt3d(i,j,k)* (pprs(i,j,k+1) - pprs(i,j,k))
         end do
      end do
   end do

!Geopotential for all levels (nlev+1), From Bottom (eta=1) to top (eta=0) (dP-thickness)
     do k=nlev, 1, -1
        do j=1,np
           do i=1,np
              gp3d(i,j,k) =   gp3d(i,j,k+1) + dgp3d(i,j,k) 
           end do
        end do
    end do
 
! In any case, dP or GP ...
   do k=1,nlev+1
        do j=1,np
           do i=1,np
               lnpr(i,j,k)= log(pr3d(i,j,k))
           end do
        end do
   end do

! Derivation of temperature and height at the "full"  levels 
!!  Based on Simmons & Burridge (MWR, 1981) + Peter's thesis
!!for GP thckness:
!
!    do k=1,nlev
!       do j=1,np
!          do i=1,np
!
!             dpr = pr3d(i,j,k+1) - pr3d(i,j,k)
!             alphak =  1.0D0 - (pr3d(i,j,k)/dpr)* (lnpr(i,j,k+1)-lnpr(i,j,k))
!             if (k == 1)  alphak = log(2.0D0)
!
!             pk = pr3d(i,j,k+1)*exp(-alphak)
!             t3d(i,j,k) =  pt3d(i,j,k)* (pk / p_0)**rdcp      !regular temp
!
!             ht3d(i,j,k) = gp3d(i,j,k+1) +  r_d* t3d(i,j,k) *   alphak 
!             peta(i,j,k) =  log(pk)
!
!             !ht3d(i,j,k) = gp3d(i,j,k+1) +  r_d* t3d(i,j,k) *    &
!             !                     (1.0D0 - (pr3d(i,j,k)/dpr)* (lnpr(i,j,k+1)-lnpr(i,j,k))  )
!             !peta(i,j,k) =  (pr3d(i,j,k+1)*lnpr(i,j,k+1) - pr3d(i,j,k)*lnpr(i,j,k))/dpr 
!             ! peta(i,j,k) =  log(pp)
!             ! peta(i,j,k) =  lnpr(i,j,k+1) - alphak
!               
!          end do
!       end do
!    end do

!! dP -case 
    do k=1,nlev
       do j=1,np
          do i=1,np
           pp = (pr3d(i,j,k) + pr3d(i,j,k+1))*0.5D0
               t3d(i,j,k)=  pt3d(i,j,k)* (pp/p_0)**rdcp
               dpr = pr3d(i,j,k+1) - pr3d(i,j,k)
               ht3d(i,j,k)= gp3d(i,j,k+1)+ r_d*t3d(i,j,k)*(1.0D0-(pr3d(i,j,k)/dpr)*(lnpr(i,j,k+1)-lnpr(i,j,k)))
                peta(i,j,k)=  log(pp)
             !  peta(i,j,k)=  (lnpr(i,j,k) + lnpr(i,j,k+1))*0.5D0
            end do
         end do
      end do
!=======================================================================================================!
 end subroutine lagrangian_surfvars 
!=======================================================================================================! 
function pt2temp(pr3d,pt3d) result(t3d)
 real (kind=real_kind),intent(in) :: pr3d(np,np,nlev),pt3d(np,np,nlev)
 real (kind=real_kind)            :: t3d(np,np,nlev)
 integer:: i,j,k    
 real (kind=real_kind):: rdcp,pp
 
 rdcp = r_d/c_p
 do k=1,nlev
 do j=1,np
 do i=1,np
    pp = (pr3d(i,j,k) + pr3d(i,j,k+1))*0.5D0
    t3d(i,j,k)=  pt3d(i,j,k)*(pp/p_0)**rdcp 
 enddo
 enddo
 enddo
end function pt2temp 
!=======================================================================================================! 
function temp2pt(pr3d,t3d) result(pt3d)
 real (kind=real_kind),intent(in) :: pr3d(np,np,nlev),t3d(np,np,nlev)
 real (kind=real_kind)            :: pt3d(np,np,nlev)
 integer:: i,j,k    
 real (kind=real_kind):: rdcp,pp
 
 rdcp = r_d/c_p
 do k=1,nlev
 do j=1,np
 do i=1,np
    pp = (pr3d(i,j,k) + pr3d(i,j,k+1))*0.5D0
    pt3d(i,j,k)=  t3d(i,j,k)*(p_0/pp)**rdcp 
 enddo
 enddo
 enddo
end function temp2pt 
!=======================================================================================================!
 subroutine lagrangian_surfvars_temp(ie,sgp,ptop,pt3d,dp3d,pr3d,gp3d,ht3d,peta,t3d)
!=======================================================================================================!
    integer, intent(in)                  :: ie
    real (kind=real_kind), intent(in)    :: sgp(np,np), ptop(np,np) 
    real (kind=real_kind), intent(in)    :: pt3d(np,np,nlev),dp3d(np,np,nlev)
    real (kind=real_kind), intent(out)   :: pr3d(np,np,nlev+1),gp3d(np,np,nlev+1)
    real (kind=real_kind), intent(out)   :: ht3d(np,np,nlev),peta(np,np,nlev), t3d(np,np,nlev)
!=======================================================================================================!
    real (kind=real_kind), dimension(np,np,nlev)   :: dgp3d, dlnp
    real (kind=real_kind), dimension(np,np,nlev+1) :: lnpr
    real (kind=real_kind), dimension(np,np,nlev)   :: pprs
    real (kind=real_kind) ::  eta(nlev), ak(nlev+1), bk(nlev+1), etk(nlev+1)
    real (kind=real_kind) ::  xx(4), yy(4)          
    real (kind=real_kind) ::  rdcp, p1,p2,pp, cprd, pp_0, dpr
    integer :: i,j,k 
!=======================================================================================================!
    rdcp = r_d / c_p
    cprd = 1.0D0 / rdcp 
    pp_0 = 1.0D0 / (p_0)**rdcp

    pr3d(:,:,1) = ptop(:,:)

 !Computed pressure at interfaces from GP thickness (for GP-thickness case)
 !Note: dp3d --> GP thickness
    do k=1, nlev
       do j=1,np
          do i=1,np
             pr3d(i,j,k+1) = p_0 * ( (pr3d(i,j,k)/p_0)**rdcp  + dp3d(i,j,k)/(c_p* pt3d(i,j,k)) )**cprd
          end do
       end do
    end do


 !Pressure at etk(k) levels (interfaces)   for dP thickness case
  !do k=1,nlev
  !   do j=1,np
  !      do i=1,np
  !        pr3d(i,j,k+1) = pr3d(i,j,k) + dp3d(i,j,k)
  !      end do
  !   end do
  !end do

  !do k=1, nlev+1
  !   do j=1,np
  !      do i=1,np
  !        pprs(i,j,k) = (pr3d(i,j,k))**rdcp 
  !      end do
  !   end do
  !end do

 !Computed geopotential thickness [d(phi) = -Cp*Theta_0*d(P^kappa)

  !do k=nlev, 1, -1
  !   do j=1,np
  !      do i=1,np
  !          dgp3d(i,j,k) =   c_p * pp_0 * pt3d(i,j,k)* (pprs(i,j,k+1) - pprs(i,j,k))
  !      end do
  !   end do
  !end do

 !Geopotential for all levels (nlev+1), From Bottom (eta=1) to top (eta=0)

    gp3d(:,:,nlev+1) = sgp(:,:)

    do k=nlev, 1, -1
       do j=1,np
          do i=1,np
           ! gp3d(i,j,k) =   gp3d(i,j,k+1) + dgp3d(i,j,k) 
             gp3d(i,j,k) =   gp3d(i,j,k+1) + dp3d(i,j,k)    !for GP-thickness case
          end do
       end do
    end do
 
    do k=1,nlev+1
        do j=1,np
           do i=1,np
               lnpr(i,j,k) = log(pr3d(i,j,k))
           end do
        end do
    end do

!!  Based on Simmons & Burridge (MWR, 1981)

    do k=1,nlev
       do j=1,np
          do i=1,np

                     pp = (pr3d(i,j,k) + pr3d(i,j,k+1))*0.5D0

             t3d(i,j,k) =  pt3d(i,j,k)* (pp / p_0)**rdcp      !regular temp

          ! for dP thickness case only
          ! ht3d(i,j,k) = gp3d(i,j,k+1) +  r_d* t3d(i,j,k)*  &
          !             (1.0D0 - (pr3d(i,j,k)/dp3d(i,j,k))* log(pr3d(i,j,k+1)/pr3d(i,j,k))  )
       
       
          !for GP thckness:

             dpr = pr3d(i,j,k+1) - pr3d(i,j,k)
             ht3d(i,j,k) = gp3d(i,j,k+1) +  r_d* t3d(i,j,k) *    &
                                  (1.0D0 - (pr3d(i,j,k)/dpr)* (lnpr(i,j,k+1)-lnpr(i,j,k))  )

            !peta(i,j,k) =  (pr3d(i,j,k+1)*lnpr(i,j,k+1) - pr3d(i,j,k)*lnpr(i,j,k))/dpr 
             peta(i,j,k) =  log(pp)
          end do
       end do
    end do
!=======================================================================================================!
 end subroutine lagrangian_surfvars_temp 
!=======================================================================================================!
!=======================================================================================================!
  subroutine etalevel_prs_gp(ie,sgp,ptop,pt3d,dp3d,pr3d,gp3d,ht3d,peta,t3d)
!=======================================================================================================!
    integer, intent(in)                  :: ie
    real (kind=real_kind), intent(in)    :: sgp(np,np), ptop(np,np) 
    real (kind=real_kind), intent(in)    :: pt3d(np,np,nlev),dp3d(np,np,nlev)
    real (kind=real_kind), intent(out)   :: pr3d(np,np,nlev+1),gp3d(np,np,nlev+1)
    real (kind=real_kind), intent(out)   :: ht3d(np,np,nlev),peta(np,np,nlev), t3d(np,np,nlev)
!=======================================================================================================!
    real (kind=real_kind), dimension(np,np,nlev)   :: dgp3d, dlnp
    real (kind=real_kind), dimension(np,np,nlev+1) :: lnpr
    real (kind=real_kind), dimension(np,np,nlev)   :: pprs
    real (kind=real_kind) ::  eta(nlev), ak(nlev+1), bk(nlev+1), etk(nlev+1)
    real (kind=real_kind) ::  xx(4), yy(4)          
    real (kind=real_kind) ::  rdcp, p1,p2,pp, cprd, pp_0, dpr
    integer :: i,j,k 
!=======================================================================================================!
    rdcp = r_d / c_p
    cprd = 1.0D0 / rdcp 
    pp_0 = 1.0D0 / (p_0)**rdcp

  !!  Call eta_levels(ak,bk,eta,etk)
    pr3d(:,:,1) = ptop(:,:)

 !Computed pressure at interfaces from GP thickness (for GP-thickness case)
 !Note: dp3d --> GP thickness
    do k=1, nlev
       do j=1,np
          do i=1,np
             pr3d(i,j,k+1) = p_0 * ( (pr3d(i,j,k)/p_0)**rdcp  + dp3d(i,j,k)/(c_p* pt3d(i,j,k)) )**cprd
          end do
       end do
    end do


 !Pressure at etk(k) levels (interfaces)   for dP thickness case
  !do k=1,nlev
  !   do j=1,np
  !      do i=1,np
  !        pr3d(i,j,k+1) = pr3d(i,j,k) + dp3d(i,j,k)
  !      end do
  !   end do
  !end do

  !do k=1, nlev+1
  !   do j=1,np
  !      do i=1,np
  !        pprs(i,j,k) = (pr3d(i,j,k))**rdcp 
  !      end do
  !   end do
  !end do

 !Computed geopotential thickness [d(phi) = -Cp*Theta_0*d(P^kappa)

  !do k=nlev, 1, -1
  !   do j=1,np
  !      do i=1,np
  !          dgp3d(i,j,k) =   c_p * pp_0 * pt3d(i,j,k)* (pprs(i,j,k+1) - pprs(i,j,k))
  !      end do
  !   end do
  !end do

 !Geopotential for all levels (nlev+1), From Bottom (eta=1) to top (eta=0)

    gp3d(:,:,nlev+1) = sgp(:,:)

    do k=nlev, 1, -1
       do j=1,np
          do i=1,np
           ! gp3d(i,j,k) =   gp3d(i,j,k+1) + dgp3d(i,j,k) 
             gp3d(i,j,k) =   gp3d(i,j,k+1) + dp3d(i,j,k)    !for GP-thickness case
          end do
       end do
    end do
 
    do k=1,nlev+1
        do j=1,np
           do i=1,np
               lnpr(i,j,k) = log(pr3d(i,j,k))
           end do
        end do
    end do

!!  Based on Simmons & Burridge (MWR, 1981)

    do k=1,nlev
       do j=1,np
          do i=1,np

                     pp = (pr3d(i,j,k) + pr3d(i,j,k+1))*0.5D0

             t3d(i,j,k) =  pt3d(i,j,k)* (pp / p_0)**rdcp      !regular temp

          ! for dP thickness case only
          ! ht3d(i,j,k) = gp3d(i,j,k+1) +  r_d* t3d(i,j,k)*  &
          !             (1.0D0 - (pr3d(i,j,k)/dp3d(i,j,k))* log(pr3d(i,j,k+1)/pr3d(i,j,k))  )
       
          !for GP thckness:

             dpr = pr3d(i,j,k+1) - pr3d(i,j,k)
             ht3d(i,j,k) = gp3d(i,j,k+1) +  r_d* t3d(i,j,k) *    &
                                  (1.0D0 - (pr3d(i,j,k)/dpr)* (lnpr(i,j,k+1)-lnpr(i,j,k))  )

            !peta(i,j,k) =  (pr3d(i,j,k+1)*lnpr(i,j,k+1) - pr3d(i,j,k)*lnpr(i,j,k))/dpr 
             peta(i,j,k) =  log(pp)
          end do
       end do
    end do
!=======================================================================================================!
 end subroutine etalevel_prs_gp 
!=======================================================================================================!

!=======================================================================================================!
subroutine etalevel_temp(ps,pt3d,tt3d)
!=======================================================================================================!

    real (kind=real_kind), dimension(np,np), intent(in)         :: ps 
    real (kind=real_kind), dimension(np,np,nlev),intent(in)     :: pt3d
    real (kind=real_kind), dimension(np,np,nlev),intent(out)    :: tt3d

    real (kind=real_kind) ::  eta(nlev), ak(nlev+1), bk(nlev+1), etk(nlev+1)
    real (kind=real_kind) ::  rdcp, peta 

    integer :: i,j,k 

           rdcp = r_d / c_p

      Call eta_levels(ak,bk,eta,etk)

 !Temparature  at eta(k) levels  from Pot.temp

    do k=1,nlev
       do j=1,np
          do i=1,np
            !peta = ((ak(k) + ak(k+1))*p_0 + (bk(k) + bk(k+1))* ps(i,j) ) *0.5D0
            !tt3d(i,j,k+1) = pt3d(i,j,k) * (peta / p_0)**rdcp
             peta =  ak(k)*p_0 + bk(k)* ps(i,j) 
             tt3d(i,j,k) = pt3d(i,j,k) * (peta / p_0)**rdcp
          end do
       end do
    end do
!=======================================================================================================!
 end subroutine etalevel_temp   
!=======================================================================================================!
!=======================================================================================================!
subroutine eta_levels(ak,bk,eta,etk)

    real (kind=real_kind), dimension(nlev+1), intent(out) ::  ak, bk, etk
    real (kind=real_kind), dimension(nlev),   intent(out) ::  eta


    real (kind=real_kind) :: a18(18+1) = &
    (/ 0.00251499D0,  0.00710361D0,  0.01904260D0, 0.04607560D0, 0.08181860D0, &
       0.07869805D0,  0.07463175D0,  0.06955308D0, 0.06339061D0, 0.05621774D0, &
       0.04815296D0,  0.03949230D0,  0.03058456D0, 0.02193336D0, 0.01403670D0, &
       0.007458598D0, 0.002646866D0, 0.00000000D0, 0.00000000D0  /)

    real (kind=real_kind) :: b18(18+1) = &
    (/ 0.000000D0,   0.000000D0,   0.000000D0,  0.000000D0,  0.000000D0, &
       0.03756984D0, 0.08652625D0, 0.1476709D0, 0.221864D0,  0.308222D0, &
       0.4053179D0,  0.509588D0,   0.6168328D0, 0.7209891D0, 0.816061D0, &
       0.8952581D0,  0.953189D0,   0.985056D0,  1.000000D0   /)
    

    real (kind=real_kind) :: a26(26+1) = &
    (/ 0.002194067D0, 0.004895209D0, 0.009882418D0, 0.01805201D0, 0.02983724D0, &
       0.04462334D0,  0.06160587D0,  0.07851243D0,  0.07731271D0, 0.07590131D0, &
       0.07424086D0,  0.07228744D0,  0.06998933D0,  0.06728574D0, 0.06410509D0, &
       0.06036322D0,  0.05596111D0,  0.05078225D0,  0.04468960D0, 0.03752191D0, &
       0.02908949D0,  0.02084739D0,  0.01334443D0,  0.00708499D0, 0.00252136D0, &
       0.00000000D0,  0.00000000D0  /)
  
    real (kind=real_kind) :: b26(26+1) = &
    (/ 0.00000000D0, 0.00000000D0, 0.00000000D0, 0.00000000D0, 0.00000000D0, &
        0.00000000D0, 0.00000000D0, 0.00000000D0, 0.01505309D0, 0.03276228D0, &
        0.05359622D0, 0.07810627D0, 0.1069411D0,  0.1408637D0,  0.1807720D0,  &
        0.2277220D0,  0.2829562D0,  0.3479364D0,  0.4243822D0,  0.5143168D0,  &
        0.6201202D0,  0.7235355D0,  0.8176768D0,  0.8962153D0,  0.9534761D0,  &
        0.9851122D0,  1.0000000D0  /)

    real (kind=real_kind) :: a49(49+1) = &
    (/  0.0022518650D0, 0.0039838900D0, 0.0067043640D0, 0.0107323100D0, 0.0163423300D0, &
	0.0236711900D0, 0.0326145600D0, 0.0427452700D0, 0.0538261000D0, 0.0651217500D0, &
	0.0756985000D0, 0.0845428300D0, 0.0839631000D0, 0.0833410300D0, 0.0826735200D0, &
	0.0819572500D0, 0.0811886600D0, 0.0803639300D0, 0.0794789500D0, 0.0785293400D0, &
	0.0775103600D0, 0.0764169500D0, 0.0752436800D0, 0.0739847000D0, 0.0726337500D0, &
	0.0711841400D0, 0.0696286300D0, 0.0679595000D0, 0.0661684600D0, 0.0642465800D0, &
	0.0621843300D0, 0.0599714400D0, 0.0575969000D0, 0.0550489200D0, 0.0523148300D0, &
	0.0493810200D0, 0.0462329200D0, 0.0428548700D0, 0.0392300600D0, 0.0353404900D0, &
	0.0311668100D0, 0.0266882500D0, 0.0218825700D0, 0.0167637100D0, 0.0120817100D0, &
	0.0079596120D0, 0.0045102970D0, 0.0018312150D0, 0.0000000000D0, 0.0000000000D0 /)

    real (kind=real_kind) :: b49(49+1) = &
    (/  0.0000000000D0, 0.0000000000D0, 0.0000000000D0, 0.0000000000D0, 0.0000000000D0, &
	0.0000000000D0, 0.0000000000D0, 0.0000000000D0, 0.0000000000D0, 0.0000000000D0, &
	0.0000000000D0, 0.0000000000D0, 0.0067551120D0, 0.0140036400D0, 0.0217816400D0, &
	0.0301277800D0, 0.0390835600D0, 0.0486935200D0, 0.0590054200D0, 0.0700705600D0, &
	0.0819439400D0, 0.0946845900D0, 0.1083559000D0, 0.1230258000D0, 0.1387673000D0, &
	0.1556586000D0, 0.1737837000D0, 0.1932327000D0, 0.2141024000D0, 0.2364965000D0, &
	0.2605264000D0, 0.2863115000D0, 0.3139801000D0, 0.3436697000D0, 0.3755280000D0, &
	0.4097133000D0, 0.4463958000D0, 0.4857576000D0, 0.5279946000D0, 0.5733168000D0, &
	0.6219495000D0, 0.6741346000D0, 0.7301315000D0, 0.7897776000D0, 0.8443334000D0, &
	0.8923650000D0, 0.9325572000D0, 0.9637744000D0, 0.9851122000D0, 1.0000000000D0 /)
    
    integer :: k
!=======================================================================================================!
    if (nlev == 18 ) then

       do k = 1, nlev+1
          ak(k)= a18(k)
          bk(k)= b18(k)
          etk(k) = ak(k) + bk(k)
       end do

    elseif (nlev == 26) then

       do k = 1, nlev+1
          ak(k)= a26(k)
          bk(k)= b26(k)
          etk(k) = ak(k) + bk(k)
       end do

    elseif (nlev == 49) then

       do k = 1, nlev+1
          ak(k)= a49(k)
          bk(k)= b49(k)
          etk(k) = ak(k) + bk(k)
       end do

    endif 

       do k = 1, nlev
          eta(k) = (etk(k) + etk(k+1)) * 0.5D0 
       end do

 end subroutine eta_levels   
!=======================================================================================================!
subroutine eta_all_levels(ak,bk,etk,am,bm,eta)

    real (kind=real_kind), dimension(nlev+1), intent(out) ::  ak, bk, etk    
    real (kind=real_kind), dimension(nlev),   intent(out) ::  am, bm, eta


    real (kind=real_kind) :: a18(18+1) = &
    (/ 0.00251499D0,  0.00710361D0,  0.01904260D0, 0.04607560D0, 0.08181860D0, &
       0.07869805D0,  0.07463175D0,  0.06955308D0, 0.06339061D0, 0.05621774D0, &
       0.04815296D0,  0.03949230D0,  0.03058456D0, 0.02193336D0, 0.01403670D0, &
       0.007458598D0, 0.002646866D0, 0.00000000D0, 0.00000000D0  /)

    real (kind=real_kind) :: b18(18+1) = &
    (/ 0.000000D0,   0.000000D0,   0.000000D0,  0.000000D0,  0.000000D0, &
       0.03756984D0, 0.08652625D0, 0.1476709D0, 0.221864D0,  0.308222D0, &
       0.4053179D0,  0.509588D0,   0.6168328D0, 0.7209891D0, 0.816061D0, &
       0.8952581D0,  0.953189D0,   0.985056D0,  1.000000D0   /)
    

    real (kind=real_kind) :: a26(26+1) = &
    (/ 0.002194067D0, 0.004895209D0, 0.009882418D0, 0.01805201D0, 0.02983724D0, &
       0.04462334D0,  0.06160587D0,  0.07851243D0,  0.07731271D0, 0.07590131D0, &
       0.07424086D0,  0.07228744D0,  0.06998933D0,  0.06728574D0, 0.06410509D0, &
       0.06036322D0,  0.05596111D0,  0.05078225D0,  0.04468960D0, 0.03752191D0, &
       0.02908949D0,  0.02084739D0,  0.01334443D0,  0.00708499D0, 0.00252136D0, &
       0.00000000D0,  0.00000000D0  /)
  
    real (kind=real_kind) :: b26(26+1) = &
    (/ 0.00000000D0, 0.00000000D0, 0.00000000D0, 0.00000000D0, 0.00000000D0, &
        0.00000000D0, 0.00000000D0, 0.00000000D0, 0.01505309D0, 0.03276228D0, &
        0.05359622D0, 0.07810627D0, 0.1069411D0,  0.1408637D0,  0.1807720D0,  &
        0.2277220D0,  0.2829562D0,  0.3479364D0,  0.4243822D0,  0.5143168D0,  &
        0.6201202D0,  0.7235355D0,  0.8176768D0,  0.8962153D0,  0.9534761D0,  &
        0.9851122D0,  1.0000000D0  /)

    real (kind=real_kind) :: a49(49+1) = &
    (/  0.0022518650D0, 0.0039838900D0, 0.0067043640D0, 0.0107323100D0, 0.0163423300D0, &
	0.0236711900D0, 0.0326145600D0, 0.0427452700D0, 0.0538261000D0, 0.0651217500D0, &
	0.0756985000D0, 0.0845428300D0, 0.0839631000D0, 0.0833410300D0, 0.0826735200D0, &
	0.0819572500D0, 0.0811886600D0, 0.0803639300D0, 0.0794789500D0, 0.0785293400D0, &
	0.0775103600D0, 0.0764169500D0, 0.0752436800D0, 0.0739847000D0, 0.0726337500D0, &
	0.0711841400D0, 0.0696286300D0, 0.0679595000D0, 0.0661684600D0, 0.0642465800D0, &
	0.0621843300D0, 0.0599714400D0, 0.0575969000D0, 0.0550489200D0, 0.0523148300D0, &
	0.0493810200D0, 0.0462329200D0, 0.0428548700D0, 0.0392300600D0, 0.0353404900D0, &
	0.0311668100D0, 0.0266882500D0, 0.0218825700D0, 0.0167637100D0, 0.0120817100D0, &
	0.0079596120D0, 0.0045102970D0, 0.0018312150D0, 0.0000000000D0, 0.0000000000D0 /)

    real (kind=real_kind) :: b49(49+1) = &
    (/  0.0000000000D0, 0.0000000000D0, 0.0000000000D0, 0.0000000000D0, 0.0000000000D0, &
	0.0000000000D0, 0.0000000000D0, 0.0000000000D0, 0.0000000000D0, 0.0000000000D0, &
	0.0000000000D0, 0.0000000000D0, 0.0067551120D0, 0.0140036400D0, 0.0217816400D0, &
	0.0301277800D0, 0.0390835600D0, 0.0486935200D0, 0.0590054200D0, 0.0700705600D0, &
	0.0819439400D0, 0.0946845900D0, 0.1083559000D0, 0.1230258000D0, 0.1387673000D0, &
	0.1556586000D0, 0.1737837000D0, 0.1932327000D0, 0.2141024000D0, 0.2364965000D0, &
	0.2605264000D0, 0.2863115000D0, 0.3139801000D0, 0.3436697000D0, 0.3755280000D0, &
	0.4097133000D0, 0.4463958000D0, 0.4857576000D0, 0.5279946000D0, 0.5733168000D0, &
	0.6219495000D0, 0.6741346000D0, 0.7301315000D0, 0.7897776000D0, 0.8443334000D0, &
	0.8923650000D0, 0.9325572000D0, 0.9637744000D0, 0.9851122000D0, 1.0000000000D0 /)
    
    integer :: k
!=======================================================================================================!
    if (nlev == 18 ) then

       do k = 1, nlev+1
          ak(k)= a18(k)
          bk(k)= b18(k)
          etk(k) = ak(k) + bk(k)
       end do

    elseif (nlev == 26) then

       do k = 1, nlev+1
          ak(k)= a26(k)
          bk(k)= b26(k)
          etk(k) = ak(k) + bk(k)
       end do

    elseif (nlev == 49) then

       do k = 1, nlev+1
          ak(k)= a49(k)
          bk(k)= b49(k)
          etk(k) = ak(k) + bk(k)
       end do

    endif 

       do k = 1, nlev
          am(k) = (ak(k) + ak(k+1)) * 0.5D0            
	  bm(k) = (bk(k) + bk(k+1)) * 0.5D0          
          eta(k)= (etk(k)+ etk(k+1))* 0.5D0 
       end do

end subroutine eta_all_levels   
!=======================================================================================================!
!=======================================================================================================!
        Function index_search(grid,n0,n1,xp) result(ii)
!=======================================================================================================!
!	Searching the position of (xp) on a "vertical" grid by Bisection
        Implicit None

          Integer, Intent(in) :: n0,n1
          real(kind=real_kind), Intent(in) :: xp
          real(kind=real_kind), Intent(in), Dimension(n0:nlev+n1):: grid
          Integer  :: ii, nm,na,nb
            na = n0
            nb = nlev+n1
             do
               if  ((nb-na) <=  1)  exit
               nm = (nb + na)/2
                if (xp  >  grid(nm)) then
                 na = nm
                else
                 nb = nm
                endif
             enddo

              ii = na

      end Function index_search
!=======================================================================================================!       
!=======================================================================================================!
!	Linear interpolation
!=======================================================================================================!
Function  Linear_Int(x1,x2,y1,y2,xin)  result(yout)

      Implicit None
      Real(Kind=real_kind), Intent(in) :: x1,x2,y1,y2 ,xin

      Real(Kind=real_kind) :: ss, yout

          if (x2 == x1)  then
               yout = y1
            else
                ss = (y2 - y1)/(x2 - x1)
                  yout = y1 + ss * (xin - x1)
          endif
!=======================================================================================================!
end Function Linear_Int
!=======================================================================================================!
!=======================================================================================================!
function intp_cubiclag(xi,yy,xval) result(yval)

   real (kind=real_kind), dimension(4), intent(in) :: xi, yy
   real (kind=real_kind), intent(in)   :: xval

  ! real (kind=real_kind)   :: xx(4),yy(4),xval

    real (kind=real_kind)   :: xx(4), w(10), yval, prod,sm1,sm2,eps

    integer :: i,j, n 
!=======================================================================================================!
        xx(:) = xi(:)

          n = 4          !For cubic interpolation  y=f(x)
        eps = 1.0D-10

        yval = 0.0D0 
 
    do j = 1, n-1
          if ((xi(j+1)-xi(j)) ==  0.0D0) then
               print*, xi(j), 'good luck'
              xx(j+1) = xi(j+1) + eps
          endif
    enddo
 
    do  j=1,n
         if (xval == xx(j)) then
            yval=yy(j)
            exit
         endif

           prod=1.0D0 
 
        do  i=1,n
           if (j /= i) then
             prod = prod*(xx(j)-xx(i))
           endif
           w(j)=1.0D0 /(prod*(xval-xx(j)))
        end do

    end do

        sm1  = 0.0D0 
        sm2  = 0.0D0 

      do  j=1,n
        sm1 = sm1 + w(j)*yy(j)
        sm2 = sm2 + w(j)
      end do

        if (sm2  ==  0.0D0) print*, 'Intp-error: CubicLag', xx(1),xx(2),xx(3),xx(4)
        yval = sm1/sm2
!=======================================================================================================!
 end function Intp_CubicLag 
!=======================================================================================================!
!=======================================================================================================!
  subroutine fluxjacobian_bcl(metinv,uvbuf,dpbuf,ptbuf,uvcomp,dp,pt,ptop,gbot,fjbcl)
!=======================================================================================================!
    Implicit None
    integer, parameter :: south=1, east=2, north=3, west=4
    real (kind=real_kind),                                  intent(in)   :: metinv(2,2,np,np)
    real (kind=real_kind), dimension(0:np+1,0:np+1,nlev),   intent(in)   :: dpbuf, ptbuf
    real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev), intent(in)   :: uvbuf
    real (kind=real_kind), dimension(np,np,2,nlev),         intent(in)   :: uvcomp
    real (kind=real_kind), dimension(np,np,nlev),           intent(in)   :: pt, dp
    real (kind=real_kind), dimension(np,np),                intent(in)   :: ptop,gbot
    real (kind=real_kind), dimension(4,nlev),               intent(out)  :: fjbcl
!=======================================================================================================!
    real (kind=real_kind), dimension(np,4,2,nlev)  :: uv_senw, uv_edge
    real (kind=real_kind), dimension(np,4,nlev+1)  :: pr_edge, pr_senw, gh_edge, gh_senw,pp_edge, pp_senw
    real (kind=real_kind), dimension(np,4,nlev)    :: dp_edge, pt_edge, uvmax,ghmax
    real (kind=real_kind), dimension(np,4,nlev)    :: dp_senw, pt_senw, dgp_senw, dgp_edge, ptg_senw
    real (kind=real_kind), dimension(np,4,nlev)    :: ph_senw, ph_edge
    real (kind=real_kind), dimension(np,4)         :: gij_edge, ptop_edge, gbot_edge
!=======================================================================================================!
    real (kind=real_kind):: rtmp,ghedge, ghsenw
    real (kind=real_kind):: rdcp, cprd, pp_0, c_pp, alfa, r_pp, temp
    Integer:: i,j,k, wall
!=======================================================================================================!
           rdcp = r_d / c_p
           cprd = 1.0D0 / rdcp
           c_pp = c_p / (p_0)**rdcp
           r_pp = r_d / (p_0)**rdcp
!=======================================================================================================!
!=======================================================================================================!
!       Boundary values for velocity & flux terms from the neighbors
!=======================================================================================================!
  do k = 1, nlev
   do j = 1, 2
   do i = 1, np
      uv_senw(i,south,j,k) = uvbuf(i,   0,j,k)
      uv_senw(i,east,j,k)  = uvbuf(np+1,i,j,k)
      uv_senw(i,north,j,k) = uvbuf(i,np+1,j,k)
      uv_senw(i,west,j,k)  = uvbuf(0,i,   j,k)

      uv_edge(i,south,j,k) = uvcomp(i, 1,j,k)
      uv_edge(i,east,j,k)  = uvcomp(np,i,j,k)
      uv_edge(i,north,j,k) = uvcomp(i,np,j,k)
      uv_edge(i,west,j,k)  = uvcomp(1,i, j,k)
   end do
   end do
  end do
!=======================================================================================================!
 do k = 1, nlev
  do i = 1, np
    dp_senw(i,south,k) = dpbuf(i,   0,k)
    dp_senw(i,east,k)  = dpbuf(np+1,i,k)
    dp_senw(i,north,k) = dpbuf(i,np+1,k)
    dp_senw(i,west,k)  = dpbuf(0,i   ,k)

    pt_senw(i,south,k) = ptbuf(i,   0,k)
    pt_senw(i,east,k)  = ptbuf(np+1,i,k)
    pt_senw(i,north,k) = ptbuf(i,np+1,k)
    pt_senw(i,west,k)  = ptbuf(0,i   ,k)

    dp_edge(i,south,k) = dp(i, 1,k)
    dp_edge(i,east,k)  = dp(np,i,k)
    dp_edge(i,north,k) = dp(i,np,k)
    dp_edge(i,west,k)  = dp(1,i ,k)

    pt_edge(i,south,k) = pt(i, 1,k)
    pt_edge(i,east,k)  = pt(np,i,k)
    pt_edge(i,north,k) = pt(i,np,k)
    pt_edge(i,west,k)  = pt(1,i ,k)
  end do
 end do

  do i = 1, np
    ptop_edge(i,south) = ptop(i, 1)
    ptop_edge(i,east)  = ptop(np,i)
    ptop_edge(i,north) = ptop(i,np)
    ptop_edge(i,west)  = ptop(1, i)

    gbot_edge(i,south) = gbot(i, 1)
    gbot_edge(i,east)  = gbot(np,i)
    gbot_edge(i,north) = gbot(i,np)
    gbot_edge(i,west)  = gbot(1, i)

    gij_edge(i,south) = sqrt(metinv(2,2,i,1))
    gij_edge(i,north) = sqrt(metinv(2,2,i,np))
    gij_edge(i,east)  = sqrt(metinv(1,1,np,i))
    gij_edge(i,west)  = sqrt(metinv(1,1,1,i))
  end do

! For pressure thickness case
!  pr_edge(:,:,1) = ptop_edge(:,:)
!  pr_senw(:,:,1) = ptop_edge(:,:)
!
! do k = 1, nlev
!   do wall = 1, 4
!     do i = 1, np
!      pr_edge(i,wall,k+1) = pr_edge(i,wall,k) + dp_edge(i,wall,k)
!      pr_senw(i,wall,k+1) = pr_senw(i,wall,k) + dp_senw(i,wall,k)
!     end do
!   end do
! end do


 !!Geopotential thickness from the element edges and the neighbouring halo regions
 !!In this case,    dp = dgp 

  gh_edge(:,:,nlev+1) = abs(gbot_edge(:,:))
  gh_senw(:,:,nlev+1) = abs(gbot_edge(:,:))

 do k = nlev, 1, -1
   do wall = 1, 4
    do i = 1, np
      gh_edge(i,wall,k) = gh_edge(i,wall,k+1) + dp_edge(i,wall,k)
      gh_senw(i,wall,k) = gh_senw(i,wall,k+1) + dp_senw(i,wall,k)
    end do
   end do
 end do

 do k = 1, nlev
    do wall = 1, 4
      do i = 1, np
          ghmax(i,wall,k) = max(abs(gh_edge(i,wall,k)),abs(gh_senw(i,wall,k)))
      end do
    end do
  end do

 do k = 1, nlev
      do i = 1, np
        uvmax(i,south,k) = max(abs(uv_edge(i,south,2,k)),abs(uv_senw(i,south,2,k)))
        uvmax(i,north,k) = max(abs(uv_edge(i,north,2,k)),abs(uv_senw(i,north,2,k)))
         uvmax(i,east,k) = max(abs(uv_edge(i,east,1,k)),abs(uv_senw(i,east,1,k)))
         uvmax(i,west,k) = max(abs(uv_edge(i,west,1,k)),abs(uv_senw(i,west,1,k)))
      end do
  end do

 do k = 1, nlev
    do wall = 1, 4
      alfa = 0.0D0
      do i = 1, np
         alfa = max(alfa,(uvmax(i,wall,k) + sqrt(ghmax(i,wall,k)*gij_edge(i,wall))  ))
      end do
        fjbcl(wall,k) = alfa
    end do
  end do

!=======================================================================================================!
!=======================================================================================================!
 end subroutine fluxjacobian_bcl
!=======================================================================================================!
end module dg3d_vertical_mod
