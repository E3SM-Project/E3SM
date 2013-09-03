#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dg3d_remap_mod
!=======================================================================================================!
! R. D. Nair 04/04/07 
!=======================================================================================================!
  use kinds
  ! ------------------------
  use physical_constants
  ! ------------------------
  use dimensions_mod
  ! ------------------------
  use derivative_mod
  ! ------------------------
  use coordinate_systems_mod
  ! ------------------------
  use quadrature_mod
  ! ------------------------
  use time_mod
  ! ------------------------
 use dg3d_vertical_mod, only : eta_levels
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
   real (kind=real_kind), private, parameter  :: r_d = 287.0D0      ! Gas const (dry)
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
  public  :: linear_remap
  public  :: parabolic_remap  
  public  :: energy_remap
  public  :: monotonic_remap 

  private  :: intp_cubiclag
  private  :: linear_int
  private  :: index_search
  private  :: ppm_remap
  private  :: ppm_coef
  private  :: density_new
  private  :: monoton_hoint 
!=======================================================================================================!
 contains
!=======================================================================================================!
!=======================================================================================================!
subroutine linear_remap(ie,sgp,ptop,pt3d,uv3d,dp3d)
!=======================================================================================================!
    integer, intent(in)                  :: ie
    real (kind=real_kind), intent(in)    :: sgp(np,np), ptop(np,np) 
    real (kind=real_kind), intent(inout) :: pt3d(np,np,nlev),uv3d(np,np,2,nlev)   
    real (kind=real_kind), intent(inout) :: dp3d(np,np,nlev)       
!=======================================================================================================!
    real (kind=real_kind), dimension(np,np,nlev+1) :: pr3d_ref, pr3d
    real (kind=real_kind), dimension(np,np,nlev) :: dgp3d
    real (kind=real_kind), dimension(0:nlev+1)   :: pgrid, logp, pref, uk,vk,tk
    real (kind=real_kind) ::  xx(4), yy(4)          
    real (kind=real_kind) ::  rdcp, p1,p2,pp, px,f1,f2, cprd
    integer :: indx(nlev)
    integer :: i,j,k , ik
    real (kind=real_kind) ::  eta(nlev), etv(nlev), etk(nlev+1)
    real (kind=real_kind) ::  ak(nlev+1),bk(nlev+1)
    real (kind=real_kind) ::  trm1,trm2,trm3,trm4, term
!=======================================================================================================!
!	Resetting progostic variables at reference coordinates
!=======================================================================================================!
           rdcp = r_d / c_p
           cprd = c_p / r_d 

      Call eta_levels(ak,bk,eta,etk)

   !Initial reference pressure


   !Pressure at interfaces from pressure-thickness (Lagrangian)

    pr3d(:,:,1) = ptop(:,:)

  ! For dP thickness 
    do k=1,nlev
       do j=1,np
          do i=1,np
            pr3d(i,j,k+1) = pr3d(i,j,k) + dp3d(i,j,k)
          end do
       end do
    end do

 !Computed pressure at interfaces from GP thickness (for GP-thickness case)
 !Note: dp3d --> GP thickness
  !  do k=1, nlev
  !     do j=1,np
  !        do i=1,np
  !           pr3d(i,j,k+1) = p_0 * ( (pr3d(i,j,k)/p_0)**rdcp  + dp3d(i,j,k)/(c_p* pt3d(i,j,k)) )**cprd
  !        end do
  !     end do
  !  end do

  !New reference Pressure at etk(k) levels (Eulerian interfaces)

   do k = 1, nlev+1
    do j=1,np
        do i=1,np
          pr3d_ref(i,j,k) = p_0 * ak(k) + bk(k)* pr3d(i,j,nlev+1)
        enddo
    enddo
   enddo

  !New Pressure thickness defined (dP-thickness case)
   do k = 1, nlev
    do j=1,np
        do i=1,np
            dp3d(i,j,k) =  pr3d_ref(i,j,k+1) - pr3d_ref(i,j,k)
        enddo
    enddo
   enddo

  !Remapping onto the new reference coordinates

     do j=1,np
          do i=1,np

            do k=1,nlev
                pgrid(k)= (pr3d(i,j,k+1) + pr3d(i,j,k))*0.5D0
                pref(k) = (pr3d_ref(i,j,k) + pr3d_ref(i,j,k+1))*0.5D0 
                uk(k) =  uv3d(i,j,1,k)
                vk(k) =  uv3d(i,j,2,k)
                tk(k) =  pt3d(i,j,k)
            end do

                pgrid(0) = pr3d(i,j,1)
                pgrid(nlev+1) = pr3d(i,j,nlev+1)

                uk(0) =  (3.0D0 *uk(1) - uk(2) )* 0.5D0
                uk(nlev+1) =  (3.0D0 *uk(nlev) - uk(nlev-1) )* 0.5D0
                vk(0) =  (3.0D0 *vk(1) - vk(2) )* 0.5D0
                vk(nlev+1) =  (3.0D0 *vk(nlev) - vk(nlev-1) )* 0.5D0

                tk(0) =  (3.0D0 *tk(1) - tk(2) )* 0.5D0
                tk(nlev+1) =  (3.0D0 *tk(nlev) - tk(nlev-1) )* 0.5D0

        ! index  search

             do k = 1, nlev
               indx(k) = index_search(pgrid,0,1,pref(k))
             enddo

             do k = 1, nlev
                px = pref(k)
                ik = indx(k)
                p1 = pgrid(ik)
                p2 = pgrid(ik+1)

                f1 = uk(ik)
                f2 = uk(ik+1)
               uv3d(i,j,1,k) = Linear_Int(p1,p2,f1,f2,px)

                f1 = vk(ik)
                f2 = vk(ik+1)
               uv3d(i,j,2,k) = Linear_Int(p1,p2,f1,f2,px)

              !if ((j==1).and.(i==3)) print*,k, pref(k), pgrid(ik+1)

                p1 = log(pgrid(ik))
                p2 = log(pgrid(ik+1))
                f1 = tk(ik)
                f2 = tk(ik+1)
                px = log(px)
               pt3d(i,j,k) = Linear_Int(p1,p2,f1,f2,px)
               
             enddo
              
       end do
    end do


  !New GP thickness on the reference grid
  ! do k=nlev, 1, -1
  !    do j=1,np
  !       do i=1,np
  !           dp3d(i,j,k) =   c_p * pt3d(i,j,k)* ((pr3d_ref(i,j,k+1)/p_0)**(rdcp) - &
  !                                                (pr3d_ref(i,j,k)/p_0)**(rdcp))
  !       end do
  !    end do
  ! end do

!=======================================================================================================!
 end subroutine linear_remap 
!=======================================================================================================!
!=======================================================================================================!
subroutine energy_remap(metinv,sgp,ptop,pt3d,uv3d,thick)

    real (kind=real_kind), intent(in)    :: metinv(2,2,np,np)
    real (kind=real_kind), intent(in)    :: sgp(np,np), ptop(np,np)
    !real (kind=real_kind), intent(inout) :: pt3d(np,np,nlev), u3d(np,np,nlev),v3d(np,np,nlev)
    real (kind=real_kind), intent(inout) :: pt3d(np,np,nlev),uv3d(np,np,2,nlev)   
    real (kind=real_kind), intent(inout) :: thick(np,np,nlev)     !thickness (dP or dGP)

    real (kind=real_kind), dimension(nlev+1,np,np) :: pr3d_ref, pr3d
    real (kind=real_kind), dimension(0:nlev+1)   :: pgrid, logp, pref, uk,vk,tk
    real (kind=real_kind), dimension(nlev,np,np) :: uz, vz, tz, dpz
    real (kind=real_kind), dimension(nlev) :: ru, rv, rt

    real (kind=real_kind), dimension(np,np,nlev) :: energy_lag, energy_eul, dlnp, dprs
    real (kind=real_kind), dimension(np,np,nlev+1) :: pr_eta, gp_eta, lnpr,prka

    integer :: indx(nlev), ipr(nlev)
    integer :: i,j,k , ik

    real (kind=real_kind) ::  trm1,trm2,trm3,trm4, term
    real (kind=real_kind) ::  xx(4), yy(4)
    real (kind=real_kind) ::  rdcp, p1,p2,pp, px,f1,f2, cprd,c_pp
    real (kind=real_kind) ::  eta(nlev), etv(nlev), etk(nlev+1)
    real (kind=real_kind) ::  ak(nlev+1),bk(nlev+1), sigma_level
    real (kind=real_kind) ::  pge(nlev+1),pgl(nlev+1)
    real (kind=real_kind) ::  dp(nlev), tza(nlev),coef(nlev,8), fo(nlev)

    real (kind=real_kind) ::  top,bot,temp, peta, e_pe,e_th,e_ke, pp_0, dgp3d ,u1,u2,v1,v2
    real (kind=real_kind) ::  sigtop, sigbot, sigk 


   ! Reseting progostic variables at reference coordinates

           rdcp = r_d / c_p
           cprd = c_p / r_d

           pp_0 = (p_0)**rdcp
           c_pp = c_p / pp_0


      Call eta_levels(ak,bk,eta,etk)

 !Pressure at etk(k) levels (interfaces)

   pr_eta(:,:,1) = ptop(:,:)
   gp_eta(:,:,nlev+1) = sgp(:,:)

  !! for dP-thickness 
    do k=1,nlev
       do j=1,np
          do i=1,np
            pr_eta(i,j,k+1) = pr_eta(i,j,k) + thick(i,j,k)
          end do
       end do
    end do

 !Computed pressure at interfaces from GP thickness (for GP-thickness case)
 !Note: thick --> GP thickness
 !  do k=1, nlev
 !    do j=1,np
 !       do i=1,np
 !          pr_eta(i,j,k+1) = p_0 * ( (pr_eta(i,j,k)/p_0)**rdcp  + thick(i,j,k)/(c_p* pt3d(i,j,k)) )**cprd
 !       end do
 !    end do
 !  end do

  do k=1,nlev
    do j=1,np
       do i=1,np
          dprs(i,j,k) = pr_eta(i,j,k+1) - pr_eta(i,j,k)
       end do
    end do
  end do

   do k=1,nlev+1
      do j=1,np
         do i=1,np
           lnpr(i,j,k) = log(pr_eta(i,j,k))
           prka(i,j,k) = (pr_eta(i,j,k))**rdcp
         end do
      end do
   end do

 ! Geopotential at interfaces


 !for pressure thickness case (dP)
    do k=nlev, 1, -1
       do j=1,np
          do i=1,np
              dgp3d =  c_p * pt3d(i,j,k)* ((pr_eta(i,j,k+1)/p_0)**(rdcp) - (pr_eta(i,j,k)/p_0)**(rdcp))
              gp_eta(i,j,k) =   gp_eta(i,j,k+1) + dgp3d
          end do
       end do
    end do

 !For GP thickness case
 !   do k=nlev, 1, -1
 !      do j=1,np
 !         do i=1,np
 !             gp_eta(i,j,k) =   gp_eta(i,j,k+1) + thick(i,j,k)
 !         end do
 !      end do
 !   end do

    !Total energy in Lagrangan coordinates:

    do j=1,np
        do i=1,np
           do k=1,nlev

             peta  = (pr_eta(i,j,k+1) + pr_eta(i,j,k))*0.5D0

             e_th = c_p * pt3d(i,j,k)/(rdcp*pp_0) * (prka(i,j,k+1) - prka(i,j,k))/ (lnpr(i,j,k+1)- lnpr(i,j,k))  

             e_pe = (pr_eta(i,j,k+1)*gp_eta(i,j,k+1) - pr_eta(i,j,k)*gp_eta(i,j,k))/ dprs(i,j,k)

             !!Note : (v1,v2) --> covariant;  (u1,u2) --> contravariant
                    v1 = uv3d(i,j,1,k)
                    v2 = uv3d(i,j,2,k)
                    u1 = metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2
                    u2 = metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2
             e_ke = 0.5D0 * (u1*v1 + u2*v2)

             energy_lag(i,j,k) = (e_th + e_pe + e_ke)

           end do
        end do
    end do


 !!!reversing the k-loop

       do j=1,np
          do i=1,np
             do k=1,nlev
               dpz(k,i,j) = dprs(i,j,k)
                !tz(k,i,j) = pt3d(i,j,k)
                tz(k,i,j) = energy_lag(i,j,k)
                uz(k,i,j) = uv3d(i,j,1,k)
                vz(k,i,j) = uv3d(i,j,2,k)
             end do
          end do
       end do

   !Pressure at interfaces from pressure-thickness (Lagrangian)

    pr3d(1,:,:) = ptop(:,:)

    do j=1,np
        do i=1,np
            do k=1,nlev
              pr3d(k+1,i,j) = pr3d(k,i,j) + dpz(k,i,j)
            end do
        end do
    end do

  !New reference Pressure at etk(k) levels (Eulerian interfaces)
          do j=1,np
           do i=1,np
              do k = 1, nlev+1
                 pr3d_ref(k,i,j) = p_0 * ak(k) + bk(k)* pr3d(nlev+1,i,j)
            enddo
            enddo
         enddo

! Normalized pressure (Sigma)  reference coordinates 
!        do j=1,np
!          do i=1,np
!               sigtop = ptop(i,j)/pr3d(nlev+1,i,j)
!               sigbot = 1.0D0 
!              do k = 1, nlev+1
!               sigk =  sigtop + (k-1)*(sigbot - sigtop )/dble(nlev) 
!               pr3d_ref(k,i,j) =  sigk* pr3d(nlev+1,i,j)
!              enddo
!          enddo
!        enddo

   !New Pressure thickness defined

    do j=1,np
        do i=1,np
          do k = 1, nlev
            dprs(i,j,k) =  pr3d_ref(k+1,i,j) - pr3d_ref(k,i,j)
          enddo
       enddo
    enddo

  !Remapping onto the new reference coordinates
    do j=1,np
          do i=1,np

            do k=1,nlev
                pgrid(k) = (pr3d(k+1,i,j) + pr3d(k,i,j))*0.5D0
                 pref(k) = (pr3d_ref(k,i,j) + pr3d_ref(k+1,i,j))*0.5D0
                   uk(k) = uz(k,i,j)
                   vk(k) = vz(k,i,j)
                   tk(k) = tz(k,i,j)
            end do

          !Pre-computation for parabolic remap

               pgl(:) = pr3d(:,i,j)
               pge(:) = pr3d_ref(:,i,j)
                dp(:) = dpz(:,i,j)
               ipr(:) = 0

            do k = 1, nlev
               ipr(k) = index_search(pgl,1,1,pge(k))
            enddo

             Call ppm_coef(dp,coef)

             Call ppm_remap(pge,pgl,uz(:,i,j),dp,ipr,coef,ru)
             Call ppm_remap(pge,pgl,vz(:,i,j),dp,ipr,coef,rv)
             Call ppm_remap(pge,pgl,tz(:,i,j),dp,ipr,coef,rt)

           ! Linear extrapolations for top/bot boundary regions

                pgrid(0) = pr3d(1,i,j)
                pgrid(nlev+1) = pr3d(nlev+1,i,j)

                uk(0) =  (3.0D0 *uk(1) - uk(2) )* 0.5D0
                uk(nlev+1) =  (3.0D0 *uk(nlev) - uk(nlev-1) )* 0.5D0
                vk(0) =  (3.0D0 *vk(1) - vk(2) )* 0.5D0
                vk(nlev+1) =  (3.0D0 *vk(nlev) - vk(nlev-1) )* 0.5D0

                tk(0) =  (3.0D0 *tk(1) - tk(2) )* 0.5D0
                tk(nlev+1) =  (3.0D0 *tk(nlev) - tk(nlev-1) )* 0.5D0

          do k = 1, nlev

              if ((k > 2).and.(k < nlev-1)) then      !parabolic remap
             !if ((k > 3).and.(k < nlev-2)) then      !parabolic remap

                uz(k,i,j) = ru(k)
                vz(k,i,j) = rv(k)
                tz(k,i,j) = rt(k)

              else                                    !Linear remap

                ik = index_search(pgrid,0,1,pref(k))
                px = pref(k)
                p1 = pgrid(ik)
                p2 = pgrid(ik+1)

                f1 = uk(ik)
                f2 = uk(ik+1)

               uz(k,i,j) = Linear_Int(p1,p2,f1,f2,px)

                f1 = vk(ik)
                f2 = vk(ik+1)

               vz(k,i,j) = Linear_Int(p1,p2,f1,f2,px)

                p1 = log(pgrid(ik))
                p2 = log(pgrid(ik+1))
                f1 = tk(ik)
                f2 = tk(ik+1)
                px = log(px)

               tz(k,i,j) = Linear_Int(p1,p2,f1,f2,px)

             endif

             enddo

       end do
    end do

 !reversing the k-loop (for remapped fields)

    do k=1,nlev
       do j=1,np
          do i=1,np
               !pt3d(i,j,k) = tz(k,i,j)
                energy_eul(i,j,k) = tz(k,i,j)
                uv3d(i,j,1,k) = uz(k,i,j)
                uv3d(i,j,2,k) = vz(k,i,j)
             end do
          end do
      end do

  !Potential temperature retrieval from the total energy
 
   pr_eta(:,:,1) = ptop(:,:)

  do k=1,nlev
      do j=1,np
         do i=1,np
           pr_eta(i,j,k+1) = pr_eta(i,j,k) + dprs(i,j,k)
           dlnp(i,j,k) = log(pr_eta(i,j,k+1) / pr_eta(i,j,k))
         end do
      end do
   end do

  do k=1,nlev+1
      do j=1,np
         do i=1,np
           prka(i,j,k) = pr_eta(i,j,k)**rdcp
         end do
      end do
   end do


   gp_eta(:,:,nlev+1) = sgp(:,:)

  do k=nlev, 1, -1
    do j=1,np
        do i=1,np

             !!Note : (v1,v2) --> covariant;  (u1,u2) --> contravariant
                v1 = uv3d(i,j,1,k)
                v2 = uv3d(i,j,2,k)
                u1 = metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2
                u2 = metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2

              e_ke = 0.5D0 * (u1*v1 + u2*v2)     !new KE

                top =  energy_eul(i,j,k) - e_ke - gp_eta(i,j,k+1)
                bot =  c_p * (1.0D0 - rdcp *(pr_eta(i,j,k) /dprs(i,j,k))*dlnp(i,j,k) )

               temp = top / bot              ! temperature retrieval

               gp_eta(i,j,k) = gp_eta(i,j,k+1) + r_d * temp * dlnp(i,j,k)

              pt3d(i,j,k) = temp*(rdcp*pp_0) * dlnp(i,j,k) / (prka(i,j,k+1) - prka(i,j,k))

           end do
    end do
  end do

 !New thickness (dP or dGP)
   do k= 1,nlev
      do j=1,np
         do i=1,np
          ! thick(i,j,k) =  gp_eta(i,j,k) -  gp_eta(i,j,k+1)    !for GP case
            thick(i,j,k) = dprs(i,j,k)                          !for dP case
         end do
      end do
   end do

 end subroutine energy_remap
!=======================================================================================================!
!=======================================================================================================!
    subroutine  ppm_remap(pge,pgl,za,dp,ipr,coef,fo)
!=======================================================================================================!

    Implicit None

    Integer, Intent(in), Dimension(nlev)  :: ipr
    real(kind=real_kind), Intent(in), Dimension(nlev) :: za, dp
    real(kind=real_kind), Intent(in), Dimension(nlev+1) :: pge, pgl
    real(kind=real_kind), Intent(in), Dimension(nlev,8) :: coef
    real(kind=real_kind), Intent(out), Dimension(nlev) :: fo

    real(kind=real_kind), Dimension(nlev) :: zv,dzb,zb,a6,a6m,czb,zbl,zbr
    real(kind=real_kind), Dimension(nlev) :: daj,slp,del,ami,bal,bar,ba6
    real(kind=real_kind), Dimension(nlev) :: zmas,cmas

    real(kind=real_kind) ::  deno, fact,pr1,pr2, dmin, bmn,   &
                              fac1,fac2, dx,term, a,b,c, pom, del6, sm

    real(kind=real_kind) :: r12, xp, c1,c2,c3,av12, diff,absd, dens

    real(kind=real_kind) ::  tol, db3

    Integer :: i,j,k,kj, mono 

              
!---------------------------
!   mono = 0, 1, 2 are for positive-definite, monotonic
!    and un-filtered  solutions respectively.
!    za  -> cell-zonal average (ni gridpoinj values)
!    zb  -> values at the border of the cell
!    coef(:,:) is the precomputed coefficients for parabola making
!    Ref: Nair & Machenhauer (MWR, 2002); CISL paper 
!---------------------------

      tol = 1.0D-12
      db3 = 3.0D0
      mono = 1 

! Non-uniform grid
    do j = 2, nlev-1
         daj(j) = coef(j,1)* ( coef(j,2)*(za(j+1)-za(j)) + coef(j,3)*(za(j)-za(j-1)) )
    enddo

! Monotonize the slope

    do j = 2, nlev-1
       dmin = 0.0D0
        pom = 1.0D0
        if (daj(j) < 0.0D0) pom = -1.0D0
        term = (za(j+1) - za(j)) * (za(j) - za(j-1))
          c1 = 2.0D0* abs(za(j) - za(j+1))
          c2 = 2.0D0* abs(za(j) - za(j-1))
          c3 = abs(daj(j))

         if (term > 0.0D0) then
           dmin = min(c1,c2,c3)
         else
           dmin = 0.0D0
         endif
          slp(j) = dmin * pom
      enddo

!   Interpolation of the zonal values onto the cell border

       do j = 2, nlev-2

           fac1 = za(j+1) - za(j)
           fac2 = dp(j+1)*coef(j,6)*slp(j) - dp(j)*coef(j,5) * slp(j+1)
        zb(j+1) = za(j) + coef(j,4)* fac1 + ( fac1*coef(j,7) + fac2 ) / coef(j,8)
       enddo

!     Monotonization of the parabola.

       do j = 3, nlev-1
         zbl(j) = zb(j)
         zbr(j) = zb(j+1)
         dzb(j) = zbr(j) - zbl(j)
         bal(j) = zbl(j)
         bar(j) = zbr(j)
      enddo

      do j = 3, nlev-1
        a6(j) =  (za(j) - (zbr(j) + zbl(j))* 0.5D0) * 6.0D0
        ba6(j) = a6(j)
      enddo

!   For the Monotonic case

    if (mono  ==  1) then  
        do j = 3, nlev-1
           pr1 = a6(j) * dzb(j)
           pr2 = dzb(j) * dzb(j)
           if (slp(j) ==   0.0D0) then
             bal(j) = za(j)
             bar(j) = za(j)
             ba6(j) = 0.0D0
           elseif (pr1 <    (-pr2)) then
             ba6(j) = db3* (zbl(j) - za(j))
             bar(j) = zbl(j) - ba6(j)
             bal(j) = zbl(j)
           elseif (pr1 >   pr2) then
             ba6(j) = db3* (zbr(j) - za(j))
             bal(j) = zbr(j) - ba6(j)
             bar(j) = zbr(j)
           endif
       enddo
    endif


!   Redefining the coeffients of  a6m

    do j = 3, nlev-1
       if (mono ==  2) then       !Non-Monotonic case
          a6m(j) = 6.0D0 * ( za(j) - (zbl(j) + zbr(j))* 0.5D0 )
          del(j) = zbr(j) - zbl(j)
       else
          a6m(j) = ba6(j)
          del(j) = bar(j) - bal(j)
       endif
    enddo

!  Conservative remapping using PPM
!  Accumulated mass along the Lagrangian cells

              sm = 0.0D0

        do j = 1, nlev
              sm = sm + za(j) * dp(j)
          cmas(j) = sm
          zmas(j) = sm     !false filling
            fo(j) = 0.0D0
        enddo

!    Remapping to the "Eulerian zones"

        do j = 3, nlev-1
       !do j = 3, nlev-2

          kj = ipr(j)
          
         ! if ((kj-1) < 1) then 
         !    print *, 'kj=', kj, j 
         !    call abort
         ! end if
          
         dx = pge(j) - pgl(kj)

            c1 = za(kj)
            c2 = del(kj)
            c3 = a6m(kj)

          av12 = c1 * 12.0D0
          del6 = c2 * 6.0D0

!       When xp is in the rage [-1/2,+1/2]

          xp = dx / dp(kj)  -  0.5D0
          diff = xp +  0.5D0
          absd = abs(1.0D0 - diff)

        if (diff < tol) then
           dens = 0.0D0
          elseif (absd <  tol) then
            dens = za(kj)
          else
           dens =  Density_new(xp,c1,c2,c3,av12,del6)
        endif

          zmas(j) =  cmas(kj-1) + dens  * dp(kj)

      enddo

  !skipping two cells from the top & bottom levels

         do j = 3, nlev-2
            fo(j)= (zmas(j+1) - zmas(j))/(pge(j+1) - pge(j))
         enddo

   end subroutine ppm_remap
!=======================================================================================================!
!=======================================================================================================!
  subroutine ppm_coef(dp,coef)
!=======================================================================================================!
! Pre-computation of parbla coefficients & edge-zone interpolations

      Implicit None

      real (kind=real_kind), dimension(nlev), intent(in)      :: dp
      real (kind=real_kind), dimension(nlev,8),intent(out)    :: coef

      real (kind=real_kind) :: term
      integer :: j

       coef = 0.0D0

      do j = 2, nlev-1
          coef(j,1) = dp(j) / (dp(j-1) + dp(j) + dp(j+1))
          coef(j,2) = (2.0D0 *dp(j-1) + dp(j) ) / (dp(j+1) + dp(j))
          coef(j,3) = (2.0D0 *dp(j+1) + dp(j) ) / (dp(j-1) + dp(j))
      enddo

!   Interpolation of the zonal values onto the cell border

        do j = 2, nlev-2
          term = 2.0D0 *dp(j+1) * dp(j) / (dp(j) + dp(j+1))
          coef(j,4) = dp(j) / (dp(j) + dp(j+1))
          coef(j,5) = (dp(j-1) + dp(j)) / (2.0D0 *dp(j) + dp(j+1))
          coef(j,6) = (dp(j+1) + dp(j+2)) / (2.0D0 *dp(j+1) + dp(j))
          coef(j,7) = term* (coef(j,5) - coef(j,6))
          coef(j,8) = dp(j-1) + dp(j) + dp(j+1) + dp(j+2)
       enddo

 end subroutine ppm_coef
!=======================================================================================================!
!=======================================================================================================!
function  density_new(xp,c1,c2,c3,av12,del6)  result(den)
!=======================================================================================================!
! Analytic Integration of the parabola
      Implicit None

      real(kind=real_kind), Intent(in) :: c1,c2,c3,xp,av12,del6
      real(kind=real_kind) :: zup,zlow,den

          zup = xp * (av12 + del6*xp + c3*(1.0D0 - xp*xp*4.0D0) )
         zlow = 1.5D0 * c2 - 6.0D0 * c1

         den = (zup - zlow) / 12.0D0

    end function density_new
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
!   This is based on  conservative PPM interpolation without total energy constraint
!=======================================================================================================!
subroutine parabolic_remap(metinv,sgp,ptop,pt3d,uv3d,thick,qt3d)

    real (kind=real_kind), intent(in)    :: metinv(2,2,np,np)
    real (kind=real_kind), intent(in)    :: sgp(np,np), ptop(np,np)
    real (kind=real_kind), intent(inout) :: pt3d(np,np,nlev),uv3d(np,np,2,nlev)   
    real (kind=real_kind), intent(inout) :: thick(np,np,nlev)     !thickness (dP or dGP)
    real (kind=real_kind), intent(inout) :: qt3d(np,np,nlev)

    real (kind=real_kind), dimension(nlev+1,np,np) :: pr3d_ref, pr3d
    real (kind=real_kind), dimension(0:nlev+1)   :: pgrid, logp, pref, uk,vk,tk, qk 
    real (kind=real_kind), dimension(nlev,np,np) :: uz, vz, tz, qz, dpz
    real (kind=real_kind), dimension(nlev) :: ru, rv, rt, rq 

    real (kind=real_kind), dimension(np,np,nlev) :: energy_lag, energy_eul, dlnp, dprs
    real (kind=real_kind), dimension(np,np,nlev+1) :: pr_eta, gp_eta, lnpr,prka

    integer :: indx(nlev), ipr(nlev)
    integer :: i,j,k , ik

    real (kind=real_kind) ::  trm1,trm2,trm3,trm4, term
    real (kind=real_kind) ::  xx(4), yy(4)
    real (kind=real_kind) ::  rdcp, p1,p2,pp, px,f1,f2, cprd,c_pp
    real (kind=real_kind) ::  eta(nlev), etv(nlev), etk(nlev+1)
    real (kind=real_kind) ::  ak(nlev+1),bk(nlev+1)
    real (kind=real_kind) ::  pge(nlev+1),pgl(nlev+1)
    real (kind=real_kind) ::  dp(nlev), tza(nlev),coef(nlev,8), fo(nlev)

    real (kind=real_kind) ::  top,bot,temp, peta, e_pe,e_th,e_ke, pp_0, dgp3d ,u1,u2,v1,v2


   ! Reseting progostic variables at reference coordinates

           rdcp = r_d / c_p
           cprd = c_p / r_d

           pp_0 = (p_0)**rdcp
           c_pp = c_p / pp_0

      Call eta_levels(ak,bk,eta,etk)

 !Pressure at etk(k) levels (interfaces)

   pr_eta(:,:,1) = ptop(:,:)

 ! for dP-thickness 
    do k=1,nlev
       do j=1,np
          do i=1,np
            pr_eta(i,j,k+1) = pr_eta(i,j,k) + thick(i,j,k)
          end do
       end do
    end do

 !Computed pressure at interfaces from GP thickness (for GP-thickness case)
 !Note: thick --> GP thickness           (on Lagrangian surface)
 ! do k=1, nlev
 !  do j=1,np
 !      do i=1,np
 !         pr_eta(i,j,k+1) = p_0 * ( (pr_eta(i,j,k)/p_0)**rdcp  + thick(i,j,k)/(c_p* pt3d(i,j,k)) )**cprd
 !      end do
 !  end do
 ! end do

  do k=1,nlev
    do j=1,np
       do i=1,np
          dprs(i,j,k) = pr_eta(i,j,k+1) - pr_eta(i,j,k)
       end do
    end do
  end do

    !Energy Remapping (Total energy in Lagrangan coordinates):

 !!!reversing the k-loop

       do j=1,np
          do i=1,np
             do k=1,nlev
               dpz(k,i,j) = dprs(i,j,k)
                tz(k,i,j) = pt3d(i,j,k)
                qz(k,i,j) = qt3d(i,j,k)
                uz(k,i,j) = uv3d(i,j,1,k)
                vz(k,i,j) = uv3d(i,j,2,k)
             end do
          end do
       end do

   !Pressure at interfaces from pressure-thickness (Lagrangian) on reversed loop 

    pr3d(1,:,:) = ptop(:,:)

    do j=1,np
        do i=1,np
            do k=1,nlev
              pr3d(k+1,i,j) = pr3d(k,i,j) + dpz(k,i,j)
            end do
        end do
    end do

  !New reference Pressure at etk(k+1/2) levels (Eulerian interfaces)

  do j=1,np
    do i=1,np
        do k = 1, nlev+1
           pr3d_ref(k,i,j) = p_0 * ak(k) + bk(k)* pr3d(nlev+1,i,j)
        enddo
    enddo
  enddo

   !New Pressure thickness defined  (Eulerian)

    do j=1,np
        do i=1,np
          do k = 1, nlev
            dprs(i,j,k) =  pr3d_ref(k+1,i,j) - pr3d_ref(k,i,j)
          enddo
       enddo
    enddo

  !Remapping onto the new reference coordinates
    do j=1,np
          do i=1,np

            do k=1,nlev
                pgrid(k) = (pr3d(k+1,i,j) + pr3d(k,i,j))*0.5D0
                 pref(k) = (pr3d_ref(k,i,j) + pr3d_ref(k+1,i,j))*0.5D0
                   uk(k) = uz(k,i,j)
                   vk(k) = vz(k,i,j)
                   tk(k) = tz(k,i,j)
                   qk(k) = qz(k,i,j)
            end do

          !Pre-computation for parabolic remap

               pgl(:) = pr3d(:,i,j)
               pge(:) = pr3d_ref(:,i,j)
                dp(:) = dpz(:,i,j)
               ipr(:) = 0

            do k = 2, nlev
               ipr(k) = index_search(pgl,1,1,pge(k))
            enddo

             Call ppm_coef(dp,coef)

             Call ppm_remap(pge,pgl,uz(:,i,j),dp,ipr,coef,ru)
             Call ppm_remap(pge,pgl,vz(:,i,j),dp,ipr,coef,rv)
             Call ppm_remap(pge,pgl,tz(:,i,j),dp,ipr,coef,rt)
             Call ppm_remap(pge,pgl,qz(:,i,j),dp,ipr,coef,rq)

           ! Linear extrapolations for top/bot boundary regions

                pgrid(0) = pr3d(1,i,j)
                pgrid(nlev+1) = pr3d(nlev+1,i,j)

                uk(0) =  (3.0D0 *uk(1) - uk(2) )* 0.5D0
                uk(nlev+1) =  (3.0D0 *uk(nlev) - uk(nlev-1) )* 0.5D0
                vk(0) =  (3.0D0 *vk(1) - vk(2) )* 0.5D0
                vk(nlev+1) =  (3.0D0 *vk(nlev) - vk(nlev-1) )* 0.5D0

                tk(0) =  (3.0D0 *tk(1) - tk(2) )* 0.5D0
                tk(nlev+1) =  (3.0D0 *tk(nlev) - tk(nlev-1) )* 0.5D0

                qk(0) =  (3.0D0 *qk(1) - qk(2) )* 0.5D0
                qk(nlev+1) =  (3.0D0 *qk(nlev) - qk(nlev-1) )* 0.5D0

          do k = 1, nlev

              if ((k > 2).and.(k < nlev-1)) then      !parabolic remap

                uz(k,i,j) = ru(k)
                vz(k,i,j) = rv(k)
                tz(k,i,j) = rt(k)
                qz(k,i,j) = rq(k)

              else                                    !Linear remap

                ik = index_search(pgrid,0,1,pref(k))
                px = pref(k)
                p1 = pgrid(ik)
                p2 = pgrid(ik+1)

                f1 = uk(ik)
                f2 = uk(ik+1)

               uz(k,i,j) = Linear_Int(p1,p2,f1,f2,px)

                f1 = vk(ik)
                f2 = vk(ik+1)

               vz(k,i,j) = Linear_Int(p1,p2,f1,f2,px)

                p1 = log(pgrid(ik))
                p2 = log(pgrid(ik+1))
                f1 = tk(ik)
                f2 = tk(ik+1)
                px = log(px)

               tz(k,i,j) = Linear_Int(p1,p2,f1,f2,px)

                f1 = qk(ik)
                f2 = qk(ik+1)
               qz(k,i,j) = Linear_Int(p1,p2,f1,f2,px)

             endif

             enddo

       end do
    end do

 !reversing the k-loop (for remapped fields)

    do k=1,nlev
       do j=1,np
          do i=1,np
                pt3d(i,j,k) = tz(k,i,j)       !PT
                qt3d(i,j,k) = qz(k,i,j)       !Q
                uv3d(i,j,1,k) = uz(k,i,j)     !u1
                uv3d(i,j,2,k) = vz(k,i,j)     !u2
               thick(i,j,k) =   dprs(i,j,k)   !dP
             end do
          end do
      end do

end subroutine parabolic_remap
!=======================================================================================================!
!=======================================================================================================!
!   This is based on  conservative PPM interpolation without total energy constraint
!=======================================================================================================!
subroutine monotonic_remap(metinv,sgp,ptop,pt3d,uv3d,thick)

    real (kind=real_kind), intent(in)    :: metinv(2,2,np,np)
    real (kind=real_kind), intent(in)    :: sgp(np,np), ptop(np,np)
    real (kind=real_kind), intent(inout) :: pt3d(np,np,nlev),uv3d(np,np,2,nlev)   
    real (kind=real_kind), intent(inout) :: thick(np,np,nlev)     !thickness (dP or dGP)

    real (kind=real_kind), dimension(nlev+1,np,np) :: pr3d_ref, pr3d
    real (kind=real_kind), dimension(0:nlev+1)   :: pgrid, logp, pref, uk,vk,tk
    real (kind=real_kind), dimension(nlev,np,np) :: uz, vz, tz, dpz
    real (kind=real_kind), dimension(nlev) :: ru, rv, rt
    real (kind=real_kind), dimension(nlev+1) ::  pge,pgl

    real (kind=real_kind), dimension(np,np,nlev) :: energy_lag, energy_eul, dlnp, dprs
    real (kind=real_kind), dimension(np,np,nlev+1) :: pr_eta, gp_eta, lnpr,prka

    integer :: indx(nlev), ipr(nlev)
    integer :: i,j,k , ik

    real (kind=real_kind) ::  trm1,trm2,trm3,trm4, term
    real (kind=real_kind) ::  xx(4), yy(4)
    real (kind=real_kind) ::  rdcp, p1,p2,pp, px,f1,f2, cprd,c_pp
    real (kind=real_kind) ::  eta(nlev), etv(nlev), etk(nlev+1)
    real (kind=real_kind) ::  ak(nlev+1),bk(nlev+1)
    real (kind=real_kind) ::  dp(nlev), tza(nlev),coef(nlev,8), fo(nlev)

    real (kind=real_kind) ::  top,bot,temp, peta, e_pe,e_th,e_ke, pp_0, dgp3d ,u1,u2,v1,v2


   ! Reseting progostic variables at reference coordinates

           rdcp = r_d / c_p
           cprd = c_p / r_d

           pp_0 = (p_0)**rdcp
           c_pp = c_p / pp_0

      Call eta_levels(ak,bk,eta,etk)

 !Pressure at etk(k) levels (interfaces)

   pr_eta(:,:,1) = ptop(:,:)

 ! for dP-thickness 
    do k=1,nlev
       do j=1,np
          do i=1,np
            pr_eta(i,j,k+1) = pr_eta(i,j,k) + thick(i,j,k)
          end do
       end do
    end do

  do k=1,nlev
    do j=1,np
       do i=1,np
          dprs(i,j,k) = pr_eta(i,j,k+1) - pr_eta(i,j,k)
       end do
    end do
  end do

 !!!reversing the k-loop

       do j=1,np
          do i=1,np
             do k=1,nlev
               dpz(k,i,j) = dprs(i,j,k)
                tz(k,i,j) = pt3d(i,j,k)
                uz(k,i,j) = uv3d(i,j,1,k)
                vz(k,i,j) = uv3d(i,j,2,k)
             end do
          end do
       end do

   !Pressure at interfaces from pressure-thickness (Lagrangian) on reversed loop 

    pr3d(1,:,:) = ptop(:,:)

    do j=1,np
        do i=1,np
            do k=1,nlev
              pr3d(k+1,i,j) = pr3d(k,i,j) + dpz(k,i,j)
            end do
        end do
    end do

  !New reference Pressure at etk(k+1/2) levels (Eulerian interfaces)

  do j=1,np
    do i=1,np
        do k = 1, nlev+1
           pr3d_ref(k,i,j) = p_0 * ak(k) + bk(k)* pr3d(nlev+1,i,j)
        enddo
    enddo
  enddo

   !New Pressure thickness defined  (Eulerian)

    do j=1,np
        do i=1,np
          do k = 1, nlev
            dprs(i,j,k) =  pr3d_ref(k+1,i,j) - pr3d_ref(k,i,j)
          enddo
       enddo
    enddo

  !Remapping onto the new reference coordinates
    do j=1,np
          do i=1,np

            do k=1,nlev
                pgrid(k) = (pr3d(k+1,i,j) + pr3d(k,i,j))*0.5D0
                 pref(k) = (pr3d_ref(k,i,j) + pr3d_ref(k+1,i,j))*0.5D0
                   uk(k) = uz(k,i,j)
                   vk(k) = vz(k,i,j)
                   tk(k) = tz(k,i,j)
            end do

          !Pre-computation for parabolic remap

              !pgl(:) = log(pr3d(:,i,j))
              !pge(:) = log(pr3d_ref(:,i,j))
              pgl(:) = pr3d(:,i,j)
              pge(:) = pr3d_ref(:,i,j)
                dp(:) = dpz(:,i,j)
               ipr(:) = 1

            do k = 2, nlev
               ipr(k) = index_search(pgl,1,1,pge(k))
            enddo

           ! Call ppm_coef(dp,coef)

           ! Call ppm_remap(pge,pgl,uz(:,i,j),dp,ipr,coef,ru)
           ! Call ppm_remap(pge,pgl,vz(:,i,j),dp,ipr,coef,rv)
           ! Call ppm_remap(pge,pgl,tz(:,i,j),dp,ipr,coef,rt)

             Call monoton_hoint(pgl,uz(:,i,j),pge,ipr,ru)
             Call monoton_hoint(pgl,vz(:,i,j),pge,ipr,rv)
             Call monoton_hoint(pgl,tz(:,i,j),pge,ipr,rt)

           ! Linear extrapolations for top/bot boundary regions

                pgrid(0) = pr3d(1,i,j)
                pgrid(nlev+1) = pr3d(nlev+1,i,j)

                uk(0) =  (3.0D0 *uk(1) - uk(2) )* 0.5D0
                uk(nlev+1) =  (3.0D0 *uk(nlev) - uk(nlev-1) )* 0.5D0
                vk(0) =  (3.0D0 *vk(1) - vk(2) )* 0.5D0
                vk(nlev+1) =  (3.0D0 *vk(nlev) - vk(nlev-1) )* 0.5D0

               tk(0) =  (3.0D0 *tk(1) - tk(2) )* 0.5D0
                tk(nlev+1) =  (3.0D0 *tk(nlev) - tk(nlev-1) )* 0.5D0

          do k = 1, nlev

              if ((k > 2).and.(k < nlev-1)) then      !parabolic remap

                uz(k,i,j) = ru(k)
                vz(k,i,j) = rv(k)
                tz(k,i,j) = rt(k)

              else                                    !Linear remap

                ik = index_search(pgrid,0,1,pref(k))
                px = pref(k)
                p1 = pgrid(ik)
                p2 = pgrid(ik+1)

                f1 = uk(ik)
                f2 = uk(ik+1)

               uz(k,i,j) = Linear_Int(p1,p2,f1,f2,px)

                f1 = vk(ik)
                f2 = vk(ik+1)

               vz(k,i,j) = Linear_Int(p1,p2,f1,f2,px)

                p1 = log(pgrid(ik))
                p2 = log(pgrid(ik+1))
                f1 = tk(ik)
                f2 = tk(ik+1)
                px = log(px)

               tz(k,i,j) = Linear_Int(p1,p2,f1,f2,px)

             endif

           ! if ((ie==88).and.(j==1).and.(i==3)) print*,k,ipr(k), pgl(k),  tz(k,i,j)

             enddo

       end do
    end do

 !reversing the k-loop (for remapped fields)

    do k=1,nlev
       do j=1,np
          do i=1,np
                pt3d(i,j,k) = tz(k,i,j)
                uv3d(i,j,1,k) = uz(k,i,j)
                uv3d(i,j,2,k) = vz(k,i,j)
               thick(i,j,k) =   dprs(i,j,k)    !!dP case
             end do
          end do
      end do


end subroutine monotonic_remap
!!=======================================================================================================!
!!R.Nair NCAR/IMAGe/2008
!!High-Order 1D monotone interpolation
!!Atmost Cubic or alteast quadraic/linear  monotone values depending on data
!!Inspired from the Akima's interpolation

       Subroutine   monoton_hoint(sgrid,fin,tgrid,inx,fout)

        Implicit None

           Real(kind=real_kind), Intent(in), dimension(nlev+1) :: sgrid, tgrid 
           Real(kind=real_kind), Intent(in), dimension(nlev) :: fin
           Integer, Intent(in), dimension(nlev)   :: inx
           Real(kind=real_kind), Intent(out), dimension(nlev)  :: fout

           Real(kind=real_kind), dimension(0:nlev+1) ::  sl,hx, yp,pl
           Real(kind=real_kind), dimension(nlev) ::  ak, bk
           Real(kind=real_kind) :: ai,bi,ci,di,xi, hr
           Integer ::  k, ik

            !fin   --> function value (known)
            !sgrid --> source grid  (non-uniform)
            !ik = inx(k) --> index
            !tgrid --> target grid,  s.t.,  sgrid(ik) < tgrid(k) < sgrid(ik+1)

           do k = 1, nlev
             hx(k) =  (sgrid(k+1) - sgrid(k))
           enddo

          !derivatives (slopes)
           do k = 1, nlev-1
             sl(k) = (fin(k+1)  - fin(k)) /hx(k)
           enddo

           do k = 2, nlev-1
             pl(k) = (sl(k-1)*hx(k) + sl(k)*hx(k-1)) /(hx(k) + hx(k-1))
           enddo

          !Monotonization
           do k = 2, nlev-1
             yp(k) = (sign(1.0D0,sl(k-1)) + sign(1.0D0,sl(k)) ) * &
                                            min(abs(sl(k-1)), abs(sl(k)),0.5D0*abs(pl(k)))
           enddo


           !Boundaries

            hr = hx(1) / (hx(1) + hx(2))
            pl(1) = sl(1) * (1.0D0 + hr) - sl(2) * hr 

              if  (pl(1)*sl(1) <= 0.0D0 ) then 
                  yp(1) = 0.0D0
              elseif (abs(pl(1)) > 2.0D0* abs(sl(1))) then 
                  yp(1) = 2.0D0 * sl(1)
              else
                  yp(1) = pl(1)
              endif

             hr = hx(nlev) / (hx(nlev) + hx(nlev-1))
             pl(nlev) = sl(nlev-1) * (1.0D0 + hr) - sl(nlev-2) * hr 

              if  (pl(nlev)*sl(nlev-1) <= 0.0D0 ) then
                  yp(nlev) = 0.0D0
              elseif (abs(pl(nlev)) > 2.0D0* abs(sl(nlev-1))) then
                  yp(nlev) = 2.0D0 * sl(nlev-1)
              else
                       yp(nlev) = pl(nlev)
              endif

          !! Cubic of the form:  fi(xo) = ai *(x-xo)**3 + bi*(x-xo)**2 + ci*(x-xo) + di

            do k = 1, nlev-1
                ak(k) = (yp(k) + yp(k+1) - 2.0D0* sl(k))/ (hx(k)*hx(k))
                bk(k) = (3.0D0*sl(k) - 2.0D0*yp(k) - yp(k+1)) / hx(k)
            enddo

            do k = 2, nlev-1
                ik = inx(k)
                ai = ak(ik)
                bi = bk(ik)
                ci = yp(ik)
                di = fin(ik)
                xi = tgrid(k) - sgrid(ik)
                fout(k) = xi * (xi*(ai * xi + bi) + ci) + di
            enddo

   End Subroutine   monoton_hoint

!=======================================================================================================!
!=======================================================================================================!
end module dg3d_remap_mod
