#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
!  This file contains the initial condititions for two baroclinic
!  instability probelms:
!
!  Polvani, Scott and Thomas, MWR 2004
!  Jablonowski and Williamson, QJR (2006) 132 
!


module baroclinic_inst_mod
!
!  This module contains the initial condititions for two baroclinic
!  instability probelms:
!
!  Polvani, Scott and Thomas, MWR 2004
!  Jablonowski and Williamson, QJR (2006) 132 
!
    ! ====================
    use kinds, only : real_kind, iulog
    ! ====================
    use physical_constants, only : omega, rearth, rgas, p0, dd_pi, Cp,g
    ! ====================
    use dimensions_mod, only : nlev,np, qsize
    ! ====================
    use quadrature_mod, only : quadrature_t, gauss
    ! ====================
    use element_mod, only : element_t
    ! ====================
    use global_norms_mod, only : global_integral
    ! ====================
    use hybrid_mod, only : hybrid_t
    ! ====================
    use control_mod, only : integration, test_case, u_perturb
    ! ====================
    use hybvcoord_mod, only : hvcoord_t 
    ! ====================
    use coordinate_systems_mod, only : spherical_polar_t
    use viscosity_mod, only : compute_zeta_C0


  implicit none
  private

  public :: binst_init_state, jw_baroclinic



contains


subroutine jw_baroclinic(elem,hybrid,hvcoord,nets,nete)
!=======================================================================================================!
    type(element_t), intent(inout) :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete
!=======================================================================================================!
    real (kind=real_kind)  :: tbar(nlev)     
    real (kind=real_kind) ::  eta(nlev), etv(nlev)
    real (kind=real_kind) ::  latc, lonc, rr,rc,  aa, lat,lon, snlat, cslat , v1,v2 
    real (kind=real_kind) ::  trm1,trm2,trm3,trm4, term
    integer :: tl,i,j,k,ie, idex
!=======================================================================================================!
   real (kind=real_kind),  parameter  :: u0 = 35.0D0        ! Zonal Mean wind 
   real (kind=real_kind),  parameter  :: t0 = 288.0D0       ! Mean temp
   real (kind=real_kind),  parameter  :: gama = 0.005D0     ! Lapse rate
!     Nair implementation was using r_d = 287.0D0, c_p = 1004.64D0
!     which is slightly different than Rgas and Cp.  should not matter
!   real (kind=real_kind),  parameter  :: r_d = 287.0D0      ! Gas const (dry)
!   real (kind=real_kind),  parameter  :: c_p = 1004.64D0    ! Cp        
!   real (kind=real_kind),  parameter  :: grv = 9.80616D0    ! Gravity   
!   real (kind=real_kind),  parameter  :: omg = 7.29212D-05  ! Omega     
!   real (kind=real_kind),  parameter  :: erad = 6.371229D06 ! Earth radius
   real (kind=real_kind),  parameter  :: ddt = 4.8D05       ! Temp  gradient 
   real (kind=real_kind),  parameter  :: eta_t  = 0.2D0     ! eta level 
   real (kind=real_kind),  parameter  :: eta_s  = 1.0D0     ! eta level 
   real (kind=real_kind),  parameter  :: eta_0  = 0.252D0   ! eta level 

   real (kind=real_kind) :: r_d,omg,grv,erad

   real(kind=real_kind), allocatable :: var3d(:,:,:,:)

   if (hybrid%masterthread) write(iulog,*) 'initializing Jablonowski and Williamson baroclinic instability test V1'

!      Call eta_levels(ak,bk,eta,etai)
!          interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
!          midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
!          etai = hyai(k)+hybi(k)
!          etam = hyam(k)+hybm(k)

! Note: for this test case, at t=0, ps0=ps=1000 mb.  
! so p = (a(k)+b(k)), OR  p = eta ps   and eta = sigma = p/ps
! (code is using mb for surface pressure?  )
!
!  
      r_d = Rgas
      omg=omega
      erad=rearth
      grv=g



       do k = 1, nlev
          eta(k) = hvcoord%etam(k)
          etv(k) = (eta(k) - eta_0) * DD_PI * 0.5D0 
       end do
       
       latc = DD_PI * (2.0D0/9.0D0) 
       lonc = DD_PI * (1.0D0/9.0D0) 


do ie=nets,nete

    ! initial velocity       
    do k=1, nlev 
       do j=1,np
          do i=1,np

          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          snlat=SIN(lat)
          cslat=COS(lat)


          aa = SIN(latc)*snlat + COS(latc)*cslat*COS(lon - lonc) 
          rc =  10.0D0  * ACOS(aa)
          v1 = u0 * (cos(etv(k)))**1.5D0 * (sin(2.0D0 * lat))**2  +  u_perturb*exp(-rc*rc) 

          elem(ie)%state%v(i,j,1,k,:)=v1
          elem(ie)%state%v(i,j,2,k,:)=0
 
          end do
       end do
    end do

 ! Layer-mean Temperature fields 

   do k = 1, nlev
    if (eta(k) <= eta_t) then 
       tbar(k) = t0 * eta(k)**(r_d*gama/grv) + ddt * (eta_t - eta(k))**5 
    else 
       tbar(k) = t0 * eta(k)**(r_d*gama/grv)
    endif
   end do


    do k=1,nlev
       do j=1,np
          do i=1,np
             lon = elem(ie)%spherep(i,j)%lon
             lat = elem(ie)%spherep(i,j)%lat

             snlat=SIN(lat)
             cslat=COS(lat)

             trm1 = 0.75D0 * (eta(k) * DD_PI*u0 /r_d) * sin(etv(k)) *sqrt(cos(etv(k)))
             trm2 = -2.0D0 *snlat**6 *(cslat**2 + 1.0D0/3.0D0) + 10.0D0/63.0D0 
             trm3 =  2.0D0 * u0* (cos(etv(k)))**1.5D0
             trm4 = (1.60D0 *cslat**3 *(snlat**2 + 2.0D0/3.0D0) - DD_PI *0.25D0)* erad*omg

             elem(ie)%state%T(i,j,k,:) = tbar(k) + trm1 *(trm2 * trm3 + trm4 )
 
          end do
       end do
    end do
 !Surface geopotential

    do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          snlat=SIN(lat)
          cslat=COS(lat)

          trm1 =  u0* ( cos((eta_s - eta_0)*DD_PI*0.5D0) )**1.5D0

          trm2 = -2.0D0 *snlat**6 *(cslat**2 + 1.0D0/3.0D0) + 10.0D0/63.0D0 
          trm3 = (1.60D0 *cslat**3 *(snlat**2 + 2.0D0/3.0D0) - DD_PI *0.25D0)* erad*omg

          elem(ie)%state%phis(i,j) = trm1 *(trm2 * trm1 + trm3) 
          elem(ie)%state%lnps(i,j,:) = LOG(p0)
          elem(ie)%state%ps_v(i,j,:) = p0
       end do
    end do
enddo


! initialize passive tracers.  just for testing.  
! first tracer:  temperature at T=0

if (qsize>=1) then
   do idex=1,qsize
   do ie=nets,nete
      !do tl=1,3
      do tl=1,1
         elem(ie)%state%Q(:,:,:,idex) = elem(ie)%state%T(:,:,:,tl)/400
      enddo
   enddo
   enddo
endif

! vorticity -2e-5 .. 2e-5
! total tracer mass is close to zero since Q is positive and negative
! so lets take truncated values to test qneg fixer
if (qsize>=2) then
   idex=2 ! prevents a compiler warning when qsize<2
   allocate(var3d(np,np,nlev,nets:nete))
   call compute_zeta_C0(var3d,elem,hybrid,nets,nete,1)
   do ie=nets,nete
      !do tl=1,3
      do tl=1,1
         do k=1,nlev
         do j=1,np
         do i=1,np
            elem(ie)%state%Q(i,j,k,idex) = var3d(i,j,k,ie)/2e-5
            if (var3d(i,j,k,ie)<0) elem(ie)%state%Q(i,j,k,idex)=0
         enddo
         enddo
         enddo
      enddo
   enddo
   deallocate(var3d)
endif


! third tracer is constant
if (qsize>=3) then
   idex=3
   do ie=nets,nete
      elem(ie)%state%Q(:,:,:,idex) = 1.0d0
   enddo
endif


! if we run a test case with 10 tracers, assume we want to reproduce the AVEC benchmarks
! using SJ's test tracers:
if (qsize==10) then
   do ie=nets,nete
      do j=1,np
      do i=1,np
         term = sin(9.*elem(ie)%spherep(i,j)%lon)*sin(9.*elem(ie)%spherep(i,j)%lat)

         do k=1,nlev
         do idex=4,qsize
            if ( term < 0. ) then
               elem(ie)%state%Q(i,j,k,idex) = 0
            else
               elem(ie)%state%Q(i,j,k,idex) = 1
            endif
         enddo
         enddo
      enddo
      enddo
   enddo
endif


!=======================================================================================================!
 end subroutine 



  subroutine binst_init_state(elem, hybrid,nets,nete,hvcoord)
    type(element_t), intent(inout) :: elem(:)
    integer, intent(in)   :: nets
    integer, intent(in)   :: nete
    type (hvcoord_t), intent(in) :: hvcoord
    type (hybrid_t), intent(in) :: hybrid

    ! Local variables

    integer, parameter :: nlat = 128

    integer, parameter :: nhl  = nlev+1, nfl = nlev

    type (quadrature_t) :: gs

    real (kind=real_kind), dimension(nlat) :: lat
    real (kind=real_kind), dimension(nlat) :: wts
#if 1
    real (kind=real_kind), parameter :: href = 7.34D0
#else
    real (kind=real_kind), parameter :: href = 7.D0
#endif 
    real (kind=real_kind), parameter :: hrgas = href/Rgas
    real (kind=real_kind), parameter :: z0 = 22.D0
    real (kind=real_kind), parameter :: z1 = 30.D0
    real (kind=real_kind), parameter :: zd = 5.D0
    real (kind=real_kind), parameter :: u0 = 50.D0
    real (kind=real_kind), parameter :: sat0 = 288.15D0

    real (kind=real_kind), dimension(nlat)   :: mu, mufac, tp, ef
    real (kind=real_kind), dimension(nlat-1) :: dmu
    real (kind=real_kind), dimension(nlat) :: ufull, dtdphi
    real (kind=real_kind),allocatable     :: t1(:,:,:)
    real (kind=real_kind), dimension(nfl) :: tfull, avg
    real (kind=real_kind) :: dlat
    real (kind=real_kind) :: z, fac
    real (kind=real_kind) :: tz, sz, fp, dudz, uf

    real (kind=real_kind), dimension(nlev+1) :: shalf, phalf
    real (kind=real_kind), dimension(nlev)   :: sfull, pfull, zfull

    real (kind=real_kind) :: v1, v2, sume, sumo
    real (kind=real_kind) :: latp,lonp
    real (kind=real_kind) :: tmp

    integer :: ie,i,j,k,l,ip,jp
    integer :: nm1 
    integer :: n0 
    integer :: np1

    integer :: nptsp,nptsv

    nm1= 1
    n0 = 2
    np1= 3

    if (hybrid%masterthread) write(iulog,*) 'initializing Jablonowski and Williamson ASP baroclinic instability test'

    do k=1,nhl
#if 1
       shalf(k) = (k-1)*1.0D0/REAL(nlev,kind=real_kind)
#else
       shalf(k) = (k-1)*0.05D0
#endif
       !      phalf(k) = shalf(k) * hvcoord%ps0 * 100.D0
    end do

    do k=1,nlev
       pfull(k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*p0
       zfull(k) = -href * log(pfull(k)/p0)
    end do

    !   sfull = ( shalf(2:nhl) + shalf(1:nhl-1) ) / 2.0D0
    !   pfull = sfull * p0 * 100.D0
    !   zfull = -href * log(sfull)

    do k=1,nfl
       if ( zfull(k) <= 11.D0 ) then
          tfull(k) = -6.5D0*zfull(k)
       else if ( zfull(k) <= 20.D0 ) then
          tfull(k) = -6.5D0*11.D0
       else if ( zfull(k) <= 32.D0 ) then
          tfull(k) = -6.5D0*11.D0+(zfull(k)-20.D0)
       else if ( zfull(k) <= 47.D0 ) then
          tfull(k) = -6.5D0*11.D0+(32.D0-20.D0)+2.8D0*(zfull(k)-32.D0)
       else if ( zfull(k) <= 51.D0 ) then
          tfull(k) = -6.5D0*11.D0+(32.D0-20.D0)+2.8D0*(47.D0-32.D0)
       else if ( zfull(k) <= 71.D0 ) then
          tfull(k) = -6.5D0*11.D0+(32.D0-20.D0)+2.8D0*(47.D0-32.D0)-2.8D0*(zfull(k)-51.D0)
       else if ( zfull(k) <= 80.D0 ) then
          tfull(k) = -6.5D0*11.D0+(32.D0-20.D0)+2.8D0*(47.D0-32.D0)-2.8D0*(71.D0-51.D0)- &
               2.D0*(zfull(k)-71.D0)
       else
#if 1
          tfull(k) = -6.5D0*11.D0+(32.D0-20.D0)+2.8D0*(47.D0-32.D0)-2.8D0*(71.D0-51.D0)- &
               2.D0*(84.852D0-71.D0)
#else
          tfull(k) = 1.0D+99       
#endif
       end if
       tfull(k) = tfull(k) + sat0
    end do

    nptsp = SIZE(elem(nets)%state%T(:,:,1,n0),1)
    nptsv = SIZE(elem(nets)%state%v(:,:,:,1,n0),1)
    allocate(t1(nptsp,nptsp,nets:nete))

    do ie=nets,nete
       elem(ie)%state%lnps(:,:,:) = LOG(p0)
       elem(ie)%state%ps_v(:,:,:) = p0
       !       elem(ie)%state%lnps(:,:,nm1)= elem(ie)%state%lnps(:,:,n0)
       !       elem(ie)%state%lnps(:,:,np1)= 0.0D0
    end do

    gs = gauss(nlat)

    do k=1,nlev
       t1(:,:,nets:nete) = 0.D0
       do ie=nets,nete

          z  = zfull(k)
          tz = tanh((z-z0)/zd)
          sz = 1.0D0/cosh((z-z0)/zd)
          fp = (-3.D0/10.D0) * (tz**2) * (sz**2) * sin(DD_PI*z/z1) + &
               (DD_PI/60.D0) * (1.D0 - tz**3) * cos(DD_PI*z/z1)

          do j=1,nptsp
             do i=1,nptsp

                lonp = elem(ie)%spherep(i,j)%lon
                latp = elem(ie)%spherep(i,j)%lat

                ! Gauss points and weights

                lat(:) = (latp+DD_PI*0.5D0)*0.5D0*gs%points(:) + (-DD_PI*0.5D0+latp)*0.5D0
                wts(:) = (latp+DD_PI*0.5D0)*0.5D0*gs%weights(:)

                do l=1,nlat
                   mu(l)   = sin(lat(l))
                   tp(l)   = tan(lat(l))
                   ef(l)   = 2.D0 * omega * mu(l)
                   !JMD The following expression was causing the IBM compiler to generate 
                   !    bad code.  So I removed it... 
                   !    mulfac(l) = (sin(DD_PI*mu(l)*mu(l)))**3
                   tmp     = sin(DD_PI*mu(l)*mu(l))
                   mufac(l) = (tmp)**3

                   fac = u0 * 0.5D0 * (1.D0-tanh((z-z0)/zd)**3) * sin(DD_PI*z/z1)
                   ufull(l) = fac * mufac(l)

                   dudz = u0 * fp * mufac(l)
                   dtdphi(l) = -hrgas*(rearth*ef(l)+2.D0*ufull(l)*tp(l))*dudz
                   if ( mu(l) < 0.0D0 ) then
                      dtdphi(l) = 0.0D0
                   end if
                end do

                tmp = 0.D0 
                do l=1,nlat
                   tmp = tmp + wts(l)*dtdphi(l)
                end do

                t1(i,j,ie) = tmp
                elem(ie)%state%T(i,j,k,n0) = t1(i,j,ie) 

             end do
          end do

       end do

       avg(k) = real(global_integral(elem, t1(:,:,nets:nete),hybrid,nptsp,nets,nete))

    end do
    if(test_case.eq."aquaplanet") then
       do k=1,nlev
          do ie=nets,nete
             elem(ie)%state%T(:,:,k,nm1)=elem(ie)%state%T(:,:,k,n0)
             elem(ie)%state%T(:,:,k,np1)=0.0D0
          end do
       end do
    else
       do k=1,nlev
          do ie=nets,nete
             do j=1,nptsp
                do i=1,nptsp

                   lonp = elem(ie)%spherep(i,j)%lon
                   latp = elem(ie)%spherep(i,j)%lat

                   elem(ie)%state%T(i,j,k,n0) = elem(ie)%state%T(i,j,k,n0) + tfull(k) - avg(k)

                   elem(ie)%state%T(i,j,k,n0)=elem(ie)%state%T(i,j,k,n0) + &
                        1.D0/cosh(3.D0*(lonp-DD_PI*0.5D0))**2 * &
                        1.D0/cosh(6.D0*(latp-DD_PI*0.25D0))**2

                end do
             end do

             elem(ie)%state%T(:,:,k,nm1)=elem(ie)%state%T(:,:,k,n0)
             elem(ie)%state%T(:,:,k,np1)=elem(ie)%state%T(:,:,k,n0)

          end do
       end do
    endif

    do ie=nets,nete
       do k=1,nlev

          z = zfull(k)

          do j=1,nptsv
             do i=1,nptsv

                lonp = elem(ie)%spherep(i,j)%lon
                latp = elem(ie)%spherep(i,j)%lat

                fac = u0 * 0.5D0 * (1.D0-tanh((z-z0)/zd)**3) * sin(DD_PI*z/z1)
                tmp = sin(DD_PI*sin(latp)*sin(latp))
                v1 = fac * (tmp)**3
                if (sin(latp) < 0) then
                   v1 = 0.0D0
                end if

                v2 = 0.0D0

#if 0
                if (( integration == "explicit" ).or.( integration == "full_imp" )) then
                   ! explicit covariant
                   elem(ie)%state%v(i,j,1,k,n0)= v1*elem(ie)%D(i,j,1,1) + v2*elem(ie)%D(i,j,2,1)
                   elem(ie)%state%v(i,j,2,k,n0)= v1*elem(ie)%D(i,j,1,2) + v2*elem(ie)%D(i,j,2,2)
                else
                   ! semi-implicit contravariant
                   elem(ie)%state%v(i,j,1,k,n0)= v1*elem(ie)%Dinv(i,j,1,1) + v2*elem(ie)%Dinv(i,j,1,2)
                   elem(ie)%state%v(i,j,2,k,n0)= v1*elem(ie)%Dinv(i,j,2,1) + v2*elem(ie)%Dinv(i,j,2,2)
                endif
#else
                elem(ie)%state%v(i,j,1,k,n0)= v1
                elem(ie)%state%v(i,j,2,k,n0)= v2
#endif
             end do
          end do

          elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,n0)
          elem(ie)%state%v(:,:,:,k,np1)=elem(ie)%state%v(:,:,:,k,n0)

       end do
    end do
    ! =======================================
    ! Initialize Surface Geopotential
    ! =======================================

    do ie=nets,nete
       elem(ie)%state%ps_v(:,:,:) =hvcoord%ps0
       !       elem(ie)%state%ps_v(:,:,nm1)=hvcoord%ps0
       !       elem(ie)%state%ps_v(:,:,np1)=0.0D0
       elem(ie)%state%phis(:,:) = 0.0D0
       elem(ie)%derived%fm = 0.0D0
    end do
    deallocate(t1)

  end subroutine binst_init_state





end module baroclinic_inst_mod



