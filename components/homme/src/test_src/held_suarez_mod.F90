#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module held_suarez_mod


  use coordinate_systems_mod, only: spherical_polar_t
  use dimensions_mod,         only: nlev,np,qsize
  use element_mod,            only: element_t
  use element_state,          only: timelevels
  use element_ops,            only: get_field, set_thermostate
  use hybrid_mod,             only: hybrid_t
  use hybvcoord_mod,          only: hvcoord_t
  use kinds,                  only: real_kind, iulog
  use physical_constants,     only: p0, kappa,g, dd_pi
  use physics_mod,            only: prim_condense
  use time_mod,               only: secpday

implicit none
private

  real (kind=real_kind), public, parameter :: sigma_b = 0.70D0
  real (kind=real_kind), public, parameter :: k_a     = 1.0D0/(40.0D0*secpday)
  real (kind=real_kind), public, parameter :: k_f     = 1.0D0/(1.0D0*secpday)
  real (kind=real_kind), public, parameter :: k_s     = 1.0D0/(4.0D0*secpday)
  real (kind=real_kind), public, parameter :: dT_y    = 60.0D0
  real (kind=real_kind), public, parameter :: dtheta_z= 10.0D0

  public :: hs_v_forcing
  public :: hs_T_forcing
  public :: hs0_init_state
  public :: hs_forcing

contains

  subroutine hs_forcing(elemin,hvcoord,nm1,nm1_Q,dt)

    type (element_t)                 :: elemin
    type (hvcoord_t)                 :: hvcoord
    integer                          :: nm1,nm1_Q  ! timelevel to use
    real (kind=real_kind)                   :: dt

    ! local
    real (kind=real_kind)                   :: pmid,r0,r1,dtf_q,dp,rdp,FQ
    real (kind=real_kind), dimension(np,np) :: psfrc 
    real (kind=real_kind)                   :: temperature(np,np,nlev)
    integer                                 :: i,j,k,q

    dtf_q = dt
    call get_field(elemin,'temperature',temperature,hvcoord,nm1,nm1_Q)
        
    do j=1,np
       do i=1,np
          psfrc(i,j) = (elemin%state%ps_v(i,j,nm1))
       end do
    end do

    elemin%derived%FT(:,:,:) = elemin%derived%FT(:,:,:) + &
         hs_T_forcing(hvcoord,psfrc(1,1),               &
         temperature,elemin%spherep,np, nlev)
    
    elemin%derived%FM(:,:,:,:)= elemin%derived%FM(:,:,:,:) + &
         hs_v_forcing(hvcoord,psfrc(1,1),& 
         elemin%state%v(1,1,1,1,nm1),np,nlev)

    if (qsize>=1) then
       ! HS with tracer  (Galewsky type forcing, with flux of  2.3e-5 kg/m^2/s
       ! MASS in kg/m^2   = < Q dp_in_Pa / g >   
       ! flux in kg/m^2/s = < FQ dp_in_Pa / g >   
       ! We want < FQ dp_in_Pa / g > = 2.3e-5  so:  FQ = 2.3e-5*g/dp_in_Pa 

       ! lowest layer thickness, in Pa
       dp = ( hvcoord%hyai(nlev+1) - hvcoord%hyai(nlev) ) + &
               ( hvcoord%hybi(nlev+1) - hvcoord%hybi(nlev) )*1000*100
       rdp = 1./ dp
       q=1
       do j=1,np
          do i=1,np
             FQ = rdp * g * 2.3E-5 * COS(elemin%spherep(i,j)%lat)**2
             elemin%derived%FQ(i,j,nlev,q) =elemin%derived%FQ(i,j,nlev,q)+FQ
          enddo
       enddo

       do j=1,np
          do i=1,np
             do k=1,nlev
                pmid = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*(elemin%state%ps_v(i,j,nm1))
                r0=elemin%state%Q(i,j,k,q)
                r1=r0
                call Prim_Condense(r1,temperature(i,j,k),pmid)
                elemin%derived%FQ(i,j,k,q) = elemin%derived%FQ(i,j,k,q) + &
                     (r1-r0)/(dtf_q)
             enddo
          enddo
       enddo
    endif

    
  end subroutine hs_forcing

  function hs_v_forcing(hvcoord,ps,v,npts,nlevels) result(hs_v_frc)

    integer, intent(in)               :: npts
    integer, intent(in)               :: nlevels
    type (hvcoord_t), intent(in)       :: hvcoord
    real (kind=real_kind), intent(in) :: ps(npts,npts)

    real (kind=real_kind), intent(in) :: v(npts,npts,2,nlevels)
    real (kind=real_kind)             :: hs_v_frc(npts,npts,2,nlevels)

    ! Local variables

    integer i,j,k
    real (kind=real_kind) :: k_v
    real (kind=real_kind) :: rps(npts,npts)
    real (kind=real_kind) :: p,etam
    do j=1,npts
       do i=1,npts
         rps(i,j)     = 1.0D0/ps(i,j)
       end do
    end do

    do k=1,nlevels
       do j=1,npts
          do i=1,npts
             p    = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
             etam      = hvcoord%hyam(k) + hvcoord%hybm(k)
             k_v = k_f*MAX(0.0_real_kind,(etam - sigma_b )/(1.0_real_kind - sigma_b))
             hs_v_frc(i,j,1,k) = -k_v*v(i,j,1,k)
             hs_v_frc(i,j,2,k) = -k_v*v(i,j,2,k)
          end do
       end do
    end do

  end function hs_v_forcing

  function hs_T_forcing(hvcoord,ps,T,sphere,npts,nlevels) result(hs_T_frc)

    integer, intent(in) :: npts
    integer, intent(in) :: nlevels

    type (hvcoord_t), intent(in)          :: hvcoord
    real (kind=real_kind), intent(in)    :: ps(npts,npts)
    real (kind=real_kind), intent(in)    :: T(npts,npts,nlevels)
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)

    real (kind=real_kind)             :: hs_T_frc(npts,npts,nlevels)

    ! Local variables

    real (kind=real_kind) :: p,logprat,pratk,Teq
    real (kind=real_kind) :: logps0,etam
    real (kind=real_kind) :: lat,snlat

    real (kind=real_kind) :: k_t(npts,npts)
    real (kind=real_kind) :: snlatsq(npts,npts)
    real (kind=real_kind) :: cslatsq(npts,npts)

    real (kind=real_kind) :: rps(npts,npts)

    real (kind=real_kind) :: rec_one_minus_sigma_b

    integer i,j,k

    logps0    = LOG(hvcoord%ps0 )

    do j=1,npts
       do i=1,npts
         snlat        = SIN(sphere(i,j)%lat)
         snlatsq(i,j) = snlat*snlat
         cslatsq(i,j) = 1.0D0 - snlatsq(i,j)
         rps(i,j)     = 1.0D0/ps(i,j)
       end do
    end do

    rec_one_minus_sigma_b = 1.0D0/(1.0D0 - sigma_b)

    do k=1,nlevels
       do j=1,npts
          do i=1,npts
             p         = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
             logprat   = LOG(p)-logps0
             pratk     = EXP(kappa*(logprat))
             etam      = hvcoord%hyam(k) + hvcoord%hybm(k)

             k_t(i,j)= k_a + (k_s-k_a)*cslatsq(i,j)*cslatsq(i,j)* &
                       MAX(0.0D0,(etam - sigma_b)/(1.0D0 - sigma_b))
             Teq     = MAX(200.0D0,(315.0D0 - dT_y*snlatsq(i,j) - dtheta_z*logprat*cslatsq(i,j))*pratk)

#if 0
             ! ======================================
             ! This is a smooth forcing 
             ! for debugging purposes only...
             ! ======================================

             k_t(i,j)= k_a 
             pratk   = EXP(0.081*(logprat))
             Teq     = (315.0D0 - dT_y*snlatsq(i,j))*pratk
#endif
             hs_T_frc(i,j,k)= -k_t(i,j)*(T(i,j,k)-Teq)
          end do
       end do
    end do
      
  end function hs_T_forcing

  subroutine hs0_init_state(elem, hybrid, hvcoord,nets,nete,Tinit)

    type(element_t),        intent(inout) :: elem(:)
    type(hybrid_t),         intent(in)    :: hybrid                   ! hybrid parallel structure
    type (hvcoord_t),       intent(in)    :: hvcoord
    integer,                intent(in)    :: nets
    integer,                intent(in)    :: nete
    real (kind=real_kind),  intent(in)    :: Tinit

    ! Local variables
    
    integer ie,i,j,k,q,tl
    integer :: nm1 
    integer :: n0 
    integer :: np1
    real (kind=real_kind) :: lat_mtn,lon_mtn,r_mtn,h_mtn,rsq,lat,lon
    real (kind=real_kind) :: temperature(np,np,nlev)

    if (hybrid%masterthread) write(iulog,*) 'initializing Held-Suarez primitive equations test'

    nm1= 1
    n0 = 2
    np1= 3

    do ie=nets,nete

       elem(ie)%state%ps_v(:,:,n0) =hvcoord%ps0
       elem(ie)%state%ps_v(:,:,nm1)=hvcoord%ps0
       elem(ie)%state%ps_v(:,:,np1)=hvcoord%ps0

       elem(ie)%state%phis(:,:)=0.0D0
#undef HS_TOPO1
#ifdef HS_TOPO1
       lat_mtn = dd_pi/6
       lon_mtn = 3*dd_pi/2
       r_mtn = dd_pi/9
       h_mtn = 4000
       do i=1,np
          do j=1,np
             lat = elem(ie)%spherev(i,j)%lat
             lon = elem(ie)%spherev(i,j)%lon
             rsq=MIN((lat-lat_mtn)**2 + (lon-lon_mtn)**2,R_mtn**2)
             elem(ie)%state%phis(i,j)=g*h_mtn*(1.0D0 - SQRT(rsq)/R_mtn)
          enddo
       enddo
#endif

       elem(ie)%state%v(:,:,:,:,n0) =0.0D0
       elem(ie)%state%v(:,:,:,:,nm1)=elem(ie)%state%v(:,:,:,:,n0)
       elem(ie)%state%v(:,:,:,:,np1)=elem(ie)%state%v(:,:,:,:,n0)

       if (qsize>=1) then
       q=1
       elem(ie)%state%Q(:,:,:,q) =0  ! moist HS tracer IC=0
       do q=2,qsize
          elem(ie)%state%Q(:,:,:,q) =temperature(:,:,:)/400
       enddo
       endif

       temperature(:,:,:)=Tinit
       do tl=1,timelevels
          call set_thermostate(elem(ie),temperature,hvcoord,n0,1)
       enddo
    end do

  end subroutine hs0_init_state


end module held_suarez_mod

