#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define ZFROMGEO
#define FORCINGSTAT
module aquaplanet
#ifdef _PRIM
  use parallel_mod, only : abortmp, iam, mpireal_t, parallel_t, global_shared_buf, global_shared_sum
  ! ====================
  use kinds, only  : real_kind
  ! ==================== 
  use physical_constants, only : Cp, omega, rearth, rgas, kappa, G, p0, rwater_vapor, DD_PI
  ! ====================
  use dimensions_mod, only : nlev,np, nelemd
  ! ====================
  use element_mod, only : element_t
  ! ====================
  use physics_mod, only : elem_physics_t, Virtual_Temperature, Virtual_Specific_Heat, Saturation_Specific_Humidity
  ! ====================
  use hybrid_mod, only : hybrid_t
  ! ====================
  use hybvcoord_mod, only : hvcoord_t
  ! ====================
  use time_mod,      only : TimeLevel_t
  ! ====================
  use reduction_mod, only : parallelmax,parallelmin
  ! ====================
  use gravity_wave_drag_mod, only : Rayleigh, ue, ve, tme, the, pre, qve
  ! ====================
  use control_mod, only : statefreq, columnpackage
  ! ====================
  use column_types_mod, only : ColumnModelMulticloud_t
  ! ====================
  use global_norms_mod, only: wrap_repro_sum
  ! ====================
  implicit none
  private
  ! JPE aquaplanet namelist
  ! Aquaplanet intial state and forcings

  real (kind=real_kind), public :: cool_ampl
  real (kind=real_kind), public :: cool_min
  real (kind=real_kind), public :: cool_max

  integer              , public :: qv_flag
  integer              , public :: qv_pert_flag
  real (kind=real_kind), public :: qv_pert_ampl
  real (kind=real_kind), public :: qv_pert_zmin
  real (kind=real_kind), public :: qv_pert_zmax
  integer              , public :: isrf_forc
  real (kind=real_kind), public :: h_dis
  real (kind=real_kind), public :: cdrag
  real (kind=real_kind), public :: chdrag
  real (kind=real_kind), public :: wstar
  real (kind=real_kind), public :: tsurf
  real (kind=real_kind), public :: qsurf
  real (kind=real_kind), public :: u0
  real (kind=real_kind), public :: zabsampl
  real (kind=real_kind), public :: zabsmid
  real (kind=real_kind), public :: zabsmin
  integer              , public :: noisef



  logical :: exp_int=.true.

  public :: aquaplanet_init_state
  public :: aquaplanet_forcing
  public :: Potential_Temperature
  public :: presc_cooling_mc

  real (kind=real_kind),allocatable :: ff_min(:,:),ff_max(:,:)

  real (kind=real_kind) :: udrag_min,udrag_max
  real (kind=real_kind) :: vdrag_min,vdrag_max
  real (kind=real_kind) :: tsflx_min,tsflx_max
  real (kind=real_kind) :: qsflx_min,qsflx_max
  ! temporary - need to fix
  integer, parameter :: max_output_streams=5
  integer, parameter :: max_output_variables=15

  integer, parameter :: varcnt2d = 4, varcnt3d=4
  integer :: ivarID(varcnt2d+varcnt3d,max_output_streams)


  real (kind=real_kind), allocatable,public :: usf(:,:,:,:)
  real (kind=real_kind), allocatable,public :: vsf(:,:,:,:)
  real (kind=real_kind), allocatable,public :: tsf(:,:,:,:)
  real (kind=real_kind), allocatable,public :: qsf(:,:,:,:)

  real (kind=real_kind), allocatable,public :: udrag(:,:,:)
  real (kind=real_kind), allocatable,public :: vdrag(:,:,:)
  real (kind=real_kind), allocatable,public :: qsflx(:,:,:)
  real (kind=real_kind), allocatable,public :: tsflx(:,:,:)

contains

  subroutine aquaplanet_forcing(dt,ie, elemin,elemin_physics, hybrid,hvcoord,nets,nete,tl,mc)
    real (kind=real_kind),intent(in) :: dt
    type (element_t), intent(inout)  :: elemin
    type (elem_physics_t), intent(inout)  :: elemin_physics
    integer, intent(in)              :: ie
    type (hybrid_t), intent(in)      :: hybrid
    type (hvcoord_t), intent(in)     :: hvcoord
    type (TimeLevel_t), intent(in)   :: tl
    integer,intent(in)               :: nets,nete

    type (ColumnModelMulticloud_t), intent(in), optional :: mc

#ifdef FORCINGSTAT
    real (kind=real_kind) :: tcool_pmin,tcool_pmax,tcool_psum
#endif

    ! local
    real (kind=real_kind) :: cooling(nlev)
    integer :: nm1,iprint 
    integer :: i, j, k

    nm1   = tl%nm1
    iprint=0
!    if ((tl%nstep==0).and.(hybrid%masterthread).and.(ie==nets)) iprint=1

    ! cooling is embedded in MC param thus do not do it here...

    if(columnpackage .ne. "multicloud")then
       if(cool_ampl.ne.0.) then
          cooling=presc_cooling(hvcoord,np,nlev,iprint) ! cooling can change in time
          
          do k=1,nlev
             do j=1,np
                do i=1,np
                   elemin%derived%FT(i,j,k,nm1) = elemin%derived%FT(i,j,k,nm1) + cooling(k)*elemin_physics%mask(i,j)
                enddo
             enddo
          enddo
       endif

#ifdef FORCINGSTAT
!       if(ie==nete .and. MODULO(tl%nstep,statefreq) == 0) then
!          tcool_pmin = ParallelMin(minval(cooling),hybrid)
!          tcool_pmax = ParallelMax(maxval(cooling),hybrid)
!          tcool_psum = ParallelSum(sum(cooling),hybrid)
!          if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
!             write (*,100) "Tcool = ",tcool_pmin,tcool_pmax,tcool_psum
!          endif
!       endif
!100    format (A9,3(E24.15))
#endif
    endif

    if(isrf_forc.eq.1) then

       ! In the case of restart need to allocate these
       if(.not. allocated(usf)) then
          allocate(usf(np,np,nlev,nelemd))
          allocate(vsf(np,np,nlev,nelemd))
          allocate(qsf(np,np,nlev,nelemd))
          allocate(tsf(np,np,nlev,nelemd))
       end if


       call srf_flux_simple(usf,vsf,tsf,qsf,&
            ie,elemin,hybrid,hvcoord,nets,nete,tl,np,nlev,nelemd,iprint)

       do k=1,nlev
          do j=1,np
             do i=1,np
                elemin%derived%FM(i,j,1,k,nm1)= elemin%derived%FM(i,j,1,k,nm1) + usf(i,j,k,ie)
                elemin%derived%FM(i,j,2,k,nm1)= elemin%derived%FM(i,j,2,k,nm1) + vsf(i,j,k,ie)
             end do
          end do
       end do

       do k=1,nlev
          do j=1,np
             do i=1,np
                elemin%derived%FT(i,j,k,nm1)  = &
                     elemin%derived%FT(i,j,k,nm1)  + tsf(i,j,k,ie)
                elemin%derived%FQ(i,j,k,1,nm1)= &
                     elemin%derived%FQ(i,j,k,1,nm1)+ qsf(i,j,k,ie)
             end do
          end do
       end do
    endif

  end subroutine aquaplanet_forcing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function presc_cooling(hvcoord,npts,nlevels,iprint) result(cooling)

    integer, intent(in)          :: npts
    integer, intent(in)          :: nlevels
    integer, intent(in)          :: iprint
    type (hvcoord_t), intent(in) :: hvcoord
    real (kind=real_kind), dimension(nlevels) :: cooling

    ! Local variables

    real (kind=real_kind), parameter :: href = 7.34D3
    real (kind=real_kind), dimension(nlevels)   :: zfull
    real (kind=real_kind) :: cool_ampl_per_sec,func23,rnlev,a,b,c
    integer :: k

    func23(a,b,c)=max(0.0_real_kind,min(1.0_real_kind,(b-a)/(b-c)))

    zfull=zisothermal(href,nlevels)                       ! isothermal z heights
    rnlev=1.0_real_kind/(24.0_real_kind*3600.0_real_kind) ! [day]   -> [sec] 
    cool_ampl_per_sec = cool_ampl * rnlev                 ! [K/day] -> [K/sec]

    do k=1,nlevels
       cooling(k) = cool_ampl_per_sec*func23(zfull(k),cool_max,cool_min)
       if (iprint.eq.1) print 1,'cooling:',k,zfull(k), &
            cool_ampl,' [K]',cool_max,cool_min,cooling(k),' [K/sec]'
    end do
    if (iprint.eq.1) print *
1   format(a8,i3,f8.1,f6.2,a4,2f8.1,f12.8,a8)

  end function presc_cooling

  function presc_cooling_mc(hvcoord,mc,npts,nlevels,iprint) result(cooling)

    integer, intent(in)          :: npts
    integer, intent(in)          :: nlevels
    integer, intent(in)          :: iprint
    type (hvcoord_t), intent(in) :: hvcoord
    real (kind=real_kind), dimension(nlevels) :: cooling
    type (ColumnModelMulticloud_t), intent(in) :: mc

    ! Local variables

    real (kind=real_kind), parameter :: href = 7.34D3
    real (kind=real_kind), dimension(nlevels)   :: zfull
    real (kind=real_kind) :: cool_ampl_per_sec,func23,rnlev,a,b,c
    integer :: k

    do k=1,nlevels
       cooling(k) = mc%Q0R1*mc%D%psitrunc(k,1) + mc%D%Q0R2*(mc%D%psitrunc(k,2)-mc%csr*mc%D%psi2avgtrunc)
    end do

    if (iprint.eq.1)then
       do k=1,nlevels
          print 1,'cooling:',k,zfull(k), &
               cool_ampl,' [K]',cool_max,cool_min,cooling(k),' [K/sec]'
       end do
    endif

    if (iprint.eq.1) print *
1   format(a8,i3,f8.1,f6.2,a4,2f8.1,f12.8,a8)

  end function presc_cooling_mc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  subroutine aquaplanet_init_state(elem, hybrid,hvcoord,nets,nete, integration,noisestart)
  subroutine aquaplanet_init_state(elem, hybrid,hvcoord,nets,nete, integration)

    type (element_t), intent(inout) :: elem(:)
    type (hybrid_t),intent(in)       :: hybrid	
    type (hvcoord_t), intent(in) :: hvcoord
    integer, intent(in)   :: nets
    integer, intent(in)   :: nete
    character(len=*)    , intent(in) :: integration 
 !   integer, intent(in)   :: noisestart
    !ccccccccccccccccccccccccccccccccccccccccccccccccccc
    !cc environmental sounding:
    !cc "0000 UTC, 1 September 1974 GARP-GATE sounding"
    !cccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer, parameter :: npin = 23

    real(kind=real_kind) :: zin(npin)
    real(kind=real_kind) ::  th(npin)

    real(kind=real_kind), parameter :: temp(npin)=(/&
         25.26,  24.13,  21.04,  18.66,  16.50,  13.41,   9.06, &
         3.73,  -1.51,  -6.97, -14.09, -22.44, -30.57, -39.60, &
         -48.69, -57.40, -65.21, -72.58, -76.71, -74.98, -74.98, &
         -74.98, -74.98 /)

    real(kind=real_kind), parameter :: press(npin)=(/ &
         1008.00, 991.25, 945.50, 893.79, 836.06, 772.82, 705.22, &
         635.05, 564.48, 495.73, 430.71, 370.78, 316.72, 268.82, &
         226.98, 190.82, 159.87, 133.55, 111.29,  92.56,  52.31, &
         22.08,   9.32 /)

    real(kind=real_kind), parameter:: vap(npin)=(/ &
         0.178E+02, 0.172E+02, 0.156E+02, 0.134E+02, 0.111E+02, &
         0.888E+01, 0.631E+01, 0.487E+01, 0.396E+01, 0.200E+01, &
         0.984E+00, 0.806E+00, 0.370E+00, 0.135E+00, 0.599E-01, &
         0.258E-01, 0.123E-01, 0.582E-02, 0.367E-02, 0.589E-02, &
         0.104E-02, 0.247E-02, 0.585E-02 /)
    real(kind=real_kind), parameter ::  uu(npin)=(/ &
         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., &
         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)
    real(kind=real_kind), parameter ::  vv(npin)=(/ &
         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., &
         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)

    real(kind=real_kind) ::  aa,bb,cc,uur,uu0,uu1,vv1,rearthi,zearth,fcorio,tmpf
    real(kind=real_kind) ::  zabscal,zabsamp
    real(kind=real_kind) ::  u0loc
    real(kind=real_kind) ::  ff_pmin,ff_pmax

    !cccccccccccccccccccccccccccccccccccccccc
    !cc environmental sounding 
    !cc interpolated to model pressure levels
    !cccccccccccccccccccccccccccccccccccccccc
    integer, parameter :: nl = nlev
    integer, parameter :: interph = 1
    real (kind=real_kind), dimension(nlev) :: zfull,pfull
    real (kind=real_kind), dimension(np,np) :: ff

    !cccccccccccccccccccccccccccccccccccccccc
    ! Local variables and parameters
    !cccccccccccccccccccccccccccccccccccccccc
    integer :: ie,ne,i,j,k,kk,iisn,ifirst,ktop
    real (kind=real_kind)            :: coe2,pr0,ampl,fxy
    real (kind=real_kind), parameter :: tt00=273.15    ! temperature in [K] of zero [C]
    real (kind=real_kind), parameter :: href = 7.34D03 ! atmosphere stadart height [m]
    !   real (kind=real_kind), parameter :: href = 6.00D03 ! atmosphere stadart height [m]

    exp_int=.true.
    if(integration.eq."semi_imp") exp_int=.false.

    do k=1,npin
       th(k) =Potential_Temperature(temp(k),press(k),tt00)
    enddo
    !   zin=zhydrostatic(press,th  ,vap,press(1),th(1)  ,vap(1),0.0D0,npin,1,1,0.0D0)
    zin=zhydrostatic(press,temp,vap,press(1),temp(1),vap(1),0.0D0,npin,0,1,tt00)
    !----------------------------------------------------

    !   zin(npin)=37500.
    if (hybrid%masterthread) then
       print*,'INPUT SOUNDING'
       print*,' k    z,        p,      t,     th,',&
            '     qv,    u,      v'
       print*,' --------------------------------',&
            '-------------------------'
       do k=npin,1,-1
          print 921,k,zin(k),press(k),temp(k)+tt00,th(k),vap(k),uu(k),vv(k)
921       format(1x,i3,f9.2,4f8.2,2f6.2)
       enddo
       print *
    end if

    !----------------------------------------------------
    ! compute environmental profiles on state pressure half 
    ! levels (nlev+1 is surface), use pressure on half leveles
    ! phalf(:)=hvcoord%hyai(:)*hvcoord%ps0+hvcoord%hybi(:)*p0
    !----------------------------------------------------
    !cc surface data:
    !     iisn=1
    !     tme1(nlev+1)=temp(iisn)+tt00
    !     the1(nlev+1)=th(iisn)
    !     qve1(nlev+1)=vap(iisn)*1.e-3
    !      ue1(nlev+1)=uu(iisn)
    !      ve1(nlev+1)=vv(iisn)
    !
    !cc higher levels - interpolate:

#if 0
    !----------------------------------------------------
    ! linear interpolation between approximate heights
    ! compute approximate height assuming isothermal profile
    !----------------------------------------------------
    zfull=zisothermal(href,nlev)
    print *,'MAX HEIGHT:',maxval(zfull)
    do k=nlev,1,-1
       do kk=2,npin
          iisn=kk-1
          if(zin(kk).ge.zfull(k)) go to 665
       enddo
       print*,'INPUT SOUNDING DOES NOT GO HIGH ENOUGH. STOP.',k,zfull(k),kk,zin(kk)
       stop 'SOUNDING'
665    continue
       coe2=(zfull(k)-zin(iisn))/(zin(iisn+1)-zin(iisn))
       tme(k)=(coe2*temp(iisn+1) + (1.-coe2)*temp(iisn))+tt00
       qve(k)=(coe2*vap(iisn+1)  + (1.-coe2)*vap(iisn))*1.e-3
       the(k)= coe2*th(iisn+1)   + (1.-coe2)*th(iisn)
       ue(k) = coe2*uu(iisn+1)   + (1.-coe2)*uu(iisn)
       ve(k) = coe2*vv(iisn+1)   + (1.-coe2)*vv(iisn)
    end do
    pre(:)=hvcoord%hyam(:)*hvcoord%ps0+hvcoord%hybm(:)*p0
    zfull=zhydrostatic(pre,tme,qve*1.D3,press(1),temp(1),vap(1),0.0D0,nlev,0,-1,0.0D0)
    if (hybrid%masterthread) then
       print*,'ENVIRONMENTAL PROFILES ON STATE LEVELS, HEIGHT INTERPOLATION'
       do k=1,nlev
          print 200,zfull(k),pre(k),tme(k),qve(k)*1.e3,ue(k),ve(k)
       enddo
    end if
#else
    !--------------------------------------------------------
    ! linear interpolation between presssure levels
    ! compute approximate initial state pressure levels 
    ! using hybrid coordinates relationship, eq. (3.a.92) 
    ! of the CCM-2 description, (NCAR/TN-382+STR), p. 24.
    ! compute full levels (nlev is first level above surface)
    !--------------------------------------------------------
    pre(:)=hvcoord%hyam(:)*hvcoord%ps0+hvcoord%hybm(:)*p0

    !
    !    levels above sounding top are set to equal the sounding top
    !
    ktop=1 
    do k=1,nlev
       if(press(npin) > pre(k)) then
          ktop=k+1
       endif
    enddo

    iisn=npin
    do k=ktop,nlev
       do while(press(iisn).lt.pre(k) .and. iisn.gt.0)
          iisn=iisn-1
       enddo
       !!       coe2=(log(pre(k))-log(press(iisn)))/(log(press(iisn+1))-log(press(iisn)))
       coe2=(pre(k)-press(iisn))/(press(iisn+1)-press(iisn))
       tme(k)=(coe2*temp(iisn+1) + (1.-coe2)*temp(iisn))+tt00
       qve(k)=(coe2*vap(iisn+1)  + (1.-coe2)*vap(iisn))*1.e-3
       the(k)= coe2*th(iisn+1)   + (1.-coe2)*th(iisn)
       ue(k) = coe2*uu(iisn+1)   + (1.-coe2)*uu(iisn)
       ve(k) = coe2*vv(iisn+1)   + (1.-coe2)*vv(iisn)
    end do
    do k=1,ktop-1
       tme(k)=tme(ktop)
       qve(k)=qve(ktop)
       the(k)=the(ktop)
       ue(k) = ue(ktop)
       ve(k) = ve(ktop)
    enddo

    !       zfull=zhydrostatic(pre,tme,qve*1.D3,press(1),temp(1),vap(1),0.0D0,nlev,0,-1,0.0D0)
    zfull=zhydrostatic(pre,tme,qve*1.D3,p0,temp(1),vap(1),0.0D0,nlev,0,-1,0.0D0)
    if (hybrid%masterthread) then
       print*,'ENVIRONMENTAL PROFILES ON STATE LEVELS, PRESSURE INTERPOLATION'
       do k=1,nlev
          print 200,zfull(k),pre(k),tme(k),qve(k)*1.e3,ue(k),ve(k)
       enddo
       print *
    end if
#endif
200 format(1x,'z,p,tme,qve,ue,ve:',f9.2,3f7.2,2f6.2)
    !----------------------------------------------------
    !      zfull=zisothermal(href,nlev)
    !      print*,'ENVIRONMENTAL PROFILES ON STATE LEVELS, ISOTHERMAL HEIGHT TEST'
    !      do k=1,nlev
    !        print 200,zfull(k),pre(k),tme1(k),qve1(k)*1.e3,ue1(k),ve1(k)
    !      enddo
    !----------------------------------------------------

    qsurf=Saturation_Specific_Humidity(p0,tsurf)  ! p [mb]

    ! convert to potential temperature

    tsurf=tsurf*th(1)/(temp(1)+tt00)

    if(hybrid%masterthread) print *,'SURF:',tsurf,qsurf,th(1),temp(1)+tt00

    if(zabsampl.gt.0) then
       zabsamp=zabsampl/(24.0D0*3600.0D0)
       zabscal= (zabsmid-zabsmin)/DD_PI

       do k=1,nlev
          Rayleigh(k)=zabsamp*(1+tanh((zfull(k)-zabsmid)/zabscal))
       enddo
       if (hybrid%masterthread) then
          do k=1,nlev
             print 203,k,zfull(k),Rayleigh(k)
          enddo
       endif
203    format(1x,'Rayleigh friction :',i3,f10.3,2x,4e14.8)
    endif

    if(qv_flag.eq.0) then
       !----------------------------------------------------
       ! apply initial soundings to full set model variables
       !----------------------------------------------------
       ff=0.0D0
       allocate(ff_min(nets:nete,1:nlev),ff_max(nets:nete,1:nlev))
       ff_min(:,:)=100.0D0
       ff_max(:,:)=-100.0D0

       fcorio=2.0D0*omega
       rearthi=1./rearth

       ifirst=1
       do ie=nets,nete

          elem(ie)%state%lnps(:,:,:) = LOG(p0)
          elem(ie)%state%phis(:,:)   = 0.0D0

          do k=1,nlev

             if(noisef.gt.0) call noise(elem(ie)%spherep, ff,zfull(k),hybrid,np,np,k,ifirst)

             do j=1,np
                do i=1,np
                   ff_min(ie,k)=min(ff_min(ie,k),ff(i,j))
                   ff_max(ie,k)=max(ff_max(ie,k),ff(i,j))
                enddo
             enddo

             if(ie==nete .and. noisef.gt.0) then
                ff_pmin = ParallelMin(ff_min(:,k),hybrid)
                ff_pmax = ParallelMax(ff_max(:,k),hybrid)
                if (hybrid%masterthread) print *,'noise:',k,ff_pmin,ff_pmax
             endif

             ! elem(ie)%state%T(:,:,k,:)   =tme(k)
             !  elem(ie)%state%Q(:,:,1,k)   =(1.0D0+ff(:,:)*1.0D-2)*qve(k)
             !elem(ie)%state%Q(:,:,k,1,:)  = qve(k) 
             elem(ie)%state%Q(:,:,k,1)  = qve(k) !Q has one level

             zearth=rearth+zfull(k)
             uu0=u0*zearth*rearthi
             uur=u0*rearthi
             aa=uur*(2.*uu0+zearth*fcorio)+uu0*fcorio
             bb=uu0*(fcorio+uu0/zearth)
             cc=bb-.5*aa

             ! quick and dirty initialization of geopotential used only for first step forcing
             elem(ie)%derived%phi(:,:,k) = g*zfull(k)


             do j=1,np
                do i=1,np
                   elem(ie)%state%T(i,j,k,1)= tme(k)*(1.+cc/g*sin(elem(ie)%spherep(i,j)%lat)**2)
                   !  elem(ie)%state%T(i,j,k,1)= elem(ie)%state%T(i,j,k,1)+2.*ff(i,j) ! 2 degree noise
                   ! uu1=uu0*cos(elem(ie)%spherep(i,j)%lat)  ! *(ff(i,j)*0.01) ! 5% noise
                   elem(ie)%state%v(i,j,1,k,1)=uu0+ff(i,j)*1.0D-2
                   elem(ie)%state%v(i,j,2,k,1)=ve(k)

                end do
             end do
          end do

          elem(ie)%state%T(:,:,:,2)  = elem(ie)%state%T(:,:,:,1)
          elem(ie)%state%T(:,:,:,3)  = elem(ie)%state%T(:,:,:,1)
          elem(ie)%state%v(:,:,:,:,2)= elem(ie)%state%v(:,:,:,:,1)
          elem(ie)%state%v(:,:,:,:,3)= elem(ie)%state%v(:,:,:,:,1)
       end do
       deallocate(ff_min,ff_max)

    else ! QV only, initial H-S/baroclinic inst
       !----------------------------------------------------
       ! apply initial soundings to humidity variable only
       !----------------------------------------------------
       do ie=nets,nete
          elem(ie)%state%lnps(:,:,:) = LOG(p0)
          elem(ie)%state%phis(:,:)   = 0.0D0
          do k=1,nlev
             !elem(ie)%state%Q(:,:,k,1,:)   =qve(k)
             elem(ie)%state%Q(:,:,k,1)   =qve(k)
             elem(ie)%state%T(:,:,k,:)   =tme(k)
          end do
       end do
    endif

    if(qv_pert_flag.eq.1) then
       if (hybrid%masterthread) then
          print *, 'Applying initial Qv perterbation'
       end if

       !----------------------------------------------------
       ! apply the initial perturbation to humidity variable
       !----------------------------------------------------
       do ie=nets,nete
          do k=1,nlev
             if(zfull(k) >= qv_pert_zmin .and. zfull(k) <= qv_pert_zmax) then
                ampl = qv_pert_ampl
             else
                ampl = 0.0D0
             endif
             if(ie==nete) then
                if (hybrid%masterthread) print *,'Q pert:',k,ampl
             endif

             do j=1,np
                do i=1,np
                   if(abs(elem(ie)%spherep(i,j)%lat)< DD_PI/2.) then
                      fxy = cos(elem(ie)%spherep(i,j)%lon) * &
                           cos(elem(ie)%spherep(i,j)%lat) * &
                           cos(elem(ie)%spherep(i,j)%lat)
                   else
                      fxy = 0.0
                   end if
                   elem(ie)%state%Q(i,j,k,1)=(1.0+ampl*fxy)*elem(ie)%state%Q(i,j,k,1)

                   !                     if(k.eq.nlev .and. elem(ie)%spherep(i,j)%lat<(-0.21).and. &
                   !                          elem(ie)%spherep(i,j)%lat>(-0.22) .and. &
                   !                          elem(ie)%spherep(i,j)%lon<(3.76) .and.  elem(ie)%spherep(i,j)%lon>(3.75)) then
                   !                        print *, __FILE__,__LINE__,iam,i,j, elem(ie)%spherep(i,j)%lat, elem(ie)%spherep(i,j)%lon,elem(ie)%state%Q(i,j,k,1) 
                   !                     end if
                end do
             end do
          end do
       end do
    endif
    do ie=nets,nete
       elem(ie)%derived%FM  = 0.
       elem(ie)%derived%FT  = 0.
       elem(ie)%derived%FQ(:,:,:,1,:)  = 0.
    end do

    if(.not. allocated(usf)) then
       ne = SIZE(elem)
       allocate(usf(np,np,nlev,ne))
       allocate(vsf(np,np,nlev,ne))
       allocate(qsf(np,np,nlev,ne))
       allocate(tsf(np,np,nlev,ne))

    end if

    if(.not. allocated(udrag)) then
       allocate(udrag(np,np,nets:nete))
       allocate(vdrag(np,np,nets:nete))
       allocate(qsflx(np,np,nets:nete))
       allocate(tsflx(np,np,nets:nete))
       udrag=0.0_real_kind
       vdrag=0.0_real_kind
       qsflx=0.0_real_kind
       tsflx=0.0_real_kind
    endif




  end subroutine aquaplanet_init_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine srf_flux_simple(uf,vf,tf,qf,ie,elemin,&
       hybrid,hvcoord,nets,nete,tl,npts,nlevels,ne,iprint)
    implicit none	 
    real (kind=real_kind),intent(out):: uf(npts,npts,nlevels,ne)
    real (kind=real_kind),intent(out):: vf(npts,npts,nlevels,ne)
    real (kind=real_kind),intent(out):: tf(npts,npts,nlevels,ne)
    real (kind=real_kind),intent(out):: qf(npts,npts,nlevels,ne)
    type (element_t)                   :: elemin
    type (hybrid_t)                    :: hybrid	
    type (hvcoord_t)                   :: hvcoord
    type (TimeLevel_t)                 :: tl
    integer,target                     :: nets,nete
    integer,intent(in)                 :: ie,ne
    integer,intent(in)                 :: npts
    integer,intent(in)                 :: nlevels
    integer,intent(in)                 :: iprint

    real (kind=real_kind) :: zfull(np,np,nlev)
    real (kind=real_kind) :: pfull(np,np,nlev)
    real (kind=real_kind) :: utmp(np,np,nlev)
    real (kind=real_kind) :: vtmp(np,np,nlev)
    real (kind=real_kind) :: ttmp(np,np,nlev)
    real (kind=real_kind) :: qtmp(np,np,nlev)
    real (kind=real_kind) :: u_srf(np,np)
    real (kind=real_kind) :: v_srf(np,np)
#ifdef FORCINGSTAT
    real (kind=real_kind) :: udrag_pmin,udrag_pmax,udrag_psum
    real (kind=real_kind) :: vdrag_pmin,vdrag_pmax,vdrag_psum
    real (kind=real_kind) :: tsflx_pmin,tsflx_pmax,tsflx_psum
    real (kind=real_kind) :: qsflx_pmin,qsflx_pmax,qsflx_psum
    real (kind=real_kind) :: usfrc_pmin,usfrc_pmax,usfrc_psum
    real (kind=real_kind) :: vsfrc_pmin,vsfrc_pmax,vsfrc_psum
    real (kind=real_kind) :: tsfrc_pmin,tsfrc_pmax,tsfrc_psum
    real (kind=real_kind) :: qsfrc_pmin,qsfrc_pmax,qsfrc_psum
#endif
    !
    !    ! Local variables
    !
    integer i,j,k,nm1, kmin
    integer iee
    real (kind=real_kind) :: p_srf, t_srf, q_srf, rec_dz
    real (kind=real_kind) :: tsrf,qsrf, wstar2
    real (kind=real_kind) :: vau,umod,uprm,Wd,Tp,Qp,expfun
    real (kind=real_kind) :: Tvirt,Cpvirt,Rho,Rhoi_dz,G_dp,ees


    tsrf=tsurf
    qsrf=qsurf
    wstar2=wstar**2

    if(iprint.eq.1) print *,'SURF:',tsrf,qsrf


    nm1 = tl%nm1 

    !---------------------------------------------------
    ! compute pressure and heights from hydrostatic eq.
    !--------------------------------------------------- 
#ifdef ZFROMGEO
    zfull = elemin%derived%phi/g
#else

    do j=1,np
       do i=1,np
          p_srf=EXP(elemin%state%lnps(i,j,nm1))
          !jpe             p_srf=p0

          do k=1,nlevels
             pfull(i,j,k)  = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*p_srf


          end do

          zfull(i,j,:)=zhydrostatic(pfull(i,j,:), &
               elemin%state%T(i,j,:,nm1),&
               elemin%state%Q(i,j,:,1)*1.D3, &
               p_srf, &
               tsrf,&
               qsrf*1.D3, &
               0.0D0,nlevels,0,-1,0.0D0)

       end do
    end do

#endif

    !---------------------------------------------------
    ! compute surface fluxes from bulk drag parametrization
    !---------------------------------------------------

    do j=1,np
       do i=1,np
          u_srf(i,j)=elemin%state%v(i,j,1,nlevels,nm1)
          v_srf(i,j)=elemin%state%v(i,j,2,nlevels,nm1)
       end do
    end do


    Wd=0_real_kind
    Tp=0_real_kind
    Qp=0_real_kind
    uprm=0_real_kind
    umod=uprm

    do j=1,np
       do i=1,np

          p_srf=hvcoord%hyam(nlevels)*hvcoord%ps0 + hvcoord%hybm(nlevels)*EXP(elemin%state%lnps(i,j,nm1))

          t_srf = Potential_Temperature( &
               elemin%state%T(i,j,nlevels,nm1), &
               p_srf,0.0D0)
          q_srf=elemin%state%Q(i,j,nlevels,1)

          !             Wd=pelem(ie)%surfc%Wd(i,j)      ! Incorporate coupling 
          !             Tp=pelem(ie)%surfc%Tprime(i,j)  ! parameters from column
          !             Qp=pelem(ie)%surfc%Qprime(i,j)  ! model scheme

          vau=Cdrag*sqrt(u_srf(i,j)**2+v_srf(i,j)**2)

          udrag(i,j,ie)=vau*u_srf(i,j)
          vdrag(i,j,ie)=vau*v_srf(i,j)

          umod=Cdrag* sqrt(u_srf(i,j)**2+v_srf(i,j)**2+wstar2+Wd**2)

          !            uprm=Cdrag*(sqrt(u_srf(i,j)**2+v_srf(i,j)**2+wstar2+Wd**2) &
          !                       -sqrt(u_srf(i,j)**2+v_srf(i,j)**2+wstar2))

          tsflx(i,j,ie)=-umod*(t_srf-tsrf)-uprm*Tp

          qsflx(i,j,ie)=-umod*(q_srf-qsrf)-uprm*Qp

       end do
    end do




    !---------------------------------------------------
    ! compute vertical distribution of surface fluxes
    !---------------------------------------------------
    kmin=nlevels
    do k=1,nlevels

       do j=1,np
          do i=1,np

             if(zfull(i,j,k).gt.h_dis) then
                expfun=0.
                uf(i,j,k,ie)=0._real_kind
                vf(i,j,k,ie)=0._real_kind
                tf(i,j,k,ie)=0._real_kind
                qf(i,j,k,ie)=0._real_kind
             else
                expfun=exp(-zfull(i,j,k)/h_dis)
                if(k.lt.kmin) kmin=k  ! This is the highest level (minimum k) at which forcing is applied
             endif
             utmp(i,j,k)=udrag(i,j,ie)*expfun
             vtmp(i,j,k)=vdrag(i,j,ie)*expfun
             ttmp(i,j,k)=tsflx(i,j,ie)*expfun
             qtmp(i,j,k)=qsflx(i,j,ie)*expfun

          enddo
       enddo
    enddo

    !---------------------------------------------------
    ! forces for the momentum, temperature and moisture
    !---------------------------------------------------
    if(kmin .eq. 1) then
       !print *, 'it appears that surface forcing is being applied at TOA'
       !call abortmp('it appears that surface forcing is being applied at TOA')
    end if
    do k=kmin,nlevels-1
       do j=1,np
          do i=1,np
             rec_dz=1._real_kind/(zfull(i,j,k-1)-zfull(i,j,k+1))
             uf(i,j,k,ie)= (utmp(i,j,k-1)-utmp(i,j,k+1))*rec_dz
             vf(i,j,k,ie)= (vtmp(i,j,k-1)-vtmp(i,j,k+1))*rec_dz
             tf(i,j,k,ie)=-(ttmp(i,j,k-1)-ttmp(i,j,k+1))*rec_dz
             qf(i,j,k,ie)=-(qtmp(i,j,k-1)-qtmp(i,j,k+1))*rec_dz
          enddo
       enddo
    enddo
    k=nlevels
    do j=1,np
       do i=1,np
          rec_dz=1._real_kind/(zfull(i,j,k-1)-zfull(i,j,k))
          uf(i,j,k,ie)= (utmp(i,j,k-1)-utmp(i,j,k))*rec_dz
          vf(i,j,k,ie)= (vtmp(i,j,k-1)-vtmp(i,j,k))*rec_dz
          tf(i,j,k,ie)=-(ttmp(i,j,k-1)-ttmp(i,j,k))*rec_dz
          qf(i,j,k,ie)=-(qtmp(i,j,k-1)-qtmp(i,j,k))*rec_dz
       end do
    end do

#ifdef FORCINGSTAT
    if(ie==nete .and. MODULO(tl%nstep,statefreq) == 0) then
       udrag_pmin = ParallelMin(minval(udrag),hybrid)
       vdrag_pmin = ParallelMin(minval(vdrag),hybrid)
       tsflx_pmin = ParallelMin(minval(tsflx),hybrid)*Cp
       qsflx_pmin = ParallelMin(minval(qsflx),hybrid)*2.53e6
       udrag_pmax = ParallelMax(maxval(udrag),hybrid)
       vdrag_pmax = ParallelMax(maxval(vdrag),hybrid)
       tsflx_pmax = ParallelMax(maxval(tsflx),hybrid)*Cp
       qsflx_pmax = ParallelMax(maxval(qsflx),hybrid)*2.53e6
       do iee = nets, nete
       global_shared_buf(iee,1:4) = 0.d0
       do j = 1, np
       do i = 1, np
         global_shared_buf(iee,1) = global_shared_buf(iee,1) + udrag(i,j,iee)
         global_shared_buf(iee,2) = global_shared_buf(iee,2) + vdrag(i,j,iee)
         global_shared_buf(iee,3) = global_shared_buf(iee,3) + tsflx(i,j,iee)
         global_shared_buf(iee,4) = global_shared_buf(iee,4) + qsflx(i,j,iee)
       enddo
       enddo
       enddo
       call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
       udrag_psum = global_shared_sum(1)
       vdrag_psum = global_shared_sum(2)
       tsflx_psum = global_shared_sum(3)*Cp
       qsflx_psum = global_shared_sum(4)*2.53d6

       if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
          write (*,100) "Udrag = ",udrag_pmin,udrag_pmax,udrag_psum
          write (*,100) "Vdrag = ",vdrag_pmin,vdrag_pmax,vdrag_psum
          write (*,100) "Tsflx = ",tsflx_pmin,tsflx_pmax,tsflx_psum
          write (*,100) "Qsflx = ",qsflx_pmin,qsflx_pmax,qsflx_psum
          print *
       endif

       usfrc_pmin = ParallelMin(minval(usf),hybrid)
       vsfrc_pmin = ParallelMin(minval(vsf),hybrid)
       tsfrc_pmin = ParallelMin(minval(tsf),hybrid)
       qsfrc_pmin = ParallelMin(minval(qsf),hybrid)
       usfrc_pmax = ParallelMax(maxval(usf),hybrid)
       vsfrc_pmax = ParallelMax(maxval(vsf),hybrid)
       tsfrc_pmax = ParallelMax(maxval(tsf),hybrid)
       qsfrc_pmax = ParallelMax(maxval(qsf),hybrid)
       do iee = nets, nete
       global_shared_buf(iee,1:4) = 0.d0
       do k = 1, nlev
       do j = 1, np
       do i = 1, np
         global_shared_buf(iee,1) = global_shared_buf(iee,1) + usf(i,j,k,iee)
         global_shared_buf(iee,2) = global_shared_buf(iee,2) + vsf(i,j,k,iee)
         global_shared_buf(iee,3) = global_shared_buf(iee,3) + tsf(i,j,k,iee)
         global_shared_buf(iee,4) = global_shared_buf(iee,4) + qsf(i,j,k,iee)
       enddo
       enddo
       enddo
       enddo
       call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
       usfrc_psum = global_shared_sum(1)
       vsfrc_psum = global_shared_sum(2)
       tsfrc_psum = global_shared_sum(3)
       qsfrc_psum = global_shared_sum(4)

       if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
          write (*,100) "Usf = ",usfrc_pmin,usfrc_pmax,usfrc_psum
          write (*,100) "Vsf = ",vsfrc_pmin,vsfrc_pmax,vsfrc_psum
          write (*,100) "Tsf = ",tsfrc_pmin,tsfrc_pmax,tsfrc_psum
          write (*,100) "Qsf = ",qsfrc_pmin,qsfrc_pmax,qsfrc_psum
          write (*,100) "TsfW = ",tsfrc_pmin*Cp, &
               tsfrc_pmax*Cp, &
               tsfrc_psum*Cp
          write (*,100) "QsfW = ",qsfrc_pmin*2.53e6, &
               qsfrc_pmax*2.53e6, &
               qsfrc_psum*2.53e6
          print *
       endif

100    format (A9,3(E24.15))
    endif
#endif

  end subroutine srf_flux_simple


  function zisothermal(href,nlevels) result(zfull)
    real (kind=real_kind), intent(in)          :: href
    integer, intent(in)                        :: nlevels
    real (kind=real_kind), dimension(nlevels)  :: zfull
    real (kind=real_kind), dimension(nlevels)  :: sfull
    !   real (kind=real_kind), dimension(nlevels)  :: pfull
    !   real (kind=real_kind), dimension(nlevels+1):: shalf
    !   real (kind=real_kind), dimension(nlevels+1):: phalf
    real (kind=real_kind)                      :: rnlev
    integer                                    :: k

    rnlev=1.0_real_kind/REAL(nlevels,kind=real_kind)
    !-> do k=1,nlevels+1
    !->    shalf(k) = (k-1)*rnlev
    !->    phalf(k) = shalf(k) * hvcoord%ps0 * 100.D0
    !-> end do
    do k=1,nlevels
       sfull(k) = (2*k-1)*rnlev/2.0_real_kind
       !->    pfull(k) = sfull(k) * p0 * 100.D0
       zfull(k) = -href * log(sfull(k)) ! isothermal assumption
    end do
  end function zisothermal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !------------------------------------------------------
  ! compute approximate level heights for hydrostaic profile,
  ! (assume temperature is constant inside levels)
  !------------------------------------------------------

  function zhydrostatic(press,th,vap,psrf,tsrf,qsrf,zsrf,nlevs,itemp,ivert,tt00) result(zin)
    integer, intent(in) :: nlevs
    integer, intent(in) :: itemp     ! th: (1)->pot temp, (0)->temp
    integer, intent(in) :: ivert     ! ivert>0 vertical index rise up
    ! ivert<0 vertical index rise down
    real (kind=real_kind), intent(in), dimension(nlevs):: press
    real (kind=real_kind), intent(in), dimension(nlevs):: th
    real (kind=real_kind), intent(in), dimension(nlevs):: vap
    real (kind=real_kind), intent(in)                  :: psrf
    real (kind=real_kind), intent(in)                  :: tsrf
    real (kind=real_kind), intent(in)                  :: qsrf
    real (kind=real_kind), intent(in)                  :: zsrf
    real (kind=real_kind), intent(in)                  :: tt00
    real (kind=real_kind), dimension(nlevs)            :: zin
    real (kind=real_kind) :: delt,tavi,deltz,tempk,tempkm
    integer :: k,km,k0,k1,k2

    !compute approximated height of sounding pressure levels:
    if(ivert.eq.1) then
       k0=1
       k1=2
       k2=nlevs 
    elseif(ivert.eq.-1) then
       k0=nlevs
       k1=nlevs-1
       k2=1
    endif
    if(press(k0).eq.psrf) then  ! surface level
       zin(k0)=zsrf
    else                        ! first level above surface
       tempk =(th(k0)+tt00)*(1.+.6e-3*vap(k0))*(itemp*(p0/press(k0))**(-kappa)+1-itemp)
       tempkm=(tsrf  +tt00)*(1.+.6e-3*qsrf   )*(itemp*(p0/psrf     )**(-kappa)+1-itemp)
       delt=tempk-tempkm
       if (delt.gt.1.e-4) then
          tavi=log(tempk/tempkm)/delt
       else
          tavi=1./tempk
       endif
       deltz=-Rgas/(tavi*g) * log(press(k0)/psrf)
       zin(k0)=zsrf+deltz ! level above surface
    endif

    if(itemp.eq.1) then ! potential temperature
       tempk =(th(k0 )+tt00)*(p0/press(k0 ))**(-kappa)*(1.+.6e-3*vap(k0 ))
       do k=k1,k2,ivert
          km=k-ivert
          tempkm=tempk
          tempk =(th(k )+tt00)*(p0/press(k ))**(-kappa)*(1.+.6e-3*vap(k ))

          delt=tempk-tempkm
          if (delt.gt.1.e-4) then
             tavi=log(tempk/tempkm)/delt
          else
             tavi=1./tempk
          endif
          deltz=-Rgas/(tavi*g) * log(press(k)/press(km))
          zin(k)=zin(km)+deltz
       end do
    else  ! temperature
       tempk = virtual_temperature(th(k0 )+tt00,vap(k0)*1.e-3)
       do k=k1,k2,ivert
          km=k-ivert
          tempkm=tempk
          tempk = virtual_temperature(th(k )+tt00,vap(k)*1.e-3)
          delt=tempk-tempkm
          if (delt.gt.1.e-4) then
             tavi=log(tempk/tempkm)/delt
          else
             tavi=1./tempk
          endif
          deltz=-Rgas/(tavi*g) * log(press(k)/press(km))
          zin(k)=zin(km)+deltz
       end do
    endif

  end function zhydrostatic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function Potential_Temperature(temp,press,tt00) result(th)
    real (kind=real_kind), intent(in) :: tt00
    real (kind=real_kind), intent(in) :: temp  ! [K] or [C]
    real (kind=real_kind), intent(in) :: press ! [mb]
    real (kind=real_kind)             :: th    ! [K]
    real (kind=real_kind)             :: pr0   ! [mb]
    integer                           :: k

    !convert from temperature (deg C or K) into potential temperature
    pr0=1./(p0)
    th=(temp+tt00)*(press*pr0)**(-kappa)

  end function Potential_Temperature


  subroutine noise(spherep, ff,z,hybrid,nx,ny,k,ifirst)
    use coordinate_systems_mod, only : spherical_polar_t

    type (spherical_polar_t), intent(in) :: spherep(:,:)
    real (kind=real_kind), intent(out), dimension(nx,ny) :: ff
    real (kind=real_kind), intent(in) :: z
    type (hybrid_t)                   :: hybrid
    double precision rand
    real (kind=real_kind) :: ha,amp,ranf
    integer i,j,k,nx,ny,nz
    integer ia,im,ic,ierr,iran,ifirst
    save iran

    amp=max(0.0_real_kind,(ha-z)/ha)
#ifdef DONTDO
#if defined(_AIX) || defined(_BGL)      
    call random_number(ff)

    do j=1,ny
       do i=1,nx
          if(i>noisef .and. j>noisef .and. i<(nx-noisef+1) .and. j<(ny-noisef+1)) then
             ff(i,j)=(ff(i,j)-0.5)*amp*cos(spherep(i,j)%lat)
          else
             ff(i,j)=0.0_real_kind
          endif
       enddo
    end do

    return
#endif
#endif
    ha=3.0D3
    ! Numerical Recipes in Fortran - Quick and Dirty Generators p.274-275
    !     im         ia         ic           overflow
    !   86436       1093      18254           2^27
    !  117128       1277      24749           2^28
    !  145800       3661      30809           2^29
    !  139968       3877      29573           2^30
    !  134456       8121      28411           2^31
    !  233280       9301      49297           2^32

    im=86436
    ia=1093
    ic=18254


    if(ifirst.eq.1) iran=1
    iran=iran+hybrid%par%rank+k*100

    amp=max(0.0_real_kind,(ha-z)/ha)
    do j=noisef,ny-noisef+1
       do i=noisef,nx-noisef+1
          iran=mod(iran*ia+ic,im)
          ranf=float(iran)/float(im)
          ff(i,j)=(ranf-0.5)*amp*cos(spherep(i,j)%lat)
       enddo
    end do

    ifirst=0
  end subroutine noise
#elif defined(_PRIMDG)
  ! temporary place holder 
  use kinds, only  : real_kind
  use physical_constants, only : Cp, omega, rearth, rgas, kappa, G, p0, rwater_vapor, DD_PI
  use dimensions_mod, only : nlev,np, nelemd
  implicit none
  private
  real (kind=real_kind), allocatable,public :: usf(:,:,:,:)
  real (kind=real_kind), allocatable,public :: vsf(:,:,:,:)
  real (kind=real_kind), allocatable,public :: tsf(:,:,:,:)
  real (kind=real_kind), allocatable,public :: qsf(:,:,:,:)
  real (kind=real_kind), allocatable,public :: udrag(:,:,:)
  real (kind=real_kind), allocatable,public :: vdrag(:,:,:)
  real (kind=real_kind), allocatable,public :: qsflx(:,:,:)
  real (kind=real_kind), allocatable,public :: tsflx(:,:,:)
  contains
#endif


  subroutine sst_profile(lat,profile, sst)
    real(kind=real_kind), intent(in) :: lat(:,:)
    real(kind=real_kind), intent(out) :: sst(:,:)
    integer, intent(in) :: profile
    real(kind=real_kind), parameter :: CtoK=273.15
    real(kind=real_kind) :: sterm
    integer :: i,j


    if(profile.eq.1) then
       do j=1,np
          do i=1,np
             if(abs(lat(i,j)) <= DD_PI*0.33333_real_kind) then
                sterm = sin(1.5_real_kind*lat(i,j))
                sst(i,j) = CtoK +  27.0_real_kind*(1.0_real_kind - sterm*sterm)
             else
                sst(i,j) = CtoK
             endif
          end do
       end do
    end if
  end subroutine sst_profile
end module aquaplanet
