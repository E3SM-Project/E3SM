!-----------------------------------------------------------------------------------!
!MODULE SPELT_MOD---------------------------------------------------------CE-for FVM!
! SPELT_MOD File for the spelt project in HOMME                                     !
! Author: Christoph Erath                                                           !
! Date: 24.August 2011                                                              !
! MAIN module to run spelt on HOMME                                                 !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module spelt_mod

  use kinds, only : real_kind, int_kind
  use dimensions_mod, only: ne, nlev, ntrac, np, ntrac_d, nc, nhe, nip, nipm, nep
  use edge_mod, only : ghostBuffertr_t,edgebuffer_t,ghostbuffer_t
  use time_mod, only : timelevel_t
  use coordinate_systems_mod, only : spherical_polar_t, cartesian2D_t
  use element_mod, only : element_t, timelevels
  use hybrid_mod, only : hybrid_t

  implicit none
  private
  type (ghostBuffertr_t)                      :: cellghostbuf
  type (ghostbuffertr_t)                      :: factorR
  type (EdgeBuffer_t)                         :: edgeveloc
    
  type, public :: spelt_struct
    ! spelt tracer mixing ratio: (kg/kg)
    real (kind=real_kind)    :: c(1-nipm:nep+nipm,1-nipm:nep+nipm,nlev,ntrac_d,timelevels) 
    real (kind=real_kind)    :: psc(1-nipm:nep+nipm,1-nipm:nep+nipm)
!-----------------------------------------------------------------------------------!   
    real (kind=real_kind)    :: vn0(np,np,2,nlev), vn12(np,np,2,nlev) 
!-----------------------------------------------------------------------------------!
    ! define the departure grid (depends on the examples, usually the velocity)
    ! they are NOT stored for the Levels
    type (spherical_polar_t) :: asphere(nep,nep)     ! Spherical coordinates
    ! equidistant mesh spacing in each direction of the reference(!) element, depends on nep       
    real(kind=real_kind)     :: dab(1:nc)   ! cell size in alpha beta
    real(kind=real_kind)     :: drefx(1-nhe:nc+nhe,1-nhe:nc+nhe), drefy(1-nhe:nc+nhe,1-nhe:nc+nhe), pref(1:nc+1)
      
    real(kind=real_kind)     :: sga(1-nipm:nep+nipm,1-nipm:nep+nipm)     
    ! next maybe used just for test reasons
    real (kind=real_kind)    :: elem_mass
      real (kind=real_kind)   :: area_sphere(1:nc,1:nc) 
    real (kind=real_kind)    :: cstart(1:nc,1:nc)    ! for error analysis
!       real (kind=real_kind)    :: maxcfl(2,nlev)   
    !should be replace once the mapping comes from the reference element
    real (kind=real_kind)    :: Dinv(2,2,np,np)     ! Map vector field on the sphere to covariant v on cube
    real (kind=real_kind)    :: Ainv(2,2,nep,nep)
    
    !these should be outside the structure, only needed for positivity preserving
    real (kind=real_kind)    :: R(1-nhe:nc+nhe,1-nhe:nc+nhe,nlev,ntrac_d)
    real (kind=real_kind)    :: fluxhigh(nc,nc,nlev,ntrac_d,4)
    
    integer                  :: Facenumber
  end type spelt_struct
  
  public :: cellghostbuf, factorR, edgeveloc 
  public :: spelt_init1,spelt_init2, spelt_init3, spelt_mcgregordss,spelt_rkdss, spelt_grid_init
  public :: spelt_run, spelt_runair
  public :: cip_coeff, cip_interpolate, metric_term, cell_search, qmsl_cell_filter, cell_minmax, cip_cell_avr
  public :: spelt_runlimit, spelt_runpos
contains


subroutine spelt_run(elem,spelt,hybrid,deriv,tstep,tl,nets,nete)

  use derivative_mod, only : derivative_t
  ! ---------------------------------------------------------------------------------
  use edge_mod, only :  ghostVpack2d, ghostVunpack2d
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
  ! ---------------------------------------------------------------------------------
  use coordinate_systems_mod, only : spherical_to_cart, cart2cubedspherexy, ref2sphere, sphere2cubedsphere
  ! ------EXTERNAL----------------
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------  
  
  implicit none
  type (element_t), intent(inout)             :: elem(:)
  type (spelt_struct), intent(inout)          :: spelt(:)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  type (derivative_t), intent(in)             :: deriv           ! derivative struct
  real (kind=real_kind), intent(in)           :: tstep
  type (TimeLevel_t), intent(in)              :: tl              ! time level struct
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
 
  integer                                     :: i,j,k,ie,itr
        
  real (kind=real_kind)                       :: ff(nip,nip)
  real (kind=real_kind)                       :: cf(nip,nip,1-nhe:nc+nhe,1-nhe:nc+nhe)
  type (spherical_polar_t)                    :: dsphere1(1:nep,1:nep), dsphere2(1:nep,1:nep)
    
  real (kind=real_kind)                       :: slval(3), fluxval(nep,nep,2), flux(4)
  type (cartesian2D_t)                        :: dref1(1:nep,1:nep), dref2(1:nep,1:nep)
  real (kind=real_kind)                       :: sg1(1:nep,1:nep),sg2(1:nep,1:nep)  
  real (kind=real_kind)                       :: minmax(1-nhe:nc+nhe,1-nhe:nc+nhe,2)
  
  real (kind=real_kind)                       :: contrauv(1:nep,1:nep,2,3)
  
  integer                                     :: icell1(1:nep,1:nep), jcell1(1:nep,1:nep)     
  integer                                     :: icell2(1:nep,1:nep), jcell2(1:nep,1:nep)
  integer                                     :: face_nodep
  
  integer                                     :: icell, jcell
  real (kind=real_kind)                       :: dx, dy, dxyi, dt6, sg, sga
  type (cartesian2D_t)                        :: alphabeta
  real (kind=real_kind)                       :: tmp
  type (spherical_polar_t)                    :: tmpsphere
  
!   call t_startf('SPELT scheme') 
  
  dt6  = tstep/ 6.0D0
  do ie=nets,nete 
    do k=1, nlev
!       call solidbody_all(spelt(ie), dsphere1,dsphere2,contrauv,k) 
      !For Benchmark test
      call boomerang_all(spelt(ie), dsphere1,dsphere2,contrauv,k,tl%nstep)
      !For SE velocities
!       call spelt_dep_from_gll(elem(ie), deriv, spelt(ie)%asphere,dsphere1,0.5D0*tstep,tl,k)         
!       call spelt_dep_from_gll(elem(ie), deriv, spelt(ie)%asphere,dsphere2,tstep,tl,k)
!       call get_contravelocities(elem(ie),spelt(ie),contrauv, k,deriv)         
      !search has not to be done for all tracers!
      
      do j=1,nep
        do i=1,nep
          call cell_search(elem(ie),spelt(ie), dsphere1(i,j),icell1(i,j), jcell1(i,j),dref1(i,j),alphabeta, face_nodep)
          sg1(i,j)=metric_term(alphabeta)
          call cell_search(elem(ie),spelt(ie), dsphere2(i,j), icell2(i,j), jcell2(i,j),dref2(i,j),alphabeta,face_nodep)
          sg2(i,j)=metric_term(alphabeta)
        end do
      end do
      
      ! search of both point on the trajectory done
      do itr=1,ntrac
        do j=1-nhe,nc+nhe
          do i=1-nhe,nc+nhe
            icell=1+(i-1)*nipm
            jcell=1+(j-1)*nipm
            ff=spelt(ie)%c(icell:icell+nipm,jcell:jcell+nipm,k,itr,tl%n0)
            minmax(i,j,:)=cell_minmax(ff)
            call cip_coeff(spelt(ie)%drefx(i,j),spelt(ie)%drefy(i,j),ff,ff(2,2),cf(:,:,i,j))
          enddo
        enddo
        do j=1,nep
          do i=1,nep  
            sga=spelt(ie)%sga(i,j)
            slval(1)=spelt(ie)%c(i,j,k,itr,tl%n0)
 
            tmp=cip_interpolate(cf(:,:,icell1(i,j),jcell1(i,j)),dref1(i,j)%x,dref1(i,j)%y) 
            tmp=qmsl_cell_filter(icell1(i,j),jcell1(i,j),minmax,tmp)
            slval(2)=(sga/sg1(i,j))*tmp

            tmp=cip_interpolate(cf(:,:,icell2(i,j),jcell2(i,j)),dref2(i,j)%x,dref2(i,j)%y) 
            tmp=qmsl_cell_filter(icell2(i,j),jcell2(i,j),minmax,tmp)
            slval(3)=(sga/sg2(i,j))*tmp
 
            spelt(ie)%c(i,j,k,itr,tl%np1)=(sga/sg2(i,j))*tmp

!            if (mod(i,2)==1) then                   ! works only for nip=3!!!
!               fluxval(i,j,1) =  dt6 * contrauv(i,j,1,1)* (slval(1) + & 
!                      4.0D0 * slval(2) + slval(3) )     
!             endif
!             if (mod(j,2)==1) then            ! works only for nip=3!!!
!               fluxval(i,j,2) =  dt6 * contrauv(i,j,2,1)* (slval(1) + & 
!                        4.0D0 * slval(2) + slval(3) )                      
!             endif
            
            if (mod(i,2)==1) then                   ! works only for nip=3!!!
               fluxval(i,j,1) =  dt6 * (contrauv(i,j,1,1)* slval(1) + & 
                      4.0D0 * contrauv(i,j,1,2)*slval(2) + contrauv(i,j,1,3)*slval(3) )  
                         
             endif
             if (mod(j,2)==1) then            ! works only for nip=3!!!
               fluxval(i,j,2) =  dt6 * (contrauv(i,j,2,1)* slval(1) + & 
                        4.0D0 * contrauv(i,j,2,2)*slval(2) + contrauv(i,j,2,3)*slval(3) )                                
             endif
          end do
        end do 
        do jcell=1,nc
          do icell=1,nc          
              i=1+(icell-1)*nipm
              j=1+(jcell-1)*nipm
              
              sga=spelt(ie)%sga(i+1,j+1)   
              dx=spelt(ie)%dab(icell)   
              dy=spelt(ie)%dab(jcell)
              
              flux(1) = dy * (fluxval(i,j,1) + 4.0D0 * fluxval(i,j+1,1) + fluxval(i,j+2,1))/6.0D0  ! west
              flux(2) = dx * (fluxval(i,j,2) + 4.0D0 * fluxval(i+1,j,2) + fluxval(i+2,j,2))/6.0D0  ! south
              flux(3) = dy * (fluxval(i+2,j,1) + 4.0D0 * fluxval(i+2,j+1,1) + fluxval(i+2,j+2,1))/6.0D0 ! east
              flux(4) = dx * (fluxval(i+2,j+2,2) + 4.0D0 * fluxval(i+1,j+2,2) + fluxval(i,j+2,2))/6.0D0 ! north
              
              spelt(ie)%c(i+1,j+1,k,itr,tl%np1) = spelt(ie)%c(i+1,j+1,k,itr,tl%n0)/sga + &
                                        (flux(1) + flux(2) - flux(3) - flux(4) ) / (spelt(ie)%area_sphere(icell,jcell))
              spelt(ie)%c(i+1,j+1,k,itr,tl%np1)=spelt(ie)%c(i+1,j+1,k,itr,tl%np1)*sga                                              
          end do
        end do
      end do      
    end do
    call ghostVpack2d(cellghostbuf,spelt(ie)%c,nipm, nep,nlev,ntrac,0, tl%np1, timelevels,elem(ie)%desc)
  end do
!-----------------------------------------------------------------------------------! 
!   call t_startf('SPELT Communication') 
  call ghost_exchangeV(hybrid,cellghostbuf,nipm,nep,ntrac)
!   call t_stopf('SPELT Communication')
!-----------------------------------------------------------------------------------!  
!   call t_startf('SPELT Unpacking')  
  do ie=nets,nete
    call ghostVunpack2d(cellghostbuf,spelt(ie)%c,nipm, nep,nlev,ntrac,0, tl%np1, timelevels,elem(ie)%desc)
  end do
!   call t_stopf('SPELT Unpacking')
!   call t_stopf('SPELT scheme') 
  
end subroutine spelt_run


! POSITIVITY PRESERVING WORKS HERE
subroutine spelt_runpos(elem,spelt,hybrid,deriv,tstep,tl,nets,nete)

  use derivative_mod, only : derivative_t
  ! ---------------------------------------------------------------------------------
  use edge_mod, only :  ghostVpack2d, ghostVunpack2d
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
  ! ---------------------------------------------------------------------------------
  use coordinate_systems_mod, only : spherical_to_cart, cart2cubedspherexy, ref2sphere, sphere2cubedsphere
  ! ------EXTERNAL----------------
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------  
  use bndry_mod, only : ghost_exchangevfull
  use edge_mod, only : ghostbuffertr_t, ghostvpackR, ghostvunpackR, initghostbuffer,&
       freeghostbuffer, freeghostbuffertr
       
  implicit none
  type (element_t), intent(inout)             :: elem(:)
  type (spelt_struct), intent(inout)          :: spelt(:)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  type (derivative_t), intent(in)             :: deriv           ! derivative struct
  real (kind=real_kind), intent(in)           :: tstep
  type (TimeLevel_t), intent(in)              :: tl              ! time level struct
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
 
  integer                                     :: i,j,k,ie,itr
        
  real (kind=real_kind)                       :: ff(nip,nip), c_low(1:nc,1:nc)
  real (kind=real_kind)                       :: cf(nip,nip,1-nhe:nc+nhe,1-nhe:nc+nhe)
  type (spherical_polar_t)                    :: dsphere1(1:nep,1:nep), dsphere2(1:nep,1:nep)
    
  real (kind=real_kind)                       :: slval(3), fluxval(nep,nep,2), flux(4)
  real (kind=real_kind)                       :: fluxlowx(nc+1,nc), fluxlowy(nc,nc+1)
  real (kind=real_kind)                       :: contrauv(1:nep,1:nep,2,3)
  
  
  type (cartesian2D_t)                        :: dref1(1:nep,1:nep), dref2(1:nep,1:nep)
  real (kind=real_kind)                       :: sg1(1:nep,1:nep),sg2(1:nep,1:nep)  
  real (kind=real_kind)                       :: minmax(1-nhe:nc+nhe,1-nhe:nc+nhe,2)
  integer                                     :: icell1(1:nep,1:nep), jcell1(1:nep,1:nep)     
  integer                                     :: icell2(1:nep,1:nep), jcell2(1:nep,1:nep)
  
  integer                                     :: icell, jcell, jx, jy
  real (kind=real_kind)                       :: dxoy, dxyi, dt6, sg, sga
  type (cartesian2D_t)                        :: alphabeta
  real (kind=real_kind)                       :: dx, dy, tmp !, tmpR(3,3)
  
  real (kind=real_kind)                       :: c_min, c_max, cmo(4)
  
  integer                                     :: face_nodep
  
!   call t_startf('SPELT scheme') 
    
  dt6  = tstep/ 6.0D0
  do ie=nets,nete 
    spelt(ie)%R=0.0D0
    do k=1, nlev
      !FOR BENCHMARK TESTS
      !       call solidbody_all(spelt(ie), dsphere1,dsphere2,contrauv,k) 
      call boomerang_all(spelt(ie), dsphere1,dsphere2,contrauv,k,tl%nstep)
      !FOR USING THE SE VELOCITIES    
!       call spelt_dep_from_gll(elem(ie), deriv, spelt(ie)%asphere,dsphere1,0.5D0*tstep,tl,k)         
!       call spelt_dep_from_gll(elem(ie), deriv, spelt(ie)%asphere,dsphere2,tstep,tl,k)
!       call get_contravelocities(elem(ie),spelt(ie),contrauv, k,deriv)    
      !search has not to be done for all tracers!
      do j=1,nep
        do i=1,nep
          call cell_search(elem(ie),spelt(ie), dsphere1(i,j),icell1(i,j), jcell1(i,j),dref1(i,j),alphabeta,face_nodep)
          sg1(i,j)=metric_term(alphabeta)
          call cell_search(elem(ie),spelt(ie), dsphere2(i,j), icell2(i,j), jcell2(i,j),dref2(i,j),alphabeta, face_nodep)
          sg2(i,j)=metric_term(alphabeta)
        end do
      end do
      ! search of both point on the trajectory done
      do itr=1,ntrac
        do j=1-nhe,nc+nhe
          do i=1-nhe,nc+nhe
            icell=1+(i-1)*nipm
            jcell=1+(j-1)*nipm
            ff=spelt(ie)%c(icell:icell+nipm,jcell:jcell+nipm,k,itr,tl%n0)
            minmax(i,j,:)=cell_minmax(ff)
            call cip_coeff(spelt(ie)%drefx(i,j),spelt(ie)%drefy(i,j),ff,ff(2,2),cf(:,:,i,j))
          enddo
        enddo
        !reconstruction coefficients caculated done
        do j=1,nep
          do i=1,nep  
            sga=spelt(ie)%sga(i,j)
            slval(1)=spelt(ie)%c(i,j,k,itr,tl%n0)
 
            tmp=cip_interpolate(cf(:,:,icell1(i,j),jcell1(i,j)),dref1(i,j)%x,dref1(i,j)%y) 
            tmp=qmsl_cell_filter(icell1(i,j),jcell1(i,j),minmax,tmp)
            slval(2)=(sga/sg1(i,j))*tmp

            tmp=cip_interpolate(cf(:,:,icell2(i,j),jcell2(i,j)),dref2(i,j)%x,dref2(i,j)%y) 
            tmp=qmsl_cell_filter(icell2(i,j),jcell2(i,j),minmax,tmp)
            slval(3)=(sga/sg2(i,j))*tmp
            spelt(ie)%c(i,j,k,itr,tl%np1)=slval(3)

!            if (mod(i,2)==1) then                   ! works only for nip=3!!!
!               fluxval(i,j,1) =  dt6 * spelt(ie)%contrau(i,j,k)* (slval(1) + & 
!                      4.0D0 * slval(2) + slval(3) )
!             endif
!             if (mod(j,2)==1) then            ! works only for nip=3!!!
!               fluxval(i,j,2) =  dt6 * spelt(ie)%contrav(i,j,k)* (slval(1) + & 
!                      4.0D0 * slval(2) + slval(3) )    
!             endif  
            if (mod(i,2)==1) then                   ! works only for nip=3!!!
               fluxval(i,j,1) =  dt6 * (contrauv(i,j,1,1)* slval(1) + & 
                      4.0D0 * contrauv(i,j,1,2)*slval(2) + contrauv(i,j,1,3)*slval(3) )  
             
             endif
             if (mod(j,2)==1) then            ! works only for nip=3!!!
               fluxval(i,j,2) =  dt6 * (contrauv(i,j,2,1)* slval(1) + & 
                        4.0D0 * contrauv(i,j,2,2)*slval(2) + contrauv(i,j,2,3)*slval(3) )                                
             endif
          end do
        end do 

!!!! only needed for filter              
        do jcell=1,nc
          j=2+(jcell-1)*nipm
          jy=1+(jcell-1)*nipm
          do icell=1,nc
            dx=spelt(ie)%dab(icell)   
            dy=spelt(ie)%dab(jcell)
            i=2+(icell-1)*nipm
            jx=1+(icell-1)*nipm
            
            sga=spelt(ie)%sga(i,j)
            !high order flux
            spelt(ie)%fluxhigh(icell,jcell,k,itr,1) = &
                 dy * (fluxval(jx,jy,1) + 4.0D0 * fluxval(jx,jy+1,1) + fluxval(jx,jy+2,1)) / 6.0D0  ! west
            spelt(ie)%fluxhigh(icell,jcell,k,itr,2) = &
                 dx * (fluxval(jx,jy,2) + 4.0D0 * fluxval(jx+1,jy,2) + fluxval(jx+2,jy,2)) / 6.0D0  ! south
            spelt(ie)%fluxhigh(icell,jcell,k,itr,3) = &
                 dy * (fluxval(jx+2,jy,1) + 4.0D0 * fluxval(jx+2,jy+1,1) + fluxval(jx+2,jy+2,1)) / 6.0D0 ! east
            spelt(ie)%fluxhigh(icell,jcell,k,itr,4) = &
                 dx * (fluxval(jx+2,jy+2,2) + 4.0D0 * fluxval(jx+1,jy+2,2) + fluxval(jx,jy+2,2)) / 6.0D0 ! north
            
            tmp=-min(0.0D0,spelt(ie)%fluxhigh(icell,jcell,k,itr,1))-min(0.0D0,spelt(ie)%fluxhigh(icell,jcell,k,itr,2)) + &
                               max(0.0D0,spelt(ie)%fluxhigh(icell,jcell,k,itr,3))+max(0.0D0,spelt(ie)%fluxhigh(icell,jcell,k,itr,4))
            spelt(ie)%R(icell,jcell,k,itr)=0.0D0  
            if (tmp>1.0D-14) then
              spelt(ie)%R(icell,jcell,k,itr)=min(1.0D0,(spelt(ie)%c(i,j,k,itr,tl%n0)/sga)*spelt(ie)%area_sphere(icell,jcell)/tmp)  
            endif
                   
          end do
        end do
      end do !ntrac
    end do   !nlev
    call ghostVpackR(factorR, spelt(ie)%R,nhe,nc,nlev,ntrac,0,elem(ie)%desc)
  end do
       ! Anti diffusive flux are computed for each cell, done!
!   call t_startf('SPELT Communication1')      
  call ghost_exchangeV(hybrid,factorR,nhe,nc,ntrac)
!   call t_stopf('SPELT Communication1') 

  do ie=nets,nete
    call ghostVunpackR(factorR, spelt(ie)%R, nhe, nc,nlev,ntrac,0,elem(ie)%desc)
    do k=1, nlev
      do itr=1,ntrac
        !!!! only needed for filter     
        do j=1,nc !nep,2
          do i=1,nc !nep,2    
              dx=spelt(ie)%dab(i)   
              dy=spelt(ie)%dab(j)        
              icell=2+(i-1)*nipm
              jcell=2+(j-1)*nipm   
              
              sga=spelt(ie)%sga(icell,jcell)
              cmo=1.0D0 
                !west
                if(spelt(ie)%fluxhigh(i,j,k,itr,1)>0.0D0) then
                  cmo(1)=spelt(ie)%R(i-1,j,k,itr)  
                else
                  cmo(1)=spelt(ie)%R(i,j,k,itr)
                endif 
                !south                         
                if(spelt(ie)%fluxhigh(i,j,k,itr,2)>0.0D0) then
                  cmo(2)=spelt(ie)%R(i,j-1,k,itr)
                else
                  cmo(2)=spelt(ie)%R(i,j,k,itr)
                endif  
                !east  
                if(spelt(ie)%fluxhigh(i,j,k,itr,3)>=0.0D0) then
                  cmo(3)=spelt(ie)%R(i,j,k,itr)  
                else
                  cmo(3)=spelt(ie)%R(i+1,j,k,itr)
                endif 
                !north      
                if(spelt(ie)%fluxhigh(i,j,k,itr,4)>=0.0D0) then
                  cmo(4)=spelt(ie)%R(i,j,k,itr)  
                else
                  cmo(4)=spelt(ie)%R(i,j+1,k,itr) 
                endif   
                spelt(ie)%c(icell,jcell,k,itr,tl%np1) = spelt(ie)%c(icell,jcell,k,itr,tl%n0)/sga - &
                     (-cmo(1)*spelt(ie)%fluxhigh(i,j,k,itr,1) - cmo(2)*spelt(ie)%fluxhigh(i,j,k,itr,2) &
                      + cmo(3)*spelt(ie)%fluxhigh(i,j,k,itr,3) + cmo(4)*spelt(ie)%fluxhigh(i,j,k,itr,4) ) &
                      /(spelt(ie)%area_sphere(i,j)) 
                spelt(ie)%c(icell,jcell,k,itr,tl%np1)=spelt(ie)%c(icell,jcell,k,itr,tl%np1)*sga          
          end do
        end do                
      end do !ntrac
    end do   !nlev
    call ghostVpack2d(cellghostbuf,spelt(ie)%c,nipm, nep,nlev,ntrac,0, tl%np1, timelevels,elem(ie)%desc)
  end do
!-----------------------------------------------------------------------------------! 
!   call t_startf('SPELT Communication2') 
  call ghost_exchangeV(hybrid,cellghostbuf,nipm,nep,ntrac)
!   call t_stopf('SPELT Communication2')
!-----------------------------------------------------------------------------------!  
!   call t_startf('SPELT Unpacking2')  
  do ie=nets,nete
    call ghostVunpack2d(cellghostbuf,spelt(ie)%c,nipm, nep,nlev,ntrac,0, tl%np1, timelevels,elem(ie)%desc)
  end do
!   call t_stopf('SPELT Unpacking2')
!   call t_stopf('SPELT scheme')
end subroutine spelt_runpos


subroutine spelt_runair(elem,spelt,hybrid,deriv,tstep,tl,nets,nete)

  use derivative_mod, only : derivative_t
  ! ---------------------------------------------------------------------------------
  use edge_mod, only :  ghostVpack2d, ghostVunpack2d
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
  ! ---------------------------------------------------------------------------------
  use coordinate_systems_mod, only : spherical_to_cart, cart2cubedspherexy, ref2sphere, sphere2cubedsphere
  ! ------EXTERNAL----------------
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------  
  
  implicit none
  type (element_t), intent(inout)             :: elem(:)
  type (spelt_struct), intent(inout)          :: spelt(:)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  type (derivative_t), intent(in)             :: deriv           ! derivative struct
  real (kind=real_kind), intent(in)           :: tstep
  type (TimeLevel_t), intent(in)              :: tl              ! time level struct
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
 
  integer                                     :: i,j,k,ie,itr
        
  real (kind=real_kind)                       :: ff(nip,nip)
  real (kind=real_kind)                       :: cf(nip,nip,1-nhe:nc+nhe,1-nhe:nc+nhe)
  type (spherical_polar_t)                    :: dsphere1(1:nep,1:nep), dsphere2(1:nep,1:nep)
    
  real (kind=real_kind)                       :: slval(3), fluxval(nep,nep,2), flux(4)
  real (kind=real_kind)                       :: slvalair(nep,nep,3), slvalone(nep,nep,3)
  
  type (cartesian2D_t)                        :: dref1(1:nep,1:nep), dref2(1:nep,1:nep)
  real (kind=real_kind)                       :: sg1(1:nep,1:nep),sg2(1:nep,1:nep)  
  real (kind=real_kind)                       :: minmax(1-nhe:nc+nhe,1-nhe:nc+nhe,2)
  integer                                     :: icell1(1:nep,1:nep), jcell1(1:nep,1:nep)     
  integer                                     :: icell2(1:nep,1:nep), jcell2(1:nep,1:nep)
  real (kind=real_kind)                       :: contrauv(1:nep,1:nep,2,3)
  
  integer                                     :: icell, jcell
  real (kind=real_kind)                       :: dx, dy, dxyi, dt6, sg, sga
  type (cartesian2D_t)                        :: alphabeta
  real (kind=real_kind)                       :: tmp
  integer                                     :: face_nodep
  
  call t_startf('SPELT scheme') 
  
  dt6  = tstep/ 6.0D0
  do ie=nets,nete 
    do k=1, nlev
!       call solidbody_all(spelt(ie), dsphere1,dsphere2,contrauv,k) 
      !For Benchmark test
!       call boomerang_all(spelt(ie), dsphere1,dsphere2,contrauv,k,tl%nstep)
      !For SE velocities
      call spelt_dep_from_gll(elem(ie), deriv, spelt(ie)%asphere,dsphere1,0.5D0*tstep,tl,k)         
      call spelt_dep_from_gll(elem(ie), deriv, spelt(ie)%asphere,dsphere2,tstep,tl,k)
      call get_contravelocities(elem(ie),spelt(ie),contrauv, k,deriv)       
      !search has not to be done for all tracers!
      do j=1,nep
        do i=1,nep
          call cell_search(elem(ie),spelt(ie), dsphere1(i,j),icell1(i,j), jcell1(i,j),dref1(i,j),alphabeta,face_nodep)
          sg1(i,j)=metric_term(alphabeta) 
          call cell_search(elem(ie),spelt(ie), dsphere2(i,j), icell2(i,j), jcell2(i,j),dref2(i,j),alphabeta,face_nodep)
          sg2(i,j)=metric_term(alphabeta)          
        end do
      end do
      
      !ONE
!       do j=1-nhe,nc+nhe
!         do i=1-nhe,nc+nhe
!           icell=1+(i-1)*nipm
!           jcell=1+(j-1)*nipm
!           ff=spelt(ie)%c(icell:icell+nipm,jcell:jcell+nipm,k,1,tl%n0)
!           minmax(i,j,:)=cell_minmax(ff)
!           call cip_coeff(spelt(ie)%drefx(i,j),spelt(ie)%drefy(i,j),ff,ff(2,2),cf(:,:,i,j))
!         enddo
!       enddo
!       do j=1,nep
!         do i=1,nep  
!           sga=spelt(ie)%sga(i,j)
!           slvalone(i,j,1)=spelt(ie)%c(i,j,k,1,tl%n0)
! 
!           tmp=cip_interpolate(cf(:,:,icell1(i,j),jcell1(i,j)),dref1(i,j)%x,dref1(i,j)%y) 
!           tmp=qmsl_cell_filter(icell1(i,j),jcell1(i,j),minmax,tmp)
!           slvalone(i,j,2)=(sga/sg1(i,j))*tmp
! 
!           tmp=cip_interpolate(cf(:,:,icell2(i,j),jcell2(i,j)),dref2(i,j)%x,dref2(i,j)%y) 
!           tmp=qmsl_cell_filter(icell2(i,j),jcell2(i,j),minmax,tmp)
!           slvalone(i,j,3)=(sga/sg2(i,j))*tmp
!           spelt(ie)%c(i,j,k,1,tl%np1)=spelt(ie)%sga(i,j)  !slvalone(i,j,3)
! 
!         end do
!       end do 
      
      ! AIR
      do j=1-nhe,nc+nhe
        do i=1-nhe,nc+nhe
          icell=1+(i-1)*nipm
          jcell=1+(j-1)*nipm
          ff=spelt(ie)%c(icell:icell+nipm,jcell:jcell+nipm,k,1,tl%n0)
          minmax(i,j,:)=cell_minmax(ff)
          call cip_coeff(spelt(ie)%drefx(i,j),spelt(ie)%drefy(i,j),ff,ff(2,2),cf(:,:,i,j))
        enddo
      enddo
      do j=1,nep
        do i=1,nep  
          sga=spelt(ie)%sga(i,j)
          slvalair(i,j,1)=spelt(ie)%c(i,j,k,1,tl%n0)

          tmp=cip_interpolate(cf(:,:,icell1(i,j),jcell1(i,j)),dref1(i,j)%x,dref1(i,j)%y) 
!           tmp=qmsl_cell_filter(icell1(i,j),jcell1(i,j),minmax,tmp)
          slvalair(i,j,2)=(sga/sg1(i,j))*tmp

          tmp=cip_interpolate(cf(:,:,icell2(i,j),jcell2(i,j)),dref2(i,j)%x,dref2(i,j)%y) 
!           tmp=qmsl_cell_filter(icell2(i,j),jcell2(i,j),minmax,tmp)
          slvalair(i,j,3)=(sga/sg2(i,j))*tmp
          spelt(ie)%c(i,j,k,1,tl%np1)=slvalair(i,j,3)

!           if (mod(i,2)==1) then                   ! works only for nip=3!!!
!              fluxval(i,j,1) =  dt6 * (contrauv(i,j,1,1)* sga*slvalair(i,j,1)/slvalone(i,j,1) + & 
!                     4.0D0 * contrauv(i,j,1,2)*sga*slvalair(i,j,2)/slvalone(i,j,2) + &
!                     contrauv(i,j,1,3)*sga*slvalair(i,j,3)/slvalone(i,j,3) )  
!            endif
!            if (mod(j,2)==1) then            ! works only for nip=3!!!
!              fluxval(i,j,2) =  dt6 * (contrauv(i,j,2,1)* sga*slvalair(i,j,1)/slvalone(i,j,1) + & 
!                       4.0D0 * contrauv(i,j,2,2)*sga*slvalair(i,j,2)/slvalone(i,j,2) + &
!                       contrauv(i,j,2,3)*sga*slvalair(i,j,3)/slvalone(i,j,3) )                                
!            endif
            if (mod(i,2)==1) then                   ! works only for nip=3!!!
               fluxval(i,j,1) =  dt6 * (contrauv(i,j,1,1)* slvalair(i,j,1)+ & 
                      4.0D0 * contrauv(i,j,1,2)*slvalair(i,j,2)+ &
                      contrauv(i,j,1,3)*slvalair(i,j,3) )  
             endif
             if (mod(j,2)==1) then            ! works only for nip=3!!!
               fluxval(i,j,2) =  dt6 * (contrauv(i,j,2,1)* slvalair(i,j,1) + & 
                        4.0D0 * contrauv(i,j,2,2)*slvalair(i,j,2) + &
                        contrauv(i,j,2,3)*slvalair(i,j,3) )                                
             endif

! 
!          if (mod(i,2)==1) then                   ! works only for nip=3!!!
!             fluxval(i,j,1) =  dt6 * spelt(ie)%contrau(i,j,k)* sga*(slvalair(i,j,1)/slvalone(i,j,1) + & 
!                    4.0D0 * slvalair(i,j,2)/slvalone(i,j,2) + slvalair(i,j,3)/slvalone(i,j,3) )
!           endif
!           if (mod(j,2)==1) then            ! works only for nip=3!!!
!             fluxval(i,j,2) =  dt6 * spelt(ie)%contrav(i,j,k)* sga*(slvalair(i,j,1)/slvalone(i,j,1) + & 
!                    4.0D0 * slvalair(i,j,2)/slvalone(i,j,2) + slvalair(i,j,3)/slvalone(i,j,3) )    
!           endif
        end do
      end do 
      
      do jcell=1,nc
        do icell=1,nc          
            i=1+(icell-1)*nipm
            j=1+(jcell-1)*nipm
            
            sga=spelt(ie)%sga(i+1,j+1)   
            dx=spelt(ie)%dab(icell)   
            dy=spelt(ie)%dab(jcell)
            
            flux(1) = dy * (fluxval(i,j,1) + 4.0D0 * fluxval(i,j+1,1) + fluxval(i,j+2,1))/6.0D0  ! west
            flux(2) = dx * (fluxval(i,j,2) + 4.0D0 * fluxval(i+1,j,2) + fluxval(i+2,j,2))/6.0D0  ! south
            flux(3) = dy * (fluxval(i+2,j,1) + 4.0D0 * fluxval(i+2,j+1,1) + fluxval(i+2,j+2,1))/6.0D0 ! east
            flux(4) = dx * (fluxval(i+2,j+2,2) + 4.0D0 * fluxval(i+1,j+2,2) + fluxval(i,j+2,2))/6.0D0 ! north
            
            spelt(ie)%c(i+1,j+1,k,1,tl%np1) = spelt(ie)%c(i+1,j+1,k,1,tl%n0)/sga + &
                                      (flux(1) + flux(2) - flux(3) - flux(4) ) / (spelt(ie)%area_sphere(icell,jcell))
            spelt(ie)%c(i+1,j+1,k,1,tl%np1)=spelt(ie)%c(i+1,j+1,k,1,tl%np1)*sga                                              
        end do
      end do
      
      ! Tracers
      do itr=2,ntrac
        do j=1-nhe,nc+nhe
          do i=1-nhe,nc+nhe
            icell=1+(i-1)*nipm
            jcell=1+(j-1)*nipm
            ff=spelt(ie)%c(icell:icell+nipm,jcell:jcell+nipm,k,itr,tl%n0)
            minmax(i,j,:)=cell_minmax(ff)
            call cip_coeff(spelt(ie)%drefx(i,j),spelt(ie)%drefy(i,j),ff,ff(2,2),cf(:,:,i,j))
          enddo
        enddo
        do j=1,nep
          do i=1,nep  
            sga=spelt(ie)%sga(i,j)
            slval(1)=spelt(ie)%c(i,j,k,itr,tl%n0)
 
            tmp=cip_interpolate(cf(:,:,icell1(i,j),jcell1(i,j)),dref1(i,j)%x,dref1(i,j)%y) 
            tmp=qmsl_cell_filter(icell1(i,j),jcell1(i,j),minmax,tmp)
            slval(2)=(sga/sg1(i,j))*tmp

            tmp=cip_interpolate(cf(:,:,icell2(i,j),jcell2(i,j)),dref2(i,j)%x,dref2(i,j)%y) 
            tmp=qmsl_cell_filter(icell2(i,j),jcell2(i,j),minmax,tmp)
            slval(3)=(sga/sg2(i,j))*tmp
            spelt(ie)%c(i,j,k,itr,tl%np1)=slval(3)
            
!            if (mod(i,2)==1) then                   ! works only for nip=3!!!
!               fluxval(i,j,1) =  dt6 * (contrauv(i,j,1,1)*slvalair(i,j,1)*slval(1)/slvalone(i,j,1) + & 
!                      4.0D0 * contrauv(i,j,1,2)*slvalair(i,j,2)*slval(2)/slvalone(i,j,2) + &
!                      contrauv(i,j,1,3)*slvalair(i,j,3)*slval(3)/slvalone(i,j,3) )
!             endif
!             if (mod(j,2)==1) then            ! works only for nip=3!!!
!               fluxval(i,j,2) =  dt6 * (contrauv(i,j,2,1)*slvalair(i,j,1)*slval(1)/slvalone(i,j,1) + & 
!                      4.0D0 * contrauv(i,j,2,2)*slvalair(i,j,2)*slval(2)/slvalone(i,j,2) + &
!                      contrauv(i,j,2,3)*slvalair(i,j,3)*slval(3)/slvalone(i,j,3) )    
!             endif
            if (mod(i,2)==1) then                   ! works only for nip=3!!!
               fluxval(i,j,1) =  dt6 * (contrauv(i,j,1,1)*slvalair(i,j,1)*slval(1)/sga + & 
                      4.0D0 * contrauv(i,j,1,2)*slvalair(i,j,2)*slval(2)/sga + &
                      contrauv(i,j,1,3)*slvalair(i,j,3)*slval(3)/sga )
             endif
             if (mod(j,2)==1) then            ! works only for nip=3!!!
               fluxval(i,j,2) =  dt6 * (contrauv(i,j,2,1)*slvalair(i,j,1)*slval(1)/sga + & 
                      4.0D0 * contrauv(i,j,2,2)*slvalair(i,j,2)*slval(2)/sga + &
                      contrauv(i,j,2,3)*slvalair(i,j,3)*slval(3)/sga )    
             endif
          end do
        end do 
        
        do jcell=1,nc
          do icell=1,nc          
              i=1+(icell-1)*nipm
              j=1+(jcell-1)*nipm
              
              sga=spelt(ie)%sga(i+1,j+1)   
              dx=spelt(ie)%dab(icell)   
              dy=spelt(ie)%dab(jcell)
              
              flux(1) = dy * (fluxval(i,j,1) + 4.0D0 * fluxval(i,j+1,1) + fluxval(i,j+2,1))/6.0D0  ! west
              flux(2) = dx * (fluxval(i,j,2) + 4.0D0 * fluxval(i+1,j,2) + fluxval(i+2,j,2))/6.0D0  ! south
              flux(3) = dy * (fluxval(i+2,j,1) + 4.0D0 * fluxval(i+2,j+1,1) + fluxval(i+2,j+2,1))/6.0D0 ! east
              flux(4) = dx * (fluxval(i+2,j+2,2) + 4.0D0 * fluxval(i+1,j+2,2) + fluxval(i,j+2,2))/6.0D0 ! north
              
              spelt(ie)%c(i+1,j+1,k,itr,tl%np1) = (spelt(ie)%c(i+1,j+1,k,itr,tl%n0)/sga)*(spelt(ie)%c(i+1,j+1,k,1,tl%n0)/sga) + &
                                        (flux(1) + flux(2) - flux(3) - flux(4) ) / (spelt(ie)%area_sphere(icell,jcell))
              spelt(ie)%c(i+1,j+1,k,itr,tl%np1)=spelt(ie)%c(i+1,j+1,k,itr,tl%np1)*sga/(spelt(ie)%c(i+1,j+1,k,1,tl%np1)/sga)
          end do
        end do
      end do
    end do
    call ghostVpack2d(cellghostbuf,spelt(ie)%c,nipm, nep,nlev,ntrac,0, tl%np1, timelevels,elem(ie)%desc)
  end do
!-----------------------------------------------------------------------------------! 
  call t_startf('SPELT Communication') 
  call ghost_exchangeV(hybrid,cellghostbuf,nipm,nep,ntrac)
  call t_stopf('SPELT Communication')
!-----------------------------------------------------------------------------------!  
  call t_startf('SPELT Unpacking')  
  do ie=nets,nete
    call ghostVunpack2d(cellghostbuf,spelt(ie)%c,nipm, nep,nlev,ntrac,0, tl%np1, timelevels,elem(ie)%desc)
  end do
  call t_stopf('SPELT Unpacking')
  call t_stopf('SPELT scheme') 
end subroutine spelt_runair


subroutine spelt_runlimit(elem,spelt,hybrid,deriv,tstep,tl,nets,nete)
! 
  use derivative_mod, only : derivative_t
  ! ---------------------------------------------------------------------------------
  use edge_mod, only :  ghostVpack2d, ghostVunpack2d
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
  ! ---------------------------------------------------------------------------------
  use coordinate_systems_mod, only : spherical_to_cart, cart2cubedspherexy, ref2sphere, sphere2cubedsphere
  ! ------EXTERNAL----------------
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------  
  use bndry_mod, only : ghost_exchangevfull
  use edge_mod, only : ghostbuffertr_t, ghostvpack, ghostvunpack, initghostbuffer,&
       freeghostbuffer
       
  implicit none
  type (element_t), intent(inout)             :: elem(:)
  type (spelt_struct), intent(inout)          :: spelt(:)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  type (derivative_t), intent(in)             :: deriv           ! derivative struct
  real (kind=real_kind), intent(in)           :: tstep
  type (TimeLevel_t), intent(in)              :: tl              ! time level struct
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
 
!   integer                                     :: i,j,k,ie,itr
!         
!   real (kind=real_kind)                       :: ff(nip,nip), c_low(1:nc,1:nc)
!   real (kind=real_kind)                       :: cf(nip,nip,1-nhe:nc+nhe,1-nhe:nc+nhe)
!   type (spherical_polar_t)                    :: dsphere1(1:nep,1:nep), dsphere2(1:nep,1:nep)
!     
!   real (kind=real_kind)                       :: slval(3), fluxval(nep,nep,2), flux(4)
!   real (kind=real_kind)                       :: fluxlowx(nc+1,nc), fluxlowy(nc,nc+1)
!   real (kind=real_kind)                       :: fluxhigh(nets:nete,nc,nc,nlev,ntrac,4)
!   
!   type (cartesian2D_t)                        :: dref1(1:nep,1:nep), dref2(1:nep,1:nep)
!   real (kind=real_kind)                       :: sg1(1:nep,1:nep),sg2(1:nep,1:nep)  
!   real (kind=real_kind)                       :: minmax(1-nhe:nc+nhe,1-nhe:nc+nhe,2)
!   integer                                     :: icell1(1:nep,1:nep), jcell1(1:nep,1:nep)     
!   integer                                     :: icell2(1:nep,1:nep), jcell2(1:nep,1:nep)
!   
!   integer                                     :: icell, jcell, jx, jy
!   real (kind=real_kind)                       :: dxoy, dxyi, dt6, sg, sga
!   type (cartesian2D_t)                        :: alphabeta
!   real (kind=real_kind)                       :: dx, dy, tmp
!   
!   real (kind=real_kind)                       :: R(nets:nete,1-nhe:nc+nhe,1-nhe:nc+nhe,nlev,ntrac,1), c_min, c_max, cmo(4)
! 
!   integer                                     :: face_nodep
!   
!   call t_startf('SPELT scheme') 
!     
!   dt6  = tstep/ 6.0D0
!   do ie=nets,nete 
!     do k=1, nlev
!       call solidbody_all(spelt(ie), dsphere1,dsphere2,k)                                                                   
! !       call boomerang_all(spelt(ie), dsphere1,dsphere2,k,tl%nstep)
! !       call spelt_dep_from_gll(elem(ie), deriv, spelt(ie)%asphere,dsphere1,0.5D0*tstep,tl,k)         
! !       call spelt_dep_from_gll(elem(ie), deriv, spelt(ie)%asphere,dsphere2,tstep,tl,k)         
!       !search has not to be done for all tracers!
!       do j=1,nep
!         do i=1,nep
! !           call solidbody(spelt(ie)%asphere(i,j), dsphere, 0.5D0) 
!           call cell_search(elem(ie),spelt(ie), dsphere1(i,j),icell1(i,j), jcell1(i,j),dref1(i,j),alphabeta,face_nodep)
!           sg1(i,j)=metric_term(alphabeta)
! !           sg1(i,j)=metric_termref(elem(ie),dref1(i,j))
! !           call solidbody(spelt(ie)%asphere(i,j), dsphere, 1.0D0)  
!           call cell_search(elem(ie),spelt(ie), dsphere2(i,j), icell2(i,j), jcell2(i,j),dref2(i,j),alphabeta, face_nodep)
!           sg2(i,j)=metric_term(alphabeta)
! !           sg2(i,j)=metric_termref(elem(ie),dref2(i,j))
!         end do
!       end do
!       ! search of both point on the trajectory done
!       do itr=1,ntrac
!         do j=1-nhe,nc+nhe
!           do i=1-nhe,nc+nhe
!             icell=1+(i-1)*nipm
!             jcell=1+(j-1)*nipm
!             ff=spelt(ie)%c(icell:icell+nipm,jcell:jcell+nipm,k,itr,tl%n0)
!             minmax(i,j,:)=cell_minmax(ff)
!             call cip_coeff(spelt(ie)%drefx(i,j),spelt(ie)%drefy(i,j),ff,ff(2,2),cf(:,:,i,j))
!           enddo
!         enddo
!         !reconstruction coefficients caculated done
!         do j=1,nep
!           do i=1,nep  
!             sga=spelt(ie)%sga(i,j)
!             slval(1)=spelt(ie)%c(i,j,k,itr,tl%n0)
!  
!             tmp=cip_interpolate(cf(:,:,icell1(i,j),jcell1(i,j)),dref1(i,j)%x,dref1(i,j)%y) 
!             tmp=qmsl_cell_filter(icell1(i,j),jcell1(i,j),minmax,tmp)
!             slval(2)=(sga/sg1(i,j))*tmp
! 
!             tmp=cip_interpolate(cf(:,:,icell2(i,j),jcell2(i,j)),dref2(i,j)%x,dref2(i,j)%y) 
!             tmp=qmsl_cell_filter(icell2(i,j),jcell2(i,j),minmax,tmp)
!             slval(3)=(sga/sg2(i,j))*tmp
!             spelt(ie)%c(i,j,k,itr,tl%np1)=slval(3)
! 
!            if (mod(i,2)==1) then                   ! works only for nip=3!!!
!               fluxval(i,j,1) =  dt6 * spelt(ie)%contrau(i,j,k)* (slval(1) + & 
!                      4.0D0 * slval(2) + slval(3) )
!             endif
!             if (mod(j,2)==1) then            ! works only for nip=3!!!
!               fluxval(i,j,2) =  dt6 * spelt(ie)%contrav(i,j,k)* (slval(1) + & 
!                      4.0D0 * slval(2) + slval(3) )    
!             endif  
! !             if (mod(i,2)==1) then                   ! works only for nip=3!!!                                                
! !                fluxval(i,j,1) =  dt6 * (spelt(ie)%contrau(i,j,k)* slval(1) + &                                               
! !                       4.0D0 * spelt(ie)%contrau1(i,j)*slval(2) + spelt(ie)%contrau2(i,j)*slval(3) )                          
! !                                                                                                                              
! !              endif                                                                                                           
! !              if (mod(j,2)==1) then            ! works only for nip=3!!!                                                      
! !                fluxval(i,j,2) =  dt6 * (spelt(ie)%contrav(i,j,k)* slval(1) + &                                               
! !                         4.0D0 * spelt(ie)%contrav1(i,j)*slval(2) + spelt(ie)%contrav2(i,j)*slval(3) )                       \
! !                                                                                                                              
! !              endif  
!           end do
!         end do 
!         !calculate low order flux
!         dxoy=abs(elem(ie)%corners(1)%x-elem(ie)%corners(2)%x)/nc
!         
!         do j=1,nc+1
!           do i=1,nc
!             icell=2+(i-1)*nipm
!             jcell=1+(j-1)*nipm
!             
!             if (spelt(ie)%contrau(jcell,icell,k)>0.0D0) then
!               fluxlowx(j,i)=dxoy*tstep*spelt(ie)%contrau(jcell,icell,k)* &
!                             spelt(ie)%c(jcell-1,icell,k,itr,tl%n0)!*spelt(ie)%sga(jcell,icell)/spelt(ie)%sga(jcell-1,icell)
!             else
!               fluxlowx(j,i)=dxoy*tstep*spelt(ie)%contrau(jcell,icell,k)* &
!                             spelt(ie)%c(jcell+1,icell,k,itr,tl%n0)!*spelt(ie)%sga(jcell,icell)/spelt(ie)%sga(jcell+1,icell)
!             endif
!             if (spelt(ie)%contrav(icell,jcell,k)>0.0D0) then
!               fluxlowy(i,j)=dxoy*tstep*spelt(ie)%contrav(icell,jcell,k)* &
!                             spelt(ie)%c(icell,jcell-1,k,itr,tl%n0)!*spelt(ie)%sga(icell,jcell)/spelt(ie)%sga(icell,jcell-1)
!             else
!               fluxlowy(i,j)=dxoy*tstep*spelt(ie)%contrav(icell,jcell,k)* &
!                             spelt(ie)%c(icell,jcell+1,k,itr,tl%n0)!*spelt(ie)%sga(icell,jcell)/spelt(ie)%sga(icell,jcell+1)
!             endif    
!             
!           end do
!         end do
! !!!! only needed for filter              
!         do j=1,nc
!           jcell=2+(j-1)*nipm
!           jy=1+(j-1)*nipm
!           do i=1,nc
!             dx=spelt(ie)%dab(i)   
!             dy=spelt(ie)%dab(j)
!             icell=2+(i-1)*nipm
!         ! for low order flux
!             c_low(i,j) = spelt(ie)%c(icell,jcell,k,itr,tl%n0) + &
!                                   (fluxlowx(i,j) + fluxlowy(i,j) - fluxlowx(i+1,j) - fluxlowy(i,j+1) ) / (dx*dy)
!                                   
!             spelt(ie)%c(icell,jcell,k,itr,tl%np1)=c_low(i,j) 
!             !high order flux
!             jx=1+(i-1)*nipm
!             fluxhigh(ie,i,j,k,itr,1) = &
!                  dx * (fluxval(jx,jy,1) + 4.0D0 * fluxval(jx,jy+1,1) + fluxval(jx,jy+2,1)) / 6.0D0  ! west
!             fluxhigh(ie,i,j,k,itr,2) = &
!                  dy * (fluxval(jx,jy,2) + 4.0D0 * fluxval(jx+1,jy,2) + fluxval(jx+2,jy,2)) / 6.0D0  ! south
!             fluxhigh(ie,i,j,k,itr,3) = &
!                  dx * (fluxval(jx+2,jy,1) + 4.0D0 * fluxval(jx+2,jy+1,1) + fluxval(jx+2,jy+2,1)) / 6.0D0 ! east
!             fluxhigh(ie,i,j,k,itr,4) = &
!                  dy * (fluxval(jx+2,jy+2,2) + 4.0D0 * fluxval(jx+1,jy+2,2) + fluxval(jx,jy+2,2)) / 6.0D0 ! north
!             
!             R(ie,i,j,k,itr,1)=-min(0.0D0,fluxhigh(ie,i,j,k,itr,1))-min(0.0D0,fluxhigh(ie,i,j,k,itr,2)) + &
!                                max(0.0D0,fluxhigh(ie,i,j,k,itr,3))+max(0.0D0,fluxhigh(ie,i,j,k,itr,4))
! !             
!             if (R(ie,i,j,k,itr,1)>0.0D0) then
!               R(ie,i,j,k,itr,1)=min(1.0D0,spelt(ie)%c(icell,jcell,k,itr,tl%n0)*dx*dy/R(ie,i,j,k,itr,1))
!             else
!               R(ie,i,j,k,itr,1)=0.0D0  
!             endif
!             
!             !overwrite low order flux by antidiffusive flux
! !             fluxlowx(i,j)=flux(1)-fluxlowx(i,j) 
! !             fluxlowy(i,j)=flux(2)-fluxlowx(i,j)         
!           end do
!         end do
! !         do i=1,nc
! !           !high order flux
! !           jx=1+(i-1)*nipm
! !           jy=nep
! !           flux(2) = dxoy * (fluxval(jx,jy,2) + 4.0D0 * fluxval(jx+1,jy,2) + fluxval(jx+2,jy,2))  
! !           fluxlowy(i,nc+1)=flux(2)-fluxlowy(i,nc+1)      
! !         end do
!       end do !ntrac
!     end do   !nlev
!     call ghostVpack(factorR, R(ie,:,:,:,:,:),nhe,nc,nlev,ntrac,0,1,1,elem(ie)%desc)
!   end do
!        ! Anti diffusive flux are computed for each cell, done!
!   call ghost_exchangeV(hybrid,factorR,nhe,nc,ntrac)
! 
!   do ie=nets,nete
!     call ghostVunpack(factorR, R(ie,:,:,:,:,:), nhe, nc,nlev,ntrac,0, 1,1,elem(ie)%desc)
!     do k=1, nlev
!       do itr=1,ntrac
!         !!!! only needed for filter     
!         do j=1,nc !nep,2
!           do i=1,nc !nep,2    
!               dx=spelt(ie)%dab(i)   
!               dy=spelt(ie)%dab(j)        
!               icell=2+(i-1)*nipm
!               jcell=2+(j-1)*nipm   
!               cmo=1.0D0 
!                 !west
!                 if(fluxhigh(ie,i,j,k,itr,1)>=0.0D0) then
!                   cmo(1)=R(ie,i-1,j,k,itr,1)  
!                 else
!                   cmo(1)=R(ie,i,j,k,itr,1)  
!                 endif 
!                 !south                         
!                 if(fluxhigh(ie,i,j,k,itr,2)>=0.0D0) then
!                   cmo(2)=R(ie,i,j-1,k,itr,1)  
!                 else
!                   cmo(2)=R(ie,i,j,k,itr,1)  
!                 endif  
!                 !east  
!                 if(fluxhigh(ie,i,j,k,itr,3)>=0.0D0) then
!                   cmo(3)=R(ie,i,j,k,itr,1)  
!                 else
!                   cmo(3)=R(ie,i+1,j,k,itr,1)  
!                 endif 
!                 !north      
!                 if(fluxhigh(ie,i,j,k,itr,4)>=0.0D0) then
!                   cmo(4)=R(ie,i,j,k,itr,1)  
!                 else
!                   cmo(4)=R(ie,i,j+1,k,itr,1)  
!                 endif   
!                 spelt(ie)%c(icell,jcell,k,itr,tl%np1) = spelt(ie)%c(icell,jcell,k,itr,tl%n0) - &
!                      (-cmo(1)*fluxhigh(ie,i,j,k,itr,1) - cmo(2)*fluxhigh(ie,i,j,k,itr,2) + &
!                      cmo(3)*fluxhigh(ie,i,j,k,itr,3) + cmo(4)*fluxhigh(ie,i,j,k,itr,4) ) /(dx*dy)
!           end do
!         end do                
!       end do !ntrac
!     end do   !nlev
!     call ghostVpack2d(cellghostbuf,spelt(ie)%c,nipm, nep,nlev,ntrac,0, tl%np1, timelevels,elem(ie)%desc)
!   end do
!   call t_stopf('SPELT scheme')
! !-----------------------------------------------------------------------------------! 
!   call t_startf('SPELT Communication') 
!   call ghost_exchangeV(hybrid,cellghostbuf,nipm,nep,ntrac)
!   call t_stopf('SPELT Communication')
! !-----------------------------------------------------------------------------------!  
!   call t_startf('SPELT Unpacking')  
!   do ie=nets,nete
!     call ghostVunpack2d(cellghostbuf,spelt(ie)%c,nipm, nep,nlev,ntrac,0, tl%np1, timelevels,elem(ie)%desc)
!   end do
!   call t_stopf('SPELT Unpacking')
end subroutine spelt_runlimit



subroutine get_contravelocities(elem, spelt, contrauv,k,deriv)
  use physical_constants, only : rearth
  use derivative_mod, only : derivative_t, interpolate_gll2spelt_points

  implicit none
  type (element_t)                         :: elem
  type (spelt_struct), intent(in)          :: spelt
  real (kind=real_kind), intent(out)       :: contrauv(nep,nep,2,3)
  integer, intent(in)                      :: k
  type (derivative_t), intent(in)          :: deriv           ! derivative struct
  
  integer                                  :: i,j
  real (kind=real_kind)                    :: vstar(np,np,2), v1, v2 
  
  vstar=elem%derived%vstar(:,:,:,k)/rearth
  do j=1,np
    do i=1,np
      v1 = spelt%Dinv(1,1,i,j)*vstar(i,j,1) + spelt%Dinv(1,2,i,j)*vstar(i,j,2)
      v2 = spelt%Dinv(2,1,i,j)*vstar(i,j,1) + spelt%Dinv(2,2,i,j)*vstar(i,j,2)
      vstar(i,j,1)=v1
      vstar(i,j,2)=v2    
    enddo
  enddo
  contrauv(:,:,1,1)=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
  contrauv(:,:,2,1)=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
  
  vstar=spelt%vn12(:,:,:,k)/rearth
  do j=1,np
    do i=1,np
      v1 = spelt%Dinv(1,1,i,j)*vstar(i,j,1) + spelt%Dinv(1,2,i,j)*vstar(i,j,2)
      v2 = spelt%Dinv(2,1,i,j)*vstar(i,j,1) + spelt%Dinv(2,2,i,j)*vstar(i,j,2)
      vstar(i,j,1)=v1
      vstar(i,j,2)=v2    
    enddo
  enddo
  contrauv(:,:,1,2)=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
  contrauv(:,:,2,2)=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
  
  vstar=spelt%vn0(:,:,:,k)/rearth
  do j=1,np
    do i=1,np
      v1 = spelt%Dinv(1,1,i,j)*vstar(i,j,1) + spelt%Dinv(1,2,i,j)*vstar(i,j,2)
      v2 = spelt%Dinv(2,1,i,j)*vstar(i,j,1) + spelt%Dinv(2,2,i,j)*vstar(i,j,2)
      vstar(i,j,1)=v1
      vstar(i,j,2)=v2    
    enddo
  enddo
  contrauv(:,:,1,3)=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
  contrauv(:,:,2,3)=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
  
end subroutine get_contravelocities

subroutine boomerang_all(spelt, dsphere1,dsphere2,contrauv,k,nstep)
  use kinds, only : real_kind
  use time_mod, only : tstep, nmax, ndays,Time_at
  use physical_constants, only : DD_PI, rearth
  use coordinate_systems_mod, only : spherical_polar_t
  use coordinate_systems_mod, only : sphere2cubedsphere

  implicit none
  type (spelt_struct), intent(inout)          :: spelt
  type (spherical_polar_t),intent(out)        :: dsphere1(nep,nep),dsphere2(nep,nep)
  real (kind=real_kind), intent(out)          :: contrauv(nep,nep,2,3)
  
  integer, intent(in)                         :: k
  integer, intent(in)                         :: nstep
  
  integer                     :: iorder, i,j,itr, iteration
  real (kind=real_kind)       :: slat, clat, slon, clon, xg, yg, zg, ca, sa, co, so, &
                                 slon2,clon2,slonh,clonh
  real (kind=real_kind)       :: tmp_time,tstart,lon,lat,tmp_dt,tmp_lm,tmp_th, ck, tt, omega
  real (kind=real_kind)       :: dplm,dpth,trm1,trm2,sslm,cwt,swt
  real (kind=real_kind)       :: dt2,dt3,dt4,dtt,udc,uexact, vexact, u,v
  
  real (kind=real_kind)       :: Ainv(2,2), time
  
  type (cartesian2D_t)        :: alphabeta
  type (cartesian2D_t)        :: dcart
  integer                     :: number, face_nodep
  
  integer                     :: nmaxaround
  
  iteration=10   ! has to be an even number here
  iorder=5
  
  nmaxaround=12*nmax/ndays          ! in 12 days around the earth
!   tmp_time=(nstep)*5.0D0/nmaxaround
!   tmp_dt = 5.0D0/nmaxaround/iteration
  tstart=(nstep+1)*5.0D0/(nmax)
  tmp_dt = 5.0D0/(nmax)/iteration
  dt2 = tmp_dt*tmp_dt/2.0D0
  dt3 = dt2*tmp_dt/3.0D0
  dt4 = dt3*tmp_dt/4.0D0
  tt = 5.0D0  !total time
  ck = 10.0D0/tt
  omega = DD_PI/tt
  
  do j=1,nep
    do i=1,nep 
      tmp_time=tstart
      tmp_lm = spelt%asphere(i,j)%lon
      tmp_th = spelt%asphere(i,j)%lat
      do itr=1,iteration
        if (itr>1) then
          tmp_time=tmp_time-tmp_dt
        end if 
     
        cwt = cos(omega*tmp_time)
        swt = sin(omega*tmp_time)

        lon   = tmp_lm
        lat   = tmp_th
        slat  = sin(lat)
        clat  = cos(lat)
    
        slon  = sin(lon-tmp_time*2*DD_PI/tt)!solid-body rotation added
        slon2 = sin(2*(lon-tmp_time*2*DD_PI/tt))!solid-body rotation added
        clon  = cos(lon-tmp_time*2*DD_PI/tt)!solid-body rotation added
        clon2 = cos(2*(lon-tmp_time*2*DD_PI/tt))!solid-body rotation added
        sslm  = (sin(0.5*(lon-tmp_time*2*DD_PI/tt)))**2
        
        uexact =  ck*slon*slon*sin(2*lat)*cos(tmp_time*omega) + clat*2*DD_PI/(tt)
        vexact =  ck*slon2*clat*cos(tmp_time*omega)
        udc    =  2*ck*slon*slon*slat*cos(tmp_time*omega) + 2*DD_PI/tt
        ! 2nd-order
        tmp_lm = lon - tmp_dt * udc 
    
        ! 3rd-order
        if (iorder>2) then
             tmp_lm = tmp_lm - 2 * dt2 * slon * (slon * swt * omega * slat - 2 * udc * cwt * clon * slat&
             - vexact * slon * cwt * clat) * ck 
        endif
        ! 4th-order
        if (iorder>3) then
             tmp_lm = tmp_lm - dt3 * (-4 * slon * swt * (2 * udc * clon * slat + vexact * slon * clat) * &
             ck * omega + 2 * cwt * (-omega ** 2 * slat + omega ** 2 * slat * clon ** 2 &
             + 4 * udc ** 2 * clon ** 2 * slat - 2 * udc ** 2 * slat + 4 * udc * vexact * &
             slon * clat * clon - vexact ** 2 * slat + vexact ** 2 * slat * clon ** 2) * ck)
        endif
        !
        ! 5th-order
        !
        if (iorder>4) then
             tmp_lm = tmp_lm + dt4 * (-2 * swt * (-slat * omega ** 2 + omega ** 2 * slat * clon ** 2 + &
             12 * udc ** 2 * clon ** 2 * slat - 6 * udc ** 2 * slat + 12 * udc * vexact * &
             slon * clon * clat - 3 * vexact ** 2 * slat + 3 * vexact ** 2 * slat * clon ** 2) &
             * ck * omega - 2 * cwt * (vexact ** 3 * clat - 3 * vexact * omega ** 2 * clat * &
             clon ** 2 - vexact ** 3 * clat * clon ** 2 - 8 * udc ** 2 * vexact * clon ** 2 * clat &
             + 6 * udc * vexact ** 2 * slat * slon * clon + 3 * vexact * omega ** 2 * clat + 6 * &
             udc ** 2 * vexact * clat + 8 * udc ** 3 * slon * clon * slat + 6 * udc * slon * &
             omega ** 2 * clon * slat) * ck)
        endif
        ! 2nd-order
        tmp_th = lat - tmp_dt * vexact 
        ! 3rd-order
        if (iorder>2) then
             tmp_th = tmp_th &
    !- 2 * dt2 * ck * (clat * swt * omega * slon * clon - 2 * udc * clat * cwt * &
    !             clon ** 2 + udc * clat * cwt + vexact * slat * cwt * slon * clon)
             + dt2 * (-ck * slon2 * clat * swt * omega + 2 * udc * ck * clon2 *&
             clat * cwt - vexact * ck * slon2 * slat * cwt) 
        endif
        ! 4th-order
        if (iorder>3) then    
            tmp_th = tmp_th - dt3 * (-4 * omega * ck * swt * (2 * udc * clat * clon ** 2 &
             - udc * clat - vexact * slat * slon * clon) - 2 * ck * cwt * (clat * omega ** 2 * &
             slon * clon + 4 * udc ** 2 * slon * clat * clon + 12 * udc * vexact * slat &
             * clon ** 2 - 2 * udc * vexact * slat + vexact ** 2 * slon * clat * clon))
        endif
        ! 5th-order
        if (iorder>4) then
          tmp_th = tmp_th + dt4 * (2 * omega * ck * swt * (clat * omega ** 2 * slon * clon - 6 * &
          udc * vexact * slat + 12 * udc ** 2 * slon * clon * clat + 12 * udc * vexact * &
          slat * clon ** 2 + 3 * vexact ** 2 * slon * clon * clat) - 2 * ck * cwt * &
          (-vexact ** 3 * slat * slon * clon - 3 * udc * clat * omega ** 2 + 6 * udc * clat &
          * omega ** 2 * clon ** 2 - 3 * vexact * slat * omega ** 2 * slon * clon + 6 * &
          udc ** 3 * clat * clon ** 2 - 12 * udc ** 2 * vexact * slat * slon * clon + 6 * &
          udc * vexact ** 2 * clon ** 2 * clat - 4 * udc * vexact ** 2 * clat - 4 * &
          udc ** 3 * clat))
        endif
        ! need value half of the way
        if(itr==5) then
          dsphere1(i,j)%lon=tmp_lm
          dsphere1(i,j)%lat=tmp_th
          dsphere1(i,j)%r=spelt%asphere(i,j)%r
          u=uexact /( 12.0D0*3600.0D0*24.0D0/5.0D0)
          v=vexact /( 12.0D0*3600.0D0*24.0D0/5.0D0)
          time=((2*nstep+1)*5.0D0/nmax)/2.0D0
          call get_boomerang_velocities(spelt%asphere(i,j), u,v, time)
!           call cube_facepoint_ne(dsphere1,ne,dcart, alphabeta, number, face_nodep)
!           call cube_facepoint_ne(dsphere1(i,j),ne,dcart, alphabeta, number, face_nodep)
!           call get_Ainv(spelt%Facenumber, dsphere1(i,j), alphabeta,Ainv)
!           spelt%contrau1(i,j)=Ainv(1,1)*u+Ainv(2,1)*v      
!           spelt%contrav1(i,j)=Ainv(1,2)*u+Ainv(2,2)*v
          contrauv(i,j,1,2)=spelt%Ainv(1,1,i,j)*u+spelt%Ainv(2,1,i,j)*v      
          contrauv(i,j,2,2)=spelt%Ainv(1,2,i,j)*u+spelt%Ainv(2,2,i,j)*v
        endif
        if(itr==1) then
          u=uexact /( 12.0D0*3600.0D0*24.0D0/5.0D0)
          v=vexact /( 12.0D0*3600.0D0*24.0D0/5.0D0)
          contrauv(i,j,1,1)=spelt%Ainv(1,1,i,j)*u+spelt%Ainv(2,1,i,j)*v      
          contrauv(i,j,2,1)=spelt%Ainv(1,2,i,j)*u+spelt%Ainv(2,2,i,j)*v
!           alphabeta=sphere2cubedsphere(spelt%asphere(i,j), spelt%Facenumber)
!           face_nodep=spelt%Facenumber
!           call cube_facepoint_ne(spelt%asphere(i,j),ne,dcart, alphabeta, number, face_nodep)
!           if(face_nodep .ne. spelt%Facenumber) then
! !              write(*,*) 'face', face_nodep, spelt%Facenumber
!              alphabeta=sphere2cubedsphere(spelt%asphere(i,j), face_nodep)
!              write(*,*) 'alpha1',metric_term(alphabeta)
!              alphabeta=sphere2cubedsphere(spelt%asphere(i,j), spelt%Facenumber)
!              write(*,*) 'alpha2',metric_term(alphabeta)
!              
!              call get_Ainv(spelt%Facenumber, spelt%asphere(i,j), alphabeta,Ainv)
! !              write(*,*) 'vor', spelt%contrau(i,j,k),spelt%contrav(i,j,k)
! !              spelt%contrau(i,j,k)=Ainv(1,1)*u+Ainv(2,1)*v      
! !              spelt%contrav(i,j,k)=Ainv(1,2)*u+Ainv(2,2)*v
! !              write(*,*) 'nach', spelt%contrau(i,j,k),spelt%contrav(i,j,k)
! !              
!           endif
!           call get_Ainv(face_nodep, spelt%asphere(i,j), alphabeta,Ainv)
!           spelt%contrau(i,j,k)=Ainv(1,1)*u+Ainv(2,1)*v      
!           spelt%contrav(i,j,k)=Ainv(1,2)*u+Ainv(2,2)*v
        endif
      end do
      dsphere2(i,j)%lon=tmp_lm
      dsphere2(i,j)%lat=tmp_th
      dsphere2(i,j)%r=spelt%asphere(i,j)%r

      u=uexact /( 12.0D0*3600.0D0*24.0D0/5.0D0)
      v=vexact /( 12.0D0*3600.0D0*24.0D0/5.0D0)
!       call cube_facepoint_ne(dsphere2(i,j),ne,dcart, alphabeta, number, face_nodep)
!       call get_Ainv(spelt%Facenumber, dsphere2(i,j), alphabeta,Ainv)
!       spelt%contrau2(i,j)=Ainv(1,1)*u+Ainv(2,1)*v      
!       spelt%contrav2(i,j)=Ainv(1,2)*u+Ainv(2,2)*v
!       time=((nstep)*5.0D0/nmax)
      call get_boomerang_velocities(spelt%asphere(i,j), u,v, ((nstep)*5.0D0/nmax))

      contrauv(i,j,1,3)=spelt%Ainv(1,1,i,j)*u+spelt%Ainv(2,1,i,j)*v      
      contrauv(i,j,2,3)=spelt%Ainv(1,2,i,j)*u+spelt%Ainv(2,2,i,j)*v
    enddo
  enddo
end subroutine boomerang_all


subroutine get_boomerang_velocities(sphere, u,v, time)
  use coordinate_systems_mod, only : cart2cubedspherexy, spherical_to_cart
  use element_mod, only : element_t
  use fvm_control_volume_mod, only: fvm_struct
  use physical_constants, only : DD_PI, rearth
  use time_mod, only : ndays
  
  implicit none
  type (spherical_polar_t), intent(in)   :: sphere
  real (kind=real_kind),intent(in)             :: time  ! time of the arrival grid
  real (kind=real_kind),intent(out)             :: u,v
  
  integer                              :: i,j
  real (kind=real_kind)       :: slat, clat, slon, clon, &
                                 slon2,clon2, tt, ck, omega, lon, lat, tmp_time
  tt = 5.0D0  !total time
  ck = 10.0D0/tt
  omega = DD_PI/tt

  !tmp_time=(nstep+1)*5.0D0/nmax
!   tmp_time=5.0D0*time/(ndays*3600*24.0D0)  ! convert from dimensional time to dimensionless
  tmp_time=time
      lon = sphere%lon
      lat = sphere%lat

      slat  = sin(lat)
      clat  = cos(lat)

      slon  = sin(lon-tmp_time*2*DD_PI/tt)!solid-body rotation added
      slon2 = sin(2*(lon-tmp_time*2*DD_PI/tt))!solid-body rotation added
      clon  = cos(lon-tmp_time*2*DD_PI/tt)!solid-body rotation added
      clon2 = cos(2*(lon-tmp_time*2*DD_PI/tt))!solid-body rotation added

      u =  ck*slon*slon*sin(2*lat)*cos(tmp_time*omega) + clat*2*DD_PI/(tt)
      v =  ck*slon2*clat*cos(tmp_time*omega)
      ! convert from radians per dimensionless time to 
      ! meters/sec 
      u=u  /( ndays*3600*24.0D0/5.0D0)
      v=v  /( ndays*3600*24.0D0/5.0D0)

end subroutine get_boomerang_velocities

subroutine get_Ainv(FaceNum, sphere, alphabeta,Ainv)
  use physical_constants, only : DD_PI
  use coordinate_systems_mod, only : spherical_polar_t
  
  integer              , intent(in)   :: FaceNum
  type (spherical_polar_t),intent(in) :: sphere
  type (cartesian2D_t), intent(in)        :: alphabeta
  
  real (kind=real_kind), intent(out)  :: Ainv(2,2)

  real (kind=real_kind)  :: lon, lat, alpha, beta, tmp
  
  lon=sphere%lon
  lat=sphere%lat
  alpha=alphabeta%x
  beta=alphabeta%y

    if (FaceNum <= 4) then
      lon=lon-(FaceNum-1)*DD_PI/2
      tmp=1.0D0/(cos(lat)*cos(lon))
      Ainv(1,1)=tmp*cos(alpha)*cos(alpha)/cos(lon)
      Ainv(2,1)=0.0D0
      Ainv(1,2)=tmp*cos(beta)*cos(beta)*tan(lat)*tan(lon)
      Ainv(2,2)=tmp*cos(beta)*cos(beta)/cos(lat)
    endif
    if (FaceNum == 5) then
      tmp=1.0D0/(sin(lat)*sin(lat))
      Ainv(1,1)=-tmp*cos(alpha)*cos(alpha)*sin(lat)*cos(lon)
      Ainv(2,1)=tmp*cos(alpha)*cos(alpha)*sin(lon)
      Ainv(1,2)=tmp*cos(beta)*cos(beta)*sin(lat)*sin(lon)
      Ainv(2,2)=tmp*cos(beta)*cos(beta)*cos(lon)
    endif
    if (FaceNum == 6) then
      tmp=1.0D0/(sin(lat)*sin(lat))
      Ainv(1,1)=tmp*cos(alpha)*cos(alpha)*sin(lat)*cos(lon)
      Ainv(2,1)=-tmp*cos(alpha)*cos(alpha)*sin(lon)
      Ainv(1,2)=tmp*cos(beta)*cos(beta)*sin(lat)*sin(lon)
      Ainv(2,2)=tmp*cos(beta)*cos(beta)*cos(lon)
    endif
end subroutine get_Ainv

! initialize global buffers shared by all threads
subroutine spelt_init1(par)
  use edge_mod, only : initghostbuffer,initEdgebuffer
  use parallel_mod, only : parallel_t, haltmp
  
  implicit none
  type (parallel_t) :: par
  
  if (nip .ne. 3) then
     if (par%masterproc) then
        print *, "NUMBER OF INTERPOLATOIN POINTS ERROR in spelt_init1"
        print *, "Only three points for a cell a supported right now!"
     endif
     call haltmp("stopping")
  end if
  if (nhe .ne. 1) then
     if (par%masterproc) then
        print *, "PARAMTER ERROR for fvm: Number of halo zone for the extended"
        print *,"element nhe has to be 1, only this is available now! STOP!"
     endif
     call haltmp("stopping")
  end if
  
  if (ntrac>ntrac_d) then
     if (par%masterproc) print *,'ntrac,ntrac_d=',ntrac,ntrac_d
     call haltmp("PARAMTER ERROR for spelt: ntrac > ntrac_d")
  endif

  call initghostbuffer(cellghostbuf,nlev,ntrac,nipm,nep) !+1 for the air_density, which comes from SE
  call initghostbuffer(factorR,nlev,ntrac,nhe,nc)    ! use the tracer entry
  call initEdgebuffer(edgeveloc,2*nlev)
end subroutine spelt_init1

! initialization that can be done in threaded regions
subroutine spelt_init2(elem,spelt,hybrid,nets,nete,tl)
  use bndry_mod, only: compute_ghost_corner_orientation 
  use derivative_mod, only : derivative_t, derivinit, v2pinit, &
                             interpolate_gll2spelt_points
  implicit none
  type (timelevel_t)                    :: tl
  type (spelt_struct)                    :: spelt(:)
  type (element_t)                      :: elem(:)
  type (hybrid_t)                       :: hybrid
  integer                               :: ie,nets,nete
  type (derivative_t)                         :: deriv           ! derivative struct
  
  call compute_ghost_corner_orientation(hybrid,elem,nets,nete)
  call spelt_grid_init(elem,spelt, nets, nete, tl)
end subroutine spelt_init2

! first communciation of tracers
subroutine spelt_init3(elem,spelt,hybrid,nets,nete,tnp0)
  use edge_mod, only :  ghostVpack2d, ghostVunpack2d
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
  use edge_mod, only :  ghostVpack2d_single, ghostVunpack2d_single,initghostbuffer,freeghostbuffertr
  
  implicit none
  
  type (element_t),intent(inout)            :: elem(:)                 
  type (spelt_struct),intent(inout)         :: spelt(:)  
                  
  type (hybrid_t),intent(in)                :: hybrid                  
                                                                            
  integer,intent(in)                        :: nets,nete,tnp0    
     
  real (kind=real_kind)                     :: ff(nip,nip)                                                                          
  integer                                   :: ie, k, itr, i, j, icell, jcell                 
  
  type (ghostBuffertr_t)                      :: buf
  
  
  ! do it only for FVM tracers, FIRST TRACER will be the AIR DENSITY   
      !first exchange of the initial values
  do ie=nets,nete
    call ghostVpack2d(cellghostbuf,spelt(ie)%c,nipm, nep,nlev,ntrac,0, tnp0, timelevels,elem(ie)%desc)
  end do
  !-----------------------------------------------------------------------------------!  
  call ghost_exchangeV(hybrid,cellghostbuf,nipm,nep,ntrac)
  !-----------------------------------------------------------------------------------!    
  do ie=nets,nete
    call ghostVunpack2d(cellghostbuf,spelt(ie)%c,nipm, nep,nlev,ntrac,0, tnp0, timelevels,elem(ie)%desc)
  end do

! initialize test example
!   do ie=nets,nete  
!     do k=1, nlev
!       do itr=1,ntrac
!         do j=1-nhe,nc+nhe
!           do i=1-nhe,nc+nhe
!             icell=1+(i-1)*nipm
!             jcell=1+(j-1)*nipm
!             ff=spelt(ie)%c(icell:icell+nipm,jcell:jcell+nipm,k,itr,tnp0)
!             ff(2,2)=cip_cell_avr(ff) 
!             spelt(ie)%c(icell+1,jcell+1,k,itr,tnp0)=ff(2,2)
!           enddo
!         enddo
!       enddo
!     enddo
!   enddo

  call initghostbuffer(buf,1,1,nipm,nep)
  do ie=nets,nete
    call ghostVpack2d_single(buf,spelt(ie)%sga,nipm, nep,elem(ie)%desc)
  end do
!-----------------------------------------------------------------------------------! 
  call ghost_exchangeV(hybrid,buf,nipm,nep,1)
!-----------------------------------------------------------------------------------!  
  do ie=nets,nete
    call ghostVunpack2d_single(buf,spelt(ie)%sga,nipm, nep,elem(ie)%desc)
  end do
  call freeghostbuffertr(buf)
end subroutine spelt_init3

subroutine spelt_grid_init(elem,spelt,nets,nete,tl)
  use coordinate_systems_mod, only : ref2sphere, sphere2cubedsphere
  use cube_mod, only : vmap
  use kinds, only : longdouble_kind
  use physical_constants, only : DD_PI
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  use quadrature_mod, only : quadrature_t, gausslobatto
  use parallel_mod, only : haltmp, abortmp
  
    
  implicit none
  type (element_t),intent(inout)            :: elem(:)                 
  type (spelt_struct),intent(inout)          :: spelt(:)
  integer, intent(in)                       :: nets, nete
  type (TimeLevel_t), intent(in)            :: tl              ! time level struct
  
  integer                                   :: cnt, pos

  integer                                   :: ie, i, j, ic,jc, iel,jel
  real (kind=longdouble_kind)               :: xref, yref, dx
  real (kind=real_kind)                     :: tmpD(2,2), detD, tmp, lat, lon, alpha, beta
  type (cartesian2D_t)                      :: alphabeta, tmpab(nep,nep)
  type (cartesian2D_t)                      :: dref  
  
  type(cartesian2d_t)                       :: cart1, cart2, cart3, cart4, cartn(1:nep,1:nep)
  real (kind=real_kind)                       :: cartx(1:nc+1), carty(1:nc+1)
  
  real (kind=real_kind)                     :: dalpha, dbeta
  real(kind=real_kind)               :: pi,pj,qi,qj
  integer                                   :: cubeboundary, nbrsface
  logical                              :: corner
  type (quadrature_t) :: gp   ! Quadrature points and weights on pressure grid

  
  do ie=nets,nete
    spelt(ie)%Facenumber=elem(ie)%FaceNum
    ! for the np grid
    if (np .ne. nc+1) then
      call haltmp("PARAMTER ERROR for SPELT, you are in gll grid point mode, use np = nc+1")
    endif
    gp=gausslobatto(np)
    spelt(ie)%pref=gp%points
    ! for the nc grid
!     dx=2.0D0/(nc)   ! equi-distant grid on reference element in both directions!
!     do j=1,nc+1
!       spelt(ie)%pref(j)=-1+(j-1)*dx
!     end do

    do j=1,nc
      spelt(ie)%dab(j)=abs(elem(ie)%cartp(j+1,1)%x-elem(ie)%cartp(j,1)%x)   ! for np grid
!       spelt(ie)%dab(j)=abs(elem(ie)%corners(1)%x-elem(ie)%corners(2)%x)/nc  ! for nc grid
      do i=1,nc  
        spelt(ie)%drefx(i,j)=abs(spelt(ie)%pref(i+1)-spelt(ie)%pref(i)) 
        spelt(ie)%drefy(i,j)=abs(spelt(ie)%pref(j+1)-spelt(ie)%pref(j))
        do jc=1,nip
          do ic=1,nip  
            xref=spelt(ie)%pref(i)+(ic-1)*spelt(ie)%drefx(i,j)/nipm  
            yref=spelt(ie)%pref(j)+(jc-1)*spelt(ie)%drefy(i,j)/nipm  
            iel=ic+(i-1)*nipm
            jel=jc+(j-1)*nipm
            !define the arrival grid in spherical coordinates
            spelt(ie)%asphere(iel,jel)=ref2sphere(xref,yref,elem(ie)%corners,elem(ie)%FaceNum) 
            alphabeta=sphere2cubedsphere(spelt(ie)%asphere(iel,jel), elem(ie)%FaceNum)
            tmpab(iel,jel)=alphabeta
            spelt(ie)%sga(iel,jel)=metric_term(alphabeta)
        
            lat=spelt(ie)%asphere(iel,jel)%lat
            lon=spelt(ie)%asphere(iel,jel)%lon
        
            alpha=alphabeta%x
            beta=alphabeta%y
            ! Caculation transformation matrix lat/lon to alpha beta for spelt grid
            ! earth radius is assumed to be 1.0D0
            call get_Ainv(elem(ie)%FaceNum, spelt(ie)%asphere(iel,jel), alphabeta,spelt(ie)%Ainv(:,:,iel,jel))
          end do    
        end do
      end do
    end do

    !this can be deleted, once the transformation is on the reference element
    do j=1,np
      do i=1,np
        call vmap(tmpD,elem(ie)%cartp(i,j)%x,elem(ie)%cartp(i,j)%y,elem(ie)%FaceNum)
        detD = tmpD(1,1)*tmpD(2,2) - tmpD(1,2)*tmpD(2,1)      

        spelt(ie)%Dinv(1,1,i,j) =  tmpD(2,2)/detD
        spelt(ie)%Dinv(1,2,i,j) = -tmpD(1,2)/detD
        spelt(ie)%Dinv(2,1,i,j) = -tmpD(2,1)/detD
        spelt(ie)%Dinv(2,2,i,j) =  tmpD(1,1)/detD
      enddo
    enddo
    

    ! this is assuming it is not a refined mesh - and that there is at most one corner element  
    corner=.False.
    cubeboundary=0
    do j=1,8
       cnt = elem(ie)%vertex%nbrs_ptr(j+1)- elem(ie)%vertex%nbrs_ptr(j)
       if (cnt > 0) then
          pos = elem(ie)%vertex%nbrs_ptr(j) !assuming just one - not a refined mesh
          nbrsface=elem(ie)%vertex%nbrs_face(pos)
          ! note that if the element lies on a corner, it will be at j=5,6,7,8
          if ((nbrsface /= elem(ie)%FaceNum) .AND. (j<5)) then
             cubeboundary=j
          endif
       else   ! corner on the cube
          if (.NOT. corner) then
             nbrsface=-1
             cubeboundary=j
             corner=.TRUE.
          else
             print *,'Error in spelt_grid_init: '
             call abortmp('Do not allow one element per face for fvm, please increase ne!')
          endif
       endif
    end do
    

     do i=1,nc
        spelt(ie)%drefx(0,i)=spelt(ie)%drefx(1,i)
        spelt(ie)%drefy(0,i)=spelt(ie)%drefy(1,i) 
        spelt(ie)%drefx(i,0)=spelt(ie)%drefx(i,1)
        spelt(ie)%drefy(i,0)=spelt(ie)%drefy(i,1)
        spelt(ie)%drefx(nc+1,i)=spelt(ie)%drefx(nc,i) 
        spelt(ie)%drefy(nc+1,i)=spelt(ie)%drefy(nc,i) 
        spelt(ie)%drefx(i,nc+1)=spelt(ie)%drefx(i,nc)
        spelt(ie)%drefy(i,nc+1)=spelt(ie)%drefy(i,nc)
     end do
     spelt(ie)%drefx(0,0)=spelt(ie)%drefx(1,1) 
     spelt(ie)%drefy(0,0)=spelt(ie)%drefy(1,1)
     spelt(ie)%drefx(nc+nhe,0)=spelt(ie)%drefx(nc,1) 
     spelt(ie)%drefy(nc+nhe,0)=spelt(ie)%drefy(nc,1)
     spelt(ie)%drefx(0,nc+nhe)=spelt(ie)%drefx(1,nc) 
     spelt(ie)%drefy(0,nc+nhe)=spelt(ie)%drefy(1,nc)
     spelt(ie)%drefx(nc+nhe,nc+nhe)=spelt(ie)%drefx(nc,nc) 
     spelt(ie)%drefy(nc+nhe,nc+nhe)=spelt(ie)%drefy(nc,nc)
          
     do j=1,nc
       do i=1,nc          
          cart1=sphere2cubedsphere(spelt(ie)%asphere(1+(i-1)*nipm,1+(j-1)*nipm), elem(ie)%FaceNum)
          cart2=sphere2cubedsphere(spelt(ie)%asphere(1+(i)*nipm,1+(j-1)*nipm), elem(ie)%FaceNum)
          cart3=sphere2cubedsphere(spelt(ie)%asphere(1+(i-1)*nipm,1+(j)*nipm), elem(ie)%FaceNum)
          cart4=sphere2cubedsphere(spelt(ie)%asphere(1+(i)*nipm,1+(j)*nipm), elem(ie)%FaceNum)
 !          cart1=cartn(1+(i-1)*nipm,1+(j-1)*nipm)
 !          cart2=cartn(1+(i)*nipm,1+(j-1)*nipm)
 !          cart3=cartn(1+(i-1)*nipm,1+(j)*nipm)
 !          cart4=cartn(1+(i)*nipm,1+(j)*nipm)

          spelt(ie)%area_sphere(i,j) =  (I_00(tan(cart4%x),tan(cart4%y)) - I_00(tan(cart3%x),tan(cart3%y)) + &
                             I_00(tan(cart1%x),tan(cart1%y)) - I_00(tan(cart2%x),tan(cart2%y)))


 !           spelt(ie)%area_sphere(i,j)=(I_00(cart4%x,cart4%y) - I_00(cart3%x,cart3%y) + &
 !                   I_00(cart1%x,cart1%y) - I_00(cart2%x,cart2%y))
 !           write(*,*) i, cartx(i+1),cart4%x
 !           spelt(ie)%area_sphere(i,j) =dalpha
 !          spelt(ie)%area_sphere(i,j) = (I_00(cartx(i+1),carty(j+1)) - I_00(cartx(i),carty(j+1)) + &
 !                                                I_00(cartx(i),carty(j)) - I_00(cartx(i+1),carty(j)))                  
       enddo
     enddo
    
  end do
end subroutine spelt_grid_init

function I_00(x,y)
  implicit none
  real (kind=real_kind)                 :: I_00
  real (kind=real_kind), intent(in)     :: x,y

  I_00 = ATAN(x*y/SQRT(1.0D0+x*x+y*y))
end function I_00
!END SUBROUTINE FVM_MCGREGOR----------------------------------------------CE-for FVM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_MCGREGORDSS-----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 26. May 2012                                             !
! DESCRIPTION: ! using McGregor AMS 1993 scheme: Economical Determination of        !
!                Departure Points for Semi-Lagrangian Models                        !
!                McGegror version with DSS every ugradv                             !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine spelt_mcgregordss(elem,spelt,nets,nete, hybrid, deriv, tstep, ordertaylor)
  use derivative_mod, only : derivative_t, ugradv_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  
  implicit none

  type (element_t), intent(inout)             :: elem(:)
  type (spelt_struct), intent(in)              :: spelt(:)

  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)

  type (derivative_t), intent(in)                             :: deriv      ! derivative struct
  real (kind=real_kind), intent(in)                           :: tstep
  integer, intent(in)                                         :: ordertaylor

  real (kind=real_kind), dimension(nets:nete,np,np,2,nlev)    :: ugradv
  real (kind=real_kind), dimension(nets:nete,np,np,2,nlev)    :: vhat
  integer                                                     :: ie, k, order
  real (kind=real_kind), dimension(np,np,2)                   :: ugradvtmp
  real (kind=real_kind)                                       :: timetaylor

    !------------------------------------------------------------------------------------
  timetaylor=1  
  do  order=1,ordertaylor
    timetaylor=-timetaylor*tstep/(order+1)  
    do ie=nets,nete
      if (order==1)then
        ugradv(ie,:,:,:,:)=elem(ie)%derived%vstar(:,:,:,:) 
        vhat(ie,:,:,:,:)=(spelt(ie)%vn0(:,:,:,:) + ugradv(ie,:,:,:,:))/2 
      endif
      do k=1,nlev
        ugradvtmp=ugradv_sphere(vhat(ie,:,:,:,k),ugradv(ie,:,:,:,k),deriv,elem(ie))
        ugradv(ie,:,:,1,k) = elem(ie)%spheremp(:,:)*ugradvtmp(:,:,1) 
        ugradv(ie,:,:,2,k) = elem(ie)%spheremp(:,:)*ugradvtmp(:,:,2) 
      enddo 
      call edgeVpack(edgeveloc,ugradv(ie,:,:,1,:),nlev,0,elem(ie)%desc)
      call edgeVpack(edgeveloc,ugradv(ie,:,:,2,:),nlev,nlev,elem(ie)%desc)
    enddo 
    call bndry_exchangeV(hybrid,edgeveloc)
    do ie=nets,nete
       call edgeVunpack(edgeveloc,ugradv(ie,:,:,1,:),nlev,0,elem(ie)%desc)
       call edgeVunpack(edgeveloc,ugradv(ie,:,:,2,:),nlev,nlev,elem(ie)%desc)
       do k=1, nlev  
         ugradv(ie,:,:,1,k)=ugradv(ie,:,:,1,k)*elem(ie)%rspheremp(:,:)
         ugradv(ie,:,:,2,k)=ugradv(ie,:,:,2,k)*elem(ie)%rspheremp(:,:)
         elem(ie)%derived%vstar(:,:,:,k)=elem(ie)%derived%vstar(:,:,:,k) + timetaylor*ugradv(ie,:,:,:,k)
       end do
    end do
  end do  

end subroutine spelt_mcgregordss
!END SUBROUTINE FVM_MCGREGORDSS-------------------------------------------CE-for FVM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE SPELT_RKDSS-----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR, 06. December 2012                                             !
! DESCRIPTION: ! create a runge kutta taylor serios mixture to calculate the departure grid                            !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine spelt_rkdss(elem,spelt,nets,nete, hybrid, deriv, tstep, ordertaylor)
  use derivative_mod, only : derivative_t, ugradv_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  
  implicit none

  type (element_t), intent(inout)                :: elem(:)
  type (spelt_struct), intent(in)              :: spelt(:)

  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)

  type (derivative_t), intent(in)                             :: deriv      ! derivative struct
  real (kind=real_kind), intent(in)                           :: tstep
  integer, intent(in)                                         :: ordertaylor

  integer                                                     :: ie, k, order
  real (kind=real_kind), dimension(np,np,2)                   :: ugradvtmp
  real (kind=real_kind)                                       :: timetaylor
!
! RK-SSP 2 stage 2nd order:
!     x*(t+1) = x(t) + U(x(t),t)                          
!     x(t+1) = x(t) +  1/2 ( U(x*(t+1),t+1) + U(x(t),t) )       
! apply taylor series:
!  U(x*(t+1),t+1) = U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1)
!
!  (x(t+1)-x(t))/dt =  1/2( U(x(t),t+1)+U(x(t),t)) + dt 1/2 U(x(t),t) gradU(x(t),t+1)  
!
! suppose dt = -tstep (we go backward)
!  (x(t-tstep)-x(t))/-tstep =  1/2( U(x(t),t-tstep)+U(x(t),t)) - tstep 1/2 U(x(t),t) gradU(x(t),t-tstep)  
!
!  x(t-tstep) = x(t)) -tstep * [ 1/2( U(x(t),t-tstep)+U(x(t),t)) - tstep 1/2 U(x(t),t) gradU(x(t),t-tstep) ]  
!
!    !------------------------------------------------------------------------------------
    do ie=nets,nete
       ! vn0 = U(x,t)
       ! vstar = U(x,t+1)
      do k=1,nlev
!         ugradvtmp(:,:,:)=ugradv_sphere(fvm(ie)%vn0(:,:,:,k),elem(ie)%derived%vstar(:,:,:,k),deriv,elem(ie))
        ugradvtmp(:,:,:)=ugradv_sphere(elem(ie)%derived%vstar(:,:,:,k),spelt(ie)%vn0(:,:,:,k),deriv,elem(ie))
        
        elem(ie)%derived%vstar(:,:,:,k) = &
             (elem(ie)%derived%vstar(:,:,:,k) + spelt(ie)%vn0(:,:,:,k))/2   - tstep*ugradvtmp(:,:,:)/2

        elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%spheremp(:,:)
        elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%spheremp(:,:)
      enddo 
      call edgeVpack(edgeveloc,elem(ie)%derived%vstar,2*nlev,0,elem(ie)%desc)
    enddo 
    call bndry_exchangeV(hybrid,edgeveloc)
    do ie=nets,nete
       call edgeVunpack(edgeveloc,elem(ie)%derived%vstar,2*nlev,0,elem(ie)%desc)
       do k=1, nlev  
          elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%rspheremp(:,:)
          elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%rspheremp(:,:)
       end do
    end do

end subroutine spelt_rkdss
!END SUBROUTINE SPELT_rkdss-------------------------------------------CE-for FVM!

subroutine cell_search(elem, spelt, dsphere, icell, jcell,dref, alphabeta, face_nodep) 

  implicit none

  type (spherical_polar_t), intent(in)     :: dsphere

  type (element_t), intent(in)             :: elem
  type (spelt_struct), intent(in)          :: spelt
  
  integer, intent(out)                     :: icell, jcell
  type (cartesian2D_t),intent(out)         :: dref
  type (cartesian2D_t), intent(out)        :: alphabeta
    
  type (cartesian2D_t)                     :: dcart
  integer                                  :: number, endi, starti, tmpi
  real (kind=real_kind)                    :: refnc(1:nc+1), dxcell, dycell, tmp

  integer                                  :: i, face_nodep      
  integer                                  :: tmp_i 
  integer                                  :: j, cnt, loc, nbr, pos


  call cube_facepoint_ne(dsphere,ne,dcart, alphabeta, number, face_nodep)
    ! Search index along "x"  (bisection method)
  starti = 1
  endi = nc+1
  do
    if  ((endi-starti) <=  1)  exit
    tmpi = (endi + starti)/2
    if (dcart%x  >  spelt%pref(tmpi)) then
      starti = tmpi
    else
      endi = tmpi
    endif
  enddo
  icell = starti
  dxcell=spelt%drefx(icell,1) !abs(spelt%pref(icell+1)-spelt%pref(icell))
  ! Search index along "y"
  starti = 1
  endi = nc+1
  do
    if  ((endi-starti) <=  1)  exit
    tmpi = (endi + starti)/2
    if (dcart%y  >  spelt%pref(tmpi)) then
      starti = tmpi
    else
      endi = tmpi
    endif
  enddo
  jcell = starti
  dycell=spelt%drefy(1,jcell) !abs(spelt%pref(jcell+1)-spelt%pref(jcell))
    
  dref%x=dcart%x-spelt%pref(icell)
  dref%y=dcart%y-spelt%pref(jcell)
!     
  if ((dref%x<-1.0D-12) .or.(dref%y<-1.0D-12) .or.(dref%x>dxcell+1.0D-12) .or. (dref%y>dycell+1.0D-12) ) then
    write(*,*) 'Something is wrong in search!'
    write(*,*) dref%x,dxcell, dref%y,dycell
    stop
  endif
    
  if ((icell<1) .or.(icell>nc) .or. (jcell<1) .or. (jcell>nc)) then
    write(*,*) 'icell, jcell,Something is wrong in search!'
    stop
  endif
  
  !NOTE (AB): the following code assumes that we are not using a refined mesh - so there is 
  ! at most one corner element 

  !need to find which of the 8 locations has the nbr match for the number
  loc = 0
  do j=1,8
     cnt = elem%vertex%nbrs_ptr(j+1)- elem%vertex%nbrs_ptr(j)
     if (cnt > 0) then
        pos = elem%vertex%nbrs_ptr(j) !assuming just one - not a refined mesh
        nbr = elem%vertex%nbrs(pos)
        if (number == nbr) then
           loc = j
           exit
        endif
     end if
  end do


  if(loc == 1) then  !west
    if ((elem%FaceNum<=4) .or. (face_nodep==elem%FaceNum)) then
      icell=icell-nc
    elseif (elem%FaceNum==5) then
      tmpi=icell
      icell=-jcell+1
      jcell=tmpi
      tmp=dref%x
      dref%x=dycell-dref%y
      dref%y=tmp
    elseif (elem%FaceNum==6) then
      tmpi=icell
      icell=jcell-nc
      jcell=nc+1-tmpi
      tmp=dref%x
      dref%x=dref%y
      dref%y=dxcell-tmp
    endif
  end if   
         
  if(loc == 2) then   !east
    if ((elem%FaceNum<=4).or. (face_nodep==elem%FaceNum)) then
      icell=icell+nc
    elseif (elem%FaceNum==5) then
      tmpi=icell
      icell=jcell+nc
      jcell=nc+1-tmpi
      tmp=dref%x
      dref%x=dref%y
      dref%y=dxcell-tmp
    elseif (elem%FaceNum==6) then
      tmpi=icell
      icell=nc+(nc+1-jcell)
      jcell=tmpi
      tmp=dref%x
      dref%x=dycell-dref%y
      dref%y=tmp
    endif
  end if   
         
  if(loc == 3) then   !south
    if ((elem%FaceNum==1) .OR. (elem%FaceNum==6) .or. (face_nodep==elem%FaceNum)) then
      jcell=jcell-nc
    elseif (elem%FaceNum==2) then
      tmpi=icell
      icell=nc+1-jcell
      jcell=tmpi-nc
      tmp=dref%x
      dref%x=dycell-dref%y
      dref%y=tmp
    elseif ((elem%FaceNum==3) .or. (elem%FaceNum==5)) then
      icell=nc+1-icell
      jcell=-jcell+1
      dref%x=dxcell-dref%x
      dref%y=dycell-dref%y
    elseif (elem%FaceNum==4) then
      tmpi=icell
      icell=jcell
      jcell=1-tmpi
      tmp=dref%x
      dref%x=dref%y
      dref%y=dxcell-tmp
    endif
  end if   
         
  if(loc == 4) then   !north
    if ((elem%FaceNum==1) .OR. (elem%FaceNum==5) .or. (face_nodep==elem%FaceNum)) then
      jcell=nc+jcell
    elseif (elem%FaceNum==2) then
      tmpi=icell
      icell=jcell
      jcell=nc+(nc+1-tmpi)
      tmp=dref%x
      dref%x=dref%y
      dref%y=dxcell-tmp
    elseif ((elem%FaceNum==3) .or. (elem%FaceNum==6)) then
      icell=nc+1-icell
      jcell=nc+(nc+1-jcell)
      dref%x=dxcell-dref%x
      dref%y=dycell-dref%y
    elseif (elem%FaceNum==4) then
      tmpi=icell
      icell=1+(nc-jcell)
      jcell=nc+tmpi
      tmp=dref%x
      dref%x=dycell-dref%y
      dref%y=tmp
    endif
  end if   

  ! cases southwest, southeast, northwest, northeast on a cube corner do not exist
  if(loc == 5) then   !southwest
    if ((face_nodep==elem%FaceNum) .or. (elem%FaceNum==1)) then
      icell=icell-nc
      jcell=jcell-nc
    elseif (elem%FaceNum==2) then
      if (face_nodep==1) then
        icell=icell-nc
        jcell=jcell-nc
      else
        tmpi=icell
        icell=1-jcell
        jcell=-nc+tmpi
        tmp=dref%x
        dref%x=dycell-dref%y
        dref%y=tmp
      endif
    elseif (elem%FaceNum==3) then
      if (face_nodep==2) then
        icell=icell-nc
        jcell=jcell-nc
      else          
        icell=1-icell
        jcell=1-jcell
        dref%x=dxcell-dref%x
        dref%y=dycell-dref%y
      endif
    elseif (elem%FaceNum==4) then
      if (face_nodep==3) then
        icell=icell-nc
        jcell=jcell-nc
      else
        tmpi=icell
        icell=-nc+jcell
        jcell=-tmpi+1
        tmp=dref%x
        dref%x=dref%y
        dref%y=dxcell-tmp
      endif
    elseif (elem%FaceNum==5) then
      if (face_nodep==4) then
        tmpi=icell
        icell=1-jcell
        jcell=tmpi-nc
        tmp=dref%x
        dref%x=dycell-dref%y
        dref%y=tmp
      else
        icell=1-icell
        jcell=1-jcell
        dref%x=dxcell-dref%x
        dref%y=dycell-dref%y
      endif
    elseif (elem%FaceNum==6) then
      if (face_nodep==1) then          
        icell=icell-nc
        jcell=jcell-nc
      else
        tmpi=icell
        icell=-nc+jcell
        jcell=-tmpi+1
        tmp=dref%x
        dref%x=dref%y
        dref%y=dxcell-tmp  
      endif
    endif
  endif   
       
  if(loc == 6) then   !southeast
    if ((elem%FaceNum==1) .or. (face_nodep==elem%FaceNum)) then
      icell=nc+icell
      jcell=jcell-nc
    elseif (elem%FaceNum==2) then
      if (face_nodep==3) then
        icell=nc+icell
        jcell=jcell-nc
      else
        tmpi=icell
        icell=nc+(nc+1-jcell)
        jcell=-nc+tmpi
        tmp=dref%x
        dref%x=dycell-dref%y
        dref%y=tmp  
      endif
    elseif (elem%FaceNum==3) then
      if (face_nodep==4) then
        icell=nc+icell
        jcell=jcell-nc
      else
        icell=nc+(nc+1-icell)
        jcell=1-jcell
        dref%x=dxcell-dref%x
        dref%y=dycell-dref%y
      endif
    elseif (elem%FaceNum==4) then
      if (face_nodep==1) then
        icell=nc+icell
        jcell=jcell-nc  
      else
        tmpi=icell
        icell=nc+jcell
        jcell=-tmpi+1
        tmp=dref%x
        dref%x=dref%y
        dref%y=dxcell-tmp
      endif
    elseif (elem%FaceNum==5) then
      if (face_nodep==2) then
        tmpi=icell
        icell=nc+jcell
        jcell=-tmpi+1
        tmp=dref%x
        dref%x=dref%y
        dref%y=dxcell-tmp
      else
        icell=nc+(nc+1-icell)
        jcell=1-jcell
        dref%x=dxcell-dref%x
        dref%y=dycell-dref%y
      endif
    elseif (elem%FaceNum==6) then
      if(face_nodep==1) then
        icell=nc+icell
        jcell=jcell-nc
      else
        tmpi=icell
        icell=nc+(nc+1-jcell)
        jcell=-nc+tmpi
        tmp=dref%x
        dref%x=dycell-dref%y
        dref%y=tmp
      endif
    endif
  endif   

  if(loc == 7) then   !northwest
    if ((elem%FaceNum==1) .or. (face_nodep==elem%FaceNum)) then
      icell=icell-nc
      jcell=nc+jcell
    elseif (elem%FaceNum==2) then
      if (face_nodep==1) then
        icell=icell-nc
        jcell=nc+jcell
      else
        tmpi=icell
        icell=jcell-nc
        jcell=nc+(nc+1-tmpi)
        tmp=dref%x
        dref%x=dref%y
        dref%y=dxcell-tmp
      endif
    elseif (elem%FaceNum==3) then
      if (face_nodep==2) then
        icell=icell-nc
        jcell=nc+jcell
      else
        icell=1-icell
        jcell=nc+(nc+1-jcell)
        dref%x=dxcell-dref%x
        dref%y=dycell-dref%y
      endif
    elseif (elem%FaceNum==4) then
      if (face_nodep==3) then
        icell=icell-nc
        jcell=nc+jcell
      else
        tmpi=icell
        icell=-jcell+1
        jcell=nc+tmpi
        tmp=dref%x
        dref%x=dycell-dref%y
        dref%y=tmp
      endif
    elseif (elem%FaceNum==5) then
      if (face_nodep==1) then
        icell=icell-nc
        jcell=nc+jcell
      else
        tmpi=icell
        icell=-jcell+1
        jcell=nc+tmpi
        tmp=dref%x
        dref%x=dycell-dref%y
        dref%y=tmp   
      endif
    elseif (elem%FaceNum==6) then
      if (face_nodep==4) then
        tmpi=icell
        icell=jcell-nc
        jcell=nc+(nc+1-tmpi)
        tmp=dref%x
        dref%x=dref%y
        dref%y=dxcell-tmp
      else
        icell=1-icell
        jcell=nc+(nc+1-jcell)
        dref%x=dxcell-dref%x
        dref%y=dycell-dref%y
      endif
    endif
  end if   
       
  if(loc == 8) then   !northeast
    if ((elem%FaceNum==1) .or. (face_nodep==elem%FaceNum)) then
      icell=nc+icell
      jcell=nc+jcell
    elseif (elem%FaceNum==2) then
      if (face_nodep==3) then
        icell=nc+icell
        jcell=nc+jcell
      else
        tmpi=icell
        icell=nc+jcell
        jcell=nc+(nc+1-tmpi)
        tmp=dref%x
        dref%x=dref%y
        dref%y=dxcell-tmp 
      endif
    elseif (elem%FaceNum==3) then
      if (face_nodep==4) then
        icell=nc+icell
        jcell=nc+jcell
      else
        icell=nc+(nc+1-icell)
        jcell=nc+(nc+1-jcell)
        dref%x=dxcell-dref%x
        dref%y=dycell-dref%y
      endif
    elseif (elem%FaceNum==4) then
      if (face_nodep==1) then
        icell=nc+icell
        jcell=nc+jcell
      else
        tmpi=icell
        icell=nc+(nc+1-jcell)
        jcell=nc+tmpi
        tmp=dref%x
        dref%x=dycell-dref%y
        dref%y=tmp
      endif
    elseif (elem%FaceNum==5) then
      if (face_nodep==1) then
        icell=nc+icell
        jcell=nc+jcell
      else
        tmpi=icell
        icell=nc+jcell
        jcell=nc+(nc+1-tmpi)
        tmp=dref%x
        dref%x=dref%y
        dref%y=dxcell-tmp
      endif
    elseif (elem%FaceNum==6) then
      if (face_nodep==2) then
        tmpi=icell
        icell=nc+(nc+1-jcell)
        jcell=nc+tmpi
        tmp=dref%x
        dref%x=dycell-dref%y
        dref%y=tmp
      else
        icell=nc+(nc+1-icell)
        jcell=nc+(nc+1-jcell)
        dref%x=dxcell-dref%x
        dref%y=dycell-dref%y
      endif
    endif
  end if


  if ((icell<0) .or.(icell>nc+1) .or. (jcell<0) .or. (jcell>nc+1)) then
    write(*,*) '2 Something is wrong in search!'
    tmp_i = elem%vertex%nbrs_ptr(9) - 1
    write(*,*) number, elem%vertex%nbrs(1:tmp_i)
    write(*,*) icell, jcell, elem%GlobalId, elem%FaceNum, face_nodep
    stop
  endif
  !if ((dref%x<-1.0D-12) .or.(dref%y<-1.0D-12) .or.(dref%x>dxcell+1.0D-12) .or. (dref%y>dxcell+1.0D-12) ) then
  !  write(*,*) '3 Something is wrong in search!'
  !   tmp_i = elem%vertex%nbrs_ptr(9) - 1

  !  write(*,*) number, elem%vertex%nbrs(1:tmp_i)
  !  write(*,*) dref
  !  stop
  !endif

end subroutine cell_search


! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_DEP_FROM_GLL----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR 14. December 2011                            !
! DESCRIPTION: calculates the deparute grid for spelt coming from the gll points    !
!                                                                                   !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine spelt_dep_from_gll(elem, deriv, asphere,dsphere,dt,tl,klev)
  use physical_constants, only : DD_PI, rearth
  use coordinate_systems_mod, only : cartesian3D_t, change_coordinates
  use derivative_mod, only : derivative_t, interpolate_gll2spelt_points

  implicit none
  type (element_t), intent(in)          :: elem
  type (derivative_t)  , intent(in)     :: deriv
  type (spherical_polar_t),intent(in)   :: asphere(nep,nep)
  type (spherical_polar_t),intent(out)  :: dsphere(nep,nep)
  type (timelevel_t),intent(in)         :: tl
  real (kind=real_kind),intent(in)      :: dt
  integer,intent(in)                    :: klev
  
  real (kind=real_kind)                 :: uxyz_gll(np,np,3),uxyz(nep,nep,3)
  real (kind=real_kind)                 :: un,ue,ur,clon,clat,slon,slat
  type(cartesian3D_t)                   :: acart
  integer                               :: i,j
  
  real (kind=real_kind)  :: vstar(np,np,2)
  
  ! convert velocity from lat/lon to cartesian 3D
  vstar = elem%derived%vstar(:,:,:,klev)
    
  do i=1,np
    do j=1,np
      clon = cos(elem%spherep(i,j)%lon)
      slon = sin(elem%spherep(i,j)%lon)
      clat = cos(elem%spherep(i,j)%lat)
      slat = sin(elem%spherep(i,j)%lat)

      ur = 0
      ue = vstar(i,j,1) 
      un = vstar(i,j,2)
      uxyz_gll(i,j,1)= clon*clat*ur - clon*slat*un - slon*ue
      uxyz_gll(i,j,2)= slon*clon*ur - slon*slat*un + clon*ue
      uxyz_gll(i,j,3)= slat          *ur + clat          *un    
    enddo
  enddo
  ! interpolate velocity to spelt nodes
  do i=1,3
    uxyz(:,:,i)=interpolate_gll2spelt_points(uxyz_gll(:,:,i),deriv)
  end do 
  ! compute departure point 
  ! crude, 1st order accurate approximation.  to be improved
  do i=1,nep
    do j=1,nep
      ! note: asphere()%r=1, so we need to convert velocity to radians/sec:
      acart = change_coordinates(asphere(i,j))  
      acart%x = acart%x - dt*uxyz(i,j,1)/rearth
      acart%y = acart%y - dt*uxyz(i,j,2)/rearth
      acart%z = acart%z - dt*uxyz(i,j,3)/rearth
      dsphere(i,j) = change_coordinates(acart)
      dsphere(i,j)%r = asphere(i,j)%r
    enddo
  enddo
end subroutine spelt_dep_from_gll
!END SUBROUTINE FVM_DEP_FROM_GLL------------------------------------------CE-for FVM!

function metric_term(alphabeta) result(sg) 
  implicit none
  type (cartesian2D_t), intent(in)      :: alphabeta
  real (kind=real_kind)                 :: sg 

  real (kind=real_kind) :: rr, ca,cb,ta,tb, cos2ab, rho2, rho3
  real (kind=real_kind) :: ta2, tb2                

  !Square root of G-determinent  [sqrt(|g_ij|) = g = R^2/(rho^3*cab)]

  rr = 1.0D0   !*sqrt(3)        !radius^2 for unit cubed-sphere
  ca = cos(alphabeta%x)
  cb = cos(alphabeta%y)
  ta = tan(alphabeta%x)
  tb = tan(alphabeta%y)

  ta2 = ta*ta 
  tb2 = tb*tb 
  rho2 = (1.0D0 + ta2 + tb2)

  cos2ab = (ca*cb)**2
  rho3 = (rho2)**(1.50D0)
  sg =  rr / (rho3 * cos2ab)

end function metric_term


function cell_minmax(fel) result(minmax)  
  implicit none
  real (kind=real_kind), intent(in) :: fel(nip,nip) 
  real (kind=real_kind) :: vmin, vmax, val, minmax(2) 
  integer :: i,j              

  vmin = 1.0D19
  vmax = - 1.0D19
  do j = 1, nip 
    do i = 1, nip 
      val = fel(i,j) 
      if (val < vmin) then
        vmin=val
      endif
      if (val > vmax) then
        vmax=val
      endif
!        vmin = min(vmin,val) 
!        vmax = max(vmax,val) 
    enddo
  enddo
  minmax(1) = vmin 
  minmax(2) = vmax 
end function cell_minmax

!  QMSL filter for given point 
function qmsl_cell_filter(ix,iy,minmax,val) result(qval) 
  implicit none
  integer, intent(in)                 :: ix,iy
  real (kind=real_kind), intent(in)   :: minmax(1-nhe:nc+nhe,1-nhe:nc+nhe,2), val
  real (kind=real_kind)               :: qval

  qval=val
  if (val < minmax(ix,iy,1)) then 
    qval = minmax(ix,iy,1) 
  endif 
  if (val > minmax(ix,iy,2)) then 
    qval = minmax(ix,iy,2) 
  endif 
end function qmsl_cell_filter

function cip_interpolate(cf,dxi,det) result(val) 
  implicit none
  real (kind=real_kind), intent(in) :: dxi,det 
  real (kind=real_kind), intent(in), dimension(nip,nip) :: cf 
  real (kind=real_kind) :: val, dxidet,d2xi,d2et, d2xid2et, d2xidet,dxid2et

     d2xi = dxi*dxi          
     d2et = det*det          
     dxidet = dxi*det          
     d2xidet = dxi* dxidet          
     dxid2et = det* dxidet          
     d2xid2et = d2et* d2xi
  
     val = cf(1,1)  + cf(2,1) *  dxi  + cf(1,2)*  det   + cf(2,2)* dxidet  &
                    + cf(3,1) *  d2xi + cf(1,3)*  d2et  + cf(3,2)* d2xidet &
                    + cf(2,3)*dxid2et + cf(3,3)* d2xid2et
end function cip_interpolate


subroutine cip_coeff(dxcell,dycell,ff,avr,cf) 
  use physical_constants, only : DD_PI

  implicit none
  real (kind=real_kind) , intent(in), dimension(nip,nip) :: ff
  real (kind=real_kind) , intent(in) :: avr 
  real (kind=real_kind) , intent(out) :: cf(nip,nip)  
  real (kind=real_kind) , intent(in)  :: dxcell, dycell
  integer ::  i,j
!   write(*,*) dxcell, dycell
!   d(1) = 1.0D0/dxcell
!   d(2) = 2.0D0/(dxcell*dxcell) 
!   d(3) = 1.0D0/dycell
!   d(4) = 2.0D0/(dycell*dycell)          
!   d(5) = 2.0D0/(dxcell*dycell)          
!   d(6) = 1.0D0/(dxcell*dxcell*dycell)   
!   d(7) = 1.0D0/(dxcell*dycell*dycell)          
!   d(8) = 3.0D0/(dxcell*dxcell*dycell*dycell)

  cf(1,1) = ff(1,1)
  cf(2,1) =  (-3.0D0*ff(1,1) + 4.0D0*ff(2,1) - ff(3,1) )/dxcell 
  cf(3,1) = 2.0D0 * (       ff(1,1) - 2.0D0*ff(2,1) + ff(3,1) )/(dxcell*dxcell) 

  cf(1,2) = 1.0D0 * (-3.0D0*ff(1,1) + 4.0D0*ff(1,2) - ff(1,3) )/dycell 
  cf(1,3) = 2.0D0 * (       ff(1,1) - 2.0D0*ff(1,2) + ff(1,3) )/(dycell*dycell) 

  cf(2,2) = 2.0D0 * ( 4.0D0*(ff(1,1) - ff(2,3) - ff(3,2))  + ff(1,3) &
                     + ff(3,1) - 8.0D0*(ff(2,1) + ff(1,2))  + 18.0D0*avr )/(dxcell*dycell) 

  cf(3,2) = 1.0D0 * (-5.0D0*(ff(1,1) + ff(3,1)) + 12.0D0*(ff(1,2) + ff(3,2))  - ff(1,3) &
                     + 8.0D0*ff(2,3)  + 16.0D0*ff(2,1) - ff(3,3) - 36.0D0*avr )/(dxcell*dxcell*dycell) 

  cf(2,3) = 1.0D0 * (-5.0D0*(ff(1,1) + ff(1,3)) + 12.0D0*(ff(2,1)+ ff(2,3)) &
                     + 16.0D0*ff(1,2) - ff(3,1) + 8.0D0*ff(3,2) - ff(3,3) - 36.0D0*avr )/(dxcell*dycell*dycell) 

  cf(3,3) = 3.0D0 * (-4.0D0*(ff(1,2) + ff(2,1) + ff(2,3) + ff(3,2))   &
                     + ff(1,1) + ff(1,3) + ff(3,1) + ff(3,3) + 12.0D0*avr )/(dxcell*dxcell*dycell*dycell) 
end  subroutine cip_coeff

function cip_cell_avr(ff)  result(qint)
  implicit none
  real (kind=real_kind), intent(in), dimension(nip,nip) :: ff
  real (kind=real_kind) :: qint
  real (kind=real_kind) :: f(3)       
  integer ::   j

!  delta = wel(1) * wel(1) * 0.25D0
  qint = 0.0D0

! cell-wise integral  by Simpson's 3pt rule 
  do j = 1, 3  
    f(j) = ff(1,j) + 4.0D0*ff(2,j) + ff(3,j) 
  enddo 
  qint = (f(1) + 4.0D0 * f(2) + f(3))/36.0D0 
end function cip_cell_avr


subroutine cube_facepoint_ne(sphere,ne,cart, cube, number, face_no)
  use coordinate_systems_mod, only : sphere2cubedsphere
  use physical_constants, only : DD_PI
  use parallel_mod,           only : abortmp
  
  implicit none

  type (spherical_polar_t), intent (in) :: sphere
  integer             , intent(in)      :: ne
  type (cartesian2D_t), intent(out)     :: cart
  type (cartesian2D_t), intent(out)     :: cube
  
  integer             , intent(out)     :: number, face_no

  real (kind=real_kind) :: xp,yp
  integer               :: ie, je
  real (kind=real_kind) :: x1,x2
  real (kind=real_kind) :: dx

  face_no = panel_index(sphere%lon,sphere%lat)
  
  cube    = sphere2cubedsphere(sphere, face_no)
  xp      = cube%x
  yp      = cube%y

  ! MNL: for uniform grids (on cube face), analytic solution is fine
  x1 = xp + 0.25D0*DD_PI
  x2 = yp + 0.25D0*DD_PI

  dx = (0.5D0*DD_PI)/ne
  ie = INT(ABS(x1)/dx)
  je = INT(ABS(x2)/dx)
  ! if we are exactly on an element edge, we can put the point in
  ! either the ie or ie+1 element, EXCEPT if ie==ne.
  if ( ABS(x1)/dx < ne  ) ie=ie+1
  if ( ABS(x2)/dx < ne  ) je=je+1
  if (ie>ne .or. je>ne) then
     write(*,*)'ERROR: ',ie,je,ne
     write(*,*)'lat,lon=',sphere%lat,sphere%lon
     write(*,*)'face no=',face_no
     write(*,*)xp,yp,dx,x1,x2,x1/dx,x2/dx
     call abortmp('cube_facepoint_ne: bad argument')
  endif

  ! bug fix MT 1/2009.  This was creating a plotting error at
  ! the row of elements in iface=2 at 50 degrees (NE=16 128x256 lat/lon grid)
  ! For point on element edge, we can have ie=2, but x1=dx
  ! but if ie>1, we must execute this statement.
  ! The only time we can skip this statement is if ie=1, but then
  ! the statement has no effect, so lets never skip it:
  !    if (x1 > dx ) then
  x1 = x1 - dble(ie-1)*dx
  !    endif

  x1 = 2.0D0*(x1/dx)-1.0D0

  !    if (x2 > dx ) then    ! removed MT 1/2009, see above
  x2 = x2 - dble(je-1)*dx
  !    endif

  x2 = 2.0D0*(x2/dx)-1.0D0

  ! coordinates within an element [-1,1]
  cart%x = x1
  cart%y = x2
  number = ie + (je-1)*ne + (face_no-1)*ne*ne
end subroutine cube_facepoint_ne

!================================================
!  (Nair) Cube face index and local coordinates
!================================================

function panel_index(lam,the) result(iface)  
implicit none
real (kind=real_kind), intent(in) :: lam,the 
real (kind=real_kind) :: xx,yy,zz, pm 
integer :: ix,iy,iz, iface 

        xx  = cos(lam)*cos(the)
        yy  = sin(lam)*cos(the)
        zz  = sin(the)

        pm = max(abs(xx),abs(yy),abs(zz))

        ix = Integer_Part(xx,pm)
        iy = Integer_Part(yy,pm)
        iz = Integer_Part(zz,pm)

! Panel index determination

       if (iz  ==  1) then
           iface = 6
       elseif (iz  == -1) then
             iface = 5
       elseif ((ix == 1).and.(iy /= 1)) then
                iface = 1
       elseif ((ix == -1).and.(iy /= -1)) then
                iface = 3
       elseif ((iy == 1).and.(ix /= -1)) then
            iface = 2
       elseif ((iy == -1).and.(ix /=  1)) then
          iface = 4
       endif

end function panel_index

function integer_part(a,b)  result(ii) 
   implicit none
   real (kind=real_kind), Intent(in) :: a,b   
   real (kind=real_kind) :: tol,a1,b1 
   Integer :: ii
  
     tol = 1.0D-14

      a1 = abs(a)
      b1 = abs((b - a1)/b)

      if (b1 < tol )  then
         if (a < 0.0D0)  then
            ii = -1 
          else
            ii =  1 
         endif 
       else
         ii = int(a/b)
      endif 

end  function integer_part


subroutine solidbody_all(spelt, dsphere1,dsphere2,contrauv,k)
  use kinds, only : real_kind
  use time_mod, only : tstep, nmax, Time_at
  use physical_constants, only : DD_PI, rearth
  use coordinate_systems_mod, only : spherical_polar_t

  implicit none
  type (spelt_struct), intent(inout)          :: spelt
  type (spherical_polar_t),intent(out)        :: dsphere1(nep,nep),dsphere2(nep,nep)
  real (kind=real_kind), intent(out)          :: contrauv(nep,nep,2,3)
  
  integer, intent(in)                         :: k
! for test case solid-body rotation on the sphere with alpha
  real (kind=real_kind)                 :: alpha, omega
  real (kind=real_kind)                 :: lap,thp,lamrot, therot,tmplondep,tmplatdep
  real (kind=real_kind)                 ::tmpflux
  integer ie,i,j
  real (kind=real_kind)                 :: lon, lat, u0, u, v


  ! set values for solid-body rotation on the sphere with alpha, this should be 
  ! outside 
  alpha=0.0D0 !DD_PI/4!-0.9*DD_PI/4.0D0 !DD_PI/4  !DD_PI/4 !1.3!0.78
!   omega=2*DD_PI/Time_at(nmax)          ! angular velocity: around the earth
  
  omega=2*DD_PI/1036800                !in 12 days around the earth
  lap=DD_PI
  thp=DD_PI/2-alpha

  u0=2*DD_PI/dble(1036800)

  do j=1,nep
    do i=1,nep
      lon = spelt%asphere(i,j)%lon
      lat = spelt%asphere(i,j)%lat
      if (alpha==0.0D0) then      ! move along the equator
        tmplondep=lon-omega*tstep/2.0D0
        if (tmplondep<0.0D0) then
          tmplondep=tmplondep+2*DD_PI
        endif
        tmplatdep=lat
      else            
        ! rotate sphere with the pole (lap,thp) with respect to the unrotated system
      call vortex_rotatedsphere(lap,thp,lon,lat,lamrot,therot)
        lamrot=lamrot-omega*tstep/2.0D0
    !     if (lamrot<0.0D0) lamrot=lamrot+2*DD_PI        
      call vortex_rotatedsphereback(lap,thp,lamrot,therot,&
                                   tmplondep,tmplatdep)                                                      
      endif
      dsphere1(i,j)%lon=tmplondep
      dsphere1(i,j)%lat=tmplatdep
      dsphere1(i,j)%r=spelt%asphere(i,j)%r

      if (alpha==0.0D0) then      ! move along the equator
        tmplondep=lon-omega*tstep
        if (tmplondep<0.0D0) then
          tmplondep=tmplondep+2*DD_PI
        endif
        tmplatdep=lat
      else            
        ! rotate sphere with the pole (lap,thp) with respect to the unrotated system
      call vortex_rotatedsphere(lap,thp,lon,lat,lamrot,therot)
        lamrot=lamrot-omega*tstep
    !     if (lamrot<0.0D0) lamrot=lamrot+2*DD_PI        
      call vortex_rotatedsphereback(lap,thp,lamrot,therot,&
                                   tmplondep,tmplatdep)                                                      
      endif
      dsphere2(i,j)%lon=tmplondep
      dsphere2(i,j)%lat=tmplatdep
      dsphere2(i,j)%r=spelt%asphere(i,j)%r

      u = u0*(cos(alpha)*cos(lat)+sin(alpha)*cos(lon)*sin(lat))
      v = -u0*sin(alpha)*sin(lon)
      ! convert from radians per dimensionless time to 
      ! meters/sec 
!       vstar(i,j,1)=u !* rearth /( 12*3600*24/5)
!       vstar(i,j,2)=v !* rearth /( 12*3600*24/5)
      ! transform to contravariant (alpha/beta velocities)
      contrauv(i,j,1,1)=spelt%Ainv(1,1,i,j)*u+spelt%Ainv(2,1,i,j)*v      
      contrauv(i,j,2,1)=spelt%Ainv(1,2,i,j)*u+spelt%Ainv(2,2,i,j)*v
      contrauv(i,j,1,2)=contrauv(i,j,1,1)
      contrauv(i,j,2,2)=contrauv(i,j,2,1)
      contrauv(i,j,1,3)=contrauv(i,j,1,1)
      contrauv(i,j,2,3)=contrauv(i,j,2,1)

    enddo
  enddo
end subroutine solidbody_all



!-----------------------------------------------------------------------------------!
! INPUT:  Transformation between a rotated system with north pole (lap,thp) and     !
!         the unrotated system (lau,thu)                                            !
! OUTPUT: coordinate system (lar,thr) of the rotated system,                        !
!         lar\in [-pi,pi]; thr\in [-pi/2,pi/2]                                      !
! Literature: *W. T. M. Verkley. The construction of barotropic modons on a         !
!              sphere. J. Atmos. Sci., 41, 2492-2504, 1984                          !
!             *R. D. Nair and C. Jablonowski. Moving Vortices on the Sphere: A Test !
!              Case for Horizontal Advection Problems. Mon. Wea. Rev, 136, 699-711, !
!              2008                                                                 !  
! Remark: if thr=pi/2, lar can be arbitrary (is pole, singularity)                  !
!-----------------------------------------------------------------------------------!

subroutine vortex_rotatedsphere(lap,thp,lau,thu,lar,thr)
  ! Rotate to new North Pole at (lap,thp)
  ! (lar,thr) are the rotated coordinates coorsponding to (lau,thu) in 
  ! the regular sphere
  use physical_constants, only : DD_PI
  
  implicit none
  real (kind=real_kind), intent(in)  :: lap,thp,lau,thu
  real (kind=real_kind), intent(out) :: lar,thr

  real (kind=real_kind) :: cost,sint,sinp,cosp
  real (kind=real_kind) ::  trm, trm1,trm2,trm3
  real (kind=real_kind) ::  pi2

  pi2=2.0D0*DD_PI

  sinp = sin(thp)
  cosp = cos(thp)
  cost = cos(thu)
  sint = sin(thu)

  trm  = cost * cos(lau- lap)
  trm1 = cost * sin(lau- lap)
  trm2 = sinp * trm  - cosp * sint
  trm3 = sinp * sint + cosp * trm

  lar = atan2(trm1,trm2)
  if (lar < 0.0D0 ) lar = lar + pi2
  if (lar > pi2)    lar = lar - pi2
  thr = asin(trm3)      

end subroutine vortex_rotatedsphere


!-----------------------------------------------------------------------------------!
! INPUT:  (BACK) Transformation between a rotated system with north pole (lap,thp)  !
!         with the coordinate system (lar, thr) to the unrotated system             !
! OUTPUT: coordinate system (lau,thu) of the unrotated system,                      !
!         rla\in [-pi,pi]; rth\in [-pi/2,pi/2]                                      !
! Literature: *W. T. M. Verkley. The construction of barotropic modons on a         !
!              sphere. J. Atmos. Sci., 41, 2492-2504, 1984                          !
!             *R. D. Nair and C. Jablonowski. Moving Vortices on the Sphere: A Test !
!              Case for Horizontal Advection Problems. Mon. Wea. Rev, 136, 699-711, !
!              2008                                                                 !
! Remark: if thr=pi/2, lar can be arbitrary (is pole, singularity)                  !
!-----------------------------------------------------------------------------------!

subroutine vortex_rotatedsphereback(lap,thp,lar,thr,lau,thu)
  use physical_constants, only : DD_PI

  implicit none
  real (kind=real_kind), intent(in)  :: lap,thp,lar,thr
  real (kind=real_kind), intent(out) :: lau,thu
  !
  real (kind=real_kind) :: cost,sint,cosp,sinp,clam,slam 
  real (kind=real_kind) :: trm, t1,t2,t3
  real (kind=real_kind) ::  pi2

  pi2=2.0D0*DD_PI

  !
  !* Back to unrotated system

  cost = cos(thr)
  sint = sin(thr)
  clam = cos(lar)
  slam = sin(lar)
  cosp = cos(thp)
  sinp = sin(thp)

  t1 = slam * cost
  t2 = sint*cosp + cost*clam*sinp
  t3 = sint*sinp - cost*clam*cosp
  lau =  lap+ atan2(t1,t2)
  if (lau < 0.0D0 )  lau = lau + pi2
  if (lau > pi2)     lau = lau - pi2
  thu =  asin(t3)

end subroutine vortex_rotatedsphereback

! subroutine filter(cf,minmax)
! 
!   implicit none
!   real (kind=real_kind), intent(in), dimension(nip,nip) :: cf
!   real (kind=real_kind), intent(in)   :: minmax(2), val
!   
! 
! end subroutine filter
! ----------------------------------------------------------------------------------!
!SUBROUTINE SLOPELIMITER_VAL----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
! DESCRIPTION: returns the value from the reconstruction (3rd order Taylor polynom) !
!              at the point (cartx,carty) -> in cube CARTESIAN coordinates          !
!                                                                                   !
! INPUT: value  ...  point value (calculated here by recons_val_cart)               !
!        cell_value ...  tracer value (in the cell center) of the cell              !
!        local_min ...  minmal value in the patch                                   !
!        local_max ...  maximal value in the patch                                  !
! INPUT/OUTPUT: min_phi ... slope limiter, inout because we go through any possible !
!                           extrema on the cell                                     !
!-----------------------------------------------------------------------------------!
subroutine slopelimiter_val(value, cell_value, local_min, local_max, min_phi)

  implicit none
  real (kind=real_kind), intent(in)    :: value, cell_value
  real (kind=real_kind), intent(in)    :: local_min, local_max
  real (kind=real_kind), intent(inout) :: min_phi

  real (kind=real_kind) :: phi 
  
  phi= 0.0D0
  ! Check against the minimum bound on the reconstruction
  if (value - cell_value > 1.0D-12 * value) then
    phi = (local_max - cell_value) / (value - cell_value)
    if (phi < min_phi) then
      min_phi = phi
    endif
  ! Check against the maximum bound on the reconstruction
  elseif (value - cell_value < -1.0D-12 * value) then
    phi = (local_min - cell_value) / (value - cell_value)    
    if(phi < min_phi) then
      min_phi = phi
    endif
  endif

end subroutine slopelimiter_val


end module spelt_mod
