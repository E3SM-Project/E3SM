!-----------------------------------------------------------------------------------!
!MODULE CWFV_MOD----------------------------------------------------------CE-for FVM!
! CWFV_MOD File for the fvm project in HOMME                                         !
! Author: Christoph Erath                                                           !
! Date: 25.January 2011                                                             !
! MAIN module to run fvm on HOMME                                                   !
! 14.November 2011: reorganisation done                                             !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module cwfv_mod

  use kinds, only : real_kind, int_kind
  use dimensions_mod, only: ne, nlev, ntrac, np, ntrac_d, nc, nhe, nip, nipm, nep
  use edge_mod, only : ghostBuffertr_t,edgebuffer_t
  use time_mod, only : timelevel_t
  use coordinate_systems_mod, only : spherical_polar_t, cartesian2D_t
  use element_mod, only : element_t, timelevels
  use hybrid_mod, only : hybrid_t

  implicit none
  private
  type (ghostBuffertr_t)                      :: cellghostbuf
  type (EdgeBuffer_t)                         :: edgeveloc
    
  type, public :: cwfv_struct
    ! fvm tracer mixing ratio: (kg/kg)
    real (kind=real_kind)    :: c(1-nipm:nep+nipm,1-nipm:nep+nipm,nlev,ntrac_d,timelevels) 
!-----------------------------------------------------------------------------------!   
    real (kind=real_kind)    :: vn0(np,np,2,nlev) 
    real (kind=real_kind)    :: contrau(1:nep,1:nep), contrav(1:nep,1:nep)
!-----------------------------------------------------------------------------------!
    ! define the departure grid (depends on the examples, usually the velocity)
    ! they are NOT stored for the Levels
    type (spherical_polar_t) :: asphere(nep,nep)     ! Spherical coordinates
    ! equidistant mesh spacing in each direction of the reference(!) element, depends on nep       
    real(kind=real_kind)     :: dx  
    real(kind=real_kind)     :: sga(1:nep,1:nep)   
    ! next maybe used just for test reasons
    real (kind=real_kind)    :: elem_mass
!       real (kind=real_kind)    :: area_sphere(1-nhe:nc+nhe,1-nhe:nc+nhe) 
    real (kind=real_kind) :: cstart(1:nc,1:nc)    ! for error analysis
!       real (kind=real_kind)    :: maxcfl(2,nlev)   
    !should be replace once the mapping comes from the reference element
    real (kind=real_kind)    :: Dinv(2,2,np,np)     ! Map vector field on the sphere to covariant v on cube

  end type cwfv_struct
  
  public :: cellghostbuf, edgeveloc, fvm_init1,fvm_init2, fvm_mcgregordsscwfv, spelt_run, fvm_grid_init
  public :: cip_coeff, cip_interpolate, metric_term, metric_termref, cell_search, qmsl_cell_filter, cell_minmax, cip_cell_avr
contains

subroutine spelt_run(elem,fvm,hybrid,deriv,tstep,tl,nets,nete)

  use derivative_mod, only : derivative_t
  ! ---------------------------------------------------------------------------------
  use edge_mod, only :  ghostVpack2d, ghostVunpack2d
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
  ! ---------------------------------------------------------------------------------
  use coordinate_systems_mod, only : spherical_to_cart, cart2cubedspherexy, sphere2cubedsphere
  ! ------EXTERNAL----------------
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------  
  
  implicit none
  type (element_t), intent(inout)             :: elem(:)
  type (cwfv_struct), intent(inout)           :: fvm(:)
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
  integer                                     :: icell1(1:nep,1:nep), jcell1(1:nep,1:nep)     
  integer                                     :: icell2(1:nep,1:nep), jcell2(1:nep,1:nep)
  
  integer                                     :: icell, jcell
  real (kind=real_kind)                       :: dxoy, dxyi, dt6, sg, sga
  type (cartesian2D_t)                        :: alphabeta
  real (kind=real_kind)                       :: tmp
  

  
  dt6  = tstep/ 6.0D0
  do ie=nets,nete 
    tmp=abs(elem(ie)%corners(1)%x-elem(ie)%corners(2)%x)/nc
    dxoy=tmp / 6.0D0
    dxyi=1.0D0/(tmp*tmp)
    do k=1, nlev
      call fvm_dep_from_gll(elem(ie), deriv, fvm(ie)%asphere,dsphere1,0.5D0*tstep,tl,k)         
      call fvm_dep_from_gll(elem(ie), deriv, fvm(ie)%asphere,dsphere2,tstep,tl,k)         
      !search has not to be done for all tracers!
      do j=1,nep
        do i=1,nep
!           call solidbody(fvm(ie)%asphere(i,j), dsphere, 0.5D0) 
          call cell_search(elem(ie), dsphere1(i,j),icell1(i,j), jcell1(i,j),dref1(i,j),alphabeta)
          sg1(i,j)=metric_term(alphabeta)
!           sg1(i,j)=metric_termref(elem(ie),dref1(i,j))
!           call solidbody(fvm(ie)%asphere(i,j), dsphere, 1.0D0)  
          call cell_search(elem(ie), dsphere2(i,j), icell2(i,j), jcell2(i,j),dref2(i,j),alphabeta)
          sg2(i,j)=metric_term(alphabeta)
!           sg2(i,j)=metric_termref(elem(ie),dref2(i,j))
          
        end do
      end do
      ! search of both point on the trajectory done
      do itr=1,ntrac
        do j=1-nhe,nc+nhe
          do i=1-nhe,nc+nhe
            icell=1+(i-1)*nipm
            jcell=1+(j-1)*nipm
            ff=fvm(ie)%c(icell:icell+nipm,jcell:jcell+nipm,k,itr,tl%n0)
            minmax(i,j,:)=cell_minmax(ff)
            call cip_coeff(ff,ff(2,2),cf(:,:,i,j))
          enddo
        enddo
        do j=1,nep
          do i=1,nep  
            sga=fvm(ie)%sga(i,j)
            slval(1)=fvm(ie)%c(i,j,k,itr,tl%n0)
 
            tmp=cip_interpolate(cf(:,:,icell1(i,j),jcell1(i,j)),dref1(i,j)%x,dref1(i,j)%y) 
            tmp=qmsl_cell_filter(icell1(i,j),jcell1(i,j),minmax,tmp)
            slval(2)=(sga/sg1(i,j))*tmp

            tmp=cip_interpolate(cf(:,:,icell2(i,j),jcell2(i,j)),dref2(i,j)%x,dref2(i,j)%y) 
            tmp=qmsl_cell_filter(icell2(i,j),jcell2(i,j),minmax,tmp)
            slval(3)=(sga/sg2(i,j))*tmp
            fvm(ie)%c(i,j,k,itr,tl%np1)=slval(3)

           if (mod(i,2)==1) then                   ! works only for nip=3!!!
              fluxval(i,j,1) =  dt6 * fvm(ie)%contrau(i,j)* (slval(1) + & 
                     4.0D0 * slval(2) + slval(3) )
            endif
            if (mod(j,2)==1) then            ! works only for nip=3!!!
              fluxval(i,j,2) =  dt6 * fvm(ie)%contrav(i,j)* (slval(1) + & 
                     4.0D0 * slval(2) + slval(3) )    
            endif
          end do
        end do 
        
        do j=1,nep,2
          do i=1,nep,2            
              flux(1) = dxoy * (fluxval(i,j,1) + 4.0D0 * fluxval(i,j+1,1) + fluxval(i,j+2,1))  ! west
              flux(2) = dxoy * (fluxval(i,j,2) + 4.0D0 * fluxval(i+1,j,2) + fluxval(i+2,j,2))  ! south
              flux(3) = dxoy * (fluxval(i+2,j,1) + 4.0D0 * fluxval(i+2,j+1,1) + fluxval(i+2,j+2,1)) ! east
              flux(4) = dxoy * (fluxval(i+2,j+2,2) + 4.0D0 * fluxval(i+1,j+2,2) + fluxval(i,j+2,2)) ! north
              fvm(ie)%c(i+1,j+1,k,itr,tl%np1) = fvm(ie)%c(i+1,j+1,k,itr,tl%n0) + &
                                        (flux(1) + flux(2) - flux(3) - flux(4) ) * dxyi
          end do
        end do 
      end do
    end do
    call ghostVpack2d(cellghostbuf,fvm(ie)%c,nipm, nep,nlev,ntrac,0, tl%np1, timelevels,elem(ie)%desc)
  end do
!-----------------------------------------------------------------------------------! 
  call t_startf('FVM Communication') 
  call ghost_exchangeV(hybrid,cellghostbuf,nipm,nep,ntrac)
  call t_stopf('FVM Communication')
!-----------------------------------------------------------------------------------!  
  call t_startf('FVM Unpacking')  
  do ie=nets,nete
    call ghostVunpack2d(cellghostbuf,fvm(ie)%c,nipm, nep,nlev,ntrac,0, tl%np1, timelevels,elem(ie)%desc)
  end do
  call t_stopf('FVM Unpacking')
end subroutine spelt_run

! initialize global buffers shared by all threads
subroutine fvm_init1(par)
  use edge_mod, only : initghostbuffer,initEdgebuffer
  use parallel_mod, only : parallel_t, haltmp
  type (parallel_t) :: par
  
  if (ntrac>ntrac_d) then
     if (par%masterproc) print *,'ntrac,ntrac_d=',ntrac,ntrac_d
     call haltmp("PARAMTER ERROR for fvm: ntrac > ntrac_d")
  endif

  call initghostbuffer(cellghostbuf,nlev,ntrac,nipm,nep) !+1 for the air_density, which comes from SE
  call initEdgebuffer(edgeveloc,2*nlev)
end subroutine fvm_init1

! initialization that can be done in threaded regions
subroutine fvm_init2(elem,fvm,hybrid,nets,nete,tl)
  use bndry_mod, only: compute_ghost_corner_orientation 

  type (timelevel_t)                    :: tl
  type (cwfv_struct)                    :: fvm(:)
  type (element_t)                      :: elem(:)
  type (hybrid_t)                       :: hybrid
  integer                               :: ie,nets,nete
  
  call compute_ghost_corner_orientation(hybrid,elem,nets,nete)
  call fvm_grid_init(elem,fvm, nets, nete, tl)
end subroutine fvm_init2


subroutine fvm_grid_init(elem,fvm,nets,nete,tl)
  use coordinate_systems_mod, only : sphere2cubedsphere
  use cube_mod, only : vmap, ref2sphere
  implicit none
  type (element_t),intent(inout)            :: elem(:)                 
  type (cwfv_struct),intent(inout)          :: fvm(:)
  integer, intent(in)                       :: nets, nete
  type (TimeLevel_t), intent(in)            :: tl              ! time level struct
  
  integer                                   :: ie, i, j
  real (kind=real_kind)                     :: xref, yref, dx, tmpD(2,2), detD
  type (cartesian2D_t)                      :: alphabeta
  type (cartesian2D_t)                      :: dref
  
  do ie=nets,nete
    dx=2.0D0/(nep-1)   ! equi-distant grid on reference element in both directions!
    fvm(ie)%dx=dx
    do j=1,nep
      yref=-1+(j-1)*dx
      do i=1,nep  
        xref=-1+(i-1)*dx
        !define the arrival grid in spherical coordinates
        fvm(ie)%asphere(i,j)=ref2sphere(xref,yref,elem(ie)%corners,elem(ie)%FaceNum)             
        alphabeta=sphere2cubedsphere(fvm(ie)%asphere(i,j), elem(ie)%FaceNum)
        fvm(ie)%sga(i,j)=metric_term(alphabeta)
        
        dref%x=xref
        dref%y=yref
!         fvm(ie)%sga(i,j)=metric_termref(elem(ie),dref)
        
      end do 
    end do
    !this can be deleted, once the transformation is on the reference element
    do j=1,np
      do i=1,np
        call vmap(tmpD,elem(ie)%cartp(i,j)%x,elem(ie)%cartp(i,j)%y,elem(ie)%FaceNum)
        detD = tmpD(1,1)*tmpD(2,2) - tmpD(1,2)*tmpD(2,1)      

        fvm(ie)%Dinv(1,1,i,j) =  tmpD(2,2)/detD
        fvm(ie)%Dinv(1,2,i,j) = -tmpD(1,2)/detD
        fvm(ie)%Dinv(2,1,i,j) = -tmpD(2,1)/detD
        fvm(ie)%Dinv(2,2,i,j) =  tmpD(1,1)/detD
      enddo
    enddo
  end do
end subroutine fvm_grid_init

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
subroutine fvm_mcgregordsscwfv(elem,fvm,nets,nete, hybrid, deriv, tstep, ordertaylor)
  use derivative_mod, only : derivative_t, ugradv_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  
  implicit none

  type (element_t), intent(inout)             :: elem(:)
  type (cwfv_struct), intent(in)              :: fvm(:)

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
        vhat(ie,:,:,:,:)=(fvm(ie)%vn0(:,:,:,:) + ugradv(ie,:,:,:,:))/2 
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

end subroutine fvm_mcgregordsscwfv
!END SUBROUTINE FVM_MCGREGORDSS-------------------------------------------CE-for FVM!

subroutine cell_search(elem, dsphere, icell, jcell,dref, alphabeta) 
  use interpolate_mod, only : cube_facepoint_ne 
  use coordinate_systems_mod, only : cube_face_number_from_sphere 
  use physical_constants, only : DD_PI

  implicit none

  type (spherical_polar_t), intent(in)     :: dsphere

  type (element_t), intent(in)             :: elem
  integer, intent(out)                     :: icell, jcell
  type (cartesian2D_t),intent(out)         :: dref
  type (cartesian2D_t), intent(out)        :: alphabeta
    
  type (cartesian2D_t)                     :: dcart
  integer                                  :: number, endi, starti, tmpi
  real (kind=real_kind)                    :: refnc(1:nc+1), dxcell, tmp

  integer                                  :: i, face_nodep      

  tmp=nc
  do i=1,nc+1
    refnc(i)= 2*(i-1)/tmp - 1
  end do
  dxcell=abs(refnc(2)-refnc(1))
  
!   dxcell=abs(elem%corners(1)%x-elem%corners(2)%x)/(nc)
!   write(*,*) dxcell
  call cube_facepoint_ne(dsphere,ne,dcart, alphabeta, number, face_nodep)
    ! Search index along "x"  (bisection method)
  starti = 1
  endi = nc+1
  do
    if  ((endi-starti) <=  1)  exit
    tmpi = (endi + starti)/2
    if (dcart%x  >  refnc(tmpi)) then
      starti = tmpi
    else
      endi = tmpi
    endif
  enddo
  icell = starti
  ! Search index along "y"
  starti = 1
  endi = nc+1
  do
    if  ((endi-starti) <=  1)  exit
    tmpi = (endi + starti)/2
    if (dcart%y  >  refnc(tmpi)) then
      starti = tmpi
    else
      endi = tmpi
    endif
  enddo
  jcell = starti
    
  dref%x=dcart%x-refnc(icell)
  dref%y=dcart%y-refnc(jcell)
!     
  if ((dref%x<-1.0D-12) .or.(dref%y<-1.0D-12) .or.(dref%x>dxcell+1.0D-12) .or. (dref%y>dxcell+1.0D-12) ) then
    write(*,*) 'Something is wrong in search!'
    write(*,*) dref%x, dref%y
    stop
  endif
    
  if ((icell<1) .or.(icell>nc) .or. (jcell<1) .or. (jcell>nc)) then
    write(*,*) 'icell, jcell,Something is wrong in search!'
    stop
  endif

  if(number==elem%vertex%nbrs(1)%n) then  !west
    if ((elem%FaceNum<=4) .or. (face_nodep==elem%FaceNum)) then
      icell=icell-nc
    elseif (elem%FaceNum==5) then
      tmpi=icell
      icell=-jcell+1
      jcell=tmpi
      tmp=dref%x
      dref%x=dxcell-dref%y
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
         
  if(number==elem%vertex%nbrs(2)%n) then   !east
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
      dref%x=dxcell-dref%y
      dref%y=tmp
    endif
  end if   
         
  if(number==elem%vertex%nbrs(3)%n) then   !south
    if ((elem%FaceNum==1) .OR. (elem%FaceNum==6) .or. (face_nodep==elem%FaceNum)) then
      jcell=jcell-nc
    elseif (elem%FaceNum==2) then
      tmpi=icell
      icell=nc+1-jcell
      jcell=tmpi-nc
      tmp=dref%x
      dref%x=dxcell-dref%y
      dref%y=tmp
    elseif ((elem%FaceNum==3) .or. (elem%FaceNum==5)) then
      icell=nc+1-icell
      jcell=-jcell+1
      dref%x=dxcell-dref%x
      dref%y=dxcell-dref%y
    elseif (elem%FaceNum==4) then
      tmpi=icell
      icell=jcell
      jcell=1-tmpi
      tmp=dref%x
      dref%x=dref%y
      dref%y=dxcell-tmp
    endif
  end if   
         
  if(number==elem%vertex%nbrs(4)%n) then   !north
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
      dref%y=dxcell-dref%y
    elseif (elem%FaceNum==4) then
      tmpi=icell
      icell=1+(nc-jcell)
      jcell=nc+tmpi
      tmp=dref%x
      dref%x=dxcell-dref%y
      dref%y=tmp
    endif
  end if   
  ! cases southwest, southeast, northwest, northeast on a cube corner do not exist
  if(number==elem%vertex%nbrs(5)%n) then   !southwest
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
        dref%x=dxcell-dref%y
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
        dref%y=dxcell-dref%y
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
        dref%x=dxcell-dref%y
        dref%y=tmp
      else
        icell=1-icell
        jcell=1-jcell
        dref%x=dxcell-dref%x
        dref%y=dxcell-dref%y
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
       
  if(number==elem%vertex%nbrs(6)%n) then   !southeast
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
        dref%x=dxcell-dref%y
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
        dref%y=dxcell-dref%y
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
        dref%y=dxcell-dref%y
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
        dref%x=dxcell-dref%y
        dref%y=tmp
      endif
    endif
  endif   

  if(number==elem%vertex%nbrs(7)%n) then   !northwest
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
        dref%y=dxcell-dref%y
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
        dref%x=dxcell-dref%y
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
        dref%x=dxcell-dref%y
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
        dref%y=dxcell-dref%y
      endif
    endif
  end if   
       
  if(number==elem%vertex%nbrs(8)%n) then   !northeast
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
        dref%y=dxcell-dref%y
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
        dref%x=dxcell-dref%y
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
        dref%x=dxcell-dref%y
        dref%y=tmp
      else
        icell=nc+(nc+1-icell)
        jcell=nc+(nc+1-jcell)
        dref%x=dxcell-dref%x
        dref%y=dxcell-dref%y
      endif
    endif
  end if
      
  if ((icell<0) .or.(icell>nc+1) .or. (jcell<0) .or. (jcell>nc+1)) then
    write(*,*) '2 Something is wrong in search!'
    write(*,*) number, elem%vertex%nbrs(1)%n, elem%vertex%nbrs(2)%n, elem%vertex%nbrs(3)%n, elem%vertex%nbrs(4)%n, elem%vertex%nbrs(5)%n, elem%vertex%nbrs(6)%n, elem%vertex%nbrs(7)%n, elem%vertex%nbrs(8)%n
    write(*,*) icell, jcell, elem%GlobalId, elem%FaceNum, face_nodep
    stop
  endif
  if ((dref%x<-1.0D-12) .or.(dref%y<-1.0D-12) .or.(dref%x>dxcell+1.0D-12) .or. (dref%y>dxcell+1.0D-12) ) then
    write(*,*) '3 Something is wrong in search!'
    write(*,*) number, elem%vertex%nbrs(1)%n, elem%vertex%nbrs(2)%n, elem%vertex%nbrs(3)%n, elem%vertex%nbrs(4)%n, elem%vertex%nbrs(5)%n, elem%vertex%nbrs(6)%n, elem%vertex%nbrs(7)%n, elem%vertex%nbrs(8)%n
    write(*,*) dref
    stop
  endif
end subroutine cell_search


! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_DEP_FROM_GLL----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR 14. December 2011                            !
! DESCRIPTION: calculates the deparute grid for fvm coming from the gll points      !
!                                                                                   !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine fvm_dep_from_gll(elem, deriv, asphere,dsphere,dt,tl,klev)
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
  ! interpolate velocity to fvm nodes
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
end subroutine fvm_dep_from_gll
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

function metric_termref(elem,ref) result(sg) 
  use cube_mod, only : dmap
  implicit none
  type (element_t), intent(in)          :: elem
  
  type (cartesian2D_t), intent(in)      :: ref
  real (kind=real_kind)                 :: sg                

  real (kind=real_kind)  :: Jp(2,2)
  real (kind=real_kind)  :: D(2,2)


  ! input (a,b) shold be a point in the reference element [-1,1]
  ! compute Jp(a,b)
  Jp(1,1) = elem%u2qmap(2,1) + elem%u2qmap(4,1)*ref%y
  Jp(1,2) = elem%u2qmap(3,1) + elem%u2qmap(4,1)*ref%x
  Jp(2,1) = elem%u2qmap(2,2) + elem%u2qmap(4,2)*ref%y
  Jp(2,2) = elem%u2qmap(3,2) + elem%u2qmap(4,2)*ref%x
  
  call Dmap(D, elem, ref%x,ref%y)
  
  
  sg=(abs(D(1,2)*D(2,1)-D(1,1)*D(2,2)))
end function metric_termref

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

! Positivity preserving QMSL filter for given point 
function qmsl_cell_filter(ix,iy,minmax,val) result(qval) 
  implicit none
  integer, intent(in)                 :: ix,iy
  real (kind=real_kind), intent(in)   :: minmax(1-nhe:nc+nhe,1-nhe:nc+nhe,2), val
  real (kind=real_kind)               :: qval

  qval=val
  if (val <= minmax(ix,iy,1)) then 
    qval = minmax(ix,iy,1) 
  endif 
  if (val >= minmax(ix,iy,2)) then 
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


subroutine cip_coeff(ff,avr,cf) 
  use physical_constants, only : DD_PI

  implicit none
  real (kind=real_kind) , intent(in), dimension(nip,nip) :: ff
  real (kind=real_kind) , intent(in) :: avr 
  real (kind=real_kind) , intent(out) :: cf(nip,nip)  
  real (kind=real_kind)  :: dxcell, dycell, d(8) 
  integer ::  i,j

  dxcell = 1.0D0/(2.0D0/nc)  
!    dxcell=1.0D0/(DD_PI/2/nc)
  dycell = dxcell               
  d(1) = dxcell
  d(2) = 2.0D0*dxcell*dxcell 
  d(3) = d(1) 
  d(4) = d(2)          
  d(5) = d(2)          
  d(6) = dxcell*dycell*dycell   
  d(7) = d(6)          
  d(8) = 3.0D0*d(6)*dxcell  

  cf(1,1) = ff(1,1)
  cf(2,1) = d(1) * (-3.0D0*ff(1,1) + 4.0D0*ff(2,1) - ff(3,1) ) 
  cf(3,1) = d(2) * (       ff(1,1) - 2.0D0*ff(2,1) + ff(3,1) ) 

  cf(1,2) = d(3) * (-3.0D0*ff(1,1) + 4.0D0*ff(1,2) - ff(1,3) ) 
  cf(1,3) = d(4) * (       ff(1,1) - 2.0D0*ff(1,2) + ff(1,3) ) 

  cf(2,2) = d(5) * ( 4.0D0*(ff(1,1) - ff(2,3) - ff(3,2))  + ff(1,3) &
                     + ff(3,1) - 8.0D0*(ff(2,1) + ff(1,2))  + 18.0D0*avr ) 

  cf(3,2) = d(6) * (-5.0D0*(ff(1,1) + ff(3,1)) + 12.0D0*(ff(1,2) + ff(3,2))  - ff(1,3) &
                     + 8.0D0*ff(2,3)  + 16.0D0*ff(2,1) - ff(3,3) - 36.0D0*avr ) 

  cf(2,3) = d(7) * (-5.0D0*(ff(1,1) + ff(1,3)) + 12.0D0*(ff(2,1)+ ff(2,3)) &
                     + 16.0D0*ff(1,2) - ff(3,1) + 8.0D0*ff(3,2) - ff(3,3) - 36.0D0*avr ) 

  cf(3,3) = d(8) * (-4.0D0*(ff(1,2) + ff(2,1) + ff(2,3) + ff(3,2))   &
                     + ff(1,1) + ff(1,3) + ff(3,1) + ff(3,3) + 12.0D0*avr ) 
end  subroutine cip_coeff

function cip_cell_avr(ff)  result(qint)
  implicit none
  real (kind=real_kind), intent(in), dimension(nip,nip) :: ff
  real (kind=real_kind) :: qint
  real (kind=real_kind) :: f(3)       
  real (kind=real_kind) :: quart, delta
  integer ::   i,j,k,l, ie,je ,ip

!  delta = wel(1) * wel(1) * 0.25D0
  qint = 0.0D0

! cell-wise integral  by Simpson's 3pt rule 
  do j = 1, 3  
    f(j) = ff(1,j) + 4.0D0*ff(2,j) + ff(3,j) 
  enddo 
  qint = (f(1) + 4.0D0 * f(2) + f(3))/36.0D0 
end function cip_cell_avr

end module cwfv_mod
