!-----------------------------------------------------------------------------------!
! MODULE FVM_MOD----------------------------------------------------------CE-for FVM!
! fvm_MOD File for the fvm project in HOMME                                         !
! Author: Christoph Erath                                                           !
! Date: 25.January 2011                                                             !
! MAIN module to run fvm on HOMME                                                   !
! 14.November 2011: reorganisation done                                             !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module fvm_bench_mod
  use kinds, only : real_kind, int_kind, longdouble_kind
  use edge_mod, only : freeghostbuffertr, ghostVpack, ghostVunpack, &
                       edgeVpack, edgeVunpack, freeedgebuffer 
  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, np, ne, nc, nhc, nhe, nlev, ntrac
  use time_mod, only : timelevel_t
  use element_mod, only : element_t, timelevels
  use hybrid_mod, only : hybrid_t

contains


subroutine cslam_run_bench(elem,fvm,red,hybrid,nets,nete,tl)
  ! ---------------------------------------------------------------------------------  
  use fvm_bsp_mod, only: fvm_bsp, get_boomerang_velocities_gll, get_solidbody_velocities_gll
  ! ---------------------------------------------------------------------------------  
  use fvm_control_volume_mod, only: fvm_struct
  ! ---------------------------------------------------------------------------------
  use fvm_mod, only: cslam_runairdensity, fvm_init1,fvm_init2, fvm_init3, fvm_mcgregor,fvm_mcgregordss, cellghostbuf, edgeveloc
  ! ---------------------------------------------------------------------------------
  use fvm_line_integrals_mod, only: compute_weights
  ! ---------------------------------------------------------------------------------  
  use fvm_filter_mod, only: monotonic_gradient_cart
  ! ---------------------------------------------------------------------------------
  use checksum_mod, only: test_ghost
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t, derivative_stag_t, derivinit, deriv_print
  ! ---------------------------------------------------------------------------------
  use reduction_mod, only : ReductionBuffer_ordered_1d_t
  ! ---------------------------------------------------------------------------------
  use thread_mod, only : nthreads
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV, bndry_exchangeV
  ! ---------------------------------------------------------------------------------
  use time_mod, only : tstep, nmax, time_at, timelevel_update, timelevel_init
  ! ---------------------------------------------------------------------------------
  use thread_mod, only : NThreads
  ! ---------------------------------------------------------------------------------
  use coordinate_systems_mod, only : spherical_polar_t,spherical_to_cart, &
                                     cart2cubedspherexy
  ! ---------------------------------------------------------------------------------
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  ! ---------------------------------------------------------------------------------
  use parallel_mod, only: global_shared_buf, global_shared_sum
  ! ---------------------------------------------------------------------------------
  use global_norms_mod, only: wrap_repro_sum
  ! ---------------------------------------------------------------------------------
  use reduction_mod, only : parallelmax, parallelmin
  ! ---------------------------------------------------------------------------------  
  use physical_constants, only : g
  ! ------EXTERNAL----------------
  use perf_mod, only : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  ! -----------------------------------------------
  
#ifdef PIO_INTERP
    use interp_movie_mod, only : interp_movie_init, interp_movie_output, interp_movie_finish
#else
    use shal_movie_mod, only : shal_movie_init, shal_movie_output, shal_movie_finish
#endif
  
  implicit none
  type (element_t), intent(inout)                :: elem(:)
  type (fvm_struct), intent(inout)             :: fvm(:)
  type (ReductionBuffer_ordered_1d_t),intent(in)    :: red   ! reduction buffer         (shared)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)

  real (kind=real_kind)                       :: massstart, mass, maxc, maxcstart,minc, mincstart, tmp, tmpref  
  real (kind=real_kind)                       :: tmp1(nets:nete), tmp2(nets:nete)
  real (kind=real_kind)                       :: l1,l2, lmax 
  
  integer                                     :: i,j,k,ie,itr, jx, jy, jdx, jdy, h
  type (TimeLevel_t)                          :: tl              ! time level struct
  type (derivative_t)                         :: deriv           ! derivative struct
  type (spherical_polar_t)                    :: tmpsphereincart   
 
  character (len=99)                          :: filename
  
  real (kind=real_kind)   , dimension(10*(nc+2*nhe)*(nc+2*nhe),6)  :: weights_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_eul_index_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_lgr_index_all
  integer (kind=int_kind)                                          :: jall
    
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons  
  real (kind=real_kind), dimension(nelemd,1-nhe:nc+nhe,1-nhe:nc+nhe) :: area    
  real (kind=real_kind)                                              :: xtmp
  real (kind=longdouble_kind)                                        :: fvm_nodes(nc+1)
  
  real (kind=real_kind), dimension(np,np,2)    :: vstar, vhat
  real (kind=real_kind)                        :: maxcflx, maxcfly  
  
  
  integer  choosetrac, chooselev   !for test reason the output
 !-----------------------------------------------------------------------------------!  
 choosetrac=4
 chooselev=1
 
  if(hybrid%masterthread) then 
    print *,"!-----------------------------------------------------------------------!"
    print *,"!  Test CASE for fvm, Christoph Erath                                 !" 
    print *,"!-----------------------------------------------------------------------!" 
  endif
     
  ! Initialize derivative structure
  ! fvm nodes are equally spaced in alpha/beta
  ! HOMME with equ-angular gnomonic projection maps alpha/beta space
  ! to the reference element via simple scale + translation
  ! thus, fvm nodes in reference element [-1,1] are a tensor product of
  ! array 'fvm_nodes(:)' computed below:
  xtmp=nc 
  do i=1,nc+1
    fvm_nodes(i)= 2*(i-1)/xtmp - 1
  end do
  call derivinit(deriv,fvm_corners=fvm_nodes)

!-----------------------------------------------------------------------------------!    
  do ie=nets,nete
    call fvm_bsp(fvm(ie),tl)
    fvm(ie)%elem_mass=0
    do j=1,nc
      do i=1,nc
        if (choosetrac==1) then   ! mass of air, code is not optimal
          fvm(ie)%elem_mass=fvm(ie)%elem_mass + &
                        fvm(ie)%area_sphere(i,j)*fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)
        else
          fvm(ie)%elem_mass=fvm(ie)%elem_mass + &
                        fvm(ie)%area_sphere(i,j)*fvm(ie)%c(i,j,chooselev,1,tl%n0)*&
                                                   fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)
        endif
        fvm(ie)%cstart(i,j)=fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)
      enddo
    enddo
    !

!     call ghostVpack(cellghostbuf, fvm(ie)%c,nhc,nc,nlev,ntrac,0,tl%n0,timelevels,elem(ie)%desc)
    ! reset the new unknown
    do k=1,nlev
      do itr=1,ntrac
        do j=1-nhc,nc+nhc
          do i=1-nhc,nc+nhc 
          fvm(ie)%c(i,j,k,itr,tl%np1)=0.0D0
          end do
        end do
      enddo
    enddo
  end do
  
  !first exchange of the initial values
  call fvm_init3(elem,fvm,hybrid,nets,nete,tl%n0)
!-----------------------------------------------------------------------------------!  
!   call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc)
!-----------------------------------------------------------------------------------!    

  do ie=nets,nete
!      call ghostVunpack(cellghostbuf, fvm(ie)%c, nhc, nc,nlev,ntrac, 0, tl%n0, timelevels,elem(ie)%desc)
    ! for the mass value
    global_shared_buf(ie,1)=0D0
    global_shared_buf(ie,1)=fvm(ie)%elem_mass
    ! for the max value on the sphere
    tmp1(ie) = MAXVAL(fvm(ie)%c(:,:,chooselev,choosetrac,tl%n0))
    tmp2(ie) = MINVAL(fvm(ie)%c(:,:,chooselev,choosetrac,tl%n0))   
  end do


!-----------------------------------------------------------------------------------!
  
  !need the buffer cellghostbuf in the time loop
  ! for mass calculation
  call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
  massstart=global_shared_sum(1)
  maxcstart = parallelmax(tmp1,hybrid)
  mincstart = parallelmin(tmp2,hybrid)
!-----------------------------------------------------------------------------------!  
!Initialize Output via geopotential (should be changed, separate output for fvm
!write first time step to IO 
#ifdef PIO_INTERP
  call interp_movie_init(elem,hybrid,nets,nete,tl=tl)    
  call interp_movie_output(elem,tl, hybrid, 0D0, deriv, nets, nete,fvm)
#else
    call shal_movie_init(elem,hybrid,fvm)
    call shal_movie_output(elem,tl, hybrid, 0D0, nets, nete,deriv,fvm)
#endif 
!-----------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------!
  if(hybrid%masterthread) then 
    print *
    print *,"Arrival grid created , interpolation points calculated, initialization done. " 
    print *
  endif
  tmp=0
  
  call t_barrierf('fvm time loop', hybrid%par%comm)
  call t_startf('fvm')
  
  !BEGIN TIME LOOP, start at 0, calculate then next step
  DO WHILE(tl%nstep<nmax)
! start old mcgregor----------------------
!     do ie=nets,nete
!       do k=1,nlev
!         vstar = get_boomerang_velocities_gll(elem(ie), time_at(tl%nstep+1))
!         vhat= (get_boomerang_velocities_gll(elem(ie), time_at(tl%nstep)) + vstar) / 2.0D0
!         ! calculate high order approximation
!         call fvm_mcgregor(elem(ie), deriv, tstep, vhat, vstar, 3)
!      
!         ! apply DSS to make vstar C0
!         elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%spheremp(:,:)*vstar(:,:,1) 
!         elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%spheremp(:,:)*vstar(:,:,2) 
!       enddo
!       call edgeVpack(edgeveloc,elem(ie)%derived%vstar(:,:,1,:),nlev,0,elem(ie)%desc)
!       call edgeVpack(edgeveloc,elem(ie)%derived%vstar(:,:,2,:),nlev,nlev,elem(ie)%desc)
!     enddo 
!     call bndry_exchangeV(hybrid,edgeveloc)
!     do ie=nets,nete
!        call edgeVunpack(edgeveloc,elem(ie)%derived%vstar(:,:,1,:),nlev,0,elem(ie)%desc)
!        call edgeVunpack(edgeveloc,elem(ie)%derived%vstar(:,:,2,:),nlev,nlev,elem(ie)%desc)
!        do k=1, nlev  
!          elem(ie)%derived%vstar(:,:,1,k)=elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%rspheremp(:,:)
!          elem(ie)%derived%vstar(:,:,2,k)=elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%rspheremp(:,:)
!        end do
!     end do
!end old mcgegor-----------------------
! ! start mcgregordss
    do ie=nets,nete
      do k=1,nlev
        elem(ie)%derived%vstar(:,:,:,k)=get_boomerang_velocities_gll(elem(ie), time_at(tl%nstep+1))
        fvm(ie)%vn0(:,:,:,k)=get_boomerang_velocities_gll(elem(ie),time_at(tl%nstep))
!         elem(ie)%derived%vstar(:,:,:,k)=get_solidbody_velocities_gll(elem(ie), time_at(tl%nstep+1))
!         fvm(ie)%vn0(:,:,:,k)=get_solidbody_velocities_gll(elem(ie),time_at(tl%nstep))
      end do
    end do
    call fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, tstep, 3)
! ! end mcgregordss   
    call cslam_runairdensity(elem,fvm,hybrid,deriv,tstep,tl,nets,nete)
    
    call TimeLevel_update(tl,"forward")
     
if (mod(tl%nstep,50)==0) then  
    do ie=nets,nete
    ! prepare data for I/O
      global_shared_buf(ie,1)=0D0  ! for mass calculation
      ! test mass, just for chooselev and choosetrac, it is not optimized yet
      do j=1,nc
        do i=1,nc   
          if (choosetrac==1) then
            global_shared_buf(ie,1)=global_shared_buf(ie,1)+fvm(ie)%area_sphere(i,j)*&
                                    fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)
          else   
            global_shared_buf(ie,1)=global_shared_buf(ie,1)+fvm(ie)%area_sphere(i,j)*&
                 fvm(ie)%c(i,j,chooselev,1,tl%n0)*fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)
          endif
        end do
      end do
      ! for the max/min value on the sphere
      tmp1(ie) = MAXVAL(fvm(ie)%c(:,:,chooselev,choosetrac,tl%n0))
      tmp2(ie) = MINVAL(fvm(ie)%c(:,:,chooselev,choosetrac,tl%n0))
    end do
!-----------------------------------------------------------------------------------!
    ! for mass calculation
    call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
    mass=global_shared_sum(1)
    maxc = parallelmax(tmp1,hybrid)
    minc = parallelmin(tmp2,hybrid)
    maxcflx = parallelmax(fvm(:)%maxcfl(1,chooselev),hybrid)
    maxcfly = parallelmax(fvm(:)%maxcfl(2,chooselev),hybrid)
    !
    if  (hybrid%masterthread) then 
      write(*,*) 'time=', time_at(tl%nstep), 'timeatmax',Time_at(nmax)
      write(*,*) 'chooselev=', chooselev, 'choosetrac=', choosetrac
      write(*,*) 'STEP',tl%nstep,'MAXSTEP',nmax, &
                 't0', tl%n0, 't1', tl%np1
               write(*,*) 'massbegin', massstart, 'massend', mass 
      write(*,*) 'rel', (mass-massstart)/massstart           
      write(*,*) 'maxvaluestart:', maxcstart, 'minvaluestart:', mincstart
      write(*,*) 'maxvalue:     ', maxc,       'minvalue:    ', minc
      write(*,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly 
      print *
    endif
endif
!-----------------------------------------------------------------------------------!  

#ifdef PIO_INTERP
    call interp_movie_output(elem, tl, hybrid, 0D0, deriv, nets, nete,fvm)
#else     
    call shal_movie_output(elem, tl, hybrid, 0D0, nets, nete,deriv,fvm)
#endif
!-----------------------------------------------------------------------------------!  
  END DO
!------------END TIME LOOP-------------END TIME LOOP--------------END TIME LOOP-----!
!-----------------------------------------------------------------------------------! 
  call t_stopf('fvm')


  call freeghostbuffertr(cellghostbuf)
  call freeedgebuffer(edgeveloc)
#ifdef PIO_INTERP
    call interp_movie_finish
#else
    call shal_movie_finish
#endif
!-----------------------------------------------------------------------------------!  
! Error analysis/ complicated, but for a first try o.k.
    do ie=nets,nete
      tmp=0.0D0
      tmpref=0.0D0
      global_shared_buf(ie,:)=0.0D0
      do j=1,nc
        do i=1,nc
          global_shared_buf(ie,1)=global_shared_buf(ie,1)+fvm(ie)%area_sphere(i,j)*abs(fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)-fvm(ie)%cstart(i,j))
          global_shared_buf(ie,2)=global_shared_buf(ie,2)+fvm(ie)%area_sphere(i,j)*abs(fvm(ie)%cstart(i,j))
          
          global_shared_buf(ie,3)=global_shared_buf(ie,3)+fvm(ie)%area_sphere(i,j)*(fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)-fvm(ie)%cstart(i,j))* &
                                            (fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)-fvm(ie)%cstart(i,j))
          global_shared_buf(ie,4)=global_shared_buf(ie,4)+fvm(ie)%area_sphere(i,j)*(fvm(ie)%cstart(i,j))*(fvm(ie)%cstart(i,j))
          
          tmp=max(tmp,abs(fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)-fvm(ie)%cstart(i,j)))
          tmpref=max(tmpref,abs(fvm(ie)%cstart(i,j)))
        end do
      end do
      tmp1(ie)=tmp
      tmp2(ie)=tmpref
    end do
    call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
    l1=global_shared_sum(1)/global_shared_sum(2)
    l2=sqrt(global_shared_sum(3)/global_shared_sum(4))
    
    lmax = parallelmax(tmp1,hybrid)/parallelmax(tmp2,hybrid)


!SUMMARY
  if(hybrid%masterthread) then 
    print *
    print *,"!-----------------------------------------------------------------------!"
    print *,"!  Test CASE for FVM, Christoph Erath                                   !" 
    print *,"!-----------------------------------------------------------------------!"
    print *  
    write(*,*) 'number of elements', 6*ne*ne*nc*nc

    print *
    write(*,*) 'chooselev=', chooselev, 'choosetrac=', choosetrac
    write(*,*) 'massbegin', massstart, 'massend', mass 
    write(*,*) 'rel', (mass-massstart)/massstart
    write(*,*) 'maxvaluestart:', maxcstart, 'minvaluestart:', mincstart
    write(*,*) 'maxvalue:     ', maxc,      'minvalue:     ', minc
    write(*,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly 
    write(*,*) "l1 = ", l1, "l2 = ", l2, "lmax = ", lmax
    write(*,*) "ne*nc = ", ne*nc, "timestep = ", tstep
  endif

  0817 format("*****ELEMENT ",I6,2x,I6,2x,I1)
end subroutine cslam_run_bench


end module fvm_bench_mod
