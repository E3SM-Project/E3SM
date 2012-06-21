!-----------------------------------------------------------------------------------!
!MODULE CWFV_MOD-------------------------------------------------------CE-for CWFV!
! CWFV_MOD File for the cwfv project in HOMME                                      !
! Author: Christoph Erath                                                           !
! Date: 13.January 2012                                                             !
! MAIN module to run CWFV on HOMME                                                  !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module cwfv_bench_mod
  use kinds, only : real_kind, int_kind, longdouble_kind
  use edge_mod, only : ghostbuffertr_t, initghostbuffer, freeghostbuffertr, &
                       ghostVpack, ghostVunpack
  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, ne, np, nlev, ntrac
  use time_mod, only : timelevel_t
  use element_mod, only : element_t, timelevels
  use hybrid_mod, only : hybrid_t
  use cwfv_mod, only: cellghostbuf, cwfv_struct


  public :: fvm_init1,fvm_init2
contains



subroutine cwfv_run_bench(elem,fvm,red,hybrid,nets,nete,tl)
  ! ---------------------------------------------------------------------------------
  use fvm_bsp_mod, only: analytical_function, set_boomerang_velocities_gll, get_boomerang_velocities_gll
  use cwfv_mod, only: fvm_mcgregordss, cwfv_run
!    
  ! ---------------------------------------------------------------------------------
  use fvm_line_integrals_mod, only: compute_weights
  ! ---------------------------------------------------------------------------------  
  use fvm_filter_mod, only: monotonic_gradient_cart
  ! ---------------------------------------------------------------------------------
  use fvm_reconstruction_mod, only: reconstruction
  ! ---------------------------------------------------------------------------------
  use checksum_mod, only: test_ghost
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t, derivative_stag_t, derivinit, deriv_print
  ! ---------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------
  use reduction_mod, only : ReductionBuffer_ordered_1d_t
  ! ---------------------------------------------------------------------------------
  use thread_mod, only : nthreads
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
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
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------
  
#ifdef PIO_INTERP
    use interp_movie_mod, only : interp_movie_init, interp_movie_output, interp_movie_finish
#else
    use shal_movie_mod, only : shal_movie_init, shal_movie_output, shal_movie_finish
#endif
  
  implicit none
  type (element_t), intent(inout)                :: elem(:)
  type (cwfv_struct), intent(inout)             :: fvm(:)
  type (ReductionBuffer_ordered_1d_t),intent(in)    :: red   ! reduction buffer         (shared)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)

  real (kind=real_kind)                       :: massstart, mass, maxc, maxcstart,minc, mincstart, tmp  
  real (kind=real_kind)                       :: tmp1(nets:nete), tmp2(nets:nete)
  
  integer                                     :: i,j,k,ie,itr, kptr, jx, jy, jdx, jdy, h
  type (TimeLevel_t)                          :: tl              ! time level struct
  type (derivative_t)                         :: deriv           ! derivative struct
  type (spherical_polar_t)                    :: tmpsphereincart   
 
  character (len=99)                          :: filename
     
  real (kind=real_kind) :: xtmp
  
  integer  choosetrac, chooselev   !for test reason the output
 !-----------------------------------------------------------------------------------!  
 choosetrac=2
 chooselev=1
 
  if(hybrid%masterthread) then 
    print *,"!-----------------------------------------------------------------------!"
    print *,"!  Test CASE for CWFV, Christoph Erath                                  !" 
    print *,"!-----------------------------------------------------------------------!" 
  endif
     
  
  ! Initialize derivative structure
  ! CWFV nodes are equally spaced in alpha/beta
  ! HOMME with equ-angular gnomonic projection maps alpha/beta space
  ! to the reference element via simple scale + translation
  ! thus, CWFV nodes in reference element [-1,1] are a tensor product of
  ! array 'fvm_nodes(:)' computed below:
  call derivinit(deriv)


!-----------------------------------------------------------------------------------!    
  kptr=0
  do ie=nets,nete
    do k=1, nlev
        fvm(ie)%c(:,:,k,1,tl%n0)=1.0D0    !density of the air
        do itr=2,ntrac
          do j=1,np
            do i=1,np               
              call analytical_function(fvm(ie)%c(i,j,k,itr,tl%n0),elem(ie)%spherep(i,j),k,itr)      
            end do
          end do
        end do 
      end do
    !
    !first exchange of the initial values
    call ghostVpack(cellghostbuf, fvm(ie)%c,np,np,nlev,ntrac,0,tl%n0,timelevels,elem(ie)%desc)
  end do
!-----------------------------------------------------------------------------------!  
  call ghost_exchangeV(hybrid,cellghostbuf,np,np)
!-----------------------------------------------------------------------------------!    
  do ie=nets,nete
    kptr=0
     call ghostVunpack(cellghostbuf, fvm(ie)%c, np, np,nlev,ntrac, 0, tl%n0, timelevels,elem(ie)%desc)
    ! for the mass value
!     global_shared_buf(ie,1)=0D0
!     global_shared_buf(ie,1)=fvm(ie)%elem_mass
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
!Initialize Output via geopotential (should be changed, separate output for CWFV
!write first time step to IO 
  do ie=nets,nete
    elem(ie)%state%p(:,:,:,tl%n0)=g*fvm(ie)%c(1:np,1:np,:,choosetrac,tl%n0)
    elem(ie)%state%ps(:,:)=0.0D0
  end do
#ifdef PIO_INTERP
  call interp_movie_init(elem,hybrid,nets,nete,tl=tl)    
  call interp_movie_output(elem,tl, hybrid, 0D0, deriv, nets, nete)
#else
    call shal_movie_init(elem,hybrid,fvm)
    call shal_movie_output(elem,tl, hybrid, 0D0, nets, nete,deriv)
#endif 
!-----------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------!
  if(hybrid%masterthread) then 
    print *
    print *,"Arrival grid created , interpolation points calculated, initialization done. " 
    print *
  endif
  tmp=0
  call t_startf('CWFV')
  
  !BEGIN TIME LOOP, start at 0, calculate then next step
!   DO WHILE(tl%nstep<nmax)
  DO WHILE(tl%nstep<nmax)
! HERE STARTS THE NEW SCHEME --------------------------------------------------------
  kptr=0
! ! start mcgregordss
    do ie=nets,nete
      do k=1,nlev
        elem(ie)%derived%vstar(:,:,:,k)=get_boomerang_velocities_gll(elem(ie), time_at(tl%nstep-1))
        fvm(ie)%vn0(:,:,:,k)=get_boomerang_velocities_gll(elem(ie),time_at(tl%nstep))
      end do
    end do
    
    call fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, tstep, 3)
    call cwfv_run(elem,fvm,hybrid,deriv,tstep,tl,nets,nete)


   do ie=nets,nete
    ! prepare data for I/O
      global_shared_buf(ie,1)=0D0  ! for mass calculation
      global_shared_buf(ie,2)=0D0  ! for max value
      ! test mass, just for chooselev and choosetrac
!       do j=1,nc
!         do i=1,nc        
!           mass=mass+fvm(ie)%area_sphere(i,j)*fvm(ie)%c(i,j,chooselev,choosetrac,tl%np1)
!           global_shared_buf(ie,1)=global_shared_buf(ie,1)+fvm(ie)%area_sphere(i,j)*fvm(ie)%c(i,j,chooselev,choosetrac,tl%np1)
!         end do
!       end do
      ! for the max/min value on the sphere
      tmp1(ie) = MAXVAL(fvm(ie)%c(:,:,chooselev,choosetrac,tl%np1))
      tmp2(ie) = MINVAL(fvm(ie)%c(:,:,chooselev,choosetrac,tl%np1))
    end do

    call TimeLevel_update(tl,"forward") 

!-----------------------------------------------------------------------------------!
    ! for mass calculation
    call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
    mass=global_shared_sum(1)
    maxc = parallelmax(tmp1,hybrid)
    minc = parallelmin(tmp2,hybrid)
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
      print *
    endif
!-----------------------------------------------------------------------------------!  
  do ie=nets,nete
    elem(ie)%state%p(:,:,:,tl%n0)=g*fvm(ie)%c(1:np,1:np,:,choosetrac,tl%n0)
    elem(ie)%state%ps(:,:)=0.0D0
  end do
#ifdef PIO_INTERP
    call interp_movie_output(elem, tl, hybrid, 0D0, deriv, nets, nete)
#else     
    call shal_movie_output(elem, tl, hybrid, 0D0, nets, nete,deriv)
#endif
!-----------------------------------------------------------------------------------!  
  END DO
!------------END TIME LOOP-------------END TIME LOOP--------------END TIME LOOP-----!
!-----------------------------------------------------------------------------------! 
  call t_stopf('CWFV')


  call freeghostbuffertr(cellghostbuf)
 
#ifdef PIO_INTERP
    call interp_movie_finish
#else
    call shal_movie_finish
#endif
!-----------------------------------------------------------------------------------!  

!SUMMARY
  if(hybrid%masterthread) then 
    print *
    print *,"!-----------------------------------------------------------------------!"
    print *,"!  Test CASE for CWFV, Christoph Erath                                  !" 
    print *,"!-----------------------------------------------------------------------!"
    print *  
    write(*,*) 'number of elements', 6*ne*ne

    print *
    write(*,*) 'chooselev=', chooselev, 'choosetrac=', choosetrac
    write(*,*) 'massbegin', massstart, 'massend', mass 
    write(*,*) 'rel', (mass-massstart)/massstart
    write(*,*) 'maxvaluestart:', maxcstart, 'minvaluestart:', mincstart
    write(*,*) 'maxvalue:     ', maxc,      'minvalue:     ', minc
  endif

  0817 format("*****ELEMENT ",I6,2x,I6,2x,I1)
  
!   do ie=nets,nete
!     write(*,*) elem(ie)%LocalId,elem(ie)%GlobalId,elem(ie)%vertex%nbrs(:)%n
!   end do

end subroutine cwfv_run_bench

end module cwfv_bench_mod
