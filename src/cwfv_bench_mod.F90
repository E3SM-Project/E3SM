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
  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, ne, nc, nhc, nhe, nlev, ntrac
  use time_mod, only : timelevel_t
  use element_mod, only : element_t, timelevels
  use cslam_control_volume_mod, only: cslam_struct
  use hybrid_mod, only : hybrid_t
  use cslam_mod, only: cellghostbuf


  public :: cslam_init1,cslam_init2
contains



subroutine cwfv_run_bench(elem,cslam,red,hybrid,nets,nete,tl)
  ! ---------------------------------------------------------------------------------
  use cslam_bsp_mod, only: cslam_bsp, set_boomerang_velocities_gll
  use cslam_mod, only: cslam_run
  ! ---------------------------------------------------------------------------------
  use cslam_line_integrals_mod, only: compute_weights
  ! ---------------------------------------------------------------------------------  
  use cslam_filter_mod, only: monotonic_gradient_cart
  ! ---------------------------------------------------------------------------------
  use cslam_reconstruction_mod, only: reconstruction
  ! ---------------------------------------------------------------------------------
  use checksum_mod, only: test_ghost
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t, derivative_stag_t, derivinit, deriv_print
  ! ---------------------------------------------------------------------------------
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
  type (cslam_struct), intent(inout)             :: cslam(:)
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
  
  real (kind=real_kind)   , dimension(10*(nc+2*nhe)*(nc+2*nhe),6)  :: weights_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_eul_index_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_lgr_index_all
  integer (kind=int_kind)                                          :: jall
    
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons  
  real (kind=real_kind), dimension(nelemd,1-nhe:nc+nhe,1-nhe:nc+nhe) :: area    
  real (kind=real_kind) :: xtmp
  real (kind=longdouble_kind) :: cslam_nodes(nc+1)
  
  integer  choosetrac, chooselev   !for test reason the output
 !-----------------------------------------------------------------------------------!  
 choosetrac=1
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
  ! array 'cslam_nodes(:)' computed below:
  xtmp=nc 
  do i=1,nc+1
    cslam_nodes(i)= 2*(i-1)/xtmp - 1
  end do
  call derivinit(deriv,cslam_corners=cslam_nodes)


!-----------------------------------------------------------------------------------!    
  kptr=0
  do ie=nets,nete
    call cslam_bsp(cslam(ie),tl)
    cslam(ie)%elem_mass=0
    do j=1,nc
      do i=1,nc
        cslam(ie)%elem_mass=cslam(ie)%elem_mass + &
                      cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,chooselev,choosetrac,tl%n0)
      enddo
    enddo
    !
    !first exchange of the initial values
    call ghostVpack(cellghostbuf, cslam(ie)%c,nhc,nc,nlev,ntrac,0,tl%n0,timelevels,elem(ie)%desc)
    ! reset the new unknown
    do k=1,nlev
      do itr=1,ntrac
        do j=1-nhc,nc+nhc
          do i=1-nhc,nc+nhc 
          cslam(ie)%c(i,j,k,itr,tl%np1)=0.0D0
          end do
        end do
      enddo
    enddo
  end do

!-----------------------------------------------------------------------------------!  
  call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc)
!-----------------------------------------------------------------------------------!    

  do ie=nets,nete
    kptr=0
     call ghostVunpack(cellghostbuf, cslam(ie)%c, nhc, nc,nlev,ntrac, 0, tl%n0, timelevels,elem(ie)%desc)
    ! for the mass value
    global_shared_buf(ie,1)=0D0
    global_shared_buf(ie,1)=cslam(ie)%elem_mass
    ! for the max value on the sphere
    tmp1(ie) = MAXVAL(cslam(ie)%c(:,:,chooselev,choosetrac,tl%n0))
    tmp2(ie) = MINVAL(cslam(ie)%c(:,:,chooselev,choosetrac,tl%n0))   
  ! BEGIN Testoutput: write data in p (interpolation) to use existing IO
  ! prepare date for I/O
    do k=1,nlev
      do j=1,nc+1
        do i=1,nc+1 
          !write it in p because of IO 
          !first only with three elements in the patch
          if ((cslam(ie)%cubeboundary==swest) .AND. (j==1) .AND. (i==1)) then
            elem(ie)%state%p(i,j,k,tl%n0)=g*(cslam(ie)%area_sphere(i,j-1)*cslam(ie)%c(i,j-1,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i-1,j)*cslam(ie)%c(i-1,j,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,k,choosetrac,tl%n0))/ &
            (cslam(ie)%area_sphere(i,j-1)+cslam(ie)%area_sphere(i-1,j)+cslam(ie)%area_sphere(i,j))
          elseif ((cslam(ie)%cubeboundary==seast) .AND. (j==1) .AND. (i==nc+1)) then
            elem(ie)%state%p(i,j,k,tl%n0)=g*(cslam(ie)%area_sphere(i-1,j-1)*cslam(ie)%c(i-1,j-1,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i-1,j)*cslam(ie)%c(i-1,j,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,k,choosetrac,tl%n0))/(cslam(ie)%area_sphere(i-1,j-1)+ &
            cslam(ie)%area_sphere(i-1,j)+cslam(ie)%area_sphere(i,j)) 
          elseif ((cslam(ie)%cubeboundary==nwest) .AND. (j==nc+1) .AND. (i==1)) then
            elem(ie)%state%p(i,j,k,tl%n0)=g*(cslam(ie)%area_sphere(i-1,j-1)*cslam(ie)%c(i-1,j-1,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i,j-1)*cslam(ie)%c(i,j-1,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,k,choosetrac,tl%n0))/(cslam(ie)%area_sphere(i-1,j-1)+ &
            cslam(ie)%area_sphere(i,j-1)+cslam(ie)%area_sphere(i,j)) 
          elseif ((cslam(ie)%cubeboundary==neast) .AND. (j==nc+1) .AND. (i==nc+1)) then   
            elem(ie)%state%p(i,j,k,tl%n0)=g*(cslam(ie)%area_sphere(i-1,j-1)*cslam(ie)%c(i-1,j-1,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i,j-1)*cslam(ie)%c(i,j-1,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i-1,j)*cslam(ie)%c(i-1,j,k,choosetrac,tl%n0))/(cslam(ie)%area_sphere(i-1,j-1)+ &
            cslam(ie)%area_sphere(i,j-1)+cslam(ie)%area_sphere(i-1,j))
          else
            elem(ie)%state%p(i,j,k,tl%n0)=g*(cslam(ie)%area_sphere(i-1,j-1)*cslam(ie)%c(i-1,j-1,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i,j-1)*cslam(ie)%c(i,j-1,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i-1,j)*cslam(ie)%c(i-1,j,k,choosetrac,tl%n0)+ &
            cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,k,choosetrac,tl%n0))/(cslam(ie)%area_sphere(i-1,j-1)+ &
            cslam(ie)%area_sphere(i,j-1)+cslam(ie)%area_sphere(i-1,j)+cslam(ie)%area_sphere(i,j))
          end if
        end do
      end do
    end do
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
#ifdef PIO_INTERP
  call interp_movie_init(elem,hybrid,nets,nete,tl=tl)    
  call interp_movie_output(elem,tl, hybrid, 0D0, deriv, nets, nete)
#else
    call shal_movie_init(elem,hybrid,cslam)
    call shal_movie_output(elem,tl, hybrid, 0D0, nets, nete,deriv,cslam)
#endif 
!-----------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------!  
#ifdef PIO_INTERP
  call interp_movie_output(elem, tl, hybrid, 0D0, deriv, nets, nete)
#else     
  call shal_movie_output(elem, tl, hybrid, 0D0, nets, nete,deriv,cslam)
#endif
!-----------------------------------------------------------------------------------!
  if(hybrid%masterthread) then 
    print *
    print *,"Arrival grid created , interpolation points calculated, initialization done. " 
    print *
  endif
  tmp=0
  call t_startf('CWFV')
  
  !BEGIN TIME LOOP, start at 0, calculate then next step
  DO WHILE(tl%nstep<nmax)

! HERE STARTS THE NEW SCHEME --------------------------------------------------------
  kptr=0
  do ie=nets, nete
     !loop through all levels
     do k=1,nlev
        !loop through all tracers
        do itr=1,ntrac
           !reset value
           do j=1-nhc,nc+nhc
              do i=1-nhc,nc+nhc 
                 cslam(ie)%c(i,j,k,itr,tl%np1)=0.0D0
              end do
           end do
           
           !!!! PUT IN NEW SCHEME HERE!!!!
           do j=1,nc
              do i=1,nc
                 cslam(ie)%c(i,j,k,itr,tl%np1)=cslam(ie)%c(i,j,k,itr,tl%n0)
              end do
           end do
        enddo  !End Tracer
     end do  !End Level
     !note write tl%np1 in buffer
     call ghostVpack(cellghostbuf, cslam(ie)%c(:,:,:,:,tl%np1),nhc,nc,nlev,ntrac,kptr,elem(ie)%desc)
  end do

  !-----------------------------------------------------------------------------------!
  !
  ! go one time step forward and do the data exchange, new values shift from tl%np1 to tl%n0
  call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc)
  !-----------------------------------------------------------------------------------!
  do ie=nets,nete
     ! note time level has changed, write buffer in tl%np1 
     call ghostVunpack(cellghostbuf, cslam(ie)%c(:,:,:,:,tl%np1), nhc, nc,nlev,ntrac, kptr, elem(ie)%desc)
  enddo

! HERE ENDS THE NEW SCHEME ----------------------------------------------------------





   do ie=nets,nete
    ! prepare data for I/O
      global_shared_buf(ie,1)=0D0  ! for mass calculation
      global_shared_buf(ie,2)=0D0  ! for max value
      do k=1, nlev
        do j=1,nc+1
          do i=1,nc+1 
            !write it in p because of IO 
            !first only with three elements in the patch
            if ((cslam(ie)%cubeboundary==swest) .AND. (j==1) .AND. (i==1)) then
              elem(ie)%state%p(i,j,k,tl%np1)=g*(cslam(ie)%area_sphere(i,j-1)*cslam(ie)%c(i,j-1,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i-1,j)*cslam(ie)%c(i-1,j,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,k,choosetrac,tl%np1))/ &
              (cslam(ie)%area_sphere(i,j-1)+cslam(ie)%area_sphere(i-1,j)+cslam(ie)%area_sphere(i,j))
            elseif ((cslam(ie)%cubeboundary==seast) .AND. (j==1) .AND. (i==nc+1)) then
              elem(ie)%state%p(i,j,k,tl%np1)=g*(cslam(ie)%area_sphere(i-1,j-1)*cslam(ie)%c(i-1,j-1,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i-1,j)*cslam(ie)%c(i-1,j,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,k,choosetrac,tl%np1))/(cslam(ie)%area_sphere(i-1,j-1)+ &
              cslam(ie)%area_sphere(i-1,j)+cslam(ie)%area_sphere(i,j)) 
            elseif ((cslam(ie)%cubeboundary==nwest) .AND. (j==nc+1) .AND. (i==1)) then
              elem(ie)%state%p(i,j,k,tl%np1)=g*(cslam(ie)%area_sphere(i-1,j-1)*cslam(ie)%c(i-1,j-1,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i,j-1)*cslam(ie)%c(i,j-1,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,k,choosetrac,tl%np1))/(cslam(ie)%area_sphere(i-1,j-1)+ &
              cslam(ie)%area_sphere(i,j-1)+cslam(ie)%area_sphere(i,j)) 
            elseif ((cslam(ie)%cubeboundary==neast) .AND. (j==nc+1) .AND. (i==nc+1)) then   
              elem(ie)%state%p(i,j,k,tl%np1)=g*(cslam(ie)%area_sphere(i-1,j-1)*cslam(ie)%c(i-1,j-1,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i,j-1)*cslam(ie)%c(i,j-1,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i-1,j)*cslam(ie)%c(i-1,j,k,choosetrac,tl%np1))/(cslam(ie)%area_sphere(i-1,j-1)+ &
              cslam(ie)%area_sphere(i,j-1)+cslam(ie)%area_sphere(i-1,j))
            else
              elem(ie)%state%p(i,j,k,tl%np1)=g*(cslam(ie)%area_sphere(i-1,j-1)*cslam(ie)%c(i-1,j-1,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i,j-1)*cslam(ie)%c(i,j-1,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i-1,j)*cslam(ie)%c(i-1,j,k,choosetrac,tl%np1)+ &
              cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,k,choosetrac,tl%np1))/(cslam(ie)%area_sphere(i-1,j-1)+ &
              cslam(ie)%area_sphere(i,j-1)+cslam(ie)%area_sphere(i-1,j)+cslam(ie)%area_sphere(i,j))
            end if 
          end do
        end do
      end do
      ! test mass, just for chooselev and choosetrac
      do j=1,nc
        do i=1,nc        
          mass=mass+cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,chooselev,choosetrac,tl%np1)
          global_shared_buf(ie,1)=global_shared_buf(ie,1)+cslam(ie)%area_sphere(i,j)*cslam(ie)%c(i,j,chooselev,choosetrac,tl%np1)
        end do
      end do
      ! for the max/min value on the sphere
      tmp1(ie) = MAXVAL(cslam(ie)%c(:,:,chooselev,choosetrac,tl%np1))
      tmp2(ie) = MINVAL(cslam(ie)%c(:,:,chooselev,choosetrac,tl%np1))
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

#ifdef PIO_INTERP
    call interp_movie_output(elem, tl, hybrid, 0D0, deriv, nets, nete)
#else     
    call shal_movie_output(elem, tl, hybrid, 0D0, nets, nete,deriv,cslam)
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
    write(*,*) 'number of elements', 6*ne*ne*nc*nc

    print *
    write(*,*) 'chooselev=', chooselev, 'choosetrac=', choosetrac
    write(*,*) 'massbegin', massstart, 'massend', mass 
    write(*,*) 'rel', (mass-massstart)/massstart
    write(*,*) 'maxvaluestart:', maxcstart, 'minvaluestart:', mincstart
    write(*,*) 'maxvalue:     ', maxc,      'minvalue:     ', minc
  endif

  0817 format("*****ELEMENT ",I6,2x,I6,2x,I1)
end subroutine cwfv_run_bench

end module cwfv_bench_mod
