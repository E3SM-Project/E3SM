!-----------------------------------------------------------------------------------!
!MODULE CSLAM_MOD-------------------------------------------------------CE-for CSLAM!
! CSLAM_MOD File for the cslam project in HOMME                                     !
! Author: Christoph Erath                                                           !
! Date: 25.January 2011                                                             !
! MAIN module to run CSLAM on HOMME                                                 !
! 14.November 2011: reorganisation done                                             !
! 7.Februar 2012: cslam_run and cslam_runair                                        !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module cslam_mod

#ifndef MESH

  use kinds, only : real_kind, int_kind, longdouble_kind
  use edge_mod, only : ghostbuffertr_t, initghostbuffer, freeghostbuffertr, &
                       ghostVpack, ghostVunpack,  edgebuffer_t, initEdgebuffer
  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, ne, nc, nhc, nhe, nlev, ntrac, np, ntrac_d
  use time_mod, only : timelevel_t
  use element_mod, only : element_t
  use cslam_control_volume_mod, only: cslam_struct
  use hybrid_mod, only : hybrid_t

  implicit none
  private
  type (ghostBuffertr_t)                      :: cellghostbuf
  type (EdgeBuffer_t)                         :: edgeveloc
  
  public :: cslam_run, cslam_runair, cslam_runairdensity
  public :: cellghostbuf, edgeveloc, cslam_init1,cslam_init2, cslam_mcgregor, cslam_mcgregordss
contains

subroutine cslam_run(elem,cslam,hybrid,deriv,tstep,tl,nets,nete)
  ! ---------------------------------------------------------------------------------
  use cslam_line_integrals_mod, only: compute_weights
  ! ---------------------------------------------------------------------------------  
  use cslam_filter_mod, only: monotonic_gradient_cart
  ! ---------------------------------------------------------------------------------
  use cslam_reconstruction_mod, only: reconstruction
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
  ! ---------------------------------------------------------------------------------
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------
   
  implicit none
  type (element_t), intent(inout)                :: elem(:)
  type (cslam_struct), intent(inout)             :: cslam(:)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
  real (kind=real_kind)                       :: tstep
  
  integer                                     :: i,j,k,ie,itr, jx, jy, jdx, jdy, h
  type (TimeLevel_t)                          :: tl              ! time level struct
  type (derivative_t)                         :: deriv           ! derivative struct
 
  real (kind=real_kind)   , dimension(10*(nc+2*nhe)*(nc+2*nhe),6)  :: weights_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_eul_index_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_lgr_index_all
  integer (kind=int_kind)                                          :: jall
    
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons
  real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer0 
  real (kind=real_kind), dimension(1:nc,1:nc)                        :: tracer1
   
  do ie=nets, nete
    do k=1,nlev
      call cslam_mesh_dep(elem(ie),deriv,cslam(ie),tstep,tl,k)
      !-Departure CSLAM Meshes, initialization done                                                               
      call compute_weights(cslam(ie),6,weights_all,weights_eul_index_all, &
             weights_lgr_index_all,jall) 
      !loop through all tracers
      do itr=1,ntrac
        tracer0=cslam(ie)%c(:,:,k,itr,tl%n0)
        call reconstruction(tracer0, cslam(ie),elem(ie)%corners(1),recons)
        call monotonic_gradient_cart(tracer0, cslam(ie),recons, elem(ie)%desc)
        tracer1=0.0D0   
        do h=1,jall
          jx  = weights_lgr_index_all(h,1)
          jy  = weights_lgr_index_all(h,2)
          jdx = weights_eul_index_all(h,1)
          jdy = weights_eul_index_all(h,2)
              
          call remap(tracer0(jdx,jdy),tracer1(jx,jy),weights_all(h,:),&
                     recons(:,jdx,jdy),cslam(ie)%spherecentroid(:,jdx,jdy))             
        end do
        ! finish scheme
        do j=1,nc
          do i=1,nc
            cslam(ie)%c(i,j,k,itr,tl%np1)=tracer1(i,j)/cslam(ie)%area_sphere(i,j)
          end do
        end do
      enddo  !End Tracer
    end do  !End Level
    !note write tl%np1 in buffer
    call ghostVpack(cellghostbuf, cslam(ie)%c(:,:,:,:,tl%np1),nhc,nc,nlev,ntrac,0,elem(ie)%desc)
  end do
  call t_startf('CSLAM Communication')
  call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc)
  call t_stopf('CSLAM Communication')
  !-----------------------------------------------------------------------------------!
  do ie=nets,nete
    call ghostVunpack(cellghostbuf, cslam(ie)%c(:,:,:,:,tl%np1), nhc, nc,nlev,ntrac, 0, elem(ie)%desc)
  enddo
end subroutine cslam_run



! use this subroutine for benchmark tests, couple airdensity with tracer concentration
subroutine cslam_runairdensity(elem,cslam,hybrid,deriv,tstep,tl,nets,nete)
  ! ---------------------------------------------------------------------------------
  use cslam_line_integrals_mod, only: compute_weights
  ! ---------------------------------------------------------------------------------  
  use cslam_filter_mod, only: monotonic_gradient_cart
  ! ---------------------------------------------------------------------------------
  use cslam_reconstruction_mod, only: reconstruction
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
  ! ---------------------------------------------------------------------------------
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------
   
  implicit none
  type (element_t), intent(inout)                :: elem(:)
  type (cslam_struct), intent(inout)             :: cslam(:)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
  real (kind=real_kind)                       :: tstep
  
  integer                                     :: i,j,k,ie,itr, jx, jy, jdx, jdy, h
  type (TimeLevel_t)                          :: tl              ! time level struct
  type (derivative_t)                         :: deriv           ! derivative struct
 
  real (kind=real_kind)   , dimension(10*(nc+2*nhe)*(nc+2*nhe),6)  :: weights_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_eul_index_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_lgr_index_all
  integer (kind=int_kind)                                          :: jall
    
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons
  real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer0 
  
  real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer_air0   
  real (kind=real_kind), dimension(1:nc,1:nc)                        :: tracer1, tracer_air1 
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons_air   
   
  do ie=nets, nete
    do k=1,nlev
      call cslam_mesh_dep(elem(ie),deriv,cslam(ie),tstep,tl,k)
      !-Departure CSLAM Meshes, initialization done                                                               
      call compute_weights(cslam(ie),6,weights_all,weights_eul_index_all, &
             weights_lgr_index_all,jall) 
             
      tracer_air0=1.0D0 !elem(ie)%derived%dp(:,:,k)       
      call reconstruction(tracer_air0, cslam(ie),elem(ie)%corners(1),recons_air)
      call monotonic_gradient_cart(tracer_air0, cslam(ie),recons_air, elem(ie)%desc)
      tracer_air1=0.0D0   
      do h=1,jall
        jx  = weights_lgr_index_all(h,1)
        jy  = weights_lgr_index_all(h,2)
        jdx = weights_eul_index_all(h,1)
        jdy = weights_eul_index_all(h,2)
        call remap(tracer_air0(jdx,jdy),tracer_air1(jx,jy),weights_all(h,:),&
                   recons_air(:,jdx,jdy),cslam(ie)%spherecentroid(:,jdx,jdy))             
      end do
      ! finish scheme
      do j=1,nc
        do i=1,nc
          tracer_air1(i,j)=tracer_air1(i,j)/cslam(ie)%area_sphere(i,j)
        end do
      end do
      !loop through all tracers
      do itr=1,ntrac
        tracer0=cslam(ie)%c(:,:,k,itr,tl%n0)
        call reconstruction(tracer0, cslam(ie),elem(ie)%corners(1),recons)
        call monotonic_gradient_cart(tracer0, cslam(ie),recons, elem(ie)%desc)
        tracer1=0.0D0   
        do h=1,jall
          jx  = weights_lgr_index_all(h,1)
          jy  = weights_lgr_index_all(h,2)
          jdx = weights_eul_index_all(h,1)
          jdy = weights_eul_index_all(h,2)
         
          call remap_air(tracer0(jdx,jdy),tracer1(jx,jy),tracer_air0(jdx,jdy),&
                       weights_all(h,:), recons(:,jdx,jdy),recons_air(:,jdx,jdy),&
                       cslam(ie)%spherecentroid(:,jdx,jdy))
        end do                     
        ! finish scheme
        do j=1,nc
          do i=1,nc
            cslam(ie)%c(i,j,k,itr,tl%np1)= &
            (tracer1(i,j)/cslam(ie)%area_sphere(i,j))/tracer_air1(i,j)
          end do
        end do            
      enddo  !End Tracer
    end do  !End Level
    !note write tl%np1 in buffer
    call ghostVpack(cellghostbuf, cslam(ie)%c(:,:,:,:,tl%np1),nhc,nc,nlev,ntrac,0,elem(ie)%desc)
  end do
  call t_startf('CSLAM Communication')
  call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc)
  call t_stopf('CSLAM Communication')
  !-----------------------------------------------------------------------------------!
  do ie=nets,nete
     call ghostVunpack(cellghostbuf, cslam(ie)%c(:,:,:,:,tl%np1), nhc, nc,nlev,ntrac, 0, elem(ie)%desc)
  enddo
end subroutine cslam_runairdensity






! use this subroutine for benchmark tests, couple airdensity with tracer concentration
subroutine cslam_runair(elem,cslam,hybrid,deriv,tstep,tl,nets,nete)
  ! ---------------------------------------------------------------------------------
  use cslam_line_integrals_mod, only: compute_weights
  ! ---------------------------------------------------------------------------------  
  use cslam_filter_mod, only: monotonic_gradient_cart
  ! ---------------------------------------------------------------------------------
  use cslam_reconstruction_mod, only: reconstruction
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
  ! ---------------------------------------------------------------------------------
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------
   
  implicit none
  type (element_t), intent(inout)                :: elem(:)
  type (cslam_struct), intent(inout)             :: cslam(:)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
  real (kind=real_kind)                       :: tstep
  
  integer                                     :: i,j,k,ie,itr, jx, jy, jdx, jdy, h
  type (TimeLevel_t)                          :: tl              ! time level struct
  type (derivative_t)                         :: deriv           ! derivative struct
 
  real (kind=real_kind)   , dimension(10*(nc+2*nhe)*(nc+2*nhe),6)  :: weights_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_eul_index_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_lgr_index_all
  integer (kind=int_kind)                                          :: jall
    
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons
  real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer0 
  
  real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer_air0   
  real (kind=real_kind), dimension(1:nc,1:nc)                        :: tracer1, tracer_air1 
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons_air   
   
  do ie=nets, nete
    do k=1,nlev
      call cslam_mesh_dep(elem(ie),deriv,cslam(ie),tstep,tl,k)
      !-Departure CSLAM Meshes, initialization done                                                               
      call compute_weights(cslam(ie),6,weights_all,weights_eul_index_all, &
             weights_lgr_index_all,jall) 
      !loop through all tracers
      do itr=1,ntrac
        tracer0=cslam(ie)%c(:,:,k,itr,tl%n0)
        call reconstruction(tracer0, cslam(ie),elem(ie)%corners(1),recons)
        call monotonic_gradient_cart(tracer0, cslam(ie),recons, elem(ie)%desc)
        tracer1=0.0D0   
        if (itr==1) then !calculation for air density  (the first tracer is supposed to be)
          recons_air=recons
          tracer_air0=tracer0
          do h=1,jall
            jx  = weights_lgr_index_all(h,1)
            jy  = weights_lgr_index_all(h,2)
            jdx = weights_eul_index_all(h,1)
            jdy = weights_eul_index_all(h,2)
                
            call remap(tracer0(jdx,jdy),tracer1(jx,jy),weights_all(h,:),&
                       recons(:,jdx,jdy),cslam(ie)%spherecentroid(:,jdx,jdy))             
          end do
          ! finish scheme
          do j=1,nc
            do i=1,nc
              tracer_air1(i,j)=tracer1(i,j)/cslam(ie)%area_sphere(i,j)
              cslam(ie)%c(i,j,k,itr,tl%np1)=tracer_air1(i,j)
            end do
          end do
        else !calculation for the other tracers
          do h=1,jall
            jx  = weights_lgr_index_all(h,1)
            jy  = weights_lgr_index_all(h,2)
            jdx = weights_eul_index_all(h,1)
            jdy = weights_eul_index_all(h,2)
           
            call remap_air(tracer0(jdx,jdy),tracer1(jx,jy),tracer_air0(jdx,jdy),&
                         weights_all(h,:), recons(:,jdx,jdy),recons_air(:,jdx,jdy),&
                         cslam(ie)%spherecentroid(:,jdx,jdy))
          end do                     
          ! finish scheme
          do j=1,nc
            do i=1,nc
              cslam(ie)%c(i,j,k,itr,tl%np1)= &
              (tracer1(i,j)/cslam(ie)%area_sphere(i,j))/tracer_air1(i,j)
            end do
          end do            
        end if
      enddo  !End Tracer
    end do  !End Level
    !note write tl%np1 in buffer
    call ghostVpack(cellghostbuf, cslam(ie)%c(:,:,:,:,tl%np1),nhc,nc,nlev,ntrac,0,elem(ie)%desc)
  end do
  call t_startf('CSLAM Communication')
  call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc)
  call t_stopf('CSLAM Communication')
  !-----------------------------------------------------------------------------------!
  do ie=nets,nete
     call ghostVunpack(cellghostbuf, cslam(ie)%c(:,:,:,:,tl%np1), nhc, nc,nlev,ntrac, 0, elem(ie)%desc)
  enddo
end subroutine cslam_runair


! do the remapping
subroutine remap(tracer0,tracer1,weights,recons,centroid)
  real (kind=real_kind), intent(in)               :: tracer0
  real (kind=real_kind), intent(inout)            :: tracer1
  real (kind=real_kind), intent(in)               :: weights(6),recons(5),centroid(5)
    
  tracer1 = tracer1+weights(1)*(&
       ! all constant terms 
       tracer0 - recons(1)*centroid(1) - recons(2)*centroid(2) &
       + recons(3)*(2.0D0*centroid(1)**2 -centroid(3)) &
       + recons(4)*(2.0D0*centroid(2)**2 -centroid(4)) &
       + recons(5)*(2.0D0*centroid(1)*centroid(2)-centroid(5))) + &
       ! linear terms
       weights(2)*&
       (recons(1) - recons(3)*2.0D0*centroid(1)- recons(5)*centroid(2)) + &
       weights(3)*&
       (recons(2) - recons(4)*2.0D0*centroid(2)- recons(5)*centroid(1)) + &
       ! quadratic terms
       weights(4)*recons(3)+&
       weights(5)*recons(4)+&
       weights(6)*recons(5)
end subroutine remap

! do remapping with air (i.e. conserve mass of air density * concentration),
! see Nair et.al 2010 in JCP: A class of deformational flow test cases for linear transport
! schemes on the sphere, Appendix B
subroutine remap_air(tracer0,tracer1,tracer_air,weights,recons, recons_air,centroid)
  real (kind=real_kind), intent(in)         :: tracer0, tracer_air
  real (kind=real_kind), intent(inout)      :: tracer1
  real (kind=real_kind), intent(in)         :: weights(6),recons(5),recons_air(5),centroid(5)
    
  tracer1 = tracer1+&
        ! air density times tracer reconstruction
        tracer_air*(weights(1)*(&      ! 1 is for air
        ! all constant terms 
        tracer0 - recons(1)*centroid(1) - recons(2)*centroid(2) &
        + recons(3)*(2.0D0*centroid(1)**2 -centroid(3)) &
        + recons(4)*(2.0D0*centroid(2)**2 -centroid(4)) &
        + recons(5)*(2.0D0*centroid(1)*centroid(2)-centroid(5))) + &
        ! linear terms
        weights(2)* &
        (recons(1)- recons(3)*2.0D0*centroid(1) - recons(5)*centroid(2)) + &
        weights(3)* &
        (recons(2) - recons(4)*2.0D0*centroid(2) - recons(5)*centroid(1)) + &
        ! quadratic terms
        weights(4)*recons(3) + &
        weights(5)*recons(4) + &
        weights(6)*recons(5)) + &

        !tracer times air reconstruction
        tracer0*(weights(1)*(&      
        ! all constant terms 
!       tracer_air &  this term cancels it out
        - recons_air(1)*centroid(1) - recons_air(2)*centroid(2) &
        + recons_air(3)*(2.0D0*centroid(1)**2 -centroid(3)) &
        + recons_air(4)*(2.0D0*centroid(2)**2 -centroid(4)) &
        + recons_air(5)*(2.0D0*centroid(1)*centroid(2)-centroid(5))) + &
        ! linear terms
        weights(2)* &
        (recons_air(1) - recons_air(3)*2.0D0*centroid(1) - recons_air(5)*centroid(2)) + &
        weights(3)* &
        (recons_air(2) - recons_air(4)*2.0D0*centroid(2) - recons_air(5)*centroid(1)) + &
        ! quadratic terms
        weights(4)*recons_air(3)+&
        weights(5)*recons_air(4)+&
        weights(6)*recons_air(5))
end subroutine remap_air


! initialize global buffers shared by all threads
subroutine cslam_init1(par)
  use parallel_mod, only : parallel_t, haltmp
  type (parallel_t) :: par
  
  if (nc<4) then
     if (par%masterproc) then
        print *, "NUMBER OF CELLS ERROR for CSLAM: Number of cells parmeter"
        print *, "parameter nc at least 4 (nc>=4), nc*nc cells per element. This is"
        print *, "needed for the cubic reconstruction, which is only implemented yet! STOP"
     endif
     call haltmp("stopping")
  end if
  if (nhe .ne. 1) then
     if (par%masterproc) then
        print *, "PARAMTER ERROR for CSLAM: Number of halo zone for the extended"
        print *,"element nhe has to be 1, only this is available now! STOP!"
     endif
     call haltmp("stopping")
  end if
  if (ntrac>ntrac_d) then
     if (par%masterproc) print *,'ntrac,ntrac_d=',ntrac,ntrac_d
     call haltmp("PARAMTER ERROR for CSLAM: ntrac > ntrac_d")
  endif

  call initghostbuffer(cellghostbuf,nlev,ntrac+1,nhc,nc) !+1 for the air_density, which comes from SE
  
  call initEdgebuffer(edgeveloc,2*nlev)
end subroutine cslam_init1



! initialization that can be done in threaded regions
subroutine cslam_init2(elem,cslam,hybrid,nets,nete,tl)
  use cslam_control_volume_mod, only: cslam_mesh_ari
  use cslam_analytic_mod, only: computexytosphere_moments
  use bndry_mod, only: compute_ghost_corner_orientation 


  type (timelevel_t) :: tl
  type (cslam_struct) :: cslam(:)
  type (element_t) :: elem(:)
  type (hybrid_t)                             :: hybrid
  integer :: ie,nets,nete

  call compute_ghost_corner_orientation(hybrid,elem,nets,nete)
  ! run some tests:
  !   call test_ghost(hybrid,elem,nets,nete)
  

  do ie=nets,nete
    call cslam_mesh_ari(elem(ie),cslam(ie),tl)
    call computexytosphere_moments(cslam(ie),elem(ie)%desc)
  enddo

end subroutine cslam_init2


! ----------------------------------------------------------------------------------!
!SUBROUTINE CSLAM_MESH_DEP----------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 30.October 2011                                          !
! DESCRIPTION: Calculates the departure mesh                                        !
!              Please select the test cases by an COMMENT/UNCOMMENT basis           !
!                                                                                   !
! CALLS: solidbody, cart2cubedspherexy, analytical_function                         !
! INTPUT/OUTPUT:                                                                    !
!        cslam   ...  structure                                                     !
! INPUT: nstep   ... actual step                                                    !
!-----------------------------------------------------------------------------------!
subroutine cslam_mesh_dep(elem, deriv, cslam, dt, tl, klev)
  use coordinate_systems_mod, only : cartesian2D_t
  use element_mod, only : element_t
  use cslam_control_volume_mod, only: cslam_struct
  use time_mod, only : timelevel_t, time_at
  use parallel_mod, only : haltmp
  use control_mod, only : test_cfldep

  use derivative_mod, only : derivative_t
#ifdef _CSLAM
  use cslam_bsp_mod, only: boomerang, solidbody
#endif
  implicit none
  type (derivative_t)  , intent(in) :: deriv
  type (cslam_struct), intent(inout)   :: cslam
  type (timelevel_t),intent(in)        :: tl
  integer,intent(in)                   :: klev
  
  type (element_t), intent(inout)      :: elem
  real (kind=real_kind)                :: time,dt, dx, dy, cflx, cfly, maxcflx, maxcfly  
  integer                              :: i,j
  type (cartesian2D_t)                 :: dcart(nc+1,nc+1)
  

! for the benchmark test, use more accurate departure point creation
#if 0
!CE: define new mesh for fvm cslam on an equal angular grid
! go from alpha,beta -> cartesian xy on cube -> lat lon on the sphere
#ifdef _CSLAM
  do j=1,nc+1
     do i=1,nc+1               
!                 call solidbody(cslam%asphere(i,j), cslam%dsphere(i,j))
        call boomerang(cslam%asphere(i,j), cslam%dsphere(i,j),tl%nstep)
     end do
  end do
#endif
#else

! for given velocities in the element structure
  call cslam_dep_from_gll(elem, deriv, cslam%asphere,cslam%dsphere,dt,tl,klev)
#endif

   if (test_cfldep) then
     call check_departurecell(cslam,klev) 
   endif 
  
end subroutine cslam_mesh_dep


! ----------------------------------------------------------------------------------!
!SUBROUTINE CSLAM_DEP_FROM_GLL------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR 14. December 2011                            !
! DESCRIPTION: calculates the deparute grid for CSLAM coming from the gll points    !
!                                                                                   !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine cslam_dep_from_gll(elem, deriv, asphere,dsphere,dt,tl,klev)
  use physical_constants, only : DD_PI, rearth
  use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t, change_coordinates
  use time_mod, only : timelevel_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, interpolate_gll2cslam_corners

  implicit none
  type (element_t), intent(in)          :: elem
  type (derivative_t)  , intent(in) :: deriv
  type (spherical_polar_t),intent(in)   :: asphere(nc+1,nc+1)
  type (spherical_polar_t),intent(out)  :: dsphere(nc+1,nc+1)
  type (timelevel_t),intent(in)        :: tl
  real (kind=real_kind),intent(in)      :: dt
  integer,intent(in)                   :: klev
  
  real (kind=real_kind)                 :: uxyz_gll(np,np,3),uxyz(nc+1,nc+1,3)
  real (kind=real_kind)                 :: un,ue,ur,clon,clat,slon,slat
  type(cartesian3D_t)                   :: acart
  integer :: i,j,k
  
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
   ! interpolate velocity to CSLAM nodes
   do i=1,3
      uxyz(:,:,i)=interpolate_gll2cslam_corners(uxyz_gll(:,:,i),deriv)
   end do 
   ! compute departure point 
   ! crude, 1st order accurate approximation.  to be improved
   do i=1,nc+1
      do j=1,nc+1
         ! note: asphere()%r=1, so we need to convert velocity to radians/sec:
         acart = change_coordinates(asphere(i,j))  
         acart%x = acart%x - dt*uxyz(i,j,1)/rearth
         acart%y = acart%y - dt*uxyz(i,j,2)/rearth
         acart%z = acart%z - dt*uxyz(i,j,3)/rearth
         dsphere(i,j) = change_coordinates(acart)
         dsphere(i,j)%r = asphere(i,j)%r
      enddo
   enddo
  
end subroutine cslam_dep_from_gll
!END SUBROUTINE CSLAM_DEP_FROM_GLL--------------------------------------CE-for CSLAM!

!SUBROUTINE CHECK_DEPARTURECELL-----------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: Check if the departure cell is degenerated or if the                 !
!              CFL number > 1.0D0                                                   !
!         IF THE CHECK FAILS, THE PROGRAM WILL STOP                                 !
!                                                                                   !
! INPUT:  cslam  ... cslam structur                                                 ! 
!         klev   ... Level (vertical)                                               !
!-----------------------------------------------------------------------------------!
subroutine check_departurecell(cslam,klev)
  use coordinate_systems_mod, only: cartesian2D_t, cart2cubedspherexy, spherical_to_cart
  use cslam_control_volume_mod, only:  cslam_struct                                         
  implicit none

  type (cslam_struct), intent(inout)   :: cslam
  integer,intent(in)                   :: klev
  
  type (cartesian2D_t)                 :: dcart(nc+1,nc+1)
  real (kind=real_kind)                :: cflx, cfly, maxcflx, maxcfly  
  integer                              :: i,j
  logical                              :: crossline
            
  ! calculate xy Cartesian on the cube of departure points on the corresponding face  
  ! the calculation of dcart is also done in cslam_line_integrals_mod.F90, one could 
  ! save the points here...
  ! BUT: this subroutine should only be used for test reasons
  do j=1,nc+1
     do i=1,nc+1               
        call cart2cubedspherexy(spherical_to_cart(cslam%dsphere(i,j)),&
             cslam%faceno,dcart(i,j))              
     end do
  end do                        
  cflx=0.0D0
  cfly=0.0D0
  cslam%maxcfl(1,klev)=cflx
  cslam%maxcfl(2,klev)=cfly
  crossline=.FALSE.
  do j=1,nc
     do i=1,nc  
      ! equidistant mesh in alpha/beta coordinates
      cflx=abs((atan(dcart(i,j)%x)-atan(cslam%acartx(i)))/cslam%dalpha)
      cfly=abs((atan(dcart(i,j)%y)-atan(cslam%acarty(j)))/cslam%dbeta)
      if(cflx>cslam%maxcfl(1,klev)) cslam%maxcfl(1,klev)=cflx
      if(cfly>cslam%maxcfl(2,klev)) cslam%maxcfl(2,klev)=cfly
    
      !one could stop immediately here, if crossline=.TRUE., but want to calculate 
      ! all CFL first
      call check_lines_cross(dcart(i,j),dcart(i+1,j),dcart(i,j+1),dcart(i+1,j+1),crossline)  
      call check_lines_cross(dcart(i,j),dcart(i,j+1),dcart(i+1,j),dcart(i+1,j+1),crossline)
     end do
  end do  
  ! nodes on the north and east boundary
  do i=1,nc+1  
    cflx=abs((atan(dcart(i,nc+1)%x)-atan(cslam%acartx(i)))/cslam%dalpha)
    cfly=abs((atan(dcart(nc+1,j)%y)-atan(cslam%acarty(j)))/cslam%dbeta)
    if(cflx>cslam%maxcfl(1,klev)) cslam%maxcfl(1,klev)=cflx
    if(cfly>cslam%maxcfl(2,klev)) cslam%maxcfl(2,klev)=cfly
  end do

  if (cslam%maxcfl(1,klev) > 1.0D0 .OR. cslam%maxcfl(2,klev) > 1.0D0) then
    write(*,*) "Error in cslam_mod.F90: CFL number too big"
    write(*,*) "CFL has to be < 1.0D0"
    write(*,*) "Choose a smaller time step!"
    write(*,*) "max CFL in this element: maxcflx", cslam%maxcfl(1,klev), "maxcfly", cslam%maxcfl(2,klev) 
    STOP "Exit program!"
  endif

  if (crossline) then
    write(*,*) "FATAL Error in cslam_mod.F90: departure cell is degenerated!"
    write(*,*) "Choose a smaller time step! max CFL in this element: maxcflx", cslam%maxcfl(1,klev), "maxcfly", cslam%maxcfl(2,klev) 
    STOP "Exit program!"
  endif
end subroutine check_departurecell


! see, e.g., http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
subroutine check_lines_cross(p1,p2,q1,q2,crossline)
  use coordinate_systems_mod, only: cartesian2D_t
  implicit none
  type(cartesian2D_t), intent (in)  :: p1,p2,q1,q2
  logical            , intent(inout)  :: crossline
  !
  ! local workspace
  !
  REAL (kind=real_kind)    :: dn,tp,tq
  
  dn = (q2%y-q1%y)*(p2%x-p1%x)-(q2%x-q1%x)*(p2%y-p1%y)

  if (abs(dn)>1.0D-12) then
    tp = ((q2%x-q1%x)*(p1%y-q1%y)-(q2%y-q1%y)*(p1%x-q1%x))/dn
    tq = ((p2%x-p1%x)*(p1%y-q1%y)-(p2%y-p1%y)*(p1%x-q1%x))/dn
    ! implement the next two lines faster!!!!!!
    if (tp>=0.0D0 .and. tp<=1.0D00 .and. &
        tq>=0.0D0 .and. tq<=1.0D00) then
        crossline=.TRUE.
    endif
  endif
end subroutine check_lines_cross

!SUBROUTINE CHECK_DEPARTURECELLEST--------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: Check if the departure cell is degenerated or if the                 !
!              CFL number > 1.0D0                                                   !
!         IF THE CHECK FAILS, THE PROGRAM WILL STOP                                 !
!                                                                                   !
! INPUT:  acartx, acarty  ... arrival grid coordinates in the element               ! 
!         dcart           ... departure grid coordinates in the element             !                                         !
!                                                                                   !
!OLD AND DIFFERENT CHECK FOR THE SHAPE OF THE DEPARTURE CELL------------------------!
subroutine check_departurecellest(acartx,acarty,dcart)
  use coordinate_systems_mod, only: cartesian2D_t
  implicit none

  real (kind=real_kind), intent(in)    :: acartx(-nhe:nc+2+nhe), acarty(-nhe:nc+2+nhe)
  type (cartesian2D_t), intent(in)     :: dcart(nc+1,nc+1)
  real (kind=real_kind)                :: dx, dy, cflx, cfly, maxcflx, maxcfly  
  integer                              :: i,j
  logical                              :: crossline
                          
  cflx=0.0D0
  cfly=0.0D0
  maxcflx=cflx
  maxcfly=cfly
  do j=1,nc
     do i=1,nc  
      !estimate (rough) CFL number
      dx=acartx(i)-acartx(i+1)
      dy=acarty(j)-acarty(j+1)

      cflx=abs((dcart(i,j)%x-acartx(i))/dx)
      cfly=abs((dcart(i,j)%y-acarty(j))/dy)
      
      if(cflx>maxcflx) maxcflx=cflx
      if(cfly>maxcfly) maxcfly=cfly
      call check_lines_cross(dcart(i,j),dcart(i+1,j),dcart(i,j+1),dcart(i+1,j+1),crossline)  
      call check_lines_cross(dcart(i,j),dcart(i,j+1),dcart(i+1,j),dcart(i+1,j+1),crossline)
     end do
  end do  
  
  do i=1,nc  
    !estimate (rough) CFL number
    dx=acartx(i)-acartx(i+1)
    dy=acarty(j)-acarty(j+1)

    cflx=abs((dcart(i,nc+1)%x-acartx(i))/dx)
    cfly=abs((dcart(nc+1,j)%y-acarty(j))/dy)
      
    if(cflx>maxcflx) maxcflx=cflx
    if(cfly>maxcfly) maxcfly=cfly
  end do
  
  if (maxcflx > 1.0D0 .OR. maxcfly > 1.0D0) then
    write(*,*) "FATAL Error in cslam_line_integrals.F90: CFL number too big"
    write(*,*) "CFL has to be < 1.0D0"
    write(*,*) "Choose a smaller time step!"
    write(*,*) "max CFL in this element: maxcflx", maxcflx, "maxcfly", maxcfly 
    STOP "Exit program!"
  endif
  
  if (crossline) then
    write(*,*) "FATAL Error in cslam_line_integrals.F90: departure cell is degenerated!"
    write(*,*) "Choose a smaller time step!"
    write(*,*) "max CFL in this element: maxcflx", maxcflx, "maxcfly", maxcfly 
    STOP "Exit program!"
  endif
end subroutine check_departurecellest


! ----------------------------------------------------------------------------------!
!SUBROUTINE CSLAM_MCGREGOR----------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 09. January 2012                                         !
! DESCRIPTION: ! using McGregor AMS 1993 scheme: Economical Determination of        !
!                Departure Points for Semi-Lagrangian Models                        !
!                                                                                   !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine cslam_mcgregor(elem, deriv, tstep, vhat, vstar,order)
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, gradient_sphere, ugradv_sphere, vorticity_sphere
  implicit none

  type (element_t), intent(in)                                :: elem
  type (derivative_t), intent(in)                             :: deriv      ! derivative struct
  real (kind=real_kind), intent(in)                           :: tstep
  real (kind=real_kind), dimension(np,np,2), intent(inout)    :: vstar
  real (kind=real_kind), dimension(np,np,2), intent(in)       :: vhat
  
  integer, intent(in)                                         :: order

  integer                                            :: i
  real (kind=real_kind), dimension(np,np,2)          :: ugradv
  real (kind=real_kind)                              :: timetaylor
    
  ugradv=vstar
  timetaylor=1
  do i=1,order
!     tmp = 0.5D0*(vgradv(:,:,1)**2 + vgradv(:,:,2)**2) 
! 
!     gradvstar = gradient_sphere(tmp,deriv,elem%Dinv)    ! scalar -> latlon vector
!     tmp = vorticity_sphere(vgradv,deriv,elem)                 ! latlon vector -> scalar 
! 
!     ! v \nabla v expressed through gradient of mean scalar and relative velocity
!     ! see (17) in Williamson et.al. 1992, JCP 102 (211-224)
!     vgradv(:,:,1)= gradvstar(:,:,1) - vstarold(:,:,2)*tmp(:,:)
!     vgradv(:,:,2)= gradvstar(:,:,2) + vstarold(:,:,1)*tmp(:,:)
!     
!     timetaylor=-timetaylor*tstep/(i+1)
!     vstar=vstar + timetaylor*vgradv
    
    ugradv=ugradv_sphere(vhat,ugradv,deriv,elem)
    timetaylor=-timetaylor*tstep/(i+1)
    
    vstar=vstar + timetaylor*ugradv  
  end do
end subroutine cslam_mcgregor
!END SUBROUTINE CSLAM_MCGREGOR------------------------------------------CE-for CSLAM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE CSLAM_MCGREGORDSS-------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 26. May 2012                                             !
! DESCRIPTION: ! using McGregor AMS 1993 scheme: Economical Determination of        !
!                Departure Points for Semi-Lagrangian Models                        !
!                McGegror version with DSS every ugradv                             !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine cslam_mcgregordss(elem,cslam,nets,nete, hybrid, deriv, tstep, ordertaylor)
  use derivative_mod, only : derivative_t, ugradv_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  
  implicit none

  type (element_t), intent(inout)                :: elem(:)
  type (cslam_struct), intent(in)                :: cslam(:)

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
        vhat(ie,:,:,:,:)=(cslam(ie)%vn0(:,:,:,:) + ugradv(ie,:,:,:,:))/2 
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

end subroutine cslam_mcgregordss
!END SUBROUTINE CSLAM_MCGREGORDSS---------------------------------------CE-for CSLAM!
#endif

end module cslam_mod
