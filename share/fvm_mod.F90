!-----------------------------------------------------------------------------------!
!MODULE FVM_MOD-----------------------------------------------------------CE-for FVM!
! FVM_MOD File for the fvm project in HOMME                                         !
! Author: Christoph Erath                                                           !
! Date: 25.January 2011                                                             !
! MAIN module to run fvm on HOMME                                                   !
! 14.November 2011: reorganisation done                                             !
! 7.Februar 2012: cslam_run and cslam_runair                                            !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module fvm_mod

  use kinds, only : real_kind, int_kind, longdouble_kind
  use edge_mod, only : ghostbuffertr_t, initghostbuffer, freeghostbuffertr, &
                       ghostVpack, ghostVunpack,  edgebuffer_t, initEdgebuffer
  use bndry_mod, only: ghost_exchangeV                     
  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, ne, nc, nhc, nhe, nlev, ntrac, np, ntrac_d
  use time_mod, only : timelevel_t
  use element_mod, only : element_t, timelevels
  use fvm_control_volume_mod, only: fvm_struct
  use hybrid_mod, only : hybrid_t

  implicit none
  private
  type (ghostBuffertr_t)                      :: cellghostbuf
  type (EdgeBuffer_t)                         :: edgeveloc
  
!   public :: cslam_run, cslam_runair, cslam_runairdensity
  public :: cslam_runairdensity
  
  public :: cellghostbuf, edgeveloc, fvm_init1,fvm_init2, fvm_mcgregor, fvm_mcgregordss
  public :: fvm_init3
contains

! subroutine cslam_run(elem,fvm,hybrid,deriv,tstep,tl,nets,nete)
!   ! ---------------------------------------------------------------------------------
!   use fvm_line_integrals_mod, only: compute_weights
!   ! ---------------------------------------------------------------------------------  
!   use fvm_filter_mod, only: monotonic_gradient_cart
!   ! ---------------------------------------------------------------------------------
!   use fvm_reconstruction_mod, only: reconstruction
!   ! ---------------------------------------------------------------------------------
!   use derivative_mod, only : derivative_t
!   ! ---------------------------------------------------------------------------------
!   use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
!   ! -----------------------------------------------
!    
!   implicit none
!   type (element_t), intent(inout)                :: elem(:)
!   type (fvm_struct), intent(inout)             :: fvm(:)
!   type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
!   integer, intent(in)                         :: nets  ! starting thread element number (private)
!   integer, intent(in)                         :: nete  ! ending thread element number   (private)
!   real (kind=real_kind)                       :: tstep
!   
!   integer                                     :: i,j,k,ie,itr, jx, jy, jdx, jdy, h
!   type (TimeLevel_t)                          :: tl              ! time level struct
!   type (derivative_t)                         :: deriv           ! derivative struct
!  
!   real (kind=real_kind)   , dimension(10*(nc+2*nhe)*(nc+2*nhe),6)  :: weights_all
!   integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_eul_index_all
!   integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_lgr_index_all
!   integer (kind=int_kind)                                          :: jall
!     
!   real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons
!   real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer0 
!   real (kind=real_kind), dimension(1:nc,1:nc)                        :: tracer1
!    
!   do ie=nets, nete
!     do k=1,nlev
!       call fvm_mesh_dep(elem(ie),deriv,fvm(ie),tstep,tl,k)
!       !-Departure fvm Meshes, initialization done                                                               
!       call compute_weights(fvm(ie),6,weights_all,weights_eul_index_all, &
!              weights_lgr_index_all,jall) 
!       !loop through all tracers
!       do itr=1,ntrac
!         tracer0=fvm(ie)%c(:,:,k,itr,tl%n0)
!         call reconstruction(tracer0, fvm(ie),recons)
!         call monotonic_gradient_cart(tracer0, fvm(ie),recons, elem(ie)%desc)
!         tracer1=0.0D0   
! !         do h=1,jall
! !           jx  = weights_lgr_index_all(h,1)
! !           jy  = weights_lgr_index_all(h,2)
! !           jdx = weights_eul_index_all(h,1)
! !           jdy = weights_eul_index_all(h,2)
! !               
! !           call cslam_remap(tracer0(jdx,jdy),tracer1(jx,jy),weights_all(h,:),&
! !                      recons(:,jdx,jdy),fvm(ie)%spherecentroid(:,jdx,jdy))             
! !         end do
!         call cslam_remap(tracer0,tracer1,weights_all, recons, &
!                    fvm(ie)%spherecentroid, weights_eul_index_all, weights_lgr_index_all, jall)
!         ! finish scheme
!         do j=1,nc
!           do i=1,nc
!             fvm(ie)%c(i,j,k,itr,tl%np1)=tracer1(i,j)/fvm(ie)%area_sphere(i,j)
!           end do
!         end do
!       enddo  !End Tracer
!     end do  !End Level
!     !note write tl%np1 in buffer
!     call ghostVpack(cellghostbuf, fvm(ie)%c,nhc,nc,nlev,ntrac,0,tl%np1,timelevels,elem(ie)%desc)
!   end do
!   call t_startf('FVM Communication')
!   call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc)
!   call t_stopf('FVM Communication')
!   !-----------------------------------------------------------------------------------!
!   do ie=nets,nete
!      call ghostVunpack(cellghostbuf, fvm(ie)%c, nhc, nc,nlev,ntrac, 0, tl%np1, timelevels,elem(ie)%desc)
!   enddo
! end subroutine cslam_run



! use this subroutine for benchmark tests, couple airdensity with tracer concentration
subroutine cslam_runairdensity(elem,fvm,hybrid,deriv,tstep,tl,nets,nete)
  ! ---------------------------------------------------------------------------------
  use fvm_line_integrals_mod, only: compute_weights
  ! ---------------------------------------------------------------------------------  
  use fvm_filter_mod, only: monotonic_gradient_cart
  ! ---------------------------------------------------------------------------------
  use fvm_reconstruction_mod, only: reconstruction
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t
  ! ---------------------------------------------------------------------------------
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------
  use edge_mod, only :  ghostBuffertr_t,ghostVpack2d_level, ghostVunpack2d_level,initghostbuffer,freeghostbuffertr
   
   
  implicit none
  type (element_t), intent(inout)                :: elem(:)
  type (fvm_struct), intent(inout)             :: fvm(:)
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
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: reconsalt
  
  real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer0 
  
  real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer_air0   
  real (kind=real_kind), dimension(1:nc,1:nc)                        :: tracer1, tracer_air1 
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons_air   
  
  type (ghostBuffertr_t)                      :: buflatlon
  
  
  call t_startf('CSLAM scheme') 
  
  do ie=nets, nete
    do k=1,nlev
      call fvm_mesh_dep(elem(ie),deriv,fvm(ie),tstep,tl,k)
    end do
  end do
  
  call initghostbuffer(buflatlon,nlev,2,2,nc+1)    ! use the tracer entry 2 for lat lon
  
  do ie=nets,nete
    call ghostVpack2d_level(buflatlon,fvm(ie)%dsphere(:,:,:)%lat,1,2, nc+1,nlev,elem(ie)%desc) !kptr = 1 for lat
    call ghostVpack2d_level(buflatlon,fvm(ie)%dsphere(:,:,:)%lon,2,2, nc+1,nlev,elem(ie)%desc) !kptr =2 for lon
  end do
!-----------------------------------------------------------------------------------! 
  call ghost_exchangeV(hybrid,buflatlon,2,nc+1)
!-----------------------------------------------------------------------------------!  
  do ie=nets,nete
    call ghostVunpack2d_level(buflatlon,fvm(ie)%dsphere(:,:,:)%lat,1,2, nc+1,nlev,elem(ie)%desc)
    call ghostVunpack2d_level(buflatlon,fvm(ie)%dsphere(:,:,:)%lon,2,2, nc+1,nlev,elem(ie)%desc)
    fvm(ie)%dsphere(:,:,:)%r=1.0D0  !!! RADIUS IS ASSUMED TO BE 1.0DO !!!!
    
  end do
  call freeghostbuffertr(buflatlon)
  
  do ie=nets, nete
    do k=1,nlev
      !-Departure fvm Meshes, initialization done                                                               
      call compute_weights(fvm(ie),6,weights_all,weights_eul_index_all, &
             weights_lgr_index_all,k,jall)    
       
      ! THE FIRST TRACER IS AIRDENSITY
      tracer_air0=fvm(ie)%c(:,:,k,1,tl%n0)       
      call reconstruction(tracer_air0, fvm(ie),recons_air)

      call monotonic_gradient_cart(tracer_air0, fvm(ie),recons_air, elem(ie)%desc)
!       recons_air(1,:,:)=0.0D0
!       recons_air(2,:,:)=0.0D0
!       recons_air(3,:,:)=0.0D0
!       recons_air(4,:,:)=0.0D0
!       recons_air(5,:,:)=0.0D0
      
      tracer_air1=0.0D0   

      call cslam_remap(tracer_air0,tracer_air1,weights_all, recons_air, &
                 fvm(ie)%spherecentroid, weights_eul_index_all, weights_lgr_index_all, jall)             
      ! finish scheme
      do j=1,nc
        do i=1,nc
          tracer_air1(i,j)=tracer_air1(i,j)/fvm(ie)%area_sphere(i,j)
          fvm(ie)%c(i,j,k,1,tl%np1)=tracer_air1(i,j)
        end do
      end do
      !loop through all tracers
      do itr=2,ntrac
        tracer0=fvm(ie)%c(:,:,k,itr,tl%n0)
        call reconstruction(tracer0, fvm(ie),recons)
!         reconsalt=recons
        call monotonic_gradient_cart(tracer0, fvm(ie),recons, elem(ie)%desc)
!         recons(1,:,:)=0.0D0
!         recons(2,:,:)=0.0D0
!         recons(3,:,:)=0.0D0
!         recons(4,:,:)=0.0D0
!         recons(5,:,:)=0.0D0
        
        tracer1=0.0D0   
        
!         call cslam_remap(tracer0,tracer1,weights_all, recons, &
!                    fvm(ie)%spherecentroid, weights_eul_index_all, weights_lgr_index_all, jall)                   
        call cslam_remap_air(tracer0,tracer1,tracer_air0,weights_all, recons,recons_air,&
                       fvm(ie)%spherecentroid,weights_eul_index_all, weights_lgr_index_all, jall)
                       
        ! finish scheme
        do j=1,nc
          do i=1,nc
            fvm(ie)%c(i,j,k,itr,tl%np1)= &
            (tracer1(i,j)/fvm(ie)%area_sphere(i,j))/tracer_air1(i,j)
!             fvm(ie)%c(i,j,k,itr,tl%np1)=tracer1(i,j)/fvm(ie)%area_sphere(i,j)
          end do
        end do            
      enddo  !End Tracer
    end do  !End Level
    !note write tl%np1 in buffer
    call ghostVpack(cellghostbuf, fvm(ie)%c,nhc,nc,nlev,ntrac,0,tl%np1,timelevels,elem(ie)%desc)
  end do
  call t_stopf('CSLAM scheme')
  call t_startf('FVM Communication')
  call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc)
  call t_stopf('FVM Communication')
  !-----------------------------------------------------------------------------------!
  call t_startf('FVM Unpack')
  do ie=nets,nete
     call ghostVunpack(cellghostbuf, fvm(ie)%c, nhc, nc,nlev,ntrac, 0, tl%np1, timelevels,elem(ie)%desc)
  enddo
  call t_stopf('FVM Unpack')
end subroutine cslam_runairdensity


! ! use this subroutine for benchmark tests, couple airdensity with tracer concentration
! subroutine cslam_runair(elem,fvm,hybrid,deriv,tstep,tl,nets,nete)
!   ! ---------------------------------------------------------------------------------
!   use fvm_line_integrals_mod, only: compute_weights
!   ! ---------------------------------------------------------------------------------  
!   use fvm_filter_mod, only: monotonic_gradient_cart
!   ! ---------------------------------------------------------------------------------
!   use fvm_reconstruction_mod, only: reconstruction
!   ! ---------------------------------------------------------------------------------
!   use derivative_mod, only : derivative_t
!   ! ---------------------------------------------------------------------------------
!   use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
!   ! -----------------------------------------------
!    
!   implicit none
!   type (element_t), intent(inout)                :: elem(:)
!   type (fvm_struct), intent(inout)             :: fvm(:)
!   type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
!   integer, intent(in)                         :: nets  ! starting thread element number (private)
!   integer, intent(in)                         :: nete  ! ending thread element number   (private)
!   real (kind=real_kind)                       :: tstep
!   
!   integer                                     :: i,j,k,ie,itr, jx, jy, jdx, jdy, h
!   type (TimeLevel_t)                          :: tl              ! time level struct
!   type (derivative_t)                         :: deriv           ! derivative struct
!  
!   real (kind=real_kind)   , dimension(10*(nc+2*nhe)*(nc+2*nhe),6)  :: weights_all
!   integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_eul_index_all
!   integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_lgr_index_all
!   integer (kind=int_kind)                                          :: jall
!     
!   real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons
!   real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer0 
!   
!   real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer_air0   
!   real (kind=real_kind), dimension(1:nc,1:nc)                        :: tracer1, tracer_air1 
!   real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons_air   
!    
!   do ie=nets, nete
!     do k=1,nlev
!       call fvm_mesh_dep(elem(ie),deriv,fvm(ie),tstep,tl,k)
!       !-Departure fvm Meshes, initialization done                                                               
!       call compute_weights(fvm(ie),6,weights_all,weights_eul_index_all, &
!              weights_lgr_index_all,jall) 
!       !loop through all tracers
!       do itr=1,ntrac
!         tracer0=fvm(ie)%c(:,:,k,itr,tl%n0)
! !         call reconstruction(tracer0, fvm(ie),recons)
! !         call monotonic_gradient_cart(tracer0, fvm(ie),recons, elem(ie)%desc)
!         recons(1,:,:)=0.0D0
!         recons(2,:,:)=0.0D0
!         recons(3,:,:)=0.0D0
!         recons(4,:,:)=0.0D0
!         recons(5,:,:)=0.0D0
!         tracer1=0.0D0   
!         if (itr==1) then !calculation for air density  (the first tracer is supposed to be)
!           recons_air=recons
!           tracer_air0=tracer0
!           call cslam_remap(tracer0,tracer1,weights_all, recons, &
!                  fvm(ie)%spherecentroid, weights_eul_index_all, weights_lgr_index_all, jall)
!           ! finish scheme
!           do j=1,nc
!             do i=1,nc
!               tracer_air1(i,j)=tracer1(i,j)/fvm(ie)%area_sphere(i,j)
!               fvm(ie)%c(i,j,k,itr,tl%np1)=tracer_air1(i,j)
!             end do
!           end do
!         else !calculation for the other tracers    
!           call cslam_remap_air(tracer0,tracer1,tracer_air0,weights_all, recons,recons_air,&
!                        fvm(ie)%spherecentroid,weights_eul_index_all, weights_lgr_index_all, jall)                
!           ! finish scheme
!           do j=1,nc
!             do i=1,nc
!               fvm(ie)%c(i,j,k,itr,tl%np1)= &
!               (tracer1(i,j)/fvm(ie)%area_sphere(i,j))/tracer_air1(i,j)
!             end do
!           end do            
!         end if
!       enddo  !End Tracer
!     end do  !End Level
!     !note write tl%np1 in buffer
!     call ghostVpack(cellghostbuf, fvm(ie)%c,nhc,nc,nlev,ntrac,0,tl%np1,timelevels,elem(ie)%desc)
!   end do
!   call t_startf('FVM Communication')
!   call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc)
!   call t_stopf('FVM Communication')
!   !-----------------------------------------------------------------------------------!
!   do ie=nets,nete
!      call ghostVunpack(cellghostbuf, fvm(ie)%c, nhc, nc,nlev,ntrac, 0, tl%np1, timelevels,elem(ie)%desc)
!   enddo
! end subroutine cslam_runair


! do the remapping
subroutine cslam_remap(tracer0,tracer1,weights,recons,centroid, &
                 weights_eul_index_all, weights_lgr_index_all, jall)
  real (kind=real_kind), intent(in)           :: tracer0(1-nhc:nc+nhc,1-nhc:nc+nhc)
  real (kind=real_kind), intent(inout)        :: tracer1(1:nc,1:nc)
  real (kind=real_kind), intent(in)           :: weights(10*(nc+2*nhe)*(nc+2*nhe),6)
  real (kind=real_kind), intent(in)           :: recons(5,1-nhe:nc+nhe,1-nhe:nc+nhe)
  real (kind=real_kind), intent(in)           :: centroid(5,1-nhe:nc+nhe,1-nhe:nc+nhe)
  integer (kind=int_kind), intent(in)        :: weights_eul_index_all(10*(nc+2*nhe)*(nc+2*nhe),2)
  integer (kind=int_kind), intent(in)        :: weights_lgr_index_all(10*(nc+2*nhe)*(nc+2*nhe),2)
  integer (kind=int_kind), intent(in)        :: jall  
  
  integer                                     :: h, jx, jy, jdx, jdy
       
    do h=1,jall
      jx  = weights_lgr_index_all(h,1)
      jy  = weights_lgr_index_all(h,2)
      jdx = weights_eul_index_all(h,1)
      jdy = weights_eul_index_all(h,2)

      tracer1(jx,jy) = tracer1(jx,jy)+weights(h,1)*(&
         ! all constant terms 
         tracer0(jdx,jdy) - recons(1,jdx,jdy)*centroid(1,jdx,jdy) &
         - recons(2,jdx,jdy)*centroid(2,jdx,jdy) &
         + recons(3,jdx,jdy)*(2.0D0*centroid(1,jdx,jdy)**2 -centroid(3,jdx,jdy)) &
         + recons(4,jdx,jdy)*(2.0D0*centroid(2,jdx,jdy)**2 -centroid(4,jdx,jdy)) &
         + recons(5,jdx,jdy)*(2.0D0*centroid(1,jdx,jdy)*centroid(2,jdx,jdy)-centroid(5,jdx,jdy))) + &
         ! linear terms
         weights(h,2)*&
         (recons(1,jdx,jdy) - recons(3,jdx,jdy)*2.0D0*centroid(1,jdx,jdy) &
          - recons(5,jdx,jdy)*centroid(2,jdx,jdy)) + &
         weights(h,3)*&
         (recons(2,jdx,jdy) - recons(4,jdx,jdy)*2.0D0*centroid(2,jdx,jdy) &
         - recons(5,jdx,jdy)*centroid(1,jdx,jdy)) + &
         ! quadratic terms
         weights(h,4)*recons(3,jdx,jdy)+&
         weights(h,5)*recons(4,jdx,jdy)+&
         weights(h,6)*recons(5,jdx,jdy)
     end do
end subroutine cslam_remap

! do remapping with air (i.e. conserve mass of air density * concentration),
! see Nair et.al 2010 in JCP: A class of deformational flow test cases for linear transport
! schemes on the sphere, Appendix B
subroutine cslam_remap_air(tracer0,tracer1,tracer_air, weights,recons, recons_air, centroid, &
                 weights_eul_index_all, weights_lgr_index_all, jall)
                 
  real (kind=real_kind), intent(in)           :: tracer0(1-nhc:nc+nhc,1-nhc:nc+nhc)
  real (kind=real_kind), intent(inout)        :: tracer1(1:nc,1:nc)
  real (kind=real_kind), intent(in)           :: tracer_air(1-nhc:nc+nhc,1-nhc:nc+nhc)
  real (kind=real_kind), intent(in)           :: recons(5,1-nhe:nc+nhe,1-nhe:nc+nhe)
  real (kind=real_kind), intent(in)           :: recons_air(5,1-nhe:nc+nhe,1-nhe:nc+nhe)
  
  real (kind=real_kind), intent(in)           :: weights(10*(nc+2*nhe)*(nc+2*nhe),6)
  real (kind=real_kind), intent(in)           :: centroid(5,1-nhe:nc+nhe,1-nhe:nc+nhe)
  integer (kind=int_kind), intent(in)         :: weights_eul_index_all(10*(nc+2*nhe)*(nc+2*nhe),2)
  integer (kind=int_kind), intent(in)         :: weights_lgr_index_all(10*(nc+2*nhe)*(nc+2*nhe),2)
  integer (kind=int_kind), intent(in)         :: jall  

  integer                                     :: h, jx, jy, jdx, jdy    
  
  do h=1,jall
    jx  = weights_lgr_index_all(h,1)
    jy  = weights_lgr_index_all(h,2)
    jdx = weights_eul_index_all(h,1)
    jdy = weights_eul_index_all(h,2)                     
    tracer1(jx,jy) = tracer1(jx,jy)+&
          ! air density times tracer reconstruction
          tracer_air(jdx,jdy)*(weights(h,1)*(&      ! 1 is for air
          ! all constant terms 
          tracer0(jdx,jdy) - recons(1,jdx,jdy)*centroid(1,jdx,jdy) &
          - recons(2,jdx,jdy)*centroid(2,jdx,jdy) &
          + recons(3,jdx,jdy)*(2.0D0*centroid(1,jdx,jdy)**2 -centroid(3,jdx,jdy)) &
          + recons(4,jdx,jdy)*(2.0D0*centroid(2,jdx,jdy)**2 -centroid(4,jdx,jdy)) &
          + recons(5,jdx,jdy)*(2.0D0*centroid(1,jdx,jdy)*centroid(2,jdx,jdy)-centroid(5,jdx,jdy))) + &
          ! linear terms
          weights(h,2)* &
          (recons(1,jdx,jdy)- recons(3,jdx,jdy)*2.0D0*centroid(1,jdx,jdy) &
           - recons(5,jdx,jdy)*centroid(2,jdx,jdy)) + &
          weights(h,3)* &
          (recons(2,jdx,jdy) - recons(4,jdx,jdy)*2.0D0*centroid(2,jdx,jdy) &
          - recons(5,jdx,jdy)*centroid(1,jdx,jdy)) + &
          ! quadratic terms
          weights(h,4)*recons(3,jdx,jdy) + &
          weights(h,5)*recons(4,jdx,jdy) + &
          weights(h,6)*recons(5,jdx,jdy)) + &

          !tracer times air reconstruction
          tracer0(jdx,jdy)*(weights(h,1)*(&      
          ! all constant terms 
  !       tracer_air &  this term cancels it out
          - recons_air(1,jdx,jdy)*centroid(1,jdx,jdy) - recons_air(2,jdx,jdy)*centroid(2,jdx,jdy) &
          + recons_air(3,jdx,jdy)*(2.0D0*centroid(1,jdx,jdy)**2 -centroid(3,jdx,jdy)) &
          + recons_air(4,jdx,jdy)*(2.0D0*centroid(2,jdx,jdy)**2 -centroid(4,jdx,jdy)) &
          + recons_air(5,jdx,jdy)*(2.0D0*centroid(1,jdx,jdy)*centroid(2,jdx,jdy)-centroid(5,jdx,jdy))) + &
          ! linear terms
          weights(h,2)* &
          (recons_air(1,jdx,jdy) - recons_air(3,jdx,jdy)*2.0D0*centroid(1,jdx,jdy) &
          - recons_air(5,jdx,jdy)*centroid(2,jdx,jdy)) + &
          weights(h,3)* &
          (recons_air(2,jdx,jdy) - recons_air(4,jdx,jdy)*2.0D0*centroid(2,jdx,jdy) &
          - recons_air(5,jdx,jdy)*centroid(1,jdx,jdy)) + &
          ! quadratic terms
          weights(h,4)*recons_air(3,jdx,jdy)+&
          weights(h,5)*recons_air(4,jdx,jdy)+&
          weights(h,6)*recons_air(5,jdx,jdy))
  end do
end subroutine cslam_remap_air


! initialize global buffers shared by all threads
subroutine fvm_init1(par)
  use parallel_mod, only : parallel_t, haltmp
  type (parallel_t) :: par
  
  if (nc<4) then
     if (par%masterproc) then
        print *, "NUMBER OF CELLS ERROR for fvm: Number of cells parmeter"
        print *, "parameter nc at least 4 (nc>=4), nc*nc cells per element. This is"
        print *, "needed for the cubic reconstruction, which is only implemented yet! STOP"
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
     call haltmp("PARAMTER ERROR for fvm: ntrac > ntrac_d")
  endif

  call initghostbuffer(cellghostbuf,nlev,ntrac,nhc,nc) !+1 for the air_density, which comes from SE
  
  call initEdgebuffer(edgeveloc,2*nlev)
end subroutine fvm_init1



! initialization that can be done in threaded regions
subroutine fvm_init2(elem,fvm,hybrid,nets,nete,tl)
  use fvm_control_volume_mod, only: fvm_mesh_ari
  use fvm_analytic_mod, only: computexytosphere_moments
  use bndry_mod, only: compute_ghost_corner_orientation 


  type (timelevel_t) :: tl
  type (fvm_struct) :: fvm(:)
  type (element_t) :: elem(:)
  type (hybrid_t)                             :: hybrid
  integer :: ie,nets,nete

  call compute_ghost_corner_orientation(hybrid,elem,nets,nete)
  ! run some tests:
  !   call test_ghost(hybrid,elem,nets,nete)
  

  do ie=nets,nete
    call fvm_mesh_ari(elem(ie),fvm(ie),tl)
    call computexytosphere_moments(fvm(ie),elem(ie)%desc)
  enddo

end subroutine fvm_init2

! first communciation of FVM tracers
subroutine fvm_init3(elem,fvm,hybrid,nets,nete,tnp0)
  
  type (element_t),intent(inout)            :: elem(:)                 
  type (fvm_struct),intent(inout)           :: fvm(:)  
                  
  type (hybrid_t),intent(in)                :: hybrid                  
                                                                            
  integer,intent(in)                        :: nets,nete,tnp0       
                                                                            
  integer                                   :: ie                 
  
  ! do it only for FVM tracers, FIRST TRACER will be the AIR DENSITY   
  do ie=nets,nete 
    call ghostVpack(cellghostbuf, fvm(ie)%c,nhc,nc,nlev,ntrac,0,tnp0,timelevels,elem(ie)%desc)
  end do
  !exchange values for the initial data
  call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc)
  !-----------------------------------------------------------------------------------!
  do ie=nets,nete
    call ghostVunpack(cellghostbuf, fvm(ie)%c, nhc, nc,nlev,ntrac, 0, tnp0, timelevels,elem(ie)%desc)
  enddo

end subroutine fvm_init3
! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_MESH_DEP--------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 30.October 2011                                          !
! DESCRIPTION: Calculates the departure mesh                                        !
!              Please select the test cases by an COMMENT/UNCOMMENT basis           !
!                                                                                   !
! CALLS: solidbody, cart2cubedspherexy, analytical_function                         !
! INTPUT/OUTPUT:                                                                    !
!        fvm   ...  structure                                                       !
! INPUT: nstep   ... actual step                                                    !
!-----------------------------------------------------------------------------------!
subroutine fvm_mesh_dep(elem, deriv, fvm, dt, tl, klev)
  use coordinate_systems_mod, only : cartesian2D_t
  use element_mod, only : element_t
  use fvm_control_volume_mod, only: fvm_struct
  use time_mod, only : timelevel_t, time_at
  use parallel_mod, only : haltmp
  use control_mod, only : test_cfldep

  use derivative_mod, only : derivative_t

#ifndef _PRIM
  use fvm_bsp_mod, only: boomerang, solidbody
#endif

  implicit none
  type (derivative_t)  , intent(in) :: deriv
  type (fvm_struct), intent(inout)   :: fvm
  type (timelevel_t),intent(in)        :: tl
  integer,intent(in)                   :: klev
  
  type (element_t), intent(inout)      :: elem
  real (kind=real_kind)                :: time,dt, dx, dy, cflx, cfly, maxcflx, maxcfly  
  integer                              :: i,j
  type (cartesian2D_t)                 :: dcart(nc+1,nc+1)
  

! for the benchmark test, use more accurate departure point creation
#if 0
!CE: define new mesh for fvm fvm on an equal angular grid
! go from alpha,beta -> cartesian xy on cube -> lat lon on the sphere
! #ifdef _FVM
  do j=1,nc+1
     do i=1,nc+1               
!                 call solidbody(fvm%asphere(i,j), fvm%dsphere(i,j))
        call boomerang(fvm%asphere(i,j), fvm%dsphere(i,j,klev),tl%nstep)
     end do
  end do
! #endif
#else

! for given velocities in the element structure
  call fvm_dep_from_gll(elem, deriv, fvm%asphere,fvm%dsphere,dt,tl,klev)
#endif

   if (test_cfldep) then
     call check_departurecell(fvm,klev) 
   endif 
  
end subroutine fvm_mesh_dep

! subroutine fvm_dep_from_gll_iter(elem, deriv, fvm, dt, tl, klev)
!   use coordinate_systems_mod, only : cartesian2D_t
!   use element_mod, only : element_t
!   use fvm_control_volume_mod, only: fvm_struct
!   use time_mod, only : timelevel_t, time_at
!   use parallel_mod, only : haltmp
!   use control_mod, only : test_cfldep
! 
!   use derivative_mod, only : derivative_t
! 
!   implicit none
!   type (element_t), intent(in)          :: elem
!   type (derivative_t)  , intent(in) :: deriv
!   type (spherical_polar_t),intent(in)   :: asphere(nc+1,nc+1)
!   type (spherical_polar_t),intent(out)  :: dsphere(nc+1,nc+1)
!   type (timelevel_t),intent(in)        :: tl
!   real (kind=real_kind),intent(in)      :: dt
!   integer,intent(in)                   :: klev
!   
!   real (kind=real_kind)                 :: uxyz_gll(np,np,3),uxyz(nc+1,nc+1,3)
!   real (kind=real_kind)                 :: un,ue,ur,clon,clat,slon,slat
!   type(cartesian3D_t)                   :: acart
!   integer :: i,j,k
!   
!   real (kind=real_kind)  :: vstar(np,np,2)
!   
!    ! convert velocity from lat/lon to cartesian 3D
!    vstar = elem%derived%vstar(:,:,:,klev)
!     
!    do i=1,np
!      do j=1,np
!         clon = cos(elem%spherep(i,j)%lon)
!         slon = sin(elem%spherep(i,j)%lon)
!         clat = cos(elem%spherep(i,j)%lat)
!         slat = sin(elem%spherep(i,j)%lat)
! 
!         ur = 0
!         ue = vstar(i,j,1) 
!         un = vstar(i,j,2)
! 
!         uxyz_gll(i,j,1)= clon*clat*ur - clon*slat*un - slon*ue
!         uxyz_gll(i,j,2)= slon*clon*ur - slon*slat*un + clon*ue
!         uxyz_gll(i,j,3)= slat          *ur + clat          *un
!     
!      enddo
!    enddo
!    ! interpolate velocity to fvm nodes
!    do i=1,3
!       uxyz(:,:,i)=interpolate_gll2fvm_corners(uxyz_gll(:,:,i),deriv)
!    end do 
!    ! compute departure point 
!    ! crude, 1st order accurate approximation.  to be improved
!    do i=1,nc+1
!       do j=1,nc+1
!          ! note: asphere()%r=1, so we need to convert velocity to radians/sec:
!          acart = change_coordinates(asphere(i,j))  
!          acart%x = acart%x - dt*uxyz(i,j,1)/rearth
!          acart%y = acart%y - dt*uxyz(i,j,2)/rearth
!          acart%z = acart%z - dt*uxyz(i,j,3)/rearth
!          dsphere(i,j) = change_coordinates(acart)
!          dsphere(i,j)%r = asphere(i,j)%r
!       enddo
!    enddo
!   
! end subroutine fvm_dep_from_gll_iter

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
  use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t, change_coordinates
  use time_mod, only : timelevel_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, interpolate_gll2fvm_corners

  implicit none
  type (element_t), intent(in)          :: elem
  type (derivative_t)  , intent(in) :: deriv
  type (spherical_polar_t),intent(in)   :: asphere(nc+1,nc+1)
  type (spherical_polar_t),intent(out)  :: dsphere(-1:nc+3,-1:nc+3,nlev)
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
   ! interpolate velocity to fvm nodes
   do i=1,3
      uxyz(:,:,i)=interpolate_gll2fvm_corners(uxyz_gll(:,:,i),deriv)
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
         dsphere(i,j,klev) = change_coordinates(acart)
         dsphere(i,j,klev)%r = asphere(i,j)%r
      enddo
   enddo
  
end subroutine fvm_dep_from_gll
!END SUBROUTINE FVM_DEP_FROM_GLL------------------------------------------CE-for FVM!

!SUBROUTINE CHECK_DEPARTURECELL-------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: Check if the departure cell is degenerated or if the                 !
!              CFL number > 1.0D0                                                   !
!         IF THE CHECK FAILS, THE PROGRAM WILL STOP                                 !
!                                                                                   !
! INPUT:  fvm  ... fvm structur                                                 ! 
!         klev   ... Level (vertical)                                               !
!-----------------------------------------------------------------------------------!
subroutine check_departurecell(fvm,klev)
  use coordinate_systems_mod, only: cartesian2D_t, cart2cubedspherexy, spherical_to_cart
  use fvm_control_volume_mod, only:  fvm_struct                                         
  implicit none

  type (fvm_struct), intent(inout)   :: fvm
  integer,intent(in)                   :: klev
  
  type (cartesian2D_t)                 :: dcart(nc+1,nc+1)
  real (kind=real_kind)                :: cflx, cfly, maxcflx, maxcfly  
  integer                              :: i,j
  logical                              :: crossline
            
  ! calculate xy Cartesian on the cube of departure points on the corresponding face  
  ! the calculation of dcart is also done in fvm_line_integrals_mod.F90, one could 
  ! save the points here...
  ! BUT: this subroutine should only be used for test reasons
  do j=1,nc+1
     do i=1,nc+1               
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(i,j,klev)),&
             fvm%faceno,dcart(i,j))              
     end do
  end do                        
  cflx=0.0D0
  cfly=0.0D0
  fvm%maxcfl(1,klev)=cflx
  fvm%maxcfl(2,klev)=cfly
  crossline=.FALSE.
  do j=1,nc
     do i=1,nc  
      ! equidistant mesh in alpha/beta coordinates
      cflx=abs((atan(dcart(i,j)%x)-atan(fvm%acartx(i)))/fvm%dalpha)
      cfly=abs((atan(dcart(i,j)%y)-atan(fvm%acarty(j)))/fvm%dbeta)
      if(cflx>fvm%maxcfl(1,klev)) fvm%maxcfl(1,klev)=cflx
      if(cfly>fvm%maxcfl(2,klev)) fvm%maxcfl(2,klev)=cfly
    
      !one could stop immediately here, if crossline=.TRUE., but want to calculate 
      ! all CFL first
      call check_lines_cross(dcart(i,j),dcart(i+1,j),dcart(i,j+1),dcart(i+1,j+1),crossline)  
      call check_lines_cross(dcart(i,j),dcart(i,j+1),dcart(i+1,j),dcart(i+1,j+1),crossline)
     end do
  end do  
  ! nodes on the north and east boundary
  do i=1,nc+1  
    cflx=abs((atan(dcart(i,nc+1)%x)-atan(fvm%acartx(i)))/fvm%dalpha)
    cfly=abs((atan(dcart(nc+1,j)%y)-atan(fvm%acarty(j)))/fvm%dbeta)
    if(cflx>fvm%maxcfl(1,klev)) fvm%maxcfl(1,klev)=cflx
    if(cfly>fvm%maxcfl(2,klev)) fvm%maxcfl(2,klev)=cfly
  end do

  if (fvm%maxcfl(1,klev) > 1.0D0 .OR. fvm%maxcfl(2,klev) > 1.0D0) then
    write(*,*) "Error in fvm_mod.F90: CFL number too big"
    write(*,*) "CFL has to be < 1.0D0"
    write(*,*) "Choose a smaller time step!"
    write(*,*) "max CFL in this element: maxcflx", fvm%maxcfl(1,klev), "maxcfly", fvm%maxcfl(2,klev) 
    STOP "Exit program!"
  endif

  if (crossline) then
    write(*,*) "FATAL Error in fvm_mod.F90: departure cell is degenerated!"
    write(*,*) "Choose a smaller time step! max CFL in this element: maxcflx", fvm%maxcfl(1,klev), "maxcfly", fvm%maxcfl(2,klev) 
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

!SUBROUTINE CHECK_DEPARTURECELLEST----------------------------------------CE-for FVM!
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
    write(*,*) "FATAL Error in fvm_line_integrals.F90: CFL number too big"
    write(*,*) "CFL has to be < 1.0D0"
    write(*,*) "Choose a smaller time step!"
    write(*,*) "max CFL in this element: maxcflx", maxcflx, "maxcfly", maxcfly 
    STOP "Exit program!"
  endif
  
  if (crossline) then
    write(*,*) "FATAL Error in fvm_line_integrals.F90: departure cell is degenerated!"
    write(*,*) "Choose a smaller time step!"
    write(*,*) "max CFL in this element: maxcflx", maxcflx, "maxcfly", maxcfly 
    STOP "Exit program!"
  endif
end subroutine check_departurecellest


! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_MCGREGOR--------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 09. January 2012                                         !
! DESCRIPTION: ! using McGregor AMS 1993 scheme: Economical Determination of        !
!                Departure Points for Semi-Lagrangian Models                        !
!                                                                                   !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine fvm_mcgregor(elem, deriv, tstep, vhat, vstar,order)
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
end subroutine fvm_mcgregor
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
subroutine fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, tstep, ordertaylor)
  use derivative_mod, only : derivative_t, ugradv_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  
  implicit none

  type (element_t), intent(inout)                :: elem(:)
  type (fvm_struct), intent(in)                :: fvm(:)

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

end subroutine fvm_mcgregordss
!END SUBROUTINE FVM_MCGREGORDSS-------------------------------------------CE-for FVM!

end module fvm_mod
