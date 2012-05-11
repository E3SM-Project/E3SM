#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module dg_tvdrk_mod
!=======================================================================================================!
! ---------------------
 use kinds, only: real_kind
! ---------------------
 use physical_constants, only: rearth , g, dd_pi
! ---------------------
 use dimensions_mod, only: np,nlev
! ---------------------
 use element_mod, only: element_t
! ---------------------
 use edge_mod, only: EdgeBuffer_t, edgevpack,edgerotate,edgevunpack,edgeDGVpack,edgeDGVunpack
! ---------------------
 use filter_mod, only:  filter_t, filter_P
! ---------------------
 use hybrid_mod, only: hybrid_t
! ---------------------
 use reduction_mod, only: reductionbuffer_ordered_1d_t
! ---------------------
! use parallel_mod
! ---------------------
 use quadrature_mod, only: quadrature_t, gauss, gausslobatto
! ---------------------
 use derivative_mod, only: derivative_t, gradient_wk
! ---------------------
 use time_mod, only: timelevel_t, smooth
! ---------------------
 use control_mod, only: filter_freq, filter_counter, topology, test_case, nu
! ---------------------
 use cg_mod, only: cg_t
! ---------------------
 use bndry_mod, only: bndry_exchangev
!----------------------    
 use dg_tests_mod, only: sw2_init,sw5_init, galewsky_init 
!----------------------------------------------
 use dg_core_mod, only: dg_sw_model,contra2co,co2contra,height2phi,phi2height, &
                        dg_diff_grads_uv, dg_diff_grads, dg_diff_flux, dgsw_uvh_rhs , &
                        cov2sphere,sphere2cov,sphere2contra,contra2sphere 
!=======================================================================================================!
!=======================================================================================================!
 public :: dg_tvdrk
 public :: dg_uvhrk
 public :: dg_tvdrk_advect
 private:: dg_stage
 public :: dg_ssprk2     !testbed for diffusion LDG etc.
!=======================================================================================================!
 contains
!=======================================================================================================! 
subroutine dg_stage(stage,alpha_rk,beta_rk)    
!=======================================================================================================!
 use kinds, only : real_kind
!=======================================================================================================! 
 integer, intent(in):: stage
 real (kind=real_kind), dimension(stage,stage), intent(inout):: alpha_rk,beta_rk
 real (kind=real_kind) :: zero,one,two,three,four
!=======================================================================================================!
 zero = 0.0D0
 one  = 1.0D0
 two  = 2.0D0
 three= 3.0D0
 four = 4.0D0
!=======================================================================================================! 
if (stage==3) then
 alpha_rk= zero; beta_rk= zero
!=======================================================================================================! 
 alpha_rk(1,1)= one; beta_rk(1,1)= one
!=======================================================================================================!  
 alpha_rk(2,1)= three/four; alpha_rk(2,2)= one/four; beta_rk(2,2)= one/four 
!=======================================================================================================! 
 alpha_rk(3,1)= one/three; alpha_rk(3,3)= two/three; beta_rk(3,3)= two/three 
!=======================================================================================================!  
endif
!=======================================================================================================! 
end subroutine dg_stage
!=======================================================================================================!
! dg_tvdrk(stage,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
! 
! 	v(i,j,2,nlev,timelevel)	=> Contravariant velocity u^j(u^1,u^2)		
!	couv(i,j,2,nlev)	=> Covariant velocity u_i(u_1,u_2)
!	ht(i,j,nlev)		=> Height field (h)
!       hs(i,j,nlev)		=> Mountain height field (hs)
!	phi(i,j,nlev)		=> Geopotential Phi= g(h+hs)
!	psi(i,j,nlev)		=> Height field (h)
!  	fcor(i,j)	 	=> Coriolis Parameter
!
!=======================================================================================================!
subroutine dg_tvdrk(elem, stage,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
!=======================================================================================================!
!use dg_tests_mod, only: sw2_init,sw5_init
!----------------------------------------------
!use dg_core_mod, only: dg_sw_model,contra2co,co2contra,height2phi,phi2height 
!=======================================================================================================!
!   In and Out variable
!=======================================================================================================!
 implicit none
 type(element_t)      , intent(inout) :: elem(:)
 integer              , intent(in) :: stage
 type (EdgeBuffer_t)  , intent(in) :: edge3
 type (derivative_t)  , intent(in) :: deriv
 type (filter_t)                   :: flt
 type (hybrid_t)      , intent(in) :: hybrid
 type (quadrature_t)               :: gll
 real (kind=real_kind), intent(in) :: dt
 real (kind=real_kind), intent(in) :: pmean
 type (TimeLevel_t)   , intent(in) :: tl
 integer              , intent(in) :: nets
 integer              , intent(in) :: nete
!=======================================================================================================!
!   Local variable
!=======================================================================================================!
 integer:: ig, ntime
 real (kind=real_kind), dimension(:,:), pointer :: fcor,rmv,mv,metdet
 real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv
!=======================================================================================================!
 real (kind=real_kind), dimension(0:np+1,0:np+1,nlev,nets:nete)    :: htbuf
 real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: uvbuf
 real (kind=real_kind), dimension(np,np,3)                         :: sw_rhs
!=======================================================================================================! 
 real (kind=real_kind), dimension(np,np,nets:nete,0:stage)  :: ht_rk   
 real (kind=real_kind), dimension(np,np,2,nets:nete,0:stage):: uv_rk
 real (kind=real_kind), dimension(stage,stage):: alpha_rk,beta_rk
 real (kind=real_kind) :: zero,one,two,three,four
 integer :: s,im1,intrk
!=======================================================================================================!
! real (kind=real_kind) :: delx,dely
! real (kind=real_kind) :: dx,dy,rdx,rdy
 real (kind=real_kind) :: lenscale,grv 
 real (kind=real_kind) :: hdt 
 real*8                :: st,et,time_adv
 integer    :: i,j,k,ie
 integer    :: kptr
 integer    :: nm1,n0,np1
 integer    :: nstep
!=======================================================================================================!
!   Time level 		  ; Time step level ; Time step size ; Lenscale
!   nm1= 1, n0= 2, np1= 3 ; nstep= nstep    ; hdt= dt= tstep ; lenscale= rearth= 6.37122 x 10^6 m
!=======================================================================================================!
 nm1   = tl%nm1
 n0    = tl%n0
 np1   = tl%np1
 nstep = tl%nstep
!=======================================================================================================! 
 hdt = dt
 grv = g
 lenscale=rearth
!=======================================================================================================!   
 zero = 0.0D0
!=======================================================================================================!
 call dg_stage(stage,alpha_rk,beta_rk)
!=======================================================================================================!
!   Starting Time step
!=======================================================================================================!

 If (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0) then
    filter_counter = filter_counter + 1
    uvbuf = zero
!=======================================================================================================!
!   Starting Nlev
!=======================================================================================================!
 k= 1
!=======================================================================================================!
!   Initializing
!=======================================================================================================!
!=======================================================================================================!
   do ie=nets,nete
      if(nstep==1)then
        if (topology == "cube" .and. test_case=="swtc2") then
           call sw2_init(elem,tl,ie,k,pmean)
        elseif (topology == "cube" .and. test_case=="swtc5") then
           call sw5_init(elem,tl,ie,k,pmean)
        elseif (topology == "cube" .and. test_case=="galewsky") then
           call galewsky_init(elem,tl,ie,k,pmean)
        endif
        elem(ie)%state%couv(:,:,:,k)= contra2co(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%met(:,:,:,:))  
      endif
       !================================================================================================!
       ! Coverting covarinat to contravariant, height to phi, etc
       !================================================================================================!       
       elem(ie)%state%psi(:,:,k)= height2phi(elem(ie)%state%ht(:,:,k),elem(ie)%metdet(:,:))
       elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))       
       ht_rk(:,:,ie,0)  = elem(ie)%state%psi(:,:,k)
       uv_rk(:,:,:,ie,0)= elem(ie)%state%couv(:,:,:,k)
   enddo
!=======================================================================================================!
!   TVD-RK Methods
!=======================================================================================================!
do intrk=1,stage
!=======================================================================================================!
im1= intrk-1   
!=======================================================================================================!
   do ie=nets,nete
       !================================================================================================!
       ! DG: Packing Contravariant velocity & height fields
       !================================================================================================!
       kptr=0*nlev
       call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
       kptr=2*nlev
       call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)

       !================================================================================================!
       ! DG: Rotating velocity (contra)
       !================================================================================================!
       kptr=0*nlev
       call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)
   enddo
!=======================================================================================================!
   ! ===================================================================================================!
   ! Insert communications here: for shared memory, just a single
   ! sync is required
   ! ===================================================================================================!

    call bndry_exchangeV(hybrid,edge3)

!=======================================================================================================!
!   TVD-RK stage
!=======================================================================================================!    
!=======================================================================================================!
   do ie=nets,nete

       !================================================================================================!
       ! Unpack the edges for uvcomp and height
       !================================================================================================!
       kptr=0*nlev
       call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
       kptr =2*nlev
       call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)

       !================================================================================================!
       ! Boundary values for velocity & height fields from the neighbors,
       ! stored in the buffer
       !================================================================================================!
       call dg_sw_model(elem(ie),deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),          	     			&
            		elem(ie)%state%v(:,:,:,k,n0),elem(ie)%state%couv(:,:,:,k),   			&
	    		elem(ie)%state%ht(:,:,k),elem(ie)%state%hs(:,:,k),sw_rhs)
       
       ht_rk(:,:,ie,intrk)  = zero
       uv_rk(:,:,:,ie,intrk)= zero
       do s=0,im1
          ht_rk(:,:,ie,intrk)  = ht_rk(:,:,ie,intrk)+							&
	  			 alpha_rk(intrk,s+1)*ht_rk(:,:,ie,s) + beta_rk(intrk,s+1)*hdt*sw_rhs(:,:,1)
       	  uv_rk(:,:,1,ie,intrk)= uv_rk(:,:,1,ie,intrk)+							&
	  			 alpha_rk(intrk,s+1)*uv_rk(:,:,1,ie,s)+ beta_rk(intrk,s+1)*hdt*sw_rhs(:,:,2)       	  
	  uv_rk(:,:,2,ie,intrk)= uv_rk(:,:,2,ie,intrk)+							&
	  			 alpha_rk(intrk,s+1)*uv_rk(:,:,2,ie,s)+ beta_rk(intrk,s+1)*hdt*sw_rhs(:,:,3)				
       enddo       
       elem(ie)%state%ht(:,:,k)= phi2height(ht_rk(:,:,ie,intrk),elem(ie)%metdet(:,:))
       elem(ie)%state%v(:,:,:,k,n0)= co2contra(uv_rk(:,:,:,ie,intrk),elem(ie)%metinv(:,:,:,:))	
       elem(ie)%state%couv(:,:,:,k)= uv_rk(:,:,:,ie,intrk)
   enddo
!=======================================================================================================!
enddo
!=======================================================================================================!
!   Ending Nlev
!=======================================================================================================!
!=======================================================================================================!
!   Ending Time Step
!=======================================================================================================!
 ENDIF

!=======================================================================================================!
end subroutine dg_tvdrk
!=======================================================================================================!
!=======================================================================================================!
! dg_tvdrk_advect(stage,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
! 
! 	v(i,j,2,nlev,timelevel)	=> Contravariant velocity u^j(u^1,u^2)		
!	couv(i,j,2,nlev)	=> Covariant velocity u_i(u_1,u_2)
!	ht(i,j,nlev)		=> Height field (h)
!	psi(i,j,nlev)		=> Height field (h)
!
!=======================================================================================================!
subroutine dg_tvdrk_advect(elem,stage,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
!=======================================================================================================!
 use dg_tests_mod, only: sw1_init  
!----------------------------------------------
 use dg_core_mod, only: dg_adv_model,contra2co,co2contra,height2phi,phi2height
!=======================================================================================================!
!   In and Out variable
!=======================================================================================================!
 implicit none
 type(element_t)      , intent(inout) :: elem(:)
 integer              , intent(in) :: stage
 type (EdgeBuffer_t)  , intent(in) :: edge3
 type (derivative_t)  , intent(in) :: deriv
 type (filter_t)                   :: flt
 type (hybrid_t)      , intent(in) :: hybrid
 type (quadrature_t)               :: gll
 real (kind=real_kind), intent(in) :: dt
 real (kind=real_kind), intent(in) :: pmean
 type (TimeLevel_t)   , intent(in) :: tl
 integer              , intent(in) :: nets
 integer              , intent(in) :: nete
!=======================================================================================================!
!   Local variable
!=======================================================================================================!
 integer :: ig, ntime
 real (kind=real_kind), dimension(:,:), pointer :: fcor,rmv,mv,metdet
 real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv
!=======================================================================================================!
 real (kind=real_kind), dimension(0:np+1,0:np+1,nlev,nets:nete)    :: htbuf
 real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: uvbuf
 real (kind=real_kind), dimension(np,np)                           :: ht_rhs
!=======================================================================================================! 
 real (kind=real_kind), dimension(np,np,nets:nete,0:stage):: ht_rk
 real (kind=real_kind), dimension(stage,stage):: alpha_rk,beta_rk
 real (kind=real_kind) :: zero,one,two,three,four
 integer :: s,im1,intrk
!=======================================================================================================!
! real (kind=real_kind) :: delx,dely
! real (kind=real_kind) :: dx,dy,rdx,rdy
 real (kind=real_kind) :: lenscale,grv 
 real (kind=real_kind) :: hdt 
 real*8                :: st,et,time_adv
 integer :: i,j,k,ie
 integer :: kptr
 integer :: nm1,n0,np1,nstep
!=======================================================================================================!
!   Time level 		  ; Time step level ; Time step size ; Lenscale
!   nm1= 1, n0= 2, np1= 3 ; nstep= nstep    ; hdt= dt= tstep ; lenscale= rearth= 6.37122 x 10^6 m
!=======================================================================================================!
 nm1   = tl%nm1
 n0    = tl%n0
 np1   = tl%np1
 nstep = tl%nstep
!=======================================================================================================! 
 hdt = dt
 lenscale=rearth
!=======================================================================================================!   
 zero = 0.0D0
!=======================================================================================================!
 call dg_stage(stage,alpha_rk,beta_rk)
!=======================================================================================================!
!   Starting Time step
!=======================================================================================================!

 If (nstep > 0 ) then
    filter_counter = filter_counter + 1
    uvbuf = zero
!=======================================================================================================!
!   Starting Nlev
!=======================================================================================================!
 k=1
!=======================================================================================================!
!   Initializing
!=======================================================================================================!
!=======================================================================================================!
    do ie=nets,nete
       if(nstep==1)then
          call sw1_init(elem,tl,ie,k,pmean)
       endif

       elem(ie)%state%psi(:,:,k)   = height2phi(elem(ie)%state%ht(:,:,k),elem(ie)%metdet(:,:))
       elem(ie)%state%couv(:,:,:,k)= contra2co(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%met(:,:,:,:))
       ht_rk(:,:,ie,0)= elem(ie)%state%psi(:,:,k)
    enddo
!=======================================================================================================!
!   TVD-RK Methods
!=======================================================================================================!
do intrk=1,stage
!=======================================================================================================!
im1= intrk-1   
!=======================================================================================================!       
    do ie=nets,nete    
       !======================================================
       ! DG: Packing Contravariant velocity & height fields
       !======================================================
       kptr=0*nlev
       call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
       kptr=2*nlev
       call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)

       !======================================================
       ! DG: Rotating velocity (contra)
       !======================================================
       kptr=0*nlev
       call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)


    enddo
!=======================================================================================================!
   ! =============================================================
   ! Insert communications here: for shared memory, just a single
   ! sync is required
   ! =============================================================
    call bndry_exchangeV(hybrid,edge3)

!=======================================================================================================!
!   TVD-RK stage 
!=======================================================================================================!
    do ie=nets,nete

       !===========================================================
       ! Unpack the edges for uvcomp and height  (nair)
       !============================================================
       kptr=0*nlev
       call edgeDGVunpack(edge3,uvbuf(0,0,1,1,ie),2*nlev,kptr,elem(ie)%desc)
       kptr =2*nlev
       call edgeDGVunpack(edge3,htbuf(0,0,1,ie),nlev,kptr,elem(ie)%desc)

       call dg_adv_model(elem(ie),deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),     				&
       		         elem(ie)%state%v(:,:,:,k,n0),elem(ie)%state%ht(:,:,k),ht_rhs)	
       
       ht_rk(:,:,ie,intrk)= zero
       do s=0,im1
          ht_rk(:,:,ie,intrk) = ht_rk(:,:,ie,intrk)+							&
	  			alpha_rk(intrk,s+1)*ht_rk(:,:,ie,s) + beta_rk(intrk,s+1)*hdt*ht_rhs(:,:)
       enddo
       elem(ie)%state%ht(:,:,k)= phi2height(ht_rk(:,:,ie,intrk),elem(ie)%metdet(:,:))
    enddo
!=======================================================================================================!
enddo   
!=======================================================================================================!
!   Ending Time Step
!=======================================================================================================!
 ENDIF
!=======================================================================================================!
end subroutine dg_tvdrk_advect


!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!!------uvh- Expt. 12/6/2011 
subroutine dg_uvhrk(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!=======================================================================================================!
!   In and Out variable
!=======================================================================================================!
   implicit none
   type(element_t)      , intent(inout) :: elem(:)
   type (EdgeBuffer_t)  , intent(in) :: edge3
   type (derivative_t)  , intent(in) :: deriv
   type (filter_t)                   :: flt
   type (hybrid_t)      , intent(in) :: hybrid
   type (quadrature_t)               :: gll
   real (kind=real_kind), intent(in) :: dt
   real (kind=real_kind), intent(in) :: pmean
   type (TimeLevel_t)   , intent(in) :: tl
   integer              , intent(in) :: nets
   integer              , intent(in) :: nete
!=======================================================================================================!
!   Local variable
!=======================================================================================================!
   integer:: ig, ntime
   real (kind=real_kind), dimension(:,:), pointer :: fcor,rmv,mv,metdet
   real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv
!=======================================================================================================!
   real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: psi_rk0
   real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: uv_rk0
   real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: dif_uv, difgr1, difgr2
   real (kind=real_kind), dimension(0:np+1,0:np+1,nlev,nets:nete)    :: htbuf
   real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: uvbuf, dubuf, dvbuf
   real (kind=real_kind), dimension(np,np,2)                         :: gr1,gr2
   real (kind=real_kind), dimension(np,np,3)                         :: sw_rhs
   real (kind=real_kind), dimension(np,np)                           :: difu, difv
!=======================================================================================================!
  real (kind=real_kind) :: zero,one,two,three,four
   real (kind=real_kind) :: delx,dely, xlat,ylon
   real (kind=real_kind) :: dx,dy,rdx,rdy
   real (kind=real_kind) :: lenscale,grv
   real (kind=real_kind) :: hdt , r2d
   real*8                :: st,et,time_adv
   integer    :: i,j,k,ie
   integer    :: kptr
   integer    :: nm1,n0,np1
   integer    :: nstep
!=======================================================================================================!
!   Time level            ; Time step level ; Time step size ; Lenscale
!   nm1= 1, n0= 2, np1= 3 ; nstep= nstep    ; hdt= dt= tstep ; lenscale= rearth= 6.37122 x 10^6 m
!=======================================================================================================!
   nm1   = tl%nm1
   n0    = tl%n0
   np1   = tl%np1
   nstep = tl%nstep
!=======================================================================================================!
   hdt = dt
   grv = g
   lenscale=rearth
!=======================================================================================================!
!   Starting Time step
!=======================================================================================================!
If (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0) then
   filter_counter = filter_counter + 1
   uvbuf = zero
!=======================================================================================================!
!   Starting Nlev
!=======================================================================================================!
     k= 1         !Always
!=======================================================================================================!
!   Initializing
!=======================================================================================================!
  do ie=nets,nete
     if(nstep==1)then
       if (topology == "cube" .and. test_case=="swtc2") then
          call sw2_init(elem,tl,ie,k,pmean)
       elseif (topology == "cube" .and. test_case=="swtc5") then
          call sw5_init(elem,tl,ie,k,pmean)
      elseif (topology == "cube" .and. test_case=="galewsky") then
         call galewsky_init(elem,tl,ie,k,pmean)
       endif
       elem(ie)%state%uv(:,:,:,k) = elem(ie)%state%v(:,:,:,k,n0)
      !if (ie == nets)  print*, " dg_uvh activated  ", "tl = ", n0 
     endif
  !================================================================================================!

      elem(ie)%state%psi(:,:,k)= height2phi(elem(ie)%state%ht(:,:,k),elem(ie)%metdet(:,:))
      uv_rk0(:,:,:,k,ie) = elem(ie)%state%uv(:,:,:,k)
      elem(ie)%state%v(:,:,:,k,n0)= sphere2contra(elem(ie)%state%uv(:,:,:,k),elem(ie)%Dinv)

      ! DG: Packing (u,v) spherical  velocity & height fields
      !======================================================
      kptr=0*nlev
      call edgeDGVpack(edge3,elem(ie)%state%uv(:,:,:,:),2*nlev,kptr,elem(ie)%desc)
      kptr=2*nlev
      call edgeDGVpack(edge3,elem(ie)%state%ht(:,:,:),nlev,kptr,elem(ie)%desc)

  enddo

  ! =============================================================
      call bndry_exchangeV(hybrid,edge3)
  ! =============================================================

!++++
 
  if (nu .ne.  0.0D0 ) then      !activate diffusion 
!=======================================================================================================!
!   LDG-Step1
!=======================================================================================================!
  do ie=nets,nete

      !============================================================
      ! Unpack the edges for uvcomp and height  (nair)
      !============================================================
      kptr=0*nlev
      call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
      kptr =2*nlev
      call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)


      !====================================================
      ! stored in the buffer
      !====================================================

      call dg_diff_grads_uv(elem(ie),deriv,uvbuf(0,0,1,k,ie),elem(ie)%state%uv(1,1,1,k),gr1,gr2)

        difgr1(:,:,:,k,ie) = gr1(:,:,:)
        difgr2(:,:,:,k,ie) = gr2(:,:,:)

  end do

   do ie = nets, nete
       kptr=3*nlev
       call edgeDGVpack(edge3,difgr1(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

       kptr=5*nlev
       call edgeDGVpack(edge3,difgr2(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
   end do

      call bndry_exchangeV(hybrid,edge3)


 do ie=nets,nete
       !============================================================
       ! Unpack the edges for uvcomp and height  (nair)
       !============================================================
       kptr=3*nlev
       call edgeDGVunpack(edge3, dubuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
       kptr=5*nlev
       call edgeDGVunpack(edge3, dvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

       call dg_diff_flux(elem(ie),deriv,dubuf(0,0,1,1,ie),difgr1(1,1,1,1,ie),difu)

       call dg_diff_flux(elem(ie),deriv,dvbuf(0,0,1,1,ie),difgr2(1,1,1,1,ie),difv)

            !dif_uv(:,:,1,k,ie) = difu(:,:)
            !dif_uv(:,:,2,k,ie) = difv(:,:)

   ! Converting to [u,v] = A^-T [u_1,u_2] 

   do j=1,np
   do i=1,np
      dif_uv(i,j,1,k,ie) = difu(i,j) *elem(ie)%Dinv(1,1,i,j) + difv(i,j) * elem(ie)%Dinv(2,1,i,j) 
      dif_uv(i,j,2,k,ie) = difu(i,j) *elem(ie)%Dinv(1,2,i,j) + difv(i,j) * elem(ie)%Dinv(2,2,i,j) 
   enddo
   enddo

  end do

  else      ! No diffusion (LDG) 

         dif_uv = 0.0D0
!++++
  do ie=nets,nete
      !============================================================
      ! Unpack the edges for uvcomp and height  (nair)
      !============================================================
      kptr=0*nlev
      call edgeDGVunpack(edge3, uvbuf(:,:,:,:,ie), 2*nlev, kptr, elem(ie)%desc)
      kptr =2*nlev
      call edgeDGVunpack(edge3, htbuf(:,:,:,ie), nlev, kptr, elem(ie)%desc)
  end do

  endif 

!=======================================================================================================!
!   SSP-RK2,  stage 1
!=======================================================================================================!
  do ie=nets,nete

      call dgsw_uvh_rhs(elem(ie),deriv,uvbuf(:,:,:,k,ie),htbuf(0,0,k,ie),    &
                        elem(ie)%state%uv(:,:,:,k),                          &
                        elem(ie)%state%ht(:,:,k),elem(ie)%state%hs(:,:,k),sw_rhs)

       psi_rk0(:,:,k,ie)  = elem(ie)%state%psi(:,:,k) + dt* sw_rhs(:,:,1)
       elem(ie)%state%uv(:,:,1,k)= uv_rk0(:,:,1,k,ie) + dt* (sw_rhs(:,:,2) + dif_uv(:,:,1,k,ie))
       elem(ie)%state%uv(:,:,2,k)= uv_rk0(:,:,2,k,ie) + dt* (sw_rhs(:,:,3) + dif_uv(:,:,2,k,ie)) 
  end do


  do ie = nets, nete

     elem(ie)%state%ht(:,:,k)    = phi2height(psi_rk0(:,:,k,ie),elem(ie)%metdet(:,:))

      kptr=0*nlev
      call edgeDGVpack(edge3,elem(ie)%state%uv(:,:,:,k),2*nlev,kptr,elem(ie)%desc)
      kptr =2*nlev
      call edgeDGVpack(edge3,elem(ie)%state%ht(:,:,:),nlev,kptr,elem(ie)%desc)
  end do

  call bndry_exchangeV(hybrid,edge3)

!============================
!   SSP-RK2:  stage 2
!============================
   do ie = nets, nete

        kptr=0*nlev
        call edgeDGVunpack(edge3, uvbuf(:,:,:,:,ie), 2*nlev, kptr, elem(ie)%desc)
        kptr=2*nlev
        call edgeDGVunpack(edge3, htbuf(:,:,:,ie), nlev, kptr, elem(ie)%desc)

        call dgsw_uvh_rhs(elem(ie),deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),    &
                          elem(ie)%state%uv(:,:,:,k),                          &
                          elem(ie)%state%ht(:,:,k),elem(ie)%state%hs(:,:,k),sw_rhs)

        psi_rk0(:,:,k,ie)   = (elem(ie)%state%psi(:,:,k)  + psi_rk0(:,:,k,ie)         + dt* sw_rhs(:,:,1))*0.5D0
        elem(ie)%state%uv(:,:,1,k) = (uv_rk0(:,:,1,k,ie) + elem(ie)%state%uv(:,:,1,k) + dt* (sw_rhs(:,:,2)  &
                                                                                      + dif_uv(:,:,1,k,ie)))*0.5D0
        elem(ie)%state%uv(:,:,2,k) = (uv_rk0(:,:,2,k,ie) + elem(ie)%state%uv(:,:,2,k) + dt* (sw_rhs(:,:,3)  & 
                                                                                      + dif_uv(:,:,1,k,ie)))*0.5D0

        elem(ie)%state%ht(:,:,k) = phi2height(psi_rk0(:,:,k,ie),elem(ie)%metdet(:,:))

       elem(ie)%state%v(:,:,:,k,n0)= sphere2contra(elem(ie)%state%uv(:,:,:,k),elem(ie)%Dinv)

   end do

ENDIF

end subroutine dg_uvhrk 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!=======================================================================================================!
subroutine dg_ssprk2(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
!=======================================================================================================!
!   In and Out variable
!=======================================================================================================!
   implicit none
   type(element_t)      , intent(inout) :: elem(:)
   type (EdgeBuffer_t)  , intent(in) :: edge3
   type (derivative_t)  , intent(in) :: deriv
   type (filter_t)                   :: flt
   type (hybrid_t)      , intent(in) :: hybrid
   type (quadrature_t)               :: gll
   real (kind=real_kind), intent(in) :: dt
   real (kind=real_kind), intent(in) :: pmean
   type (TimeLevel_t)   , intent(in) :: tl
   integer              , intent(in) :: nets
   integer              , intent(in) :: nete
!=======================================================================================================!
!   Local variable
!=======================================================================================================!
   integer:: ig, ntime
   real (kind=real_kind), dimension(:,:), pointer :: fcor,rmv,mv,metdet
   real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv
!=======================================================================================================!
   real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: psi_rk
   real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: cuv_rk,cuv_tmp,contrauv 
   real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: dif_uv, difgr1, difgr2
   real (kind=real_kind), dimension(0:np+1,0:np+1,nlev,nets:nete)    :: htbuf
   real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: uvbuf, dubuf, dvbuf
   real (kind=real_kind), dimension(np,np,2)                         :: gr1,gr2
   real (kind=real_kind), dimension(np,np,3)                         :: sw_rhs
   real (kind=real_kind), dimension(np,np)                           :: difu, difv
!=======================================================================================================!
  real (kind=real_kind) :: zero,one,two,three,four
   real (kind=real_kind) :: delx,dely, xlat,ylon
   real (kind=real_kind) :: dx,dy,rdx,rdy
   real (kind=real_kind) :: lenscale,grv
   real (kind=real_kind) :: hdt , r2d
   real*8                :: st,et,time_adv
   integer    :: i,j,k,ie
   integer    :: kptr
   integer    :: nm1,n0,np1
   integer    :: nstep
!=======================================================================================================!
!   Time level            ; Time step level ; Time step size ; Lenscale
!   nm1= 1, n0= 2, np1= 3 ; nstep= nstep    ; hdt= dt= tstep ; lenscale= rearth= 6.37122 x 10^6 m
!=======================================================================================================!
   nm1   = tl%nm1
   n0    = tl%n0
   np1   = tl%np1
   nstep = tl%nstep
!=======================================================================================================!
   hdt = dt
   grv = g
   lenscale=rearth
!=======================================================================================================!
   zero = 0.0D0
   one  = 1.0D0
   two  = 2.0D0
   three= 3.0D0
   four = 4.0D0
   r2d = 180.0D0 / dd_pi
!=======================================================================================================!
!   Starting Time step
!=======================================================================================================!
If (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0) then
   filter_counter = filter_counter + 1
   uvbuf = zero
!=======================================================================================================!
!   Starting Nlev
!=======================================================================================================!
k= 1         !Always
!=======================================================================================================!
!   Initializing
!=======================================================================================================!
!=======================================================================================================!
  do ie=nets,nete
     if(nstep==1)then
       if (topology == "cube" .and. test_case=="swtc2") then
          call sw2_init(elem,tl,ie,k,pmean)
       elseif (topology == "cube" .and. test_case=="swtc5") then
          call sw5_init(elem,tl,ie,k,pmean)
       !elseif (topology == "cube" .and. test_case=="swtc8") then
        ! call sw8_init(elem,tl,ie,k,pmean)
       endif
      !elem(ie)%state%couv(:,:,:,k)= contra2co(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%met(:,:,:,:))
       elem(ie)%state%couv(:,:,:,k)= sphere2cov(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%D)
     endif
      !================================================================================================!
      ! Coverting covarinat to contravariant, height to phi, etc
      !================================================================================================!
      elem(ie)%state%psi(:,:,k)= height2phi(elem(ie)%state%ht(:,:,k),elem(ie)%metdet(:,:))
      cuv_tmp(:,:,:,k,ie)= elem(ie)%state%couv(:,:,:,k)
      elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))
      !contrauv(:,:,:,k,ie)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))

      !======================================================
      ! DG: Packing Contravariant velocity & height fields
      !======================================================
      kptr=0*nlev
      call edgeDGVpack(edge3,elem(ie)%state%v(1,1,1,1,n0),2*nlev,kptr,elem(ie)%desc)
      !call edgeDGVpack(edge3,contrauv(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
      kptr=2*nlev
      call edgeDGVpack(edge3,elem(ie)%state%ht(1,1,1),nlev,kptr,elem(ie)%desc)
      !======================================================
      ! DG: Rotating velocity (contra)
      !======================================================
      kptr=0*nlev
      call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)

  enddo

  ! =============================================================
  ! Insert communications here: for shared memory, just a single
  ! sync is required
  ! =============================================================

   call bndry_exchangeV(hybrid,edge3)

!=======================================================================================================!
!   LDG-Step1
!=======================================================================================================!
  do ie=nets,nete

      !============================================================
      ! Unpack the edges for uvcomp and height  (nair)
      !============================================================
      kptr=0*nlev
      call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
      kptr =2*nlev
      call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)


      !====================================================
      ! stored in the buffer
      !====================================================

      call dg_diff_grads(elem(ie),deriv,uvbuf(0,0,1,k,ie),elem(ie)%state%v(1,1,1,k,n0), &
    ! call dg_diff_grads(elem(ie),deriv,uvbuf(0,0,1,k,ie),contrauv(1,1,1,k,ie), &
                         elem(ie)%state%couv(1,1,1,k), gr1,gr2)

        difgr1(:,:,:,k,ie) = gr1(:,:,:)
        difgr2(:,:,:,k,ie) = gr2(:,:,:)
  end do

   do ie = nets, nete
       kptr=3*nlev
       call edgeDGVpack(edge3,difgr1(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
       kptr=3*nlev
       call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)
       kptr=5*nlev
       call edgeDGVpack(edge3,difgr2(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
       kptr=5*nlev
       call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)
   end do

      call bndry_exchangeV(hybrid,edge3)

       ! dif_uv = 0.0D0

   do ie=nets,nete
       !============================================================
       ! Unpack the edges for uvcomp and height  (nair)
       !============================================================
       kptr=3*nlev
       call edgeDGVunpack(edge3, dubuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
       kptr=5*nlev
       call edgeDGVunpack(edge3, dvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

       call dg_diff_flux(elem(ie),deriv,dubuf(0,0,1,1,ie),difgr1(1,1,1,1,ie),difu)

       call dg_diff_flux(elem(ie),deriv,dvbuf(0,0,1,1,ie),difgr2(1,1,1,1,ie),difv)

             dif_uv(:,:,1,k,ie) = difu(:,:)
             dif_uv(:,:,2,k,ie) = difv(:,:)

   end do

!=======================================================================================================!
!   SSP-RK2,  stage 1
!=======================================================================================================!
  do ie=nets,nete

      !====================================================
      ! Boundary values for velocity & height fields from the neighbors,
      ! stored in the buffer
      !====================================================
      !Note:- Reuse of the buffer values
      call dg_sw_model(elem(ie),deriv,uvbuf(0,0,1,k,ie),htbuf(0,0,k,ie),                                     &
                       elem(ie)%state%v(1,1,1,k,n0),elem(ie)%state%couv(1,1,1,k),                      &
                       !contrauv(1,1,1,k,ie),elem(ie)%state%couv(1,1,1,k),                      &
                       elem(ie)%state%ht(1,1,k),elem(ie)%state%hs(1,1,k),sw_rhs)

      psi_rk(:,:,k,ie)  = elem(ie)%state%psi(:,:,k)  + hdt*sw_rhs(:,:,1)
      cuv_rk(:,:,1,k,ie)= elem(ie)%state%couv(:,:,1,k)+ hdt* (sw_rhs(:,:,2) + dif_uv(:,:,1,k,ie))
      cuv_rk(:,:,2,k,ie)= elem(ie)%state%couv(:,:,2,k)+ hdt* (sw_rhs(:,:,3) + dif_uv(:,:,2,k,ie))
  end do


  do ie = nets, nete

      elem(ie)%state%couv(:,:,:,k)= cuv_rk(:,:,:,k,ie)
      elem(ie)%state%ht(:,:,k)    = phi2height(psi_rk(:,:,k,ie),elem(ie)%metdet(:,:))
      elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))
      !contrauv(:,:,:,k,ie)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))

      !======================================================
      ! DG: Packing  contravariant vectors,  ht-field  (nair)
      !======================================================
      kptr=0*nlev
       call edgeDGVpack(edge3,elem(ie)%state%v(1,1,1,1,n0),2*nlev,kptr,elem(ie)%desc)
      !call edgeDGVpack(edge3,contrauv(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
      kptr =2*nlev
      call edgeDGVpack(edge3,elem(ie)%state%ht(1,1,1),nlev,kptr,elem(ie)%desc)
      !======================================================
      ! DG: Rotating vectors (nair)
      !======================================================
      kptr=0*nlev
      call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)

  end do

  ! =============================================================
  ! Insert communications here: for shared memory, just a single
  ! sync is required
  ! =============================================================

  call bndry_exchangeV(hybrid,edge3)

!=======================================================================================================!
!   SSP-RK stage 2
!=======================================================================================================!
  do ie = nets, nete

        !===========================================================
        ! Unpack the edges for  flux vectors and scalar fields  (nair)
        !============================================================
        kptr=0*nlev
        call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
        kptr=2*nlev
        call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
        !====================================================
        !   DG: S-E-N-W  boundary  flux operations
        !====================================================
        Call dg_sw_model(elem(ie),deriv,uvbuf(0,0,1,k,ie),htbuf(0,0,k,ie),                           &
                         elem(ie)%state%v(1,1,1,k,n0),elem(ie)%state%couv(1,1,1,k),                  &
                         !contrauv(1,1,1,k,ie),elem(ie)%state%couv(1,1,1,k),                  &
                         elem(ie)%state%ht(1,1,k),elem(ie)%state%hs(1,1,k),sw_rhs)

        psi_rk(:,:,k,ie)  = (elem(ie)%state%psi(:,:,k)  + psi_rk(:,:,k,ie)  + hdt*sw_rhs(:,:,1))*0.5D0
        cuv_rk(:,:,1,k,ie)= (cuv_tmp(:,:,1,k,ie)+ cuv_rk(:,:,1,k,ie)+ hdt* (sw_rhs(:,:,2) + &
                                                                            dif_uv(:,:,1,k,ie)) )*0.5D0
        cuv_rk(:,:,2,k,ie)= (cuv_tmp(:,:,2,k,ie)+ cuv_rk(:,:,2,k,ie)+ hdt* (sw_rhs(:,:,3) + &
                                                                            dif_uv(:,:,2,k,ie)) )*0.5D0

        elem(ie)%state%couv(:,:,:,k)= cuv_rk(:,:,:,k,ie)
        elem(ie)%state%ht(:,:,k) = phi2height(psi_rk(:,:,k,ie),elem(ie)%metdet(:,:))

        elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))
       !contrauv(:,:,:,k,ie)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))
       !elem(ie)%state%v(:,:,:,k,n0)= cov2sphere(cuv_rk(:,:,:,k,ie),elem(ie)%Dinv)


   end do
!=======================================================================================================!
!   Ending Nlev
!=======================================================================================================!
!=======================================================================================================!
!   Ending Time Step
!=======================================================================================================!
ENDIF
!=======================================================================================================!
end subroutine dg_ssprk2
!=======================================================================================================!

!=======================================================================================================!
end module dg_tvdrk_mod
