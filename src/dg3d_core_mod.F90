#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!=======================================================================================!
module dg3d_core_mod
!=======================================================================================!
 use kinds,           only : real_kind
 use dimensions_mod,  only : np, nlev
 use derivative_mod,  only : derivative_t, gradient_wk
 use physical_constants, only : rearth, rrearth, g
 use dg_core_mod, only : sphere2cov
!=======================================================================================!
 implicit none
 private
!=======================================================================================! 
 public :: dg3d_gradient_mass   
 public :: dg3d_gradient_mom    
 public :: dg3d_source_term
 public :: dg3d_laxfred_flux 
 public :: dg3d_fluxjacobian 
 public :: check_fj 
 public :: real_vorticity 
!=======================================================================================!
 public :: dg3d_fnum_flux
 public :: dg3d_signalwave
!=======================================================================================!
 integer, private:: riemanntype= 0
!=======================================================================================!
 real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv,D,Dinv
 real (kind=real_kind), dimension(:,:),     pointer :: fcor,mvi,metdet,rmetdet
 real (kind=real_kind), dimension(np,np) :: mmx,mmxi
!=======================================================================================!
 real (kind=real_kind), parameter, private:: zero = 0.0D0
 real (kind=real_kind), parameter, private:: one  = 1.0D0  
 real (kind=real_kind), parameter, private:: scale= 1.0D0       
!=======================================================================================!
 contains
 !=======================================================================================!
!  Gradient Mass					   				!    
!  Dvv_twt ->  transpose of derivax *G-weights						!
!  Mvv_twt ->  Diadgonal matrix with G-weights						!
!=======================================================================================!
function  real_vorticity(uv,sg,D,deriv) result(vor)
!==========================================================================================
    type (derivative_t)         :: deriv
    real (kind=real_kind), intent(in) :: D(2,2,np,np)
    real (kind=real_kind), dimension(np,np), intent(in) :: sg
    real (kind=real_kind), dimension(np,np,2), intent(in) :: uv
    real (kind=real_kind), dimension(np,np,2) :: couv
    real (kind=real_kind), dimension(np,np) :: vor
    real(kind=real_kind) ::  dvdx00, dudy00
    real (kind=real_kind):: term,delm
    integer:: i,j,k,l

    couv(:,:,:) = sphere2cov(uv,D)

    do j=1,np
       do l=1,np
          dudy00=0.0D0
          dvdx00=0.0D0
          do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l)* couv(i,j,2)
             dudy00 = dudy00 + deriv%Dvv(i,l)* couv(j,i,1)
          enddo
          vor(l,j) = dvdx00
          deriv%vvtemp(j,l) = dudy00
       enddo
    enddo

    do j=1,np
       do i=1,np
           vor(i,j)=(vor(i,j)-deriv%vvtemp(i,j))*rrearth /sg(i,j)
       end do
    end do
end function  real_vorticity
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function dg3d_gradient_mass(deriv,uvflx)  result(gradf)
!=======================================================================================!
    type (derivative_t):: deriv
    real(kind=real_kind), intent(in) :: uvflx(np,np,2)
    real(kind=real_kind)             :: gradf(np,np)
    integer:: i,j,l
    logical, parameter :: UseUnroll = .TRUE.
    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11
!=======================================================================================!
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2
          sumx00=zero
          sumx01=zero
          sumx10=zero
          sumx11=zero

          sumy00=zero
          sumy01=zero
          sumy10=zero
          sumy11=zero

          do i=1,np
             sumx00  = sumx00  + deriv%Dvv_twt(i,l  ) * uvflx(i,  j,1)
             sumx01  = sumx01  + deriv%Dvv_twt(i,l+1) * uvflx(i,  j,1)
             sumx10  = sumx10  + deriv%Dvv_twt(i,l  ) * uvflx(i,j+1,1)
             sumx11  = sumx11  + deriv%Dvv_twt(i,l+1) * uvflx(i,j+1,1)

             sumy00  = sumy00  + deriv%Mvv_twt(i,  l) * uvflx(i,  j,2)
             sumy01  = sumy01  + deriv%Mvv_twt(i,l+1) * uvflx(i,  j,2)
             sumy10  = sumy10  + deriv%Mvv_twt(i,  l) * uvflx(i,j+1,2)
             sumy11  = sumy11  + deriv%Mvv_twt(i,l+1) * uvflx(i,j+1,2)
          end do

          deriv%vvtempt(j  ,l  ,1) = sumx00
          deriv%vvtempt(j  ,l+1,1) = sumx01
          deriv%vvtempt(j+1,l  ,1) = sumx10
          deriv%vvtempt(j+1,l+1,1) = sumx11

          deriv%vvtempt(j  ,l  ,2) = sumy00
          deriv%vvtempt(j  ,l+1,2) = sumy01
          deriv%vvtempt(j+1,l  ,2) = sumy10
          deriv%vvtempt(j+1,l+1,2) = sumy11

       end do
    end do


    do j=1,np,2
       do i=1,np,2
          sumx00=zero
          sumx01=zero
          sumx10=zero
          sumx11=zero

          sumy00=zero
          sumy01=zero
          sumy10=zero
          sumy11=zero

          do l=1,np
             sumx00 = sumx00 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i  ,1)
             sumx01 = sumx01 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i  ,1)
             sumx10 = sumx10 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i+1,1)
             sumx11 = sumx11 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i+1,1)

             sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i  ,2)
             sumy01 = sumy01 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i  ,2)
             sumy10 = sumy10 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i+1,2)
             sumy11 = sumy11 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i+1,2)
          end do

          gradf(i,j)    = (sumx00 + sumy00) * rrearth
          gradf(i,j+1)  = (sumx01 + sumy01) * rrearth
          gradf(i+1,j)  = (sumx10 + sumy10) * rrearth
          gradf(i+1,j+1)= (sumx11 + sumy11) * rrearth
       end do
    end do
else
    do j=1,np
       do l=1,np
          sumx00=zero
	  sumy00=zero
          do i=1,np
             sumx00  = sumx00  + deriv%Dvv_twt(i,l) * uvflx(i,j,1)
             sumy00  = sumy00  + deriv%Mvv_twt(i,l) * uvflx(i,j,2)
	  enddo
          deriv%vvtempt(j,l,1) = sumx00
          deriv%vvtempt(j,l,2) = sumy00
      enddo
   enddo
    do j=1,np
       do i=1,np
          sumx00=zero
	  sumy00=zero
          do l=1,np
             sumx00 = sumx00 +  deriv%Mvv_twt(l,j)*deriv%vvtempt(l,i,1)
             sumy00 = sumy00 +  deriv%Dvv_twt(l,j)*deriv%vvtempt(l,i,2)
	  enddo
          gradf(i,j) = (sumx00 + sumy00) * rrearth
       enddo
    enddo
endif

!=======================================================================================!
end function  dg3d_gradient_mass
!=======================================================================================!
!  Gradient Momentum					   				!    
!  Dvv_twt ->  transpose of derivax *G-weights						!
!  Mvv_twt ->  Diadgonal matrix with G-weights						!
!=======================================================================================!
subroutine dg3d_gradient_mom(deriv,energy,gradu1,gradu2)
!=======================================================================================!
    type (derivative_t) :: deriv
    real(kind=real_kind), intent(in) :: energy(np,np)
    real(kind=real_kind), intent(out):: gradu1(np,np), gradu2(np,np)
    integer:: i,j,l    
    logical, parameter :: UseUnroll = .TRUE.
    real(kind=real_kind):: sumx00,sumx01
    real(kind=real_kind):: sumy00,sumy01
    real(kind=real_kind):: sumx10,sumx11
    real(kind=real_kind):: sumy10,sumy11
!=======================================================================================!
    !Grad-u 
    !Grad-v 

if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2
          sumx00=zero
          sumx01=zero
          sumx10=zero
          sumx11=zero

          sumy00=zero
          sumy01=zero
          sumy10=zero
          sumy11=zero

          do i=1,np
             sumx00 = sumx00 + deriv%Dvv_twt(i,l)  * energy(i,j)
             sumx01 = sumx01 + deriv%Dvv_twt(i,l+1)* energy(i,j)
             sumx10 = sumx10 + deriv%Dvv_twt(i,l)  * energy(i,j+1)
             sumx11 = sumx11 + deriv%Dvv_twt(i,l+1)* energy(i,j+1)

             sumy00 = sumy00 + deriv%Mvv_twt(i,l)  * energy(i,j)
             sumy01 = sumy01 + deriv%Mvv_twt(i,l+1)* energy(i,j)
             sumy10 = sumy10 + deriv%Mvv_twt(i,l)  * energy(i,j+1)
             sumy11 = sumy11 + deriv%Mvv_twt(i,l+1)* energy(i,j+1)
          end do

          deriv%vvtempt(j  ,l  ,1) = sumx00
          deriv%vvtempt(j  ,l+1,1) = sumx01
          deriv%vvtempt(j+1,l  ,1) = sumx10
          deriv%vvtempt(j+1,l+1,1) = sumx11

          deriv%vvtempt(j  ,l  ,2) = sumy00
          deriv%vvtempt(j  ,l+1,2) = sumy01
          deriv%vvtempt(j+1,l  ,2) = sumy10
          deriv%vvtempt(j+1,l+1,2) = sumy11

       end do
    end do




    do j=1,np,2
       do i=1,np,2
          sumx00=zero
          sumx01=zero
          sumx10=zero
          sumx11=zero

          sumy00=zero
          sumy01=zero
          sumy10=zero
          sumy11=zero

          do l=1,np
             sumx00 = sumx00 +  deriv%Mvv_twt(l,j)*deriv%vvtempt(l,i,1)
             sumx01 = sumx01 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i,1)
             sumx10 = sumx10 +  deriv%Mvv_twt(l,j)*deriv%vvtempt(l,i+1,1)
             sumx11 = sumx11 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i+1,1)

             sumy00 = sumy00 +  deriv%Dvv_twt(l,j)*deriv%vvtempt(l,i,2)
             sumy01 = sumy01 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i,2)
             sumy10 = sumy10 +  deriv%Dvv_twt(l,j)*deriv%vvtempt(l,i+1,2)
             sumy11 = sumy11 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i+1,2)
          end do
          gradu1(i,j)    = rrearth*sumx00 
          gradu1(i,j+1)  = rrearth*sumx01 
          gradu1(i+1,j)  = rrearth*sumx10 
          gradu1(i+1,j+1)= rrearth*sumx11 

          gradu2(i,j)    = rrearth*sumy00
          gradu2(i,j+1)  = rrearth*sumy01
          gradu2(i+1,j)  = rrearth*sumy10
          gradu2(i+1,j+1)= rrearth*sumy11
       end do
    end do
else
    do j=1,np
       do l=1,np
          sumx00=zero
          sumy00=zero
          do i=1,np
             sumx00 = sumx00 + deriv%Dvv_twt(i,l) * energy(i,j)
             sumy00 = sumy00 + deriv%Mvv_twt(i,l) * energy(i,j)
	  enddo
          deriv%vvtempt(j,l,1) = sumx00
          deriv%vvtempt(j,l,2) = sumy00
      enddo
    enddo
    do j=1,np
       do i=1,np
          sumx00=zero
	  sumy00=zero
          do l=1,np
             sumx00 = sumx00 +  deriv%Mvv_twt(l,j)*deriv%vvtempt(l,i,1)
             sumy00 = sumy00 +  deriv%Dvv_twt(l,j)*deriv%vvtempt(l,i,2)
	  enddo
          gradu1(i,j) = rrearth*sumx00
          gradu2(i,j) = rrearth*sumy00
       enddo
    enddo
endif
!=======================================================================================!
end subroutine  dg3d_gradient_mom
!=======================================================================================================!
!=======================================================================================================!
subroutine dg3d_source_term(mv,rmv,deriv,gcori,contrauv,couv,pgrad,force,source)
!==========================================================================================
    Implicit None
    real (kind=real_kind),intent(in)    :: mv(np,np)          ! mass matrix on velocity grid
    real (kind=real_kind), intent(in)    :: rmv(np,np)         ! inverse mass matrix on velocity grid
    type (derivative_t):: deriv
    real (kind=real_kind), dimension(np,np), intent(in) :: gcori
    real (kind=real_kind), dimension(np,np,2), intent(in) :: contrauv, couv, pgrad  
    real (kind=real_kind), dimension(np,np,4), intent(in)::  force     
    real (kind=real_kind), dimension(np,np,4), intent(out):: source     
    real (kind=real_kind), dimension(np,np) :: vor, vort
    real(kind=real_kind) ::  dvdx00,dvdx01
    real(kind=real_kind) ::  dvdx10,dvdx11
    real(kind=real_kind) ::  dudy00,dudy01
    real(kind=real_kind) ::  dudy10,dudy11
    real (kind=real_kind):: term,delm
    integer:: i,j,k,l    
    logical, parameter :: UseUnroll = .TRUE.    
!=======================================================================================================!    
!  mv  => elem(ie)%mv   
!  rmv => elem(ie)%rmv
!=======================================================================================================!
  mmx(:,:) = mv(:,:)
  do j = 1, np
  do i = 1, np
      mmxi(i,j) = one / mmx(i,j) 
  enddo 
  enddo 
!=======================================================================================================!
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2

          dudy00=zero
          dudy01=zero
          dudy10=zero
          dudy11=zero

          dvdx00=zero
          dvdx01=zero
          dvdx10=zero
          dvdx11=zero

          do i=1,np

             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* couv(i  ,j,2)
             dvdx01 = dvdx01 + deriv%Dvv(i,l+1)* couv(i  ,j,2)
             dvdx10 = dvdx10 + deriv%Dvv(i,l  )* couv(i,j+1,2)
             dvdx11 = dvdx11 + deriv%Dvv(i,l+1)* couv(i,j+1,2)

             dudy00 = dudy00 + deriv%Dvv(i,l  )* couv(j  ,i,1)
             dudy01 = dudy01 + deriv%Dvv(i,l+1)* couv(j  ,i,1)
             dudy10 = dudy10 + deriv%Dvv(i,l  )* couv(j+1,i,1)
             dudy11 = dudy11 + deriv%Dvv(i,l+1)* couv(j+1,i,1)

          end do

          vor(l  ,j  ) = dvdx00
          vor(l+1,j  ) = dvdx01
          vor(l  ,j+1) = dvdx10
          vor(l+1,j+1) = dvdx11

          deriv%vvtemp(j  ,l  ) = dudy00
          deriv%vvtemp(j  ,l+1) = dudy01
          deriv%vvtemp(j+1,l  ) = dudy10
          deriv%vvtemp(j+1,l+1) = dudy11

       end do
    end do
else
    do j=1,np
       do l=1,np
          dudy00=zero
	  dvdx00=zero
          do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l)* couv(i,j,2)
             dudy00 = dudy00 + deriv%Dvv(i,l)* couv(j,i,1)
	  enddo
          vor(l,j) = dvdx00
          deriv%vvtemp(j,l) = dudy00
       enddo
    enddo
endif


    do j=1,np
       do i=1,np
          vor(i,j)=(vor(i,j)-deriv%vvtemp(i,j))*rrearth 
       end do
    end do

    source(:,:,:)= zero
    do j = 1, np
       do i = 1, np

          term =  (vor(i,j) +  gcori(i,j))

          source(i,j,1) = ( contrauv(i,j,2) * term - pgrad(i,j,1) + force(i,j,1) ) * mmx(i,j) 
          source(i,j,2) = (-contrauv(i,j,1) * term - pgrad(i,j,2) + force(i,j,2) ) * mmx(i,j)       
	  source(i,j,3) = 0.0D0
          source(i,j,4) = 0.0D0  
          !source(i,j,4) = force(i,j,4) * mmx(i,j) 
       enddo
    enddo
!=======================================================================================================! 
end subroutine dg3d_source_term
!=======================================================================================================! 
!=======================================================================================================!
!  Signal Wave								   				!
!=======================================================================================================!
subroutine dg3d_signalwave(contrauv,contrauv_halo,gh11,gh22,gh11_halo,gh22_halo,fjmax,lambda,SL,SR)
!=======================================================================================! 
 Implicit None 
 real (kind=real_kind),dimension(np,np,2),intent(in):: contrauv
 real (kind=real_kind), dimension(np,4,2),intent(in):: contrauv_halo   
 real (kind=real_kind),dimension(np,np),  intent(in):: gh11,gh22 
 real (kind=real_kind),dimension(np,4),   intent(in):: gh11_halo,gh22_halo
 real (kind=real_kind),dimension(4),intent(out)   :: fjmax  
 real (kind=real_kind),dimension(4,np),intent(out):: lambda,SL,SR
 integer, parameter:: south=1,east=2,north=3,west=4
 integer:: i,j,wall
!=======================================================================================================! 
!=======================================================================================================!
!=======================================================================================================!
    SL= zero
    SR= zero
    do i=1, np
       SL(south,i)= abs(contrauv(i,1,2))         + sqrt(gh22(i,1)) 
       SR(south,i)= abs(contrauv_halo(i,south,2))+ sqrt(gh22_halo(i,south))
    !=====================================================================================!       
       SL(north,i)= abs(contrauv_halo(i,north,2))+ sqrt(gh22_halo(i,north))
       SR(north,i)= abs(contrauv(i,np,2))        + sqrt(gh22(i,np))
    enddo   
    
    do j=1, np
       SL(west,j)= abs(contrauv_halo(j,west,1))+ sqrt(gh11_halo(j,west))
       SR(west,j)= abs(contrauv(1,j,1))        + sqrt(gh11(1,j))
    !=====================================================================================!  
       SL(east,j)= abs(contrauv(np,j,1))       + sqrt(gh11(np,j))
       SR(east,j)= abs(contrauv_halo(j,east,1))+ sqrt(gh11_halo(j,east))
    enddo   
!=======================================================================================!    
    lambda= zero
    do wall=1,4
    do i=1,np
       lambda(wall,i)= max( SL(wall,i), SR(wall,i) )
    enddo   
    enddo
    do wall=1,4    
    fjmax(wall)= zero
    do i=1,np 
       fjmax(wall)= max(lambda(wall,i),fjmax(wall))
    enddo
    enddo 
!=======================================================================================!    
    SR(:,:) = lambda(:,:)
    SL(:,:) =-lambda(:,:)
!=======================================================================================!
!=======================================================================================!
end subroutine dg3d_signalwave 
!=======================================================================================!
!=======================================================================================!
subroutine dg3d_fnum_flux(deriv,dx,dy,fjmax,lambda,SL,SR,vec,vec_halo,flx,flx_halo,numflux)
!=======================================================================================! 
 Implicit None
 type (derivative_t):: deriv
 real (kind=real_kind),intent(in):: dx,dy
 real (kind=real_kind),dimension(4),intent(in)   :: fjmax  
 real (kind=real_kind),dimension(4,np),intent(in):: lambda,SL,SR
 real (kind=real_kind),dimension(np,np,4),  intent(in):: vec 
 real (kind=real_kind),dimension(np,4,4),   intent(in):: vec_halo 
 real (kind=real_kind),dimension(np,np,2,4),intent(in):: flx
 real (kind=real_kind),dimension(np,4,2,4), intent(in):: flx_halo
 real (kind=real_kind),dimension(np,np,4), intent(out):: numflux
!=======================================================================================================! 
 integer, parameter:: south=1,east=2,north=3,west=4, neqn=4 
 integer:: i,j,k,eqn,wall 
 real (kind=real_kind),dimension(4,4,np):: fx_left,fx_right,fy_left,fy_right
 real (kind=real_kind),dimension(4,4,np):: vec_left,vec_right,flux_left,flux_right,flux_senw
 real (kind=real_kind),dimension(4,np):: nx,ny,del,edgemass
!=======================================================================================================!
 nx= zero
 ny= zero
 nx(south,:)= zero
 ny(south,:)=-one
 nx(east,:) = one
 ny(east,:) = zero      
 nx(north,:)= zero
 ny(north,:)= one     
 nx(west,:) =-one
 ny(west,:) = zero 
!=======================================================================================!
 del(south,:)= dx
 del(east,:) = dy
 del(north,:)= dx
 del(west,:) = dy    
!=======================================================================================!
 edgemass= zero
 do wall=1,4
 do i=1,np
   edgemass(wall,i) = deriv%Mvv_twt(i,i)
 enddo 
 enddo
!=======================================================================================!
 vec_left = zero
 vec_right= zero
 do eqn=1,neqn
   vec_left(eqn,south,:)= vec(:,1,eqn)
   vec_left(eqn,east,:) = vec(np,:,eqn)
   vec_left(eqn,north,:)= vec(:,np,eqn)
   vec_left(eqn,west,:) = vec(1,:,eqn)
!=======================================================================================!
   vec_right(eqn,south,:)= vec_halo(:,south,eqn) 
   vec_right(eqn,east,:) = vec_halo(:,east,eqn) 
   vec_right(eqn,north,:)= vec_halo(:,north,eqn)
   vec_right(eqn,west,:) = vec_halo(:,west,eqn)
 enddo 
!=======================================================================================!
 fx_left = zero
 fx_right= zero
 do eqn=1,neqn
   fx_left(eqn,south,:)= flx(:,1,1,eqn)
   fx_left(eqn,east,:) = flx(np,:,1,eqn)
   fx_left(eqn,north,:)= flx(:,np,1,eqn)
   fx_left(eqn,west,:) = flx(1,:,1,eqn)
!=======================================================================================!
   fx_right(eqn,south,:)= flx_halo(:,south,1,eqn)
   fx_right(eqn,east,:) = flx_halo(:,east,1,eqn) 
   fx_right(eqn,north,:)= flx_halo(:,north,1,eqn)
   fx_right(eqn,west,:) = flx_halo(:,west,1,eqn)
 enddo  
!=======================================================================================!
 fy_left = zero
 fy_right= zero
 do eqn=1,neqn
   fy_left(eqn,south,:)= flx(:,1,2,eqn)
   fy_left(eqn,east,:) = flx(np,:,2,eqn)
   fy_left(eqn,north,:)= flx(:,np,2,eqn)
   fy_left(eqn,west,:) = flx(1,:,2,eqn)
!=======================================================================================!
   fy_right(eqn,south,:)= flx_halo(:,south,2,eqn)
   fy_right(eqn,east,:) = flx_halo(:,east,2,eqn) 
   fy_right(eqn,north,:)= flx_halo(:,north,2,eqn)
   fy_right(eqn,west,:) = flx_halo(:,west,2,eqn)
 enddo 
!=======================================================================================! 
 flux_left = zero
 flux_right= zero
 do eqn=1,neqn
 do wall=1,4
   flux_left(eqn,wall,:) =  fx_left(eqn,wall,:)*nx(wall,:)+fy_left(eqn,wall,:)*ny(wall,:)
   flux_right(eqn,wall,:)= fx_right(eqn,wall,:)*nx(wall,:)+fy_right(eqn,wall,:)*ny(wall,:)
 enddo 
 enddo  
!=======================================================================================!
!  For SW-system: numerical flux (phi,u1,u2) order (3 equations)			!
!=======================================================================================!   
 flux_senw(:,:,:)= zero
 do eqn =1,neqn
 do wall=1,4
 do j=1,np
 if (riemanntype==0) then
   flux_senw(eqn,wall,j)= 0.5D0*(flux_left(eqn,wall,j)+flux_right(eqn,wall,j) 		&
   			 -scale*fjmax(wall)*(vec_right(eqn,wall,j) - vec_left(eqn,wall,j))) 
 elseif (riemanntype==1) then
   flux_senw(eqn,wall,j)= 0.5D0*( flux_left(eqn,wall,j)+flux_right(eqn,wall,j) 		&
   			 -scale*lambda(wall,j)*(vec_right(eqn,wall,j) - vec_left(eqn,wall,j)))
 endif			 
 enddo
 enddo 
 enddo
!=======================================================================================!
!  For SW-system: numerical flux (phi,u1,u2) order (3 equations)			!
!=======================================================================================!   
 numflux(:,:,:)= zero
 do eqn=1,neqn
 do j=1,np
       numflux(j,1,eqn) = numflux(j,1,eqn) +edgemass(south,j)*flux_senw(eqn,south,j)*del(south,j)       
       numflux(np,j,eqn)= numflux(np,j,eqn)+edgemass(east,j)*flux_senw(eqn,east,j)*del(east,j)
       numflux(j,np,eqn)= numflux(j,np,eqn)+edgemass(north,j)*flux_senw(eqn,north,j)*del(north,j)       
       numflux(1,j,eqn) = numflux(1,j,eqn) +edgemass(west,j)*flux_senw(eqn,west,j)*del(west,j)
 enddo
 enddo
!=======================================================================================!
!=======================================================================================!
end subroutine dg3d_fnum_flux 
!=======================================================================================!    
!=======================================================================================!
subroutine dg3d_laxfred_flux(numeqn,deriv,fjmax,si,si_senw,uvflx,uvflx_senw,fluxout)

   integer, parameter :: south=1, east=2, north=3, west=4
   integer, intent(in):: numeqn
   type (derivative_t)                                  :: deriv
   real (kind=real_kind), dimension(4),    intent(in)   :: fjmax
   real (kind=real_kind), dimension(np,np,2,numeqn),intent(in) :: uvflx
   real (kind=real_kind), dimension(np,np,numeqn),  intent(in) :: si
   real (kind=real_kind), dimension(np,4,2,numeqn), intent(in) :: uvflx_senw
   real (kind=real_kind), dimension(np,4,numeqn),   intent(in) :: si_senw

   real (kind=real_kind), dimension(np,np,numeqn), intent(out) :: fluxout

   real (kind=real_kind), dimension(np,np) :: mij
   real (kind=real_kind), dimension(np)   :: lf_south,lf_north,lf_east,lf_west
   real(kind=real_kind) ::  fj(4), ul,ur , left, right, f_left, f_right, s1,s2

   integer i,j,k,eqn

      mij(:,:) = 0.0D0
      mij(1,1) = 1.0D0
    mij(np,np) = 1.0D0

        do j = 1, 4
           fj(j) = fjmax(j) *1.00D0
        enddo

    !For 3D-SW-system  (u1,u2,dp,pt,qt) order   (5 equations)

    fluxout(:,:,:) = 0.0D0

      do eqn = 1, numeqn

           ! East & West   LF flux  (fjmax <- max of flux Jacobian)

         do j = 1, np
                  f_left  = uvflx_senw(j,west,1,eqn)
                  f_right = uvflx(1,j,1,eqn)
                    left  = si_senw(j,west,eqn)
                   right  = si(1,j,eqn)
               lf_west(j) = 0.5D0 *(f_left + f_right - fj(west)*(right - left))

                  f_left  = uvflx(np,j,1,eqn)
                  f_right = uvflx_senw(j,east,1,eqn)
                    left  = si(np,j,eqn)
                   right  = si_senw(j,east,eqn)
               lf_east(j) = 0.5D0 *(f_left + f_right - fj(east)*(right - left))
          end do

           ! North& South  LF flux  (fjmax <- max of flux Jacobian)

         do i = 1, np
                   f_left = uvflx_senw(i,south,2,eqn)
                  f_right = uvflx(i,1,2,eqn)
                     left = si_senw(i,south,eqn)
                   right  = si(i,1,eqn)
              lf_south(i) = 0.5D0 *(f_left + f_right - fj(south)*(right - left))

                  f_left  = uvflx(i,np,2,eqn)
                  f_right = uvflx_senw(i,north,2,eqn)
                    left  = si(i,np,eqn)
                    right = si_senw(i,north,eqn)
              lf_north(i) = 0.5D0 *(f_left + f_right - fj(north)*(right - left))
         end do

        !Flux integral along the element boundary

        do j = 1, np
        do i = 1, np
            s1 = (lf_east(j) *mij(i,np) - lf_west(j) *mij(i,1) )* deriv%Mvv_twt(j,j)
            s2 = (lf_north(i)*mij(j,np) - lf_south(i)*mij(j,1) )* deriv%Mvv_twt(i,i)
            fluxout(i,j,eqn) = (s1  + s2)*rrearth
        end do
        end do

    end do

 end subroutine  dg3d_laxfred_flux

!=======================================================================================!   
function dg3d_fluxjacobian(contrauv,contrauv_halo,gh11,gh22,gh11_halo,gh22_halo) result(fjmax)
!=======================================================================================!
 Implicit None
 real (kind=real_kind),dimension(np,np,2),intent(in):: contrauv
 real (kind=real_kind), dimension(np,4,2),intent(in):: contrauv_halo
 real (kind=real_kind),dimension(np,np),  intent(in):: gh11,gh22
 real (kind=real_kind),dimension(np,4),   intent(in):: gh11_halo,gh22_halo
 real (kind=real_kind),dimension(4)   :: fjmax
 integer, parameter:: south=1,east=2,north=3,west=4

 integer:: i,j,wall
 real (kind=real_kind):: alfa1,alfa2, ul,ur
!=======================================================================================================!
! Max of the eigenvalues of the flux-jacobian from "left" and "right" edges
!=======================================================================================================!

     alfa1 = 0.0D0
     alfa2 = 0.0D0
   do i = 1, np
        ul = abs(contrauv(i,1,2))          + sqrt(gh22(i,1))
        ur = abs(contrauv_halo(i,south,2)) + sqrt(gh22_halo(i,south))
      alfa1= max(alfa1,ul,ur)
 
        ul = abs(contrauv_halo(i,north,2)) + sqrt(gh22_halo(i,north))
        ur = abs(contrauv(i,np,2))         + sqrt(gh22(i,np))
     alfa2= max(alfa2,ul,ur)
   enddo

     fjmax(south) = alfa1
     fjmax(north) = alfa2

     alfa1 = 0.0D0
     alfa2 = 0.0D0
   do j = 1, np
        ul = abs(contrauv_halo(j,west,1))  + sqrt(gh11_halo(j,west))
        ur = abs(contrauv(1,j,1))          + sqrt(gh11(1,j))
      alfa1= max(alfa1,ul,ur)

        ul = abs(contrauv(np,j,1))         + sqrt(gh11(np,j))
        ur = abs(contrauv_halo(j,east,1))  + sqrt(gh11_halo(j,east))
      alfa2= max(alfa2,ul,ur)
   enddo

     fjmax(west) = alfa1
     fjmax(east) = alfa2

!=======================================================================================================!
 end function  dg3d_fluxjacobian
!=======================================================================================================!
function check_fj(contrauv,contrauv_halo,gh11,gh22,gh11_halo,gh22_halo) result(fjmax)
!=======================================================================================!
 Implicit None
 real (kind=real_kind),dimension(np,np,2),intent(in):: contrauv
 real (kind=real_kind), dimension(np,4,2),intent(in):: contrauv_halo
 real (kind=real_kind),dimension(np,np),  intent(in):: gh11,gh22
 real (kind=real_kind),dimension(np,4),   intent(in):: gh11_halo,gh22_halo
 real (kind=real_kind),dimension(4)   :: fjmax
 integer, parameter:: south=1,east=2,north=3,west=4

 integer:: i,j,wall
 real (kind=real_kind):: alfa1,alfa2, ul,ur, sr
!=======================================================================================================!
! Max of the eigenvalues of the flux-jacobian from "left" and "right" edges
!=======================================================================================================!

     alfa1 = 0.0D0
     alfa2 = 0.0D0
   do i = 1, np
        sr = sqrt(gh22(i,1))
        ul = abs(contrauv(i,1,2))          + sr
        ur = abs(contrauv_halo(i,south,2)) + sr
      alfa1= max(alfa1,ul,ur)
 
        sr = sqrt(gh22(i,np))
        ul = abs(contrauv_halo(i,north,2)) + sr
        ur = abs(contrauv(i,np,2))         + sr
     alfa2= max(alfa2,ul,ur)
   enddo

     fjmax(south) = alfa1
     fjmax(north) = alfa2

     alfa1 = 0.0D0
     alfa2 = 0.0D0
   do j = 1, np
        sr = sqrt(gh11(1,j))
        ul = abs(contrauv_halo(j,west,1))  + sr
        ur = abs(contrauv(1,j,1))          + sr
      alfa1= max(alfa1,ul,ur)

        sr = sqrt(gh11(np,j))
        ul = abs(contrauv(np,j,1))         + sr
        ur = abs(contrauv_halo(j,east,1))  + sr
      alfa2= max(alfa2,ul,ur)
   enddo

     fjmax(west) = alfa1
     fjmax(east) = alfa2

!=======================================================================================================!
 end function  check_fj
!=======================================================================================================!
!=======================================================================================================!
!=======================================================================================================!
end module dg3d_core_mod


