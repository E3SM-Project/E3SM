#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dg_flux_mod
!=======================================================================================!
 use kinds,           only : real_kind
 use dimensions_mod,  only : np, nlev
 use derivative_mod,  only : derivative_t
!=======================================================================================!
 implicit none
 private
 public :: fnum_flux 
 public :: signalwave
!=======================================================================================================!
 integer, public :: riemanntype
 real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv
 real (kind=real_kind), dimension(:,:),     pointer :: rmv,mv,mvi,metdet,rmetdet
 real (kind=real_kind), public, dimension(np,np) :: mmxi            
!=======================================================================================!
 contains
!=======================================================================================================!
!  Shallow Water Flux Evaluation					   				!
!=======================================================================================================!  
!  Primitive Variable Terms: U=[sqrt(det(G(i,j)))*h, u_1, u_2]						!
!		si(i,j,1)= phi(i,j)									!
!		si(i,j,2)= couv(i,j,1)									!
!		si(i,j,3)= couv(i,j,2)									!
!=======================================================================================================!
!  Flux1 term=> fxy(i,j,1,1)= sqrt(det(G(i,j)))*(ht-hill)*u^1(contra)					!
!		fxy(i,j,1,2)= energy(i,j)								!
!		fxy(i,j,1,3)= 0										!
!=======================================================================================================!
!  Flux2 term=>	fxy(i,j,2,1)= sqrt(det(G(i,j)))*(ht-hill)*u^2(contra)					!
!		fxy(i,j,2,2)= 0										!
!		fxy(i,j,2,3)= energy(i,j)								!
!=======================================================================================================!
subroutine fnum_flux(deriv,dx,dy,fjmax,lambda,SL,SR,vec,vec_halo,flx,flx_halo,numflux)
!=======================================================================================! 
 Implicit None 
 type (derivative_t), intent(in) :: deriv
 real (kind=real_kind),intent(in):: dx,dy
 real (kind=real_kind),dimension(4),intent(in)   :: fjmax  
 real (kind=real_kind),dimension(4,np),intent(in):: lambda,SL,SR
 real (kind=real_kind),dimension(np,np,3),  intent(in):: vec 
 real (kind=real_kind),dimension(np,4,3),   intent(in):: vec_halo 
 real (kind=real_kind),dimension(np,np,2,3),intent(in):: flx
 real (kind=real_kind),dimension(np,4,2,3), intent(in):: flx_halo
 real (kind=real_kind),dimension(np,np,3), intent(out):: numflux 
 integer, parameter:: south=1,east=2,north=3,west=4
 integer:: i,j,k,eqn,wall 
 real (kind=real_kind),dimension(3,4,np):: fx_left,fx_right,fy_left,fy_right
 real (kind=real_kind),dimension(3,4,np):: vec_left,vec_right,flux_left,flux_right,flux_senw
 real (kind=real_kind),dimension(4,np):: nx,ny,del,edgemass
!=======================================================================================================!
 nx= 0.0D0
 ny= 0.0D0
 nx(south,:)= 0.0D0
 ny(south,:)=-1.0D0
 nx(east,:) = 1.0D0
 ny(east,:) = 0.0D0      
 nx(north,:)= 0.0D0
 ny(north,:)= 1.0D0     
 nx(west,:) =-1.0D0
 ny(west,:) = 0.0D0 
!=======================================================================================!
 del(south,:)= dx
 del(east,:) = dy
 del(north,:)= dx
 del(west,:) = dy    
!=======================================================================================!
 edgemass= 0.0D0
 do wall=1,4
 do i=1,np
   edgemass(wall,i) = deriv%Mvv_twt(i,i)
 enddo 
 enddo
!=======================================================================================!
 vec_left = 0.0D0
 vec_right= 0.0D0
 do eqn=1,3
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
 fx_left = 0.0D0
 fx_right= 0.0D0
 do eqn=1,3
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
 fy_left = 0.0D0
 fy_right= 0.0D0
 do eqn=1,3
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
 flux_left = 0.0D0
 flux_right= 0.0D0
 do eqn=1,3
 do wall=1,4
   flux_left(eqn,wall,:) =  fx_left(eqn,wall,:)*nx(wall,:)+fy_left(eqn,wall,:)*ny(wall,:)
   flux_right(eqn,wall,:)= fx_right(eqn,wall,:)*nx(wall,:)+fy_right(eqn,wall,:)*ny(wall,:)
 enddo 
 enddo  
!=======================================================================================!
!  For SW-system: numerical flux (phi,u1,u2) order (3 equations)			!
!=======================================================================================!   
 flux_senw(:,:,:)= 0.0D0
 do eqn =1,3
 do wall=1,4
 do j=1,np
 if (riemanntype==0) then
   flux_senw(eqn,wall,j)= 0.5D0*(flux_left(eqn,wall,j)+flux_right(eqn,wall,j) 		&
   			 -fjmax(wall)*(vec_right(eqn,wall,j) - vec_left(eqn,wall,j))) 
 elseif (riemanntype==1) then
   flux_senw(eqn,wall,j)= 0.5D0*( flux_left(eqn,wall,j)+flux_right(eqn,wall,j) 		&
   			 -lambda(wall,j)*(vec_right(eqn,wall,j) - vec_left(eqn,wall,j)))
 endif			 
 enddo
 enddo 
 enddo
!=======================================================================================!
!  For SW-system: numerical flux (phi,u1,u2) order (3 equations)			!
!=======================================================================================!   
 numflux(:,:,:)= 0.0D0
 do eqn = 1, 3
 do j=1, np
       numflux(j,1,eqn) = numflux(j,1,eqn) +edgemass(south,j)*flux_senw(eqn,south,j)*del(south,j)       
       numflux(np,j,eqn)= numflux(np,j,eqn)+edgemass(east,j)*flux_senw(eqn,east,j)*del(east,j)
       numflux(j,np,eqn)= numflux(j,np,eqn)+edgemass(north,j)*flux_senw(eqn,north,j)*del(north,j)       
       numflux(1,j,eqn) = numflux(1,j,eqn) +edgemass(west,j)*flux_senw(eqn,west,j)*del(west,j)
 enddo
 enddo
!=======================================================================================!
 return
!=======================================================================================!
end subroutine fnum_flux 

!=======================================================================================================!
!  Signal Wave								   				!
!=======================================================================================================!
subroutine signalwave(contrauv,contrauv_halo,gh11,gh22,gh11_halo,gh22_halo,fjmax,lambda,SL,SR)
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
! real (kind=real_kind):: alfa1,alfa2
!=======================================================================================================!
!    alfa1 = 0.0D0
!    alfa2 = 0.0D0
!    do i = 1, np
!       ul = abs(contrauv(i,1,2))      + sqrt(gh22(i,1,2))
!       ur = abs(contrauv_halo(i,south,2)) + sqrt(gh22_halo(i,south))
!       alfa1= max(alfa1,ul,ur)
!       
!       ul = abs(contrauv_halo(i,north,2)) + sqrt(gh22_halo(i,north))
!       ur = abs(contrauv(i,np,2))     + sqrt(gh22(i,np,2))
!       alfa2= max(alfa2,ul,ur)
!    enddo
!    fjmax(south) = alfa1
!    fjmax(north) = alfa2
!    alfa1 = 0.0D0
!    alfa2 = 0.0D0
!    do j = 1, np
!       ul = abs(contrauv_halo(j,west,1))  + sqrt(gh11_halo(j,west))
!       ur = abs(contrauv(1,j,1))      + sqrt(gh11(1,j,1))
!       alfa1= max(alfa1,ul,ur)
!       ul = abs(contrauv(np,j,1))     + sqrt(gh11(np,j,1))
!       ur = abs(contrauv_halo(j,east,1))  + sqrt(gh11_halo(j,east))
!       alfa2= max(alfa2,ul,ur)
!    enddo
!    fjmax(west) = alfa1
!    fjmax(east) = alfa2 
!=======================================================================================================!
    SL= 0.0D0
    SR= 0.0D0
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
    lambda= 0.0D0
    do wall=1,4
    do i=1,np
       lambda(wall,i)= max( SL(wall,i), SR(wall,i) )
    enddo   
    enddo
    fjmax= 0.0D0
    do wall=1,4
    do i=1,np 
       fjmax(wall)= max(lambda(wall,i),fjmax(wall))
    enddo
    enddo 
!=======================================================================================!    
    SR(:,:)= lambda(:,:)
    SL(:,:)=-lambda(:,:)
!=======================================================================================!
 return
!=======================================================================================!
end subroutine signalwave    
!=======================================================================================!
!=======================================================================================!
end module dg_flux_mod


