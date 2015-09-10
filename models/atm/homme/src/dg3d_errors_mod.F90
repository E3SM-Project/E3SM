#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dg3d_errors_mod
!=======================================================================================================!
!=======================================================================================================!
  ! ---------------------
  use cube_mod
  ! ------------------------
  use kinds
  ! ------------------------
  use physical_constants
  ! ------------------------
  use dimensions_mod
  ! ------------------------
  use derivative_mod
  ! ------------------------
  use edge_mod
  ! ------------------------
  use bndry_mod
  ! ------------------------
  use coordinate_systems_mod
  ! ------------------------
  use quadrature_mod
  ! ------------------------
  use ref_state_mod
  ! ------------------------
  use global_norms_mod
  ! ------------------------
  use time_mod
  ! ------------------------
  use element_mod, only : element_t
  ! ------------------------
  !use state_mod
  ! ------------------------
  use hybrid_mod
  ! ------------------------
  use reduction_mod, only : red_max, red_sum, pmax_mt
  ! ------------------------  
  use dg3d_vertical_mod, only : eta_levels, linear_int, index_search
  ! ------------------------
  use interpolate_mod  
  ! ------------------------
  use parallel_mod, only : syncmp, global_shared_buf, global_shared_sum
!=======================================================================================================!
!=======================================================================================================!
implicit none
private
!=======================================================================================================! 
  real (kind=real_kind), private :: delta  = 1.0D-12 
!=======================================================================================================! 
  public  :: jw_bcl_errors
  public  :: slice850
  public  :: slice_at_prlevel
!=======================================================================================================! 
  private :: l2_norm  
  private :: sphere_integral 
!=======================================================================================================!    
  public  :: cube2polar  
  public  :: jw_bcl_zonal
!=======================================================================================================!
 contains
!=======================================================================================================!
!	jw_bcl_errors:
!=======================================================================================================!
subroutine jw_bcl_errors(elem,iounit,tl,hybrid,nets,nete)
    type(element_t) , intent(in) :: elem(:)
    integer   :: iounit
    type (TimeLevel_t)   , intent(in) :: tl         ! model time struct
    type (hybrid_t)      , intent(in) :: hybrid     
    integer, intent(in)   :: nets
    integer, intent(in)   :: nete

    ! Local variables 
    real (kind=real_kind) :: time_tmp,sum_1,sum_2,sum_del
    real (kind=real_kind) :: eta(nlev), etv(nlev), etk(nlev+1)
    real (kind=real_kind) :: ak(nlev+1), bk(nlev+1)
    real (kind=real_kind) :: l2_1(nlev),l2_2(nlev),del_etak(nlev)
    real (kind=real_kind) :: u_true(np,np,nlev,nets:nete)
    real (kind=real_kind) :: u_time(np,np,nlev,nets:nete)    

    integer ie,k
!=======================================================================================================!
    if (tl%nstep == 0) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          open(iounit+0,file='./errors/l2norm_1.dat',form="formatted")          
	  open(iounit+1,file='./errors/l2norm_2.dat',form="formatted")
       end if
    end if
!=======================================================================================================!
    do ie=nets,nete   
    do k=1,nlev 
       u_true(:,:,k,ie)= elem(ie)%state%uv0(:,:,1,k)
       u_time(:,:,k,ie)= elem(ie)%state%uv(:,:,1,k)
    end do
    end do
!=======================================================================================================!
    do k=1,nlev     
       l2_1(k)= l2_norm(elem,1,u_time(:,:,k,nets:nete),u_true(:,:,k,nets:nete),hybrid,np,nets,nete)
       l2_2(k)= l2_norm(elem,2,u_time(:,:,k,nets:nete),u_true(:,:,k,nets:nete),hybrid,np,nets,nete)
    end do
    
    call eta_levels(ak,bk,eta,etk)
    do k=1,nlev     
       del_etak(k)= abs(etk(k+1)-etk(k))
    enddo
    sum_1= 0.0D0    
    sum_2= 0.0D0
    sum_del= 0.0D0
    do k=1,nlev     
       sum_1= sum_1 + l2_1(k)*del_etak(k)       
       sum_2= sum_2 + l2_2(k)*del_etak(k)
       sum_del= sum_del + del_etak(k)
    enddo
    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       write(iounit+0,30)Time_at(tl%nstep)/secpday,SQRT(sum_1/sum_del)       
       write(iounit+1,30)Time_at(tl%nstep)/secpday,SQRT(sum_2/sum_del)
    end if
30  format(f11.6,4x,e13.6)
!=======================================================================================================!
end subroutine jw_bcl_errors
!=======================================================================================================!
function l2_norm(elem, l2_type,h,ht,hybrid,npts,nets,nete) result(l2)
   type(element_t), intent(in) :: elem(:)
   integer  , intent(in) :: l2_type,npts,nets,nete
   real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)  ! computed soln
   real (kind=real_kind), intent(in) :: ht(npts,npts,nets:nete) ! true soln
   type (hybrid_t)      , intent(in) :: hybrid
   real (kind=real_kind) :: l2   

   ! Local variables

   real (kind=real_kind) :: dh2(npts,npts,nets:nete)
   real (kind=real_kind) :: ht2(npts,npts,nets:nete)
   real (kind=real_kind) :: dh2_int
   real (kind=real_kind) :: ht2_int
   integer i,j,ie

   do ie=nets,nete
        do j=1,npts
           do i=1,npts
  	      dh2(i,j,ie)=(h(i,j,ie)-ht(i,j,ie))**2
  	      ht2(i,j,ie)= ht(i,j,ie)**2
           end do
        end do
   end do

   dh2_int = sphere_integral(elem, dh2(:,:,nets:nete),hybrid,npts,nets,nete)
   ht2_int = sphere_integral(elem, ht2(:,:,nets:nete),hybrid,npts,nets,nete)


!   l2 = SQRT(dh2_int)/SQRT(ht2_int)
if (l2_type==1) then
   l2 = dh2_int
else if (l2_type==2) then
   l2 = dh2_int/ht2_int
end if

end function l2_norm  
!=======================================================================================================!
function sphere_integral(elem,h,hybrid,npts,nets,nete) result(I_sphere)
    type(element_t), intent(in) :: elem(:)
    integer  , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind) :: I_sphere

    real (kind=real_kind) :: I_priv
    real (kind=real_kind) :: I_shared
    common /gblintcom/I_shared

    ! Local variables

    integer :: ie,j,i
    real(kind=real_kind) :: I_tmp(1)

    real (kind=real_kind) :: da
    
    I_shared = 0.0D0
    if (npts==np) then
       do ie=nets,nete
          I_priv   = 0.0D0
          do j=1,np
          do i=1,np
           da = elem(ie)%mp(i,j)*elem(ie)%metdet(i,j)
           I_priv = I_priv + da*h(i,j,ie)
          end do
          end do
          global_shared_buf(ie,1) = I_priv
       end do
    else if (npts==np) then
       do ie=nets,nete
          I_priv   = 0.0D0
          do j=1,np
          do i=1,np
           da = elem(ie)%mp(i,j)*elem(ie)%metdet(i,j)
           I_priv = I_priv + da*h(i,j,ie)
          end do
          end do
          global_shared_buf(ie,1) = I_priv
       end do
    end if
    call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
    red_sum%buf(1) = global_shared_sum(1)

    I_sphere = red_sum%buf(1)/(4.0D0*DD_PI)

end function sphere_integral
!=======================================================================================================!
function slice850(ie,var3d) result(var850)
 implicit none
 integer, intent(in)  :: ie
 real (kind=real_kind), intent(in):: var3d(np,np,nlev)   
 real (kind=real_kind):: var850(np,np)
 real (kind=real_kind):: t14(np,np),t15(np,np)
 real (kind=real_kind):: t22(np,np),t23(np,np)
 real (kind=real_kind):: t43(np,np),t44(np,np)
 real (kind=real_kind):: p850,p14,p15,p22,p23,p43,p44
!=======================================================================================================!
if (nlev==18) then 
 p14 = 796.93134D0
 p15 = 877.88709D0
 p850= 850.D0
 t14 = 0.0D0
 t15 = 0.0D0
 var850= 0.0D0
 t14(:,:) =  var3d(:,:,14)
 t15(:,:) =  var3d(:,:,15)
 var850(:,:)= t14(:,:) + (p850-p14)*(t15(:,:)-t14(:,:))/(p15-p14)  
!=======================================================================================================!
elseif (nlev==26) then 
 p22 = 798.13911D0
 p23 = 878.65064D0
 p850= 850.D0
 t22 = 0.0D0
 t23 = 0.0D0
 var850= 0.0D0
 t22(:,:) =  var3d(:,:,22)
 t23(:,:) =  var3d(:,:,23)
 var850(:,:)= t22(:,:) + (p850-p22)*(t23(:,:)-t22(:,:))/(p23-p22)
!=======================================================================================================!
elseif (nlev==49) then 
 p43 = 817.2279823575D0
 p44 = 867.7626102075D0
 p850= 850.D0
 t43 = 0.0D0
 t44 = 0.0D0
 var850= 0.0D0
 !t43(:,:) =  var3d(:,:,43)
 !t44(:,:) =  var3d(:,:,44)
 var850(:,:)= t43(:,:) + (p850-p43)*(t44(:,:)-t43(:,:))/(p44-p43)
endif  
!=======================================================================================================!
end function slice850
!=======================================================================================================!
! Vertical interpolation at a given pressure level
!=======================================================================================================!
function slice_at_prlevel(plevel,pr3d,var3d) result(varout)
 implicit none
 real (kind=real_kind), intent(in):: plevel         
 real (kind=real_kind), intent(in):: pr3d(np,np,nlev+1)   
 real (kind=real_kind), intent(in):: var3d(np,np,nlev)   
 real (kind=real_kind):: varout(np,np)
 real (kind=real_kind):: pgrid(0:nlev+1) 
 real (kind=real_kind):: var(nlev) 
 real (kind=real_kind):: p1,p2,f1,f2,px, hpa
 integer :: i,j,k,ik 
!=======================================================================================================!
!Pressure converted  to hPa , "plevel" is assumed to be in hPa units 

       hpa = 0.01D0 

     do j=1,np
          do i=1,np
                pgrid(0) = pr3d(i,j,1) *hpa 
                pgrid(nlev+1) = pr3d(i,j,nlev+1) *hpa 

            do k=1,nlev
                pgrid(k)= (pr3d(i,j,k+1) + pr3d(i,j,k))*0.5D0 *hpa 
                var(k) =  var3d(i,j,k)
            end do

        ! The  index  search

             ik = index_search(pgrid,0,1,plevel)

             do k = 1, nlev
             !   p1 = pgrid(ik)
             !   p2 = pgrid(ik+1)
             !   px = plevel 
             !   f1 = var(ik)
             !   f2 = var(ik+1)
             !   varout(i,j) = linear_int(p1,p2,f1,f2,px)

          !vertical log/linear option

                 p1 = log(pgrid(ik))
                 p2 = log(pgrid(ik+1))
                 px = log(plevel) 
                 f1 = var(ik)
                 f2 = var(ik+1)
                varout(i,j) = linear_int(p1,p2,f1,f2,px)
               
           end do
       end do
    end do

!=======================================================================================================!
end function slice_at_prlevel
!=======================================================================================================!
!=======================================================================================================!
!	Polar Grid Output
!=======================================================================================================!
subroutine cube2polar(elem,tl,hybrid)
    implicit none    
    type (element_t), intent(in) :: elem(:)
    type (TimeLevel_t), intent(in) :: tl 
    type (hybrid_t)   , intent(in) :: hybrid   
    integer, parameter :: iunit=88    
    integer:: ie,ierr,i,j,k,klon,klat,klev
    type (quadrature_t)  :: gll
    type (interpolate_t) :: interp
    real (kind=real_kind):: lat,lon
    real (kind=real_kind):: local_cube(np,np,nlev)
    real (kind=real_kind):: global_cube(np*ne,np*ne,6,nlev)
    real (kind=real_kind),dimension(:,:,:),allocatable:: global_pg 
!=======================================================================================================!
    gll=gausslobatto(np)
    !call interpolate_create(gll,interp)
    
    if (tl%nstep == 0 .and. hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       print *,'GLL coordinate'    
       do i=1,np
          print *,i,interp%glp(i)
       end do
       print *,'GLL intpolant'
       do j=1,np
          print *,j,(interp%Imat(i,j),i=1,np)
       enddo	  
    endif       
!=======================================================================================================!    
    call syncmp(hybrid%par)
    do ie=1,nelemd
    do k= 1,nlev
       local_cube(:,:,k)= elem(ie)%state%uv(:,:,1,k)
    enddo
    ierr= cube_assemble(global_cube,local_cube,elem(ie),hybrid%par,nelemd,nelem,ie)
    enddo
!=======================================================================================================!     
    call syncmp(hybrid%par) 
    klon= (np-1)*ne*4 + 1
    klat= (np-1)*ne*2 + 1
    klev= 1 
    allocate(global_pg(klon,klat,nlev))

    call syncmp(hybrid%par)       
    if (tl%nstep == 0 .and. hybrid%par%masterproc .and. hybrid%ithr==0) then
      print *,'dimension (lon,lat):',klon,',',klat
      open (iunit,file='./data/uv.dat',form='formatted') 
      do j=1,klat
      lat= 0.5D0*DD_PI - (j-0.5D0)*(DD_PI/klat) + delta
        do i=1,klon
          lon= (i-1)*((2.0D0*DD_PI)/klon) + delta
          write(iunit,10)lon,lat,global_pg(i,j,klev)
        end do
      end do    
      close(iunit)
    endif    

! IO format statements..
 10 format(e21.15,1x,e21.15,1x,e21.15)
 20 format("")       
 deallocate(global_pg)
!=======================================================================================================!
end subroutine cube2polar
!=======================================================================================================!
!	jw_bcl_zonal:
!=======================================================================================================!
subroutine jw_bcl_zonal(elem,tl,hybrid)
    implicit none    
    type (element_t), intent(in) :: elem(:)
    type (TimeLevel_t), intent(in) :: tl 
    type (hybrid_t)   , intent(in) :: hybrid   
    integer, parameter :: iunit=88    
    integer:: ie,ierr,i,j,k,klon,klat,klev
    type (quadrature_t)  :: gll
    type (interpolate_t) :: interp
    real (kind=real_kind):: lat,lon,sum_global,sum_1,sum_2,sum_eta
    real (kind=real_kind):: local_cube(np,np,nlev)
    real (kind=real_kind):: global_cube(np*ne,np*ne,6,nlev)
    real (kind=real_kind),dimension(:,:,:),allocatable:: global_true,global_time
    real (kind=real_kind),dimension(:,:),allocatable  :: zonal_true,zonal_time    
    real (kind=real_kind) :: eta(nlev), etv(nlev), etk(nlev+1)
    real (kind=real_kind) :: ak(nlev+1), bk(nlev+1),del_etak(nlev)
!=======================================================================================================!
    gll=gausslobatto(np)
    !call interpolate_create(gll, interp)
    
    klon= (np-1)*ne*4 + 1
    klat= (np-1)*ne*2 + 1    
    if (tl%nstep == 0 .and. hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       print *,'GLL coordinate'    
       do i=1,np
          print *,i,interp%glp(i)
       end do
       print *,'GLL intpolant'
       do j=1,np
          print *,j,(interp%Imat(i,j),i=1,np)
       enddo    
    endif	  
    allocate(global_true(klon,klat,nlev))    
    allocate(global_time(klon,klat,nlev)) 
    allocate(zonal_true(klat,nlev))         
    allocate(zonal_time(klat,nlev))
!=======================================================================================================!        
    call syncmp(hybrid%par)   
    do ie=1,nelemd
    do k= 1,nlev
       local_cube(:,:,k)= elem(ie)%state%uv0(:,:,1,k)
    enddo
    ierr= cube_assemble(global_cube,local_cube,elem(ie),hybrid%par,nelemd,nelem,ie)
    enddo
    do k=1,nlev    
    do j=1,klat
       sum_global= 0.0D0
       do i=1,klon
          sum_global= sum_global + global_true(i,j,k)
       enddo
       zonal_true(j,k)= sum_global/klon
    enddo
    enddo 

    call syncmp(hybrid%par) 
    do ie=1,nelemd
    do k= 1,nlev
       local_cube(:,:,k)= elem(ie)%state%uv(:,:,1,k)
    enddo
    ierr= cube_assemble(global_cube,local_cube,elem(ie),hybrid%par,nelemd,nelem,ie)
    enddo
    do k=1,nlev    
    do j=1,klat
       sum_global= 0.0D0
       do i=1,klon
          sum_global= sum_global + global_time(i,j,k)
       enddo
       zonal_time(j,k)= sum_global/klon
    enddo
    enddo     
!=======================================================================================================!
    if (tl%nstep == 0 .and. hybrid%par%masterproc .and. hybrid%ithr==0) then
      open(iunit+0,file='./errors/l2uniform_1.dat',form="formatted")          
      open(iunit+1,file='./errors/l2uniform_2.dat',form="formatted")
    end if
   
    call eta_levels(ak,bk,eta,etk)
    sum_eta=0.0D0
    do k=1,nlev     
       del_etak(k)= abs(etk(k+1)-etk(k))
       sum_eta= sum_eta + del_etak(k)
    enddo
    
    sum_1=0.0D0
    sum_2=0.0D0  
    do k=1,nlev
    do j=1,klat    
       do i=1,klon    
          sum_1= sum_1 + ABS(global_time(i,j,k) - zonal_time(j,k))**2.0D0*del_etak(k)/(4.0D0*DD_PI)
       enddo
       sum_2= sum_2 + ABS(zonal_time(j,k) - zonal_true(j,k))**2.0D0*del_etak(k)/2.0D0
    enddo
    enddo    
    if (hybrid%par%masterproc .and. hybrid%ithr==0) then
       write(iunit+0,30)Time_at(tl%nstep)/secpday,SQRT(sum_1/sum_eta)       
       write(iunit+1,30)Time_at(tl%nstep)/secpday,SQRT(sum_2/sum_eta)
    end if
30  format(f11.6,4x,e13.6)    
!=======================================================================================================!
    deallocate(global_true)    
    deallocate(global_time) 
    deallocate(zonal_true)         
    deallocate(zonal_time)
!=======================================================================================================!
end subroutine jw_bcl_zonal
!=======================================================================================================! 
end module dg3d_errors_mod
