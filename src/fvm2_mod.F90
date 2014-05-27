!-----------------------------------------------------------------------------------!
! fvm2_MOD File for the fvm project in HOMME, FLUX VERSION                      !
! Author: Christoph Erath                                                           !
! Date: 25.January 2010                                                             !
! MAIN module to run fvm on HOMME                                                 !
! 14.June 2011: reorganisation done                                                 !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module fvm_mod
  use kinds, only : real_kind
  
contains

subroutine cslam_run(elem,cellbuffer,pointsbuffer,red,par,ithr,nets,nete)
  ! ---------------------------------------------------------------------------------
  use dimensions_mod, only: nc
  ! ---------------------------------------------------------------------------------
  use fvm_control_volume_mod, only: fvm_mesh, fvm_struct
  ! ---------------------------------------------------------------------------------
  use fvm_flux_mod, only: simplified_flux
  ! ---------------------------------------------------------------------------------
  use fvm_ghostcell_mod, only: ghostcellpack, ghostcellunpack, ghostpointspack, ghostpointsunpack
  ! ---------------------------------------------------------------------------------
  use fvm_line_integrals_mod, only: my_compute_weights
  ! ---------------------------------------------------------------------------------
  use edge_mod, only : EdgeBuffer_t
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: bndry_exchangeV 
  ! ---------------------------------------------------------------------------------
  use hybrid_mod, only : hybrid_t, hybrid_create
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t, derivative_stag_t, derivinit, deriv_print
  ! ---------------------------------------------------------------------------------
  use element_mod, only : element_t
  ! ---------------------------------------------------------------------------------
  use parallel_mod, only : parallel_t
  ! ---------------------------------------------------------------------------------
  use reduction_mod, only : ReductionBuffer_ordered_1d_t
  ! ---------------------------------------------------------------------------------
  use thread_mod, only : nthreads
  ! ---------------------------------------------------------------------------------
  use hybrid_mod, only : hybrid_t, hybrid_create
  ! ---------------------------------------------------------------------------------
  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, ne
  ! ---------------------------------------------------------------------------------
  use time_mod, only : tstep, nmax, timelevel_t, time_at, timelevel_update, timelevel_init
  ! ---------------------------------------------------------------------------------
  use thread_mod, only : NThreads
  ! ---------------------------------------------------------------------------------
  use coordinate_systems_mod, only : spherical_polar_t,spherical_to_cart, cart2cubedspherexy
  ! ---------------------------------------------------------------------------------
  
#ifdef PIO_INTERP
    use interp_movie_mod, only : interp_movie_init, interp_movie_output, interp_movie_finish
#else
    use shal_movie_mod, only : shal_movie_init, shal_movie_output, shal_movie_finish
#endif
  
  implicit none

  type (element_t), intent(inout),pointer     :: elem(:)
  type (EdgeBuffer_t), intent(in)             :: cellbuffer, pointsbuffer ! buffer entity             (shared)
  type (ReductionBuffer_ordered_1d_t),intent(in)    :: red   ! reduction buffer         (shared)
  type (parallel_t), intent(in)               :: par   ! distributed parallel structure (shared)
  integer, intent(in)                         :: ithr  ! thread number                  (private)
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)

  type (fvm_struct), pointer :: fvm(:)

  real (kind=real_kind)                       :: spherearea, massstart, mass, maxc, maxcstart  
  integer                                     :: i,j,k,ie, kptr,time
  type (hybrid_t)                             :: hybrid
  type (TimeLevel_t)                          :: tl              ! time level struct
  type (derivative_t)                         :: deriv           ! derivative struct
  type (spherical_polar_t)                    :: tmpsphereincart   
 
  character (len=99)                          :: filename
 !-----------------------------------------------------------------------------------!  
   if(par%masterproc) then 
     print *,"!-----------------------------------------------------------------------!"
     print *,"!  SOLID-BODY ROTATION Test CASE for simplified fvm, Christoph Erath  !" 
     print *,"!-----------------------------------------------------------------------!" 
   endif
     
  ! for communication of threads (exchange data)  
  hybrid = hybrid_create(par,ithr,NThreads)  
  ! initialize time management
  call timeLevel_init(tl)
  
  spherearea=0
  massstart=0
  
  !-Create new fvm Mesh
  allocate(fvm(nelemd))
  call fvm_mesh(elem,fvm,tl,nets,nete)
  !-Arrival and departure fvm Meshes, initialization done  
  if(par%masterproc) then 
    print *
    print *,"Arrival and departure grid created, initialization done. " 
    print *
  endif
!-----------------------------------------------------------------------------------!  
  
  ! Initialize derivative structure
  call derivinit(deriv)
    
!  massstart=massstart/spherearea
!-----------------------------------------------------------------------------------!  
kptr=0
!pre exchange of ghost cells
do ie=nets, nete
  call ghostcellpack(cellbuffer, elem(ie)%state%c(:,:,:,tl%n0),nlev,kptr,elem(ie)%desc) 
  ! reset the new unknown
    do j=1-nc,nc+nc
      do i=1-nc,nc+nc 
        elem(ie)%state%c(i,j,1,tl%np1)=0.0D0
      end do
    end do
    spherearea=spherearea+fvm(ie)%elem_areasphere
    massstart=massstart+fvm(ie)%elem_mass
end do    
call bndry_exchangeV(hybrid,cellbuffer)

do ie=nets, nete
  call ghostcellunpack(cellbuffer, elem(ie)%state%c(:,:,:,tl%n0),nlev,kptr,elem(ie)%desc)
end do
!-----------------------------------------------------------------------------------!  
!Initialize Output via geopotential (should be changed, separate output for fvm
!write first time step to IO 
#ifdef PIO_INTERP
	  call interp_movie_init(elem,hybrid,nets,nete,tl=tl)    
    call interp_movie_output(elem,tl, hybrid, 0D0, nets, nete)
#else
	  call shal_movie_init(elem,hybrid)
    call shal_movie_output(elem,tl, hybrid, 0D0, nets, nete,deriv)
#endif 
!-----------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------!  
! BEGIN Testoutput: write data in p (interpolation) to use existing IO
maxcstart=0
  do ie=nets, nete
    do j=1,nc+1
      do i=1,nc+1 
         !write it in p because of IO 
         elem(ie)%state%p(i,j,1,tl%n0)=(elem(ie)%state%c(i-1,j-1,1,tl%n0)+ &
          elem(ie)%state%c(i,j-1,1,tl%n0)+elem(ie)%state%c(i-1,j,1,tl%n0)+elem(ie)%state%c(i,j,1,tl%n0))/4
      end do
    end do
  end do
  
!-----------------------------------------------------------------------------------!  
#ifdef PIO_INTERP
       call interp_movie_output(elem, tl, hybrid, 0D0, nets, nete)
#else     
       call shal_movie_output(elem, tl, hybrid, 0D0, nets, nete,deriv)
#endif
!-----------------------------------------------------------------------------------!
!BEGIN TIME LOOP, start at 0, calculate then next step
do while(tl%nstep<0)  
  ! calculate with simplified flux fvm
  call simplified_flux(elem(ie)%state%c(:,:,:,:),fvm,tl, nets, nete)
  ! go one time step forward and do the data exchange, new values shift from tl%np1 to tl%n0
  call TimeLevel_update(tl,"forward") 
  ! calculated values are now in tl%n0  
  do ie=nets, nete
    call ghostcellpack(cellbuffer, elem(ie)%state%c(:,:,:,tl%n0),nlev,kptr,elem(ie)%desc) 
    ! reset the new unknown
      do j=1-nhc,nc+nhc
        do i=1-nhc,nc+nhc 
          elem(ie)%state%c(i,j,1,tl%np1)=0.0D0
        end do
      end do
  end do    
  call bndry_exchangeV(hybrid,cellbuffer)
  
  do ie=nets, nete
    call ghostcellunpack(cellbuffer, elem(ie)%state%c(:,:,:,tl%n0),nlev,kptr,elem(ie)%desc)
      ! calculate fluxes for all elements, write in ghostcells fluxes, which are not
      ! needed just because of conservation of fluxes
  end do
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! BEGIN Testoutput: Just for test reasons
  mass=0D0
  do ie=nets, nete
    !print *, ""
    !write (*,0817) elem(ie)%LocalId,elem(ie)%GlobalId,elem(ie)%FaceNum
    do j=1,nc+1
      do i=1,nc+1 
         if ((i<=nc) .AND. (j<=nc)) then  
       !    write (*,*) elem(ie)%state%c(i,j,1,tl%n0),fvm(ie)%area_sphere(i,j)      
       !    print *
           mass=mass+fvm(ie)%area_sphere(i,j)*elem(ie)%state%c(i,j,1,tl%n0)
         endif
         !write it in p because of IO 
         elem(ie)%state%p(i,j,1,tl%n0)=(elem(ie)%state%c(i-1,j-1,1,tl%n0)+ &
         elem(ie)%state%c(i,j-1,1,tl%n0)+elem(ie)%state%c(i-1,j,1,tl%n0)+elem(ie)%state%c(i,j,1,tl%n0))/4
      end do
    end do
  end do 
! END Testoutput  
0817 format("Element ",I4,2x,I4,2x,I1)
 
  if  (par%masterproc) then 
    time=Time_at(tl%nstep)
    write(*,*) 'time=', time, 'timeatmax',Time_at(nmax),'step',tl%nstep,'maxstep',nmax, 't0', tl%n0, 't1', tl%np1
    write(*,*) 'massbegin', massstart, 'massend', mass, (mass-massstart)/massstart
 endif  
! test mass, just with one processor right
 ! mass=mass/spherearea

!   if (abs((mass-massstart)/massstart)>1D-14) then
!     print *, "FATAL Error in mass conservation, value greater than 1D-14"
!     if(par%masterproc) then 
!       time=Time_at(tl%nstep)
!       write(*,*) 'time=', time, 'timeatmax',Time_at(nmax),'step',tl%nstep,'maxstep',nmax
!       print *
!     endif
!     STOP "Exit program!"
!   endif
!-----------------------------------------------------------------------------------!  
#ifdef PIO_INTERP
       call interp_movie_output(elem, tl, hybrid, 0D0, nets, nete)
#else     
       call shal_movie_output(elem, tl, hybrid, 0D0, nets, nete,deriv)
#endif
!-----------------------------------------------------------------------------------!  
end do
!END TIME LOOP
!-----------------------------------------------------------------------------------!  
#ifdef PIO_INTERP
    call interp_movie_finish
#else
    call shal_movie_finish
#endif
!-----------------------------------------------------------------------------------!  


!SUMMARY
if(par%masterproc) then 
  print *
  print *,"!-----------------------------------------------------------------------!"
  print *,"!  SOLID-BODY ROTATION Test CASE for simplified fvm, Christoph Erath  !" 
  print *,"!-----------------------------------------------------------------------!"
  print *  
endif

maxc=0
mass=0
maxcstart=0
do ie=nets, nete
!   print *, ""
!   write (*,0817) elem(ie)%LocalId,elem(ie)%GlobalId,elem(ie)%FaceNum
  do j=1,nc
    do i=1,nc 
!          write (*,*) elem(ie)%state%c(i,j,1,tl%n0),fvm(ie)%area_sphere(i,j)
         if (maxc<elem(ie)%state%c(i,j,1,tl%n0)) then
           maxc=elem(ie)%state%c(i,j,1,tl%n0)   
         endif 
         mass=mass+fvm(ie)%area_sphere(i,j)*elem(ie)%state%c(i,j,1,tl%n0)      
    end do
  end do
  if (maxcstart<fvm(ie)%maxc) then
    maxcstart=fvm(ie)%maxc   
  endif
end do
print *
write(*,*) 'maxvaluestart:', maxcstart, 'maxvalue:', maxc, 'massstart',&
            massstart, 'mass', mass, 'massrel', (mass-massstart)/massstart
write(*,*) 'number of elements', 6*ne*ne*nc*nc


!-END TESTENVIRONMENT__--------------------------------------------------------------
end subroutine cslam_run
end module fvm_mod
