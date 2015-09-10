#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#undef  _CUBE_VEL
#define _SPHERE_VEL

module dg_movie_mod
#ifndef PIO_INTERP
  ! ---------------------
  use kinds, only : real_kind
  ! ---------------------
  use dimensions_mod, only : ne, np, nelem, nelemd, nlev, nelemdmax, &
       GlobalUniqueCols, npsq
  ! ---------------------
  use hybrid_mod, only : hybrid_t
  ! ---------------------
#ifdef _MPI
  use parallel_mod, only : syncmp, haltmp, iam, abortmp, mpireal_t, mpi_max, mpi_sum, mpiinteger_t
#else
  use parallel_mod, only : syncmp, haltmp, iam, abortmp, mpireal_t
#endif
  ! ---------------------
  use time_mod, only : timelevel_t , tstep, ndays, time_at, secpday 
  ! ---------------------
  use cube_mod, only : cube_assemble
  ! ---------------------
  use control_mod, only : test_case, runtype
  ! ---------------------
  use element_mod
  ! ---------------------
  use coordinate_systems_mod, only : cartesian2d_t, spherical_polar_t
  ! ---------------------
  use common_io_mod, only : &
       output_start_time,   &
       output_end_time,     &
       output_frequency,    &
       max_output_variables,&
       max_output_streams,  &
       varname_len,         &
       nfsizekind, &
       nf_handle, &
       nf_selectedvar, &
       nf_double, &
       nf_int, &
       get_current_varnames

  use netcdf_io_mod, only : &
       nf_output_init_begin,&
       nf_output_register_dims, &
       nf_output_register_variables,&
       nf_variable_attributes, &
       nf_output_init_complete,  &
       nf_put_var, &
       nf_advance_frame, &
       nf_close_all, &
       nf_get_frame


  ! ---------------------
  use dof_mod, only : UniqueNcolsP, Uniquepoints, UniqueCoords, CreateUniqueIndex
  ! ---------------------
  ! ------------------------
#ifdef _PRIMDG
  use interpolate_mod   
#endif
  ! ------------------------
  implicit none

  !=======================================================================================================!
  private
  public :: dg_movie_init
  public :: dg_movie_output
  public :: dg_movie_finish  
  public :: setvarnames
  !=======================================================================================================!

#ifdef _PRIMDG  
  integer, parameter :: varcnt =11
  integer, parameter :: maxdims=4
  character*(*),parameter::dimnames(maxdims)=(/'ncolp','nlev ','nelem','time '/)  
  integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,4,0,0,  &
       1,2,4,0,  &
       1,2,4,0,  &
       1,2,4,0,  &
       1,0,0,0,  &
       1,0,0,0,  &
       4,0,0,0,	 &
       1,2,4,0,	 &
       1,2,4,0,	 &
       1,2,4,0,	 &
       1,2,4,0/),&
       shape=(/maxdims,varcnt/))
  character*(*),parameter::varnames(varcnt)=(/'ps   ','geop ','u    ','v    ',&
       'lonp ','latp ','time ',	      &
       'T    ','p    ','zeta ','dgs  '/)
  integer, parameter :: vartype(varcnt)=(/nf_double,nf_double,nf_double,nf_double,&
       nf_double,nf_double,nf_double,&
       nf_double,nf_double,nf_double,nf_double/)
  logical, parameter :: varrequired(varcnt)=(/.false.,.false.,.false.,.false.,&
       .true.,.true.,.true.,&
       .false.,.false.,.false.,.false./)
#else
  integer, parameter :: varcnt =7
  integer, parameter :: maxdims=4
  character*(*),parameter::dimnames(maxdims)=(/'ncolp','nlev ','nelem','time '/)  
  integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,4,0,0,  &
       1,2,4,0,  &
       1,2,4,0,  &
       1,2,4,0,  &
       1,0,0,0,  &
       1,0,0,0,  &
       4,0,0,0/),&
       shape=(/maxdims,varcnt/))
  character*(*),parameter::varnames(varcnt)=(/'ps   ','geop ','u    ','v    ',&
       'lonp ','latp ','time '/)
  integer, parameter :: vartype(varcnt)=(/nf_double,nf_double,nf_double,nf_double,&
       nf_double, nf_double, nf_double/)
  logical, parameter :: varrequired(varcnt)=(/.false.,.false.,.false.,.false.,&
       .true.,.true.,.true./)
#endif

  ! local size of variable block for output
  type(nf_handle), target :: ncdf(max_output_streams)

  !=======================================================================================================!
  !=======================================================================================================!
contains
  !=======================================================================================================!
  subroutine dg_movie_init(elem, hybrid)
    type (element_t), intent(in) :: elem(:)
    type (hybrid_t), intent(in)     :: hybrid
    ! Local variables
    integer ie,i,j,k,ios
    integer :: v1(4), vstart, gstart
    integer(kind=nfsizekind) :: start(1), count(1)
    real (kind=real_kind), dimension(npsq) :: latp,lonp
    integer :: istartp, istartv, nxypg, nxyvg, ierr, elemstart
    integer :: dimsize(maxdims)
    integer, allocatable :: startg(:,:), sendbuf(:)

    call nf_output_init_begin(ncdf,hybrid%par%masterproc,hybrid%par%nprocs,hybrid%par%rank, &
         hybrid%par%comm,test_case,runtype)

    nxyp=0
    nxyv=0
    do ie=1,nelemd
       nxyp=nxyp+elem(ie)%idxP%NumUniquePts
       nxyv=nxyv+elem(ie)%idxV%NumUniquePts
    enddo
    dimsize = (/nxyp,nlev,nelem,0/)

    call nf_output_register_dims(ncdf, maxdims, dimnames, dimsize)
    call nf_output_register_variables(ncdf,varcnt,varnames,vardims,vartype,varrequired)
    call nf_variable_attributes(ncdf, 'latp', 'column latitude','radians')
    call nf_variable_attributes(ncdf, 'lonp', 'column longitude','radians')
    call nf_variable_attributes(ncdf, 'time', 'Model elapsed time','days')
    call nf_output_init_complete(ncdf)

    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0) ) then
          do ie=1,nelemdmax
             if(ie<=nelemd) then
                start(1)=elem(ie)%idxP%UniquPtOffset
                count(1)=elem(ie)%idxP%NumUniquePts
                call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep,latp, lonp)

             else
                count=0
             end if
             !
             ! Because nf_put_var is collective (when using pnetcdf) we need to assure that it is called
             ! the same number of times on each proc
             !
             call nf_put_var(ncdf(ios),lonp,start, count, name='lonp')
             call nf_put_var(ncdf(ios),latp,start, count, name='latp')
          enddo
       end if
    end do

  end subroutine dg_movie_init

  subroutine dg_movie_output(elem, tl,hybrid, phimean)
    use physical_constants, only : g
    use time_mod, only : Timelevel_t, time_at
    type (element_t), intent(in) :: elem(:)
    type (TimeLevel_t), intent(in)  :: tl
    type (hybrid_t), intent(in)     :: hybrid
    real (kind=real_kind), intent(in) :: phimean
    real (kind=real_kind) :: varptmp(np,np,nlev)
    integer :: ios, ierr, istat(4)
    integer :: vcntv2d, vcntv3d, vcntp3d,vcntp2d, ie, k
    character(len=varname_len), pointer :: output_varnames(:)
    real(kind=real_kind),parameter :: dayspersec=1./(3600.*24.)
    integer(kind=nfsizekind) :: start(3), count(3), start2d(2),count2d(2)
    integer :: ncnt 

    real (kind=real_kind)  :: varp2d(npsq)
    real (kind=real_kind)  :: varp3d(npsq,nlev)

    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0) .and. (output_start_time(ios) .le. tl%nstep) .and. &
            (output_end_time(ios) .ge. tl%nstep) .and. MODULO(tl%nstep,output_frequency(ios)) .eq. 0) then

          output_varnames => get_current_varnames(ios)
          start2d(1)=1
          start2d(2)=nf_get_frame(ncdf(ios))
          count2d(1)=0
          count2d(2)=1
          start(1)=1
          start(2)=1
          count(1)=0
          count(2)=nlev
          start(3)=nf_get_frame(ncdf(ios))
          count(3)=1

          if(nf_selectedvar('ps', output_varnames)) then
             do ie=1,nelemdmax
                if(ie<=nelemd) then
                   start2d(1) = elem(ie)%idxP%UniquPtOffset
                   count2d(1) = elem(ie)%idxP%NumUniquePts
                   ncnt = count2d(1)
#ifdef _PRIMDG
                   call UniquePoints(elem(ie)%idxP,elem(ie)%state%pr3d(:,:,nlev+1),varp2D)
#else
                   call UniquePoints(elem(ie)%idxP,elem(ie)%state%ps,varp2D)
#endif
                else
                   ncnt=1
                   count2d=0
                end if
                call nf_put_var(ncdf(ios),varp2d(1:ncnt),start2d,count2d,name='ps')
             enddo
          endif
#ifdef _PRIMDG          
	  if(nf_selectedvar('dgs', output_varnames)) then
             do ie=1,nelemdmax
                if(ie<=nelemd) then
                   start(1) = elem(ie)%idxP%UniquPtOffset
                   count(1) = elem(ie)%idxP%NumUniquePts
                   count(2) = 3
                   count(3) = 1
                   ncnt = count(1)                   
		   do k=1,3
                      varptmp(:,:,k) = elem(ie)%state%out3d(:,:,k)
                   end do
                   call UniquePoints(elem(ie)%idxP,3,varptmp,varp3D)		
		else
                   ncnt=1
                   count=0                
		end if
                call nf_put_var(ncdf(ios),varp3d(1:ncnt,:),start,count,name='dgs')             
	     enddo
	  endif
          if(nf_selectedvar('p', output_varnames)) then             
             do ie=1,nelemdmax
                if(ie<=nelemd) then
                   start(1) = elem(ie)%idxP%UniquPtOffset
                   count(1) = elem(ie)%idxP%NumUniquePts
                   count(2) = nlev
                   count(3) = 1
                   ncnt = count(1)		   
		   do k=1,nlev
                      varptmp(:,:,k) = elem(ie)%state%pr3d(:,:,k)
                   end do
                   call UniquePoints(elem(ie)%idxP,nlev,varptmp,varp3D)
		else
                   ncnt=1
                   count=0
                end if
                call nf_put_var(ncdf(ios),varp3d(1:ncnt,:),start,count,name='p')             
	     enddo
	  endif
          if(nf_selectedvar('T', output_varnames)) then
             do ie=1,nelemdmax
                if(ie<=nelemd) then
                   start(1) = elem(ie)%idxP%UniquPtOffset
                   count(1) = elem(ie)%idxP%NumUniquePts
                   count(2) = nlev
                   count(3) = 1
                   ncnt = count(1)		   
		   do k=1,nlev
                      varptmp(:,:,k) = elem(ie)%state%t3d(:,:,k)
                   end do
                   call UniquePoints(elem(ie)%idxP,nlev,varptmp,varp3D)
		else
                   ncnt=1
                   count=0                
		end if
                call nf_put_var(ncdf(ios),varp3d(1:ncnt,:),start,count,name='T')             
	     enddo
	  endif
	  if(nf_selectedvar('zeta', output_varnames)) then             
             do ie=1,nelemdmax
                if(ie<=nelemd) then
                   start(1) = elem(ie)%idxP%UniquPtOffset
                   count(1) = elem(ie)%idxP%NumUniquePts
                   count(2) = nlev
                   count(3) = 1
                   ncnt = count(1)		   
		   do k=1,nlev
                      varptmp(:,:,k) = elem(ie)%state%zeta(:,:,k)
                   end do
                   call UniquePoints(elem(ie)%idxP,nlev,varptmp,varp3D)
		else
                   ncnt=1
                   count=0
                end if
                call nf_put_var(ncdf(ios),varp3d(1:ncnt,:),start,count,name='zeta')             
	     enddo
	  endif
#endif

          if(nf_selectedvar('geop', output_varnames)) then
             start(1)=1
             count(1)=0
             do ie=1,nelemdmax
                if(ie<=nelemd) then
                   start(1) = elem(ie)%idxP%UniquPtOffset
                   count(1) = elem(ie)%idxP%NumUniquePts
                   count(2) = nlev
                   count(3) = 1
                   ncnt = count(1)
#ifdef _SWDG
		   call UniquePoints(elem(ie)%idxP,nlev,elem(ie)%state%ht,varp3D)
#else
                   do k=1,nlev
                      varptmp(:,:,k) = (elem(ie)%state%p(:,:,k,tl%n0) + elem(ie)%state%ps + phimean)/g
                   end do
                   call UniquePoints(elem(ie)%idxP,nlev,varptmp,varp3D)
#endif
                else
                   ncnt=1
                   count=0
                end if

                call nf_put_var(ncdf(ios),varp3d(1:ncnt,:),start,count,name='geop')
             enddo
          endif

          if(nf_selectedvar('u', output_varnames)) then
             do ie=1,nelemdmax
                if(ie<=nelemd) then
                   do k=1,nlev
                      varvTMP(:,:,k) = elem(ie)%D(1,1,:,:)*elem(ie)%state%v(:,:,1,k,tl%n0)+ &
                           elem(ie)%D(1,2,:,:)*elem(ie)%state%v(:,:,2,k,tl%n0)
                   end do
                   start(1)=elem(ie)%idxV%UniquPtOffset
                   count(1)=elem(ie)%idxV%NumUniquePts
                   count(2)=nlev
                   count(3)=1
                   ncnt=1
                   call UniquePoints(elem(ie)%idxV,nlev,varvtmp,varv3D)
                else
                   ncnt=1
                   count=0
                endif
                call nf_put_var(ncdf(ios), varv3D(1:ncnt,:),start,count,name='u')
             enddo
          endif

          if(nf_selectedvar('v', output_varnames)) then
             do ie=1,nelemdmax
                if(ie<=nelemd) then
                   do k=1,nlev
                      varvTMP(:,:,k) = elem(ie)%D(:,:,2,1)*elem(ie)%state%v(:,:,1,k,tl%n0)+ &
                           elem(ie)%D(:,:,2,2)*elem(ie)%state%v(:,:,2,k,tl%n0)
                   end do
                   start(1)=elem(ie)%idxV%UniquPtOffset
                   count(1)=elem(ie)%idxV%NumUniquePts
                   count(2)=nlev
                   count(3)=1
                   ncnt=1
                   call UniquePoints(elem(ie)%idxV,nlev,varvtmp,varv3D)
                else
                   ncnt=1
                   count=0
                endif
                call nf_put_var(ncdf(ios), varv3D(1:ncnt,:),start,count,name='v')
             enddo
          endif

          count(3) = 1
          call nf_put_var(ncdf(ios),real(dayspersec*time_at(tl%nstep),kind=real_kind),&
               start(3:3),count(3:3),name='time')
          call nf_advance_frame(ncdf(ios))

       end if


    end do
  end subroutine dg_movie_output

  subroutine dg_movie_finish
    integer :: istat(4)

    call nf_close_all(ncdf)
    istat=0
  end subroutine dg_movie_finish

  ! 
  ! Called by control_mod to set the list of variables available in this model
  !
  subroutine setvarnames(nlvarnames)
    character*(*), intent(out) :: nlvarnames(:)
    nlvarnames(1:varcnt) = varnames
  end subroutine setvarnames
#endif
  !=======================================================================================================!
end module dg_movie_mod


