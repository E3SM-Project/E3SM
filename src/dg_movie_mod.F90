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
  use control_mod, only : test_case, runtype, accumstart, accumstop, accumfreq
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
#ifdef _PRIMDG
  public :: dg_accum_output  
  public :: dg_accum_latlon
#endif
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
                      varvTMP(:,:,k) = elem(ie)%D(2,1,:,:)*elem(ie)%state%v(:,:,1,k,tl%n0)+ &
                           elem(ie)%D(2,2,:,:)*elem(ie)%state%v(:,:,2,k,tl%n0)
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

  !=======================================================================================================!
  ! ===============================================
  ! dg_accum_output:
  !
  ! Accumulator routine...
  ! ===============================================
#ifdef _PRIMDG
  subroutine dg_accum_output(elem, tl, hybrid)
    type (element_t), target :: elem(:)
    type (TimeLevel_t)  :: tl
    type (hybrid_t)     :: hybrid

    ! Local variables

    type (cartesian2D_t), dimension(:,:), pointer :: cart

    integer ie,i,j,k,l,kadv,nx1,ny1
    integer ierr
    integer npts
    integer:: naccum, simday		     ! current day of simulation
    integer:: accumday= accumfreq*secpday 

    real (kind=real_kind), dimension(:,:,:), allocatable :: x,y
    real (kind=real_kind), dimension(:,:,:), allocatable :: T
    real (kind=real_kind), dimension(:,:,:), allocatable :: u

    real (kind=real_kind) :: v1,v2,epsi=1.0D-06
    real (kind=real_kind), dimension(:,:,:,:), allocatable :: gbl_fld
    real (kind=real_kind), dimension(:,:,:,:), allocatable :: gbl_x
    real (kind=real_kind), dimension(:,:,:,:), allocatable :: gbl_y

    logical               :: accum_on, accum_done
    real (kind=real_kind) :: rnaccump1
    character(len=2)      :: quadtype
    character(len=6)      :: chartag
    character(len=6)      :: qtag
    character(len=80)     :: filename
    !=======================================================================================================!
    quadtype = "gl"

    !$OMP BARRIER
    if (hybrid%ithr==0) then

       ! Take care of accumulator logic...

       if ( ABS(Time_at(tl%nstep)/secpday-accumstart)<epsi ) then
          if (hybrid%par%masterproc) then
             write(*,'(a)')'Accumulation for held-suarez started ============================'
          endif
          naccum=0
          do ie=1,nelemd
             elem(ie)%accum%T=0.0D0
             elem(ie)%accum%u=0.0D0
          end do
       end if

       if ( ABS(Time_at(tl%nstep)/secpday-accumstop)<epsi ) then          
          if (hybrid%par%masterproc) then
             write(*,'(a)')'Accumulation for held-suarez finished ============================'
          endif
          accum_done=.true.
       else
          accum_done=.false.
       end if

       if( Time_at(tl%nstep) > accumstart*secpday  .and. &
            Time_at(tl%nstep) < accumstop*secpday   .and. &
            MODULO(tl%nstep,accumfreq)==0 ) then          
          if (hybrid%par%masterproc) then
             write(*,'(a,1x,i8,a,f10.2,a)')'Accumulation for held-suarez at:',tl%nstep,', time=',Time_at(tl%nstep)/secpday,' days'
          endif
          accum_on=.true.
          rnaccump1=1.0D0/REAL(naccum+1,kind=real_kind)
       else
          accum_on=.false.
       end if

       if (ndays > 0) then
          simday=NINT(Time_at(tl%nstep)/secpday)
       else
          simday=tl%nstep
       end if


       !=======================================================================================================! 
       ! ===================================
       ! output T average...
       ! ===================================

       if (accum_on) then
          do ie=1,nelemd      
             elem(ie)%accum%T= rnaccump1*( elem(ie)%state%t3d + naccum*elem(ie)%accum%T)
          end do
          call syncmp(hybrid%par)
       end if
       !=======================================================================================================! 
       if (accum_done) then

          ! =======================================
          ! construct lon,lat for Temperature grid
          ! =======================================

          npts = SIZE(elem(1)%state%t3d,1)

          allocate(x(npts,npts,1))
          allocate(y(npts,npts,1))
          allocate(gbl_x(npts*ne,npts*ne,6,1))
          allocate(gbl_y(npts*ne,npts*ne,6,1))
          allocate(gbl_fld(npts*ne,npts*ne,6,nlev))

          do ie=1,nelemd
             if (npts==np) then
                cart => elem(ie)%cartp
             end if
             x(:,:,1)=cart(:,:)%x
             ierr=cube_assemble(gbl_x,x,elem(ie),hybrid%par,nelemd,nelem,ie)    
             if (ierr /= 0) call haltmp("error in prim_movie_output, halting...")   
          end do
          call syncmp(hybrid%par)

          do ie=1,nelemd
             if (npts==np) then
                cart => elem(ie)%cartp
             end if
             y(:,:,1)=cart(:,:)%y
             ierr=cube_assemble(gbl_y,y,elem(ie),hybrid%par,nelemd,nelem,ie)       
          end do
          call syncmp(hybrid%par)

          ! =======================================
          ! output accumulated average temperature
          ! =======================================

          !          write(chartag,'(i6)')simday
          !          filename=TRIM(ADJUSTL(output_dir))//TRIM(output_prefix)//TRIM(ADJUSTL(test_case))//".Tavg."//TRIM(ADJUSTL(chartag))

          do ie=1,nelemd      
             ierr=cube_assemble(gbl_fld,elem(ie)%accum%T,elem(ie),hybrid%par,nelemd,nelem,ie)
             if (ierr /= 0) call haltmp("error in prim_movie_output, halting...")
          end do
          call syncmp(hybrid%par)

          !          if (hybrid%par%masterproc) then
          !             call movie_output(ne,npts,quadtype,gbl_fld,gbl_x(:,:,:,1),gbl_y(:,:,:,1),filename)
          !          end if
          !          call syncmp(hybrid%par)
          !=======================================================================================================! 
          ! 3D writing: 												!
          !=======================================================================================================!             
          nx1 = ne*npts
          ny1 = ne*npts      
          if (hybrid%par%masterproc) then
             open (unit=37,file='./data/accum_temp.out')
             write(37,9)((((gbl_fld(i,j,k,l),i=1,nx1),j=1,ny1),k=1,6),l=1,nlev)
             close(37)
9            format(1x,5E16.6)
          end if
          call syncmp(hybrid%par)
          !=======================================================================================================!
          deallocate(x)
          deallocate(y)
          deallocate(gbl_fld)
          deallocate(gbl_x)
          deallocate(gbl_y)

       end if
       !=======================================================================================================! 
       if (accum_on) then

          ! ========================
          ! output u velocity
          ! ========================

          npts=SIZE(elem(1)%state%v(:,:,1,1,tl%n0),1)

          allocate(u(npts,npts,nlev))

          do ie=1,nelemd
             do k=1,nlev
                do j=1,npts
                   do i=1,npts
                      v1 = elem(ie)%state%v(i,j,1,k,tl%n0)
                      v2 = elem(ie)%state%v(i,j,2,k,tl%n0)

#if defined(_CUBE_VEL)
                      u(i,j,k)= v1
#else
                      u(i,j,k)= v1*elem(ie)%D(1,1,i,j) + v2*elem(ie)%D(1,2,i,j)
#endif
                   end do
                end do
             end do
             elem(ie)%accum%u= rnaccump1*(u + naccum*elem(ie)%accum%u)
          end do

          deallocate(u)

          call syncmp(hybrid%par)

       end if
       !=======================================================================================================! 
       if (accum_done) then

          ! ===================================
          ! construct lon,lat on velocity grid
          ! ===================================

          npts=SIZE(elem(1)%accum%u,1)

          allocate(x(npts,npts,1))
          allocate(y(npts,npts,1))
          allocate(gbl_x(npts*ne,npts*ne,6,1))
          allocate(gbl_y(npts*ne,npts*ne,6,1))
          allocate(gbl_fld(npts*ne,npts*ne,6,nlev))

          do ie=1,nelemd
             if (npts==np) then
                cart => elem(ie)%cartp
             end if
             x(:,:,1)=cart(:,:)%x
             ierr=cube_assemble(gbl_x,x,elem(ie),hybrid%par,nelemd,nelem,ie)       
             if (ierr /= 0) call haltmp("error in prim_movie_output, halting...")
          end do
          call syncmp(hybrid%par)

          do ie=1,nelemd
             if (npts==np) then
                cart => elem(ie)%cartp
             end if

             y(:,:,1)=cart(:,:)%y
             ierr=cube_assemble(gbl_y,y,elem(ie),hybrid%par,nelemd,nelem,ie)       

             if (ierr /= 0) call haltmp("error in prim_movie_output, halting...")
          end do
          call syncmp(hybrid%par)

          !          write(chartag,'(i6)') simday
          !          filename =TRIM(ADJUSTL(output_dir))//TRIM(output_prefix)//TRIM(ADJUSTL(test_case))//".Uavg."//TRIM(ADJUSTL(chartag))

          do ie=1,nelemd
             ierr=cube_assemble(gbl_fld,elem(ie)%accum%u,elem(ie),hybrid%par,nelemd,nelem,ie)
             if (ierr /= 0) call haltmp("error in prim_movie_output, halting...")
          end do
          call syncmp(hybrid%par)

          !          if (hybrid%par%masterproc) then
          !             call movie_output(ne,npts,quadtype,gbl_fld,gbl_x(:,:,:,1),gbl_y(:,:,:,1),filename)
          !          end if
          !          call syncmp(hybrid%par)
          !=======================================================================================================! 
          ! 3D writing: 												!
          !=======================================================================================================!        
          nx1 = ne*npts
          ny1 = ne*npts         
          if (hybrid%par%masterproc) then
             open (unit=37,file='./data/accum_vel.out')
             write(37,19)((((gbl_fld(i,j,k,l),i=1,nx1),j=1,ny1),k=1,6),l=1,nlev)
             close(37)
19           format(1x,5E16.6)
          end if
          call syncmp(hybrid%par)
          !=======================================================================================================!          
          deallocate(x)
          deallocate(y)
          deallocate(gbl_fld)
          deallocate(gbl_x)
          deallocate(gbl_y)
       end if

       !=======================================================================================================!
       if (accum_on) then
          naccum=naccum+1
       end if

    end if
    !$OMP BARRIER
  end subroutine dg_accum_output

  !=======================================================================================================!
  !   dg_accum_latlon
  !=======================================================================================================!
  subroutine dg_accum_latlon(elem,tl,hybrid)    
    use quadrature_mod, only : quadrature_t, gausslobatto
    implicit none    
    type (element_t) :: elem(:)
    type (TimeLevel_t)  :: tl
    type (hybrid_t)     :: hybrid
    ! Local variables  
    integer, parameter :: iunit=88     
    type (quadrature_t)  :: gll
    type (interpolate_t) :: interp 
    character(len=80):: filename   
    logical:: accum_on,accum_done
    integer:: i,j,k,ie,ierr,klon,klat,klev,naccum,simday
    integer:: accumday= accumfreq*secpday 
    real (kind=real_kind):: lat,lon,sum_global 
    real (kind=real_kind):: rnaccump1
    real (kind=real_kind):: v1,v2,epsi=1.0D-06    
    real (kind=real_kind):: local_cube(np,np,nlev)
    real (kind=real_kind):: global_cube(np*ne,np*ne,6,nlev)
    real (kind=real_kind), dimension(:,:,:), allocatable:: global_T,global_u
    real (kind=real_kind), dimension(:,:), allocatable  :: zonal_T,zonal_u 
    !=======================================================================================================!
    if (hybrid%ithr==0) then
       if ( ABS(Time_at(tl%nstep)/secpday-accumstart)<epsi ) then
          if (hybrid%par%masterproc) then
             write(*,'(a)')'Accumulation for held-suarez started ============================'
          endif
          naccum=0
          do ie=1,nelemd
             elem(ie)%accum%T=0.0D0
             elem(ie)%accum%u=0.0D0
          end do
       end if

       if ( ABS(Time_at(tl%nstep)/secpday-accumstop)<epsi ) then          
          if (hybrid%par%masterproc) then
             write(*,'(a)')'Accumulation for held-suarez finished ============================'
          endif
          accum_done=.true.
       else
          accum_done=.false.
       end if

       if( Time_at(tl%nstep) > accumstart*secpday  .and. &
            Time_at(tl%nstep) < accumstop*secpday   .and. &
            MODULO(tl%nstep,accumfreq)==0 ) then          
          if (hybrid%par%masterproc) then
             write(*,'(a,1x,i8,a,f10.2,a)')'Accumulation for held-suarez at:',tl%nstep,', time=',Time_at(tl%nstep)/secpday,' days'
          endif
          accum_on=.true.
          rnaccump1=1.0D0/REAL(naccum+1,kind=real_kind)
       else
          accum_on=.false.
       end if

       if (ndays > 0) then
          simday=NINT(Time_at(tl%nstep)/secpday)
       else
          simday=tl%nstep
       end if
       !=======================================================================================================!
       !   T and u average...
       !=======================================================================================================! 
       if (accum_on) then
          do ie=1,nelemd
             do k=1,nlev
                elem(ie)%accum%u(:,:,k)= rnaccump1*(elem(ie)%state%uv(:,:,1,k) + naccum*elem(ie)%accum%u(:,:,k))          
                elem(ie)%accum%T(:,:,k)= rnaccump1*(elem(ie)%state%t3d(:,:,k)+ naccum*elem(ie)%accum%T(:,:,k))
             end do
          end do
          call syncmp(hybrid%par)
       end if

       if (accum_on) then
          naccum=naccum+1
       end if
       !=======================================================================================================!
       !   lat-lon Output   
       !=======================================================================================================!
       if (accum_done) then
          gll=gausslobatto(np)
          call interpolate_create(gll,interp)    

          klon= (np-1)*ne*4 + 1
          klat= (np-1)*ne*2 + 1
          if (hybrid%par%masterproc .and. hybrid%ithr==0) then
             print *,'GLL coordinate'    
             do i=1,np
                print *,i,interp%glp(i)
             end do
             print *,'GLL intpolant'
             do j=1,np
                print *,j,(interp%Imat(i,j),i=1,np)
             enddo
             print *,'Global Longitude:',klon
             print *,'Global Latitude :',klat           
             print *,'Vertical Layers :',nlev                
          endif
          allocate(global_T(klon,klat,nlev))    
          allocate(global_u(klon,klat,nlev)) 
          allocate(zonal_T(klat,nlev))         
          allocate(zonal_u(klat,nlev))    

          call syncmp(hybrid%par)   
          do ie=1,nelemd
             do k= 1,nlev
                local_cube(:,:,k)= elem(ie)%accum%u(:,:,k)
             enddo
             ierr= cube_assemble(global_cube,local_cube,elem(ie),hybrid%par,nelemd,nelem,ie)
          enddo
          do k=1,nlev
             call polar_grid_assemble(global_cube(:,:,:,k),global_u(:,:,k),interp)  
          enddo
          do k=1,nlev    
             do j=1,klat
                sum_global= 0.0D0
                do i=1,klon
                   sum_global= sum_global + global_u(i,j,k)
                enddo
                zonal_u(j,k)= sum_global/klon
             enddo
          enddo

          call syncmp(hybrid%par)   
          do ie=1,nelemd
             do k= 1,nlev
                local_cube(:,:,k)= elem(ie)%accum%T(:,:,k)
             enddo
             ierr= cube_assemble(global_cube,local_cube,elem(ie),hybrid%par,nelemd,nelem,ie)
          enddo
          do k=1,nlev
             call polar_grid_assemble(global_cube(:,:,:,k),global_T(:,:,k),interp)  
          enddo
          do k=1,nlev    
             do j=1,klat
                sum_global= 0.0D0
                do i=1,klon
                   sum_global= sum_global + global_T(i,j,k)
                enddo
                zonal_T(j,k)= sum_global/klon
             enddo
          enddo
       end if
       !=======================================================================================================!
       if (accum_done) then
          if (hybrid%par%masterproc .and. hybrid%ithr==0) then
             open(iunit+0,file='./movies/zonal_T.dat',form="formatted",status='unknown')          
             open(iunit+1,file='./movies/zonal_u.dat',form="formatted",status='unknown')       
             do k=1,nlev    
                do j=1,klat
                   write(iunit+0,10)zonal_T(j,k) 
                   write(iunit+1,10)zonal_u(j,k)         
                enddo
             enddo
10           format(e22.15)   
             close(iunit+0)
             close(iunit+1)      
          end if

          call syncmp(hybrid%par)
          deallocate(global_T)    
          deallocate(global_u) 
          deallocate(zonal_T)         
          deallocate(zonal_u) 
       end if
    end if
    !=======================================================================================================!
  end subroutine dg_accum_latlon
#endif
#endif
  !=======================================================================================================!
end module dg_movie_mod


