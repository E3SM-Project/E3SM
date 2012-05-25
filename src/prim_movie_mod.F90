#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_movie_mod
#ifndef PIO_INTERP
  ! ---------------------
  use kinds, only : real_kind, longdouble_kind
  ! ---------------------
  use dimensions_mod, only :  nlev, nelem, nelemd, np, ne, nelemdmax, GlobalUniqueCols, nlevp, qsize, ntrac, nc
  ! ---------------------
  use hybvcoord_mod, only :  hvcoord_t 
  ! ---------------------
#ifdef _MPI
  use parallel_mod, only : syncmp, iam, mpireal_t, mpi_max, mpi_sum, mpiinteger_t, parallel_t, haltmp, abortmp
#else
  use parallel_mod, only : syncmp, iam, mpireal_t, parallel_t
#endif
  ! ---------------------
  use time_mod, only : Timelevel_t, tstep, ndays, time_at, secpday, nendstep,nmax
  ! ---------------------
  use element_mod, only : element_t
  use cslam_control_volume_mod, only : cslam_struct
  ! ---------------------
  use cube_mod, only : cube_assemble
  ! ---------------------
!  use prim_state_mod, only : naccum
  ! ---------------------
  use control_mod, only : test_case, runtype, accumstart, &
       accumstop, accumfreq, restartfreq, &
       integration, columnpackage, hypervis_power
  ! ---------------------
  use common_io_mod, only : &
       output_start_time,   &
       output_end_time,     &
       output_frequency,    &
       output_dir,          &
       max_output_variables,&
       max_output_streams,  &
       varname_len,         &
       nf_handle,           &
       get_current_varnames, &
       nfsizekind,             &
       nf_selectedvar
       
  use surfaces_mod, only : cvlist, InitControlVolumesData, InitControlVolumes

  use netcdf_io_mod, only:  nf_output_init_begin,& 
       nf_global_attribute, &
       nf_output_init_complete,  &
       nf_output_register_variables,&
       nf_put_var, &
       nf_close_all, &
       nf_output_register_dims, &
       nf_advance_frame, &
       nf_variable_attributes, &
       nf_get_frame

  ! ---------------------
  use coordinate_systems_mod, only : cartesian2D_t, spherical_polar_t, cartesian3D_t, spherical_to_cart
  ! ---------------------
  use physical_constants, only : g, kappa, p0, dd_pi
  ! ---------------------
  use aquaplanet_io_mod, only : aq_movie_init, aq_movie_output, aq_set_varnames
  ! ---------------------
  use physics_io_mod, only : physics_movie_init, physics_movie_output, physics_set_varnames
  ! ---------------------
  use dof_mod, only : UniquePoints, UniqueCoords, UniqueNcolsP, createmetadata
  ! ---------------------
  use pio, only  : io_desc_t !_EXTERNAL


    use hybrid_mod, only : hybrid_t, hybrid_create
    use edge_mod, only : EdgeBuffer_t

    use common_movie_mod, only: varrequired, vartype, varnames, varcnt, vardims, &
	dimnames, maxdims

  implicit none
  private
  public :: prim_movie_output, prim_movie_init, &
       prim_movie_finish,  nextoutputstep

  type(nf_handle), save :: ncdf(max_output_streams)
  integer, private :: nxyp
  integer(kind=nfsizekind) :: piostart2d, piocount2d, piocount3d(2), piostart3d(2)

contains

  subroutine GetDOF(elem, gcols, nz, compdof)

    type(element_t), intent(in) :: elem(:)
    integer, intent(in) :: gcols, nz
    integer, intent(out) :: compdof(:)
    integer :: k, i, ie, icnt

    icnt=0
    do k=1,nz
       do ie=1,nelemd
          do i=1,elem(ie)%idxp%NumUniquePts
             icnt=icnt+1
             compDOF(icnt)=elem(ie)%idxp%UniquePtOffset+i-1+(k-1)*GCols
          end do
       end do
    end do
  end subroutine GetDOF

  subroutine prim_movie_init(elem, par, hvcoord,tl)
    use hybvcoord_mod, only : hvcoord_t
    use parallel_mod, only : abortmp
    use pio, only : PIO_InitDecomp, pio_setdebuglevel, pio_int, pio_double, pio_closefile !_EXTERNAL
    use netcdf_io_mod, only : iodesc2d, iodesc3d, iodesc2d_nc, iodesc3d_nc, iodesc3d_subelem, iodesct, pio_subsystem 
    use common_io_mod, only : num_io_procs, num_agg, io_stride
#ifdef _CSLAM
    use cslam_control_volume_mod, only : cslam_mesh_ari
#endif
    type (element_t), intent(in) :: elem(:)
    type (parallel_t), intent(in)     :: par
    type (hvcoord_t), intent(in) :: hvcoord
    type(timelevel_t) :: tl
    ! Local variables
    type (cslam_struct) :: cslam_tmp
    type (hybrid_t) :: hybrid
    real (kind=real_kind),allocatable, dimension(:) :: latp,lonp
    integer :: ie, v1(4), i, ios, istartP
    integer,dimension(maxdims) :: dimsize
    integer :: st, en, icnt, kmax
    integer :: j,jj,cc,ii,k, iorank,base, global_nc, global_nsub
    integer(kind=nfsizekind) :: start(2), count(2)
    integer, allocatable :: compDOF(:)
    integer, allocatable :: dof(:)
    type(io_desc_t) :: iodescv, iodescvp1
    integer,allocatable  :: subelement_corners(:,:)
    real(kind=real_kind),allocatable  :: var1(:,:),var2(:,:)
    character(len=varname_len), pointer :: output_varnames(:)

    real (kind=real_kind) :: vartmp(np,np,nlev)
    real (kind=real_kind),allocatable :: var3d(:,:)


#ifdef _MPI
    integer :: ierr
#endif

    num_agg = 1
    call PIO_setDebugLevel(0)
#ifdef _AIX
!    call unbind()
#endif


    call nf_output_init_begin(ncdf,par%masterproc,par%nprocs,par%rank, &
         par%comm,test_case,runtype)


    nxyp=0
    do ie=1,nelemd
      nxyp=nxyp+elem(ie)%idxp%NumUniquePts
    enddo
    global_nc=nc*nc*nelem  ! total number of physics points
    global_nsub=(np-1)*(np-1)*nelem  ! total number of subelements
    dimsize = (/GlobalUniqueCols,nlev,nlevp,nelemd,0,global_nc,global_nsub/)
    call nf_output_register_dims(ncdf,maxdims, dimnames, dimsize)


    allocate(compdof(nxyp*nlev), latp(nxyp),lonp(nxyp))
    
    ! Create the DOF arrays for GLL points
    iorank=pio_subsystem%io_rank

    call getDOF(elem, GlobalUniqueCols, 1, compdof)
    call PIO_initDecomp(pio_subsystem, pio_double,(/GlobalUniqueCols/),&
         compDOF(1:nxyp),IOdesc2D)

    call getDOF(elem, GlobalUniqueCols, nlev, compdof)
    call PIO_initDecomp(pio_subsystem, pio_double,(/GlobalUniqueCols,nlev/),&
         compDOF,IOdesc3D)

    


! trivial case for vertical variables
    if(par%masterproc) then
       do k=1,nlevp
          compdof(k)=k
       end do
    else
       compdof=0
    end if
    call pio_initdecomp(pio_subsystem, pio_double, (/nlev/), compdof(1:nlev), iodescv)
    call pio_initdecomp(pio_subsystem, pio_double, (/nlevp/), compdof(1:nlevp), iodescvp1)
    
! this is a trivial case for the time variable
    if(iorank==0) then
       compdof(1)=1
    else		    	
       compdof(1)=0
    end if
    start=-1
    count=-1

    call PIO_initDecomp(pio_subsystem,pio_double,(/1/),&
         compDOF(1:1),IOdescT)

    deallocate(compdof)

! the CSLAM grid
    allocate(dof(nc*nc*nelemd*nlev))
    jj=0
    do cc=0,nlev-1
       do ie=1,nelemd
          base = ((elem(ie)%globalid-1)+cc*nelem)*(nc*nc)
          ii=0
          do j=1,nc
             do i=1,nc
                ii=ii+1
                jj=jj+1
                dof(jj) = base+ii
             end do
          end do
       end do
    end do
    call pio_initdecomp(pio_subsystem, pio_double, (/global_nc,nlev/), dof, iodesc3d_nc)
    call PIO_initDecomp(pio_subsystem, pio_double,(/global_nc/),&
         dof(1:(nelemd*nc*nc)),IOdesc2D_nc)
    deallocate(dof)

! the GLL based element subgrid
    allocate(dof((np-1)*(np-1)*nelemd*nlev))
    jj=0
    do cc=0,nlev-1
       do ie=1,nelemd
          base = ((elem(ie)%globalid-1)+cc*nelem)*(np-1)*(np-1)
          ii=0
          do j=1,np-1
             do i=1,np-1
                ii=ii+1
                jj=jj+1
                dof(jj) = base+ii
             end do
          end do
       end do
    end do
    call pio_initdecomp(pio_subsystem, pio_int, (/global_nsub,nlev/), dof, iodesc3d_subelem)
    deallocate(dof)


    if (par%masterproc) print *,'registering NETCDF variables'
    call nf_output_register_variables(ncdf,varcnt,varnames,vardims,vartype,varrequired)
    call nf_global_attribute(ncdf, 'np', np)
    call nf_global_attribute(ncdf, 'ne', ne)

    call nf_variable_attributes(ncdf, 'ps', 'surface pressure','pascals','coordinates','lat lon')
    call nf_variable_attributes(ncdf, 'area', 'area weights','radians^2','coordinates','lat lon')
    call nf_variable_attributes(ncdf, 'u', 'longitudinal wind component','meters/second')
    call nf_variable_attributes(ncdf, 'v', 'latitudinal wind component','meters/second')
    call nf_variable_attributes(ncdf, 'T', 'Temperature','degrees kelvin')
    call nf_variable_attributes(ncdf, 'lat', 'column latitude','degrees_north')
    call nf_variable_attributes(ncdf, 'lon', 'column longitude','degrees_east')
    call nf_variable_attributes(ncdf, 'phys_lat', 'column latitude','degrees_north')
    call nf_variable_attributes(ncdf, 'phys_lon', 'column longitude','degrees_east')
    call nf_variable_attributes(ncdf, 'phys_cv_lat', 'column control volume latitude','degrees_north')
    call nf_variable_attributes(ncdf, 'phys_cv_lon', 'column control volume longitude','degrees_east')
    call nf_variable_attributes(ncdf, 'time', 'Model elapsed time','days')
    call nf_variable_attributes(ncdf, 'lev' ,'hybrid level at midpoints' ,'level','positive','down') !,'formula_terms','a: hyam b: hybm p0: P0 ps: PS')
    call nf_variable_attributes(ncdf, 'ilev','hybrid level at interfaces','level','positive','down') !,'formula_terms','a: hyai b: hybi p0: P0 ps: PS')
    call nf_variable_attributes(ncdf, 'hyam','hybrid A coefficiet at layer midpoints' ,'dimensionless') 
    call nf_variable_attributes(ncdf, 'hybm','hybrid B coefficiet at layer midpoints' ,'dimensionless') 
    call nf_variable_attributes(ncdf, 'hyai','hybrid A coefficiet at layer interfaces' ,'dimensionless') 
    call nf_variable_attributes(ncdf, 'hybi','hybrid B coefficiet at layer interfaces' ,'dimensionless') 


    if(test_case.eq.'aquaplanet') then
       call aq_movie_init(ncdf)
    end if
    if(columnpackage.ne.'none') then
       call physics_movie_init(ncdf)
    end if

    call nf_output_init_complete(ncdf)


    do ios=1,max_output_streams
       output_varnames=>get_current_varnames(ios)
       if( (nf_selectedvar('cv_lat', output_varnames)) .or.  &
            (nf_selectedvar('cv_lon', output_varnames)) ) then
          if (.not. allocated(cvlist)) then
             if (par%masterproc) print *,'computing GLL dual grid for  control volumes:'
             call InitControlVolumesData(nelemd)
             ! single thread
             hybrid = hybrid_create(par,0,1)
             call InitControlVolumes(elem,hybrid,1,nelemd)
             if (par%masterproc) print *,'done.'
          endif
       endif

       if( (nf_selectedvar('corners', output_varnames))) then
          if (par%masterproc) print *,'writing subelement metadata'
          allocate(subelement_corners((np-1)*(np-1)*nelemd,nlev))
          subelement_corners=0
          call createmetadata(par, elem, subelement_corners(:,1:4))
          call nf_put_var(ncdf(ios),subelement_corners,start(1:1),count(1:1),&
               name='corners',iodescin=iodesc3d_subelem)
          deallocate(subelement_corners)
       endif
    enddo



    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0) ) then
          st=1
          if (par%masterproc) print *,'writing coordinates ios=',ios
          do ie=1,nelemdmax
            ! if (par%masterproc .and. mod(ie,1).eq.0 ) print *,'ie=',ie
	    if(ie<=nelemd) then
               en=st+elem(ie)%idxp%NumUniquePts-1
               call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep,latp(st:en), lonp(st:en)) 
               st=en+1
            end if
          enddo

          latp=latp*180/dd_pi
          lonp=lonp*180/dd_pi
          call nf_put_var(ncdf(ios),latp,start(1:1),count(1:1),name='lat', iodescin=iodesc2d)
          call nf_put_var(ncdf(ios),lonp,start(1:1),count(1:1),name='lon', iodescin=iodesc2d)

          start(1)=1
          count(1)=nlev
          call nf_put_var(ncdf(ios),hvcoord%hyam,start(1:1),count(1:1),name='hyam', iodescin=iodescv)
          call nf_put_var(ncdf(ios),hvcoord%hybm,start(1:1),count(1:1),name='hybm', iodescin=iodescv)
          call nf_put_var(ncdf(ios),hvcoord%etam,start(1:1),count(1:1),name='lev', iodescin=iodescv)
          start(1)=1
          count(1)=nlevp
          call nf_put_var(ncdf(ios),hvcoord%hyai,start(1:1),count(1:1),name='hyai',iodescin=iodescvp1)
          call nf_put_var(ncdf(ios),hvcoord%hybi,start(1:1),count(1:1),name='hybi',iodescin=iodescvp1)
          call nf_put_var(ncdf(ios),hvcoord%etai,start(1:1),count(1:1),name='ilev',iodescin=iodescvp1)


          ! check if reqested 




          output_varnames=>get_current_varnames(ios)
          kmax = 0
          if(nf_selectedvar('cv_lon',output_varnames).or. &
               nf_selectedvar('cv_lat',output_varnames).or. &
               nf_selectedvar('faceno',output_varnames)) then
             do ie=1,nelemd
                kmax = MAX(kmax,MAXVAL(cvlist(ie)%nvert))
             enddo
             if ( nlev < kmax ) call abortmp('cv output requires nlev >= max number of vertex')
          endif
          
          if(nf_selectedvar('cv_lon', output_varnames)) then
             if (par%masterproc) print *,'writing cv_lon...'
             allocate(var3d(nxyp,nlev))
             st=1
             do ie=1,nelemd
                en=st+elem(ie)%idxp%NumUniquePts-1
                vartmp = 0
                do k=1,kmax
                   vartmp(:,:,k)=cvlist(ie)%vert_latlon(k,:,:)%lon*180/dd_pi
                enddo
                call UniquePoints(elem(ie)%idxp, nlev, vartmp, var3d(st:en,:))
                st=en+1
             enddo
             Call nf_put_var(ncdf(ios),var3d,start, count, name='cv_lon')
             deallocate(var3d)
          end if
          if(nf_selectedvar('cv_lat', output_varnames)) then
             if (par%masterproc) print *,'writing cv_lat...'
             allocate(var3d(nxyp,nlev))
             st=1
             do ie=1,nelemd
                en=st+elem(ie)%idxp%NumUniquePts-1
                vartmp = 0
                do k=1,kmax
                   vartmp(:,:,k)=cvlist(ie)%vert_latlon(k,:,:)%lat*180/dd_pi
                enddo
                call UniquePoints(elem(ie)%idxp,nlev,vartmp,var3d(st:en,:))
                st=en+1
             enddo
             call nf_put_var(ncdf(ios),var3d,start, count, name='cv_lat')
             deallocate(var3d)
          end if
          if(nf_selectedvar('faceno', output_varnames)) then
             if (par%masterproc) print *,'writing face_no...'
             allocate(var3d(nxyp,nlev))
             ! also output face_no
             st=1
             do ie=1,nelemd
                en=st+elem(ie)%idxp%NumUniquePts-1
                vartmp = 0
                do k=1,kmax
                   vartmp(:,:,k)=cvlist(ie)%face_no(k,:,:)
                enddo
                call UniquePoints(elem(ie)%idxp,nlev,vartmp,var3d(st:en,:))
                st=en+1
             enddo
             call nf_put_var(ncdf(ios),var3d,start, count, name='faceno')
             deallocate(var3d)
          end if

#ifdef _CSLAM          
          allocate(var1(nc*nc*nelemd,nlev))
          allocate(var2(nc*nc*nelemd,nlev))
          if( (nf_selectedvar('phys_lat', output_varnames))) then
             jj=0
             do ie=1,nelemd
                ! cslam grid for element ie.
                ! (if ntrac=0, cslam grid was never computed, so we cant use cslam()%
                call cslam_mesh_ari(elem(ie),cslam_tmp,tl)
                ii=0
                do j=1,nc
                   do i=1,nc
                      ii=ii+1
                      jj=jj+1
                      ! take center of mass of elem(ie)%
                      var1(jj,1) = cslam_tmp%centersphere(i,j)%lat
                      var2(jj,1) = cslam_tmp%centersphere(i,j)%lon
                   end do
                end do
             end do

             var1=var1*180/dd_pi
             var2=var2*180/dd_pi
             if (par%masterproc) print *,'writing phys_lat'
             call nf_put_var(ncdf(ios),var1(:,1),start(1:1),count(1:1),&
                  name='phys_lat',iodescin=iodesc2d_nc)
             if (par%masterproc) print *,'writing phys_lon'
             call nf_put_var(ncdf(ios),var2(:,1),start(1:1),count(1:1),&
                  name='phys_lon',iodescin=iodesc2d_nc)
          endif

          if( (nf_selectedvar('phys_area', output_varnames))) then
             jj=0
             do ie=1,nelemd

                ! cslam grid for element ie.
                ! (if ntrac=0, cslam grid was never computed, so we cant use cslam()%
                call cslam_mesh_ari(elem(ie),cslam_tmp,tl)

                ii=0
                do j=1,nc
                   do i=1,nc
                      ii=ii+1
                      jj=jj+1
                      var1(jj,1)= cslam_tmp%area_sphere(i,j)
                   end do
                end do
             end do
             call nf_put_var(ncdf(ios),var1(:,1),start(1:1),count(1:1),&
                  name='phys_area',iodescin=iodesc2d_nc)
          endif

          if( (nf_selectedvar('phys_cv_lat', output_varnames))) then
             if (par%masterproc) print *,'writing physics grid cv metadata'
             var1=0
             var2=0
             jj=0
             do ie=1,nelemd

                ! cslam grid for element ie.
                ! (if ntrac=0, cslam grid was never computed, so we cant use cslam()%
                call cslam_mesh_ari(elem(ie),cslam_tmp,tl)

                ii=0
                do j=1,nc
                   do i=1,nc
                      ii=ii+1
                      jj=jj+1
                      ! take center of mass of elem(ie)%
                      var1(jj,1) = cslam_tmp%asphere(i,j)%lat
                      var2(jj,1) = cslam_tmp%asphere(i,j)%lon
                      var1(jj,2) = cslam_tmp%asphere(i+1,j)%lat
                      var2(jj,2) = cslam_tmp%asphere(i+1,j)%lon
                      var1(jj,3) = cslam_tmp%asphere(i+1,j+1)%lat
                      var2(jj,3) = cslam_tmp%asphere(i+1,j+1)%lon
                      var1(jj,4) = cslam_tmp%asphere(i,j+1)%lat
                      var2(jj,4) = cslam_tmp%asphere(i,j+1)%lon
                   end do
                end do
             end do
             var1=var1*180/dd_pi
             var2=var2*180/dd_pi
             call nf_put_var(ncdf(ios),var1(:,:),start(1:1),count(1:1),&
                  name='phys_cv_lat',iodescin=iodesc3d_nc)
             call nf_put_var(ncdf(ios),var2(:,:),start(1:1),count(1:1),&
                  name='phys_cv_lon',iodescin=iodesc3d_nc)
                
          endif
          deallocate(var1)
          deallocate(var2)
#else
          if( (nf_selectedvar('phys_lat', output_varnames))) then
             if (par%masterproc) print *,'WARNING: compile with -D_CSLAM to output CSLAM grid coordinate data'
          endif
#endif
          if (par%masterproc) print *,'done writing coordinates ios=',ios
       end if
    end do

    
    deallocate(latp,lonp)
#ifdef _AIX
!    call rebind()
#endif

  end subroutine prim_movie_init
  subroutine prim_movie_finish
! ncdf is a module global
    call nf_close_all(ncdf)
   
   
  end subroutine prim_movie_finish
!
! This function returns the next step number in which an output (either restart or movie) 
! needs to be written.
!
  integer function nextoutputstep(tl)
    type(timelevel_t), intent(in) :: tl
    integer :: ios, nstep(max_output_streams)

    nstep(:) = nEndStep
    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0)) then
          if ((output_start_time(ios) .le. tl%nstep) .and. &
               (output_end_time(ios) .ge. tl%nstep)) then
             nstep(ios)=nstep(ios)+ output_frequency(ios) - &
                  MODULO(nstep(ios),output_frequency(ios))
          end if
       end if
    end do
    nextoutputstep=minval(nstep)
    if(restartfreq>0) then
       nextoutputstep=min(nextoutputstep,tl%nstep+restartfreq-MODULO(tl%nstep,restartfreq))    
    end if
 end function nextoutputstep

  subroutine prim_movie_output(elem, tl, hvcoord, hybrid, nets,nete, cslam)
    use piolib_mod, only : Pio_SetDebugLevel !_EXTERNAL
    use perf_mod, only : t_startf, t_stopf !_EXTERNAL
    use viscosity_mod, only : compute_zeta_C0
    use netcdf_io_mod, only : iodesc3d_nc

    type (element_t)    :: elem(:)
    type (cslam_struct), optional   :: cslam(:)
    type (TimeLevel_t)  :: tl
    type (hvcoord_t)    :: hvcoord
    type (hybrid_t)      , intent(in) :: hybrid
    integer              :: nets,nete

    real*8              :: st_write, et_write, dt_write, dt_write_global
    character(len=varname_len), pointer :: output_varnames(:)
    integer :: ie,ios, i, j, k,jj
    real (kind=real_kind) :: pfull, pr0
    real(kind=real_kind),parameter :: dayspersec=1d0/(3600.*24.)
    real (kind=real_kind) :: vartmp(np,np,nlev),arealocal(np,np)
    real (kind=real_kind) :: var2d(nxyp), var3d(nxyp,nlev), ke(np,np,nlev)
    real (kind=real_kind) :: temp3d(np,np,nlev,nets:nete)
    real (kind=real_kind) :: varphys(nc*nc*nelemd,nlev)

    integer :: st, en, kmax, qindex
    character(len=2) :: vname

    integer(kind=nfsizekind) :: start(3), count(3), start2d(2),count2d(2)
    integer :: ncnt
    call t_startf('prim_movie_output:pio')
#ifdef _AIX
!    call unbind()
#endif


    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0)) then
          if ((output_start_time(ios) .le. tl%nstep) .and. &
               (output_end_time(ios) .ge. tl%nstep) .and. &
               MODULO(tl%nstep,output_frequency(ios)) .eq. 0) then
             output_varnames=>get_current_varnames(ios)
             start2d(1)=piostart2d
             start2d(2)=nf_get_frame(ncdf(ios))
             count2d(1)=piocount2d
             count2d(2)=1

             count(1:2)=piocount3d
             start(1:2)=piostart3d
             start(3)=nf_get_frame(ncdf(ios))
             count(3)=1

             if(nf_selectedvar('ps', output_varnames)) then
                if (hybrid%masterthread) print *,'writing ps...'
                st=1
                do ie=1,nelemd
                   vartmp(:,:,1) = exp(elem(ie)%state%lnps(:,:,tl%n0))
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxP,vartmp(:,:,1),var2d(st:en))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var2d,start2d,count2d,name='ps')
             endif

             if(nf_selectedvar('hypervis', output_varnames)) then
                if (hybrid%masterthread) print *,'writing hypervis...'
                st=1
                do ie=1,nelemd
                   vartmp(:,:,1) = elem(ie)%variable_hyperviscosity(:,:)
                   ! scale back to a length scale
                   if (hypervis_power /= 0 ) vartmp(:,:,1)=vartmp(:,:,1)**(2d0/hypervis_power)
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxP,vartmp(:,:,1),var2d(st:en))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var2d,start2d,count2d,name='hypervis')
             endif


             if(nf_selectedvar('geos', output_varnames)) then
                if (hybrid%masterthread) print *,'writing geos...'
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxP,elem(ie)%state%phis,var2d(st:en))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var2d,start2d,count2d,name='geos')
             endif

             if(nf_selectedvar('area', output_varnames)) then
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   arealocal=1/elem(ie)%rspheremp(:,:)         
                   call UniquePoints(elem(ie)%idxP,arealocal,var2d(st:en))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var2d,start2d,count2d,name='area')
             endif


             if(nf_selectedvar('zeta', output_varnames)) then
                if (hybrid%masterthread) print *,'writing zeta...'
                ! velocities are on sphere for primitive equations
                call compute_zeta_C0(temp3d,elem,hybrid,nets,nete,tl%n0)

                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxp,nlev,temp3d(:,:,:,ie),var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='zeta')
             end if


             if(nf_selectedvar('T', output_varnames)) then
                if (hybrid%masterthread) print *,'writing T...'
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxP,nlev,elem(ie)%state%T(:,:,:,tl%n0),var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='T')
             end if


             if(nf_selectedvar('Th', output_varnames)) then
                pr0=1./(p0)
                st=1
                do ie=1,nelemd
                   do k=1,nlev
                      do j=1,np
                         do i=1,np
                            pfull = hvcoord%hyam(k)*hvcoord%ps0  &
                                 + hvcoord%hybm(k)*exp(elem(ie)%state%lnps(i,j,tl%n0))
                            varTMP(i,j,k)=elem(ie)%state%T(i,j,k,tl%n0)* &
                                 (pfull*pr0)**(-kappa)
                         end do
                      end do
                   end do
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxP,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                end do
                call nf_put_var(ncdf(ios),var3d,start, count, name='Th')
             end if
             if(nf_selectedvar('u', output_varnames)) then
                if (hybrid%masterthread) print *,'writing u...'
                st=1
                do ie=1,nelemd
                   do k=1,nlev
                      vartmp(:,:,k) = elem(ie)%state%v(:,:,1,k,tl%n0)
                   end do
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxv,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='u')
             end if

             if(nf_selectedvar('v', output_varnames)) then
                if (hybrid%masterthread) print *,'writing v...'
                st=1
                do ie=1,nelemd
                   do k=1,nlev
                      vartmp(:,:,k) = elem(ie)%state%v(:,:,2,k,tl%n0)
                   end do
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxv,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='v')
             end if

             if(nf_selectedvar('ke', output_varnames)) then
                st=1
                do ie=1,nelemd
                   do k=1,nlev 
                      ke(:,:,k) = (elem(ie)%state%v(:,:,1,k,tl%n0)**2 + &
        	       elem(ie)%state%v(:,:,2,k,tl%n0)**2 )/2
                   enddo
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxp,nlev,ke, var3d(st:en,:))
                   st=en+1
                end do
                call nf_put_var(ncdf(ios),var3d,start, count, name='ke')
             end if

             do qindex=1,min(qsize,4)
                write(vname,'(a1,i1)') 'Q',qindex
                if (qindex==1) vname='Q'
                if(nf_selectedvar(vname, output_varnames)) then
                   if (hybrid%masterthread) print *,'writing ',vname
                   st=1
                   do ie=1,nelemd
                      en=st+elem(ie)%idxp%NumUniquePts-1
                      call UniquePoints(elem(ie)%idxP,nlev,elem(ie)%state%Q(:,:,:,qindex,tl%n0), var3d(st:en,:))
                      st=en+1
                   end do
                   call nf_put_var(ncdf(ios),var3d,start, count, name=vname)
                end if
             enddo

             if (present(cslam)) then
             do qindex=1,min(ntrac,4)
                write(vname,'(a1,i1)') 'C',qindex
                if (qindex==1) vname='C'
                if(nf_selectedvar(vname, output_varnames)) then

                   
                   if (hybrid%masterthread) print *,'writing ',vname
                   do k=1,nlev
                      jj=0
                      do ie=1,nelemd
                         do j=1,nc
                            do i=1,nc
                               jj=jj+1
                               varphys(jj,k)= cslam(ie)%c(i,j,k,qindex,tl%n0)
                            end do
                         end do
                      end do
                   end do
                   call nf_put_var(ncdf(ios),varphys,start, count, name=vname,iodescin=iodesc3d_nc)
                endif
             enddo
             endif

             if(nf_selectedvar('geo', output_varnames)) then
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxP,nlev,elem(ie)%derived%phi,var3d(st:en,:))
                   st=en+1
                end do
                call nf_put_var(ncdf(ios),var3d,start, count, name='geo')
             end if

             if(nf_selectedvar('omega', output_varnames)) then
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxp,nlev,elem(ie)%derived%omega_p,var3d(st:en,:))
                   st=en+1
                end do
                call nf_put_var(ncdf(ios),var3d,start, count, name='omega')
             end if
             if(test_case.eq.'aquaplanet') then
                call aq_movie_output(ncdf(ios), elem, output_varnames,nxyp, nlev)
             end if
             if(columnpackage.ne.'none') then
                call physics_movie_output(ncdf(ios),elem, output_varnames, nxyp)
             end if
!             call PIO_SetDebugLevel(3)
             call nf_put_var(ncdf(ios),real(dayspersec*time_at(tl%nstep),kind=real_kind),&
                  start(3:3),count(3:3),name='time')
             call nf_advance_frame(ncdf(ios))
!             call PIO_SetDebugLevel(0)

          end if
       end if
    end do
#ifdef _AIX
!    call rebind()
#endif
    call t_stopf('prim_movie_output:pio')
 end subroutine prim_movie_output

#endif
end module prim_movie_mod

