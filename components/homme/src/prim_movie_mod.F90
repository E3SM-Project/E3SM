! Sept 2019 O. Guba Add w_i, mu_i, geo_i, pnh to native output

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_movie_mod
#ifndef HOMME_WITHOUT_PIOLIBRARY
#ifndef PIO_INTERP
  use kinds, only : real_kind, longdouble_kind
  use dimensions_mod, only :  nlev, nelem, nelemd, np, ne, nelemdmax, GlobalUniqueCols, nlevp, qsize
  use hybvcoord_mod, only :  hvcoord_t
#ifdef _MPI
  use parallel_mod, only : syncmp, iam, mpireal_t, mpi_max, mpi_sum, mpiinteger_t, parallel_t, haltmp, abortmp
#else
  use parallel_mod, only : syncmp, iam, mpireal_t, parallel_t
#endif
  use time_mod, only : Timelevel_t, tstep, ndays, time_at, secpday, nendstep,nmax, Timelevel_Qdp
  use element_mod, only : element_t

  use cube_mod, only : cube_assemble
  use control_mod, only : test_case, runtype, geometry, &
       restartfreq, &
       integration, hypervis_power, qsplit
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
       nf_selectedvar,       &
       PIOFS

  use surfaces_mod, only : cvlist, InitControlVolumesData, InitControlVolumes

  use netcdf_io_mod, only:  nf_output_init_begin,&
       nf_global_attribute, &
       nf_output_init_complete,  &
       nf_output_register_variables,&
       nf_put_var => nf_put_var_netcdf, &
       nf_close_all, &
       nf_output_register_dims, &
       nf_advance_frame, &
       nf_variable_attributes, &
       nf_get_frame

  use coordinate_systems_mod, only : cartesian2D_t, spherical_polar_t, cartesian3D_t, spherical_to_cart
  use physical_constants, only : g, kappa, p0, dd_pi
  use dof_mod, only : UniquePoints, UniqueCoords, UniqueNcolsP, createmetadata
#ifndef HOMME_WITHOUT_PIOLIBRARY
  use pio, only  : io_desc_t, pio_iotask_rank !_EXTERNAL
#endif

    use hybrid_mod, only : hybrid_t, hybrid_create
    use edgetype_mod, only : EdgeBuffer_t

    use common_movie_mod, only: varrequired, vartype, varnames, varcnt, vardims, &
                                dimnames, maxdims

  implicit none
  private
  save
  public :: prim_movie_output, prim_movie_init, &
       prim_movie_finish,  nextoutputstep

  type(nf_handle) :: ncdf(max_output_streams)
  integer, private :: nxyp

contains

  subroutine GetDOF(elem, gcols, nz, compdof)

    type(element_t), intent(in) :: elem(:)
    integer, intent(in) :: gcols, nz
    integer*8, intent(out) :: compdof(:)
    integer :: k, i, ie
    integer*8 :: icnt,gcols8

    icnt=0
    do k=1,nz
       do ie=1,nelemd
          do i=1,elem(ie)%idxp%NumUniquePts
             icnt=icnt+1
             if (icnt < 0 ) call abortmp('ERROR: UniquePts integer overflow')
             Gcols8=Gcols ! force I*8 calcluation
             compDOF(icnt)=elem(ie)%idxp%UniquePtOffset+i-1+(k-1)*GCols8
             if (compDOF(icnt) < 0 ) call abortmp('ERROR: compDOF integer overflow')
          end do
       end do
    end do
  end subroutine GetDOF

  subroutine prim_movie_init(elem, par, hvcoord,tl)
    use hybvcoord_mod, only : hvcoord_t
    use parallel_mod, only : abortmp
    use pio, only : PIO_InitDecomp, pio_setdebuglevel, pio_int, pio_double, pio_closefile !_EXTERNAL
    use netcdf_io_mod, only : iodesc2d, iodesc3d, iodesc3d_subelem, iodesct, iodesc3dp1
    use common_io_mod, only : num_io_procs, num_agg, io_stride
    use reduction_mod, only : parallelmax
    type (element_t), intent(in) :: elem(:)
    type (parallel_t), intent(in)     :: par
    type (hvcoord_t), intent(in) :: hvcoord
    type(timelevel_t) :: tl
    ! Local variables
    type (hybrid_t) :: hybrid
    real (kind=real_kind),allocatable, dimension(:) :: latp,lonp
    integer :: ie, v1(4), i, ios, istartP
    integer,dimension(maxdims) :: dimsize
    integer :: st, en, icnt, kmax,kmax2
    integer :: j,cc,k, iorank, global_nsub
    integer*8 :: ii,jj,base
    integer(kind=nfsizekind) :: start(2), count(2)
    integer*8, allocatable :: compDOF(:)
    integer*8, allocatable :: compDOFp1(:)
    type(io_desc_t) :: iodescv, iodescvp1
    integer,allocatable  :: subelement_corners(:,:)
    real(kind=real_kind),allocatable  :: var1(:,:),var2(:,:)
    character(len=varname_len), pointer :: output_varnames(:)

    real (kind=real_kind) :: vartmp(np,np,nlev)
    real (kind=real_kind),allocatable :: var3d(:,:),var2d(:)

    real (kind=real_kind) :: latlon_adj_factor
#ifdef _MPI
    integer :: ierr
#endif

    call PIO_setDebugLevel(0)

    if (geometry=="sphere") then
      latlon_adj_factor = 180.0D0/dd_pi
    else if (geometry=="plane") then
      latlon_adj_factor = 1.0D0
    end if

    call nf_output_init_begin(ncdf,par%masterproc,par%nprocs,par%rank, &
         par%comm,test_case,runtype)


    nxyp=0
    do ie=1,nelemd
      nxyp=nxyp+elem(ie)%idxp%NumUniquePts
    enddo
    global_nsub=(np-1)*(np-1)*nelem  ! total number of subelements
    dimsize = (/GlobalUniqueCols,nlev,nlevp,nelem,0,global_nsub/)
    call nf_output_register_dims(ncdf,maxdims, dimnames, dimsize)


    allocate(compdof(nxyp*nlevp), latp(nxyp),lonp(nxyp))

    ! Create the DOF arrays for GLL points
    iorank=pio_iotask_rank(PIOFS)

    if (par%masterproc) print *,'compDOF for 2d'
    call getDOF(elem, GlobalUniqueCols, 1, compdof)
    call PIO_initDecomp(PIOFS, pio_double,(/GlobalUniqueCols/),&
         compDOF(1:nxyp),IOdesc2D)

    if (par%masterproc) print *,'compDOF for 3d nlev'
    call getDOF(elem, GlobalUniqueCols, nlev, compdof)
    call PIO_initDecomp(PIOFS, pio_double,(/GlobalUniqueCols,nlev/),&
         compDOF(1:nxyp*nlev),IOdesc3D)

    if (par%masterproc) print *,'compDOF for 3d nlevp'
    call getDOF(elem, GlobalUniqueCols, nlevp, compdof)
    call PIO_initDecomp(PIOFS,pio_double,(/GlobalUniqueCols,nlevp/),&
         compDOF,iodesc3dp1)

! trivial case for vertical variables
    if(par%masterproc) then
       do k=1,nlevp
          compdof(k)=k
       end do
    else
       compdof=0
    end if
    call pio_initdecomp(PIOFS, pio_double, (/nlev/), compdof(1:nlev), iodescv)
    call pio_initdecomp(PIOFS, pio_double, (/nlevp/), compdof(1:nlevp), iodescvp1)

! this is a trivial case for the time variable
    if(iorank==0) then
       compdof(1)=1
    else
       compdof(1)=0
    end if
    start=-1
    count=-1

    call PIO_initDecomp(PIOFS,pio_double,(/1/),&
         compDOF(1:1),IOdescT)

    deallocate(compdof)


! the GLL based element subgrid
    if (par%masterproc) print *,'compDOF for subcell '
    allocate(compdof((np-1)*(np-1)*nelemd*nlev))
    jj=0
    do cc=0,nlev-1
       do ie=1,nelemd
          base = (INT(elem(ie)%globalid-1,8)+INT(cc,8)*nelem)*(np-1)*(np-1)
          ii=0
          do j=1,np-1
             do i=1,np-1
                ii=ii+1
                jj=jj+1
                compdof(jj) = base+ii
                if (compdof(jj)<0) call abortmp('ERROR: compDOF subcell integer overflow')
             end do
          end do
       end do
    end do
    call pio_initdecomp(PIOFS, pio_int, (/global_nsub,nlev/), compdof, iodesc3d_subelem)
    deallocate(compdof)


    if (par%masterproc) print *,'registering NETCDF variables'
    call nf_output_register_variables(ncdf,varcnt,varnames,vardims,vartype,varrequired)
    call nf_global_attribute(ncdf, 'np', np)
    call nf_global_attribute(ncdf, 'ne', ne)

  if (geometry=="sphere") then
    call nf_variable_attributes(ncdf, 'ps', 'surface pressure','Pa','coordinates','lat lon')
    call nf_variable_attributes(ncdf, 'area', 'area weights','radians^2','coordinates','lat lon')
    call nf_variable_attributes(ncdf, 'u', 'longitudinal wind component','meters/second')
    call nf_variable_attributes(ncdf, 'v', 'latitudinal wind component','meters/second')
  else if (geometry=="plane") then
    call nf_variable_attributes(ncdf, 'ps', 'surface pressure','Pa','coordinates','x y')
    call nf_variable_attributes(ncdf, 'area', 'area weights','m^2','coordinates','x y')
    call nf_variable_attributes(ncdf, 'u', 'x-dir wind component','meters/second')
    call nf_variable_attributes(ncdf, 'v', 'y-dir wind component','meters/second')
  end if

    call nf_variable_attributes(ncdf, 'T', 'Temperature','degrees kelvin')
    call nf_variable_attributes(ncdf, 'Th','potential temperature \theta','degrees kelvin')
    call nf_variable_attributes(ncdf, 'w', 'vertical wind component','meters/second')
    call nf_variable_attributes(ncdf, 'w_i',  'vertical wind component on interfaces','meters/second')
    call nf_variable_attributes(ncdf, 'mu_i', 'mu=dp/d\pi on interfaces','dimensionless')
    call nf_variable_attributes(ncdf, 'geo_i','geopotential on interfaces','meters')
    call nf_variable_attributes(ncdf, 'pnh',  'total pressure','Pa')
#ifdef _PRIM
    call nf_variable_attributes(ncdf, 'geos', 'surface geopotential','m^2/s^2')
    call nf_variable_attributes(ncdf, 'PHIS', 'surface geopotential','m^2/s^2')
    call nf_variable_attributes(ncdf, 'precl','Precipitation rate','meters of water/s')
#endif
  if (geometry=="sphere") then
    call nf_variable_attributes(ncdf, 'lat', 'column latitude','degrees_north')
    call nf_variable_attributes(ncdf, 'lon', 'column longitude','degrees_east')
  else if (geometry=="plane") then
    call nf_variable_attributes(ncdf, 'lat', 'column y','m')
    call nf_variable_attributes(ncdf, 'lon', 'column x','m')
  end if

    call nf_variable_attributes(ncdf, 'time', 'Model elapsed time','days')
    call nf_variable_attributes(ncdf, 'lev' ,'hybrid level at midpoints' ,'level','positive','down') !,'formula_terms','a: hyam b: hybm p0: P0 ps: PS')
    call nf_variable_attributes(ncdf, 'ilev','hybrid level at interfaces','level','positive','down') !,'formula_terms','a: hyai b: hybi p0: P0 ps: PS')
    call nf_variable_attributes(ncdf, 'hyam','hybrid A coefficiet at layer midpoints' ,'dimensionless')
    call nf_variable_attributes(ncdf, 'hybm','hybrid B coefficiet at layer midpoints' ,'dimensionless')
    call nf_variable_attributes(ncdf, 'hyai','hybrid A coefficiet at layer interfaces' ,'dimensionless')
    call nf_variable_attributes(ncdf, 'hybi','hybrid B coefficiet at layer interfaces' ,'dimensionless')
    call nf_output_init_complete(ncdf)

    do ios=1,max_output_streams
       output_varnames=>get_current_varnames(ios)
       if( (nf_selectedvar('cv_lat', output_varnames)) .or.  &
            (nf_selectedvar('cv_lon', output_varnames)) ) then
          if (.not. allocated(cvlist)) then
             if (par%masterproc) print *,'computing GLL dual grid for  control volumes:'
             call InitControlVolumesData(par,elem,nelemd)
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

          ! MAKE THIS AN ADJUSTMENT FACTOR!
          latp=latp*latlon_adj_factor
          lonp=lonp*latlon_adj_factor
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
             kmax2=0
             do ie=1,nelemd
                kmax2 = MAX(kmax2,MAXVAL(cvlist(ie)%nvert))
             enddo
             kmax = ParallelMax(kmax2,hybrid)
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
                   vartmp(:,:,k)=cvlist(ie)%vert_latlon(k,:,:)%lon*latlon_adj_factor
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
                   vartmp(:,:,k)=cvlist(ie)%vert_latlon(k,:,:)%lat*latlon_adj_factor
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

          if(nf_selectedvar('geos', output_varnames)) then
             allocate(var2d(nxyp))
             if (par%masterproc) print *,'writing geos...'
             st=1
             do ie=1,nelemd
                en=st+elem(ie)%idxp%NumUniquePts-1
                call UniquePoints(elem(ie)%idxP,elem(ie)%state%phis,var2d(st:en))
                st=en+1
             enddo
             call nf_put_var(ncdf(ios),var2d,start,count,name='geos')
             deallocate(var2d)
          endif
          if(nf_selectedvar('PHIS', output_varnames)) then
             allocate(var2d(nxyp))
             if (par%masterproc) print *,'writing geos as PHIS...'
             st=1
             do ie=1,nelemd
                en=st+elem(ie)%idxp%NumUniquePts-1
                call UniquePoints(elem(ie)%idxP,elem(ie)%state%phis,var2d(st:en))
                st=en+1
             enddo
             call nf_put_var(ncdf(ios),var2d,start,count,name='PHIS')
             deallocate(var2d)
          endif

          if (par%masterproc) print *,'done writing coordinates ios=',ios
       end if
    end do
    deallocate(latp,lonp)
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

  subroutine prim_movie_output(elem, tl, hvcoord, par)
    use piolib_mod, only : Pio_SetDebugLevel !_EXTERNAL
    use pio, only : pio_syncfile !_EXTERNAL
    use perf_mod, only : t_startf, t_stopf !_EXTERNAL
    use viscosity_mod, only : compute_zeta_C0
    use element_ops, only : get_field, get_field_i
    use dcmip16_wrapper, only: precl
    use netcdf_io_mod, only : iodesc3dp1

    type (element_t)    :: elem(:)

    type (TimeLevel_t)  :: tl
    type (hvcoord_t)    :: hvcoord
    type (parallel_t)   :: par

    real*8              :: st_write, et_write, dt_write, dt_write_global
    character(len=varname_len), pointer :: output_varnames(:)
    integer :: ie,ios, i, j, k,jj
    real (kind=real_kind) :: pfull, pr0
    real(kind=real_kind),parameter :: dayspersec=1d0/(3600.*24.)
    real (kind=real_kind) :: vartmp(np,np,nlev), vartmp1(np,np,nlevp), arealocal(np,np)
    real (kind=real_kind) :: var2d(nxyp), var3d(nxyp,nlev), var3dp1(nxyp,nlevp), ke(np,np,nlev)
    real (kind=real_kind) :: temp3d(np,np,nlev,nelemd)

    integer :: st, en, kmax, qindex, n0, n0_Q
    character(len=2) :: vname

    integer(kind=nfsizekind) :: start(3), count(3), start2d(2),count2d(2), &
                                startp1(3), countp1(3)
    integer :: ncnt
    call t_startf('prim_movie_output:pio')

    n0=tl%n0
    call TimeLevel_Qdp( tl, qsplit, n0_Q)

    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0)) then
          if ((output_start_time(ios) .le. tl%nstep) .and. &
               (output_end_time(ios) .ge. tl%nstep) .and. &
               MODULO(tl%nstep,output_frequency(ios)) .eq. 0) then
             output_varnames=>get_current_varnames(ios)
             start2d(1)=0
             start2d(2)=nf_get_frame(ncdf(ios))
             count2d(1)=0
             count2d(2)=1

             count(1:2)=0
             start(1:2)=0
             start(3)=nf_get_frame(ncdf(ios))
             count(3)=1

             countp1(1:2)=0
             startp1(1:2)=0
             startp1(3)=start(3)
             countp1(3)=1

             if(nf_selectedvar('ps', output_varnames)) then
                if (par%masterproc) print *,'writing ps...'
                st=1
                do ie=1,nelemd
                   vartmp(:,:,1) = (elem(ie)%state%ps_v(:,:,n0))
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxP,vartmp(:,:,1),var2d(st:en))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var2d,start2d,count2d,name='ps')
             endif

             if(nf_selectedvar('precl', output_varnames)) then
                if (par%masterproc) print *,'writing precl...'
                st=1
                do ie=1,nelemd
                   vartmp(:,:,1) = precl(:,:,ie)
                   !vartmp(:,:,1) = elem(ie)%state%Q(:,:,(nlev*2)/3,1)  ! hack for movies
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxP,vartmp(:,:,1),var2d(st:en))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var2d,start2d,count2d,name='precl')
             endif

             if(nf_selectedvar('hypervis', output_varnames)) then
                if (par%masterproc) print *,'writing hypervis...'
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
                if (par%masterproc) print *,'writing zeta...'
                ! velocities are on sphere for primitive equations
                call compute_zeta_C0(temp3d,elem,par,n0)

                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxp,nlev,temp3d(:,:,:,ie),var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='zeta')
             end if


             if(nf_selectedvar('T', output_varnames)) then
                if (par%masterproc) print *,'writing T...'
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call get_field(elem(ie),'temperature',vartmp,hvcoord,n0,n0_Q)
                   call UniquePoints(elem(ie)%idxP,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='T')
             end if


             if(nf_selectedvar('Th', output_varnames)) then
                pr0=1./(p0)
                st=1
                do ie=1,nelemd
                   call get_field(elem(ie),'pottemp',vartmp,hvcoord,n0,n0_Q)
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxP,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                end do
                call nf_put_var(ncdf(ios),var3d,start, count, name='Th')
             end if

            if(nf_selectedvar('rho', output_varnames)) then
                if (par%masterproc) print *,'writing rho...'
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call get_field(elem(ie),'rho',vartmp,hvcoord,n0,n0_Q)
                   call UniquePoints(elem(ie)%idxP,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='T')
             end if

             if(nf_selectedvar('u', output_varnames)) then
                if (par%masterproc) print *,'writing u...'
                st=1
                do ie=1,nelemd
                   do k=1,nlev
                      vartmp(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)
                   end do
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxp,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='u')
             end if

             if(nf_selectedvar('v', output_varnames)) then
                if (par%masterproc) print *,'writing v...'
                st=1
                do ie=1,nelemd
                   do k=1,nlev
                      vartmp(:,:,k) = elem(ie)%state%v(:,:,2,k,n0)
                   end do
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxp,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='v')
             end if

             if(nf_selectedvar('ke', output_varnames)) then
                st=1
                do ie=1,nelemd
                   do k=1,nlev
                      ke(:,:,k) = (elem(ie)%state%v(:,:,1,k,n0)**2 + &
                      elem(ie)%state%v(:,:,2,k,n0)**2 )/2
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
                   if (par%masterproc) print *,'writing ',vname
                   st=1
                   do ie=1,nelemd
                      en=st+elem(ie)%idxp%NumUniquePts-1
                      call UniquePoints(elem(ie)%idxP,nlev,elem(ie)%state%Q(:,:,:,qindex), var3d(st:en,:))
                      st=en+1
                   end do
                   call nf_put_var(ncdf(ios),var3d,start, count, name=vname)
                end if
             enddo

             if(nf_selectedvar('geo', output_varnames)) then
                if (par%masterproc) print *,'writing geo...'
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call get_field(elem(ie),'geo',vartmp,hvcoord,n0,n0_Q)
                   call UniquePoints(elem(ie)%idxP,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='geo')
             end if

             if(nf_selectedvar('geo_i', output_varnames)) then
                if (par%masterproc) print *,'writing geo_i...'
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call get_field_i(elem(ie),'geo_i',vartmp1,hvcoord,n0)
                   call UniquePoints(elem(ie)%idxP,nlevp,vartmp1,var3dp1(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3dp1,startp1,countp1,name='geo_i',iodescin=iodesc3dp1)
             end if

             if(nf_selectedvar('w', output_varnames)) then
                if (par%masterproc) print *,'writing w...'
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call get_field(elem(ie),'w',vartmp,hvcoord,n0,n0_Q)
                   call UniquePoints(elem(ie)%idxP,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start, count, name='w')
             end if

             if(nf_selectedvar('w_i', output_varnames)) then
                if (par%masterproc) print *,'writing w_i...'
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call get_field_i(elem(ie),'w_i',vartmp1,hvcoord,n0)
                   call UniquePoints(elem(ie)%idxP,nlevp,vartmp1,var3dp1(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3dp1,startp1, countp1, name='w_i',iodescin=iodesc3dp1)
             end if


             if(nf_selectedvar('mu_i', output_varnames)) then
                if (par%masterproc) print *,'writing mu_i...'
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call get_field_i(elem(ie),'mu_i',vartmp1,hvcoord,n0)
                   call UniquePoints(elem(ie)%idxP,nlevp,vartmp1,var3dp1(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3dp1,startp1, countp1, name='mu_i',iodescin=iodesc3dp1)
             end if

             if(nf_selectedvar('pnh', output_varnames)) then
                if (par%masterproc) print *,'writing pnh...'
                st=1
                do ie=1,nelemd
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call get_field(elem(ie),'pnh',vartmp,hvcoord,n0,n0_Q)
                   call UniquePoints(elem(ie)%idxP,nlev,vartmp,var3d(st:en,:))
                   st=en+1
                enddo
                call nf_put_var(ncdf(ios),var3d,start,count,name='pnh')
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


             if(nf_selectedvar('dp3d', output_varnames)) then
                st=1
                do ie=1,nelemd
                   do k=1,nlev
                      ke(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,n0) !reutilize ke
                   enddo
                   en=st+elem(ie)%idxp%NumUniquePts-1
                   call UniquePoints(elem(ie)%idxp,nlev,ke, var3d(st:en,:))
                   st=en+1
                end do
                call nf_put_var(ncdf(ios),var3d,start, count, name='dp3d')
             end if


             if (par%masterproc) print *,'writing time...'
             call nf_put_var(ncdf(ios),real(dayspersec*time_at(tl%nstep),kind=real_kind),&
                  start(3:3),count(3:3),name='time')
             call nf_advance_frame(ncdf(ios))
             call pio_syncfile(ncdf(ios)%fileid)
             if (par%masterproc) print *,'finished I/O sync'
          end if
       end if
    end do
    call t_stopf('prim_movie_output:pio')
 end subroutine prim_movie_output

#endif
!for PIOINTERP
#endif
!for WITHOUT_PIOLIB
end module prim_movie_mod
