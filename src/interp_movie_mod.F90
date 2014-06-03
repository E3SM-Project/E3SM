#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module interp_movie_mod
  use kinds, only : real_kind
  use dimensions_mod, only :  nlev, nelemd, np, ne, qsize, ntrac, nc
  use interpolate_mod, only : interpolate_t, setup_latlon_interp, interpdata_t, &
       get_interp_parameter, get_interp_lat, get_interp_lon, interpolate_scalar, interpolate_vector, &
       set_interp_parameter, interpol_phys_latlon
#if defined(_SPELT)
  use interpolate_mod, only : interpol_spelt_latlon
#endif
  use pio_io_mod, only : & 
       nf_output_init_begin,&
       nf_output_init_complete,  &
       nf_output_register_variables,&
       nf_put_var, &
       nf_close_all, &
       nf_output_register_dims, &
       nf_selectedvar, &
       nf_advance_frame, &
       nf_handle, &
       get_current_varnames, &
       nf_variable_attributes, &
       nfsizekind,             &
       nf_get_frame,           &
       PIO_double,              &
       nf_init_decomp, &
       get_varindex

  use control_mod, only : test_case, runtype, accumstart, &
       accumstop, accumfreq, restartfreq, &
       integration, columnpackage, kmass, nu
  use common_io_mod, only:  &
       output_start_time,   & 	
       output_end_time,     &
       output_frequency,    &
       output_dir,          &
       max_output_variables,&
       max_output_streams,  &
       varname_len,         &
       nf_addrequiredvar,   &
       num_io_procs,        &
       PIOFS
  use fvm_control_volume_mod, only : fvm_struct
  use spelt_mod, only : spelt_struct

  implicit none
#undef V_IS_LATLON
#if defined(_PRIM) || defined(_PRIMDG)
#define V_IS_LATLON
  integer, parameter :: varcnt = 43
  integer, parameter :: maxdims =  5
  character*(*), parameter :: varnames(varcnt)=(/'ps       ', &
                                                 'geos     ', &
                                                 'zeta     ', &
                                                 'dp3d     ', &
                                                 'div      ', &
                                                 'T        ', &
                                                 'Th       ', &
                                                 'u        ', &
                                                 'v        ', &
                                                 'ke       ', &
                                                 'Q        ', &
                                                 'Q2       ', &
                                                 'Q3       ', &
                                                 'Q4       ', &
                                                 'Q5       ', &
                                                 'psC      ', &
                                                 'C        ', &
                                                 'C2       ', &
                                                 'C3       ', &
                                                 'C4       ', &
                                                 'C5       ', &                                                 
                                                 'geo      ', &
                                                 'omega    ', &
                                                 'FU       ', &
                                                 'FV       ', &
                                                 'DIFFU    ', &
                                                 'DIFFV    ', &
                                                 'DIFFT    ', &
                                                 'hypervis ', &
                                                 'max_dx   ', &
                                                 'min_dx   ', &
                                                 'CONVU    ', &
                                                 'CONVV    ', &
                                                 'lat      ', &
                                                 'lon      ', &
                                                 'gw       ', &
                                                 'lev      ', &
                                                 'ilev     ', &
                                                 'hyam     ', &
                                                 'hybm     ', &
                                                 'hyai     ', &
                                                 'hybi     ', &
                                                 'time     '/)
  integer, parameter :: vartype(varcnt)=(/PIO_double,PIO_double,PIO_double,PIO_double,&
                                          PIO_double,&
                                          PIO_double,PIO_double,PIO_double,PIO_double,&
                                          PIO_double,PIO_double,PIO_double,PIO_double,&
                                          PIO_double,PIO_double,PIO_double,&
                                          PIO_double,PIO_double,&
                                          PIO_double,PIO_double,PIO_double,PIO_double,&
                                          PIO_double,PIO_double,PIO_double,PIO_double,&
                                          PIO_double,&
                                          PIO_double,&
                                          PIO_double,PIO_double,&
                                          PIO_double,&
                                          PIO_double,PIO_double,PIO_double,PIO_double,&
                                          PIO_double,PIO_double,PIO_double,&
                                          PIO_double,PIO_double,&
                                          PIO_double,PIO_double,&
                                          PIO_double/)
  logical, parameter :: varrequired(varcnt)=(/.false.,.false.,.false.,.false.,.false.,&
                                              .false.,&   
                                              .false.,.false.,.false.,.false.,.false.,&
                                              .false.,.false.,.false.,.false.,.false.,&
                                              .false.,.false.,.false.,.false.,.false.,&
                                              .false.,.false.,.false.,.false.,.false.,&
                                              .false.,.false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,.false.,.true. ,.true. ,&
                                              .true.,.true. ,.true. ,&   ! gw,lev,ilev
                                              .true. ,.true. ,&   ! hy arrays
                                              .true. ,.true. ,&   ! hy arrays
                                              .true./)

  integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,2,5,0,0,  &
       1,2,0,0,0,  &   ! geos
       1,2,3,5,0,  &   ! zeta
       1,2,3,5,0,  &   ! dp3d
       1,2,3,5,0,  &   ! div
       1,2,3,5,0,  &   ! T
       1,2,3,5,0,  &   ! Th
       1,2,3,5,0,  &   ! u
       1,2,3,5,0,  &   ! v
       1,2,3,5,0,  &   ! ke
       1,2,3,5,0,  &   ! Q
       1,2,3,5,0,  &   ! Q2
       1,2,3,5,0,  &   ! Q3
       1,2,3,5,0,  &   ! Q4
       1,2,3,5,0,  &   ! Q5
       1,2,5,0,0,  &   ! psC
       1,2,3,5,0,  &   ! C
       1,2,3,5,0,  &   ! C2
       1,2,3,5,0,  &   ! C3
       1,2,3,5,0,  &   ! C4
       1,2,3,5,0,  &   ! C5
       1,2,3,5,0,  &   ! geo
       1,2,3,5,0,  &   ! omega
       1,2,3,5,0,  &   ! FU
       1,2,3,5,0,  &   ! FV
       1,2,3,5,0,  &   ! DIFFU
       1,2,3,5,0,  &   ! DIFFV
       1,2,3,5,0,  &   ! DIFFT
       1,2,0,0,0,  &   ! hypervis
       1,2,0,0,0,  &   ! max dx
       1,2,0,0,0,  &   ! min dx
       1,2,3,5,0,  &   ! CONVU
       1,2,3,5,0,  &   ! CONVV
       2,0,0,0,0,  &   ! lat
       1,0,0,0,0,  &   ! lon
       2,0,0,0,0,  &   ! gw
       3,0,0,0,0,  &   ! lev
       4,0,0,0,0,  &   ! ilev
       3,0,0,0,0,  &   ! hyam
       3,0,0,0,0,  &   ! hybm
       4,0,0,0,0,  &   ! hyai
       4,0,0,0,0,  &   ! hybi
       5,0,0,0,0 /),&
       shape=(/maxdims,varcnt/))

  character*(*),parameter::dimnames(maxdims)=(/'lon ','lat ','lev ','ilev','time'/)  
#else
  integer, parameter :: varcnt = 19
  integer, parameter :: maxdims=4
  character*(*),parameter::dimnames(maxdims)=(/'lon ','lat ','lev ','time'/)  
  integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,2,4,0,  &
                                                               1,2,3,4,  &
                                                               1,2,3,4,  &
                                                               1,2,3,4,  &
                                                               1,2,3,4,  &
                                                               1,0,0,0,  &
                                                               2,0,0,0,  &
                                                               2,0,0,0,  &
                                                               4,0,0,0, &
                                                               1,2,0,0, &
                                                               1,2,0,0, &
                                                               1,2,0,0, &
                                                               1,2,4,0, &
                                                               1,2,3,4,  &
                                                               1,2,3,4,  &
                                                               1,2,3,4,  &
                                                               1,2,3,4,  &
                                                               1,2,3,4,  &
                                                               1,2,3,4/),&
                                                               shape=(/maxdims,varcnt/))
  character*(*),parameter::varnames(varcnt)=(/'ps      ','geop    ','u       ', &
                                              'v       ','zeta    ','lon     ', &
                                              'lat     ','gw      ','time    ', &
                                              'hypervis','max_dx  ','min_dx  ', &
                                              'psC     ','C       ','C1      ','C2      ', &
                                              'C3      ','C4      ',            &
                                              'div     '/)
  integer, parameter :: vartype(varcnt)=(/PIO_double,PIO_double,PIO_double,PIO_double, &
                                          PIO_double,PIO_double,PIO_double,PIO_double, &
                                          PIO_double,PIO_double,PIO_double,PIO_double, &
                                          PIO_double, PIO_double,PIO_double,PIO_double, &
                                          PIO_double, PIO_double, PIO_double/)
  logical, parameter :: varrequired(varcnt)=(/.false.,.false.,.false.,.false.,&
                                              .false.,.true.,.true.,.true.,.true.,&
                                              .false.,.false.,.false.,.false., &
                                              .false.,.false.,.false.,.false.,.false.,.false./)

#endif
  type(interpolate_t) :: interp
  type(interpdata_t), allocatable :: interpdata(:)

  integer(kind=nfsizekind) :: start2d(3), count2d(3), start3d(4), count3d(4)
  type(nf_handle), save :: ncdf(max_output_streams)


contains
  subroutine interp_movie_init(elem,hybrid,nets,nete,hvcoord,tl)
    use time_mod, only : timelevel_t
    use hybrid_mod, only : hybrid_t
    use element_mod, only: element_t
    use pio, only : pio_setdebuglevel, PIO_Put_att, pio_put_var, pio_global ! _EXTERNAL
    use parallel_mod, only : parallel_t, haltmp, syncmp
    use interpolate_mod, only : get_interp_lat, get_interp_lon, get_interp_gweight
#if defined(_PRIM) || defined(_PRIMDG)
	use hybvcoord_mod, only : hvcoord_t
    use aquaplanet_io_mod, only : aq_movie_init
    use physics_io_mod, only : physics_movie_init
#endif

    type (TimeLevel_t), intent(in)         :: tl     ! time level struct
    type(element_t) :: elem(:)
    
    type(hybrid_t), target  :: hybrid
#if defined(_PRIM) || defined(_PRIMDG)
    type(hvcoord_t), intent(in), optional :: hvcoord
#else
    ! ignored
    integer, optional :: hvcoord
#endif
    integer, intent(in) :: nets, nete
    integer :: dimsize(maxdims)   
    integer, pointer :: ldof2d(:),ldof3d(:), iodof2d(:), iodof3d(:)
    integer, pointer :: latdof(:), londof(:), iodoflon(:), iodoflat(:)

    integer :: icnt, i, j, k, lcount, iorank, nlat, nlon, tdof(1), tiodof(1), ios, ie
    type(parallel_t), pointer :: par
    integer(kind=nfsizekind) :: start1d(1), count1d(1)
    real(kind=real_kind), allocatable :: lat(:), lon(:), gw(:)
    real(kind=real_kind), allocatable :: lev(:),ilev(:)
    integer :: varid, vindex
    integer :: ierr
    character(len=9)     :: charnum
    character(len=90)    :: hname

    allocate(interpdata(nelemd))

    call setup_latlon_interp(elem,interpdata, hybrid%par)

    lcount = sum(interpdata(nets:nete)%n_interp)
    par => hybrid%par

    nlat = get_interp_parameter('nlat')
    nlon = get_interp_parameter('nlon')


    if (runtype==0) then
       hname = test_case
    else
       ! restart runs (1) or (2)
       write(charnum,'(i9.9)') tl%nstep
       hname = trim(ADJUSTL(test_case)) // "-" // charnum //"-"
    endif

    call PIO_setDebugLevel(0)
    call nf_output_init_begin(ncdf,par%masterproc,par%nprocs,par%rank, &
         par%comm,hname,runtype)
#if defined(_PRIM) || defined(_PRIMDG)
    dimsize=(/nlon,nlat,nlev,nlev+1,0/)
#else
    dimsize=(/nlon,nlat,nlev,0/)
#endif
    call nf_output_register_dims(ncdf, maxdims, dimnames, dimsize)

    iorank=piofs%io_rank

    ! Create the DOF arrays
    allocate(ldof2d(lcount))
    allocate(ldof3d(lcount*nlev))
    icnt=0
    do ie=nets,nete
       do i=1,interpdata(ie)%n_interp
          icnt=icnt+1
          ldof2d(icnt)=interpdata(ie)%ilon(i)+(interpdata(ie)%ilat(i)-1)*nlon
       end do
    end do
    icnt=0
    do k=1,nlev
       do ie=nets,nete
          do i=1,interpdata(ie)%n_interp
             icnt=icnt+1
             ldof3d(icnt)=interpdata(ie)%ilon(i)+(interpdata(ie)%ilat(i)-1)*nlon+(k-1)*nlat*nlon
          end do
       end do
    end do
    call getiodof(2, (/nlon,nlat/), iorank, iodof2d, start2d(1:2), count2d(1:2))

    call nf_init_decomp(ncdf, (/1,2/), ldof2d, iodof2d,start2d(1:2),count2d(1:2))

    call getiodof(3, (/nlon,nlat,nlev/), iorank, iodof3d, start3d(1:3), count3d(1:3))

    call nf_init_decomp(ncdf, (/1,2,3/), ldof3d, iodof3d,start3d(1:3),count3d(1:3))

    deallocate(iodof2d, iodof3d, ldof2d,ldof3d)

    call nf_output_register_variables(ncdf,varcnt,varnames,vardims,vartype,varrequired)
    do ios = 1, max_output_streams
       if((output_frequency(ios) .gt. 0) ) then
          i=PIO_Put_att(ncdf(ios)%fileid,pio_global,'np',np)
          i=PIO_Put_att(ncdf(ios)%fileid,pio_global,'ne',ne)
       endif
    enddo
    !call nf_global_attribute(ncdf, 'np', np)
    !call nf_global_attribute(ncdf, 'ne', ne)

    call nf_variable_attributes(ncdf, 'ps', 'surface pressure','Pa')
    call nf_variable_attributes(ncdf, 'u', 'longitudinal wind component','meters/second')
    call nf_variable_attributes(ncdf, 'v', 'latitudinal wind component','meters/second')
    call nf_variable_attributes(ncdf, 'zeta', 'Relative vorticity','1/s')
#if defined(_PRIM) || defined(_PRIMDG)
    call nf_variable_attributes(ncdf, 'geo', 'Geopotential','m^2/s^2')
    call nf_variable_attributes(ncdf, 'geos', 'Surface geopotential','m^2/s^2')
    call nf_variable_attributes(ncdf, 'T', 'Temperature','degrees kelvin')
    call nf_variable_attributes(ncdf, 'dp3d', 'delta p','Pa')
    call nf_variable_attributes(ncdf, 'Q', 'concentration','kg/kg')
    call nf_variable_attributes(ncdf, 'Q2', 'concentration','kg/kg')
    call nf_variable_attributes(ncdf, 'Q3', 'concentration','kg/kg')
    call nf_variable_attributes(ncdf, 'Q4', 'concentration','kg/kg')
    call nf_variable_attributes(ncdf, 'Q5', 'concentration','kg/kg')
    call nf_variable_attributes(ncdf, 'lev' ,'hybrid level at midpoints' ,'level','positive','down') !,'formula_terms','a: hyam b: hybm p0: P0 ps: PS')
    call nf_variable_attributes(ncdf, 'ilev','hybrid level at interfaces','level','positive','down') !,'formula_terms','a: hyai b: hybi p0: P0 ps: PS')
    call nf_variable_attributes(ncdf, 'hyam','hybrid A coefficiet at layer midpoints' ,'dimensionless') 
    call nf_variable_attributes(ncdf, 'hybm','hybrid B coefficiet at layer midpoints' ,'dimensionless') 
    call nf_variable_attributes(ncdf, 'hyai','hybrid A coefficiet at layer interfaces' ,'dimensionless') 
    call nf_variable_attributes(ncdf, 'hybi','hybrid B coefficiet at layer interfaces' ,'dimensionless') 
#endif
    call nf_variable_attributes(ncdf, 'gw', 'gauss weights','dimensionless')
    call nf_variable_attributes(ncdf, 'lat', 'column latitude','degrees_north')
    call nf_variable_attributes(ncdf, 'lon', 'column longitude','degrees_east')
    call nf_variable_attributes(ncdf, 'time', 'Model elapsed time','days')
    call nf_variable_attributes(ncdf, 'psC', 'surface pressure','Pa')
    call nf_variable_attributes(ncdf, 'C', 'concentration','kg/kg')
    call nf_variable_attributes(ncdf, 'C2', 'concentration','kg/kg')
    call nf_variable_attributes(ncdf, 'C3', 'concentration','kg/kg')
    call nf_variable_attributes(ncdf, 'C4', 'concentration','kg/kg')
    call nf_variable_attributes(ncdf, 'C5', 'concentration','kg/kg')

#if defined(_PRIM) || defined(_PRIMDG)
    if(test_case.eq.'aquaplanet') then
       call aq_movie_init(ncdf)
    end if
    if(columnpackage.ne.'none') then
       call physics_movie_init(ncdf)
    end if
#endif
    call nf_output_init_complete(ncdf)
    allocate(lon(nlon), lat(nlat), gw(nlat))
    allocate(lev(nlev), ilev(nlev+1))
    lon = get_interp_lon()
    lat = get_interp_lat()


    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0) ) then
          if(iorank==0) print *,"writing coordinates to ios=",ios

        vindex = get_varindex('lon',ncdf(ios)%varlist)
        varid = ncdf(ios)%varlist(vindex)%vardesc%varid
	ierr = pio_put_var(ncdf(ios)%FileID,varid, lon)
        vindex = get_varindex('lat',ncdf(ios)%varlist)
        varid = ncdf(ios)%varlist(vindex)%vardesc%varid
	ierr = pio_put_var(ncdf(ios)%FileID,varid, lat)

        ! output gauss weights
        gw = get_interp_gweight()
        vindex = get_varindex('gw',ncdf(ios)%varlist)
        varid = ncdf(ios)%varlist(vindex)%vardesc%varid
        ierr = pio_put_var(ncdf(ios)%FileID,varid, gw)


#if defined(_PRIM) || defined(_PRIMDG)
          if (present(hvcoord)) then
             vindex = get_varindex('lev',ncdf(ios)%varlist)
             varid = ncdf(ios)%varlist(vindex)%vardesc%varid
             ierr = pio_put_var(ncdf(ios)%FileID,varid, hvcoord%etam)

             vindex = get_varindex('ilev',ncdf(ios)%varlist)
             varid = ncdf(ios)%varlist(vindex)%vardesc%varid
             ierr = pio_put_var(ncdf(ios)%FileID,varid, hvcoord%etai)

             vindex = get_varindex('hyam',ncdf(ios)%varlist)
             varid = ncdf(ios)%varlist(vindex)%vardesc%varid
             ierr = pio_put_var(ncdf(ios)%FileID,varid, hvcoord%hyam)

             vindex = get_varindex('hybm',ncdf(ios)%varlist)
             varid = ncdf(ios)%varlist(vindex)%vardesc%varid
             ierr = pio_put_var(ncdf(ios)%FileID,varid, hvcoord%hybm)

             vindex = get_varindex('hyai',ncdf(ios)%varlist)
             varid = ncdf(ios)%varlist(vindex)%vardesc%varid
             ierr = pio_put_var(ncdf(ios)%FileID,varid, hvcoord%hyai)

             vindex = get_varindex('hybi',ncdf(ios)%varlist)
             varid = ncdf(ios)%varlist(vindex)%vardesc%varid
             ierr = pio_put_var(ncdf(ios)%FileID,varid, hvcoord%hybi)
          end if
#endif
       endif
    end do
    deallocate(lat,lon,gw)
    deallocate(lev,ilev)
    call syncmp(par)
  end subroutine interp_movie_init



  subroutine interp_movie_finish
    call nf_close_all(ncdf)
  end subroutine interp_movie_finish



  subroutine interp_movie_output(elem, tl, hybrid, phimean, nets,nete, fvm, hvcoord)

    use kinds, only : int_kind, real_kind
    use element_mod, only : element_t
    use time_mod, only : Timelevel_t, tstep, ndays, time_at, secpday, nendstep,nmax
    use parallel_mod, only : parallel_t, abortmp
#if defined(_PRIM) 
    use hybvcoord_mod, only :  hvcoord_t 
    use aquaplanet_io_mod, only : aq_movie_output
    use physics_io_mod, only : physics_movie_output
#elif defined _PRIMDG
    use hybvcoord_mod, only :  hvcoord_t
#endif
    use physical_constants, only : omega, g, rrearth, dd_pi, kappa, p0
    use derivative_mod, only : derivinit, derivative_t, vorticity, laplace_sphere_wk
    use hybrid_mod, only : hybrid_t
    use pio, only : pio_setdebuglevel, pio_syncfile ! _EXTERNAL

    use viscosity_mod, only : compute_zeta_C0, make_c0, compute_zeta_c0_contra,&
                              compute_div_c0,compute_div_c0_contra
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    ! ---------------------    
    type (element_t),target    :: elem(:)
    
#if defined(_SPELT)
    type (spelt_struct), optional   :: fvm(:)
#else
    type (fvm_struct), optional   :: fvm(:)    
#endif
    
    type (TimeLevel_t)  :: tl
    type (parallel_t)     :: par
#if defined(_PRIM)
    type (hvcoord_t)    :: hvcoord
#elif defined(_PRIMDG)
    type (hvcoord_t)    :: hvcoord
#else
    integer,optional    :: hvcoord
#endif
    real (kind=real_kind), intent(in) :: phimean

    type (hybrid_t)      , intent(in) :: hybrid
    integer              :: nets,nete

    character(len=varname_len), pointer :: output_varnames(:)
    integer :: ie,ios, i, j, k
    real (kind=real_kind) :: pfull, pr0
    real(kind=real_kind),parameter :: dayspersec=1d0/(3600.*24.)
    real(kind=real_kind), allocatable :: datall(:,:), var3d(:,:,:,:)
    real(kind=real_kind), allocatable :: varvtmp(:,:,:,:), ulatlon(:,:,:,:,:)
    
    integer :: st, en

    integer :: ierr

    integer :: ncnt,n0,n0q,itype,qindex,cindex
    character(len=2) :: vname

    real (kind=real_kind) :: vco(np,np,2),ke(np,np,nlev)
    real (kind=real_kind) :: v1,v2

    type (derivative_t)  :: deriv


    call t_startf('interp_movie_output')
    if (hybrid%NThreads /= 1) &
         call abortmp('Error: interp_movie_output can only be called outside threaded region')

    n0 = tl%n0
    call derivinit(deriv)

!    if (0==piofs%io_rank) write(*,'(a,i4,a,i1)') &
!         "lat/lon interp movie output: ios=",ios," interpolation type=",&
!         get_interp_parameter("itype")
     
    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0)) then
          if ((output_start_time(ios) .le. tl%nstep) .and. &
               (output_end_time(ios) .ge. tl%nstep) .and. &
               MODULO(tl%nstep,output_frequency(ios)) .eq. 0) then

             ncnt = sum(interpdata(nets:nete)%n_interp)  ! ncnt not defined if output disabled
             output_varnames=>get_current_varnames(ios)

             start2d(3)=nf_get_frame(ncdf(ios))
             count2d(3)=1
             start3d(4)=nf_get_frame(ncdf(ios))
             count3d(4)=1


             if(nf_selectedvar('ps', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing ps...'
                st=1
                allocate (datall(ncnt,1))

                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
#ifdef _PRIM
                   call interpolate_scalar(interpdata(ie),elem(ie)%state%ps_v(:,:,n0), &
                        np, datall(st:en,1))
#elif defined _PRIMDG
                   call interpolate_scalar(interpdata(ie),elem(ie)%state%pr3d(:,:,nlev+1), &
                        np, datall(st:en,1))
#else
                   call interpolate_scalar(interpdata(ie),elem(ie)%state%ps, &
                        np, datall(st:en,1))
#endif
                   st=st+interpdata(ie)%n_interp
                enddo
	        
#ifdef _PRIM
                if (p0 < 2000)  then  ! convert to Pa, if using mb
                   datall(:,1) = 100*(datall(:,1)) 
                endif
#endif
                call nf_put_var(ncdf(ios),datall(:,1),start2d,count2d,name='ps')
                deallocate(datall)
             endif


             if(nf_selectedvar('zeta', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing zeta...'
                allocate(datall(ncnt,nlev))
                allocate(var3d(np,np,nlev,nets:nete))
#ifdef V_IS_LATLON
                ! velocities are on sphere for primitive equations
                call compute_zeta_C0(var3d,elem,hybrid,nets,nete,n0)
#else
                call compute_zeta_C0_contra(var3d,elem,hybrid,nets,nete,n0)
#endif
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie), var3d(:,:,:,ie), &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='zeta')
                deallocate(datall, var3d)
             end if

             if(nf_selectedvar('div', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing div...'
                allocate(datall(ncnt,nlev))
                allocate(var3d(np,np,nlev,nets:nete))
#ifdef V_IS_LATLON
                ! velocities are on sphere for primitive equations
                call compute_div_C0(var3d,elem,hybrid,nets,nete,n0)
#else
                call compute_div_C0_contra(var3d,elem,hybrid,nets,nete,n0)
#endif
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie), var3d(:,:,:,ie), &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='div')
                deallocate(datall, var3d)
             end if

             if(nf_selectedvar('u', output_varnames) .or. &
                  nf_selectedvar('v', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing u,v...'
                allocate(var3d(ncnt,nlev,2,1))
#ifdef V_IS_LATLON
                st=1                
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_vector(interpdata(ie), elem(ie), &
                                        elem(ie)%state%v(:,:,:,:,n0), np, nlev, var3d(st:en,:,:,1), 0)
                   st=st+interpdata(ie)%n_interp
                enddo
#else
                allocate(varvtmp(np,np,2,nlev))
                ! I = D^-1 D 
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   do k=1,nlev
                       varvtmp(:,:,1,k) = elem(ie)%D(1,1,:,:)*elem(ie)%state%v(:,:,1,k,n0) + &
                                          elem(ie)%D(1,2,:,:)*elem(ie)%state%v(:,:,2,k,n0)
                       varvtmp(:,:,2,k) = elem(ie)%D(2,1,:,:)*elem(ie)%state%v(:,:,1,k,n0) + &
                                          elem(ie)%D(2,2,:,:)*elem(ie)%state%v(:,:,2,k,n0)
                   end do
                   call interpolate_vector(interpdata(ie), elem(ie), &
                                        varvtmp, np, nlev, var3d(st:en,:,:,1), 0)
                   st=st+interpdata(ie)%n_interp
                enddo
                deallocate(varvtmp)

#endif
                if(nf_selectedvar('u', output_varnames)) then
                   call nf_put_var(ncdf(ios),var3d(:,:,1,1),start3d, count3d, name='u')
                end if

                if(nf_selectedvar('v', output_varnames)) then
                   call nf_put_var(ncdf(ios),var3d(:,:,2,1),start3d, count3d, name='v')
                end if
                deallocate(var3d)
             end if

#if defined(_PRIMDG) 
             if(nf_selectedvar('T', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing real temperature Y (not potential)...'
                allocate(datall(ncnt,nlev))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie), elem(ie)%state%T(:,:,:,n0), &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='T')
                deallocate(datall)
             end if
             if(nf_selectedvar('Q', output_varnames) .and. qsize>=1) then
                if (hybrid%par%masterproc) print *,'writing DG moist Q...'
                allocate(datall(ncnt,nlev))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   ! Not sure if this interpolation is the right one, HOMME-SE uses something different
                   ! but for now we get what we can to the output file.
                   call interpolate_scalar(interpdata(ie), elem(ie)%state%Q(:,:,:,n0), &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='Q')
                 deallocate(datall)
             end if

#endif

#if defined(_FVM) 
           if(nf_selectedvar('psC', output_varnames)) then
              if (hybrid%par%masterproc) print *,'writing psC...'
              st=1
              allocate(datall(ncnt,1))
              do ie=nets,nete
                 en=st+interpdata(ie)%n_interp-1
                 call interpol_phys_latlon(interpdata(ie),fvm(ie)%psc, &
                                    fvm(ie),elem(ie)%corners,elem(ie)%desc,datall(st:en,1))
                 st=st+interpdata(ie)%n_interp
              enddo
        
#ifdef _PRIM
              if (p0 < 2000)  then  ! convert to Pa, if using mb
                 datall(:,1) = 100*(datall(:,1)) 
              endif
#endif
              call nf_put_var(ncdf(ios),datall(:,1),start2d,count2d,name='psC')
              deallocate(datall)
           endif          
           
            do cindex=1,min(ntrac,5)  ! allow a maximum output of 5 tracers
               write(vname,'(a1,i1)') 'C',cindex
               if (cindex==1) vname='C'

               if(nf_selectedvar(vname, output_varnames)) then
                  if (hybrid%par%masterproc) print *,'writing ',vname
                  allocate(datall(ncnt,nlev))
                  st=1
                  do ie=nets,nete
                     en=st+interpdata(ie)%n_interp-1
                     do k=1,nlev                       
                       call interpol_phys_latlon(interpdata(ie),fvm(ie)%c(:,:,k,cindex,n0), &
                                          fvm(ie),elem(ie)%corners,elem(ie)%desc,datall(st:en,k))
                     end do
                     st=st+interpdata(ie)%n_interp
                  enddo                  
                  call nf_put_var(ncdf(ios),datall,start3d, count3d, name=vname)                  
                  deallocate(datall)                  
               end if
            enddo
#endif

#if defined(_SPELT) 
           if(nf_selectedvar('psC', output_varnames)) then
              if (hybrid%par%masterproc) print *,'writing for SPELT: psC...'
              st=1
              allocate(datall(ncnt,1))
              do ie=nets,nete
                 en=st+interpdata(ie)%n_interp-1
                 call interpol_spelt_latlon(interpdata(ie),fvm(ie)%psc, &
                                    fvm(ie),elem(ie)%corners,datall(st:en,1))
                 st=st+interpdata(ie)%n_interp
              enddo
        
#ifdef _PRIM
              if (p0 < 2000)  then  ! convert to Pa, if using mb
                 datall(:,1) = 100*(datall(:,1)) 
              endif
#endif
              call nf_put_var(ncdf(ios),datall(:,1),start2d,count2d,name='psC')
              deallocate(datall)
           endif


            do cindex=1,min(ntrac,5)  ! allow a maximum output of 5 tracers
               write(vname,'(a1,i1)') 'C for SPELT: ',cindex
               if (cindex==1) vname='C'

               if(nf_selectedvar(vname, output_varnames)) then
                  if (hybrid%par%masterproc) print *,'writing for SPELT: ',vname
                  allocate(datall(ncnt,nlev))
                  st=1


                  do ie=nets,nete
                     en=st+interpdata(ie)%n_interp-1
                     do k=1,nlev                       
                       call interpol_spelt_latlon(interpdata(ie),fvm(ie)%c(:,:,k,cindex,n0), &
                                          fvm(ie),elem(ie)%corners,datall(st:en,k))                                          
                     end do
                     st=st+interpdata(ie)%n_interp
                  enddo                  
                  call nf_put_var(ncdf(ios),datall,start3d, count3d, name=vname)                  
                  deallocate(datall)                  
               end if
            enddo
#endif



             if(nf_selectedvar('geop', output_varnames)) then
                allocate(datall(ncnt,nlev),var3d(np,np,nlev,1))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   do k=1,nlev
#ifdef _PRIM
                      var3d(:,:,k,1) = 0  ! need to compute PHI from hydrostatic relation
#elif defined _SWDG
                      var3d(:,:,k,1) = elem(ie)%state%ht(:,:,k)
#elif defined _PRIMDG
                      var3d(:,:,k,1) = (elem(ie)%state%p(:,:,k,n0) + elem(ie)%state%phis(:,:) + phimean)/g 
#else
                      if(test_case.eq.'vortex') then
                         var3d(:,:,k,1) = elem(ie)%state%p(:,:,k,n0)
                      elseif(test_case.eq.'swirl') then
                         var3d(:,:,k,1) = elem(ie)%state%p(:,:,k,n0)
                      else
                         var3d(:,:,k,1) = (elem(ie)%state%p(:,:,k,n0) + elem(ie)%state%ps(:,:) + phimean)/g
                      end if
                      if (kmass.ne.-1) then
                         ! p(:,:,kmass) = is the density, 
                         ! other levels are tracers.  Output concentration:
                         if(k.ne.kmass) &
                              var3d(:,:,k,1)=var3d(:,:,k,1)/elem(ie)%state%p(:,:,kmass,n0)
                      endif

#endif
                      call interpolate_scalar(interpdata(ie), var3d(:,:,k,1), &
                           np, datall(st:en,k))
                   end do
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='geop')
                deallocate(datall,var3d)
             end if



#if 0
             ! DEBUG code to output laplace_sphere_wk of surface pressure:
             if(nf_selectedvar('hypervis', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing hypervis...'
                allocate(datall(ncnt,nlev), var3d(np,np,nlev,nets:nete))

                do ie=nets,nete
                   do k=1,nlev
                      var3d(:,:,k,ie) = 0
#ifdef _PRIM
                      var3d(:,:,k,ie) = elem(ie)%state%ps_v(:,:,n0)
#elif defined _SWDG
                      var3d(:,:,k,ie) = 0  ! set this to surface pressure
#elif defined _PRIMDG
                      var3d(:,:,k,ie) = 0  ! set this to surface pressure
#else
                      var3d(:,:,k,ie) = elem(ie)%state%p(:,:,k,n0)
#endif
                      var3d(:,:,k,ie)=laplace_sphere_wk(var3d(:,:,k,ie),&
                           deriv,elem(ie),var_coef=.true.)
                      ! laplace_sphere_wk returns weak form with mass matrix
                      ! already applied.  remove mass matrix, since make_C0
                      ! routine below will also multiply by mass matrix before DSS
                      var3d(:,:,k,ie)=var3d(:,:,k,ie)/elem(ie)%spheremp(:,:)

                      ! viscosity coefficient
                      ! (normally we apply this operator twice, but here only 
                      ! once, so use sqrt(nu)
                      var3d(:,:,k,ie)=var3d(:,:,k,ie)*sqrt(nu)
                   enddo
                enddo
                call make_C0(var3d,elem,hybrid,nets,nete)
                print *,'min/max hypervis: ',minval(var3d),maxval(var3d)
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   do k=1,nlev
                      call interpolate_scalar(interpdata(ie), var3d(:,:,k,ie), &
                           np, datall(st:en,k))
                   end do
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='hypervis')
                deallocate(datall,var3d)
             end if
#else
             ! End hyperviscosity
             ! MNL: output hyperviscosity
             if(nf_selectedvar('hypervis', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing hypervis...'
                allocate(datall(ncnt,nlev),var3d(np,np,nlev,1))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   do k=1,nlev
                       var3d(:,:,k,1) = nu*elem(ie)%variable_hyperviscosity(:,:)**2

                      call interpolate_scalar(interpdata(ie), var3d(:,:,k,1), &
                           np, datall(st:en,k))
                   end do
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='hypervis')
                deallocate(datall,var3d)
             end if
             ! End hyperviscosity
#endif

             ! MNL: output hypervis length scale 
             if(nf_selectedvar('max_dx', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing max_dx...'
                allocate(datall(ncnt,nlev),var3d(np,np,nlev,1))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   do k=1,nlev
                       var3d(:,:,k,1) = elem(ie)%dx_long

                      call interpolate_scalar(interpdata(ie), var3d(:,:,k,1), &
                           np, datall(st:en,k))
                   end do
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='max_dx')
                deallocate(datall,var3d)
             end if

             ! MNL: output CFL length scale 
             if(nf_selectedvar('min_dx', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing min_dx...'
                allocate(datall(ncnt,nlev),var3d(np,np,nlev,1))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   do k=1,nlev
                       var3d(:,:,k,1) = elem(ie)%dx_short 

                      call interpolate_scalar(interpdata(ie), var3d(:,:,k,1), &
                           np, datall(st:en,k))
                   end do
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='min_dx')
                deallocate(datall,var3d)
             end if
             ! End dx outputs


#if defined(_PRIM) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            large block of _PRIM only I/O
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             if(nf_selectedvar('geos', output_varnames)) then
                if (nf_get_frame(ncdf(ios))==1) then
                st=1
                allocate (datall(ncnt,1))
                
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie),elem(ie)%state%phis(:,:), &
                        np, datall(st:en,1))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall(:,1),start2d,count2d,name='geos')
                deallocate(datall)
                endif
             endif

             if(nf_selectedvar('T', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing T...'
                allocate(datall(ncnt,nlev))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie), elem(ie)%state%T(:,:,:,n0), &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='T')
                deallocate(datall)
             end if

             if(nf_selectedvar('dp3d', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing dp3d...'
                allocate(datall(ncnt,nlev))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie), elem(ie)%state%dp3d(:,:,:,n0), &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='dp3d')
                deallocate(datall)
             end if

             if(nf_selectedvar('ke', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing ke...'
                allocate(datall(ncnt,nlev))
                st=1
                do ie=nets,nete
                   do k=1,nlev 
                      ke(:,:,k) = (elem(ie)%state%v(:,:,1,k,n0)**2 + &
        	       elem(ie)%state%v(:,:,2,k,n0)**2 )/2
                   enddo
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie), ke, &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='ke')
                deallocate(datall)
             end if

             if(nf_selectedvar('Th', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing Th...'
                pr0=1./(p0)
                st=1
                allocate(datall(ncnt,nlev),var3d(np,np,nlev,1))
                do ie=nets,nete
                   do k=1,nlev
                      do j=1,np
                         do i=1,np
                            pfull = hvcoord%hyam(k)*hvcoord%ps0  &
                                 + hvcoord%hybm(k)*exp(elem(ie)%state%lnps(i,j,n0))
                            var3d(i,j,k,1)=elem(ie)%state%T(i,j,k,n0)* &
                                 (pfull*pr0)**(-kappa)
                         end do
                      end do
                   end do
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie), var3d(:,:,:,1), &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                end do
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='Th')
                deallocate(datall,var3d)
             end if


             do qindex=1,min(qsize,5)
                write(vname,'(a1,i1)') 'Q',qindex
                if (qindex==1) vname='Q'

                if(nf_selectedvar(vname, output_varnames)) then
                   if (hybrid%par%masterproc) print *,'writing ',vname
                   ! switch to bilinear interpolation for tracers
                   itype=get_interp_parameter("itype")
                   call set_interp_parameter("itype",1)
                   allocate(datall(ncnt,nlev))
                   st=1
                   do ie=nets,nete
                      en=st+interpdata(ie)%n_interp-1
                      call interpolate_scalar(interpdata(ie), elem(ie)%state%Q(:,:,:,qindex), &
                           np, nlev, datall(st:en,:))
                      st=st+interpdata(ie)%n_interp
                   enddo
                   call nf_put_var(ncdf(ios),datall,start3d, count3d, name=vname)
                   deallocate(datall)
                   call set_interp_parameter("itype",itype)
                end if
             enddo


             if(nf_selectedvar('geo', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing geo...'
                allocate(datall(ncnt,nlev))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie), elem(ie)%derived%phi, &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='geo')
                deallocate(datall)
             end if

             if(nf_selectedvar('omega', output_varnames)) then
                if (hybrid%par%masterproc) print *,'writing omega...'
                allocate(datall(ncnt,nlev), var3d(np,np,nlev,nets:nete))
                do ie=nets,nete
                   do k=1,nlev
                      do j=1,np
                         do i=1,np
                            pfull = hvcoord%hyam(k)*hvcoord%ps0  &
                                 + hvcoord%hybm(k)*exp(elem(ie)%state%lnps(i,j,n0))
                            var3d(i,j,k,ie)=elem(ie)%derived%omega_p(i,j,k)*pfull
                         end do
                      end do
                   end do
                end do
                call make_C0(var3d,elem,hybrid,nets,nete)
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie), var3d(:,:,:,ie), &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='omega')
                deallocate(datall,var3d)
             end if

             ! note: this is kind of a hack: forcing is computed during the
             ! timestep, from timelevel nm1 and stored in FM(nm1). after
             ! the timestep is over, nm1 data will be in FM(np1)
             if(nf_selectedvar('FU', output_varnames) .or. &
                  nf_selectedvar('FV', output_varnames)) then
                allocate(var3d(ncnt,2,nlev,1))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_vector(interpdata(ie), elem(ie),  &
                        elem(ie)%derived%FM(:,:,:,:,tl%np1), np, nlev, var3d(st:en,:,:,1), 0)
                   st=st+interpdata(ie)%n_interp
                enddo

                if(nf_selectedvar('FU', output_varnames)) then
                   call nf_put_var(ncdf(ios),var3d(:,1,:,1),start3d, count3d, name='FU')
                end if

                if(nf_selectedvar('FV', output_varnames)) then
                   call nf_put_var(ncdf(ios),var3d(:,2,:,1),start3d, count3d, name='FV')
                end if
                deallocate(var3d)
             end if

#ifdef ENERGY_DIAGNOSTICS
             if(nf_selectedvar('DIFFU', output_varnames) .or. &
                  nf_selectedvar('DIFFV', output_varnames)) then
                allocate(var3d(ncnt,2,nlev,1))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_vector(interpdata(ie), elem(ie), &
                        elem(ie)%accum%DIFF(:,:,:,:), np, nlev, var3d(st:en,:,:,1), 0)
                   st=st+interpdata(ie)%n_interp
                enddo

                if(nf_selectedvar('DIFFU', output_varnames)) then
                   call nf_put_var(ncdf(ios),var3d(:,1,:,1),start3d, count3d, name='DIFFU')
                end if

                if(nf_selectedvar('DIFFV', output_varnames)) then
                   call nf_put_var(ncdf(ios),var3d(:,2,:,1),start3d, count3d, name='DIFFV')
                end if
                deallocate(var3d)
             end if
             if(nf_selectedvar('CONVU', output_varnames) .or. &
                  nf_selectedvar('CONVV', output_varnames)) then
                allocate(var3d(ncnt,2,nlev,1))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_vector(interpdata(ie), elem(ie), &
                        elem(ie)%accum%CONV(:,:,:,:), np, nlev, var3d(st:en,:,:,1), 0)
                   st=st+interpdata(ie)%n_interp
                enddo

                if(nf_selectedvar('CONVU', output_varnames)) then
                   call nf_put_var(ncdf(ios),var3d(:,1,:,1),start3d, count3d, name='CONVU')
                end if

                if(nf_selectedvar('CONVV', output_varnames)) then
                   call nf_put_var(ncdf(ios),var3d(:,2,:,1),start3d, count3d, name='CONVV')
                end if
                deallocate(var3d)
             end if
             if(nf_selectedvar('DIFFT', output_varnames)) then
                allocate(datall(ncnt,nlev))
                st=1
                do ie=nets,nete
                   en=st+interpdata(ie)%n_interp-1
                   call interpolate_scalar(interpdata(ie), elem(ie)%accum%DIFFT, &
                        np, nlev, datall(st:en,:))
                   st=st+interpdata(ie)%n_interp
                enddo
                call nf_put_var(ncdf(ios),datall,start3d, count3d, name='DIFFT')
                deallocate(datall)
             end if
#endif

             ! note: HOMME now compiles both interp_movie_mod.F90 and native grid output modules 
             ! (prim_movie_mod.F90 or shal_movie_mod.F90) into a single executable.
             ! the intent is to make PIO_INTERP a run time variable and remove the #ifdef
             !
             ! However, these two routines are still conditionally compiled for either PIO or PIO_INTERP
             ! and hence must be protected by an #ifdef:
#ifdef PIO_INTERP
             if(test_case.eq.'aquaplanet') then
                call aq_movie_output(ncdf(ios), elem, interpdata, output_varnames, ncnt, nlev)
             end if
             if(columnpackage.ne.'none') then
                call physics_movie_output(ncdf(ios),elem, interpdata, output_varnames, ncnt)
             end if
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            end _PRIM only I/O
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif




             if(piofs%io_rank==0) then
                count2d(3:3)=1
             else
                count2d(3:3)=0
             end if
             if (hybrid%par%masterproc) print *,'writing time...'
             call nf_put_var(ncdf(ios),real(dayspersec*time_at(tl%nstep),kind=real_kind),&
                  start2d(3:3),count2d(3:3),name='time')
             call nf_advance_frame(ncdf(ios))
             call pio_syncfile(ncdf(ios)%fileid)
             if (hybrid%par%masterproc) print *,'finished I/O sync'
          end if ! if output needs to be written
       end if ! output stream is enabled
    end do  ! do ios=1,max_output_streams
    call t_stopf('interp_movie_output')
  end subroutine interp_movie_output





  subroutine GetIODOF(ndims, gdims, iorank, iodof, start, count)
    integer, intent(in) :: ndims
    integer, intent(in) :: gdims(ndims)
    integer, intent(in) :: iorank
    integer(kind=nfsizekind), intent(out) :: start(ndims), count(ndims)
    integer, pointer :: iodof(:) ! gcc4.2 didn't like intent(out)

    integer :: nzrank, nxrank, nx, k, i, j, icnt




    if(iorank>=0) then
       nx=num_io_procs
       count(ndims)=max(1,gdims(ndims)/num_io_procs)
       nx=max(1,num_io_procs/gdims(ndims))
       nzrank=iorank/nx

       k=num_io_procs-gdims(ndims)*nx

       if(iorank>num_io_procs-k-1.and.k>0) then
          nzrank=gdims(ndims)-1
       end if

       start(ndims)=nzrank*count(ndims)+1
       if(gdims(ndims)>num_io_procs) then
          k=gdims(ndims)-num_io_procs*count(ndims)
          if(k>iorank) then
             count(ndims)=count(ndims)+1
          end if
          if(k>=iorank) then
             start(ndims)=start(ndims)+iorank
          else
             start(ndims)=start(ndims)+k
          end if

       end if

       if(k>0 .and.nzrank==gdims(ndims)-1 ) then
          nx=nx+k
       end if
       nxrank=iorank

       do i=ndims-1,1,-1
          !           print *, nxrank, nx

          nxrank=mod(nxrank,nx)
          count(i)=gdims(i)/nx
          k=gdims(i)-count(i)*nx
          if(nxrank<k) then
             count(i)=count(i)+1
             start(i)=count(i)*nxrank+1
          else
             start(i)=count(i)*nxrank+k+1
          end if
          nx=max(1,num_io_procs/gdims(i))

       end do

       icnt=0
       if(ndims.eq.1) then
          allocate(iodof(count(1)))
          do i=1,count(1)
             iodof(i)=start(1)+i-1
          end do
       else if(ndims.eq.2) then
          allocate(iodof(count(1)*count(2)))

          do j=1,count(2)
             do i=1,count(1)
                icnt=icnt+1
                iodof(icnt)=start(1)+gdims(1)*(start(2)-1)+icnt-1
             end do
          end do

       else
          allocate(iodof(count(1)*count(2)*count(3)))

          do k=1,count(3)
             do j=1,count(2)
                do i=1,count(1)
                   icnt=icnt+1
                   iodof(icnt)=start(1)+gdims(1)*(start(2)-1)+ &
                        gdims(1)*gdims(2)*(start(3)-1)+icnt-1
                end do
             end do
          end do
       end if

    else	
       allocate(iodof(1))
       iodof=-1
    end if

  end subroutine GetIODOF






end module interp_movie_mod
