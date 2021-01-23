#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
! This is code to interpolate from the cubed sphere grid to a lat lon grid
! using the cubed sphere basis functions.  
!
! September 2007
! Jim Edwards 




module interpolate_driver_mod
#ifndef HOMME_WITHOUT_PIOLIBRARY
  use pio, only : file_desc_t, var_desc_t , io_desc_t, & ! _EXTERNAL
   pio_get_att, pio_setdebuglevel, pio_closefile, &
   pio_put_att, pio_global, pio_put_var, pio_iotask_rank, &
   pio_read_darray, pio_setframe, pio_get_var, &
   PIO_OFFSET_KIND, pio_char, &
   pio_inq_varid, pio_inq_attname, pio_copy_att, pio_inq_varnatts,&
   pio_inq_att,&
   PIO_MAX_NAME, pio_double

  use pio_nf_utils, only : copy_pio_var ! _EXTERNAL
  use pio_io_mod, only : nfsizekind ! _EXTERNAL
  use common_io_mod, only : nf_handle, max_output_streams, MAX_INFILES,&
              infilenames,piofs, output_type
  use parallel_mod, only : abortmp, syncmp, haltmp
  implicit none
  private
!#include "pnetcdf.inc"
#endif
  public :: interpolate_driver, pio_read_phis, &
       read_gll_topo_file, read_physgrid_topo_file, write_physgrid_topo_file, &
       write_physgrid_smoothed_phis_file

#ifndef HOMME_WITHOUT_PIOLIBRARY
  integer :: nlat, nlon

  type dim_t
     integer :: len
     character(len=PIO_MAX_NAME) :: name
  end type dim_t

  type var_t
     integer, pointer :: ndims(:)
     integer, pointer :: vtype(:)
     integer, pointer :: dimids(:,:)
     character*(PIO_MAX_NAME), pointer :: name(:)
     type(var_desc_t), pointer :: vardesc(:)
     logical, pointer :: timedependent(:)
     logical, pointer :: decomposed(:)
  end type var_t

  type file_t
     type(File_Desc_t) :: FileID
     
     type(var_t) :: vars
     type(dim_t), pointer :: dims(:)
     integer :: unlimid
     integer :: natts
  end type file_t

  integer, parameter :: maxdims = 4
  integer :: ncoldimid=-1

  type(dim_t) :: ewdim(maxdims), nsdim(maxdims), zdim(maxdims), ncoldim(maxdims)
  type(dim_t) :: ewdimi(maxdims), nsdimi(maxdims), zdimi(maxdims)

  logical :: init_done = .false.

  character*(pio_max_name) ewnames (maxdims)
  character*(pio_max_name) nsnames (maxdims)
  character*(pio_max_name) znames (maxdims)
  character*(pio_max_name) ncolnames (maxdims)
  character(len=6) :: ncoldimname='ncol'   ! change to ncol_d if detected
  integer :: newnames, nnsnames, nznames, nncolnames


  type(nf_handle), private, target, save :: ncdf(max_output_streams)
!  type(io_desc_t), pointer :: iodesc3d, iodesc2d, iodesc3dp1
  type(io_desc_t) , save:: iodesc3d, iodesc2d, iodesc3dp1

#endif

contains
  subroutine interpolate_driver(elem,hybrid)
    use hybrid_mod, only : hybrid_t
    use element_mod, only : element_t
#ifndef HOMME_WITHOUT_PIOLIBRARY
    use dimensions_mod, only : ne, nelem, np, nlev
    use common_io_mod, only : varnames=>output_varnames1, nf_handle
    !  use interpolate_mod
    use dof_mod
    use interpolate_mod, only : interpdata_t
    !  use domain_mod, only : domain1d_t, decompose
    !  use thread_mod, only : omp_get_thread_num
    !  use reduction_mod, only : reductionbuffer_ordered_1d_t
#endif
    implicit none

    type(element_t) :: elem(:)
    type(hybrid_t), intent(in) :: hybrid

#ifndef HOMME_WITHOUT_PIOLIBRARY
    integer, parameter :: maxvars=30
    type(file_t) :: infile
    !  type(element_t), allocatable :: elem(:)
    !  type (domain1d_t), allocatable:: dom_mt(:)
    !  type (parallel_t) :: par
    !  type (ReductionBuffer_ordered_1d_t) :: red

    type(interpdata_t), pointer :: interpdata(:)
    type(nf_handle), pointer :: outfile
    integer :: ret, ie
    integer, allocatable :: gdof(:,:,:)
    character(len=160), pointer :: infilename 
    integer :: i
!
! For each infilename in the namelist, read and interpolate to an outputfile.    
! The infile is expected to have np and ne (the cubed sphere horizontal dimensions) 
! defined as global attributes
!
    do i=1,MAX_INFILES
       if(len(trim(infilenames(i)))>0) then
          infilename=>infilenames(i)
          call infile_initialize(elem, hybrid%par,infilename, varnames, infile)
          call outfile_initialize(elem, hybrid, infile, outfile, interpdata, infilename)
          if (hybrid%par%masterproc) print *,'interpolating...'
          call interpolate_vars(infile, outfile, elem, hybrid%par, interpdata)
          call pio_closefile(infile%fileid)
          call pio_closefile(outfile%fileid)
          call free_infile(infile)
          deallocate(interpdata)
       end if
    end do
#endif
  end subroutine interpolate_driver

#ifndef HOMME_WITHOUT_PIOLIBRARY
  subroutine free_infile(infile)
    type(file_t), intent(inout) :: infile

    deallocate(infile%vars%vardesc)
    nullify(infile%vars%vardesc)
    deallocate(infile%vars%vtype)
    nullify(infile%vars%vtype)
    deallocate(infile%vars%name)
    nullify(infile%vars%name)
    deallocate(infile%vars%ndims)
    nullify(infile%vars%ndims)
    deallocate(infile%vars%timedependent)
    nullify(infile%vars%timedependent)

    deallocate(infile%vars%decomposed)
    nullify(infile%vars%decomposed)
    deallocate(infile%vars%dimids)
    nullify(infile%vars%dimids)

  end subroutine free_infile

  subroutine infile_initialize(elem, par, infilename, varnames, infile)
    use pio ! _EXTERNAL
    use common_io_mod, only : io_stride, num_io_procs, num_agg
    use parallel_mod, only : parallel_t, mpireal_t
    use interp_movie_mod, only : getiodof
    use element_mod, only: element_t
    use dimensions_mod, only : ne, np, nlev,nelem
    type(element_t) :: elem(:)
    type(parallel_t) :: par
    character(len=*), intent(in) :: infilename
    character(len=*), intent(inout) :: varnames(:)
    type(file_t), intent(out) :: infile
    type(io_desc_t), pointer :: iodesc
    character(len=80), allocatable :: varnames_list(:)
    
    integer :: ndims, nvars,  varcnt
    integer :: ret, i, ncid, n
    integer :: varid, ncols, lev
    integer :: ne_file, np_file, nlev_file
    
    integer*8, pointer :: ldof(:)
    integer(kind=PIO_OFFSET_KIND) :: start(3), count(3)
    integer :: iorank, dimcnt
    character(len=80) :: name

    if (par%masterproc) print *,'initializing input file: ',trim(infilename)

!    call pio_setdebuglevel(1)
    if(output_type.eq.'netcdf') then
       ret = PIO_OpenFile(PIOFS, InFile%FileID, iotype_netcdf, infilename)
    else
       ret = PIO_OpenFile(PIOFS, InFile%FileID, iotype_pnetcdf, infilename)
    end if

    ret = PIO_inquire(InFile%FileID, nDimensions=ndims, nVariables=nvars, &
         nAttributes=infile%natts, unlimitedDimId=InFile%unlimid)

!    print *, __FILE__,__LINE__,trim(infilename), ndims, nvars, infile%natts, infile%unlimid

    allocate(infile%dims(ndims))

    nlev_file=-1

    ! is this an output file with "ncol", or "ncol_d"?
    ! set ncoldimname before calling "is_ncolfile"
    do i=1,ndims
       ret = PIO_inq_dimname(Infile%fileid, i, infile%dims(i)%name)
       if(infile%dims(i)%name.eq.'ncol_d') ncoldimname='ncol_d'
    end do
    if (par%masterproc) print *,'interpolating dimension = ',ncoldimname

    ncols=0
    do i=1,ndims
       ret = PIO_inq_dimname(Infile%fileid, i, infile%dims(i)%name)
       ret = PIO_inq_dimlen(Infile%fileid, i, infile%dims(i)%len)
       if(is_ncoldim(infile%dims(i)%name)) ncols=infile%dims(i)%len
       if(infile%dims(i)%name.eq.'lev') nlev_file = infile%dims(i)%len
    end do
    if(trim(varnames(1)) == 'all') then
       varcnt=nvars
    else
       varcnt=0
       allocate(varnames_list(size(varnames)))
       do i=1,size(varnames)
          if(varnames(i) == '') exit
          ret = PIO_inq_varid(InFile%FileID, varnames(i), varid)
          if (ret==PIO_NOERR) then
             varcnt=varcnt+1
             varnames_list(varcnt)=varnames(i)
          else
             if (par%masterproc) print *,'skipping requested variable=',trim(varnames(i))
          endif
       end do
       do i=1,varcnt
          varnames(i)=varnames_list(i)
       enddo
       deallocate(varnames_list)
    end if

    if(is_ncolfile(infile%dims)) then
       np_file=-1
       ne_file=-1
       do i=1,infile%natts
          ret = pio_inq_attname(infile%FileID, PIO_GLOBAL, i, name)
          if (trim(name)=='np') ret=PIO_get_att(infile%FileID, PIO_GLOBAL, 'np', np_file)
          if (trim(name)=='ne') ret=PIO_get_att(infile%FileID, PIO_GLOBAL, 'ne', ne_file)
       enddo

       ! abort if file attribute doesn't match model resolution
       ! some files, like topo files, dont have these attributes.
       if (ncols/=2+nelem*(np-1)**2) then
          print *,'nelem,np,ncols',nelem,np,ncols
          call abortmp('File resolution incoorect: ncols <> 2+nelem*(np-1)^2')
       endif

       if(ne_file/=ne .and. ne_file/=-1) then
          print *,'ne, ne_file',ne,ne_file
          call abortmp('The variable ne in the namelist must be the same as that of the file.')
       end if
       if(np_file/=np .and. np_file/=-1) then
          print *,'np, np_file',np,np_file
          call abortmp('The variable dimensions_mod::np must be the same as that of the file, you will need to recompile.')
       end if
    else
       call abortmp('The input file is missing required ncol dimensions')
    end if

! define the decompositions...
    iorank=pio_iotask_rank(piofs)
    call getcompdof(ldof, elem, 1)
    call PIO_initDecomp(piofs,pio_double,(/ncols/),ldof, iodesc2d)
    deallocate(ldof)

    call getcompdof(ldof, elem, nlev)
    call PIO_initDecomp(piofs,pio_double,(/ncols,nlev/), &
         ldof, iodesc3d)

    deallocate(ldof)
    

    call getcompdof(ldof, elem, nlev+1)
    call PIO_initDecomp(piofs,pio_double,(/ncols,nlev+1/), &
         ldof, iodesc3dp1)

    deallocate(ldof)

    allocate(infile%vars%vardesc(varcnt))   
    allocate(infile%vars%vtype(varcnt))
    allocate(infile%vars%name(varcnt))   
    allocate(infile%vars%ndims(varcnt))   
    allocate(infile%vars%timedependent(varcnt))   
    allocate(infile%vars%decomposed(varcnt))   
    allocate(infile%vars%dimids(ndims,varcnt))   
    

    infile%vars%timedependent(1:varcnt)=.false.
    infile%vars%decomposed(1:varcnt)=.false.

! Get the name and dimension information of each variable in the input file,
! determine if the variable needs interpolation and if it has a time dimension
    do i=1,varcnt
       if(trim(varnames(1)) == 'all') then
          infile%vars%vardesc(i)%varid=i
          ret = PIO_inq_varname(InFile%FileID, infile%vars%vardesc(i), infile%vars%name(i))
       else
          if (par%masterproc) print *,'inquiring for variable: ',varnames(i)
          ret = PIO_inq_varid(InFile%FileID, varnames(i), infile%vars%vardesc(i))
          infile%vars%name(i)=varnames(i)
       end if

       ret = PIO_inq_vartype(infile%FileID, infile%vars%vardesc(i),  &
            infile%vars%vtype(i))
       ret = PIO_inq_varndims(infile%FileID, infile%vars%vardesc(i), &
            infile%vars%ndims(i))
       ret = PIO_inq_vardimid(infile%FileID, infile%vars%vardesc(i), &
            infile%vars%dimids(1:infile%vars%ndims(i),i))
       lev=1
       do n=1,infile%vars%ndims(i)
          if(infile%dims(infile%vars%dimids(n,i))%name.eq.'lev') then
             ! if we are reading a variable with levels, verifty nlev == nlev_file
             if(nlev_file/=nlev .and. nlev_file/=-1) then
                print *,'nlev, nlev_file',nlev,nlev_file
                call abortmp('Error: dimensions_mod::nlev does not match file nlev')
             end if
             lev=nlev
          endif
          if(infile%dims(infile%vars%dimids(n,i))%name.eq.'ilev') then
             ! if we are reading a variable with levels, verifty nlev == nlev_file
             if(nlev_file/=nlev .and. nlev_file/=-1) then
                print *,'nlev, nlev_file',nlev,nlev_file
                call abortmp('Error: dimensions_mod::nlev does not match file nlev')
             end if
             lev=nlev+1
          endif
          if(infile%dims(infile%vars%dimids(n,i))%name.eq.'time') infile%vars%timedependent(i)=.true.
          if(is_ncoldim(infile%dims(infile%vars%dimids(n,i))%name)) infile%vars%decomposed(i)=.true.
       end do
    end do
    if (par%masterproc) print *,'input file: nlev,ne=',nlev,ne
  end subroutine infile_initialize

! Open the output file and go through all of the netcdf definition steps.
! variables that had the ncols dimension in the input file have lon,lat dimensions in the 
! output file.

  subroutine outfile_initialize(elem, hybrid, infile, outfile, interpdata, infilename )
    use kinds, only : real_kind
    use dimensions_mod, only : nelemd, nelem, nlev
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t
    use interpolate_mod, only : interpdata_t, setup_latlon_interp, &
         get_interp_lat, get_interp_lon, get_interp_gweight, get_interp_parameter

    use parallel_mod, only : abortmp
    use common_io_mod, only : output_frequency, output_start_time, &
         output_end_time
    use pio_io_mod, only : nf_output_init_begin, nf_output_init_complete, & ! _EXTERNAL
         nf_output_register_dims, nf_output_register_variables, get_varindex
    implicit none
    type(element_t) :: elem(:)
    type(hybrid_t) :: hybrid
    type(file_t) :: infile

    type(interpdata_t), pointer :: interpdata(:)
    type(nf_handle), pointer :: outfile

    integer :: ncoldimid, londimid, latdimid, i, j, k
    integer :: ofdims, nvars, ierr
    character*(pio_max_name), allocatable :: dimnames(:)
    integer , allocatable :: vardims(:,:), dimsize(:), otype(:)
    logical, allocatable :: varrequired(:)

    real(kind=real_kind), allocatable :: lon(:), lat(:), gw(:), lev(:), ilev(:)
    integer(kind=nfsizekind) :: count1d(1), start1d(1)
    integer :: ios, iorank, ndnt, vindex
    character(len=*) :: infilename


    output_frequency(2:max_output_streams)=0
    output_frequency(1)=1
    output_start_time(1)=0
    output_end_time(1)=1

    allocate(interpdata(nelemd))
    call setup_latlon_interp(elem, interpdata, hybrid%par)
    nlat = get_interp_parameter('nlat')
    nlon = get_interp_parameter('nlon')

    if(is_ncolfile(infile%dims)) then
       ofdims = size(infile%dims) + 1
    else
       call abortmp('this branch is not yet supported')
    end if
    allocate(dimnames(ofdims))
    allocate(dimsize(ofdims))
    nvars = size(infile%vars%name)


    allocate(varrequired(nvars))
    allocate(vardims(ofdims,nvars))
    allocate(otype(nvars))
    varrequired=.true.

    do i=1,ofdims
       if(i==ofdims) then
          dimsize(i)=nlat
          dimnames(i)='lat'
          latdimid=i
       else if(is_ncoldim(infile%dims(i)%name)) then
          dimnames(i)='lon'
          londimid=i
          ncoldimid=i
          dimsize(i)=nlon
       else 
          dimnames(i)=infile%dims(i)%name
          if(i==infile%unlimid) then
             dimsize(i)=0
          else
             dimsize(i)=infile%dims(i)%len
          end if
       end if
    end do

    if (hybrid%masterthread) print *,'Examining input variables:'
    vardims=0
    do i=1,nvars
       k=1
       otype(i)=infile%vars%vtype(i)
       do j=1,infile%vars%ndims(i)
         if(infile%vars%dimids(j,i).eq.ncoldimid) then
            vardims(k,i) = londimid
            vardims(k+1,i) = latdimid
            k=k+2
          else
             vardims(k,i)=infile%vars%dimids(j,i)
             k=k+1
          end if
       end do
       if (hybrid%masterthread) then
          print *,'var=',trim(infile%vars%name(i)),' remap=',infile%vars%decomposed(i)
       endif

       ! if the file happens to have a variable 'lat','lon' or 'gw', rename so it doesnt conflict
       ! with the interpolated output
       if (trim(infile%vars%name(i))=='lat')  infile%vars%name(i)='lat_orig'
       if (trim(infile%vars%name(i))=='lon')  infile%vars%name(i)='lon_orig'
       if (trim(infile%vars%name(i))=='gw')  infile%vars%name(i)='gw_orig'
    end do

    ndnt = index(infilename,".nc")

    call nf_output_init_begin(ncdf,hybrid%par%masterproc,hybrid%par%nprocs,hybrid%par%rank, &
         hybrid%par%comm, infilename(1:ndnt)//'interp'   ,0)
    
    call nf_output_register_dims(ncdf, ofdims, dimnames, dimsize)
    if (hybrid%par%masterproc) print *,'creating PIO decompositions'
    call create_output_decomps(ncdf, interpdata, nlon, nlat)
    if (hybrid%par%masterproc) then
      print *,'registering output variables -- '//&
           'if you get a seg fault, try increasing the number of processors used to help memory management'
    end if
    call nf_output_register_variables(ncdf, nvars, infile%vars%name, &
         vardims, otype, varrequired)

    vardims(1,1)=latdimid
    call nf_output_register_variables(ncdf, 1, (/'gw'/), vardims(1:1,:), &
         (/pio_double/), (/.true./))
    ! if the user did not specify lat,lon, be sure to add them:
    call nf_output_register_variables(ncdf, 1, (/'lat'/), vardims(1:1,:), &
         (/pio_double/), (/.true./))
    vardims(1,1)=londimid
    call nf_output_register_variables(ncdf, 1, (/'lon'/), vardims(1:1,:), &
         (/pio_double/), (/.true./))


    deallocate(dimnames, dimsize, varrequired, vardims, otype)

    outfile=>ncdf(1)
    iorank=pio_iotask_rank(piofs)
    
    if (hybrid%par%masterproc) print *,'copying attributes'
    call copy_attributes(infile, outfile)

    if (hybrid%par%masterproc) print *,'output file attribute: ',&
	'interpolated from file '//trim(infilename)
    ierr= PIO_Put_att(outfile%fileid, PIO_GLOBAL, 'history',&
         'interpolated from file '//trim(infilename))

    call nf_output_init_complete(ncdf)

! write the new lat and lon arrays, and possibly gauss weights as well

    if(iorank.eq.0) then
       allocate(lon(nlon), lat(nlat), gw(nlat))
       lon = get_interp_lon()
       lat = get_interp_lat()
    else
       allocate(lon(1), lat(1), gw(1))
    end if

    vindex = get_varindex('lon', outfile%varlist)
    if(vindex>=0) then
       ierr = pio_put_var(outfile%fileid, outfile%varlist(vindex)%ivarid, lon)   
    endif
    vindex = get_varindex('lat', outfile%varlist)
    if(vindex>=0) then
       ierr = pio_put_var(outfile%fileid, outfile%varlist(vindex)%ivarid, lat)   
    end if

    ! output latitude weights
    if(iorank.eq.0) then
       gw = get_interp_gweight()
    end if
    vindex = get_varindex('gw', outfile%varlist)
    if(vindex>=0) then
       ierr = pio_put_var(outfile%fileid, outfile%varlist(vindex)%ivarid, gw)   
    end if

    deallocate(lon, lat, gw)
  end subroutine outfile_initialize

! interpolate all variables from the input file to the output file.   
! variables which are not decomposed are simply copied.
!
  subroutine interpolate_vars(infile, outfile, elem, par, interpdata)
    use dof_mod, only : putuniquepoints
    use kinds, only : real_kind
    use pio_io_mod, only : nf_handle, nf_put_var => nf_put_var_pio! _EXTERNAL
    use interpolate_mod, only : interpdata_t, interpolate_scalar, interpolate_vector, &
         get_interp_parameter,var_is_vector_vvar,var_is_vector_uvar, replace_vec_by_vordiv
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
    use edgetype_mod, only : edgebuffer_t
    use dimensions_mod, only : nelemd, nlev, np
    use parallel_mod, only : parallel_t, syncmp
    use bndry_mod, only : bndry_exchangeV
    use common_io_mod, only : readystate, nf_variable
    use viscosity_mod, only : compute_zeta_c0, compute_div_c0

    type(element_t), intent(inout) :: elem(:)
    type(parallel_t),intent(in) :: par
    type(file_t), intent(inout) :: infile
    type(nf_handle), intent(inout) :: outfile
    type(interpdata_t), intent(in) :: interpdata(:)

    type(edgeBuffer_t) :: edge    
    real(kind=real_kind), pointer :: array(:,:), varray(:,:,:)
    real(kind=real_kind), pointer :: farray(:)	, fvarray(:,:), ftmp(:,:), fvtmp(:,:,:)
    real(kind=real_kind), pointer :: zeta(:,:,:,:),div(:,:,:,:)

    integer(kind=nfsizekind) :: start(3), count(3)
    integer :: starti(3), counti(3), index(1)
    integer :: len, ii, iu, iv, k, it, nd, offset, iuvar
    integer :: nvars, i, j, n, ncnt_out, ncnt_in, st, en, ncols, ierr, ie, lev, ntimes, ndims
    character*(PIO_MAX_NAME) :: tmpname
    integer(kind=nfsizekind) :: start2d(3), count2d(3), start3d(4), count3d(4)

                   


    if(outfile%state .ne. readystate) then
       print *,__FILE__,__LINE__,outfile%state
       call abortmp('outfile not in readystate')
    end if

    nvars = size(infile%vars%name)

    ncnt_in = sum(elem(1:nelemd)%idxp%numUniquePts)
    ncnt_out = sum(interpdata(1:nelemd)%n_interp)

    call initedgebuffer(par,edge,elem,2*nlev)


    VARLOOP: do i=nvars,1,-1
       if(trim(infile%vars%name(i)).eq.'lat' .or. trim(infile%vars%name(i)).eq.'lon' .or. &
            trim(infile%vars%name(i)).eq.'lat_d' .or. trim(infile%vars%name(i)).eq.'lon_d'    ) then
          ! do nothing
       else 
          lev=1
          do n=1,infile%vars%ndims(i)
             if(infile%dims(infile%vars%dimids(n,i))%name.eq.'lev') lev=nlev
             if(infile%dims(infile%vars%dimids(n,i))%name.eq.'ilev') lev=nlev+1
          end do
          ndims = infile%vars%ndims(i) 

          ntimes=1
          if(infile%vars%timedependent(i)) then
             ntimes = get_dimlen(infile%dims,'time')
          end if
          if(infile%vars%decomposed(i)) then
             start3d=1
             count3d(1)=nlon
             count3d(2)=nlat
             count3d(3)=lev
             count3d(4)=1
             if(par%masterproc) print *,'considering: ',trim(infile%vars%name(i))

             if (var_is_vector_Vvar(infile%vars%name(i))/=0) then
                if(par%masterproc) print *,'skipping V component- will be interpolated with U component'
                cycle VARLOOP
             endif

             iuvar=var_is_vector_Uvar(infile%vars%name(i))  ! index in vector_uvar()
             iv=0                                           ! index in infile%vars%name() 
             if (iuvar>0) then
                ! i = index of u componnet, now find iv = index of v component
                VECVARLOOP: do iv=1,nvars
                   if(iuvar==var_is_vector_Vvar(infile%vars%name(iv))) exit VECVARLOOP
                end do VECVARLOOP
                if (iv>nvars) then
                   if(par%masterproc) print *,'Error finding V component. skipping...'
                   cycle VARLOOP
                endif
             endif

             ! 3 cases:
             ! iv=0         interpolation variable is not a vector
             ! iv>nvars     interpolation variable is vector, but V component NOT FOUND
             ! 0<iv<=nvars  U and V components identified
             if(iv>0) then
                if(par%masterproc) print *,'Interpolating vector (',trim(infile%vars%name(i)), ',', &
                     trim(infile%vars%name(iv)),')'
                if(par%masterproc) then
                   if (get_interp_parameter('itype')==0) print *,'using native spectral element interpolation'
                   if (get_interp_parameter('itype')==1) print *,'using bilinear interpolation'
                endif
                do n=1,ntimes
                   if(par%masterproc) print *, 'interpolating for timelevel ',n,'/',ntimes
                   if(infile%vars%timedependent(i)) then
                      call pio_setframe(infile%FileID, infile%vars%vardesc(i),int(n,kind=PIO_OFFSET_KIND))
                      call pio_setframe(infile%FileID, infile%vars%vardesc(iv),int(n,kind=PIO_OFFSET_KIND))
                   end if
                   allocate(fvarray(ncnt_in*lev,2))
                   fvarray=0.0d0
                   if(lev==1) then
                      call pio_read_darray(infile%FileID, infile%vars%vardesc(i), iodesc2d, fvarray(:,1), ierr)
                      call pio_read_darray(infile%FileID, infile%vars%vardesc(iv), iodesc2d, fvarray(:,2), ierr)
                   else if(lev==nlev) then
                      call pio_read_darray(infile%FileID, infile%vars%vardesc(i),  iodesc3d, fvarray(:,1), ierr)
                      call pio_read_darray(infile%FileID, infile%vars%vardesc(iv), iodesc3d, fvarray(:,2), ierr)
                   else
                      call pio_read_darray(infile%FileID, infile%vars%vardesc(i),  iodesc3dp1, fvarray(:,1), ierr)
                      call pio_read_darray(infile%FileID, infile%vars%vardesc(iv), iodesc3dp1, fvarray(:,2), ierr)
                   end if
                   offset=0
                   do ie=1,nelemd
                      elem(ie)%state%v(:,:,:,:,1) = 0.0d0
	              allocate(fvtmp(elem(ie)%idxP%NumUniquePts, 2, lev))
                      do k=1,lev
                         do ii=1,elem(ie)%idxP%NumUniquePts
                            !it=(elem(ie)%idxp%UniquePtOffset+ii+(k-1)*ncnt_in)-1
                            it=(offset+ii+(k-1)*ncnt_in)
                            fvtmp(ii,1,k) = fvarray(it,1)
                            fvtmp(ii,2,k) = fvarray(it,2)
                         end do
                      end do
                      offset = offset+elem(ie)%idxP%NumUniquePts
                      call putUniquePoints(elem(ie)%idxP, 2, lev, &
                           fvtmp,elem(ie)%state%v(:,:,:,:,1))
                      deallocate(fvtmp)
                      call edgevpack(edge, elem(ie)%state%v(:,:,:,:,1),2*lev,0,ie)
                   end do
                   deallocate(fvarray)

                   call bndry_exchangeV(par, edge)
                   do ie=1,nelemd
                      call edgeVunpack(edge, elem(ie)%state%v(:,:,:,:,1),2*lev,0,ie)
                   enddo

                   ! hack to get native vorticity/divergence
                   if ( replace_vec_by_vordiv(iuvar) ) then
                      if(par%masterproc) print *,'VORDIV flag set. Overwriting (',&
                           trim(infile%vars%name(i)), ',',trim(infile%vars%name(iv)),')',&
                           ' with (vor,div)'
                      allocate(zeta(np,np,nlev,nelemd))
                      allocate(div(np,np,nlev,nelemd))
                      call compute_zeta_c0(zeta,elem,par,1)
                      call compute_div_c0(div,elem,par,1)    
                      ! overwrite velocity with result:
                      do ie=1,nelemd
                         elem(ie)%state%v(:,:,1,:,1)=zeta(:,:,:,ie)
                         elem(ie)%state%v(:,:,2,:,1)=div(:,:,:,ie)
                      enddo
                      deallocate(zeta)
                      deallocate(div)
                   endif

                   st = 1
                   allocate(varray(ncnt_out,lev,2))
                   do ie=1,nelemd
                      en=st+interpdata(ie)%n_interp-1
                      call interpolate_vector(interpdata(ie), elem(ie), &
                           elem(ie)%state%V(:,:,:,1:lev,1), lev, varray(st:en,:,:), 0)
                      st=st+interpdata(ie)%n_interp
                   end do
                   nd=2
                   if(lev>1) nd=nd+1
                   if(infile%vars%timedependent(i)) then
                      nd=nd+1
                      start3d(nd)=n
                      count3d(nd)=1
                   end if
                   if(lev.eq.1) then
                      call nf_put_var(outfile, varray(:,1,1), start3d(1:nd), count3d(1:nd), name=infile%vars%name(i))
                      call nf_put_var(outfile, varray(:,1,2), start3d(1:nd), count3d(1:nd), name=infile%vars%name(iv))
                   else
                      call nf_put_var(outfile, varray(:,:,1), start3d(1:nd), count3d(1:nd), name=infile%vars%name(i))
                      call nf_put_var(outfile, varray(:,:,2), start3d(1:nd), count3d(1:nd), name=infile%vars%name(iv))
                   end if
                   
                   deallocate(varray)
                end do

             else 
                if(par%masterproc) print *,'Interpolating ',trim(infile%vars%name(i)),' nlev=',lev
                if(par%masterproc) then
                   if (get_interp_parameter('itype')==0) print *,'using native spectral element interpolation'
                   if (get_interp_parameter('itype')==1) print *,'using bilinear interpolation'
                endif
! Interpolate a scalar field.   
                start = 1
                count = 1
                count(1) = ncnt_in
                if(lev>1) then
                   count(2)=lev
                endif
               
                do n=1,ntimes
                   if(par%masterproc) print *, 'interpolating for timelevel ',n,'/',ntimes
                   if(infile%vars%timedependent(i)) then
                      call pio_setframe(infile%FileID, infile%vars%vardesc(i),int(n,kind=PIO_OFFSET_KIND))
                   end if
                   allocate(farray(ncnt_in*lev))
                   farray = 1.0e-37
                   if(lev==1) then
                      call pio_read_darray(infile%FileID, infile%vars%vardesc(i), iodesc2d, farray, ierr)
                   else if(lev==nlev) then
                      call pio_read_darray(infile%FileID, infile%vars%vardesc(i), iodesc3d, farray, ierr)
                   else
                      call pio_read_darray(infile%FileID, infile%vars%vardesc(i), iodesc3dp1, farray, ierr)
                   end if
                   offset=0
                   allocate(ftmp(np*np, lev))
                   do ie=1,nelemd
                      do k=1,lev
                         do ii=1,elem(ie)%idxP%NumUniquePts
                            !iv=(elem(ie)%idxp%UniquePtOffset+ii+(k-1)*ncnt_in)-1
                            iv=(offset+ii+(k-1)*ncnt_in)
                            ftmp(ii,k) = farray(iv)
                         end do
                      end do
                      offset = offset+elem(ie)%idxP%NumUniquePts
                      elem(ie)%state%Q(:,:,:,1) = 0.0d0
                      call putUniquePoints(elem(ie)%idxP, lev, ftmp, elem(ie)%state%Q(:,:,1:lev,1))
                      call edgevpack(edge, elem(ie)%state%Q(:,:,:,1),lev,0,ie)
                   end do
                   deallocate(ftmp)

                   call bndry_exchangeV(par, edge)
                   st = 1
                   deallocate(farray)
                   allocate(array(ncnt_out,lev))
                   array=0
                   do ie=1,nelemd
                      en=st+interpdata(ie)%n_interp-1
                      call edgeVunpack(edge, elem(ie)%state%Q(:,:,:,1),lev,0,ie)
                      
                      call interpolate_scalar(interpdata(ie), elem(ie)%state%Q(:,:,1:lev,1), &
                           np, lev, array(st:en,:))

                      st=st+interpdata(ie)%n_interp
                   end do
                   
                   nd=2
                   if(lev>1) nd=nd+1
                   if(infile%vars%timedependent(i)) then
                      nd=nd+1
                      start3d(nd)=n
                      count3d(nd)=1
                   end if
                   if(lev.eq.1) then
                      call nf_put_var(outfile, array(:,1), start3d(1:nd), count3d(1:nd), name=infile%vars%name(i))
                   else
                      call nf_put_var(outfile, array, start3d(1:nd), count3d(1:nd), name=infile%vars%name(i))
                   end if
                   deallocate(array)
                end do
             end if

          else
             if(par%masterproc) print *,'Copying ',trim(infile%vars%name(i))
             ! copy non-decomposed data
             len =1
             if( outfile%varlist(i)%vtype .eq. PIO_Char) then
                if (infile%vars%ndims(i)>=2) then
                   ! SCORPIO chokes on these.
                   if(par%masterproc) print *,"SKIPPING. cant copy 2D char array"
                else   
                   do n=2,infile%vars%ndims(i)
                      len = len*infile%dims(infile%vars%dimids(n,i))%len
                   end do
                   call copy_pio_var(infile%FileID, Outfile%fileid, infile%vars%vardesc(i), &
                        outfile%varlist(i)%vardesc, infile%dims(infile%vars%dimids(1,i))%len, len)
                endif
             else
                if(infile%vars%ndims(i)<=1) then
                   do n=1,infile%vars%ndims(i)
                      len = len*infile%dims(infile%vars%dimids(n,i))%len
                   end do
                   if(infile%vars%timedependent(i)) then
                      ! bug in copy_pio_var - it will not copy the unlim dimension
                      ! workaround: (assumes variable can be read into a real*8)
                      allocate(farray(len))
                      ierr = pio_get_var(infile%FileID, infile%vars%vardesc(i), farray)
                      do j=1,len
                         index(1)=j
                         ierr = pio_put_var(Outfile%fileid, outfile%varlist(i)%vardesc,index, farray(j))
                      enddo
                      deallocate(farray)
                   else
                      call copy_pio_var(infile%FileID, Outfile%fileid, infile%vars%vardesc(i), &
                           outfile%varlist(i)%vardesc, len)
                   endif
                else
                   do n=1,infile%vars%ndims(i)
                      counti(n) = infile%dims(infile%vars%dimids(n,i))%len
                   end do
                   if (infile%vars%ndims(i) <= 2 ) then
                      call copy_pio_var(infile%FileID, Outfile%fileid, infile%vars%vardesc(i), &
                           outfile%varlist(i)%vardesc, counti(1:infile%vars%ndims(i)))
                   else
                      ! copy_pio_var only copies variables with up to 2 dimensions
                      if(par%masterproc) print *,"SKIPPING copy. ndims>2.  ndims=",ndims
                   endif
                end if
             end if
          end if
       end if
       call syncmp(par)
    end do VARLOOP

    call freeedgebuffer(edge)
    
  end subroutine interpolate_vars
#endif
! read a variable from a file
!
! if we ever need to read something other than PHIS, this routine should
! be replaced with a more general routine to read any field
  subroutine pio_read_phis(elem, par, varname)
    use element_mod, only : element_t
    use parallel_mod, only : parallel_t, syncmp
#ifndef HOMME_WITHOUT_PIOLIBRARY
    use dof_mod, only : putuniquepoints
    use kinds, only : real_kind
    use edge_mod, only : edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
    use edgetype_mod, only : edgebuffer_t
    use dimensions_mod, only : nelemd, nlev, np
    use bndry_mod, only : bndry_exchangeV
    use common_io_mod, only : varname_len,infilenames
#endif
    type(element_t), intent(inout) :: elem(:)
    type(parallel_t),intent(in) :: par
    character(*), intent(in), optional :: varname
#ifndef HOMME_WITHOUT_PIOLIBRARY
    ! local
    character(len=varname_len), dimension(1) :: varnames
    type(file_t)     :: infile
    type(edgeBuffer_t) :: edge    
    real(kind=real_kind), allocatable :: farray(:)
    real(kind=real_kind), allocatable :: ftmp(:,:)

    integer :: ii,k,ie,ilev,iv, ierr,offset
    integer :: ncnt_in

    ilev=1
    call initedgebuffer(par,edge,elem,ilev)
    ncnt_in = sum(elem(1:nelemd)%idxp%numUniquePts)

    varnames(1)="PHIS"
    if (present(varname)) varnames(1) = trim(varname)
    call infile_initialize(elem, par,infilenames(1), varnames, infile)


    allocate(farray(ncnt_in))
    farray = 1.0e-37
    call pio_read_darray(infile%FileID, infile%vars%vardesc(1), iodesc2d, farray, ierr)
    !call pio_read_darray(infile%FileID, infile%vars%vardesc(i), iodesc3d, farray, ierr)
    !call pio_read_darray(infile%FileID, infile%vars%vardesc(i), iodesc3dp1, farray, ierr)

    offset=0
    do ie=1,nelemd
       allocate(ftmp(elem(ie)%idxP%NumUniquePts, ilev))
       do k=1,ilev
          do ii=1,elem(ie)%idxP%NumUniquePts
             iv=(offset+ii+(k-1)*ncnt_in)
             ftmp(ii,k) = farray(iv)
          end do
       end do
       offset = offset+elem(ie)%idxP%NumUniquePts
       elem(ie)%state%phis(:,:)=0
       call putUniquePoints(elem(ie)%idxP, ftmp(:,1), elem(ie)%state%phis(:,:))
       call edgevpack(edge, elem(ie)%state%phis(:,:),ilev,0,ie)
       deallocate(ftmp)
    end do
    call bndry_exchangeV(par, edge)
    do ie=1,nelemd
       call edgeVunpack(edge, elem(ie)%state%phis(:,:),ilev,0,ie)
    end do


    deallocate(farray)
    call freeedgebuffer(edge)

    call pio_closefile(infile%fileid)
    call free_infile(infile)
#endif
  end subroutine pio_read_phis
  
!
! Create the pio decomps for the output file.
!
#ifndef HOMME_WITHOUT_PIOLIBRARY
  subroutine create_output_decomps(ncdf,interpdata,nlon,nlat)
    use pio_io_mod ! _EXTERNAL
    use dimensions_mod
    use interpolate_mod
    use interp_movie_mod, only : getiodof
    implicit none
    type(nf_handle) :: ncdf(:)
    type(interpdata_t):: interpdata(:)
    integer, intent(in) :: nlon, nlat
    integer*8, pointer :: ldof2d(:),ldof3d(:)

    integer :: icnt, ie, i, lcount
    integer :: k, iorank
    integer(kind=nfsizekind) :: start1d(1), count1d(1)
    integer :: londim, latdim, timedim, levdim, ilevdim
    integer(kind=nfsizekind) :: start2d(3), count2d(3), start3d(4), count3d(4)
    integer*8 :: nlon8,nlat8

    ! files might only have 2D data, or only lev data
    levdim=0
    ilevdim=0
    do i=1,size(ncdf(1)%dimlist)
       if(ncdf(1)%dimlist(i)%dimname.eq.'lon') londim=ncdf(1)%dimlist(i)%dimID
       if(ncdf(1)%dimlist(i)%dimname.eq.'lat') latdim=ncdf(1)%dimlist(i)%dimID
       if(ncdf(1)%dimlist(i)%dimname.eq.'lev') levdim=ncdf(1)%dimlist(i)%dimID
       if(ncdf(1)%dimlist(i)%dimname.eq.'time') timedim=ncdf(1)%dimlist(i)%dimID
       if(ncdf(1)%dimlist(i)%dimname.eq.'ilev') ilevdim=ncdf(1)%dimlist(i)%dimID
    end do

    lcount = sum(interpdata(1:nelemd)%n_interp)
    iorank = pio_iotask_rank(piofs)

    nlon8=nlon ! force all calculations involving these dims to i*8
    nlat8=nlat 

    ! Create the DOF arrays
    allocate(ldof2d(lcount))
    icnt=0
    do ie=1,nelemd
       do i=1,interpdata(ie)%n_interp
          icnt=icnt+1
          ldof2d(icnt)=interpdata(ie)%ilon(i)+(interpdata(ie)%ilat(i)-1)*nlon8
       end do
    end do
    call getiodof(2, (/nlon,nlat/), iorank, start2d(1:2), count2d(1:2))
    call nf_init_decomp(ncdf, (/londim,latdim/), ldof2d, start2d(1:2),count2d(1:2))
    deallocate(ldof2d)

    if (levdim>0) then
    allocate(ldof3d(lcount*nlev))
    icnt=0
    do k=1,nlev
       do ie=1,nelemd
          do i=1,interpdata(ie)%n_interp
             icnt=icnt+1
             ldof3d(icnt)=interpdata(ie)%ilon(i)+(interpdata(ie)%ilat(i)-1)*nlon8+(k-1)*nlat8*nlon8
          end do
       end do
    end do
    call getiodof(3, (/nlon,nlat,nlev/), iorank, start3d(1:3), count3d(1:3))
    call nf_init_decomp(ncdf, (/londim,latdim,levdim/), ldof3d, start3d(1:3),count3d(1:3))
    deallocate(ldof3d)
    endif

    if (ilevdim>0) then
    allocate(ldof3d(lcount*(nlev+1)))
    icnt=0
    do k=1,nlev+1
       do ie=1,nelemd
          do i=1,interpdata(ie)%n_interp
             icnt=icnt+1
             ldof3d(icnt)=interpdata(ie)%ilon(i)+(interpdata(ie)%ilat(i)-1)*nlon8+(k-1)*nlat8*nlon8
          end do
       end do
    end do
    call getiodof(3, (/nlon,nlat,nlev+1/), iorank, start3d(1:3), count3d(1:3))
    call nf_init_decomp(ncdf, (/londim,latdim,ilevdim/), ldof3d(1:icnt), start3d(1:3),count3d(1:3))
    deallocate(ldof3d)
    endif

  end subroutine create_output_decomps
  logical function is_ncoldim (name)
    character*(*) name

    integer n

    if (.not. init_done) call init_dims
    is_ncoldim = .false.

    do n=1,nncolnames
       if (name == ncolnames(n)) then
          is_ncoldim = .true.
       end if
    end do

    return
  end function is_ncoldim
  logical function is_ncolfile(dims)
    type(dim_t), intent(in) :: dims(:)
    integer :: i

    is_ncolfile=.false.
    do i=1,size(dims)
       is_ncolfile=is_ncoldim(dims(i)%name)
       if(is_ncolfile) exit
    end do
  end function is_ncolfile
  !-------------------------------------------------------------------------------



  !-------------------------------------------------------------------------------
  subroutine init_dims

    integer n

    ewnames(:) = ' '
    nsnames(:) = ' '
    znames(:)  = ' '
    ncolnames(:)  = ' '

    ewnames(1) = 'lon'
    ewnames(2) = 'slon'
    nsnames(1) = 'lat'
    nsnames(2) = 'slat'
    znames(1)  = 'lev'
    znames(2)  = 'ilev'
    ncolnames(1) = ncoldimname

    do n=1,size (ewnames)
       if (ewnames(n) == ' ') exit
    end do
    newnames = n - 1

    do n=1,size (nsnames)
       if (nsnames(n) == ' ') exit
    end do
    nnsnames = n - 1

    do n=1,size (znames)
       if (znames(n) == ' ') exit
    end do
    nznames = n - 1

    do n=1,size (ncolnames)
       if (ncolnames(n) == ' ') exit
    end do
    nncolnames = n - 1

    ewdim(:)%name = ' '
    nsdim(:)%name = ' '
    zdim(:)%name = ' '
    ncoldim(:)%name = ' '

    ewdimi(:)%name = ' '
    nsdimi(:)%name = ' '
    zdimi(:)%name = ' '

    init_done = .true.

    return
  end subroutine init_dims
  !-------------------------------------------------------------------------------
  logical function is_ewdim (name)
    character*(*) name

    integer n

    if (.not. init_done) call init_dims
    is_ewdim = .false.

    do n=1,newnames
       if (name == ewnames(n)) then
          is_ewdim = .true.
       end if
    end do

    return
  end function is_ewdim


  !-------------------------------------------------------------------------------

  logical function is_nsdim (name)
    character*(*) name

    integer n

    if (.not. init_done) call init_dims
    is_nsdim = .false.

    do n=1,nnsnames
       if (name == nsnames(n)) then
          is_nsdim = .true.
       end if
    end do

    return
  end function is_nsdim
  !-------------------------------------------------------------------------------
  logical function is_zdim (name)
    character*(*) name

    integer n

    if (.not. init_done) call init_dims
    is_zdim = .false.

    do n=1,nznames
       if (name == znames(n)) then
          is_zdim = .true.
       end if
    end do

    return
  end function is_zdim
  !-------------------------------------------------------------------------------

  integer function get_dimlen (arr, dimname)
    type(dim_t) :: arr(:)
    character(len=*) dimname

    integer n

    do n=1,size (arr,1)
       if (trim(arr(n)%name) == trim(dimname)) then
          get_dimlen = arr(n)%len
          return
       end if
    end do

    write(6,*) 'WARNING: get_dimlen: dimname ',trim(dimname), ' not found'
    get_dimlen=-1

  end function get_dimlen

  subroutine copy_attributes(infile, outfile)
    use common_io_mod, only : nf_handle

    type(file_t) :: infile
    type(nf_handle) :: outfile

    integer :: i, j, ierr, natts, vid

    character(len=80) :: name

! copy global attributes
!	print *,'copying global  attributes'
    do i=1,infile%natts
!	print *,i,infile%natts
       ierr = pio_inq_attname(infile%FileID, PIO_GLOBAL, i, name)
!	print *,'name=',trim(name)
       ierr = PIO_copy_att(infile%FileID, PIO_GLOBAL, name, outfile%FileID, PIO_GLOBAL)
!	print *,'copied=',trim(name)
    end do
       
! copy variable attributes
    do j=1,size(infile%vars%ndims)
       natts=0
       vid = infile%vars%vardesc(j)%varid
       ierr = pio_inq_varnatts(infile%FileID,vid ,natts)
       do i=1,natts
          ierr = pio_inq_attname(infile%FileID, vid, i, name)
          ierr = PIO_copy_att(infile%FileID, vid, name, outfile%FileID, j)
       end do
    end do

  end subroutine copy_attributes

  subroutine getCompdof(Compdof, elem, lev)
    use dimensions_mod, only : nelemd, Gcols=>GlobalUniqueCols
    use element_mod, only : element_t
    use kinds, only: long_kind
    
    type(element_t), intent(in) :: elem(:)
    integer*8, pointer :: Compdof(:)
    integer, intent(in) :: lev

    integer :: nxyp, icnt, i, ie, k

    nxyp=0
    do ie=1,nelemd
      nxyp=nxyp+elem(ie)%idxp%NumUniquePts
    enddo
    allocate(Compdof(nxyp*lev))
    
    icnt=0
    do k=1,lev
       do ie=1,nelemd
          do i=1,elem(ie)%idxp%NumUniquePts
             icnt=icnt+1
             ! force calculation to be done integer*8
             compDOF(icnt)=elem(ie)%idxp%UniquePtOffset+i-1+(k-1)*int(GCols,kind=long_kind)
             if (compDOF(icnt)<0) then
                print *,lev,nxyp,GCols
                call abortmp('getCompdof(): compDOF intetger overflow')
             endif
          end do
       end do
    end do
	
  end subroutine getcompdof
#endif

  ! ----------------------------------------------------------------------------
  ! GLL-physgrid topography file utilites

  subroutine read_gll_topo_file(filename, elem, par, fields, fieldnames)
    ! fields(:np,:np,:nelemd,i) is field i in the list
    !     PHIS, SGH, SGH30, LANDM_COSLAT, LANDFRAC

    use element_mod, only: element_t
    use parallel_mod, only: parallel_t
    use kinds, only: real_kind
    use dimensions_mod, only: nelemd, nlev, np, npsq
    use common_io_mod, only: varname_len
#ifndef HOMME_WITHOUT_PIOLIBRARY
    use dof_mod, only: putuniquepoints
    use edge_mod, only: edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
    use edgetype_mod, only: edgebuffer_t
    use bndry_mod, only: bndry_exchangeV
#endif

    character(len=*), intent(in) :: filename
    type(element_t), intent(inout) :: elem(:)
    type(parallel_t),intent(in) :: par
    real(kind=real_kind), intent(out), dimension(np,np,nelemd,5) :: fields
    character(len=varname_len), intent(out) :: fieldnames(5)

#ifndef HOMME_WITHOUT_PIOLIBRARY
    type(file_t) :: infile
    type(edgeBuffer_t) :: edge    
    real(kind=real_kind), allocatable :: farray(:)
    real(kind=real_kind) :: ftmp(npsq)
    real(kind=real_kind), pointer :: arr3(:,:,:)

    integer :: ii,ie,ilev,iv,ierr,offset,vari,ncnt_in,nlyr

    fieldnames(1) = 'PHIS'
    fieldnames(2) = 'SGH'
    fieldnames(3) = 'SGH30'
    fieldnames(4) = 'LANDM_COSLAT'
    fieldnames(5) = 'LANDFRAC'

    ilev = 1
    nlyr = ilev*size(fieldnames)

    call initedgebuffer(par, edge, elem, nlyr)
    ncnt_in = sum(elem(1:nelemd)%idxp%numUniquePts)

    call infile_initialize(elem, par, trim(filename), fieldnames, infile)

    allocate(farray(ncnt_in))
    do vari = 1,nlyr
       farray = 1.0e-37
       call pio_read_darray(infile%FileID, infile%vars%vardesc(vari), iodesc2d, farray, ierr)

       offset = 0
       do ie = 1,nelemd
          do ii = 1,elem(ie)%idxP%NumUniquePts
             iv = offset + ii
             ftmp(ii) = farray(iv)
          end do
          offset = offset+elem(ie)%idxP%NumUniquePts
          fields(:,:,ie,vari) = 0
          call putUniquePoints(elem(ie)%idxP, ftmp(:elem(ie)%idxP%NumUniquePts), fields(:,:,ie,vari))
          call edgevpack(edge, fields(:,:,ie,vari), 1, vari-1, ie)
       end do
    end do

    call bndry_exchangeV(par, edge)

    do ie = 1,nelemd
       do vari = 1,nlyr
          call edgeVunpack(edge, fields(:,:,ie,vari), 1, vari-1, ie)
       end do
    end do

    deallocate(farray)
    call freeedgebuffer(edge)

    call pio_closefile(infile%fileid)
    call free_infile(infile)
#endif
  end subroutine read_gll_topo_file

  subroutine read_physgrid_topo_file(infilename, elem, par, fieldnames, nphys, pg_fields, stat)
    use element_mod, only: element_t
    use parallel_mod, only: parallel_t
    use kinds, only: real_kind
    use common_io_mod, only: varname_len, io_stride, num_io_procs, num_agg
#ifndef HOMME_WITHOUT_PIOLIBRARY
    use dimensions_mod, only: nelemd, nlev, np, npsq, nelem
    use control_mod, only: max_string_len
    use pio, only: pio_init, pio_openfile, pio_rearr_box, pio_inquire, pio_inq_dimname, &
         pio_inq_dimlen, pio_initdecomp
#endif

    character(len=*), intent(in) :: infilename
    type(element_t), intent(in) :: elem(:)
    type(parallel_t), intent(in) :: par
    character(len=varname_len), intent(in) :: fieldnames(:)
    real(kind=real_kind), intent(out) :: pg_fields(:,:,:)
    integer, intent(out) :: nphys, stat

#ifndef HOMME_WITHOUT_PIOLIBRARY
    type(file_desc_t) :: fileid
    type(var_desc_t) :: vardesc
    integer :: ndims, ncol, nvars, natts, nfield, fldi, iotype, nf2, ie
    character(len=pio_max_name) :: dimname
    integer*8, pointer :: dof(:)
    real(real_kind), allocatable :: raw(:)

    iotype = get_iotype()
    stat = pio_openfile(piofs, fileid, iotype, infilename)
    stat = pio_inquire(fileid, ndimensions=ndims, nvariables=nvars)
    if (ndims /= 1) then
       if (par%masterproc) print *, 'read_physgrid_topo expects input file to have 1 dim'
       stat = -1; return
    end if

    stat = pio_inq_dimname(fileid, 1, dimname)
    if (dimname /= 'ncol') then
       if (par%masterproc) print *, 'read_physgrid_topo expects dimname "ncol"'
       stat = -1; return
    end if
    stat = pio_inq_dimlen(fileid, 1, ncol)

    nphys = nint(sqrt(real(ncol/nelem, real_kind)))
    if (nphys*nphys*nelem /= ncol) then
       if (par%masterproc) then
          print *, 'read_physgrid_topo has inconsistent nelem, ncol, nphys:', &
               nelem, ncol, nphys
       end if
       stat = -1; return
    end if
    if (nphys > np) then
       if (par%masterproc) print *, 'read_physgrid_topo has nphys > np:', nphys, np
       stat = -1; return
    end if
    nf2 = nphys*nphys

    call make_physgrid_dof(elem, nphys, dof)
    call pio_initdecomp(piofs, pio_double, (/ncol/), dof, iodesc2d)
    deallocate(dof)

    allocate(raw(nf2*nelemd))
    nfield = size(fieldnames)
    do fldi = 1,nfield
       stat = pio_inq_varid(fileid, fieldnames(fldi), vardesc)
       if (stat /= 0) then
          if (par%masterproc) print *, 'read_physgrid_topo: could not find var', fieldnames(fldi)
          return
       end if
       call pio_read_darray(fileid, vardesc, iodesc2d, raw, stat)
       if (stat /= 0) then
          if (par%masterproc) print *, 'read_physgrid_topo: pio_read_darray returned stat', stat
          return
       end if
       do ie = 1,nelemd
          pg_fields(:nf2,ie,fldi) = raw(nf2*(ie-1)+1 : nf2*ie)
       end do
    end do
    deallocate(raw)

    call pio_closefile(fileid)
    stat = 0
#endif
  end subroutine read_physgrid_topo_file

  subroutine write_physgrid_topo_file(infilename, outfilenameprefix, elem, par, &
       gll_fields, pg_fields, latlon, fieldnames, nphys, history)
    ! gll_fields and fieldnames are as output from pio_read_gll_topo_file.

    use element_mod, only: element_t
    use parallel_mod, only: parallel_t
    use kinds, only: real_kind
    use dimensions_mod, only: nelemd, nlev, np, npsq, nelem
    use common_io_mod, only: varname_len
#ifndef HOMME_WITHOUT_PIOLIBRARY
    use pio_io_mod, only: nf_output_init_complete, nf_output_register_variables, nf_put_var_pio
    use control_mod, only: max_string_len
#endif

    integer, parameter :: nvar = 8, nvar_old = 5

    character(len=*), intent(in) :: infilename, outfilenameprefix, history
    type(element_t), intent(in) :: elem(:)
    type(parallel_t), intent(in) :: par
    real(kind=real_kind), intent(in) :: &
         gll_fields(np, np,      nelemd, nvar-1), &
         pg_fields (nphys*nphys, nelemd, nvar-1), &
         latlon    (nphys*nphys, nelemd, 2) ! (:,:,1) is lat, (:,:,2) is lon
    character(len=varname_len), intent(in) :: fieldnames(nvar_old)
    integer, intent(in) :: nphys

#ifndef HOMME_WITHOUT_PIOLIBRARY
    character(len=varname_len) :: varnames(nvar), name
    character(len=max_string_len) :: save_state(2)
    integer :: nf2, i, j, k, n, vardims(1,nvar), vartypes(nvar)
    integer(kind=nfsizekind) :: unused(1)
    logical :: varreqs(nvar)
    type(file_t) :: infile

    nf2 = nphys*nphys
    call set_output_vars(save_state)
    call physgrid_topo_begin_write(elem, par, outfilenameprefix, nphys)

    ! variables
    do i = 1,nvar_old
       varnames(i) = fieldnames(i)
       vardims(1,i) = 1
    end do
    varnames(nvar_old+1) = 'lat'
    vardims(1,nvar_old+1) = 1
    varnames(nvar_old+2) = 'lon'
    vardims(1,nvar_old+2) = 1
    varnames(nvar) = 'PHIS_d'
    vardims(1,nvar) = 2
    varreqs = .true.
    vartypes = pio_double
    call nf_output_register_variables(ncdf, nvar, varnames, vardims, vartypes, varreqs)

    ! Copy variable and global attributes from the GLL source topography file.
    call infile_initialize(elem, par, trim(infilename), varnames(1:nvar-1), infile)
    call copy_attributes(infile, ncdf(1))
    ! Copy PHIS attributes to PHIS_d.
    k = infile%vars%vardesc(1)%varid ! PHIS id
    j = pio_inq_varnatts(infile%fileid, k, n) ! n atts
    do i = 1,n
       j = pio_inq_attname(infile%fileid, k, i, name) ! att name
       j = pio_copy_att(infile%fileid, k, name, ncdf(1)%fileid, nvar) ! PHIS_d has id nvar
    end do
    call pio_closefile(infile%fileid)
    call free_infile(infile)
    j = pio_put_att(ncdf(1)%fileid, pio_global, 'history', history)
    
    call nf_output_init_complete(ncdf)

    ! Write physgrid topo fields.
    do i = 1,nvar_old
       call nf_put_var_pio(ncdf(1), reshape(pg_fields(:,:,i), (/nf2*nelemd/)), &
            unused, unused, ncdf(1)%varlist(i))
    end do
    ! Write lat-lon.
    do i = nvar_old+1,nvar-1
       call nf_put_var_pio(ncdf(1), reshape(latlon(:,:,i-nvar_old), (/nf2*nelemd/)), &
            unused, unused, ncdf(1)%varlist(i))
    end do
    ! Write GLL field PHIS_d.
    call write_gll_field(elem, ncdf(1), gll_fields(:,:,:,1), ncdf(1)%varlist(nvar))

    call pio_closefile(ncdf(1)%fileid)
    call restore_output_vars(save_state)
#endif
  end subroutine write_physgrid_topo_file

  subroutine write_physgrid_smoothed_phis_file(outfilenameprefix, elem, par, &
       gll_fields, pg_fields, nphys, history, output_latlon)
    use element_mod, only: element_t
    use parallel_mod, only: parallel_t
    use kinds, only: real_kind
    use dimensions_mod, only: nelemd, nlev, np, npsq, nelem
#ifndef HOMME_WITHOUT_PIOLIBRARY
    use common_io_mod, only: varname_len
    use pio_io_mod, only: nf_output_init_complete, nf_output_register_variables, nf_put_var_pio
    use control_mod, only: max_string_len
#endif

    integer, parameter :: max_nvar = 6

    character(len=*), intent(in) :: outfilenameprefix, history
    type(element_t), intent(in) :: elem(:)
    type(parallel_t), intent(in) :: par
    real(kind=real_kind), intent(in) :: gll_fields(:,:,:,:), pg_fields(:,:,:)
    integer, intent(in) :: nphys
    logical, optional, intent(in) :: output_latlon

#ifndef HOMME_WITHOUT_PIOLIBRARY
    character(len=varname_len) :: varnames(max_nvar), name
    character(len=max_string_len) :: save_state(2)
    integer :: i, nf2, stat, nvar, vardims(1,max_nvar), vartypes(max_nvar)
    integer(kind=nfsizekind) :: unused(1)
    logical :: varreqs(max_nvar), outll

    outll = .false.
    if (present(output_latlon)) outll = output_latlon
    if (outll .and. (size(gll_fields,4) < 3 .or. size(pg_fields,3) < 3)) then
       if (par%masterproc) then
          print *, 'write_physgrid_smoothed_phis_file: gll_fields and pg_fields must&
               & have 3 fields if output_latlon=true; setting output_latlon=false'
       end if
       outll = .false.
    end if

    nvar = 2
    if (outll) nvar = 6

    nf2 = nphys*nphys
    call set_output_vars(save_state)
    call physgrid_topo_begin_write(elem, par, outfilenameprefix, nphys)

    ! variables
    varnames(1) = 'PHIS'  ; vardims(1,1) = 1
    varnames(2) = 'PHIS_d'; vardims(1,2) = 2
    if (outll) then
       varnames(3) = 'lat'  ; vardims(1,3) = 1
       varnames(4) = 'lat_d'; vardims(1,4) = 2
       varnames(5) = 'lon'  ; vardims(1,5) = 1
       varnames(6) = 'lon_d'; vardims(1,6) = 2
    end if
    varreqs = .true.
    vartypes = pio_double
    call nf_output_register_variables(ncdf, nvar, varnames, vardims, vartypes, varreqs)
    stat = pio_put_att(ncdf(1)%fileid, pio_global, 'history', history)
    call nf_output_init_complete(ncdf)

    call nf_put_var_pio(ncdf(1), reshape(pg_fields(:nf2,:,1), (/nf2*nelemd/)), &
         unused, unused, ncdf(1)%varlist(1))
    call write_gll_field(elem, ncdf(1), gll_fields(:,:,:,1), ncdf(1)%varlist(2))
    if (outll) then
       do i = 1,2
          call nf_put_var_pio(ncdf(1), reshape(pg_fields(:nf2,:,i+1), (/nf2*nelemd/)), &
               unused, unused, ncdf(1)%varlist(1 + 2*i))
          call write_gll_field(elem, ncdf(1), gll_fields(:,:,:,i+1), ncdf(1)%varlist(2 + 2*i))
       end do
    end if

    call pio_closefile(ncdf(1)%fileid)
    call restore_output_vars(save_state)
#endif
  end subroutine write_physgrid_smoothed_phis_file

  ! ------------------------------------
  ! Utils

#ifndef HOMME_WITHOUT_PIOLIBRARY
  function get_iotype() result(iotype)
    use pio, only: iotype_netcdf, iotype_pnetcdf

    integer :: iotype

    if (output_type == 'netcdf') then
       iotype = iotype_netcdf
    else
       iotype = iotype_pnetcdf
    end if
  end function get_iotype

  subroutine make_physgrid_dof(elem, nphys, dof)
    ! Caller must deallocate dof when done.

    use element_mod, only: element_t
    use dimensions_mod, only: nelemd

    type(element_t), intent(in) :: elem(:)
    integer, intent(in) :: nphys
    integer*8, intent(out), pointer :: dof(:)
    
    integer :: nf2, ie, j

    nf2 = nphys*nphys
    allocate(dof(nelemd*nf2))
    do ie = 1,nelemd
       do j = 1,nf2
          dof(nf2*(ie-1) + j) = nf2*(elem(ie)%globalid-1) + j
       end do
    end do
  end subroutine make_physgrid_dof

  subroutine set_output_vars(save_state)
    use common_io_mod, only: output_frequency, output_start_time, output_end_time, &
         output_dir, output_prefix
    use control_mod, only: max_string_len
    
    character(len=max_string_len), intent(out) :: save_state(2)

    ! Save global output-file state.
    save_state(1) = output_dir
    save_state(2) = output_prefix

    ! We get a spurious '1' at the end of the output filename, but other than
    ! that, we get exactly what we want.
    output_prefix = ''
    output_dir = ''
    output_frequency(2:max_output_streams) = 0
    output_frequency(1) = 1
    output_start_time(1) = 0
    output_end_time(1) = 1
  end subroutine set_output_vars

  subroutine restore_output_vars(save_state)
    use common_io_mod, only: output_frequency, output_start_time, output_end_time, &
         output_dir, output_prefix
    use control_mod, only: max_string_len
    
    character(len=max_string_len), intent(in) :: save_state(2)

    output_dir = save_state(1)
    output_prefix = save_state(2)
  end subroutine restore_output_vars

  subroutine physgrid_topo_begin_write(elem, par, outfilenameprefix, nphys)
    use element_mod, only: element_t
    use parallel_mod, only: parallel_t
    use dimensions_mod, only: np, nelem
    use pio_io_mod, only: nf_output_init_begin, nf_output_register_dims, nf_init_decomp
    use common_io_mod, only: varname_len

    integer, parameter :: ndim = 2

    type(element_t), intent(in) :: elem(:)
    type(parallel_t), intent(in) :: par
    character(len=*), intent(in) :: outfilenameprefix
    integer, intent(in) :: nphys

    character(len=varname_len) :: dimnames(ndim)
    integer :: dimsizes(ndim), nf2
    integer*8 :: itmp(1)
    integer(kind=nfsizekind) :: unused(1)
    integer*8, pointer :: dof(:)

    call nf_output_init_begin(ncdf, par%masterproc, par%nprocs, par%rank, par%comm, &
         outfilenameprefix, 0)

    ! dimensions
    nf2 = nphys*nphys
    dimnames(1) = 'ncol'
    dimsizes(1) = nelem*nf2
    dimnames(2) = 'ncol_d'
    ! Euler's formula applied to the GLL grid:
    !   v - e + f = 2 => v = e - f + 2 = 2 f - f + 2 = (np-1)^2 nelem + 2
    dimsizes(2) = (np-1)**2*nelem + 2
    call nf_output_register_dims(ncdf, ndim, dimnames, dimsizes)

    ! physgrid decomp
    call make_physgrid_dof(elem, nphys, dof)
    call nf_init_decomp(ncdf, (/1/), dof, &
         unused, unused) ! these args are unused
    deallocate(dof)
    ! GLL decomp
    call getcompdof(dof, elem, 1)
    call nf_init_decomp(ncdf, (/2/), dof, unused, unused)
    deallocate(dof)
  end subroutine physgrid_topo_begin_write
  
  subroutine write_gll_field(elem, ncdf, phis, nfvar)
    use element_mod, only: element_t
    use common_io_mod, only: nf_variable
    use dof_mod, only: UniquePoints
    use dimensions_mod, only: nelemd
    use kinds, only: real_kind
    use pio_io_mod, only: nf_put_var_pio

    type(element_t), intent(in) :: elem(:)
    type(nf_handle), intent(inout) :: ncdf
    real(real_kind), intent(in) :: phis(:,:,:)
    type(nf_variable), intent(in) :: nfvar

    integer(kind=nfsizekind) :: unused(1)
    real(kind=real_kind), allocatable :: gll_unique(:)
    integer :: k, i

    allocate(gll_unique(sum(elem%idxp%NumUniquePts)))
    k = 1
    do i = 1,nelemd
       call UniquePoints(elem(i)%idxP, phis(:,:,i), &
            gll_unique(k : k + elem(i)%idxp%NumUniquePts - 1))
       k = k + elem(i)%idxp%NumUniquePts
    end do
    call nf_put_var_pio(ncdf, gll_unique, unused, unused, nfvar)
    deallocate(gll_unique)
  end subroutine write_gll_field
#endif

end module interpolate_driver_mod
