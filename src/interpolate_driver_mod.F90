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
#ifdef PIO_INTERP
  use pio, only : file_desc_t, var_desc_t , io_desc_t, & ! _EXTERNAL
   pio_get_att, pio_setdebuglevel, pio_closefile, &
   pio_put_att, pio_global, pio_put_var, &
   pio_read_darray, pio_setframe, pio_get_var, &
   pio_Offset, pio_char, &
   pio_inq_varid, pio_inq_attname, pio_copy_att, pio_inq_varnatts,&
   PIO_MAX_NAME, pio_double

  use pio_nf_utils, only : copy_pio_var ! _EXTERNAL
  use pio_io_mod, only : nfsizekind ! _EXTERNAL
  use common_io_mod, only : nf_handle, max_output_streams, MAX_INFILES,&
              infilenames,piofs, output_type
  use parallel_mod, only : abortmp, syncmp, haltmp
  implicit none
  private
!#include "pnetcdf.inc"

  public :: interpolate_driver

  integer :: nlat, nlon

  type dim_t
     integer :: len
     character(len=PIO_MAX_NAME) :: name
  end type dim_t

  type var_t
     integer, pointer :: ndims(:)
!     integer, pointer :: type(:)
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

  logical init_done
  data init_done /.false./

  character*(pio_max_name) ewnames (maxdims)
  character*(pio_max_name) nsnames (maxdims)
  character*(pio_max_name) znames (maxdims)
  character*(pio_max_name) ncolnames (maxdims)

  integer :: newnames, nnsnames, nznames, nncolnames


  type(nf_handle), private, target, save :: ncdf(max_output_streams)
!  type(io_desc_t), pointer :: iodesc3d, iodesc2d, iodesc3dp1
  type(io_desc_t) , save:: iodesc3d, iodesc2d, iodesc3dp1


contains
  subroutine interpolate_driver(elem,hybrid)
    use dimensions_mod, only : ne, nelem, np, nlev
    use common_io_mod, only : varnames=>output_varnames1, nf_handle
    use element_mod, only : element_t
    !  use interpolate_mod
    use dof_mod
    use hybrid_mod, only : hybrid_t
    use interpolate_mod, only : interpdata_t
    !  use domain_mod, only : domain1d_t, decompose
    !  use thread_mod, only : omp_get_thread_num
    !  use reduction_mod, only : reductionbuffer_ordered_1d_t
    implicit none

    integer, parameter :: maxvars=30
    type(element_t) :: elem(:)
    type(file_t) :: infile
    !  type(element_t), allocatable :: elem(:)
    type(hybrid_t), intent(in) :: hybrid
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
          if (hybrid%par%masterproc) print *,'input file: ',trim(infilename)
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
  end subroutine interpolate_driver

  subroutine free_infile(infile)
    type(file_t), intent(inout) :: infile

    deallocate(infile%vars%vardesc)
    nullify(infile%vars%vardesc)
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
    use dimensions_mod, only : ne, np, nlev
    type(element_t) :: elem(:)
    type(parallel_t) :: par
    character(len=*), intent(in) :: infilename
    character(len=*), intent(inout) :: varnames(:)
    type(file_t), intent(out) :: infile
    type(io_desc_t), pointer :: iodesc

    integer :: ndims, nvars,  varcnt
    integer :: ret, i, ncid, n
    integer :: varid, ncols, lev
    integer :: ne_file, np_file, nlev_file

    integer, pointer :: ldof(:)
    integer(kind=PIO_Offset) :: start(3), count(3)
    integer :: iorank, dimcnt

    call PIO_Init(par%rank, par%comm, num_io_procs, num_agg, &
         io_stride, pio_rearr_box, PIOFS)

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

    do i=1,ndims
       ret = PIO_inq_dimname(Infile%fileid, i, infile%dims(i)%name)
       ret = PIO_inq_dimlen(Infile%fileid, i, infile%dims(i)%len)
       if(infile%dims(i)%name.eq.'ncol') ncols=infile%dims(i)%len
       if(infile%dims(i)%name.eq.'lev') nlev_file = infile%dims(i)%len
    end do
    if(trim(varnames(1)) == 'all') then
       varcnt=nvars
    else
       varcnt=0
       do i=1,size(varnames)
          if(varnames(i) == '') exit
	  varcnt=varcnt+1
       end do
    end if

    if(is_ncolfile(infile%dims)) then
       ret = PIO_get_att(infile%FileID, PIO_GLOBAL, 'np', np_file)
       ret = PIO_get_att(infile%FileID, PIO_GLOBAL, 'ne', ne_file)
       
       nlev_file=get_dimlen(infile%dims,"lev")
       
       if(ne_file/=ne) then
          print *,'ne, ne_file',ne,ne_file
          call abortmp('The variable ne in the namelist must be the same as that of the file.')
       end if
       if(nlev_file/=nlev) then
          print *,'nlev, nlev_file',nlev,nlev_file
          call abortmp('The variable nlev in Params.inc must be the same as that of the file, you will need to recompile.')
       end if
       if(np_file/=np) then
          print *,'np, np_file',np,np_file
          call abortmp('The variable np in Params.inc must be the same as that of the file, you will need to recompile.')
       end if
    else
       call abortmp('The input file is missing required dimensions')
    end if

! define the decompositions...
    iorank=piofs%io_rank
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
            infile%vars%vardesc(i)%type)
       ret = PIO_inq_varndims(infile%FileID, infile%vars%vardesc(i), &
            infile%vars%ndims(i))
       ret = PIO_inq_vardimid(infile%FileID, infile%vars%vardesc(i), &
            infile%vars%dimids(1:infile%vars%ndims(i),i))
       lev=1
       do n=1,infile%vars%ndims(i)
          if(infile%dims(infile%vars%dimids(n,i))%name.eq.'lev') lev=nlev
          if(infile%dims(infile%vars%dimids(n,i))%name.eq.'ilev') lev=nlev+1
          if(infile%dims(infile%vars%dimids(n,i))%name.eq.'time') infile%vars%timedependent(i)=.true.
          if(infile%dims(infile%vars%dimids(n,i))%name.eq.'ncol') infile%vars%decomposed(i)=.true.

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
    logical lonregistered, latregistered
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

    lonregistered=.false.
    latregistered=.false.
    vardims=0
    do i=1,nvars
       k=1
       otype(i)=infile%vars%vardesc(i)%type
       do j=1,infile%vars%ndims(i)
          if(infile%vars%dimids(j,i).eq.ncoldimid) then
             if(infile%vars%name(i).eq.'lon') then
                vardims(k,i) = londimid
                lonregistered=.true.
             else if(infile%vars%name(i).eq.'lat') then
                vardims(k,i) = latdimid
                latregistered=.true.
             else
                vardims(k,i) = londimid
                vardims(k+1,i) = latdimid
                k=k+1
             end if
             k=k+1
          else
             vardims(k,i)=infile%vars%dimids(j,i)
             k=k+1
          end if
       end do
    end do

    ndnt = index(infilename,".nc")

    call nf_output_init_begin(ncdf,hybrid%par%masterproc,hybrid%par%nprocs,hybrid%par%rank, &
         hybrid%par%comm, infilename(1:ndnt)//'interp'   ,0)
    
    call nf_output_register_dims(ncdf, ofdims, dimnames, dimsize)
    call create_output_decomps(ncdf, interpdata, nlon, nlat)
    if (hybrid%par%masterproc) print *,'registering output variables -- if you get a seg fault, try increasing the number of processors used to help memory management'
    call nf_output_register_variables(ncdf, nvars, infile%vars%name, &
         vardims, otype, varrequired)

    vardims(1,1)=latdimid
    call nf_output_register_variables(ncdf, 1, (/'gw'/), vardims(1:1,:), &
         (/pio_double/), (/.true./))
    ! if the user did not specify lat,lon, be sure to add them:
    if (.not. latregistered) then
       call nf_output_register_variables(ncdf, 1, (/'lat'/), vardims(1:1,:), &
            (/pio_double/), (/.true./))
    endif
    if (.not. lonregistered) then
       vardims(1,1)=londimid
       call nf_output_register_variables(ncdf, 1, (/'lon'/), vardims(1:1,:), &
            (/pio_double/), (/.true./))
    endif


    deallocate(dimnames, dimsize, varrequired, vardims, otype)

    outfile=>ncdf(1)
    iorank=piofs%io_rank
    
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
    use pio_io_mod, only : nf_handle, nf_put_var ! _EXTERNAL
    use interpolate_mod, only : interpdata_t, interpolate_scalar, interpolate_vector, &
         get_interp_parameter,var_is_vector_vvar,var_is_vector_uvar
    use element_mod, only : element_t
    use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
    use dimensions_mod, only : nelemd, nlev, np
    use parallel_mod, only : parallel_t, syncmp
    use bndry_mod, only : bndry_exchangeV
    use common_io_mod, only : readystate, nf_variable

    type(element_t), intent(inout) :: elem(:)
    type(parallel_t),intent(in) :: par
    type(file_t), intent(inout) :: infile
    type(nf_handle), intent(inout) :: outfile
    type(interpdata_t), intent(in) :: interpdata(:)

    type(edgeBuffer_t) :: edge    
    real(kind=real_kind), pointer :: array(:,:), varray(:,:,:)
    real(kind=real_kind), pointer :: farray(:)	, fvarray(:,:), ftmp(:,:), fvtmp(:,:,:)

    integer(kind=nfsizekind) :: start(3), count(3)
    integer :: starti(3), counti(3)
    integer :: len, ii, iu, iv, k, it, nd, offset, iuvar
    integer :: nvars, i, n, ncnt_out, ncnt_in, st, en, ncols, ierr, ie, lev, ntimes, ndims
    character*(PIO_MAX_NAME) :: tmpname
    integer(kind=nfsizekind) :: start2d(3), count2d(3), start3d(4), count3d(4)

    if(outfile%state .ne. readystate) then
       print *,__FILE__,__LINE__,outfile%state
       call abortmp('outfile not in readystate')
    end if

    nvars = size(infile%vars%name)

    ncnt_in = sum(elem(1:nelemd)%idxp%numUniquePts)
    ncnt_out = sum(interpdata(1:nelemd)%n_interp)

!    call initedgebuffer(edge,2*nlev)
    call initedgebuffer(edge,4*nlev)


    VARLOOP: do i=nvars,1,-1
       if(trim(infile%vars%name(i)).eq.'lat' .or. trim(infile%vars%name(i)).eq.'lon') then
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
                      call pio_setframe(infile%vars%vardesc(i),int(n,kind=PIO_Offset))
                      call pio_setframe(infile%vars%vardesc(iv),int(n,kind=PIO_Offset))
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
                      call edgevpack(edge, elem(ie)%state%v(:,:,:,:,1),2*lev,0,elem(ie)%desc)
                   end do
                   call bndry_exchangeV(par, edge)
                   st = 1
                   deallocate(fvarray)
                   allocate(varray(ncnt_out,lev,2))
                   do ie=1,nelemd
                      en=st+interpdata(ie)%n_interp-1
                      call edgeVunpack(edge, elem(ie)%state%v(:,:,:,:,1),2*lev,0,elem(ie)%desc)
                      call interpolate_vector(interpdata(ie), elem(ie), &
                           elem(ie)%state%V(:,:,:,1:lev,1),np, lev, varray(st:en,:,:), 0)
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
! Ignore these
!             else if(scan(infile%vars%name(i),'V').eq.1) then  ! ignores everything that starts with V
!             else if(infile%vars%name(i).eq.'lon') then    ! already ignored above
!             else if(infile%vars%name(i).eq.'lat') then    ! already ignored above
! 
             else 
                if(par%masterproc) print *,'Interpolating ',trim(infile%vars%name(i))
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
                      call pio_setframe(infile%vars%vardesc(i),int(n,kind=PIO_Offset))
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
                   do ie=1,nelemd
                      allocate(ftmp(elem(ie)%idxP%NumUniquePts, lev))
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
                      call edgevpack(edge, elem(ie)%state%Q(:,:,1:lev,1),lev,0,elem(ie)%desc)
                      deallocate(ftmp)
                   end do

                   call bndry_exchangeV(par, edge)
                   st = 1
                   deallocate(farray)
                   allocate(array(ncnt_out,lev))
                   array=0
                   do ie=1,nelemd
                      en=st+interpdata(ie)%n_interp-1
                      call edgeVunpack(edge, elem(ie)%state%Q(:,:,1:lev,1),lev,0,elem(ie)%desc)
                      
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
             if( outfile%varlist(i)%vardesc%type .eq. PIO_Char) then
                do n=2,infile%vars%ndims(i)
                   len = len*infile%dims(infile%vars%dimids(n,i))%len
                end do
                call copy_pio_var(infile%FileID, Outfile%fileid, infile%vars%vardesc(i), &
                     outfile%varlist(i)%vardesc, infile%dims(infile%vars%dimids(1,i))%len, len)
             else
                if(infile%vars%ndims(i)<=1) then
                   do n=1,infile%vars%ndims(i)
                      len = len*infile%dims(infile%vars%dimids(n,i))%len
                   end do
                   call copy_pio_var(infile%FileID, Outfile%fileid, infile%vars%vardesc(i), &
                        outfile%varlist(i)%vardesc, len)
                else
                   do n=1,infile%vars%ndims(i)
                      counti(n) = infile%dims(infile%vars%dimids(n,i))%len
                   end do
                   call copy_pio_var(infile%FileID, Outfile%fileid, infile%vars%vardesc(i), &
                        outfile%varlist(i)%vardesc, counti(1:infile%vars%ndims(i)))
                end if
             end if
          end if
       end if
       call syncmp(par)
    end do VARLOOP

    call freeedgebuffer(edge)
    
  end subroutine interpolate_vars
!
! Create the pio decomps for the output file.
!
  subroutine create_output_decomps(ncdf,interpdata,nlon,nlat)
    use pio_io_mod ! _EXTERNAL
    use dimensions_mod
    use interpolate_mod
    use interp_movie_mod, only : getiodof
    implicit none
    type(nf_handle) :: ncdf(:)
    type(interpdata_t):: interpdata(:)
    integer, intent(in) :: nlon, nlat
    integer, pointer :: ldof2d(:),ldof3d(:), iodof2d(:), iodof3d(:)
    integer, pointer :: latdof(:), londof(:)

    integer :: icnt, ie, i, lcount, tdof(1), tiodof(1)
    integer :: k, iorank
    integer(kind=nfsizekind) :: start1d(1), count1d(1)
    integer :: londim, latdim, timedim, levdim, ilevdim
    integer(kind=nfsizekind) :: start2d(3), count2d(3), start3d(4), count3d(4)

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
    iorank = piofs%io_rank

    ! Create the DOF arrays
    allocate(ldof2d(lcount))
    icnt=0
    do ie=1,nelemd
       do i=1,interpdata(ie)%n_interp
          icnt=icnt+1
          ldof2d(icnt)=interpdata(ie)%ilon(i)+(interpdata(ie)%ilat(i)-1)*nlon
       end do
    end do
    call getiodof(2, (/nlon,nlat/), iorank, iodof2d, start2d(1:2), count2d(1:2))
    call nf_init_decomp(ncdf, (/londim,latdim/), ldof2d, iodof2d,start2d(1:2),count2d(1:2))
    deallocate(iodof2d, ldof2d)

    if (levdim>0) then
    allocate(ldof3d(lcount*nlev))
    icnt=0
    do k=1,nlev
       do ie=1,nelemd
          do i=1,interpdata(ie)%n_interp
             icnt=icnt+1
             ldof3d(icnt)=interpdata(ie)%ilon(i)+(interpdata(ie)%ilat(i)-1)*nlon+(k-1)*nlat*nlon
          end do
       end do
    end do
    call getiodof(3, (/nlon,nlat,nlev/), iorank, iodof3d, start3d(1:3), count3d(1:3))
    call nf_init_decomp(ncdf, (/londim,latdim,levdim/), ldof3d, iodof3d,start3d(1:3),count3d(1:3))
    deallocate(iodof3d, ldof3d)
    endif

    if (ilevdim>0) then
    allocate(ldof3d(lcount*(nlev+1)))
    icnt=0
    do k=1,nlev+1
       do ie=1,nelemd
          do i=1,interpdata(ie)%n_interp
             icnt=icnt+1
             ldof3d(icnt)=interpdata(ie)%ilon(i)+(interpdata(ie)%ilat(i)-1)*nlon+(k-1)*nlat*nlon
          end do
       end do
    end do
    call getiodof(3, (/nlon,nlat,nlev+1/), iorank, iodof3d, start3d(1:3), count3d(1:3))
    call nf_init_decomp(ncdf, (/londim,latdim,ilevdim/), ldof3d(1:icnt), iodof3d,start3d(1:3),count3d(1:3))
    deallocate(iodof3d, ldof3d)
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
    ncolnames(1) = 'ncol'

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

    type(element_t), intent(in) :: elem(:)
    integer, pointer :: Compdof(:)
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
             compDOF(icnt)=elem(ie)%idxp%UniquePtOffset+i-1+(k-1)*GCols
          end do
       end do
    end do
	
  end subroutine getcompdof

    
#endif

end module interpolate_driver_mod





