#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module pio_io_mod
#if defined PIO || defined PIO_INTERP
  use kinds, only : int_kind, real_kind
  use pio ! _EXTERNAL
  use pio_types ! _EXTERNAL
  use pio_kinds, only: nfsizekind=>PIO_OffSet ! _EXTERNAL
#ifdef _MPI
  !HOMME Specific:  abortmp, mpireal_t, mpiinteger_t
  use parallel_mod, only : abortmp,  mpi_info_null, mpireal_t,  mpiinteger_t, &
       mpi_integer, mpi_sum
#else
  !HOMME Specific: abortmp, mpireal_t,mpiinteger_t
  use parallel_mod, only : abortmp, mpireal_t, mpiinteger_t
!  use kinds, only : nfsizekind=>int_kind
#endif
  use common_io_mod, only : varname_len, output_start_time, output_end_time, &
       output_frequency, max_output_streams, output_varnames1, output_varnames2, &
       output_varnames3, output_varnames4, output_varnames5, max_output_streams, &
       num_io_procs, num_agg, io_stride, output_dir, output_prefix, max_output_variables, &
       nf_selectedvar, get_varindex, get_current_varnames, output_type, &
       nf_dim, nf_variable, nf_handle, beginstate, dimsstate, varsstate, readystate, nf_decomp, piofs
  use control_mod, only : MAX_STRING_LEN         !HOMME Specific: MAX_STRING_LEN
  implicit none  
!  private   
!commented out the above private line because gfortran has trouble:
!Error: The component 'decomplist' is a PRIVATE type and cannot be a component of 'nf_handle', which is PUBLIC at (1)

  interface nf_put_var
     module procedure nf_put0DI
     module procedure nf_put0DR
     module procedure nf_put1DI
     module procedure nf_put1DR
     module procedure nf_put2DR
     module procedure nf_put3DR
  end interface

  ! public interfaces  
  public :: nf_output_init_begin, nf_output_init_complete, nf_output_register_variables, &
       nf_put_var, nf_close, nf_close_all, &
       nf_output_register_dims, nf_advance_frame, nf_handle, nf_variable, get_variable_handle, &
       nf_variable_attributes, nf_get_frame, nfsizekind, nf_dim, nf_init_decomp
  public :: unlim_dim

!#include <pnetcdf.inc> ! _EXTERNAL

!  integer(kind=nfsizekind),parameter :: unlim_dim=nf_unlimited
  integer(kind=nfsizekind),parameter :: unlim_dim=PIO_unlimited

contains

  !
  !  Advances the time (unlimited dimension) level of the associated file
  !

  subroutine nf_advance_frame(ncdf)
    type(nf_handle), intent(inout) :: ncdf
    ncdf%iframe=ncdf%iframe+1
    !
    ! This call flushs the i/o buffers after a write, it may affect performance
    ! and can be removed if nessesary
    ! intentionally blank for pio
  end subroutine nf_advance_frame

  integer function nf_get_frame(ncdf)
    type(nf_handle), intent(in) :: ncdf
    nf_get_frame=ncdf%iframe
  end function nf_get_frame

  !
  !  These are all internal variations of nf_put_var which writes a timeslice of 
  !  a registered variable to the file associated with nf_handle ncdf
  !
  subroutine nf_put0DI(ncdf, var, start, count,varid, name)
    type(nf_handle), intent(in) :: ncdf
    integer, intent(in) :: var
    integer(kind=nfsizekind), intent(in) :: start(:), count(:)
    type(nf_variable), intent(in), optional, target :: varid
    character*(*), intent(in),optional :: name
    integer :: vindex, extent, vartmp(1)
    type (nf_variable), pointer :: varptr
    !$OMP SINGLE

    if(ncdf%state /= readystate) then
       print *,__FILE__,__LINE__,ncdf%state, readystate       
       call abortmp('sanity check failed, bad file handle')
    end if

    if(present(varid)) then
       varptr=>varid
    else if(present(name)) then
       vindex = get_varindex(name, ncdf%varlist)    
       varptr => ncdf%varlist(vindex)
    else
       call abortmp('one of the optional arguments name and varid must be provided')   
    end if

    if(varptr%ndims>0) then
       !       varptr%lstart(1)=ncdf%iframe
       extent=count(1)
    else
       extent=1
    end if

    call abortmp('function not supported')

    !$OMP END SINGLE
  end subroutine nf_put0DI

  subroutine nf_put0DR(ncdf, var, start, count, varid, name)
    type(nf_handle), intent(inout) :: ncdf
    real(kind=real_kind), intent(in) :: var
    integer(kind=nfsizekind), intent(in) :: start(:), count(:)
    type(nf_variable), intent(in), optional, target :: varid
    character*(*), intent(in),optional :: name
    type (nf_variable), pointer :: varptr
    integer :: ierr, vindex

    !$OMP SINGLE
    if(ncdf%state /= readystate) then
       print *,__FILE__,__LINE__,ncdf%state, readystate 
       call abortmp('sanity check failed, bad file handle')
    end if
    if(present(varid)) then
       varptr=>varid
    else if(present(name)) then
       vindex = get_varindex(name, ncdf%varlist)    
       varptr => ncdf%varlist(vindex)
    else
       call abortmp('one of the optional arguments name and varid must be provided')   
    end if

    if(varptr%timedependent) then
       ierr = pio_put_var(ncdf%FileID, varptr%vardesc, int(start(1:1)), var)	
    else
       ierr = pio_put_var(ncdf%FileID, varptr%vardesc, var)	
    end if
    !$OMP END SINGLE


  end subroutine nf_put0DR

  subroutine nf_put1DI(ncdf, var, start, count, varid, name)
    type(nf_handle), intent(inout) :: ncdf
    integer, intent(in) :: var(:)
    integer(kind=nfsizekind), intent(in) :: start(:), count(:)
    type(nf_variable), intent(in), optional, target :: varid
    character*(*), intent(in),optional :: name
    integer :: vindex, extent, ierr
    type (nf_variable), pointer :: varptr
    type(io_desc_t), pointer :: iodesc
    integer :: id
    !$OMP SINGLE

    if(ncdf%state /= readystate) then
       print *,__FILE__,__LINE__,ncdf%state, readystate       
       call abortmp('sanity check failed, bad file handle')
    end if
    if(present(varid)) then
       varptr=>varid
    else if(present(name)) then
       vindex = get_varindex(name, ncdf%varlist)    
       varptr => ncdf%varlist(vindex)
    else
       call abortmp('one of the optional arguments name and varid must be provided')   
    end if
    if(varptr%timedependent) call PIO_SetFrame(varptr%vardesc,start(2))
    nullify(iodesc)
    do id=1,size(ncdf%decomplist)
       if(varptr%dimsid==ncdf%decomplist(id)%dimsid) then
          iodesc=>ncdf%decomplist(id)%iodesc
          exit
       end if
    end do
    if(.not. associated(iodesc)) then
       call abortmp('no iodesc found for variable')
    end if
    call PIO_write_darray(ncdf%FileID, varptr%vardesc, iodesc, var, ierr)

    !$OMP END SINGLE
  end subroutine nf_put1DI

  subroutine nf_put1DR(ncdf, var, start, count, varid, name)
    type(nf_handle), intent(inout) :: ncdf
    real(kind=real_kind), intent(in) :: var(:)
    integer(kind=nfsizekind), intent(in) :: start(:), count(:)
    type(nf_variable), intent(in), optional, target :: varid
    character*(*), intent(in),optional :: name
    integer :: vindex, extent, ierr, ndims
    type (nf_variable), pointer :: varptr
    type (io_desc_t), pointer :: iodesc
    integer :: id

    !$OMP SINGLE

    if(ncdf%state /= readystate) then
       print *,__FILE__,__LINE__,ncdf%state, readystate       
       call abortmp('sanity check failed, bad file handle')
    end if
    if(present(varid)) then
       varptr=>varid
    else if(present(name)) then
       vindex = get_varindex(name, ncdf%varlist)    
       varptr => ncdf%varlist(vindex)
    else
       call abortmp('one of the optional arguments name and varid must be provided')   
    end if
    ndims=size(start)	
    if(varptr%timedependent) then
       call PIO_SetFrame(varptr%vardesc,start(ndims))
    end if

    nullify(iodesc)
    do id=1,size(ncdf%decomplist)
       if(varptr%dimsid==ncdf%decomplist(id)%dimsid) then
          iodesc=>ncdf%decomplist(id)%iodesc
          exit
       end if
    end do
    if(.not. associated(iodesc)) then
       call abortmp('no iodesc found for variable')
    end if
    call PIO_write_darray(ncdf%FileID, varptr%vardesc, iodesc, var, ierr)

    !$OMP END SINGLE
  end subroutine nf_put1DR
  subroutine nf_put2DR(ncdf, var, start, count, varid, name)
    type(nf_handle), intent(inout) :: ncdf
    real(kind=real_kind), intent(in) :: var(:,:)
    integer(kind=nfsizekind), intent(in) :: start(:), count(:)
    type(nf_variable), intent(in), optional, target :: varid
    character*(*), intent(in),optional :: name
    integer :: vindex, extent, ierr, msize, ndims
    type (nf_variable), pointer :: varptr
    type (io_desc_t), pointer :: iodesc
    integer :: id

    !$OMP SINGLE
    if(ncdf%state /= readystate) then
       print *,__FILE__,__LINE__,ncdf%state, readystate       
       call abortmp('sanity check failed, bad file handle')
    end if
    if(present(varid)) then
       varptr=>varid
    else if(present(name)) then
       vindex = get_varindex(name, ncdf%varlist)    
       varptr => ncdf%varlist(vindex)
    else
       call abortmp('one of the optional arguments name and varid must be provided')   
    end if

    if(varptr%timedependent) then
       ndims=size(start)
       call PIO_SetFrame(varptr%vardesc,start(ndims))
    endif
    msize = size(var)
    nullify(iodesc)
    do id=1,size(ncdf%decomplist)
       if(varptr%dimsid==ncdf%decomplist(id)%dimsid) then
          iodesc=>ncdf%decomplist(id)%iodesc
          exit
       end if
    end do
    if(.not. associated(iodesc)) then
       call abortmp('no iodesc found for variable')
    end if
    call PIO_write_darray(ncdf%FileID, varptr%vardesc, iodesc, reshape(var,(/msize/)), ierr)

    !$OMP END SINGLE
  end subroutine nf_put2DR


  subroutine nf_put3DR(ncdf, var, start, count, varid, name)
    type(nf_handle), intent(in) :: ncdf
    real(kind=real_kind), intent(in) :: var(:,:,:)
    integer(kind=nfsizekind), intent(in) :: start(:), count(:)
    type(nf_variable), intent(in), optional, target :: varid
    character*(*), intent(in),optional :: name
    integer :: vindex, extent
    type (nf_variable), pointer :: varptr

    !$OMP SINGLE
    if(ncdf%state /= readystate) then
       print *,__FILE__,__LINE__,ncdf%state, readystate       
       call abortmp('sanity check failed, bad file handle')
    end if
    if(present(varid)) then
       varptr=>varid
    else if(present(name)) then
       vindex = get_varindex(name, ncdf%varlist)    
       varptr => ncdf%varlist(vindex)
    else
       call abortmp('one of the optional arguments name and varid must be provided')   
    end if
    call abortmp('function not supported')
    !$OMP END SINGLE
  end subroutine nf_put3DR

  !
  ! returns the nf_variable associated with name in the nf_handle ncdf
  !
  function get_variable_handle(ncdf,name)
    type(nf_handle), intent(in) :: ncdf
    character*(*) :: name
    type(nf_variable), pointer :: get_variable_handle
    integer :: vindex

    vindex = get_varindex(name, ncdf%varlist)    

    if(vindex > 0) &
         get_variable_handle => ncdf%varlist(vindex)
    ! if not?

  end function get_variable_handle

  integer function decompdimsid(ncdf, dims) result(dimsid)
    type(nf_handle), intent(in) :: ncdf
    integer, intent(in) :: dims(:)
    integer :: ndims
    integer :: id
    ndims = size(dims)
    dimsid=0
    do id=1,ndims
       dimsid=dimsid+ncdf%dimlist(dims(id))%dimID*10**(id-1)
    end do

  end function decompdimsid


  subroutine nf_init_decomp(ncdf_list, dims, ldof,iodof,start,count)
    type(nf_handle), intent(inout), target :: ncdf_list(:)
    integer, intent(in) :: dims(:), ldof(:), iodof(:)
    integer(kind=nfsizekind) :: start(:), count(:)
    
    type(nf_handle), pointer :: ncdf
    type(nf_decomp) :: decomp
    type(nf_decomp), allocatable :: tmpdecomplist(:)
    integer, allocatable :: dimsize(:)
    integer :: ndims, ios, id, olddecompcnt, ierr

    ndims = size(dims)
    allocate(dimsize(ndims))

    do ios=1,max_output_streams
       ncdf =>ncdf_list(ios)
       decomp%dimsid=0

       if(ncdf%state > 0) then
          decomp%dimsid=decompdimsid(ncdf, dims)

          do id=1,ndims
             dimsize(id)=ncdf%dimlist(dims(id))%dsize
          end do
          ! special handling for data dependent only on time
          if(ndims==1.and.dimsize(1)==0) then
             dimsize(1)=1
          end if

          call PIO_initDecomp(piofs,pio_double,dimsize,ldof, decomp%iodesc)

          olddecompcnt=0
          if(associated(ncdf%decomplist)) then
             ! reallocate decomplist at the larger size
             olddecompcnt=size(ncdf%decomplist)
             allocate(tmpdecomplist(olddecompcnt),stat=ierr)
             if(ierr.ne.0) then 
                print *, __LINE__, ierr,olddecompcnt
                call abortmp('allocate error in pio_io_mod')
             end if
             tmpdecomplist = ncdf%decomplist
             deallocate(ncdf%decomplist)
             allocate(ncdf%decomplist(olddecompcnt+1),stat=ierr)
             if(ierr.ne.0) then 
                print *, __LINE__, ierr,olddecompcnt+1
                call abortmp('allocate error in pio_io_mod')
             end if
             
             ncdf%decomplist(1:olddecompcnt)=tmpdecomplist
             deallocate(tmpdecomplist)
          else
             allocate(ncdf%decomplist(1),stat=ierr)
             if(ierr.ne.0) then 
                call abortmp('allocate error in pio_io_mod')
             end if
          endif
          ncdf%decomplist(olddecompcnt+1)=decomp
       end if
    end do



  end subroutine nf_init_decomp

  !
  !  Variables that could be written to the output file must be registered with a call to 
  !  nf_output_register_variables, this call must occur after the call to nf_output_register_dims but 
  !  before the call to nf_output_init_complete.  Only variables that are listed as output_variables
  !  in the analysis_nl namelist are actually registered.
  !
  subroutine nf_output_register_variables(ncdf_list,varcnt, varname, dims,vartype_in,varrequired_in,lsize)
    type(nf_handle), intent(inout), target :: ncdf_list(:)
    integer, intent(in) :: varcnt
    character*(*), intent(in) :: varname(varcnt)

    integer, intent(in) :: dims(:,:)
    integer, intent(in), target, optional :: vartype_in(:)
    logical, intent(in), target, optional :: varrequired_in(:)
    integer, intent(in), optional :: lsize  ! local size of the decomposed dimension

    character(len=varname_len), pointer :: output_varnames(:)
    integer :: hDimID
    integer :: ierr, ios
    integer :: ivarid, i, ii, j, oldvarcnt, newvarcnt
    type(nf_handle), pointer :: ncdf
    type(nf_variable), pointer :: varptr, tmpvarlist(:)
    integer :: dimcnt, id, dimsid
    integer, allocatable :: vardims(:)
    logical, pointer :: varrequired(:)
    integer, pointer :: vartype(:)
    logical :: time_dep_var
    integer :: maxdims

    !$OMP SINGLE
    if(present(varrequired_in)) then
       varrequired=>varrequired_in
    else
       allocate(varrequired(varcnt))
       varrequired=.false.
    end if
    if(present(vartype_in)) then
       vartype=>vartype_in
    else
       allocate(vartype(varcnt))
       vartype=PIO_double
    end if

    do ios=1,max_output_streams
       ncdf=>ncdf_list(ios)
       if(ncdf%state == dimsstate .or. ncdf%state == varsstate) then

          output_varnames => get_current_varnames(ios)

          ncdf%state = varsstate
          maxdims = size(ncdf%dimlist)
          allocate(vardims(maxdims),stat=ierr)

          if(ierr.ne.0) then 
             print *, __LINE__, ierr, maxdims
             call abortmp('allocate error in pio_io_mod')
          end if

          !  Count the variables actually being written to this file
          if(associated(ncdf%varlist)) then
             oldvarcnt=size(ncdf%varlist)
          else
             oldvarcnt=0
          end if
          newvarcnt=oldvarcnt
          do i=1,varcnt
             if(varrequired(i).or.nf_selectedvar(varname(i),output_varnames)) then
                newvarcnt=newvarcnt+1
             endif
          end do

          if(newvarcnt.ne.oldvarcnt) then
             if(associated(ncdf%varlist)) then
                ! reallocate varlist at the larger size
                allocate(tmpvarlist(oldvarcnt),stat=ierr)
                if(ierr.ne.0) then 
                   print *, __LINE__, ierr,oldvarcnt
                   call abortmp('allocate error in pio_io_mod')
                end if
                tmpvarlist = ncdf%varlist
                deallocate(ncdf%varlist)
                allocate(ncdf%varlist(newvarcnt),stat=ierr)
                if(ierr.ne.0) then 
                   print *, __LINE__, ierr,newvarcnt
                   call abortmp('allocate error in pio_io_mod')
                end if

                ncdf%varlist(1:oldvarcnt)=tmpvarlist
                deallocate(tmpvarlist)
             else
                allocate(ncdf%varlist(newvarcnt),stat=ierr)
                if(ierr.ne.0) then 
                   print *, __LINE__, ierr,newvarcnt
                   call abortmp('allocate error in pio_io_mod')
                end if
             endif

             ii=oldvarcnt+1
             do i=1,varcnt
                time_dep_var=.false.
                vardims=unlim_dim
                dimcnt=0
                do j=1,size(dims,1)
                   if(dims(j,i) > 0) then
                      vardims(j)=ncdf%dimlist(dims(j,i))%dimID
                      if(ncdf%dimlist(dims(j,i))%dsize==0) time_dep_var=.true.
                      dimcnt=dimcnt+1
                   endif
                end do
                if(varrequired(i).or.nf_selectedvar(varname(i),output_varnames)) then
                   varptr=>ncdf%varlist(ii)
                   varptr%timedependent=time_dep_var
                   ii=ii+1

                   if(ierr.ne.0) then 
                      print *, __LINE__, ierr,newvarcnt
                      call abortmp('allocate error in pio_io_mod')
                   end if
                   varptr%ndims=dimcnt
                   varptr%varname=varname(i)
                   varptr%required = varrequired(i)

                   ierr = PIO_def_var(ncdf%FileID, varname(i), vartype(i), vardims(1:dimcnt), &
                        varptr%varDesc)
	           varptr%ivarid=varptr%vardesc%varid

	           if(varptr%timedependent .and. (dimcnt>1)) then
                      varptr%dimsid = decompdimsid(ncdf, vardims(1:dimcnt-1))
                   else
                      varptr%dimsid = decompdimsid(ncdf, vardims(1:dimcnt))
                   end if
                   
                end if
             end do
          end if
          deallocate(vardims)
       end if
    end do
    !$OMP END SINGLE
  end subroutine nf_output_register_variables


  !
  !  Add the standard text attributes long_name and units to a variable 
  !

  subroutine nf_variable_attributes(ncdf_list,varname,long_name,units,otherattname,otherattval)
    implicit none
    type(nf_handle), intent(in),target :: ncdf_list(:)
    character*(*),intent(in) :: varname
    character*(*),intent(in),optional :: long_name
    character*(*),intent(in),optional :: units
    character*(*),intent(in),optional :: otherattname
    character*(*),intent(in),optional :: otherattval
    type(nf_variable), pointer :: var
    type(nf_handle), pointer :: ncdf
    integer :: vindex, ios, ierr

    !$OMP SINGLE
    do ios = 1, max_output_streams
       ncdf => ncdf_list(ios)
       if(ncdf%state==varsstate) then
          vindex = get_varindex(varname, ncdf%varlist)
          if(vindex > 0) then
             var => ncdf%varlist(vindex)

             if(present(long_name)) then
                ierr = PIO_put_att(ncdf%FileID,var%vardesc%varid,'long_name',trim(long_name))
             end if
             if(present(units)) then
                ierr = PIO_put_att(ncdf%FileID,var%vardesc%varid,'units',trim(units))
             end if
             if(present(otherattname)) then
                ierr = PIO_put_att(ncdf%FileID,var%vardesc%varid,otherattname,trim(otherattval))
             end if
          end if
       end if
    end do
    !$OMP END SINGLE

  end subroutine nf_variable_attributes

  !
  ! Open files for each nf_handle in the array ncdf  
  !
   !
  subroutine nf_output_init_begin(ncdf, masterproc,nprocs,rank,comm,file_prefix,runtype)
    implicit none
    type(nf_handle), intent(out) :: ncdf(:)
    logical, intent(in) :: masterproc
    integer, intent(in) :: nprocs
    integer, intent(in) :: rank
    integer, intent(in) :: comm

    integer, intent(in) :: runtype
    character*(*),intent(in) :: file_prefix


    integer :: ios

    !
    ! Loop through output_streams, identify which will be used and open files for them
    !
    !$OMP SINGLE
    call PIO_Init(rank,comm,num_io_procs,num_agg,io_stride,PIO_rearr_box,PIOFS)
    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0) .and. (output_start_time(ios) .lt. output_end_time(ios))) then 
          ncdf(ios)%iframe=1

          call nf_open_file(masterproc,nprocs,comm,rank,ios,&
               output_prefix,file_prefix,runtype,ncdf(ios)%ncFileID, ncdf(ios)%FileID)
          ncdf(ios)%state = beginstate
       else
          ncdf(ios)%ncFileID=-1
          ncdf(ios)%state = -1
       end if
       nullify(ncdf(ios)%varlist)
       nullify(ncdf(ios)%dimlist)
       nullify(ncdf(ios)%decomplist)
    end do

    !$OMP END SINGLE 
  end subroutine nf_output_init_begin

  subroutine nf_output_register_dims(ncdf_list,dimcnt, dimname, dsize)
    implicit none
    type(nf_handle), intent(inout), target :: ncdf_list(:)
    integer, intent(in) :: dimcnt
    character*(*), intent(in) :: dimname(dimcnt)
    integer, intent(in) :: dsize(dimcnt)
    type(nf_dim), pointer :: tmpdimlist(:)
    type(nf_handle), pointer :: ncdf

    integer :: i, ierr, ios
    integer :: olddimcnt, newdimcnt
    !$OMP SINGLE

    do ios=1,max_output_streams
       ncdf=> ncdf_list(ios)
       if(ncdf%state == dimsstate .or. ncdf%state == beginstate) then
       
          !
          ! We should make sure the user is not trying to define the same dimension again
          !
          ncdf%state = dimsstate
          if(associated(ncdf%dimlist)) then
             olddimcnt = size(ncdf%dimlist)
             newdimcnt = olddimcnt+dimcnt
             ! reallocate dimlist at the larger size
             allocate(tmpdimlist(olddimcnt),stat=ierr)
             if(ierr.ne.0) then 
                print *, __LINE__, ierr,olddimcnt
                call abortmp('allocate error in pio_io_mod')
             end if
             
             tmpdimlist = ncdf%dimlist
             deallocate(ncdf%dimlist)
             allocate(ncdf%dimlist(newdimcnt),stat=ierr)
             if(ierr.ne.0) then 
                print *, __LINE__, ierr,newdimcnt
                call abortmp('allocate error in pio_io_mod')
             end if
             ncdf%dimlist(1:olddimcnt)=tmpdimlist
             deallocate(tmpdimlist)
          else
             allocate(ncdf%dimlist(dimcnt),stat=ierr)
             if(ierr.ne.0) then 
                print *, __LINE__, ierr,dimcnt
                call abortmp('allocate error in pio_io_mod')
             end if
          endif
          do i=1,dimcnt
             ierr = PIO_def_dim(ncdf%FileID, dimname(i), dsize(i), &
                  ncdf%dimlist(i)%DimID)
             ncdf%dimlist(i)%dsize=dsize(i)
             if(dsize(i).eq.0) then
                ncdf%timeDimID=ncdf%dimlist(i)%DimID
             end if
             ncdf%dimlist(i)%dimname=dimname(i)
          end do
       end if
    end do
    !$OMP END SINGLE

  end subroutine nf_output_register_dims



  subroutine nf_output_init_complete(ncdf_list )
    type(nf_handle), intent(in), target :: ncdf_list(:)
    type(nf_handle), pointer :: ncdf
    integer :: ios, ierr

    !$OMP SINGLE
    do ios=1,max_output_streams
       ncdf=>ncdf_list(ios)
       if(ncdf%state == varsstate) then
          ncdf%state = readystate
          ierr = PIO_enddef(ncdf%FileID)
       end if
    end do
    ! ==========================================
    ! Deallocate the temporary coordinate arrays 
    ! =========================================
    !$OMP END SINGLE
  end subroutine nf_output_init_complete

  subroutine check(status, line)
    integer, intent ( in) :: status
    integer, intent(in), optional :: line

    if(status /= PIO_noerr) then
#ifdef AIX
       call xl__trbk()
#endif
       call abortmp("pio_io_mod Error")
    end if
  end subroutine check

  ! close the files and deallocate the lists
  subroutine nf_close_all(ncdf_list)
    type(nf_handle), intent(inout) :: ncdf_list(:)
    integer :: ios

    do ios=1,max_output_streams
       call nf_close(ncdf_list(ios))
    end do

  end subroutine nf_close_all

  subroutine nf_close(ncdf)
    type(nf_handle), intent(inout) :: ncdf
    ! Close netCDF file

    !$OMP SINGLE

    if(ncdf%state>0) then              
       call PIO_closeFile(ncdf%FileID)
       ncdf%state=-1
       if(associated(ncdf%varlist)) then
          deallocate(ncdf%varlist)
       end if
       if(associated(ncdf%dimlist)) then
          deallocate(ncdf%dimlist)
       end if
    end if

    !$OMP END SINGLE
  end subroutine nf_close

  subroutine nf_open_file(masterproc,nprocs,comm,iam,ios,output_prefix,file_prefix,runtype,ncFileID, FileID)
    use pio_types, only : file_desc_t, PIO_64BIT_OFFSET ! _EXTERNAL
    logical, intent(in) :: masterproc
    integer, intent(in) :: nprocs
    integer, intent(in) :: comm
    integer, intent(in) :: iam, ios
    character*(*),intent(in) :: file_prefix, output_prefix
    integer, intent(in) :: runtype
    integer, intent(out) :: ncFileID
    type(File_desc_t), intent(out) :: FileID
    character(len=MAX_STRING_LEN) :: filename,charnum

    integer :: ierr

    write(filename,'(a,i1.1,a)') &
         TRIM(ADJUSTL(output_dir))//TRIM(output_prefix)//TRIM(ADJUSTL(file_prefix)),ios,".nc"

    if(runtype==0) then
       if(output_type.eq.'netcdf') then
          if(masterproc) print *, 'Opening file ',trim(filename), ' using netcdf'
          ierr = PIO_CreateFile(PIOFS, FileID, iotype_netcdf  ,trim(filename), PIO_64BIT_OFFSET)
       else
          if(masterproc) print *, 'Opening file ',trim(filename), ' using pnetcdf'
          ierr = PIO_CreateFile(PIOFS, FileID, iotype_pnetcdf  ,trim(filename), PIO_64BIT_OFFSET)
       end if
    else
       ! this code is broken, as our NETCDF code cannot append to existing files:
       !if(masterproc) print *, 'Appending file ', filename
       !ierr = PIO_OpenFile(FileID, trim(filename), 1 ),__LINE__)
       !ierr = PIO_Redef(FileID))

       ! for now, just assume the calling program created a new, different filename:
       if(masterproc) print *, 'Opening file ',trim(filename),' using pnetcdf'
       ierr = PIO_CreateFile(PIOFS, FileID, iotype_pnetcdf, trim(filename), PIO_64BIT_OFFSET)
    end if
    ncFileID=fileid%fh
  end subroutine nf_open_file
#endif
end module pio_io_mod
