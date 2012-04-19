#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
module netcdf_io_mod
  use kinds, only : int_kind, real_kind
  !HOMME Specific: abortmp, mpireal_t,mpiinteger_t
  use parallel_mod, only : abortmp, mpireal_t, mpiinteger_t

  use pio !_EXTERNAL
  use pio_types ! _EXTERNAL	
  use control_mod, only : MAX_STRING_LEN         !HOMME Specific: MAX_STRING_LEN
  use common_io_mod, only : varname_len, output_start_time, output_end_time, &
       output_frequency, max_output_streams, output_varnames1, output_varnames2, &
       output_varnames3, output_varnames4, output_varnames5, max_output_streams, &
       num_io_procs, num_agg, io_stride, output_dir, output_prefix, max_output_variables, &
       nf_selectedvar, get_varindex, get_current_varnames, get_varindex, &
       nf_dim, nf_variable, nf_handle, beginstate, dimsstate, varsstate, readystate, &
       nfsizekind, unlim_dim, output_type
  use common_io_mod, only : nf_double, nf_noerr

  implicit none  
  private
  ! end of analysis_nl namelist variables
  type(io_desc_t), public :: iodesc3d, iodesc2d, iodesct, iodesc2d_nc, iodesc3d_nc, iodesc3d_subelem
  type(iosystem_desc_t),public :: pio_subsystem


  interface nf_put_var
     module procedure nf_put0DR
     module procedure nf_put1DI
     module procedure nf_put1DR
     module procedure nf_put2DR
     module procedure nf_put2DI
  end interface

  ! public interfaces  
  public :: nf_output_init_begin, nf_output_init_complete, nf_output_register_variables, &
       nf_put_var, nf_close, nf_close_all,  &
       nf_output_register_dims, nf_advance_frame,  &
       nf_variable_attributes, nf_get_frame, nf_global_attribute

  !
  ! output file information is contained in the linked list nf_handle, that list in turn contains
  !  a linked list of output variable information for that file.


contains


  !
  !  Advances the time (unlimited dimension) level of the associated file
  !


  subroutine nf_advance_frame(ncdf)
    type(nf_handle), intent(inout) :: ncdf
    !$OMP SINGLE
    ncdf%iframe=ncdf%iframe+1
    !
    ! This call flushs the i/o buffers after a write, it may affect performance
    ! and can be removed if nessesary
    ! Frames are advanced on a per variable basis in pio
    !$OMP END SINGLE
  end subroutine nf_advance_frame

  integer function nf_get_frame(ncdf)
    type(nf_handle), intent(in) :: ncdf
    nf_get_frame=ncdf%iframe
  end function nf_get_frame



  subroutine nf_put0DR(ncdf, var, start, count, varid, name)
    type(nf_handle), intent(inout) :: ncdf
    real(kind=real_kind), intent(in) :: var
    integer(kind=nfsizekind), intent(in) :: start(:), count(:)
    type(nf_variable), intent(in), optional, target :: varid
    character*(*), intent(in),optional :: name
    integer :: vindex, extent, ierr
    type (nf_variable), pointer :: varptr
    real(kind=real_kind) :: vartmp(1)

    !$OMP SINGLE
    if(ncdf%state /= readystate) then
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
    ierr = PIO_Put_var(ncdf%FileID, varptr%vardesc%varid, int(start), &
         int(count), (/var/))
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
    !$OMP SINGLE

    if(ncdf%state /= readystate) then
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
    call PIO_write_darray(ncdf%FileID, varptr%vardesc, iodesct, var, ierr)

    !$OMP END SINGLE
  end subroutine nf_put1DI
  subroutine nf_put1DR(ncdf, var, start, count, varid, name, iodescin)
    type(nf_handle), intent(inout) :: ncdf
    real(kind=real_kind), intent(in) :: var(:)
    integer(kind=nfsizekind), intent(in) :: start(:), count(:)
    type(nf_variable), intent(in), optional, target :: varid
    character*(*), intent(in),optional :: name
    type(io_desc_t), intent(inout), optional :: iodescin

    integer :: vindex, extent, ierr
    type (nf_variable), pointer :: varptr

    !$OMP SINGLE
    if(ncdf%state /= readystate) then
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

!    if(associated(varptr%vardesc%iodesc)) then
       if(varptr%timedependent) call PIO_SetFrame(varptr%vardesc,start(2))
       if(present(iodescin)) then
          call PIO_write_darray(ncdf%FileID, varptr%vardesc, iodescin, var, ierr)
       else	
          call PIO_write_darray(ncdf%FileID, varptr%vardesc, iodesc2d, var, ierr)
       end if
    

    !$OMP END SINGLE
  end subroutine nf_put1DR


  subroutine nf_put2DI(ncdf, var, start, count, varid, name, iodescin)
    type(nf_handle), intent(inout) :: ncdf
    integer, intent(in) :: var(:,:)
    integer(kind=nfsizekind), intent(in) :: start(:), count(:)
    type(nf_variable), intent(in), optional, target :: varid
    character*(*), intent(in),optional :: name
    type(io_desc_t), intent(inout), optional :: iodescin
    integer :: vindex, extent, ierr, msize
    type (nf_variable), pointer :: varptr
    !$OMP SINGLE
    if(ncdf%state /= readystate) then
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
    if(varptr%timedependent) call PIO_SetFrame(varptr%vardesc,start(3))
    msize = size(var)
    if(present(iodescin)) then
       call PIO_write_darray(ncdf%FileID, varptr%vardesc, iodescin, reshape(var,(/msize/)), ierr)
    else
       call PIO_write_darray(ncdf%FileID, varptr%vardesc, iodesc3d, reshape(var,(/msize/)), ierr)
    endif
    !$OMP END SINGLE
  end subroutine nf_put2DI



  subroutine nf_put2DR(ncdf, var, start, count, varid, name, iodescin)
    type(nf_handle), intent(inout) :: ncdf
    real(kind=real_kind), intent(in) :: var(:,:)
    integer(kind=nfsizekind), intent(in) :: start(:), count(:)
    type(nf_variable), intent(in), optional, target :: varid
    character*(*), intent(in),optional :: name
    type(io_desc_t), intent(inout), optional :: iodescin
    integer :: vindex, extent, ierr, msize
    type (nf_variable), pointer :: varptr
    !$OMP SINGLE
    if(ncdf%state /= readystate) then
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
    if(varptr%timedependent) call PIO_SetFrame(varptr%vardesc,start(3))
    msize = size(var)
    if(present(iodescin)) then
       call PIO_write_darray(ncdf%FileID, varptr%vardesc, iodescin, reshape(var,(/msize/)), ierr)
    else
       call PIO_write_darray(ncdf%FileID, varptr%vardesc, iodesc3d, reshape(var,(/msize/)), ierr)
    endif
    !$OMP END SINGLE
  end subroutine nf_put2DR







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
    integer :: dimcnt, maxdims
    integer, allocatable :: vardims(:)
    logical, pointer :: varrequired(:)
    integer, pointer :: vartype(:)
    logical :: time_dep_var

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
       vartype=nf_double
    end if
    
    do ios=1,max_output_streams
       ncdf=>ncdf_list(ios)
       if(ncdf%ncFileID>=0) then
          output_varnames => get_current_varnames(ios)

          if(ncdf%state /= dimsstate .and. ncdf%state /= varsstate) then
             call abortmp('sanity check failed, file handle in wrong state')
          end if
          ncdf%state = varsstate
          maxdims = size(ncdf%dimlist)
          allocate(vardims(maxdims),stat=ierr)
          
          if(ierr.ne.0) then 
             print *, __LINE__, ierr, maxdims
             call abortmp('allocate error in netcdf_io_mod')
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
                   call abortmp('allocate error in netcdf_io_mod')
                end if
                tmpvarlist = ncdf%varlist
                deallocate(ncdf%varlist)
                allocate(ncdf%varlist(newvarcnt),stat=ierr)
                if(ierr.ne.0) then 
                   print *, __LINE__, ierr,newvarcnt
                   call abortmp('allocate error in netcdf_io_mod')
                end if

                ncdf%varlist(1:oldvarcnt)=tmpvarlist
                deallocate(tmpvarlist)
             else
                allocate(ncdf%varlist(newvarcnt),stat=ierr)
                if(ierr.ne.0) then 
                   print *, __LINE__, ierr,newvarcnt
                   call abortmp('allocate error in netcdf_io_mod')
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
                      if(ncdf%dimlist(dims(j,i))%dsize==0) then
                         time_dep_var=.true.
                      endif
                      dimcnt=dimcnt+1
                   endif
                end do

                if(varrequired(i).or.nf_selectedvar(varname(i),output_varnames)) then
                   varptr=>ncdf%varlist(ii)
                   varptr%timedependent=time_dep_var
                   ii=ii+1

                   if(ierr.ne.0) then 
                      print *, __LINE__, ierr,newvarcnt
                      call abortmp('allocate error in netcdf_io_mod')
                   end if
                   varptr%ndims=dimcnt
                   varptr%varname=varname(i)
                   varptr%required = varrequired(i)
                   !print *,__FILE__,__LINE__,varname(i),vartype(i),vardims(1:dimcnt)

                   ierr = PIO_def_var(ncdf%FileID, varname(i), vartype(i), vardims(1:dimcnt), &
                        varptr%varDesc)
                   varptr%ivarid=ivarid
                end if
             end do
          end if
          deallocate(vardims)
       end if
    end do
    !$OMP END SINGLE
  end subroutine nf_output_register_variables

  subroutine nf_global_attribute(ncdf_list,attname, attval)
    implicit none
    type(nf_handle), target, intent(in) :: ncdf_list(:)
    character(*), intent(in) :: attname
    integer, intent(in) :: attval
    type(nf_handle), pointer :: ncdf   
    integer :: ios, ierr
    integer(kind=nfsizekind), parameter :: one=1


    do ios = 1, max_output_streams
       ncdf => ncdf_list(ios)
       if(ncdf%state==varsstate) then
          ierr = PIO_Put_att(ncdf%fileid,pio_global,attname,attval)
       endif
    enddo

  end subroutine nf_global_attribute

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

!	print *, __FILE__,__LINE__,varname, var%ivarid, long_name
             if(present(long_name)) then
                ierr = PIO_put_att(ncdf%FileID,var%varDesc%varid,'long_name',trim(long_name))
             end if
             if(present(units)) then
                ierr = PIO_put_att(ncdf%FileID,var%varDesc%varid,'units',trim(units))
             end if
             if(present(otherattname)) then
                ierr = PIO_put_att(ncdf%FileID,var%varDesc%varid,otherattname,trim(otherattval))
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
    use parallel_mod, only : haltmp
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
    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0) .and. (output_start_time(ios) .lt. output_end_time(ios))) then 
          ncdf(ios)%iframe=1

          call PIO_Init(rank,comm,num_io_procs,num_agg,io_stride,&
		PIO_REARR_BOX,pio_subsystem)

          call nf_open_file(masterproc,nprocs,comm,rank,ios,&
               output_prefix,file_prefix,runtype,ncdf(ios)%ncFileID, ncdf(ios)%FileID)

          ncdf(ios)%state = beginstate
       else
          ncdf(ios)%ncFileID=-1
          ncdf(ios)%state = -1
       end if
       nullify(ncdf(ios)%varlist)
       nullify(ncdf(ios)%dimlist)
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
       if((ncdf%ncFileID) >= 0) then

          !
          ! We should make sure the user is not trying to define the same dimension again
          !
          if(ncdf%state /= dimsstate .and. ncdf%state /= beginstate) then
             call abortmp('sanity check failed, file handle in wrong state ')
          end if
          ncdf%state = dimsstate
          if(associated(ncdf%dimlist)) then
             olddimcnt = size(ncdf%dimlist)
             newdimcnt = olddimcnt+dimcnt
             ! reallocate dimlist at the larger size
             allocate(tmpdimlist(olddimcnt),stat=ierr)
             if(ierr.ne.0) then 
                print *, __LINE__, ierr,olddimcnt
                call abortmp('allocate error in netcdf_io_mod')
             end if

             tmpdimlist = ncdf%dimlist
             deallocate(ncdf%dimlist)
             allocate(ncdf%dimlist(newdimcnt),stat=ierr)
             if(ierr.ne.0) then 
                print *, __LINE__, ierr,newdimcnt
                call abortmp('allocate error in netcdf_io_mod')
             end if
             ncdf%dimlist(1:olddimcnt)=tmpdimlist
             deallocate(tmpdimlist)
          else
             allocate(ncdf%dimlist(dimcnt),stat=ierr)
             if(ierr.ne.0) then 
                print *, __LINE__, ierr,dimcnt
                call abortmp('allocate error in netcdf_io_mod')
             end if
          endif

          do i=1,dimcnt

             ierr = PIO_def_dim(ncdf%FileID, dimname(i), dsize(i), &
                  ncdf%dimlist(i)%DimID)
             ncdf%dimlist(i)%dsize=dsize(i)

             ncdf%dimlist(i)%dimname=dimname(i)
             if(dsize(i).eq.0) then
                ncdf%timeDimID=ncdf%dimlist(i)%DimID
             end if
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
       if((ncdf%ncfileID)>=0) then
          if(ncdf%state /= varsstate) then
             call abortmp('sanity check failed, file handle in wrong state ')
          end if
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


    if(status /= nf_noerr) then
       print *, status
!       if(present(line)) then
!          print *, trim(nfmpi_strerror(status)),' at line ',line, ' of file ',__FILE__
!       else
!          print *, trim(nfmpi_strerror(status))
!       end if
#ifdef AIX
       call xl__trbk()
#endif
       call abortmp("netcdf_io_mod Error")
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
       deallocate(ncdf%varlist)
       deallocate(ncdf%dimlist)
    end if
    !$OMP END SINGLE
  end subroutine nf_close

  subroutine nf_open_file(masterproc,nprocs,comm,iam,ios,output_prefix,file_prefix,runtype,ncFileID, FileID)
      logical, intent(in) :: masterproc
      integer, intent(in) :: nprocs
      integer, intent(in) :: comm
      integer, intent(in) :: iam, ios
      character*(*),intent(in) :: file_prefix, output_prefix
      integer, intent(in) :: runtype
      integer, intent(out) :: ncFileID
      type(File_desc_t), intent(out) :: FileID
      character(len=MAX_STRING_LEN) :: filename
      integer :: output_info,ierr

      write(filename,'(a,i1.1,a)') &
           TRIM(ADJUSTL(output_dir))//TRIM(output_prefix)//TRIM(ADJUSTL(file_prefix)),ios,".nc"

      if(output_type.eq.'netcdf') then
         ierr = PIO_CreateFile(pio_subsystem, FileID, iotype_netcdf, filename, PIO_64BIT_OFFSET)
      else
         ierr = PIO_CreateFile(pio_subsystem, FileID, iotype_pnetcdf, filename, PIO_64BIT_OFFSET)
      endif

      if(masterproc) print *, 'Opening file ',trim(filename), fileid%fh, output_type
    end subroutine nf_open_file
  end module netcdf_io_mod


