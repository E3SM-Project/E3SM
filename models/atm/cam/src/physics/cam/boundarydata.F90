#define _FILE 'physics/cam1/boundarydata.F90 '
module boundarydata
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use spmd_utils,      only: masterproc
  use ppgrid,          only: pcols, pver, begchunk, endchunk
  use physics_types,   only: physics_state
  use abortutils,      only: endrun
#if ( defined SPMD )
  use mpishorthand,    only: mpicom, mpir8, mpiint
#endif
  use netcdf
  use error_messages,  only: handle_ncerr
  use cam_logfile,     only: iulog
  implicit none
  private
  integer, parameter :: ptrtim=12, ptrlon=1

  type boundarydata_type
     integer           :: ncid
     integer           :: fieldcnt
     integer           :: nm
     integer           :: np
     integer           :: latsiz
     integer           :: levsiz
     integer           :: ncolsiz
     integer           :: timsiz
     integer           :: vertextrap
     logical           :: iszonal, isncol
     integer           :: ndims
     integer           :: thistimedim
     integer           :: psid
     integer  :: map(4)
     integer  :: dimids(4)
     integer,  pointer :: dataid(:) => null()
     integer,  pointer :: columnmap(:) => null()
     integer, pointer  :: start(:,:) => null()
     integer, pointer  :: count(:,:) => null()
     real(r8), pointer :: lat(:) => null()
     real(r8), pointer :: zi(:) => null()
     real(r8), pointer :: pin(:) => null()
     real(r8), pointer :: cdates(:) => null()
     real(r8), pointer :: fields(:,:,:,:,:) => null()
     real(r8), pointer :: datainst(:,:,:,:) => null()
     real(r8), pointer :: hybi(:) => null()
     real(r8), pointer :: ps(:,:,:) => null()
  end type boundarydata_type

  public boundarydata_init
  public boundarydata_update
  public boundarydata_vert_interp
  public boundarydata_type

contains
  subroutine boundarydata_init(bndyfilename,phys_state,fieldnames,fieldcnt,bndydata,vertextrap)
    implicit none    
    character(len=*),intent(in) :: bndyfilename
    type(physics_state), intent(in):: phys_state(begchunk:endchunk)
    integer,intent(in) :: fieldcnt
    character(len=*), intent(in) :: fieldnames(fieldcnt)
    type(boundarydata_type),intent(out) :: bndydata 
    integer,intent(in), optional :: vertextrap ! if 0 set values outside output grid to 0
                                               ! if 1 set to boundary value
                                               ! if 2 set to cyclic boundaries 
                                               ! if 3 leave on input data grid and extrapolate later
    real(r8), pointer :: datain(:,:,:,:,:)
    integer :: lchnk

    bndydata%fieldcnt=fieldcnt
    if(present(vertextrap)) then
       bndydata%vertextrap=vertextrap
    else
       bndydata%vertextrap=0
    end if
    nullify(bndydata%fields)

    call boundarydata_read(phys_state,bndyfilename,fieldcnt,fieldnames,bndydata,datain)

    if(bndydata%iszonal) then
       call boundarydata_interpolate(phys_state,datain,bndydata)

       allocate(bndydata%datainst(size(bndydata%fields,1),size(bndydata%fields,2), &
            begchunk:endchunk,bndydata%fieldcnt))
    
       deallocate(datain)    
    end if
  end subroutine boundarydata_init

  subroutine boundarydata_update(phys_state, bndydata, update_out)
    use interpolate_data,only : get_timeinterp_factors
    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
    type(boundarydata_type), intent(inout) :: bndydata
    logical, intent(out), optional :: update_out
    real(r8) :: cdate
    integer :: nm, np, lchnk, j, k, fld, cols, cole, ncol, ndims
    real(r8) :: fact1, fact2
    real(r8), allocatable :: datain(:,:,:,:,:)
    logical :: update
    integer :: latspan
    integer :: kmax
    integer :: count(4), start(4), ierr

    
    call get_data_bounding_date_indices(bndydata%cdates,bndydata%nm,bndydata%np,cdate,update)
    if(present(update_out)) update_out=update
    nm= bndydata%nm
    np= bndydata%np

    call get_timeinterp_factors(.true., np, bndydata%cdates(nm), bndydata%cdates(np), &
         cdate, fact1, fact2, _FILE)

    if(size(bndydata%fields,5).eq.2) then 
       nm=1
       np=2
       if(update) then ! we need to read in the next month and interpolate
          if(bndydata%isncol) then
             bndydata%fields(:,:,:,:,nm)=bndydata%fields(:,:,:,:,np)
             do lchnk=begchunk,endchunk
                ncol=phys_state(lchnk)%ncol
                cols=1
                cole=cols+bndydata%count(cols,lchnk)-1
                do while(cole<=ncol)
                   
                   if(bndydata%levsiz==1) then               
                      ndims=2
                      start=(/bndydata%start(cols,lchnk),bndydata%np,-1,-1/)
                      count=(/bndydata%count(cols,lchnk),1,-1,-1/)
                   else
                      ndims=3
                      start=(/bndydata%start(cols,lchnk),bndydata%levsiz,bndydata%np,-1/)
                      count=(/bndydata%count(cols,lchnk),1,1,-1/)
                   end if
                   do fld=1,bndydata%fieldcnt
                      call handle_ncerr( nf90_get_var(bndydata%ncid, bndydata%dataid(fld) , &
                           bndydata%fields(cols:cole,:,lchnk,fld,np),  &
                           start(1:ndims), count(1:ndims)),&
                           _FILE,__LINE__)

                   end do
                   if(cols==ncol) exit
                   cols=cols+bndydata%count(cols,lchnk)
                   cole=cols+bndydata%count(cols,lchnk)-1
                end do
             end do

          else
             allocate(datain(ptrlon,bndydata%levsiz,bndydata%latsiz,1,bndydata%fieldcnt))
             if(masterproc) then
                count(1)=ptrlon
                count(2)=bndydata%levsiz
                count(3)=bndydata%latsiz
                count(4)=1
                start(1)=1
                start(2)=1
                start(3)=1
                start(4)=bndydata%np
                write(iulog,*) 'boundarydata reading data for month: ',bndydata%np
                do fld=1,bndydata%fieldcnt
                   call handle_ncerr( nf90_get_var(bndydata%ncid, bndydata%dataid(fld), &
                        datain(:,:,:,:,fld), start, count),_FILE,__LINE__)
                end do
             end if
#ifdef SPMD
             call mpibcast (datain, bndydata%levsiz*bndydata%latsiz*1*bndydata%fieldcnt, mpir8, 0, mpicom, ierr)
#endif          
             bndydata%fields(:,:,:,:,nm) = bndydata%fields(:,:,:,:,np) 
             call boundarydata_interpolate(phys_state,datain,bndydata)
             deallocate(datain)
          end if
       end if
    end if
    kmax = size(bndydata%fields,2)
   
    do fld=1,bndydata%fieldcnt
       do lchnk=begchunk,endchunk
          if(bndydata%isncol) then
             latspan = phys_state(lchnk)%ncol
          else
             latspan=phys_state(lchnk)%ulatcnt
          end if
          do k=1,kmax
             do j=1,latspan
                bndydata%datainst(j,k,lchnk,fld)=bndydata%fields(j,k,lchnk,fld,nm)*fact1 + &
                     bndydata%fields(j,k,lchnk,fld,np)*fact2
             end do
          end do
       end do
    end do
  end subroutine boundarydata_update


  subroutine boundarydata_read(phys_state,bndyfilename,fieldcnt,fieldnames,bndydata,datain)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Do initial read of time-variant boundary dataset, containing
    ! 12 monthly fields as a function of latitude and pressure.  Determine the two
    ! consecutive months between which the current date lies.
    ! 
    ! Method: 
    ! 
    ! Author: NCAR CMS
    !-----------------------------------------------------------------------
    use ioFileMod, only : getfil


    implicit none    
    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
    character(len=*),intent(in) :: bndyfilename
    integer,intent(in) :: fieldcnt
    character(len=*), intent(in) :: fieldnames(fieldcnt)
    type(boundarydata_type), intent(inout) :: bndydata
    real(r8), pointer :: datain(:,:,:,:,:)  ! 
    !
    ! Local variables
    !
    integer :: londimid
    integer :: latdimid
    integer :: ncoldimid
    integer :: levdimid
    integer :: ilevdimid
    integer :: timdimid
    integer :: ndims
    integer :: dimlen
    integer :: ilevsiz
    integer :: ncolsiz
    character(len=nf90_max_name) :: dimname


!    integer ::  ncid              ! netcdf id for file
    integer ::  dateid            ! netcdf id for date variable
    integer ::  secid             ! netcdf id for seconds variable
    integer ::  lonid             ! netcdf id for longitude variable
    integer ::  ncolid             ! netcdf id for longitude variable
    integer ::  latid             ! netcdf id for latitude variable
    integer ::  levid             ! netcdf id for level variable
    integer ::  timid             ! netcdf id for time variable
    integer :: hybid

    integer ::  dataid  ! netcdf id for data fields

    integer ::  lonsiz            ! size of longitude dimension on tracer dataset
    integer ::  levsiz            ! size of level dimension on tracer dataset
    integer ::  latsiz            ! size of latitude dimension on tracer dataset

    integer ::  j,n,k,nt,id          ! indices
    integer ::  ki,ko,ji,jo       ! indices
    integer :: date_tr(ptrtim), sec_tr(ptrtim)

    integer :: dimids(4), start(4), count(4)
    integer, pointer :: columnmap(:)
    real(r8) :: calday        ! current calendar day
    real(r8), pointer :: pin(:)
    real(r8), allocatable :: tmp_ps(:,:), tmp_fld(:,:,:)
    integer :: mincid,maxcid
    real(r8), allocatable, target:: lati(:)
    integer :: cols, cole
    integer ::  ierr, dimcnt
    integer ::  i, ncol, lchnk
    character(len=256) :: locfn    ! netcdf local filename to open 

    !
    !-----------------------------------------------------------------------
    !
    ! SPMD: Master reads dataset and does broadcast.  All subsequent interpolation is
    !       done in every process.  This is not required, one could remove this conditional
    !       and read the dataset independently on each task.
    !
    if(masterproc) then
       write(iulog,*)'boundarydata_read: Reading from: ', trim(bndyfilename)
#ifndef USE_MASTERPROC
    end if
#endif
    call getfil(bndyfilename, locfn)
    call handle_ncerr( nf90_open(locfn, 0, bndydata%ncid),&
         _FILE,__LINE__)

    !       write(iulog,*)'boundarydata_read: NCOPN returns id ',bndydata%ncid,' for file ',trim(locfn)
    !
    !------------------------------------------------------------------------
    ! Read tracer data
    !------------------------------------------------------------------------
    !
    ! Get dimension info
    !
    nullify(columnmap)
    nullify(pin)
    call handle_ncerr( nf90_inquire(bndydata%ncid, bndydata%ndims, unlimiteddimid=timdimid), &
         _FILE,__LINE__)
    ncolsiz=-1
    levsiz=-1
    lonsiz=-1
    latsiz=-1
    do i=1,bndydata%ndims
       call handle_ncerr( nf90_inquire_dimension(bndydata%ncid, i, dimname, dimlen),&
            _FILE,__LINE__)
       if (dimname(1:3).eq.'lat') then
          latdimid=i
          latsiz=dimlen
       else if (dimname(1:3) .eq.'lon') then
          londimid=i
          lonsiz=dimlen
       else if (dimname(1:4) .eq. 'ncol') then
          ncoldimid=i
          ncolsiz=dimlen
          bndydata%isncol=.true. 
       else if (dimname(1:3) .eq. 'lev') then
          levdimid=i
          levsiz=dimlen
       else if (dimname(1:4) .eq. 'ilev') then
          ilevdimid=i
          ilevsiz=dimlen
       else if (dimname(1:4) .eq. 'time') then
          if(timdimid/=i) then
             timdimid=i
          end if
          bndydata%timsiz=dimlen
       else
          write(iulog,*) 'Warning: do not know how to handle dimension ',&
               trim(dimname), ' in boundarydata.F90:313'
       end if
    end do

    if (latsiz>0 .and. lonsiz<=1) then
       bndydata%iszonal=.true.
       allocate(bndydata%lat(latsiz))
    end if
    if(bndydata%isncol) then
!       allocate (columnmap(ncolsiz))
!       call handle_ncerr( nf90_inq_varid(bndydata%ncid, 'ncol'    , ncolid),&
!            _FILE,__LINE__)
!       call handle_ncerr( nf90_get_var(bndydata%ncid,ncolid,columnmap), &
!            _FILE,__LINE__)
       if(levsiz>0) then
          allocate(bndydata%fields(pcols,levsiz,begchunk:endchunk,fieldcnt,2))

          ierr = nf90_inq_varid(bndydata%ncid, 'PS', bndydata%psid)
          if(ierr.eq.NF90_NOERR) then
             allocate(bndydata%ps(pcols,begchunk:endchunk,2))
             allocate(bndydata%hybi(levsiz+1))
             call handle_ncerr(nf90_inq_varid(bndydata%ncid,'hybi',hybid),&
                  _FILE,__LINE__)
             call handle_ncerr( nf90_get_var(bndydata%ncid, hybid, bndydata%hybi ),&
                  _FILE,__LINE__)
          else 
             call endrun('Did not recognize a vertical coordinate variable')
          end if
       else
          levsiz=1
          allocate(bndydata%fields(pcols,1,begchunk:endchunk,fieldcnt,2))
       end if
    else
       allocate(datain(lonsiz,levsiz,latsiz,2,fieldcnt))
       !
       ! Check dimension info
       !
       if (lonsiz/=ptrlon) then
          call endrun ('BOUNDARYDATA_READ: longitude dependence not implemented')
       endif

       if (bndydata%timsiz /= ptrtim) then
          write(iulog,*)'BOUNDARYDATA_READ: timsiz=',bndydata%timsiz,' must = ptrtim=',ptrtim
          call endrun
       end if
       if( bndydata%vertextrap.lt.3) then
          allocate(pin(levsiz))
       else
          allocate(bndydata%pin(levsiz))
          pin => bndydata%pin
       end if
       allocate(bndydata%lat(latsiz))


       allocate(datain(ptrlon,levsiz,latsiz,2,fieldcnt))

       call handle_ncerr( nf90_inq_varid(bndydata%ncid, 'lat'    , latid),&
            _FILE,__LINE__)
    end if
    !
    ! Determine necessary dimension and variable id's
    !
    allocate(bndydata%cdates(bndydata%timsiz))

    call handle_ncerr(nf90_inq_varid(bndydata%ncid, 'date'   , dateid), &
         _FILE,__LINE__)
    call handle_ncerr( nf90_get_var(bndydata%ncid, dateid, date_tr),&
         _FILE,__LINE__)
    ierr = nf90_inq_varid(bndydata%ncid, 'datesec', secid)
    if(ierr==NF90_NOERR) then
       call handle_ncerr( nf90_get_var(bndydata%ncid, secid , sec_tr),&
            _FILE,__LINE__)
    else
       sec_tr=0
    end if

    if (mod(date_tr(1),10000)/100 /= 1) then
       call endrun ('(boundarydata_read): error when cycling data: 1st month must be 1')
    end if
    if (mod(date_tr(ptrtim),10000)/100 /= 12) then
       call endrun ('(boundarydata_read): error when cycling data: last month must be 12')
    end if
    !
    !    return the calander dates of the file data
    !
    do n=1,ptrtim
       call bnddyi(date_tr(n), sec_tr(n), bndydata%cdates(n))
    end do
!    else
!       call handle_ncerr( nf90_inq_varid(bndydata%ncid, 'time', dateid),&
!            _FILE,__LINE__)
!
!       call handle_ncerr( nf90_get_var(bndydata%ncid, dateid, bndydata%cdates),&
!            _FILE,__LINE__)
!
!    end if
#ifdef USE_MASTERPROC
 else
    allocate(bndydata%cdates(ptrtim))
 end if
#ifdef SPMD
 call mpibcast (bndydata%cdates, ptrtim, mpir8, 0, mpicom, ierr)
#endif
#endif
 bndydata%nm=12
 bndydata%np=1
 call get_data_bounding_date_indices(bndydata%cdates,bndydata%nm,bndydata%np)
#ifdef USE_MASTERPROC
 if(masterproc) then
#endif
    !
    ! Obtain entire date and sec variables. Assume that will always
    ! cycle over 12 month data.
    !
    !
    ! Obtain input data latitude and level arrays.
    !
    if(bndydata%iszonal) then
       call handle_ncerr( nf90_get_var(bndydata%ncid, latid, bndydata%lat),&
            _FILE,__LINE__)
       ierr = nf90_inq_varid(bndydata%ncid, 'lev'    , levid)
       call handle_ncerr( nf90_get_var(bndydata%ncid, levid, pin ),&
            _FILE,__LINE__)
    end if

    allocate(bndydata%dataid(fieldcnt))
    if(masterproc) then
       write(iulog,*) 'boundarydata reading data for months: ',bndydata%nm,bndydata%np
    end if
    do i=1,fieldcnt
       call handle_ncerr( nf90_inq_varid(bndydata%ncid, fieldnames(i)   , bndydata%dataid(i)),&
            _FILE,__LINE__)
    end do
    if(bndydata%isncol) then
       allocate(bndydata%start(pcols,begchunk:endchunk), &
            bndydata%count(pcols,begchunk:endchunk))

!
!  For i/o efficiency we read in a block of data which includes the data needed on this 
!  processor but which may in fact include data not needed here.  physics cids are just the
!  offset into the file.
! 

       bndydata%start=-1
       bndydata%count=1
       mincid=2147483647
       maxcid=-1
       do lchnk=begchunk,endchunk
          ncol=phys_state(lchnk)%ncol
          i=minval(phys_state(lchnk)%cid(1:ncol)) 
          if(i < mincid) mincid = i
          i=maxval(phys_state(lchnk)%cid(1:ncol))
          if(i > maxcid) maxcid = i
       end do

       allocate(tmp_ps(mincid:maxcid,2))
       start=(/mincid,bndydata%nm,1,-1/)
       if(bndydata%np>bndydata%nm) then
          count=(/maxcid-mincid+1,2,-1,-1/)
       else
          count=(/maxcid-mincid+1,1,-1,-1/)
       end if
       if(associated(bndydata%ps) ) then
          call handle_ncerr( nf90_get_var(bndydata%ncid, bndydata%psid , &
               tmp_ps(:,1:count(2)), start(1:2), &
               count(1:2)),&
               _FILE,__LINE__)
          if(bndydata%np<bndydata%nm) then
             start(2)=bndydata%np
             call handle_ncerr( nf90_get_var(bndydata%ncid, bndydata%psid , &
                  tmp_ps(:,2:2), start(1:2), &
                  count(1:2)),&
                  _FILE,__LINE__)
          end if

          do lchnk=begchunk,endchunk
             do n=1,phys_state(lchnk)%ncol
                bndydata%ps(n ,lchnk,:) = tmp_ps(phys_state(lchnk)%cid(n),:)
             end do
          end do
          deallocate(tmp_ps)

       end if

       if(levsiz>1) then
          dimcnt=3
       else
          dimcnt=2
       end if
       start(2)=1
       count(2)=levsiz
                
       if(bndydata%np>bndydata%nm) then
          count(dimcnt)=2
       else
          count(dimcnt)=1
       end if
       start(dimcnt)=bndydata%nm

       allocate(tmp_fld(mincid:maxcid,count(2),2))

       do i=1,fieldcnt
          call handle_ncerr( nf90_get_var(bndydata%ncid, bndydata%dataid(i) , &
               tmp_fld(:,:,1:count(dimcnt)), &
               start(1:dimcnt), count(1:dimcnt)),&
               _FILE,__LINE__)

          do lchnk=begchunk,endchunk
             do n=1,phys_state(lchnk)%ncol
                bndydata%fields(n,:,lchnk,i,1:count(dimcnt)) = tmp_fld(phys_state(lchnk)%cid(n),:,:)
             end do
          end do
       end do
       if(bndydata%np<bndydata%nm) then
          start(dimcnt)=bndydata%np
          do i=1,fieldcnt
             call handle_ncerr( nf90_get_var(bndydata%ncid, bndydata%dataid(i) , &
                  tmp_fld(:,:,2:2), &
                  start(1:dimcnt), count(1:dimcnt)),&
                  _FILE,__LINE__)
             do lchnk=begchunk,endchunk
                do n=1,phys_state(lchnk)%ncol
                   bndydata%fields(n,:,lchnk,i,1:count(dimcnt)) = tmp_fld(phys_state(lchnk)%cid(n),:,:)
                end do
             end do
          end do
       end if
       deallocate(tmp_fld)
!       deallocate(columnmap)
    else
       !
       ! get the dimension orientation info from the first variable assume but verify that
       ! all variables requested have the same orientation
       !  
       allocate(bndydata%start(4,1),bndydata%count(4,1))
       call handle_ncerr( nf90_inquire_variable(bndydata%ncid,bndydata%dataid(1), &
            ndims=bndydata%ndims,dimids=bndydata%dimids),_FILE,__LINE__)

       bndydata%start=1
       do id=1,bndydata%ndims
          if(bndydata%dimids(id).eq.londimid) then
             bndydata%map(id)=1
             bndydata%count(id,1)=lonsiz
          else if(bndydata%dimids(id).eq.levdimid) then
             bndydata%map(id)=lonsiz
             bndydata%count(id,1)=levsiz
          else if(bndydata%dimids(id).eq.latdimid) then
             bndydata%map(id)=lonsiz
             if(any(bndydata%dimids.eq.levdimid)) bndydata%map(id)=bndydata%map(id)*levsiz
             bndydata%count(id,1)=latsiz
          else if(bndydata%dimids(id).eq.timdimid) then
             bndydata%map(id)=lonsiz*latsiz
             if(any(bndydata%dimids.eq.levdimid)) bndydata%map(id)=bndydata%map(id)*levsiz
             bndydata%count(id,1)=2
             bndydata%start(id,1)=bndydata%nm
             bndydata%thistimedim=id
          else
             write(iulog,*) __LINE__,fieldnames(i),id,bndydata%dimids(id),londimid, &
                  levdimid,latdimid,timdimid
             call endrun(_FILE)
          end if
       end do

       do i=1,fieldcnt
          call handle_ncerr( nf90_inq_varid(bndydata%ncid, fieldnames(i), bndydata%dataid(i)),&
               _FILE,__LINE__)

          call handle_ncerr( nf90_inquire_variable(bndydata%ncid,bndydata%dataid(i), &
               ndims=ndims,dimids=dimids),_FILE,__LINE__)
          if(ndims/=bndydata%ndims .or. dimids(1)/=bndydata%dimids(1).or.&
               dimids(2)/=bndydata%dimids(2) .or. dimids(3)/=bndydata%dimids(3)) then
             call endrun('Variable dims or order does not match')
          end if

          if(bndydata%np .gt. bndydata%nm) then
             call handle_ncerr( nf90_get_var(bndydata%ncid, bndydata%dataid(i), &
                  datain(:,:,:,:,i),bndydata%start(:,1), bndydata%count(:,1), &
                  map=bndydata%map),_FILE,__LINE__)
          else
             bndydata%count(bndydata%thistimedim,1)=1
             call handle_ncerr( nf90_get_var(bndydata%ncid, bndydata%dataid(i), &
                  datain(:,:,:,1:1,i), bndydata%start(:,1), bndydata%count(:,1), &
                  map=bndydata%map), _FILE,__LINE__)

             bndydata%start(bndydata%thistimedim,1)=bndydata%np
             call handle_ncerr( nf90_get_var(bndydata%ncid, bndydata%dataid(i), &
                  datain(:,:,:,2:2,i), bndydata%start(:,1), bndydata%count(:,1), &
                  map=bndydata%map), _FILE,__LINE__)

          endif
       end do

    end if
#ifdef USE_MASTERPROC
 end if
#ifdef SPMD 
 call mpibcast (levsiz, 1, mpiint, 0, mpicom, ierr)
 call mpibcast (latsiz, 1, mpiint, 0, mpicom, ierr)
#endif
#endif
 bndydata%latsiz=latsiz
 bndydata%levsiz=levsiz
#ifdef USE_MASTERPROC
 if(.not. masterproc) then
    if( bndydata%vertextrap.lt.3) then
       allocate(pin(levsiz))
    else
       allocate(bndydata%pin(levsiz))
       pin => bndydata%pin
    end if
    allocate(bndydata%lat(latsiz))
    allocate(datain(ptrlon,levsiz,latsiz,2,fieldcnt))
 endif
#ifdef SPMD
 call mpibcast (bndydata%lat, latsiz, mpir8, 0, mpicom, ierr)
 call mpibcast (pin, levsiz, mpir8, 0, mpicom, ierr)
 call mpibcast (datain, levsiz*latsiz*2*fieldcnt, mpir8, 0, mpicom, ierr)

#endif
#endif
 ! Convert input pressure from millibars to pascals.
 if(associated(pin)) then
    pin=pin*100._r8
    if(bndydata%vertextrap.lt.3) then
       allocate(bndydata%zi(levsiz))
    !
    !
    ! Convert input pressure levels to height (m).
    !
       do k=1,levsiz
          bndydata%zi(k) = 7.0e3_r8 * log (1.0e5_r8 / pin(k))
       end do
       deallocate(pin)
    end if
 end if
end subroutine boundarydata_read

  subroutine boundarydata_interpolate(phys_state, datain, bndydata)
    use ref_pres, only : pref_mid
    use interpolate_data,only : interp_type, lininterp_init, &
         lininterp_finish, lininterp
    use physconst, only: pi

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
    real(r8),intent(in)                  :: datain(:,:,:,:,:)
    type(boundarydata_type), intent(inout) :: bndydata
    type(interp_type) :: interp_wgts, lev_wgts
      
    integer :: k, lchnk, nt, j, fcnt
    real(r8) :: zo(pver)
    real(r8) :: lato(pcols)
    integer :: ulatcnt
    integer :: maxlatcnt
    integer :: timesize, tvalout

    !------------------------------------------------------------------------
    ! Interpolate tracer data to model grid
    !------------------------------------------------------------------------
    !
    ! Loop over all input times.
    !

    timesize=2

    maxlatcnt=1
    do lchnk=begchunk,endchunk
       maxlatcnt=max(maxlatcnt,phys_state(lchnk)%ulatcnt)
    end do
    if(bndydata%vertextrap.lt.3) then
       !
       ! Convert approximate cam pressure levels to height (m).
       !
       do k=1,pver
          zo (k) = 7.0e3_r8 * log (1.0e5_r8 / pref_mid(k))
       end do

       call lininterp_init(bndydata%zi,size(bndydata%zi),zo,pver,bndydata%vertextrap,lev_wgts)
       if(.not. associated(bndydata%fields)) then
          allocate(bndydata%fields(maxlatcnt,pver,begchunk:endchunk,bndydata%fieldcnt,timesize))
          bndydata%fields=0_r8
       end if
    else
       if(.not. associated(bndydata%fields)) then
          allocate(bndydata%fields(maxlatcnt,bndydata%levsiz,begchunk:endchunk,bndydata%fieldcnt,timesize))
          bndydata%fields=0_r8
       end if
    endif
    do lchnk=begchunk,endchunk
       ulatcnt=phys_state(lchnk)%ulatcnt

       !
       ! Convert cam model latitudes to degrees.
       ! Input model latitudes already in degrees.
       !
       do j=1,ulatcnt
          lato(j) = phys_state(lchnk)%ulat(j)*180._r8/pi          
       end do

       call lininterp_init(bndydata%lat,size(bndydata%lat),lato(1:ulatcnt),ulatcnt,1,interp_wgts)
       timesize =  size(datain,4)
       do fcnt=1,bndydata%fieldcnt
          do nt = 1,timesize
             if(timesize.gt.1) then
                tvalout=nt
             else
                tvalout=2
             end if
            if(bndydata%vertextrap.lt.3) then
                call lininterp(transpose(datain(1,:,:,nt,fcnt)),bndydata%latsiz,bndydata%levsiz, &
                     bndydata%fields(1:ulatcnt,:,lchnk,fcnt,tvalout), ulatcnt, pver, interp_wgts, lev_wgts) 
             else
                do k=1,bndydata%levsiz
                   call lininterp(datain(1,k,:,nt,fcnt),bndydata%latsiz, &
                        bndydata%fields(1:ulatcnt,k,lchnk,fcnt,tvalout), ulatcnt, interp_wgts)
                end do
             end if
          end do
       end do                 ! end loop over time samples
       call lininterp_finish(interp_wgts)
    end do
    if(bndydata%vertextrap.lt.3) &
         call lininterp_finish(lev_wgts)

    return
  end subroutine boundarydata_interpolate

  subroutine get_data_bounding_date_indices(cdates,nm,np, cdayout, update)
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
         is_perpetual
    real(r8), intent(in) :: cdates(ptrtim)
    real(r8), intent(out),optional :: cdayout
    logical, intent(out) ,optional :: update
    integer, intent(inout) :: nm, np
    integer :: n, np1
    real(r8) :: calday
    integer :: yr, mon, day   ! components of a date
    integer :: ncdate         ! current date in integer format [yyyymmdd]
    integer :: ncsec          ! current time of day [seconds]

    calday = get_curr_calday()
    if(present(cdayout)) cdayout=calday
    if(present(update)) update=.false. ! initialize output variable

    if(min(nm,np) .ge. 1 .and. max(nm,np) .le. 12 .and.  &
         calday>cdates(nm) .and. calday<=cdates(np)) return
    if((nm==12 .and. np==1) .and. (calday <= cdates(np) .or. &
         calday > cdates(nm))) return

    if(present(update)) update=.true.

    if(calday <= cdates(1) .or. calday > cdates(12)) then
       nm=12
       np=1
    else
       nm=0
       do n=1,ptrtim-1
          if(calday>cdates(n) .and. calday<=cdates(n+1)) then
             nm=n
             np=n+1
          end if
       end do
       if(nm .eq. 0) then
          if ( is_perpetual() ) then
             call get_perp_date(yr, mon, day, ncsec)
          else
             call get_curr_date(yr, mon, day, ncsec)
          end if
          ncdate = yr*10000 + mon*100 + day

          write(iulog,*)'model date:', ncdate, ncsec,'boundary data dates:', cdates
          call endrun('BOUNDARYDATA_READ: Failed to find dates bracketing dates')
       end if
    end if

  end subroutine get_data_bounding_date_indices


  !================================================================================================
  subroutine boundarydata_vert_interp(lchnk, ncol, levsiz, fldcnt, pin, pmid, datain, dataout)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Interpolate ozone from current time-interpolated values to model levels
    ! 
    ! Method: Use pressure values to determine interpolation levels
    ! 
    ! Author: Bruce Briegleb
    ! 
    !--------------------------------------------------------------------------
    implicit none
    ! Arguments
    !
    integer, intent(in) :: lchnk               ! chunk identifier
    integer, intent(in) :: ncol                ! number of atmospheric columns
    integer, intent(in) :: levsiz
    integer, intent(in) :: fldcnt
    real(r8), intent(in) :: pin(levsiz)
    real(r8), intent(in) :: pmid(pcols,pver)   ! level pressures (mks)
    real(r8), intent(in) :: datain(pcols,levsiz,fldcnt)
    real(r8), intent(out) :: dataout(pcols,pver,fldcnt)    ! ozone mass mixing ratio
    !
    ! local storage
    !

    integer ::  i                   ! longitude index
    integer ::  k, kk, kkstart      ! level indices
    integer ::  kupper(pcols)       ! Level indices for interpolation
    integer ::  kount               ! Counter
    integer ::  fld
    real(r8) dpu                ! upper level pressure difference
    real(r8) dpl                ! lower level pressure difference
    !--------------------------------------------------------------------------
    !
    ! Initialize index array
    !
    do i=1,ncol
       kupper(i) = 1
    end do

    do k=1,pver
       !
       ! Top level we need to start looking is the top level for the previous k
       ! for all longitude points
       !
       kkstart = levsiz
       do i=1,ncol
          kkstart = min0(kkstart,kupper(i))
       end do
       kount = 0
       !
       ! Store level indices for interpolation
       !
       do kk=kkstart,levsiz-1
          do i=1,ncol
             if (pin(kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(kk+1)) then
                kupper(i) = kk
                kount = kount + 1
             end if
          end do
          !
          ! If all indices for this level have been found, do the interpolation and
          ! go to the next level
          !
          if (kount.eq.ncol) then
             do fld=1,fldcnt
                do i=1,ncol
                   dpu = pmid(i,k) - pin(kupper(i))
                   dpl = pin(kupper(i)+1) - pmid(i,k)
                   dataout(i,k,fld) = (datain(i,kupper(i),fld )*dpl + &
                        datain(i,kupper(i)+1,fld)*dpu)/(dpl + dpu)

                end do
             end do
             goto 35
          end if
       end do
       !
       ! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
       ! must extrapolate from the bottom or top data level for at least some
       ! of the longitude points.
       !
       do fld=1,fldcnt
          do i=1,ncol
             if (pmid(i,k) .lt. pin(1)) then
                dataout(i,k,fld) = datain(i,1,fld)*pmid(i,k)/pin(1)
             else if (pmid(i,k) .gt. pin(levsiz)) then
                dataout(i,k,fld) = datain(i,levsiz,fld)
             else
                dpu = pmid(i,k) - pin(kupper(i))
                dpl = pin(kupper(i)+1) - pmid(i,k)
                dataout(i,k,fld) = (datain(i,kupper(i),fld )*dpl + &
                     datain(i,kupper(i)+1,fld)*dpu)/(dpl + dpu)
             end if
          end do
       end do
       if (kount.gt.ncol) then
          call endrun ('ozone_data_vert_interp: Bad ozone data: non-monotonicity suspected')
       end if
35     continue
    end do
  end subroutine boundarydata_vert_interp
#if 0  
  subroutine ncol_read_bracket(cid,columnmap,start,count,ncol)
    integer, intent(in) :: cid(:), columnmap(:), ncol
    integer, intent(out) :: start(:), count(:)

    integer :: i, j, tcol
    
    tcol = size(columnmap)
    count=1
    do i=1,ncol
#if 1
	start(i)=cid(i)
	count(i)=1
#else
       do j=1,tcol
          if(columnmap(j).eq.cid(i)) then
             start(i)=j
             count(i)=1
             exit
          end if
       end do
#endif
    end do
    do i=1,ncol-1
       do j=1,ncol-i
          if(columnmap(start(i+j)).eq.columnmap(start(i)+j)) then
             count(i)=count(i)+1
          else
             exit
          end if
       end do
    end do
    write(iulog,*) __LINE__,cid(1),cid(ncol),minval(cid(1:ncol)),maxval(cid(1:ncol))

  end subroutine ncol_read_bracket
#endif
end module boundarydata
