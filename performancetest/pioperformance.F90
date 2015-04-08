program pioperformance
#ifndef NO_MPIMOD
  use mpi
#endif
  use perf_mod, only : t_initf, t_finalizef
  use pio, only : pio_iotype_netcdf, pio_iotype_pnetcdf, pio_iotype_netcdf4p, &
       pio_iotype_netcdf4c, pio_rearr_subset, pio_rearr_box
  implicit none
#ifdef NO_MPIMOD
#include <mpif.h>
#endif  
  integer, parameter :: max_io_task_array_size=64, max_decomp_files=64


  integer :: ierr, mype, npe, i
  logical :: Mastertask
  character(len=256) :: decompfile(max_decomp_files)
  character(len=8) :: pio_typenames(4)
  integer :: piotypes(4), niotypes
  integer :: rearrangers(2)
  integer, parameter :: max_nvars=12
  integer :: niotasks(max_io_task_array_size)
  integer :: nv, nframes, nvars(max_nvars)
  namelist /pioperf/ decompfile, pio_typenames, rearrangers, niotasks, nframes, nvars
#ifdef BGQTRY
  external :: print_memusage
#endif
  !
  ! Initialize MPI
  !
  call MPI_Init(ierr)
  call CheckMPIreturn(__LINE__,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)
  call CheckMPIreturn(__LINE__,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, npe,  ierr)
  call CheckMPIreturn(__LINE__,ierr)
  if(mype==0) then
     Mastertask=.true.
  else
     Mastertask=.false.
  endif
#ifdef BGQTRY
  call print_memusage()
#endif
  nvars = 0
  niotasks = -1 ! loop over all possible values
  rearrangers = 0
  nframes = 5
  decompfile = ' '
  pio_typenames = ' '
  piotypes = -1
  if(mype==0) then
     open(unit=12,file='pioperf.nl',status='old')
     read(12,pioperf)
     close(12)

     do i=1,4
        if(pio_typenames(i) .eq. 'netcdf') then
           piotypes(i) = PIO_IOTYPE_NETCDF
        else if(pio_typenames(i) .eq. 'netcdf4p') then
           piotypes(i) = PIO_IOTYPE_NETCDF4P
        else if(pio_typenames(i) .eq. 'netcdf4c') then
           piotypes(i) = PIO_IOTYPE_NETCDF4C
        else if(pio_typenames(i) .eq. 'pnetcdf') then
           piotypes(i) = PIO_IOTYPE_PNETCDF
        else
           exit
        endif
     enddo

  endif

  call MPI_Bcast(decompfile,256*max_decomp_files,MPI_CHARACTER,0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(piotypes,4, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(rearrangers, 2, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(niotasks, max_io_task_array_size, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nframes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nvars, max_nvars, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

  call t_initf('pioperf.nl', LogPrint=.false., mpicom=MPI_COMM_WORLD, MasterTask=MasterTask)
  niotypes = 0
  do i=1,4
     if (piotypes(i) > -1) niotypes = niotypes+1
  enddo
  if(rearrangers(1)==0) then
    rearrangers(1)=1
    rearrangers(2)=2
  endif  

  do i=1,max_decomp_files
     if(len_trim(decompfile(i))==0) exit
     if(mype == 0) print *, decompfile(i)
     do nv=1,max_nvars
       if(nvars(nv)>0) then
         call pioperformancetest(decompfile(i), piotypes(1:niotypes), mype, npe, rearrangers, niotasks, nframes, nvars(nv))
     endif
   enddo
  enddo
  call t_finalizef()

  call MPI_Finalize(ierr)
contains

  subroutine pioperformancetest(filename, piotypes, mype, npe_base, rearrangers, niotasks,nframes, nvars)
    use pio
    use pio_support, only : pio_readdof
    use perf_mod
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mype, npe_base
    integer, intent(in) :: piotypes(:)
    integer, intent(in) :: rearrangers(:)
    integer, intent(inout) :: niotasks(:)
    integer, intent(in) :: nframes 
    integer, intent(in) :: nvars

    integer(kind=PIO_Offset_kind), pointer :: compmap(:)
    integer :: ntasks
    integer :: comm
    integer :: npe
    integer :: color
    integer(kind=PIO_Offset_kind) :: maplen, gmaplen
    integer :: ndims
    integer, pointer :: gdims(:)
    character(len=20) :: fname
    type(var_desc_t) :: vari(nvars), varr(nvars), vard(nvars)
    type(iosystem_desc_t) :: iosystem
    integer :: stride, n
    integer, allocatable :: ifld(:,:), ifld_in(:,:)
    real, allocatable :: rfld(:)
    double precision, allocatable :: dfld(:,:), dfld_in(:,:)
    type(file_desc_t) :: File
    type(io_desc_t) :: iodesc_i4, iodesc_r4, iodesc_r8
    integer :: ierr
    integer(pio_offset_kind) :: frame=1
    integer :: iotype, rearr, rearrtype
    integer :: j, k, errorcnt
    character(len=8) :: varname
    double precision :: wall(2), sys(2), usr(2)
    integer :: niomin, niomax
    integer :: nv, mode
    integer,  parameter :: c0 = -1
    double precision, parameter :: cd0 = 1.0e30
    nullify(compmap)



    call pio_readdof(filename, ndims, gdims, compmap, MPI_COMM_WORLD)
    maplen = size(compmap)

!    color = 0
!    if(maplen>0) then
       color = 1
!    endif

    call MPI_Comm_split(MPI_COMM_WORLD, color, mype, comm, ierr)

    call MPI_Comm_size(comm, npe,  ierr)
    call CheckMPIreturn(__LINE__,ierr)
    niomin=1
    niomax=min(npe,max_io_task_array_size)
    if(niotasks(1)<=0) then
       do j=1,min(max_io_task_array_size, npe)
          niotasks(j)=npe-j+1
       enddo
    endif

    if(mype < npe) then

       call MPI_ALLREDUCE(maplen,gmaplen,1,MPI_INTEGER8,MPI_SUM,comm,ierr)

!       if(gmaplen /= product(gdims)) then
!          print *,__FILE__,__LINE__,gmaplen,gdims
!       endif
    
       allocate(ifld(maplen,nvars))
       allocate(ifld_in(maplen,nvars))
       allocate(rfld(maplen))
       allocate(dfld(maplen,nvars))
       allocate(dfld_in(maplen,nvars))

!       ifld = mype
       rfld = mype
       dfld = mype
       do nv=1,nvars
          do j=1,maplen
!             ifld(j,nv) = nv*100000000 +mype*1000000 + compmap(j)
             ifld(j,nv) = compmap(j)
             dfld(j,nv) = ifld(j,nv)/1000000.0
!             if(nv==2) then
!                ifld(j+(nv-1)*maplen) = -(mype*1000000 + compmap(j))
!             else
!                ifld(j+(nv-1)*maplen) = mype*1000000 + compmap(j)
!             endif
          enddo
       enddo

#ifdef BGQTRY
  call print_memusage()
#endif

       do k=1,size(piotypes)
          iotype = piotypes(k)
          call MPI_Barrier(comm,ierr)
          if(mype==0) then
             print *,'iotype=',piotypes(k)
          endif
          if(iotype==PIO_IOTYPE_PNETCDF) then
             mode = PIO_64BIT_DATA
          else
             mode = 0
          endif
          do rearrtype=1,2
             rearr = rearrangers(rearrtype)
             if(rearr /= PIO_REARR_SUBSET .and. rearr /= PIO_REARR_BOX) exit

             do n=niomin,niomax
                ntasks = niotasks(n)

                if(ntasks<=0 .or. ntasks>npe) exit
                stride = max(1,npe/ntasks)

                call pio_init(mype, comm, ntasks, 0, stride, PIO_REARR_SUBSET, iosystem)
                   
                write(fname, '(a,i1,a,i4.4,a,i1,a)') 'pioperf.',rearr,'-',ntasks,'-',iotype,'.nc'
		
                ierr =  PIO_CreateFile(iosystem, File, iotype, trim(fname), mode)

                call WriteMetadata(File, gdims, vari, varr, vard)
                call MPI_Barrier(comm,ierr)
                call t_stampf(wall(1), usr(1), sys(1))
                call PIO_InitDecomp(iosystem, PIO_INT, gdims, compmap, iodesc_i4, rearr=rearr)
!                call PIO_InitDecomp(iosystem, PIO_REAL, gdims, compmap, iodesc_r4, rearr=rearr)
!                call PIO_InitDecomp(iosystem, PIO_DOUBLE, gdims, compmap, iodesc_r8, rearr=rearr)
!                print *,__FILE__,__LINE__,minval(dfld),maxval(dfld),minloc(dfld),maxloc(dfld)

                do frame=1,nframes
!                   if(mype==0) print *,__FILE__,__LINE__,frame
                   do nv=1,nvars   
                      call PIO_setframe(File, vari(nv), frame)
                      call pio_write_darray(File, vari(nv), iodesc_i4, ifld(:,nv)    , ierr, fillval= -nv)
                   enddo
! multiversion  
!                 call pio_write_darray(File, vari, iodesc_i4, ifld, ierr)

!!$
!!$                   
!!$                   call PIO_setframe(File, varr, frame)
!!$                   call pio_write_darray(File, varr, iodesc_r4, rfld, ierr)
!!$                   
!!$                   
!!$                   call PIO_setframe(File, vard, frame)
!!$                   
!!$                   call pio_write_darray(File, vard, iodesc_r8, dfld, ierr)
!!$
!!$
                enddo
                call pio_closefile(File)
!                print *,__FILE__,__LINE__

                call MPI_Barrier(comm,ierr)

                call t_stampf(wall(2), usr(2), sys(2))
                wall(1) = wall(2)-wall(1)
                call MPI_Reduce(wall(1), wall(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
                if(mype==0) then
                   ! print out performance in MB/s
                   print *, 'write ',rearr, ntasks, nvars, nvars*nframes*gmaplen*4.0/(1048576.0*wall(2))
#ifdef BGQTRY
  call print_memusage()
#endif
                end if
! Now the Read
                ierr = PIO_OpenFile(iosystem, File, iotype, trim(fname), mode=PIO_NOWRITE);
                do nv=1,nvars
                   write(varname,'(a,i4.4)') 'vari',nv
                   ierr =  pio_inq_varid(File, varname, vari(nv))
                enddo
                call MPI_Barrier(comm,ierr)
                call t_stampf(wall(1), usr(1), sys(1))
                
                do frame=1,nframes                   
                   do nv=1,nvars
                      call PIO_setframe(File, vari(nv), frame)
                      call pio_read_darray(File, vari(nv), iodesc_i4, ifld_in(:,nv), ierr)
                   enddo
                enddo
                
                call pio_closefile(File)
                call MPI_Barrier(comm,ierr)
                call t_stampf(wall(2), usr(2), sys(2))
                wall(1) = wall(2)-wall(1)
                call MPI_Reduce(wall(1), wall(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
                errorcnt = 0
                do nv=1,nvars
                   do j=1,maplen
                      if(ifld(j,nv) /= ifld_in(j,nv) .and. compmap(j) /= 0) then
                         if(errorcnt < 10) then
                            print *,__LINE__,mype,j,nv,ifld(j,nv),ifld_in(j,nv),compmap(j)
                         endif
                         errorcnt = errorcnt+1
                      endif
                   enddo
                enddo
                j = errorcnt
                call MPI_Reduce(j, errorcnt, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
                
                if(mype==0) then
                   if(errorcnt > 0) then
                      print *,'ERROR: INPUT/OUTPUT data mismatch ',errorcnt
                   endif
                   print *, 'read ',rearr, ntasks,nvars, nvars*nframes*gmaplen*4.0/(1048576.0*wall(2))
#ifdef BGQTRY 
  call print_memusage()
#endif
                end if
                
                call PIO_freedecomp(iosystem, iodesc_r4)
                call PIO_freedecomp(iosystem, iodesc_r8)
                call PIO_freedecomp(iosystem, iodesc_i4)
                
                call pio_finalize(iosystem, ierr)
             enddo
          enddo
       enddo
!       deallocate(compmap)
       deallocate(ifld)
       deallocate(ifld_in)
       deallocate(rfld)
       deallocate(dfld)
       deallocate(dfld_in)
    endif

  end subroutine pioperformancetest

  subroutine WriteMetadata(File, gdims, vari, varr, vard)
    use pio
    type(file_desc_t) :: File
    integer, intent(in) :: gdims(:)
    type(var_desc_t),intent(out) :: vari(:), varr(:), vard(:)
    integer :: ndims
    character(len=12) :: dimname
    character(len=8) :: varname
    integer, allocatable :: dimid(:)
    integer :: i, iostat, nv
    integer :: nvars

    nvars = size(vari)

    ndims = size(gdims)
    
    allocate(dimid(ndims+1))

    do i=1,ndims
       write(dimname,'(a,i6.6)') 'dim',i  
       iostat = PIO_def_dim(File, trim(dimname), int(gdims(i),pio_offset_kind), dimid(i))
    enddo
    iostat = PIO_def_dim(File, 'time', PIO_UNLIMITED, dimid(ndims+1))

    do nv=1,nvars
       write(varname,'(a,i4.4)') 'vari',nv
       iostat = PIO_def_var(File, varname, PIO_INT, dimid, vari(nv))
       write(varname,'(a,i4.4)') 'varr',nv
       iostat = PIO_def_var(File, varname, PIO_REAL, dimid, varr(nv))
       write(varname,'(a,i4.4)') 'vard',nv
       iostat = PIO_def_var(File, varname, PIO_DOUBLE, dimid, vard(nv))
    enddo

    iostat = PIO_enddef(File)

  end subroutine WriteMetadata


!=============================================
!  CheckMPIreturn:
!
!      Check and prints an error message
!  if an error occured in a MPI subroutine.
!=============================================
  subroutine CheckMPIreturn(line,errcode)
    use MPI
    implicit none
    integer, intent(in) :: errcode
    integer, intent(in) :: line
    character(len=MPI_MAX_ERROR_STRING) :: errorstring
    
    integer :: errorlen
    
    integer :: ierr
    
    if (errcode .ne. MPI_SUCCESS) then
       call MPI_Error_String(errcode,errorstring,errorlen,ierr)
       write(*,*) errorstring(1:errorlen)
    end if
  end subroutine CheckMPIreturn

end program pioperformance
