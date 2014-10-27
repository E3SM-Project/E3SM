program pioperformance
  use mpi
  use perf_mod, only : t_initf, t_finalizef
  implicit none
  integer :: ierr, mype, npe
  character(len=30) :: decompfile

!  decompfile = "piodecomp30tasks01dims06.dat"
!  decompfile = "piodecomp30tasks02dims05.dat"
  decompfile = "piodecomp04tasks02dims03.dat"
!  decompfile = "simple.dat"
  !
  ! Initialize MPI
  !
  call MPI_Init(ierr)
  call CheckMPIreturn(__LINE__,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)
  call CheckMPIreturn(__LINE__,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, npe,  ierr)
  call CheckMPIreturn(__LINE__,ierr)

  call t_initf('pioperf.nl', LogPrint=.false.)

  call pioperformancetest(decompfile, mype, npe)

  call t_finalizef()

  call MPI_Finalize(ierr)
contains

  subroutine pioperformancetest(filename, mype, npe_base)
    use pio
    use pio_support, only : pio_readdof
    use perf_mod
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mype, npe_base
    integer(kind=PIO_Offset_kind), pointer :: compmap(:)
    integer :: ntasks
    integer :: comm
    integer :: npe
    integer :: color
    integer(kind=PIO_Offset_kind) :: maplen, gmaplen
    integer :: ndims
    integer, pointer :: gdims(:)
    character(len=20) :: fname
    type(var_desc_t) :: vari, varr, vard
    type(iosystem_desc_t) :: iosystem
    integer :: stride, niotasks
    integer, allocatable :: ifld(:)
    real, allocatable :: rfld(:)
    double precision, allocatable :: dfld(:)
    type(file_desc_t) :: File
    type(io_desc_t) :: iodesc_i4, iodesc_r4, iodesc_r8
    integer :: ierr
    integer(pio_offset_kind) :: frame=1
    integer :: iotype, rearr, rearrtype
    integer :: j
    integer :: nframes = 1
    double precision :: wall(2), sys(2), usr(2)
    

    nullify(compmap)

    call pio_readdof(filename, ndims, gdims, compmap, MPI_COMM_WORLD)
    maplen = size(compmap)

    color = 0
    if(maplen>0) then
       color = 1
    endif

    call MPI_Comm_split(MPI_COMM_WORLD, color, mype, comm, ierr)

    call MPI_Comm_size(comm, npe,  ierr)
    call CheckMPIreturn(__LINE__,ierr)


    if(mype < npe) then

       call MPI_ALLREDUCE(maplen,gmaplen,1,MPI_OFFSET,MPI_SUM,comm,ierr)

       if(gmaplen /= product(gdims)) then
          print *,__FILE__,__LINE__,gmaplen,gdims
       endif
    
       allocate(ifld(maplen))
       allocate(rfld(maplen))
       allocate(dfld(maplen))

!       ifld = mype
       rfld = mype
       dfld = mype
       iotype = PIO_IOTYPE_PNETCDF

       do j=1,maplen
          ifld(j) = mype*1000000 + compmap(j)
       enddo
       rearr = PIO_REARR_SUBSET
       do rearrtype=1,2

          do niotasks=1,npe
             stride = max(1,npe/niotasks)
             if(real(stride) == real(npe)/real(niotasks)) then
                call pio_init(mype, comm, niotasks, 0, stride, PIO_REARR_SUBSET, iosystem)
                
                write(fname, '(a,i1,a,i4.4,a)') 'pioperf.',rearrtype,'-',niotasks,'.nc'
                ierr =  PIO_CreateFile(iosystem, File, iotype, trim(fname))
                
                call WriteMetadata(File, gdims, vari, varr, vard)
                call t_stampf(wall(1), usr(1), sys(1))
                call PIO_InitDecomp(iosystem, PIO_INT, gdims, compmap, iodesc_i4, rearr=rearr)
                call PIO_InitDecomp(iosystem, PIO_REAL, gdims, compmap, iodesc_r4, rearr=rearr)
                call PIO_InitDecomp(iosystem, PIO_DOUBLE, gdims, compmap, iodesc_r8, rearr=rearr)

                do frame=1,nframes
                   
                   call PIO_setframe(File, vari, frame)
                   print *,__FILE__,__LINE__
                   call pio_write_darray(File, vari, iodesc_i4, ifld, ierr)
                   print *,__FILE__,__LINE__
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
                call PIO_freedecomp(iosystem, iodesc_r4)
                call PIO_freedecomp(iosystem, iodesc_r8)
                
                call PIO_freedecomp(iosystem, iodesc_i4)
                call pio_closefile(File)
                call t_stampf(wall(2), usr(2), sys(2))
                wall(1) = wall(2)-wall(1)
                call MPI_Reduce(wall(1), wall(2), 1, MPI_DOUBLE, MPI_MAX, 0, comm, ierr)
                if(mype==0) then
                   print *, rearrtype, niotasks, wall(2)
                end if
             
                call pio_finalize(iosystem, ierr)
             endif
          enddo
          rearr = PIO_REARR_BOX
       enddo

       deallocate(compmap)
       deallocate(ifld)
    endif

  end subroutine pioperformancetest

  subroutine WriteMetadata(File, gdims, vari, varr, vard)
    use pio
    type(file_desc_t) :: File
    integer, intent(in) :: gdims(:)
    type(var_desc_t),intent(out) :: vari, varr, vard
    integer :: ndims
    character(len=12) :: dimname
    integer, allocatable :: dimid(:)
    integer :: i, iostat

    ndims = size(gdims)
    
    allocate(dimid(ndims+1))

    do i=1,ndims
       write(dimname,'(a,i6.6)') 'dim',i  
       iostat = PIO_def_dim(File, trim(dimname), int(gdims(i),pio_offset_kind), dimid(i))
    enddo
    iostat = PIO_def_dim(File, 'time', PIO_UNLIMITED, dimid(ndims+1))

    iostat = PIO_def_var(File, 'vari', PIO_INT, dimid, vari)
    
    iostat = PIO_def_var(File, 'varr', PIO_REAL, dimid, varr)
    
    iostat = PIO_def_var(File, 'vard', PIO_DOUBLE, dimid, vard)

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
