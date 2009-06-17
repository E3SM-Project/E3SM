#define _FILE_ "piolib_mod.F90"
! Node of interest for TESTMEM (-1 means all)
#define NODE 0

#define DEBUG_REARR 0
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module piolib_mod

  !BOP
  ! !MODULE: piolib_mod
  ! !DESCRIPTION:
  !  This module provides a parallel input/output interface
  !  for writing fields.
  !
  ! !REVISION HISTORY:
  !  CVS:$Id: piolib_mod.F90,v 1.25 2007/02/06 21:37:14 dennis Exp $
  !  CVS:$Name:  $

  ! !USES:

  !--------------
  use pio_kinds
  !--------------
  use pio_types
  !--------------
  use alloc_mod
  !--------------
  use pio_support, only : piodie, Debug, DebugIO, CheckMPIReturn
  !
  use nf_mod, only : create_nf, open_nf,close_nf, sync_nf
  use pionfread_mod, only : read_nf
  use pionfwrite_mod, only : write_nf

  use iompi_mod
  use rearrange
#ifdef TIMING
  use perf_mod, only : t_startf, t_stopf     ! _EXTERNAL
#endif

  implicit none
  private
  save

  include 'mpif.h'    ! _EXTERNAL

  ! !PUBLIC MEMBER FUNCTIONS:

  public ::  PIO_init,      &
       PIO_initDecomp,    &
       PIO_OpenFile,      &
       PIO_SyncFile,      &
       PIO_CreateFile,    &
       PIO_CloseFile,     &
       PIO_setIOtype,     &
       PIO_NumToRead,     &
       PIO_NumToWrite,    &
       PIO_setframe,      &
       PIO_AdvanceFrame,  &
       PIO_SetDebugLevel, &
       PIO_SetErrorHandling, &
       PIO_get_local_array_size, &
       PIO_FreeDecomp

#ifdef DOTHIS
  !
  !  Attribute functions 
  !
  public ::  PIO_inquire, &
       PIO_put_var,       &
       PIO_get_var,       &
       PIO_ReDef,         &
       PIO_inq_attname, &
       PIO_inq_varid, &
       PIO_inq_varname, &
       PIO_inq_vartype, &
       PIO_inq_varndims, &
       PIO_inq_vardimid, &
       PIO_inq_varnatts, &
       PIO_inq_att,   &
       PIO_inq_attlen,&
       PIO_inq_dimid, &
       PIO_inq_dimname, &
       PIO_inq_dimlen, &
       PIO_put_att,   &
       PIO_get_att,   &
       PIO_copy_att

  !
  !  Define functions
  !
  public :: PIO_def_dim, &
       PIO_def_var, &
       PIO_enddef
#endif
  public :: PIO_setNUMAgg,     &
       PIO_dupIODesc
#ifdef MEMCHK
  integer :: lastrss=0
#endif

  !EOP
  !BOC
  !-----------------------------------------------------------------------
  !
  !  module variables
  !
  !-----------------------------------------------------------------------

  interface PIO_OpenFile
     module procedure OpenFile
  end interface

  interface PIO_SyncFile
     module procedure SyncFile
  end interface

  interface PIO_CreateFile
     module procedure CreateFile
  end interface

  interface PIO_SetFrame
     module procedure SetFrame
  end interface

  interface PIO_AdvanceFrame
     module procedure AdvanceFrame
  end interface

  interface PIO_CloseFile
     module procedure CloseFile
  end interface

  interface PIO_FreeDecomp
     module procedure FreeDecomp_ios
     module procedure FreeDecomp_file
  end interface

  interface PIO_init
     module procedure init
  end interface

  interface PIO_initDecomp
     module procedure initDecomp_1dof_nf, &
          initDecomp_1dof_nf_box, &
          initDecomp_1dof_bin,& 
          initDecomp_2dof_nf, &
          initDecomp_2dof_bin
  end interface

  interface PIO_dupIOdesc
     module procedure dupIOdesc
  end interface

  interface PIO_setIOtype 
     module procedure setIOtype
  end interface

  interface PIO_NumToRead
     module procedure numToRead
  end interface

  interface PIO_NumToWrite
     module procedure numToWrite
  end interface

  interface PIO_setNUMAgg
     module procedure setNUMAgg
  end interface

  interface PIO_SetDebugLevel
     module procedure SetDebugLevel
  end interface

  interface PIO_SetErrorHandling
     module procedure SetErrorHandlingF
     module procedure SetErrorHandlingI
  end interface


  !EOC
  !***********************************************************************


contains
  integer function pio_get_local_array_size(iodesc)
    type(io_desc_t), intent(in) :: iodesc
    pio_get_local_array_size = iodesc%compsize
  end function pio_get_local_array_size


  !************************************
  ! AdvanceFrame
  !

  subroutine AdvanceFrame(vardesc)
    type(Var_Desc_t), intent(inout) :: vardesc
    vardesc%rec=vardesc%rec+1
  end subroutine AdvanceFrame



  !************************************
  ! SetFrame
  !

  subroutine SetFrame(vardesc,frame)
    type(Var_Desc_t), intent(inout) :: vardesc
    integer(kind=PIO_OffSet), intent(in) :: frame
    vardesc%rec=frame
  end subroutine SetFrame



  !************************************
  ! SetDebugLevel
  !

  subroutine SetDebugLevel(level)
    integer(i4), intent(in) :: level	
    if(level.eq.0) then
       Debug=.FALSE.
       DebugIO=.FALSE.
    else if(level.eq.1) then
       Debug=.TRUE.
       DebugIO=.FALSE.
    else if(level.eq.2) then
       Debug=.FALSE.
       DebugIO=.TRUE.
    else if(level.ge.3) then
       Debug=.TRUE.
       DebugIO=.TRUE.
    end if
  end subroutine SetDebugLevel


  subroutine SetErrorHandlingF(File, method)
    type(file_desc_t), intent(inout) :: File
    integer, intent(in) :: method

    call SetErrorHandlingI(File%iosystem, method)
  end subroutine SetErrorHandlingF

  subroutine SetErrorHandlingI(IOsystem, method)
    type(iosystem_desc_t), intent(inout) :: IOsystem
    integer, intent(in) :: method

    iosystem%error_handling=method

    if(method > PIO_INTERNAL_ERROR .or. method < PIO_RETURN_ERROR) then
       call piodie(_FILE_,__LINE__,'invalid error handling method requested')
    end if
  end subroutine SetErrorHandlingI

  !************************************
  ! initDecomp_2dof_BIN
  !************************************
  subroutine initDecomp_2dof_BIN(Iosystem,basepioTYPE,dims,lenBLOCKS,compDOF,ioDofr,iodofw,IOdesc)
    type (Iosystem_desc_t), intent(in) :: Iosystem
    integer(i4), intent(in)           :: basepioTYPE
    integer(i4)                       :: baseTYPE
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenBLOCKS
    integer (i4), intent(in)          :: compDOF(:)   ! Global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: ioDofR(:)     ! Global degrees of freedom for IO decomposition 
    integer (i4), intent(in)          :: ioDofW(:)     ! Global degrees of freedom for IO decomposition 
    type (IO_desc_t), intent(out)     :: IOdesc

    integer(kind=PIO_Offset) :: start(1), count(1)

    integer (i4) :: i,ndims,glength,n_iotasks
    logical :: UseRearranger
    integer (i4) ::  nDispR,nDispW
    integer (i4) :: lengthR, lengthW
    integer (i4), allocatable :: displaceR(:),displaceW(:)

    baseTYPE=pio_type_to_mpi_type(basepioTYPE)

    !-------------------------------------------
    ! For testing purposes Set the IOMap
    ! (DecompMap_t) to something basic for
    ! testing.
    !-------------------------------------------
    UseRearranger = Iosystem%UseRearranger

    !---------------------
    ! number of dimensions
    !---------------------
    ndims = SIZE(dims)
    !---------------------
    ! total global size
    !---------------------
    gLength=1
    do i=1,ndims
       gLength = gLength*dims(i)
    enddo
    lengthR = SIZE(ioDofR);
    lengthW = SIZE(ioDofW)
    nDispW=SIZE(ioDofW)/lenBLOCKS; nDispR=SIZE(ioDofR)/lenBLOCKS

    allocate(displaceR(nDispR),displaceW(nDispW))

    !--------------------------------------------
    ! calculate MPI data structure displacements
    !--------------------------------------------
    !DBG    print *,'PIO_initDecomp: before call to CalcDisplace'
    call CalcDisplace(lenBLOCKS,ioDofR,displaceR)
    call CalcDisplace(lenBLOCKS,ioDofW,displaceW)

    n_iotasks = Iosystem%num_iotasks

    IOdesc%glen = gLength

    if(Debug) print *,'IAM: ',Iosystem%io_rank,'initDecomp: UseRearranger: ',UseRearranger

    !---------------------------------------------
    !  The setup for the MPI-IO type information
    !---------------------------------------------
    if(Iosystem%IOproc) then
       !-----------------------------------------------
       ! setup the data structure for the read operation
       !-----------------------------------------------
       IODesc%Read%n_elemTYPE = SIZE(displaceR)
       IODesc%Read%n_words    = IODesc%Read%n_elemTYPE*lenBLOCKS
       call GenIndexedBlock(lenBLOCKS,baseTYPE,IODesc%Read%elemTYPE,IODesc%Read%fileTYPE,displaceR)

       !-------------------------------------------------
       ! setup the data structure for the write operation
       !-------------------------------------------------
       IODesc%Write%n_elemTYPE = SIZE(displaceW)
       IODesc%Write%n_words    = IODesc%Write%n_elemTYPE*lenBLOCKS
       call GenIndexedBlock(lenBLOCKS,baseTYPE,IODesc%Write%elemTYPE,IODesc%Write%fileTYPE,displaceW)

       if(Debug) print *,'initDecomp: At the end of subroutine'
       !       if(IODesc%Read%n_elemTYPE == 0 .and. IODesc%Write%n_elemTYPE == 0) IoSystem%IOproc = .FALSE.
    endif

    deallocate(displaceR,displaceW)


  end subroutine initDecomp_2dof_BIN


  !************************************
  ! initDecomp_1dof_BIN
  !

  subroutine initDecomp_1dof_BIN(Iosystem,baseTYPE,dims,lenBLOCKS,compDOF,ioDofr,IOdesc)
    type (Iosystem_desc_t), intent(in) :: Iosystem
    integer(i4), intent(in)           :: baseTYPE
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenBLOCKS
    integer (i4), intent(in)          :: compDOF(:)   ! Global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: ioDofR(:)     ! Global degrees of freedom for IO decomposition 
    type (IO_desc_t), intent(out)     :: IOdesc

    integer(kind=PIO_Offset) :: start(1), count(1)
    ! these are not used in the binary interface

    start(1)=-1
    count(1)=-1
    call initDecomp_1dof_nf(Iosystem,baseTYPE,dims,lenBLOCKS,compDOF,ioDofr,start, count, IOdesc)
  end subroutine initDecomp_1dof_BIN



  !************************************
  ! initDecomp_2dof_nf
  !

  subroutine initDecomp_2dof_nf(Iosystem,baseTYPE,dims,lenBLOCKS,compDOF,ioDofr,iodofw,start, count, IOdesc)
    type (Iosystem_desc_t), intent(in) :: Iosystem
    integer(i4), intent(in)           :: baseTYPE
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenBLOCKS
    integer (i4), intent(in)          :: compDOF(:)   ! Global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: ioDofR(:)     ! Global degrees of freedom for IO decomposition 
    integer (i4), intent(in)          :: ioDofW(:)     ! Global degrees of freedom for IO decomposition 

    type (IO_desc_t), intent(out)     :: IOdesc

    integer(kind=PIO_Offset), intent(in) :: start(:), count(:)
    type (io_desc_t) :: tmp


    call initDecomp_1dof_nf(Iosystem, baseType, dims, lenblocks, compdof, iodofR, start, count, iodesc)

    call initDecomp_1dof_nf(Iosystem, baseType, dims, lenblocks, compdof, iodofW, start, count, tmp)

    call dupiodesc2(iodesc%write,tmp%write)

    if(Debug) then
       print *, _FILE_,__LINE__,iodesc%read%filetype,iodesc%read%elemtype,&
            iodesc%read%n_elemtype,iodesc%read%n_words   
       print *, _FILE_,__LINE__,iodesc%write%filetype,iodesc%write%elemtype,&
            iodesc%write%n_elemtype,iodesc%write%n_words
    end if

  end subroutine initDecomp_2dof_nf



  !************************************
  ! initDecomp_1dof_nf
  !

  subroutine initDecomp_1dof_nf(Iosystem,basepioTYPE,dims,lenblocks,compDOF,ioDof,start, count, IOdesc)
    type (Iosystem_desc_t), intent(in) :: Iosystem
    integer(i4), intent(in)           :: basepioTYPE
    integer(i4)                       :: baseTYPE
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in) :: lenBLOCKS
    integer (i4), intent(in)          :: compDOF(:)   ! Global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: ioDof(:)     ! Global degrees of freedom for IO decomposition 
    type (IO_desc_t), intent(out)     :: IOdesc

    integer(kind=PIO_Offset), intent(in) :: start(:), count(:)

    integer(i4) :: length,n_iotasks
    integer(i4) :: ndims

    integer (i4), pointer :: displace(:)  ! The displacements for the MPI data structure (Read)

    integer(i4) :: prev
    integer(i4) :: gLength    ! global length in words
    integer(i4) :: ii,i,dis,ierr
    integer(i4),pointer, dimension(:) :: blocklen,disp
    logical(log_kind) ::  UseRearranger
    logical, parameter :: Check = .TRUE.
    integer(i4) :: nDisp
#ifdef MEMCHK
    integer :: msize, rss, mshare, mtext, mstack
#endif

    baseTYPE=pio_type_to_mpi_type(basepioTYPE)

    !-------------------------------------------
    ! For testing purposes Set the IOMap
    ! (DecompMap_t) to something basic for
    ! testing.
    !-------------------------------------------
#ifdef TIMING
    call t_startf("pio_initDecomp")
#endif
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
    UseRearranger = Iosystem%UseRearranger
    !---------------------
    ! number of dimensions
    !---------------------
    ndims = SIZE(dims)
    !---------------------
    ! total global size
    !---------------------
    gLength=1
    do i=1,ndims
       gLength = gLength*dims(i)
    enddo

    nDisp=SIZE(ioDof)/lenBLOCKS

    call alloc_check(displace,nDisp)

#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif

    call alloc_check(IODesc%start,ndims)
    call alloc_check(IODesc%count,ndims)
    IODesc%start(:) = start(:)
    IODesc%count(:) = count(:)
    !--------------------------------------------
    ! calculate MPI data structure displacements 
    !--------------------------------------------
    if(Debug) print *,'PIO_initDecomp: CalcDisplace',ndisp,size(iodof),lenblocks
    call CalcDisplace(lenBLOCKS,ioDof,displace)
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif

    n_iotasks = Iosystem%num_iotasks
    length = size(iodof)
    !
    !   This facilitates the use of seperate read and write descripters. 
    !
    IOdesc%IOmap%start  = Iosystem%io_rank*length
    IOdesc%IOmap%length = length
    IOdesc%glen = glength

    if(Debug) print *,'IAM: ',Iosystem%io_rank,'initDecomp: UseRearranger: ',UseRearranger, glength
    if(UseRearranger) then 
       call rearrange_create(Iosystem,compDOF,ioDOF,ioDesc)
    endif
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif


    !---------------------------------------------
    !  The setup for the MPI-IO type information 
    !---------------------------------------------
    if(Iosystem%IOproc) then 
       !-----------------------------------------------
       ! setup the data structure for the io operation 
       !-----------------------------------------------
       IODesc%Write%n_elemTYPE = SIZE(displace)
       IODesc%Write%n_words    = IODesc%Write%n_elemTYPE*lenBLOCKS
       call GenIndexedBlock(lenBLOCKS,baseTYPE,IODesc%Write%elemTYPE,IODesc%Write%fileTYPE,displace)

       if(Debug) print *,'initDecomp: At the end of subroutine',IODesc%Write%n_elemTYPE,IODesc%Write%n_words
    endif
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
    call dupiodesc2(iodesc%write,iodesc%read)
    if(Debug) then
       print *, _FILE_,__LINE__,iodesc%read%filetype,iodesc%read%elemtype,&
            iodesc%read%n_elemtype,iodesc%read%n_words   
       print *, _FILE_,__LINE__,iodesc%write%filetype,iodesc%write%elemtype,&
            iodesc%write%n_elemtype,iodesc%write%n_words
    end if
    call dealloc_check(displace)

#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
#ifdef TIMING
    call t_stopf("pio_initDecomp")
#endif
  end subroutine initDecomp_1dof_nf




  subroutine initDecomp_1dof_nf_box(Iosystem,basepioTYPE,dims,compDOF, IOdesc, start, count, method)
    type (Iosystem_desc_t), intent(inout) :: Iosystem
    integer(i4), intent(in)           :: basepioTYPE
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: compDOF(:)   ! Global degrees of freedom for computational decomposition
    integer (kind=PIO_OFFSET), optional :: start(:), count(:)    
    integer (i4), optional, intent(in) :: method
    type (IO_desc_t), intent(out)     :: IOdesc

    integer(i4) :: length,n_iotasks
    integer(i4) :: ndims
    integer (i4)                       :: lenBLOCKS
    integer(i4)                       :: baseTYPE

    integer (i4), pointer :: displace(:)  ! The displacements for the MPI data structure (Read)

    integer(i4) :: prev
    integer(i4) :: gLength    ! global length in words
    integer(i4) :: ii,i,dis,ierr
    integer(i4),pointer, dimension(:) :: blocklen,disp
    logical(log_kind) ::  UseRearranger
    logical, parameter :: Check = .TRUE.
    integer(i4) :: nDisp
    integer(i4) :: iosize               ! rml
    integer(i4) :: mthd
#ifdef MEMCHK
    integer :: msize, rss, mshare, mtext, mstack
#endif

    integer ierror

    if(minval(dims)<=0) then
       print *,_FILE_,__LINE__,dims
       call piodie(_FILE_,__LINE__,'bad value in dims argument')
    end if

#ifdef TIMING
    call t_startf("pio_initDecomp_box")
#endif
    if (Iosystem%comp_rank == 0 .and. Debug) &
         print *,Iosystem%comp_rank,': Invoking initDecomp_1dof_nf_box'

    baseTYPE=pio_type_to_mpi_type(basepioTYPE)

    !-------------------------------------------
    ! For testing purposes Set the IOMap
    ! (DecompMap_t) to something basic for
    ! testing.
    !-------------------------------------------
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
    UseRearranger = Iosystem%UseRearranger
    !---------------------
    ! number of dimensions
    !---------------------
    ndims = SIZE(dims)
    !---------------------
    ! total global size
    !---------------------

    gLength=1
    do i=1,ndims
       gLength = gLength*dims(i)
    enddo


    ! remember count() is only defined on io procs
    call alloc_check(IODESC%start,ndims)
    call alloc_check(IODESC%count,ndims)
    iodesc%baseTYPE=baseTYPE

    iodesc%compsize=size(compdof)

    iodesc%start=0
    iodesc%count=0
    if (Iosystem%IOproc) then

       if(present(method)) then
          mthd=method
       else
          mthd=2
       end if
       if(present(start) .and. present(count)) then
          iodesc%start = start
          iodesc%count = count
       else if(present(start) .or. present(count)) then
          call piodie( _FILE_,__LINE__, &
               'Both optional parameters start and count must be provided')
       else
          call getiostartandcount(ndims, dims, Iosystem%num_iotasks, Iosystem%io_rank, IODESC%start, IODESC%count)

       end if
       iosize=1
       do i=1,ndims
          iosize=iosize*IODESC%count(i)
       end do
       call MPI_ALLREDUCE(iosize, iodesc%maxiobuflen, 1, mpi_integer, mpi_sum, iosystem%io_comm, ierr)
       call CheckMPIreturn('MPI_ALLREDUCE in initdecomp',ierr)


       lenBLOCKS=IODESC%count(1)


       if(lenblocks>0) then
          nDisp=iosize/lenBLOCKS
       else
          ndisp=0
       end if
       call alloc_check(displace,nDisp)

       !--------------------------------------------
       ! calculate MPI data structure displacements 
       !--------------------------------------------

       if(Debug) print *,'PIO_initDecomp: CalcDisplace', &
            ndisp,iosize,lenblocks, iodesc%start, iodesc%count
       call CalcDisplace_box(dims,IODESC%start,IODESC%count,ndims,displace)

       n_iotasks = Iosystem%num_iotasks
       length = iosize                      ! rml

       !
       !   This facilitates the use of seperate read and write descripters. 
       !

       IOdesc%IOmap%start  = Iosystem%io_rank*length
       IOdesc%IOmap%length = length
       IOdesc%glen = glength
    endif


#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
    if(Debug) print *,'IAM: ',Iosystem%io_rank, &
         'initDecomp: UseRearranger: ',UseRearranger, glength

    if(UseRearranger) then 
       call rearrange_create( Iosystem,compDOF,dims,ndims,ioDesc)
    endif

#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif



    !---------------------------------------------
    !  The setup for the MPI-IO type information 
    !---------------------------------------------
    if(Iosystem%IOproc) then 
       !-----------------------------------------------
       ! setup the data structure for the io operation 
       !-----------------------------------------------
       IODesc%Write%n_elemTYPE = SIZE(displace)
       IODesc%Write%n_words    = IODesc%Write%n_elemTYPE*lenBLOCKS
       call GenIndexedBlock(lenBLOCKS,baseTYPE,IODesc%Write%elemTYPE,IODesc%Write%fileTYPE,displace)

       if(Debug) print *,'initDecomp: At the end of subroutine',IODesc%Write%n_elemTYPE,IODesc%Write%n_words, &
            iodesc%write%elemtype

       IODesc%Read%n_elemTYPE = SIZE(displace)
       IODesc%Read%n_words    = IODesc%Read%n_elemTYPE*lenBLOCKS
       call GenIndexedBlock(lenBLOCKS,baseTYPE,IODesc%Read%elemTYPE,IODesc%Read%fileTYPE,displace)
    else
       iodesc%write%n_elemtype=0
       iodesc%write%n_words=0
       iodesc%write%elemType=-1
       iodesc%write%fileType=-1
    endif

#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif

    call dupiodesc2(iodesc%write,iodesc%read)

    if(Debug) then
       print *, _FILE_,__LINE__,iodesc%read%filetype,iodesc%read%elemtype,&
            iodesc%read%n_elemtype,iodesc%read%n_words   
       print *, _FILE_,__LINE__,iodesc%write%filetype,iodesc%write%elemtype,&
            iodesc%write%n_elemtype,iodesc%write%n_words
    end if

    if (Iosystem%IOproc) then
       call dealloc_check(displace)
    endif

#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
#ifdef TIMING
    call t_stopf("pio_initDecomp_box")
#endif

  end subroutine initDecomp_1dof_nf_box



  subroutine getiostartandcount(ndims, gdims, num_io_procs, iorank, start, count)

    implicit none

    real, parameter :: eps = 1.0e-4
    integer, intent(in) :: ndims, gdims(ndims), num_io_procs

    integer(kind=PIO_Offset), intent(out) :: start(ndims), count(ndims)   ! Start and count arrays 


    integer :: n,m
    integer :: use_io_procs,iorank, sdims, cnt
    logical :: done
    integer :: xiam, xpes, ps, pe, ds, de
    integer :: size,tsize
    integer, parameter :: minblocksize=16        ! minimum block size on a task

    !    integer(kind=PIO_Offset), allocatable :: start(:,:)
    !    integer(kind=PIO_Offset), allocatable :: count(:,:)


    tsize=1
    do n=1,ndims
       tsize=tsize*gdims(n)
    end do

    use_io_procs=num_io_procs
    do while(tsize/minblocksize < use_io_procs .and. use_io_procs>1)
       use_io_procs=use_io_procs-1
    end do

    start = 1
    count = 0


    if(iorank>=use_io_procs) return 

    !-----------------

    cnt = 1
    sdims = ndims
    do while (cnt < use_io_procs .and. sdims > 0)
       cnt = cnt * gdims(sdims)
       sdims = sdims - 1
    enddo
    if (sdims < 0) then
       write(6,*) 'error in sdims ', sdims
       stop
    endif

    do m = 1,sdims
       start(m) = 1
       count(m) = gdims(m)
    enddo

    xpes = use_io_procs
    xiam = iorank   ! goes from 0 to xpes-1
    do m = ndims, sdims+1, -1

       ds = int(float(gdims(m)*(xiam  ))/float(xpes) + 1)
       de = int(float(gdims(m)*(xiam+1))/float(xpes))

       start(m) = ds
       count(m) = max(de-ds+1,1)

       !     write(6,*) 'tcx1 ',n,m,start(n,m),count(n,m),iorank

       ps = int(float(xpes*(start(m)-1))/float(gdims(m)) - eps + 1)
       pe = int(float(xpes*(start(m)  ))/float(gdims(m)) - eps)

       xpes = pe - ps + 1
       xiam = xiam - ps

       !     write(6,*) 'tcx2 ',n,m,ps,pe,xpes,xiam

    enddo



  end subroutine getiostartandcount






  !************************************
  ! dupiodesc2
  !

  subroutine dupiodesc2(src, dest)
    type(io_desc2_t), intent(in) :: src
    type(io_desc2_t), intent(out) :: dest

    dest%fileType = src%fileType
    dest%elemType = src%elemType
    dest%n_elemType = src%n_elemType
    dest%n_words = src%n_words
  end subroutine dupiodesc2



  !************************************
  ! GenIndexedBlock
  !
  ! given input lenBLOCKS, baseTYPE, and displacement
  ! create two MPI types: 
  !   elemTYPE - a single block of baseTYPE repeated lenBLOCKS times
  !   fileTYPE - elemTYPE repeated at each entry in displacement()
  !              (i.e. size(displacement) entries)
  !


  subroutine GenIndexedBlock(lenBLOCKS,baseTYPE,elemTYPE,fileTYPE,displace)
#ifdef NO_MPI2
    use pio_support, only : MPI_TYPE_CREATE_INDEXED_BLOCK
#endif
    integer(i4), intent(in) :: lenBLOCKS     ! length of blocks
    integer(i4), intent(in) :: baseTYPE      ! base MPI type 
    integer(i4), intent(inout) :: elemTYPE   ! elementary MPI type
    integer(i4), intent(inout) :: fileTYPE   ! file MPII type 
    integer(i4), intent(in) :: displace(:)   ! MPI displacement in the array

    integer(i4) :: numBLOCKS,i,ierr,prev
    logical, parameter :: Check = .TRUE.

    integer:: nints, nadds, ndtypes, comb, lbasetype

    if (lenBLOCKS < 1) then
       return
       call piodie( _FILE_,__LINE__, &
            'GenIndexedBlock: lenBLOCKS=',lenBLOCKS,' (expected >=1)')
    endif

    numBLOCKS = SIZE(displace)

    !tcx - allow empty iodofs
    if (numBLOCKS > 0) then

       if(Debug) print *,'GenIndexedBlock: Setting up MPI-IO data structure stuff, numBLOCKS, lenBLOCKS',&
            numBLOCKS,lenBLOCKS
       prev = displace(1)
       do i=2,numBLOCKS
          if(prev > displace(i)) then
             print *,'GenIndexedBlock: error detected: non-monotonic increasing displace detected!'
          endif
          prev = displace(i)
       enddo

    endif
    select case(basetype)
    case (PIO_DOUBLE)
       lbasetype=MPI_REAL8
    case (PIO_REAL  )
       lbasetype=MPI_REAL4
    case (PIO_INT)
       lbasetype=MPI_INTEGER
    case (PIO_CHAR)
       lbasetype=MPI_CHARACTER
    case default
       lbasetype=basetype
    end select


#ifdef _MPISERIAL
    ! When compiling w/mpiserial for snetcdf output, these fields are not used
    elemTYPE=0
    fileTYPE=0
    ! _MPISERIAL
#else

    call MPI_Type_contiguous(lenBLOCKS,lbaseTYPE,elemTYPE,ierr)
    if(Check) call CheckMPIreturn('GenIndexedBlock: after call to Type_contiguous: ',ierr)
    call MPI_Type_commit(elemTYPE,ierr)
    if(Check) call CheckMPIreturn('GenIndexedBlock: after call to Type_commit: ',ierr)
    call MPI_Type_create_indexed_block(numBLOCKS,1,displace,elemTYPE,fileTYPE,ierr)
    if(Check) call CheckMPIreturn('GenIndexedBlock: after call to Type_create_indexed_block: ',ierr)
    call MPI_Type_commit(fileTYPE,ierr)
    if(Check) call CheckMPIreturn('GenIndexedBlock: after call to Type_commit: ',ierr)
    call mpi_type_get_envelope(elemtype, nints, nadds, ndtypes, comb, ierr)

    ! _MPISERIAL
#endif

  end subroutine GenIndexedBlock



  !************************************
  ! init
  !

  subroutine init(comp_rank, comp_comm, num_iotasks, num_aggregator, stride,  rearr, IOsystem,base)

    integer(i4), intent(in) :: comp_rank
    integer(i4), intent(in) :: comp_comm
    integer(i4), intent(in) :: num_iotasks 
    integer(i4), intent(in) :: num_aggregator
    integer(i4), intent(in) :: stride
    !    integer(i4), intent(in) :: ioType
    integer(i4), intent(in) :: rearr
    type (IOsystem_desc_t), intent(out)  :: IOsystem  ! IO descriptor to initalize
    integer(i4), intent(in),optional :: base

    integer(i4) :: n_iotasks
    integer(i4) :: length
    integer(i4) :: ngseg,io_rank,i,lbase, io_comm,ierr 
    integer(i4) :: lstride
    integer(i4), pointer :: ioranks(:)

    character(len=5) :: cb_nodes
    logical(log_kind), parameter :: Check = .TRUE.
    logical(log_kind) :: usehints

    integer(i4) :: MPI_GROUP_WORLD, MPI_GROUP_IO

#ifdef TIMING
    call t_startf("pio_init")
#endif
    IOsystem%error_handling = PIO_INTERNAL_ERROR
    IOsystem%comp_comm = comp_comm
    IOsystem%comp_rank = comp_rank
    IOsystem%rearr = rearr

    call MPI_COMM_size(comp_comm,IOsystem%num_tasks,ierr)
    if(Check) call CheckMPIreturn('init: after call to COMM_size: ',ierr)
    ! ---------------------------------------
    ! need some more error checking code for 
    ! setting of number of IO nodes
    ! ---------------------------------------

    n_iotasks=num_iotasks

    if (n_iotasks>IOsystem%num_tasks) then
       n_iotasks=IOsystem%num_tasks
       if (IOsystem%comp_rank==0) then
          print *,'***WARNING, reducing io tasks to ',n_iotasks, &
               ' because there are not enough processors'
       endif
    endif

    lbase = 0
    ! unless you are using all procs, shift off the masterproc
    if(n_iotasks<IOsystem%num_tasks) then
       lbase=1
    end if
    if (present(base)) then
       if(base>0 .and. base<IOsystem%num_tasks) lbase = base
    endif
    IOsystem%IOMaster=lbase

    IOsystem%num_iotasks = n_iotasks
    if(Debug) print *,'init: IOsystem%num_tasks,n_iotasks,num_aggregator: ',IOsystem%num_tasks,n_iotasks,num_aggregator

    ! --------------------------
    ! select which nodes are IO
    ! nodes and set IOProc
    ! --------------------------
    lstride = stride

    if (lbase+(n_iotasks-1)*lstride >= IOsystem%num_tasks) then
       print *,_FILE_,__LINE__,lbase,n_iotasks,lstride,IOsystem%num_tasks
       call piodie(_FILE_,__LINE__,'Not enough procs for the stride')
    endif

    call alloc_check(ioranks,n_iotasks,'init:n_ioranks')
    IOsystem%IOProc = .FALSE.
    do i=1,n_iotasks
       ioranks(i)=(lbase + (i-1)*lstride)

       if (ioranks(i)>=IOsystem%num_tasks) then
          call piodie( _FILE_,__LINE__, &
               'tried to assign io processor beyond max rank ',&
               ioranks(i), &
               ' num_tasks=',IOsystem%num_tasks )
       endif

       if(comp_rank == ioranks(i))  IOsystem%IOProc = .TRUE.
    enddo

    ! rml stride looks like it was never set before
    IOsystem%io_stride=lstride       ! rml fix

    IOsystem%io_rank=0

    !-----------------------
    ! Setup IO_comm and IO_rank
    !-----------------------
    IOsystem%io_rank=-1
    call MPI_COMM_GROUP(comp_comm,MPI_GROUP_WORLD,ierr)
    if(Check) call CheckMPIreturn('init: after call to COMM_GROUP: ',ierr)

    call MPI_GROUP_INCL(MPI_GROUP_WORLD,n_iotasks,ioranks,MPI_GROUP_IO,ierr)
    if(Check) call CheckMPIreturn('init: after call to GROUP_RANGE_INCL: ',ierr)

    call dealloc_check(ioranks)


    call MPI_COMM_CREATE(comp_comm,MPI_GROUP_IO,IOsystem%io_comm,ierr)
    if(Check) call CheckMPIreturn('init: after call to COMM_CREATE: ',ierr)

    if(IOsystem%IOproc) call MPI_COMM_rank(IOsystem%io_comm,IOsystem%io_rank,ierr)
    if(Check) call CheckMPIreturn('init: after call to COMM_rank: ',ierr)


    if(Debug) print *,'init: IAM: ',comp_rank,'IO processor: ',IOsystem%IOproc, 'IO rank ',IOsystem%io_rank



    if(Debug) print *,'init: IAM: ',comp_rank,' before allocate(status): n_iotasks: ',n_iotasks

    usehints = .FALSE.

    if (rearr == PIO_rearr_none) then
       IOsystem%UseRearranger= .FALSE.
    else
       IOsystem%UseRearranger= .TRUE.
    endif


#if defined(USEMPIIO) || defined(_PNETCDF)

    call MPI_info_create(IOsystem%info,ierr)
    ! Turn on MPI-IO aggregation 
    if(num_aggregator .gt. 0) then 
       call setNUMAgg(IOsystem,num_aggregator)  ! let mpi-io do aggregation
       usehints = .TRUE.
    endif

#ifdef PIO_GPFS_HINTS
    call MPI_info_set(IOsystem%info,"IBM_largeblock_io","true",ierr)
    usehints = .TRUE.
#endif
#ifdef PIO_LUSTRE_HINTS
    call MPI_info_set(IOsystem%info,'romio_ds_read','disable',ierr)
    call MPI_info_set(IOsystem%info,'romio_ds_write','disable',ierr)
    usehints = .TRUE.
#endif

    if(.not. usehints) IOsystem%info = MPI_INFO_NULL

#endif

    if(Debug) print *,'IAM: ',IOsystem%io_rank,__LINE__,'init: UseRearranger: ',IOsystem%UseRearranger

    !---------------------------------
    ! initialize the rearranger system 
    !---------------------------------

    if (IOsystem%UseRearranger) then
       call rearrange_init(IOsystem)
    endif
#ifdef TIMING
    call t_stopf("pio_init")
#endif

  end subroutine init



  !************************************
  ! setNUMAgg
  !

  subroutine setNUMAgg(iosystem,numAgg)
    type (iosystem_desc_t), intent(inout) :: iosystem
    integer(i4), intent(in)  :: numAgg

    character(len=5) :: cb_nodes
    integer(i4) :: ierr
    logical(log_kind), parameter ::  Check = .true.
#if defined(USEMPIIO) || defined(_PNETCDF)
    write(cb_nodes,('(i5)')) numAgg
    call MPI_info_create(iosystem%info,ierr)    
    if(Check) call CheckMPIreturn('setnumagg: after call to MPI_info_create: ',ierr)

#ifdef _BGL
    call MPI_info_set(iosystem%info,'bgl_cb_nodes',TRIM(ADJUSTL(cb_nodes)),ierr)
#else
    call MPI_info_set(iosystem%info,'cb_nodes',TRIM(ADJUSTL(cb_nodes)),ierr)
#endif
#endif
  end subroutine setNUMAgg





  !=============================================
  !  dupIOdesc:
  !
  !   Duplicate the IO descriptor
  !
  !=============================================


  ! rml: possible problem here wrt dubbing the box rearranger
  ! data, as well as maybe the mct rearranger???

  subroutine dupIOdesc(src,dest)

    integer :: n
    type (IO_desc_t), intent(in) :: src
    type (IO_desc_t), intent(inout) :: dest


    dest%glen        =  src%glen
    if(associated(src%start)) then
       n = SIZE(src%start)
       allocate(dest%start(n))
       dest%start(:)       =  src%start(:)
    endif

    if(associated(src%count)) then
       n = SIZE(src%count)
       allocate(dest%count(n))
       dest%count(:)       =  src%count(:)
    endif

    !DBG    print *,'before dupiodesc2'
    call dupiodesc2(src%Read, dest%Read)
    call dupiodesc2(src%Write, dest%Write)
    !DBG    print *,'after dupiodesc2'
#ifdef _USEMCT
    !DBG    print *, _FILE_,__LINE__, 'This function is broken'
    dest%lsize_comp = src%lsize_comp
    dest%lsize_io = src%lsize_io

    if(associated(src%compDOF_index)) then
       n = SIZE(src%compDOF_index)
       allocate(dest%compDOF_index(n))
       dest%compDOF_index(:) = src%compDOF_index(:)
    endif

#endif
    dest%baseTYPE = src%baseTYPE

    if(associated(src%dest_ioproc)) then 
       n = SIZE(src%dest_ioproc)
       allocate(dest%dest_ioproc(n))
       dest%dest_ioproc(:) = src%dest_ioproc(:)
    endif

    if(associated(src%dest_ioindex)) then 
       n = SIZE(src%dest_ioindex)
       allocate(dest%dest_ioindex(n))
       dest%dest_ioindex(:) = src%dest_ioindex(:)
    endif

    if(associated(src%rfrom)) then 
       n = SIZE(src%rfrom)
       allocate(dest%rfrom(n))
       dest%rfrom(:) = src%rfrom(:)
    endif

    if(associated(src%rtype)) then 
       n = SIZE(src%rtype)
       allocate(dest%rtype(n))
       dest%rtype(:) = src%rtype(:)
    endif

    if(associated(src%scount)) then 
       n = SIZE(src%scount)
       allocate(dest%scount(n))
       dest%scount(:) = src%scount(:)
    endif

    if(associated(src%stype)) then 
       n = SIZE(src%stype)
       allocate(dest%stype(n))
       dest%stype(:) = src%stype(:)
    endif

    call copy_decompMap(src%IOmap,dest%IOmap)
    call copy_decompMap(src%COMPmap,dest%COMPmap)

    dest%compsize = src%compsize


  end subroutine dupIOdesc



  !=============================================
  !  copy_decompMAP:
  !
  !   Copy decompMAP_t data structures
  !
  !=============================================

  subroutine copy_decompMap(src,dest)

    type (DecompMap_t), intent(in) :: src
    type (DecompMap_t), intent(inout) :: dest


    dest%start    = src%start
    dest%length   = src%length

  end subroutine copy_decompMap



  !=============================================
  !  setIOType:
  !   
  !   Set the desired type of IO to perform 
  !
  !=============================================

  subroutine setIOType(File,iotype,rearr)

    type (file_desc_t), intent(inout) :: File
    integer(i4), intent(in) :: iotype 
    integer(i4), intent(in) :: rearr

    file%iotype = iotype
    file%Iosystem%rearr = rearr

  end subroutine setIOType



  !*****************************************
  ! numToRead_io
  !
  !*****************************************

  integer function numToRead(ioDesc) result(num)

    type (IO_desc_t) :: ioDesc

    num = ioDesc%Read%n_words

  end function numToRead



  !*****************************************
  ! numToWrite
  !
  !*****************************************

  integer function numToWrite(ioDesc) result(num)

    type (IO_desc_t) :: ioDesc

    num = ioDesc%Write%n_words

  end function numToWrite




  !=============================
  ! CreateFile:
  !
  !  Open a disk file
  ! ============================

  integer function CreateFile(IOSystem, File,iotype, fname, amode_in) result(ierr)
    type (iosystem_desc_t), intent(inout), target :: IOsystem
    type (File_desc_t), intent(out) :: File
    integer, intent(in) :: iotype
    character(len=*), intent(in)  :: fname
    integer, optional, intent(in) :: amode_in
    ! ===================
    !  Local variables
    ! ===================
    integer    :: amode
    logical, parameter :: Check = .TRUE.
    character(len=9) :: rd_buffer

#ifdef TIMING
    call t_startf("pio_createfile")
#endif

    if(Debug) print *,'CreateFile: {comp,io}_rank:',IOsystem%comp_rank,IOsystem%io_rank,'IO proc: ',IOsystem%IOproc
    ierr=PIO_noerr

    if(present(amode_in)) then
       amode = amode_in
    else	
       amode = 0
    end if

    File%iosystem => iosystem

    !--------------------------------
    ! Set some IOtype specific stuff
    !--------------------------------

    if(iosystem%num_iotasks.eq.1.and.iotype.eq.iotype_pnetcdf) then	
#if defined(_NETCDF)
       File%iotype=iotype_netcdf
#else
       File%iotype = iotype 
#endif
    else
       File%iotype = iotype 
    end if
#if defined(USEMPIIO) || defined(_PNETCDF)
    if ( (File%iotype==iotype_pbinary .or. File%iotype==iotype_direct_pbinary) &
         .and. (.not. IOsystem%UseRearranger) ) then
       write(rd_buffer,('(i9)')) 16*1024*1024
       call MPI_info_set(IOsystem%info,'cb_buffer_size',TRIM(ADJUSTL(rd_buffer)),ierr)
    endif
#endif


    select case(iotype)
    case(iotype_pbinary, iotype_direct_pbinary)
       if(present(amode_in)) then
          print *, 'Warning, the mode argument is currently ignored for binary file operations'
       end if
       ierr = create_mpiio(File,fname)
    case( iotype_pnetcdf, iotype_netcdf)
       ierr = create_nf(File,fname, amode)	
       if(debug .and. iosystem%io_rank==0)print *,_FILE_,__LINE__,' open: ', fname, file%fh
    case(iotype_binary)
       print *,'CreateFile: IO type not supported'
    end select


    if(Debug) print *, _FILE_,__LINE__,'CreateFile complete', amode
#ifdef TIMING
    call t_stopf("pio_createfile")
#endif
  end function CreateFile



  !=============================
  ! OpenFile:
  !
  !  Open a disk file
  ! ============================

  integer function OpenFile(IOsystem, File, iotype, fname,mode) result(ierr)
    type (IOsystem_desc_t), intent(inout), target :: IOsystem
    type (File_desc_t), intent(out) :: File
    integer, intent(in) :: iotype
    character(len=*), intent(in)  :: fname
    integer, optional, intent(in) :: mode
    ! ===================
    !  Local variables
    ! ================
    integer    :: amode
    logical, parameter :: Check = .TRUE.
    character(len=9) :: rd_buffer

#ifdef TIMING
    call t_startf("pio_openfile")
#endif

    if(Debug) print *,'OpenFile: {comp,io}_rank:',Iosystem%comp_rank,Iosystem%io_rank,'IO proc: ',Iosystem%IOproc

    ierr=PIO_noerr

    File%iosystem => iosystem

    !--------------------------------
    ! Set some IOtype specific stuff
    !--------------------------------

    if(iosystem%num_iotasks.eq.1.and.iotype.eq.iotype_pnetcdf) then	
#if defined(_NETCDF)
       File%iotype=iotype_netcdf
#else
       File%iotype = iotype 
#endif       
    else
       File%iotype = iotype 
    end if
#if defined(USEMPIIO) || defined(_PNETCDF)
    if ( (File%iotype==iotype_pbinary .or. File%iotype==iotype_direct_pbinary) &
         .and. (.not. IOsystem%UseRearranger) ) then
       write(rd_buffer,('(i9)')) 16*1024*1024
       call MPI_info_set(IOsystem%info,'cb_buffer_size',TRIM(ADJUSTL(rd_buffer)),ierr)
    endif
#endif
    select case(iotype)
    case(iotype_pbinary, iotype_direct_pbinary)
       if(present(mode)) then
          print *, 'Warning, the mode argument is currently ignored for binary file operations'
       end if
       ierr = open_mpiio(File,fname)
    case( iotype_pnetcdf, iotype_netcdf)
       if(present(mode)) then
          ierr = open_nf(File,fname,mode)
       else
          ierr = open_nf(File,fname)
       end if
       if(debug .and. iosystem%io_rank==0)print *,_FILE_,__LINE__,' open: ', fname, file%fh
    case(iotype_binary)   ! appears to be a no-op

    end select

#ifdef TIMING
    call t_stopf("pio_openfile")
#endif

  end function OpenFile



  !=============================
  ! Sync a netcdf file
  ! ============================
  subroutine SyncFile(File)
    implicit none
    type (File_desc_t) :: File
    integer :: ierr

    select case(File%iotype)
    case( iotype_pnetcdf, iotype_netcdf)
       ierr = sync_nf(File)
    case(iotype_pbinary, iotype_direct_pbinary)
    case(iotype_binary) 
    end select
  end subroutine SyncFile



  !=============================
  ! Free a decomp
  !
  !  free all allocated mct storage
  ! ============================
  subroutine freeDecomp_ios(IOS,ioDesc)
    implicit none
    type (IOsystem_desc_t) :: IOS
    type (IO_desc_t) :: ioDesc

    call rearrange_free(IOS,ioDesc)
    if(associated(IODESC%start)) then
       call dealloc_check(IODESC%start,'iodesc%start')
       nullify(IODESC%start)
    end if

    if(associated(IODESC%count)) then
       call dealloc_check(IODESC%count,'iodesc%count')    
       nullify(IODESC%count)
    end if
  end subroutine freeDecomp_ios

  !=============================
  ! Free a decomp
  !
  !  free all allocated mct storage
  ! ============================
  subroutine freeDecomp_file(File,ioDesc)
    implicit none
    type (File_desc_t) :: File
    type (IO_desc_t) :: ioDesc

    call freedecomp_ios(File%iosystem, iodesc)

  end subroutine freeDecomp_file



  !=============================
  ! CloseFile:
  !
  !  Close a disk file
  ! ============================

  subroutine CloseFile(File)

    type (File_desc_t),intent(inout)   :: File

    integer :: ierr
    integer :: iotype 
    logical, parameter :: Check = .TRUE.

#ifdef TIMING
    call t_startf("pio_closefile")
#endif
    if(debug) print *,_FILE_,__LINE__,'close: ',File%fh
    iotype = File%iotype 
    select case(iotype)
    case(iotype_pbinary, iotype_direct_pbinary)
       ierr = close_mpiio(File)
    case( iotype_pnetcdf, iotype_netcdf)
       ierr = close_nf(File)
    case(iotype_binary)
       print *,'CloseFile: IO type not supported'
    end select
#ifdef TIMING
    call t_stopf("pio_closefile")
#endif


  end subroutine CloseFile


  !******************************
  ! read_ascii
  !

  subroutine read_ascii(rank,IOBUF,size)

    integer, intent(in) :: rank
    real (r8), dimension(:) :: IOBUF
    integer, intent(in) :: size

    character(len=80) filename
    integer lun
    integer ios
    integer i

    lun=10+rank
    write(filename,"('fort.',I2)" ) lun
    write(6,*) 'Filename is:', filename

    open(lun,FILE=filename,STATUS='OLD',IOSTAT=ios)
    if (ios /= 0) then
       write(6,*) rank,': Could not open ascii file: ',filename
    endif

    do i=1,size
       read(UNIT=lun,FMT=*,IOSTAT=ios) IOBUF(i)
       if (ios /= 0) then
          write (6,*) rank,': Error reading item ',i,' of ',size
          call abort
       endif

    end do



100 continue

    close(lun)

  end subroutine read_ascii



  !*****************************
  ! CalcDisplace
  !

  subroutine CalcDisplace(bsize,dof,displace)

    integer(i4), intent(in) :: bsize    ! length of contigious blocks of numbers
    integer(i4), intent(in) :: dof(:)   ! degree of freedom on which to setup the displacement array
    integer(i4), intent(inout) :: displace(:)  ! array of MPI displacments

    integer :: numBLocks,lenBlocks,i,ii,dis

    numBLOCKS = SIZE(displace)
    lenBLOCKS = bsize
    do i=1,numBLOCKS
       ii = (i-1)*lenBLOCKS+1
       dis = DOF(ii)-1
       dis = dis/lenBLOCKS   
       displace(i) = dis
    enddo
    do i=1,numBLOCKS-1	
       if(displace(i+1) .lt. displace(i)) then
          print *,'CalcDisplace: Error with displacement arrays',i,displace(i:i+1),numblocks,size(dof),dof(numblocks)
          call piodie( _FILE_,__LINE__)
       endif
    enddo

  end subroutine CalcDisplace


  !*********************************
  !
  ! new version
  !

  subroutine CalcDisplace_box(gsize,start,count,ndim,displace)

    integer(i4),intent(in) :: gsize(:)   ! global size of output domain
    integer(kind=PIO_Offset),intent(in) :: start(:), count(:)
    integer(i4), intent(in) :: ndim
    integer(i4),intent(inout) :: displace(:)  ! MPI displacments

    !!

    integer ndisp
    integer(i4) :: gstride(ndim)
    integer i,j
    integer iosize
    integer(i4) :: myloc(ndim)
    integer(i4) :: ub(ndim)
    integer idim
    logical done
    integer gindex

    gstride(1)=gsize(1)
    do i=2,ndim
       gstride(i)=gsize(i)*gstride(i-1)
    end do

    iosize=1
    do i=2,ndim
       iosize=iosize*count(i)
    end do

    ndisp=size(displace)

    if (iosize<1 .or. ndisp<1) return


    if (ndisp/=iosize) then
       call piodie(_FILE_,__LINE__,'ndisp=',ndisp,' /= iosize=',iosize)
    endif

    do i=1,ndim
       ub(i)=start(i)+count(i)-1
    end do

    ! skip x dimension (start of each block)
    ! generate displacement for every 1,y,z
    !  i.e. loop over y,z,...
    !       compute corresponding global index
    !       divide by lenBLOCKS

    displace(1)=1
    myloc=start

    do i=1,iosize
       ! go from myloc() to 1-based global index
       gindex=myloc(1)
       do j=2,ndim
          gindex=gindex+(myloc(j)-1)*gstride(j-1)
       end do

       ! rml
       ! following original but is that right???
       ! seems like the 'if' is erroneous

       gindex=gindex-1

       gindex=gindex/count(1)    ! gindex/lenblock

       displace(i)=gindex

       ! increment myloc to next position


       idim=2                    ! dimension to increment
       done=.false.

       if (i<iosize) then
          do while (.not. done)
             if (myloc(idim)<ub(idim)) then
                myloc(idim)=myloc(idim)+1
                done=.true.
             else
                myloc(idim)=start(idim)
                idim=idim+1
                if (idim>ndim) call piodie(_FILE_,__LINE__,'dim overflow')
             endif
          end do
       endif

    end do

    do i=2,ndim
       if (myloc(i) /= ub(i)) then
          print *,'myloc=',myloc
          print *,'ub=',ub
          call piodie( _FILE_,__LINE__,'myloc/=ub')
       endif
    end do


    ! check for strictly increasing

    do i=1,nDisp-1	
       if(displace(i) .gt. displace(i+1)) then
          call piodie(_FILE_,__LINE__,'displace is not increasing')
       endif
    enddo

  end subroutine CalcDisplace_box



end module piolib_mod

  !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
