#ifdef HAVE_CONFIG_H
#include "config.h"
#endif



module restart_io_mod 
   !------------------
   use kinds , only : int_kind, real_kind
   !------------------
   use dimensions_mod, only : nelem, ne, np, nlev, nelemd
   !------------------
   use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
#ifdef _MPI
   use parallel_mod, only : iam, mpiinteger_t, mpireal_t, mpi_status_size, &
      mpi_success, mpi_max, syncmp, haltmp, &
      MPIinteger_t,parallel_t, &
      MPI_BYTE, abortmp
   use parallel_mod, only : abortmp
#ifdef _PRESTART
   use parallel_mod, only :  mpi_offset_kind, mpi_address_kind, mpi_mode_wronly,&
      mpi_mode_create, mpi_info_null, mpi_mode_rdonly, mpi_order_fortran
#endif
#else
   use parallel_mod, only : iam, mpiinteger_t, mpireal_t,parallel_t,haltmp
#endif
   !------------------
   use time_mod, only : timelevel_t, nendstep, nmax
   !------------------
   use element_mod, only: elem_state_t, element_t
   !------------------
   use control_mod, only : restartdir, restartfile, columnpackage
   !------------------
#ifdef _PRIM
   use physics_types_mod, only : physics_state_t, pelem
#endif
   !------------------
   use schedule_mod, only : Schedule
   !------------------
   implicit none

! these should be accessed through parallel_mod, above:
! intel fortran will complain if they appear in both places
!#ifdef _MPI
!#include <mpif.h> ! _EXTERNAL
!#endif


   private 

   integer,public               :: nwordsRestartBuffer_t
   logical,public               :: collective_io_write = .false.
   logical,public               :: collective_io_read  = .false.

   !===============================================
   ! Define the RestartHeader Structure 
   !===============================================
! 
!  The sizeof this header structure is hardcoded in createrestartheader below.  If you change it 
!  here you must change it there as well.
!

   type, public :: RestartHeader_t
      integer(kind=int_kind)    :: size          !   1
      integer(kind=int_kind)    :: version       !   2
      integer(kind=int_kind)    :: ne            !   3
      integer(kind=int_kind)    :: np            !   4
      integer(kind=int_kind)    :: nlev          !   5
      integer(kind=int_kind)    :: number        !   6
      integer(kind=int_kind)    :: ElemRecLength !   7
      type (TimeLevel_t)        :: TimeLevel     !  12 
   end type

   integer, parameter :: RESTART_HDR_CNT = 12 


   ! MT 2010:  note: when using COLLECTIVE IO, the code we have now, 
   ! the offset between elements in 
   ! RestartBuffer must match the offset used in the restart file.
   ! This means for non-Emanual physics restarts, puffer cannot be 
   ! included in the RestartBuffer_t struct below.  
   !
   ! I tried to include the pelem() data in the restart file for all
   ! cases, including non-Emanual physics like Held-Suarez run, 
   ! but the code was crashing when creating the first MPI type, maybe because
   ! pelem() has not been initialized. 
   !
   ! best solution would be to modify the code to allow a larger offset
   ! between elements in RestartBuffer than in the file
   ! for now, disable restart when running Emanual physics.  
   ! 
   ! Emanual restart can be re-enabled by turning of COLLECTIVE IO 
   ! and adding the puffer struct back in below (change PRIMXXX back to PRIM)
   !
   type, public :: RestartBuffer_t
      type (elem_state_t)        :: buffer
#ifdef _PRIMXXX
      type (physics_state_t)     :: puffer
#endif
   end type


   type, public :: StateDesc_t
      integer(kind=int_kind)         :: numComponents
      integer(kind=int_kind)         :: ptr
      integer(kind=MPI_OFFSET_KIND),pointer :: disp(:)
      integer(kind=int_kind),pointer :: blklen(:)
      integer(kind=int_kind),pointer :: kind(:)
      integer(kind=int_kind),pointer :: type(:)
      integer(kind=int_kind)         :: nwords
      integer(kind=int_kind)         :: padding
   end type 

   type, public :: File_elem_t
      integer(kind=int_kind)    :: fh
      integer(kind=int_kind)    :: ElemState   ! MPI datastructure single element state
      integer(kind=int_kind)    :: State       ! MPI datastructure all elemnts on a process
      integer(kind=int_kind)    :: GlobalState ! MPI datastructure global view of State
      character(len=80)         :: fname       ! Name of File
      type (parallel_t)         :: par
   end type 
   type, public :: File_face_t
     integer(kind=int_kind)     :: fh
     integer(kind=int_kind)     :: face      ! MPI datastructure for a cube face
     integer(kind=int_kind)     :: buffer    ! MPI datastructure for a collection of faces
     integer(kind=int_kind)     :: GlobalFace ! MPI datastructure gloal view of face
     type (parallel_t)          :: par
   end type 

   ! =========================================
   !  Some variables used by all MPI routines 
   ! =========================================
   integer                         :: errorcode,errorlen,ierr
   character(len=80)               :: errorstring
   
   ! ====================================================
   !  Routines for Restart files 
   ! ====================================================
   public :: AddStateField
   public :: CreateStateDescriptor,PrintStateDescriptor
   public :: ConstructElementFile
#if defined(_MPI) && defined(_PRESTART)
   public :: PrintTypeInfo
#endif

   type (RestartBuffer_t),allocatable,target,public   :: RestartBuffer(:)


   integer,parameter              :: RestartVersion=4
   ! =========================================
   !  Some variables used by all MPI routines
   ! =========================================

#if defined(_MPI) && defined(_PRESTART)
   integer(kind=MPI_OFFSET_KIND)   :: offset,nbytes
   integer :: STATUS(MPI_STATUS_SIZE)
#else
   integer(kind=int_kind)          :: offset,nbytes
#endif

   type (File_elem_t),public          :: RestFile

   ! ====================================================
   !  Routines for Restart files
   ! ====================================================
   public :: ReadRestart,WriteRestart
contains 
!=============================
! WriteState
!
!  Writes a field out to disk
! ============================
   subroutine WriteState(File,variable,recl)
     implicit none
     type (File_elem_t),intent(in)                      :: File
     type (RestartBuffer_t),target,intent(in)      :: variable(:)
     integer, intent(in) :: recl

     !==================
     ! Local variables
     !==================
     integer           :: ie, ig
     integer           :: amode, info
#if defined(_MPI) && defined(_PRESTART)

     integer (kind=int_kind)	:: type_ablock, type_fview
     integer                    :: array_disp(nelemd)
     integer                    :: array_blen(nelemd)
     integer                    :: isiz
     integer*4			:: acount
     integer(kind=MPI_OFFSET_KIND) :: disp, iext,lb
#endif
     call t_startf('WriteState')

#if defined(_MPI) && defined(_PRESTART)
if (COLLECTIVE_IO_WRITE) then
     call mpi_type_contiguous( recl, MPI_BYTE, type_ablock, ierr )
     call mpi_type_commit( type_ablock, ierr )


     do ie=1,nelemd
#ifdef _PREDICT
       ig = Schedule(iam)%Local2Global(ie)
#else
       ig = Schedule(1)%Local2Global(ie)
#endif
       array_disp(ie) = ig-1
       array_blen(ie) = 1

     enddo

     disp = array_disp(1)
     do ie=1,nelemd
       array_disp(ie) = array_disp(ie) - disp
     enddo

     call mpi_type_create_indexed_block( nelemd, 1, array_disp, type_ablock, type_fview, ierr )
     call mpi_type_commit( type_fview, ierr )
     call mpi_type_get_extent( type_fview, iext, lb, ierr )
     call mpi_type_size( type_fview, isiz, ierr )
#if 0
     print *, 'debug: array_disp(1)=', array_disp(1), ', recl=', recl, ', iext=' &
          , iext, ', isiz=', isiz
#endif
endif


     amode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
     info  = MPI_INFO_NULL
     call MPI_file_open(File%par%comm,File%fname,amode,info,File%fh,ierr)

     if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        call abortmp(errorstring)
     endif

if ( COLLECTIVE_IO_WRITE ) then
     disp = disp * recl
     call MPI_File_set_view( File%fh, disp, type_ablock, type_fview, "native", MPI_INFO_NULL, ierr )
     call mpi_file_write_all( File%fh, variable, nelemd, type_ablock, status, ierr )
     if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        call abortmp(errorstring)
     endif

else

     do ie=1,nelemd
#ifdef _PREDICT
       ig = Schedule(iam)%Local2Global(ie)
#else
       ig = Schedule(1)%Local2Global(ie)
#endif
       ! offset = (ig-1)*recl   this will overflow integer*4 at high-res
       offset = ig-1
       offset = offset*recl
       call MPI_file_write_at(File%fh,offset,variable(ie),recl,MPI_BYTE,status,ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,__FILE__,__LINE__, &
                        ierr,errorstring(1:errorlen)
       endif
     enddo
endif

     call CloseFile(File%fh)
#else
     if(iam == 1 ) then
       open(unit=56,file=File%fname,status='UNKNOWN',form='UNFORMATTED', &
	     recl=recl,ACCESS='DIRECT')
     endif
#ifdef _MPI
     call syncmp(File%par)
     if (iam /= 1) then 
       open(unit=56,file=File%fname,status='OLD',form='UNFORMATTED', &
             recl=recl,ACCESS='DIRECT')

     endif   
#endif
     do ie=1,nelemd
#ifdef _PREDICT
       ig = Schedule(iam)%Local2Global(ie)
#else
       ig = Schedule(1)%Local2Global(ie)
#endif
!	print *, __FILE__,__LINE__,ie,ig,sizeof(variable(ie))

        if(columnpackage == "emanuel")then
   	  write(56,rec=ig) variable(ie)
        else
	  write(56,rec=ig) variable(ie)%buffer
	endif
     enddo
     close(56)
#endif
     call t_stopf('WriteState')

   end subroutine WriteState
!=============================
! ReadState:
!
!  Reads a field out to disk
! ============================
   subroutine ReadState(File,variable,recl)

     type (File_elem_t),intent(in)                      :: File
     type (RestartBuffer_t),target,intent(out)     :: variable(:)
     integer, intent(in)  :: recl

     !==================
     ! Local variables
     !==================
     integer           :: ie, ig, info, amode

#if defined(_MPI) && defined(_PRESTART)
     integer (kind=int_kind)    :: type_ablock, type_fview
     integer                    :: array_disp(nelemd)
     integer                    :: array_blen(nelemd)
     integer*4                  :: acount, isiz
     integer (kind=mpi_offset_kind) :: disp, iext,lb
#endif

     call t_startf('ReadState')

#if defined(_MPI) && defined(_PRESTART)
if ( COLLECTIVE_IO_READ ) then
     call mpi_type_contiguous( recl, MPI_BYTE, type_ablock, ierr )
     call mpi_type_commit( type_ablock, ierr )

     do ie=1,nelemd
#ifdef _PREDICT
       ig = Schedule(iam)%Local2Global(ie)
#else
       ig = Schedule(1)%Local2Global(ie)
#endif
       array_disp(ie) = ig-1
       array_blen(ie) = 1

     enddo

     disp = array_disp(1)
     do ie=1,nelemd
       array_disp(ie) = array_disp(ie) - disp
     enddo

!    call mpi_type_indexed( nelemd, array_blen, array_disp, type_ablock, type_fview, ierr )
     call mpi_type_create_indexed_block( nelemd, 1, array_disp, type_ablock, type_fview, ierr )
     call mpi_type_commit( type_fview, ierr )
     call mpi_type_get_extent( type_fview, iext,lb, ierr )
     call mpi_type_size( type_fview, isiz, ierr )
#if 0
     print *, 'debug: array_disp(1)=', array_disp(1), ', recl=', recl, ', iext=', iext, ', isiz=', isiz
#endif
endif

     amode = MPI_MODE_RDONLY     
     info  = MPI_INFO_NULL
     call MPI_file_open(File%par%comm,File%fname,amode,info,File%fh,ierr)
     if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,__FILE__,__LINE__,errorstring(1:errorlen)
     endif

if ( COLLECTIVE_IO_READ ) then
     disp = disp * recl
     call MPI_File_set_view( File%fh, disp, type_ablock, type_fview, "native", MPI_INFO_NULL, ierr )
     call mpi_file_read_all( File%fh, variable, nelemd, type_ablock, status, ierr )
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,__FILE__,__LINE__, &
                        ierr,errorstring(1:errorlen), offset
       endif
     
else
     do ie=1,nelemd
#ifdef _PREDICT
       ig = Schedule(iam)%Local2Global(ie)
#else
       ig = Schedule(1)%Local2Global(ie)
#endif
       ! offset = (ig-1)*recl   this will overflow integer*4 at high-res
       offset = ig-1
       offset = offset*recl
       call MPI_file_read_at(File%fh,offset,variable(ie),recl,MPI_BYTE,status,ierr)

       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'ReadState:  After MPI_file_read_at(variable(i)):  ', &
                        ierr,errorstring(1:errorlen)
       endif
     enddo
endif
     call CloseFile(File%fh)
#else
     open(unit=56,file=File%fname,status='OLD',form='UNFORMATTED', &
	     recl=recl,ACCESS='DIRECT')
     do ie=1,nelemd
#ifdef _PREDICT
       ig = Schedule(iam)%Local2Global(ie)
#else
       ig = Schedule(1)%Local2Global(ie)
#endif
         if(columnpackage == "emanuel")then
           read(56,rec=ig) variable(ie)
	 else
           read(56,rec=ig) variable(ie)%buffer
         endif
     enddo
     close(56)
#endif
     call t_stopf('ReadState')

   end subroutine ReadState
!==========================================================
! CreateRestartHeader:
!
!  Create the restart header to write to the restart file.
!==========================================================
   subroutine CreateRestartHeader(header, TimeLevel)
    implicit none

    type (RestartHeader_t),intent(inout) :: header
    type (TimeLevel_t), intent(in) :: TimeLevel
    if (0==ne) call abortmp('Error in CreateRestartHeader: ne is zero')

    header%version   = RestartVersion
    header%ne        = ne
    header%np        = np
    header%nlev      = nlev
    header%number    = 0 ! reserved for future use
    header%TimeLevel = TimeLevel
!
! one integer marks the element number
!
    header%ElemRecLength   = nwordsRestartBuffer_t*real_kind
!
! The intrinsic sizeof is not (yet) a part of the F90 standard.  The code is left as 
! documentation of the hardcoded value RESTART_HDR_CNT.
! 
! only used by old, serial restart code
#ifndef _PRESTART
#ifdef _AIX
    if(RESTART_HDR_CNT*int_kind .ne. sizeof(header)) then
	call haltmp('bad restart header size')
    endif
#endif
#endif
    header%size       = RESTART_HDR_CNT*int_kind

   end subroutine CreateRestartHeader
!======================================================
! CheckRestartHeader:
!
!     Check to make sure the restart file you are
!       trying to use matches the models resolution
!======================================================
   function CheckRestartHeader(header) result(ierr)
    implicit none

    integer                              :: ierr
    type (RestartHeader_t),intent(inout) :: header


    if (0==ne) call abortmp('Error in CreateRestartHeader: ne is zero')
    ierr = 0
    if( header%version .ne. RestartVersion) then
      write(*,*) 'CheckRestartHeader: Restart tape VERSION does not match model (', &
                header%version,' != ',RestartVersion,')'
      ierr = -1
    endif

    if( header%ne .ne. ne) then
      write(*,*) 'CheckRestartHeader: Restart tape NE does not match model (', &
                header%ne,' != ',ne,')'
      ierr = -1
    endif

    if( header%np .ne. np) then
      write(*,*) 'CheckRestartHeader: Restart tape NP does not match model (', &
                header%np,' != ',np,')'
      ierr = -1
    endif

    if( header%nlev .ne. nlev) then
      write(*,*) 'CheckRestartHeader: Restart tape NLEV does not match model (', &
                header%nlev,' != ',nlev,')'
      ierr = -1
    endif

   end function CheckRestartHeader
!==============================================================================
!==============================================================================
   subroutine WriteRestartHeader(File,Header)

     type (File_elem_t), intent(in)          :: File
     type (RestartHeader_t),intent(in)  :: Header

     integer                   :: info,len,amode
     character(len=10)         :: datarep
     character(len=84)               :: headername

     headername = trim(File%fname)//'.hdr'
     
     ! =============================================
     !  Load the information into the Restart Header
     ! =============================================
     offset = 0

#if defined(_MPI) && defined(_PRESTART)
     amode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
     info  = MPI_INFO_NULL
     call MPI_file_open(File%par%comm,headername,amode,info,File%fh,ierr)
     if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,__FILE__,__LINE__,errorstring(1:errorlen)
	call abortmp(headername)
     endif

     info  = MPI_INFO_NULL
     datarep = 'native'

     call MPI_File_set_view(File%fh,offset,MPIinteger_t,MPIinteger_t,datarep,info,ierr)
     if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'WriteRestartHeader: After MPI_file_set_view: ',errorstring(1:errorlen)
     endif

     if(File%par%masterproc) then
        ! =========================================================
        !  If you are the master process, write the restart header
        ! =========================================================
        len = Header%size/int_kind
        call MPI_file_write(File%fh,Header,len,MPIinteger_t,status,ierr)
        if(ierr .ne. MPI_SUCCESS) then
           errorcode=ierr
           call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
           print *,'WriteRestartHeader: After MPI_file_write_at: ',errorstring(1:errorlen)
        endif
     endif
     call syncmp(File%par)
     call CloseFile(File%fh)
#else
    if(iam.eq.1) then
      open(unit=57,file=headername,status='NEW',form='UNFORMATTED', &
	     recl=(int_kind*RESTART_HDR_CNT),ACCESS='DIRECT')

      write(unit=57,rec=1) Header
      close(57)
    endif
#endif

   end subroutine WriteRestartHeader
!=============================================================================================
!
!=============================================================================================
   subroutine ReadRestartHeader(File,Header,TimeLevel)

     type (File_elem_t),intent(in)  :: File
     type (RestartHeader_t),intent(out) :: Header
     type (timelevel_t), intent(out) :: TimeLevel
     integer                         :: info,len,amode
     character(len=10)               :: datarep
     character(len=84)               :: headername

     headername = trim(File%fname)//'.hdr'

#if defined(_MPI) && defined(_PRESTART)
     amode = MPI_MODE_RDONLY     
     info  = MPI_INFO_NULL
     call MPI_file_open(File%par%comm,headername,amode,info,File%fh,ierr)
     if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,__FILE__,__LINE__,errorstring(1:errorlen)
     endif



     info  = MPI_INFO_NULL

     ! =============================================
     !  Load the information into the Restart Header
     ! =============================================
     offset = 0
     datarep = 'native'
     call MPI_File_set_view(File%fh,offset,MPIinteger_t,MPIinteger_t,datarep,info,ierr)
     if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'ReadRestartHeader: After MPI_file_set_view: ',errorstring(1:errorlen)
     endif

     ! ======================================================================
     !  Read the first two words off the header, so see how large things are
     ! ======================================================================
     len = 2
     call MPI_file_read_at(File%fh,offset,Header,len,MPIinteger_t,status,ierr)

     ! =========================================================
     !  Now, read the complete restart header
     ! =========================================================
     len = (Header%size/int_kind) 
     call MPI_file_read_at(File%fh,offset,Header,len,MPIinteger_t,status,ierr)
     if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'ReadRestartHeader: After MPI_file_write_at: ',errorstring(1:errorlen)
     endif
     ierr=CheckRestartHeader(Header)
     if(ierr .ne. 0) then
        call haltmp("ReadRestartHeader: Error with restart file")
     endif
     offset = Header%size 
!     call MPI_File_set_view(File%fh,offset,File%ElemState,File%GlobalState,datarep,info,ierr)
!     if(ierr .ne. MPI_SUCCESS) then
!        errorcode=ierr
!        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
!        print *,'ReadRestartHeader: After MPI_file_set_view: ',errorstring(1:errorlen)
!     endif
     call CloseFile(File%fh)
#else
    if(iam.eq.1) then
       open(unit=57,file=headername,status='OLD',form='UNFORMATTED', &
	     recl=(int_kind*RESTART_HDR_CNT),ACCESS='DIRECT')

       read(57,rec=1) Header 
       close(57)
    endif
#ifdef _MPI
    call MPI_bcast(Header,RESTART_HDR_CNT,MPIinteger_t,File%par%root,File%par%comm,ierr)
#endif
    
#endif

    TimeLevel = Header%TimeLevel
    nEndStep = nmax + TimeLevel%nstep

   end subroutine ReadRestartHeader
!========================================================================================
!
!========================================================================================
    subroutine WriteRestart(elem, ithr,nets,nete,tl)
    type(element_t), intent(in)      :: elem(:)       
    integer,            intent(in)         :: ithr,nete,nets
    type (TimeLevel_t), intent(in)         :: tl    ! time level struct

    integer                                :: ie
    character(len=186)                      :: fname
    character(len=9)                       :: charnum
    type (RestartHeader_t)                 :: RestartHeader
 
    logical, parameter                     :: Debug = .FALSE.


    do ie=nets,nete
       RestartBuffer(ie)%buffer = elem(ie)%state
#ifdef _PRIMXXX
       if(columnpackage == "emanuel") then
         RestartBuffer(ie)%puffer = pelem(ie)%state
       endif
#endif
    enddo

!$OMP BARRIER
!$OMP MASTER
    write(charnum,'(i9.9)') tl%nstep
    if(iam .eq. 1) then 
       if(Debug) print *,'WriteRestart: restnum is: ',charnum
    endif

    fname = TRIM(ADJUSTL(restartdir))//"/R"//TRIM(ADJUSTL(charnum))

    RestFile%fname=fname

    call CreateRestartHeader(RestartHeader,tl)
    call WriteRestartHeader(RestFile,RestartHeader)

    call WriteState(RestFile,RestartBuffer,RestartHeader%ElemRecLength)
!$OMP END MASTER

!DBG    print *,'WriteRestart: point #10'
    end subroutine WriteRestart

 
    subroutine ReadRestart(elem, ithr,nets,nete,tl)
    implicit none    
    type (element_t),   intent(inout)     :: elem(:)
    integer,            intent(in)         :: ithr,nets,nete
    type (TimeLevel_t), intent(out)            :: tl     ! time level struct
    type (RestartHeader_t)                 :: RestartHeader
    integer :: ie
!   MT: adding barrier here, since we will change tl, so need to make
!       sure all threads are done using it
!$OMP BARRIER
!$OMP MASTER
         ! ===================================
         ! Read the restart info off the file
         ! ===================================
    RestFile%fname=restartfile
    call ReadRestartHeader(RestFile,RestartHeader,tl)
    call ReadState(RestFile,RestartBuffer,RestartHeader%ElemRecLength)
!$OMP END MASTER
!$OMP BARRIER
    do ie=nets,nete
       elem(ie)%state = RestartBuffer(ie)%buffer
#ifdef _PRIMXXX
       if(columnpackage == "emanuel")then
          pelem(ie)%state = RestartBuffer(ie)%puffer
       endif
#endif
    enddo


    end subroutine ReadRestart

!=============================
! CloseFile:
!
!  Close a disk file
! ============================
   subroutine CloseFile(fh)

   integer,intent(in)   :: fh

#if defined(_MPI) && defined(_PRESTART)
     call MPI_file_close(fh,ierr)
     if(ierr .ne. MPI_SUCCESS) then
       errorcode=ierr
       call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
       print *,'CloseFile: error with MPI_file_close: ',errorstring(1:errorlen)
     endif
#else
     close(unit=fh)
#endif

   end subroutine CloseFile
! =========================================================
! ConstructElementFile
!
!  Initalizes MPI I-O  File structures 
! =========================================================
    subroutine ConstructElementFile(descriptor,File,ierr)

    type (StateDesc_t),intent(in) :: descriptor
    type (File_elem_t),intent(out) :: File
    integer (kind=int_kind),intent(out)  :: ierr

    integer,parameter            :: ndims = 1
    integer                      :: sizes(ndims), &
                                    subsizes(ndims), &
                                    displace(ndims)

    integer,allocatable          :: elem_blklengths(:), &
                                    elem_buftype(:), &
				    elem_filetype(:)
    integer(kind=MPI_ADDRESS_KIND),allocatable          :: elem_disp(:)
    integer (kind=int_kind)      :: ie,ig
 

#if defined(_MPI) && defined(_PRESTART)
    if (collective_io_read .or. collective_io_write ) then
       ! we need to setup collective MPI IO structs and data types
    else
       return
    endif

    call PrintStateDescriptor(descriptor)

    call mpi_type_contiguous( descriptor%nwords*real_kind, MPI_BYTE, File%ElemState, ierr )
!    call MPI_type_create_Struct(descriptor%NumComponents,descriptor%blklen, &
!	 descriptor%disp,descriptor%type,File%ElemState,ierr)
    if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'ConstructStateDescriptor:  After MPI_type_Struct for [ElemState] ', &
                        ierr,errorstring(1:errorlen)
    endif

    call MPI_Type_Commit(File%ElemState,ierr)
    if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'ConstructStateDescriptor:  After MPI_type_Commit for [ElemState] ', &
                        ierr,errorstring(1:errorlen)
    endif
#ifdef DEBUG_IO
    print *,'ConstructStateDescriptor: After creation of [ElemState]'
    call PrintTypeInfo(File%ElemState)
#endif
    ! =========================
    !  Allocate working arrays 
    ! =========================
    allocate(elem_blklengths(nelemd))
    allocate(elem_disp(nelemd))
    allocate(elem_buftype(nelemd))
    allocate(elem_filetype(nelemd))

    ! ========================================================
    !  Setup the buffer structure  
    ! ========================================================
    sizes(1)   = 1
    sizes(1)   = nelemd
    subsizes(1) = 1
    subsizes(1) = 1
    displace(1)  = 0

    do ie=1,nelemd
       displace(1) = ie - 1
       call MPI_Type_Create_Subarray(ndims,sizes,subsizes, &
             displace,MPI_ORDER_FORTRAN,File%ElemState,elem_buftype(ie),ierr)
       elem_blklengths(ie)    = 1
       elem_disp(ie)          = 0
       if (ierr.ne.MPI_SUCCESS) call abortmp('Error in MPI_Type_create_subarray')
    enddo


    call MPI_Type_create_Struct(nelemd,elem_blklengths, &
             elem_disp,elem_buftype,File%State,ierr)
    if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'ConstructStateDescriptor:  After MPI_type_Struct for [State] ', &
                        ierr,errorstring(1:errorlen)
        call abortmp(errorstring)
    endif
    call MPI_Type_Commit(File%State,ierr)
    if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'ConstructStateDescriptor:  After MPI_type_Commit for [ElemState] ', &
                        ierr,errorstring(1:errorlen)
    endif
#ifdef DEBUG_IO
    call PrintTypeInfo(File%State)
#endif
    ! ========================================================
    !  Setup the file structure  
    ! ========================================================
    sizes(1)    = 1
    sizes(1)    = nelem
    subsizes(1) = 1
    subsizes(1) = 1
    displace(1) = 0

    do ie=1,nelemd
#ifdef _PREDICT
       ig = Schedule(iam)%Local2Global(ie)
#else
       ig = Schedule(1)%Local2Global(ie)
#endif
       displace(1) = ig - 1
       call MPI_Type_Create_Subarray(ndims,sizes,subsizes,displace, &
             MPI_ORDER_FORTRAN,File%ElemState,elem_filetype(ie),ierr)
       elem_blklengths(ie)  = 1
       elem_disp(ie)        = 0
    enddo

    call MPI_type_create_Struct(nelemd,elem_blklengths, elem_disp, &
             elem_filetype, File%GlobalState,ierr)
    if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'ConstructStateDescriptor:  After MPI_type_Struct for [GlobalState] ', &
                        ierr,errorstring(1:errorlen)
    endif
    call MPI_Type_Commit(File%GlobalState,ierr)
    if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'ConstructStateDescriptor:  After MPI_type_Commit for [GlobalState] ', &
                        ierr,errorstring(1:errorlen)
    endif
#ifdef DEBUG_IO
    call PrintTypeInfo(File%GlobalState)
#endif
    ! ============================
    !  Deallocate working arrays 
    ! ============================
    deallocate(elem_blklengths)
    deallocate(elem_disp)
    deallocate(elem_buftype)
    deallocate(elem_filetype)
#endif

   end subroutine ConstructElementFile

!==========================================================================
!  AddStateField:
!     Adds a statefield to the descriptor
!==========================================================================
    subroutine AddStateField(desc,len,type)

    type (StateDesc_t),intent(inout)   :: desc
    integer (kind=int_kind),intent(in) :: len,type

    integer  :: i
    
    i = desc%ptr
	
    desc%blklen(i)=len
    desc%type(i)=type

    if(type .eq. MPIinteger_t) desc%kind(i)=4
    if(type .eq. MPIReal_t)    desc%kind(i)=8

    if (i==1) then
       desc%disp(i) = 0  ! should this be 0 or 1?
    else
       desc%disp(i) = desc%disp(i-1) + desc%kind(i-1)*desc%blklen(i-1)
    endif

    desc%nwords = desc%nwords + len   ! record the number of words 
    desc%ptr=i+1

 
    end subroutine AddStateField
#if defined(_MPI) && defined(_PRESTART)
!=========================================================================
!  PrintTypeInfo:
!    Prints out information about an MPI datatype 
!=========================================================================
  subroutine PrintTypeInfo(type)

    integer (kind=int_kind) :: type 
    integer (kind=int_kind) :: size
    integer (kind=MPI_ADDRESS_KIND) :: lb,extent
        
    print *,'PrintTypeInfo: MPI_datatype # is ',type
    call MPI_Type_Size(type,size,ierr)
    if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'PrintTypeInfo:  After MPI_type_Size', &
                        ierr,errorstring(1:errorlen)
    else
        print *,'PrintTypeInfo: SIZE in bytes of datatype is ',size
    endif

    call MPI_Type_Get_Extent(type,lb,extent,ierr)
    if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'PrintTypeInfo:  After MPI_Type_Get_Extent', &
                        ierr,errorstring(1:errorlen)
    else
        print *,'PrintTypeInfo: LB and EXTENT are  ',lb,extent
    endif

  end subroutine PrintTypeInfo
#endif
!=========================================================================
!  PrintStateDescriptor
!     Prints out the state descriptor
!=========================================================================
    subroutine PrintStateDescriptor(desc)
	type (StateDesc_t), intent(in) :: desc
     integer (kind=int_kind) :: numComp 
     integer (kind=int_kind) :: i

   if(iam .eq. 1) then 
       print *,'PrintStateDescriptor:' 
       numComp = desc%NumComponents
       do i=1,numComp
	  print 100,i,desc%blklen(i),desc%disp(i)
! ,desc%type(i)
       enddo
    endif 
100 format('component ',i3,' length=',i5,' displacement=',i7)
    end subroutine PrintStateDescriptor
!==========================================================================
!  CreateStateDescriptor:
!     Creates a descriptor to make the formation of Restart buffers easier
!===========================================================================
    function CreateStateDescriptor(numComp) result(desc)
    
    integer (kind=int_kind), intent(in) :: numComp
    type (StateDesc_t) :: desc
  
    integer (kind=int_kind) :: len
	
    len = numComp
    allocate(desc%blklen(len))
    allocate(desc%disp(len))
    allocate(desc%kind(len))
    allocate(desc%type(len))

    desc%blklen(:)=0
    desc%disp(:)=0
    desc%kind(:)=0
    desc%type(:)=0

    desc%nwords = 0
    desc%ptr=1
    desc%numComponents = numComp

    end function CreateStateDescriptor
!===============================================================


end module restart_io_mod
