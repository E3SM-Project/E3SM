#define _FILE_ "piolib_mod.f90"
#define debug_rearr 0
#ifdef BGP
#define BGx
#endif
#ifdef BGL
#define BGx
#endif

module piolib_mod
  !--------------
  use pio_kinds
  !--------------
  use pio_types, only : file_desc_t, iosystem_desc_t, var_desc_t, io_desc_t, &
	pio_iotype_pbinary, pio_iotype_binary, pio_iotype_direct_pbinary, &
	pio_iotype_netcdf, pio_iotype_pnetcdf, pio_iotype_netcdf4p, pio_iotype_netcdf4c, &
        pio_noerr
  !--------------
  use alloc_mod
  !--------------
  use pio_support, only : piodie, debug, debugio, debugasync, checkmpireturn
  !
  use ionf_mod, only : create_nf, open_nf,close_nf, sync_nf
  use pionfread_mod, only : read_nf
  use pionfwrite_mod, only : write_nf

  use pio_mpi_utils, only : PIO_type_to_mpi_type 
  use iompi_mod
  use rearrange
#ifdef TIMING
  use perf_mod, only : t_startf, t_stopf     ! _EXTERNAL
#endif
  use pio_msg_mod

  implicit none
  private
  save

  include 'mpif.h'    ! _EXTERNAL

  ! !public member functions:

  public :: PIO_init,     &
       PIO_finalize,      &
       PIO_initdecomp,    &
       PIO_openfile,      &
       PIO_syncfile,      &
       PIO_createfile,    &
       PIO_closefile,     &
       PIO_setiotype,     &
       PIO_numtoread,     &
       PIO_numtowrite,    &
       PIO_setframe,      &
       PIO_advanceframe,  &
       PIO_setdebuglevel, &
       PIO_seterrorhandling, &
       PIO_get_local_array_size, &
       PIO_freedecomp,     &
       PIO_dupiodesc,     &
       PIO_getnumiotasks, &
       PIO_set_hint,      &
       PIO_FILE_IS_OPEN

#ifdef MEMCHK
!> this is an internal variable for memory leak debugging 
!! it is used when macro memchk is defined and it causes each task to print the 
!! memory resident set size anytime it changes within pio.
!<
  integer :: lastrss=0
#endif

  !eop
  !boc
  !-----------------------------------------------------------------------
  !
  !  module variables
  !
  !-----------------------------------------------------------------------
!> 
!! @defgroup PIO_openfile PIO_openfile
!< 
  interface PIO_openfile
     module procedure PIO_openfile
  end interface

!> 
!! @defgroup PIO_syncfile PIO_syncfile
!<
  interface PIO_syncfile
     module procedure syncfile
  end interface

!> 
!! @defgroup PIO_createfile PIO_createfile
!<
  interface PIO_createfile
     module procedure createfile
  end interface

!> 
!! @defgroup PIO_setframe PIO_setframe
!! @brief sets the unlimited dimension for netcdf file for record number for binary files
!<
  interface PIO_setframe
     module procedure setframe
  end interface

!> 
!! @defgroup PIO_advanceframe PIO_advanceframe
!<
  interface PIO_advanceframe
     module procedure advanceframe
  end interface

!> 
!! @defgroup PIO_closefile PIO_closefile
!<
  interface PIO_closefile
     module procedure closefile
  end interface


!> 
!! @defgroup PIO_freedecomp PIO_freedecomp
!! free memory associated with a io descriptor
!<
  interface PIO_freedecomp
     module procedure freedecomp_ios
     module procedure freedecomp_file
  end interface

!> 
!! @defgroup PIO_init PIO_init
!! initializes the pio subsystem
!<
  interface PIO_init
     module procedure init_intracom
     module procedure init_intercom
     
  end interface

!> 
!! @defgroup PIO_finalize PIO_finalize
!! Shuts down and cleans up any memory associated with the pio library.
!<
  interface PIO_finalize
     module procedure finalize
  end interface

!>
!! @defgroup PIO_initdecomp PIO_initdecomp
!! @brief PIO_initdecomp is an overload interface the models decomposition to pio.
!<
  interface PIO_initdecomp
     module procedure PIO_initdecomp_dof  ! previous name: initdecomop_1dof_nf_box
     module procedure PIO_initdecomp_bc
     module procedure PIO_initdecomp_dof_dof
     module procedure initdecomp_1dof_nf
     module procedure initdecomp_1dof_bin
     module procedure initdecomp_2dof_nf
     module procedure initdecomp_2dof_bin
  end interface

!> 
!! @defgroup PIO_dupiodesc PIO_dupiodesc
!! duplicates an eisting io descriptor
!<
  interface PIO_dupiodesc
     module procedure dupiodesc
  end interface

!> 
!! @defgroup PIO_setiotype PIO_setiotype
!!  sets the io type used by pio
!<
  interface PIO_setiotype 
     module procedure setiotype
  end interface

!> 
!! @defgroup PIO_numtoread PIO_numtoread
!! returns the total number of words to read
!<
  interface PIO_numtoread
     module procedure numtoread
  end interface

!> 
!! @defgroup PIO_numtowrite PIO_numtowrite
!! returns the total number of words to write
!<
  interface PIO_numtowrite
     module procedure numtowrite
  end interface


!> 
!! @defgroup PIO_getnumiotasks PIO_getnumiotasks
!!  returns the actual number of IO-tasks used.  PIO 
!!  will reset the total number of IO-tasks if certain 
!!  conditions are meet
!<
  interface PIO_getnumiotasks
     module procedure getnumiotasks
  end interface

!> 
!!  @defgroup PIO_setdebuglevel PIO_setdebuglevel
!!  sets the level of debug information that pio will generate.
!<
  interface PIO_setdebuglevel
     module procedure setdebuglevel
  end interface

!> 
!!  @defgroup PIO_seterrorhandling PIO_seterrorhandling
!!  sets the form of error handling for pio.
!!
!! By default pio handles errors internally by printing a string
!! describing the error and calling mpi_abort.  Application
!! developers can change this behavior for calls to the underlying netcdf
!! libraries with a call to PIO_seterrorhandling. For example if a
!! developer wanted to see if an input netcdf format file contained the variable
!! 'u' they might write the following
!! @verbinclude errorhandle
!<
  interface PIO_seterrorhandling
     module procedure seterrorhandlingf
     module procedure seterrorhandlingi
  end interface

!>
!! @defgroup PIO_get_local_array_size PIO_get_local_array_size
!<

  !eoc
  !***********************************************************************

  private :: nextlarger

contains

  logical function PIO_FILE_IS_OPEN(File)
    type(file_desc_t), intent(in) :: file
    pio_file_is_open = file%file_is_open
  end function PIO_FILE_IS_OPEN


!> 
!! @public 
!! @ingroup PIO_get_local_array_size
!! @brief This function returns the expected local size of an array associated with iodesc
!! @details
!! @param iodesc @copydoc io_desc_t
!<
  integer function PIO_get_local_array_size(iodesc)
    type(io_desc_t), intent(in) :: iodesc   
    PIO_get_local_array_size = iodesc%compsize
  end function PIO_get_local_array_size

!> 
!! @public 
!! @ingroup PIO_advanceframe
!! @brief advances the record dimension of a variable in a netcdf format file 
!!  or the block address in a binary file
!! @details
!! @param vardesc @copydoc var_desc_t
!<
  subroutine advanceframe(vardesc)
    type(var_desc_t), intent(inout) :: vardesc
    vardesc%rec=vardesc%rec+1
  end subroutine advanceframe

!> 
!! @public 
!! @ingroup PIO_setframe 
!! @brief sets the record dimension of a variable in a netcdf format file 
!! or the block address in a binary file
!! @details
!! @param vardesc @copydoc var_desc_t
!! @param frame   : frame number to set
!<
  subroutine setframe(vardesc,frame)
    type(var_desc_t), intent(inout) :: vardesc
    integer(kind=PIO_offset), intent(in) :: frame
    vardesc%rec=frame
  end subroutine setframe

!>  
!! @public
!! @ingroup PIO_setdebuglevel
!! @brief sets the level of debug information output to stdout by pio 
!! @details
!! @param level : default value is 0, allowed values 0-3
!<
  subroutine setdebuglevel(level)
    integer(i4), intent(in) :: level	
    if(level.eq.0) then
       debug=.false.
       debugio=.false.
       debugasync=.false.
    else if(level.eq.1) then
       debug=.true.
       debugio=.false.
       debugasync=.false.
    else if(level.eq.2) then
       debug=.false.
       debugio=.true.
       debugasync=.false.
    else if(level.eq.3) then
       debug=.true.
       debugio=.true.
       debugasync=.false.
    else if(level.eq.4) then
       debug=.false.
       debugio=.false.
       debugasync=.true.
    else if(level.eq.5) then
       debug=.true.
       debugio=.false.
       debugasync=.true.
    else if(level.ge.6) then
       debug=.true.
       debugio=.true.
       debugasync=.true.
    
    end if
  end subroutine setdebuglevel

!>
!! @ingroup PIO_seterrorhandling
!! @public
!! @brief set the pio error handling method for a file
!!
!! @param file @copydoc file_desc_t
!! @param method :
!! @copydoc PIO_error_method
!<
  subroutine seterrorhandlingf(file, method)
    type(file_desc_t), intent(inout) :: file
    integer, intent(in) :: method

    call seterrorhandlingi(file%iosystem, method)
  end subroutine seterrorhandlingf

!>
!! @ingroup PIO_seterrorhandling 
!! @public
!! @brief set the pio error handling method for the iosystem
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param method :
!! @copydoc PIO_error_method
!<
  subroutine seterrorhandlingi(ios, method)
    use pio_types, only : pio_internal_error, pio_return_error
    use pio_msg_mod, only : pio_msg_seterrorhandling
    type(iosystem_desc_t), intent(inout) :: ios
    integer, intent(in) :: method
    integer :: msg, ierr

    if(ios%async_interface .and. .not. ios%ioproc ) then
       msg=PIO_MSG_SETERRORHANDLING
       if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
       call MPI_BCAST(method,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , ierr)
    end if
    if(Debugasync) print *,__FILE__,__LINE__,method
    ios%error_handling=method

    if(method > PIO_internal_error .or. method < PIO_return_error) then
       call piodie(_FILE_,__LINE__,'invalid error handling method requested')
    end if
  end subroutine seterrorhandlingi

!> 
!! @public 
!! @ingroup PIO_initdecomp
!! @brief Implements the @ref decomp_bc for PIO_initdecomp
!! @details  This provides the ability to describe a computational 
!! decomposition in PIO that has a block-cyclic form.  That is 
!! something that can be described using start and count arrays.
!! Optional parameters for this subroutine allows for the specification
!! of io decomposition using iostart and iocount arrays.  If iostart
!! and iocount arrays are not specified by the user, and rearrangement 
!! is turned on then PIO will calculate a suitable IO decomposition
!! @param iosystem @copydoc iosystem_desc_t
!! @param basepiotype @copydoc use_PIO_kinds
!! @param dims An array of the global length of each dimesion of the variable(s)
!! @param compstart The start index into the block-cyclic computational decomposition
!! @param compcount The count for the block-cyclic computational decomposition
!! @param iodesc @copydoc iodesc_generate
!! @param iostart   The start index for the block-cyclic io decomposition
!! @param iocount   The count for the block-cyclic io decomposition
!<
  subroutine PIO_initdecomp_bc(iosystem,basepiotype,dims,compstart,compcount,iodesc,iostart,iocount)
    type (iosystem_desc_t), intent(inout) :: iosystem
    integer(i4), intent(in)               :: basepiotype
    integer(i4), intent(in)               :: dims(:)
    integer (kind=PIO_OFFSET)             :: compstart(:)  
    integer (kind=PIO_OFFSET)             :: compcount(:)    
    type (IO_desc_t), intent(out)         :: iodesc
    integer (kind=PIO_OFFSET),optional    :: iostart(:)  
    integer (kind=PIO_OFFSET),optional    :: iocount(:)    

!    character(len=*), parameter :: '::PIO_initdecomp_bc'

    call piodie(__FILE__,__LINE__,'subroutine not yet implemented')

  end subroutine PIO_initdecomp_bc

!>
!! @public
!! @ingroup PIO_initdecomp
!! @brief Implements the @ref decomp_dof for PIO_initdecomp
!! @details  This provides the ability to describe a computational
!! decomposition in PIO using degrees of freedom method. This is  
!! a decomposition that can not be easily described using a start  
!! and count metehod (see @ref decomp_dof).  This subroutine also 
!! requires the user to specify the IO decomposition using the 
!! degree of freedom method.  This version of the subroutine 
!! is most suitable for those who want complete control over 
!! the actions of PIO.
!! @param iosystem @copydoc iosystem_desc_t
!! @param basepiotype @copydoc use_PIO_kinds
!! @param dims An array of the global length of each dimesion of the variable(s)
!! @param compdof Mapping of the storage order for the computatinal decomposition to its memory order
!! @param iodesc @copydoc iodesc_generate
!! @param iodof Mapping of the storage order for the IO decomposition its memory order
!<
  subroutine PIO_initdecomp_dof_dof(iosystem,basepiotype,dims,compdof,iodesc,iodof)
    type (iosystem_desc_t), intent(inout)          :: iosystem
    integer(i4), intent(in)                        :: basepiotype
    integer(i4), intent(in)                        :: dims(:)
    integer(i4), intent(in)                        :: compdof(:)
    type (IO_desc_t), intent(out)                   :: iodesc
    integer(i4), intent(in)                        :: iodof(:)

!    character(len=*), parameter :: subName=modName//'::PIO_initdecomp_dof_dof'

!    call piodie(subname,__LINE__,'subroutine not yet implemented')

  end subroutine PIO_initdecomp_dof_dof

!> 
!! @public 
!! @ingroup PIO_initdecomp
!! @brief A deprecated interface to the PIO_initdecomp method.
!! @details
!! @deprecated
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param basepiotype : the type of variable(s) associated with this iodesc.  
!! @copydoc PIO_kinds
!! @param dims : an array of the global length of each dimesion of the variable(s)
!! @param lenblocks :
!! @param compdof : mapping of the storage order of the variable to its memory order
!! @param iodofr :
!! @param iodofw :
!! @param iodesc @copydoc iodesc_generate
!<
  subroutine initdecomp_2dof_bin(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodofw,iodesc)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4)                       :: basetype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenblocks
    integer (i4), intent(in)          :: compdof(:)   !> global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: iodofr(:)     !> global degrees of freedom for io decomposition 
    integer (i4), intent(in)          :: iodofw(:)     !> global degrees of freedom for io decomposition 
    type (io_desc_t), intent(out)     :: iodesc

    integer(kind=PIO_offset) :: start(1), count(1)

    integer (i4) :: i,ndims,glength,n_iotasks
    logical :: userearranger
    integer (i4) ::  ndispr,ndispw
    integer (i4) :: lengthr, lengthw
    integer (i4), pointer :: displacer(:),displacew(:)

    basetype=PIO_type_to_mpi_type(basepiotype)

    !-------------------------------------------
    ! for testing purposes set the iomap
    ! (decompmap_t) to something basic for
    ! testing.
    !-------------------------------------------
    userearranger = iosystem%userearranger

    !---------------------
    ! number of dimensions
    !---------------------
    ndims = size(dims)
    !---------------------
    ! total global size
    !---------------------
    glength=1
    do i=1,ndims
       glength = glength*dims(i)
    enddo
    lengthr = size(iodofr);
    lengthw = size(iodofw)
    if(lenblocks>0) then
       ndispw=size(iodofw)/lenblocks 
       ndispr=size(iodofr)/lenblocks
    else
       ndispw=0
       ndispr=0
    end if
    call alloc_check(displacer,ndispr)
    call alloc_check(displacew,ndispw)

    !--------------------------------------------
    ! calculate mpi data structure displacements
    !--------------------------------------------
    !dbg    print *,'PIO_initdecomp: before call to calcdisplace'
    if(lenblocks>0) then
       call calcdisplace(lenblocks,iodofr,displacer)
       call calcdisplace(lenblocks,iodofw,displacew)
    end if
    n_iotasks = iosystem%num_iotasks

    iodesc%glen = glength

    if(debug) print *,'iam: ',iosystem%io_rank,'initdecomp: userearranger: ',userearranger

    !---------------------------------------------
    !  the setup for the mpi-io type information
    !---------------------------------------------
    if(iosystem%ioproc) then
       !-----------------------------------------------
       ! setup the data structure for the read operation
       !-----------------------------------------------
       iodesc%read%n_elemtype = ndispr
       iodesc%read%n_words    = iodesc%read%n_elemtype*lenblocks
       call genindexedblock(lenblocks,basetype,iodesc%read%elemtype,iodesc%read%filetype,displacer)

       !-------------------------------------------------
       ! setup the data structure for the write operation
       !-------------------------------------------------
       iodesc%write%n_elemtype = ndispw
       iodesc%write%n_words    = iodesc%write%n_elemtype*lenblocks

       call genindexedblock(lenblocks,basetype,iodesc%write%elemtype,iodesc%write%filetype,displacew)

       if(debug) print *,'initdecomp: at the end of subroutine'
       !       if(iodesc%read%n_elemtype == 0 .and. iodesc%write%n_elemtype == 0) iosystem%ioproc = .false.
    endif

    deallocate(displacer,displacew)


  end subroutine initdecomp_2dof_bin


!> 
!! @public 
!! @ingroup PIO_initdecomp
!! @brief A deprecated interface to the PIO_initdecomp method.
!! @details
!! @deprecated
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param basepiotype : the type of variable(s) associated with this iodesc.  
!! @copydoc PIO_kinds
!! @param dims : an array of the global length of each dimesion of the variable(s)
!! @param lenblocks : 
!! @param compdof : mapping of the storage order of the variable to its memory order
!! @param iodofr : 
!! @param iodesc @copydoc iodesc_generate
!<
  subroutine initdecomp_1dof_bin(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodesc)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenblocks
    integer (i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: iodofr(:)     ! global degrees of freedom for io decomposition 
    type (io_desc_t), intent(out)     :: iodesc

    integer(kind=PIO_offset) :: start(1), count(1)
    ! these are not used in the binary interface

    start(1)=-1
    count(1)=-1
    call initdecomp_1dof_nf(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,start, count, iodesc)
  end subroutine initdecomp_1dof_bin

!> 
!! @public 
!! @ingroup PIO_initdecomp
!! @brief A deprecated interface to the PIO_initdecomp method.
!! @details
!! @deprecated
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param basepiotype : the type of variable(s) associated with this iodesc.
!! @copydoc PIO_kinds
!! @param dims : an array of the global length of each dimesion of the variable(s)
!! @param lenblocks : 
!! @param compdof : mapping of the storage order of the variable to its memory order
!! @param iodofr : 
!! @param iodofw :
!! @param start : used with count to give a block description of the shape of the data
!! @param count : 
!! @param iodesc @copydoc iodesc_generate
!<
  subroutine initdecomp_2dof_nf(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodofw,start, count, iodesc)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenblocks
    integer (i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: iodofr(:)     ! global degrees of freedom for io decomposition 
    integer (i4), intent(in)          :: iodofw(:)     ! global degrees of freedom for io decomposition 

    type (io_desc_t), intent(out)     :: iodesc

    integer(kind=PIO_offset), intent(in) :: start(:), count(:)
    type (io_desc_t) :: tmp


    call initdecomp_1dof_nf(iosystem, basepiotype, dims, lenblocks, compdof, iodofr, start, count, iodesc)

    call initdecomp_1dof_nf(iosystem, basepiotype, dims, lenblocks, compdof, iodofw, start, count, tmp)

    call dupiodesc2(iodesc%write,tmp%write)

    if(debug) then
       print *, _FILE_,__LINE__,iodesc%read%filetype,iodesc%read%elemtype,&
            iodesc%read%n_elemtype,iodesc%read%n_words   
       print *, _FILE_,__LINE__,iodesc%write%filetype,iodesc%write%elemtype,&
            iodesc%write%n_elemtype,iodesc%write%n_words
    end if

  end subroutine initdecomp_2dof_nf

!> 
!! @public 
!! @ingroup PIO_initdecomp
!! @brief A deprecated interface to the PIO_initdecomp method.
!! @details
!! @deprecated
!! @param iosystem : a defined PIO system descriptor, see pio_types
!! @param basepiotype : The type of variable(s) associated with this iodesc.  
!! @copydoc PIO_kinds
!! @param dims : an array of the global length of each dimesion of the variable(s)
!! @param lenblocks : 
!! @param compdof : mapping of the storage order of the variable to its memory order
!! @param iodof : 
!! @param start :
!! @param count :
!! @param iodesc @copydoc iodesc_generate
!<
  subroutine initdecomp_1dof_nf(iosystem,basepiotype,dims,lenblocks,compdof,iodof,start, count, iodesc)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in) :: lenblocks
    integer (i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: iodof(:)     ! global degrees of freedom for io decomposition 
    type (io_desc_t), intent(out)     :: iodesc
    integer :: piotype
    integer(kind=PIO_offset), intent(in) :: start(:), count(:)

    integer(i4) :: length,n_iotasks
    integer(i4) :: ndims

    integer (i4), pointer :: displace(:)  ! the displacements for the mpi data structure (read)

    integer(i4) :: prev
    integer(i4) :: glength    ! global length in words
    integer(i4) :: ii,i,dis,ierr
    integer(i4),pointer, dimension(:) :: blocklen,disp
    logical(log_kind) ::  userearranger
    logical, parameter :: check = .true.
    integer(i4) :: ndisp
#ifdef MEMCHK
    integer :: msize, rss, mshare, mtext, mstack
#endif

    piotype=PIO_type_to_mpi_type(basepiotype)

    !-------------------------------------------
    ! for testing purposes set the iomap
    ! (decompmap_t) to something basic for
    ! testing.
    !-------------------------------------------
#ifdef TIMING
    call t_startf("PIO_initdecomp")
#endif
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
    userearranger = iosystem%userearranger
    !---------------------
    ! number of dimensions
    !---------------------
    ndims = size(dims)
    !---------------------
    ! total global size
    !---------------------
    glength=1
    do i=1,ndims
       glength = glength*dims(i)
    enddo
    if(lenblocks>0) then
       ndisp=size(iodof)/lenblocks
    else
       ndisp=0
    end if
    call alloc_check(displace,ndisp)

#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif

    call alloc_check(iodesc%start,ndims)
    call alloc_check(iodesc%count,ndims)
    iodesc%start(1:size(start)) = start(:)
    iodesc%count(1:size(count)) = count(:)
    !--------------------------------------------
    ! calculate mpi data structure displacements 
    !--------------------------------------------
    if(lenblocks>0) then
       if(debug) print *,'PIO_initdecomp: calcdisplace',ndisp,size(iodof),lenblocks
       call calcdisplace(lenblocks,iodof,displace)
    end if
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif

    n_iotasks = iosystem%num_iotasks
    length = size(iodof)
    !
    !   this facilitates the use of seperate read and write descripters. 
    !
    iodesc%iomap%start  = iosystem%io_rank*length
    iodesc%iomap%length = length
    iodesc%glen = glength

    if(debug) print *,'iam: ',iosystem%io_rank,'initdecomp: userearranger: ',userearranger, glength
    if(userearranger) then 
       call rearrange_create(iosystem,compdof,iodof,iodesc)
    endif
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif


    !---------------------------------------------
    !  the setup for the mpi-io type information 
    !---------------------------------------------
    if(iosystem%ioproc) then 
       !-----------------------------------------------
       ! setup the data structure for the io operation 
       !-----------------------------------------------
       iodesc%write%n_elemtype = ndisp
       iodesc%write%n_words    = iodesc%write%n_elemtype*lenblocks

       call genindexedblock(lenblocks,piotype,iodesc%write%elemtype,iodesc%write%filetype,displace)

       if(debug) print *,'initdecomp: at the end of subroutine',iodesc%write%n_elemtype,iodesc%write%n_words
    endif
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
    call dupiodesc2(iodesc%write,iodesc%read)
    if(debug) then
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
    call t_stopf("PIO_initdecomp")
#endif
  end subroutine initdecomp_1dof_nf

!>
!! @public
!! @ingroup PIO_initdecomp
!! @brief Implements the @ref decomp_dof for PIO_initdecomp (previous name: \b initdecomp_1dof_nf_box)
!! @details  This provides the ability to describe a computational
!! decomposition in PIO using degrees of freedom method. This is
!! a decomposition that can not be easily described using a start
!! and count metehod (see @ref decomp_dof).
!! Optional parameters for this subroutine allows for the specififcation of
!! io decomposition using iostart and iocount arrays.  If iostart
!! and iocount arrays are not specified by the user, and rearrangement
!! is turned on then PIO will calculate an suitable IO decomposition.
!! Note that this subroutine was previously called \em initdecomp_1dof_nf_box
!! @param iosystem @copydoc iosystem_desc_t
!! @param basepiotype @copydoc use_PIO_kinds
!! @param dims An array of the global length of each dimesion of the variable(s)
!! @param compdof Mapping of the storage order for the computatinal decomposition to its memory order
!! @param iodesc @copydoc iodesc_generate
!! @param iostart   The start index for the block-cyclic io decomposition
!! @param iocount   The count for the block-cyclic io decomposition
!<
  subroutine PIO_initdecomp_dof(iosystem,basepiotype,dims,compdof, iodesc, iostart, iocount)
    type (iosystem_desc_t), intent(inout) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (kind=PIO_offset), optional :: iostart(:), iocount(:)
    type (io_desc_t), intent(out)     :: iodesc

    integer(i4) :: length,n_iotasks
    integer(i4) :: ndims
    integer (i4)                       :: lenblocks
    integer(i4)                       ::  piotype

    integer (i4), pointer :: displace(:)  ! the displacements for the mpi data structure (read)

    integer(i4) :: prev
    integer(i4) :: glength    ! global length in words
    integer(i4) :: ii,i,dis,ierr
    integer(i4),pointer, dimension(:) :: blocklen,disp
    logical(log_kind) ::  userearranger
    logical, parameter :: check = .true.
    integer(i4) :: ndisp
    integer(i4) :: iosize               ! rml
    integer(i4) :: msg
    logical :: is_async=.false.
#ifdef MEMCHK
    integer :: msize, rss, mshare, mtext, mstack
#endif

    integer ierror


#ifdef TIMING
    call t_startf("PIO_initdecomp_dof")
#endif
    if(iosystem%async_interface .and. .not. iosystem%ioproc) then
       msg = PIO_MSG_INITDECOMP_DOF
       is_async=.true.
       if(DebugAsync) print*,__FILE__,__LINE__, iosystem%ioranks
       if(iosystem%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, iosystem%ioroot, 1, iosystem%union_comm, ierr)
       end if
       if(DebugAsync) print*,__FILE__,__LINE__, ierr, iosystem%ioroot, iosystem%comp_rank

       call mpi_bcast(basepiotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
       if(DebugAsync) print*,__FILE__,__LINE__

       call mpi_bcast(size(dims), 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
       call mpi_bcast(dims, size(dims), mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

       if(DebugAsync) print*,__FILE__,__LINE__
       call mpi_bcast(iodesc%async_id, 1, mpi_integer, iosystem%iomaster, iosystem%intercomm, ierr)  
       if(DebugAsync) print*,__FILE__,__LINE__, iodesc%async_id
    endif

    if(minval(dims)<=0) then
       print *,_FILE_,__LINE__,dims
       call piodie(_FILE_,__LINE__,'bad value in dims argument')
    end if

    if (iosystem%comp_rank == 0 .and. debug) &
         print *,iosystem%comp_rank,': invoking PIO_initdecomp_dof'

    if(DebugAsync) print*,__FILE__,__LINE__
    piotype=PIO_type_to_mpi_type(basepiotype)
       
    !-------------------------------------------
    ! for testing purposes set the iomap
    ! (decompmap_t) to something basic for
    ! testing.
    !-------------------------------------------
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
    userearranger = iosystem%userearranger
    !---------------------
    ! number of dimensions
    !---------------------
    ndims = size(dims)
    !---------------------
    ! total global size
    !---------------------
    
    glength=1
    do i=1,ndims
       glength = glength*dims(i)
    enddo
       

    ! remember iocount() is only defined on io procs
    call alloc_check(iodesc%start,ndims)
    call alloc_check(iodesc%count,ndims)
    iodesc%basetype=piotype
       
    iodesc%compsize=size(compdof)
       
    iodesc%start=0
    iodesc%count=0
    if(DebugAsync) print*,__FILE__,__LINE__, iosystem%num_tasks, iosystem%num_iotasks, iosystem%io_rank, iosystem%io_comm, iosystem%ioranks

    if (iosystem%ioproc) then
       if(present(iostart) .and. present(iocount)) then
          iodesc%start = iostart
          iodesc%count = iocount
       else if(present(iostart) .or. present(iocount)) then
          call piodie( _FILE_,__LINE__, &
               'both optional parameters start and count must be provided')
       else
          call getiostartandcount(iosystem%num_tasks, ndims, dims, iosystem%num_iotasks, iosystem%io_rank, iosystem%io_comm, iodesc%start, iodesc%count)
       end if
       iosize=1
       do i=1,ndims
          iosize=iosize*iodesc%count(i)
       end do
       call mpi_allreduce(iosize, iodesc%maxiobuflen, 1, mpi_integer, mpi_max, iosystem%io_comm, ierr)
       call checkmpireturn('mpi_allreduce in initdecomp',ierr)

       if(debug) print *,'IAM: ',iosystem%comp_rank,' after getiostartandcount: count is: ',iodesc%count


       lenblocks=iodesc%count(1)


       if(lenblocks>0) then
          ndisp=iosize/lenblocks
       else
          ndisp=0
       end if
       call alloc_check(displace,ndisp)
       
       !--------------------------------------------
       ! calculate mpi data structure displacements 
       !--------------------------------------------
      
       if(debug) print *,'PIO_initdecomp: calcdisplace', &
            ndisp,iosize,lenblocks, iodesc%start, iodesc%count
       call calcdisplace_box(dims,iodesc%start,iodesc%count,ndims,displace)
          
       n_iotasks = iosystem%num_iotasks
       length = iosize                      ! rml

       !
       !   this facilitates the use of seperate read and write descripters. 
       !

       iodesc%iomap%start  = iosystem%io_rank*length
       iodesc%iomap%length = length
       iodesc%glen = glength
    endif
    if(DebugAsync) print*,__FILE__,__LINE__

#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
    if(debug) print *,__FILE__,__LINE__,'iam: ',iosystem%io_rank, &
         'initdecomp: userearranger: ',userearranger, glength

    if(userearranger) then 
       if(DebugAsync) print*,__FILE__,__LINE__
       call rearrange_create( iosystem,compdof,dims,ndims,iodesc)
    endif
    if(DebugAsync) print*,__FILE__,__LINE__
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif
    
    !---------------------------------------------
    !  the setup for the mpi-io type information 
    !---------------------------------------------
    if(iosystem%ioproc) then 
       !-----------------------------------------------
       ! setup the data structure for the io operation 
       !-----------------------------------------------
       iodesc%write%n_elemtype = ndisp
       iodesc%write%n_words    = iodesc%write%n_elemtype*lenblocks
       call genindexedblock(lenblocks,piotype,iodesc%write%elemtype,iodesc%write%filetype,displace)
       
       if(debug) print *,__FILE__,__LINE__,iodesc%write%n_elemtype,iodesc%write%n_words,iodesc%write%elemtype,iodesc%write%filetype
    else
       iodesc%write%n_elemtype=0
       iodesc%write%n_words=0
       iodesc%write%elemtype = mpi_datatype_null
       iodesc%write%filetype = mpi_datatype_null
    endif
    
#ifdef MEMCHK	
    call get_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,_FILE_,__LINE__,'mem=',rss
    end if
#endif

    call dupiodesc2(iodesc%write,iodesc%read)
    
!    if(debug) then
!       print*, _FILE_,__LINE__,iodesc%read%filetype,iodesc%read%elemtype,&
!            iodesc%read%n_elemtype,iodesc%read%n_words   
!       print *, _FILE_,__LINE__,iodesc%write%filetype,iodesc%write%elemtype,&
!            iodesc%write%n_elemtype,iodesc%write%n_words
!    end if
    
    if (iosystem%ioproc) then
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
    call t_stopf("PIO_initdecomp_dof")
#endif

  end subroutine PIO_initdecomp_dof

  subroutine getiostartandcount(ntasks, ndims, gdims, num_io_procs, iorank, iocomm, start, count)

    implicit none

    integer, intent(in) :: ntasks
    integer, intent(in) :: ndims, gdims(ndims), num_io_procs

    integer, intent(in) :: iocomm
    integer(kind=PIO_offset), intent(out) :: start(ndims), count(ndims)   ! start and count arrays 

    integer :: ierr
    integer :: n,m
    integer :: use_io_procs,iorank, sdims, cnt
    logical :: done
    integer :: xiam, xpes, ps, pe, ds, de, ns, pe1, ps1
    integer :: size,tsize
    integer :: it
    real(r8) :: fanfactor, fanlimit
    real(r8) :: rtmp
    integer :: testvalue
    integer, parameter :: minblocksize=16        ! minimum block size on a task
    integer, parameter :: maxit=1               ! maximum number of times to iterate on the fanin/out limiter  (Probably want a better solution)
    integer,allocatable :: pes_per_dim(:), step(:)
    integer,allocatable :: bsize(:),nblocks(:),fblocks(:)


    tsize=1
    do n=1,ndims
       tsize=tsize*gdims(n)
    end do
    allocate(pes_per_dim(ndims))
    allocate(step(ndims))
    allocate(bsize(ndims))
    allocate(nblocks(ndims))
    allocate(fblocks(ndims))

    use_io_procs=num_io_procs
    do while(tsize/minblocksize < use_io_procs .and. use_io_procs>1)
       use_io_procs=use_io_procs-1
    end do
    if(Debug) print *,'iorank: ',iorank,' getiostartandcount: use_io_procs: ',use_io_procs

    start = 1
    count = 0


!!    if(iorank>=use_io_procs) return 

    !-----------------

    cnt = 1
    sdims = ndims
    do while (cnt < use_io_procs .and. sdims > 0)
       !       cnt = cnt * ceiling(dble(gdims(sdims))/dble(minblocksize))
       cnt = cnt * gdims(sdims)
       sdims = sdims - 1
    enddo
    if (sdims < 0) then
       call piodie( _FILE_,__LINE__, &
            'error in sdims',sdims)
    endif

    do m = 1,sdims
       start(m) = 1
       count(m) = gdims(m)
       bsize(m) = gdims(m)
       nblocks(m) = 1
       fblocks(m) = 1
    enddo
    if(Debug) print *,'iorank: ',iorank, ' getiostartandcount: sdims: ',sdims
    if(Debug) print *,'iorank: ',iorank, ' getiostartandcount: count: ',count

    fanlimit  = 50.00
    fanfactor = fanlimit + 1.0  !we want at least one trip through the do while loop  

    it = 0
    step(:) = 1
    do while (fanfactor > fanlimit .and. it < maxit ) 
       xpes = use_io_procs
       xiam = iorank   ! goes from 0 to xpes-1
       do m = ndims, sdims+1, -1
!          if(xpes >= gdims(m)) then
          if(xpes > gdims(m)) then
             ps = -1
             ns = 1
             do while (ps < 0 .and. ns <= gdims(m) )
                ps1 = int((dble(xpes)*dble((ns-1)*step(m)))/dble(gdims(m)))
                pe1 = int((dble(xpes)*dble(ns*step(m)  ))/dble(gdims(m))) - 1
                if (xiam >= ps1 .and. xiam <= pe1) then
                   ps = ps1
                   pe = pe1
                   start(m) = (ns-1)*step(m) + 1
                   count(m) = step(m)
                end if
                ns = ns+1
             end do
             xpes = pe - ps + 1
             xiam = xiam - ps
             !          write(6,*) 'tcx1 ',iorank,m,start(m),count(m)
             step(m)=nextlarger(step(m),gdims(m))
             if(step(m) == gdims(m)) fanlimit = fanlimit + 10.0
          else
             if (m /= sdims+1) then
                call piodie( _FILE_,__LINE__, &
                     'm /= sdims+1',ival1=m,ival2=sdims)
             endif
             ds = int((dble(gdims(m))*dble(xiam  ))/dble(xpes)) + 1
             de = int((dble(gdims(m))*dble(xiam+1))/dble(xpes))
             start(m) = ds
             count(m) = de-ds+1
          end if

          if (start(m) < 1 .or. count(m) < 1) then
             print *, 'start =',start, ' count=',count
             call piodie( _FILE_,__LINE__, &
                  'start or count failed to converge')
          endif

       enddo
       do m=1,ndims
          pes_per_dim(m) = gdims(m)/count(m)
       enddo
       ! -----------------------------------------------
       ! note this caculation assumes that the the first 
       ! two horizontal dimensions are decomposed,
       ! -----------------------------------------------
       if(ndims==1) then
          fanfactor = ntasks/pes_per_dim(1)
       else
          fanfactor = ntasks/(pes_per_dim(1)*pes_per_dim(2))
       end if
       call mpi_allreduce(fanfactor,rtmp,1,MPI_REAL8,MPI_MAX,iocomm,ierr)
       fanfactor=rtmp
       if(Debug) print *,'iorank: ',iorank,'getiostartandcount: pes_per_dim is: ',pes_per_dim
       if(Debug) print *,'iorank: ',iorank,' getiostartandcount: fan factor is: ',fanfactor
       it=it+1
    enddo
    deallocate(step)
    deallocate(pes_per_dim)
    deallocate(bsize,nblocks,fblocks)
    !   stop 'end of getiostartandcount'

    if(iorank>=use_io_procs) then 
	start = 1
        count = 0 
    endif

  end subroutine getiostartandcount

  integer function nextlarger(current,value) result(res)

  integer :: current,value,rem
  !-----------------------------------------
  ! This function finds a value res that is
  !  larger than current that divides value evenly
  !-----------------------------------------
  res=current
  rem = 1
  do while ((rem .ne. 0) .and. (res .lt. value))
     res = res + 1
     rem = MOD(value,res)
!     print *,'res,rem: ',res,rem
  enddo

  end function nextlarger

  !************************************
  ! dupiodesc2
  !

  subroutine dupiodesc2(src, dest)
    use pio_types, only : io_desc2_t
    type(io_desc2_t), intent(in) :: src
    type(io_desc2_t), intent(out) :: dest

    dest%filetype = src%filetype
    dest%elemtype = src%elemtype
    dest%n_elemtype = src%n_elemtype
    dest%n_words = src%n_words
  end subroutine dupiodesc2



  !************************************
  ! genindexedblock
  !
  ! given input lenblocks, basetype, and displacement
  ! create two mpi types: 
  !   elemtype - a single block of basetype repeated lenblocks times
  !   filetype - elemtype repeated at each entry in displacement()
  !              (i.e. size(displacement) entries)
  !


  subroutine genindexedblock(lenblocks,basetype,elemtype,filetype,displace)
    use pio_types, only : pio_double, pio_int, pio_real, pio_char
#ifdef NO_MPI2
    use pio_support, only : mpi_type_create_indexed_block
#endif
    integer(i4), intent(in) :: lenblocks     ! length of blocks
    integer(i4), intent(in) :: basetype      ! base mpi type 
    integer(i4), intent(inout) :: elemtype   ! elementary mpi type
    integer(i4), intent(inout) :: filetype   ! file mpi type 
    integer(i4), intent(in) :: displace(:)   ! mpi displacement in the array

    integer(i4) :: numblocks,i,ierr,prev
    logical, parameter :: check = .true.

    integer:: nints, nadds, ndtypes, comb, lbasetype

    numblocks = size(displace)

    !tcx - allow empty iodofs
    if (numblocks > 0) then
       prev = displace(1)
       do i=2,numblocks
          if(prev > displace(i)) then
             print *,'genindexedblock: error detected: non-monotonic increasing displace detected!'
          endif
          prev = displace(i)
       enddo

    endif
    select case(basetype)
    case (PIO_double)
       lbasetype=mpi_real8
    case (PIO_real  )
       lbasetype=mpi_real4
    case (PIO_int)
       lbasetype=mpi_integer
    case (PIO_char)
       lbasetype=mpi_character
    case default
       lbasetype=basetype
    end select


#ifdef _MPISERIAL
    ! when compiling w/mpiserial for snetcdf output, these fields are not used
    elemtype=0
    filetype=0
    ! _MPISERIAL
#else
    if(lenblocks<1) then
       elemtype = lbasetype
       filetype = lbasetype
    else
       call mpi_type_contiguous(lenblocks,lbasetype,elemtype,ierr)
       if(check) call checkmpireturn('genindexedblock: after call to type_contiguous: ',ierr)
       call mpi_type_commit(elemtype,ierr)
       if(check) call checkmpireturn('genindexedblock: after call to type_commit: ',ierr)
       call mpi_type_create_indexed_block(numblocks,1,displace,elemtype,filetype,ierr)
       if(check) call checkmpireturn('genindexedblock: after call to type_create_indexed_block: ',ierr)
       call mpi_type_commit(filetype,ierr)
       if(check) call checkmpireturn('genindexedblock: after call to type_commit: ',ierr)
       call mpi_type_get_envelope(elemtype, nints, nadds, ndtypes, comb, ierr)
    end if
    ! _MPISERIAL
#endif

  end subroutine genindexedblock
!> 
!! @public
!! @ingroup PIO_init
!! @brief initialize the pio subsystem. 
!! @details  This is a collective call.  Input parameters are read on comp_rank=0
!!   values on other tasks are ignored.  This variation of PIO_init locates the IO tasks on a subset 
!!   of the compute tasks.
!! @param comp_rank mpi rank of each participating task,
!! @param comp_comm the mpi communicator which defines the collective.
!! @param num_iotasks the number of iotasks to define.
!! @param num_aggregator the mpi aggregator count
!! @param stride the stride in the mpi rank between io tasks.
!! @param rearr @copydoc PIO_rearr_method
!! @param iosystem a derived type which can be used in subsequent pio operations (defined in PIO_types).
!! @param base @em optional argument can be used to offset the first io task - default base is task 1.
!<
  subroutine init_intracom(comp_rank, comp_comm, num_iotasks, num_aggregator, stride,  rearr, iosystem,base)
    use pio_types, only : pio_internal_error, pio_rearr_none
    integer(i4), intent(in) :: comp_rank
    integer(i4), intent(in) :: comp_comm
    integer(i4), intent(in) :: num_iotasks 
    integer(i4), intent(in) :: num_aggregator
    integer(i4), intent(in) :: stride
    integer(i4), intent(in) :: rearr
    type (iosystem_desc_t), intent(out)  :: iosystem  ! io descriptor to initalize
    integer(i4), intent(in),optional :: base

    integer(i4) :: n_iotasks
    integer(i4) :: length
    integer(i4) :: ngseg,io_rank,i,lbase, io_comm,ierr 
    integer(i4) :: lstride, itmp
    integer(i4), pointer :: iotmp(:),iotmp2(:)

    integer :: mpi_comm_io, intercomm

    character(len=5) :: cb_nodes
    logical(log_kind), parameter :: check = .true.
    logical :: async_setup = .false.

    integer(i4) :: j

    integer(i4) :: mpi_group_world, mpi_group_io, mpi_group_compute

    integer(i4) :: iotask
    integer(i4) :: rearrFlag

#ifdef TIMING
    call t_startf("PIO_init")
#endif

    iosystem%error_handling = PIO_internal_error
    iosystem%union_comm = comp_comm
    iosystem%comp_comm = comp_comm
    iosystem%comp_rank = comp_rank
    iosystem%intercomm = MPI_COMM_NULL
    iosystem%my_comm = comp_comm
    iosystem%async_interface = .false.
#ifndef _MPISERIAL
    iosystem%info = mpi_info_null
#endif
    call mpi_comm_size(comp_comm,iosystem%num_tasks,ierr)

    iosystem%num_comptasks = iosystem%num_tasks
    iosystem%union_rank = comp_rank
    iosystem%rearr = rearr

    if(check) call checkmpireturn('init: after call to comm_size: ',ierr)
    ! ---------------------------------------
    ! need some more error checking code for 
    ! setting of number of io nodes
    ! ---------------------------------------

    n_iotasks=num_iotasks

    if (n_iotasks>iosystem%num_tasks) then
       n_iotasks=iosystem%num_tasks
       if (iosystem%comp_rank==0) then
          print *,'***warning, reducing io tasks to ',n_iotasks, &
               ' because there are not enough processors'
       endif
    endif

    lbase = 0
    ! unless you are using all procs, shift off the masterproc
    if(n_iotasks<iosystem%num_tasks) then
       lbase=1
    end if
    if (present(base)) then
       if(base>=0 .and. base<iosystem%num_tasks) lbase = base
    endif

    if(debug) print *,'init: iosystem%num_tasks,n_iotasks,num_aggregator: ',iosystem%num_tasks,n_iotasks,num_aggregator

    ! --------------------------
    ! select which nodes are io
    ! nodes and set ioproc
    ! --------------------------
    lstride = stride
    ! Check sanity of input arguments

    call mpi_bcast(iosystem%rearr, 1, mpi_integer, 0, iosystem%comp_comm, ierr)
    call mpi_bcast(n_iotasks, 1, mpi_integer, 0, iosystem%comp_comm, ierr)
    call mpi_bcast(lstride, 1, mpi_integer, 0, iosystem%comp_comm, ierr)
    call mpi_bcast(lbase, 1, mpi_integer, 0, iosystem%comp_comm, ierr)

    if (lbase+(n_iotasks-1)*lstride >= iosystem%num_tasks) then
       print *,_FILE_,__LINE__,lbase,n_iotasks,lstride,iosystem%num_tasks
       call piodie(_FILE_,__LINE__,'not enough procs for the stride')
    endif

    iosystem%ioproc = .false.

#ifdef BGx

    call alloc_check(iotmp,iosystem%num_tasks,'init:num_tasks')
    call alloc_check(iotmp2,iosystem%num_tasks,'init:num_tasks')
    !---------------------------------------------------
    ! Note for Blue Gene n_iotasks get overwritten in 
    ! determineiotasks   
    !
    ! Entry: it is the number of IO-clients per IO-node
    ! Exit:  is is the total number of IO-tasks
    !---------------------------------------------------
    if (iosystem%rearr == PIO_rearr_none) then
       rearrFlag = 0
    else
       rearrFlag = 1
    endif
    call determineiotasks(comp_comm,n_iotasks,lbase,lstride,rearrFlag,iotask)


    iotmp(:)=0
    if(iotask==1) then 
       iosystem%ioproc = .true.
       iotmp(comp_rank+1) = 1
    endif
    iotmp2(:)=0 
    call MPI_allreduce(iotmp,iotmp2,iosystem%num_tasks,MPI_INTEGER,MPI_SUM,comp_comm,ierr)

    n_iotasks=SUM(iotmp2)
    iosystem%num_iotasks =n_iotasks
    call alloc_check(iosystem%ioranks,n_iotasks,'init:n_ioranks')
    j=1
    iosystem%iomaster=-1
    do i=1,iosystem%num_tasks
       if(iotmp2(i) == 1) then 
          iosystem%ioranks(j) = i-1
	  j=j+1
	  if(iosystem%iomaster<0) iosystem%iomaster=i-1
       endif
    enddo
    call dealloc_check(iotmp)
    call dealloc_check(iotmp2)
    iotask = 0
    if(iosystem%ioproc) then 
       iotask = 1
    endif
    call identity(comp_comm,iotask)

#else

    iosystem%num_iotasks = n_iotasks
    call alloc_check(iosystem%ioranks,n_iotasks,'init:n_ioranks')

    do i=1,n_iotasks
       iosystem%ioranks(i)=(lbase + (i-1)*lstride)

       if (iosystem%ioranks(i)>=iosystem%num_tasks) then
          call piodie( _FILE_,__LINE__, &
               'tried to assign io processor beyond max rank ',&
               iosystem%ioranks(i), &
               ' num_tasks=',iosystem%num_tasks )
       endif

       if(comp_rank == iosystem%ioranks(i))  iosystem%ioproc = .true.
    enddo
#endif

    iosystem%iomaster = iosystem%ioranks(1)
    iosystem%ioroot = iosystem%ioranks(1)
    if(debug) print *,'init: iam: ',comp_rank,'io processor: ',iosystem%ioproc, 'io rank ',&
         iosystem%io_rank, iosystem%iomaster		  


    if(debug) print *,'init: iam: ',comp_rank,' before allocate(status): n_iotasks: ',n_iotasks

    if (iosystem%rearr == PIO_rearr_none) then
       iosystem%userearranger= .false.
    else
       iosystem%userearranger= .true.
    endif

#if defined(USEMPIIO) || defined(_PNETCDF) || defined(_NETCDF4)
#ifndef _MPISERIAL
    call mpi_info_create(iosystem%info,ierr)
#endif

    if(debug) print *,'iam: ',iosystem%io_rank,__LINE__,'init: userearranger: ',iosystem%userearranger

    !---------------------------------
    ! initialize the rearranger system 
    !---------------------------------

    if (iosystem%userearranger) then
       call rearrange_init(iosystem)
    endif

    iosystem%io_rank=-1
    call mpi_comm_group(comp_comm,mpi_group_world,ierr)
    if(check) call checkmpireturn('init: after call to comm_group: ',ierr)

    call mpi_group_incl(mpi_group_world,n_iotasks,iosystem%ioranks,mpi_group_io,ierr)
    if(check) call checkmpireturn('init: after call to group_range_incl: ',ierr)

    if(DebugAsync) print *,__FILE__,__LINE__,'n: ',n_iotasks, ' r: ',iosystem%ioranks, ' g: ',mpi_group_io

    !-----------------------
    ! setup io_comm and io_rank
    !-----------------------

    call mpi_comm_create(comp_comm,mpi_group_io,iosystem%io_comm,ierr)
    if(check) call checkmpireturn('init: after call to comm_create: ',ierr)
    
    if(iosystem%ioproc) call mpi_comm_rank(iosystem%io_comm,iosystem%io_rank,ierr)
    if(check) call checkmpireturn('init: after call to comm_rank: ',ierr)
    ! turn on mpi-io aggregation 
    !DBG    print *,'PIO_init: before call to setnumagg'
    itmp = num_aggregator
    call mpi_bcast(itmp, 1, mpi_integer, 0, iosystem%comp_comm, ierr)

    if(itmp .gt. 0) then 
       write(cb_nodes,('(i5)')) itmp
#ifdef BGx
       call PIO_set_hint(iosystem,"bgl_nodes_pset",trim(adjustl(cb_nodes)))
#else
       call PIO_set_hint(iosystem,"cb_nodes",trim(adjustl(cb_nodes)))
#endif       
    endif

#ifdef PIO_GPFS_HINTS
    call PIO_set_hint(iosystem,"ibm_largeblock_io","true")
#endif
#ifdef PIO_LUSTRE_HINTS
    call PIO_set_hint(iosystem, 'romio_ds_read','disable') 
    call PIO_set_hint(iosystem,'romio_ds_write','disable') 
#endif
#endif

#ifdef TIMING
    call t_stopf("PIO_init")
#endif
  end subroutine init_intracom


!> 
!! @public
!! @ingroup PIO_init
!! @brief Initialize the pio subsystem.
!! @details  This is a collective call.  Input parameters are read on comp_rank=0
!!   values on other tasks are ignored.  This variation of PIO_init sets up a distinct set of tasks
!!   to handle IO, these tasks do not return from this call.  Instead they go to an internal loop 
!!   and wait to receive further instructions from the computational tasks 
!! @param component_count The number of computational components to associate with this IO component
!! @param peer_comm  The communicator from which all other communicator arguments are derived
!! @param comp_comms The computational communicator for each of the computational components
!! @param io_comm    The io communicator 
!! @param iosystem a derived type which can be used in subsequent pio operations (defined in PIO_types).
!<
  subroutine init_intercom(component_count, peer_comm, comp_comms, io_comm, iosystem)
    use pio_types, only : pio_internal_error, pio_rearr_box
    integer, intent(in) :: component_count
    integer, intent(in) :: peer_comm
    integer, intent(in) :: comp_comms(component_count)   !  The compute communicator
    integer, intent(in) :: io_comm     !  The io communicator

    type (iosystem_desc_t), intent(out)  :: iosystem(component_count)  ! io descriptor to initalize

    integer :: ierr
    logical :: is_inter
    logical, parameter :: check=.true.

    integer :: i, j, iam, io_leader, comp_leader
    integer(i4), pointer :: iotmp(:)

#ifdef TIMING
    call t_startf("PIO_init")
#endif
#if defined(NO_MPI2) || defined(_MPISERIAL)
    call piodie( _FILE_,__LINE__,'The PIO async interface requires an MPI2 complient MPI library')
#else 
    do i=1,component_count
       iosystem(i)%error_handling = PIO_internal_error
       iosystem(i)%comp_comm = comp_comms(i)
       iosystem(i)%io_comm = io_comm
       iosystem(i)%info = mpi_info_null
       iosystem(i)%comp_rank= -1
       iosystem(i)%io_rank  = -1
       iosystem(i)%async_interface = .true.
       iosystem(i)%comproot = MPI_PROC_NULL
       iosystem(i)%ioroot = MPI_PROC_NULL
       iosystem(i)%compmaster= MPI_PROC_NULL
       iosystem(i)%iomaster = MPI_PROC_NULL 


       if(io_comm/=MPI_COMM_NULL) then
          ! Find the rank of the io leader in peer_comm
          call mpi_comm_rank(io_comm,iosystem(i)%io_rank, ierr)
          if(iosystem(i)%io_rank==0) then 
             call mpi_comm_rank(peer_comm, iam, ierr)
          else
             iam = -1
          end if
          call mpi_allreduce(iam, io_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
          ! Find the rank of the comp leader in peer_comm
          iam = -1
          call mpi_allreduce(iam, comp_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
          ! create the intercomm
          call mpi_intercomm_create(io_comm, 0, peer_comm, comp_leader, i, iosystem(i)%intercomm, ierr)
          ! create the union_comm
          call mpi_intercomm_merge(iosystem(i)%intercomm, .true., iosystem(i)%union_comm, ierr)
       else
          ! Find the rank of the io leader in peer_comm
          iam = -1
          call mpi_allreduce(iam, io_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
          ! Find the rank of the comp leader in peer_comm
          iosystem(i)%comp_rank = -1
          if(comp_comms(i)/=MPI_COMM_NULL) then
             call mpi_comm_rank(comp_comms(i),iosystem(i)%comp_rank, ierr)          
             if(iosystem(i)%comp_rank==0) then
                call mpi_comm_rank(peer_comm, iam, ierr)
             else
                iam=-1
             end if
          end if
          call mpi_allreduce(iam, comp_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
          
          ! create the intercomm
          call mpi_intercomm_create(comp_comms(i), 0, peer_comm, io_leader, i, iosystem(i)%intercomm, ierr)
          ! create the union comm
          call mpi_intercomm_merge(iosystem(i)%intercomm, .false., iosystem(i)%union_comm, ierr)
       end if
       if(Debugasync) print *,__FILE__,__LINE__,i, iosystem(i)%intercomm, iosystem(i)%union_comm

       if(iosystem(i)%union_comm /= MPI_COMM_NULL) then
          call mpi_comm_rank(iosystem(i)%union_comm, iosystem(i)%union_rank, ierr)
          if(check) call checkmpireturn('init: after call to comm_rank: ',ierr)
          call mpi_comm_size(iosystem(i)%union_comm, iosystem(i)%num_tasks, ierr)
          if(check) call checkmpireturn('init: after call to comm_size: ',ierr)

             
          if(io_comm /= MPI_COMM_NULL) then
             call mpi_comm_size(io_comm, iosystem(i)%num_iotasks, ierr)
             if(check) call checkmpireturn('init: after call to comm_size: ',ierr)

             if(iosystem(i)%io_rank==0) then
                iosystem(i)%iomaster = MPI_ROOT
                iosystem(i)%ioroot = iosystem(i)%union_rank
             end if
             iosystem(i)%ioproc = .true.
             iosystem(i)%compmaster = 0

             call pio_msg_handler_init(io_comm, iosystem(i)%io_rank)
          end if


          if(comp_comms(i) /= MPI_COMM_NULL) then
             call mpi_comm_size(comp_comms(i), iosystem(i)%num_comptasks, ierr)
             if(check) call checkmpireturn('init: after call to comm_size: ',ierr)

             iosystem(i)%iomaster = 0
             iosystem(i)%ioproc = .false.
             if(iosystem(i)%comp_rank==0) then
                iosystem(i)%compmaster = MPI_ROOT
                iosystem(i)%comproot = iosystem(i)%union_rank
             end if

          end if

          iosystem(i)%userearranger = .true.
          iosystem(i)%rearr = PIO_rearr_box
          
          if(Debugasync) print *,__FILE__,__LINE__
          
          call MPI_allreduce(iosystem(i)%comproot, j, 1, MPI_INTEGER, MPI_MAX,iosystem(i)%union_comm,ierr)
          iosystem%comproot=j
          call MPI_allreduce(iosystem(i)%ioroot, j, 1, MPI_INTEGER, MPI_MAX,iosystem(i)%union_comm,ierr)
          iosystem%ioroot=j

          if(Debugasync) print *,__FILE__,__LINE__, i, iosystem(i)%comproot, iosystem(i)%ioroot

          if(io_comm/=MPI_COMM_NULL) then
             call mpi_bcast(iosystem(i)%num_comptasks, 1, mpi_integer, iosystem(i)%compmaster,iosystem(i)%intercomm, ierr)

             call mpi_bcast(iosystem(i)%num_iotasks, 1, mpi_integer, iosystem(i)%iomaster, iosystem(i)%intercomm, ierr)

             call alloc_check(iotmp,iosystem(i)%num_iotasks,'init:iotmp')
             iotmp(:) = 0
             iotmp( iosystem(i)%io_rank+1)=iosystem(i)%union_rank

          end if
          if(comp_comms(i)/=MPI_COMM_NULL) then
             call mpi_bcast(iosystem(i)%num_comptasks, 1, mpi_integer, iosystem(i)%compmaster, iosystem(i)%intercomm, ierr)

             call mpi_bcast(iosystem(i)%num_iotasks, 1, mpi_integer, iosystem(i)%iomaster, iosystem(i)%intercomm, ierr)

             call alloc_check(iotmp,iosystem(i)%num_iotasks,'init:iotmp')
             iotmp(:)=0

          end if

          iosystem(i)%my_comm = iosystem(i)%intercomm

          call alloc_check(iosystem(i)%ioranks, iosystem(i)%num_iotasks,'init:n_ioranks')
          if(Debugasync) print *,__FILE__,__LINE__,iotmp
          call MPI_allreduce(iotmp,iosystem(i)%ioranks,iosystem(i)%num_iotasks,MPI_INTEGER,MPI_MAX,iosystem(i)%union_comm,ierr)
          
          if(Debugasync) print *,__FILE__,__LINE__,iosystem(i)%ioranks
          call dealloc_check(iotmp)
          
          !---------------------------------
          ! initialize the rearranger system 
          !---------------------------------
          if (iosystem(i)%userearranger) then
             call rearrange_init(iosystem(i))
          endif
       end if
    
#if defined(USEMPIIO) || defined(_PNETCDF) || defined(_NETCDF4)
#ifndef _MPISERIAL
       ! turn on mpi-io aggregation 
       !DBG    print *,'PIO_init: before call to setnumagg'
       itmp = num_aggregator
       call mpi_bcast(itmp, 1, mpi_integer, 0, iosystem%union_comm, ierr)
       call mpi_info_create(iosystem(i)%info,ierr)
       if(itmp .gt. 0) then 
          write(cb_nodes,('(i5)')) itmp
#ifdef BGx
          call PIO_set_hint(iosystem(i),"bgl_nodes_pset",trim(adjustl(cb_nodes)))
#else
          call PIO_set_hint(iosystem(i),"cb_nodes",trim(adjustl(cb_nodes)))
#endif       
       endif

#ifdef PIO_GPFS_HINTS
       call PIO_set_hint(iosystem(i),"ibm_largeblock_io","true")
#endif
#ifdef PIO_LUSTRE_HINTS
       call PIO_set_hint(iosystem(i), 'romio_ds_read','disable') 
       call PIO_set_hint(iosystem(i),'romio_ds_write','disable') 
#endif
#endif
#endif
    end do

    if(DebugAsync) print*,__FILE__,__LINE__, iosystem(1)%ioranks



    ! This routine does not return
    if(io_comm /= MPI_COMM_NULL) call pio_msg_handler(component_count,iosystem) 
    
    if(DebugAsync) print*,__FILE__,__LINE__, iosystem(1)%ioranks
#ifdef TIMING
    call t_stopf("PIO_init")
#endif
#endif
  end subroutine init_intercom

!>
!! @public
!! @defgroup PIO_recommend_iotasks
!! @brief Recommend a subset of tasks in comm to use as IO tasks
!! @details  This subroutine will give PIO's best recommendation for the number and
!!    location of iotasks for a given system there is no requirement to follow this recommendation.
!!    Using the recommendation requires that PIO_BOX_RERRANGE be used
!! @param A communicator of mpi tasks to choose from
!! @param miniotasks \em optional The minimum number of IO tasks the caller desires
!! @param maxiotasks \em optional The maximum number of IO tasks the caller desires
!! @param iotask if true pio recommends that this task be used as an iotask

  subroutine pio_recommend_iotasks(comm, ioproc, numiotasks, miniotasks, maxiotasks )
    integer, intent(in) :: comm
    logical, intent(out) :: ioproc
    integer, intent(out) :: numiotasks
    integer, optional, intent(in) :: miniotasks, maxiotasks

    integer :: num_tasks, ierr, iotask, iotasks, iam

    integer(i4), pointer :: iotmp(:),iotmp2(:)

    call mpi_comm_size(comm,num_tasks,ierr)
    call mpi_comm_rank(comm,iam,ierr)

#ifdef BGx    
    call alloc_check(iotmp,num_tasks,'init:num_tasks')
    call alloc_check(iotmp2,num_tasks,'init:num_tasks')
    !---------------------------------------------------
    ! Note for Blue Gene n_iotasks get overwritten in 
    ! determineiotasks   
    !
    ! Entry: it is the number of IO-clients per IO-node
    ! Exit:  is is the total number of IO-tasks
    !---------------------------------------------------

    numiotasks=-(miniotasks+maxiotasks)/2
    call determineiotasks(comm,numiotasks,1,0,1,iotask)

    iotmp(:)=0
    if(iotask==1) then 
       ioproc = .true.
       iotmp(iam+1) = 1
    endif
    iotmp2(:)=0 
    call MPI_allreduce(iotmp,iotmp2,num_tasks,MPI_INTEGER,MPI_SUM,comm,ierr)

    numiotasks=SUM(iotmp2)

    call dealloc_check(iotmp)
    call dealloc_check(iotmp2)

    call identity(comm,iotask)
#endif


  end subroutine pio_recommend_iotasks


!> 
!! @public
!! @defgroup PIO_set_hint  PIO_set_hint
!! @brief set file system hints using mpi_info_set
!! @details This is a collective call which expects the following parameters:
!! @param iosystem @copydoc io_desc_t
!! @param hint  the string name of the hint to define
!! @param hintval  the string value to set the hint to
!! @retval ierr @copydoc  error_return
!<
  subroutine PIO_set_hint(iosystem, hint, hintval)
    type (iosystem_desc_t), intent(out)  :: iosystem  ! io descriptor to initalize
    character(len=*), intent(in) :: hint, hintval
    
    integer :: ierr
#if defined(USEMPIIO) || defined(_PNETCDF) || defined(_NETCDF4)
    if(iosystem%ioproc) then
#ifndef _MPISERIAL
       call mpi_info_set(iosystem%info,hint,hintval,ierr)
       if(iosystem%io_rank==0 || Debug) print *,'Setting mpi info: ',hint,'=',hintval
#endif
       call checkmpireturn('PIO_set_hint',ierr)
    end if
#endif
  end subroutine PIO_set_hint


!> 
!! @public
!! @ingroup PIO_finalize 
!! @brief finalizes the pio subsystem.
!! @details This is a collective call which expects the following parameters
!! @param iosystem : @copydoc io_desc_t
!! @retval ierr @copydoc  error_return
!<
  subroutine finalize(iosystem,ierr)
     type (iosystem_desc_t), intent(inout) :: iosystem 
     integer(i4), intent(out) :: ierr
     
     integer :: msg

     if(iosystem%async_interface .and. iosystem%comp_rank==0) then
       msg = PIO_MSG_EXIT
       call mpi_send(msg, 1, mpi_integer, iosystem%ioroot, 1, iosystem%union_comm, ierr)
     end if

#ifndef _MPISERIAL
     if(iosystem%info .ne. mpi_info_null) then 
        call mpi_info_free(iosystem%info,ierr) 
     endif
     if(iosystem%io_comm .ne. mpi_comm_null) then 
        call mpi_comm_free(iosystem%io_comm,ierr)
     endif
#endif
     ierr = 0

  end subroutine finalize


!>
!! @public
!! @ingroup PIO_getnumiotasks
!! @brief This returns the number of IO-tasks that PIO is using 
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param numiotasks : the number of IO-tasks
!<
   subroutine getnumiotasks(iosystem,numiotasks)
       type (iosystem_desc_t), intent(in) :: iosystem
       integer(i4), intent(out) :: numiotasks
       numiotasks = iosystem%num_iotasks
   end subroutine getnumiotasks




  !=============================================
  !  dupiodesc:
  !
  !   duplicate the io descriptor
  !
  !=============================================


  ! rml: possible problem here wrt dubbing the box rearranger
  ! data, as well as maybe the mct rearranger???

!> 
!! @public 
!! @ingroup PIO_dupiodesc
!! @brief duplicates an existing io descriptor
!! @details
!! @param src :  an io description handle returned from @ref PIO_initdecomp (see PIO_types)
!! @param dest : the newly created io descriptor with the same characteristcs as src.
!<
  subroutine dupiodesc(src,dest)

    integer :: n
    type (io_desc_t), intent(in) :: src
    type (io_desc_t), intent(inout) :: dest


    dest%glen        =  src%glen
    if(associated(src%start)) then
       n = size(src%start)
       allocate(dest%start(n))
       dest%start(:)       =  src%start(:)
    endif

    if(associated(src%count)) then
       n = size(src%count)
       allocate(dest%count(n))
       dest%count(:)       =  src%count(:)
    endif

    !dbg    print *,'before dupiodesc2'
    call dupiodesc2(src%read, dest%read)
    call dupiodesc2(src%write, dest%write)
    !dbg    print *,'after dupiodesc2'

    dest%basetype = src%basetype

    if(associated(src%dest_ioproc)) then 
       n = size(src%dest_ioproc)
       allocate(dest%dest_ioproc(n))
       dest%dest_ioproc(:) = src%dest_ioproc(:)
    endif

    if(associated(src%dest_ioindex)) then 
       n = size(src%dest_ioindex)
       allocate(dest%dest_ioindex(n))
       dest%dest_ioindex(:) = src%dest_ioindex(:)
    endif

    if(associated(src%rfrom)) then 
       n = size(src%rfrom)
       allocate(dest%rfrom(n))
       dest%rfrom(:) = src%rfrom(:)
    endif

    if(associated(src%rtype)) then 
       n = size(src%rtype)
       allocate(dest%rtype(n))
       dest%rtype(:) = src%rtype(:)
    endif

    if(associated(src%scount)) then 
       n = size(src%scount)
       allocate(dest%scount(n))
       dest%scount(:) = src%scount(:)
    endif

    if(associated(src%stype)) then 
       n = size(src%stype)
       allocate(dest%stype(n))
       dest%stype(:) = src%stype(:)
    endif

    call copy_decompmap(src%iomap,dest%iomap)
    call copy_decompmap(src%compmap,dest%compmap)

    dest%compsize = src%compsize


  end subroutine dupiodesc

  !=============================================
  !  copy_decompmap:
  !
  !   copy decompmap_t data structures
  !
  !=============================================

  subroutine copy_decompmap(src,dest)
    use pio_types, only : decompmap_t
    type (decompmap_t), intent(in) :: src
    type (decompmap_t), intent(inout) :: dest


    dest%start    = src%start
    dest%length   = src%length

  end subroutine copy_decompmap

!> 
!! @public 
!! @ingroup PIO_setiotype
!! @brief sets the desired type of io to perform
!! @details
!! @param file @copydoc file_desc_t
!! @param iotype : @copydoc PIO_iotype
!! @param rearr : @copydoc PIO_rearr_method
!<
  subroutine setiotype(file,iotype,rearr)

    type (file_desc_t), intent(inout) :: file
    integer(i4), intent(in) :: iotype 
    integer(i4), intent(in) :: rearr

    file%iotype = iotype
    file%iosystem%rearr = rearr

  end subroutine setiotype

!>
!! @public
!! @ingroup PIO_numtoread
!! @brief returns the global number of words to read for this io descriptor
!! @details
!! @param iodesc : @copydoc io_desc_t
!! @retval num   :  the number of words to read 
!<
  integer function numtoread(iodesc) result(num)

    type (io_desc_t) :: iodesc

    num = iodesc%read%n_words

  end function numtoread

!>
!! @public
!! @ingroup PIO_numtowrite
!! @brief returns the global number of words to write for this io descriptor
!! @details
!! @param iodesc : @copydoc io_desc_t
!<
  integer function numtowrite(iodesc) result(num)

    type (io_desc_t) :: iodesc

    num = iodesc%write%n_words

  end function numtowrite

!> 
!! @public
!! @ingroup PIO_createfile 
!! @brief create a file using pio
!! @details  Input parameters are read on comp task 0 and ignored elsewhere
!! @param iosystem : a defined pio system descriptor created by a call to @ref PIO_init (see PIO_types)
!! @param file	:  the returned file descriptor
!! @param iotype : @copydoc PIO_iotype
!! @param fname : the name of the file to open
!! @param amode_in : the creation mode flag. the following flags are available: PIO_clobber, PIO_noclobber. 
!! @retval ierr @copydoc error_return
!<
  integer function createfile(iosystem, file,iotype, fname, amode_in) result(ierr)
    type (iosystem_desc_t), intent(inout), target :: iosystem
    type (file_desc_t), intent(out) :: file
    integer, intent(in) :: iotype
    character(len=*), intent(in)  :: fname
    integer, optional, intent(in) :: amode_in
    
    ! ===================
    !  local variables
    ! ===================
    logical :: iscallback
    integer    :: amode
    integer :: msg
    logical, parameter :: check = .true.
    character(len=9) :: rd_buffer
    character(len=char_len)  :: myfname
#ifdef TIMING
    call t_startf("PIO_createfile")
#endif

    if(debug.or.debugasync) print *,'createfile: {comp,io}_rank:',iosystem%comp_rank,iosystem%io_rank,'io proc: ',iosystem%ioproc, iosystem%async_interface, iotype
    ierr=PIO_noerr
    

    if(present(amode_in)) then
       amode = amode_in
    else	
       amode = 0
    end if

    file%iotype = iotype 
    myfname = fname

    if(.not. (iosystem%async_interface .and. iosystem%ioproc)) then
       call mpi_bcast(amode, 1, MPI_INTEGER, 0, iosystem%comp_comm, ierr)
       call mpi_bcast(file%iotype, 1, MPI_INTEGER, 0, iosystem%comp_comm, ierr)

       if(len(fname) > char_len) then
          print *,'Length of filename exceeds compile time max, increase char_len in pio_kinds and recompile', len(fname), char_len
          call piodie( _FILE_,__LINE__)
       end if

       call mpi_bcast(myfname, len(fname), mpi_character, 0, iosystem%comp_comm, ierr)
    end if

    file%iosystem => iosystem

    !--------------------------------
    ! set some iotype specific stuff
    !--------------------------------

#if defined(USEMPIIO) 
    if ( (file%iotype==pio_iotype_pbinary .or. file%iotype==pio_iotype_direct_pbinary) &
         .and. (.not. iosystem%userearranger) ) then
       write(rd_buffer,('(i9)')) 16*1024*1024
       call PIO_set_hint(iosystem, "cb_buffer_size",trim(adjustl(rd_buffer)))
    endif
#endif

#ifndef _NETCDF4
    if(file%iotype==pio_iotype_netcdf4p .or. file%iotype==pio_iotype_netcdf4c) then
       print *, 'WARNING: PIO was not built with NETCDF 4 support changing iotype to netcdf'
       file%iotype = pio_iotype_netcdf
    end if
#endif
    if(iosystem%async_interface .and. .not. iosystem%ioproc) then
       msg = PIO_MSG_CREATE_FILE
       if(iosystem%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, iosystem%ioroot, 1, iosystem%union_comm, ierr)
       end if

       call mpi_bcast(myfname, char_len, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
       call mpi_bcast(iotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
       call mpi_bcast(amode, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

    end if
    select case(iotype)
    case(pio_iotype_pbinary, pio_iotype_direct_pbinary)
       if(present(amode_in)) then
          print *, 'warning, the mode argument is currently ignored for binary file operations'
       end if
       ierr = create_mpiio(file,myfname)
    case( pio_iotype_pnetcdf, pio_iotype_netcdf, pio_iotype_netcdf4p, pio_iotype_netcdf4c)
       ierr = create_nf(file,myfname, amode)	
       if(debug .and. iosystem%io_rank==0)print *,_FILE_,__LINE__,' open: ', myfname, file%fh
    case(pio_iotype_binary)
       print *,'createfile: io type not supported'
    end select
    if(ierr==0) file%file_is_open=.true.

    if(debug .and. file%iosystem%io_rank==0) print *,_FILE_,__LINE__,'open: ',file%fh, myfname

#ifdef TIMING
    call t_stopf("PIO_createfile")
#endif
  end function createfile

!> 
!! @public
!! @ingroup PIO_openfile 
!! @brief open an existing file using pio
!! @details  Input parameters are read on comp task 0 and ignored elsewhere.
!! @param iosystem : a defined pio system descriptor created by a call to @ref PIO_init (see PIO_types)
!! @param file	:  the returned file descriptor
!! @param iotype : @copydoc PIO_iotype
!! @param fname : the name of the file to open
!! @param mode : a zero value (or PIO_nowrite) specifies the default
!! behavior: open the dataset with read-only access, buffering and
!! caching accesses for efficiency otherwise, the creation mode is
!! PIO_write. setting the PIO_write flag opens the dataset with
!! read-write access. ("writing" means any kind of change to the dataset,
!! including appending or changing data, adding or renaming dimensions,
!! variables, and attributes, or deleting attributes.) 
!! @retval ierr @copydoc error_return
!<
  integer function PIO_openfile(iosystem, file, iotype, fname,mode) result(ierr)
    type (iosystem_desc_t), intent(inout), target :: iosystem
    type (file_desc_t), intent(out) :: file
    integer, intent(in) :: iotype
    character(len=*), intent(in)  :: fname
    integer, optional, intent(in) :: mode

    ! ===================
    !  local variables
    ! ================
    integer    :: amode, msg
    logical, parameter :: check = .true.
    character(len=9) :: rd_buffer
    character(len=char_len) :: myfname

#ifdef TIMING
    call t_startf("PIO_openfile")
#endif



    if(Debug .or. Debugasync) print *,'PIO_openfile: {comp,io}_rank:',iosystem%comp_rank,iosystem%io_rank,'io proc: ',iosystem%ioproc
    ierr=PIO_noerr

    file%iosystem => iosystem

    if(present(mode)) then
       amode = mode
    else	
       amode = 0
    end if
    !--------------------------------
    ! set some iotype specific stuff
    !--------------------------------

    if(iosystem%num_iotasks.eq.1.and.iotype.eq.pio_iotype_pnetcdf) then	
#if defined(_NETCDF)
       file%iotype=pio_iotype_netcdf
#else
       file%iotype = iotype 
#endif       
    else
       file%iotype = iotype 
    end if


    myfname = fname

#if defined(USEMPIIO)
    if ( (file%iotype==pio_iotype_pbinary .or. file%iotype==pio_iotype_direct_pbinary) &
         .and. (.not. iosystem%userearranger) ) then
       write(rd_buffer,('(i9)')) 16*1024*1024
       call PIO_set_hint(iosystem, "cb_buffer_size",trim(adjustl(rd_buffer)))
    endif
#endif
#ifndef _NETCDF4
    if(file%iotype==pio_iotype_netcdf4p .or. file%iotype==pio_iotype_netcdf4c) then
       print *, 'WARNING: PIO was not built with NETCDF 4 support changing iotype to netcdf'
       file%iotype = pio_iotype_netcdf
    end if
#endif
    if(.not. (iosystem%ioproc .and. iosystem%async_interface)) then
       call mpi_bcast(amode, 1, MPI_INTEGER, 0, iosystem%comp_comm, ierr)
       call mpi_bcast(file%iotype, 1, MPI_INTEGER, 0, iosystem%comp_comm, ierr)
       if(len(fname) > char_len) then
          print *,'Length of filename exceeds compile time max, increase char_len in pio_kinds and recompile'
          call piodie( _FILE_,__LINE__)
       end if

       call mpi_bcast(myfname, len(fname), mpi_character, 0, iosystem%comp_comm, ierr)
    end if

    if(iosystem%async_interface .and. .not. iosystem%ioproc) then
       msg = PIO_MSG_OPEN_FILE
       if(iosystem%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, iosystem%ioroot, 1, iosystem%union_comm, ierr)
       end if
       
       call mpi_bcast(myfname, char_len, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
       call mpi_bcast(iotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
       call mpi_bcast(amode, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
    end if

    select case(iotype)
    case(pio_iotype_pbinary, pio_iotype_direct_pbinary)
       if(amode /=0) then
          print *, 'warning, the mode argument is currently ignored for binary file operations'
       end if
       ierr = open_mpiio(file,myfname)
    case( pio_iotype_pnetcdf, pio_iotype_netcdf, pio_iotype_netcdf4c, pio_iotype_netcdf4p)
       ierr = open_nf(file,myfname,amode)
       if(debug .and. iosystem%io_rank==0)print *,_FILE_,__LINE__,' open: ', myfname, file%fh
    case(pio_iotype_binary)   ! appears to be a no-op
       
    end select
    if(Debug .and. file%iosystem%io_rank==0) print *,_FILE_,__LINE__,'open: ',file%fh, myfname
    if(ierr==0) file%file_is_open=.true.
#ifdef TIMING
    call t_stopf("PIO_openfile")
#endif

  end function PIO_openfile

!> 
!! @public 
!! @ingroup PIO_syncfile 
!! @brief synchronizing a file forces all writes to complete before the subroutine returns. 
!!
!! @param file @copydoc file_desc_t
!<
  subroutine syncfile(file)
    implicit none
    type (file_desc_t), target :: file
    integer :: ierr, msg
    type(iosystem_desc_t), pointer :: ios
     
 
    ios => file%iosystem
    if(ios%async_interface .and. .not. ios%ioproc) then
       msg = PIO_MSG_SYNC_FILE
       if(ios%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
       end if
      
       call mpi_bcast(file%fh, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
    end if

    select case(file%iotype)
    case( pio_iotype_pnetcdf, pio_iotype_netcdf)
       ierr = sync_nf(file)
    case(pio_iotype_pbinary, pio_iotype_direct_pbinary)
    case(pio_iotype_binary) 
    end select
  end subroutine syncfile
!> 
!! @public 
!! @ingroup PIO_freedecomp
!! @brief free all allocated storage associated with this decomposition
!! @details
!! @param ios :  a defined pio system descriptor created by call to @ref PIO_init (see PIO_types)
!! @param iodesc @copydoc io_desc_t
!<
  subroutine freedecomp_ios(ios,iodesc)
    implicit none
    type (iosystem_desc_t) :: ios
    type (io_desc_t) :: iodesc
    integer :: ierr


    call rearrange_free(ios,iodesc)

#ifndef _MPISERIAL
    if(ios%ioproc) then
!       if(debug) print *,__FILE__,__LINE__,iodesc%write%n_elemtype,iodesc%write%n_words,iodesc%write%elemtype,iodesc%write%filetype

       if((iodesc%read%filetype .ne. mpi_datatype_null)  &
	  .and. (iodesc%read%filetype .ne. iodesc%write%filetype) .and. &
	  iodesc%read%n_words>0) then 
          call mpi_type_free(iodesc%read%filetype,ierr)
          call checkmpireturn('freedecomp mpi_type_free: ',ierr)
          call mpi_type_free(iodesc%read%elemtype,ierr)
          call checkmpireturn('freedecomp mpi_type_free: ',ierr)
          iodesc%read%filetype=mpi_datatype_null
       endif
       if(iodesc%write%filetype .ne. mpi_datatype_null .and. &
	  iodesc%write%n_words>0) then 
          call mpi_type_free(iodesc%write%filetype,ierr)
          call checkmpireturn('freedecomp mpi_type_free: ',ierr)
          call mpi_type_free(iodesc%write%elemtype,ierr)
          call checkmpireturn('freedecomp mpi_type_free: ',ierr)
          iodesc%write%filetype=mpi_datatype_null
       endif
   
    end if
#endif

    if(associated(iodesc%start)) then
       call dealloc_check(iodesc%start,'iodesc%start')
       nullify(iodesc%start)
    end if

    if(associated(iodesc%count)) then
       call dealloc_check(iodesc%count,'iodesc%count')    
       nullify(iodesc%count)
    end if
  end subroutine freedecomp_ios
!>
!! @public 
!! @ingroup PIO_freedecomp
!! @brief free all allocated storage associated with this decomposition
!! @details
!! @param file @copydoc file_desc_t
!! @param iodesc : @copydoc io_desc_t
!! @retval ierr @copydoc error_return
!<
  subroutine freedecomp_file(file,iodesc)
    implicit none
    type (file_desc_t) :: file
    type (io_desc_t) :: iodesc

    call freedecomp_ios(file%iosystem, iodesc)

  end subroutine freedecomp_file

!> 
!! @public
!! @ingroup PIO_closefile
!! @brief close a disk file
!! @details
!! @param file @copydoc file_desc_t
!< 
  subroutine closefile(file)

    type (file_desc_t),intent(inout)   :: file

    integer :: ierr, msg
    integer :: iotype 
    logical, parameter :: check = .true.

#ifdef TIMING
    call t_startf("PIO_closefile")
#endif
    if(file%iosystem%async_interface .and. .not. file%iosystem%ioproc) then
       msg = PIO_MSG_CLOSE_FILE
       if(file%iosystem%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, file%iosystem%ioroot, 1, file%iosystem%union_comm, ierr)
       end if
       call mpi_bcast(file%fh, 1, mpi_integer, file%iosystem%compmaster, file%iosystem%intercomm, ierr)
    end if

    if(debug .and. file%iosystem%io_rank==0) print *,_FILE_,__LINE__,'close: ',file%fh
    iotype = file%iotype 
    select case(iotype)
    case(pio_iotype_pbinary, pio_iotype_direct_pbinary)
       ierr = close_mpiio(file)
    case( pio_iotype_pnetcdf, pio_iotype_netcdf, pio_iotype_netcdf4p, pio_iotype_netcdf4c)
       ierr = close_nf(file)
    case(pio_iotype_binary)
       print *,'closefile: io type not supported'
    end select
    if(ierr==0) file%file_is_open=.false.

#ifdef TIMING
    call t_stopf("PIO_closefile")
#endif


  end subroutine closefile


  !******************************
  ! read_ascii
  !

  subroutine read_ascii(rank,iobuf,size)

    integer, intent(in) :: rank
    real (r8), dimension(:) :: iobuf
    integer, intent(in) :: size

    character(len=80) filename
    integer lun
    integer ios
    integer i

    lun=10+rank
    write(filename,"('fort.',i2)" ) lun
    write(6,*) 'filename is:', filename

    open(lun,file=filename,status='old',iostat=ios)
    if (ios /= 0) then
       write(6,*) rank,': could not open ascii file: ',filename
    endif

    do i=1,size
       read(unit=lun,fmt=*,iostat=ios) iobuf(i)
       if (ios /= 0) then
          write (6,*) rank,': error reading item ',i,' of ',size
          call abort
       endif

    end do

    close(lun)

  end subroutine read_ascii



  !*****************************
  ! calcdisplace
  !

  subroutine calcdisplace(bsize,dof,displace)

    integer(i4), intent(in) :: bsize    ! length of contigious blocks of numbers
    integer(i4), intent(in) :: dof(:)   ! degree of freedom on which to setup the displacement array
    integer(i4), intent(inout) :: displace(:)  ! array of mpi displacments

    integer :: numblocks,lenblocks,i,ii,dis

    numblocks = size(displace)
    lenblocks = bsize
    do i=1,numblocks
       ii = (i-1)*lenblocks+1
       dis = dof(ii)-1
       dis = dis/lenblocks   
       displace(i) = dis
    enddo
    do i=1,numblocks-1	
       if(displace(i+1) .lt. displace(i)) then
          print *,'calcdisplace: error with displacement arrays',i,displace(i:i+1),numblocks,size(dof),dof(numblocks)
          call piodie( _FILE_,__LINE__)
       endif
    enddo

  end subroutine calcdisplace


  subroutine calcdisplace_box(gsize,start,count,ndim,displace)

    integer(i4),intent(in) :: gsize(:)   ! global size of output domain
    integer(kind=PIO_offset),intent(in) :: start(:), count(:)
    integer(i4), intent(in) :: ndim
    integer(i4),intent(inout) :: displace(:)  ! mpi displacments

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

    iosize=min(int(count(1)),1)
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
    !       divide by lenblocks

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

    do i=1,ndisp-1	
       if(displace(i) .gt. displace(i+1)) then
          call piodie(_FILE_,__LINE__,'displace is not increasing')
       endif
    enddo

  end subroutine calcdisplace_box


end module piolib_mod

  !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
