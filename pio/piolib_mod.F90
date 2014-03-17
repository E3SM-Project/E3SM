#define __PIO_FILE__ "piolib_mod.f90"
#define debug_rearr 0
#ifdef BGP
#define BGx
#endif
#ifdef BGL
#define BGx
#endif
#ifdef BGQ
#define BGx
#endif
!>
!! @file 
!! @brief Initialization Routines for PIO
!! 
!! $Revision$
!! $LastChangedDate$
!<
module piolib_mod
  use iso_c_binding
  !--------------
  use pio_kinds
  !--------------
  use pio_types, only : file_desc_t, iosystem_desc_t, var_desc_t, io_desc_t, &
	pio_iotype_pbinary, pio_iotype_binary, pio_iotype_direct_pbinary, &
	pio_iotype_netcdf, pio_iotype_pnetcdf, pio_iotype_netcdf4p, pio_iotype_netcdf4c, &
        pio_noerr
  !--------------
  use pio_support, only : piodie, debug, debugio, debugasync, checkmpireturn
  !

  use pio_mpi_utils, only : PIO_type_to_mpi_type 
#ifdef TIMING
  use perf_mod, only : t_startf, t_stopf     ! _EXTERNAL
#endif
#ifndef NO_MPIMOD
  use mpi    ! _EXTERNAL
#endif
  implicit none
  private
#ifdef NO_MPIMOD
  include 'mpif.h'    ! _EXTERNAL
#endif
  ! !public member functions:

  public :: PIO_init,     &
       PIO_finalize,      &
       PIO_initdecomp,    &
       PIO_openfile,      &
       PIO_syncfile,      &
       PIO_createfile,    &
       PIO_closefile,     &
       PIO_setframe,      &
       PIO_advanceframe,  &
       PIO_setdebuglevel, &
       PIO_seterrorhandling, &
       PIO_get_local_array_size, &
       PIO_freedecomp,     &
       PIO_getnumiotasks, &
       PIO_set_hint,      &
       PIO_FILE_IS_OPEN, &
       PIO_deletefile, &
       PIO_get_numiotasks, &
       PIO_get_iorank


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
     module procedure PIO_initdecomp_dof_i4  ! previous name: initdecomop_1dof_nf_box
     module procedure PIO_initdecomp_dof_i8  ! previous name: initdecomop_1dof_nf_box
     module procedure initdecomp_1dof_nf_i4
     module procedure initdecomp_1dof_nf_i8
     module procedure initdecomp_1dof_bin_i4
     module procedure initdecomp_1dof_bin_i8
     module procedure initdecomp_2dof_nf_i4
     module procedure initdecomp_2dof_nf_i8
     module procedure initdecomp_2dof_bin_i4
     module procedure initdecomp_2dof_bin_i8
!     module procedure PIO_initdecomp_bc
!     module procedure PIO_initdecomp_dof_dof
  end interface

!> 
!! @defgroup PIO_dupiodesc PIO_dupiodesc
!! duplicates an eisting io descriptor
!<
!  interface PIO_dupiodesc
!     module procedure dupiodesc
!  end interface



!> 
!! @defgroup PIO_numtoread PIO_numtoread
!! returns the total number of words to read
!<
!  interface PIO_numtoread
!     module procedure numtoread
!  end interface

!> 
!! @defgroup PIO_numtowrite PIO_numtowrite
!! returns the total number of words to write
!<
!  interface PIO_numtowrite
!     module procedure numtowrite
!  end interface


!> 
!! @defgroup PIO_getnumiotasks PIO_getnumiotasks
!!  returns the actual number of IO-tasks used.  PIO 
!!  will reset the total number of IO-tasks if certain 
!!  conditions are meet
!<
  interface PIO_get_numiotasks
     module procedure getnumiotasks
  end interface
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

contains
!> 
!! @public 
!! @ingroup PIO_file_is_open
!! @brief This logical function indicates if a file is open.
!! @details
!! @param File @copydoc file_desc_t
!<
  logical function PIO_FILE_IS_OPEN(File)
    type(file_desc_t), intent(in) :: file
    interface
       integer(C_INT) function PIOc_File_is_Open(ncid) &
            bind(C,NAME="PIOc_File_is_Open")
         use iso_c_binding
         implicit none
         integer(c_int), value :: ncid
       end function PIOc_File_is_Open
    end interface

    if(PIOc_File_is_Open(file%fh)==1) then
       PIO_FILE_IS_OPEN = .true.
    else
       PIO_FILE_IS_OPEN = .false.
    endif


  end function PIO_FILE_IS_OPEN


!> 
!! @public 
!! @ingroup PIO_get_local_array_size
!! @brief This function returns the expected local size of an array associated with iodesc
!! @details
!! @param iodesc 
!! @copydoc io_desc_t
!<
  integer function PIO_get_local_array_size(iodesc)
    type(io_desc_t), intent(in) :: iodesc   
    interface
       integer(C_INT) function PIOc_get_local_array_size(ioid) &
            bind(C,NAME="PIOc_get_local_array_size")
         use iso_c_binding
         implicit none
         integer(C_INT), value :: ioid
       end function PIOc_get_local_array_size
    end interface
    PIO_get_local_array_size = PIOc_get_local_array_size(iodesc%ioid)
  end function PIO_get_local_array_size

!> 
!! @public 
!! @ingroup PIO_advanceframe
!! @brief advances the record dimension of a variable in a netcdf format file 
!!  or the block address in a binary file
!! @details
!! @param[in,out] vardesc @copybrief var_desc_t 
!<
  subroutine advanceframe(file, vardesc)
    type(file_desc_t), intent(in) :: file
    type(var_desc_t), intent(inout) :: vardesc
    integer ierr;
    interface
       integer(C_INT) function PIOc_advanceframe(fileid, varid) &
            bind(C,NAME="PIOc_advanceframe")
         use iso_c_binding
         implicit none
         integer(C_INT), value :: fileid
         integer(C_INT), value :: varid
       end function PIOc_advanceframe
    end interface
    ierr = PIOc_advanceframe(file%fh, vardesc%varid-1) 
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
  subroutine setframe(file, vardesc,frame)
    type(file_desc_t) :: file
    type(var_desc_t), intent(inout) :: vardesc
    integer(PIO_OFFSET_KIND), intent(in) :: frame
    integer :: ierr, iframe
    interface
       integer(C_INT) function PIOc_setframe(ncid, varid, frame) &
            bind(C,NAME="PIOc_setframe")
         use iso_c_binding
         implicit none
         integer(C_INT), value :: ncid
         integer(C_INT), value :: varid
         integer(C_INT), value :: frame
       end function PIOc_setframe
    end interface
    iframe = frame-1
    ierr = PIOc_setframe(file%fh, vardesc%varid-1, iframe)
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
    type(iosystem_desc_t), intent(inout) :: ios
    integer, intent(in) :: method

    interface
       integer(c_int) function PIOc_Set_IOSystem_Error_Handling(ios, method) &
            bind(C,name="PIOc_Set_IOSystem_Error_Handling")
         use iso_c_binding
         integer(c_int), value :: ios
         integer(c_int), value :: method
       end function PIOc_Set_IOSystem_Error_Handling
    end interface
    integer :: ierr

    ierr = PIOc_Set_IOSystem_Error_Handling(ios%iosysid, method)

  end subroutine seterrorhandlingi

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
  subroutine initdecomp_2dof_bin_i4(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodofw,iodesc)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4)                       :: basetype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenblocks
    integer (i4), intent(in)          :: compdof(:)   !> global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: iodofr(:)     !> global degrees of freedom for io decomposition 
    integer (i4), intent(in)          :: iodofw(:)     !> global degrees of freedom for io decomposition 
    type (io_desc_t), intent(inout)     :: iodesc


    call initdecomp_2dof_bin_i8(iosystem,basepiotype,dims,lenblocks,int(compdof,PIO_OFFSET_KIND),int(iodofr,PIO_OFFSET_KIND), &
         int(iodofw,PIO_OFFSET_KIND),iodesc)


  end subroutine initdecomp_2dof_bin_i4
  subroutine initdecomp_2dof_bin_i8(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodofw,iodesc)
!    use calcdisplace_mod, only : calcdisplace
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4)                       :: basetype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenblocks
    integer (PIO_OFFSET_KIND), intent(in)          :: compdof(:)   !> global degrees of freedom for computational decomposition
    integer (PIO_OFFSET_KIND), intent(in)          :: iodofr(:)     !> global degrees of freedom for io decomposition 
    integer (PIO_OFFSET_KIND), intent(in)          :: iodofw(:)     !> global degrees of freedom for io decomposition 
    type (io_desc_t), intent(inout)     :: iodesc
#ifdef DOTHIS
    integer(PIO_OFFSET_KIND) :: start(1), count(1)

    integer (i4) :: i,ndims,n_iotasks
    integer(PIO_OFFSET_KIND) glength
    logical :: userearranger
    integer (PIO_OFFSET_KIND) ::  ndispr,ndispw
    integer (PIO_OFFSET_KIND) :: lengthr, lengthw
    integer (PIO_OFFSET_KIND), pointer :: displacer(:),displacew(:)


    nullify(iodesc%start)
    nullify(iodesc%count)

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
    glength= product(int(dims,PIO_OFFSET_KIND))
    if(glength > int(huge(i),PIO_OFFSET_KIND)) then
       call piodie( __PIO_FILE__,__LINE__, &
            'requested array size too large for this interface ')       
    endif



    lengthr = size(iodofr);
    lengthw = size(iodofw)
    if(lenblocks>0) then
       ndispw=size(iodofw)/lenblocks 
       ndispr=size(iodofr)/lenblocks
    else
       ndispw=0
       ndispr=0
    end if
    call alloc_check(displacer,int(ndispr))
    call alloc_check(displacew,int(ndispw))

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
       call genindexedblock(lenblocks,basetype,iodesc%read%elemtype,iodesc%read%filetype,int(displacer))

       !-------------------------------------------------
       ! setup the data structure for the write operation
       !-------------------------------------------------
       iodesc%write%n_elemtype = ndispw
       iodesc%write%n_words    = iodesc%write%n_elemtype*lenblocks

       call genindexedblock(lenblocks,basetype,iodesc%write%elemtype,iodesc%write%filetype,int(displacew))

       if(debug) print *,'initdecomp: at the end of subroutine'
       !       if(iodesc%read%n_elemtype == 0 .and. iodesc%write%n_elemtype == 0) iosystem%ioproc = .false.
    endif

    deallocate(displacer,displacew)
#endif

  end subroutine initdecomp_2dof_bin_i8


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
  subroutine initdecomp_1dof_bin_i8(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodesc)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer(i4), intent(in)          :: lenblocks
    integer(PIO_OFFSET_KIND), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer(PIO_OFFSET_KIND), intent(in)          :: iodofr(:)     ! global degrees of freedom for io decomposition 
    type (io_desc_t), intent(inout)     :: iodesc

    integer(PIO_OFFSET_KIND) :: start(1), count(1)
    ! these are not used in the binary interface

    start(1)=-1
    count(1)=-1
    call initdecomp_1dof_nf_i8(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,start, count, iodesc)
  end subroutine initdecomp_1dof_bin_i8

  subroutine initdecomp_1dof_bin_i4(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodesc)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenblocks
    integer (i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: iodofr(:)     ! global degrees of freedom for io decomposition 
    type (io_desc_t), intent(inout)     :: iodesc

    integer(PIO_OFFSET_KIND) :: start(1), count(1)
    ! these are not used in the binary interface

    start(1)=-1
    count(1)=-1
    call initdecomp_1dof_nf_i8(iosystem,basepiotype,dims,lenblocks, &
         int(compdof,PIO_OFFSET_KIND),int(iodofr,PIO_OFFSET_KIND),start, count, iodesc)
  end subroutine initdecomp_1dof_bin_i4

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
  subroutine initdecomp_2dof_nf_i4(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodofw,start, count, iodesc)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenblocks
    integer (i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: iodofr(:)     ! global degrees of freedom for io decomposition 
    integer (i4), intent(in)          :: iodofw(:)     ! global degrees of freedom for io decomposition 

    type (io_desc_t), intent(inout)     :: iodesc

    integer(PIO_OFFSET_KIND), intent(in) :: start(:), count(:)
    type (io_desc_t) :: tmp


    call pio_initdecomp(iosystem, basepiotype,dims,lenblocks,int(compdof,PIO_OFFSET_KIND),int(iodofr,PIO_OFFSET_KIND), &
         int(iodofw,PIO_OFFSET_KIND),start,count,iodesc)

  end subroutine initdecomp_2dof_nf_i4

  subroutine initdecomp_2dof_nf_i8(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodofw,start, count, iodesc)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenblocks
    integer (PIO_OFFSET_KIND), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (PIO_OFFSET_KIND), intent(in)          :: iodofr(:)     ! global degrees of freedom for io decomposition 
    integer (PIO_OFFSET_KIND), intent(in)          :: iodofw(:)     ! global degrees of freedom for io decomposition 

    type (io_desc_t), intent(inout)     :: iodesc

    integer(PIO_OFFSET_KIND), intent(in) :: start(:), count(:)
    type (io_desc_t) :: tmp
    integer :: ierr

    call initdecomp_1dof_nf_i8(iosystem, basepiotype, dims, lenblocks, compdof, iodofr, start, count, iodesc)

    call initdecomp_1dof_nf_i8(iosystem, basepiotype, dims, lenblocks, compdof, iodofw, start, count, tmp)
    call mpi_abort(mpi_comm_world, 0, ierr)
!    call dupiodesc2(iodesc%write,tmp%write)

  end subroutine initdecomp_2dof_nf_i8

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
  subroutine initdecomp_1dof_nf_i4(iosystem,basepiotype,dims,lenblocks,compdof,iodof,start, count, iodesc)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in) :: lenblocks
    integer (i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (i4), intent(in)          :: iodof(:)     ! global degrees of freedom for io decomposition 
    type (io_desc_t), intent(inout)     :: iodesc
    integer :: piotype	
    integer(PIO_OFFSET_KIND), intent(in) :: start(:), count(:)

    call initdecomp_1dof_nf_i8(iosystem, basepiotype,dims,lenblocks,int(compdof,PIO_OFFSET_KIND),int(iodof,PIO_OFFSET_KIND),&
         start,count,iodesc)

  end subroutine initdecomp_1dof_nf_i4
  subroutine initdecomp_1dof_nf_i8(iosystem,basepiotype,dims,lenblocks,compdof,iodof,start, count, iodesc)
!    use calcdisplace_mod, only : calcdisplace
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in) :: lenblocks
    integer (PIO_OFFSET_KIND), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (PIO_OFFSET_KIND), intent(in)          :: iodof(:)     ! global degrees of freedom for io decomposition 
    type (io_desc_t), intent(inout)     :: iodesc
    integer :: piotype
    integer(PIO_OFFSET_KIND), intent(in) :: start(:), count(:)
#ifdef DOTHIS
    integer(i4) :: length,n_iotasks
    integer(i4) :: ndims

    integer (PIO_OFFSET_KIND), pointer :: displace(:)  ! the displacements for the mpi data structure (read)

    integer(i4) :: prev
    integer(PIO_OFFSET_KIND) :: glength    ! global length in words
    integer(i4) :: ii,i,dis,ierr
    integer(i4),pointer, dimension(:) :: blocklen,disp
    logical(log_kind) ::  userearranger
    logical, parameter :: check = .true.
    integer(PIO_OFFSET_KIND) :: ndisp
#ifdef MEMCHK
    integer :: msize, rss, mshare, mtext, mstack
#endif
    nullify(iodesc%start)
    nullify(iodesc%count)

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
    call GPTLget_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,__PIO_FILE__,__LINE__,'mem=',rss
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
    glength= product(int(dims,PIO_OFFSET_KIND))
    if(glength > huge(ndisp)) then
       print *,__FILE__,__LINE__,dims,glength
       call piodie( __PIO_FILE__,__LINE__, &
            'requested array size too large for this interface ')       
    endif

    if(lenblocks>0) then
       ndisp=size(iodof)/lenblocks
    else
       ndisp=0
    end if
    call alloc_check(displace,int(ndisp))

#ifdef MEMCHK	
    call GPTLget_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,__PIO_FILE__,__LINE__,'mem=',rss
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
    call GPTLget_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,__PIO_FILE__,__LINE__,'mem=',rss
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
             call piodie( __PIO_FILE__,__LINE__, &
                  'this interface does not use rearranger')
       
    endif
#ifdef MEMCHK	
    call GPTLget_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,__PIO_FILE__,__LINE__,'mem=',rss
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

       call genindexedblock(lenblocks,piotype,iodesc%write%elemtype,iodesc%write%filetype,int(displace))

       

       if(debug) print *,'initdecomp: at the end of subroutine',iodesc%write%n_elemtype,iodesc%write%n_words
    endif
#ifdef MEMCHK	
    call GPTLget_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,__PIO_FILE__,__LINE__,'mem=',rss
    end if
#endif
    call dupiodesc2(iodesc%write,iodesc%read)
    if(debug) then
       print *, __PIO_FILE__,__LINE__,iodesc%read%filetype,iodesc%read%elemtype,&
            iodesc%read%n_elemtype,iodesc%read%n_words   
       print *, __PIO_FILE__,__LINE__,iodesc%write%filetype,iodesc%write%elemtype,&
            iodesc%write%n_elemtype,iodesc%write%n_words
    end if
    call dealloc_check(displace)

#ifdef MEMCHK	
    call GPTLget_memusage(msize, rss, mshare, mtext, mstack)
    if(rss>lastrss) then
       lastrss=rss
       print *,__PIO_FILE__,__LINE__,'mem=',rss
    end if
#endif
#ifdef TIMING
    call t_stopf("PIO_initdecomp")
#endif
#endif
  end subroutine initdecomp_1dof_nf_i8

!>
!! @public
!! @ingroup PIO_initdecomp
!! @brief Implements the @ref decomp_dof for PIO_initdecomp (previous name: \b initdecomp_1dof_nf_box)
!! @details  This provides the ability to describe a computational
!! decomposition in PIO using degrees of freedom method. This is
!! a decomposition that can not be easily described using a start
!! and count method (see @ref decomp_dof).
!! Optional parameters for this subroutine allows for the specififcation of
!! io decomposition using iostart and iocount arrays.  If iostart
!! and iocount arrays are not specified by the user, and rearrangement
!! is turned on then PIO will calculate an suitable IO decomposition.
!! Note that this subroutine was previously called \em initdecomp_1dof_nf_box
!! @param iosystem @copydoc iosystem_desc_t
!! @param basepiotype @copydoc use_PIO_kinds
!! @param dims An array of the global length of each dimesion of the variable(s)
!! @param compdof Mapping of the storage order for the computational decomposition to its memory order
!! @param iodesc @copydoc iodesc_generate
!! @param iostart   The start index for the block-cyclic io decomposition
!! @param iocount   The count for the block-cyclic io decomposition
!<
  subroutine PIO_initdecomp_dof_i4(iosystem,basepiotype,dims,compdof, iodesc, iostart, iocount)
    type (iosystem_desc_t), intent(inout) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (PIO_OFFSET_KIND), optional :: iostart(:), iocount(:)
    type (io_desc_t), intent(inout)     :: iodesc
    integer(PIO_OFFSET_KIND), pointer :: internal_compdof(:)
    integer(i4), intent(in)           :: dims(:)

    allocate(internal_compdof(size(compdof)))
    internal_compdof = int(compdof,PIO_OFFSET_KIND)
    
    if(present(iostart) .and. present(iocount) ) then
       call pio_initdecomp_dof_i8(iosystem, basepiotype, dims, internal_compdof, iodesc, iostart, iocount)
    else 
       call pio_initdecomp_dof_i8(iosystem, basepiotype, dims, internal_compdof, iodesc)
    endif
    deallocate(internal_compdof)

  end subroutine PIO_initdecomp_dof_i4


  subroutine PIO_initdecomp_dof_i8(iosystem,basepiotype,dims,compdof, iodesc, iostart, iocount)
    type (iosystem_desc_t), intent(inout) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (PIO_OFFSET_KIND), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (PIO_OFFSET_KIND), optional :: iostart(:), iocount(:)
    type (io_desc_t), intent(inout)     :: iodesc
    integer(c_int) :: ndims
    integer(c_int), dimension(:), allocatable, target :: cdims, cstart, ccount
    integer(PIO_OFFSET_KIND), allocatable :: ccompmap(:)
    interface
       integer(C_INT) function PIOc_InitDecomp(iosysid,basetype,ndims,dims, &
            maplen, compmap, ioidp, iostart, iocount)  &
            bind(C,name="PIOc_InitDecomp")
         use iso_c_binding
         integer(C_INT), value :: iosysid
         integer(C_INT), value :: basetype
         integer(C_INT), value :: ndims
         integer(C_INT) :: dims(*)
         integer(C_INT), value :: maplen
         integer(C_SIZE_T) :: compmap(*)
         integer(C_INT) :: ioidp
         type(C_PTR), value :: iostart
         type(C_PTR), value :: iocount
       end function PIOc_InitDecomp
    end interface
    integer :: ierr, i, maplen
    
#ifdef TIMING
    call t_startf("PIO_initdecomp_dof")
#endif
    call mpi_barrier(MPI_COMM_WORLD,ierr)

    ndims = size(dims)
    allocate(cdims(ndims))
    do i=1,ndims
       cdims(i) = dims(ndims-i+1)
    end do
    maplen = size(compdof)
    allocate(ccompmap(maplen))
    ccompmap = compdof - 1

    if(present(iostart) .and. present(iocount)) then
       allocate(cstart(ndims), ccount(ndims))
       do i=1,ndims
          cstart(i) = iostart(ndims-i+1)
          ccount(i) = iocount(ndims-i+1)
       end do
       ierr = PIOc_InitDecomp(iosystem%iosysid, basepiotype, ndims, cdims, &
            maplen, ccompmap, iodesc%ioid, C_LOC(cstart), C_LOC(ccount))
       deallocate(cstart, ccount)
    else
       ierr = PIOc_InitDecomp(iosystem%iosysid, basepiotype, ndims, cdims, &
            maplen, ccompmap, iodesc%ioid, C_NULL_PTR, C_NULL_PTR)
    end if
    deallocate(cdims)
    deallocate(ccompmap)
#ifdef TIMING
    call t_stopf("PIO_initdecomp_dof")
#endif

  end subroutine PIO_initdecomp_dof_i8


  !************************************
  ! dupiodesc2
  !

!  subroutine dupiodesc2(src, dest)
!    use pio_types, only : io_desc2_t
!    type(io_desc2_t), intent(in) :: src
!    type(io_desc2_t), intent(out) :: dest

!    dest%filetype = src%filetype
!    dest%elemtype = src%elemtype
!    dest%n_elemtype = src%n_elemtype
!    dest%n_words = src%n_words
!  end subroutine dupiodesc2



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

    integer(i4) :: numblocks,i,ierr, prev

    logical, parameter :: check = .true.

    integer:: nints, nadds, ndtypes, comb, lbasetype

    numblocks = size(displace)

    !tcx - allow empty displace array
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
       elemtype = mpi_datatype_null
       filetype = mpi_datatype_null
       call mpi_type_contiguous(lenblocks,lbasetype,elemtype,ierr)
       if(check) call checkmpireturn('genindexedblock: after call to type_contiguous: ',ierr)
       call mpi_type_commit(elemtype,ierr)
       if(check) call checkmpireturn('genindexedblock: after call to type_commit: ',ierr)
       if(numblocks>0) then
         call mpi_type_create_indexed_block(numblocks,1,displace,elemtype,filetype,ierr)
         if(check) call checkmpireturn('genindexedblock: after call to type_create_indexed_block: ',ierr)
         call mpi_type_commit(filetype,ierr)
         if(check) call checkmpireturn('genindexedblock: after call to type_commit: ',ierr)
         if(debug) then
            call mpi_type_get_envelope(filetype, nints, nadds, ndtypes, comb, ierr)
            print *,__FILE__,__LINE__,nints,nadds,ndtypes,comb,ierr
         endif
       endif

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
    integer :: lbase
    integer :: ierr
    interface
       integer(c_int) function PIOc_Init_Intracomm_from_F90(f90_comp_comm, num_iotasks, stride,base,iosysidp) &
            bind(C,name="PIOc_Init_Intracomm_from_F90")
         use iso_c_binding
         integer(C_INT), value :: f90_comp_comm
         integer(C_INT), value :: num_iotasks
         integer(C_INT), value :: stride
         integer(C_INT), value :: base
         integer(C_INT) :: iosysidp
       end function PIOc_Init_Intracomm_from_F90
    end interface

#ifdef TIMING
    call t_startf("PIO_init")
#endif
    lbase=0
    if(present(base)) lbase=base
    ierr = PIOc_Init_Intracomm_from_F90(comp_comm,num_iotasks,stride,lbase,iosystem%iosysid)

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
#ifdef DOTHIS
    integer :: ierr
    logical :: is_inter
    logical, parameter :: check=.true.
  
    integer :: i, j, iam, io_leader, comp_leader
    integer(i4), pointer :: iotmp(:)
    character(len=5) :: cb_nodes
    integer :: itmp
    
#ifdef TIMING
    call t_startf("PIO_init")
#endif
#if defined(NO_MPI2) || defined(_MPISERIAL)
    call piodie( __PIO_FILE__,__LINE__, &
     'The PIO async interface requires an MPI2 complient MPI library')
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
       iosystem(i)%numOST = PIO_num_OST


       if(io_comm/=MPI_COMM_NULL) then
          ! Find the rank of the io leader in peer_comm
          call mpi_comm_rank(io_comm,iosystem(i)%io_rank, ierr)
          if(iosystem(i)%io_rank==0) then 
             call mpi_comm_rank(peer_comm, iam, ierr)
          else
             iam = -1
          end if
          call mpi_allreduce(iam, io_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
          call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,__FILE__,__LINE__)
          ! Find the rank of the comp leader in peer_comm
          iam = -1
          call mpi_allreduce(iam, comp_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
          call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,__FILE__,__LINE__)
          ! create the intercomm
          call mpi_intercomm_create(io_comm, 0, peer_comm, comp_leader, i, iosystem(i)%intercomm, ierr)
          ! create the union_comm
          call mpi_intercomm_merge(iosystem(i)%intercomm, .true., iosystem(i)%union_comm, ierr)
       else
          ! Find the rank of the io leader in peer_comm
          iam = -1
          call mpi_allreduce(iam, io_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
          call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,__FILE__,__LINE__)

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
          call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,__FILE__,__LINE__)

          ! create the intercomm
          call mpi_intercomm_create(comp_comms(i), 0, peer_comm, io_leader, i, iosystem(i)%intercomm, ierr)
          ! create the union comm
          call mpi_intercomm_merge(iosystem(i)%intercomm, .false., iosystem(i)%union_comm, ierr)
       end if
       if(Debugasync) print *,__PIO_FILE__,__LINE__,i, iosystem(i)%intercomm, iosystem(i)%union_comm

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
          
          if(Debugasync) print *,__PIO_FILE__,__LINE__
          
          call MPI_allreduce(iosystem(i)%comproot, j, 1, MPI_INTEGER, MPI_MAX,iosystem(i)%union_comm,ierr)
          call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,__FILE__,__LINE__)
          
          iosystem%comproot=j
          call MPI_allreduce(iosystem(i)%ioroot, j, 1, MPI_INTEGER, MPI_MAX,iosystem(i)%union_comm,ierr)
          call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,__FILE__,__LINE__)

          iosystem%ioroot=j

          if(Debugasync) print *,__PIO_FILE__,__LINE__, i, iosystem(i)%comproot, iosystem(i)%ioroot

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
          if(Debugasync) print *,__PIO_FILE__,__LINE__,iotmp
          call MPI_allreduce(iotmp,iosystem(i)%ioranks,iosystem(i)%num_iotasks,MPI_INTEGER,MPI_MAX,iosystem(i)%union_comm,ierr)
          call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,__FILE__,__LINE__)

          if(Debugasync) print *,__PIO_FILE__,__LINE__,iosystem(i)%ioranks
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
       call mpi_info_create(iosystem(i)%info,ierr)
       ! turn on mpi-io aggregation 
       !DBG    print *,'PIO_init: before call to setnumagg'
!       itmp = num_aggregator
!       call mpi_bcast(itmp, 1, mpi_integer, 0, iosystem%union_comm, ierr)
!       if(itmp .gt. 0) then 
!          write(cb_nodes,('(i5)')) itmp
!#ifdef BGx
!          call PIO_set_hint(iosystem(i),"bgl_nodes_pset",trim(adjustl(cb_nodes)))
!#else
!          call PIO_set_hint(iosystem(i),"cb_nodes",trim(adjustl(cb_nodes)))
!#endif       
!       endif

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

    if(DebugAsync) print*,__PIO_FILE__,__LINE__, iosystem(1)%ioranks


    iosystem%num_aiotasks = iosystem%num_iotasks
    iosystem%numost = PIO_NUM_OST

    ! This routine does not return
    if(io_comm /= MPI_COMM_NULL) call pio_msg_handler(component_count,iosystem) 
    
    if(DebugAsync) print*,__PIO_FILE__,__LINE__, iosystem(1)%ioranks
#ifdef TIMING
    call t_stopf("PIO_init")
#endif
#endif
#endif
  end subroutine init_intercom

!>
!! @public
!! @defgroup PIO_recommend_iotasks PIO_recommend_iotasks
!! @brief Recommend a subset of tasks in comm to use as IO tasks
!! @details  This subroutine will give PIO's best recommendation for the number and
!!    location of iotasks for a given system there is no requirement to follow this recommendation.
!!    Using the recommendation requires that PIO_BOX_RERRANGE be used
!! @param A communicator of mpi tasks to choose from
!! @param miniotasks \em optional The minimum number of IO tasks the caller desires
!! @param maxiotasks \em optional The maximum number of IO tasks the caller desires
!! @param iotask if true pio recommends that this task be used as an iotask
!<

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
    call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,__FILE__,__LINE__)

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
    type (iosystem_desc_t), intent(inout)  :: iosystem  ! io descriptor to initalize
    character(len=*), intent(in) :: hint, hintval
    integer :: ierr
    
    interface
       integer(C_INT) function PIOc_set_hint(iosysid, key, val) &
            bind(C,name="PIOc_set_hint")
         use iso_c_binding
         integer(C_INT), intent(in), value :: iosysid
         character(C_CHAR), intent(in) :: key
         character(C_CHAR), intent(in) :: val
       end function PIOc_set_hint
    end interface


    ierr = PIOc_set_hint(iosystem%iosysid, hint, hintval)


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
    interface
       integer(C_INT) function PIOc_finalize(iosysid) &
            bind(C,name="PIOc_finalize")
         use iso_c_binding
         integer(C_INT), intent(in), value :: iosysid
       end function PIOc_finalize
    end interface

    ierr = PIOc_finalize(iosystem%iosysid)
     
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
       integer :: ierr
       interface
          integer(C_INT) function PIOc_get_numiotasks(iosysid,numiotasks) &
               bind(C,name="PIOc_get_numiotasks")
            use iso_c_binding
            integer(C_INT), intent(in), value :: iosysid
            integer(C_INT), intent(out) :: numiotasks
          end function PIOc_get_numiotasks
       end interface
       ierr = PIOc_get_numiotasks(iosystem%iosysid, numiotasks)

   end subroutine getnumiotasks

   subroutine PIO_get_iorank(iosystem, iorank)
     type(iosystem_desc_t), intent(in) :: iosystem
     integer, intent(out) :: iorank
       integer :: ierr
       interface
          integer(C_INT) function PIOc_get_iorank(iosysid,iorank) &
               bind(C,name="PIOc_get_iorank")
            use iso_c_binding
            integer(C_INT), intent(in), value :: iosysid
            integer(C_INT), intent(out) :: iorank
          end function PIOc_get_iorank
       end interface

       ierr = PIOc_get_iorank(iosystem%iosysid, iorank)
     end subroutine PIO_get_iorank

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

#ifdef WHAT_REASON_FOR_THIS_IN_THE_API
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
#endif
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
    integer :: mode
    interface
       integer(C_INT) function PIOc_createfile(iosysid, fh, iotype, fname,mode) &
         bind(C,NAME='PIOc_createfile')
         use iso_c_binding
         implicit none
         integer(c_int), value :: iosysid
         integer(c_int) :: fh
         integer(c_int) :: iotype
         character(kind=c_char) :: fname
         integer(c_int), value :: mode
       end function PIOc_createfile
    end interface
#ifdef TIMING
    call t_startf("PIO_createfile")
#endif
    mode = 0
    if(present(amode_in)) mode = amode_in
    ierr = PIOc_createfile(iosystem%iosysid, file%fh, iotype, trim(fname)//C_NULL_CHAR, mode)
    file%iosystem => iosystem
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
!! @param iotype : @copybrief PIO_iotype
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
  integer function PIO_openfile(iosystem, file, iotype, fname,mode, CheckMPI) result(ierr)
    type (iosystem_desc_t), intent(inout), target :: iosystem
    type (file_desc_t), intent(out) :: file
    integer, intent(in) :: iotype
    character(len=*), intent(in)  :: fname
    integer, optional, intent(in) :: mode
    logical, optional, intent(in) :: CheckMPI

    interface
       integer(C_INT) function PIOc_openfile(iosysid, fh, iotype, fname,mode, CheckMPI) &
         bind(C,NAME='PIOc_openfile')
         use iso_c_binding
         implicit none
         integer(c_int), value :: iosysid
         integer(c_int) :: fh
         integer(c_int) :: iotype
         character(kind=c_char) :: fname
         integer(c_int), value :: mode
         logical(c_bool), value :: CheckMPI
       end function PIOc_openfile
    end interface
    logical(c_bool) :: iCheckMPI=.true.
    integer :: imode=0
#ifdef TIMING
    call t_startf("PIO_openfile")
#endif
    if(present(Checkmpi)) icheckmpi=logical(Checkmpi,c_bool)
    if(present(mode)) imode = mode
    ierr = PIOc_openfile( iosystem%iosysid, file%fh, iotype, &
         trim(fname)//C_NULL_CHAR, imode, iCheckMPI)

    file%iosystem => iosystem
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
    integer :: ierr
    interface 
       integer(C_INT) function PIOc_sync(ncid) &
            bind(C,name="PIOc_sync")
         use iso_c_binding
         integer(C_INT), intent(in), value :: ncid
       end function PIOc_sync
    end interface

    ierr = PIOc_sync(file%fh)

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
    interface 
       integer(C_INT) function PIOc_freedecomp(iosysid, ioid) &
            bind(C,name="PIOc_freedecomp")
         use iso_c_binding
         integer(C_INT), intent(in), value :: iosysid, ioid
       end function PIOc_freedecomp
    end interface

    ierr = PIOc_freedecomp(ios%iosysid, iodesc%ioid)

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
    type(file_desc_t) :: file
    integer :: ierr
    interface
       integer(c_int) function PIOc_closefile(ncid) &
            bind(C,name="PIOc_closefile")
         use iso_c_binding
         integer(C_INT), value :: ncid
       end function PIOc_closefile
    end interface
#ifdef TIMING
    call t_startf("PIO_closefile")
#endif
    ierr = PIOc_closefile(file%fh)

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
#ifndef CPRNAG
          call abort
#else
          stop
#endif
       endif

    end do

    close(lun)

  end subroutine read_ascii

  subroutine pio_deletefile(ios, fname)
    type(iosystem_desc_t) :: ios
    character(len=*) :: fname
    integer :: ierr
    interface
       integer(c_int) function PIOc_deletefile(iosid, fname) &
            bind(C,name="PIOc_deletefile")
         use iso_c_binding
         integer(C_INT), value :: iosid
         character(kind=c_char) :: fname         
       end function PIOc_deletefile
    end interface

    ierr = PIOc_deletefile(ios%iosysid, trim(fname)//C_NULL_CHAR)

  end subroutine pio_deletefile





end module piolib_mod

  !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
