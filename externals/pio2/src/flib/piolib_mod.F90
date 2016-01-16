#define __PIO_FILE__ "piolib_mod.f90"
#define debug_rearr 0

!>
!! @file 
!! @brief Initialization Routines for PIO
!! 
!<
module piolib_mod
  use iso_c_binding
  !--------------
  use pio_kinds
  !--------------
  use pio_types, only : file_desc_t, iosystem_desc_t, var_desc_t, io_desc_t, &
	pio_iotype_netcdf, pio_iotype_pnetcdf, pio_iotype_netcdf4p, pio_iotype_netcdf4c, &
        pio_noerr, pio_rearr_subset
  !--------------
  use pio_support, only : piodie, debug, debugio, debugasync, checkmpireturn
  !


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
       PIO_iotype_available

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
!! @brief sets the unlimited dimension for netcdf file
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
!! @details initdecomp_1dof_bin_i8, initdecomp_1dof_nf_i4, initdecomp_2dof_bin_i4,
!! and initdecomp_2dof_nf_i4 are all depreciated, but supported for backwards 
!! compatibility. 
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
     module procedure PIO_initdecomp_bc
!     module procedure PIO_initdecomp_dof_dof
  end interface

!> 

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

!!$#ifdef __GFORTRAN__
!!$    pure function fptr ( inArr ) result ( ptr )
!!$        integer (PIO_OFFSET_KIND), dimension(:), target, intent(in) :: inArr
!!$        integer (PIO_OFFSET_KIND), target :: ptr
!!$        ptr = inArr(1)
!!$    end function fptr
!!$#elif CPRNAG
!!$! no-op -- nothing here for nag.
!!$#else
#define fptr(arg) arg
!!$#endif

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
    PIO_FILE_IS_OPEN = .false.
    if(associated(file%iosystem)) then
      if(PIOc_File_is_Open(file%fh)==1) then
        PIO_FILE_IS_OPEN = .true.
      endif
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
!! @param level : default value is 0, allowed values 0-6
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
  subroutine seterrorhandlingf(file, method, oldmethod)
    type(file_desc_t), intent(inout) :: file
    integer, intent(in) :: method
    integer, intent(out), optional :: oldmethod
    call seterrorhandlingi(file%iosystem, method, oldmethod)
  end subroutine seterrorhandlingf

!>
!! @ingroup PIO_seterrorhandling 
!! @public
!! @brief set the pio error handling method for the iosystem
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param method :
!! @copydoc PIO_error_method
!<
  subroutine seterrorhandlingi(ios, method, oldmethod)
    type(iosystem_desc_t), intent(inout) :: ios
    integer, intent(in) :: method
    integer, intent(out), optional :: oldmethod

    interface
       integer(c_int) function PIOc_Set_IOSystem_Error_Handling(ios, method) &
            bind(C,name="PIOc_Set_IOSystem_Error_Handling")
         use iso_c_binding
         integer(c_int), value :: ios
         integer(c_int), value :: method
       end function PIOc_Set_IOSystem_Error_Handling
    end interface
    integer(c_int) ::  loldmethod

    loldmethod = PIOc_Set_IOSystem_Error_Handling(ios%iosysid, method)
    if(present(oldmethod)) oldmethod = loldmethod


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
!<
  subroutine PIO_initdecomp_bc(iosystem,basepiotype,dims,compstart,compcount,iodesc)
    type (iosystem_desc_t), intent(inout) :: iosystem
    integer(i4), intent(in)               :: basepiotype
    integer(i4), intent(in)               :: dims(:)
    integer (kind=PIO_OFFSET_KIND)             :: compstart(:)  
    integer (kind=PIO_OFFSET_KIND)             :: compcount(:)    
    type (IO_desc_t), intent(out)         :: iodesc

    interface
       integer(C_INT) function PIOc_InitDecomp_bc(iosysid, basetype, ndims, dims, compstart, compcount, ioidp) &
            bind(C,name="PIOc_InitDecomp_bc")
         use iso_c_binding
         integer(C_INT), value :: iosysid
         integer(C_INT), value :: basetype
         integer(C_INT), value :: ndims
         integer(C_INT) :: dims(*)
         integer(C_INT) :: ioidp
         integer(C_SIZE_T) :: compstart(*)
         integer(C_SIZE_T) :: compcount(*)
       end function PIOc_InitDecomp_bc
    end interface
    integer :: i, ndims
    integer, allocatable ::  cdims(:)
    integer(PIO_Offset_kind), allocatable :: cstart(:), ccount(:)
    integer :: ierr

    ndims = size(dims)

    allocate(cstart(ndims), ccount(ndims), cdims(ndims))

    do i=1,ndims
       cdims(i) = dims(ndims-i+1)
       cstart(i)  = compstart(ndims-i+1)-1
       cstart(i)  = compcount(ndims-i+1)
    end do

    ierr = PIOc_InitDecomp_bc(iosystem%iosysid, basepiotype, ndims, cdims, &
         cstart, ccount, iodesc%ioid)


    deallocate(cstart, ccount, cdims)


  end subroutine PIO_initdecomp_bc





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
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in)          :: lenblocks
    integer (PIO_OFFSET_KIND), intent(in)          :: compdof(:)   !> global degrees of freedom for computational decomposition
    integer (PIO_OFFSET_KIND), intent(in)          :: iodofr(:)     !> global degrees of freedom for io decomposition 
    integer (PIO_OFFSET_KIND), intent(in)          :: iodofw(:)     !> global degrees of freedom for io decomposition 
    type (io_desc_t), intent(inout)     :: iodesc




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
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (i4), intent(in) :: lenblocks
    integer (PIO_OFFSET_KIND), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (PIO_OFFSET_KIND), intent(in)          :: iodof(:)     ! global degrees of freedom for io decomposition 
    type (io_desc_t), intent(inout)     :: iodesc
    integer :: piotype
    integer(PIO_OFFSET_KIND), intent(in) :: start(:), count(:)




    if(any(iodof/=compdof)) then
       call piodie( __PIO_FILE__,__LINE__, &
            'Not sure what to do here')
    else
       call PIO_initdecomp_dof_i8(iosystem,basepiotype,dims,compdof, iodesc,PIO_REARR_SUBSET, start,count)
    endif


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
  subroutine PIO_initdecomp_dof_i4(iosystem,basepiotype,dims,compdof, iodesc, rearr, iostart, iocount)
    type (iosystem_desc_t), intent(inout) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer, optional, target :: rearr
    integer (PIO_OFFSET_KIND), optional :: iostart(:), iocount(:)
    type (io_desc_t), intent(inout)     :: iodesc
    integer(PIO_OFFSET_KIND), pointer :: internal_compdof(:)
    integer(i4), intent(in)           :: dims(:)

    allocate(internal_compdof(size(compdof)))
    internal_compdof = int(compdof,PIO_OFFSET_KIND)
    
    if(present(iostart) .and. present(iocount) ) then
       call pio_initdecomp_dof_i8(iosystem, basepiotype, dims, internal_compdof, iodesc, &
            PIO_REARR_SUBSET, iostart, iocount)
    else 
       call pio_initdecomp_dof_i8(iosystem, basepiotype, dims, internal_compdof, iodesc, rearr)
    endif
    deallocate(internal_compdof)

  end subroutine PIO_initdecomp_dof_i4


  subroutine PIO_initdecomp_internal(iosystem,basepiotype,dims,maplen, compdof, iodesc, rearr, iostart, iocount)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer, intent(in) :: maplen
    integer (PIO_OFFSET_KIND), intent(in) :: compdof(maplen)   ! global degrees of freedom for computational decomposition
    integer, optional, target :: rearr
    integer (PIO_OFFSET_KIND), optional :: iostart(:), iocount(:)
    type (io_desc_t), intent(inout)     :: iodesc

    integer(c_int) :: ndims
    integer(c_int), dimension(:), allocatable, target :: cdims
    integer(PIO_OFFSET_KIND), dimension(:), allocatable, target :: cstart, ccount

    type(C_PTR) :: crearr
    interface
       integer(C_INT) function PIOc_InitDecomp(iosysid,basetype,ndims,dims, &
            maplen, compmap, ioidp, rearr, iostart, iocount)  &
            bind(C,name="PIOc_InitDecomp")
         use iso_c_binding
         integer(C_INT), value :: iosysid
         integer(C_INT), value :: basetype
         integer(C_INT), value :: ndims
         integer(C_INT) :: dims(*)
         integer(C_INT), value :: maplen
         integer(C_SIZE_T) :: compmap(*)
         integer(C_INT) :: ioidp
         type(C_PTR), value :: rearr
         type(C_PTR), value :: iostart
         type(C_PTR), value :: iocount
       end function PIOc_InitDecomp
    end interface
    integer :: ierr,i 
    
    ndims = size(dims)
    allocate(cdims(ndims))
    do i=1,ndims
       cdims(i) = dims(ndims-i+1)
    end do

    if(present(rearr)) then
       crearr = C_LOC(rearr)
    else
       crearr = C_NULL_PTR
    endif

    if(present(iostart) .and. present(iocount)) then
       allocate(cstart(ndims), ccount(ndims))
       do i=1,ndims
          cstart(i) = iostart(ndims-i+1)-1
          ccount(i) = iocount(ndims-i+1)
       end do

       ierr = PIOc_InitDecomp(iosystem%iosysid, basepiotype, ndims, cdims, &
            maplen, compdof, iodesc%ioid, crearr, C_LOC(cstart), C_LOC(ccount))
       deallocate(cstart, ccount)
    else
        ierr = PIOc_InitDecomp(iosystem%iosysid, basepiotype, ndims, cdims, &
            maplen, compdof, iodesc%ioid, crearr, C_NULL_PTR, C_NULL_PTR)
    end if

    deallocate(cdims)


  end subroutine PIO_initdecomp_internal


  subroutine PIO_initdecomp_dof_i8(iosystem,basepiotype,dims,compdof, iodesc, rearr, iostart, iocount)
    type (iosystem_desc_t), intent(in) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (PIO_OFFSET_KIND), intent(in) :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer, optional, target :: rearr
    integer (PIO_OFFSET_KIND), optional :: iostart(:), iocount(:)
    type (io_desc_t), intent(inout)     :: iodesc
    integer :: maplen

#ifdef TIMING
    call t_startf("PIO:initdecomp_dof")
#endif

    maplen = size(compdof)

    call PIO_initdecomp_internal(iosystem, basepiotype, dims, maplen, compdof, iodesc, rearr, iostart,iocount)


#ifdef TIMING
    call t_stopf("PIO:initdecomp_dof")
#endif

  end subroutine PIO_initdecomp_dof_i8

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
    use pio_types, only : pio_internal_error
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
       integer(c_int) function PIOc_Init_Intracomm_from_F90(f90_comp_comm, num_iotasks, stride,base,rearr,iosysidp) &
            bind(C,name="PIOc_Init_Intracomm_from_F90")
         use iso_c_binding
         integer(C_INT), value :: f90_comp_comm
         integer(C_INT), value :: num_iotasks
         integer(C_INT), value :: stride
         integer(C_INT), value :: base
         integer(C_INT), value :: rearr
         integer(C_INT) :: iosysidp
       end function PIOc_Init_Intracomm_from_F90
    end interface

#ifdef TIMING
    call t_startf("PIO:init")
#endif
    lbase=0
    if(present(base)) lbase=base
    ierr = PIOc_Init_Intracomm_from_F90(comp_comm,num_iotasks,stride,lbase,rearr,iosystem%iosysid)

    call CheckMPIReturn("Bad Initialization in PIO_Init_Intracomm:  ", ierr,__FILE__,__LINE__)
#ifdef TIMING
    call t_stopf("PIO:init")
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
    call t_startf("PIO:init")
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
       call mpi_info_create(iosystem(i)%info,ierr)
       ! turn on mpi-io aggregation 
       !DBG    print *,'PIO_init: before call to setnumagg'
!       itmp = num_aggregator
!       call mpi_bcast(itmp, 1, mpi_integer, 0, iosystem%union_comm, ierr)
!       if(itmp .gt. 0) then 
!          write(cb_nodes,('(i5)')) itmp
!#ifdef BGQ
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
    end do

    if(DebugAsync) print*,__PIO_FILE__,__LINE__, iosystem(1)%ioranks


    iosystem%num_aiotasks = iosystem%num_iotasks
    iosystem%numost = PIO_NUM_OST

    ! This routine does not return
    if(io_comm /= MPI_COMM_NULL) call pio_msg_handler(component_count,iosystem) 
    
    if(DebugAsync) print*,__PIO_FILE__,__LINE__, iosystem(1)%ioranks
#ifdef TIMING
    call t_stopf("PIO:init")
#endif
#endif
#endif
  end subroutine init_intercom


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
    if(iosystem%iosysid /= -1) then
       ierr = PIOc_finalize(iosystem%iosysid)
    endif
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

   logical function pio_iotype_available( iotype) result(available)
     integer, intent(in) :: iotype
     interface
        integer(C_INT) function PIOc_iotype_available(iotype) &
             bind(C,name="PIOc_iotype_available")
          use iso_c_binding
          integer(C_INT), intent(in), value :: iotype
        end function PIOc_iotype_available
     end interface
     available= (PIOc_iotype_available(iotype) == 1)

   end function pio_iotype_available


!> 
!! @public
!! @ingroup PIO_createfile 
!! @brief  Create a NetCDF or PNetCDF file using PIO.
!! @details  Input parameters are read on comp task 0 and ignored elsewhere
!! @param iosystem : A defined pio system descriptor created by a call to @ref PIO_init (see PIO_types)
!! @param file	:  The returned file descriptor
!! @param iotype : @copydoc PIO_iotype
!! @param fname : The name of the file to open
!! @param amode_in : The NetCDF creation mode flag. the following flags are available: 
!! (1) zero value or NC_NOWRITE is default and opens the file with read-only access. 
!! (2) NC_WRITE for read-write access. 
!! (3) NC_SHARE is used for NetCDF classic, and dangerous with this application. 
!! (4) NC_WRITE|NC_SHARE
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
         character(kind=c_char) :: fname(*)
         integer(c_int), value :: mode
       end function PIOc_createfile
    end interface
    character, allocatable :: cfname(:)
    integer :: i, nl
#ifdef TIMING
    call t_startf("PIO:createfile")
#endif
    mode = 0
    if(present(amode_in)) mode = amode_in
    nl = len_trim(fname)
    allocate(cfname(nl+1))
    do i=1,nl
       cfname(i) = fname(i:i)
    enddo
    cfname(nl+1)=C_NULL_CHAR
    ierr = PIOc_createfile(iosystem%iosysid, file%fh, iotype, cfname, mode)
    deallocate(cfname)
    file%iosystem => iosystem
#ifdef TIMING
    call t_stopf("PIO:createfile")
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
  integer function PIO_openfile(iosystem, file, iotype, fname,mode) result(ierr)

!    use ifcore, only: tracebackqq
    type (iosystem_desc_t), intent(inout), target :: iosystem
    type (file_desc_t), intent(out) :: file
    integer, intent(in) :: iotype
    character(len=*), intent(in)  :: fname
    integer, optional, intent(in) :: mode
    integer :: iorank
    interface
       integer(C_INT) function PIOc_openfile(iosysid, fh, iotype, fname,mode) &
         bind(C,NAME='PIOc_openfile')
         use iso_c_binding
         implicit none
         integer(c_int), value :: iosysid
         integer(c_int) :: fh
         integer(c_int) :: iotype
         character(kind=c_char) :: fname(*)
         integer(c_int), value :: mode
       end function PIOc_openfile
    end interface
    integer :: imode=0, i, nl
    character, allocatable :: cfname(:)
#ifdef TIMING
    call t_startf("PIO:openfile")
#endif
    if(present(mode)) imode = mode
    nl = len_trim(fname)
    allocate(cfname(nl+1))
    do i=1,nl
       cfname(i) = fname(i:i)
    enddo
    cfname(nl+1)=C_NULL_CHAR
    ierr = PIOc_openfile( iosystem%iosysid, file%fh, iotype, cfname, imode)
    deallocate(cfname)
    file%iosystem => iosystem

#ifdef TIMING
    call t_stopf("PIO:openfile")
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
  
    call syncfile(file)     

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
    call t_startf("PIO:closefile")
#endif
    ierr = PIOc_closefile(file%fh)
    nullify(file%iosystem)
#ifdef TIMING
    call t_stopf("PIO:closefile")
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

!>
!! @public 
!! @ingroup PIO_deletefile
!! @brief Delete a file 
!! @details
!! @param ios : a pio system handle
!! @param fname : a filename
!<
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
