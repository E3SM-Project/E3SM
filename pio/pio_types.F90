#define __PIO_FILE__ "pio_types.F90"
module pio_types
    use pio_kinds

#ifdef _NETCDF
     use netcdf                                  ! _EXTERNAL
#endif
    implicit none
    private 

    !-------------------------------------------
    !  data structure to describe decomposition
    !-------------------------------------------
    type, public :: DecompMap_t
#ifdef SEQUENCE
	sequence
#endif
        integer(i4) :: start
        integer(i4) :: length
    end type

    !------------------------------------
    !  a file descriptor data structure
    !------------------------------------
!>
!! @public
!! @struct iosystem_desc_t iosystem_desc_t
!! @brief A defined PIO system descriptor created by @ref PIO_init (see pio_types)
!<
    type, public :: IOSystem_desc_t
#ifdef SEQUENCE
	sequence
#endif
        
        integer(i4)              :: union_comm         ! The intracomm union of comp and io communicators (for async only)
        integer(i4)              :: IO_comm            ! The IO communicator
        integer(i4)              :: comp_comm          ! The Compute communicator
        integer(i4)              :: intercomm          ! the intercomm (may be MPI_COMM_NULL)
        
        integer(i4)              :: my_comm            ! either comp_comm or intercomm
        integer(i4)              :: num_tasks          !  number of tasks
        integer(i4)              :: num_iotasks        ! total number of IO tasks
        integer(i4)              :: num_aiotasks       ! number of actual IO tasks
        integer(i4)              :: num_comptasks

        integer(i4)              :: union_rank
        integer(i4)              :: comp_rank          ! the computational rank
        integer(i4)              :: io_rank            ! the io rank if io_rank = -1 not an IO processor
!
        integer(i4)              :: Info               ! MPI-IO info structure
        integer(i4)              :: numOST
        
! rank of the io and comp roots in the intercomm
        integer(i4)              :: IOMaster           ! The intercom of the io_rank 0
        integer(i4)              :: compMaster           ! The intercom of the comp_rank 0

! rank of the io and comp roots in the union_comm
        integer(i4)              :: IOroot           ! The union_rank of the io_rank 0
        integer(i4)              :: comproot           ! The union_rank of the comp_rank 0

        logical(log_kind)        :: IOproc             ! .true. if an IO processor
        logical(log_kind)        :: UseRearranger      ! .true. if data rearrangement is necessary
        logical(log_kind)        :: async_interface    ! .true. if using the async interface model
        integer(i4)              :: rearr         ! type of rearranger
                                                  ! e.g. rearr_{none,box}
	integer(i4)              :: error_handling ! how pio handles errors
        integer(i4),pointer      :: ioranks(:)         ! the computational ranks for the IO tasks

	! This holds the IODESC
    end type

    type iosystem_list_t
       type(iosystem_desc_t), pointer :: this_iosystem
    end type iosystem_list_t


    integer, parameter :: MAX_IO_SYSTEMS=6
    type(iosystem_list_t) :: iosystems(MAX_IO_SYSTEMS)


     
!> 
!! @public
!! @struct file_desc_t file_desc_t
!! @brief File descriptor returned by \ref PIO_openfile or \ref PIO_createfile (see pio_types)
!>
    type, public :: File_desc_t
       type(iosystem_desc_t), pointer :: iosystem
       integer(i4) :: fh
       integer(kind=PIO_OFFSET) :: offset             ! offset into file
       integer(i4)              :: iotype             ! Type of IO to perform see parameter statement below     
       logical                  :: file_is_open = .false.
    end type File_desc_t


    !------------------------------------------------------
    !  data structure to describe a data movement operator
    !------------------------------------------------------
    type, public :: IO_desc2_t
#ifdef SEQUENCE
	sequence
#endif
        integer(i4)         ::  fileTYPE   ! MPI data types for file
        integer(i4)         ::  elemTYPE
        integer(i4)         :: n_elemTYPE
        integer(i4)         :: n_words
    end type IO_desc2_t

!>
!! @private
!! @defgroup iodesc_generate io descriptors, generating
!! @brief The io descriptor structure in defined in this subroutine 
!! and subsequently used in @ref PIO_read_darray, @ref PIO_write_darray,
!! @ref PIO_put_var, @ref PIO_get_var calls (see pio_types).
!<

!>
!! @public 
!! @struct io_desc_t io_desc_t
!! @brief  An io descriptor handle that is generated in @ref PIO_initdecomp 
!! (see pio_types)
!<
   type, public :: IO_desc_t
#ifdef SEQUENCE
	sequence
#endif
        type(IO_desc2_t)    :: Read
        type(IO_desc2_t)    :: Write
	integer(kind=PIO_Offset), pointer :: start(:) => NULL()
	integer(kind=PIO_Offset), pointer :: count(:) => NULL()

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! fields for box-based rearranger
        ! should put this in its own derived type later

        integer :: baseTYPE

        integer, pointer :: dest_ioproc(:)=> NULL()     ! for each dof
        integer, pointer :: dest_ioindex(:)=> NULL()    ! for each dof


        ! Values needed only on io procs
        integer,pointer :: rfrom(:)=> NULL()   ! rfrom(nrecvs)= rank of ith sender
        integer,pointer :: rtype(:)=> NULL()   ! rtype(nrecvs)=mpi types for receives

        
        ! needed on all procs
        integer,pointer :: scount(:)=> NULL()  ! scount(num_iotasks)= # sends to ith ioproc
        integer,pointer :: stype(:)=> NULL()   ! stype(num_iotasks)=mpi type for sends

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(i4) :: async_id


        type (DecompMap_t)  :: IOmap      ! IO decomposition map
        type (DecompMap_t)  :: COMPmap    ! Computational decomposition map
        integer :: nrecvs                      ! valid for io procs
        integer(i4)         :: glen       ! global length of array in words
	integer(i4)         :: compsize   ! size of expected comp buffer
        integer(i4)         :: maxiobuflen   ! size of largest iobuffer
    end type

!>
!! @public
!! @struct var_desc_t var_desc_t
!! @brief A variable descriptor returned from @ref PIO_def_var (see pio_types) 
!<
    type, public :: Var_desc_t
#ifdef SEQUENCE
	sequence
#endif	
        integer(i4)     :: varID
        integer(i4)     :: rec   ! This is a record number or pointer into the unlim dimension of the	    
                                 ! netcdf file
	integer(i4)     :: type
        integer(i4)     :: ndims ! number of dimensions as defined on the netcdf file.
    end type 

!>
!! @defgroup PIO_iotype PIO_iotype
!! @public
!! @brief An integer parameter which controls the iotype
!! @details
!!   - PIO_iotype_pbinary : Use MPI-IO to read/write C like binary file
!!   - PIO_iotype_direct_pbinary: Use MPI-IO to read/write direct access binary files
!!   - PIO_iotype_binary : serial read/write of binary files using 'base_node'
!!   - PIO_iotype_pnetcdf : parallel read/write of pNetCDF files (netcdf3)
!!   - PIO_iotype_netcdf : serial read/write of NetCDF files using 'base_node' (netcdf3)
!!   - PIO_iotype_netcdf4c : parallel read/serial write of NetCDF4 (HDF5) files with data compression
!!   - PIO_iotype_netcdf4p : parallel read/write of NETCDF4 (HDF5) files 
!>
    integer(i4), public, parameter ::  &
        PIO_iotype_pbinary = 1, &! use MPI-IO with data types to read/write C like binary files
        PIO_iotype_direct_pbinary = 2,  & !use MPI-IO with data types to read/write direct access binary files
        PIO_iotype_binary  = 4, &   ! serial read/write of binary files using 'base_node'
        PIO_iotype_pnetcdf = 5, &   ! parallel read/write of pNetCDF files
        PIO_iotype_netcdf  = 6, &   ! serial read/write of NetCDF file using 'base_node'
        PIO_iotype_netcdf4c = 7, &  ! netcdf4 (hdf5 format) file opened for compression (serial write access only)   
        PIO_iotype_netcdf4p = 8     ! netcdf4 (hdf5 format) file opened in parallel (all netcdf4 files for read will be opened this way)

    integer(i4), public, parameter ::                       &
        iotype_pbinary = PIO_iotype_pbinary,                &
        iotype_direct_pbinary = PIO_iotype_direct_pbinary,  &
        iotype_binary  = PIO_iotype_binary,                 &
        iotype_pnetcdf = PIO_iotype_pnetcdf,                &
        iotype_netcdf  = PIO_iotype_netcdf


!>
!! @defgroup PIO_rearr_method PIO_rearr_method
!! @public 
!! @brief The three choices to control rearrangement are:
!! @details
!!  - PIO_rearr_none : Do not use any form of rearrangement
!!  - PIO_rearr_box : Use a PIO internal box rearrangement
!>
    integer(i4), public, parameter :: PIO_rearr_none = 0
    integer(i4), public, parameter :: PIO_rearr_box =  1

!> 
!! @public
!! @defgroup PIO_error_method error handling method 
!! @brief  The three types of error handling methods are: 
!! @details
!!  - PIO_INTERNAL_ERROR  : abort on error from any task
!!  - PIO_BCAST_ERROR     : broadcast an error from io_rank 0 to all tasks in comm
!!  - PIO_RETURN_ERROR    : do nothing - allow the user to handle it
!<
  integer(i4), public :: PIO_INTERNAL_ERROR = -51
  integer(i4), public :: PIO_BCAST_ERROR = -52
  integer(i4), public :: PIO_RETURN_ERROR = -53

!>
!! @public 
!! @defgroup error_return error return codes
!! @brief : The error return code; ierr != PIO_noerr indicates
!! an error. (see @ref PIO_seterrorhandling ) 
!> 

!>
!! @struct use_PIO_kinds
!! @brief The type of variable(s) associated with this iodesc.
!! @copydoc PIO_kinds
!<

!>
!! @public 
!! @defgroup PIO_kinds PIO_kinds
!! @brief The base types supported by PIO are:
!! @details
!!  - PIO_double : 8-byte reals or double precision 
!!  - PIO_real : 4-byte reals
!!  - PIO_int :  4-byte integers
!!  - PIO_char : character
!<
#ifdef _PNETCDF
#include <pnetcdf.inc>   /* _EXTERNAL */
   integer, public, parameter :: PIO_global = nf_global
   integer, public, parameter :: PIO_unlimited = nf_unlimited
   integer, public, parameter :: PIO_double = nf_double
   integer, public, parameter :: PIO_real   = nf_real
   integer, public, parameter :: PIO_int    = nf_int
   integer, public, parameter :: PIO_char   = nf_char
   integer, public, parameter :: PIO_noerr  = nf_noerr
   integer, public, parameter :: PIO_WRITE  = nf_write
   integer, public, parameter :: PIO_nowrite  = nf_nowrite
   integer, public, parameter :: PIO_CLOBBER = nf_clobber	
   integer, public, parameter :: PIO_NOCLOBBER = nf_NOclobber	
   integer, public, parameter :: PIO_NOFILL = nf_nofill
   integer, public, parameter :: PIO_MAX_NAME = nf_max_name
   integer, public, parameter :: PIO_MAX_VAR_DIMS = nf_max_var_dims
   integer, public, parameter :: PIO_64BIT_OFFSET = nf_64bit_offset
   integer, public, parameter :: PIO_num_OST =  16
#else
#ifdef _NETCDF
   integer, public, parameter :: PIO_global = nf90_global
   integer, public, parameter :: PIO_unlimited = nf90_unlimited
   integer, public, parameter :: PIO_double = nf90_double
   integer, public, parameter :: PIO_real   = nf90_real
   integer, public, parameter :: PIO_int    = nf90_int
   integer, public, parameter :: PIO_char   = nf90_char
   integer, public, parameter :: PIO_noerr  = nf90_noerr
   integer, public, parameter :: PIO_WRITE  = nf90_write
   integer, public, parameter :: PIO_nowrite = nf90_nowrite
   integer, public, parameter :: PIO_CLOBBER = nf90_clobber	
   integer, public, parameter :: PIO_NOCLOBBER = nf90_NOclobber	
   integer, public, parameter :: PIO_NOFILL = nf90_nofill
   integer, public, parameter :: PIO_MAX_NAME = nf90_max_name
   integer, public, parameter :: PIO_MAX_VAR_DIMS = nf90_max_var_dims
   integer, public, parameter :: PIO_64BIT_OFFSET = nf90_64bit_offset
   integer, public, parameter :: PIO_num_OST =  16
#else
   integer, public, parameter :: PIO_global = 0
   integer, public, parameter :: PIO_double = 6
   integer, public, parameter :: PIO_real   = 5
   integer, public, parameter :: PIO_int    = 4
   integer, public, parameter :: PIO_char   = 2
   integer, public, parameter :: PIO_noerr  = 0
   integer, public, parameter :: PIO_MAX_NAME = 25
   integer, public, parameter :: PIO_MAX_VAR_DIMS = 6
   integer, public, parameter :: PIO_CLOBBER = 10
   integer, public, parameter :: PIO_NOCLOBBER = 11
   integer, public, parameter :: PIO_WRITE = 20
   integer, public, parameter :: PIO_NOWRITE = 21
   integer, public, parameter :: PIO_64BIT_OFFSET = 0
   integer, public, parameter :: PIO_num_OST =  16
#endif
#endif
! should be defined in the mct mpiserial library.
#ifdef _MPISERIAL
   integer, public, parameter :: MPI_DATATYPE_NULL=0
#endif

end module pio_types
