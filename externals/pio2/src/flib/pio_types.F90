!>
!! @file
!! @brief Derived datatypes and constants for PIO Fortran API
!!
!<
module pio_types
    use pio_kinds
    use iso_c_binding
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
!! @struct iosystem_desc_t
!! @brief A defined PIO system descriptor created by @ref PIO_init (see pio_types)
!<
    type, public :: IOSystem_desc_t
        integer(kind=c_int) :: iosysid = -1
     end type IOSystem_desc_t

!>
!! @private
!! @struct io_data_list
!! @brief Linked list of buffers for pnetcdf non-blocking interface
!>
!    type, public :: io_data_list
!       integer :: request
!       real(r4), pointer :: data_real(:) => null()
!       integer(i4), pointer :: data_int(:) => null()
!       real(r8), pointer :: data_double(:) => null()
!       type(io_data_list), pointer :: next => null()
!    end type io_data_list


!>
!! @public
!! @struct file_desc_t
!! @brief File descriptor returned by \ref PIO_openfile or \ref PIO_createfile (see pio_types)
!!
!>
    type, public :: File_desc_t
       integer(kind=c_int) :: fh
       type(iosystem_desc_t), pointer :: iosystem => null()
    end type File_desc_t


!>
!! @private
!! @defgroup iodesc_generate io descriptors, generating
!! @brief The io descriptor structure in defined in this subroutine
!! and subsequently used in @ref PIO_read_darray, @ref PIO_write_darray,
!! @ref PIO_put_var, @ref PIO_get_var calls (see pio_types).
!<

!>
!! @public
!! @struct io_desc_t
!! @brief  An io descriptor handle that is generated in @ref PIO_initdecomp
!! (see pio_types)
!<
   type, public :: io_desc_t
#ifdef SEQUENCE
	sequence
#endif
        integer(i4)         :: ioid
    end type

!>
!! @public
!! @struct var_desc_t
!! @brief A variable descriptor returned from @ref PIO_def_var (see pio_types)
!<
    type, public :: Var_desc_t
#ifdef SEQUENCE
       sequence
#endif
       integer(i4) :: varID
       integer(i4) :: ncid
    end type Var_desc_t

!>
!! @defgroup PIO_iotype PIO_iotype
!! @public
!! @brief An integer parameter which controls the iotype
!! @details
!!   - PIO_iotype_pnetcdf : parallel read/write of pNetCDF files (netcdf3)
!!   - PIO_iotype_netcdf : serial read/write of NetCDF files using 'base_node' (netcdf3)
!!   - PIO_iotype_netcdf4c : parallel read/serial write of NetCDF4 (HDF5) files with data compression
!!   - PIO_iotype_netcdf4p : parallel read/write of NETCDF4 (HDF5) files
!>
    integer(i4), public, parameter ::  &
        PIO_iotype_pnetcdf = 1, &   ! parallel read/write of pNetCDF files
        PIO_iotype_netcdf  = 2, &   ! serial read/write of NetCDF file using 'base_node'
        PIO_iotype_netcdf4c = 3, &  ! netcdf4 (hdf5 format) file opened for compression (serial write access only)
        PIO_iotype_netcdf4p = 4    ! netcdf4 (hdf5 format) file opened in parallel (all netcdf4 files for read will be opened this way)


! These are for backward compatability and should not be used or expanded upon
    integer(i4), public, parameter ::                       &
        iotype_pnetcdf = PIO_iotype_pnetcdf,                &
        iotype_netcdf  = PIO_iotype_netcdf


!>
!! @defgroup PIO_rearr_method PIO_rearr_method
!! @public
!! @brief The three choices to control rearrangement are:
!! @details
!!  - PIO_rearr_none : Do not use any form of rearrangement
!!  - PIO_rearr_box : Use a PIO internal box rearrangement
!! -  PIO_rearr_subset : Use a PIO internal subsetting rearrangement
!>

    integer(i4), public, parameter :: PIO_rearr_box =  1
    integer(i4), public, parameter :: PIO_rearr_subset =  2

!>
!! @public
!! @defgroup PIO_error_method error_methods
!! @details
!! The three types of error handling methods are:
!!  - PIO_INTERNAL_ERROR  : abort on error from any task
!!  - PIO_BCAST_ERROR     : broadcast an error from io_rank 0 to all tasks in comm
!!  - PIO_RETURN_ERROR    : do nothing - allow the user to handle it
!<
  integer(i4), public, parameter :: PIO_INTERNAL_ERROR = -51
  integer(i4), public, parameter :: PIO_BCAST_ERROR = -52
  integer(i4), public, parameter :: PIO_RETURN_ERROR = -53

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
   integer, public, parameter :: PIO_64BIT_DATA = nf_64bit_data
   integer, public, parameter :: PIO_FILL_INT = nf_fill_int;
   real, public, parameter :: PIO_FILL_FLOAT = nf_fill_float;
   double precision, public, parameter :: PIO_FILL_DOUBLE = nf_fill_double;

#else
#ifdef _NETCDF
#include <netcdf.inc>   /* _EXTERNAL */
   integer, public, parameter :: PIO_global = nf_global
   integer, public, parameter :: PIO_unlimited = nf_unlimited
   integer, public, parameter :: PIO_double = nf_double
   integer, public, parameter :: PIO_real   = nf_real
   integer, public, parameter :: PIO_int    = nf_int
   integer, public, parameter :: PIO_char   = nf_char
   integer, public, parameter :: PIO_noerr  = nf_noerr
   integer, public, parameter :: PIO_WRITE  = nf_write
   integer, public, parameter :: PIO_nowrite = nf_nowrite
   integer, public, parameter :: PIO_CLOBBER = nf_clobber
   integer, public, parameter :: PIO_NOCLOBBER = nf_NOclobber
   integer, public, parameter :: PIO_NOFILL = nf_nofill
   integer, public, parameter :: PIO_MAX_NAME = nf_max_name
   integer, public, parameter :: PIO_MAX_VAR_DIMS = nf_max_var_dims
   integer, public, parameter :: PIO_64BIT_OFFSET = nf_64bit_offset
   integer, public, parameter :: PIO_64BIT_DATA = 0
   integer, public, parameter :: PIO_FILL_INT = nf_fill_int;
   real, public, parameter :: PIO_FILL_FLOAT = nf_fill_float;
   double precision, public, parameter :: PIO_FILL_DOUBLE = nf_fill_double;
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
   integer, public, parameter :: PIO_64BIT_DATA = 0
   integer, public, parameter :: PIO_FILL_INT = -2147483647;
   real, public, parameter :: PIO_FILL_FLOAT =  9.9692099683868690e+36;
   double precision, public, parameter :: PIO_FILL_DOUBLE = 9.9692099683868690e+36;
#endif
#endif
   integer, public, parameter :: PIO_num_OST =  16

   type, public :: PIO_rearr_comm_fc_opt_t
      logical :: enable_hs            ! Enable handshake?
      logical :: enable_isend         ! Enable isends?
      integer :: max_pend_req         ! Maximum pending requests
    end type PIO_rearr_comm_fc_opt_t

    type, public :: PIO_rearr_opt_t
      integer                         :: comm_type
      integer                         :: fcd       ! Flow control direction
      type(PIO_rearr_comm_fc_opt_t)   :: comm_fc_opts_comp2io
      type(PIO_rearr_comm_fc_opt_t)   :: comm_fc_opts_io2comp
    end type PIO_rearr_opt_t

end module pio_types
