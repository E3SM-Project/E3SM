#define _FILE_ "pio_types.F90"
module pio_types
    use pio_kinds
#ifdef _USEMCT
! needed for MCT
    use m_GlobalSegMap, only: GlobalSegMap      ! _EXTERNAL
    use m_Rearranger, only: Rearranger          ! _EXTERNAL
#endif

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
    type, public :: IOSystem_desc_t
#ifdef SEQUENCE
	sequence
#endif
        integer(i4)              :: IO_comm            ! The IO communicator
        integer(i4)              :: comp_comm          ! The Compute communicator
        integer(i4)              :: num_tasks        ! total number of tasks
        integer(i4)              :: num_iotasks        ! total number of IO tasks
        integer(i4)              :: io_stride          ! stride between IO tasks
        integer(i4)              :: comp_rank          ! the computational rank
        integer(i4)              :: io_rank            ! the io rank if io_rank = -1 not an IO processor
!
        integer(i4)              :: Info               ! MPI-IO info structure
        integer(i4)              :: IOMaster           ! The comp_rank of the io_rank 0
        logical(log_kind)        :: IOproc             ! .true. if an IO processor
        logical(log_kind)        :: UseRearranger      ! .true. if data rearrangement is necessary
        integer(i4)              :: rearr         ! type of rearranger
                                                  ! e.g. rearr_{none,mct,box}
	integer(i4)              :: error_handling ! how pio handles errors

	! This holds the IODESC
    end type
    
    type, public :: File_desc_t
       type(iosystem_desc_t), pointer :: iosystem
       integer(i4) :: fh
       integer(kind=PIO_OFFSET) :: offset             ! offset into file
       integer(i4)              :: iotype             ! Type of IO to perform see parameter statement below    
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


    !------------------------------------------------------
    ! IO_desc_t
    !------------------------------------------------------


    type, public :: IO_desc_t
#ifdef SEQUENCE
	sequence
#endif
        type(IO_desc2_t)    :: Read
        type(IO_desc2_t)    :: Write
	integer(kind=PIO_Offset), pointer :: start(:)
	integer(kind=PIO_Offset), pointer :: count(:)
#ifdef _USEMCT
        ! MCT GlobalSegMaps defined on comp_comm

	integer(i4) :: lsize_comp      ! local size of GSMap for comp layout
        integer(i4) :: lsize_io        ! local size of GSMap for IO layout

        type (Rearranger) :: rearr_comp_to_io   ! mct rearranger comp->io
        type (Rearranger) :: rearr_io_to_comp   ! mct rearranger io->comp

        integer(i4), pointer :: compDOF_index(:)  ! permutation array
#endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! fields for box-based rearranger
        ! should put this in its own derived type later

        integer :: baseTYPE

        integer, pointer :: dest_ioproc(:)     ! for each dof
        integer, pointer :: dest_ioindex(:)    ! for each dof


        ! Values needed only on io procs
        integer,pointer :: rfrom(:)   ! rfrom(nrecvs)= rank of ith sender
        integer,pointer :: rtype(:)   ! rtype(nrecvs)=mpi types for receives

        
        ! needed on all procs
        integer,pointer :: scount(:)  ! scount(num_iotasks)= # sends to ith ioproc
        integer,pointer :: stype(:)   ! stype(num_iotasks)=mpi type for sends

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        


        type (DecompMap_t)  :: IOmap      ! IO decomposition map
        type (DecompMap_t)  :: COMPmap    ! Computational decomposition map
        integer :: nrecvs                      ! valid for io procs
        integer(i4)         :: glen       ! global length of array in words
	integer(i4)         :: compsize   ! size of expected comp buffer
        integer(i4)         :: maxiobuflen   ! size of largest iobuffer
    end type

    type, public :: Var_desc_t
#ifdef SEQUENCE
	sequence
#endif	
        integer(i4)     :: varID
        integer(i4)     :: rec   ! This is a record number or pointer into the unlim dimension of the	    
	integer(i4)     :: type
                                 ! netcdf file
    end type 

! !PUBLIC parameters:
    integer(i4), public, parameter ::  &
        iotype_pbinary = 1, &! use MPI-IO with data types to read/write C like binary files
        iotype_direct_pbinary = 2,  & !use MPI-IO with data types to read/write direct access binary files
        iotype_binary  = 4, &   ! serial read/write of binary files using 'base_node'
        iotype_pnetcdf = 5, &   ! parallel read/write of pNetCDF files
        iotype_netcdf  = 6      ! serial read/write of NetCDF file using 'base_node'


! Constants to define rearranger type

    integer(i4), public, parameter :: PIO_rearr_none = 0
    integer(i4), public, parameter :: PIO_rearr_mct =  1
    integer(i4), public, parameter :: PIO_rearr_box =  2

!> constants for error handling methods
!!  Three choices for error handling:
!!  1: abort on error from any task           PIO_INTERNAL_ERROR
!!  2: broadcast an error from io_rank 0      PIO_BCAST_ERROR
!!  3: do nothing - allow the user to handle it PIO_RETURN_ERROR
!<
  integer(i4), public :: PIO_INTERNAL_ERROR = -51
  integer(i4), public :: PIO_BCAST_ERROR = -52
  integer(i4), public :: PIO_RETURN_ERROR = -53




! pnetcdf/netcdf constants should be the same, but
! have to pick one...


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
   integer, public, parameter :: PIO_64bit_offset = 0
#endif
#endif
! should be defined in the mct mpiserial library.
#ifdef _MPISERIAL
   integer, public, parameter :: MPI_DATATYPE_NULL=0
#endif

   public:: pio_type_to_mpi_type


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! pio_type_to_mpi_type


  integer function pio_type_to_mpi_type(ptype) result(mtype)
    use pio_support
    implicit none
    include 'mpif.h'            ! _EXTERNAL
    integer, intent(in):: ptype

    select case(ptype)
      case (PIO_double)
         mtype=MPI_REAL8
      case (PIO_real)
         mtype=MPI_REAL4
      case (PIO_int)
         mtype=MPI_INTEGER
      case (PIO_char)
         mtype=MPI_CHARACTER
      case default
         call piodie( _FILE_,__LINE__, &
                      'Could not convert pio type=',ptype,' to an mpi type')
    end select

  end function pio_type_to_mpi_type



end module pio_types
