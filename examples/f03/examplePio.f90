!> @file
!! A simple Fortran example for the ParallelIO Library.
module pioExample

    use pio, only : PIO_init, PIO_rearr_subset, iosystem_desc_t, file_desc_t
    use pio, only : PIO_finalize, PIO_noerr, PIO_iotype_netcdf, PIO_createfile
    use pio, only : PIO_int,var_desc_t, PIO_redef, PIO_def_dim, PIO_def_var, PIO_enddef
    use pio, only : PIO_closefile, io_desc_t, PIO_initdecomp, PIO_write_darray
    use pio, only : PIO_freedecomp, PIO_clobber, PIO_read_darray, PIO_syncfile, PIO_OFFSET_KIND
    use pio, only : PIO_nowrite, PIO_openfile

    implicit none
    save
    private

    include 'mpif.h'

    integer, parameter :: LEN = 16  ! total length of the array we are using.  This is then divided among MPI processes
    integer, parameter :: VAL = 42  ! value used for array that will be written to netcdf file
    integer, parameter :: ERR_CODE = 99

    !> @brief A class to hold example code and data.
    !! This class contains the data and functions to execute the
    !! example.
    type, public :: pioExampleClass

        integer :: myRank
        integer :: ntasks
        integer :: niotasks
        integer :: stride        ! stride in the mpi rank between io tasks
        integer :: numAggregator
        integer :: optBase

        type(iosystem_desc_t) :: pioIoSystem
        type(file_desc_t)     :: pioFileDesc
        type(var_desc_t)      :: pioVar
        type(io_desc_t)       :: iodescNCells

        integer               :: iotype
        integer               :: pioDimId
        integer               :: ista
        integer               :: isto
        integer               :: arrIdxPerPe
        integer, dimension(1) :: dimLen

        integer, allocatable  :: dataBuffer(:)
        integer, allocatable  :: readBuffer(:)
        integer, allocatable  :: compdof(:)

        character(len=255) :: fileName

    contains

        procedure,  public  :: init
        procedure,  public  :: createDecomp
        procedure,  public  :: createFile
        procedure,  public  :: defineVar
        procedure,  public  :: writeVar
        procedure,  public  :: readVar
        procedure,  public  :: closeFile
        procedure,  public  :: cleanUp
        procedure,  private :: errorHandle

    end type pioExampleClass

contains

    !> @brief Initialize MPI, ParallelIO, and example data.
    !! Initialize the MPI and ParallelIO libraries. Also allocate
    !! memory to write and read the sample data to the netCDF file.
    subroutine init(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: ierr,i

        !
        ! initialize MPI
        !

        call MPI_Init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, this%myRank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, this%ntasks , ierr)

        !
        ! set up PIO for rest of example
        !

        this%stride        = 1
        this%numAggregator = 0
        this%optBase       = 1
        this%iotype        = PIO_iotype_netcdf
        this%fileName      = "examplePio_f90.nc"
        this%dimLen(1)     = LEN

        this%niotasks = this%ntasks ! keep things simple - 1 iotask per MPI process

        !write(*,*) 'this%niotasks ',this%niotasks

        call PIO_init(this%myRank,      & ! MPI rank
            MPI_COMM_WORLD,             & ! MPI communicator
            this%niotasks,              & ! Number of iotasks (ntasks/stride)
            this%numAggregator,         & ! number of aggregators to use
            this%stride,                & ! stride
            PIO_rearr_subset,           & ! do not use any form of rearrangement
            this%pioIoSystem,           & ! iosystem
            base=this%optBase)            ! base (optional argument)

        ! 
        ! set up some data that we will write to a netcdf file
        !

        this%arrIdxPerPe = LEN / this%ntasks

        if (this%arrIdxPerPe < 1) then
            call this%errorHandle("Not enough work to distribute among pes", ERR_CODE)
        endif

        this%ista = this%myRank * this%arrIdxPerPe + 1
        this%isto = this%ista + (this%arrIdxPerPe - 1)

        allocate(this%compdof(this%ista:this%isto))
        allocate(this%dataBuffer(this%ista:this%isto))
        allocate(this%readBuffer(this%ista:this%isto))

        this%compdof(this%ista:this%isto) = (/(i, i=this%ista,this%isto, 1)/)
        this%dataBuffer(this%ista:this%isto) = this%myRank + VAL
        this%readBuffer(this%ista:this%isto) = 0

    end subroutine init

    !> @fn createDecomp
    !! @brief Create the decomposition for the example.
    !! This subroutine creates the decomposition for the example.
    subroutine createDecomp(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer(PIO_OFFSET_KIND) :: start(1)
        integer(PIO_OFFSET_KIND) :: count(1)

        start(1) = this%ista
        count(1) = this%arrIdxPerPe

        call PIO_initdecomp(this%pioIoSystem, PIO_int, this%dimLen, this%compdof(this%ista:this%isto), &
            this%iodescNCells)

    end subroutine createDecomp

    !> @brief Create netCDF output file.
    !! This subroutine creates the netCDF output file for the example.
    subroutine createFile(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: retVal

        retVal = PIO_createfile(this%pioIoSystem, this%pioFileDesc, this%iotype, trim(this%fileName), PIO_clobber)
        call this%errorHandle("Could not create "//trim(this%fileName), retVal)

    end subroutine createFile

    !> @brief Define the netCDF metadata.
    !! This subroutine defines the netCDF dimension and variable used
    !! in the output file.
    subroutine defineVar(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: retVal

        retVal = PIO_def_dim(this%pioFileDesc, 'x', this%dimLen(1) , this%pioDimId)
        call this%errorHandle("Could not define dimension x", retVal)

        retVal = PIO_def_var(this%pioFileDesc, 'foo', PIO_int, (/this%pioDimId/), this%pioVar)
        call this%errorHandle("Could not define variable foo", retVal)

        retVal = PIO_enddef(this%pioFileDesc)
        call this%errorHandle("Could not end define mode", retVal)

    end subroutine defineVar

    !> @brief Write the sample data to the output file.
    !! This subroutine writes the sample data array to the netCDF
    !! output file.
    subroutine writeVar(this)

    implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: retVal

        call PIO_write_darray(this%pioFileDesc, this%pioVar, this%iodescNCells, this%dataBuffer(this%ista:this%isto), retVal)
        call this%errorHandle("Could not write foo", retVal)
        call PIO_syncfile(this%pioFileDesc)

    end subroutine writeVar

    !> @brief Read the sample data from the output file.
    !! This subroutine reads the sample data array from the netCDF
    !! output file.
    subroutine readVar(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: retVal

        call PIO_read_darray(this%pioFileDesc, this%pioVar, this%iodescNCells,  this%readBuffer, retVal)
        call this%errorHandle("Could not read foo", retVal)

    end subroutine readVar

    !> @brief Close the netCDF output file.
    !! This subroutine closes the output file used by this example.
    subroutine closeFile(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        call PIO_closefile(this%pioFileDesc)

    end subroutine closeFile

    !> @brief Clean up resources.
    !! This subroutine cleans up resources used in the example. The
    !! ParallelIO and MPI libraries are finalized, and memory
    !! allocated in this example program is freed.
    subroutine cleanUp(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: ierr

        deallocate(this%compdof)
        deallocate(this%dataBuffer)
        deallocate(this%readBuffer)

        call PIO_freedecomp(this%pioIoSystem, this%iodescNCells)
        call PIO_finalize(this%pioIoSystem, ierr)
        call MPI_Finalize(ierr)

    end subroutine cleanUp

    !> @brief Handle errors.
    !! This subroutine is called if there is an error.
    subroutine errorHandle(this, errMsg, retVal)

        implicit none

        class(pioExampleClass), intent(inout) :: this
        character(len=*),       intent(in)    :: errMsg
        integer,                intent(in)    :: retVal

        if (retVal .ne. PIO_NOERR) then
            write(*,*) retVal,errMsg
            call PIO_closefile(this%pioFileDesc)
            call mpi_abort(MPI_COMM_WORLD,0,retVal)
        end if

    end subroutine errorHandle

end module pioExample

!> @brief Main execution of example code.
!! This is an example program for the ParallelIO library.
!!
!! This program creates a netCDF output file with the ParallelIO
!! library, then writes and reads some data to and from the file.
!!
!! This example does the following:
!!
!! - initialization initializes the MPI library, initializes the
!!   ParallelIO library with ::PIO_int(). Then allocate memory for a
!!   data array of sample data to write, and an array to read the data
!!   back into. Also allocate an array to hold decomposition
!!   information.
!!
!! - creation of decomposition by calling ::PIO_initdecomp().
!!
program main

    use pioExample, only : pioExampleClass

    implicit none

    type(pioExampleClass) :: pioExInst

    call pioExInst%init()
    call pioExInst%createDecomp()
    call pioExInst%createFile()
    call pioExInst%defineVar()
    call pioExInst%writeVar()
    call pioExInst%readVar()
    call pioExInst%closeFile()
    call pioExInst%cleanUp()

end program main
