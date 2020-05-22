!> @file
!! A simple Fortran example for the ParallelIO Library.
module pioExample

    use piolib_mod, only : PIO_init, PIO_finalize, PIO_createfile, PIO_closefile, &
        PIO_initdecomp, PIO_freedecomp, PIO_syncfile, PIO_openfile, PIO_setframe
    use pio_types,  only : iosystem_desc_t, file_desc_t, &
        pio_noerr, PIO_iotype_netcdf, var_desc_t, io_desc_t, PIO_int, &
        pio_clobber, PIO_nowrite, PIO_unlimited
    use pio_kinds,  only : PIO_OFFSET_KIND, i4 
    use pio_nf,     only :  PIO_redef, PIO_def_dim, PIO_def_var, PIO_enddef
    use piodarray,  only : PIO_write_darray, PIO_read_darray 

    implicit none

    private

    include 'mpif.h'

    integer(i4), public, parameter :: PIO_rearr_subset =  2
    !> @brief Length of the data array we are using.  This is then
    !! divided among MPI processes.
    integer, parameter :: LEN = 16

    !> @brief Value used for array that will be written to netcdf file.
    integer, parameter :: VAL = 42  

    !> @brief Error code if anything goes wrong.
    integer, parameter :: ERR_CODE = 99

    !> mpi communication group to be used by module
    integer :: mpicom
    !> @brief A class to hold example code and data.
    !! This class contains the data and functions to execute the
    !! example.
    type, public :: pioExampleClass

        !> @brief Rank of processor running the code.
        integer :: myRank

        !> @brief Number of processors participating in MPI communicator.
        integer :: ntasks

        !> @brief Number of processors performing I/O.
        integer :: niotasks

        !> @brief Stride in the mpi rank between io tasks.
        integer :: stride        

        !> @brief Number of aggregator.
        integer :: numAggregator

        !> @brief Start index of I/O processors.
        integer :: optBase

        !> @brief The ParallelIO system set up by @ref PIO_init.
        type(iosystem_desc_t) :: pioIoSystem

        !> @brief Contains data identifying the file.
        type(file_desc_t)     :: pioFileDesc

        !> @brief The netCDF variable ID.
        type(var_desc_t)      :: pioVar

        !> @brief An io descriptor handle that is generated in @ref PIO_initdecomp.
        type(io_desc_t)       :: iodescNCells

        !> @brief Specifies the flavor of netCDF output.
        integer               :: iotype

        !> @brief The netCDF dimension ID.
        integer               :: pioDimId(2)

        !> @brief 1-based index of start of this processors data in full data array.
        integer               :: ista

        !> @brief Size of data array for this processor.
        integer               :: isto

        !> @brief Number of elements handled by each processor.
        integer               :: arrIdxPerPe

        !> @brief The length of the dimension of the netCDF variable.
        integer, dimension(2) :: dimLen

        !> @brief Buffer to hold sample data that is written to netCDF file.
        integer, allocatable  :: dataBuffer(:)

        !> @brief Buffer to read data into.
        integer, allocatable  :: readBuffer(:)

        !> @brief Array describing the decomposition of the data.
        integer, allocatable  :: compdof(:)

        !> @brief Name of the sample netCDF file written by this example.
        character(len=255) :: fileName

    contains

        !> @brief Initialize MPI, ParallelIO, and example data.
        !! Initialize the MPI and ParallelIO libraries. Also allocate
        !! memory to write and read the sample data to the netCDF file.
        procedure,  public  :: init

        !> @brief Create the decomposition for the example.
        !! This subroutine creates the decomposition for the example.
        procedure,  public  :: createDecomp

        !> @brief Create netCDF output file.
        !! This subroutine creates the netCDF output file for the example.
        procedure,  public  :: createFile

        !> @brief Define the netCDF metadata.
        !! This subroutine defines the netCDF dimension and variable used
        !! in the output file.
        procedure,  public  :: defineVar

        !> @brief Write the sample data to the output file.
        !! This subroutine writes the sample data array to the netCDF
        !! output file.
        procedure,  public  :: writeVar

        !> @brief Read the sample data from the output file.
        !! This subroutine reads the sample data array from the netCDF
        !! output file.
        procedure,  public  :: readVar

        !> @brief Close the netCDF output file.
        !! This subroutine closes the output file used by this example.
        procedure,  public  :: closeFile

        !> @brief Clean up resources.
        !! This subroutine cleans up resources used in the example. The
        !! ParallelIO and MPI libraries are finalized, and memory
        !! allocated in this example program is freed.
        procedure,  public  :: cleanUp

        !> @brief Handle errors.
        !! This subroutine is called if there is an error.
        procedure,  private :: errorHandle

    end type pioExampleClass

contains

    !> @brief Initialize MPI, ParallelIO, and example data.
    !! Initialize the MPI and ParallelIO libraries. Also allocate
    !! memory to write and read the sample data to the netCDF file.
    subroutine init(this,mpicom_in)

        implicit none

        class(pioExampleClass), intent(inout) :: this
        integer,intent(IN)            :: mpicom_in              ! mpi communicator

        integer :: ierr,i

        !
        ! initialize MPI
        !

!        call MPI_Init(ierr)
        mpicom = mpicom_in
        call MPI_Comm_rank(mpicom, this%myRank, ierr)
        call MPI_Comm_size(mpicom, this%ntasks , ierr)

        !
        ! set up PIO for rest of example
        !

        this%stride        = 1
        this%numAggregator = 0
        this%optBase       = 1
        this%iotype        = PIO_iotype_netcdf
        this%fileName      = "examplePio_f90.nc"
        this%dimLen(1)     = LEN
        this%dimLen(2)     = 4

        this%niotasks = this%ntasks ! keep things simple - 1 iotask per MPI process

        write(*,*) 'this%myRank ',this%myRank
        write(*,*) 'this%niotasks ',this%niotasks
        write(*,*) 'this%numAggregator ',this%numAggregator
        write(*,*) 'this%stride ',this%stride
        write(*,*) 'this%pioIoSystem ',this%pioIoSystem
        write(*,*) 'this%iotype ',this%iotype
        write(*,*) 'this%optBase ',this%optBase
        write(*,*) 'mpicom ',mpicom

        call PIO_init(this%myRank,      & ! MPI rank
            mpicom,             & ! MPI communicator
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
        if (this%myRank == this%ntasks-1) this%isto = LEN

        allocate(this%compdof(this%ista:this%isto))
        allocate(this%dataBuffer(this%ista:this%isto))
        allocate(this%readBuffer(this%ista:this%isto))

        this%compdof(this%ista:this%isto) = (/(i, i=this%ista,this%isto, 1)/)
        this%dataBuffer(this%ista:this%isto) = this%myRank + VAL
        this%readBuffer(this%ista:this%isto) = 0

    end subroutine init

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

    subroutine createFile(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: retVal

        retVal = PIO_createfile(this%pioIoSystem, this%pioFileDesc, this%iotype, trim(this%fileName), PIO_clobber)
        call this%errorHandle("Could not create "//trim(this%fileName), retVal)

    end subroutine createFile

    subroutine defineVar(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: retVal

        retVal = PIO_def_dim(this%pioFileDesc, 'x', this%dimLen(1) , this%pioDimId(1))
        call this%errorHandle("Could not define dimension x", retVal)

        retVal = PIO_def_dim(this%pioFileDesc, 't', pio_unlimited , this%pioDimId(2))
        call this%errorHandle("Could not define dimension x", retVal)

        retVal = PIO_def_var(this%pioFileDesc, 'foo', PIO_int, (/this%pioDimId/), this%pioVar)
        call this%errorHandle("Could not define variable foo", retVal)

        retVal = PIO_enddef(this%pioFileDesc)
        call this%errorHandle("Could not end define mode", retVal)

    end subroutine defineVar

    subroutine writeVar(this,step)

    implicit none

        class(pioExampleClass), intent(inout) :: this
!        integer, intent(in), optional :: step

        integer(kind=pio_offset_kind), intent(in), optional :: step
        integer :: retVal

        if (present(step)) then
          call PIO_setframe(this%pioFileDesc,this%pioVar,step)
          this%dataBuffer(this%ista:this%isto) = this%dataBuffer(this%ista:this%isto) + 100
        end if
        call PIO_write_darray(this%pioFileDesc, this%pioVar, this%iodescNCells, this%dataBuffer(this%ista:this%isto), retVal)
        call this%errorHandle("Could not write foo", retVal)
        call PIO_syncfile(this%pioFileDesc)

    end subroutine writeVar

    subroutine readVar(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: retVal

        call PIO_read_darray(this%pioFileDesc, this%pioVar, this%iodescNCells,  this%readBuffer, retVal)
        call this%errorHandle("Could not read foo", retVal)

    end subroutine readVar

    subroutine closeFile(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        call PIO_closefile(this%pioFileDesc)

    end subroutine closeFile

    subroutine cleanUp(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: ierr

        deallocate(this%compdof)
        deallocate(this%dataBuffer)
        deallocate(this%readBuffer)

        call PIO_freedecomp(this%pioIoSystem, this%iodescNCells)
        call PIO_finalize(this%pioIoSystem, ierr)

    end subroutine cleanUp

    subroutine errorHandle(this, errMsg, retVal)

        implicit none

        class(pioExampleClass), intent(inout) :: this
        character(len=*),       intent(in)    :: errMsg
        integer,                intent(in)    :: retVal

        if (retVal .ne. PIO_NOERR) then
            write(*,*) retVal,errMsg
            call PIO_closefile(this%pioFileDesc)
            call mpi_abort(mpicom,0,retVal)
        end if

    end subroutine errorHandle

end module pioExample
