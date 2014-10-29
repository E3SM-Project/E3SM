module pioExample

    ! 
    ! simple example of using PIO with netcdf
    !

    use pio, only : PIO_init, PIO_rearr_none, iosystem_desc_t, file_desc_t
    use pio, only : PIO_finalize, PIO_noerr, PIO_iotype_netcdf, PIO_createfile
    use pio, only : PIO_int,var_desc_t, PIO_redef, PIO_def_dim, PIO_def_var, PIO_enddef
    use pio, only : PIO_closefile, io_desc_t, PIO_initdecomp, PIO_write_darray
    use pio, only : PIO_freedecomp, PIO_clobber, PIO_readvar

    implicit none
    save
    private

    include 'mpif.h'

    integer, parameter :: LEN = 1024
    integer, parameter :: VAL = 42

    type, public :: pioExampleClass

        integer :: myRank
        integer :: ntasks
        integer :: niotasks
        integer :: stride        ! stride in the mpi rank between io tasks
        integer :: numAggregator
        integer :: optBase

        type(iosystem_desc_t) :: pioIoSystem
        type(file_desc_t)     :: pioFileDesc
        integer               :: iotype
        integer               :: pioDimId
        type(var_desc_t)      :: pioVar
        type(io_desc_t)       :: iodescNCells
        integer               :: dataBuffer(LEN)
        integer               :: readBuffer(LEN)
        integer               :: compdof(LEN)
        integer, dimension(1) :: dimLen

        character(len=255) :: fileName

    contains
        procedure, non_overridable, public :: init
        procedure,                  public :: createDecomp
        procedure,                  public :: createFile
        procedure,                  public :: defineVar
        procedure,                  public :: writeVar
        procedure,                  public :: readVar
        procedure,                  public :: closeFile
        procedure,                  public :: delete
        procedure,                  private :: errorHandle
    end type pioExampleClass

    !
    ! Fortran way of getting a user-defined ctor in the
    ! <inst> = <newClass> form.
    !
    interface pioExampleClass
        module procedure newPioExampleClass
    end interface

contains

    function newPioExampleClass()

        implicit none

        type(pioExampleClass) :: newPioExampleClass

    end function newPioExampleClass

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
        this%compdof       = (/ (i, i = 1,LEN) /)

        this%niotasks = this%ntasks/this%stride

        call PIO_init(this%myRank,      & ! MPI rank
            MPI_COMM_WORLD,             & ! MPI communicator
            this%niotasks,              & ! Number of iotasks (ntasks/stride)
            this%numAggregator,         & ! number of aggregators to use
            this%stride,                & ! stride
            PIO_rearr_none,             & ! do not use any form of rearrangement
            this%pioIoSystem,           & ! iosystem
            base=this%optBase)            ! base (optional argument)

        ! 
        ! set up some data that we will write and read from a netcdf file
        !

        this%dataBuffer = this%compdof * VAL

    end subroutine init

    subroutine createDecomp(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        call PIO_initdecomp(this%pioIoSystem, PIO_int, this%dimLen, this%compdof, this%iodescNCells)

    end subroutine createDecomp

    subroutine createFile(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: retVal
        character(len=255) :: errMsg

        retVal = PIO_createfile(this%pioIoSystem, this%pioFileDesc, this%iotype, trim(this%fileName),PIO_clobber)
        call this%errorHandle("Could not create "//trim(this%fileName), retVal)

    end subroutine createFile

    subroutine defineVar(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: retVal
        character(len=255) :: errMsg

        retVal = PIO_def_dim(this%pioFileDesc, 'x', this%dimLen(1) , this%pioDimId)
        call this%errorHandle("Could not define dimension x", retVal)

        retVal = PIO_def_var(this%pioFileDesc, 'foo', PIO_int, (/this%pioDimId/), this%pioVar)
        call this%errorHandle("Could not define variable foo", retVal)

        retVal = PIO_enddef(this%pioFileDesc)
        call this%errorHandle("Could not end define mode", retVal)

    end subroutine defineVar

    subroutine writeVar(this)

    implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: retVal

        call PIO_write_darray(this%pioFileDesc, this%pioVar, this%iodescNCells, this%dataBuffer, retVal, -1)
        call this%errorHandle("Could not write foo", retVal)

    end subroutine writeVar

    subroutine readVar(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: ierr
        integer :: retVal

        call PIO_read_darray(this%pioFileDesc, this%pioVar, this%iodescNCells,  this%readBuffer, retVal)

        write(*,*) 'After PIO read, elements 1:10 ',this%readBuffer(1:10)

    end subroutine readVar

    subroutine closeFile(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: ierr

        call PIO_closefile(this%pioFileDesc)

    end subroutine closeFile

    subroutine delete(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: ierr

        call PIO_freedecomp(this%pioIoSystem, this%iodescNCells)
        call PIO_finalize(this%pioIoSystem, ierr)
        call MPI_Finalize(ierr)

    end subroutine delete

    subroutine errorHandle(this, errMsg, retVal)

        implicit none

        class(pioExampleClass), intent(inout) :: this
        character(len=*), intent(in) :: errMsg
        integer, intent(in) :: retVal

        if (retVal .ne. PIO_NOERR) then
            write(*,*) retVal,errMsg
            call PIO_closefile(this%pioFileDesc)
            call mpi_abort(MPI_COMM_WORLD,0,retVal)
        end if

    end subroutine errorHandle


end module pioExample

program main

    use pioExample, only : pioExampleClass

    implicit none

    type(pioExampleClass) :: pioExInst

    pioExInst = pioExampleClass()
    call pioExInst%init()
    call pioExInst%createDecomp()
    call pioExInst%createFile()
    call pioExInst%defineVar()
    call pioExInst%writeVar()
    call pioExInst%readVar()
    call pioExInst%closeFile()
    call pioExInst%delete()

end program main
