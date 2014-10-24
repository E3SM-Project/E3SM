module pioExample

    ! 
    ! simple example of using PIO with netcdf
    !

    use pio, only : PIO_init, PIO_rearr_none, iosystem_desc_t, file_desc_t
    use pio, only : PIO_finalize

    implicit none
    save
    private

    include 'mpif.h'

    type, public :: pioExampleClass

        integer :: myRank
        integer :: ntasks
        integer :: niotasks
        integer :: stride        ! stride in the mpi rank between io tasks
        integer :: numAggregator
        integer :: optBase

        type(iosystem_desc_t) :: pio_iosystem
        type(file_desc_t)     :: pio_file

    contains
        procedure, non_overridable, public :: init
        procedure,                  public :: createAndWrite
        procedure,                  public :: readAndClose
        procedure,                  public :: delete
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

        write(*,*) ' pioExample::new  - ctor '

    end function newPioExampleClass

    subroutine init(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: ierr

        this%stride = 1
        this%numAggregator = 0
        this%optBase = 1

        call MPI_Init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, this%myRank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, this%ntasks , ierr)

        write(*,*) ' pioExample::init ',this%myRank,this%ntasks, this%stride

        this%niotasks = this%ntasks/this%stride

        !   call PIO_init
        call PIO_init(this%myRank,      & ! MPI rank
            MPI_COMM_WORLD,             & ! MPI communicator
            this%niotasks,              & ! Number of iotasks (ntasks/stride)
            this%numAggregator,         & ! number of aggregators to use
            this%stride,                & ! stride
            PIO_rearr_none,             & ! do not use any form of rearrangement
            this%pio_iosystem,          & ! iosystem
            base=this%optBase)            ! base (optional argument)

    end subroutine init

    subroutine createAndWrite(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: ierr



    end subroutine createAndWrite

    subroutine readAndClose(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: ierr



    end subroutine readAndClose

    subroutine delete(this)

        implicit none

        class(pioExampleClass), intent(inout) :: this

        integer :: ierr

        call PIO_finalize(this%pio_iosystem, ierr)
        call MPI_Finalize(ierr)

        write(*,*) ' pioExample::delete - dtor '

end subroutine delete

end module pioExample

program main

    use pioExample, only : pioExampleClass

    implicit none

    type(pioExampleClass) :: pioExInst

    pioExInst = pioExampleClass()
    call pioExInst%init()
    call pioExInst%delete()
    call pioExInst%createAndWrite()
    call pioExInst%readAndClose()

end program main
