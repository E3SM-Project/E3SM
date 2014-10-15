module pioExample

    implicit none
    save
    private

    type, public :: pioExampleClass

    contains
        procedure, non_overridable, public :: init
        !procedure, public :: write
        !procedure, public :: read
        !procedure, public :: delete
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

        write(*,*) ' pioExample::init  - init '

    end subroutine init

end module pioExample

program main

    use pioExample, only : pioExampleClass

    implicit none

    type(pioExampleClass) :: pioExInst

    pioExInst = pioExampleClass()
    call pioExInst%init()

end program main