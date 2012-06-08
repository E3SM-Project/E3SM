#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!------------------------------------------------------------------------------
! High Order Method Modeling Environment (HOMME)
!------------------------------------------------------------------------------
!
! MODULE: Manager
!
!
! DESCRIPTION:
!> @brief The Manager controls global shared resources in applications using HOMME.
!!
!! This module is designed after the <b>singleton</b> pattern for Object Oriented programming.
!! See \a Design \a Patterns by Gamma, Helm, Johnson and Vlissides.
!! The idea is to avoid passing resources from high level routines in HOMME to lower
!! level routines via subroutines parameters. \n
!! For example, if one has an array of 'cell'
!! which is utilized everywhere in HOMME, rather than forcing the arguments of the routines
!! to take a 'cell' argument and passed forward like this
!! \code{.F90}
!! subroutine higher(..., cell my_cell, ...)
!!    ...some code  doing something...
!!          \! Pass along the argument cell
!!          call lower(.., my_cell, ...)
!! end subroutine higher
!! \endcode
!! the function at the higher level can avoid the parameter my_cell all together and instead the function
!! at the lower level, the one that really needs the resource obtains the resource like this
!! \code{.F90}
!! subroutine lower(...)
!!   use Manager, only: cell_get
!!    ...some code  doing something...
!!          \! get the cell resource
!!          call cell_get(cell)
!! end subroutine lower
!! \endcode
!! Notice the concept of a resource is an entity in the program that is unique and that can be obtained with
!! a handle. In Fortran 90 we implement this idea via pointers to the resource, so a caller must create
!! a pointer outside, nullify it, and call the appropiatte subroutine to get its
!! own pointer assigned to the resource. For example,
!! \code{.F90}
!! subroutine InitColumnModel(elem, cm,hvcoord,hybrid,tl,nets,nete,runtype)
!!
!!    use Manager
!!    ...
!!    nullify(elem_physics)
!!    call element_physics_get(elem_physics)
!!    ...
!! end subroutine InitColumnModel
!! \endcode
!! Important
!! 1. Clearly this is not thread safe by default unless you use a mutex to control the resource. \n
!! 2. The naming convention for the subroutines to "get" the resources is \<resource name\>\_get. \n
!! 3. If a pointer is not null when using a get routine, the program will abort because we can't override
!!    a valid address.
!! \n
!! \n
!> @author Jose Garcia (jgarcia@ucar.edu). NCAR
!------------------------------------------------------------------------------

module Manager

    use physics_mod, only : elem_physics_t
    use control_mod, only : physics
    implicit none

    private
    !> @brief This is the physics resource.
    type (elem_physics_t), pointer  :: my_elem_physics(:)

    ! --------------------------------- Private members ------------------------------------------
    private:: element_physics_set
    ! --------------------------------- Public members ------------------------------------------
    public:: ManagerInit
    public:: element_physics_get


contains


    !> Initialize the physics resource to the total number of elements (dimensions_mod::nelemd) for the corresponding MPI rank.
    subroutine element_physics_set()

        use dimensions_mod, only : nelemd

        allocate(my_elem_physics(nelemd))

    end subroutine element_physics_set

    !> Assign to the given pointer the address of the elem_physics resource.
    subroutine element_physics_get(elem_physics)
        use parallel_mod, only: abortmp
        use physics_mod,  only : elem_physics_t

        type (elem_physics_t), dimension(:), pointer, intent(out)  :: elem_physics

        if  (physics == 0) then
            call abortmp('The variable physics is set to zero, yet you wish to get a handle to the physics resource. Not possible, aborting!')
        endif
        if (associated(elem_physics)) then
            call abortmp('Pointer <elem_physics> is currently associated, setting it will cause a memory leak, aborting.')
        endif
        elem_physics => my_elem_physics

    end subroutine element_physics_get

    !> Initialize the manager so resources are created and allocated properly.
    subroutine ManagerInit()
        use parallel_mod, only: abortmp

        if (physics == 0) then
            print *, 'No physics package prescribed, do not attempt to get access to the physics resource, it is uninitialized'
        else if (physics == 1) then
            ! for now only the multicloud package works.
            call element_physics_set()
        else
            call abortmp('Physics package values can be either 0 or 1, aborting!')
        endif
    end subroutine ManagerInit

end module Manager
