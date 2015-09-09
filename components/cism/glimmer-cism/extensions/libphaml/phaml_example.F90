!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   phaml_example.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


module phaml_example  
contains
!need to specify which phaml module
!phaml_example is 1
    !------------------------------------------------------------------------------
    !SUBROUTINE: phaml_init
    !ARGUMENTS: cism_model, phaml_solution
    !DESCRIPTION:
    ! This subroutine initializes the user_mod variables before the phaml run.
    ! Technically, the phaml solution isn't needed, but in the future there might
    ! be other things which need to be initialized.  It is important to call 
    ! the user_close subroutine when finished.
    !------------------------------------------------------------------------------
    subroutine phaml_init(cism_model,phaml_solution)
        use phaml
        use phaml_user_mod
        use glide_types
        implicit none
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        
        call user_init(get_ewn(cism_model),get_nsn(cism_model), &
        get_dew(cism_model),get_dns(cism_model),1)
    end subroutine phaml_init

!-------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !SUBROUTINE: phaml_setup
    !ARGUMENTS: cism_model, phaml_solution
    !DESCRIPTION:
    ! This subroutine is if a nonlinear type of solution is desired where you want
    ! to initialize the problem without solving or closing the solution.
    !------------------------------------------------------------------------------
    subroutine phaml_setup(cism_model,phaml_solution)
        use phaml
        use phaml_user_mod
        use phaml_support
        use glide_types
        implicit none
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        
        
        !initialise phaml and the variables needed for solution
        call phaml_init(cism_model,phaml_solution)
        !generate the mesh file
        call make_ice_poly_file(cism_model)
        !set the mesh and create the phaml solution
        call phaml_create(phaml_solution,nproc=1,update_umod=.true., &
            triangle_files="mesh.1", &
            !spawn_form=DEBUG_SLAVE, &
            draw_grid_who=NO_ONE)
            
        !set the initial conditions array    
        uphaml = cism_model%phaml%uphaml
        
        !set the initial conditions
        !update usermod must be called at least twice before iconds
        !can be set with phaml_solve_pde
        call update_usermod(phaml_solution)
        !this sets the initial conditions
        call phaml_solve_pde(phaml_solution,            &
                     !max_vert=32000,                   &
                     task=SET_INITIAL,                  &
                     refterm=KEEP_NELEM,                &
                     max_eq=3200,                       &
                     error_estimator=INITIAL_CONDITION, &
                     print_header_who=NO_ONE,           &
                     print_trailer_who=NO_ONE,          &
                     degree=3,                          &
                     draw_grid_when=NEVER)
                     
        !necessary here???????????????
        !put the solution in the uphaml variable                 
        call phaml_getsolution(phaml_solution, cism_model%phaml%uphaml)           
    end subroutine phaml_setup

!-------------------------------------------------------------------------
    

    !------------------------------------------------------------------------------
    !SUBROUTINE: phaml_lin_evolve
    !ARGUMENTS: cism_model, phaml_solution
    !DESCRIPTION:
    ! This subroutine is if a linear type of solution is desired where you have
    ! initialized the problem but have not solved it yet. or closing the solution.
    !------------------------------------------------------------------------------
    subroutine phaml_lin_evolve(cism_model,phaml_solution)
        use phaml
        use phaml_user_mod
        use glide_types
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        
        call phaml_copy_soln_to_old(phaml_solution)
        call phaml_solve_pde(phaml_solution,            &
                         draw_grid_when=NEVER,          &
                         pause_at_start=.false.,        &
                         pause_after_phases=.false.,    &
                         task = SOLVE_ONLY,             &   !these two options make it
                         max_refsolveloop= 1,           &   !solve just one iteration
                         print_error_who=MASTER,        & 
                         mg_cycles=20,                  &
                         mg_tol=MG_NO_TOL,              &
                         print_header_who=NO_ONE,       &
                         print_trailer_who=NO_ONE)
                         
    end subroutine phaml_lin_evolve

!-------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !SUBROUTINE: phaml_nonlin_evolve
    !ARGUMENTS: cism_model, phaml_solution
    !DESCRIPTION:
    ! This subroutine is for a nonlinear type of solution where you have
    ! initialized the problem but have not solved it yet. or closing the solution.
    !------------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine phaml_nonlin_evolve(cism_model,phaml_solution)
        use phaml
        use phaml_user_mod
        use glide_types
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        
        call phaml_copy_soln_to_old(phaml_solution)
        call phaml_solve_pde(phaml_solution,            &
                         draw_grid_when=NEVER,          &
                         pause_at_start=.false.,        &
                         pause_after_phases=.false.,    &
                         task = SOLVE_ONLY,             &   !these two options make it
                         max_refsolveloop= 1,           &   !solve just one iteration
                         print_error_who=MASTER,        & 
                         mg_cycles=20,                  &
                         mg_tol=MG_NO_TOL,              &
                         print_header_who=NO_ONE,       &
                         print_trailer_who=NO_ONE)
                         
    end subroutine phaml_nonlin_evolve

!-------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !SUBROUTINE: phaml_evolve
    !ARGUMENTS: cism_model, phaml_solution
    !DESCRIPTION:
    ! Once the init function has been called this function will create the mesh,
    ! initialize the problem, solve it, then store the data in the phaml custom
    ! type in the CISM model.  After this is finished the phaml_close subroutine
    ! must be called.
    !------------------------------------------------------------------------------
    subroutine phaml_evolve(cism_model,phaml_solution)
        use phaml
        use phaml_user_mod
        use phaml_support
        use glide_types
        implicit none
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        !initialise phaml and the variables needed for solution
        !call phaml_init(cism_model,phaml_solution)
        !generate the mesh file
        call make_ice_poly_file(cism_model)
        
        !set the mesh and create the phaml solutiony
        call phaml_create(phaml_solution,nproc=1,update_umod=.true., &
            triangle_files="mesh.1", &
            spawn_form=DEBUG_SLAVE, &!DEBUG_BOTH, &
            draw_grid_who=NO_ONE)
            
        !set the initial conditions array    
        uphaml = cism_model%phaml%uphaml
        !set the initial conditions
        !update usermod must be called at least twice before iconds
        !can be set with phaml_solve_pde
        call update_usermod(phaml_solution)
        write(*,*) 'phaml_solve'
        !solve the pde
        call phaml_solve_pde(phaml_solution,                      &
                         task=SET_INITIAL,          &
                         refterm=KEEP_NELEM,        &
                         max_eq=3200,              &
                         error_estimator=INITIAL_CONDITION, &
                         print_header_who=NO_ONE,   &
                         print_trailer_who=NO_ONE,  &
                         degree=3,                  &
                         draw_grid_when=NEVER)
        call phaml_solve_pde(phaml_solution,            &
                         draw_grid_when=NEVER,          &
                         pause_at_start=.false.,        &
                         pause_after_phases=.false.,    &
                         print_error_who=MASTER, & 
                         task = SOLVE_ONLY,            &   !these two options make it
                         max_refsolveloop= 1,          &   !solve just one iteration
                         mg_cycles=20,                  &
                         mg_tol=MG_NO_TOL,              &
                         print_header_who=NO_ONE,       &
                         print_trailer_who=NO_ONE)
        write(*,*) 'done'                         
        !put the solution in the uphaml variable                 
        call phaml_getsolution(phaml_solution, cism_model%phaml%uphaml)
        !close the phaml session
        !call phaml_close(phaml_solution)
    end subroutine phaml_evolve

!-------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------
    !SUBROUTINE: phaml_getsolution
    !ARGUMENTS: phaml_solution, uout
    !DESCRIPTION:
    ! In order to get the solution of the pde an x and y vector must be passed to
    ! the phaml_evaluate function and the corresponding solution is returned in a u
    ! vector.  This function takes care of creating those arrays and reshaping them 
    ! to a 2d array which it returns in uout.
    !------------------------------------------------------------------------------
    subroutine phaml_getsolution(phaml_solution, uout)
        use phaml
        use phaml_user_mod
        implicit none
        type(phaml_solution_type) :: phaml_solution
        real(my_real), intent(out), dimension(:,:) :: uout
        real(my_real), allocatable, dimension(:) :: x,y,u
        
        !the globals are set in init and in phaml_user_mod
        allocate(x(gewn*gnsn))
        allocate(y(gewn*gnsn))
        allocate(u(gewn*gnsn))
        !this sets up the locations of our grid points to pass to phaml_evaluate  
        call get_xyarrays(x,y)
        !get the initial solution in u
        call phaml_evaluate(phaml_solution,x,y,u)
        !write solution to uphaml in the model
        call reshape_array_to_two(uout,u)
        !free memory
        deallocate(x)
        deallocate(y)
        deallocate(u)
    end subroutine phaml_getsolution
!-------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------
    !SUBROUTINE: phaml_close
    !ARGUMENTS: phaml_solution
    !DESCRIPTION:
    ! This subroutine simply destroys the phaml solution instance as well as 
    ! deallocating all of the user data.
    !------------------------------------------------------------------------------
    subroutine phaml_close(phaml_solution)
        use phaml
        use phaml_user_mod
        implicit none
        type(phaml_solution_type) :: phaml_solution
        
        !destroy phaml session
        call phaml_destroy(phaml_solution)
        call user_close()
    end subroutine phaml_close

!-------------------------------------------------------------------------


end module phaml_example
