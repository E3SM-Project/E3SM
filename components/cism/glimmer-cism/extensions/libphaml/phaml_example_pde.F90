
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   phaml_example_pde.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module phaml_example_pde
    !----------------------------------------------------
    ! This file contains the user supplied external subroutines that define
    ! the PDE(s) to be solved, and other required external subroutines.
    !   pdecoefs bconds boundary_point boundary_npiece boundary_param iconds trues
    !   truexs trueys update_usermod phaml_integral_kernel
    !----------------------------------------------------
    use phaml
contains    
    !-------------------------------------------------------------------------
    
    subroutine example_pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
        !----------------------------------------------------
        ! This subroutine returns the coefficient and right hand side of the PDE
        ! at the point (x,y)
        !
        ! The PDE is
        !
        !    -( cxx(x,y) * u  )  -( cyy(x,y) * u  ) + c(x,y) * u = rs(x,y)
        !                   x  x                y  y
        !
        ! For eigenvalue problems, the right hand side is lambda * u * rs(x,y)
        !
        ! cxy, cx and cy are not currently used and should not be set.  They are
        ! for future expansion.
        !
        ! NOTE: BE CAREFUL TO GET THE SIGNS RIGHT
        ! e.g. cxx=cyy=1 means rs=-(uxx+uyy)
        !
        !----------------------------------------------------
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                      c(:,:),rs(:)
        !----------------------------------------------------
        ! Begin executable code
        
        cxx(1,1) = 1.0_my_real
        cyy(1,1) = 1.0_my_real
        c(1,1) = 0.0_my_real
        rs(1) = 0.0_my_real
        
        cxy=0; cx=0; cy=0
    
    end subroutine example_pdecoefs
    
    !-------------------------------------------------------------------------
    subroutine example_bconds(x,y,bmark,itype,c,rs)
        !----------------------------------------------------
        ! This subroutine returns the boundary conditions at the point (x,y).
        !
        ! Each boundary condition is either
        !
        !    u = rs(x,y) or  u  + c(x,y)*u = rs(x,y)
        !                     n
        !
        ! itype indicates which type of condition applies to each point, using
        ! symbolic constants from module phaml.  It is DIRICHLET for the first
        ! condition, NATURAL for the second condition with c==0, and MIXED for
        ! the second condition with c/=0.
        !
        !----------------------------------------------------
        use phaml_user_mod
        use glide_mask
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        integer, intent(in) :: bmark
        integer, intent(out) :: itype(:)
        real(my_real), intent(out) :: c(:,:),rs(:)
        integer :: ew,ns
        real :: middle
        !----------------------------------------------------
        ! Non-module procedures used are:
        
        interface
        
           function trues(x,y,comp,eigen) ! real (my_real)
           use phaml
           real (my_real), intent(in) :: x,y
           integer, intent(in) :: comp,eigen
           real (my_real) :: trues
           end function trues
        
        end interface
        
        !----------------------------------------------------
        ! Begin executable code
        ! Dirichlet boundary conditions
        itype = DIRICHLET    
        rs(1) = 0.0_my_real
        c = 0.0_my_real
        !if(is_grounding_line(bmark)) then
        !    rs(1) = 1.0_my_real
        !end if
        middle = REAL(gdew*(gewn*0.5))
        if (abs(x-y)< gdew*.2 .and. &
            (x .gt. middle) .and. &
            (y .gt. middle)) then
            rs=1.0_my_real
        end if

    end subroutine example_bconds
    
    !-------------------------------------------------------------------------
    function example_iconds(x,y,comp,eigen)
        !----------------------------------------------------
        ! This routine returns the initial condition for a time dependent problem.
        ! It can also be used for the initial guess for a nonlinear problem, or
        ! systems of equations solved by successive substitution.
        ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
        ! eigenvectors from an eigenvalue problem, and is ignored in this example.
        ! For problems where there are no initial conditions, it is a dummy.
        !----------------------------------------------------
        use phaml_user_mod
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        integer, intent(in) :: comp,eigen
        real(my_real) :: ret_value
        real(my_real) :: example_iconds
        integer :: ew, ns
        
        !maps an x, y to nearest integer r,c value for array lookup
        !in general, DO NOT use these functions in this way
        ew = getew(x)
        ns = getns(y)
        !----------------------------------------------------
        ! Begin executable code
        ret_value = 0.0
        if (uphaml(ew,ns) .gt. 0.0) then
            ret_value = uphaml(ew,ns)
        end if
        example_iconds = ret_value
    
    end function example_iconds
    
    !-------------------------------------------------------------------------
    function example_trues(x,y,comp,eigen) ! real (my_real)
        !----------------------------------------------------
        ! This is the true solution of the differential equation, if known.
        ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
        ! eigenvectors from an eigenvalue problem, and is ignored in this example.
        !----------------------------------------------------
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        integer, intent(in) :: comp,eigen
        real (my_real) :: example_trues
        !----------------------------------------------------
        ! Begin executable code
        
        !trues = x**2 + y**2
        example_trues = 0.0_my_real !return 0 since we don't know
    
    end function example_trues
    
    !-------------------------------------------------------------------------
    function example_truexs(x,y,comp,eigen) ! real (my_real)
        !----------------------------------------------------
        ! This is the x derivative of the true solution of the differential
        ! equation, if known.
        ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
        ! eigenvectors from an eigenvalue problem, and is ignored in this example.
        !----------------------------------------------------
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        integer, intent(in) :: comp,eigen
        real (my_real) :: example_truexs
        !----------------------------------------------------
        ! Begin executable code
        
        example_truexs = 0.0_my_real
    
    end function example_truexs
    
    !-------------------------------------------------------------------------
    function example_trueys(x,y,comp,eigen) ! real (my_real)
        !----------------------------------------------------
        ! This is the y derivative of the true solution of the differential
        ! equation, if known.
        ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
        ! eigenvectors from an eigenvalue problem, and is ignored in this example.
        !----------------------------------------------------
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        integer, intent(in) :: comp,eigen
        real (my_real) :: example_trueys
        !----------------------------------------------------
        ! Begin executable code
        
        example_trueys = 0.0_my_real
    
    end function example_trueys
    
    !-------------------------------------------------------------------------
    subroutine example_boundary_point(ipiece,s,x,y)
        !----------------------------------------------------
        ! This routine defines the boundary of the domain.  It returns the point
        ! (x,y) with parameter s on piece ipiece.
        ! If boundary_npiece <= 0 it is a dummy and the domain is given by triangle
        ! data files.
        !----------------------------------------------------
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        integer, intent(in) :: ipiece
        real(my_real), intent(in) :: s
        real(my_real), intent(out) :: x,y
        !----------------------------------------------------
        ! Begin executable code
        ! Dummy version
        
        x = 0.0_my_real
        y = 0.0_my_real
    
    end subroutine example_boundary_point
    
    !-------------------------------------------------------------------------
    function example_boundary_npiece(hole)
        !----------------------------------------------------
        ! This routine gives the number of pieces in the boundary definition.
        ! If boundary_npiece <= 0 the domain is given by triangle data files.
        !----------------------------------------------------
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        integer, intent(in) :: hole
        integer :: example_boundary_npiece
        !----------------------------------------------------
        ! Begin executable code
        !write(*,*) 'boundary_npiece'
        example_boundary_npiece = 0
    
    end function example_boundary_npiece
    
    !-------------------------------------------------------------------------
    subroutine example_boundary_param(start,finish)
        !----------------------------------------------------
        ! This routine gives the range of parameter values for each piece of the
        ! boundary.
        ! If boundary_npiece <= 0 it is a dummy and the domain is given by triangle
        ! data files.
        !----------------------------------------------------
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(out) :: start(:), finish(:)
        !----------------------------------------------------
        ! Begin executable code
        
        ! Dummy version
        
        start = 0.0_my_real
        finish = 0.0_my_real
    
    end subroutine example_boundary_param
    
    
    !-------------------------------------------------------------------------
    function example_phaml_integral_kernel(kernel,x,y)
        integer, intent(in) :: kernel
        real(my_real), intent(in) :: x,y
        real(my_real) :: example_phaml_integral_kernel
        
        ! Identity function
        
        example_phaml_integral_kernel = 1.0
    
    end function example_phaml_integral_kernel
    
    !-------------------------------------------------------------------------
    function example_regularity(x,y)
        real(my_real), intent(in) :: x(3),y(3)
        real(my_real) :: example_regularity
        
        ! Dummy version, assume infinitely differentiable everywhere.
        
        example_regularity = 0.0_my_real
    
    end function example_regularity


    subroutine example_update_usermod(phaml_solution) 
        !----------------------------------------------------
        ! This routine updates the module variables on the slave processes by
        ! sending them from the master process
        !----------------------------------------------------
        use phaml_user_mod
        implicit none
        type(phaml_solution_type), intent(in) :: phaml_solution
        integer :: iparam(6)
        real(my_real),allocatable,dimension(:) :: rparam
     
        
        !the slaves need the globals passed first so that the allocate
        !for rparam has the values for the size.  This means after
        !phaml_create another call to update_usermod must be made
        
        allocate(rparam(gnsn*gewn*num_arrays))
        iparam(1) = gnsn
        iparam(2) = gewn
        iparam(3) = gdns
        iparam(4) = gdew
        iparam(5) = num_arrays
        iparam(6) = modnum
        !if more arrays are needed reshape and combine
        call reshape_array_to_one(uphaml,rparam)
        !Call the routine that performs the actual exchange.   
        call master_to_slaves(phaml_solution,iparam,rparam)
        gnsn = iparam(1)
        gewn = iparam(2)
        gdns = iparam(3)
        gdew = iparam(4)
        num_arrays = iparam(5)
        modnum = iparam(6)
        if(num_arrays .eq. 1) then
            call reshape_array_to_two(uphaml,rparam)
        end if
        !an else should be added here to handle passing more than one array
        !using the concat_arrays, split_arrays function in combination with the
        !combine and separate functions.
        deallocate(rparam)
    end subroutine example_update_usermod

end module phaml_example_pde
