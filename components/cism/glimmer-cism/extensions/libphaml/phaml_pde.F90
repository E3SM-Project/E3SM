
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   phaml_pde.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!----------------------------------------------------
! This file contains the user supplied external subroutines that define
! the PDE(s) to be solved, and other required external subroutines.
!   pdecoefs bconds boundary_point boundary_npiece boundary_param iconds trues
!   truexs trueys update_usermod phaml_integral_kernel
!----------------------------------------------------

!-------------------------------------------------------------------------
    
!> This subroutine returns the coefficient and right hand side of the PDE
!! at the point (x,y) 

subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
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
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    implicit none
    !----------------------------------------------------
    ! Dummy arguments
    real(my_real), intent(in) :: x,y
    real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                  c(:,:),rs(:)
                          
    !----------------------------------------------------
    ! Begin executable code
    if(modnum .eq. 1) then 
        call example_pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
    else
        cxx(1,1) = 1.0_my_real
        cyy(1,1) = 1.0_my_real
        c(1,1) = 0.0_my_real
        rs(1) = 0.0_my_real
    end if
    
    cxy=0; cx=0; cy=0

end subroutine pdecoefs

!-------------------------------------------------------------------------

!> This subroutine returns the boundary conditions at the point (x,y).
subroutine bconds(x,y,bmark,itype,c,rs)
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
    use phaml
    use phaml_example_pde
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
    if(modnum .eq. 1) then 
        call example_bconds(x,y,bmark,itype,c,rs)
    else
        itype = DIRICHLET
        rs(1) = 0.0_my_real
        c = 0.0_my_real
    end if

end subroutine bconds

!-------------------------------------------------------------------------
    
!> This routine returns the initial condition for a time dependent problem.
function iconds(x,y,comp,eigen)
    !----------------------------------------------------
    ! This routine returns the initial condition for a time dependent problem.
    ! It can also be used for the initial guess for a nonlinear problem, or
    ! systems of equations solved by successive substitution.
    ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
    ! eigenvectors from an eigenvalue problem, and is ignored in this example.
    ! For problems where there are no initial conditions, it is a dummy.
    !----------------------------------------------------
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    implicit none
    !----------------------------------------------------
    ! Dummy arguments
    real(my_real), intent(in) :: x,y
    integer, intent(in) :: comp,eigen
    real(my_real) :: ret_value
    real(my_real) :: iconds

    if(modnum .eq. 1) then
        ret_value = example_iconds(x,y,comp,eigen)
    else
        ret_value = 0.0
    end if
    iconds = ret_value
end function iconds

!-------------------------------------------------------------------------

!> This is the true solution of the differential equation, if known.
function trues(x,y,comp,eigen) ! real (my_real)
    !----------------------------------------------------
    ! This is the true solution of the differential equation, if known.
    ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
    ! eigenvectors from an eigenvalue problem, and is ignored in this example.
    !----------------------------------------------------
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    implicit none
    !----------------------------------------------------
    ! Dummy arguments
    real(my_real), intent(in) :: x,y
    integer, intent(in) :: comp,eigen
    real (my_real) :: trues
    !----------------------------------------------------
    ! Begin executable code
    
    if(modnum .eq. 1) then 
        trues = example_trues(x,y,comp,eigen)
    else
        trues = huge(0.0_my_real) !return 0 since we don't know
    end if
end function trues

!-------------------------------------------------------------------------
    
!> This is the x derivative of the true solution of the differential equation, if known

function truexs(x,y,comp,eigen) ! real (my_real)
    !----------------------------------------------------
    ! This is the x derivative of the true solution of the differential
    ! equation, if known.
    ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
    ! eigenvectors from an eigenvalue problem, and is ignored in this example.
    !----------------------------------------------------
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    implicit none
    !----------------------------------------------------
    ! Dummy arguments
    real(my_real), intent(in) :: x,y
    integer, intent(in) :: comp,eigen
    real (my_real) :: truexs
    !----------------------------------------------------
    ! Begin executable code
    if(modnum .eq. 1) then 
        truexs = example_truexs(x,y,comp,eigen)
    else
        truexs = huge(0.0_my_real)
    end if

end function truexs

!-------------------------------------------------------------------------

!> This is the y derivative of the true solution of the differential equation, if known

function trueys(x,y,comp,eigen) ! real (my_real)
    !----------------------------------------------------
    ! This is the y derivative of the true solution of the differential
    ! equation, if known.
    ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
    ! eigenvectors from an eigenvalue problem, and is ignored in this example.
    !----------------------------------------------------
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    implicit none
    !----------------------------------------------------
    ! Dummy arguments
    real(my_real), intent(in) :: x,y
    integer, intent(in) :: comp,eigen
    real (my_real) :: trueys
    !----------------------------------------------------
    ! Begin executable code
    if(modnum .eq. 1) then 
        trueys = example_trueys(x,y,comp,eigen)
    else
        trueys = huge(0.0_my_real)
    end if
end function trueys

!-------------------------------------------------------------------------
!> This routine defines the boundary of the domain.
subroutine boundary_point(ipiece,s,x,y)
    !----------------------------------------------------
    ! This routine defines the boundary of the domain.  It returns the point
    ! (x,y) with parameter s on piece ipiece.
    ! If boundary_npiece <= 0 it is a dummy and the domain is given by triangle
    ! data files.
    !----------------------------------------------------
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    implicit none
    !----------------------------------------------------
    ! Dummy arguments
    integer, intent(in) :: ipiece
    real(my_real), intent(in) :: s
    real(my_real), intent(out) :: x,y
    !----------------------------------------------------
    ! Begin executable code
    ! Dummy version
    if(modnum .eq. 1) then 
        call example_boundary_point(ipiece,s,x,y)
    else
        x = 0.0_my_real
        y = 0.0_my_real
    end if
end subroutine boundary_point

!-------------------------------------------------------------------------

!> This routine gives the number of pieces in the boundary definition.    

function boundary_npiece(hole)
    !----------------------------------------------------
    ! This routine gives the number of pieces in the boundary definition.
    ! If boundary_npiece <= 0 the domain is given by triangle data files.
    !----------------------------------------------------
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    implicit none
    !----------------------------------------------------
    ! Dummy arguments
    integer, intent(in) :: hole
    integer :: boundary_npiece
    !----------------------------------------------------
    ! Begin executable code
    if(modnum .eq. 1) then 
        boundary_npiece = example_boundary_npiece(hole)
    else
        boundary_npiece = 0
    end if
end function boundary_npiece

!-------------------------------------------------------------------------
!> This routine gives the range of parameter values for each piece of the boundary.
subroutine boundary_param(start,finish)
    !----------------------------------------------------
    ! This routine gives the range of parameter values for each piece of the
    ! boundary.
    ! If boundary_npiece <= 0 it is a dummy and the domain is given by triangle
    ! data files.
    !----------------------------------------------------
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    implicit none
    !----------------------------------------------------
    ! Dummy arguments
    real(my_real), intent(out) :: start(:), finish(:)
    !----------------------------------------------------
    ! Begin executable code
    
    ! Dummy version
    if(modnum .eq. 1) then 
        call example_boundary_param(start,finish)
    else
        start = 0.0_my_real
        finish = 0.0_my_real
    end if

end subroutine boundary_param


!-------------------------------------------------------------------------
!> This is the identity function that PHAML requires.
function phaml_integral_kernel(kernel,x,y)
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    integer, intent(in) :: kernel
    real(my_real), intent(in) :: x,y
    real(my_real) :: phaml_integral_kernel
    
    ! Identity function
    if(modnum .eq. 1) then 
        phaml_integral_kernel = example_phaml_integral_kernel(kernel,x,y)
    else
        phaml_integral_kernel = 1.0
    end if

end function phaml_integral_kernel

!-------------------------------------------------------------------------
!> Provides the \emph{a priori} knowledge about the singular nature of the solution if applicable.
function regularity(x,y)
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    real(my_real), intent(in) :: x(3),y(3)
    real(my_real) :: regularity
    
    ! Dummy version, assume infinitely differentiable everywhere.
    if(modnum .eq. 1) then 
        regularity = example_regularity(x,y)
    else
        regularity = huge(0.0_my_real)
    end if
end function regularity

!-------------------------------------------------------------------------
!> This routine updates the module variables on the slave processes by sending them from the master process.
    
subroutine update_usermod(phaml_solution) 
    !----------------------------------------------------
    ! This routine updates the module variables on the slave processes by
    ! sending them from the master process
    !----------------------------------------------------
    use phaml
    use phaml_example_pde
    use phaml_user_mod
    implicit none
    type(phaml_solution_type), intent(in) :: phaml_solution
    integer :: iparam(6)
    real(my_real),allocatable,dimension(:) :: rparam
    logical, save :: first_call = .true.

    !the slaves need the globals passed first so that the allocate
    !for rparam has the values for the size.  This means after
    !phaml_create another call to update_usermod must be made
    
    !this also means the modnum won't be known in the slaves so the first call
    !variables must be the same for all modules
    if (first_call) then
        allocate(rparam(1))
        iparam(1) = gnsn
        iparam(2) = gewn
        iparam(3) = gdns
        iparam(4) = gdew
        iparam(5) = num_arrays
        iparam(6) = modnum
        call master_to_slaves(phaml_solution,iparam,rparam)
        gnsn = iparam(1)
        gewn = iparam(2)
        gdns = iparam(3)
        gdew = iparam(4)
        num_arrays = iparam(5)
        modnum = iparam(6)
        deallocate(rparam)
        call array_init()
        first_call = .false.
    else
        if(modnum .eq. 1) then 
            call example_update_usermod(phaml_solution) 
        else
            allocate(rparam(gnsn*gewn*num_arrays))
            iparam(1) = gnsn
            iparam(2) = gewn
            iparam(3) = gdns
            iparam(4) = gdew
            iparam(5) = num_arrays
            iparam(6) = modnum
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
        end if
    end if
end subroutine update_usermod

