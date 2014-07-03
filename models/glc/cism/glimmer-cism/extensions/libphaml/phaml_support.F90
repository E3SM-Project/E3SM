!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   phaml_support.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module phaml_support
!    use phaml
!    use phaml_user_mod
!    implicit none
contains
    !------------------------------------------------------------------------------
    !SUBROUTINE: is_ice_edge
    !ARGUMENTS: model (glimmer), ew, ns, ret
    !DESCRIPTION:
    !   ew, ns is a point in the model. The function checks the four adjacent
    !   points on the mask.  if (ew,ns) has ice, but one of the adjacent points
    !   does not have ice, then the function returns true.  This is a way to find 
    !   a pseudo grounding line/ice edge.
    !------------------------------------------------------------------------------
    
    subroutine is_ice_edge(model,ew,ns,ret)
        use glide_types
        use glide_mask
        use phaml_user_mod
        implicit none
        
        type(glide_global_type) :: model
        integer, intent(in) :: ew,ns
        logical, intent(out) :: ret
        logical :: m,m1,m2,m3,m4
        
        !m (ew,ns)
        m = .false.
        if(has_ice(model%geometry%thkmask(ew,ns)) .eqv. .true.) then
            m = .true.
        end if
        
        !m1 (ew,ns+1)
        m1 = .true.
        if(ns .eq. gnsn) then
            m1 = .false.
        else if(has_ice(model%geometry%thkmask(ew,ns+1)) .eqv. .false.) then
            m1 = .false.
        end if
        
        !m2 (ew,ns-1)
        m2 = .true.
        if(ns .eq. 1) then
            m2 = .false.
        else if(has_ice(model%geometry%thkmask(ew,ns-1)) .eqv. .false.) then
            m2 = .false.
        end if
        
        !m3 (ew+1,ns)
        m3 = .true.
        if(ew .eq. gewn) then
            m3 = .false.
        else if(has_ice(model%geometry%thkmask(ew+1,ns)) .eqv. .false.) then
            m3 = .false.
        end if
        
        !m4 (ew-1,ns)        
        m4 = .true.
        if(ew .eq. 1) then
            m4 = .false.
        else if(has_ice(model%geometry%thkmask(ew-1,ns)) .eqv. .false.) then
            m4 = .false.
        end if
        ret = .false.
        if(m .eqv. .true.) then
            if((m1 .eqv. .false.) .or. &
                (m2 .eqv. .false.) .or. &
                (m3 .eqv. .false.) .or. &
                (m4 .eqv. .false.)) then
                ret = .true.    
            end if
        end if
        !is_ice_edge = ret
    end subroutine is_ice_edge
    !------------------------------------------------------------------------------
    !SUBROUTINE: get_bmark
    !ARGUMENTS: model (glimmer), ew1, ns1, ew2, ns2, bmark
    !DESCRIPTION:
    !   ew1, ns1 is a point in the model, and ew2,ns2 is an adjacent point such 
    !   that there is an edge between the two.  This function determines if that
    !   edge is on the grounding line or within the boundaries and sets the bmark
    !   to the appropriate value.
    !------------------------------------------------------------------------------
    subroutine get_bmark(model,ew1,ns1,ew2,ns2,bmark)
        use glide_types
        use glide_mask
        implicit none
        
        type(glide_global_type) :: model
        integer, intent(in) :: ns1,ew1,ns2,ew2
        integer, intent(out) :: bmark
        logical :: ret1, ret2
        ret1 = .false.
        ret2 = .false.
        bmark=0 ! default value
        if ((has_ice(model%geometry%thkmask(ew1,ns1)) .eqv. .true.) .and. &
            (has_ice(model%geometry%thkmask(ew2,ns2)) .eqv. .true.)) then 
            bmark=model%geometry%thkmask(ew1,ns1) 
        end if
        if ((is_grounding_line(model%geometry%thkmask(ew1,ns1)) .eqv. .true.) .and. &
            (is_grounding_line(model%geometry%thkmask(ew2,ns2)) .eqv. .true.)) then 
            bmark=model%geometry%thkmask(ew1,ns1)
        end if
        call is_ice_edge(model,ew1,ns1,ret1)
        call is_ice_edge(model,ew2,ns2,ret2)
        if (ret1 .and. ret2) then  
            bmark=glide_mask_grounding_line !model%geometry%thkmask(ew1,ns1)
        end if
    end subroutine get_bmark
    
    !------------------------------------------------------------------------------
    !SUBROUTINE: make_ice_poly_file
    !ARGUMENTS: model (glimmer)
    !DESCRIPTION:
    ! This subroutine creates a poly file for an ice sheet by using the CISM mask 
    ! variable to add only nodes which have ice.
    ! bmark is set to be the mask value at that node.
    !
    ! Running triangle is compiler dependent.  That call may need to be changed, 
    ! or may not be supported at all, on some compilers.
    !------------------------------------------------------------------------------
    subroutine make_ice_poly_file(model)
        use glide
        use glide_mask
        implicit none
        
        integer :: nx, ny, node1ew,node1ns,node2ew,node2ns
        real :: dx,dy   
        integer :: i, j, count,bmark,status,tot_nodes
        logical :: ret
        type(glide_global_type) :: model
        !node_table addressed x,y, stores node number  and is -1 if not included, also bmark
        integer, dimension(:,:,:), allocatable :: node_table
        integer, dimension(:,:), allocatable :: edge_table   
        integer, external :: system
        character string*2
    
        
        !grid variables
        nx = get_ewn(model)
        ny = get_nsn(model)
        dx = get_dew(model)
        dy = get_dns(model)!model%numerics%dns
        
        allocate(node_table(ny,nx,2))
        !declares the max number of edges needed if grid was square 1 is x, 2 is y
        count = max(nx,ny) 
        allocate(edge_table((4*count-2)*(count-1),3))
        
        !write a .poly file with the grid of nodes, boundary edges, and bmark
        
        open(unit=21,file="mesh.poly",status="replace")
        !count the nodes
        tot_nodes = 0
        count = 0
        do i=0,nx-1
            do j=0,ny-1
                bmark=0
                node1ew = i+1
                node1ns = j+1
                if (has_ice(model%geometry%thkmask(node1ew,node1ns)) .eqv. .true.) then
                    bmark=model%geometry%thkmask(node1ew,node1ns)
                end if
                if (is_grounding_line(model%geometry%thkmask(node1ew,node1ns)) .eqv. .true.) then 
                    bmark=model%geometry%thkmask(node1ew,node1ns)
                end if
                ret = .false.
                call is_ice_edge(model,node1ew,node1ns,ret)
                if (ret .eqv. .true.) then 
                    bmark=glide_mask_grounding_line!model%geometry%thkmask(node1ew,node1ns)!5
                end if
                if (bmark > 0) then
                    count = count + 1
                    node_table(node1ns, node1ew,1) = count
                    node_table(node1ns, node1ew,2) = bmark
                    tot_nodes = tot_nodes + 1
                else
                    node_table(node1ns, node1ew,1) = -1
                end if
            end do
        end do
        write(*,*) 'Writing Vertices'
        write(21,*) tot_nodes, 2, 0, 1
        do i=1,nx
            do j=1,ny
                if (node_table(j,i,1) .ge. 0) then
                    write(21,*) node_table(j,i,1), dx*(i-1), dy*(j-1), node_table(j,i,2)
                end if
            end do
        end do
        
        write(*,*) 'Writing Edges'
        count = 0
        !write edges connecting columns (vertical)
        do i=0,nx-1
            do j=0,ny-2
                bmark=0
                node1ns = j + 1
                node1ew = i + 1
                node2ns = j + 2
                node2ew = i + 1
                if ((node_table(node1ns,node1ew,1).ge. 0) .and. &
                    (node_table(node2ns,node2ew,1).ge. 0)) then
                    call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                    count = count + 1
                    edge_table(count,1) = node_table(node1ns,node1ew,1)
                    edge_table(count,2) = node_table(node2ns,node2ew,1)
                    edge_table(count,3) = bmark
                end if
            end do
        end do
        !write edges connecting rows (horizontal)
        do i=0,nx-2
            do j=0,ny-1
                bmark=0
                node1ns = j + 1
                node1ew = i + 1
                node2ns = j + 1
                node2ew = i + 2
                if ((node_table(node1ns,node1ew,1).ge. 0) .and. &
                    (node_table(node2ns,node2ew,1).ge. 0)) then
                    call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                    count = count + 1
                    edge_table(count,1) = node_table(node1ns,node1ew,1)
                    edge_table(count,2) = node_table(node2ns,node2ew,1)
                    edge_table(count,3) = bmark
                end if
            end do
        end do
        !write forward cross edges    
        do i=0,nx-2
            do j=0,ny-2
                bmark=0
                node1ns = j + 1
                node1ew = i + 1
                node2ns = j + 2
                node2ew = i + 2
                if ((node_table(node1ns,node1ew,1).ge. 0) .and. &
                    (node_table(node2ns,node2ew,1).ge. 0)) then
                    call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                    count = count + 1
                    edge_table(count,1) = node_table(node1ns,node1ew,1)
                    edge_table(count,2) = node_table(node2ns,node2ew,1)
                    edge_table(count,3) = bmark
                end if
            end do
        end do    
        !write backwards cross edges    
        do i=0,nx-2
            do j=0,ny-2
                bmark=0
                node1ns = j + 1
                node1ew = i + 2
                node2ns = j + 2
                node2ew = i + 1
                if ((node_table(node1ns,node1ew,1).ge. 0) .and. &
                    (node_table(node2ns,node2ew,1).ge. 0)) then
                    if ((node_table(j+1,i+1,1).lt. 0) .or. &
                        (node_table(j+2,i+2,1).lt. 0)) then
                        call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                        count = count + 1
                        edge_table(count,1) = node_table(node1ns,node1ew,1)
                        edge_table(count,2) = node_table(node2ns,node2ew,1)
                        edge_table(count,3) = bmark
                    end if
                end if
            end do
        end do
        write(21,*) count, 1
        do i=1,count
            write(21,*) i, edge_table(i,1), edge_table(i,2), edge_table(i,3)
        end do
        write(21,*) 0
        
        close(unit=21)
        
        !deallocate the dynamic arrays
        deallocate( node_table )
        deallocate( edge_table )
        ! run triangle
        ! NOTE: THIS IS COMPILER DEPENDENT.  You may need to change this statement.
        status = System("triangle -neqj mesh.poly")
    
    end subroutine make_ice_poly_file

    
    
    !------------------------------------------------------------------------------
    !SUBROUTINE: make_full_poly_file
    !ARGUMENTS: model (glimmer)
    !DESCRIPTION:
    ! This subroutine creates a node file for an nx X ny grid on the rectangular
    ! domain (xmin,xmax) X (ymin,ymax), and runs triangle to create the
    ! required .node, .ele, .neigh and .edge files.
    !
    ! bmark is set to be the mask value at that node
    !
    ! Running triangle is compiler dependent.  That call may need to be changed, 
    ! or may not be supported at all on some compilers.
    !------------------------------------------------------------------------------
    subroutine make_full_poly_file(model)
        use glide
        use glide_mask
        implicit none   
        integer :: nx, ny, node1ew,node1ns,node2ew,node2ns
        real :: dx,dy,xmin,xmax,ymin,ymax    
        integer :: i, j, count, bmark,status
        type(glide_global_type) :: model
        integer, dimension(:), allocatable :: xvals,yvals
        integer, dimension(:,:),pointer :: mask
        integer, external :: system
        
        !this is to change the grid type later.  not current used.
    
        !grid variables
        nx = get_ewn(model)
        ny = get_nsn(model)
        xmin = 0
        dx = get_dew(model)
        ymin = 0
        dy = get_dns(model)!model%numerics%dns
        !for testing
        xmax = nx*dx
        ymax = ny*dy
        
        allocate( xvals(nx) )
        allocate( yvals(ny) )
        
        do i=0,nx-1
           xvals(i+1) = dx*i
        end do
        
        do i=0,ny-1
           yvals(i+1) = dy*i
        end do
        
        ! write a .poly file with the grid of nodes, boundary edges, and bmark
        
        open(unit=21,file="mesh.poly",status="replace")
        
    
        write(*,*) 'Writing Vertices'
        write(21,*) (nx)*(ny), 2, 0, 1
        count = 0
        do i=0,nx-1
            do j=0,ny-1
                bmark=0
                node1ew = j+1
                node1ns = i+1
                bmark=model%geometry%thkmask(node1ew,node1ns)
                count = count + 1
                write(21,*) count, xvals(node1ns), yvals(node1ew), bmark
            end do
        end do
        
        write(*,*) 'Writing Edges'
        count = 0
        !number of edges = (nx-1)*ny + (ny-1)*nx + (ny-1)*(nx-1)
        write(21,*) (nx-1)*ny + (ny-1)*nx + (ny-1)*(nx-1), 1
        !write edges connecting columns (vertical)
        do i=0,nx-1
            do j=0,ny-2
                bmark=0
                node1ns = j + 1
                node1ew = i + 1
                node2ns = j + 2
                node2ew = i + 1
                call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                count = count + 1
                write(21,*) count,i*nx+j+1,i*nx+j+2,bmark
            end do
        end do
        !write edges connecting rows (horizontal)
        do i=0,nx-2
            do j=0,ny-1
                bmark=0
                node1ns = j + 1
                node1ew = i + 1
                node2ns = j + 1
                node2ew = i + 2
                call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                count = count + 1
                write(21,*) count,i*nx+j+1,(i+1)*nx+j+1,bmark
            end do
        end do
        !write cross edges    
        do i=0,nx-2
            do j=0,ny-2
                bmark=0
                node1ns = j + 1
                node1ew = i + 1
                node2ns = j + 2
                node2ew = i + 2
                call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                count = count + 1
                write(21,*) count,i*nx+j+1,(i+1)*nx+j+2,bmark
            end do
        end do    
        write(21,*) 0
        close(unit=21)
        
        !deallocate the dynamic arrays
        deallocate( xvals )
        deallocate( yvals )
        ! run triangle
        ! NOTE: THIS IS COMPILER DEPENDENT.  You may need to change this statement.
        status = System("triangle -neqj mesh.poly")
    
    end subroutine make_full_poly_file
!---------------------------------------------------------------------------
end module phaml_support





