!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   felix_dycore_interface.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


module felix_dycore_interface

   use glimmer_physcon,  only : scyr
   use glimmer_paramets, only : vel0, tau0, vis0
   use glide_types
   use glimmer_log
   use parallel
   use glissade_grid_operators, only: glissade_stagger 
   !use glimmer_to_dycore

   implicit none
   private


   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: felix_velo_init, &
             felix_velo_driver

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************


contains


!***********************************************************************
!
!  routine felix_velo_init
!
!> \brief   Initializes the external Albany-FELIX velocity solver
!> \author  Irina Kalashnikova
!> \date    13 September 2013
!> \version SVN:$Id$
!> \details
!>  This routine initializes the external Albany-FELIX ice velocity solver.
!
!-----------------------------------------------------------------------

   subroutine felix_velo_init(model)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(glide_global_type),intent(inout) :: model

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------


      if (this_rank == 0) print *, 'DEBUG: Inside felix_velo_init.'

      ! === First do any preparations needed on the CISM side (if any)


      ! === Now call the external Albany code for any init that it needs to do
      !call gtd_init_dycore(model,dycore_model_index)
      ! Doug - does this interface still make sense here?
      ! Doug - what needs to change (if anything) if the code is compiled without
      !        external Felix libraries?  Do we need a stub module?
      ! Doug - We might need to do some rearranging to make sure the call 
      !        to gtd_init_dycore_interface happens in the right place 
      !        (presumably in simple_glide/simple_felix/cism_driver).
      !        (I think I see how to do this, but will wait for now.)

   !--------------------------------------------------------------------
   end subroutine felix_velo_init




!***********************************************************************
!
!  routine felix_velo_driver
!
!> \brief   Makes preparations and calls the external Albany-FELIX velocity solver
!> \author  Irina Kalashnikova
!> \date    13 September 2013
!> \version SVN:$Id$
!> \details
!>  This routine makes preparations and calls the external
!>  Albany-FELIX velocity solver.
!
!-----------------------------------------------------------------------

   subroutine felix_velo_driver(model)

      use glimmer_global, only : dp
      use glimmer_physcon, only: gn, scyr
      use glimmer_paramets, only: thk0, len0, vel0, vis0
      use glimmer_log
      use glide_types
      use glide_mask
      
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      
      type(glide_global_type),intent(inout) :: model

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      

      call get_parallel_finite_element_mesh_data(model%general%ewn,      model%general%nsn ,&
                                                 model%general%upn, &
                                                 model%numerics%sigma, &
                                                 nhalo, &
                                                 len0 * model%numerics%dew, &
                                                 len0 * model%numerics%dns, &
                                                 thk0 * model%geometry%thck, &
                                                 thk0 * model%geometry%usrf, &
                                                 thk0 * model%geometry%topg,&
                                                 thk0 * model%numerics%thklim, &
                                                 (tau0 / vel0 / scyr) *model%velocity%beta, &
                                                 (vis0*scyr) *model%temper%flwa)


      !IK, 10/24/13, notes to self:
      !To use constant flwa = 1e-16, set flow_law = 0 in input (config) file
      !To use beta field from .nc file, set which_ho_babc = 5 in input (config)
      !file; to use no-slip, set which_ho_babc = 4
     

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------


      if (this_rank == 0) print *, 'DEBUG: Inside felix_velo_driver.'

      ! === First do any preparations needed on the CISM side


      ! === Now call the external Albany code
      !call gtd_run_dycore(dycore_model_index,cur_time,time_inc)
      ! Doug - does this interface still make sense here?
      ! Doug - what needs to change (if anything) if the code is compiled without
      !        external Felix libraries?  Do we need a stub module?


   !--------------------------------------------------------------------
   end subroutine felix_velo_driver



!***********************************************************************
! Private subroutines:
!***********************************************************************



!***********************************************************************
!
!  routine get_parallel_finite_element_mesh_data
!
!> \author  Irina Kalashnikova
!> \date    18 September 2013
!> \version SVN:$Id$
!> \details
!
! Naming convention:
! - cells, vertices are in 2D
! - elements, nodes are in 3D
!
! The function get_parallel_finite_element_mesh_data creates a parallel mesh of
! a given geometry using the data.  In particular, global node and element IDs
! are created, and an offset it added to the x and y coordinates on
! multi-processor runs.  The following are data that would be needed in
! Albany/FELIX (so these would need to be passed through an interface b/w the 2
! codes): 
!
! xyz_at_nodes:                       
! Double array of size (nx-1)*(ny-1)*nz x 3.  It gives the x, y
! and z coordinates of all the nodes on each  processor.  
! Note: Right now this array consists of the full mesh, in
! particular, non-active nodes have not been removed.  This is OK for
! Albany/FELIX -- the non-active nodes will not be assembled as they will not
! appear in global_element_conn_active, the element connectivity array.  We could
! remove the non-active nodes at some point to avoid passing stuff b/w the
! codes that isn't needed.  
! Note 2: the nodes need to be converted to km prior to being
! passed to Albany/FELIX b/c Albany/FELIX works with meshes in km. 
!
! global_node_id_owned_map: 
! Integer array of size (nx-1)*(ny-1)*nz x 1.  This is
! effectively a map from local node IDs to global node
! IDs.  It is 1-based, consistent w/ Fortran numbering (so the first node
! is node number 1).
!
! global_element_conn_active: 
! Dynamically allocated integer array of size nCellsActive*(nz-1) x 8 where nCellsActive is
! the number of active elements (with ice) in 2D.  This array is the element
! connectivity array.  The 8 columns of this array give the element
! connectivity (node #s defining a given element), 1-based.  
! Note: The global element numbering in global_element_conn_active will be
! non-contiguous and there will be some element #s
! missing (e.g., if elements 1 and 2 are not active, they will not appear in
! global_element_conn_active).  This is OK for Albany/FELIX.  Also some of the nodes
! (the non-active ones) will not appear in the connectivity array.  This
! is OK too. 
! 
! global_element_id_active_owned_map: 
! Dynamically allocated integer array of size
! nCellsActive*(nz-1) x 1 where nCellsActive is the number of active
! elements (with ice) in 2D.  This is a map from local element IDs to global element IDs.
! It is 1-based.  Only active elements are included. 
!
! global_basal_face_conn_active: 
! Dynamically allocated integer array of size nCellsActive x 5 where nCellsActive is the number of
! active elements (with ice) in 2D.  This array is the basal face connectivity
! array.  The first column gives the global number of the element to
! which the face belongs.  The next 4 columns give the face connectivity (node #s
! defining the face of the element), again 1-based.
! Note: Same comment as for global_element_conn_active.
!
! global_basal_face_id_active_owned_map: 
! Dynamically allocated integer array of size
! nCellsActive x 1 where nCellsActive is the number of
! active elements (with ice) in 2D.  This is a map from local
! face IDs to global face IDs.  It is 1-based.  Only active
! faces are included.
!
! surf_height_at_nodes: 
! Double array of size (nx-1)*(ny-1)*nz.  This is effectively
! stagusrf extended to 3D (we need it as a 3D data structure in Albany/FELIX). 
! Note: Like the xyz_at_nodes array this would be defined at all the
! nodes in the original mesh, so it would include non-active
! nodes.  We can change this at some point if we don't want to pass extra stuff b/w
! codes.  Note also that this needs to be converted to km as Albany/FELIX uses meshes in km not
! meters. 
!
!
!-----------------------------------------------------------------------

   subroutine get_parallel_finite_element_mesh_data(nx,         ny,           &
                                                    nz,         sigma,        &
                                                    nhalo,                    &
                                                    dx,         dy,           &
                                                    thck,       usrf,         &
                                                    topg,                     &
                                                    thklim,     beta,         &
                                                    flwa)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

       integer, intent(in) ::   &
       nx, ny,               &  ! number of grid cells in each direction
       nz,                   &  ! number of vertical levels where velocity is computed
                                ! (same as model%general%upn)
       nhalo                    ! number of rows/columns of halo cells

       real(dp), dimension(:), intent(in) :: &
       sigma

       real(dp), intent(in) ::  &
       dx,  dy                  ! grid cell length and width (m)
                                ! assumed to have the same value for each grid
                                ! cell

       real(dp), dimension(:,:), intent(in) ::  &
       thck,                 &  ! ice thickness (m)
       usrf,                 &  ! upper surface elevation (m)
       topg                     ! elevation of topography (m)


       real(dp), intent(in) ::   &
       thklim                   ! minimum ice thickness for active cells (m)

       real(dp), dimension(:,:), intent(in) ::  &
       beta                     ! basal traction parameter

       real(dp), dimension(:,:,:), intent(in) ::  &
       flwa                     ! flow factor parameter


      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !IK, 9/8/13: xyz_at_nodes will need to be passed to Albany/FELIX
      !These are divided by 1000 to convert from meters to km, as Albany/FELIX
      !takes meshes in km.  
      !TO DO: make xyz_at_nodes have intent(out) 
      real(dp), dimension((nx-2*nhalo+1)*(ny-2*nhalo+1)*nz,3) :: &
      xyz_at_nodes       ! x, y and z coordinates of each vertex

      !IK, 9/8/13: GlobalNodeID_3D will need to be passed to Albany/FELIX 
      !TO DO: make GlobalNodeID_3D  have intent(out) 
      integer, dimension((nx-2*nhalo+1)*(ny-2*nhalo+1)*nz) ::  &
      global_node_id_owned_map    !This is effectively a map from local -> global IDs
                                  !for the full 3D mesh        

      !IK, 9/8/13: surf_height_at_nodes will need to be passed to Albany/FELIX 
      !These values are divided by 1000 to convert from meters to km, as
      !Albany/FELIX
      !takes meshes in km.  
      !TO DO: make surf_height_at_nodes  have intent(out) 
      real(dp), dimension((nx-2*nhalo+1)*(ny-2*nhalo+1)*nz) ::  &
      surf_height_at_nodes                  !This is an extension of
                                            !stagusrf to 3D        

      !IK, 9/8/13: global_element_conn_active will need to be passed to
      !Albany/FELIX 
      !TO DO: make global_element_conn_active  have intent(out) 
      integer, dimension(:, :), allocatable ::   &
       global_element_conn_active         !Like global_element_conn but first column is
                                      !removed and only active elements (cells)
                                      !are included 


      !IK, 9/12/13: global_element_id_active_owned_map will need to be passed
      !to Albany/FELIX 
      !TO DO: make global_element_id_active_owned_map  have intent(out) 
      integer, dimension(:, :), allocatable ::   &
       global_element_id_active_owned_map         !First column of global_element_conn but and only 
                                               !active elements (cells) are
                                               !included.
                                               !This is effectively a map from
                                               !local -> global IDs for the
                                               !elements

      !TO DO: make global_basal_face_conn_active  have intent(out) 
      integer, dimension(:, :), allocatable ::   &
       global_basal_face_conn_active              !Like global_basal_face_conn but only active
                                      !elements (cells) are included 

      !IK, 9/12/13: global_basal_face_id_active_owned_map will need to be passed to
      !Albany/FELIX 
      !TO DO: make global_basal_face_id_active_owned_map  have intent(out) 
      integer, dimension(:, :), allocatable ::   &
       global_basal_face_id_active_owned_map         !First column of global_basal_face_conn but and only 
                                          !active elements (cells) are included.
                                          !This is effectively a map from
                                          !local -> global IDs for the
                                          !elements

      !IK, 10/24/13: beta_at_nodes will need to be passed to Albany/FELIX 
      !These values are divided by 1000 to convert from meters to kPa a m^(-1),
      !as
      !Albany/FELIX takes meshes in km (so beta needs to be converted to the
      !appropriate units).  
      !TO DO: make beta_at_nodes  have intent(out) 
      real(dp), dimension((nx-2*nhalo+1)*(ny-2*nhalo+1)*nz) ::  &
      beta_at_nodes                         !This is an extension of
                                            !beta to 3D        

      !IK, 10/24/13: flwa_at_active_cells will need to be passed to Albany/FELIX 
      !This is the value of the flow factor at the elements
      !These values are multilied by 1.0e12 to convert to Albany/FELIX units
      !TO DO: make flwa_at_active_cells  have intent(out) 
      real(dp), dimension(:, :), allocatable ::  &
      flwa_at_active_cells                  !This is essentially flwa in 
                                            !vector form and at only the active
                                            !cells      


      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      
      integer, dimension(nx,ny) ::    &
      imask                  ! = 1 where ice is present, else = 0

      logical, dimension(nx,ny) ::     &
      active_cell            ! true for active cells (thck > thklim and border locally owned vertices)


      real(dp), dimension(nx-1,ny-1) ::   &
      stagusrf,            & ! upper surface averaged to vertices
      stagthck               ! ice thickness averaged to vertices

      real(dp), dimension((nx-2*nhalo+1)*(ny-2*nhalo+1),2) :: &
      xy_at_vertices       ! x and y coordinates of each vertex

      integer, dimension((nx-2*nhalo+1)*(ny-2*nhalo+1)) ::  &
      global_vertex_id_owned_map !global IDs of 2D mesh 

      integer, dimension((nx-2*nhalo)*(ny-2*nhalo)) ::  &
      global_cell_id_owned_map   !global IDs of cells in 2D mesh 

      logical, dimension((nx-2*nhalo)*(ny-2*nhalo)) ::  &
      active_cell_vector          !This is like active_cell except in vector
                                  !form 

      logical, dimension((nx-2*nhalo)*(ny-2*nhalo)*(nz-1)) ::  &
      active_cell_vector3D           !This is an extension of
                                     !active_cell_vertex to 3D        

      integer, dimension((nx-2*nhalo)*(ny-2*nhalo)*(nz-1), 9) ::  &
      global_element_conn    !First column is effectively a map from local -> global IDs
                             !Remaining 8 columns give element connectivity (with
                             !global node #s)             

      integer, dimension((nx-2*nhalo)*(ny-2*nhalo), 6) ::  &
      global_basal_face_conn          !First column is effectively a map from local ->global IDs for basal faces
                                      !Second column gives global # of element to which this
                                      !boundary face belongs
                                      !Next 4 columns give the connectivity for the boundary
                                      !face   

      real(dp), dimension((nx-2*nhalo)*(ny-2*nhalo)*(nz-1)) ::  &
      flwa_at_cells                 !This is a vector form of flwa, defined at
                                    !all the elements


      integer :: i, j, k, l
      real(dp) :: x, y !x and y coordinates of vertices
      integer :: gnx, gny !for temporary calculation of global vertex/cell # in x  global vertex/cell # in y
      integer :: nNodes2D, nNodesProc2D !total # virtices, # vertices on this proc (in 2D)
      integer :: nEles2D, nElesProc2D, nElesProc3D !total # cells (in 2D), # cells on this proc (in 2D), # elements on this proc (in 3D)
      integer :: nodes_x !total # nodes in x 
      integer :: x_GID, y_GID, z_GID, x_GIDplus1, y_GIDplus1, z_GIDplus1, elem_GID, xy_plane !for creating element numbering 
      integer :: nCellsActive !# active cells (with ice) in 2D  


     !--------------------------------------------------------------------
     ! TO DO (IK, 9/18/13): 
     ! - Make stuff that needs to be passed to Albany/FELIX an out argument of
     ! this function 
     !--------------------------------------------------------------------

     !IK, 9/9/13: printing for debug 
     !print *, 'In glissade_velo_higher_data! IK'
     !print *, 'Proc #: ', this_rank
     !print *, 'nx: ', nx 
     !print *, 'ny: ', ny
     !print *, 'dx: ', dx 
     !print *, 'dy: ', dy
     !print *, 'ewlb: ', ewlb
     !print *, 'nslb: ', nslb
     !print *, 'global_ewn: ', global_ewn
     !print *, 'global_nsn:', global_nsn
     !print *, 'nhalo:', nhalo
     !print *, 'nz:', nz

     !---------------------------------------------------------------------------------------
     ! Creation of global node numbering of vertices/nodes (IK, 9/8/13)
     !---------------------------------------------------------------------------------------

     !IK, 9/8/13: first, create global vertices for 2D mesh to be extruded as 3D
     !mesh
     k = 1
     do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
         x = (ewlb+1)*dx + i*dx !xVertex(i,j) 
         y = (nslb+1)*dy + j*dy !yVertex(i,j) 
         gnx = ewlb + 1 + i - nhalo + 1
         gny = nslb + 1 + j - nhalo + 1
         global_vertex_id_owned_map(k) = gnx + (global_ewn+1)*(gny - 1)
         xy_at_vertices(k,1) = x/1000.0 !divide by 1000 to convert to km for
                                    !Albany/FELIX  
         xy_at_vertices(k,2) = y/1000.0 !divide by 1000 to convert to km for 
                                    !Albany/FELIX
         k = k + 1
       enddo
     enddo

     !IK, 9/8/13: now, create global nodes for 3D mesh obtained by extruding 3D
     !mesh in z-direction 
     !------------------------------------------------------------------------------
     ! Compute masks:
     ! mask = 1 where dynamically active ice is present, 0 elsewhere
     !------------------------------------------------------------------------------

     do j = 1, ny
        do i = 1, nx
           if (thck(i,j) > thklim) then
              imask(i,j) = 1
           else
              imask(i,j) = 0
           endif
        enddo
     enddo

     !------------------------------------------------------------------------------
     ! Compute ice thickness and upper surface on staggered grid
     ! (requires that thck and usrf are up to date in halo cells)
     !------------------------------------------------------------------------------

     call glissade_stagger(nx,      ny,         &
                           thck,    stagthck,   &
                           imask,   stagger_margin_in = 1)

     call glissade_stagger(nx,      ny,         &
                           usrf,    stagusrf,   &
                           imask,   stagger_margin_in = 1)

     !------------------------------------------------------------------------------

     nNodes2D = (global_ewn + 1)*(global_nsn + 1)
     nNodesProc2D = (nx - 2*nhalo+1)*(ny- 2*nhalo+1)
     do l = 1, nz  !loop over vertical layers
        global_node_id_owned_map((l-1)*nNodesProc2D + 1:nNodesProc2D*l) = global_vertex_id_owned_map + nNodes2D*(l - 1)
        xyz_at_nodes((l-1)*nNodesProc2D + 1:nNodesProc2D*l, 1:2) = xy_at_vertices
     enddo
     !IK, 9/8/13: set z-coordinate of mesh 
     k = 1
     do l = 1, nz
       do j = nhalo, ny-nhalo
         do i = nhalo, nx-nhalo
      !  do j = 1+nhalo, ny-nhalo+1
      !  do i = 1+nhalo, nx-nhalo+1
           !divide by 1000 to convert to km for Albany/FELIX
           xyz_at_nodes(k,3) = (stagusrf(i,j) - sigma(l)*stagthck(i,j))/1000.0
           surf_height_at_nodes(k) = stagusrf(i,j)/1000.0 
           beta_at_nodes(k) = beta(i,j)/1000.0;
           k = k + 1
         enddo
       enddo
     enddo

     !IK, 9/12/13: printing output for debugging/checking node
     !numbering/coordinates 
     if (this_rank == 0) then
       do l=1, (nx-2*nhalo+1)*(ny-2*nhalo+1)*nz
         print *, 'x, y, z: ',  xyz_at_nodes(l,1), xyz_at_nodes(l,2), xyz_at_nodes(l,3)
         print *, 'global node: ', global_node_id_owned_map(l)
         print *, 'sh: ', surf_height_at_nodes(l)
         print *, 'beta: ', beta_at_nodes(l)
       enddo
      endif

     ! Identify the active cells.
     ! Include all cells that border locally owned vertices and contain ice.

     nCellsActive = 0 !start counter keeping track of how many active cells
                      !there are on each processor  

     active_cell(:,:) = .false.

     do j = 1+nhalo, ny-nhalo
     do i = 1+nhalo, nx-nhalo
       if (thck(i,j) > thklim) then
          active_cell(i,j) = .true.
          nCellsActive = nCellsActive + 1
       endif
     enddo
     enddo

     !IK, 10/24/13: populate flwa_at_cells array from flwa array, and change
     !units to Albany/FELIX units
     k = 1
     do l = 1, nz-1
       do j = 1+nhalo, ny-nhalo
         do i = 1+nhalo, nx-nhalo
           flwa_at_cells(k) = flwa(l,i,j)*(1.0E12) !scale flwa by 1e12 to get units
                                                   !consistent with those in
                                                   !Albany/FELIX
           k = k + 1;
         enddo
       enddo
     enddo



     !--------------------------------------------------------------------------
     ! Creation of hexahedral mesh and global numbering of elements (IK, 9/8/13)
     !--------------------------------------------------------------------------

     nEles2D = global_ewn*global_nsn
     nElesProc2D = (nx - 2*nhalo)*(ny - 2*nhalo)
     k = 1 !local cell number
     do j = 1+nhalo, ny-nhalo
       do i = 1+nhalo, nx-nhalo
         gnx = ewlb + 1 + i - nhalo
         gny = nslb + 1 + j - nhalo
         global_cell_id_owned_map(k) = gnx + global_ewn*(gny - 1)
         active_cell_vector(k) = active_cell(i,j)
         k = k + 1
       enddo
     enddo
     do l = 1, nz - 1 !loop over vertical layers
       global_element_conn((l-1)*nElesProc2D + 1:nElesProc2D*l, 1) = global_cell_id_owned_map+ nEles2D*(l-1)
       active_cell_vector3D((l-1)*nElesProc2D + 1:nElesProc2D*l) = active_cell_vector
     enddo

     nodes_x = global_ewn + 1 !# nodes in x 
     nElesProc3D = (nx - 2*nhalo)*(ny - 2*nhalo)*(nz - 1) !number of elements on proc
     k = 1 ! counter for incrementing boundary faces  
     do i = 1, nElesProc3D
       elem_GID = global_element_conn(i, 1) - 1
       z_GID = elem_GID/nEles2D !mesh column number
       xy_plane = mod(elem_GID, nEles2D)
       x_GID = mod(xy_plane, global_ewn) !mesh column number 
       y_GID = xy_plane/(global_ewn) !mesh row number
       x_GIDplus1 = x_GID + 1
       y_GIDplus1 = y_GID + 1
       z_GIDplus1 = z_GID + 1
       ! find and mark boundary faces on basal boundary 
       if (z_GIDplus1 == nz - 1) then
         global_basal_face_conn(:, 1) = global_cell_id_owned_map
         global_basal_face_conn(k, 2) = global_element_conn(i,1)
         !IK, 9/8/13: below the +1 is added to make the connectivity 1-based
         !like in Fortran -- the node numbering has been created with this
         !convention    
         global_basal_face_conn(k, 3) = x_GID      + nodes_x*y_GID      + nNodes2D*z_GIDplus1 + 1
         global_basal_face_conn(k, 4) = x_GIDplus1 + nodes_x*y_GID      + nNodes2D*z_GIDplus1 + 1
         global_basal_face_conn(k, 5) = x_GIDplus1 + nodes_x*y_GIDplus1 + nNodes2D*z_GIDplus1 + 1
         global_basal_face_conn(k, 6) = x_GID      + nodes_x*y_GIDplus1 + nNodes2D*z_GIDplus1 + 1
         k = k + 1
       endif
       !IK, 9/8/13: below the +1 is added to make the connectivity 1-based
       !like in Fortran -- the node numbering has been created with this
       !convention    
       global_element_conn(i, 2) = x_GID      + nodes_x*y_GID      + nNodes2D*z_GIDplus1 + 1
       global_element_conn(i, 3) = x_GIDplus1 + nodes_x*y_GID      + nNodes2D*z_GIDplus1 + 1
       global_element_conn(i, 4) = x_GIDplus1 + nodes_x*y_GIDplus1 + nNodes2D*z_GIDplus1 + 1
       global_element_conn(i, 5) = x_GID      + nodes_x*y_GIDplus1 + nNodes2D*z_GIDplus1 + 1
       global_element_conn(i, 6) = x_GID      + nodes_x*y_GID      + nNodes2D*z_GID + 1
       global_element_conn(i, 7) = x_GIDplus1 + nodes_x*y_GID      + nNodes2D*z_GID + 1
       global_element_conn(i, 8) = x_GIDplus1 + nodes_x*y_GIDplus1 + nNodes2D*z_GID + 1
       global_element_conn(i, 9) = x_GID      + nodes_x*y_GIDplus1 + nNodes2D*z_GID + 1
     enddo

     !dynamically allocate arrays that depent on # active cells
     allocate(global_element_conn_active(nCellsActive*(nz-1), 8))
     allocate(global_element_id_active_owned_map(nCellsActive*(nz-1),1))
     allocate(global_basal_face_conn_active(nCellsActive, 5))
     allocate(global_basal_face_id_active_owned_map(nCellsActive,1))
     allocate(flwa_at_active_cells(nCellsActive*(nz-1), 1))
     !IK, 9/9/13: do dynamically-allocated arrays need to be deallocated/deleted?  
     k = 1
     do i = 1, nElesProc3D
       if (active_cell_vector3D(i)) then
         global_element_conn_active(k, 1:8) = global_element_conn(i,2:9)
         global_element_id_active_owned_map(k,1) = global_element_conn(i,1)
         flwa_at_active_cells(k,1) = flwa_at_cells(i)
         k = k + 1
       endif
     enddo
     k = 1
     do i = 1, nElesProc2D
       if (active_cell_vector(i)) then
         global_basal_face_conn_active(k, :) = global_basal_face_conn(i, 2:6)
         global_basal_face_id_active_owned_map(k,1) = global_basal_face_conn(i, 1)
         k = k + 1
       endif
     enddo

    !IK, 9/12/13: printing output for debugging/checking element numbering 
    if (this_rank == 0) then
     do l=1, nCellsActive*(nz-1)
       print *, 'element connectivity active: ', global_element_conn_active(l,1:8)
       print *, 'global element #: ', global_element_id_active_owned_map(l,1)
       print *, 'flwa: ', flwa_at_active_cells(l,1)
     enddo
     endif

    !IK, 9/12/13: printing output for debugging/checking basal face numbering 
    if (this_rank == 0) then
    do l=1, nCellsActive
      print *, 'face connectivity active: ', global_basal_face_conn_active(l,1:5)
      print *, 'global face #: ', global_basal_face_id_active_owned_map(l,1)
    enddo
    endif


   !--------------------------------------------------------------------
   end subroutine get_parallel_finite_element_mesh_data


end module felix_dycore_interface
