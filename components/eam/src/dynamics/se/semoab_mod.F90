#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module semoab_mod    
#ifdef HAVE_MOAB
  use iso_c_binding
  use kinds, only : real_kind, iulog, long_kind, int_kind
!  use edge_mod, only : ghostbuffertr_t, initghostbufferTR, freeghostbuffertr, &
!       ghostVpack, ghostVunpack,  edgebuffer_t, initEdgebuffer

  use dimensions_mod, only: nelem, ne, np, nelemd, nlev
  use element_mod, only : element_t
  use parallel_mod, only : parallel_t

  use m_MergeSorts,     only: IndexSet, IndexSort
 
  use cam_grid_support, only:  iMap
  use cam_abortutils,   only : endrun
  use edgetype_mod,           only: edgedescriptor_t
  use gridgraph_mod,          only: gridvertex_t

  use seq_comm_mct,  only: MHID, MHFID !  app id on moab side, for homme moab coarse and fine mesh
  use seq_comm_mct,  only: MHPGID      !  app id on moab side, for PGx style mesh, uniform from se
  use seq_comm_mct,  only: atm_pg_active ! turn it on when PG style mesh active

  use dyn_grid,      only: fv_nphys ! phys grid mesh will be replicated too
  use gllfvremap_mod,         only:  gfr_f_get_corner_latlon

  use control_mod, only :  west, east, south, north  ! 1, 2, 3, 4
  implicit none

  save

  integer local_map(np,np) !  what is the index of gll point (i,j) in a local moabconn(start: start+(np-1)*(np-1)*4-1)
  integer, allocatable :: moabconn(:) ! will have the connectivity in terms of local index in verts
  
contains

  integer function search_in(intarr, leng, value)
     integer, intent(in) :: leng
     integer, intent(in) :: intarr(leng)
     integer, intent(in) :: value

     ! binary search, as the array is ordered
     integer current, left, right
     left = 1
     right = leng

     search_in = 0
     if ( (value .gt. intarr(leng) ) .or.  ( value .lt. intarr(1) ) ) goto 10

     if ( value .eq. intarr(1) ) then
        search_in = 1
        goto 10
     endif

     if ( value .eq. intarr(leng) ) then
        search_in = leng
        goto 10
     endif

     do while (left < right )
        current = (right+left)/2
        if ( intarr(current) .eq. value ) then
            search_in = current
            goto 10
        else if ( intarr(current) .lt. value ) then
            left = current
        else
            right = current
        endif
        if ( left .eq. right -1) goto 10
     enddo
10  continue
    if ( intarr(right) .eq. value) search_in = right

  end function search_in

  subroutine create_moab_meshes(par, elem)

    use ISO_C_BINDING
    use iMOAB, only: iMOAB_CreateVertices, iMOAB_WriteMesh, iMOAB_CreateElements, &
      iMOAB_ResolveSharedEntities, iMOAB_UpdateMeshInfo, iMOAB_DefineTagStorage, &
      iMOAB_SetIntTagStorage, iMOAB_ReduceTagsMax, iMOAB_GetIntTagStorage
    use coordinate_systems_mod, only :  cartesian3D_t,  spherical_to_cart, spherical_polar_t

    type (element_t), intent(inout) :: elem(:)
    type (parallel_t)      , intent(in) :: par

    integer ierr, i, j, ie, iv, block_ID, k, numvals
    integer icol, irow, je, linx ! local indices in fine el connect

    real(kind=real_kind), allocatable :: moab_vert_coords(:)

    integer moab_dim_cquads, ix, idx, nverts, nverts_c ! used for indexing in loops; nverts will have the number of local vertices

    integer nelemd2  ! do not confuse this with dimensions_mod::nelemd

    integer(kind=long_kind), dimension(:), allocatable :: gdofv
    !  this will be moab vertex handle locally
    integer, dimension(:), allocatable :: moabvh
    integer, dimension(:), allocatable :: indx  !  this will be ordered

    integer, dimension(:), allocatable :: vdone, elemids, vgids, gdofel
    integer, dimension(:), allocatable :: vdone_c, moabconn_c, moabvh_c
    integer  currentval, dimcoord, dimen, num_el, mbtype, nve

    character*100 outfile, wopts, localmeshfile, lnum, tagname, newtagg
    integer  tagtype, numco, tag_sto_len, ent_type, tagindex
    type (cartesian3D_t)             :: cart
    integer  igcol, ii, neigh

    integer nedges_c, nverts_pg, nelem_pg, edge_index, j1
    integer, dimension(:), allocatable :: local_cell_gids, indx_cell
    integer, dimension(:,:), allocatable :: elem_edge, edge
    integer, dimension(:), allocatable :: vdone_pg, moabconn_pg

    integer nat_edge_order(4)
    integer  internal_edges, boundary_edges, reverse_edges
    integer  edge_verts(4) ! local per coarse element ! nverts_c < edge_verts <= nverts_c + edge_index
    integer  middle_vertex ! nverts_c + edge_index  < middle_vertex <= verts_pg
    type (spherical_polar_t)  ::  current_2d_vertex
    logical  pos_edge ! when looping over edges , use gdof for marking !!

    nat_edge_order = (/south, east, north, west/)

    ! for np=4,
    !      28, 32, 36, 35
    !      25, 29, 33, 34
    ! j |  13, 17, 21, 22
    !      1,  5,  9,  10
    !(1,1)     i->

    ! character*100 outfile, wopts, localmeshfile, lnum, tagname
    ! integer  tagtype, numco, tag_sto_len, ent_type, tagindex
    do j=1,np-1
      do i =1, np-1
        ix = (j-1)*(np-1)+i-1
        local_map(i,j) = ix*4 + 1
      enddo
    enddo
    do j=1, np-1
      i = j
      local_map(np, j) = ((np-1)*j-1)*4 + 2
      local_map(i, np) = ( (np-1)*(np-2)+i-1)*4 + 4
    enddo
    local_map(np, np) = ((np-1)*(np-1)-1)*4 + 3

    nelemd2 = (nelemd)*(np-1)*(np-1)
    moab_dim_cquads = (nelemd)*4*(np-1)*(np-1)

    if(par%masterproc) then
      write (iulog, *) " MOAB: semoab_mod module: create_moab_mesh_fine;  on processor " , par%rank ," nelemd: " , nelemd
    endif

    if ( nelemd > 0 ) then
      allocate(gdofv(moab_dim_cquads))
      allocate(elemids(nelemd2))
    endif

    k=0 !   will be the index for element global dofs
    do ie=1,nelemd
      do j=1,np-1
        do i=1,np-1
          ix = (ie-1)*(np-1)*(np-1)+(j-1)*(np-1)+i-1
          gdofv(ix*4+1) = elem(ie)%gdofP(i,j)
          gdofv(ix*4+2) = elem(ie)%gdofP(i+1,j)
          gdofv(ix*4+3) = elem(ie)%gdofP(i+1,j+1)
          gdofv(ix*4+4) = elem(ie)%gdofP(i,j+1)
          elemids(ix+1) = (elem(ie)%GlobalId-1)*(np-1)*(np-1)+(j-1)*(np-1)+i
        enddo
      enddo
    enddo

!     order according to global dofs

    if ( nelemd > 0 ) then
      allocate(moabvh(moab_dim_cquads))
      allocate(indx(moab_dim_cquads))

      allocate(moabconn(moab_dim_cquads))
      call IndexSet(moab_dim_cquads, indx)
      call IndexSort(moab_dim_cquads, indx, gdofv, descend=.false.)
!     after sort, gdofv( indx(i)) < gdofv( indx(i+1) )
    endif 
    idx=0
    currentval = 0
    if ( nelemd > 0 ) then
      currentval = gdofv( indx(1))
      idx = 1
    endif
     
    do ix=1,moab_dim_cquads
      if (gdofv(indx(ix)) .ne. currentval ) then
        idx=idx+1
        currentval = gdofv(indx(ix))
      endif
      moabvh(ix) = idx  ! the vertex in connectivity array will be at this local index
      ! this will be the moab connectivity
      moabconn(indx(ix)) = idx
    enddo

    nverts = idx
    if(par%masterproc) then
      write (iulog, *) " MOAB: there are ", nverts, " local vertices on master task ", currentval, " is the max local gdof"
    endif
    if ( nelemd > 0 ) then
      allocate(moab_vert_coords(3*nverts) )
      allocate(vdone(nverts))
      vdone = 0;
    endif
    if ( nelemd > 0 ) currentval = gdofv( indx(1)) ! start over to identify coordinates of the vertices

    do ix=1,moab_dim_cquads
      idx = indx(ix)   ! index in initial array, vertices in all fine quads
      k = (idx-1)/(4*(np-1)*(np-1))  ! index of coarse quad, locally, starts at 0
      ie = 1 + k  ! this is the element number; starts at nets=1
      je = ( idx -1 -k*(np-1)*(np-1)*4 ) / 4 + 1 ! local fine quad in coarse, 1 to (np-1) ^ 2
      irow = (je-1)/(np-1)+1
      icol = je - (np-1)*(irow-1)
      linx = idx - k*(np-1)*(np-1)*4 -(je-1)*4  ! this should be 1, 2, 3, 4
      if( linx == 1) then
        j = irow
        i = icol
      else if (linx == 2) then
        j = irow
        i = icol + 1
      else if (linx == 3) then
        j = irow + 1
        i = icol + 1
      else ! linx == 4
        j = irow + 1
        i = icol
      endif

      iv = moabvh(ix)
      if (vdone(iv) .eq. 0) then
        cart = spherical_to_cart (elem(ie)%spherep(i,j) )
        moab_vert_coords ( 3*(iv-1)+1 ) = cart%x
        moab_vert_coords ( 3*(iv-1)+2 ) = cart%y
        moab_vert_coords ( 3*(iv-1)+3 ) = cart%z
        vdone(iv) = gdofv(indx(ix)) ! this will be now our tag used for resolving shared entities ! convert to int, from long int
      endif

     enddo

     dimcoord = nverts*3
     dimen = 3
     if ( nelemd > 0 ) then
       ierr = iMOAB_CreateVertices(MHFID, dimcoord, dimen, moab_vert_coords)
       if (ierr > 0 )  &
         call endrun('Error: fail to create MOAB vertices ')
     endif
     !!num_el = nelemd2
     mbtype = 3 !  quadrilateral
     nve = 4;
     block_ID = 200 ! this will be for coarse mesh
     
     if ( nelemd > 0 ) then
       ierr = iMOAB_CreateElements( MHFID, nelemd2, mbtype, nve, moabconn, block_ID );
       if (ierr > 0 )  &
         call endrun('Error: fail to create MOAB elements')
     endif
      ! nverts: num vertices; vdone will store now the markers used in global resolve
      ! for this particular problem, markers are the global dofs at corner nodes
! set the global id for vertices
!   first, retrieve the tag
     tagname='GDOF'//C_NULL_CHAR
     tagtype = 0  ! dense, integer
     numco = 1
     ierr = iMOAB_DefineTagStorage(MHFID, tagname, tagtype, numco,  tagindex )
     if (ierr > 0 )  &
       call endrun('Error: fail to retrieve global id tag')
     ! now set the values
     ent_type = 0 ! vertex type
     if ( nverts > 0 ) then
       ierr = iMOAB_SetIntTagStorage ( MHFID, tagname, nverts , ent_type, vdone)
       if (ierr > 0 )  &
         call endrun('Error: fail to set marker id tag for vertices')
     endif 

     ierr = iMOAB_ResolveSharedEntities( MHFID, nverts, vdone );
     if (ierr > 0 )  &
       call endrun('Error: fail to resolve shared entities')

     if ( nelemd > 0) then
       vdone = -1 !  reuse vdone for the new tag, GLOBAL_ID (actual tag that we want to store global dof )
     endif
! use element offset for actual global dofs
     ! tagtype = 0  ! dense, integer
     ! numco = 1
     newtagg='GLOBAL_ID'//C_NULL_CHAR
     ierr = iMOAB_DefineTagStorage(MHFID, newtagg, tagtype, numco,  tagindex )
     if (ierr > 0 )  &
       call endrun('Error: fail to create new GDOF tag')
     do ie=1,nelemd
       do ii=1,elem(ie)%idxp%NumUniquePts
         i=elem(ie)%idxp%ia(ii)
         j=elem(ie)%idxp%ja(ii)
         igcol = elem(ie)%idxp%UniquePtoffset+ii-1
         ix = local_map(i,j)
         idx = moabconn((ie-1)*(np-1)*(np-1)*4 + ix) ! should
         vdone ( idx ) = igcol
       end do
     end do
     ! now set the values
     ent_type = 0 ! vertex type
     if ( nverts > 0 ) then
       ierr = iMOAB_SetIntTagStorage ( MHFID, newtagg, nverts , ent_type, vdone)
       if (ierr > 0 )  &
         call endrun('Error: fail to set global dof tag for vertices')
     endif

     ierr = iMOAB_ReduceTagsMax ( MHFID, tagindex, ent_type)
     if (ierr > 0 )  &
       call endrun('Error: fail to reduce max tag')

     ! set global id tag for elements
     ent_type = 1 ! now set the global id tag on elements
     if ( nelemd2 > 0 ) then
       ierr = iMOAB_SetIntTagStorage ( MHFID, newtagg, nelemd2 , ent_type, elemids)
       if (ierr > 0 )  &
         call endrun('Error: fail to set global id tag for elements')
     endif

!  now, after reduction, we can get the actual global ids for each vertex in the fine mesh
!  before, some vertices that were owned in MOAB but not owned in CAM did not have the right global ID tag
!  so vdone will be now correct on every task (no -1 anymore )
     ent_type = 0 ! vertex type
     if ( nverts > 0 ) then
       allocate(vgids(nverts))
       ierr = iMOAB_GetIntTagStorage ( MHFID, newtagg, nverts , ent_type, vgids)
       if (ierr > 0 )  &
         call endrun('Error: fail to retrieve GLOBAL ID on each task')
     endif
     ierr = iMOAB_UpdateMeshInfo(MHFID)
     if (ierr > 0 )  &
       call endrun('Error: fail to update mesh info')
#ifdef MOABDEBUG
!     write out the mesh file to disk, in parallel
     outfile = 'wholeFineATM.h5m'//C_NULL_CHAR
     wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
     ierr = iMOAB_WriteMesh(MHFID, outfile, wopts)
     if (ierr > 0 )  &
       call endrun('Error: fail to write the mesh file')
#endif


!    now create the coarse mesh, but the global dofs will come from fine mesh, after solving
    ! nelemd2 = nelemd
    moab_dim_cquads = (nelemd)*4

     if ( nelemd > 0 ) then
       allocate(gdofel(nelemd*np*np))
     endif
     k=0 !   will be the index for element global dofs
     do ie=1,nelemd
       ix = ie-1
       !
       gdofv(ix*4+1) = elem(ie)%gdofP(1,1)
       gdofv(ix*4+2) = elem(ie)%gdofP(np,1)
       gdofv(ix*4+3) = elem(ie)%gdofP(np,np)
       gdofv(ix*4+4) = elem(ie)%gdofP(1,np)
       elemids(ix+1) = elem(ie)%GlobalId
     enddo
! now original order

!     order according to global dofs
!     allocate(indx(moab_dim_cquads))
     if ( nelemd > 0 ) then
       call IndexSet(moab_dim_cquads, indx)
       call IndexSort(moab_dim_cquads, indx, gdofv, descend=.false.)
!      after sort, gdofv( indx(i)) < gdofv( indx(i+1) )

       allocate(moabvh_c(moab_dim_cquads))

       allocate(moabconn_c(moab_dim_cquads))
     endif
     idx = 0
     if ( nelemd > 0 ) then
       idx=1
       currentval = gdofv( indx(1))
     endif
     do ix=1,moab_dim_cquads
        if (gdofv(indx(ix)) .ne. currentval ) then
          idx=idx+1
          currentval = gdofv(indx(ix))
        endif
        moabvh_c(ix) = idx  ! the vertex in connectivity array will be at this local index
        ! this will be the moab connectivity
        moabconn_c(indx(ix)) = idx
     enddo
     nverts_c = idx
     if(par%masterproc) then
       write (iulog, *) " MOAB: there are ", nverts_c, " local vertices on master task, coarse mesh"
     endif
!     allocate(moab_vert_coords(3*idx) )
     if ( nelemd > 0 ) then 
       allocate(vdone_c(nverts_c))
       vdone_c = 0;
       currentval = gdofv( indx(1)) ! start over to identify coordinates of the vertices
     endif

     do ix=1,moab_dim_cquads
        i = indx(ix)   ! index in initial array
        ie = 1+ (i-1)/4 ! this is the element number
        j = i - ( i-1)/4*4 ! local index of vertex in element i
        iv = moabvh_c(ix)
        if (vdone_c(iv) .eq. 0) then
            moab_vert_coords ( 3*(iv-1)+1 ) = elem(ie)%corners3d(j)%x
            moab_vert_coords ( 3*(iv-1)+2 ) = elem(ie)%corners3d(j)%y
            moab_vert_coords ( 3*(iv-1)+3 ) = elem(ie)%corners3d(j)%z
            vdone_c(iv) = gdofv(indx(ix)) ! this will be now our tag used for resolving shared entities
        endif

     enddo

     dimcoord = nverts_c*3
     dimen = 3
     if ( nverts_c > 0 ) then
       ierr = iMOAB_CreateVertices(MHID, dimcoord, dimen, moab_vert_coords)
       if (ierr > 0 )  &
         call endrun('Error: fail to create MOAB vertices ')
     endif
     ! num_el = nelemd
     mbtype = 3 !  quadrilateral
     nve = 4;
     block_ID = 100 ! this will be for coarse mesh

     if ( nelemd > 0 ) then 
       ierr = iMOAB_CreateElements( MHID, nelemd, mbtype, nve, moabconn_c, block_ID );
       if (ierr > 0 )  &
         call endrun('Error: fail to create MOAB elements')
     endif
      ! idx: num vertices; vdone will store now the markers used in global resolve
      ! for this particular problem, markers are the global dofs at corner nodes
! set the global id for vertices
!   first, retrieve the tag
     tagname='GDOFV'//C_NULL_CHAR
     tagtype = 0  ! dense, integer
     numco = 1
     ierr = iMOAB_DefineTagStorage(MHID, tagname, tagtype, numco,  tagindex )
     if (ierr > 0 )  &
       call endrun('Error: fail to retrieve GDOFV id tag')
     ierr = iMOAB_DefineTagStorage(MHID, newtagg, tagtype, numco,  tagindex )
     if (ierr > 0 )  &
       call endrun('Error: fail to retrieve GLOBAL_ID tag on coarse mesh')
     ! now set the values
     ent_type = 0 ! vertex type
     if ( nverts_c > 0 ) then
       ierr = iMOAB_SetIntTagStorage ( MHID, tagname, nverts_c , ent_type, vdone_c)
       if (ierr > 0 )  &
         call endrun('Error: fail to set GDOFV tag for vertices')
     endif
      ! set global id tag for coarse elements, too; they will start at nets=1, end at nete=nelemd
     ent_type = 1 ! now set the global id tag on elements
     if ( nelemd > 0 ) then 
       ierr = iMOAB_SetIntTagStorage ( MHID, newtagg, nelemd , ent_type, elemids)
       if (ierr > 0 )  &
         call endrun('Error: fail to set global id tag for vertices')
     endif 

     ierr = iMOAB_ResolveSharedEntities( MHID, idx, vdone_c );
     if (ierr > 0 )  &
       call endrun('Error: fail to resolve shared entities')

!   global dofs are the GLL points are set for each element
     tagname='GLOBAL_DOFS'//C_NULL_CHAR
     tagtype = 0  ! dense, integer
     numco = np*np !  usually, it is 16; each element will have the dofs in order
     ierr = iMOAB_DefineTagStorage(MHID, tagname, tagtype, numco,  tagindex )
     if (ierr > 0 )  &
       call endrun('Error: fail to create global DOFS tag')
     ! now set the values
     ! set global dofs tag for coarse elements, too; they will start at nets=1, end at nete=nelemd
     ent_type = 1 ! now set the global id tag on elements
     numvals = nelemd*np*np ! input is the total number of values
     ! form gdofel from vgids
     do ie=1, nelemd
       ix = (ie-1)*np*np ! ie: index in coarse element
       je = (ie-1) * 4 * (np-1) * (np -1) !  index in moabconn array
       ! vgids are global ids for fine vertices (1,nverts)
       iv = 1
       do j=1,np
         do i=1,np
           k = local_map(i,j)
           gdofel(ix+iv) = vgids( moabconn( je + k ) )
           iv = iv + 1
         enddo
       enddo
       !  extract global ids
       vdone_c( moabconn_c( (ie-1)*4+1) ) = vgids ( moabconn(je+1 ))
       vdone_c( moabconn_c( (ie-1)*4+2) ) = vgids ( moabconn(je+ 4*(np-2)+2 )) ! valid for np = 4, 10
       vdone_c( moabconn_c( (ie-1)*4+3) ) = vgids ( moabconn(je+ 4*((np-1)*(np-1)-1) + 3 )) ! for np = 4,  35
       vdone_c( moabconn_c( (ie-1)*4+4) ) = vgids ( moabconn(je+ 4*(np-2)*(np-1) + 4 ))  !  28 for np = 4
     enddo
     if ( nelemd > 0 ) then
       ierr = iMOAB_SetIntTagStorage ( MHID, tagname, numvals, ent_type, gdofel)
       if (ierr > 0 )  &
         call endrun('Error: fail to set globalDOFs tag for coarse elements')
     endif 

      ! set the global ids for coarse vertices the same as corresponding fine vertices
     ent_type = 0 ! vertex type
     if ( nverts_c > 0 ) then
       ierr = iMOAB_SetIntTagStorage ( MHID, newtagg, nverts_c , ent_type, vdone_c)
       if (ierr > 0 )  &
         call endrun('Error: fail to set GLOBAL_DOFS tag values')
     endif 

     ierr = iMOAB_UpdateMeshInfo(MHID)
     if (ierr > 0 )  &
       call endrun('Error: fail to update mesh info')
#ifdef MOABDEBUG
!    write out the mesh file to disk, in parallel
     outfile = 'wholeATM.h5m'//C_NULL_CHAR
     wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
     ierr = iMOAB_WriteMesh(MHID, outfile, wopts)
     if (ierr > 0 )  &
       call endrun('Error: fail to write the mesh file')
#endif

     if (fv_nphys > 0 ) then
       ! create FV mesh, base on PGx
       atm_pg_active = .true. ! from now on, we will migrate / send tags for FV / ATM_PG2 mesh
       ! first count the number of edges in the coarse mesh;
       ! use euler: v-m+f = 2 => m = v + f - 2
       nedges_c = nverts_c + nelemd - 1 ! could be more, if unconnected regions ?
       if ( nedges_c < 0 ) nedges_c = 0 ! it cannot be negative
       internal_edges = 0
       boundary_edges = 0
       reverse_edges = 0
       nelem_pg =  fv_nphys * fv_nphys * nelemd ! each coarse cell is divided in fv_nphys x fv_nphys subcells
       !
       ! there are new vertices on each coarse edge (fv_phys - 1) , and (fv_nphys - 1) * (fv_nphys - 1)
       !   new vertices on each coarse cell
       if ( nelemd > 0 ) then 
         allocate (local_cell_gids(nelemd))
         allocate (indx_cell(nelemd))
         allocate (edge(2,nedges_c)) !
       endif
       do ie=1, nelemd !
         local_cell_gids(ie) = elem(ie)%GlobalID
       enddo
       if ( nelemd > 0 ) then
         call IndexSet(nelemd, indx_cell)
         call IndexSort(nelemd, indx_cell, local_cell_gids, descend=.false.)
       ! print *, ' local_cell_gids ', local_cell_gids
       ! print *, ' indx_cell ', indx_cell
         allocate( elem_edge (4, nelemd) )
       ! print *, '------------------------------- '
       ! print *, "RANK:", par%rank
       endif
       edge_index = 0
       do ie=1, nelemd !
           ! we need to check if neighbor is with id smaller; that means it was already created ?
           ! print *, '-------------- '
           ! print *, ' elem ', ie, elem(ie)%desc%actual_neigh_edges, elem(ie)%vertex%number, elem(ie)%GlobalID
           ! print *, '   nodes ', ( moabconn_c( (ie-1)*4+j), j=1,4 )
           ! print *, '   ids ', (vdone_c( moabconn_c( (ie-1)*4+j) ), j=1,4)
           ! print *, '   neigh: ', (elem(ie)%desc%globalID(j), j=1,4)
           ! print *, '   neigh order ', ( elem(ie)%desc%globalID(nat_edge_order(j)),j = 1,4 )
           k = elem(ie)%GlobalID ! current id
           do j = 1,4
               ix = j+1
               if (ix .eq. 5) ix = 1 ! next vertex in connectivity array
               neigh = elem(ie)%desc%globalID(nat_edge_order(j)) ! id neighbor
               idx = search_in(local_cell_gids, nelemd, neigh) ! index in local cells
               ! print *, '   ', j, 'neigh:', neigh, ' index' , idx

               if ( idx .gt. 0 ) then
                   ! a local edge is interior

                   if (k < neigh) then ! form the edge, increment edge index
                       edge_index = edge_index + 1
                       edge(1, edge_index) = moabconn_c(4*(ie-1) + j) ! first vertex
                       edge(2, edge_index) = moabconn_c(4*(ie-1) + ix) ! second vertex index
                       elem_edge(j, ie) = edge_index
                       internal_edges = internal_edges + 1
                       ! print *, ' edge:', edge_index, edge(1, edge_index), edge(2, edge_index), 'verts:' , &
                       !   vdone_c(edge(1, edge_index)), vdone_c(edge(2, edge_index)), ' element ', ie, ' intedge:', internal_edges

                   else
                       ! find the edge in the other list elem(idx)%globalID( nat_edge_order(j) )
                       do j1 = 1,4
                           if ( elem(idx)%desc%globalID( nat_edge_order(j1) ) .eq. k ) then
                               elem_edge(j, ie) = - elem_edge(j1, idx) ! inverse oriented
                               reverse_edges = reverse_edges + 1
                               ! print *, ' negative edge: ', elem_edge(j, ie), edge(1, -elem_edge(j, ie)), edge(2, -elem_edge(j, ie)), &
                               ! 'verts:', vdone_c(edge(1, -elem_edge(j, ie))), vdone_c(edge(2, -elem_edge(j, ie))), 'indx neg', reverse_edges
                           endif
                       enddo

                   endif
               else ! idx is 0, so it means the edge is on the boundary, form it
                   edge_index = edge_index + 1
                   edge(1, edge_index) = moabconn_c(4*(ie-1) + j) ! first vertex
                   edge(2, edge_index) = moabconn_c(4*(ie-1) + ix) ! second vertex index
                   elem_edge(j, ie) = edge_index
                   boundary_edges = boundary_edges + 1
                   ! print *, ' edge:', edge_index, edge(1, edge_index), edge(2, edge_index), 'verts:' , &
                   !       vdone_c(edge(1, edge_index)), vdone_c(edge(2, edge_index)), ' element ', ie, &
                   !       ' bedge:', boundary_edges
               endif
           enddo
       enddo
       ! show off
       nverts_pg = nverts_c + (fv_nphys - 1) * edge_index + (fv_nphys - 1) * (fv_nphys - 1) * nelemd
       ! print *, " MOAB: there are ", nverts_pg, " local vertices on master task on pg mesh ", edge_index , " local coarse edges ", &
       !  boundary_edges , ' boundary edges '
       if(par%masterproc) then
         write (iulog, *) " MOAB: there are ", nverts_pg, " local vertices on master task on pg mesh ", edge_index , " local coarse edges "
       endif
       !print *, '\n ELEMENTS: '
       !do ie=1,nelemd
       !    print *, ie, elem(ie)%GlobalID, ' local nodes:', ( moabconn_c( (ie-1)*4+j), j=1,4 ), ' edges:', (elem_edge(j, ie), j=1,4)
       !enddo
       !print *, '\n EDGES:'
       !do ie=1,edge_index
       !    print *, ie, (edge(j, ie), j=1,2)
       !enddo
       ! now generate phys grid, uniform FV type mesh;
       ! 2 cases: fv_nphys is 1 or 2; when 2, we need new nodes; will use the same id as
       ! the gdof on edge is used, with the smaller id chosen, among
       if ( nelemd > 0 ) then
         allocate(moabconn_pg(4*nelem_pg)) ! connectivity
         ! reuse moab_vert_coords for coordinates of pg mesh
         ! the first nverts_c coords are the same as coarse mesh; but we do create new
         allocate(vdone_pg(nverts_pg))
       endif
       do iv = 1, nverts_c
          vdone_pg(iv) = vdone_c(iv) ! also the coordinates will be the same !!
       enddo

       ! copy the coordinates from the middle
       j1 = 0 ! index in edge vertices; increase only for positive edges
       !  still need some
       if (fv_nphys .eq. 2) then
          current_2d_vertex%r = 1.
          do ie = 1,nelemd
              ix = (ie-1)*np*np ! ie: index in coarse element
              do j=1,4
                 idx = elem_edge(j, ie) !
                 if (idx .gt. 0) then ! increment edges, add vertex !
                     j1 = j1 + 1 ! index in moab_vert_coords for edges ! nverts_c + j1 for vertex edges !
                     ! current_2d_vertex%lat, current_2d_vertex%lon
                     pos_edge = .true.
                     iv = nverts_c + j1
                     edge_verts(j) = iv !  to form the local connectivity array
                     if ( vdone_c(edge(1, idx)) .gt.  vdone_c(edge(2, idx)) ) pos_edge = .false.
                     if (j .eq. 1) then
                         call gfr_f_get_corner_latlon(ie, 1, 1, 2, current_2d_vertex%lat, current_2d_vertex%lon)
                         if (pos_edge) then
                             vdone_pg (iv) = gdofel(ix + 2) ! elem(ie)%gdofP(2,1) ! gdofel(ix+ (j-1)*np + i)
                         else
                             vdone_pg (iv) = gdofel(ix + np - 1) !elem(ie)%gdofP(np-1,1) !
                         endif
                     else if (j .eq. 2) then
                         call gfr_f_get_corner_latlon(ie, 2, 1, 3, current_2d_vertex%lat, current_2d_vertex%lon)
                         if (pos_edge) then
                             vdone_pg (iv) = gdofel(ix + (2 - 1) * np + np)!elem(ie)%gdofP(np,2) ! ! gdofel(ix+ (j-1)*np + i)
                         else
                             vdone_pg (iv) = gdofel(ix + (np - 2) * np + np)!elem(ie)%gdofP(np,np - 1) !
                         endif
                     else if (j .eq. 3) then
                         call gfr_f_get_corner_latlon(ie, 2, 2, 4, current_2d_vertex%lat, current_2d_vertex%lon)
                         if (pos_edge) then
                             vdone_pg (iv) = gdofel(ix+ (np - 1) * np + np - 1)!elem(ie)%gdofP(np-1,np) !
                         else
                             vdone_pg (iv) = gdofel(ix+ (np-1)*np + 2) !elem(ie)%gdofP(2,np) !
                         endif
                     else ! if (j .eq. 4)
                         call gfr_f_get_corner_latlon(ie, 1, 2, 1, current_2d_vertex%lat, current_2d_vertex%lon)
                         if (pos_edge) then
                             vdone_pg (iv) = gdofel(ix+ (np - 2)*np + 1) !elem(ie)%gdofP(1,np-1) !
                         else
                             vdone_pg (iv) = gdofel(ix+ ( 2 - 1 )*np + 1) ! elem(ie)%gdofP(1,2) !
                         endif
                     endif
                     ! create the 3d vertex !
                     cart = spherical_to_cart (current_2d_vertex )
                     moab_vert_coords ( 3*(iv-1)+1 ) = cart%x
                     moab_vert_coords ( 3*(iv-1)+2 ) = cart%y
                     moab_vert_coords ( 3*(iv-1)+3 ) = cart%z
                     ! print *, 'ie, j, iv, vdone_pg(iv): ', ie, j, iv, vdone_pg(iv)
                 else ! the vertex was already created, but we need the index for connectivity of local fv cells
                     edge_verts(j) = nverts_c + ( -idx ) ! idx is index of edge (negative for already created)
                 endif

              enddo ! do j=1,4
              ! create the middle vertex too, in the center
              call gfr_f_get_corner_latlon(ie, 1, 1, 3, current_2d_vertex%lat, current_2d_vertex%lon)
              iv = nverts_c + edge_index + ie ! middle vertices are after corners, and edge vertices
              middle_vertex = iv
              vdone_pg (middle_vertex) = gdofel(ix+ np + 2)!elem(ie)%gdofP(2,2) ! first in the interior, not on edges!
              ! print *, 'ie, middle = iv, vdone_pg(iv): ', ie, iv, vdone_pg(iv)
              cart = spherical_to_cart (current_2d_vertex )
              moab_vert_coords ( 3*(iv-1)+1 ) = cart%x
              moab_vert_coords ( 3*(iv-1)+2 ) = cart%y
              moab_vert_coords ( 3*(iv-1)+3 ) = cart%z

              ! now form the local 2x2 cells, one by one; set the global id tag too!
              idx = (ie-1)*4
              ix = idx * 4 !
              ! first
              moabconn_pg(ix + 1) = moabconn_c(4*(ie-1)+1)
              moabconn_pg(ix + 2) = edge_verts(1)
              moabconn_pg(ix + 3) = middle_vertex
              moabconn_pg(ix + 4) = edge_verts(4)
              elemids(idx+1) = (elem(ie)%GlobalId-1)*4+1
              ! second
              moabconn_pg(ix + 4 + 1) = edge_verts(1)
              moabconn_pg(ix + 4 + 2) = moabconn_c(4*(ie-1)+2)
              moabconn_pg(ix + 4 + 3) = edge_verts(2)
              moabconn_pg(ix + 4 + 4) = middle_vertex
              elemids(idx+2) = (elem(ie)%GlobalId-1)*4+2
              ! third
              moabconn_pg(ix + 8 + 1) = edge_verts(4)
              moabconn_pg(ix + 8 + 2) = middle_vertex
              moabconn_pg(ix + 8 + 3) = edge_verts(3)
              moabconn_pg(ix + 8 + 4) = moabconn_c(4*(ie-1)+4)
              elemids(idx+3) = (elem(ie)%GlobalId-1)*4+3
              ! fourth
              moabconn_pg(ix + 12 + 1) = middle_vertex
              moabconn_pg(ix + 12 + 2) = edge_verts(2)
              moabconn_pg(ix + 12 + 3) = moabconn_c(4*(ie-1)+3)
              moabconn_pg(ix + 12 + 4) = edge_verts(3)
              elemids(idx+4) = (elem(ie)%GlobalId-1)*4+4

          enddo
          ! now copy from coarse for pg mesh

          dimcoord = nverts_pg*3
          dimen = 3
          if ( nverts_pg > 0 ) then
            ierr = iMOAB_CreateVertices(MHPGID, dimcoord, dimen, moab_vert_coords)
            if (ierr > 0 )  &
              call endrun('Error: fail to create MOAB vertices ')
          endif 
         ! num_el = nelem_pg  *
          mbtype = 3 !  quadrilateral
          nve = 4;
          block_ID = 300 ! this will be for pg mesh

          if ( nelem_pg > 0 ) then 
            ierr = iMOAB_CreateElements( MHPGID, nelem_pg, mbtype, nve, moabconn_pg, block_ID );
            if (ierr > 0 )  &
              call endrun('Error: fail to create MOAB elements')
          endif
          tagname='GLOBAL_ID'//C_NULL_CHAR
          tagtype = 0  ! dense, integer
          numco = 1
          ierr = iMOAB_DefineTagStorage(MHPGID, tagname, tagtype, numco,  tagindex )
          if (ierr > 0 )  &
             call endrun('Error: fail to retrieve GLOBAL id tag')

          ! now set the values
          ent_type = 0 ! vertex type
          if ( nverts_pg > 0 ) then 
            ierr = iMOAB_SetIntTagStorage ( MHPGID, tagname, nverts_pg , ent_type, vdone_pg)
            if (ierr > 0 )  &
              call endrun('Error: fail to set global id tag for vertices')
          endif
          ! set global id tag for pg2 elements, too; they will start at nets=1, end at nete=nelemd*4
          ent_type = 1 ! now set the global id tag on elements
          if ( nelem_pg > 0 ) then
            ierr = iMOAB_SetIntTagStorage ( MHPGID, tagname, nelem_pg , ent_type, elemids)
            if (ierr > 0 )  &
              call endrun('Error: fail to set global id tag for edges')
          endif 

          ierr = iMOAB_ResolveSharedEntities( MHPGID, nverts_pg, vdone_pg );
          if (ierr > 0 )  &
             call endrun('Error: fail to resolve shared ents for pg2 mesh')

          ierr = iMOAB_UpdateMeshInfo(MHPGID)
          if (ierr > 0 )  &
             call endrun('Error: fail to update mesh info for pg2 mesh')
#ifdef MOABDEBUG
    !     write out the mesh file to disk, in parallel
          outfile = 'wholeATM_PG2.h5m'//C_NULL_CHAR
          wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
          ierr = iMOAB_WriteMesh(MHPGID, outfile, wopts)
          if (ierr > 0 )  &
            call endrun('Error: fail to write the mesh file')
#endif
       endif  ! only valid for pg == 2
       if ( nelemd > 0 ) then 
         deallocate (local_cell_gids)
         deallocate (indx_cell)
         deallocate (edge) !
         deallocate(moabconn_pg) ! connectivity
         deallocate(vdone_pg)
       endif

     endif

     ! deallocate
     if ( nelemd > 0 ) then
       deallocate(moabvh)
       deallocate(moabconn) ! do not keep it anymore, we are not setting another tag on fine mesh
       deallocate(vdone)
       deallocate(gdofel)
       deallocate(indx)
       deallocate(elemids)
       deallocate(gdofv)
       deallocate(moabvh_c)
       deallocate(moabconn_c)
       deallocate(vdone_c)
     endif
!    end copy

  end subroutine create_moab_meshes
  
#endif
end module semoab_mod
