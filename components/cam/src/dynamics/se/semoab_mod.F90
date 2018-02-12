#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module semoab_mod    

  use kinds, only : real_kind, iulog
!  use edge_mod, only : ghostbuffertr_t, initghostbufferTR, freeghostbuffertr, &
!       ghostVpack, ghostVunpack,  edgebuffer_t, initEdgebuffer

  use dimensions_mod, only: nelem, ne, np
  use element_mod, only : element_t
  use parallel_mod, only : parallel_t

  use m_MergeSorts,     only: IndexSet, IndexSort
 
  use cam_grid_support, only:  iMap
  use cam_abortutils,   only : endrun

  use seq_comm_mct,  only: MHID, MHFID !  app id on moab side, for homme moab coarse and fine mesh

  implicit none

#include "moab/MOABConfig.h"
  
contains

  subroutine create_moab_mesh_coarse(par, elem, nets, nete)

    use ISO_C_BINDING

    type (element_t), intent(inout) :: elem(:)

    type (parallel_t)      , intent(in) :: par

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete

    integer ierr, i, j, ie, iv, block_ID, k, numvals

    real(kind=real_kind), allocatable, target :: moab_vert_coords(:)
    INTEGER (C_INT) ,allocatable , target :: moab_corner_quads(:)
    integer moab_dim_cquads, ix, idx ! used for indexing in loops; idx will have the number of local vertices

    integer nelemd  ! do not confuse this with dimensions_mod::nelemd

! do we really need this?
    integer , external :: iMOAB_CreateVertices, iMOAB_WriteMesh, iMOAB_CreateElements, &
        iMOAB_ResolveSharedEntities, iMOAB_UpdateMeshInfo, iMOAB_DefineTagStorage, &
        iMOAB_SetIntTagStorage

    integer(iMap), dimension(:), allocatable :: gdofv
    integer, dimension(:), allocatable :: indx  !  this will be ordered 
    
    !  this will be moab vertex handle locally
    integer, target, allocatable :: moabvh(:), moabconn(:), vdone(:), elemids(:), gdofel(:)
    integer  currentval, dimcoord, dimen, num_el, mbtype, nve

    character*100 outfile, wopts, localmeshfile, lnum, tagname
    integer  tagtype, numco, tag_sto_len, ent_type, tagindex


     nelemd = nete-nets+1
     moab_dim_cquads = (nete-nets+1)*4

     if(par%masterproc) then
       write (iulog, *) " MOAB: semoab_mod module: create_moab_mesh_coarse;  on processor " , par%rank ," elements: " ,  nets, nete
     endif
!
!     allocate(moab_fvm_coords(moab_dims_fvc))
     allocate(moab_corner_quads(moab_dim_cquads))
     allocate(elemids(nelemd))
     allocate(gdofel(nelemd*np*np))
     k=0 !   will be the index for element global dofs
     do ie=nets,nete
       ix = ie-nets
       moab_corner_quads(ix*4+1) = elem(ie)%gdofP(1,1)
       moab_corner_quads(ix*4+2) = elem(ie)%gdofP(np,1)
       moab_corner_quads(ix*4+3) = elem(ie)%gdofP(np,np)
       moab_corner_quads(ix*4+4) = elem(ie)%gdofP(1,np)
       elemids(ix+1) = elem(ie)%GlobalId
       do i=1,np
         do j=1,np
           k = k+1
           gdofel(k)=elem(ie)%gdofP(i,j)
         enddo
       enddo
     enddo
! now original order
     allocate(gdofv(moab_dim_cquads))
     do ix=1,moab_dim_cquads
        gdofv(ix) = moab_corner_quads(ix)
     enddo

!     order according to global dofs
     allocate(indx(moab_dim_cquads))
     call IndexSet(moab_dim_cquads, indx)
     call IndexSort(moab_dim_cquads, indx, gdofv, descend=.false.)
!      after sort, gdofv( indx(i)) < gdofv( indx(i+1) ) 

     allocate(moabvh(moab_dim_cquads))

     allocate(moabconn(moab_dim_cquads))
     idx=1
     currentval = gdofv( indx(1))
     do ix=1,moab_dim_cquads
        if (gdofv(indx(ix)) .ne. currentval ) then
          idx=idx+1
          currentval = gdofv(indx(ix))
        endif
        moabvh(ix) = idx  ! the vertex in connectivity array will be at this local index
        ! this will be the moab connectivity
        moabconn(indx(ix)) = idx
     enddo
     if(par%masterproc) then
       write (iulog, *) " MOAB: there are ", idx, " local vertices on master task"
     endif
     allocate(moab_vert_coords(3*idx) )
     allocate(vdone(idx))
     vdone = 0;
     currentval = gdofv( indx(1)) ! start over to identify coordinates of the vertices

     do ix=1,moab_dim_cquads
        i = indx(ix)   ! index in initial array
        ie = nets+ (i-1)/4 ! this is the element number
        j = i - ( i-1)/4*4 ! local index of vertex in element i
        iv = moabvh(ix)
        if (vdone(iv) .eq. 0) then
            moab_vert_coords ( 3*(iv-1)+1 ) = elem(ie)%corners3d(j)%x
            moab_vert_coords ( 3*(iv-1)+2 ) = elem(ie)%corners3d(j)%y
            moab_vert_coords ( 3*(iv-1)+3 ) = elem(ie)%corners3d(j)%z
            vdone(iv) = gdofv(indx(ix)) ! this will be now our tag used for resolving shared entities
        endif

     enddo

      dimcoord = idx*3
      dimen = 3
      ierr = iMOAB_CreateVertices(MHID, dimcoord, dimen, moab_vert_coords)
      if (ierr > 0 )  &
        call endrun('Error: fail to create MOAB vertices ')

      num_el = nete-nets+1
      mbtype = 3 !  quadrilateral
      nve = 4;
      block_ID = 100 ! this will be for coarse mesh

      ierr = iMOAB_CreateElements( MHID, num_el, mbtype, nve, moabconn, block_ID );
      if (ierr > 0 )  &
        call endrun('Error: fail to create MOAB elements')
      ! idx: num vertices; vdone will store now the markers used in global resolve
      ! for this particular problem, markers are the global dofs at corner nodes
! set the global id for vertices
!   first, retrieve the tag
      tagname='GLOBAL_ID'//CHAR(0)
      tagtype = 0  ! dense, integer
      numco = 1
      ierr = iMOAB_DefineTagStorage(MHID, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to retrieve global id tag')
      ! now set the values
      ent_type = 0 ! vertex type
      ierr = iMOAB_SetIntTagStorage ( MHID, tagname, idx , ent_type, vdone)
      if (ierr > 0 )  &
        call endrun('Error: fail to set global id tag for vertices')
      ! set global id tag for coarse elements, too; they will start at nets, end at nete
      ent_type = 1 ! now set the global id tag on elements
      ierr = iMOAB_SetIntTagStorage ( MHID, tagname, nelemd , ent_type, elemids)
      if (ierr > 0 )  &
        call endrun('Error: fail to set global id tag for vertices')

      ierr = iMOAB_ResolveSharedEntities( MHID, idx, vdone );
      if (ierr > 0 )  &
        call endrun('Error: fail to resolve shared entities')

!   global dofs are the GLL points are set for each element
      tagname='GLOBAL_DOFS'//CHAR(0)
      tagtype = 0  ! dense, integer
      numco = np*np !  usually, it is 16; each element will have the dofs in order
      ierr = iMOAB_DefineTagStorage(MHID, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create global DOFS tag')
      ! now set the values
      ! set global dofs tag for coarse elements, too; they will start at nets, end at nete
      ent_type = 1 ! now set the global id tag on elements
      numvals = nelemd*np*np ! input is the total number of values
      ierr = iMOAB_SetIntTagStorage ( MHID, tagname, numvals, ent_type, gdofel)
      if (ierr > 0 )  &
        call endrun('Error: fail to set global id tag for vertices')
! write in serial, on each task, before ghosting
      if (par%rank .lt. 5) then
        write(lnum,"(I0.2)")par%rank
        localmeshfile = 'owned_'//trim(lnum)// '.h5m' // CHAR(0)
        wopts = CHAR(0)
        ierr = iMOAB_WriteMesh(MHID, localmeshfile, wopts)
        if (ierr > 0 )  &
          call endrun('Error: fail to write local mesh file')
      endif

      ierr = iMOAB_UpdateMeshInfo(MHID)
      if (ierr > 0 )  &
        call endrun('Error: fail to update mesh info')
!     write out the mesh file to disk, in parallel
      outfile = 'wholeATM.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(MHID, outfile, wopts)
      if (ierr > 0 )  &
        call endrun('Error: fail to write the mesh file')

     ! deallocate
     deallocate(moab_corner_quads)
     deallocate(moabvh)
     deallocate(moabconn)
     deallocate(vdone)
     deallocate(gdofel)
     deallocate(indx)
     deallocate(elemids)

  end subroutine create_moab_mesh_coarse
  
  subroutine create_moab_mesh_fine(par, elem, nets, nete)

    use ISO_C_BINDING
    use coordinate_systems_mod, only :  cartesian3D_t,  spherical_to_cart
    type (element_t), intent(inout) :: elem(:)

    type (parallel_t)      , intent(in) :: par

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete

    integer ierr, i, j, ie, iv, block_ID, k, numvals
    integer icol, irow, je, linx ! local indices in fine el connect

    real(kind=real_kind), allocatable, target :: moab_vert_coords(:)
    INTEGER (C_INT) ,allocatable , target :: moab_corner_quads(:)
    integer moab_dim_cquads, ix, idx, nverts ! used for indexing in loops; nverts will have the number of local vertices

    integer nelemd  ! do not confuse this with dimensions_mod::nelemd

! do we really need this?
    integer , external :: iMOAB_CreateVertices, iMOAB_WriteMesh, iMOAB_CreateElements, &
        iMOAB_ResolveSharedEntities, iMOAB_UpdateMeshInfo, iMOAB_DefineTagStorage, &
        iMOAB_SetIntTagStorage, iMOAB_ReduceTagsMax

    integer(iMap), dimension(:), allocatable :: gdofv
    integer, dimension(:), allocatable :: indx  !  this will be ordered

    !  this will be moab vertex handle locally
    integer, target, allocatable :: moabvh(:), moabconn(:), vdone(:), elemids(:)
    integer  currentval, dimcoord, dimen, num_el, mbtype, nve

    character*100 outfile, wopts, localmeshfile, lnum, tagname, newtagg
    integer  tagtype, numco, tag_sto_len, ent_type, tagindex
    type (cartesian3D_t)             :: cart
    integer  igcol, ii
    integer local_map(np,np) !  what is the index of gll point (i,j) in a local moabconn(start: start+(np-1)*(np-1)*4-1)
    ! for np=4,
    !      28, 32, 36, 35
    !      25, 29, 33, 34
    ! j |  13, 17, 21, 22
    !      1,  5,  9,  10
    !(1,1)     i->

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

     nelemd = (nete-nets+1)*(np-1)*(np-1)
     moab_dim_cquads = (nete-nets+1)*4*(np-1)*(np-1)

     if(par%masterproc) then
       write (iulog, *) " MOAB: semoab_mod module: create_moab_mesh_fine;  on processor " , par%rank ," elements: " ,  nets, nete
     endif
!
!     allocate(moab_fvm_coords(moab_dims_fvc))
     allocate(moab_corner_quads(moab_dim_cquads))
     allocate(elemids(nelemd))

     k=0 !   will be the index for element global dofs
     do ie=nets,nete
       do j=1,np-1
         do i=1,np-1
           ix = (ie-nets)*(np-1)*(np-1)+(j-1)*(np-1)+i-1
           moab_corner_quads(ix*4+1) = elem(ie)%gdofP(i,j)
           moab_corner_quads(ix*4+2) = elem(ie)%gdofP(i+1,j)
           moab_corner_quads(ix*4+3) = elem(ie)%gdofP(i+1,j+1)
           moab_corner_quads(ix*4+4) = elem(ie)%gdofP(i,j+1)
           elemids(ix+1) = (elem(ie)%GlobalId-1)*(np-1)*(np-1)+(j-1)*(np-1)+i
         enddo
       enddo

     enddo
! now original order
     allocate(gdofv(moab_dim_cquads))
     do ix=1,moab_dim_cquads
        gdofv(ix) = moab_corner_quads(ix)
     enddo

!     order according to global dofs
     allocate(indx(moab_dim_cquads))
     call IndexSet(moab_dim_cquads, indx)
     call IndexSort(moab_dim_cquads, indx, gdofv, descend=.false.)
!      after sort, gdofv( indx(i)) < gdofv( indx(i+1) )
     allocate(moabvh(moab_dim_cquads))

     allocate(moabconn(moab_dim_cquads))
     idx=1
     currentval = gdofv( indx(1))
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
     allocate(moab_vert_coords(3*nverts) )
     allocate(vdone(nverts))
     vdone = 0;
     currentval = gdofv( indx(1)) ! start over to identify coordinates of the vertices

     do ix=1,moab_dim_cquads
        idx = indx(ix)   ! index in initial array, vertices in all fine quads
        k = (idx-1)/(4*(np-1)*(np-1))  ! index of coarse quad, locally, starts at 0
        ie = nets + k  ! this is the element number; starts at nets
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
            vdone(iv) = gdofv(indx(ix)) ! this will be now our tag used for resolving shared entities
        endif

     enddo

      dimcoord = nverts*3
      dimen = 3
      ierr = iMOAB_CreateVertices(MHFID, dimcoord, dimen, moab_vert_coords)
      if (ierr > 0 )  &
        call endrun('Error: fail to create MOAB vertices ')

      num_el = nelemd
      mbtype = 3 !  quadrilateral
      nve = 4;
      block_ID = 200 ! this will be for coarse mesh

      ierr = iMOAB_CreateElements( MHFID, num_el, mbtype, nve, moabconn, block_ID );
      if (ierr > 0 )  &
        call endrun('Error: fail to create MOAB elements')
      ! nverts: num vertices; vdone will store now the markers used in global resolve
      ! for this particular problem, markers are the global dofs at corner nodes
! set the global id for vertices
!   first, retrieve the tag
      tagname='GDOF'//CHAR(0)
      tagtype = 0  ! dense, integer
      numco = 1
      ierr = iMOAB_DefineTagStorage(MHFID, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to retrieve global id tag')
      ! now set the values
      ent_type = 0 ! vertex type
      ierr = iMOAB_SetIntTagStorage ( MHFID, tagname, nverts , ent_type, vdone)
      if (ierr > 0 )  &
        call endrun('Error: fail to set marker id tag for vertices')

      ierr = iMOAB_ResolveSharedEntities( MHFID, nverts, vdone );
      if (ierr > 0 )  &
        call endrun('Error: fail to resolve shared entities')

      vdone = -1 !  reuse vdone for the new tag, GDOF
! use element offset for actual global dofs
      ! tagtype = 0  ! dense, integer
      ! numco = 1
      newtagg='GLOBAL_ID'//CHAR(0)
      ierr = iMOAB_DefineTagStorage(MHFID, newtagg, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create new GDOF tag')
      do ie=nets,nete
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
      ierr = iMOAB_SetIntTagStorage ( MHFID, newtagg, nverts , ent_type, vdone)
      if (ierr > 0 )  &
        call endrun('Error: fail to set global dof tag for vertices')

      ierr = iMOAB_ReduceTagsMax ( MHFID, tagindex, ent_type)
      if (ierr > 0 )  &
        call endrun('Error: fail to reduce max tag')

      ! set global id tag for elements
      ent_type = 1 ! now set the global id tag on elements
      ierr = iMOAB_SetIntTagStorage ( MHFID, newtagg, nelemd , ent_type, elemids)
      if (ierr > 0 )  &
        call endrun('Error: fail to set global id tag for elements')

! write in serial, on each task, before ghosting
      if (par%rank .lt. 4) then
        write(lnum,"(I0.2)")par%rank
        localmeshfile = 'fineh_'//trim(lnum)// '.h5m' // CHAR(0)
        wopts = CHAR(0)
        ierr = iMOAB_WriteMesh(MHFID, localmeshfile, wopts)
        if (ierr > 0 )  &
          call endrun('Error: fail to write local mesh file')
      endif

      ierr = iMOAB_UpdateMeshInfo(MHFID)
      if (ierr > 0 )  &
        call endrun('Error: fail to update mesh info')
!     write out the mesh file to disk, in parallel
      outfile = 'wholeFineATM.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(MHFID, outfile, wopts)
      if (ierr > 0 )  &
        call endrun('Error: fail to write the mesh file')

     ! deallocate
     deallocate(moab_corner_quads)
     deallocate(moabvh)
     deallocate(moabconn)
     deallocate(vdone)
     deallocate(indx)
     deallocate(elemids)

  end subroutine create_moab_mesh_fine

end module semoab_mod


