#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module semoab_mod    
#ifdef HAVE_MOAB
  use kinds, only : real_kind, iulog, long_kind, int_kind
!  use edge_mod, only : ghostbuffertr_t, initghostbufferTR, freeghostbuffertr, &
!       ghostVpack, ghostVunpack,  edgebuffer_t, initEdgebuffer

  use dimensions_mod, only: nelem, ne, np, nlev
  use element_mod, only : element_t
  use parallel_mod, only : parallel_t

  use m_MergeSorts,     only: IndexSet, IndexSort
 
  use cam_grid_support, only:  iMap
  use cam_abortutils,   only : endrun

  use seq_comm_mct,  only: MHID, MHFID !  app id on moab side, for homme moab coarse and fine mesh

  implicit none

  save

  integer local_map(np,np) !  what is the index of gll point (i,j) in a local moabconn(start: start+(np-1)*(np-1)*4-1)
  integer, allocatable :: moabconn(:) ! will have the connectivity in terms of local index in verts
  integer  ::                     num_calls_export
  
contains

  subroutine create_moab_mesh_fine(par, elem, nets, nete)

    use ISO_C_BINDING
    use coordinate_systems_mod, only :  cartesian3D_t,  spherical_to_cart
    type (element_t), intent(inout) :: elem(:)

    type (parallel_t)      , intent(in) :: par

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete

    integer ierr, i, j, ie, iv, block_ID, k, numvals
    integer icol, irow, je, linx ! local indices in fine el connect

    real(kind=real_kind), allocatable :: moab_vert_coords(:)

    integer moab_dim_cquads, ix, idx, nverts, nverts_c ! used for indexing in loops; nverts will have the number of local vertices

    integer nelemd2  ! do not confuse this with dimensions_mod::nelemd

! do we really need this?
    integer , external :: iMOAB_CreateVertices, iMOAB_WriteMesh, iMOAB_CreateElements, &
        iMOAB_ResolveSharedEntities, iMOAB_UpdateMeshInfo, iMOAB_DefineTagStorage, &
        iMOAB_SetIntTagStorage, iMOAB_ReduceTagsMax, iMOAB_GetIntTagStorage

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
    integer  igcol, ii

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

     nelemd2 = (nete-nets+1)*(np-1)*(np-1)
     moab_dim_cquads = (nete-nets+1)*4*(np-1)*(np-1)

     if(par%masterproc) then
       write (iulog, *) " MOAB: semoab_mod module: create_moab_mesh_fine;  on processor " , par%rank ," elements: " ,  nets, nete
     endif

     allocate(gdofv(moab_dim_cquads))
     allocate(elemids(nelemd2))

     k=0 !   will be the index for element global dofs
     do ie=nets,nete
       do j=1,np-1
         do i=1,np-1
           ix = (ie-nets)*(np-1)*(np-1)+(j-1)*(np-1)+i-1
           gdofv(ix*4+1) = elem(ie)%gdofP(i,j)
           gdofv(ix*4+2) = elem(ie)%gdofP(i+1,j)
           gdofv(ix*4+3) = elem(ie)%gdofP(i+1,j+1)
           gdofv(ix*4+4) = elem(ie)%gdofP(i,j+1)
           elemids(ix+1) = (elem(ie)%GlobalId-1)*(np-1)*(np-1)+(j-1)*(np-1)+i
         enddo
       enddo
     enddo

!     order according to global dofs

     allocate(moabvh(moab_dim_cquads))
     allocate(indx(moab_dim_cquads))

     allocate(moabconn(moab_dim_cquads))
     call IndexSet(moab_dim_cquads, indx)
     call IndexSort(moab_dim_cquads, indx, gdofv, descend=.false.)
!      after sort, gdofv( indx(i)) < gdofv( indx(i+1) )

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
            vdone(iv) = gdofv(indx(ix)) ! this will be now our tag used for resolving shared entities ! convert to int, from long int
        endif

     enddo

      dimcoord = nverts*3
      dimen = 3
      ierr = iMOAB_CreateVertices(MHFID, dimcoord, dimen, moab_vert_coords)
      if (ierr > 0 )  &
        call endrun('Error: fail to create MOAB vertices ')

      num_el = nelemd2
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

      vdone = -1 !  reuse vdone for the new tag, GLOBAL_ID (actual tag that we want to store global dof )
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
      ierr = iMOAB_SetIntTagStorage ( MHFID, newtagg, nelemd2 , ent_type, elemids)
      if (ierr > 0 )  &
        call endrun('Error: fail to set global id tag for elements')

!  now, after reduction, we can get the actual global ids for each vertex in the fine mesh
!  before, some vertices that were owned in MOAB but not owned in CAM did not have the right global ID tag
!  so vdone will be now correct on every task (no -1 anymore )
      ent_type = 0 ! vertex type
      allocate(vgids(nverts))
      ierr = iMOAB_GetIntTagStorage ( MHFID, newtagg, nverts , ent_type, vgids)
      if (ierr > 0 )  &
        call endrun('Error: fail to retrieve GLOBAL ID on each task')
#ifdef MOABDEBUG
! write in serial, on each task, before ghosting
      if (par%rank .lt. 4) then
        write(lnum,"(I0.2)")par%rank
        localmeshfile = 'fineh_'//trim(lnum)// '.h5m' // CHAR(0)
        wopts = CHAR(0)
        ierr = iMOAB_WriteMesh(MHFID, localmeshfile, wopts)
        if (ierr > 0 )  &
          call endrun('Error: fail to write local mesh file')
      endif
#endif
      ierr = iMOAB_UpdateMeshInfo(MHFID)
      if (ierr > 0 )  &
        call endrun('Error: fail to update mesh info')
#ifdef MOABDEBUG
!     write out the mesh file to disk, in parallel
      outfile = 'wholeFineATM.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(MHFID, outfile, wopts)
      if (ierr > 0 )  &
        call endrun('Error: fail to write the mesh file')
#endif

     ! deallocate
!    deallocate(moabvh)
!     deallocate(moabconn)
!     deallocate(vdone)
!     deallocate(indx)
!     deallocate(elemids)




!    now create the coarse mesh, but the global dofs will come from fine mesh, after solving
     nelemd2 = nete-nets+1
     moab_dim_cquads = (nete-nets+1)*4

     allocate(gdofel(nelemd2*np*np))
     k=0 !   will be the index for element global dofs
     do ie=nets,nete
       ix = ie-nets
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
     call IndexSet(moab_dim_cquads, indx)
     call IndexSort(moab_dim_cquads, indx, gdofv, descend=.false.)
!      after sort, gdofv( indx(i)) < gdofv( indx(i+1) )

     allocate(moabvh_c(moab_dim_cquads))

     allocate(moabconn_c(moab_dim_cquads))
     idx=1
     currentval = gdofv( indx(1))
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
     allocate(vdone_c(nverts_c))
     vdone_c = 0;
     currentval = gdofv( indx(1)) ! start over to identify coordinates of the vertices

     do ix=1,moab_dim_cquads
        i = indx(ix)   ! index in initial array
        ie = nets+ (i-1)/4 ! this is the element number
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
      ierr = iMOAB_CreateVertices(MHID, dimcoord, dimen, moab_vert_coords)
      if (ierr > 0 )  &
        call endrun('Error: fail to create MOAB vertices ')

      num_el = nete-nets+1
      mbtype = 3 !  quadrilateral
      nve = 4;
      block_ID = 100 ! this will be for coarse mesh

      ierr = iMOAB_CreateElements( MHID, num_el, mbtype, nve, moabconn_c, block_ID );
      if (ierr > 0 )  &
        call endrun('Error: fail to create MOAB elements')
      ! idx: num vertices; vdone will store now the markers used in global resolve
      ! for this particular problem, markers are the global dofs at corner nodes
! set the global id for vertices
!   first, retrieve the tag
      tagname='GDOFV'//CHAR(0)
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
      ierr = iMOAB_SetIntTagStorage ( MHID, tagname, nverts_c , ent_type, vdone_c)
      if (ierr > 0 )  &
        call endrun('Error: fail to set GDOFV tag for vertices')
      ! set global id tag for coarse elements, too; they will start at nets, end at nete
      ent_type = 1 ! now set the global id tag on elements
      ierr = iMOAB_SetIntTagStorage ( MHID, newtagg, nelemd2 , ent_type, elemids)
      if (ierr > 0 )  &
        call endrun('Error: fail to set global id tag for vertices')

      ierr = iMOAB_ResolveSharedEntities( MHID, idx, vdone_c );
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
      numvals = nelemd2*np*np ! input is the total number of values
      ! form gdofel from vgids
      do ie=1, nelemd2
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
      ierr = iMOAB_SetIntTagStorage ( MHID, tagname, numvals, ent_type, gdofel)
      if (ierr > 0 )  &
        call endrun('Error: fail to set globalDOFs tag for coarse elements')


      ! set the global ids for coarse vertices the same as corresponding fine vertices
      ent_type = 0 ! vertex type
      ierr = iMOAB_SetIntTagStorage ( MHID, newtagg, nverts_c , ent_type, vdone_c)

      ! create a new tag, for transfer example ; will use it now for temperature on the surface
      !  (bottom atm to surface of ocean)
      tagname='a2oTbot'//CHAR(0) !  atm to ocean temp bottom tag
      tagtype = 1  ! dense, double
      numco = np*np !  usually, it is 16; each element will have the same order as dofs
      ierr = iMOAB_DefineTagStorage(MHID, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create atm to ocean temp bottom tag')

      tagname='a2oUbot'//CHAR(0) !  atm to ocean U bottom tag
      ierr = iMOAB_DefineTagStorage(MHID, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create atm to ocean U velocity bottom tag')

      tagname='a2oVbot'//CHAR(0) !  atm to ocean V bottom tag
      ierr = iMOAB_DefineTagStorage(MHID, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create atm to ocean V velocity bottom tag')


      ! create a new tag, for transfer example ; will use it now for temperature on the surface
      !  (bottom atm to surface of ocean); for debugging, use it on fine mesh
      tagname='a2o_T'//CHAR(0) !  atm to ocean tag
      tagtype = 1  ! dense, double
      numco = 1 !  usually, it is 1; one value per gdof
      ierr = iMOAB_DefineTagStorage(MHFID, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create atm to ocean tag')

#ifdef MOABDEBUG
! write in serial, on each task, before ghosting
      if (par%rank .lt. 5) then
        write(lnum,"(I0.2)")par%rank
        localmeshfile = 'owned_'//trim(lnum)// '.h5m' // CHAR(0)
        wopts = CHAR(0)
        ierr = iMOAB_WriteMesh(MHID, localmeshfile, wopts)
        if (ierr > 0 )  &
          call endrun('Error: fail to write local mesh file')
      endif
#endif
      ierr = iMOAB_UpdateMeshInfo(MHID)
      if (ierr > 0 )  &
        call endrun('Error: fail to update mesh info')
#ifdef MOABDEBUG
!     write out the mesh file to disk, in parallel
      outfile = 'wholeATM.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(MHID, outfile, wopts)
      if (ierr > 0 )  &
        call endrun('Error: fail to write the mesh file')
#endif

     ! initialize
     num_calls_export = 0

     ! deallocate
     deallocate(moabvh)
!     deallocate(moabconn) keep it , it is useful to set the tag on fine mesh
     deallocate(vdone)
     deallocate(gdofel)
     deallocate(indx)
     deallocate(elemids)
     deallocate(gdofv)
     deallocate(moabvh_c)
     deallocate(moabconn_c)
     deallocate(vdone_c)
!    end copy

  end subroutine create_moab_mesh_fine

  subroutine moab_export_data(elem)

    type(element_t),    pointer :: elem(:)

    integer num_elem, ierr
    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
    integer, external :: iMOAB_GetMeshInfo, iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
    integer :: size_tag_array, nvalperelem, ie, i, j, je, ix, ent_type, idx

    real(kind=real_kind), allocatable :: valuesTag(:)
    character*100 outfile, wopts, tagname, lnum

    ! count number of calls
    num_calls_export = num_calls_export + 1

    ierr  = iMOAB_GetMeshInfo ( MHID, nvert, nvise, nbl, nsurf, nvisBC );
    ! find out the number of local elements in moab mesh
    num_elem = nvise(1)
    ! now print the temperature from the state, and set it
    nvalperelem = np*np
    size_tag_array = nvalperelem*num_elem
    !print *, 'num_elem  = ', num_elem
    !print *, ((local_map(i,j), i=1,np), j=1,np)
    !print *, (moabconn(i), i=1,np*np)
    ! now load the values on both tags
    allocate(valuesTag(size_tag_array))  ! will use the same array for vertex array

    do ie=1,num_elem
      do j=1,np
        do i=1,np
          valuesTag ( (ie-1)*np*np+(j-1)*np + i ) = elem(ie)%state%T(i,j,nlev,1) ! time level 1?
        enddo
      enddo
    enddo
    ! set the tag
    tagname='a2oTbot'//CHAR(0) !  atm to ocean tag for temperature
    ent_type = 1 ! element type
    ierr = iMOAB_SetDoubleTagStorage ( MHID, tagname, size_tag_array, ent_type, valuesTag)
    if (ierr > 0 )  &
      call endrun('Error: fail to set a2oTbot tag for coarse elements')

    ! loop now for U velocity ( a2oUbot tag)
    do ie=1,num_elem
      do j=1,np
        do i=1,np
          valuesTag ( (ie-1)*np*np+(j-1)*np + i ) = elem(ie)%state%v(i,j,1,nlev,1) ! time level 1, U comp
        enddo
      enddo
    enddo
    ! set the tag
    tagname='a2oUbot'//CHAR(0) !  atm to ocean tag for U velocity
    ent_type = 1 ! element type
    ierr = iMOAB_SetDoubleTagStorage ( MHID, tagname, size_tag_array, ent_type, valuesTag)
    if (ierr > 0 )  &
      call endrun('Error: fail to set a2oUbot tag for coarse elements')

    ! loop now for V velocity ( a2oVbot tag)
    do ie=1,num_elem
      do j=1,np
        do i=1,np
          valuesTag ( (ie-1)*np*np+(j-1)*np + i ) = elem(ie)%state%v(i,j,2,nlev,1) ! time level 1, V comp
        enddo
      enddo
    enddo
    ! set the tag
    tagname='a2oVbot'//CHAR(0) !  atm to ocean tag for V velocity
    ent_type = 1 ! element type
    ierr = iMOAB_SetDoubleTagStorage ( MHID, tagname, size_tag_array, ent_type, valuesTag)
    if (ierr > 0 )  &
      call endrun('Error: fail to set a2oVbot tag for coarse elements')


#ifdef MOABDEBUG
    !     write out the mesh file to disk, in parallel
    write(lnum,"(I0.2)")num_calls_export
    outfile = 'wholeATM_T_'//trim(lnum)// '.h5m' // CHAR(0)
    wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
    ierr = iMOAB_WriteMesh(MHID, outfile, wopts)
    if (ierr > 0 )  &
      call endrun('Error: fail to write the mesh file')
#endif

    ! for debugging, set the tag on fine mesh too (for visu)
    do ie=1,num_elem
      je = (ie-1)*(np-1)*(np-1)*4
      do j=1,np
        do i= 1,np
         ix = local_map(i,j)
         idx = moabconn( je + ix ) !
         valuesTag ( idx ) = elem(ie)%state%T(i,j,nlev,1)
        end do
      end do
    end do

    tagname='a2o_T'//CHAR(0) !  atm to ocean tag, on fine mesh
    ierr  = iMOAB_GetMeshInfo ( MHFID, nvert, nvise, nbl, nsurf, nvisBC );
    ent_type = 0 ! vertex type
    ierr = iMOAB_SetDoubleTagStorage ( MHFID, tagname, nvert(1), ent_type, valuesTag)
    if (ierr > 0 )  &
      call endrun('Error: fail to set a2o_T tag for fine vertices')

#ifdef MOABDEBUG
    !     write out the mesh file to disk, in parallel

    outfile = 'wholeFineATM_T_'//trim(lnum)// '.h5m' // CHAR(0)

    ierr = iMOAB_WriteMesh(MHFID, outfile, wopts)
    if (ierr > 0 )  &
      call endrun('Error: fail to write the fine mesh file, with a temperature on it')
#endif

    deallocate(valuesTag)
  end subroutine moab_export_data
#endif
end module semoab_mod
