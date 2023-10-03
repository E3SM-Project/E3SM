module mpas_moabmesh
!  use, intrinsic :: ISO_C_BINDING
#ifdef HAVE_MOAB
  use mpas_log     
  use mpas_derived_types, only: dm_info, domain_type
  use mpas_field_routines
  use mpas_sort
  use mpas_stream_manager
  use mpas_pool_routines

  !use mpas_vector_operations
#include "moab/MOABConfig.h"
  implicit none

  contains

  SUBROUTINE errorout(ierr, message)
  integer ierr
  character*(*) message
  if (ierr.ne.0) then
    print *, message
    call exit (1)
  end if
  return
  end subroutine
  
  subroutine init_moab_mpas(domain, ext_comp_id, pidmoab)
    use iMOAB, only : iMOAB_RegisterApplication, &
     iMOAB_CreateVertices, iMOAB_CreateElements, &
     iMOAB_ResolveSharedEntities, iMOAB_DetermineGhostEntities, &
     iMOAB_DefineTagStorage,  iMOAB_SetIntTagStorage , &
     iMOAB_UpdateMeshInfo, iMOAB_SetDoubleTagStorage
    use seq_flds_mod, only: seq_flds_dom_fields

    type (domain_type), intent(inout) :: domain
    integer , intent(in) :: ext_comp_id
    integer , Intent(inout) :: pidmoab

    type (block_type), pointer :: block
    type (mpas_pool_type), pointer :: meshPool
    integer, pointer :: nCells, nVertices, maxEdges
    integer  :: pid, nblocks
    integer, dimension(:,:), pointer :: verticesOnCell
    integer, dimension(:), pointer :: nEdgesOnCell
    integer, dimension(:), pointer :: indexToVertexID, indexToCellID
    real(kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, areaCell
    logical, pointer :: on_a_sphere, is_periodic
    real(kind=RKIND), pointer :: x_period, y_period
    integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve

    integer :: c_comm, i1, j1, ic, lastvertex,n
    character*12 appname
    integer :: ierr, num_verts_in_cells
    real(kind=RKIND), allocatable, target :: moab_vert_coords(:), data(:)
    integer, allocatable, target :: indexUsed(:), invMap(:), localIds(:)
    integer dimcoord, dimen, mbtype, block_ID, proc_id
    integer  ,allocatable , target :: all_connects(:)
    character*100  tagname, lnum
    integer tagtype, numco, tag_sto_len, ent_type, tagindex, currentVertex
    real (kind=RKIND), pointer :: sphere_radius

    c_comm = domain % dminfo % comm
    write(lnum,"(I0.2)")ext_comp_id
    appname = 'MPAS_MB_'//trim(lnum)//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, c_comm, ext_comp_id, pid)
    pidmoab = pid ! this is exported, need for send to work
    call errorout(ierr, 'fail to register MPAS_MOAB mesh')
    proc_id = domain % dminfo % my_proc_id
    call mpas_log_write('MOAB MPAS app pid:  $i  task  $i ', intArgs=(/pid, proc_id/) )

!   blocks should be merged if there is more than one block per task
    nblocks = 0
    block => domain % blocklist
    do while (associated(block)) !{{{
       nblocks = nblocks + 1
        ! allocate scratch memory 
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
       call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
       call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
       call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
       call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
       call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
       call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
       call mpas_pool_get_config(meshPool, 'x_period', x_period)
       call mpas_pool_get_config(meshPool, 'y_period', y_period)
       call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
       call mpas_pool_get_array(meshPool, 'indexToVertexID', indexToVertexID)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
       call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
       call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
!       call mpas_pool_get_array(meshPool, 'xCell', xCell)
!       call mpas_pool_get_array(meshPool, 'yCell', yCell)
!       call mpas_pool_get_array(meshPool, 'zCell', zCell)

       call mpas_log_write(' MOAB instance: number of vertices:: $i number of cells:: $i solve: v:$i c:$i', intArgs=(/nVertices, nCells, nVerticesSolve, nCellsSolve/) )
        !!
       allocate(indexUsed(nVertices), invMap(nVertices) ) ! conservative, invMap should be smaller
       indexUsed = 0
       invMap = 0
    ! fill now connectivity array, nCellsSolve; fist pad to max nc
       num_verts_in_cells = nCellsSolve * maxEdges
       allocate(all_connects(num_verts_in_cells))
!    collect all vertices, and also pad 
       j1 = 0
       do ic=1, nCellsSolve
         do i1 = 1, nEdgesOnCell(ic)
           j1 = j1 + 1
           all_connects(j1) = verticesOnCell( i1, ic)
           indexUsed(all_connects(j1)) = 1
         enddo
         lastvertex = verticesOnCell( nEdgesOnCell (ic), ic)
         ! pad the rest with the last vertex
         do i1 =   nEdgesOnCell (ic) + 1, maxEdges
           j1 = j1 + 1
           all_connects(j1) = lastvertex ! repeat the last vertex (pad)
         enddo
         ! call mpas_log_write('cell: $i  v:: $i $i $i $i $i $i', intArgs=(/ic, all_connects(j1-5), all_connects(j1-4), all_connects(j1-3), all_connects(j1-2), all_connects(j1-1), all_connects(j1)/) )
       enddo

       currentVertex = 0
       do i1 = 1, nVertices
         if (indexUsed(i1) > 0) then
           currentVertex = currentVertex + 1
           indexUsed(i1) = currentVertex
           invMap(currentVertex) = i1
         endif
       enddo
       !! convert all_connects to indexUsed
       do i1 = 1, num_verts_in_cells
         all_connects(i1) = indexUsed( all_connects(i1) )
       enddo
       allocate(moab_vert_coords(3*currentVertex))
       do i1 =1, currentVertex
         moab_vert_coords(3*i1-2) = xVertex(invMap(i1))/sphere_radius
         moab_vert_coords(3*i1-1) = yVertex(invMap(i1))/sphere_radius
         moab_vert_coords(3*i1  ) = zVertex(invMap(i1))/sphere_radius
         ! call mpas_log_write('i:: $i coords:: $r $r $r $r', intArgs=(/i1/), realArgs=(/moab_vert_coords(3*i1-2),moab_vert_coords(3*i1-1), moab_vert_coords(3*i1)/) )
       enddo
       dimcoord = 3*currentVertex
       dimen = 3
       ierr = iMOAB_CreateVertices(pid, dimcoord, dimen, moab_vert_coords)
       call errorout(ierr, 'fail to create vertices')
       call mpas_log_write(' MOAB instance: created $i vertices on local proc $i ',intArgs=(/currentVertex, proc_id/))
!  so we know we used only currentvertex vertices in the pool (the rest are in halo)
       mbtype = 4 !  polygon

       block_ID = 100*ext_comp_id + proc_id + nblocks ! we should have only one block right now

       ierr = iMOAB_CreateElements( pid, nCellsSolve, mbtype, maxEdges, all_connects, block_ID );
       call errorout(ierr, 'fail to create polygons')
 
!      set the global id for vertices
!   first, retrieve the tag
       tagname='GLOBAL_ID'//C_NULL_CHAR
       tagtype = 0  ! dense, integer
       numco = 1
       ierr = iMOAB_DefineTagStorage(pid, tagname, tagtype, numco,  tagindex )
       call errorout(ierr, 'fail to get global id tag')
! now set the values
       ent_type = 0 ! vertex type
       allocate(localIds(currentVertex))
       do i1 = 1, currentVertex
         localIds(i1) = indexToVertexID (invMap(i1))
       enddo
       ierr = iMOAB_SetIntTagStorage ( pid, tagname, currentVertex , ent_type, localIds )
       call errorout(ierr, 'fail to set global id tag for vertices')

       ! set global id tag for elements
       ent_type = 1 ! now set the global id tag on elements
       ierr = iMOAB_SetIntTagStorage ( pid, tagname, nCellsSolve, ent_type, indexToCellID)
       call errorout(ierr, 'fail to set global id tag for polygons')
       ! get next block
       ierr = iMOAB_ResolveSharedEntities( pid, currentVertex, localIds );
       call errorout(ierr, 'fail to resolve shared entities')

       allocate(data(nCellsSolve))
       n=0
       do ic=1, nCellsSolve
         n= n+1
         data(n)=areaCell(ic) / (sphere_radius * sphere_radius)
       enddo

       ! add domain tags
       tagname=trim(seq_flds_dom_fields)//C_NULL_CHAR
       tagtype = 1
       ierr = iMOAB_DefineTagStorage(pid, tagname, tagtype, numco,  tagindex )
       if (ierr > 0 )  &
          call errorout(ierr, 'Error: fail to define domain tags ')

       tagname='area'//C_NULL_CHAR
       ierr = iMOAB_SetDoubleTagStorage ( pid, tagname, nCellsSolve , ent_type, data)
       if (ierr > 0 )  &
          call errorout(ierr,'Error: fail to set area tag ')

       ! set it to the same area, but change it eventually with the area from mapping routing (tempestremap for MOAB)
       ! comment it 
       ! tagname='aream'//C_NULL_CHAR
       ! ierr = iMOAB_SetDoubleTagStorage ( pid, tagname, nCellsSolve , ent_type, data)
       ! if (ierr > 0 )  &
       !    call errorout(ierr,'Error: fail to set aream tag ')

       deallocate (moab_vert_coords)
       deallocate (data)
       deallocate (all_connects)
       deallocate (indexUsed)
       deallocate (invMap)
       deallocate (localIds)
       block => block % next

    end do !}}}
    

    if (nblocks .ne. 1) then
      call errorout(1, 'more than one block per task')
    endif 
    ierr = iMOAB_UpdateMeshInfo(pid)
    call errorout(ierr, 'fail to update mesh info')

    
    
  end subroutine init_moab_mpas
#endif
end module mpas_moabmesh
