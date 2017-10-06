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

  implicit none

#include "moab/MOABConfig.h"
  
  integer, public :: MHID !   app id on moab side, for homme moab coarse mesh
  
  

contains

  subroutine create_moab_mesh_coarse(par, elem, nets, nete)

    use ISO_C_BINDING

    type (element_t), intent(inout) :: elem(:)

    type (parallel_t)      , intent(in) :: par

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete

    integer ierr, i, j, ie, iv, block_ID, num_layers, dimgh, bridge

    real(kind=real_kind), allocatable, target :: moab_vert_coords(:)
    INTEGER (C_INT) ,allocatable , target :: moab_corner_quads(:)
    integer moab_dim_cquads, ix, idx ! used for indexing in loops; idx will have the number of local vertices

    integer nelemd

! do we really need this?
    integer , external :: iMOAB_CreateVertices, iMOAB_WriteMesh, iMOAB_CreateElements, &
        iMOAB_ResolveSharedEntities, iMOAB_DetermineGhostEntities

    integer(iMap), dimension(:), allocatable :: gdofv
    integer, dimension(:), allocatable :: indx  !  this will be ordered 
    
    !  this will be moab vertex handle locally
    integer, target, allocatable :: moabvh(:), moabconn(:), vdone(:)
    integer  currentval, dimcoord, dimen, num_el, mbtype, nve

    character*100 outfile, wopts, localmeshfile, lnum


     nelemd = nete-nets+1
     moab_dim_cquads = (nete-nets+1)*4

     if(par%masterproc) then
       write (iulog, *) " MOAB: semoab_mod module: create_moab_mesh_coarse;  on processor " , par%rank ," elements: " ,  nets, nete
     endif
!
!     allocate(moab_fvm_coords(moab_dims_fvc))
     allocate(moab_corner_quads(moab_dim_cquads))
     do ie=nets,nete
       ix = ie-nets
       moab_corner_quads(ix*4+1) = elem(ie)%gdofP(1,1)
       moab_corner_quads(ix*4+2) = elem(ie)%gdofP(np,1)
       moab_corner_quads(ix*4+3) = elem(ie)%gdofP(np,np)
       moab_corner_quads(ix*4+4) = elem(ie)%gdofP(1,np)
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
      ierr = iMOAB_ResolveSharedEntities( MHID, idx, vdone );
      if (ierr > 0 )  &
        call endrun('Error: fail to resolve shared entities')

! write in serial, on each task, before ghosting
      if (par%rank .lt. 10) then
        write(lnum,"(I0.2)")par%rank
        localmeshfile = 'owned_'//trim(lnum)// '.h5m' // CHAR(0)
        wopts = CHAR(0)
        ierr = iMOAB_WriteMesh(MHID, localmeshfile, wopts)
        if (ierr > 0 )  &
          call endrun('Error: fail to write local mesh file')
      endif

!     write out the mesh file to disk, in parallel
      outfile = 'whole.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(MHID, outfile, wopts)
      if (ierr > 0 )  &
        call endrun('Error: fail to write the mesh file')

      ! (iMOAB_AppID MHID, int * ghost_dim, int *num_ghost_layers, int * bridge_dim )
      dimgh = 2 ! will ghost quads, topological dim 2
      bridge = 0 ! use vertex as bridge
      num_layers = 1 ! so far, one layer only
      ierr = iMOAB_DetermineGhostEntities( MHID, dimgh, num_layers, bridge)
      if (ierr > 0 )  &
        call endrun('Error: fail to determine ghosts')

      ! write in serial, on each task
      if (par%rank .lt. 10) then
        write(lnum,"(I0.2)")par%rank
        localmeshfile = 'localmesh_'//trim(lnum)// '.h5m' // CHAR(0)
        wopts = CHAR(0)
        ierr = iMOAB_WriteMesh(MHID, localmeshfile, wopts)
        if (ierr > 0 )  &
          call endrun('Error: fail to write local mesh file')
      endif
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)

     ! deallocate
     deallocate(moab_corner_quads)
     deallocate(moabvh)
     deallocate(moabconn)
     deallocate(vdone)
     deallocate(indx)

  end subroutine create_moab_mesh_coarse
  

end module semoab_mod


