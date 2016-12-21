#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module moabmesh_mod    
!  use, intrinsic :: ISO_C_BINDING

  use kinds, only : real_kind, iulog
!  use edge_mod, only : ghostbuffertr_t, initghostbufferTR, freeghostbuffertr, &
!       ghostVpack, ghostVunpack,  edgebuffer_t, initEdgebuffer

  use dimensions_mod, only: nelem, ne, np
!  use time_mod, only : timelevel_t
  use element_mod, only : element_t
!  use fvm_control_volume_mod, only: fvm_struct
!  use fvm_mod, only: cellghostbuf, edgeveloc
!  use hybrid_mod, only : hybrid_t
  use parallel_mod, only : parallel_t

!  use control_mod, only: tracer_transport_type, TRACERTRANSPORT_MOAB , TRACERTRANSPORT_MOAB_LINEAR

  implicit none
! used to sort the dofs
    type groupsort_t
        integer :: order    ! original order of unsorted data
        integer :: value       ! values to be sorted by (dof global in this case)
    end type groupsort_t

#include "moab/MOABConfig.h"
!#include "moab/iMOAB.h"

!  public :: moabmesh

contains

!  qsort in fortran / reinvent the wheel
recursive subroutine MBQSort(a,na)

    ! DUMMY ARGUMENTS
    integer, intent(in) :: na
    type (groupsort_t), dimension(na), intent(in out) :: a

    ! LOCAL VARIABLES
    integer :: left, right
    real :: random
    integer :: pivot
    type (groupsort_t) :: temp
    integer :: marker

    if (nA > 1) then

        call random_number(random)
        pivot = A(int(random*real(nA-1))+1)%value   ! random pivot (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1

        do while (left < right)
            right = right - 1
            do while (A(right)%value > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left)%value < pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
            end if
        end do

        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if

        call MBQSort(A(:marker-1),marker-1)
        call MBQSort(A(marker:),nA-marker+1)

    end if

end subroutine MBQSort



  SUBROUTINE errorout(ierr, message)
  integer ierr
  character*(*) message
  if (ierr.ne.0) then
    print *, message
    call exit (1)
  end if
  return
  end subroutine

  subroutine create_moab_mesh(par, elem, nets, nete)

    use ISO_C_BINDING

    type (element_t), intent(inout) :: elem(:)

    type (parallel_t)      , intent(in) :: par

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete

    integer ierr, i, j, ie, pid, iv, block_ID, num_layers, dimgh, bridge

    real(kind=real_kind), allocatable, target :: moab_vert_coords(:)
    INTEGER (C_INT) ,allocatable , target :: moab_corner_quads(:)
    integer moab_dim_cquads, ix, idx ! used for indexing in loops; idx will have the number of local vertices

    integer nelemd
    character*12 appname
! do we really need this?
    integer , external :: iMOAB_InitializeFortran, iMOAB_RegisterFortranApplication, &
        iMOAB_CreateVertices, iMOAB_WriteMesh, iMOAB_CreateElements, &
        iMOAB_ResolveSharedEntities, iMOAB_DetermineGhostEntities

    type(groupsort_t), target, allocatable :: vertdof(:)
    !  this will be moab vertex handle locally
    integer, target, allocatable :: moabvh(:), moabconn(:), vdone(:)
    ! nve number of verts per element (4 for coarse quad)
    integer  currentval, dimcoord, dimen, num_el, mbtype, nve

    character*100 outfile, wopts, localmeshfile, lnum

 !  hybrid%par

     nelemd = nete-nets+1
 !    moab_dims_fvc = (nc+1)*(nc+1)*(nete-nets+1)*3
     moab_dim_cquads = (nete-nets+1)*4

!    write (iulog, *) " on processor " , par%rank ," elements: " ,  nets, nete
!     allocate(moab_fvm_coords(moab_dims_fvc))
     allocate(moab_corner_quads(moab_dim_cquads))
     do ie=nets,nete
       ix = ie-nets
!       do j=1,4
!         moab_corner_quads(ix*4+j) = elem(ie)%node_numbers(j)
!       enddo
!      write(iulog,*)( elem(ie)%node_numbers(j), j=1,4)
       moab_corner_quads(ix*4+1) = elem(ie)%gdofP(1,1)
       moab_corner_quads(ix*4+2) = elem(ie)%gdofP(np,1)
       moab_corner_quads(ix*4+3) = elem(ie)%gdofP(np,np)
       moab_corner_quads(ix*4+4) = elem(ie)%gdofP(1,np)
       ! some coordinates?
!       if (par%masterproc) then
!
!      write (iulog, *) "element: " , elem(ie)%vertex%number
!      do j=1,4
!        write (iulog, *) moab_corner_quads(ix*4+j), elem(ie)%corners3d(j)%x, &
!            elem(ie)%corners3d(j)%y, elem(ie)%corners3d(j)%z
!      enddo
!       end if
     enddo
! now original order
     allocate(vertdof(moab_dim_cquads))
     do ix=1,moab_dim_cquads
        vertdof(ix)%order=ix
        vertdof(ix)%value=moab_corner_quads(ix)
!        write (iulog, *) vertdof(ix)%order, vertdof(ix)%value
     enddo

     call MBQSort(vertdof,moab_dim_cquads)

!     write (iulog, *) 'after sorting'

     allocate(moabvh(moab_dim_cquads))

     allocate(moabconn(moab_dim_cquads))
     idx=1
     currentval = vertdof(1)%value
     do ix=1,moab_dim_cquads
        if (vertdof(ix)%value .ne. currentval ) then
          idx=idx+1
          currentval = vertdof(ix)%value
        endif
        moabvh(ix) = idx
        ! this will be the moab connectivity
        moabconn(vertdof(ix)%order) = idx
!        write (iulog, *) vertdof(ix)%order, vertdof(ix)%value, moabvh(ix)
     enddo
!     write (iulog, *) " there are ", idx, " local vertices "
     allocate(moab_vert_coords(3*idx) )
     allocate(vdone(idx))
     vdone = 0;
     currentval = vertdof(1)%value ! start over to identify coordinates of the vertices

     do ix=1,moab_dim_cquads
        i = vertdof(ix)%order  ! index in initial array
        ie = nets+ (i-1)/4 ! this is the element number
        j = i - ( i-1)/4*4
        iv = moabvh(ix)
        if (vdone(iv) .eq. 0) then
            moab_vert_coords ( 3*(iv-1)+1 ) = elem(ie)%corners3d(j)%x
            moab_vert_coords ( 3*(iv-1)+2 ) = elem(ie)%corners3d(j)%y
            moab_vert_coords ( 3*(iv-1)+3 ) = elem(ie)%corners3d(j)%z
            vdone(iv) = vertdof(ix)%value ! this will be now our tag used for resolving shared entities
!            write (iulog, *)iv,  moab_vert_coords ( 3*(iv-1)+1 ) , moab_vert_coords ( 3*(iv-1)+2 ), &
!                 moab_vert_coords ( 3*(iv-1)+3 )
        endif

     enddo

!     do ie=nets,nete
!       ix = ie-nets
!       write (iulog, *) ix+1,  moabconn(ix*4+1), moabconn(ix*4+2), moabconn(ix*4+3), moabconn(ix*4+4)
!     enddo

      ierr = iMOAB_InitializeFortran()
      call errorout(ierr, 'fail to initialize iMOAB')

      appname = 'HM_COARSE'//CHAR(0)

      ierr = iMOAB_RegisterFortranApplication(appname, par%comm, pid)
      call errorout(ierr, 'fail to register homme spectral mesh')

      dimcoord = idx*3
      dimen = 3
      ierr = iMOAB_CreateVertices(pid, dimcoord, dimen, moab_vert_coords)
      call errorout(ierr, 'fail to create vertices')

      num_el = nete-nets+1
      mbtype = 3 !  quadrilateral
      nve = 4;
      block_ID = 100 ! this will be for coarse mesh

      ierr = iMOAB_CreateElements( pid, num_el, mbtype, nve, moabconn, block_ID );
      call errorout(ierr, 'fail to create elements')
      ! idx: num vertices; vdone will store now the markers used in global resolve
      ! for this particular problem, markers are the global dofs at corner nodes
      ierr = iMOAB_ResolveSharedEntities( pid, idx, vdone );
      call errorout(ierr, 'fail to resolve shared entities')

! write in serial, on each task, before ghosting
      if (par%rank .lt. 10) then
        write(lnum,"(I0.2)")par%rank
        localmeshfile = 'owned_'//trim(lnum)// '.h5m' // CHAR(0)
        wopts = CHAR(0)
        ierr = iMOAB_WriteMesh(pid, localmeshfile, wopts)
        call errorout(ierr, 'fail to write local mesh file')
      endif

!     write out the mesh file to disk, in parallel
      outfile = 'whole.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(pid, outfile, wopts)
      call errorout(ierr, 'fail to write the mesh file')

      ! (iMOAB_AppID pid, int * ghost_dim, int *num_ghost_layers, int * bridge_dim )
      dimgh = 2 ! will ghost quads, topological dim 2
      bridge = 0 ! use vertex as bridge
      num_layers = 1 ! so far, one layer only
      ierr = iMOAB_DetermineGhostEntities( pid, dimgh, num_layers, bridge)
      call errorout(ierr, 'fail to determine ghosts')

      ! write in serial, on each task
      if (par%rank .lt. 10) then
        write(lnum,"(I0.2)")par%rank
        localmeshfile = 'localmesh_'//trim(lnum)// '.h5m' // CHAR(0)
        wopts = CHAR(0)
        ierr = iMOAB_WriteMesh(pid, localmeshfile, wopts)
        call errorout(ierr, 'fail to write local mesh file')
      endif
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)

     ! deallocate
     deallocate(moab_corner_quads)
     deallocate(vertdof)
     deallocate(moabvh)
     deallocate(moabconn)
     deallocate(vdone)
!    call exit(1)

  end subroutine create_moab_mesh
  

end module moabmesh_mod


