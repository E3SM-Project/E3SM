module scalable_grid_init_mod
  use metagraph_mod, only: MetaVertex_t
  use parallel_mod, only: parallel_t, abortmp
  use dimensions_mod, only: npart
  use gridgraph_mod, only: GridVertex_t, GridEdge_t

  implicit none
  private

  public :: &

    ! (Almost) scalably generate MetaVertex. SFC generation is still unscalable,
    ! but ne must be > 10K before it becomes a concern. Also provide GridVertex
    ! and GridEdge. sgi_finalize will deallocate GridVertex, GridEdge, and
    ! MetaVertex, in addition to internal data.
    sgi_init_grid, &

    ! Deallocate data used in grid initialization. This is separate from
    ! sgi_init_grid because there are external clients of GridVertex and
    ! gid2igv. Both GridVertex and MetaVertex are deallocated.
    sgi_finalize, &

    ! Check each field (and subfields) of a MetaVertex against those of
    ! another. If the two are not exactly the same, output a message. This is
    ! used for correctness checking the output of sgi_init_grid.
    sgi_check, &

    ! Map a cell global ID to the entry in the GridVertex_t array GridVertex. An
    ! invalid gid input, i.e. one that does not exist in gm%gvid, is permitted;
    ! in that case, the output igv is invalid and so GridVertex(igv) /= gid.
    sgi_gid2igv

  type, public :: GridManager_t
     integer :: rank, phase, ne
     integer, allocatable :: sfcfacemesh(:,:), rank2sfc(:), gvid(:)
     logical(kind=1), allocatable :: owned_or_used(:)
     type (GridVertex_t), pointer :: gv(:) => null()
     type (GridEdge_t), pointer :: ge(:) => null()
     type (MetaVertex_t), pointer :: mv => null()
  end type GridManager_t

  type (GridManager_t), target :: sgi_gm

contains

  subroutine sgi_init_grid(par, GridVertex, GridEdge, MetaVertex)
    use dimensions_mod, only: nelem, ne
    use metagraph_mod, only: initMetaGraph

    type (parallel_t), intent(in) :: par
    type (GridVertex_t), pointer, intent(out) :: GridVertex(:)
    type (GridEdge_t), pointer, intent(out) :: GridEdge(:)
    type (MetaVertex_t), intent(out), target :: MetaVertex
    type (GridManager_t), pointer :: gm
    integer, allocatable :: sfctest(:)
    integer :: ie, i, j, face, id, sfc, nelemd, nelemdi, rank
    logical, parameter :: dbg = .false.

    gm => sgi_gm
    gm%mv => MetaVertex

    allocate(gm%rank2sfc(npart+1))
    call sgi_genspacepart(nelem, npart, gm%rank2sfc)
    allocate(gm%sfcfacemesh(ne,ne))
    call sgi_genspacecurve(ne, gm%sfcfacemesh)

    if (dbg) then
       ! sgi_genspacepart
       if (gm%rank2sfc(npart+1) /= nelem) then
          print *, 'SGI> nelem',nelem,'rank2sfc',gm%rank2sfc
       end if
       nelemd = gm%rank2sfc(2) - gm%rank2sfc(1)
       do i = 3, npart+1
          nelemdi = gm%rank2sfc(i) - gm%rank2sfc(i-1)
          if (nelemdi > nelemd) then
             print *, 'SGI> nelem',nelem,'rank2sfc',gm%rank2sfc
             exit
          end if
          nelemd = nelemdi
       end do
       ! sfc2rank
       do sfc = 0, nelem-1
          rank = sfc2rank(gm%rank2sfc, sfc)
          if (sfc < gm%rank2sfc(rank+1) .or. sfc >= gm%rank2sfc(rank+2)) then
             print *, 'SGI> sfc, rank, rank2sfc',sfc,rank,gm%rank2sfc
             exit
          end if
       end do
       ! u<->s and ->sfc
       allocate(sfctest(nelem))
       sfctest = 0
       do ie = 1, nelem
          call u2si(ie,i,j,face)
          id = s2ui(i,j,face)
          if (.not. (id == ie .and. &
               i >= 1 .and. i <= ne .and. &
               j >= 1 .and. j <= ne .and. &
               face >=1 .and. face <= 6)) then
             print *, 'SGI> u<->s:',ie,id,i,j,face
          end if
          sfc = u2sfc(gm%sfcfacemesh, id)
          if (.not. (sfc >= 0 .and. sfc < nelem)) then
             print *, 'SGI> u2sfc:',id,sfc
          end if
          sfctest(sfc+1) = sfctest(sfc+1) + 1
       end do
       do ie = 1, nelem
          if (sfctest(ie) .ne. 1) then
             print *, 'SGI> sfctest:',ie,sfctest(ie)
          end if
       end do
       deallocate(sfctest)
    end if

    gm%rank = par%rank
    call sgi_CubeTopology_phase1(gm)

    if (dbg) then
       do i = 1, size(gm%gvid)
          if (gid2igv(gm, gm%gvid(i)) /= i) then
             print *, 'SGI> igv, ret, gid, gvid', &
                  i,gid2igv(gm, gm%gvid(i)),gm%gvid(i),gm%gvid
             exit
          end if
       end do
    end if

    call sgi_CubeTopology_phase2(gm)
    call initMetaGraph(gm%rank + 1, MetaVertex, gm%gv, gm%ge)
    GridVertex => gm%gv
    GridEdge => gm%ge
  end subroutine sgi_init_grid

  subroutine sgi_finalize()
    use gridgraph_mod, only: deallocate_gridvertex_nbrs
    use metagraph_mod, only: destroyMetaGraph

    type (GridManager_t), pointer :: gm
    integer :: i

    gm => sgi_gm

    deallocate(gm%gvid)
    do i = 1, size(gm%gv)
       call deallocate_gridvertex_nbrs(gm%gv(i))
    end do
    deallocate(gm%gv)
    deallocate(gm%ge)
    call destroyMetaGraph(gm%mv)
  end subroutine sgi_finalize

  subroutine sgi_check(mv, mvo)
    use metagraph_mod, only: PrintMetaVertex, MetaEdge_t

    type (MetaVertex_t), intent(in) :: mv, mvo
    type (GridVertex_t), pointer :: gv, gvo
    type (MetaEdge_t), pointer :: me, meo
    type (GridManager_t), pointer :: gm
    integer :: i, npi, j

    gm => sgi_gm

    if (gm%rank == 0) print *, 'SGI> sgi_check'
    if (mv%number /= mvo%number) print *, 'SGI> number disagrees'
    if (mv%nmembers /= mvo%nmembers) print *, 'SGI> nmembers disagrees'
    if (mv%nedges /= mvo%nedges) print *, 'SGI> nedges disagrees'
    do i = 1, mv%nmembers
       gv => mv%members(i)
       gvo => mvo%members(i)
       ! This seems unused
       !if (gv%face_number /= gvo%face_number) print *, 'SGI> GV face_number disagrees'
       if (gv%number /= gvo%number) print *, 'SGI> GV number disagrees'
       if (gv%processor_number /= gvo%processor_number) &
            print *, 'SGI> GV processor_number disagrees'
       if (gv%SpaceCurve /= gvo%SpaceCurve) print *, 'SGI> GV SpaceCurve disagrees'
       do j = 1, size(gv%nbrs_ptr)
          if (gv%nbrs_ptr(j) /= gvo%nbrs_ptr(j)) print *, 'SGI> GV nbrs_ptr disagrees'
       end do
       do npi = 1, size(gv%nbrs_ptr)
          do j = gv%nbrs_ptr(npi), gv%nbrs_ptr(npi+1)-1
             if (gv%nbrs(j) /= gvo%nbrs(j)) print *, 'SGI> GV nbrs disagrees'
             if (gv%nbrs_face(j) /= gvo%nbrs_face(j)) print *, 'SGI> GV nbrs_face disagrees'
             if (gv%nbrs_wgt(j) /= gvo%nbrs_wgt(j)) print *, 'SGI> GV nbrs_wgt disagrees'
             if (gv%nbrs_wgt_ghost(j) /= gvo%nbrs_wgt_ghost(j)) &
                  print *, 'SGI> GV nbrs_wgt_ghost disagrees'
          end do
       end do
    end do
    do i = 1, mv%nedges
       me => mv%edges(i)
       meo => mvo%edges(i)
       if (me%number /= meo%number) print *, 'SGI> ME number disagrees'
       if (me%nmembers /= meo%nmembers) print *, 'SGI> ME nmembers disagrees'
       if (me%type /= meo%type) print *, 'SGI> ME type disagrees'
       if (me%wgtP /= meo%wgtP) print *, 'SGI> ME wgtP disagrees'
       if (me%wgtP_ghost /= meo%wgtP_ghost) print *, 'SGI> ME wgtP_ghost disagrees'
       if (me%wgtS /= meo%wgtS) print *, 'SGI> ME wgtS disagrees'
       if (me%HeadVertex /= meo%HeadVertex) print *, 'SGI> ME HeadVertex disagrees'
       if (me%TailVertex /= meo%TailVertex) print *, 'SGI> ME TailVertex disagrees'
       do j = 1, me%nmembers
          if (me%edgeptrP(j) /= meo%edgeptrP(j)) print *, 'SGI> ME edgeptrP disagrees'
          if (me%edgeptrS(j) /= meo%edgeptrS(j)) print *, 'SGI> ME edgeptrS disagrees'
          if (me%edgeptrP_ghost(j) /= meo%edgeptrP_ghost(j)) &
               print *, 'SGI> ME edgeptrP_ghost disagrees'
          if (me%members(j)%head_face /= meo%members(j)%head_face) &
               print *, 'SGI> ME GE head_face disagrees'
          if (me%members(j)%tail_face /= meo%members(j)%tail_face) &
               print *, 'SGI> ME GE tail_face disagrees'
          if (me%members(j)%head_dir /= meo%members(j)%head_dir) &
               print *, 'SGI> ME GE head_dir disagrees'
          if (me%members(j)%tail_dir /= meo%members(j)%tail_dir) &
               print *, 'SGI> ME GE tail_dir disagrees'
          if (me%members(j)%reverse .neqv. meo%members(j)%reverse) &
               print *, 'SGI> ME GE reverse disagrees'
          if (me%members(j)%head%number /= meo%members(j)%head%number) &
               print *, 'SGI> ME GE head disagrees'
          if (me%members(j)%tail%number /= meo%members(j)%tail%number) &
               print *, 'SGI> ME GE tail disagrees'
       end do
    end do
  end subroutine sgi_check

  ! Map structured (i,j,face) triple to global ID.
  function s2ui(i, j, face) result (id)
    use dimensions_mod, only: ne

    integer, intent(in) :: i, j, face
    integer :: id

    id = ((face-1)*ne + (j-1))*ne + i
  end function s2ui

  ! Map global ID to structured (i,j,face) triple.
  subroutine u2si(id, i, j, face)
    use dimensions_mod, only: ne

    integer, intent(in) :: id
    integer, intent(out) :: i, j, face
    integer :: nesq

    nesq = ne*ne
    face = (id-1)/nesq + 1
    j = modulo(id-1, nesq)/ne + 1
    i = modulo(id-1, ne) + 1
  end subroutine u2si

  ! Map structured (i,j,face) triple to SFC index.
  function s2sfc(sfcfacemesh, i, j, face) result(sfc)
    integer, intent(in) :: sfcfacemesh(:,:)
    integer, intent(in) :: i, j, face
    integer :: sfc, offset, ne, nesq

    ne = size(sfcfacemesh,1)
    nesq = ne*ne
    select case (face)
    case (1); sfc =          sfcfacemesh(i     , ne-j+1)
    case (2); sfc =   nesq + sfcfacemesh(i     , ne-j+1)
    case (3); sfc = 5*nesq + sfcfacemesh(i     , j     )
    case (4); sfc = 3*nesq + sfcfacemesh(ne-j+1, i     )
    case (5); sfc = 4*nesq + sfcfacemesh(i     , j     )
    case (6); sfc = 2*nesq + sfcfacemesh(ne-i+1, ne-j+1)
    end select
  end function s2sfc

  ! Map global element ID to SFC index.
  function u2sfc(sfcfacemesh, id) result(sfc)
    integer, intent(in) :: sfcfacemesh(:,:)
    integer, intent(in) :: id
    integer :: i, j, face, sfc

    call u2si(id, i, j, face)
    sfc = s2sfc(sfcfacemesh, i, j, face)
  end function u2sfc

  ! This routine to generate the space-filling curve (SFC) is not yet
  ! scalable. It allocates 4*ne^2 bytes.
  subroutine sgi_genspacecurve(ne, Mesh)
    use spacecurve_mod, only: IsFactorable, genspacecurve

    integer, intent(in) :: ne
    integer, intent(out) :: Mesh(ne,ne)
    integer, allocatable :: Mesh2(:,:), Mesh2_map(:,:,:), sfcij(:,:)
    integer :: i, i2, j, j2, k, ne2, sfc_index

    if(IsFactorable(ne)) then
       call GenspaceCurve(Mesh)
       !      call PrintCurve(Mesh) 
    else
       ! find the smallest ne2 which is a power of 2 and ne2>ne
       ne2=2**ceiling( log(real(ne))/log(2d0) )
       if (ne2<ne) call abortmp('Fatal SFC error')

       allocate(Mesh2(ne2,ne2))
       allocate(Mesh2_map(ne2,ne2,2))
       allocate(sfcij(0:ne2*ne2,2))

       call GenspaceCurve(Mesh2)  ! SFC partition for ne2

       ! associate every element on the ne x ne mesh (Mesh)
       ! with its closest element on the ne2 x ne2 mesh (Mesh2)
       ! Store this as a map from Mesh2 -> Mesh in Mesh2_map.
       ! elements in Mesh2 which are not mapped get assigned a value of 0
       Mesh2_map=0
       do j=1,ne
          do i=1,ne
             ! map this element to an (i2,j2) element
             ! [ (i-.5)/ne , (j-.5)/ne ]  = [ (i2-.5)/ne2 , (j2-.5)/ne2 ]
             i2=nint( ((i-.5)/ne)*ne2 + .5 )
             j2=nint( ((j-.5)/ne)*ne2 + .5 )
             if (i2<1) i2=1
             if (i2>ne2) i2=ne2
             if (j2<1) j2=1
             if (j2>ne2) j2=ne2
             Mesh2_map(i2,j2,1)=i
             Mesh2_map(i2,j2,2)=j
          enddo
       enddo

       ! create a reverse index array for Mesh2
       ! k = Mesh2(i,j) 
       ! (i,j) = (sfcij(k,1),sfci(k,2)) 
       do j=1,ne2
          do i=1,ne2
             k=Mesh2(i,j)
             sfcij(k,1)=i
             sfcij(k,2)=j
          enddo
       enddo

       ! generate a SFC for Mesh with the same ordering as the 
       ! elements in Mesh2 which map to Mesh.
       sfc_index=0
       do k=0,ne2*ne2-1
          i2=sfcij(k,1)
          j2=sfcij(k,2)
          i=Mesh2_map(i2,j2,1)
          j=Mesh2_map(i2,j2,2)
          if (i/=0) then
             ! (i2,j2) element maps to (i,j) element
             Mesh(i,j)=sfc_index
             sfc_index=sfc_index+1
          endif
       enddo
#if 0
       print *,'SFC Mapping to non powers of 2,3 used.  Mesh:'  
       do j=1,ne
          write(*,'(99i3)') (Mesh(i,j),i=1,ne)
       enddo
       call PrintCurve(Mesh2) 
#endif
       deallocate(Mesh2)
       deallocate(Mesh2_map)
       deallocate(sfcij)
    endif
  end subroutine sgi_genspacecurve

  ! Generate the SFC partitioning based on just nelem and npart. On output,
  ! rank2sfc(rank:rank+1) contains the SFC index inclusive-lower and
  ! exclusive-upper indices for part 'rank'.
  subroutine sgi_genspacepart(nelem, npart, rank2sfc)
    integer, intent(in) :: nelem, npart
    integer, intent(out) :: rank2sfc(npart+1)
    integer :: nelemd, ipart, extra, s1

    nelemd = nelem/npart
    ! every cpu gets nelemd elements, but the first 'extra' get nelemd+1
    extra = modulo(nelem,npart)
    s1 = extra*(nelemd+1)

    ! split curve into two curves:
    ! 1 ... s1  s2 ... nelem
    !
    !  s1 = extra*(nelemd+1)         (count be 0)
    !  s2 = s1+1 
    !
    ! First region gets nelemd+1 elements per Processor
    ! Second region gets nelemd elements per Processor

    rank2sfc(1) = 0
    do ipart = 1, extra
       rank2sfc(ipart+1) = ipart*(nelemd+1)
    end do
    do ipart = extra+1, npart
       rank2sfc(ipart+1) = s1 + (ipart - extra)*nelemd
    end do
  end subroutine sgi_genspacepart

  ! Map a SFC index to the rank that owns the corresponding element.
  function sfc2rank(rank2sfc, sfc) result(rank)
    integer, intent(in) :: rank2sfc(:), sfc
    integer :: npart, lo, hi, mid, rank
    
    npart = size(rank2sfc) - 1
    lo = 1
    hi = npart+1
    if (sfc >= rank2sfc(npart+1)) then
       print *, 'sfc2rank: rank2sfc, sfc',rank2sfc,sfc
       call abortmp('sfc2rank: sfc input is invalid')
    end if
    do while (hi > lo + 1)
       mid = (lo + hi)/2
       if (sfc >= rank2sfc(mid)) then
          lo = mid
       else
          hi = mid
       end if
    end do
    rank = lo - 1
  end function sfc2rank

  ! sgi_CubeTopology_phase1/2 are scalable replacements for
  ! cube_mod::CubeTopology. The GridVertex array (here, gm%gv) is constructed to
  ! have entries only for owned and remote-used elements.
  subroutine sgi_CubeTopology_phase1(gm)
    use gridgraph_mod, only: allocate_gridvertex_nbrs

    type (GridManager_t), intent(inout) :: gm
    integer, allocatable :: gvid(:)
    type (GridVertex_t), pointer :: gv
    integer :: id, ne, nelem, sfc, i, j, k, id_nbr, n_owned_or_used
    logical :: owned

    gm%phase = 1

    ne = size(gm%sfcfacemesh, 1)
    nelem = 6*ne*ne

    gm%ne = ne
    allocate(gm%owned_or_used(nelem))
    gm%owned_or_used = .false.

    do id = 1, nelem
       sfc = u2sfc(gm%sfcfacemesh, id)
       owned = sfc >= gm%rank2sfc(gm%rank+1) .and. sfc < gm%rank2sfc(gm%rank+2)

       if (.not. owned) cycle
       gm%owned_or_used(id) = .true.

       call sgi_CubeTopology_impl(gm, id, -1)
    end do

    n_owned_or_used = 0
    do id = 1, nelem
       if (gm%owned_or_used(id)) n_owned_or_used = n_owned_or_used + 1
    end do
    allocate(gm%gvid(n_owned_or_used))
    i = 1
    do id = 1, nelem
       if (gm%owned_or_used(id)) then
          gm%gvid(i) = id
          i = i + 1
       end if
    end do
    deallocate(gm%owned_or_used)
    allocate(gm%gv(n_owned_or_used))
    do i = 1, n_owned_or_used
       call allocate_gridvertex_nbrs(gm%gv(i))
       gv => gm%gv(i)
       gv%number = gm%gvid(i)
       gv%SpaceCurve = u2sfc(gm%sfcfacemesh, gv%number)
       gv%processor_number = sfc2rank(gm%rank2sfc, gv%SpaceCurve) + 1
       gv%face_number = 0
       gv%nbrs = -1
       gv%nbrs_wgt = 0
       gv%nbrs_ptr = 0
       gv%nbrs_wgt_ghost = 1
    end do
    deallocate(gm%sfcfacemesh)
    deallocate(gm%rank2sfc)
  end subroutine sgi_CubeTopology_phase1
  
  subroutine sgi_CubeTopology_phase2(gm)
    use gridgraph_mod, only: allocate_gridvertex_nbrs, deallocate_gridvertex_nbrs
    use cube_mod, only: CubeSetupEdgeIndex

    type (GridManager_t), intent(inout) :: gm
    integer :: igv, ngv, i, j, k, id_nbr, ll, loc
    type (GridVertex_t), pointer :: gv
    type (GridVertex_t) :: gv_tmp

    gm%phase = 2
    ngv = size(gm%gv)

    do igv = 1, ngv
       call sgi_CubeTopology_impl(gm, gm%gv(igv)%number, igv)
    end do

    call allocate_gridvertex_nbrs(gv_tmp)
    do igv = 1, ngv
       gv => gm%gv(igv)
       do i = 1, size(gv%nbrs)
          gv_tmp%nbrs(i) = gv%nbrs(i)
          gv_tmp%nbrs_face(i) = gv%nbrs_face(i)
          gv_tmp%nbrs_wgt(i) = gv%nbrs_wgt(i)
          gv_tmp%nbrs_wgt_ghost(i) = gv%nbrs_wgt_ghost(i)
       end do
       gv%nbrs_ptr(1) = 1
       do ll = 1,8
          loc = gv%nbrs_ptr(ll)
          if (gv_tmp%nbrs(ll) /= -1) then
             gv%nbrs(loc) = gv_tmp%nbrs(ll)
             gv%nbrs_face(loc) = gv_tmp%nbrs_face(ll)
             gv%nbrs_wgt(loc) = gv_tmp%nbrs_wgt(ll)
             gv%nbrs_wgt_ghost(loc) = gv_tmp%nbrs_wgt_ghost(ll)
             gv%nbrs_ptr(ll+1) = gv%nbrs_ptr(ll)+1
          else
             gv%nbrs_ptr(ll+1) = gv%nbrs_ptr(ll)
          end if
       end do       
    end do
    call deallocate_gridvertex_nbrs(gv_tmp)

    call sgi_initgridedge(gm)
    do i = 1, size(gm%ge)
       call CubeSetupEdgeIndex(gm%ge(i))
    end do
  end subroutine sgi_CubeTopology_phase2

  ! Set data in a single GridVertex according to gm%phase = 1 or 2.
  subroutine gm_set(gm, igv, i, j, face, dir)
    use dimensions_mod, only: np
    use control_mod, only: north, south, east, west, neast, seast, swest, nwest

    type (GridManager_t), intent(inout) :: gm
    integer, intent(in) :: igv, i, j, face, dir
    integer, parameter :: EdgeWgtP = np, CornerWgt = 1
    integer :: id, wgt

    id = s2ui(i,j,face)
    if (gm%phase == 1) then
       gm%owned_or_used(id) = .true.
    else
       gm%gv(igv)%nbrs(dir) = id
       gm%gv(igv)%nbrs_face(dir) = face
       if (dir == north .or. dir == south .or. dir == east .or. dir == west) then
          wgt = EdgeWgtP
       else
          wgt = CornerWgt
       end if
       gm%gv(igv)%nbrs_wgt(dir) = wgt
    end if
  end subroutine gm_set

  ! Called by both phases. gm_set does distinct things in each phase.
  subroutine sgi_CubeTopology_impl(gm, id, igv)
    use control_mod, only: north, south, east, west, neast, seast, swest, nwest

    type (GridManager_t), intent(inout) :: gm
    integer, intent(in) :: id, igv
    integer :: ne, i, j, k, rev

    ne = gm%ne

    call u2si(id, i, j, k)

    if (j >= 2 .and. i >= 2) then
       ! setup SOUTH, WEST, SW neighbors
       call gm_set(gm, igv, i-1, j, k, west)
       call gm_set(gm, igv, i, j-1, k, south)
       call gm_set(gm, igv, i-1, j-1, k, swest)
    end if
    if (j < ne .and. i < ne) then
       ! setup EAST, NORTH, NE neighbors
       call gm_set(gm, igv, i+1, j, k, east)
       call gm_set(gm, igv, i, j+1, k, north)
       call gm_set(gm, igv, i+1, j+1, k, neast)
    end if
    if (j > 1 .and. i < ne) then
       ! Setup the remaining SOUTH, EAST, and SE neighbors
       call gm_set(gm, igv, i, j-1, k, south)
       call gm_set(gm, igv, i+1, j, k, east)
       call gm_set(gm, igv, i+1, j-1, k, seast)
    end if
    if (j < ne .and. i > 1) then
       ! Setup the remaining NORTH, WEST, and NW neighbors
       call gm_set(gm, igv, i, j+1, k, north)
       call gm_set(gm, igv, i-1, j, k, west)
       call gm_set(gm, igv, i-1, j+1, k, nwest)
    end if
    if (k < 5) then ! west/east "belt" edges
       if (i == 1) then
          call gm_set(gm, igv, ne, j, modulo(2+k, 4)+1, west)
          if (j /= 1)  call gm_set(gm, igv, ne, j-1, modulo(2+k, 4)+1, swest)
          if (j /= ne) call gm_set(gm, igv, ne, j+1, modulo(2+k, 4)+1, nwest)
       else if (i == ne) then
          call gm_set(gm, igv, 1, j, modulo(k, 4)+1, east)
          if (j /= 1)  call gm_set(gm, igv, 1, j-1, modulo(k, 4)+1, seast)
          if (j /= ne) call gm_set(gm, igv, 1, j+1, modulo(k, 4)+1, neast)
       endif
    end if
    if (j == 1 .and. k == 1) then ! south edge of 1
       call gm_set(gm, igv, i, ne, 5, south)
       if (i /= 1)  call gm_set(gm, igv, i-1, ne, 5, swest)
       if (i /= ne) call gm_set(gm, igv, i+1, ne, 5, seast)
    end if
    if (j == ne .and. k == 5) then ! north edge of 5
       call gm_set(gm, igv, i, 1, 1, north)
       if (i /= 1)  call gm_set(gm, igv, i-1, 1, 1, nwest)
       if (i /= ne) call gm_set(gm, igv, i+1, 1, 1, neast)
    end if
    if (j == 1 .and. k == 2) then ! south edge of 2
       rev = ne+1-i
       call gm_set(gm, igv, ne, rev, 5, south)
       if (i /= 1)  call gm_set(gm, igv, ne, rev+1, 5, swest)
       if (i /= ne) call gm_set(gm, igv, ne, rev-1, 5, seast)
    end if
    if (i == ne .and. k == 5) then ! east edge of 5
       rev = ne+1-j
       call gm_set(gm, igv, rev, 1, 2, east)
       if (j /= 1)  call gm_set(gm, igv, rev+1, 1, 2, seast)
       if (j /= ne) call gm_set(gm, igv, rev-1, 1, 2, neast)
    end if
    if (j == 1 .and. k == 3) then ! south edge of 3
       rev = ne+1-i
       call gm_set(gm, igv, rev, 1, 5, south)
       if (i /= 1)  call gm_set(gm, igv, rev+1, 1, 5, swest)
       if (i /= ne) call gm_set(gm, igv, rev-1, 1, 5, seast)
    end if
    if (j == 1 .and. k == 5) then ! south edge of 5
       rev = ne+1-i
       call gm_set(gm, igv, rev, 1, 3, south)
       if (i /= 1)  call gm_set(gm, igv, rev+1, 1, 3, swest)
       if (i /= ne) call gm_set(gm, igv, rev-1, 1, 3, seast)
    end if
    if (j == 1 .and. k == 4) then ! south edge of 4
       call gm_set(gm, igv, 1, i, 5, south)
       if (i /= 1)  call gm_set(gm, igv, 1, i-1, 5, swest)
       if (i /= ne) call gm_set(gm, igv, 1, i+1, 5, seast)
    end if
    if (i == 1 .and. k == 5) then ! west edge of 5
       call gm_set(gm, igv, j, 1, 4, west)
       if (j /= 1)  call gm_set(gm, igv, j-1, 1, 4, swest)
       if (j /= ne) call gm_set(gm, igv, j+1, 1, 4, nwest)
    end if
    if (j == ne .and. k == 1) then ! north edge of 1
       call gm_set(gm, igv, i, 1, 6, north)
       if (i /= 1)  call gm_set(gm, igv, i-1, 1, 6, nwest)
       if (i /= ne) call gm_set(gm, igv, i+1, 1, 6, neast)
    end if
    if (j == 1 .and. k == 6) then ! south edge of 6
       call gm_set(gm, igv, i, ne, 1, south)
       if (i /= 1)  call gm_set(gm, igv, i-1, ne, 1, swest)
       if (i /= ne) call gm_set(gm, igv, i+1, ne, 1, seast)
    end if
    if (j == ne .and. k == 2) then ! north edge of 2
       call gm_set(gm, igv, ne, i, 6, north)
       if (i /= 1)  call gm_set(gm, igv, ne, i-1, 6, nwest)
       if (i /= ne) call gm_set(gm, igv, ne, i+1, 6, neast)
    end if
    if (i == ne .and. k == 6) then ! east edge of 6
       call gm_set(gm, igv, j, ne, 2, east)
       if (j /= 1)  call gm_set(gm, igv, j-1, ne, 2, seast)
       if (j /= ne) call gm_set(gm, igv, j+1, ne, 2, neast)
    end if
    if (j == ne .and. k == 3) then ! north edge of 3
       rev = ne+1-i
       call gm_set(gm, igv, rev, ne, 6, north)
       if (i /= 1)  call gm_set(gm, igv, rev+1, ne, 6, nwest)
       if (i /= ne) call gm_set(gm, igv, rev-1, ne, 6, neast)
    end if
    if (j == ne .and. k == 6) then ! north edge of 6
       rev = ne+1-i
       call gm_set(gm, igv, rev, ne, 3, north)
       if (i /= 1)  call gm_set(gm, igv, rev+1, ne, 3, nwest)
       if (i /= ne) call gm_set(gm, igv, rev-1, ne, 3, neast)
    end if
    if (j == ne .and. k == 4) then ! north edge of 4
       rev = ne+1-i
       call gm_set(gm, igv, 1, rev, 6, north)
       if (i /= 1)  call gm_set(gm, igv, 1, rev+1, 6, nwest)
       if (i /= ne) call gm_set(gm, igv, 1, rev-1, 6, neast)
    end if
    if (i == 1 .and. k == 6) then ! west edge of 6
       rev = ne+1-j
       call gm_set(gm, igv, rev, ne, 4, west)
       if (j /= 1)  call gm_set(gm, igv, rev+1, ne, 4, swest)
       if (j /= ne) call gm_set(gm, igv, rev-1, ne, 4, nwest)
    end if
  end subroutine sgi_CubeTopology_impl

  function gid2igv(gm, gid) result(igv)
    type (GridManager_t), intent(in) :: gm
    integer, intent(in) :: gid
    integer :: lo, hi, igv

    lo = 1
    hi = size(gm%gvid)
    do while (hi > lo)
       igv = (lo + hi)/2
       if (gm%gvid(igv) == gid) return
       if (gid > gm%gvid(igv)) then
          lo = igv+1
       else
          hi = igv-1
       end if
    end do
    igv = lo
  end function gid2igv

  ! Wrapper for external use.
  function sgi_gid2igv(gid) result (igv)
    integer, intent(in) :: gid
    integer :: igv

    igv = gid2igv(sgi_gm, gid)
  end function sgi_gid2igv

  ! Scalable replacement for GridGraph_mod::initgridedge. Only edges belonging
  ! to owned or remote-used elements are instantiated. Only owned and
  ! remote-used elements need to be represented in GridVertex.
  subroutine sgi_initgridedge(gm)
    use parallel_mod, only : abortmp
    use dimensions_mod, only : max_corner_elem
    use kinds, only: iulog
    use gridgraph_mod, only: num_neighbors

    type (GridManager_t), intent(inout) :: gm

    integer :: i,j,k,iptr,m,n,wgtV,wgtP,gid,igv
    integer :: nelem,nelem_edge,inbr,igvnbr
    integer :: mynbr_cnt, cnt, mystart, start
    logical :: owned_or_used
    logical, parameter :: Verbose = .false.

    nelem = size(gm%gv)

    ! Count the number of relevant edges.
    nelem_edge = 0
    do j = 1, nelem
       do i = 1, num_neighbors
          mynbr_cnt = gm%gv(j)%nbrs_ptr(i+1) - gm%gv(j)%nbrs_ptr(i) ! length of neighbor location
          mystart = gm%gv(j)%nbrs_ptr(i)
          do m=0, mynbr_cnt-1
             if (gm%gv(j)%nbrs_wgt(mystart + m) > 0) then ! want a non-0 weight
                owned_or_used = gm%gv(j)%processor_number == gm%rank + 1 ! want only owned or used elements
                if (.not. owned_or_used) then
                   gid = gm%gv(j)%nbrs(mystart + m)
                   igv = gid2igv(gm, gid)
                   owned_or_used = gm%gvid(igv) == gid
                end if
                if (owned_or_used) then
                   nelem_edge = nelem_edge + 1
                end if
             end if
          end do
       end do
    end do

    ! Allocate enough room for just these relevant edges.
    allocate(gm%ge(nelem_edge))
    gm%ge(:)%reverse = .false.

    ! Fill in the relevant edges.
    iptr = 1
    do j = 1, nelem
       do i = 1, num_neighbors    
          mynbr_cnt = gm%gv(j)%nbrs_ptr(i+1) - gm%gv(j)%nbrs_ptr(i)
          mystart = gm%gv(j)%nbrs_ptr(i)
          do m = 0, mynbr_cnt-1
             if (gm%gv(j)%nbrs_wgt(mystart + m) > 0) then
                owned_or_used = gm%gv(j)%processor_number == gm%rank + 1
                if (.not. owned_or_used) then
                   gid = gm%gv(j)%nbrs(mystart + m)
                   igv = gid2igv(gm, gid)
                   owned_or_used = gm%gvid(igv) == gid
                end if
                if (.not. owned_or_used) cycle

                gm%ge(iptr)%tail      => gm%gv(j)
                gm%ge(iptr)%tail_face =  mystart + m ! needs to be mystart + m (location in array)
                gm%ge(iptr)%tail_dir  =  i*max_corner_elem + m ! conversion needed for setcycle
                inbr                  =  gm%gv(j)%nbrs(mystart+m)
                igvnbr                =  gid2igv(gm, inbr)
                gm%ge(iptr)%head      => gm%gv(igvnbr)

                ! Determine which "face" of the neighbor element the edge links
                ! (i.e. the "head_face")
                do k = 1, num_neighbors
                   cnt = gm%gv(igvnbr)%nbrs_ptr(k+1) - gm%gv(igvnbr)%nbrs_ptr(k)  
                   start = gm%gv(igvnbr)%nbrs_ptr(k)
                   do n = 0, cnt-1
                      if (gm%gv(igvnbr)%nbrs(start+n) == gm%gv(j)%number) then
                         gm%ge(iptr)%head_face = start+n ! needs to be start + n (location in array)
                         gm%ge(iptr)%head_dir = k*max_corner_elem+n ! conversion (un-done in setcycle)
                      endif
                   enddo
                enddo
                iptr=iptr+1
             end if
          end do ! m loop
       end do ! i loop
    end do ! j loop
    if (nelem_edge+1 /= iptr) &
         call abortmp('Error in initgridedge: Number of edges less than expected.')

    if (Verbose) then
       print *
       write(iulog,*)"element edge tail,head list: (TEST)"
       do i=1,nelem_edge
          write(iulog,*)gm%ge(i)%tail%number,gm%ge(i)%head%number
       end do
       print *
       write(iulog,*)"element edge tail_face, head_face list: (TEST)"
       do i=1,nelem_edge
          write(iulog,*)gm%ge(i)%tail_face,gm%ge(i)%head_face
       end do
    end if
  end subroutine sgi_initgridedge
end module scalable_grid_init_mod
